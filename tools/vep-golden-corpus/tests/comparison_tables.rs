// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tool's copies of what the format comparator states in two files of its own, read by
//! their paths in `tools/vep-formats-compare/`: the class names of `data/classes.tsv` and
//! the cases of `tests/data/builder_parity/cases.tsv`.
//!
//! The class file is one line a class under a header line: a cell of a table's column is
//! empty where the class is not in the table and else its place there, counted from 0. Its
//! column `CORRECTED_IN` is no table: a cell there is the first Ensembl VEP release that
//! corrects the defect the class sets aside. The
//! case table is one line a case of three tab-separated cells: the item, its arguments as a
//! JSON list and its answer as a JSON value, a tuple a list, a set a sorted list and an
//! absent answer `null`.

mod common;

use common::RECODER_CLASSES;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use vep_golden_corpus::classify::hgvs::notation::{exonic_run, intron_offset_text};
use vep_golden_corpus::classify::introns as ie;
use vep_golden_corpus::classify::{
    allele_number, downstream, fields, hgvs, loftee, organelle, protein, sequence_edges, spliceai,
};
use vep_golden_corpus::constants::CORRECTED_BEHAVIOUR_CLASSES;
use vep_golden_corpus::json::{loads, Object, Value};

/// A file of the format comparator, by its path inside that tool.
fn comparator_file(name: &str) -> String {
    let path = Path::new(common::ROOT)
        .with_file_name("vep-formats-compare")
        .join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The class file: its header and its rows, each of as many cells as the header and no two
/// of one class.
struct ClassFile {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl ClassFile {
    fn read() -> ClassFile {
        let text = comparator_file("data/classes.tsv");
        let mut lines = text
            .lines()
            .map(|line| line.split('\t').map(str::to_string).collect::<Vec<_>>());
        let header = lines.next().expect("a header line");
        let rows: Vec<Vec<String>> = lines.collect();
        assert_eq!(header[0], "class");
        for row in &rows {
            assert_eq!(row.len(), header.len(), "{}", row[0]);
        }
        let classes: BTreeSet<&str> = rows.iter().map(|row| row[0].as_str()).collect();
        assert_eq!(classes.len(), rows.len(), "a class has one row");
        ClassFile { header, rows }
    }

    /// The classes of the column `table` in the order of their places, which count from 0
    /// without a gap.
    fn column(&self, table: &str) -> Vec<&str> {
        let at = self.header.iter().position(|name| name == table);
        let at = at.unwrap_or_else(|| panic!("the class file has no column {table}"));
        let mut placed: Vec<(usize, &str)> = self
            .rows
            .iter()
            .filter(|row| !row[at].is_empty())
            .map(|row| (row[at].parse().expect("a place"), row[0].as_str()))
            .collect();
        placed.sort_unstable();
        let places: Vec<usize> = placed.iter().map(|(place, _)| *place).collect();
        assert_eq!(places, (0..placed.len()).collect::<Vec<_>>(), "{table}");
        placed.into_iter().map(|(_, class)| class).collect()
    }

    /// Whether the file has a row of `class`, and whether that row has a place in a table.
    fn row(&self, class: &str) -> Option<bool> {
        let row = self.rows.iter().find(|row| row[0] == class)?;
        let places = self.header.iter().zip(row).skip(1);
        let mut tables =
            places.filter(|(name, _)| !name.contains('.') && name.as_str() != CORRECTED_IN);
        Some(tables.any(|(_, place)| !place.is_empty()))
    }

    /// The classes with a cell in the column `CORRECTED_IN`, each with the release the cell
    /// names as its two numbers, in the file's order.
    fn corrected_in(&self) -> Vec<(&str, (u32, u32))> {
        let at = self.header.iter().position(|name| name == CORRECTED_IN);
        let at = at.unwrap_or_else(|| panic!("the class file has no column {CORRECTED_IN}"));
        self.rows
            .iter()
            .filter(|row| !row[at].is_empty())
            .map(|row| {
                let (release, patch) = row[at].split_once('.').expect("a release and a patch");
                let number = |text: &str| text.parse().expect("a whole number");
                (row[0].as_str(), (number(release), number(patch)))
            })
            .collect()
    }
}

/// The column of the class file that names the first Ensembl VEP release that corrects a
/// class's defect.
const CORRECTED_IN: &str = "CORRECTED_IN";

/// The Ensembl VEP release the tool's corpora are written against.
const REFERENCE_RELEASE: (u32, u32) = (116, 3);

/// The classes the functions of the field dispatch (`fields::FIELD_CLASSES`) name.
const DISPATCH_CLASSES: [&str; 27] = [
    fields::MINIMAL_PLUGIN_SEES_THE_REJOINED_RECORD,
    fields::MOTIF_SCORE_REVERSE_STRAND_MNV,
    fields::CELL_TYPE_COMMA_NAMES,
    hgvs::REFSEQ_EDITED_MODEL_HGVSC_OFFSET,
    hgvs::HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT,
    hgvs::REFSEQ_EDITED_MODEL_HGVS_OFFSET,
    hgvs::SHIFTED_ALLELE_ROTATION_BY_STRAND,
    organelle::ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE,
    hgvs::HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE,
    hgvs::HGVS_SHIFTED_DELETION_CLAMPED_SPAN,
    hgvs::HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST,
    hgvs::HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE,
    hgvs::HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE,
    protein::HGVSP_STOP_CODON_PROTEIN_COMPARISON,
    protein::HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE,
    downstream::PLUGIN_HGVS_SHIFTED_CDS_COORDINATES,
    downstream::DOWNSTREAM_UTR_FRAME_PAST_CODING_END,
    downstream::DOWNSTREAM_START_PAST_OWN_CODON,
    ie::HGVSP_INTRON_END_INSERTION,
    spliceai::SPLICEAI_MULTI_ENTRY_CELL,
    ie::HGVSP_INTRON_START_INSERTION,
    ie::MAXENTSCAN_OVERLAPPING_WINDOWS,
    protein::HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE,
    allele_number::ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT,
    loftee::LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE,
    protein::HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE,
    sequence_edges::HGVS_SHIFT_SEQUENCE_EDGE_WINDOW,
];

#[test]
fn the_class_lists_are_the_columns_of_the_format_comparators_class_file() {
    let file = ClassFile::read();
    // The corrected behaviours: the rows with a place in one of four tables, but for the
    // classes whose defect Ensembl VEP corrects at the tool's reference release or before.
    let in_a_table: BTreeSet<&str> = [
        "REGISTRY_CLASSES",
        "FIELD_CLASSES",
        "PLUGIN_WARNING_CLASSES",
        "MANIFEST_CLASSES",
    ]
    .iter()
    .flat_map(|table| file.column(table))
    .collect();
    let corrected_by_ensembl_vep: BTreeSet<&str> = file
        .corrected_in()
        .into_iter()
        .filter(|(_, release)| *release <= REFERENCE_RELEASE)
        .map(|(class, _)| class)
        .collect();
    assert_eq!(
        corrected_by_ensembl_vep.iter().copied().collect::<Vec<_>>(),
        [
            "frequency_filter_allele_key_order",
            "frequency_filter_zero_value"
        ]
    );
    assert!(corrected_by_ensembl_vep.is_subset(&in_a_table));
    let corrected: BTreeSet<&str> = in_a_table
        .difference(&corrected_by_ensembl_vep)
        .copied()
        .collect();
    let own: BTreeSet<&str> = CORRECTED_BEHAVIOUR_CLASSES.into_iter().collect();
    assert_eq!(own.len(), CORRECTED_BEHAVIOUR_CLASSES.len());
    assert_eq!(own, corrected);
    // The classes the `masks` of a recoder corpus may carry, in the column's order.
    assert_eq!(file.column("RECODER_CORPUS_MASKS"), RECODER_CLASSES);
    // The field dispatch is functions in the order `classify` asks them, which is no
    // column's order. The classes its functions name are field classes of the comparison,
    // but for four the comparison's table lacks: three are no row of the file, and one is
    // a row without a place in any table.
    let field_classes = file.column("FIELD_CLASSES");
    let outside: Vec<&str> = DISPATCH_CLASSES
        .into_iter()
        .filter(|class| !field_classes.contains(class))
        .collect();
    assert_eq!(
        outside,
        [
            fields::MINIMAL_PLUGIN_SEES_THE_REJOINED_RECORD,
            fields::MOTIF_SCORE_REVERSE_STRAND_MNV,
            fields::CELL_TYPE_COMMA_NAMES,
            hgvs::SHIFTED_ALLELE_ROTATION_BY_STRAND,
        ]
    );
    let rows: Vec<Option<bool>> = outside.iter().map(|class| file.row(class)).collect();
    assert_eq!(rows, [None, Some(false), None, None]);
    assert_eq!(
        DISPATCH_CLASSES.into_iter().collect::<BTreeSet<_>>().len(),
        DISPATCH_CLASSES.len()
    );
    // The two classes asked of the format comparator are field classes there, and the
    // member an entry names its second class under is a row that is no class of a table.
    assert!(field_classes.contains(&spliceai::SPLICEAI_MULTI_ENTRY_CELL));
    assert!(field_classes.contains(&sequence_edges::HGVS_SHIFT_SEQUENCE_EDGE_WINDOW));
    assert_eq!(file.row(fields::SECOND_EXPECTED_DIVERGENCE), Some(false));
}

/// The items of the case table this tool keeps a copy of, and the cases of each.
const HELD: [(&str, usize); 24] = [
    ("CLIPPED_NOTATION_RE.pattern", 1),
    ("coordinate_intron_offset", 14),
    ("exonic_run", 196),
    ("SpliceWindows.holds", 12),
    ("INTRON_END_COPY_RE.pattern", 1),
    ("INTRON_END_DUPLICATION_RE.pattern", 1),
    ("INTRON_END_INSERTION_RE.pattern", 1),
    ("INTRON_END_PAIRS", 1),
    ("SHORT_INTRON_BASES", 1),
    ("transcript_intron_end_gaps", 18),
    ("intron_end_insertion_form", 13),
    ("copy_of_the_intron_last_base", 13),
    ("duplication_of_the_intron_last_bases", 13),
    ("insertion_that_cannot_end_the_intron", 13),
    ("insertion_that_ends_in_an_intron_end_pair", 13),
    ("INTRON_START_INSERTION_RE.pattern", 1),
    ("coding_transcript_exons", 18),
    ("transcript_intron_start_gaps", 18),
    ("insertion_after_the_exon_last_base", 10),
    ("held_insertion_gaps", 10),
    ("reverse_complement", 4),
    ("notation_inserted_bases", 8),
    ("lone_insertion_gaps", 15),
    ("notation_duplicated_positions", 7),
];

/// The items of the case table this tool keeps no copy of. `CLIPPED_NOTATION_RE.flags` is
/// the interpreter's number for the flags of a compiled pattern: the tool holds the
/// pattern's text and reads its `\d` as the ten ASCII digits, which no flag of that number
/// states.
const NOT_HELD: [(&str, usize); 1] = [("CLIPPED_NOTATION_RE.flags", 1)];

fn text(value: &Value) -> &str {
    value.as_str().expect("a text")
}

fn whole(value: &Value) -> i128 {
    i128::from(value.as_i64().expect("a whole number"))
}

fn object(value: &Value) -> &Object {
    value.as_object().expect("an object")
}

fn number(position: i128) -> Value {
    Value::from(i64::try_from(position).expect("a number of 64 bits"))
}

fn numbers(positions: impl IntoIterator<Item = i128>) -> Value {
    Value::from(positions.into_iter().map(number).collect::<Vec<_>>())
}

fn optional(text: Option<impl Into<String>>) -> Value {
    text.map_or(Value::Null, |text| Value::from(text.into()))
}

/// The positions of `first` to `last` that lie in two windows of a transcript on `strand`
/// whose introns are `introns`, `[first, last]` pairs in the transcript's order: the
/// transcript's exons are the bases between its introns and a hundred bases either side of
/// them.
fn held_positions(strand: i128, introns: &Value, first: i128, last: i128) -> Value {
    let mut introns: Vec<(i128, i128)> = introns
        .as_list()
        .expect("the introns")
        .iter()
        .map(|intron| {
            let ends = intron.as_list().expect("an intron");
            (whole(&ends[0]), whole(&ends[1]))
        })
        .collect();
    introns.sort_unstable();
    let mut starts = vec![introns[0].0 - 100];
    starts.extend(introns.iter().map(|intron| intron.1 + 1));
    let mut ends: Vec<i128> = introns.iter().map(|intron| intron.0 - 1).collect();
    ends.push(introns[introns.len() - 1].1 + 100);
    let exons: Vec<String> = starts
        .iter()
        .zip(&ends)
        .map(|(start, end)| format!(r#"{{"start": {start}, "end": {end}}}"#))
        .collect();
    let transcript = format!(r#"{{"strand": {strand}, "exons": [{}]}}"#, exons.join(", "));
    let transcript = loads(&transcript).expect("a transcript");
    let held = ie::transcript_two_window_positions(object(&transcript));
    numbers(
        held.into_iter()
            .filter(|position| (first..=last).contains(position)),
    )
}

/// The tool's answer for `item` over `arguments`.
fn answer(item: &str, arguments: &[Value]) -> Value {
    let body = || text(&arguments[0]);
    match item {
        "CLIPPED_NOTATION_RE.pattern" => Value::from(hgvs::HGVSC_CLIPPED_BODY),
        "INTRON_END_COPY_RE.pattern" => Value::from(ie::INTRON_END_COPY_BODY),
        "INTRON_END_DUPLICATION_RE.pattern" => Value::from(ie::INTRON_END_DUPLICATION_BODY),
        "INTRON_END_INSERTION_RE.pattern" => Value::from(ie::INTRON_END_INSERTION_BODY),
        "INTRON_START_INSERTION_RE.pattern" => Value::from(ie::INTRON_START_INSERTION_BODY),
        "INTRON_END_PAIRS" => {
            let mut pairs: Vec<String> = ie::INTRON_END_PAIRS
                .iter()
                .map(|pair| pair.iter().map(|base| char::from(*base)).collect())
                .collect();
            pairs.sort_unstable();
            Value::strings(&pairs)
        }
        "SHORT_INTRON_BASES" => number(ie::SHORT_INTRON_BASES),
        "coordinate_intron_offset" => optional(intron_offset_text(body())),
        "exonic_run" => exonic_run(body(), text(&arguments[1])).map_or(Value::Null, number),
        "SpliceWindows.holds" => held_positions(
            whole(&arguments[1]),
            &arguments[2],
            whole(&arguments[3]),
            whole(&arguments[4]),
        ),
        "transcript_intron_end_gaps" => {
            numbers(ie::transcript_intron_end_gaps(object(&arguments[0])))
        }
        "transcript_intron_start_gaps" => {
            numbers(ie::transcript_intron_start_gaps(object(&arguments[0])))
        }
        "coding_transcript_exons" => match ie::coding_transcript_exons(object(&arguments[0])) {
            Some(coding) => Value::from(vec![
                Value::from(coding.forward),
                Value::from(
                    coding
                        .exons
                        .iter()
                        .map(|(start, end)| numbers([*start, *end]))
                        .collect::<Vec<_>>(),
                ),
                number(coding.low),
                number(coding.high),
            ]),
            None => Value::Null,
        },
        "intron_end_insertion_form" => optional(ie::intron_end_insertion_form(body())),
        "copy_of_the_intron_last_base" => Value::from(ie::copy_of_the_intron_last_base(body())),
        "duplication_of_the_intron_last_bases" => {
            Value::from(ie::duplication_of_the_intron_last_bases(body()))
        }
        "insertion_that_cannot_end_the_intron" => {
            Value::from(ie::insertion_that_cannot_end_the_intron(body()))
        }
        "insertion_that_ends_in_an_intron_end_pair" => {
            Value::from(ie::insertion_that_ends_in_an_intron_end_pair(body()))
        }
        "insertion_after_the_exon_last_base" => {
            Value::from(ie::insertion_after_the_exon_last_base(body()))
        }
        "held_insertion_gaps" => numbers(ie::held_insertion_gaps(
            whole(&arguments[0]),
            whole(&arguments[1]),
            text(&arguments[2]),
            text(&arguments[3]),
            text(&arguments[4]),
        )),
        "reverse_complement" => Value::from(ie::reverse_complement(body())),
        "notation_inserted_bases" => {
            let pattern = match body() {
                "INTRON_END_INSERTION_RE" => &*ie::INTRON_END_INSERTION,
                "INTRON_START_INSERTION_RE" => &*ie::INTRON_START_INSERTION,
                other => panic!("no pattern {other}"),
            };
            Value::from(ie::notation_inserted_bases(pattern, text(&arguments[1])))
        }
        "lone_insertion_gaps" => numbers(ie::lone_insertion_gaps(
            whole(&arguments[0]),
            whole(&arguments[1]),
            text(&arguments[2]),
            text(&arguments[3]),
            arguments.get(4).map_or(0, whole),
            arguments.get(5).map_or(0, whole),
        )),
        "notation_duplicated_positions" => number(ie::notation_duplicated_positions(body())),
        other => panic!("the tool keeps no copy of {other}"),
    }
}

#[test]
fn the_items_kept_beside_the_format_comparator_answer_its_case_table() {
    let table = comparator_file("tests/data/builder_parity/cases.tsv");
    let not_held: BTreeMap<&str, usize> = NOT_HELD.into_iter().collect();
    let mut cases: BTreeMap<&str, usize> = BTreeMap::new();
    for (at, line) in table.lines().enumerate() {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), 3, "line {}", at + 1);
        *cases.entry(cells[0]).or_default() += 1;
        if not_held.contains_key(cells[0]) {
            continue;
        }
        let arguments = loads(cells[1]).expect("the arguments are JSON");
        let stated = loads(cells[2]).expect("the answer is JSON");
        let arguments = arguments.as_list().expect("the arguments are a list");
        assert_eq!(
            answer(cells[0], arguments),
            stated,
            "line {}: {} {}",
            at + 1,
            cells[0],
            cells[1]
        );
    }
    // Every item of the table is held or named as not held, with the count of its cases.
    let listed: BTreeMap<&str, usize> = HELD.into_iter().chain(NOT_HELD).collect();
    assert_eq!(listed.len(), HELD.len() + NOT_HELD.len());
    assert_eq!(cases, listed);
    assert_eq!(cases.values().sum::<usize>(), 403);
}
