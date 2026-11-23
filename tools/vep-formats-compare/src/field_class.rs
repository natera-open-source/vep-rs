// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `field-class` command: whether a field class accepts each of a list of differing
//! fields, for a program that classifies the differences of two outputs itself and asks
//! the comparison for one class.
//!
//! A question names a row by its Location, Allele and Feature, the differing field and its
//! two values, and carries the fields of the reference's row and of the candidate's row the
//! class reads. The answer is the class's predicate over the two rows as entries of those
//! members alone, with no record and no consequence terms: `rules::field_class_accepts`
//! over entries without a key, and for `hgvs_shift_sequence_edge_window`
//! `SequenceEdges::accepts` over entries keyed as default-format rows are.
//!
//! Another class is offered by a row of `OFFERED` with the fields of a row it reads and the
//! options its command line names, its name in `CLASSES` and among the choices `USAGE` and
//! `HELP` print, and, where it reads a side table `side` does not open, an option for the
//! table's file that `side` reads.
//!
//! The command reads no reference release: no class of `OFFERED` is one an Ensembl VEP
//! release corrects (`classes::corrected_under`), so each is applied under every release.

use crate::argparse::{
    path_text, standard_error, Action, ArgumentParser, Exit, Kind, Taken, Value as Given,
};
use crate::classes::{HGVS_SHIFT_SEQUENCE_EDGE_WINDOW, SPLICEAI_MULTI_ENTRY_CELL};
use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::input::read_text;
use crate::report::strerror;
use crate::rules::sequence_edges::SequenceEdges;
use crate::rules::spliceai::SpliceAiTable;
use crate::rules::{field_class_accepts, RuleContext, Side};
use crate::COMPONENT;
use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};

/// The first argument of a command line that is this command's.
pub const COMMAND: &str = "field-class";

/// The name the usage and the usage errors give the command.
pub const PROG: &str = "vep-formats-compare field-class";

/// The usage block and the whole help, laid out as the comparison's own at 80 columns.
const USAGE: &str = "\
usage: vep-formats-compare field-class [-h]
                                       --class {spliceai_multi_entry_cell,hgvs_shift_sequence_edge_window}
                                       [--spliceai-table FILE]
                                       [--reference-fasta FILE]
                                       [--json-cache DIR] --questions IN
                                       --answers OUT
";
const HELP: &str = "\
usage: vep-formats-compare field-class [-h]
                                       --class {spliceai_multi_entry_cell,hgvs_shift_sequence_edge_window}
                                       [--spliceai-table FILE]
                                       [--reference-fasta FILE]
                                       [--json-cache DIR] --questions IN
                                       --answers OUT

Answers, for each asked field of a row, whether a field class accepts the
difference of its two values.

options:
  -h, --help            show this help message and exit
  --class {spliceai_multi_entry_cell,hgvs_shift_sequence_edge_window}
                        the field class asked about
  --spliceai-table FILE
                        a score file the run's SpliceAI plugin read, its tabix
                        index beside it (repeatable: the snv= file and the
                        indel= file), which spliceai_multi_entry_cell reads
                        and a command line of that class names
  --reference-fasta FILE
                        the run's reference sequence, an uncompressed FASTA
                        beside its .fai index, which
                        hgvs_shift_sequence_edge_window reads and a command
                        line of that class names
  --json-cache DIR      the candidate JSON cache: its transcripts name the
                        gene of a row that states no SYMBOL, and
                        hgvs_shift_sequence_edge_window reads their spans and
                        exons, so a command line of that class names it
  --questions IN        the questions, one a line, tab-separated: Location,
                        Allele, Feature, the field's name, the reference's
                        value and the candidate's, then the fields of the
                        reference's row the class reads and those of the
                        candidate's row, - for a field a row lacks
  --answers OUT         takes one line a question, in their order: 1 where the
                        class accepts the field, else 0
";

/// The columns of a question before the fields of its two rows: Location, Allele, Feature,
/// the field's name, the reference's value and the candidate's.
pub const QUESTION_COLUMNS: usize = 6;

/// The cell of a question for a field its row lacks.
pub const LACKING: &str = "-";

/// The fields of a row a question of `spliceai_multi_entry_cell` carries, in the order of
/// their columns: the whole value, its nine members, the cutoff flag and the gene symbol.
pub const SPLICEAI_ROW_FIELDS: &[&str] = &[
    "SpliceAI_pred",
    "SpliceAI_pred_SYMBOL",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL",
    "SpliceAI_pred_DP_DG",
    "SpliceAI_pred_DP_DL",
    "SpliceAI_cutoff",
    "SYMBOL",
];

/// The cell of a question that states its row's Feature_type, which a row's key ends in
/// and no field of the row holds.
pub const FEATURE_TYPE: &str = "Feature_type";

/// The fields of a row a question of `hgvs_shift_sequence_edge_window` carries, in the order
/// of their columns: the row's Feature_type and strand, its notation, its protein notation
/// and its offset, then the length of its shift and its three position columns.
pub const SEQUENCE_EDGE_ROW_FIELDS: &[&str] = &[
    FEATURE_TYPE,
    "STRAND",
    "HGVSc",
    "HGVSp",
    "HGVS_OFFSET",
    "SHIFT_LENGTH",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
];

/// A class the command answers for, the fields of a row its questions carry and the options
/// a command line of the class names, in the order the usage gives them.
#[derive(Debug, Clone, Copy)]
pub struct Offered {
    pub class: &'static str,
    pub row_fields: &'static [&'static str],
    pub names: &'static [&'static str],
}

impl Offered {
    /// The columns of a question of the class: the six that name the field and one a field
    /// of each row.
    pub fn columns(&self) -> usize {
        QUESTION_COLUMNS + 2 * self.row_fields.len()
    }
}

/// The classes the command answers for.
pub const OFFERED: &[Offered] = &[
    Offered {
        class: SPLICEAI_MULTI_ENTRY_CELL,
        row_fields: SPLICEAI_ROW_FIELDS,
        names: &["--spliceai-table"],
    },
    Offered {
        class: HGVS_SHIFT_SEQUENCE_EDGE_WINDOW,
        row_fields: SEQUENCE_EDGE_ROW_FIELDS,
        names: &["--reference-fasta", "--json-cache"],
    },
];

/// The names of `OFFERED`, the choices of `--class`.
const CLASSES: &[&str] = &[SPLICEAI_MULTI_ENTRY_CELL, HGVS_SHIFT_SEQUENCE_EDGE_WINDOW];

/// The parsed command line.
#[derive(Debug)]
pub struct Namespace {
    pub class: &'static Offered,
    pub spliceai_table: Vec<PathBuf>,
    pub reference_fasta: Option<PathBuf>,
    pub json_cache: Option<PathBuf>,
    pub questions: PathBuf,
    pub answers: PathBuf,
}

/// The command's options, in the order its usage names them.
pub fn build_parser() -> ArgumentParser {
    let mut ap = ArgumentParser::new(PROG, USAGE, HELP);
    ap.add_argument("--class", Action::Store, Kind::Text, Some(CLASSES));
    ap.add_argument("--spliceai-table", Action::Append, Kind::Path, None);
    ap.add_argument("--reference-fasta", Action::Store, Kind::Path, None);
    ap.add_argument("--json-cache", Action::Store, Kind::Path, None);
    ap.add_argument("--questions", Action::Store, Kind::Path, None);
    ap.add_argument("--answers", Action::Store, Kind::Path, None);
    ap
}

/// The command line `argv`, the command's name left out, as `parser` reads it, or how the
/// parser ends the run: an option given twice keeps its last value, `--spliceai-table` keeps
/// every value in order, and a command line that lacks one of `--class`, `--questions` and
/// `--answers`, or an option its class names (`Offered::names`), is refused.
pub fn parse_args(
    parser: &ArgumentParser,
    argv: &[OsString],
) -> std::result::Result<Namespace, Exit> {
    let (mut class, mut questions, mut answers) = (None, None, None);
    let (mut spliceai_table, mut reference_fasta, mut json_cache) = (Vec::new(), None, None);
    for Taken { argument, value } in parser.parse_args(argv)? {
        match (parser.dest(argument), value) {
            ("class", Given::Text(text)) => {
                class = OFFERED.iter().find(|offered| offered.class == text);
            }
            ("spliceai_table", Given::Path(path)) => spliceai_table.push(path),
            ("reference_fasta", Given::Path(path)) => reference_fasta = Some(path),
            ("json_cache", Given::Path(path)) => json_cache = Some(path),
            ("questions", Given::Path(path)) => questions = Some(path),
            ("answers", Given::Path(path)) => answers = Some(path),
            _ => {}
        }
    }
    let named = |option: &str| class.is_some_and(|offered| offered.names.contains(&option));
    let lacking: Vec<&str> = [
        ("--class", class.is_none()),
        (
            "--spliceai-table",
            named("--spliceai-table") && spliceai_table.is_empty(),
        ),
        (
            "--reference-fasta",
            named("--reference-fasta") && reference_fasta.is_none(),
        ),
        (
            "--json-cache",
            named("--json-cache") && json_cache.is_none(),
        ),
        ("--questions", questions.is_none()),
        ("--answers", answers.is_none()),
    ]
    .iter()
    .filter(|(_, lacks)| *lacks)
    .map(|(option, _)| *option)
    .collect();
    match (class, questions, answers) {
        (Some(class), Some(questions), Some(answers)) if lacking.is_empty() => Ok(Namespace {
            class,
            spliceai_table,
            reference_fasta,
            json_cache,
            questions,
            answers,
        }),
        _ => Err(Exit::Error(format!(
            "the following arguments are required: {}",
            lacking.join(", ")
        ))),
    }
}

/// The side tables the class reads, of those the command line names: the reference sequence
/// and the candidate cache's transcripts for `hgvs_shift_sequence_edge_window`, else the
/// score files, with the candidate cache that names a row's gene.
fn side(args: &Namespace) -> Result<Side> {
    if args.class.class == HGVS_SHIFT_SEQUENCE_EDGE_WINDOW {
        let edges = match (&args.reference_fasta, &args.json_cache) {
            (Some(fasta), Some(cache)) => Some(SequenceEdges::new(fasta, cache)?),
            _ => None,
        };
        return Ok(Side {
            sequence_edges: edges,
            ..Side::default()
        });
    }
    let tables: Vec<&Path> = args.spliceai_table.iter().map(PathBuf::as_path).collect();
    Ok(Side {
        spliceai_table: Some(SpliceAiTable::new(&tables, args.json_cache.as_deref())?),
        ..Side::default()
    })
}

/// The row of a question as an entry: its Location, Allele and Feature, and of `names` the
/// fields whose cell is not `LACKING`, each with its cell as its value. The cell of
/// `FEATURE_TYPE` is no field: where `names` holds it the entry carries the key a
/// default-format row has, its Location, Allele, Feature and that cell.
fn row(location: &str, allele: &str, feature: &str, names: &[&str], cells: &[&str]) -> Entry {
    let cells = || names.iter().zip(cells);
    let kind = cells().find(|(name, _)| **name == FEATURE_TYPE);
    Entry {
        key: kind.map_or_else(String::new, |(_, kind)| {
            format!("{location}|{allele}|{feature}|{kind}")
        }),
        location: Some(location.to_string()),
        allele: allele.to_string(),
        feature: feature.to_string(),
        fields: cells()
            .filter(|(name, cell)| **name != FEATURE_TYPE && **cell != LACKING)
            .map(|(name, cell)| (name.to_string(), cell.to_string()))
            .collect(),
        ..Entry::default()
    }
}

/// The answers to the questions of `text`, one line a question in their order: `1` where
/// `offered`'s class accepts the question's field, else `0`. A line of another number of
/// columns than the class's questions have is refused; `path` names the questions' file in
/// that refusal. `hgvs_shift_sequence_edge_window` is asked of the side's reference sequence
/// and transcripts, as a comparison asks it, with the row's Location for its span.
pub fn answers(offered: &Offered, text: &str, path: &Path, side: &Side) -> Result<String> {
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side,
    };
    let fields = offered.row_fields.len();
    let mut out = String::new();
    for (at, line) in text.split_terminator('\n').enumerate() {
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() != offered.columns() {
            return Err(Error::comparator(format!(
                "--questions: {} line {} has {} columns, a question of {} has {}",
                path_text(path),
                at + 1,
                cells.len(),
                offered.class,
                offered.columns()
            )));
        }
        let (asked, rows) = cells.split_at(QUESTION_COLUMNS);
        let (reference, candidate) = rows.split_at(fields);
        let e = row(asked[0], asked[1], asked[2], offered.row_fields, reference);
        let a = row(asked[0], asked[1], asked[2], offered.row_fields, candidate);
        let accepts = if offered.class == HGVS_SHIFT_SEQUENCE_EDGE_WINDOW {
            let edges = side.sequence_edges.as_ref();
            edges.is_some_and(|edges| edges.accepts(&e, Some(&a), None, asked[3]))
        } else {
            field_class_accepts(
                offered.class,
                &e,
                asked[3],
                asked[4],
                asked[5],
                Some(&a),
                &context,
            )?
        };
        out.push_str(if accepts { "1\n" } else { "0\n" });
    }
    Ok(out)
}

/// Opens the side tables, answers every question and writes the answers' file, which is
/// opened once every question is answered.
fn run(args: &Namespace) -> Result<()> {
    let side = side(args)?;
    let questions = read_text(&args.questions)?;
    let text = answers(args.class, &questions, &args.questions, &side)?;
    std::fs::write(&args.answers, text).map_err(|error| {
        Error::comparator(format!(
            "cannot write {}: {}",
            path_text(&args.answers),
            strerror(&error)
        ))
    })
}

/// Runs the command over `argv`, the command's name left out, and returns the exit code: 0
/// once every question is answered and the answers are written; 2 for a refused command
/// line, a side table, a reference sequence, a cache or a questions' file that cannot be
/// read, a question of another number of columns and an answers' file that cannot be
/// written, each with nothing written to the answers' file.
pub fn main(argv: &[OsString]) -> i32 {
    let parser = build_parser();
    let args = match parse_args(&parser, argv) {
        Ok(args) => args,
        Err(Exit::Help) => {
            let _ = std::io::stdout().write_all(parser.format_help().as_bytes());
            return 0;
        }
        Err(Exit::Error(message)) => {
            let text = standard_error(&parser.error_text(&message));
            let _ = std::io::stderr().write_all(text.as_bytes());
            return 2;
        }
    };
    match run(&args) {
        Ok(()) => 0,
        Err(error) => {
            let text = standard_error(&format!("ERROR: [{COMPONENT}] {error}\n"));
            let _ = std::io::stderr().write_all(text.as_bytes());
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::spliceai::{SPLICEAI_CUTOFF_FIELD, SPLICEAI_MEMBERS, SPLICEAI_VALUE_FIELD};

    #[test]
    fn the_row_fields_of_the_spliceai_class_are_the_value_its_members_the_flag_and_the_symbol() {
        let mut names = vec![SPLICEAI_VALUE_FIELD.to_string()];
        names.extend(
            SPLICEAI_MEMBERS
                .iter()
                .map(|member| format!("{SPLICEAI_VALUE_FIELD}_{member}")),
        );
        names.push(SPLICEAI_CUTOFF_FIELD.to_string());
        names.push("SYMBOL".to_string());
        assert_eq!(SPLICEAI_ROW_FIELDS, names);
        assert_eq!(OFFERED[0].columns(), 30);
    }

    #[test]
    fn every_offered_class_is_a_choice_and_a_field_class() {
        let offered: Vec<&str> = OFFERED.iter().map(|offered| offered.class).collect();
        assert_eq!(offered, CLASSES);
        for class in CLASSES {
            assert!(crate::classes::tables().in_field_classes(class), "{class}");
            assert!(USAGE.contains(class) && HELP.contains(class), "{class}");
            // No release corrects an offered class, so the command takes no reference
            // release.
            assert_eq!(
                crate::classes::tables().corrected_in(class),
                None,
                "{class}"
            );
        }
        assert!(HELP.starts_with(USAGE));
    }

    #[test]
    fn a_row_holds_the_fields_its_question_writes() {
        let cells = ["-", "", "0.10"];
        let e = row("21:100", "C", "ENST1", &["a", "b", "c"], &cells);
        assert_eq!(e.key, "");
        assert_eq!(e.record, None);
        assert_eq!(e.location.as_deref(), Some("21:100"));
        assert_eq!((e.allele.as_str(), e.feature.as_str()), ("C", "ENST1"));
        assert_eq!(e.consequence, "");
        let held: Vec<(&str, &str)> = e
            .fields
            .iter()
            .map(|(name, value)| (name.as_str(), value.as_str()))
            .collect();
        assert_eq!(held, [("b", ""), ("c", "0.10")]);
    }

    #[test]
    fn a_row_that_states_its_feature_type_carries_the_key_of_a_default_format_row() {
        let names = [FEATURE_TYPE, "STRAND", "HGVSc"];
        let e = row("edge:10", "-", "ENST1", &names, &["Transcript", "-1", "-"]);
        assert_eq!(e.key, "edge:10|-|ENST1|Transcript");
        let held: Vec<(&str, &str)> = e
            .fields
            .iter()
            .map(|(name, value)| (name.as_str(), value.as_str()))
            .collect();
        assert_eq!(held, [("STRAND", "-1")]);
    }

    #[test]
    fn a_question_of_the_sequence_edge_class_carries_nine_cells_of_each_row() {
        let offered = OFFERED[1];
        assert_eq!(offered.class, HGVS_SHIFT_SEQUENCE_EDGE_WINDOW);
        assert_eq!(
            offered.row_fields,
            [
                "Feature_type",
                "STRAND",
                "HGVSc",
                "HGVSp",
                "HGVS_OFFSET",
                "SHIFT_LENGTH",
                "cDNA_position",
                "CDS_position",
                "Protein_position"
            ]
        );
        assert_eq!(offered.columns(), 24);
        assert_eq!(offered.names, ["--reference-fasta", "--json-cache"]);
        assert_eq!(OFFERED[0].names, ["--spliceai-table"]);
    }
}
