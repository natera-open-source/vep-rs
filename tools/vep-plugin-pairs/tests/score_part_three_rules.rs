// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tool's own classes consulted by the `score` program end to end, on a planted
//! registry: the UTRAnnotator window class on a row whose record the run's input names
//! under an identifier, the class of the shifted CDS coordinates on a record of two
//! alternates, the two classes of a plugin value that follows a term pair, and the LOFTEE
//! class of a GERP walk without the conservation database. The two classes of a
//! mitochondrial record are consulted under the tool's own registry, each under its row of
//! the tree's class table.

mod support;

use serde_json::{json, Value};
use support::{
    allele, counts, feature, feature_residuals, message, raw_counts, row_with, try_score,
    PlantedRegistry, Run, Scratch,
};
use vep_plugin_pairs::registry::BuiltInRegistry;
use vep_plugin_pairs::rules::{
    BuiltIn, Rules, LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE,
    PLUGIN_HGVS_SHIFTED_CDS_COORDINATES, PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN,
    PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR, UTRANNOTATOR_GAINED_START_WINDOW,
};

/// The tree's class table.
const COMMITTED_CLASS_TABLE: &str =
    include_str!("../../../scripts/concordance/plugin_class_table.tsv");

/// The report of a run of one engine under the tool's own classes.
fn score(
    scratch: &Scratch,
    registry: &PlantedRegistry,
    reference: &str,
    engine: &str,
    run: &Run,
) -> Value {
    try_score(
        scratch,
        reference,
        &[("engine", engine)],
        run,
        registry,
        &BuiltIn,
    )
    .unwrap_or_else(|stop| panic!("the run was stopped: {}", message(&stop)))
}

/// The committed table's one row of `klass`, as a table of its own.
fn committed_row(klass: &str) -> String {
    let declared: Vec<&str> = COMMITTED_CLASS_TABLE
        .lines()
        .filter(|line| line.starts_with(&format!("{klass}\t")))
        .collect();
    assert_eq!(declared.len(), 1, "{klass}");
    format!("{}\n", declared[0])
}

const UTR_START_LOST: &str = "5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;5UTR_annotation=CapDistanceToStart=63:type=uORF";
const UTR_START_GAINED: &str = "5UTR_consequence=5_prime_UTR_premature_start_codon_gain_variant;5UTR_annotation=CapDistanceToStart=64:type=uORF";

/// An insertion's row under the identifier `rs1` with its two UTRAnnotator fields.
fn inserted(allele: &str, fields: &str) -> String {
    row_with("1:100-101", allele, "ENST1", fields, "64-65", "rs1")
}

/// One input of a grid: what it shows, its records, or no input at all, and the keys the
/// window class sets aside under it.
type GridCase<'a> = (&'a str, Option<&'a [[&'a str; 5]]>, i64);

/// The window class under the committed table's own row, on a registry of UTRAnnotator's
/// two fields, and an insertion's reference and engine rows under the identifier `rs1`:
/// the reference reports a start codon lost, the engine one gained.
struct WindowRows {
    scratch: Scratch,
    registry: PlantedRegistry,
    table: String,
    reference: String,
    engine: String,
}

impl WindowRows {
    fn new(tag: &str) -> WindowRows {
        let scratch = Scratch::new(tag);
        let reference = scratch.output("reference.txt", &[inserted("CAT", UTR_START_LOST)]);
        let engine = scratch.output("engine.txt", &[inserted("CAT", UTR_START_GAINED)]);
        WindowRows {
            scratch,
            registry: PlantedRegistry::new(&[(
                "UTRAnnotator",
                &["5UTR_annotation", "5UTR_consequence"],
            )]),
            table: committed_row(UTRANNOTATOR_GAINED_START_WINDOW),
            reference,
            engine,
        }
    }

    /// A run under the class's row with the input of `records`, or with none; the run's
    /// files stand under `name`.
    fn run(&self, name: &str, records: Option<&[[&str; 5]]>) -> Run {
        let run = Run::new()
            .table(&self.table)
            .plugins(&["UTRAnnotator"])
            .residuals(&self.scratch.path(&format!("{name}/residuals")));
        match records {
            Some(records) => {
                run.run_input(&self.scratch.vcf(&format!("{name}/input.vcf"), records))
            }
            None => run,
        }
    }

    /// The report of `reference` against `engine` under `run`.
    fn scored(&self, reference: &str, engine: &str, run: &Run) -> Value {
        score(&self.scratch, &self.registry, reference, engine, run)
    }

    /// Scores the insertion's rows under each input of `grid`: the row's two fields are
    /// set aside, or all four values stay residual rows.
    fn holds(&self, grid: &[GridCase]) {
        for (number, (what, records, keys)) in grid.iter().enumerate() {
            let report = self.scored(
                &self.reference,
                &self.engine,
                &self.run(&format!("grid{number}"), *records),
            );
            let item = feature(&report, "UTRAnnotator");
            assert_eq!(
                item["masked"],
                json!({UTRANNOTATOR_GAINED_START_WINDOW: keys}),
                "{what}"
            );
            assert_eq!(item["residual_rows"], 4 - 2 * keys, "{what}");
        }
    }
}

const OWN: [&str; 5] = ["1", "100", "rs1", "C", "CCAT"];
const OTHER_ALLELE: [&str; 5] = ["1", "100", "rs1", "C", "CCAG"];
const ELSEWHERE: [&str; 5] = ["1", "150", "rs1", "C", "CCAT"];
const JOINED: [&str; 5] = ["1", "100", "rs1;rs2", "C", "CCAT"];

/// An insertion's row under an identifier is handed the one input line that writes it at
/// its Location with its allele, where several lines carry the identifier, and its two
/// fields are set aside as where one line alone carries it. Two lines that both write the
/// row, lines of which none does and a run without the input leave both fields residuals,
/// and two records that share the identifier and the Location are each set aside by their
/// own line.
#[test]
fn row_under_an_identifier_several_input_lines_carry_is_handed_the_record_that_writes_it() {
    let rows = WindowRows::new("several-lines");
    rows.holds(&[
        ("one line", Some(&[OWN]), 2),
        (
            "two lines at one place with different alleles",
            Some(&[OTHER_ALLELE, OWN]),
            2,
        ),
        ("two lines at different places", Some(&[ELSEWHERE, OWN]), 2),
        ("two lines that both write the row", Some(&[OWN, OWN]), 0),
        (
            "a line whose two alternates write the row in two ways, beside the insertion's line",
            Some(&[["1", "100", "rs1", "CA", "CAT,CCATA"], OWN]),
            0,
        ),
        (
            "two lines of which none writes the row",
            Some(&[OTHER_ALLELE, ELSEWHERE]),
            0,
        ),
        ("no input", None, 0),
    ]);
    // Two records under one identifier at one Location: each row reads its own line.
    let reference = rows.scratch.output(
        "two-records/reference.txt",
        &[
            inserted("CAT", UTR_START_LOST),
            inserted("CAG", UTR_START_LOST),
        ],
    );
    let engine = rows.scratch.output(
        "two-records/engine.txt",
        &[
            inserted("CAT", UTR_START_GAINED),
            inserted("CAG", UTR_START_GAINED),
        ],
    );
    let both: &[[&str; 5]] = &[OWN, OTHER_ALLELE];
    let report = rows.scored(&reference, &engine, &rows.run("two-records", Some(both)));
    let item = feature(&report, "UTRAnnotator");
    assert_eq!(counts(item), [4, 4, 0, 0, 0, 0]);
    assert_eq!(item["masked"], json!({UTRANNOTATOR_GAINED_START_WINDOW: 4}));
    assert!(
        feature_residuals(&rows.scratch.path("two-records/residuals"), "UTRAnnotator").is_empty()
    );
    // The scores are the same on one thread and on four.
    for threads in [1, 4] {
        let run = rows
            .run(&format!("two-records-{threads}"), Some(both))
            .threads(threads);
        let again = rows.scored(&reference, &engine, &run);
        assert_eq!(again["scores"], report["scores"], "{threads} threads");
    }
}

/// An insertion's row prints the first identifier of its record's ID column, so the input
/// line that joins two identifiers, or five, is its record, and its two fields are set
/// aside as under a line of that identifier alone. A line that carries the row's
/// identifier as a later element, or behind a first element that is a dot, empty or `0`,
/// is no record of the row, and both fields stay residuals.
#[test]
fn row_under_the_first_identifier_of_a_several_identifier_line_is_handed_its_record() {
    let rows = WindowRows::new("first-identifier");
    let grid = [
        ("one identifier", "rs1", 2),
        ("two identifiers", "rs1;rs2", 2),
        ("five identifiers", "rs1;rs2;rs3;rs4;rs5", 2),
        ("the identifier as a later element", "rs2;rs1", 0),
        ("a first element that is a dot", ".;rs1", 0),
        ("a first element that is empty", ";rs1", 0),
        ("a first element that is zero", "0;rs1", 0),
    ];
    let inputs = grid.map(|(_, name, _)| [["1", "100", name, "C", "CCAT"]]);
    let cases: Vec<GridCase> = grid
        .iter()
        .zip(&inputs)
        .map(|((what, _, keys), records)| (*what, Some(&records[..]), *keys))
        .collect();
    rows.holds(&cases);
    // The line of five identifiers is the row's record on one thread and on four alike.
    for threads in [1, 4] {
        let run = rows
            .run(&format!("threads{threads}"), Some(&inputs[2]))
            .threads(threads);
        let again = rows.scored(&rows.reference, &rows.engine, &run);
        let item = feature(&again, "UTRAnnotator");
        assert_eq!(
            item["masked"],
            json!({UTRANNOTATOR_GAINED_START_WINDOW: 2}),
            "{threads} threads"
        );
        assert_eq!(item["residual_rows"], 0, "{threads} threads");
    }
}

/// A line of one identifier and a line that joins that identifier with another stand under
/// one name, as two joined lines with one first identifier do. A row under the name is
/// handed the one of them that writes it at its Location with its allele, whichever of the
/// two it is and whether the other stands elsewhere or spells another allele; where both
/// write the row nothing is handed and both fields stay residuals.
#[test]
fn lines_that_share_a_first_identifier_are_told_apart_by_place_and_allele() {
    let rows = WindowRows::new("shared-first-identifier");
    rows.holds(&[
        (
            "the joined line writes the row, the plain line stands elsewhere",
            Some(&[["1", "150", "rs1", "CA", "C"], JOINED]),
            2,
        ),
        (
            "the joined line writes the row, the plain line spells another allele",
            Some(&[OTHER_ALLELE, JOINED]),
            2,
        ),
        (
            "the plain line writes the row, the joined line stands elsewhere",
            Some(&[OWN, ["1", "150", "rs1;rs2", "CA", "C"]]),
            2,
        ),
        (
            "the plain line writes the row, the joined line spells another allele",
            Some(&[OWN, ["1", "100", "rs1;rs2", "C", "CCAG"]]),
            2,
        ),
        ("both lines write the row", Some(&[OWN, JOINED]), 0),
        (
            "two joined lines of which one writes the row",
            Some(&[["1", "150", "rs1;rs3", "CA", "C"], JOINED]),
            2,
        ),
        (
            "two joined lines that both write the row",
            Some(&[["1", "100", "rs1;rs3", "C", "CCAT"], JOINED]),
            0,
        ),
    ]);
}

/// A made-up coding sequence whose three bases after coding base 4 repeat it, the 3' UTR
/// behind it, and the proteins the Downstream plugin of each engine writes for the two
/// alleles of the record: one base and two bases inserted after coding base 4.
const CODING: &str = "ATGAAAACTCTTGCTGTCCGTGGTTAA";
const THREE_PRIME_UTR: &str = "GGGCCCAAATAGCCC";
const REFERENCE_FIRST_ALLELE: &str = "NSCCPWLRAQIA";
const REFERENCE_LATER_ALLELE: &str = "KLLLSVVKGPNS";
const ENGINE_FIRST_ALLELE: &str = "KNSCCPWLRAQIA";
const ENGINE_LATER_ALLELE: &str = "KKLLLSVVKGPNS";

/// A location with its allele, then the cDNA and the CDS positions of a row there.
type Place = [&'static str; 4];

/// A row of the transcript `ENST1` at `place` with the Consequence `terms`.
fn coding_row(place: Place, terms: &str, extra: &str, uploaded: &str) -> String {
    let [location, allele, cdna, cds] = place;
    format!(
        "{uploaded}\t{location}\t{allele}\tENSG1\tENST1\tTranscript\t{terms}\t{cdna}\t{cds}\t-\t-\t-\t-\t{extra}\n"
    )
}

/// The later allele's reference row prints the CDS pair three bases on, the engine's the
/// record's own, and `plugin_hgvs_shifted_cds_coordinates` takes the row's value only as a
/// later allele's. The class is handed the record's rows with the input's line, which
/// states the alleles in order, so both alleles' values are set aside; without the input
/// the later allele's stays.
#[test]
fn later_alleles_downstream_value_is_set_aside_where_the_input_states_the_records_alleles() {
    let shifted = PLUGIN_HGVS_SHIFTED_CDS_COORDINATES;
    let scratch = Scratch::new("later-allele");
    let registry = PlantedRegistry::new(&[("Downstream", &["DownstreamProtein"])]);
    let transcript = json!({
        "stable_id": "ENST1",
        "strand": 1,
        "cdna_coding_start": 1,
        "cdna_coding_end": CODING.len(),
        "exons": [{"start": 1001, "end": 1042}],
        "variation_effect_feature_cache": {
            "translateable_seq": CODING,
            "three_prime_utr": THREE_PRIME_UTR,
            "codon_table": 1,
            "mapper": {
                "pairs": [
                    {"from_start": 1, "from_end": 42, "to_start": 1001, "to_end": 1042, "ori": 1}
                ]
            },
        },
    });
    scratch.file(
        "cache/transcripts/21/1-1000000.json",
        json!([transcript]).to_string(),
    );
    let frameshift = |allele: &'static str, cds: &'static str, protein: &str| {
        coding_row(
            ["21:1004-1005", allele, "-", cds],
            "frameshift_variant",
            &format!(
                "IMPACT=HIGH;STRAND=1;HGVS_OFFSET=3;DownstreamProtein={protein};ProteinLengthChange=6"
            ),
            "v",
        )
    };
    let reference = scratch.output(
        "reference.txt",
        &[
            frameshift("A", "4-5", REFERENCE_FIRST_ALLELE),
            frameshift("AA", "7-8", REFERENCE_LATER_ALLELE),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            frameshift("A", "4-5", ENGINE_FIRST_ALLELE),
            frameshift("AA", "4-5", ENGINE_LATER_ALLELE),
        ],
    );
    let run = Run::new()
        .table(&format!(
            "{shifted}\tDownstream\tDownstreamProtein\tmet\t-\n"
        ))
        .plugins(&["Downstream"])
        .cache(&scratch.path("cache"));
    let without = score(&scratch, &registry, &reference, &engine, &run);
    assert_eq!(
        feature(&without, "Downstream")["masked"],
        json!({shifted: 1})
    );
    let residuals = scratch.path("residuals");
    let run_input = scratch.vcf("input.vcf", &[["21", "1004", "v", "A", "AA,AAA"]]);
    let report = score(
        &scratch,
        &registry,
        &reference,
        &engine,
        &run.residuals(&residuals).run_input(&run_input),
    );
    assert_eq!(
        feature(&report, "Downstream")["masked"],
        json!({shifted: 2})
    );
    assert!(feature_residuals(&residuals, "Downstream").is_empty());
}

/// The plugins whose values follow a row's terms, each with its fields.
const TERM_FOLLOWING_PLUGINS: [(&str, &[&str]); 3] = [
    ("Downstream", &["DownstreamProtein", "ProteinLengthChange"]),
    ("NMD", &["NMD"]),
    ("LoFTEE", &["LoF", "LoF_filter", "LoF_flags", "LoF_info"]),
];
/// The plugins the committed table's row of `plugin_value_follows_stop_codon_term_pair`
/// names: those of `TERM_FOLLOWING_PLUGINS` and the four whose scores follow a row's terms.
const COMMITTED_ROW_PLUGINS: [(&str, &[&str]); 7] = [
    TERM_FOLLOWING_PLUGINS[0],
    TERM_FOLLOWING_PLUGINS[1],
    TERM_FOLLOWING_PLUGINS[2],
    ("REVEL", &["REVEL"]),
    ("AlphaMissense", &["am_class", "am_pathogenicity"]),
    ("ClinPred", &["ClinPred"]),
    ("dbNSFP", &["dbNSFP_SIFT_score"]),
];
const ESCAPING: &str = "NMD=NMD_escaping_variant";
const UNTESTED_STOP: Place = ["19:13208046-13208050", "-", "1705-?", "1630-?"];
const KEPT_STOP: Place = ["19:13208146-13208147", "TTGA", "7054-7055", "6977-6978"];
const TESTED_STOP: Place = ["2:190970701-190970707", "-", "961-967", "264-270"];
const NO_PAIR: Place = ["7:500-502", "-", "100-102", "40-42"];

/// What a class leaves of one source: the counts, the keys set aside and the residual
/// lines.
type Left<'a> = (&'a str, [i64; 6], i64, Vec<&'a str>);

/// Holds the report's entry and the residual file of each source to what `klass` leaves.
fn assert_left(report: &Value, residuals: &str, klass: &str, expected: &[Left], what: &str) {
    for (source, counted, keys, left) in expected {
        let item = feature(report, source);
        assert_eq!(counts(item), *counted, "{what} {source}");
        assert_eq!(item["masked"], json!({klass: keys}), "{what} {source}");
        assert_eq!(
            feature_residuals(residuals, source),
            *left,
            "{what} {source}"
        );
    }
}

/// A Downstream, NMD or LOFTEE value one side alone writes is set aside on a row whose two
/// term sets are a stop codon pair the class reads (the untested stop codon, the insertion
/// that keeps the stop codon, the three-base test), on whichever side's terms the plugin
/// runs on. A `LoF_flags` the table's row leaves out and the same values on a row of no
/// stop codon pair stay residuals.
#[test]
fn plugin_value_that_follows_a_stop_codon_term_pair_is_set_aside() {
    let klass = PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR;
    let scratch = Scratch::new("stop-codon-pair");
    let registry = PlantedRegistry::new(&TERM_FOLLOWING_PLUGINS);
    let reference = scratch.output(
        "reference.txt",
        &[
            coding_row(UNTESTED_STOP, "stop_retained_variant", "IMPACT=LOW", "v"),
            coding_row(
                KEPT_STOP,
                "stop_gained,frameshift_variant",
                &format!(
                    "IMPACT=HIGH;DownstreamProtein=MK;ProteinLengthChange=-3;{ESCAPING};LoF=LC;LoF_info=PERCENTILE:0.5"
                ),
                "v",
            ),
            coding_row(TESTED_STOP, "stop_retained_variant", "IMPACT=LOW", "v"),
            coding_row(NO_PAIR, "inframe_deletion", "IMPACT=MODERATE", "v"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            coding_row(
                UNTESTED_STOP,
                "frameshift_variant",
                &format!(
                    "IMPACT=HIGH;DownstreamProtein=MKV;ProteinLengthChange=-12;{ESCAPING};LoF=HC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.97"
                ),
                "v",
            ),
            coding_row(
                KEPT_STOP,
                "inframe_insertion,stop_retained_variant",
                "IMPACT=MODERATE",
                "v",
            ),
            coding_row(
                TESTED_STOP,
                "frameshift_variant",
                &format!("IMPACT=HIGH;{ESCAPING}"),
                "v",
            ),
            coding_row(
                NO_PAIR,
                "frameshift_variant",
                &format!("IMPACT=HIGH;DownstreamProtein=MR;{ESCAPING}"),
                "v",
            ),
        ],
    );
    let table = format!(
        "{klass}\tDownstream;NMD;LoFTEE\tDownstreamProtein,ProteinLengthChange;NMD;LoF,LoF_info,LoF_filter\tmet\t-\n"
    );
    let residuals = scratch.path("residuals");
    let report = score(
        &scratch,
        &registry,
        &reference,
        &engine,
        &Run::new().table(&table).residuals(&residuals),
    );
    assert_left(
        &report,
        &residuals,
        klass,
        &[
            (
                "Downstream",
                [2, 3, 0, 0, 1, 0],
                4,
                vec!["engine\t7:500-502\t-\tENST1\tDownstreamProtein\tMR"],
            ),
            (
                "NMD",
                [1, 3, 0, 0, 1, 0],
                3,
                vec!["engine\t7:500-502\t-\tENST1\tNMD\tNMD_escaping_variant"],
            ),
            (
                "LoFTEE",
                [2, 4, 0, 0, 1, 0],
                5,
                vec!["engine\t19:13208046-13208050\t-\tENST1\tLoF_flags\tSINGLE_EXON"],
            ),
        ],
        "the table's row",
    );
}

const TWO_BASES_AT_STOP: Place = ["19:13208246-13208247", "CG", "7054-7055", "6977-6978"];
const THREE_BASES_AT_STOP: Place = ["19:13208346-13208347", "CGT", "7054-7055", "6977-6978"];

/// Under an identifier the class reads an insertion that keeps the stop codon from the lines
/// of the input that carry the identifier: an NMD value the reference alone writes is set
/// aside on a row whose line inserts two bases and stays a residual on a row whose line
/// replaces two bases by three. Scored without the input each row is read by its allele's
/// length, and the two change places.
#[test]
fn plugin_value_on_an_insertion_that_keeps_the_stop_codon_follows_the_lines_of_the_input() {
    let klass = PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR;
    assert!(BuiltIn.record_reference(klass));
    let scratch = Scratch::new("stop-codon-pair-records");
    let registry = PlantedRegistry::new(&TERM_FOLLOWING_PLUGINS);
    let places = [(TWO_BASES_AT_STOP, "rs91"), (THREE_BASES_AT_STOP, "rs92")];
    let rows = |terms: &str, extra: &str| {
        places.map(|(place, name)| coding_row(place, terms, extra, name))
    };
    let reference = scratch.output(
        "reference.txt",
        &rows(
            "stop_gained,frameshift_variant",
            &format!("IMPACT=HIGH;{ESCAPING}"),
        ),
    );
    let engine = scratch.output(
        "engine.txt",
        &rows("inframe_insertion,stop_retained_variant", "IMPACT=MODERATE"),
    );
    let run_input = scratch.vcf(
        "input.vcf",
        &[
            ["19", "13208246", "rs91", "A", "ACG"],
            ["19", "13208346", "rs92", "AC", "CGT"],
        ],
    );
    let table = format!(
        "{klass}\tDownstream;NMD;LoFTEE\tDownstreamProtein,ProteinLengthChange;NMD;LoF,LoF_info,LoF_filter\tmet\t-\n"
    );
    let left = |place: Place| {
        format!(
            "reference\t{}\t{}\tENST1\tNMD\tNMD_escaping_variant",
            place[0], place[1]
        )
    };
    for (name, given, residual) in [
        ("given", Some(&run_input), left(THREE_BASES_AT_STOP)),
        ("without", None, left(TWO_BASES_AT_STOP)),
    ] {
        let residuals = scratch.path(&format!("residuals-{name}"));
        let mut run = Run::new().table(&table).residuals(&residuals);
        if let Some(input) = given {
            run = run.run_input(input);
        }
        let report = score(&scratch, &registry, &reference, &engine, &run);
        assert_eq!(
            feature(&report, "NMD")["masked"],
            json!({klass: 1}),
            "{name}"
        );
        assert_eq!(feature_residuals(&residuals, "NMD"), [residual], "{name}");
    }
}

/// Under a class table whose row names `LoF_flags`, a flag list one side alone writes is
/// set aside on a row whose two term sets are a stop codon pair the class reads, on
/// whichever side's terms LOFTEE runs on, and the same list on a row of no stop codon pair
/// stays a residual. Under a row that leaves the field out every flag list stays a
/// residual.
#[test]
fn loftee_flags_that_follow_a_stop_codon_term_pair_are_set_aside_where_the_class_table_names_the_field(
) {
    let klass = PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR;
    let scratch = Scratch::new("stop-codon-pair-flags");
    let registry = PlantedRegistry::new(&TERM_FOLLOWING_PLUGINS[2..]);
    let reference = scratch.output(
        "reference.txt",
        &[
            coding_row(UNTESTED_STOP, "stop_retained_variant", "IMPACT=LOW", "v"),
            coding_row(
                KEPT_STOP,
                "stop_gained,frameshift_variant",
                "IMPACT=HIGH;LoF=HC;LoF_flags=SINGLE_EXON",
                "v",
            ),
            coding_row(NO_PAIR, "inframe_deletion", "IMPACT=MODERATE", "v"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            coding_row(
                UNTESTED_STOP,
                "frameshift_variant",
                "IMPACT=HIGH;LoF=HC;LoF_flags=PHYLOCSF_WEAK",
                "v",
            ),
            coding_row(
                KEPT_STOP,
                "inframe_insertion,stop_retained_variant",
                "IMPACT=MODERATE",
                "v",
            ),
            coding_row(
                NO_PAIR,
                "frameshift_variant",
                "IMPACT=HIGH;LoF=HC;LoF_flags=SINGLE_EXON",
                "v",
            ),
        ],
    );
    let off_pair = [
        "engine\t7:500-502\t-\tENST1\tLoF\tHC",
        "engine\t7:500-502\t-\tENST1\tLoF_flags\tSINGLE_EXON",
    ];
    let on_pair = [
        "engine\t19:13208046-13208050\t-\tENST1\tLoF_flags\tPHYLOCSF_WEAK",
        "reference\t19:13208146-13208147\tTTGA\tENST1\tLoF_flags\tSINGLE_EXON",
    ];
    let expected = [
        (
            "LoF,LoF_info,LoF_filter,LoF_flags",
            [2, 4, 0, 0, 2, 0],
            4,
            off_pair.to_vec(),
        ),
        (
            "LoF,LoF_info,LoF_filter",
            [2, 4, 0, 1, 3, 0],
            2,
            [on_pair, off_pair].concat(),
        ),
    ];
    for (number, (fields, counted, keys, left)) in expected.into_iter().enumerate() {
        let residuals = scratch.path(&format!("{number}/residuals"));
        let report = score(
            &scratch,
            &registry,
            &reference,
            &engine,
            &Run::new()
                .table(&format!("{klass}\tLoFTEE\t{fields}\tmet\t-\n"))
                .residuals(&residuals),
        );
        assert_left(
            &report,
            &residuals,
            klass,
            &[("LoFTEE", counted, keys, left)],
            fields,
        );
    }
}

/// A row of the transcript `ENST1` for an insertion at `location`: its residues and the
/// alternate half of its Codons cell, which is the allele.
fn inserted_codons_row(
    location: &str,
    allele: &str,
    terms: &str,
    residues: &str,
    extra: &str,
) -> String {
    format!(
        "v\t{location}\t{allele}\tENSG1\tENST1\tTranscript\t{terms}\t775-776\t423-424\t141-142\t{residues}\t-/{allele}\t-\t{extra}\n"
    )
}

const SETTLED_INFRAME: [&str; 2] = ["21:10567700-10567701", "GAATAGNNNNCA"];
const SETTLED_STOP_ALONE: [&str; 2] = ["21:10567800-10567801", "TAGNNN"];
const SETTLED_FRAMESHIFT: [&str; 2] = ["21:10567900-10567901", "GAATAGNNNN"];
const STOP_AFTER_UNKNOWN: [&str; 2] = ["21:10568000-10568001", "GAANNNTAGCAT"];
const OTHER_RESIDUES: [&str; 2] = ["21:10568100-10568101", "GAATAGNNNNCA"];

/// Under the committed table's row, an NMD or LOFTEE value the engine alone writes is set
/// aside on a row whose alternate allele settles a stop codon before its unknown bases: the
/// reference reads no peptide and writes `coding_sequence_variant`, the engine `stop_gained`,
/// so the plugins run on the engine's row alone. A value both sides write beside
/// `frameshift_variant` and the same values on a row whose known bases settle no stop
/// codon, or whose residues are not the codons' translation, stay residuals.
#[test]
fn plugin_value_beside_a_stop_codon_settled_before_unknown_bases_is_set_aside() {
    let klass = PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR;
    let table = committed_row(klass);
    assert!(table.contains(
        "; stop_codon_settled_before_unknown_bases, read from the two rows alone with \
         --json-cache and without it;"
    ));
    let scratch = Scratch::new("stop-codon-settled");
    let registry = PlantedRegistry::new(&COMMITTED_ROW_PLUGINS);
    let unknown = |place: [&str; 2]| {
        inserted_codons_row(
            place[0],
            place[1],
            "coding_sequence_variant",
            "-",
            "IMPACT=MODIFIER",
        )
    };
    let truncating = "LoF=HC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.97";
    let gained = |place: [&str; 2], terms: &str, residues: &str| {
        inserted_codons_row(
            place[0],
            place[1],
            terms,
            residues,
            &format!("IMPACT=HIGH;{ESCAPING};{truncating}"),
        )
    };
    let reference = scratch.output(
        "reference.txt",
        &[
            unknown(SETTLED_INFRAME),
            unknown(SETTLED_STOP_ALONE),
            inserted_codons_row(
                SETTLED_FRAMESHIFT[0],
                SETTLED_FRAMESHIFT[1],
                "frameshift_variant",
                "-",
                &format!(
                    "IMPACT=HIGH;DownstreamProtein=MK;ProteinLengthChange=-3;{ESCAPING};LoF=HC;LoF_info=PERCENTILE:0.5"
                ),
            ),
            unknown(STOP_AFTER_UNKNOWN),
            unknown(OTHER_RESIDUES),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            gained(SETTLED_INFRAME, "stop_gained,inframe_insertion", "-/E*XX"),
            gained(SETTLED_STOP_ALONE, "stop_gained", "-/*X"),
            inserted_codons_row(
                SETTLED_FRAMESHIFT[0],
                SETTLED_FRAMESHIFT[1],
                "stop_gained,frameshift_variant",
                "-/E*XX",
                &format!(
                    "IMPACT=HIGH;DownstreamProtein=MK;ProteinLengthChange=-3;{ESCAPING};LoF=HC;LoF_info=PERCENTILE:0.6"
                ),
            ),
            gained(STOP_AFTER_UNKNOWN, "stop_gained,inframe_insertion", "-/EX*H"),
            gained(OTHER_RESIDUES, "stop_gained,inframe_insertion", "-/E*XK"),
        ],
    );
    let residuals = scratch.path("residuals");
    let report = score(
        &scratch,
        &registry,
        &reference,
        &engine,
        &Run::new().table(&table).residuals(&residuals),
    );
    let left = |place: [&str; 2], field: &str, value: &str| {
        format!(
            "engine\t{}\t{}\tENST1\t{field}\t{value}",
            place[0], place[1]
        )
    };
    let off_the_shape = |field: &str, value: &str| {
        [STOP_AFTER_UNKNOWN, OTHER_RESIDUES].map(|place| left(place, field, value))
    };
    let mut loftee = [
        vec![
            left(SETTLED_FRAMESHIFT, "LoF_info", "PERCENTILE:0.6"),
            format!(
                "reference\t{}\t{}\tENST1\tLoF_info\tPERCENTILE:0.5",
                SETTLED_FRAMESHIFT[0], SETTLED_FRAMESHIFT[1]
            ),
        ],
        off_the_shape("LoF", "HC").to_vec(),
        off_the_shape("LoF_filter", "END_TRUNC").to_vec(),
        off_the_shape("LoF_flags", "SINGLE_EXON").to_vec(),
        off_the_shape("LoF_info", "PERCENTILE:0.97").to_vec(),
    ]
    .concat();
    loftee.sort();
    let mut nmd = off_the_shape("NMD", "NMD_escaping_variant").to_vec();
    nmd.sort();
    for (source, counted, keys, expected) in [
        ("Downstream", [2, 2, 2, 2, 2, 2], 0, Vec::new()),
        ("NMD", [1, 5, 1, 1, 3, 1], 2, nmd),
        ("LoFTEE", [2, 18, 1, 2, 10, 1], 8, loftee),
    ] {
        let item = feature(&report, source);
        assert_eq!(counts(item), counted, "{source}");
        assert_eq!(item["masked"], json!({klass: keys}), "{source}");
        let mut written = feature_residuals(&residuals, source);
        written.sort();
        assert_eq!(written, expected, "{source}");
    }
}

const RETAINED_FRAMESHIFT: [&str; 2] = ["21:10568200-10568201", "GAATAGNNNN"];
const RETAINED_INFRAME: [&str; 2] = ["21:10568300-10568301", "GAATAGNNNNCA"];
const RETAINED_OTHER_RESIDUES: [&str; 2] = ["21:10568400-10568401", "GAATAGNNNN"];
const RETAINED_TWO_TERMS: [&str; 2] = ["21:10568500-10568501", "GAATAGNNNN"];

/// Under the committed table's row, a Downstream, NMD or LOFTEE value the engine alone
/// writes is set aside on a row whose alternate allele settles a stop codon before its
/// unknown bases and on which the reference writes `stop_retained_variant` for a stop codon
/// test it does not make: the engine writes `stop_gained` with `frameshift_variant` or an
/// in-frame term, so the plugins run on its row alone. The same values on a row whose
/// residues are not the codons' translation, or whose terms hold two terms beside the
/// gained stop, stay residuals.
#[test]
fn plugin_value_beside_a_gained_stop_the_reference_calls_stop_retained_is_set_aside() {
    let klass = PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR;
    let table = committed_row(klass);
    let scratch = Scratch::new("stop-codon-settled-retained");
    let registry = PlantedRegistry::new(&COMMITTED_ROW_PLUGINS);
    let retained = |place: [&str; 2]| {
        inserted_codons_row(
            place[0],
            place[1],
            "stop_retained_variant",
            "-",
            "IMPACT=LOW",
        )
    };
    let downstream = "DownstreamProtein=MK;ProteinLengthChange=-3";
    let truncating = "LoF=HC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.97";
    let shifted = |place: [&str; 2], terms: &str, residues: &str| {
        inserted_codons_row(
            place[0],
            place[1],
            terms,
            residues,
            &format!("IMPACT=HIGH;{downstream};{ESCAPING};{truncating}"),
        )
    };
    let reference = scratch.output(
        "reference.txt",
        &[
            retained(RETAINED_FRAMESHIFT),
            retained(RETAINED_INFRAME),
            retained(RETAINED_OTHER_RESIDUES),
            retained(RETAINED_TWO_TERMS),
        ],
    );
    let frameshift = "stop_gained,frameshift_variant";
    let engine = scratch.output(
        "engine.txt",
        &[
            shifted(RETAINED_FRAMESHIFT, frameshift, "-/E*XX"),
            inserted_codons_row(
                RETAINED_INFRAME[0],
                RETAINED_INFRAME[1],
                "stop_gained,inframe_insertion",
                "-/E*XX",
                &format!("IMPACT=HIGH;{ESCAPING};{truncating}"),
            ),
            shifted(RETAINED_OTHER_RESIDUES, frameshift, "-/E*XK"),
            shifted(
                RETAINED_TWO_TERMS,
                "stop_gained,frameshift_variant,inframe_insertion",
                "-/E*XX",
            ),
        ],
    );
    let residuals = scratch.path("residuals");
    let report = score(
        &scratch,
        &registry,
        &reference,
        &engine,
        &Run::new().table(&table).residuals(&residuals),
    );
    let off_the_shape = |fields: &[(&str, &str)]| {
        let mut left: Vec<String> = Vec::new();
        for place in [RETAINED_OTHER_RESIDUES, RETAINED_TWO_TERMS] {
            for (field, value) in fields {
                left.push(format!(
                    "engine\t{}\t{}\tENST1\t{field}\t{value}",
                    place[0], place[1]
                ));
            }
        }
        left.sort();
        left
    };
    for (source, counted, keys, expected) in [
        (
            "Downstream",
            [0, 6, 0, 0, 4, 0],
            2,
            off_the_shape(&[("DownstreamProtein", "MK"), ("ProteinLengthChange", "-3")]),
        ),
        (
            "NMD",
            [0, 4, 0, 0, 2, 0],
            2,
            off_the_shape(&[("NMD", "NMD_escaping_variant")]),
        ),
        (
            "LoFTEE",
            [0, 16, 0, 0, 8, 0],
            8,
            off_the_shape(&[
                ("LoF", "HC"),
                ("LoF_filter", "END_TRUNC"),
                ("LoF_flags", "SINGLE_EXON"),
                ("LoF_info", "PERCENTILE:0.97"),
            ]),
        ),
    ] {
        let item = feature(&report, source);
        assert_eq!(counts(item), counted, "{source}");
        assert_eq!(item["masked"], json!({klass: keys}), "{source}");
        let mut written = feature_residuals(&residuals, source);
        written.sort();
        assert_eq!(written, expected, "{source}");
    }
}

/// A made-up record of a transcript of one exon whose peptide holds a stop at residue 10 and
/// reads on past it, in the cells both engines write: the record's name, its Location and
/// allele, the three position columns, the residues, the Codons cell, the HGVSc after the
/// transcript and the HGVS_OFFSET pair of a row that writes one.
type HeldStopRecord = [&'static str; 10];

/// The cells of one engine's row of such a record: the Consequence, the IMPACT, the HGVSp
/// after the protein and the plugins' fields.
type HeldStopSide = [&'static str; 4];

/// A record with the reference's row and the engine's.
type HeldStopPair = (HeldStopRecord, HeldStopSide, HeldStopSide);

const HELD_STOP_SHIFTED: &str = "frameshift_variant";
const HELD_STOP_OFFSET: &str = "HGVS_OFFSET=1";

/// The seven made-up records. The first is the insertion written one base before the stop
/// codon and the last a base lost from the closing stop codon: the controls, whose two rows
/// are equal. The five between them are rows of `frameshift_inside_a_held_stop_codon`.
const HELD_STOP_PAIRS: [HeldStopPair; 7] = [
    (
        [
            "written_before_the_stop_codon",
            "1:133-134",
            "T",
            "33-34",
            "27-28",
            "9-10",
            "-/X",
            "-/T",
            "c.28dup",
            HELD_STOP_OFFSET,
        ],
        [
            HELD_STOP_SHIFTED,
            "HIGH",
            "p.Ter10LeuextTer16",
            "DownstreamProtein=LSCCCCCCCCCYPKLA;NMD=NMD_escaping_variant;ProteinLengthChange=2",
        ],
        [
            HELD_STOP_SHIFTED,
            "HIGH",
            "p.Ter10LeuextTer16",
            "DownstreamProtein=LSCCCCCCCCCYPKLA;NMD=NMD_escaping_variant;ProteinLengthChange=2",
        ],
    ),
    (
        [
            "written_inside_the_stop_codon",
            "1:134-135",
            "T",
            "34-35",
            "28-29",
            "10",
            "*/LX",
            "tga/tTga",
            "c.28dup",
            "",
        ],
        [
            "coding_sequence_variant",
            "MODIFIER",
            "p.Ter10delinsLeuTer",
            "",
        ],
        [
            HELD_STOP_SHIFTED,
            "HIGH",
            "p.Ter10LeuextTer16",
            "DownstreamProtein=LSCCCCCCCCCYPKLA;NMD=NMD_escaping_variant;ProteinLengthChange=2",
        ],
    ),
    (
        [
            "a_base_lost_inside_the_stop_codon",
            "1:136",
            "-",
            "36",
            "30",
            "10",
            "*/X",
            "tgA/tg",
            "c.30del",
            "",
        ],
        ["coding_sequence_variant", "MODIFIER", "p.Ter10%3D", ""],
        [
            HELD_STOP_SHIFTED,
            "HIGH",
            "p.Ter10TrpextTer12",
            "DownstreamProtein=WLLLLLLLLLLP;NMD=NMD_escaping_variant;ProteinLengthChange=-2",
        ],
    ),
    (
        [
            "five_bases_lost_from_inside_the_stop_codon",
            "1:135-139",
            "-",
            "35-39",
            "29-33",
            "10-11",
            "*A/X",
            "tGAGCT/t",
            "c.30_34del",
            HELD_STOP_OFFSET,
        ],
        [
            "inframe_deletion",
            "MODERATE",
            "p.Ter10_Ala12delinsCysTer",
            "",
        ],
        [
            HELD_STOP_SHIFTED,
            "HIGH",
            "p.Ter10CysextTer14",
            "DownstreamProtein=CCCCCCCCCYPKLA;NMD=NMD_escaping_variant;ProteinLengthChange=0",
        ],
    ),
    (
        [
            "four_bases_added_and_the_stop_kept",
            "1:135-136",
            "AATG",
            "35-36",
            "29-30",
            "10",
            "*/*MX",
            "tga/tgAATGa",
            "c.30_31insATGA",
            HELD_STOP_OFFSET,
        ],
        [
            "inframe_insertion,coding_sequence_variant",
            "MODERATE",
            "p.Ter10_Ala11insMetTer",
            "",
        ],
        [
            HELD_STOP_SHIFTED,
            "HIGH",
            "p.Ala11MetfsTer17",
            "DownstreamProtein=;NMD=NMD_escaping_variant;ProteinLengthChange=-14",
        ],
    ),
    (
        [
            "a_base_lost_and_the_stop_kept",
            "1:135",
            "-",
            "35",
            "29",
            "10",
            "*/X",
            "tGa/ta",
            "c.29del",
            "",
        ],
        ["coding_sequence_variant", "MODIFIER", "p.Ter10%3D", ""],
        [
            HELD_STOP_SHIFTED,
            "HIGH",
            "p.Ala11LeufsTer12",
            "DownstreamProtein=;NMD=NMD_escaping_variant;ProteinLengthChange=-14",
        ],
    ),
    (
        [
            "a_base_lost_inside_the_closing_stop_codon",
            "1:178",
            "-",
            "78",
            "72",
            "24",
            "*/X",
            "taG/ta",
            "c.72del",
            "",
        ],
        ["stop_lost", "HIGH", "p.Ter24%3D", ""],
        ["stop_lost", "HIGH", "p.Ter24%3D", ""],
    ),
];

/// The made-up transcript's row of `record` as `side` writes it, with the pairs of
/// `further` at the end of its Extra column.
fn held_stop_row(record: HeldStopRecord, side: HeldStopSide, further: &str) -> String {
    let [name, location, allele, cdna, cds, protein, residues, codons, hgvsc, offset] = record;
    let [terms, impact, hgvsp, plugins] = side;
    let mut extra = format!(
        "IMPACT={impact};STRAND=1;HGVSc=ENST00000000001.1:{hgvsc};HGVSp=ENSP00000000001.1:{hgvsp}"
    );
    for pairs in [offset, plugins, further] {
        if !pairs.is_empty() {
            extra.push(';');
            extra.push_str(pairs);
        }
    }
    format!(
        "{name}\t{location}\t{allele}\tENSG00000000001\tENST00000000001\tTranscript\t{terms}\t{cdna}\t{cds}\t{protein}\t{residues}\t{codons}\t-\t{extra}\n"
    )
}

/// The residual lines of the populated fields of `fields` that the engine's rows of `pairs`
/// write, sorted.
fn held_stop_engine_lines(pairs: &[HeldStopPair], fields: &[&str]) -> Vec<String> {
    let mut lines = Vec::new();
    for (record, _, side) in pairs {
        for (field, value) in side[3].split(';').filter_map(|pair| pair.split_once('=')) {
            if fields.contains(&field) && !value.is_empty() {
                lines.push(format!(
                    "engine\t{}\t{}\tENST00000000001\t{field}\t{value}",
                    record[1], record[2]
                ));
            }
        }
    }
    lines.sort();
    lines
}

/// Under the committed table's row, every Downstream and NMD value the engine alone writes
/// on the five made-up rows of `frameshift_inside_a_held_stop_codon` is set aside, and both
/// projections score exactly: the reference writes `coding_sequence_variant`, an in-frame
/// term or both where the engine writes `frameshift_variant`, so the plugins run on the
/// engine's row alone. Under a table that declares no class the same values are residuals.
/// The two controls pair as text under either table.
#[test]
fn plugin_values_beside_a_frameshift_inside_a_held_stop_codon_are_set_aside() {
    let klass = PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR;
    let table = committed_row(klass);
    assert!(table.contains(
        "; frameshift_inside_a_held_stop_codon, read from the two rows alone with \
         --json-cache and without it;"
    ));
    let scratch = Scratch::new("held-stop");
    let registry = PlantedRegistry::new(&COMMITTED_ROW_PLUGINS);
    let reference = scratch.output(
        "reference.txt",
        &HELD_STOP_PAIRS.map(|(record, side, _)| held_stop_row(record, side, "")),
    );
    let engine = scratch.output(
        "engine.txt",
        &HELD_STOP_PAIRS.map(|(record, _, side)| held_stop_row(record, side, "")),
    );
    let scored = |name: &str, table: &str| {
        let residuals = scratch.path(&format!("{name}/residuals"));
        let run = Run::new().table(table).residuals(&residuals);
        let report = score(&scratch, &registry, &reference, &engine, &run);
        (report, residuals)
    };

    let (report, residuals) = scored("with", &table);
    for (source, counted, keys) in [
        ("Downstream", [2, 10, 2, 2, 2, 2], 8),
        ("NMD", [1, 6, 1, 1, 1, 1], 5),
    ] {
        let item = feature(&report, source);
        assert_eq!(counts(item), counted, "{source}");
        assert_eq!(item["masked"], json!({klass: keys}), "{source}");
        assert_eq!(item["residual_rows"], 0, "{source}");
        assert!(feature_residuals(&residuals, source).is_empty(), "{source}");
        assert_eq!(counts(allele(&report, source)), counted, "{source}");
    }

    let (report, residuals) = scored("without", "#class\tplugin\tfield names\n");
    let class_pairs = &HELD_STOP_PAIRS[1..6];
    for (plugin, counted) in [(0, [2, 10, 2, 2, 10, 2]), (1, [1, 6, 1, 1, 6, 1])] {
        let (source, fields) = TERM_FOLLOWING_PLUGINS[plugin];
        let item = feature(&report, source);
        assert_eq!(counts(item), counted, "{source}");
        assert_eq!(item["masked"], json!({}), "{source}");
        let left = held_stop_engine_lines(class_pairs, fields);
        assert_eq!(left.len() as i64, counted[1] - counted[2], "{source}");
        assert_eq!(item["residual_rows"], left.len(), "{source}");
        let mut written = feature_residuals(&residuals, source);
        written.sort();
        assert_eq!(written, left, "{source}");
    }
}

const LOFTEE_VALUES: &str =
    "LoF=HC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.97";

/// Under the committed table's row, the four LOFTEE fields the engine alone writes on the
/// five made-up rows of `frameshift_inside_a_held_stop_codon` are set aside with the
/// Downstream and NMD values. A further plugin value that differs on such a row stays a
/// residual: a field of a plugin the row of the table does not name, and a
/// `ProteinLengthChange` the reference writes as well, with another value.
#[test]
fn further_plugin_value_beside_a_frameshift_inside_a_held_stop_codon_stays_a_residual() {
    let klass = PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR;
    let scratch = Scratch::new("held-stop-further");
    let mut plugins = COMMITTED_ROW_PLUGINS.to_vec();
    plugins.push(("Planted", &["PLANT_label"]));
    let registry = PlantedRegistry::new(&plugins);
    let class_pairs = &HELD_STOP_PAIRS[1..6];
    let [inside, lost] = [class_pairs[0].0, class_pairs[1].0];
    let reference: Vec<String> = class_pairs
        .iter()
        .map(|(record, side, _)| {
            let further = if *record == lost {
                "ProteinLengthChange=5"
            } else {
                ""
            };
            held_stop_row(*record, *side, further)
        })
        .collect();
    let engine: Vec<String> = class_pairs
        .iter()
        .map(|(record, _, side)| {
            let further = if *record == inside {
                format!("{LOFTEE_VALUES};PLANT_label=benign")
            } else {
                LOFTEE_VALUES.to_string()
            };
            held_stop_row(*record, *side, &further)
        })
        .collect();
    let residuals = scratch.path("residuals");
    let report = score(
        &scratch,
        &registry,
        &scratch.output("reference.txt", &reference),
        &scratch.output("engine.txt", &engine),
        &Run::new()
            .table(&committed_row(klass))
            .residuals(&residuals),
    );
    let line = |side: &str, record: HeldStopRecord, field: &str, value: &str| {
        format!(
            "{side}\t{}\t{}\tENST00000000001\t{field}\t{value}",
            record[1], record[2]
        )
    };
    let written_twice = [
        line("engine", lost, "ProteinLengthChange", "-2"),
        line("reference", lost, "ProteinLengthChange", "5"),
    ];
    assert_left(
        &report,
        &residuals,
        klass,
        &[
            (
                "Downstream",
                [1, 8, 0, 1, 1, 0],
                7,
                written_twice.iter().map(String::as_str).collect(),
            ),
            ("NMD", [0, 5, 0, 0, 0, 0], 5, Vec::new()),
            ("LoFTEE", [0, 20, 0, 0, 0, 0], 20, Vec::new()),
        ],
        "the table's row",
    );
    let planted = feature(&report, "Planted");
    assert_eq!(counts(planted), [0, 1, 0, 0, 1, 0]);
    assert_eq!(
        feature_residuals(&residuals, "Planted"),
        [line("engine", inside, "PLANT_label", "benign")]
    );
}

const SPLICE_SPAN_FIELDS: &str =
    "DownstreamProtein,ProteinLengthChange;NMD;LoF,LoF_info,LoF_filter,LoF_flags";
const IN_FRAME: [&str; 2] = [
    "splice_donor_variant,protein_altering_variant",
    "protein_altering_variant,splice_region_variant",
];
const FRAME_SHIFTED: [&str; 2] = [
    "splice_acceptor_variant,frameshift_variant",
    "frameshift_variant",
];
const TAKEN: Place = ["7:100998-101000", "ACGTAC", "998-1000", "898-900"];
const TAKEN_SHIFTED: Place = ["7:200998-201000", "ACGTA", "498-500", "398-400"];
const TWO_NEIGHBOURS: Place = ["7:300999-301000", "ACGTA", "999-1000", "899-900"];
const ACROSS: Place = ["7:400998-401000", "ACGTAC", "999-?", "899-?"];
const BOTH_WRITE: Place = ["7:500998-501000", "ACGTAC", "298-300", "198-200"];
const NAMED: Place = ["7:600999-601000", "CGTAC", "999-1000", "899-900"];
const SPELT: Place = ["7:700999-701000", "CGTAC", "999-1000", "899-900"];
const SPELT_ALLELES: &str = "7_700998_ACG/ACGTAC";

/// A made-up reference and engine output of seven rows, and the input VCF of the run.
/// Three are rows of a splice term past the record span read without their record, three
/// positions in the transcript's sequence under a longer allele and the identifier `v`,
/// which no line of the input carries: one on which the reference alone holds the donor
/// term, its row's only truncating term, and writes NMD and four LOFTEE fields; one on
/// which it holds the acceptor term beside frameshift_variant, and the engine, whose row
/// holds no splice term, alone writes the two Downstream fields; and one on which both
/// sides write a `LoF_info`, with two values. The values of the first two stand as well,
/// under `v`, on a row of two neighbouring positions and on a row whose span crosses an
/// exon's end. Two more rows of two positions inside an exon are the class's with their
/// record alone, the reference writing NMD and `LoF` beside its donor term on each: one
/// under the identifier `rs90`, whose record the input holds, and one whose first column
/// spells the alleles.
fn splice_span_outputs(scratch: &Scratch) -> (String, String, String) {
    let escaped = format!("IMPACT=HIGH;{ESCAPING};LoF=HC");
    let reference = scratch.output(
        "reference.txt",
        &[
            coding_row(
                TAKEN,
                IN_FRAME[0],
                &format!(
                    "IMPACT=HIGH;{ESCAPING};LoF=HC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.97"
                ),
                "v",
            ),
            coding_row(TAKEN_SHIFTED, FRAME_SHIFTED[0], &escaped, "v"),
            coding_row(TWO_NEIGHBOURS, FRAME_SHIFTED[0], &escaped, "v"),
            coding_row(ACROSS, IN_FRAME[0], &escaped, "v"),
            coding_row(
                BOTH_WRITE,
                IN_FRAME[0],
                "IMPACT=HIGH;LoF=HC;LoF_info=PERCENTILE:0.9",
                "v",
            ),
            coding_row(NAMED, IN_FRAME[0], &escaped, "rs90"),
            coding_row(SPELT, IN_FRAME[0], &escaped, SPELT_ALLELES),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            coding_row(TAKEN, IN_FRAME[1], "IMPACT=MODERATE", "v"),
            coding_row(
                TAKEN_SHIFTED,
                FRAME_SHIFTED[1],
                &format!(
                    "IMPACT=HIGH;DownstreamProtein=MKV;ProteinLengthChange=-12;{ESCAPING};LoF=HC"
                ),
                "v",
            ),
            coding_row(
                TWO_NEIGHBOURS,
                FRAME_SHIFTED[1],
                &format!(
                    "IMPACT=HIGH;DownstreamProtein=MKL;ProteinLengthChange=-3;{ESCAPING};LoF=HC"
                ),
                "v",
            ),
            coding_row(ACROSS, IN_FRAME[1], "IMPACT=MODERATE", "v"),
            coding_row(
                BOTH_WRITE,
                IN_FRAME[1],
                "IMPACT=MODERATE;LoF_info=PERCENTILE:0.5",
                "v",
            ),
            coding_row(NAMED, IN_FRAME[1], "IMPACT=MODERATE", "rs90"),
            coding_row(SPELT, IN_FRAME[1], "IMPACT=MODERATE", SPELT_ALLELES),
        ],
    );
    let run_input = scratch.vcf(
        "input.vcf",
        &[
            ["7", "600998", "rs90", "ACG", "ACGTAC"],
            ["7", "700998", ".", "ACG", "ACGTAC"],
        ],
    );
    (reference, engine, run_input)
}

/// On a row of a splice term past the record span an NMD or LOFTEE value the reference
/// alone writes beside its donor term and a Downstream value the engine alone writes
/// beside frameshift_variant and no splice term are set aside. The class reads a row with
/// its record where the row's first column spells the alleles and, under an identifier,
/// where the lines of the input that carry the identifier give the row by position; scored
/// without the input, the row under the identifier is read without its record and its
/// values stay residuals. The same values on a row of two neighbouring positions and on a
/// row whose span crosses an exon's end, under an identifier no input line carries, and a
/// `LoF_info` both sides write with two values stay residuals either way.
#[test]
fn plugin_value_that_follows_a_splice_term_past_the_record_span_is_set_aside() {
    let klass = PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN;
    assert_eq!(
        klass,
        "plugin_value_follows_splice_term_past_the_record_span"
    );
    assert!(BuiltIn.record_reference(klass));
    let scratch = Scratch::new("splice-span");
    let registry = PlantedRegistry::new(&TERM_FOLLOWING_PLUGINS);
    let (reference, engine, run_input) = splice_span_outputs(&scratch);
    let table = format!("{klass}\tDownstream;NMD;LoFTEE\t{SPLICE_SPAN_FIELDS}\tmet\t-\n");
    let downstream: Left = (
        "Downstream",
        [0, 4, 0, 0, 2, 0],
        2,
        vec![
            "engine\t7:300999-301000\tACGTA\tENST1\tDownstreamProtein\tMKL",
            "engine\t7:300999-301000\tACGTA\tENST1\tProteinLengthChange\t-3",
        ],
    );
    let nmd = vec!["reference\t7:400998-401000\tACGTAC\tENST1\tNMD\tNMD_escaping_variant"];
    let loftee = vec![
        "reference\t7:400998-401000\tACGTAC\tENST1\tLoF\tHC",
        "engine\t7:500998-501000\tACGTAC\tENST1\tLoF_info\tPERCENTILE:0.5",
        "reference\t7:500998-501000\tACGTAC\tENST1\tLoF_info\tPERCENTILE:0.9",
    ];
    // Scored without the input, the row under the identifier `rs90` is read without its
    // record: its NMD and LOFTEE values stay residuals.
    let unread_nmd = "reference\t7:600999-601000\tCGTAC\tENST1\tNMD\tNMD_escaping_variant";
    let unread_loftee = "reference\t7:600999-601000\tCGTAC\tENST1\tLoF\tHC";
    let cases: [(&str, Option<&str>, [Left; 3]); 2] = [
        (
            "given",
            Some(&run_input),
            [
                downstream.clone(),
                ("NMD", [6, 2, 2, 3, 2, 2], 3, nmd.clone()),
                ("LoFTEE", [11, 3, 2, 4, 3, 2], 7, loftee.clone()),
            ],
        ),
        (
            "without",
            None,
            [
                downstream,
                (
                    "NMD",
                    [6, 2, 2, 4, 2, 2],
                    2,
                    [nmd, vec![unread_nmd]].concat(),
                ),
                (
                    "LoFTEE",
                    [11, 3, 2, 5, 3, 2],
                    6,
                    [loftee, vec![unread_loftee]].concat(),
                ),
            ],
        ),
    ];
    for (name, given, scores) in cases {
        let residuals = scratch.path(&format!("residuals-{name}"));
        let mut run = Run::new().table(&table).residuals(&residuals);
        if let Some(input) = given {
            run = run.run_input(input);
        }
        let report = score(&scratch, &registry, &reference, &engine, &run);
        assert_left(&report, &residuals, klass, &scores, name);
    }
}

const MITOCHONDRIAL_SPELLING: &str = "plugin_contig_mitochondrial_spelling";
const DBNSFP_NAMED_M: &str = "dbnsfp_mitochondrial_record_named_m";

/// The row of the record `<contig>` 3460 G A on one transcript under an identifier, its
/// Extra text ending with `fields` when there are any.
fn mitochondrial(location: &str, fields: &str) -> String {
    let extra = if fields.is_empty() {
        "IMPACT=MODERATE".to_string()
    } else {
        format!("IMPACT=MODERATE;{fields}")
    };
    row_with(location, "A", "ENST00000361390", &extra, "-", "rs199476118")
}

/// The report of `reference` against `engine` under the tool's own registry and classes, in
/// a run that loaded `plugins`, read `table` and gave `words` as its command line.
fn mitochondrial_report(
    scratch: &Scratch,
    table: &str,
    plugins: &[&str],
    words: Option<&[&str]>,
    reference: &str,
    engine: &str,
) -> Value {
    let mut run = Run::new().table(table).plugins(plugins);
    if let Some(words) = words {
        run = run.run_command(words);
    }
    try_score(
        scratch,
        reference,
        &[("engine", engine)],
        &run,
        &BuiltInRegistry,
        &BuiltIn,
    )
    .unwrap_or_else(|stop| panic!("the run was stopped: {}", message(&stop)))
}

/// Under the tree's own row of the class, a field of a tabix-backed plugin that the engine
/// alone fills on a record spelt MT, chrM or chrMT is set aside whatever assembly the run's
/// command line names: the scorer reads the record's spelling and which output fills the
/// field, never the plugin's data file. A gnomADMt field, which the row gives no class, a
/// record spelt M and a field the reference alone fills stay residual rows.
#[test]
fn committed_mitochondrial_spelling_row_takes_a_tabix_plugins_field_the_engine_alone_fills() {
    let scratch = Scratch::new("mitochondrial-spelling-row");
    let table = committed_row(MITOCHONDRIAL_SPELLING);
    let plugins = ["ClinPred", "SpliceAI", "gnomADc", "gnomADMt"];
    let filled = "ClinPred=0.5;SpliceAI_pred=GENE1|0.10|0.20|0.30|0.40|-4|-3|2|1;gnomAD_mean_cov=31.2;gnomAD_hap_AF_hom=0.01";
    let scored = |name: &str, location: &str, fields: &str| {
        scratch.output(name, &[mitochondrial(location, fields)])
    };
    let masked = |report: &Value, source: &str| feature(report, source)["masked"].clone();
    let residual = |report: &Value, source: &str| feature(report, source)["residual_rows"].clone();
    for (run, words) in [
        ("grch37", Some(&["--offline", "--assembly", "GRCh37"][..])),
        ("grch38", Some(&["--offline", "--assembly", "GRCh38"][..])),
        ("unnamed", None),
    ] {
        for contig in ["MT", "chrM", "chrMT"] {
            let location = format!("{contig}:3460");
            let bare = scored(&format!("{run}/{contig}/bare.txt"), &location, "");
            let full = scored(&format!("{run}/{contig}/full.txt"), &location, filled);
            let report = mitochondrial_report(&scratch, &table, &plugins, words, &bare, &full);
            let left: Vec<&str> = ["ClinPred", "SpliceAI", "gnomADc"]
                .into_iter()
                .filter(|source| {
                    masked(&report, source) != json!({MITOCHONDRIAL_SPELLING: 1})
                        || residual(&report, source) != 0
                })
                .collect();
            assert_eq!(left, Vec::<&str>::new(), "{run} {contig}");
            assert_eq!(masked(&report, "gnomADMt"), json!({}), "{run} {contig}");
            assert_eq!(residual(&report, "gnomADMt"), 1, "{run} {contig}");
            // The same fields the reference alone fills.
            let report = mitochondrial_report(&scratch, &table, &plugins, words, &full, &bare);
            for source in ["ClinPred", "SpliceAI", "gnomADc"] {
                let what = format!("{run} {contig} {source}, the reference alone");
                assert_eq!(
                    masked(&report, source),
                    json!({MITOCHONDRIAL_SPELLING: 0}),
                    "{what}"
                );
                assert_eq!(residual(&report, source), 1, "{what}");
            }
        }
        let bare = scored(&format!("{run}/M/bare.txt"), "M:3460", "");
        let full = scored(&format!("{run}/M/full.txt"), "M:3460", filled);
        let report = mitochondrial_report(&scratch, &table, &plugins, words, &bare, &full);
        for source in ["ClinPred", "SpliceAI", "gnomADc"] {
            assert_eq!(residual(&report, source), 1, "{run} M {source}");
        }
    }
}

/// A pair the dbNSFP class leaves: what it shows, the run's command line or none, the
/// record's Location, the reference's fields, the engine's, and the residual rows of dbNSFP.
type StayingRow<'a> = (
    &'a str,
    Option<&'a [&'a str]>,
    &'a str,
    &'a str,
    &'a str,
    i64,
);

/// Under the tree's own row of the class, a dbNSFP field the reference alone fills on a
/// record spelt MT or chrMT is set aside in a run whose command line names GRCh37:
/// dbNSFP.pm:313 names such a record M to the table and :335-338 match the row at its
/// position, so against a made-up table naming the contig M whose rows count on another
/// sequence the reference writes the row's value and the engine nothing. The field is
/// dbNSFP's by the registry and the run's plugins, under its own name or the plugin's
/// prefix. It stays a residual row in a run that names another assembly or none, on a
/// record spelt chrM or M and on another chromosome, where both outputs fill it or the
/// engine alone does, and where it is another plugin's.
#[test]
fn committed_dbnsfp_named_m_row_takes_a_dbnsfp_field_the_reference_alone_fills_on_grch37() {
    let scratch = Scratch::new("dbnsfp-named-m-row");
    let table = committed_row(DBNSFP_NAMED_M);
    let plugins = ["dbNSFP", "ClinPred"];
    let table_word = "dbNSFP,/data/dbNSFP4.9a_grch37_synthetic.gz,SIFT_score,Polyphen2_HDIV_score";
    let scored = |name: &str, location: &str, fields: &str| {
        scratch.output(name, &[mitochondrial(location, fields)])
    };
    let masked = |report: &Value, source: &str| feature(report, source)["masked"].clone();
    let residual = |report: &Value, source: &str| feature(report, source)["residual_rows"].clone();
    let filled = "SIFT_score=0.01;Polyphen2_HDIV_score=0.9;ClinPred=0.5";
    let grch37: [&[&str]; 4] = [
        &["--offline", "--assembly", "GRCh37", "--plugin", table_word],
        &["--assembly=hg19", "--plugin", table_word],
        &["-a", "37", "--plugin", table_word],
        &["--assembly", "GRCh38", "--grch37", "--plugin", table_word],
    ];
    for (number, words) in grch37.iter().enumerate() {
        for contig in ["MT", "chrMT"] {
            let location = format!("{contig}:3460");
            let full = scored(
                &format!("taken{number}/{contig}/full.txt"),
                &location,
                filled,
            );
            let bare = scored(&format!("taken{number}/{contig}/bare.txt"), &location, "");
            let report =
                mitochondrial_report(&scratch, &table, &plugins, Some(words), &full, &bare);
            let what = format!("{words:?} {contig}");
            // One count a field: the row's two dbNSFP fields.
            assert_eq!(
                masked(&report, "dbNSFP"),
                json!({DBNSFP_NAMED_M: 2}),
                "{what}"
            );
            assert_eq!(residual(&report, "dbNSFP"), 0, "{what}");
            assert_eq!(masked(&report, "ClinPred"), json!({}), "{what}");
            assert_eq!(residual(&report, "ClinPred"), 1, "{what}");
        }
    }
    let prefixed = scored("prefixed/full.txt", "MT:3460", "dbNSFP_SIFT_score=0.01");
    let bare = scored("prefixed/bare.txt", "MT:3460", "");
    let report = mitochondrial_report(
        &scratch,
        &table,
        &plugins,
        Some(grch37[0]),
        &prefixed,
        &bare,
    );
    assert_eq!(masked(&report, "dbNSFP"), json!({DBNSFP_NAMED_M: 1}));
    assert_eq!(residual(&report, "dbNSFP"), 0);
    let left: [StayingRow; 9] = [
        (
            "GRCh38",
            Some(&["--assembly", "GRCh38", "--plugin", table_word]),
            "MT:3460",
            filled,
            "",
            2,
        ),
        ("no command line", None, "MT:3460", filled, "", 2),
        (
            "a command line that names no assembly",
            Some(&["--offline", "--plugin", table_word]),
            "chrMT:3460",
            filled,
            "",
            2,
        ),
        ("chrM", Some(grch37[0]), "chrM:3460", filled, "", 2),
        ("M", Some(grch37[0]), "M:3460", filled, "", 2),
        ("21", Some(grch37[0]), "21:3460", filled, "", 2),
        (
            "both fill",
            Some(grch37[0]),
            "MT:3460",
            "SIFT_score=0.01",
            "SIFT_score=0.02",
            2,
        ),
        (
            "the engine alone",
            Some(grch37[0]),
            "MT:3460",
            "",
            filled,
            2,
        ),
        (
            "an empty value of the reference",
            Some(grch37[0]),
            "MT:3460",
            "SIFT_score=",
            "SIFT_score=0.01",
            1,
        ),
    ];
    for (number, (what, words, location, reference, engine, residuals)) in left.iter().enumerate() {
        let reference = scored(&format!("left{number}/reference.txt"), location, reference);
        let engine = scored(&format!("left{number}/engine.txt"), location, engine);
        let report = mitochondrial_report(&scratch, &table, &plugins, *words, &reference, &engine);
        assert_eq!(
            masked(&report, "dbNSFP"),
            json!({DBNSFP_NAMED_M: 0}),
            "{what}"
        );
        assert_eq!(residual(&report, "dbNSFP"), *residuals, "{what}");
    }
}

const WALK_STOP_GAIN: Place = ["21:37442543", "T", "305", "130"];
const WALK_FRAMESHIFT: Place = ["21:31588081", "-", "734", "246"];
const WALK_DONOR: Place = ["21:37442703", "T", "-", "-"];
const WALK_CALLED: Place = ["21:37443313", "T", "535", "360"];
const WALKED: &str = "IMPACT=HIGH;LoF=HC;LoF_info=GERP_DIST:705,BP_DIST:702";
const WALKED_FILTERED: &str =
    "IMPACT=HIGH;LoF=LC;LoF_filter=END_TRUNC;LoF_info=GERP_DIST:150,BP_DIST:148";

/// Under the committed table's own row, which reads `not met`, the LOFTEE fields the engine
/// alone writes on a stop-gained or frameshift row whose reference row carries none of them
/// are set aside in a run that names the class and listed as not applied in a run that
/// names none. The same fields on a splice-site row and a `LoF_info` beside a `LoF` the
/// reference writes stay residuals either way. A tuple file of the same values, which
/// states no terms, is read by the reference's rows and leaves the same counts.
#[test]
fn loftee_fields_written_after_a_gerp_walk_are_set_aside_where_the_run_names_the_class() {
    let klass = LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE;
    assert_eq!(klass, "loftee_gerp_walk_without_conservation_database");
    let table = committed_row(klass);
    let cells: Vec<&str> = table.trim_end().split('\t').collect();
    assert_eq!(
        cells[1..4],
        ["LoFTEE", "LoF,LoF_info,LoF_filter,LoF_flags", "not met"]
    );
    let scratch = Scratch::new("gerp-walk");
    let registry = PlantedRegistry::new(&TERM_FOLLOWING_PLUGINS[2..]);
    let (stop_gain, frameshift, donor) =
        ("stop_gained", "frameshift_variant", "splice_donor_variant");
    let reference = scratch.output(
        "reference.txt",
        &[
            coding_row(WALK_STOP_GAIN, stop_gain, "IMPACT=HIGH", "v"),
            coding_row(WALK_FRAMESHIFT, frameshift, "IMPACT=HIGH", "v"),
            coding_row(WALK_DONOR, donor, "IMPACT=HIGH", "v"),
            coding_row(WALK_CALLED, stop_gain, "IMPACT=HIGH;LoF=HC", "v"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            coding_row(WALK_STOP_GAIN, stop_gain, WALKED, "v"),
            coding_row(WALK_FRAMESHIFT, frameshift, WALKED_FILTERED, "v"),
            coding_row(WALK_DONOR, donor, WALKED, "v"),
            coding_row(WALK_CALLED, stop_gain, WALKED, "v"),
        ],
    );
    let tuples: Vec<[&str; 5]> = [
        (WALK_STOP_GAIN, WALKED),
        (WALK_FRAMESHIFT, WALKED_FILTERED),
        (WALK_DONOR, WALKED),
        (WALK_CALLED, WALKED),
    ]
    .iter()
    .flat_map(|(place, extra)| {
        extra.split(';').skip(1).map(move |pair| {
            let (field, value) = pair.split_once('=').unwrap();
            [place[0], place[1], "ENST1", field, value]
        })
    })
    .collect();
    let tuple_file = scratch.tuples("engine.tsv", &tuples);
    let left = vec![
        "engine\t21:37442703\tT\tENST1\tLoF\tHC",
        "engine\t21:37442703\tT\tENST1\tLoF_info\tGERP_DIST:705,BP_DIST:702",
        "engine\t21:37443313\tT\tENST1\tLoF_info\tGERP_DIST:705,BP_DIST:702",
    ];
    for (name, written) in [("rows", &engine), ("tuples", &tuple_file)] {
        let residuals = scratch.path(&format!("residuals-{name}"));
        let named = score(
            &scratch,
            &registry,
            &reference,
            written,
            &Run::new()
                .table(&table)
                .met_alone()
                .applied(&[klass])
                .residuals(&residuals),
        );
        assert_left(
            &named,
            &residuals,
            klass,
            &[("LoFTEE", [1, 9, 1, 1, 4, 1], 5, left.clone())],
            name,
        );
    }
    let unnamed = score(
        &scratch,
        &registry,
        &reference,
        &engine,
        &Run::new().table(&table).met_alone(),
    );
    let item = feature(&unnamed, "LoFTEE");
    assert_eq!(raw_counts(item), [1, 9, 1]);
    assert_eq!(item["masked"], json!({}));
    assert_eq!(item["not_applied"], json!({klass: 5}));
    assert_eq!(item["residual_rows"], 8);
}
