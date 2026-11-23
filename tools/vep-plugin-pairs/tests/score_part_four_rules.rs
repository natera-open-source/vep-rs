// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `score` program under the tool's own classes on a planted registry: the
//! sibling-allele classes and the engine's rows they are handed, the REVEL table, the
//! plugin values behind a splice term past the record span, the MaxEntScan windows of the
//! candidate cache, and the committed class table on the MaxEntScan fixture corpora.

mod support;

use flate2::read::MultiGzDecoder;
use flate2::write::DeflateEncoder;
use flate2::{Compression, Crc};
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use support::{
    allele, counts, entry, feature, feature_residuals, gzip, message, on_any_thread_count, row,
    run_main, sha256, try_score, PlantedRegistry, Run, Scratch, ERROR,
};
use vep_plugin_pairs::registry::BuiltInRegistry;
use vep_plugin_pairs::rules::revel::RevelTable;
use vep_plugin_pairs::rules::{
    parse_default, BuiltIn, Consult, RecordRows, Row, Rules, ABSENT,
    ALPHAMISSENSE_SIBLING_ALLELE_ROWS, MAXENTSCAN_OVERLAPPING_WINDOWS,
    PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN, REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS,
    REVEL_SIBLING_ALLELE_ROWS, SPLICEVAULT_SIBLING_ALLELE_ROWS,
};

/// The report of `reference` against the one engine `engine` under `registry` and the
/// tool's own classes; panics on a run that was stopped.
fn scored(
    scratch: &Scratch,
    reference: &str,
    engine: &str,
    run: &Run,
    registry: &PlantedRegistry,
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

/// What a class answers over two outputs where it is handed its rows directly.
struct Direct {
    /// The fields the class accepts.
    accepted: i64,
    /// The rows on which it accepts a field.
    rows: i64,
    /// The fields whose two values differ.
    differing: i64,
}

/// The answers of `klass` over the rows both outputs write under one key, on each field of
/// `names` whose two values differ. The class is handed the engine's row, the reference's
/// rows of the row's record and, where it is a sibling-allele class, the engine's rows of
/// it, as each whole output states them beside the input's `lines`.
fn consulted_directly(
    klass: &str,
    names: &[&str],
    outputs: (&str, &str),
    lines: &[&str],
    revel_table: Option<&RevelTable>,
) -> Direct {
    let parsed = |path: &str| parse_default(&fs::read_to_string(path).unwrap()).unwrap();
    let (expected, actual) = (parsed(outputs.0), parsed(outputs.1));
    let (of_reference, of_engine) = (
        RecordRows::new(&expected, lines),
        RecordRows::new(&actual, lines),
    );
    let declared = BTreeSet::from([klass.to_string()]);
    let mut direct = Direct {
        accepted: 0,
        rows: 0,
        differing: 0,
    };
    for e in &expected {
        let Some(a) = actual.iter().find(|row| row.key == e.key) else {
            continue;
        };
        let record_rows = of_reference.of(e).unwrap();
        let candidate_rows = match BuiltIn.sibling_allele(klass) {
            true => of_engine.of(a).unwrap(),
            false => None,
        };
        let mut taken = false;
        for name in names {
            let value = |row: &Row| row.fields.get(name).unwrap_or(ABSENT).to_string();
            let (expected_value, actual_value) = (value(e), value(a));
            if expected_value == actual_value {
                continue;
            }
            direct.differing += 1;
            let accepted = BuiltIn
                .accepts(&Consult {
                    klass,
                    e,
                    name,
                    expected: &expected_value,
                    actual: &actual_value,
                    a: Some(a),
                    record_rows: record_rows.as_deref(),
                    candidate_rows: candidate_rows.as_deref(),
                    declared: &declared,
                    assembly: None,
                    source: None,
                    gwas_pairs: None,
                    revel_table,
                    utrannotator_sweep: None,
                    downstream_sequences: None,
                    splice_windows: None,
                    spliceai_table: None,
                    exon_spans: None,
                })
                .unwrap();
            if accepted {
                direct.accepted += 1;
                taken = true;
            }
        }
        direct.rows += i64::from(taken);
    }
    direct
}

/// A sibling-allele class with its plugin, the fields it reads, the first allele's values
/// and the second allele's.
struct SiblingCase {
    klass: &'static str,
    plugin: &'static str,
    names: &'static [&'static str],
    copied: &'static [&'static str],
    own: &'static [&'static str],
}

const SIBLING_ALLELE_CASES: [SiblingCase; 3] = [
    SiblingCase {
        klass: ALPHAMISSENSE_SIBLING_ALLELE_ROWS,
        plugin: "AlphaMissense",
        names: &["am_class", "am_pathogenicity"],
        copied: &["likely_benign", "0.1079"],
        own: &["likely_pathogenic", "0.9"],
    },
    SiblingCase {
        klass: REVEL_SIBLING_ALLELE_ROWS,
        plugin: "REVEL",
        names: &["REVEL"],
        copied: &["0.025"],
        own: &["0.027"],
    },
    SiblingCase {
        klass: SPLICEVAULT_SIBLING_ALLELE_ROWS,
        plugin: "SpliceVault",
        names: &["SpliceVault_site_type", "SpliceVault_site_pos"],
        copied: &["donor", "5"],
        own: &["acceptor", "7"],
    },
];
const SPELT_RECORD: &str = "3_7000200_A/C/T";
const NAMED_RECORD: [&str; 5] = ["3", "7000200", "rs77", "A", "C,T"];

/// The two rows of a made-up record of two alternates, C and T, on one transcript, both
/// giving one amino-acid change: `first` and `second` are the Extra of each allele's row.
fn two_alternate_rows(first: &str, second: &str, uploaded: &str, location: &str) -> Vec<String> {
    [("C", first), ("T", second)]
        .iter()
        .map(|(allele, extra)| {
            format!(
                "{uploaded}\t{location}\t{allele}\tENSG1\tENST1\tTranscript\tmissense_variant\t-\t-\t-\tM/L\t-\t-\t{extra}\n"
            )
        })
        .collect()
}

/// `name=value` of each field, `;` between them.
fn plugin_fields(names: &[&str], values: &[&str]) -> String {
    let pairs: Vec<String> = names
        .iter()
        .zip(values)
        .map(|(name, value)| format!("{name}={value}"))
        .collect();
    pairs.join(";")
}

/// An Extra that opens with the row's impact and writes each field.
fn plugin_extra(names: &[&str], values: &[&str]) -> String {
    format!("IMPACT=MODERATE;{}", plugin_fields(names, values))
}

#[test]
fn sibling_allele_class_is_handed_the_engines_rows_of_the_same_record() {
    for case in &SIBLING_ALLELE_CASES {
        let SiblingCase {
            klass,
            plugin,
            names,
            copied,
            own,
        } = *case;
        let scratch = Scratch::new("sibling-rows");
        let registry = PlantedRegistry::new(&[(plugin, names)]);
        let table = format!("{klass}\t{plugin}\t{}\tmet\t-\n", names.join(","));
        let run_input = scratch.vcf("input.vcf", &[NAMED_RECORD]);
        let held = names.len() as i64;
        let (every, half) = (2 * held, held);
        // The reference writes the first allele's values on both rows, the engine each
        // allele's own: the second allele's values are its sibling's only where the
        // engine's row of the first allele carries them too.
        for (label, uploaded, first, with_input, counted, aside) in [
            (
                "spelt",
                SPELT_RECORD,
                copied,
                false,
                [every, every, half, half, half, half],
                1,
            ),
            (
                "first row differs",
                SPELT_RECORD,
                own,
                false,
                [every, every, 0, every, every, 0],
                0,
            ),
            (
                "named",
                NAMED_RECORD[2],
                copied,
                true,
                [every, every, half, half, half, half],
                1,
            ),
            (
                "named, no input",
                NAMED_RECORD[2],
                copied,
                false,
                [every, every, half, every, every, half],
                0,
            ),
        ] {
            let both = plugin_extra(names, copied);
            let reference = scratch.output(
                "reference.txt",
                &two_alternate_rows(&both, &both, uploaded, "3:7000200"),
            );
            let engine = scratch.output(
                "engine.txt",
                &two_alternate_rows(
                    &plugin_extra(names, first),
                    &plugin_extra(names, own),
                    uploaded,
                    "3:7000200",
                ),
            );
            let run = Run::new().table(&table).plugins(&[plugin]);
            let run = if with_input {
                run.run_input(&run_input)
            } else {
                run
            };
            let report = scored(&scratch, &reference, &engine, &run, &registry);
            for projection in ["feature", "allele"] {
                let item = entry(&report, plugin, "engine", projection);
                assert_eq!(counts(item), counted, "{plugin} {label} {projection}");
                assert_eq!(
                    item["masked"],
                    json!({klass: aside * held}),
                    "{plugin} {label} {projection}"
                );
                assert_eq!(
                    item["residual_rows"],
                    counted[3] + counted[4] - 2 * counted[5],
                    "{plugin} {label} {projection}"
                );
            }
            // Each verdict is the class's own where it is handed the rows of the record as
            // the two whole outputs state them.
            let line = format!("{}\t.\t.\t.", NAMED_RECORD.join("\t"));
            let lines: &[&str] = if with_input { &[&line] } else { &[] };
            let direct = consulted_directly(klass, names, (&reference, &engine), lines, None);
            let item = feature(&report, plugin);
            assert_eq!(
                item["masked"],
                json!({klass: direct.accepted}),
                "{plugin} {label}"
            );
            assert_eq!(direct.rows, aside, "{plugin} {label}");
            assert_eq!(
                item["residual_rows"] == 0,
                direct.accepted == direct.differing,
                "{plugin} {label}"
            );
        }
    }
}

#[test]
fn engines_rows_that_pair_as_equal_text_are_handed_to_the_sibling_allele_classes() {
    let scratch = Scratch::new("sibling-paired");
    let rows: Vec<(&str, &[&str])> = SIBLING_ALLELE_CASES
        .iter()
        .map(|case| (case.plugin, case.names))
        .collect();
    let registry = PlantedRegistry::new(&rows);
    let table: String = SIBLING_ALLELE_CASES
        .iter()
        .map(|case| {
            format!(
                "{}\t{}\t{}\tmet\t-\n",
                case.klass,
                case.plugin,
                case.names.join(",")
            )
        })
        .collect();
    let joined = |values: fn(&SiblingCase) -> &'static [&'static str]| -> String {
        let fields: Vec<String> = SIBLING_ALLELE_CASES
            .iter()
            .map(|case| plugin_fields(case.names, values(case)))
            .collect();
        fields.join(";")
    };
    let (copied, own) = (joined(|case| case.copied), joined(|case| case.own));
    let named = ["3", "7000500", "rs78", "A", "C,T"];
    // Three records of two alternates: one whose first column spells its alleles, one
    // under an identifier, and one on which the engine's first allele's row differs too.
    let records = [
        (SPELT_RECORD, "3:7000200"),
        (named[2], "3:7000500"),
        ("3_7000900_A/C/T", "3:7000900"),
    ];
    let written = |first_of_the_last: &str, second: &str| -> Vec<String> {
        records
            .iter()
            .flat_map(|(uploaded, location)| {
                let first = match *location {
                    "3:7000900" => first_of_the_last,
                    _ => &copied,
                };
                two_alternate_rows(first, second, uploaded, location)
            })
            .collect()
    };
    let reference = scratch.output("reference.txt", &written(&copied, &copied));
    let engine = scratch.output("engine.txt", &written(&own, &own));
    let run_input = scratch.vcf("input.vcf", &[named]);
    let (report, files) = on_any_thread_count(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().table(&table).run_input(&run_input),
        &registry,
        &BuiltIn,
    );
    for case in &SIBLING_ALLELE_CASES {
        let held = case.names.len() as i64;
        let item = feature(&report, case.plugin);
        assert_eq!(
            counts(item),
            [6 * held, 6 * held, 2 * held, 4 * held, 4 * held, 2 * held],
            "{}",
            case.plugin
        );
        assert_eq!(
            item["masked"],
            json!({case.klass: 2 * held}),
            "{}",
            case.plugin
        );
        let left = String::from_utf8(files[&format!("{}.engine.feature.tsv", case.plugin)].clone())
            .unwrap();
        let left: Vec<&str> = left.lines().collect();
        assert_eq!(left.len() as i64, 4 * held, "{}", case.plugin);
        assert_eq!(item["residual_rows"], 4 * held, "{}", case.plugin);
        assert!(
            left.iter()
                .all(|line| line.split('\t').nth(1) == Some("3:7000900")),
            "{}",
            case.plugin
        );
    }
}

const REVEL_NOTICE: &str =
    "no --revel-table: a class of REVEL that reads the table sets nothing aside";

/// The class table of the two classes of REVEL.
fn revel_classes_table() -> String {
    format!(
        "{REVEL_SIBLING_ALLELE_ROWS}\tREVEL\tREVEL\tmet\t-\n{REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS}\tREVEL\tREVEL\tmet\t-\n"
    )
}

/// One BGZF block that holds `data`.
fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut deflater = DeflateEncoder::new(Vec::new(), Compression::new(6));
    deflater.write_all(data).unwrap();
    let body = deflater.finish().unwrap();
    let mut crc = Crc::new();
    crc.update(data);
    let mut block = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00".to_vec();
    block.extend(((body.len() + 25) as u16).to_le_bytes());
    block.extend(body);
    block.extend(crc.sum().to_le_bytes());
    block.extend((data.len() as u32).to_le_bytes());
    block
}

/// A made-up REVEL table on chromosome 3 in the directory `directory` of the scratch,
/// tabix-indexed on its GRCh38 column and each line a BGZF block of its own. A row is a
/// position, an alternate allele, the reference and the alternate amino acid, a score and
/// the transcript ids; the reference allele is A.
fn revel_table_file(scratch: &Scratch, directory: &str, rows: &[[&str; 6]]) -> String {
    let header = bgzf_block(
        b"#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid\n",
    );
    let body: Vec<u8> = rows
        .iter()
        .flat_map(|[position, alternate, aaref, aaalt, score, ids]| {
            bgzf_block(
                format!("3\t999\t{position}\tA\t{alternate}\t{aaref}\t{aaalt}\t{score}\t{ids}\n")
                    .as_bytes(),
            )
        })
        .collect();
    let mut table = header.clone();
    table.extend(&body);
    table.extend(bgzf_block(b""));
    let path = scratch.file(&format!("{directory}/revel.tsv.gz"), &table);
    // One sequence, `3`, with one bin of one chunk that spans every row.
    let mut index = b"TBI\x01".to_vec();
    for number in [1, 0, 1, 3, 3, i32::from(b'#'), 0, 2] {
        index.extend(number.to_le_bytes());
    }
    index.extend(b"3\0");
    index.extend(1i32.to_le_bytes());
    index.extend(0u32.to_le_bytes());
    index.extend(1i32.to_le_bytes());
    index.extend(((header.len() as u64) << 16).to_le_bytes());
    index.extend((((header.len() + body.len()) as u64) << 16).to_le_bytes());
    index.extend(0i32.to_le_bytes());
    scratch.file(&format!("{directory}/revel.tsv.gz.tbi"), gzip(&index));
    path
}

/// A missense row giving M/L, with `score` as its REVEL or, under an empty one, without
/// the field.
fn revel_row(uploaded: &str, location: &str, allele: &str, score: &str) -> String {
    let field = match score {
        "" => String::new(),
        score => format!(";REVEL={score}"),
    };
    format!(
        "{uploaded}\t{location}\t{allele}\tENSG1\tENST1\tTranscript\tmissense_variant\t-\t-\t-\tM/L\t-\t-\tIMPACT=MODERATE{field}\n"
    )
}

fn revel_registry() -> PlantedRegistry {
    PlantedRegistry::new(&[("REVEL", &["REVEL"])])
}

#[test]
fn revel_table_is_handed_to_the_classes_of_revel_which_set_nothing_aside_without_it() {
    let scratch = Scratch::new("revel-table");
    let registry = revel_registry();
    // The table's row at each position lists the transcript; three name another reference
    // amino acid than the rows' and the last the rows' own.
    let table = revel_table_file(
        &scratch,
        "table",
        &[
            ["7000200", "T", "V", "L", "0.300", "ENST1"],
            ["7000500", "T", "V", "L", "0.310", "ENST1"],
            ["7000900", "T", "V", "L", "0.320", "ENST1"],
            ["7001200", "T", "M", "L", "0.330", "ENST1"],
        ],
    );
    let records = [
        ("3_7000200_A/T", "3:7000200", "0.300"),
        ("rs78", "3:7000500", "0.310"),
        ("rs79", "3:7000900", "0.320"),
        ("3_7001200_A/T", "3:7001200", "0.330"),
    ];
    let reference_rows: Vec<String> = records
        .iter()
        .map(|(uploaded, location, score)| revel_row(uploaded, location, "T", score))
        .collect();
    let engine_rows: Vec<String> = records
        .iter()
        .map(|(uploaded, location, _)| revel_row(uploaded, location, "T", ""))
        .collect();
    let reference = scratch.output("reference.txt", &reference_rows);
    let engine = scratch.output("engine.txt", &engine_rows);
    let run_input = scratch.vcf("input.vcf", &[["3", "7000500", "rs78", "A", "T"]]);
    let run = Run::new().table(&revel_classes_table()).plugins(&["REVEL"]);
    let without = scored(
        &scratch,
        &reference,
        &engine,
        &run.clone().run_input(&run_input),
        &registry,
    );
    assert_eq!(without["revel_table"], Value::Null);
    for projection in ["feature", "allele"] {
        let item = entry(&without, "REVEL", "engine", projection);
        assert_eq!(counts(item), [4, 0, 0, 4, 0, 0], "{projection}");
        assert_eq!(item["residual_rows"], 4, "{projection}");
        assert_eq!(
            item["masked"],
            json!({REVEL_SIBLING_ALLELE_ROWS: 0, REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS: 0})
        );
    }
    let report = scored(
        &scratch,
        &reference,
        &engine,
        &run.clone().run_input(&run_input).revel(&table, false),
        &registry,
    );
    assert_eq!(
        report["revel_table"],
        json!({"path": table, "sha256": sha256(&fs::read(&table).unwrap()), "no_match": false})
    );
    for projection in ["feature", "allele"] {
        let item = entry(&report, "REVEL", "engine", projection);
        assert_eq!(counts(item), [4, 0, 0, 1, 0, 0], "{projection}");
        assert_eq!(item["residual_rows"], 1, "{projection}");
        assert_eq!(
            item["masked"],
            json!({REVEL_SIBLING_ALLELE_ROWS: 0, REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS: 3})
        );
    }
    let read =
        RevelTable::new(Path::new(&table), false).unwrap_or_else(|refusal| panic!("{}", refusal.0));
    let direct = consulted_directly(
        REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS,
        &["REVEL"],
        (&reference, &engine),
        &["3\t7000500\trs78\tA\tT\t.\t.\t."],
        Some(&read),
    );
    assert_eq!(
        (direct.accepted, direct.differing - direct.accepted),
        (3, 1)
    );
    // A table whose rows list another transcript explains the values only for a run that
    // set the plugin's no_match.
    let unlisted_rows: Vec<[&str; 6]> = records[..3]
        .iter()
        .map(|(_, location, score)| {
            let position = location.split_once(':').unwrap().1;
            [position, "T", "V", "L", *score, "ENST9"]
        })
        .collect();
    let unlisted = revel_table_file(&scratch, "unlisted", &unlisted_rows);
    let matched = scored(
        &scratch,
        &reference,
        &engine,
        &run.clone().revel(&unlisted, false),
        &registry,
    );
    assert_eq!(
        feature(&matched, "REVEL")["masked"][REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS],
        0
    );
    let unmatched = scored(
        &scratch,
        &reference,
        &engine,
        &run.clone().revel(&unlisted, true),
        &registry,
    );
    assert_eq!(unmatched["revel_table"]["no_match"], true);
    assert_eq!(
        feature(&unmatched, "REVEL")["masked"][REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS],
        3
    );
}

#[test]
fn score_of_another_alternates_row_is_set_aside_under_the_table_on_any_thread_count() {
    let scratch = Scratch::new("revel-alternates");
    // The table's row of C names the T row's alternate amino acid beside another
    // reference amino acid and comes first, so the reference writes its score on the T
    // row; at the last position the row of C names the rows' own reference amino acid.
    let table = revel_table_file(
        &scratch,
        "table",
        &[
            ["7000200", "C", "V", "L", "0.016", "ENST1"],
            ["7000200", "T", "I", "L", "0.029", "ENST1"],
            ["7000500", "C", "V", "L", "0.017", "ENST1"],
            ["7000500", "T", "I", "L", "0.028", "ENST1"],
            ["7000900", "C", "M", "L", "0.018", "ENST1"],
        ],
    );
    let named = ["3", "7000500", "rs78", "A", "C,T"];
    let records = [
        (SPELT_RECORD, "3:7000200", "0.016"),
        (named[2], "3:7000500", "0.017"),
        ("3_7000900_A/C/T", "3:7000900", "0.018"),
    ];
    let rows = |with_scores: bool| -> Vec<String> {
        records
            .iter()
            .flat_map(|(uploaded, location, score)| {
                ["C", "T"].map(|allele| {
                    let score = if with_scores && allele == "T" {
                        score
                    } else {
                        ""
                    };
                    revel_row(uploaded, location, allele, score)
                })
            })
            .collect()
    };
    let reference = scratch.output("reference.txt", &rows(true));
    let engine = scratch.output("engine.txt", &rows(false));
    let run_input = scratch.vcf("input.vcf", &[named]);
    let (report, files) = on_any_thread_count(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new()
            .table(&revel_classes_table())
            .plugins(&["REVEL"])
            .run_input(&run_input)
            .revel(&table, false),
        &revel_registry(),
        &BuiltIn,
    );
    let item = feature(&report, "REVEL");
    assert_eq!(counts(item), [3, 0, 0, 1, 0, 0]);
    assert_eq!(
        item["masked"],
        json!({REVEL_SIBLING_ALLELE_ROWS: 0, REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS: 2})
    );
    assert_eq!(
        files["REVEL.engine.feature.tsv"],
        b"reference\t3:7000900\tT\tENST1\tREVEL\t0.018\n"
    );
}

#[test]
fn summary_says_that_a_class_of_revel_sets_nothing_aside_without_the_table() {
    let scratch = Scratch::new("revel-summary");
    let registry = revel_registry();
    let table = revel_table_file(
        &scratch,
        "table",
        &[["7000200", "T", "V", "L", "0.300", "ENST1"]],
    );
    let reference = scratch.output(
        "reference.txt",
        &[revel_row("3_7000200_A/T", "3:7000200", "T", "0.300")],
    );
    let engine = scratch.output(
        "engine.txt",
        &[revel_row("3_7000200_A/T", "3:7000200", "T", "")],
    );
    let classes = scratch.file("classes.tsv", revel_classes_table());
    let run = |name: &str, more: &[&str]| {
        let engine = format!("engine={engine}");
        let (work, report, residuals) = (
            scratch.path(&format!("work-{name}")),
            scratch.path("report.json"),
            scratch.path(name),
        );
        let mut line = vec![
            "--reference",
            &reference,
            "--engine",
            &engine,
            "--plugins",
            "REVEL",
            "--class-table",
            &classes,
            "--work-dir",
            &work,
            "--sort-memory",
            "16M",
            "--report",
            &report,
            "--residuals",
            &residuals,
        ];
        line.extend(more);
        run_main(&line, &registry, &BuiltIn)
    };
    let other = REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS;

    let (code, lines, error) = run("without", &[]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert_eq!(lines.last().unwrap(), REVEL_NOTICE);
    assert!(lines.contains(&format!(
        "REVEL\tengine\tfeature\tclass {other} set aside 0 (none) [met, -]"
    )));

    let (code, lines, error) = run("with", &["--revel-table", &table]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(!lines.iter().any(|line| line == REVEL_NOTICE));
    assert!(lines.contains(&format!(
        "REVEL\tengine\tfeature\tclass {other} set aside 1 [met, -]"
    )));

    let (code, _, error) = run("no-match", &["--revel-table", &table, "--revel-no-match"]);
    assert_eq!((code, error.as_str()), (0, ""));
    let report: Value =
        serde_json::from_slice(&fs::read(scratch.path("report.json")).unwrap()).unwrap();
    assert_eq!(report["revel_table"]["no_match"], true);

    let nowhere = scratch.path("nowhere.tsv.gz");
    let (code, lines, error) = run("missing", &["--revel-table", &nowhere]);
    assert_eq!((code, lines), (2, vec![]));
    assert_eq!(error, format!("{ERROR}not a file: {nowhere}\n"));

    fs::remove_file(format!("{table}.tbi")).unwrap();
    let (code, lines, error) = run("unindexed", &["--revel-table", &table]);
    assert_eq!((code, lines), (2, vec![]));
    assert!(
        error.starts_with(&format!(
            "{ERROR}--revel-table: cannot read the tabix index of {table}: none of {table}.csi, "
        )),
        "{error}"
    );
}

/// A row of the transcript ENST1 with its consequence terms and its cDNA and CDS positions.
fn coding_row(
    place: (&str, &str),
    terms: &str,
    positions: (&str, &str),
    extra: &str,
    uploaded: &str,
) -> String {
    let (location, allele) = place;
    let (cdna, cds) = positions;
    format!(
        "{uploaded}\t{location}\t{allele}\tENSG1\tENST1\tTranscript\t{terms}\t{cdna}\t{cds}\t-\t-\t-\t-\t{extra}\n"
    )
}

const SPLICE_SPAN_PLUGINS: [(&str, &[&str]); 3] = [
    ("Downstream", &["DownstreamProtein", "ProteinLengthChange"]),
    ("NMD", &["NMD"]),
    ("LoFTEE", &["LoF", "LoF_filter", "LoF_flags", "LoF_info"]),
];

/// A reference and an engine output of seven rows and the input VCF of their run. Three
/// are rows of positions in the transcript's sequence under a longer allele and the
/// identifier `v`, which no line of the input carries: one on which the reference alone
/// holds the donor term and writes NMD and four LOFTEE fields, one on which it holds the
/// acceptor term beside frameshift_variant and the engine alone writes the two Downstream
/// fields, and one on which both sides write a `LoF_info`, with two values. The values of
/// the first two stand as well, under `v`, on a row of two neighbouring positions and on a
/// row whose span crosses an exon's end. Two more rows of two positions inside an exon are
/// read with their record alone: one under the identifier `rs90`, whose record the input
/// holds, and one whose first column spells the alleles.
fn splice_span_outputs(scratch: &Scratch) -> (String, String, String) {
    let in_frame = [
        "splice_donor_variant,protein_altering_variant",
        "protein_altering_variant,splice_region_variant",
    ];
    let shifted = [
        "splice_acceptor_variant,frameshift_variant",
        "frameshift_variant",
    ];
    let taken = (("7:100998-101000", "ACGTAC"), ("998-1000", "898-900"));
    let taken_shifted = (("7:200998-201000", "ACGTA"), ("498-500", "398-400"));
    let two_neighbours = (("7:300999-301000", "ACGTA"), ("999-1000", "899-900"));
    let across = (("7:400998-401000", "ACGTAC"), ("999-?", "899-?"));
    let both_write = (("7:500998-501000", "ACGTAC"), ("298-300", "198-200"));
    let named = (("7:600999-601000", "CGTAC"), ("999-1000", "899-900"));
    let spelt = (("7:700999-701000", "CGTAC"), ("999-1000", "899-900"));
    let spelt_alleles = "7_700998_ACG/ACGTAC";
    let escaping = "NMD=NMD_escaping_variant";
    let reference = scratch.output(
        "reference.txt",
        &[
            coding_row(
                taken.0,
                in_frame[0],
                taken.1,
                &format!(
                    "IMPACT=HIGH;{escaping};LoF=HC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.97"
                ),
                "v",
            ),
            coding_row(
                taken_shifted.0,
                shifted[0],
                taken_shifted.1,
                &format!("IMPACT=HIGH;{escaping};LoF=HC"),
                "v",
            ),
            coding_row(
                two_neighbours.0,
                shifted[0],
                two_neighbours.1,
                &format!("IMPACT=HIGH;{escaping};LoF=HC"),
                "v",
            ),
            coding_row(
                across.0,
                in_frame[0],
                across.1,
                &format!("IMPACT=HIGH;{escaping};LoF=HC"),
                "v",
            ),
            coding_row(
                both_write.0,
                in_frame[0],
                both_write.1,
                "IMPACT=HIGH;LoF=HC;LoF_info=PERCENTILE:0.9",
                "v",
            ),
            coding_row(
                named.0,
                in_frame[0],
                named.1,
                &format!("IMPACT=HIGH;{escaping};LoF=HC"),
                "rs90",
            ),
            coding_row(
                spelt.0,
                in_frame[0],
                spelt.1,
                &format!("IMPACT=HIGH;{escaping};LoF=HC"),
                spelt_alleles,
            ),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            coding_row(taken.0, in_frame[1], taken.1, "IMPACT=MODERATE", "v"),
            coding_row(
                taken_shifted.0,
                shifted[1],
                taken_shifted.1,
                &format!(
                    "IMPACT=HIGH;DownstreamProtein=MKV;ProteinLengthChange=-12;{escaping};LoF=HC"
                ),
                "v",
            ),
            coding_row(
                two_neighbours.0,
                shifted[1],
                two_neighbours.1,
                &format!(
                    "IMPACT=HIGH;DownstreamProtein=MKL;ProteinLengthChange=-3;{escaping};LoF=HC"
                ),
                "v",
            ),
            coding_row(across.0, in_frame[1], across.1, "IMPACT=MODERATE", "v"),
            coding_row(
                both_write.0,
                in_frame[1],
                both_write.1,
                "IMPACT=MODERATE;LoF_info=PERCENTILE:0.5",
                "v",
            ),
            coding_row(named.0, in_frame[1], named.1, "IMPACT=MODERATE", "rs90"),
            coding_row(
                spelt.0,
                in_frame[1],
                spelt.1,
                "IMPACT=MODERATE",
                spelt_alleles,
            ),
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

/// What the class leaves of one plugin's values: the counts, the keys set aside and the
/// residual lines.
type SpliceSpanScore = (&'static str, [i64; 6], i64, Vec<&'static str>);

/// Per plugin, what the class leaves of the outputs of `splice_span_outputs` scored with
/// their input or, where the row under the identifier `rs90` is read without its record
/// and its NMD and LOFTEE values stay residuals, without it.
fn splice_span_scores(with_input: bool) -> Vec<SpliceSpanScore> {
    let across_nmd = "reference\t7:400998-401000\tACGTAC\tENST1\tNMD\tNMD_escaping_variant";
    let loftee = [
        "reference\t7:400998-401000\tACGTAC\tENST1\tLoF\tHC",
        "engine\t7:500998-501000\tACGTAC\tENST1\tLoF_info\tPERCENTILE:0.5",
        "reference\t7:500998-501000\tACGTAC\tENST1\tLoF_info\tPERCENTILE:0.9",
    ];
    let downstream = (
        "Downstream",
        [0, 4, 0, 0, 2, 0],
        2,
        vec![
            "engine\t7:300999-301000\tACGTA\tENST1\tDownstreamProtein\tMKL",
            "engine\t7:300999-301000\tACGTA\tENST1\tProteinLengthChange\t-3",
        ],
    );
    if with_input {
        return vec![
            downstream,
            ("NMD", [6, 2, 2, 3, 2, 2], 3, vec![across_nmd]),
            ("LoFTEE", [11, 3, 2, 4, 3, 2], 7, loftee.to_vec()),
        ];
    }
    let mut left = loftee.to_vec();
    left.push("reference\t7:600999-601000\tCGTAC\tENST1\tLoF\tHC");
    vec![
        downstream,
        (
            "NMD",
            [6, 2, 2, 4, 2, 2],
            2,
            vec![
                across_nmd,
                "reference\t7:600999-601000\tCGTAC\tENST1\tNMD\tNMD_escaping_variant",
            ],
        ),
        ("LoFTEE", [11, 3, 2, 5, 3, 2], 6, left),
    ]
}

#[test]
fn rows_of_a_splice_term_past_the_record_span_are_read_alike_on_any_thread_count() {
    let scratch = Scratch::new("splice-span");
    let klass = PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN;
    let registry = PlantedRegistry::new(&SPLICE_SPAN_PLUGINS);
    let table = format!(
        "{klass}\tDownstream;NMD;LoFTEE\tDownstreamProtein,ProteinLengthChange;NMD;LoF,LoF_info,LoF_filter,LoF_flags\tmet\t-\n"
    );
    let (reference, engine, run_input) = splice_span_outputs(&scratch);
    for with_input in [true, false] {
        let run = Run::new().table(&table);
        let run = if with_input {
            run.run_input(&run_input)
        } else {
            run
        };
        let (report, files) = on_any_thread_count(
            &scratch,
            &reference,
            &[("engine", &engine)],
            &run,
            &registry,
            &BuiltIn,
        );
        for (source, counted, keys, left) in splice_span_scores(with_input) {
            let item = feature(&report, source);
            assert_eq!(counts(item), counted, "{source}, input {with_input}");
            assert_eq!(
                item["masked"],
                json!({klass: keys}),
                "{source}, input {with_input}"
            );
            let lines =
                String::from_utf8(files[&format!("{source}.engine.feature.tsv")].clone()).unwrap();
            assert_eq!(
                lines.lines().collect::<Vec<_>>(),
                left,
                "{source}, input {with_input}"
            );
        }
    }
}

/// The class table's row of the MaxEntScan class, with two further columns.
fn maxentscan_table() -> String {
    format!(
        "{MAXENTSCAN_OVERLAPPING_WINDOWS}\tMaxEntScan\tMaxEntScan_*\tmet\t-\t-\tsets nothing aside without --json-cache\n"
    )
}

const MAXENTSCAN_FIELDS: [&str; 5] = [
    "MaxEntScan_ref",
    "MaxEntScan_alt",
    "MaxEntScan_diff",
    "MaxEntScan_ref_seq",
    "MaxEntScan_alt_seq",
];

fn maxentscan_registry() -> PlantedRegistry {
    PlantedRegistry::new(&[("MaxEntScan", &MAXENTSCAN_FIELDS)])
}

/// A reference, an engine and a candidate cache of chromosome 1 whose ENST1 (forward) and
/// ENST2 (reverse) have an exon of three bases, 201 to 203, between an intron from 101 to
/// 200 and one from 204 to 300, so that 201 to 203 lie in the acceptor window of one
/// intron and the donor window of the next on both, 199 in ENST1's acceptor window alone
/// and 150 in no window. The two sides write the scores of two different windows, and
/// under `verbose` the window texts, on both transcripts at 202 and on ENST1 at 199, and
/// agree at 150.
fn maxentscan_inputs(scratch: &Scratch, verbose: bool) -> (String, String, String) {
    let exons = json!([
        {"start": 1, "end": 100},
        {"start": 201, "end": 203},
        {"start": 301, "end": 400},
    ]);
    let transcripts = json!([
        {"stable_id": "ENST1", "strand": 1, "exons": exons},
        {"stable_id": "ENST2", "strand": -1, "exons": exons},
    ]);
    scratch.file(
        "cache/transcripts/1/1-1000000.json",
        serde_json::to_string(&transcripts).unwrap(),
    );
    // A donor window's nine bases and an acceptor window's twenty-three.
    let (donor_text, acceptor_text) = ("AC".repeat(4) + "G", "AC".repeat(11) + "G");
    let extra = |scores: [&str; 3], window: &str| -> String {
        let mut fields: Vec<String> = MAXENTSCAN_FIELDS[..3]
            .iter()
            .zip(scores)
            .map(|(name, value)| format!("{name}={value}"))
            .collect();
        if verbose {
            let changed = format!("{}T", &window[..window.len() - 1]);
            fields.push(format!("MaxEntScan_ref_seq={window}"));
            fields.push(format!("MaxEntScan_alt_seq={changed}"));
        }
        fields.join(";")
    };
    let donor = extra(["8.100", "2.300", "5.800"], &donor_text);
    let acceptor = extra(["6.400", "1.100", "5.300"], &acceptor_text);
    let places = [("1:202", "ENST1"), ("1:202", "ENST2"), ("1:199", "ENST1")];
    let agreed = row("1:150", "G", "ENST1", &donor);
    let side = |fields: &str| -> Vec<String> {
        places
            .iter()
            .map(|(location, feature)| row(location, "G", feature, fields))
            .chain([agreed.clone()])
            .collect()
    };
    (
        scratch.output("reference.txt", &side(&donor)),
        scratch.output("engine.txt", &side(&acceptor)),
        scratch.path("cache"),
    )
}

#[test]
fn maxentscan_overlapping_windows_is_applied_where_the_cache_places_the_record_in_two_windows() {
    let klass = MAXENTSCAN_OVERLAPPING_WINDOWS;
    for verbose in [false, true] {
        let scratch = Scratch::new("maxentscan-windows");
        let registry = maxentscan_registry();
        let (reference, engine, cache) = maxentscan_inputs(&scratch, verbose);
        let fields: i64 = if verbose { 5 } else { 3 };
        let run = Run::new()
            .table(&maxentscan_table())
            .plugins(&["MaxEntScan"]);
        let residuals = scratch.path("with");
        let report = scored(
            &scratch,
            &reference,
            &engine,
            &run.clone().cache(&cache).residuals(&residuals),
            &registry,
        );
        assert_eq!(report["json_cache"], json!({"path": cache}));
        let (by_feature, by_allele) = (
            feature(&report, "MaxEntScan"),
            allele(&report, "MaxEntScan"),
        );
        assert_eq!(by_feature["masked"], json!({klass: 2 * fields}));
        assert_eq!(by_allele["masked"], json!({klass: fields}));
        assert_eq!(
            counts(by_feature),
            [
                4 * fields,
                4 * fields,
                fields,
                2 * fields,
                2 * fields,
                fields
            ]
        );
        let left = feature_residuals(&residuals, "MaxEntScan");
        assert_eq!(left.len() as i64, 2 * fields);
        let locations: BTreeSet<&str> = left
            .iter()
            .map(|line| line.split('\t').nth(1).unwrap())
            .collect();
        assert_eq!(locations, BTreeSet::from(["1:199"]));
        let without = scored(&scratch, &reference, &engine, &run, &registry);
        assert_eq!(without["json_cache"], Value::Null);
        assert_eq!(feature(&without, "MaxEntScan")["masked"], json!({klass: 0}));
        assert_eq!(allele(&without, "MaxEntScan")["masked"], json!({klass: 0}));
        assert_eq!(
            counts(feature(&without, "MaxEntScan")),
            [
                4 * fields,
                4 * fields,
                fields,
                4 * fields,
                4 * fields,
                fields
            ]
        );
    }
}

#[test]
fn summary_says_that_the_class_of_maxentscan_sets_nothing_aside_without_the_candidate_cache() {
    let scratch = Scratch::new("maxentscan-summary");
    let registry = maxentscan_registry();
    let (reference, engine, cache) = maxentscan_inputs(&scratch, true);
    let table = scratch.file("classes.tsv", maxentscan_table());
    let notice =
        "no --json-cache: a class of MaxEntScan that reads the candidate cache sets nothing aside";
    let run = |name: &str, more: &[&str]| {
        let engine = format!("engine={engine}");
        let (work, report) = (
            scratch.path(&format!("work-{name}")),
            scratch.path("report.json"),
        );
        let mut line = vec![
            "--reference",
            &reference,
            "--engine",
            &engine,
            "--plugins",
            "MaxEntScan",
            "--class-table",
            &table,
            "--work-dir",
            &work,
            "--sort-memory",
            "16M",
            "--report",
            &report,
        ];
        line.extend(more);
        run_main(&line, &registry, &BuiltIn)
    };
    let further = "[met, -, -, sets nothing aside without --json-cache]";
    let klass = MAXENTSCAN_OVERLAPPING_WINDOWS;

    let (code, lines, error) = run("without", &[]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert_eq!(lines.last().unwrap(), notice);
    assert!(lines.contains(&format!(
        "MaxEntScan\tengine\tfeature\tclass {klass} set aside 0 (none) {further}"
    )));

    let (code, lines, error) = run("with", &["--json-cache", &cache]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(!lines.iter().any(|line| line == notice));
    assert!(lines.contains(&format!(
        "MaxEntScan\tengine\tfeature\tclass {klass} set aside 10 {further}"
    )));
    assert!(lines.contains(&format!(
        "MaxEntScan\tengine\tallele\tclass {klass} set aside 5 {further}"
    )));
}

/// The keys under which the class sets a field aside on each run of the two MaxEntScan
/// fixture corpora (`base` the manifest's own run), in the feature projection and in the
/// allele projection: GRCh38 carries one row in two windows and GRCh37 two rows of one
/// record, each with three scores and, under verbose, two window texts.
const MAXENTSCAN_CORPUS_COUNTS: [(&str, &str, i64, i64); 6] = [
    ("GRCh37-plugin-maxentscan-fixture", "base", 6, 3),
    ("GRCh37-plugin-maxentscan-fixture", "swa-ncss", 6, 3),
    (
        "GRCh37-plugin-maxentscan-fixture",
        "swa-ncss-verbose",
        10,
        5,
    ),
    ("GRCh38-plugin-maxentscan-fixture", "base", 3, 3),
    ("GRCh38-plugin-maxentscan-fixture", "swa-ncss", 3, 3),
    ("GRCh38-plugin-maxentscan-fixture", "swa-ncss-verbose", 5, 5),
];

/// The manifest entry of a corpus run: the manifest itself for `base`, else its flag set
/// of that name.
fn corpus_target<'a>(manifest: &'a Value, run: &str) -> &'a Value {
    if run == "base" {
        return manifest;
    }
    manifest["flag_sets"]
        .as_array()
        .unwrap()
        .iter()
        .find(|item| item["name"] == run)
        .unwrap()
}

/// A corpus run's committed default-format output of the reference, unpacked, and the
/// candidate output its manifest entry documents: the same rows with the engine's value in
/// place of the reference's in every field the entry's `field_divergences` name. With
/// `changed_score` the candidate also writes another `MaxEntScan_ref` on the first row the
/// entry documents no field of.
fn documented_candidate(
    corpus: &Path,
    run: &str,
    target: &Value,
    mut changed_score: bool,
) -> (String, String) {
    let expected_dir = match run {
        "base" => "expected",
        _ => target["expected_dir"].as_str().unwrap(),
    };
    let mut reference = String::new();
    MultiGzDecoder::new(fs::File::open(corpus.join(expected_dir).join("default.txt.gz")).unwrap())
        .read_to_string(&mut reference)
        .unwrap();
    let divergences = target["field_divergences"].as_array().unwrap();
    let mut candidate = String::new();
    for line in reference.split_inclusive('\n') {
        let mut columns: Vec<String> = line
            .trim_end_matches('\n')
            .split('\t')
            .map(str::to_string)
            .collect();
        if !line.starts_with('#') {
            let text = |item: &Value, key: &str| item[key].as_str().unwrap().to_string();
            // Per field of the row, the reference's value and the engine's.
            let mut fields: Vec<(String, String, String)> = Vec::new();
            for item in divergences.iter().filter(|item| {
                item["location"] == columns[1]
                    && item["allele"] == columns[2]
                    && item["feature"] == columns[4]
            }) {
                let stated = (text(item, "field"), text(item, "vep"), text(item, "vep_rs"));
                match fields.iter_mut().find(|(field, _, _)| *field == stated.0) {
                    Some(held) => *held = stated,
                    None => fields.push(stated),
                }
            }
            let pairs: Vec<(&str, &str, &str)> = columns[13]
                .split(';')
                .map(|item| match item.split_once('=') {
                    Some((name, value)) => (name, "=", value),
                    None => (item, "", ""),
                })
                .collect();
            if fields.is_empty() && changed_score {
                if let Some((_, _, score)) =
                    pairs.iter().find(|(name, _, _)| *name == "MaxEntScan_ref")
                {
                    fields.push((
                        "MaxEntScan_ref".to_string(),
                        score.to_string(),
                        "99.999".to_string(),
                    ));
                    changed_score = false;
                }
            }
            for (name, reference_value, _) in &fields {
                assert!(
                    pairs.contains(&(name.as_str(), "=", reference_value.as_str())),
                    "{} {name}",
                    columns[1]
                );
            }
            let written: Vec<String> = pairs
                .iter()
                .map(|(name, sign, value)| {
                    match fields.iter().find(|(field, _, _)| field == name) {
                        Some((_, _, candidate_value)) => format!("{name}={candidate_value}"),
                        None => format!("{name}{sign}{value}"),
                    }
                })
                .collect();
            columns[13] = written.join(";");
        }
        candidate.push_str(&columns.join("\t"));
        candidate.push('\n');
    }
    (reference, candidate)
}

#[test]
fn committed_class_table_applies_maxentscan_overlapping_windows_on_the_fixture_corpora_at_their_count(
) {
    let klass = MAXENTSCAN_OVERLAPPING_WINDOWS;
    let tree = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let committed = tree.join("scripts/concordance/plugin_class_table.tsv");
    for (corpus, run, by_feature, by_allele) in MAXENTSCAN_CORPUS_COUNTS {
        let scratch = Scratch::new("maxentscan-corpora");
        let directory = tree.join("tests/golden/116").join(corpus);
        let manifest: Value =
            serde_json::from_slice(&fs::read(directory.join("manifest.json")).unwrap()).unwrap();
        let target = corpus_target(&manifest, run);
        let cache = directory.join("json_cache");
        let entries = |name: &str, with_cache: bool, changed_score: bool| -> (Value, Value) {
            let (reference, candidate) =
                documented_candidate(&directory, run, target, changed_score);
            let reference = scratch.file(&format!("{name}/reference.txt"), reference);
            let candidate = scratch.file(&format!("{name}/candidate.txt"), candidate);
            let options = Run::new()
                .table_file(committed.to_str().unwrap())
                .plugins(&["MaxEntScan"])
                .met_alone();
            let options = if with_cache {
                options.cache(cache.to_str().unwrap())
            } else {
                options
            };
            let report = try_score(
                &scratch,
                &reference,
                &[("vep-rs", &candidate)],
                &options,
                &BuiltInRegistry,
                &BuiltIn,
            )
            .unwrap_or_else(|stop| panic!("{corpus} {run} {name}: {}", message(&stop)));
            (
                entry(&report, "MaxEntScan", "vep-rs", "feature").clone(),
                entry(&report, "MaxEntScan", "vep-rs", "allele").clone(),
            )
        };
        let count = |item: &Value, key: &str| item[key].as_i64().unwrap();

        let (feature, allele) = entries("with", true, false);
        assert_eq!(
            target["field_divergence_summary"],
            json!({klass: by_feature}),
            "{corpus} {run}"
        );
        assert_eq!(
            target["field_divergences"].as_array().unwrap().len() as i64,
            by_feature,
            "{corpus} {run}"
        );
        assert_eq!(
            feature["masked"],
            json!({klass: by_feature}),
            "{corpus} {run}"
        );
        assert_eq!(
            allele["masked"],
            json!({klass: by_allele}),
            "{corpus} {run}"
        );
        for item in [&feature, &allele] {
            assert_eq!(item["adj_perl"], item["adj_rust"], "{corpus} {run}");
            assert_eq!(item["adj_perl"], item["adj_intersection"], "{corpus} {run}");
            assert_eq!(item["residual_rows"], 0, "{corpus} {run}");
            assert_eq!(item["not_applied"], json!({}), "{corpus} {run}");
            assert_eq!(item["perl"], item["rust"], "{corpus} {run}");
            assert!(
                count(item, "rust") > count(item, "intersection"),
                "{corpus} {run}"
            );
        }
        assert_eq!(
            count(&feature, "perl") - count(&feature, "intersection"),
            by_feature,
            "{corpus} {run}"
        );

        let (without, _) = entries("without", false, false);
        assert_eq!(without["masked"], json!({klass: 0}), "{corpus} {run}");
        assert_eq!(
            count(&without, "adj_perl") - count(&without, "adj_intersection"),
            by_feature,
            "{corpus} {run}"
        );

        let (changed, _) = entries("changed", true, true);
        assert_eq!(
            changed["masked"],
            json!({klass: by_feature}),
            "{corpus} {run}"
        );
        assert_eq!(
            count(&changed, "adj_perl") - count(&changed, "adj_intersection"),
            1,
            "{corpus} {run}"
        );
        assert!(count(&changed, "residual_rows") > 0, "{corpus} {run}");
    }
}
