// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field --spliceai-table` on the made-up runs of
//! `tests/fixtures/spliceai/`: a SpliceAI field on a row whose score record holds several
//! entries is an entry of `field_divergences` under `spliceai_multi_entry_cell`, with the
//! entries that state the row's value under the names of the VCF and JSON formats behind
//! it, every other differing field is an `unexplained_residual`, and each manifest is
//! written to the bytes of the fixture's `expected.<n>.json`. The score files are written
//! by the tests, each record line a BGZF block of its own.

mod classify_support;
mod common;
mod recorded_runs;
mod score_files;

use common::refused;
use recorded_runs::{Call, Entry, Recorded, UNEXPLAINED};
use score_files::{
    bgzf_block, record_line, score_file, write_score_file, DELETION, GENEA, GENEB, HEADER,
    SUBSTITUTION,
};
use std::collections::BTreeMap;

const CLASS: &str = "spliceai_multi_entry_cell";

/// The nine fields the VCF format writes a SpliceAI value as.
const MEMBER_FIELDS: [&str; 9] = [
    "SpliceAI_pred_SYMBOL",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL",
    "SpliceAI_pred_DP_DG",
    "SpliceAI_pred_DP_DL",
];

/// The class of each entry of a compared field by its Feature and field.
fn compared<'a>(entries: &'a [Entry], fields: &[&str]) -> BTreeMap<(&'a str, &'a str), &'a str> {
    entries
        .iter()
        .filter(|entry| fields.contains(&entry.field.as_str()))
        .map(|entry| {
            (
                (entry.feature.as_str(), entry.field.as_str()),
                entry.class.as_str(),
            )
        })
        .collect()
}

/// The two values of each entry of a field the run does not compare, by its Feature and
/// field; every entry of the run is of the class.
fn derived<'a>(
    entries: &'a [Entry],
    fields: &[&str],
) -> BTreeMap<(&'a str, &'a str), (&'a str, &'a str)> {
    assert!(entries.iter().all(|entry| entry.class == CLASS));
    entries
        .iter()
        .filter(|entry| !fields.contains(&entry.field.as_str()))
        .map(|entry| {
            (
                (entry.feature.as_str(), entry.field.as_str()),
                (entry.vep.as_str(), entry.vep_rs.as_str()),
            )
        })
        .collect()
}

#[test]
fn classify_names_the_spliceai_multi_entry_cell_class_under_the_score_files() {
    // A made-up score record of two entries, GENEA's then GENEB's: Ensembl VEP writes the
    // record's whole value on the GENEA transcript's row and nothing on the GENEB
    // transcript's, vep-rs each gene's own entry and a flag on both. Under
    // `--spliceai-table` the three differing fields are `spliceai_multi_entry_cell`, on rows
    // that state their gene as `SYMBOL` and on rows whose gene the corpus's cache names;
    // without the files, under a file that scores each gene in a record of its own, and on
    // rows whose gene nothing names, each is an `unexplained_residual`. So is a vep-rs value
    // that is the other gene's entry, with the flag beside it, and a first-entry row on
    // which Ensembl VEP writes no value.
    let run = Recorded::new("spliceai/record_of_two_entries");
    score_file(&run.root().join("joined"), DELETION, true);
    score_file(&run.root().join("single"), DELETION, false);
    let fields = ["SpliceAI_pred", "SpliceAI_cutoff"];
    let joined = ["joined/indel.synthetic.vcf.gz"];
    let single = ["single/indel.synthetic.vcf.gz"];
    let keys = [
        ("ENST1", "SpliceAI_pred"),
        ("ENST2", "SpliceAI_pred"),
        ("ENST2", "SpliceAI_cutoff"),
    ];
    let all = |class: &'static str| -> BTreeMap<(&str, &str), &str> {
        keys.iter().map(|key| (*key, class)).collect()
    };
    let classes =
        |outputs: (&str, &str), tables: &[&str], summary: (usize, &str), expected: &str| {
            let entries = run.holds(&Call {
                directory: "",
                outputs,
                fields: &fields,
                tables,
                summary,
                expected,
            });
            let found: Vec<((String, String), String)> = compared(&entries, &fields)
                .into_iter()
                .map(|((feature, field), class)| {
                    ((feature.to_string(), field.to_string()), class.to_string())
                })
                .collect();
            found
        };
    let owned = |expected: BTreeMap<(&str, &str), &str>| -> Vec<((String, String), String)> {
        expected
            .into_iter()
            .map(|((feature, field), class)| {
                ((feature.to_string(), field.to_string()), class.to_string())
            })
            .collect()
    };
    let rows = ("vep.txt", "rs.txt");
    let unexplained = (3, "{'unexplained_residual': 3}");
    let taken = (15, "{'spliceai_multi_entry_cell': 15}");
    assert_eq!(
        classes(rows, &[], unexplained, "expected.1.json"),
        owned(all(UNEXPLAINED))
    );
    assert_eq!(
        classes(rows, &joined, taken, "expected.2.json"),
        owned(all(CLASS))
    );
    assert_eq!(
        classes(rows, &single, unexplained, "expected.3.json"),
        owned(all(UNEXPLAINED))
    );
    let other = (
        5,
        "{'spliceai_multi_entry_cell': 3, 'unexplained_residual': 2}",
    );
    let mut expected = all(UNEXPLAINED);
    expected.insert(keys[0], CLASS);
    assert_eq!(
        classes(("vep.txt", "other.txt"), &joined, other, "expected.4.json"),
        owned(expected)
    );
    let silent = (
        14,
        "{'unexplained_residual': 2, 'spliceai_multi_entry_cell': 12}",
    );
    let mut expected = all(CLASS);
    expected.insert(("ENST1", "SpliceAI_pred"), UNEXPLAINED);
    expected.insert(("ENST1", "SpliceAI_cutoff"), UNEXPLAINED);
    assert_eq!(
        classes(("silent.txt", "rs.txt"), &joined, silent, "expected.5.json"),
        owned(expected)
    );
    // Rows that state no SYMBOL: nothing names their gene until the corpus holds a cache
    // that does.
    let bare = ("bare-vep.txt", "bare-rs.txt");
    assert_eq!(
        classes(bare, &joined, unexplained, "expected.6.json"),
        owned(all(UNEXPLAINED))
    );
    run.write(
        "c/json_cache/transcripts/21/1-1000000.json",
        br#"[{"stable_id": "ENST1", "gene_symbol": "GENEA"}, {"stable_id": "ENST2", "gene_symbol": "GENEB"}]"#,
    );
    assert_eq!(
        classes(bare, &joined, taken, "expected.7.json"),
        owned(all(CLASS))
    );
}

#[test]
fn classify_names_the_spliceai_multi_entry_cell_class_on_the_rows_of_a_single_base_record() {
    // A made-up score record of a substitution, one base for one base, of two entries,
    // GENEA's then GENEB's, whose rows stand at the one position with the alternate base as
    // their allele: Ensembl VEP writes the record's whole value on the GENEA transcript's
    // row and nothing on the GENEB transcript's, vep-rs each gene's own entry. Under
    // `--spliceai-table` both differing fields are `spliceai_multi_entry_cell`; without the
    // files, and under a file that scores each gene in a record of its own, each is an
    // `unexplained_residual`. So is a vep-rs value that is the other gene's entry, and an
    // Ensembl VEP value that is the first entry joined with itself.
    let run = Recorded::new("spliceai/single_base_record");
    score_file(&run.root().join("joined"), SUBSTITUTION, true);
    score_file(&run.root().join("single"), SUBSTITUTION, false);
    let fields = ["SpliceAI_pred"];
    let joined = ["joined/snv.synthetic.vcf.gz"];
    let single = ["single/snv.synthetic.vcf.gz"];
    let held = |outputs: (&str, &str),
                tables: &[&str],
                summary: (usize, &str),
                expected: &str,
                classes: [&str; 2]| {
        let entries = run.holds(&Call {
            directory: "",
            outputs,
            fields: &fields,
            tables,
            summary,
            expected,
        });
        let found = compared(&entries, &fields);
        let features = [("ENST1", "SpliceAI_pred"), ("ENST2", "SpliceAI_pred")];
        assert_eq!(found, features.into_iter().zip(classes).collect());
    };
    let rows = ("vep.txt", "rs.txt");
    let unexplained = (2, "{'unexplained_residual': 2}");
    held(rows, &[], unexplained, "expected.1.json", [UNEXPLAINED; 2]);
    let taken = (14, "{'spliceai_multi_entry_cell': 14}");
    held(rows, &joined, taken, "expected.2.json", [CLASS; 2]);
    held(
        rows,
        &single,
        unexplained,
        "expected.3.json",
        [UNEXPLAINED; 2],
    );
    let other = (
        4,
        "{'spliceai_multi_entry_cell': 3, 'unexplained_residual': 1}",
    );
    held(
        ("vep.txt", "other.txt"),
        &joined,
        other,
        "expected.4.json",
        [CLASS, UNEXPLAINED],
    );
    let twice = (
        12,
        "{'unexplained_residual': 1, 'spliceai_multi_entry_cell': 11}",
    );
    held(
        ("twice.txt", "rs.txt"),
        &joined,
        twice,
        "expected.5.json",
        [UNEXPLAINED, CLASS],
    );
}

#[test]
fn classify_names_the_spliceai_multi_entry_cell_class_on_a_row_whose_gene_two_entries_name() {
    // A made-up score record of a substitution whose two entries both name GENEA, with
    // other scores: both engines write the record's value on the GENEA transcript's row in
    // the default format, so no compared field differs, while the nine members differ,
    // Ensembl VEP's cut from the value and vep-rs's each its two values joined by a comma.
    // Under `--spliceai-table` the row holds the nine member fields and the JSON object
    // under `spliceai_multi_entry_cell`; without the files it holds nothing; and a vep-rs
    // value that is the later entry alone is an `unexplained_residual`.
    let run = Recorded::new("spliceai/entries_of_one_gene");
    let other = "GENEA|0.31|0.00|0.44|0.09|3|-2|11|-30";
    write_score_file(
        &run.root().join("table"),
        "snv.synthetic.vcf.gz",
        &[record_line(SUBSTITUTION, &[GENEA, other])],
    );
    let fields = ["SpliceAI_pred"];
    let table = ["table/snv.synthetic.vcf.gz"];
    let held = |outputs: (&str, &str), tables: &[&str], summary: (usize, &str), expected: &str| {
        run.holds(&Call {
            directory: "",
            outputs,
            fields: &fields,
            tables,
            summary,
            expected,
        })
    };
    let alike = ("vep.txt", "rs.txt");
    let taken = (10, "{'spliceai_multi_entry_cell': 10}");
    let entries = held(alike, &table, taken, "expected.1.json");
    let values = derived(&entries, &fields);
    assert_eq!(values.len(), 10);
    assert_eq!(
        values[&("ENST1", "SpliceAI_pred_SYMBOL")],
        ("GENEA", "GENEA,GENEA")
    );
    assert_eq!(
        values[&("ENST1", "SpliceAI_pred_DS_AG")],
        ("0.15", "0.15,0.31")
    );
    assert_eq!(
        values[&("ENST1", "SpliceAI_pred_DP_DL")],
        ("-9,C", "-9,-30")
    );
    assert_eq!(
        values[&("ENST1", "spliceai")],
        (
            r#"{"DP_AG":-6,"DP_AL":-41,"DP_DG":-17,"DP_DL":"-9,C","DS_AG":0.15,"DS_AL":0.08,"DS_DG":0.06,"DS_DL":0.17,"SYMBOL":"GENEA"}"#,
            r#"{"DP_AG":"-6,3","DP_AL":"-41,-2","DP_DG":"-17,11","DP_DL":"-9,-30","DS_AG":"0.15,0.31","DS_AL":"0.08,0.00","DS_DG":"0.06,0.44","DS_DL":"0.17,0.09","SYMBOL":"GENEA,GENEA"}"#
        )
    );
    assert!(held(alike, &[], (0, "{}"), "expected.2.json").is_empty());
    let later = held(
        ("vep.txt", "later.txt"),
        &table,
        (1, "{'unexplained_residual': 1}"),
        "expected.3.json",
    );
    assert_eq!(
        compared(&later, &fields),
        [(("ENST1", "SpliceAI_pred"), UNEXPLAINED)]
            .into_iter()
            .collect()
    );
}

#[test]
fn classify_names_the_spliceai_multi_entry_cell_class_under_every_record_line_at_the_allele() {
    // A made-up score file whose record line of two entries, GENEA's then GENEB's, stands
    // more than once at the substitution's allele, or beside a line of the same entries in
    // the other order or of GENEA's entry alone: Ensembl VEP writes a row the whole value of
    // the last line its gene leads, or nothing where it leads none, vep-rs the gene's entry.
    // Under `--spliceai-table` every `SpliceAI_pred` the two fill differently is
    // `spliceai_multi_entry_cell`; where GENEA's line of one entry stands last both write
    // that entry on its row, which differs in nothing. Where the lines give a gene entries
    // that differ in text, vep-rs writes the entry of the last line that names the gene, and
    // its rows are `spliceai_multi_entry_cell` under two and three lines of the same genes
    // with other scores and under two lines of other scores that each gene leads one of.
    // The row of a gene no entry names, and a row on which vep-rs's value is an earlier
    // line's entry, are each an `unexplained_residual`, and so is every differing row
    // without the files.
    let changed = "GENEA|0.31|0.00|0.44|0.09|3|-2|11|-30";
    let other = "GENEB|0.15|0.08|0.06|0.17|-6|-41|-17|-9";
    let third = "GENEA|0.00|1.00|0.10|0.05|-25|-23|37|17";
    let last = "GENEB|0.31|0.00|0.44|0.09|3|-2|11|-30";
    let both: &[&str] = &[GENEA, GENEB];
    let all = (14, "{'spliceai_multi_entry_cell': 14}");
    let none = (2, "{'unexplained_residual': 2}");
    let two = [("ENST1", CLASS), ("ENST2", CLASS)];
    // The directory of a case, its record lines, the classes of its rows under the score
    // file, and the field divergences counted under the file and without it.
    #[allow(clippy::type_complexity)]
    let cases: [(
        &str,
        &[&[&str]],
        &[(&str, &str)],
        (usize, &str),
        (usize, &str),
    ); 11] = [
        ("one-line", &[both], &two, all, none),
        ("twice", &[both, both], &two, all, none),
        ("three-times", &[both, both, both], &two, all, none),
        (
            "other-order",
            &[both, &[GENEB, GENEA]],
            &two,
            (6, "{'spliceai_multi_entry_cell': 6}"),
            none,
        ),
        (
            "one-entry-last",
            &[both, &[GENEA]],
            &[("ENST2", CLASS)],
            (11, "{'spliceai_multi_entry_cell': 11}"),
            (1, "{'unexplained_residual': 1}"),
        ),
        ("one-entry-first", &[&[GENEA], both], &two, all, none),
        (
            "gene-unnamed",
            &[both, both],
            &[("ENST1", CLASS), ("ENST2", CLASS), ("ENST3", UNEXPLAINED)],
            (
                15,
                "{'spliceai_multi_entry_cell': 14, 'unexplained_residual': 1}",
            ),
            (3, "{'unexplained_residual': 3}"),
        ),
        ("other-scores", &[both, &[changed, other]], &two, all, none),
        (
            "other-scores-three-times",
            &[both, &[changed, other], &[third, last]],
            &two,
            all,
            none,
        ),
        (
            "other-scores-other-order",
            &[both, &[other, changed]],
            &two,
            (13, "{'spliceai_multi_entry_cell': 13}"),
            none,
        ),
        (
            "earlier-line",
            &[both, &[changed, other]],
            &[("ENST1", UNEXPLAINED), ("ENST2", CLASS)],
            (
                12,
                "{'unexplained_residual': 1, 'spliceai_multi_entry_cell': 11}",
            ),
            none,
        ),
    ];
    let run = Recorded::new("spliceai/record_lines");
    let fields = ["SpliceAI_pred"];
    for (place, (directory, lines, classes, with_file, without)) in cases.into_iter().enumerate() {
        let lines: Vec<String> = lines
            .iter()
            .map(|entries| record_line(SUBSTITUTION, entries))
            .collect();
        let scores = run.root().join(directory).join("scores");
        write_score_file(&scores, "snv.synthetic.vcf.gz", &lines);
        let table = format!("{directory}/scores/snv.synthetic.vcf.gz");
        let held = |tables: &[&str], summary: (usize, &str), call: usize| {
            let entries = run.holds(&Call {
                directory,
                outputs: ("vep.txt", "rs.txt"),
                fields: &fields,
                tables,
                summary,
                expected: &format!("expected.{call}.json"),
            });
            let found: Vec<(String, String)> = compared(&entries, &fields)
                .into_iter()
                .map(|((feature, _), class)| (feature.to_string(), class.to_string()))
                .collect();
            found
        };
        let owned = |classes: &[(&str, &str)]| -> Vec<(String, String)> {
            classes
                .iter()
                .map(|(feature, class)| (feature.to_string(), class.to_string()))
                .collect()
        };
        assert_eq!(
            held(&[&table], with_file, 2 * place + 1),
            owned(classes),
            "{directory}"
        );
        let unexplained: Vec<(&str, &str)> = classes
            .iter()
            .map(|(feature, _)| (*feature, UNEXPLAINED))
            .collect();
        assert_eq!(
            held(&[], without, 2 * place + 2),
            owned(&unexplained),
            "{directory}"
        );
    }
}

#[test]
fn classify_writes_the_member_fields_of_a_row_of_the_spliceai_class_for_the_vcf_format() {
    // The VCF format writes a SpliceAI value as nine fields whatever the run's options, and
    // the golden harness reads a documented field by its format's name. For a row written
    // under `spliceai_multi_entry_cell` in a run whose default format writes the one field,
    // `classify` also writes each member field whose two parts differ: the last position
    // alone on the first entry's gene's row, where Ensembl VEP's ninth part runs on into the
    // second entry's allele, and all nine on the second entry's gene's row. A run that
    // compares the nine fields itself gains none.
    let run = Recorded::new("spliceai/member_fields");
    score_file(&run.root().join("joined"), DELETION, true);
    score_file(&run.root().join("split/joined"), DELETION, true);
    let fields = ["SpliceAI_pred"];
    let entries = run.holds(&Call {
        directory: "",
        outputs: ("vep.txt", "rs.txt"),
        fields: &fields,
        tables: &["joined/indel.synthetic.vcf.gz"],
        summary: (14, "{'spliceai_multi_entry_cell': 14}"),
        expected: "expected.1.json",
    });
    let mut members = derived(&entries, &fields);
    assert!(members.remove(&("ENST1", "spliceai")).is_some());
    assert!(members.remove(&("ENST2", "spliceai")).is_some());
    let mut expected: BTreeMap<(&str, &str), (&str, &str)> = MEMBER_FIELDS
        .into_iter()
        .zip(GENEB.split('|'))
        .map(|(field, part)| (("ENST2", field), ("-", part)))
        .collect();
    expected.insert(("ENST1", "SpliceAI_pred_DP_DL"), ("-9,G", "-9"));
    assert_eq!(members, expected);
    let entries = run.holds(&Call {
        directory: "split",
        outputs: ("vep.txt", "rs.txt"),
        fields: &MEMBER_FIELDS,
        tables: &["split/joined/indel.synthetic.vcf.gz"],
        summary: (12, "{'spliceai_multi_entry_cell': 12}"),
        expected: "expected.2.json",
    });
    let objects: Vec<(&str, &str)> = derived(&entries, &MEMBER_FIELDS).into_keys().collect();
    assert_eq!(objects, [("ENST1", "spliceai"), ("ENST2", "spliceai")]);
}

#[test]
fn classify_writes_the_json_object_of_a_row_of_the_spliceai_class_as_the_harness_prints_one() {
    // The JSON format writes a SpliceAI value as one object, which the golden harness prints
    // with its members in the byte order of their names and no spaces. For a row written
    // under `spliceai_multi_entry_cell`, `classify` also writes that object for both
    // engines: a number for each member whose text is one, written by its value as both
    // engines write it (`0.00` as `0`, `1.00` as `1`, `0.10` as `0.1`), Ensembl VEP's last
    // position a text where it runs on into the second entry's allele, the cutoff flag
    // among the members of a row that carries one, and no object for a row without the
    // value.
    let run = Recorded::new("spliceai/json_object");
    score_file(&run.root().join("joined"), DELETION, true);
    let fields = ["SpliceAI_pred", "SpliceAI_cutoff"];
    let entries = run.holds(&Call {
        directory: "",
        outputs: ("vep.txt", "rs.txt"),
        fields: &fields,
        tables: &["joined/indel.synthetic.vcf.gz"],
        summary: (15, "{'spliceai_multi_entry_cell': 15}"),
        expected: "expected.1.json",
    });
    let derived = derived(&entries, &fields);
    let own = r#"{"DP_AG":-6,"DP_AL":-41,"DP_DG":-17,"DP_DL":-9,"DS_AG":0.15,"DS_AL":0.08,"DS_DG":0.06,"DS_DL":0.17,"SYMBOL":"GENEA","SpliceAI_cutoff":"FAIL"}"#;
    let run_on = own.replace(r#""DP_DL":-9"#, r#""DP_DL":"-9,G""#);
    assert_eq!(derived[&("ENST1", "spliceai")], (run_on.as_str(), own));
    assert_eq!(
        derived[&("ENST2", "spliceai")],
        (
            "-",
            r#"{"DP_AG":-25,"DP_AL":-23,"DP_DG":37,"DP_DL":17,"DS_AG":0,"DS_AL":1,"DS_DG":0.1,"DS_DL":0.05,"SYMBOL":"GENEB","SpliceAI_cutoff":"PASS"}"#
        )
    );
}

#[test]
fn classify_reads_the_spliceai_names_of_the_vcf_and_json_formats_from_their_outputs() {
    // A run that names the VCF and JSON outputs of both engines reads the nine fields and
    // the object from them and derives nothing: it writes the entries the derivation writes
    // from the default format alone (the fixture's expected manifest), each under the class
    // by the format comparator's answer, the member fields of a row in the byte order of
    // their names. A JSON text that holds a tab cannot stand in a line of the comparator's
    // questions: its field is asked of no class and stays an `unexplained_residual`.
    let run = Recorded::new("spliceai/json_object");
    score_file(&run.root().join("joined"), DELETION, true);
    let fields = ["SpliceAI_pred", "SpliceAI_cutoff"];
    let call = Call {
        directory: "",
        outputs: ("vep.txt", "rs.txt"),
        fields: &fields,
        tables: &["joined/indel.synthetic.vcf.gz"],
        summary: (15, "{'spliceai_multi_entry_cell': 15}"),
        expected: "expected.1.json",
    };
    let derived = run.holds(&call);
    let record = "21\t5000\tdel_two_genes\tGCA\tG\t.\t.";
    let vcf = |first: &str, second: &str| {
        let members = MEMBER_FIELDS.map(|field| field.strip_prefix("SpliceAI_pred_").unwrap());
        let mut names: Vec<String> = members
            .iter()
            .map(|member| format!("SpliceAI_pred_{member}"))
            .collect();
        names.sort();
        format!(
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SYMBOL|Feature_type|Feature|SpliceAI_cutoff|{}\">\n\
             {record}\tCSQ=-|intron_variant|GENEA|Transcript|ENST1|{first},-|intron_variant|GENEB|Transcript|ENST2|{second}\n",
            names.join("|")
        )
    };
    let own = "FAIL|-6|-41|-17|-9|0.15|0.08|0.06|0.17|GENEA";
    run.write(
        "vep.vcf",
        vcf(&own.replace("|-9|", "|-9&G|"), "|||||||||").as_bytes(),
    );
    run.write(
        "rs.vcf",
        vcf(own, "PASS|-25|-23|37|17|0.00|1.00|0.10|0.05|GENEB").as_bytes(),
    );
    let json = |first: &str, second: &str| {
        let object = |feature: &str, symbol: &str, more: &str| {
            format!(
                r#"{{"variant_allele": "-", "transcript_id": "{feature}", "consequence_terms": ["intron_variant"], "gene_symbol": "{symbol}"{more}}}"#
            )
        };
        format!(
            r#"{{"input": "21\t5000\tdel_two_genes\tGCA\tG\t.\t.\t.", "transcript_consequences": [{}, {}]}}"#,
            object("ENST1", "GENEA", first),
            object("ENST2", "GENEB", second)
        ) + "\n"
    };
    let own = r#", "spliceai": {"SYMBOL": "GENEA", "DS_AG": 0.15, "DS_AL": 0.08, "DS_DG": 0.06, "DS_DL": 0.17, "DP_AG": -6, "DP_AL": -41, "DP_DG": -17, "DP_DL": -9, "SpliceAI_cutoff": "FAIL"}"#;
    let run_on = own.replace(r#""DP_DL": -9"#, r#""DP_DL": "-9,G""#);
    let run_on = format!(r#", "note": "a\tb"{run_on}"#);
    let own = format!(r#", "note": "a\tc"{own}"#);
    run.write("vep.json", json(&run_on, "").as_bytes());
    let other = r#", "spliceai": {"SYMBOL": "GENEB", "DS_AG": 0, "DS_AL": 1, "DS_DG": 0.10, "DS_DL": 0.05, "DP_AG": -25, "DP_AL": -23, "DP_DG": 37, "DP_DL": 17, "SpliceAI_cutoff": "PASS"}"#;
    run.write("rs.json", json(&own, other).as_bytes());
    let mut arguments = run.arguments(&call);
    for (option, output) in [
        ("--vep-vcf", "vep.vcf"),
        ("--vep-rs-vcf", "rs.vcf"),
        ("--vep-json", "vep.json"),
        ("--vep-rs-json", "rs.json"),
    ] {
        arguments.extend([option.to_string(), run.path(output)]);
    }
    let ran = common::run(&arguments);
    assert_eq!((ran.code, ran.stderr.as_str()), (0, ""));
    let counted = "{'spliceai_multi_entry_cell': 15, 'unexplained_residual': 1}";
    assert!(ran
        .stdout
        .ends_with(&format!("; field divergences 16 ({counted})\n")));
    // The derivation's entries of a row, its compared fields first, then the member fields
    // by name, then the fields of the JSON format by name, the note among them.
    let place = |entry: &Entry| {
        let group = match entry.field.as_str() {
            field if fields.contains(&field) => 0,
            "spliceai" | "note" => 2,
            _ => 1,
        };
        let name = if group == 0 {
            String::new()
        } else {
            entry.field.clone()
        };
        (entry.feature.clone(), group, name)
    };
    let mut expected = derived;
    let first = expected[0].clone();
    assert_eq!(first.feature, "ENST1");
    expected.push(Entry {
        field: "note".to_string(),
        vep: "a\tb".to_string(),
        vep_rs: "a\tc".to_string(),
        class: UNEXPLAINED.to_string(),
        ..first
    });
    expected.sort_by_key(place);
    assert_eq!(recorded_runs::entries(&run.manifest("")), expected);
}

#[test]
fn classify_ends_on_a_score_file_it_cannot_read_with_nothing_written() {
    let run = Recorded::new("spliceai/record_of_two_entries");
    let before = run.manifest("");
    let fields = ["SpliceAI_pred", "SpliceAI_cutoff"];
    let asked = Call {
        directory: "",
        outputs: ("vep.txt", "rs.txt"),
        fields: &fields,
        tables: &["joined/indel.synthetic.vcf.gz"],
        summary: (0, "{}"),
        expected: "expected.2.json",
    };
    let own = |detail: String| {
        refused(
            1,
            &format!("ERROR: [vep-golden-corpus] --spliceai-table: {detail}"),
        )
    };
    // A file without its index ends the run before a row is read.
    let file = score_file(&run.root().join("joined"), DELETION, true);
    let index = format!("{file}.tbi");
    let kept = std::fs::read(&index).unwrap();
    std::fs::remove_file(&index).unwrap();
    let missing = format!(
        "{index} is not a readable tabix index ([Errno 2] No such file or directory: '{index}')"
    );
    assert_eq!(run.classify(&asked), own(missing));
    assert_eq!(run.manifest(""), before);
    // A block that does not inflate, behind the record the rows stand at, ends the run at
    // the first row that asks for it.
    std::fs::write(&index, kept).unwrap();
    let lines = [
        record_line(DELETION, &[GENEA, GENEB]),
        record_line(("9000", "GCA", "G"), &[GENEA, GENEB]),
    ];
    let file = write_score_file(&run.root().join("joined"), "indel.synthetic.vcf.gz", &lines);
    let header = bgzf_block((HEADER.join("\n") + "\n").as_bytes());
    let first = bgzf_block(format!("{}\n", lines[0]).as_bytes());
    let mut second = bgzf_block(format!("{}\n", lines[1]).as_bytes());
    let size = second.len();
    second.truncate(18);
    second.resize(size, 0xff);
    std::fs::write(&file, [header, first, second, bgzf_block(b"")].concat()).unwrap();
    assert_eq!(
        run.classify(&asked),
        own(format!("{file} holds a block that does not inflate"))
    );
    assert_eq!(run.manifest(""), before);
}
