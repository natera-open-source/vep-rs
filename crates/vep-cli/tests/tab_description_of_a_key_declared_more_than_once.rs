// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tab format's description line of a key that more than one plugin or
//! `--custom` source declares, in the headers the binary writes. Ensembl VEP
//! 116.3 lists such a key once for each declaration
//! (`OutputFactory/Tab.pm:195`) and reads the text of each of its description
//! lines from one hash of the declared pairs, the plugins' and then the
//! sources' (`:149`, `:153`), so every line carries the text of the
//! declaration listed last. The cell under the key holds the value of the
//! plugin loaded last among those that wrote the row, and a source's record
//! names where no plugin wrote it (`OutputFactory.pm:2239-2241`). The key is
//! one tab column whose one line carries the text of every declaration, in the
//! order their values stand on a row: the plugins from the last loaded to the
//! first, then the source. Under `--most_severe` and `--summary` no plugin runs
//! on a row (`OutputFactory.pm:387`), and the line is the one text Ensembl VEP
//! prints. The rows held here are the ones Ensembl VEP writes, its second
//! listing of the key left out, and the default and VCF headers are Ensembl
//! VEP's own. The one listing and its text are the documented difference
//! `duplicate_output_column_collapsed`.

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The text the REVEL plugin declares for its key.
const REVEL_TEXT: &str = "Rare Exome Variant Ensemble Learner";

/// The text a source declares for its own key: its file as a header prints it.
const SOURCE_TEXT: &str = "[PATH]/site.vcf.gz";

/// The text the CADD plugin declares for `CADD_PHRED`.
const CADD_TEXT: &str = "PHRED-like scaled CADD score. CADD is only available here for non-commercial use. See CADD website for more information.";

/// The tab format's column line up to the columns a plugin or a source adds.
const TAB_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS";

/// The six rows Ensembl VEP 116.3 writes for line 4 of the REVEL corpus's
/// input under the REVEL plugin and a source named `REVEL`, its second `REVEL`
/// column left out: the plugin's score on the one row its table names and the
/// source's record name on the other five.
const SCORED_RECORD_ROWS: [&str; 6] = [
    "4637665\t21:26965413\tC\tENSG00000154736\tENST00000284987\tTranscript\tmissense_variant\t1676\t979\t327\tS/G\tAgc/Ggc\t-\tMODERATE\t-\t-1\t-\t-\t0.176",
    "4637665\t21:26965413\tC\tENSG00000154736\tENST00000652031\tTranscript\tmissense_variant,NMD_transcript_variant\t161\t163\t55\tS/G\tAgc/Ggc\t-\tMODERATE\t-\t-1\tcds_start_NF\t-\tsite_1",
    "4637665\t21:26965413\tC\tENSG00000154736\tENST00000970346\tTranscript\tmissense_variant\t2229\t979\t327\tS/G\tAgc/Ggc\t-\tMODERATE\t-\t-1\t-\t-\tsite_1",
    "4637665\t21:26965413\tC\tENSG00000154736\tENST00001108182\tTranscript\tmissense_variant\t2229\t979\t327\tS/G\tAgc/Ggc\t-\tMODERATE\t-\t-1\t-\t-\tsite_1",
    "4637665\t21:26965413\tC\tENSG00000154736\tENST00001108183\tTranscript\tmissense_variant\t2229\t979\t327\tS/G\tAgc/Ggc\t-\tMODERATE\t-\t-1\t-\t-\tsite_1",
    "4637665\t21:26965413\tC\tENSG00000154736\tENST00001131497\tTranscript\tmissense_variant\t2229\t979\t327\tS/G\tAgc/Ggc\t-\tMODERATE\t-\t-1\t-\t-\tsite_1",
];

/// The two rows Ensembl VEP 116.3 writes for line 9 of the custom corpus's
/// input under the same pair, its second `REVEL` column left out: the record
/// overlaps no transcript, the plugin writes neither row and both hold the
/// source's record name.
const UNSCORED_RECORD_ROWS: [&str; 2] = [
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2",
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2",
];

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The fixture of two plugins that write one key.
fn shared_key_fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/plugin-shared-key")
}

/// Line `line` of the input of the corpus `name`, alone under its header
/// lines, written under `dir`; `opening` is the record's chromosome and
/// position, which the line must open with.
fn record_file(dir: &Path, name: &str, line: usize, opening: &str) -> PathBuf {
    let input = fs::read_to_string(corpus(name).join("variants.vcf")).unwrap();
    let record = input.lines().nth(line - 1).expect("the record's line");
    assert!(
        record.starts_with(opening),
        "line {line} of the input of {name} is the record"
    );
    let mut text: String = input
        .lines()
        .filter(|header| header.starts_with('#'))
        .map(|header| format!("{header}\n"))
        .collect();
    text.push_str(record);
    text.push('\n');
    let path = dir.join("record.vcf");
    fs::write(&path, text).unwrap();
    path
}

/// A site file written under `dir` whose one record is `record`.
fn site_file(dir: &Path, record: &str) -> PathBuf {
    common::indexed_vcf(
        dir,
        "site.vcf.gz",
        &["##fileformat=VCFv4.2", "##contig=<ID=21>"],
        record,
    )
}

/// The REVEL plugin over the REVEL corpus's table.
fn revel_plugin() -> String {
    format!(
        "REVEL,file={}",
        corpus("GRCh38-plugin-revel")
            .join("plugin_data/revel.tsv.gz")
            .display()
    )
}

/// The source named `REVEL`, the plugin's key, over `file`, whose records the
/// variant overlaps.
fn source_named_revel(file: &Path) -> String {
    format!(
        "file={},short_name=REVEL,format=vcf,type=overlap",
        file.display()
    )
}

/// The output of a run over `input` on the cache of the corpus `cache`, under
/// `flags`; `label` names the run's output file under `dir`.
fn output(dir: &Path, label: &str, input: &Path, cache: &str, flags: &[&str]) -> String {
    let out = dir.join(format!("{label}.out"));
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus(cache).join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(flags);
    let run = cmd.output().unwrap();
    assert!(
        run.status.success(),
        "vep {label}:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    fs::read_to_string(&out).unwrap()
}

/// The header lines of `text` that open with one of `openings`, in the
/// output's order.
fn header_lines(text: &str, openings: &[&str]) -> Vec<String> {
    text.lines()
        .filter(|line| openings.iter().any(|opening| line.starts_with(opening)))
        .map(String::from)
        .collect()
}

/// The columns the tab output of `text` lists after `FLAGS`.
fn added_columns(text: &str) -> String {
    let line = text
        .lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .expect("the column line");
    line.strip_prefix(TAB_COLUMNS)
        .unwrap_or_else(|| panic!("the columns before the added ones: {line}"))
        .to_string()
}

/// The rows of `text`.
fn rows(text: &str) -> Vec<&str> {
    text.lines().filter(|line| !line.starts_with('#')).collect()
}

/// A plugin's key that a `--custom` short name repeats, on a record the plugin
/// annotates. Ensembl VEP's tab output has two `REVEL` columns, and both of
/// its `## REVEL` lines carry the source's file, also over the row whose cells
/// hold the plugin's score. The key is one column whose line carries the
/// plugin's text and then the source's, and under a `--fields` list that names
/// the key, where Ensembl VEP lists it once with the source's text, the line
/// is the same.
#[test]
fn a_plugins_key_that_a_custom_short_name_repeats_is_described_by_both_texts() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let input = record_file(dir, "GRCh38-plugin-revel", 4, "21\t26965413\t");
    let site = site_file(dir, "21\t26965413\tsite_1\tT\tC\t.\tPASS\t.");
    let plugin = revel_plugin();
    let source = source_named_revel(&site);
    let pair = ["--plugin", &plugin, "--custom", &source];
    let both_texts = format!("## REVEL : {REVEL_TEXT}; {SOURCE_TEXT}");

    let tab = output(
        dir,
        "tab",
        &input,
        "GRCh38-plugin-revel",
        &[&pair[..], &["--tab"]].concat(),
    );
    assert_eq!(header_lines(&tab, &["## REVEL "]), [both_texts.as_str()]);
    assert_eq!(added_columns(&tab), "\tSOURCE\tREVEL");
    assert_eq!(rows(&tab), SCORED_RECORD_ROWS);

    let listed = output(
        dir,
        "fields",
        &input,
        "GRCh38-plugin-revel",
        &[&pair[..], &["--tab", "--fields", "Feature,REVEL"]].concat(),
    );
    assert_eq!(header_lines(&listed, &["## REVEL "]), [both_texts.as_str()]);
    let cells: Vec<String> = SCORED_RECORD_ROWS
        .iter()
        .map(|row| {
            let columns: Vec<&str> = row.split('\t').collect();
            format!("{}\t{}", columns[4], columns[18])
        })
        .collect();
    assert_eq!(rows(&listed), cells);
}

/// The same pair on a record that overlaps no transcript, so that the plugin
/// writes neither of its rows: both hold the source's record name under the
/// one line of both texts, where both of Ensembl VEP's lines carry the
/// source's file.
#[test]
fn a_record_the_plugin_writes_no_row_of_keeps_the_line_of_both_texts() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let input = record_file(dir, "GRCh38-custom", 9, "21\t5380280\t");
    let site = site_file(dir, "21\t5380280\texact_2\tA\tC\t.\tPASS\t.");
    let plugin = revel_plugin();
    let source = source_named_revel(&site);
    let tab = output(
        dir,
        "tab",
        &input,
        "GRCh38-everything",
        &["--plugin", &plugin, "--custom", &source, "--tab"],
    );
    assert_eq!(
        header_lines(&tab, &["## REVEL "]),
        [format!("## REVEL : {REVEL_TEXT}; {SOURCE_TEXT}")]
    );
    assert_eq!(added_columns(&tab), "\tSOURCE\tREVEL");
    assert_eq!(rows(&tab), UNSCORED_RECORD_ROWS);
}

/// The headers that print a line for each declaration, which both engines
/// write alike on the scored record: the default format's two key lines and
/// the VCF header's plugin line and `##INFO` line, the plugin's before the
/// source's.
#[test]
fn the_default_and_vcf_headers_keep_a_line_for_each_declaration() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let input = record_file(dir, "GRCh38-plugin-revel", 4, "21\t26965413\t");
    let site = site_file(dir, "21\t26965413\tsite_1\tT\tC\t.\tPASS\t.");
    let plugin = revel_plugin();
    let source = source_named_revel(&site);
    let pair = ["--plugin", &plugin, "--custom", &source];

    let default = output(dir, "default", &input, "GRCh38-plugin-revel", &pair);
    assert_eq!(
        header_lines(&default, &["## REVEL "]),
        [
            format!("## REVEL : {REVEL_TEXT}"),
            format!("## REVEL : {SOURCE_TEXT}")
        ]
    );

    let vcf = output(
        dir,
        "vcf",
        &input,
        "GRCh38-plugin-revel",
        &[&pair[..], &["--vcf"]].concat(),
    );
    assert_eq!(
        header_lines(&vcf, &["##REVEL=", "##INFO=<ID=REVEL,"]),
        [
            format!("##REVEL={REVEL_TEXT}"),
            format!("##INFO=<ID=REVEL,Number=.,Type=String,Description=\"{SOURCE_TEXT}\">")
        ]
    );
}

/// A key two plugins write, on the fixture of TranscriptAnnotator under
/// `prefix=CADD` with a column headed `PHRED` and the CADD plugin. Ensembl
/// VEP's tab output has two `CADD_PHRED` columns and both of its lines carry
/// the text of the plugin loaded later, whose value stands on a row both
/// write. The key is one column whose line carries the later plugin's text
/// and then the earlier one's, in either load order.
#[test]
fn a_key_two_plugins_write_is_described_by_the_later_plugin_first() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let fixture = shared_key_fixture();
    let table = fixture.join("transcript_annotator_PHRED.tsv.gz");
    let table_text = format!("column from {}", table.display());
    let annotator = format!("TranscriptAnnotator,file={},prefix=CADD", table.display());
    let cadd = format!(
        "CADD,{}",
        fixture
            .join("CADD_shared_key_GRCh38_synthetic.tsv.gz")
            .display()
    );
    let input = fixture.join("variants.vcf");
    let lines = |label: &str, plugins: [&str; 2]| {
        let tab = output(
            dir,
            label,
            &input,
            "GRCh38-plugin-transcript-annotator",
            &["--plugin", plugins[0], "--plugin", plugins[1], "--tab"],
        );
        assert_eq!(added_columns(&tab), "\tCADD_PHRED\tCADD_RAW", "{label}");
        header_lines(&tab, &["## CADD_PHRED "])
    };
    assert_eq!(
        lines("table-then-cadd", [&annotator, &cadd]),
        [format!("## CADD_PHRED : {CADD_TEXT}; {table_text}")]
    );
    assert_eq!(
        lines("cadd-then-table", [&cadd, &annotator]),
        [format!("## CADD_PHRED : {table_text}; {CADD_TEXT}")]
    );
}

/// The same fixture with the CADD plugin, loaded last, over a table that
/// holds no row at the record, the SNV table of the CADD fixture corpus. CADD
/// writes no row, so the one cell that holds a value is the table's, and both
/// of Ensembl VEP's lines carry CADD's text. The line names both plugins as
/// it does where CADD scores the record: a header is written before any row.
#[test]
fn a_plugin_that_writes_no_row_of_the_record_is_still_named_on_the_line() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = shared_key_fixture();
    let table = fixture.join("transcript_annotator_PHRED.tsv.gz");
    let annotator = format!("TranscriptAnnotator,file={},prefix=CADD", table.display());
    let cadd = format!(
        "CADD,{}",
        corpus("GRCh38-plugin-cadd-fixture")
            .join("plugin_data/whole_genome_SNVs_synthetic.tsv.gz")
            .display()
    );
    let tab = output(
        tmp.path(),
        "table-then-cadd-with-no-row",
        &fixture.join("variants.vcf"),
        "GRCh38-plugin-transcript-annotator",
        &["--plugin", &annotator, "--plugin", &cadd, "--tab"],
    );
    assert_eq!(
        header_lines(&tab, &["## CADD_PHRED "]),
        [format!(
            "## CADD_PHRED : {CADD_TEXT}; column from {}",
            table.display()
        )]
    );
    assert_eq!(added_columns(&tab), "\tCADD_PHRED\tCADD_RAW");
    let cells: Vec<(&str, &str, &str)> = rows(&tab)
        .into_iter()
        .map(|row| {
            let columns: Vec<&str> = row.split('\t').collect();
            (columns[4], columns[17], columns[18])
        })
        .collect();
    assert_eq!(
        cells,
        [
            ("ENST00000427445", "-", "-"),
            ("ENST00000612746", "-", "-"),
            ("ENST00000612957", "-", "-"),
            ("ENST00000618007", "7.5", "-"),
            ("ENST00000622113", "-", "-"),
            ("ENST00001091160", "-", "-"),
            ("ENST00001091161", "-", "-"),
            ("ENST00001091162", "-", "-"),
            ("ENST00001131865", "-", "-"),
        ]
    );
}

/// The two flags under which no plugin runs on a row: `--most_severe` and
/// `--summary` write the record's hash alone (`OutputFactory.pm:387`), so the
/// one cell under the key holds the source's record name and the plugin's text
/// describes no cell of the run. The line is the source's alone, the text of
/// both of Ensembl VEP's lines, over the row Ensembl VEP writes.
#[test]
fn a_run_of_summary_rows_describes_the_key_by_the_source_alone() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let input = record_file(dir, "GRCh38-plugin-revel", 4, "21\t26965413\t");
    let site = site_file(dir, "21\t26965413\tsite_1\tT\tC\t.\tPASS\t.");
    let plugin = revel_plugin();
    let source = source_named_revel(&site);
    for (flag, terms) in [
        ("--most_severe", "missense_variant"),
        ("--summary", "missense_variant,NMD_transcript_variant"),
    ] {
        let tab = output(
            dir,
            flag.trim_start_matches('-'),
            &input,
            "GRCh38-plugin-revel",
            &["--plugin", &plugin, "--custom", &source, "--tab", flag],
        );
        assert_eq!(
            header_lines(&tab, &["## REVEL "]),
            [format!("## REVEL : {SOURCE_TEXT}")],
            "{flag}"
        );
        assert_eq!(added_columns(&tab), "\tSOURCE\tREVEL", "{flag}");
        assert_eq!(
            rows(&tab),
            [format!(
                "4637665\t21:26965413\t-\t-\t-\t-\t{terms}{}\tsite_1",
                "\t-".repeat(11)
            )],
            "{flag}"
        );
    }
}

/// Two plugins that write one key under `--most_severe`: neither runs on the
/// record's one row, whose cell is empty, and the line keeps the text of the
/// plugin loaded last, which both of Ensembl VEP's lines carry.
#[test]
fn two_plugins_under_a_summary_flag_keep_the_text_of_the_plugin_loaded_last() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let fixture = shared_key_fixture();
    let table = fixture.join("transcript_annotator_PHRED.tsv.gz");
    let table_text = format!("column from {}", table.display());
    let annotator = format!("TranscriptAnnotator,file={},prefix=CADD", table.display());
    let cadd = format!(
        "CADD,{}",
        fixture
            .join("CADD_shared_key_GRCh38_synthetic.tsv.gz")
            .display()
    );
    let input = fixture.join("variants.vcf");
    let lines = |label: &str, plugins: [&str; 2]| {
        let tab = output(
            dir,
            label,
            &input,
            "GRCh38-plugin-transcript-annotator",
            &[
                "--plugin",
                plugins[0],
                "--plugin",
                plugins[1],
                "--tab",
                "--most_severe",
            ],
        );
        assert_eq!(added_columns(&tab), "\tCADD_PHRED\tCADD_RAW", "{label}");
        assert_eq!(
            rows(&tab),
            [format!(
                "shared_key\t21:10521702\t-\t-\t-\t-\tsplice_region_variant{}",
                "\t-".repeat(12)
            )],
            "{label}"
        );
        header_lines(&tab, &["## CADD_PHRED "])
    };
    assert_eq!(
        lines("table-then-cadd", [&annotator, &cadd]),
        [format!("## CADD_PHRED : {CADD_TEXT}")]
    );
    assert_eq!(
        lines("cadd-then-table", [&cadd, &annotator]),
        [format!("## CADD_PHRED : {table_text}")]
    );
}
