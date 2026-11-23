// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary end to end: directories of inputs and outputs written into a scratch
//! directory, the command line run on them, and the counts, report text, stderr lines and
//! exit codes it gives; then the committed fixtures under `tests/fixtures`, whose reports,
//! stdout, stderr and exit code are held to their goldens byte for byte at three thread
//! counts.

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::{json, Value};
use vep_sv_compare::TempDir;

/// The two header lines of every input VCF written here.
const VCF_HEAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The column line of an output of seven columns.
const SEVEN_COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";

/// The four report files of a comparison.
const REPORT_FILES: [&str; 4] = [
    "concordance_report.json",
    "concordance_report.md",
    "discordant.tsv",
    "discordant_adjusted.tsv",
];

/// The first seven columns of an output row: the first column, Location, Allele, Gene,
/// Feature, Feature_type and Consequence.
type Row = [String; 7];

/// A Transcript row of gene `G`.
fn row(uploaded: &str, location: &str, allele: &str, feature: &str, consequence: &str) -> Row {
    [
        uploaded,
        location,
        allele,
        "G",
        feature,
        "Transcript",
        consequence,
    ]
    .map(str::to_string)
}

/// An output of seventeen columns: a `##` line, the column line, and every row with ten
/// columns of `-` after its seven.
fn vep_output(rows: &[Row]) -> String {
    let mut text = String::from(
        "## fastVEP output\n\
         #Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t\
         Consequence\tcDNA_position\tCDS_position\tProtein_position\t\
         Amino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\n",
    );
    for row in rows {
        text.push_str(&row.join("\t"));
        text.push_str(&"\t-".repeat(10));
        text.push('\n');
    }
    text
}

/// A scratch directory holding the directories of one comparison: `in`, `perl`, `rust`,
/// `out` and, where a test writes one, `cache`.
struct Fixture {
    dir: TempDir,
}

impl Fixture {
    fn new() -> Fixture {
        let fixture = Fixture {
            dir: TempDir::new("vep_sv_compare_cli_").expect("a scratch directory"),
        };
        for directory in ["in", "perl", "rust", "out"] {
            std::fs::create_dir(fixture.path(directory)).expect("a directory of the fixture");
        }
        fixture
    }

    fn path(&self, name: &str) -> PathBuf {
        self.dir.path.join(name)
    }

    /// Writes a file, a name ending in `.gz` gzipped, making the directories above it.
    fn write(&self, name: &str, text: &str) {
        let path = self.path(name);
        std::fs::create_dir_all(path.parent().expect("a file under the fixture"))
            .expect("the directory of a fixture file");
        let bytes = if name.ends_with(".gz") {
            let mut encoder =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            encoder.write_all(text.as_bytes()).expect("gzip in memory");
            encoder.finish().expect("gzip in memory")
        } else {
            text.as_bytes().to_vec()
        };
        std::fs::write(path, bytes).expect("a fixture file");
    }

    /// Writes an input VCF of the given data lines under `in`.
    fn write_vcf(&self, name: &str, records: &[&str]) {
        let mut text = String::from(VCF_HEAD);
        for record in records {
            text.push_str(record);
            text.push('\n');
        }
        self.write(&format!("in/{name}"), &text);
    }

    /// Writes one output of seventeen columns under each of `perl` and `rust`.
    fn write_outputs(&self, name: &str, perl: &[Row], rust: &[Row]) {
        self.write(&format!("perl/{name}.txt"), &vep_output(perl));
        self.write(&format!("rust/{name}.txt"), &vep_output(rust));
    }

    /// Writes a cache of one shard a chromosome that holds the given transcript ids.
    fn write_cache(&self, transcripts_by_chromosome: &[(&str, &[&str])]) {
        for (chromosome, ids) in transcripts_by_chromosome {
            let records: Vec<String> = ids
                .iter()
                .map(|id| format!("{{\"stable_id\": \"{id}\"}}"))
                .collect();
            self.write(
                &format!("cache/transcripts/{chromosome}/0-1.json"),
                &format!("[{}]", records.join(", ")),
            );
        }
    }

    fn read(&self, name: &str) -> String {
        read(&self.path(name))
    }

    /// The arguments of a comparison of this fixture's directories.
    fn arguments(&self, release: &str, assembly: &str) -> Vec<OsString> {
        let mut arguments = strings(&["--reference-release", release, "--assembly", assembly]);
        for (option, directory) in [
            ("--input-dir", "in"),
            ("--perl-dir", "perl"),
            ("--rust-dir", "rust"),
            ("--output-dir", "out"),
        ] {
            arguments.push(option.into());
            arguments.push(self.path(directory).into());
        }
        arguments
    }

    /// Runs the comparison of this fixture's directories with `extra` after its arguments.
    fn compare(&self, release: &str, assembly: &str, extra: &[OsString]) -> Finished {
        let mut arguments = self.arguments(release, assembly);
        arguments.extend_from_slice(extra);
        run(&arguments)
    }

    /// The JSON report of a comparison of this fixture's directories that exits 0.
    fn scored(&self, release: &str, assembly: &str, extra: &[OsString]) -> Value {
        let finished = self.compare(release, assembly, extra);
        assert_eq!(finished.code, Some(0), "stderr: {}", finished.stderr);
        self.report()
    }

    /// The JSON report under `out`.
    fn report(&self) -> Value {
        report(&self.path("out"))
    }

    /// The lines of a file under `out`.
    fn lines(&self, name: &str) -> Vec<String> {
        self.read(&format!("out/{name}"))
            .lines()
            .map(str::to_string)
            .collect()
    }
}

fn strings(texts: &[&str]) -> Vec<OsString> {
    texts.iter().map(OsString::from).collect()
}

fn read(path: &Path) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The JSON report of an output directory.
fn report(output_dir: &Path) -> Value {
    serde_json::from_str(&read(&output_dir.join("concordance_report.json")))
        .expect("the report is JSON")
}

/// How a run of the binary ended.
struct Finished {
    code: Option<i32>,
    stdout: String,
    stderr: String,
}

fn tool() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vep-sv-compare"))
}

fn finish(mut command: Command) -> Finished {
    let output = command.output().expect("the binary runs");
    Finished {
        code: output.status.code(),
        stdout: String::from_utf8(output.stdout).expect("stdout is UTF-8"),
        stderr: String::from_utf8(output.stderr).expect("stderr is UTF-8"),
    }
}

fn run(arguments: &[OsString]) -> Finished {
    let mut command = tool();
    command.args(arguments);
    finish(command)
}

/// The value under the keys of `path`; a key the report lacks fails the test.
fn at<'a>(value: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(value, |value, key| {
        value
            .get(key)
            .unwrap_or_else(|| panic!("the report has no key {key} under {path:?}"))
    })
}

fn count(value: &Value, path: &[&str]) -> i64 {
    at(value, path)
        .as_i64()
        .unwrap_or_else(|| panic!("{path:?} is no integer"))
}

fn ratio(value: &Value, path: &[&str]) -> f64 {
    at(value, path)
        .as_f64()
        .unwrap_or_else(|| panic!("{path:?} is no number"))
}

fn text<'a>(value: &'a Value, path: &[&str]) -> &'a str {
    at(value, path)
        .as_str()
        .unwrap_or_else(|| panic!("{path:?} is no text"))
}

fn object<'a>(value: &'a Value, path: &[&str]) -> &'a serde_json::Map<String, Value> {
    at(value, path)
        .as_object()
        .unwrap_or_else(|| panic!("{path:?} is no object"))
}

/// The names of an array of texts, in its order.
fn names<'a>(value: &'a Value, path: &[&str]) -> Vec<&'a str> {
    at(value, path)
        .as_array()
        .unwrap_or_else(|| panic!("{path:?} is no array"))
        .iter()
        .map(|name| name.as_str().expect("a name"))
        .collect()
}

mod aggregate_basis_end_to_end_tests {
    use super::*;

    /// Two files that both hold one tuple on each side; the first also holds a tuple of the
    /// reference alone and a tuple of the scored engine alone. File-summed, each side holds
    /// three tuples and the intersection two; as sets, two and one.
    fn report() -> Value {
        let fixture = Fixture::new();
        let truncation = |uploaded: &str, location: &str, gene: &str, feature: &str| -> Row {
            [
                uploaded,
                location,
                "<DEL>",
                gene,
                feature,
                "Transcript",
                "feature_truncation",
            ]
            .map(str::to_string)
        };
        let dup = truncation("dupvar", "21:500", "G0", "ENST_DUP");
        let perl_only = truncation("pv", "21:600", "G1", "ENST_P");
        let rust_only: Row = [
            "rv",
            "21:700",
            "<DUP>",
            "G2",
            "ENST_R",
            "Transcript",
            "feature_elongation",
        ]
        .map(str::to_string);
        for (name, extra_perl, extra_rust) in [
            ("fileA", vec![perl_only], vec![rust_only]),
            ("fileB", vec![], vec![]),
        ] {
            fixture.write_vcf(
                &format!("{name}.vcf"),
                &["21\t500\t.\tG\t<DEL>\t.\tPASS\tEND=600"],
            );
            let with = |extra: Vec<Row>| [vec![dup.clone()], extra].concat();
            fixture.write_outputs(name, &with(extra_perl), &with(extra_rust));
        }
        fixture.scored("115.2", "grch37", &[])
    }

    #[test]
    fn test_overall_uses_file_summed_counts() {
        let report = report();
        assert_eq!(count(&report, &["overall", "perl_tuples"]), 3);
        assert_eq!(count(&report, &["overall", "rust_tuples"]), 3);
        assert_eq!(count(&report, &["overall", "intersection"]), 2);
    }

    #[test]
    fn test_overall_reports_the_union_alongside() {
        let report = report();
        assert_eq!(count(&report, &["overall", "perl_tuples_union"]), 2);
        assert_eq!(count(&report, &["overall", "rust_tuples_union"]), 2);
        assert_eq!(count(&report, &["overall", "intersection_union"]), 1);
        assert_eq!(count(&report, &["overall", "cross_file_duplicate_perl"]), 1);
        assert_eq!(count(&report, &["overall", "cross_file_duplicate_rust"]), 1);
    }

    #[test]
    fn test_adjusted_uses_the_same_basis_as_raw() {
        let report = report();
        let counts = |row: &str| {
            ["perl_tuples", "rust_tuples", "intersection"].map(|key| count(&report, &[row, key]))
        };
        assert_eq!(counts("adjusted"), counts("overall"));
    }

    #[test]
    fn test_adjusted_f1_equals_raw_f1_when_nothing_is_excluded() {
        let report = report();
        let adjusted = object(&report, &["adjusted"]);
        assert_eq!(count(&report, &["adjusted", "excluded_divergences"]), 0);
        assert_eq!(count(&report, &["adjusted", "excluded_cnvtr_rust"]), 0);
        assert!(!adjusted.contains_key("excluded_csq_rust"));
        assert!(!adjusted.contains_key("excluded_csq_perl"));
        let difference = ratio(&report, &["adjusted", "f1"]) - ratio(&report, &["overall", "f1"]);
        assert_eq!((difference.abs() * 1e6).round(), 0.0);
    }

    #[test]
    fn test_adjusted_also_reports_its_union_counts() {
        let report = report();
        let adjusted = object(&report, &["adjusted"]);
        for key in [
            "perl_tuples_union",
            "rust_tuples_union",
            "intersection_union",
            "cross_file_duplicate_perl",
            "cross_file_duplicate_rust",
        ] {
            assert!(adjusted.contains_key(key), "adjusted row is missing {key}");
        }
        assert_eq!(
            count(&report, &["adjusted", "cross_file_duplicate_perl"]),
            1
        );
    }

    #[test]
    fn test_per_file_rows_are_unchanged_by_the_aggregate_override() {
        let report = report();
        for section in ["per_file", "per_file_adjusted"] {
            let rows = object(&report, &[section]);
            assert_eq!(rows.keys().collect::<Vec<_>>(), ["fileA", "fileB"]);
            for (name, row) in rows {
                let expected = if name == "fileA" { 2 } else { 1 };
                assert_eq!(count(row, &["perl_tuples"]), expected, "{section}/{name}");
                assert!(
                    row.get("perl_tuples_union").is_none(),
                    "{section}/{name} carries a union count"
                );
            }
        }
    }
}

mod fastvep_engine_end_to_end_tests {
    use super::*;

    /// One input of a deletion and a mobile-element insertion: the reference's output under
    /// the allele classes, the other output under the symbolic alleles as written.
    fn fixture() -> Fixture {
        let fixture = Fixture::new();
        fixture.write_vcf(
            "05_symbolic_del_ins.vcf",
            &[
                "21\t1000\t.\tG\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000",
                "21\t3000\t.\tG\t<INS:ME:ALU>\t.\tPASS\tSVTYPE=INS",
            ],
        );
        let rows = |deletion: &str, insertion: &str| {
            [
                row(
                    "v",
                    "21:1001-2000",
                    deletion,
                    "ENST_A",
                    "feature_truncation",
                ),
                row("v", "21:3001", insertion, "ENST_B", "intron_variant"),
            ]
        };
        fixture.write_outputs(
            "05_symbolic_del_ins",
            &rows("deletion", "Alu_insertion"),
            &rows("<DEL>", "<INS:ME:ALU>"),
        );
        fixture
    }

    #[test]
    fn test_scored_engine_fastvep_makes_raw_tokens_match() {
        let report = fixture().scored("115.2", "grch38", &strings(&["--scored-engine", "fastvep"]));
        assert_eq!(count(&report, &["overall", "intersection"]), 2);
        assert_eq!(ratio(&report, &["overall", "f1"]), 1.0);
        assert_eq!(text(&report, &["adjusted", "scored_engine"]), "fastvep");
    }

    #[test]
    fn test_default_engine_leaves_raw_tokens_as_misses() {
        let report = fixture().scored("115.2", "grch38", &[]);
        assert_eq!(count(&report, &["overall", "intersection"]), 0);
        assert_eq!(ratio(&report, &["overall", "f1"]), 0.0);
    }

    #[test]
    fn test_the_perl_side_is_parsed_without_the_mapping() {
        let fixture = fixture();
        // The two directories named again, swapped: the last of a repeated option stands.
        let extra = [
            "--perl-dir".into(),
            fixture.path("rust").into(),
            "--rust-dir".into(),
            fixture.path("perl").into(),
            "--scored-engine".into(),
            "fastvep".into(),
        ];
        let report = fixture.scored("115.2", "grch38", &extra);
        assert_eq!(count(&report, &["overall", "intersection"]), 0);
    }
}

/// The row both outputs of the tandem-repeat file hold.
fn shared_row() -> Row {
    row(
        "v0",
        "21:100-200",
        "deletion",
        "ENST_S",
        "feature_truncation",
    )
}

/// The Location of a deletion past the size cap of the reference.
const GIANT: &str = "21:1000000-21000000";

/// The scored engine's rows on the deletion past the size cap: four transcripts.
fn giant_rust_rows() -> Vec<Row> {
    ["P", "Q", "R", "S"]
        .iter()
        .map(|c| {
            row(
                "g1",
                GIANT,
                "deletion",
                &format!("ENST_{c}"),
                "feature_truncation",
            )
        })
        .collect()
}

mod transcript_selection_own_count_end_to_end_tests {
    use super::*;

    /// One file under each scored-side layer: a deletion past the size cap on which the
    /// scored engine names three transcripts the reference lacks, and a tandem-repeat gain
    /// pair differing only in representation.
    fn fixture() -> Fixture {
        let fixture = Fixture::new();
        fixture.write_vcf(
            "07_cnv_repeat.vcf",
            &["21\t1984\t.\tA\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=16"],
        );
        fixture.write_vcf(
            "09_breakends.vcf",
            &["21\t1000000\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=21000000"],
        );
        fixture.write_outputs(
            "07_cnv_repeat",
            &[
                shared_row(),
                row("v1", "21:2000-2001", "ATAT", "ENST_A", "intron_variant"),
            ],
            &[
                shared_row(),
                row(
                    "v1",
                    "21:1985-2000",
                    "tandem_repeat",
                    "ENST_A",
                    "intron_variant,feature_elongation",
                ),
            ],
        );
        fixture.write_outputs(
            "09_breakends",
            &[row("g1", GIANT, "deletion", "ENST_P", "feature_truncation")],
            &giant_rust_rows(),
        );
        fixture
    }

    #[test]
    fn test_own_count_sits_beside_the_union() {
        let report = fixture().scored("115.2", "grch37", &[]);
        let adjusted = |key: &str| count(&report, &["adjusted", key]);
        assert_eq!(adjusted("excluded_transcript_selection_rust"), 3);
        assert_eq!(adjusted("excluded_transcript_selection_perl"), 0);
        assert_eq!(adjusted("excluded_cnvtr_rust"), 1);
        assert_eq!(adjusted("transcript_selection_cnvtr_overlap_rust"), 0);
        assert_eq!(
            adjusted("excluded_rust_extra"),
            adjusted("excluded_transcript_selection_rust") + adjusted("excluded_cnvtr_rust")
                - adjusted("transcript_selection_cnvtr_overlap_rust")
        );
    }

    #[test]
    fn test_per_file_split_names_only_the_files_the_layer_touched() {
        let report = fixture().scored("115.2", "grch37", &[]);
        assert_eq!(
            at(
                &report,
                &["adjusted", "excluded_transcript_selection_by_file"]
            ),
            &json!({"09_breakends": {"rust": 3, "perl": 0}})
        );
    }

    #[test]
    fn test_own_count_keys_are_present_when_the_layer_excludes_nothing() {
        let fixture = fixture();
        fixture.scored("115.2", "grch37", &[]);
        fixture.write_outputs("09_breakends", &[shared_row()], &[shared_row()]);
        let report = fixture.scored("115.2", "grch37", &[]);
        let adjusted = |key: &str| at(&report, &["adjusted", key]);
        assert_eq!(adjusted("excluded_transcript_selection_rust"), &json!(0));
        assert_eq!(adjusted("excluded_transcript_selection_perl"), &json!(0));
        assert_eq!(
            adjusted("excluded_transcript_selection_by_file"),
            &json!({})
        );
    }
}

/// The reference's literal reading and the scored engine's symbolic reading of two
/// tandem-repeat records: a gain that differs only in representation, and a loss whose
/// literal reading carries a start term the symbolic one lacks.
fn tandem_repeat_rows() -> (Vec<Row>, Vec<Row>) {
    let tr = "tandem_repeat";
    (
        vec![
            row("v1", "21:2000-2001", "ATAT", "ENST_A", "intron_variant"),
            row(
                "v2",
                "21:2997-3000",
                "-",
                "ENST_B",
                "intron_variant,start_lost",
            ),
        ],
        vec![
            row(
                "v1",
                "21:1985-2000",
                tr,
                "ENST_A",
                "intron_variant,feature_elongation",
            ),
            row(
                "v2",
                "21:2985-3000",
                tr,
                "ENST_B",
                "intron_variant,feature_truncation",
            ),
        ],
    )
}

/// The input of the two tandem-repeat records.
fn write_tandem_repeat_vcf(fixture: &Fixture) {
    fixture.write_vcf(
        "07_cnv_repeat.vcf",
        &[
            "21\t1984\t.\tA\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=16",
            "21\t2984\t.\tA\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=16",
        ],
    );
}

mod cnv_tr_total_report_end_to_end_tests {
    use super::*;

    /// The two tandem-repeat pairs beside one tuple both sides hold.
    fn fixture() -> Fixture {
        let fixture = Fixture::new();
        write_tandem_repeat_vcf(&fixture);
        let (perl, rust) = tandem_repeat_rows();
        fixture.write_outputs(
            "07_cnv_repeat",
            &[vec![shared_row()], perl].concat(),
            &[vec![shared_row()], rust].concat(),
        );
        fixture
    }

    #[test]
    fn test_json_carries_the_class_total_beside_the_excluded_count() {
        let report = fixture().scored("115.2", "grch37", &[]);
        assert_eq!(count(&report, &["adjusted", "excluded_cnvtr_rust"]), 1);
        assert_eq!(count(&report, &["adjusted", "excluded_cnvtr_perl"]), 1);
        assert_eq!(
            count(&report, &["adjusted", "cnvtr_swap_pairs_total_rust"]),
            2
        );
        assert_eq!(
            count(&report, &["adjusted", "cnvtr_swap_pairs_total_perl"]),
            2
        );
    }

    #[test]
    fn test_json_carries_the_per_file_split() {
        let report = fixture().scored("115.2", "grch37", &[]);
        assert_eq!(
            at(&report, &["adjusted", "cnvtr_swap_pairs_total_by_file"]),
            &json!({"07_cnv_repeat": {"rust": 2, "perl": 2}})
        );
    }

    #[test]
    fn test_markdown_caption_states_n_of_m() {
        let fixture = fixture();
        fixture.scored("115.2", "grch37", &[]);
        let markdown = fixture.read("out/concordance_report.md");
        assert!(markdown.contains("(1 of 2 pairs GRCH37, 07_cnv_repeat only)"));
    }

    #[test]
    fn test_total_keys_are_present_even_when_zero() {
        let fixture = fixture();
        fixture.scored("115.2", "grch37", &[]);
        fixture.write_outputs("07_cnv_repeat", &[shared_row()], &[shared_row()]);
        let report = fixture.scored("115.2", "grch37", &[]);
        let adjusted = |key: &str| at(&report, &["adjusted", key]);
        assert_eq!(adjusted("cnvtr_swap_pairs_total_rust"), &json!(0));
        assert_eq!(adjusted("cnvtr_swap_pairs_total_perl"), &json!(0));
        assert_eq!(adjusted("cnvtr_swap_pairs_total_by_file"), &json!({}));
    }
}

mod adjusted_discordant_tsv_end_to_end_tests {
    use super::*;

    /// The Location of a deletion under the size cap.
    const SMALL: &str = "21:5000-6000";

    /// What one comparison of the fixture leaves: the fixture, its report and the lines of
    /// the two tables of one-sided tuples.
    struct Scored {
        fixture: Fixture,
        report: Value,
        raw_lines: Vec<String>,
        adj_lines: Vec<String>,
    }

    /// One tuple under each layer and three outside every layer. On the deletion past the
    /// size cap the scored engine names three transcripts the reference lacks and the
    /// reference names one transcript of another chromosome; on the deletion under the cap
    /// the reference names a transcript the scored engine lacks, which stays. Of the two
    /// tandem-repeat pairs the gain is set aside and the loss stays on both sides.
    fn scored() -> Scored {
        let fixture = Fixture::new();
        fixture.write_cache(&[
            (
                "21",
                &[
                    "ENST_P", "ENST_Q", "ENST_R", "ENST_S", "ENST_GAP", "ENST_A", "ENST_B",
                ],
            ),
            ("22", &["ENST_CHR22"]),
        ]);
        write_tandem_repeat_vcf(&fixture);
        fixture.write_vcf(
            "09_breakends.vcf",
            &[
                "21\t1000000\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=21000000",
                "21\t5000\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=6000",
            ],
        );
        let (perl, rust) = tandem_repeat_rows();
        fixture.write_outputs("07_cnv_repeat", &perl, &rust);
        let deletion = |uploaded: &str, location: &str, feature: &str| {
            row(
                uploaded,
                location,
                "deletion",
                feature,
                "feature_truncation",
            )
        };
        fixture.write_outputs(
            "09_breakends",
            &[
                deletion("g1", GIANT, "ENST_P"),
                deletion("g1", GIANT, "ENST_CHR22"),
                deletion("g2", SMALL, "ENST_GAP"),
            ],
            &giant_rust_rows(),
        );
        let extra = ["--vep-rs-cache".into(), fixture.path("cache").into()];
        let report = fixture.scored("115.2", "grch37", &extra);
        let raw_lines = fixture.lines("discordant.tsv");
        let adj_lines = fixture.lines("discordant_adjusted.tsv");
        Scored {
            fixture,
            report,
            raw_lines,
            adj_lines,
        }
    }

    /// The direction, Location, Allele and Feature of every data line.
    fn keys(lines: &[String]) -> BTreeSet<[String; 4]> {
        lines[1..]
            .iter()
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                [columns[2], columns[3], columns[4], columns[5]].map(str::to_string)
            })
            .collect()
    }

    fn key(direction: &str, location: &str, allele: &str, feature: &str) -> [String; 4] {
        [direction, location, allele, feature].map(str::to_string)
    }

    #[test]
    fn test_every_layer_fired_on_the_fixture() {
        let report = scored().report;
        let adjusted = |key: &str| count(&report, &["adjusted", key]);
        assert_eq!(adjusted("excluded_transcript_selection_rust"), 3);
        assert_eq!(adjusted("excluded_cross_chromosome_perl"), 1);
        assert_eq!(
            (
                adjusted("excluded_cnvtr_rust"),
                adjusted("excluded_cnvtr_perl")
            ),
            (1, 1)
        );
    }

    #[test]
    fn test_same_header_as_the_raw_file() {
        let scored = scored();
        assert_eq!(scored.adj_lines[0], scored.raw_lines[0]);
    }

    #[test]
    fn test_row_count_is_the_adjusted_one_sided_total() {
        let scored = scored();
        let one_sided = |row: &Value| count(row, &["only_perl"]) + count(row, &["only_rust"]);
        let rows = scored.adj_lines.len() as i64 - 1;
        assert_eq!(rows, one_sided(at(&scored.report, &["adjusted"])));
        assert_eq!(
            rows,
            object(&scored.report, &["per_file_adjusted"])
                .values()
                .map(one_sided)
                .sum::<i64>()
        );
    }

    #[test]
    fn test_raw_file_still_carries_every_one_sided_tuple() {
        let scored = scored();
        assert_eq!(
            scored.raw_lines.len() as i64 - 1,
            count(&scored.report, &["overall", "only_perl"])
                + count(&scored.report, &["overall", "only_rust"])
        );
        assert!(scored.raw_lines.len() > scored.adj_lines.len());
    }

    #[test]
    fn test_rows_are_a_subset_of_the_raw_rows_in_raw_order() {
        let scored = scored();
        let raw_index: BTreeMap<&str, usize> = scored
            .raw_lines
            .iter()
            .enumerate()
            .map(|(index, line)| (line.as_str(), index))
            .collect();
        let positions: Vec<usize> = scored
            .adj_lines
            .iter()
            .map(|line| raw_index[line.as_str()])
            .collect();
        assert!(positions.is_sorted(), "{positions:?}");
    }

    #[test]
    fn test_exactly_the_unmasked_tuples_survive() {
        let scored = scored();
        let survivors = BTreeSet::from([
            key("only_perl", SMALL, "deletion", "ENST_GAP"),
            key("only_perl", "21:2997-3000", "-", "ENST_B"),
            key("only_rust", "21:2985-3000", "tandem_repeat", "ENST_B"),
        ]);
        assert_eq!(keys(&scored.adj_lines), survivors);
        let mut masked: BTreeSet<[String; 4]> = ["Q", "R", "S"]
            .iter()
            .map(|c| key("only_rust", GIANT, "deletion", &format!("ENST_{c}")))
            .collect();
        masked.extend([
            key("only_perl", GIANT, "deletion", "ENST_CHR22"),
            key("only_perl", "21:2000-2001", "ATAT", "ENST_A"),
            key("only_rust", "21:1985-2000", "tandem_repeat", "ENST_A"),
        ]);
        let all: BTreeSet<[String; 4]> = survivors.union(&masked).cloned().collect();
        assert_eq!(keys(&scored.raw_lines), all);
    }

    #[test]
    fn test_file_is_written_empty_with_a_header_when_nothing_survives() {
        let Scored {
            fixture, raw_lines, ..
        } = scored();
        fixture.write_outputs("07_cnv_repeat", &[], &[]);
        let rows = [row("g1", GIANT, "deletion", "ENST_P", "feature_truncation")];
        fixture.write_outputs("09_breakends", &rows, &rows);
        fixture.scored("115.2", "grch37", &[]);
        assert_eq!(
            fixture.lines("discordant_adjusted.tsv"),
            [raw_lines[0].clone()]
        );
    }
}

/// The directory of the committed fixtures.
fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// Copies a directory with everything under it.
fn copy_tree(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).expect("the directory copied to");
    for entry in std::fs::read_dir(from).expect("the directory copied from") {
        let entry = entry.expect("an entry of the directory copied from");
        let target = to.join(entry.file_name());
        if entry.file_type().expect("the entry's type").is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), target).expect("a file copied");
        }
    }
}

/// The number of line feeds in a text.
fn line_feeds(text: &str) -> usize {
    text.matches('\n').count()
}

mod reference_release_fixture_tests {
    use super::*;
    use vep_sv_compare::masks::{MASKS_116_2, MASKS_RETIRED_IN_116_2};

    /// The committed fixture: `perl/` is a reference output of the 115.2 shape, with rows
    /// for a `<CPX>` and a `<NON_REF>` record, and `perl_116/` the same output as 116.2
    /// writes it, the two skipped records dropped.
    fn committed() -> PathBuf {
        fixtures().join("reference_release")
    }

    /// Runs the comparison of the committed fixture against the reference output under
    /// `perl`, with the fixture's cache: the output directory and how the run ended.
    fn run_release(release: &str, perl: &str) -> (TempDir, Finished) {
        let committed = committed();
        let out = TempDir::new("vep_sv_compare_cli_").expect("a scratch directory");
        let mut arguments = strings(&["--reference-release", release, "--assembly", "grch37"]);
        for (option, path) in [
            ("--input-dir", committed.join("in")),
            ("--perl-dir", committed.join(perl)),
            ("--rust-dir", committed.join("rust")),
            ("--vep-rs-cache", committed.join("cache")),
            ("--output-dir", out.path.clone()),
        ] {
            arguments.push(option.into());
            arguments.push(path.into());
        }
        let finished = run(&arguments);
        (out, finished)
    }

    /// A fixture of its own inputs and outputs beside a copy of the committed cache.
    fn with_committed_cache() -> (Fixture, [OsString; 2]) {
        let fixture = Fixture::new();
        copy_tree(&committed().join("cache"), &fixture.path("cache"));
        let extra = ["--vep-rs-cache".into(), fixture.path("cache").into()];
        (fixture, extra)
    }

    #[test]
    fn test_115_2_reproduces_the_published_report_byte_for_byte() {
        let (out, finished) = run_release("115.2", "perl");
        assert_eq!(finished.code, Some(0), "stderr: {}", finished.stderr);
        for name in REPORT_FILES {
            assert_eq!(
                read(&out.path.join(name)),
                read(&committed().join("expected_115_2").join(name)),
                "{name}"
            );
        }
    }

    #[test]
    fn test_116_2_sets_aside_the_skipped_records_whole_and_leaves_the_retired_classes_charged() {
        let (out, finished) = run_release("116.2", "perl_116");
        assert_eq!(finished.code, Some(0), "stderr: {}", finished.stderr);
        let report = report(&out.path);
        let adjusted = at(&report, &["adjusted"]);
        assert_eq!(text(adjusted, &["reference_release"]), "116.2");
        assert_eq!(names(adjusted, &["masks_applied"]), MASKS_116_2);
        assert_eq!(
            object(adjusted, &["masks_not_applied"])
                .keys()
                .map(String::as_str)
                .collect::<BTreeSet<_>>(),
            MASKS_RETIRED_IN_116_2
                .iter()
                .map(|(mask, _)| *mask)
                .collect::<BTreeSet<_>>()
        );
        // The `<CPX>` the reference dropped (2 rows), the gVCF block (1), the `<CPX>` with
        // no reference row (2) and the incomplete deletion (1): four records, six rows.
        assert_eq!(
            (
                count(adjusted, &["excluded_reference_skipped_rust"]),
                count(adjusted, &["reference_skipped_records"])
            ),
            (6, 4)
        );
        assert_eq!(
            at(adjusted, &["excluded_reference_skipped_by_file"]),
            &json!({"sv": {"rust": 6, "records": 4}})
        );
        for key in [
            "excluded_cross_chromosome_perl",
            "excluded_cross_chromosome_rust",
            "excluded_transcript_selection_rust",
            "excluded_non_ref_batch_perl",
        ] {
            assert_eq!(count(adjusted, &[key]), 0, "{key}");
        }
        assert_eq!(at(adjusted, &["breakend_masks_check_ran"]), &json!(true));
        assert_eq!(
            object(adjusted, &["divergence_classes"])
                .keys()
                .map(String::as_str)
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([
                "reference_skipped_record",
                "bnd_own_chromosome_rows_lost",
                "bnd_synonym_mate_rows_lost",
                "fasta_named_slice_intergenic",
                "bnd_bare_reference_allele_row",
                "bnd_mate_32bit_buffer_intergenic",
                "gvcf_variant_site_dropped",
                "dot_breakend_without_svtype_pair",
                "missing_alt_annotated",
                "non_integer_end_annotated",
                "cnv_tr_empty_expansion_dropped",
                "star_alone_svtype_annotated",
                "non_numeric_pos_numified",
                "cnv_tr_rn_zero_allele_dropped",
            ])
        );
        for key in [
            "excluded_non_integer_end_rust",
            "excluded_star_alone_svtype_perl",
            "excluded_non_numeric_pos_perl",
            "excluded_cnv_tr_rn_zero_rust",
        ] {
            assert_eq!(count(adjusted, &[key]), 0, "{key}");
        }
        assert_eq!(
            (
                count(adjusted, &["bnd_mate_32bit_buffer_records"]),
                count(adjusted, &["bnd_mate_32bit_buffer_reference_tuples"])
            ),
            (0, 0)
        );
        assert_eq!(count(adjusted, &["bnd_bare_reference_allele_rows"]), 0);
        let table = read(&out.path.join("discordant_adjusted.tsv"));
        let surviving: Vec<&str> = table.lines().skip(1).collect();
        assert_eq!(surviving.len(), 3);
        assert!(!surviving.iter().any(|line| {
            line.contains("\t21:50001-60000\t")
                || line.contains("\t21:30000001\t")
                || line.contains("\t21:70001-70100\t")
        }));
        assert!(
            surviving.iter().any(|line| line.contains("ENST22_X")),
            "a row on a transcript of another chromosome is charged under 116.2"
        );
        assert!(
            surviving.iter().any(|line| line.contains("ENST21_G3")),
            "the surplus row on a kept record is charged"
        );
        let markdown = read(&out.path.join("concordance_report.md"));
        assert!(markdown.contains(
            "- **Reference**: Ensembl VEP 116.2; masks applied: `filter_cnv_tr_expansion_divergences`"
        ));
        assert!(markdown.contains(
            "- **116.2 classes**: `reference_skipped_record` 6 vep-rs rows on 4 records set aside;"
        ));
    }

    #[test]
    fn test_116_2_refuses_a_reference_that_kept_a_record_on_a_skip_route() {
        let (out, finished) = run_release("116.2", "perl");
        assert_eq!(finished.code, Some(1));
        assert!(finished.stderr.contains(
            "ERROR: [compare_sv_concordance] skipped-record cross-check: cpx_kept N><CPX> at 21:45716816 is on the unsupported_type route"
        ));
        assert!(finished
            .stderr
            .contains("gvcf_block N><NON_REF> at 21:70001 is on the unsupported_type route"));
        assert!(
            !out.path.join("concordance_report.json").exists(),
            "nothing is scored over a failed cross-check"
        );
    }

    #[test]
    fn test_116_2_sets_aside_the_reference_allele_rows_of_a_bare_breakpoint() {
        let (fixture, cache) = with_committed_cache();
        fixture.write_vcf(
            "bare.vcf",
            &["21\t26032000\tbare_bnd\tN\t<BND>\t.\tPASS\tSVTYPE=BND"],
        );
        let own = "bare_bnd\t21:26032001\tchromosome_breakpoint\tG\tENST21_OWN\tTranscript\tintron_variant\n";
        let near = "bare_bnd\t21:26032001\tchromosome_breakpoint\tG\tENST21_G3\tTranscript\tupstream_gene_variant\n";
        fixture.write(
            "perl/bare.txt",
            &[SEVEN_COLUMNS, own, own, near, near].concat(),
        );
        fixture.write("rust/bare.txt", &[SEVEN_COLUMNS, own, near].concat());
        let report = fixture.scored("116.2", "grch37", &cache);
        assert_eq!(
            (
                count(&report, &["overall", "perl_tuples"]),
                count(&report, &["overall", "rust_tuples"])
            ),
            (2, 2)
        );
        assert_eq!(ratio(&report, &["overall", "f1"]), 1.0);
        let adjusted = at(&report, &["adjusted"]);
        assert_eq!(ratio(adjusted, &["f1"]), 1.0);
        assert_eq!(count(adjusted, &["bnd_bare_reference_allele_rows"]), 2);
        assert_eq!(
            at(adjusted, &["bnd_bare_reference_allele_rows_by_file"]),
            &json!({"bare": 2})
        );
        assert!(
            names(adjusted, &["masks_applied"]).contains(&"filter_bnd_bare_reference_allele_rows")
        );
        let class = at(
            adjusted,
            &["divergence_classes", "bnd_bare_reference_allele_row"],
        );
        assert_eq!(
            (at(class, &["excludes"]), at(class, &["excluded"])),
            (&json!(true), &json!(2))
        );
        assert_eq!(
            count(adjusted, &["excluded_perl_extra"]),
            0,
            "the rows left the tuple sets before either figure"
        );
        assert_eq!(
            (
                count(adjusted, &["perl_tuples"]),
                count(adjusted, &["intersection"])
            ),
            (2, 2)
        );
        assert!(fixture
            .read("out/concordance_report.md")
            .contains("`bnd_bare_reference_allele_row` 2 reference-allele rows set aside"));
    }

    #[test]
    fn test_116_2_exits_0_on_a_dot_breakend_whose_point_lies_inside_the_reference_span() {
        let (fixture, cache) = with_committed_cache();
        let features = [
            (
                "ENSG00000186971",
                "ENST00000334068",
                "upstream_gene_variant",
            ),
            (
                "ENSG00000240432",
                "ENST00000390690",
                "coding_sequence_variant",
            ),
        ];
        for (vid, reference, alt, span) in [
            ("d3eq_trail", "NAC", "NA.", "21:30425650-30425652"),
            ("d3eq_lead", "NAC", ".AC", "21:30425650-30425652"),
            ("d4eq_trail", "NACG", "NAC.", "21:30425650-30425653"),
        ] {
            fixture.write_vcf(
                &format!("{vid}.vcf"),
                &[&format!("21\t30425650\t{vid}\t{reference}\t{alt}\t.\t.\t.")],
            );
            for (side, location) in [("perl", span), ("rust", "21:30425651")] {
                let rows: String = features
                    .iter()
                    .map(|(gene, feature, consequence)| {
                        format!(
                            "{vid}\t{location}\tchromosome_breakpoint\t{gene}\t{feature}\tTranscript\t{consequence}\n"
                        )
                    })
                    .collect();
                fixture.write(
                    &format!("{side}/{vid}.txt"),
                    &[SEVEN_COLUMNS, &rows].concat(),
                );
            }
        }
        let report = fixture.scored("116.2", "grch38", &cache);
        let counts = |row: &Value| {
            ["perl_tuples", "rust_tuples", "intersection"].map(|key| count(row, &[key]))
        };
        assert_eq!(counts(at(&report, &["overall"])), [6, 6, 0]);
        let per_file = object(&report, &["per_file"]);
        assert_eq!(per_file.len(), 3);
        for (name, row) in per_file {
            assert_eq!(counts(row), [2, 2, 0], "{name}");
        }
        assert!(
            names(&report, &["adjusted", "masks_applied"]).contains(&"filter_dot_breakend_pairs")
        );
        assert_eq!(
            at(
                &report,
                &[
                    "adjusted",
                    "divergence_classes",
                    "dot_breakend_without_svtype_pair",
                    "excludes"
                ]
            ),
            &json!(true)
        );
        assert_eq!(
            line_feeds(&fixture.read("out/discordant.tsv")),
            13,
            "the header and the twelve one-sided rows"
        );
        assert_eq!(
            line_feeds(&fixture.read("out/discordant_adjusted.tsv")),
            1,
            "the header alone: every row paired"
        );
    }

    #[test]
    fn test_the_115_2_report_carries_no_release_block() {
        let (out, _) = run_release("115.2", "perl");
        let report = report(&out.path);
        let adjusted = object(&report, &["adjusted"]);
        assert!(!adjusted.contains_key("reference_release"));
        assert!(!adjusted.contains_key("masks_applied"));
        assert!(!read(&out.path.join("concordance_report.md")).contains("Reference"));
    }
}

mod skip_route_cross_check_tests {
    use super::*;
    use vep_sv_compare::cache::CacheAuthority;
    use vep_sv_compare::masks::skipped::{check_skip_routes, skipped_record_rows};
    use vep_sv_compare::masks::{apply_masks, MaskContext, MAX_SV_SIZE};
    use vep_sv_compare::normalize::{
        normalize_allele, normalize_consequence_set, normalize_location,
    };
    use vep_sv_compare::strings::OwnedSides;
    use vep_sv_compare::vcf::{classify_variant, VcfVariant};

    /// An input record of chromosome 21 as the reader of an input VCF builds it.
    fn record(pos: i128, reference: &str, alt: &str, info: &str, vid: &str) -> VcfVariant {
        VcfVariant {
            chrom: "21".to_string(),
            pos,
            vid: vid.to_string(),
            reference: reference.to_string(),
            alt: alt.to_string(),
            variant_type: classify_variant(reference, alt).to_string(),
            location: format!("21:{pos}"),
            info: info.to_string(),
            pos_text: None,
        }
    }

    /// A Transcript tuple, its columns normalised as an output row's are.
    fn tuple(location: &str, allele: &str, feature: &str, consequence: &str) -> [String; 5] {
        [
            normalize_location(location),
            normalize_allele(allele, "vep-rs"),
            feature.to_string(),
            "Transcript".to_string(),
            normalize_consequence_set(consequence),
        ]
    }

    fn texts(tuple: &[String; 5]) -> [&str; 5] {
        tuple.each_ref().map(String::as_str)
    }

    fn rows<'a>(tuples: &[&'a [String; 5]]) -> Vec<[&'a str; 5]> {
        tuples.iter().map(|tuple| texts(tuple)).collect()
    }

    fn sides(perl: &[&[String; 5]], rust: &[&[String; 5]]) -> OwnedSides {
        OwnedSides::from_rows(&rows(perl), &rows(rust))
    }

    /// Two kept records, a `<DEL>` both sides wrote and a `<DUP>` the reference wrote
    /// nothing for: the scored engine's row on the `<DUP>` is on no skip route. For vep-rs
    /// that stops the comparison; for another scored engine it is counted in the report and
    /// left in both figures, the run exiting 0 with the row scored as a miss. A reference
    /// that kept a record on a skip route stops every engine.
    #[test]
    fn test_a_foreign_engine_counts_the_skip_residue_instead_of_stopping() {
        let records = [
            record(1000, "N", "<DEL>", "SVTYPE=DEL;END=2000", "del"),
            record(5000, "N", "<DUP>", "SVTYPE=DUP;END=6000", "dup"),
        ];
        let deletion = tuple("21:1001-2000", "deletion", "ENST_D", "feature_truncation");
        let residue = tuple("21:5001-6000", "duplication", "ENST_U", "intron_variant");
        let owned = sides(&[&deletion], &[&deletion, &residue]);
        let residue_set = owned.set(&[texts(&residue)]);
        let both = owned.sides();
        assert_eq!(
            check_skip_routes(&both, &records, MAX_SV_SIZE, "vep-rs").len(),
            1
        );
        assert!(check_skip_routes(&both, &records, MAX_SV_SIZE, "fastvep").is_empty());
        let (excluded, confirmed, problems, counted) =
            skipped_record_rows(&both, &records, MAX_SV_SIZE, "fastvep");
        assert!(excluded.is_empty());
        assert_eq!(confirmed, 0);
        assert!(problems.is_empty());
        assert_eq!(counted, residue_set);
        let (excluded, confirmed, problems, counted) =
            skipped_record_rows(&both, &records, MAX_SV_SIZE, "vep-rs");
        assert!(excluded.is_empty());
        assert_eq!(confirmed, 0);
        assert!(counted.is_empty());
        assert_eq!(
            problems,
            ["21:5001-6000 duplication: the reference wrote nothing at this Location and no input record on a skip route starts there"]
        );
        let cpx = [record(
            45716815,
            "N",
            "<CPX>",
            "END=45730272;SVTYPE=CPX",
            "cpx",
        )];
        let kept = tuple("21:45716816-45730272", "CPX", "ENST_C", "intron_variant");
        let kept_both = sides(&[&kept], &[&kept]);
        assert_eq!(
            check_skip_routes(&kept_both.sides(), &cpx, MAX_SV_SIZE, "fastvep").len(),
            1,
            "a reference that kept a skip-route record stops every engine"
        );
        let cache = CacheAuthority::default();
        let masks = apply_masks(
            &both,
            &MaskContext {
                reference_release: "116.2",
                cache: &cache,
                scored_engine: "fastvep",
                input_variants: Some(&records),
                max_sv_size: MAX_SV_SIZE,
                perl_paths: &[],
                rust_paths: &[],
                input_by_output: &BTreeMap::new(),
                coding_ends: None,
            },
        )
        .expect("the masks of 116.2 run");
        assert!(masks.skipped_problems.is_empty());
        assert_eq!(masks.scored_engine_skip_residue, residue_set);
        assert!(masks.all_rust().is_empty());

        let fixture = Fixture::new();
        fixture.write_vcf(
            "two.vcf",
            &[
                "21\t1000\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000",
                "21\t5000\tdup\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=6000",
            ],
        );
        let deletion = "del\t21:1001-2000\tdeletion\tG\tENST_D\tTranscript\tfeature_truncation\n";
        let duplication = "dup\t21:5001-6000\tduplication\tG\tENST_U\tTranscript\tintron_variant\n";
        fixture.write("perl/two.txt", &[SEVEN_COLUMNS, deletion].concat());
        fixture.write(
            "rust/two.txt",
            &[SEVEN_COLUMNS, deletion, duplication].concat(),
        );
        let stopped = fixture.compare("116.2", "grch38", &[]);
        assert_eq!(stopped.code, Some(1));
        assert!(!fixture.path("out/concordance_report.json").exists());
        let report = fixture.scored("116.2", "grch38", &strings(&["--scored-engine", "fastvep"]));
        assert_eq!(
            ["perl_tuples", "rust_tuples", "intersection"]
                .map(|key| count(&report, &["overall", key])),
            [1, 2, 1]
        );
        let adjusted = at(&report, &["adjusted"]);
        assert_eq!(text(adjusted, &["scored_engine"]), "fastvep");
        assert_eq!(
            count(adjusted, &["scored_engine_skip_route_residue_rows"]),
            1
        );
        assert_eq!(
            at(adjusted, &["scored_engine_skip_route_residue_by_file"]),
            &json!({"two": 1})
        );
        assert_eq!(
            ["rust_tuples", "only_rust", "excluded_rust_extra"].map(|key| count(adjusted, &[key])),
            [2, 1, 0],
            "counted, not set aside"
        );
        assert_eq!(
            (ratio(&report, &["overall", "f1"]), ratio(adjusted, &["f1"])),
            (0.666667, 0.666667)
        );
        assert!(fixture.read("out/concordance_report.md").contains(
            "`scored_engine_skip_route_residue` 1 scored-engine rows at Locations the reference wrote nothing for counted, none set aside"
        ));
        assert_eq!(
            line_feeds(&fixture.read("out/discordant_adjusted.tsv")),
            2,
            "the header and the counted row"
        );
    }
}

mod input_record_hand_over_tests {
    use super::*;

    const BUCKET: &str = "splice_term_past_the_record_span";
    const LOCATION: &str = "7:100999-101000";

    /// The column line of an output of fourteen columns.
    const ROW_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

    /// Two bases replaced by five, and five bases inserted after the first of the two: both
    /// records give their rows the same Location and the same allele.
    const REPLACED: &str = "7\t100999\t1993671\tAC\tGTTCA\t.\tPASS\t.\n";
    const INSERTED: &str = "7\t100999\t1993671\tA\tAGTTCA\t.\tPASS\t.\n";

    /// One row of fourteen columns under an identifier, at [`LOCATION`].
    fn output_row(consequence: &str, cdna: &str, cds: &str, feature: &str) -> String {
        format!(
            "1993671\t{LOCATION}\tGTTCA\tENSG1\t{feature}\tTranscript\t{consequence}\t{cdna}\t{cds}\t-\t-\t-\t-\tIMPACT=LOW;STRAND=-1\n"
        )
    }

    /// The report and the lines of `discordant_adjusted.tsv` of a directory of one input
    /// and two outputs: the two rows of one pair on one transcript, the reference's and
    /// the scored engine's, beside a row both engines write.
    fn scored(record: &str, input_name: &str) -> (Value, Vec<String>) {
        let fixture = Fixture::new();
        fixture.write(&format!("in/{input_name}"), &[VCF_HEAD, record].concat());
        let matched = output_row("intron_variant", "-", "-", "ENST00000000009");
        for (side, terms) in [
            (
                "perl",
                "splice_donor_variant,frameshift_variant,splice_donor_region_variant,intron_variant",
            ),
            ("rust", "frameshift_variant,splice_region_variant"),
        ] {
            let row = output_row(terms, "999-1000", "899-900", "ENST00000635895");
            fixture.write(
                &format!("{side}/one.txt"),
                &[ROW_COLUMNS, &matched, &row].concat(),
            );
        }
        let report = fixture.scored("116.2", "grch38", &[]);
        (report, fixture.lines("discordant_adjusted.tsv"))
    }

    fn assert_set_aside(report: &Value, adjusted_lines: &[String]) {
        let adjusted = at(report, &["adjusted"]);
        assert_eq!(
            count(adjusted, &["excluded_registry_swap_by_bucket", BUCKET]),
            1
        );
        assert_eq!(adjusted_lines.len(), 1, "the header alone");
        assert_eq!(
            (
                count(adjusted, &["only_perl"]),
                count(adjusted, &["only_rust"]),
                ratio(adjusted, &["f1"])
            ),
            (0, 0, 1.0)
        );
        assert_eq!(ratio(report, &["per_file_adjusted", "one", "f1"]), 1.0);
        assert_eq!(
            (
                count(report, &["overall", "only_perl"]),
                count(report, &["overall", "only_rust"])
            ),
            (1, 1)
        );
    }

    #[test]
    fn test_the_entry_point_reads_the_record_from_the_input_folder() {
        let (report, adjusted_lines) = scored(REPLACED, "one.vcf");
        assert_set_aside(&report, &adjusted_lines);
    }

    #[test]
    fn test_a_gzipped_input_reads_as_a_plain_one() {
        let (report, adjusted_lines) = scored(REPLACED, "one.vcf.gz");
        assert_set_aside(&report, &adjusted_lines);
    }

    #[test]
    fn test_an_insertion_s_pair_stays_a_difference_with_its_record_handed() {
        let (report, adjusted_lines) = scored(INSERTED, "one.vcf");
        let adjusted = at(&report, &["adjusted"]);
        assert_eq!(
            count(adjusted, &["excluded_registry_swap_by_bucket", BUCKET]),
            0
        );
        assert_eq!(
            (
                count(adjusted, &["only_perl"]),
                count(adjusted, &["only_rust"]),
                adjusted_lines.len()
            ),
            (1, 1, 3)
        );
    }
}

mod goldens {
    use super::*;
    use vep_sv_compare::cli::usage;
    use vep_sv_compare::threads::THREADS_ENV;

    /// The comparison of a committed fixture's `in`, `rust` and `cache` directories and the
    /// reference output under `perl`, with `options` after them, writing its reports to
    /// `out`.
    fn comparison(
        fixture: &Path,
        release: &str,
        perl: &str,
        options: &[&str],
        out: &Path,
    ) -> Command {
        let mut command = tool();
        command.args(["--reference-release", release]);
        of_fixture(command, fixture, perl, options, out)
    }

    /// `command` with the directories of a committed fixture, the reference output under
    /// `perl` and `options` after them, writing its reports to `out`.
    fn of_fixture(
        mut command: Command,
        fixture: &Path,
        perl: &str,
        options: &[&str],
        out: &Path,
    ) -> Command {
        command.args(["--assembly", "grch37"]);
        command.arg("--input-dir").arg(fixture.join("in"));
        command.arg("--perl-dir").arg(fixture.join(perl));
        command.arg("--rust-dir").arg(fixture.join("rust"));
        command.arg("--vep-rs-cache").arg(fixture.join("cache"));
        command.arg("--output-dir").arg(out);
        command.args(options);
        command
    }

    /// Holds one comparison of a committed fixture to the goldens under its `expected`
    /// directory, on one thread, on two and on the host's default: the exit code, stdout,
    /// stderr and the four report files byte for byte, a report file the golden lacks being
    /// one the run does not write. The comparison reads the fixture's `in`, `rust` and
    /// `cache` directories and the reference output under `perl`, with `options` after
    /// them. A golden names the output directory `{{OUTPUT_DIR}}`.
    fn assert_goldens(fixture: &str, expected: &str, release: &str, perl: &str, options: &[&str]) {
        let fixture = fixtures().join(fixture);
        let expected = fixture.join(expected);
        let golden = |name: &str| std::fs::read_to_string(expected.join(name)).ok();
        for threads in [Some("1"), Some("2"), None] {
            let out = TempDir::new("vep_sv_compare_cli_").expect("a scratch directory");
            let mut command = comparison(&fixture, release, perl, options, &out.path);
            match threads {
                Some(count) => command.env(THREADS_ENV, count),
                None => command.env_remove(THREADS_ENV),
            };
            let finished = finish(command);
            let placed = |text: String| {
                text.replace(
                    "{{OUTPUT_DIR}}",
                    out.path.to_str().expect("a UTF-8 scratch path"),
                )
            };
            let run = format!("{}, threads {threads:?}", expected.display());
            assert_eq!(
                finished.code.map(|code| format!("{code}\n")),
                golden("exit_code.txt"),
                "exit code of {run}"
            );
            assert_eq!(
                Some(finished.stdout),
                golden("stdout.txt").map(placed),
                "stdout of {run}"
            );
            assert_eq!(
                finished.stderr,
                golden("stderr.txt").map(placed).unwrap_or_default(),
                "stderr of {run}"
            );
            for name in REPORT_FILES {
                assert_eq!(
                    std::fs::read_to_string(out.path.join(name)).ok(),
                    golden(name),
                    "{name} of {run}"
                );
            }
        }
    }

    // `reference_release`: one input whose reference output is committed in the 115.2
    // shape (`perl`) and in the 116.2 shape (`perl_116`).

    #[test]
    fn reference_release_under_115_2() {
        assert_goldens("reference_release", "expected_115_2", "115.2", "perl", &[]);
    }

    #[test]
    fn reference_release_under_116_2() {
        assert_goldens(
            "reference_release",
            "expected_116_2",
            "116.2",
            "perl_116",
            &[],
        );
    }

    #[test]
    fn reference_release_in_the_115_2_shape_is_refused_under_116_2() {
        assert_goldens(
            "reference_release",
            "expected_116_2_refused",
            "116.2",
            "perl",
            &[],
        );
    }

    // `release_classes`: four inputs with a record for every class of the 116.2 reference
    // that `reference_release` does not reach, a pair under each pair mask, a tuple two
    // files hold and rows no mask takes.

    #[test]
    fn release_classes_under_115_2() {
        assert_goldens("release_classes", "expected_115_2", "115.2", "perl", &[]);
    }

    #[test]
    fn release_classes_under_116_2() {
        assert_goldens("release_classes", "expected_116_2", "116.2", "perl", &[]);
    }

    /// The start pair of `02_one_sided` lies on `ENST21_B` at `21:38500001-38500672`. The
    /// cache's shard `0-0.json` of sequence 21 holds the record of that transcript whose span
    /// holds the Location and whose coding sequence opens on a start codon; it is filed
    /// before the shard `0-1.json`, whose record of the transcript keeps the span the other
    /// masks read. Beside a cache without that shard no record of the transcript holds the
    /// Location: the pair's two rows stay differences and no other row moves.
    #[test]
    fn release_classes_leave_the_start_pair_a_difference_beside_a_cache_without_its_record() {
        let fixture = fixtures().join("release_classes");
        let scratch = TempDir::new("vep_sv_compare_cli_").expect("a scratch directory");
        let without = scratch.path.join("cache");
        for shard in ["21/0-1.json", "22/0-1.json"] {
            let to = without.join("transcripts").join(shard);
            std::fs::create_dir_all(to.parent().expect("a sequence directory")).unwrap();
            std::fs::copy(fixture.join("cache/transcripts").join(shard), to).unwrap();
        }
        // The rows of `discordant_adjusted.tsv` at the pair's Location, and the others.
        let surviving = |cache: &Path| -> (Vec<String>, Vec<String>) {
            let out = TempDir::new("vep_sv_compare_cli_").expect("a scratch directory");
            let mut command = tool();
            command.args(["--reference-release", "116.2", "--assembly", "grch37"]);
            command.arg("--input-dir").arg(fixture.join("in"));
            command.arg("--perl-dir").arg(fixture.join("perl"));
            command.arg("--rust-dir").arg(fixture.join("rust"));
            command.arg("--vep-rs-cache").arg(cache);
            command.arg("--output-dir").arg(&out.path);
            command.env(THREADS_ENV, "1");
            let finished = finish(command);
            assert_eq!(finished.code, Some(0), "stderr: {}", finished.stderr);
            read(&out.path.join("discordant_adjusted.tsv"))
                .lines()
                .map(str::to_string)
                .partition(|line| line.contains("\t21:38500001-38500672\t"))
        };
        let (set_aside, others) = surviving(&fixture.join("cache"));
        assert_eq!(set_aside, Vec::<String>::new());
        let (left, others_without) = surviving(&without);
        let columns = |line: &String| -> Vec<String> {
            let cells: Vec<&str> = line.split('\t').collect();
            [cells[0], cells[2], cells[5]].map(str::to_string).to_vec()
        };
        assert_eq!(
            left.iter().map(columns).collect::<Vec<_>>(),
            [
                ["02_one_sided", "only_perl", "ENST21_B"],
                ["02_one_sided", "only_rust", "ENST21_B"]
            ]
        );
        assert_eq!(others_without, others);
    }

    #[test]
    fn release_classes_without_a_size_cap_are_refused_under_116_2() {
        assert_goldens(
            "release_classes",
            "expected_116_2_no_size_cap",
            "116.2",
            "perl",
            &["--reference-max-sv-size", "-1"],
        );
    }

    #[test]
    fn release_classes_scored_as_another_engine_under_116_2() {
        assert_goldens(
            "release_classes",
            "expected_116_2_other_engine",
            "116.2",
            "perl",
            &[
                "--scored-engine",
                "fastvep",
                "--reference-max-sv-size",
                "-1",
            ],
        );
    }

    /// Ensembl VEP 116.3 is scored under the rules of 116.2: on each comparison the 116.2
    /// goldens hold, a run naming 116.3 ends as the run naming 116.2 does and writes its
    /// bytes, but for the release a scored run states on stdout and in two reports. A
    /// release outside the list is a usage error that names the three listed.
    #[test]
    fn a_run_naming_116_3_writes_what_a_run_naming_116_2_writes_and_116_1_is_refused() {
        // Each place an output states the release: stdout is output 0, stderr 1, the reports 2 to 5.
        let statements = [
            (0, "Reference release RELEASE: "),
            (2, "\"reference_release\": \"RELEASE\""),
            (3, "- **Reference**: Ensembl VEP RELEASE; "),
            (3, "- **RELEASE classes**: "),
        ];
        let cases: [(&str, &str, &[&str], i32); 5] = [
            ("reference_release", "perl_116", &[], 0),
            ("reference_release", "perl", &[], 1),
            ("release_classes", "perl", &[], 0),
            (
                "release_classes",
                "perl",
                &["--reference-max-sv-size", "-1"],
                1,
            ),
            (
                "release_classes",
                "perl",
                &[
                    "--scored-engine",
                    "fastvep",
                    "--reference-max-sv-size",
                    "-1",
                ],
                0,
            ),
        ];
        for (fixture, perl, options, code) in cases {
            let case = format!("{fixture} against {perl} with {options:?}");
            let out = TempDir::new("vep_sv_compare_cli_").expect("a scratch directory");
            let outputs_under = |release: &str| -> Vec<Option<String>> {
                let mut command =
                    comparison(&fixtures().join(fixture), release, perl, options, &out.path);
                command.env_remove(THREADS_ENV);
                let finished = finish(command);
                assert_eq!(finished.code, Some(code), "{case} under {release}");
                let mut outputs = vec![Some(finished.stdout), Some(finished.stderr)];
                for name in REPORT_FILES {
                    let path = out.path.join(name);
                    outputs.push(std::fs::read_to_string(&path).ok());
                    if path.exists() {
                        std::fs::remove_file(&path).expect("a report file is removed");
                    }
                }
                outputs
            };
            let mut expected = outputs_under("116.2");
            let written = outputs_under("116.3");
            assert_eq!(
                expected[2..].iter().flatten().count(),
                if code == 0 { REPORT_FILES.len() } else { 0 },
                "the reports of {case}"
            );
            for (index, statement) in statements {
                let under_116_2 = statement.replace("RELEASE", "116.2");
                let stated = expected[index]
                    .as_deref()
                    .map_or(0, |output| output.matches(&under_116_2).count());
                assert_eq!(
                    stated,
                    usize::from(code == 0),
                    "`{under_116_2}` in output {index} of {case}"
                );
                expected[index] = expected[index].take().map(|output| {
                    output.replace(&under_116_2, &statement.replace("RELEASE", "116.3"))
                });
            }
            for (index, (expected, written)) in expected.iter().zip(&written).enumerate() {
                assert!(
                    expected == written,
                    "output {index} of {case} under 116.3 differs from 116.2's in more than the release it states"
                );
            }
        }

        let refused = run(&strings(&["--reference-release", "116.1"]));
        assert_eq!(refused.code, Some(2));
        assert_eq!(refused.stdout, "");
        assert!(usage().contains("[--reference-release {115.2,116.2,116.3}]\n"));
        assert_eq!(
            refused.stderr,
            format!(
                "{}vep-sv-compare: error: argument --reference-release: invalid choice: '116.1' \
                 (choose from '115.2', '116.2', '116.3')\n",
                usage()
            )
        );
    }

    /// A run that names no release is scored against Ensembl VEP 116.3: on a comparison
    /// that writes its reports and on one that is refused, it ends as the run naming 116.3
    /// does and writes its bytes, stdout, stderr and the four reports, and a scored run
    /// states that release.
    #[test]
    fn a_run_naming_no_release_writes_what_a_run_naming_116_3_writes() {
        let cases: [(&str, &str, i32); 3] = [
            ("reference_release", "perl_116", 0),
            ("reference_release", "perl", 1),
            ("release_classes", "perl", 0),
        ];
        for (fixture, perl, code) in cases {
            let case = format!("{fixture} against {perl}");
            let out = TempDir::new("vep_sv_compare_cli_").expect("a scratch directory");
            let outputs_of = |mut command: Command, run: &str| -> Vec<Option<String>> {
                command.env_remove(THREADS_ENV);
                let finished = finish(command);
                assert_eq!(finished.code, Some(code), "{case}, {run}");
                let mut outputs = vec![Some(finished.stdout), Some(finished.stderr)];
                for name in REPORT_FILES {
                    let path = out.path.join(name);
                    outputs.push(std::fs::read_to_string(&path).ok());
                    if path.exists() {
                        std::fs::remove_file(&path).expect("a report file is removed");
                    }
                }
                outputs
            };
            let fixture = fixtures().join(fixture);
            let named = outputs_of(
                comparison(&fixture, "116.3", perl, &[], &out.path),
                "naming 116.3",
            );
            let unnamed = outputs_of(
                of_fixture(tool(), &fixture, perl, &[], &out.path),
                "naming no release",
            );
            let stated = named[0].as_deref().map_or(0, |stdout| {
                stdout.matches("Reference release 116.3: ").count()
            });
            assert_eq!(
                stated,
                usize::from(code == 0),
                "the release stdout states, {case}"
            );
            for (index, (named, unnamed)) in named.iter().zip(&unnamed).enumerate() {
                assert!(
                    named == unnamed,
                    "output {index} of {case} naming no release differs from the run naming 116.3"
                );
            }
        }
    }
}
