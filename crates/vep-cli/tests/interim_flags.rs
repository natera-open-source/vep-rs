// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The flags of Ensembl VEP features vep-rs does not implement, through the `vep`
//! binary: `--custom` and the statistics flags are accepted and change nothing,
//! `--regulatory` is accepted, warns and declares its columns empty, `--refseq`
//! and `--merged` select the suffixed species directory of the cache root, and
//! `--stats_html`, `--all_refseq` and `--bam` are refused as unknown arguments.
//! Each test quotes the documentation sentence it proves and cites the Ensembl
//! VEP source (release/116.2) that sentence describes.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn corpus_cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

/// The committed `--custom` source: a bgzipped BED with its tabix index, one
/// region around each record of [`RECORDS`].
fn custom_fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/custom/regions.bed.gz")
}

/// Records of the GRCh38 release 116 corpus on both of its chromosomes, each
/// inside one region of the custom fixture.
const RECORDS: &str = "\
21\t30282169\tstop_ins4_twostops_rev\tC\tCTATT\t.\t.\t.
21\t30348409\tstop_del2_lost_rev\tCTA\tC\t.\t.\t.
21\t30430756\tstop_ins5_newstop_fwd\tT\tTCCTAA\t.\t.\t.
21\t34716478\tpenult_ins3_stop_fwd\tA\tAATA\t.\t.\t.
22\t20061573\tpartial_codon_del3_fwd\tGAAG\tG\t.\t.\t.
";

/// The columns `--regulatory` declares, the `regulatory` entry of `@FLAG_FIELDS`
/// (`Constants.pm:136`): `BIOTYPE`, filled only under `--biotype`
/// (`OutputFactory.pm:1580`), then the motif columns.
const REGULATORY_COLUMNS: [&str; 6] = [
    "BIOTYPE",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
];

/// The columns that stay empty under `--regulatory` even when `--biotype`, which
/// `--everything` and `--vcf` switch on, fills `BIOTYPE`.
const MOTIF_COLUMNS: [&str; 5] = {
    let [_biotype, motif_columns @ ..] = REGULATORY_COLUMNS;
    motif_columns
};

const FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

fn write_input(dir: &Path) -> PathBuf {
    let path = dir.join("in.vcf");
    fs::write(
        &path,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{RECORDS}"),
    )
    .unwrap();
    path
}

/// Runs the binary on `in.vcf` under `dir`, writing `<dir>/<label>.<format>`,
/// with `cache` naming the cache and `extra` the flags under test. Nothing is
/// asserted about the exit status.
fn run(dir: &Path, label: &str, format: &str, cache: &[&str], extra: &[&str]) -> Output {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(dir.join("in.vcf"))
        .arg("-o")
        .arg(dir.join(format!("{label}.{format}")))
        .args(cache)
        .args(["--assembly", "GRCh38", "--force_overwrite"])
        .args(extra)
        .stdin(Stdio::null());
    if format != "default" {
        cmd.arg(format!("--{format}"));
    }
    cmd.output().unwrap()
}

/// [`run`] on the corpus cache, asserting success and at least one data line
/// per record of [`RECORDS`]; the output text and the standard error.
fn annotate(dir: &Path, label: &str, format: &str, extra: &[&str]) -> (String, String) {
    let cache = corpus_cache();
    let output = run(
        dir,
        label,
        format,
        &["--json_cache", cache.to_str().unwrap()],
        extra,
    );
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{extra:?} ({format}): {stderr}");
    let text = fs::read_to_string(dir.join(format!("{label}.{format}"))).unwrap();
    assert!(
        data_lines(&text).len() >= RECORDS.lines().count(),
        "{extra:?} ({format}): the run annotates"
    );
    (text, stderr)
}

/// The output without the header lines that differ between any two runs: the
/// command line, which names the flags under test, the time of the run and the
/// cache path.
fn without_volatile_headers(text: &str) -> Vec<&str> {
    const VOLATILE: [&str; 5] = [
        "## VEP command-line:",
        "## Output produced at ",
        "## Using cache in ",
        "##VEP-command-line=",
        "##VEP=",
    ];
    text.lines()
        .filter(|line| !VOLATILE.iter().any(|prefix| line.starts_with(prefix)))
        .collect()
}

fn data_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|line| !line.starts_with('#')).collect()
}

/// Asserts the run with `flags` writes what the run without them writes, in
/// every output format, apart from the volatile headers.
fn assert_identical_output(dir: &Path, label: &str, flags: &[&str]) {
    for format in FORMATS {
        let (without, _) = annotate(dir, "plain", format, &[]);
        let (with, stderr) = annotate(dir, label, format, flags);
        assert_eq!(
            without_volatile_headers(&with),
            without_volatile_headers(&without),
            "{flags:?} changes the {format} output"
        );
        assert!(
            !stderr.contains("WARNING") && !stderr.contains("ERROR"),
            "{flags:?} ({format}): {stderr}"
        );
    }
}

/// `--custom` is accepted and the source is never opened: the output of every
/// format is the output without the flag.
///
/// Documented: "`--custom <file>,<name>,<format>,...` is accepted and the source
/// is not read, so the run proceeds and no custom field appears" (What vep-rs
/// does not do, Custom annotation sources) and "Its `##INFO=<ID=<short>_<field>,...>`
/// lines for `--custom` sources are not written, since custom annotation sources
/// are not supported" (VCF format). Ensembl VEP reads the string as
/// `file,short_name,format,type,coords,fields...` or as `key=value` pairs
/// (`AnnotationSourceAdaptor.pm:353-364`), builds one `File` source per string
/// (`:455`; `File.pm:143-197`), opens the file through the tabix parser
/// (`File/BED.pm:82`; ensembl-io `TabixParser.pm:52-63`, which dies unless the
/// file is bgzipped and indexed) and writes the source's values under its short
/// name on every overlapping record (`OutputFactory.pm:1004-1012`), with one
/// `##INFO` line per source in VCF (`OutputFactory/VCF.pm:276`).
#[test]
fn custom_is_accepted_and_changes_no_output() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let fixture = custom_fixture();
    let bytes = fs::read(&fixture).unwrap();
    assert!(
        bytes.starts_with(&[0x1f, 0x8b, 0x08, 0x04]),
        "the fixture is a BGZF file"
    );
    assert!(
        PathBuf::from(format!("{}.tbi", fixture.display())).is_file(),
        "the fixture carries its tabix index"
    );
    let source = fixture.to_str().unwrap();
    let positional = format!("{source},regions,bed,overlap");
    let keyed = format!("file={source},short_name=regions,format=bed,type=overlap");
    let absent = dir.path().join("absent.bed.gz");
    let absent = format!("{},regions,bed,overlap", absent.display());
    for (label, value) in [
        ("positional", positional.as_str()),
        ("keyed", keyed.as_str()),
        ("absent", absent.as_str()),
    ] {
        assert_identical_output(dir.path(), label, &["--custom", value]);
        let vcf = fs::read_to_string(dir.path().join(format!("{label}.vcf"))).unwrap();
        assert!(
            !vcf.lines()
                .any(|line| line.starts_with("##INFO=<ID=regions")),
            "{label}: no ##INFO line for the source"
        );
        assert!(
            !without_volatile_headers(&vcf)
                .iter()
                .any(|line| line.contains("region_")),
            "{label}: no value of the source reaches the output"
        );
    }
}

/// A cache without regulatory data warns and writes empty motif columns.
#[test]
fn missing_regulatory_data_warns_and_writes_empty_motif_columns() {
    const NOTICE: &str = "--regulatory: this cache holds no regulatory data; no RegulatoryFeature or MotifFeature row is written";
    const EVERYTHING_NOTICE: &str = NOTICE;
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());

    let (text, stderr) = annotate(dir.path(), "everything", "tab", &["--everything"]);
    assert_eq!(
        stderr.matches(EVERYTHING_NOTICE).count(),
        1,
        "everything: one warning on standard error:\n{stderr}"
    );
    let warnings = fs::read_to_string(dir.path().join("everything.tab_warnings.txt")).unwrap();
    assert!(
        warnings.contains(EVERYTHING_NOTICE),
        "everything: {warnings}"
    );
    let header: Vec<&str> = text
        .lines()
        .find(|line| line.starts_with('#') && !line.starts_with("##"))
        .expect("a tab header line")
        .split('\t')
        .collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|cell| *cell == name)
            .unwrap_or_else(|| panic!("everything: the tab header declares {name}"))
    };
    let feature_type = column("Feature_type");
    let motif_indexes = MOTIF_COLUMNS.map(column);
    for row in data_lines(&text) {
        let cells: Vec<&str> = row.split('\t').collect();
        assert!(
            !["RegulatoryFeature", "MotifFeature"].contains(&cells[feature_type]),
            "everything: no regulatory row: {row}"
        );
        for &index in &motif_indexes {
            assert_eq!(
                cells[index], "-",
                "everything: {} is empty on {row}",
                header[index]
            );
        }
    }

    let descriptions: Vec<String> = REGULATORY_COLUMNS
        .iter()
        .map(|column| format!("## {column} : "))
        .collect();
    for format in FORMATS {
        let (without, _) = annotate(dir.path(), "plain", format, &[]);
        let (with, stderr) = annotate(dir.path(), "regulatory", format, &["--regulatory"]);
        assert_eq!(
            stderr.matches(NOTICE).count(),
            1,
            "{format}: one warning on standard error:\n{stderr}"
        );
        let warnings = dir.path().join(format!("regulatory.{format}_warnings.txt"));
        let warnings = fs::read_to_string(&warnings).unwrap();
        assert_eq!(
            warnings.lines().count(),
            1,
            "{format}: one warnings-file line:\n{warnings}"
        );
        assert!(warnings.contains(NOTICE), "{format}: {warnings}");
        assert!(
            !data_lines(&with)
                .iter()
                .any(|line| line.contains("RegulatoryFeature") || line.contains("MotifFeature")),
            "{format}: no regulatory row"
        );
        let without = without_volatile_headers(&without);
        let with = without_volatile_headers(&with);
        match format {
            "json" => assert_eq!(with, without, "{format}: the records are unchanged"),
            "default" | "tab" => {
                for description in &descriptions {
                    assert_eq!(
                        with.iter()
                            .filter(|line| line.starts_with(description.as_str()))
                            .count(),
                        1,
                        "{format}: {description} declared once"
                    );
                }
                let other_lines: Vec<&str> = with
                    .iter()
                    .copied()
                    .filter(|line| !descriptions.iter().any(|d| line.starts_with(d.as_str())))
                    .collect();
                if format == "default" {
                    // The Extra column omits an empty field, so the descriptions
                    // are the whole difference.
                    assert_eq!(other_lines, without, "{format}: nothing else changes");
                } else {
                    let column_headers = REGULATORY_COLUMNS
                        .iter()
                        .map(|column| format!("\t{column}"))
                        .collect::<String>();
                    let empty_cells = "\t-".repeat(REGULATORY_COLUMNS.len());
                    let stripped: Vec<&str> = other_lines
                        .iter()
                        .map(|line| {
                            if line.starts_with("#Uploaded_variation") {
                                line.strip_suffix(column_headers.as_str())
                                    .unwrap_or_else(|| {
                                        panic!("{format}: six column names on {line}")
                                    })
                            } else if line.starts_with('#') {
                                line
                            } else {
                                line.strip_suffix(empty_cells.as_str()).unwrap_or_else(|| {
                                    panic!("{format}: six empty cells on {line}")
                                })
                            }
                        })
                        .collect();
                    assert_eq!(stripped, without, "{format}: nothing else changes");
                }
            }
            _ => {
                let fields = MOTIF_COLUMNS
                    .iter()
                    .map(|column| format!("|{column}"))
                    .collect::<String>();
                let stripped: Vec<String> = with
                    .iter()
                    .map(|line| {
                        if line.starts_with("##INFO=<ID=CSQ,") {
                            let head = line
                                .strip_suffix(&format!("{fields}\">"))
                                .unwrap_or_else(|| panic!("vcf: five field names on {line}"));
                            format!("{head}\">")
                        } else if line.starts_with('#') {
                            line.to_string()
                        } else {
                            let (record, csq) = line.split_once("CSQ=").unwrap();
                            let entries: Vec<&str> = csq
                                .split(',')
                                .map(|entry| {
                                    entry.strip_suffix("|||||").unwrap_or_else(|| {
                                        panic!("vcf: five empty fields on {entry}")
                                    })
                                })
                                .collect();
                            format!("{record}CSQ={}", entries.join(","))
                        }
                    })
                    .collect();
                assert_eq!(stripped, without, "{format}: nothing else changes");
            }
        }
    }
}

/// `--no_stats`, `--stats_file`, `--sf` and `--stats_text` are accepted and no
/// statistics file is written, with or without them.
///
/// Documented: "`--no_stats`, `--stats_file` (`--sf`) and `--stats_text` are
/// accepted and change nothing: no statistics file is written, with or without
/// them" (What vep-rs does not do, The statistics report). Ensembl VEP writes the
/// report as HTML unless `--stats_text` or `--stats_html` chooses
/// (`Runner.pm:787-789`), beside the output as `<output>_summary.html` or `.txt`,
/// or at `--stats_file` with the extension the type gives
/// (`BaseRunner.pm:254-271`), and skips it under `--no_stats` (`:300-315`).
#[test]
fn statistics_flags_are_accepted_and_write_no_summary() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let named = dir.path().join("stats.html");
    let named = named.to_str().unwrap();
    let bare = dir.path().join("report");
    let bare = bare.to_str().unwrap();
    let (reference, _) = annotate(dir.path(), "plain", "default", &[]);
    let runs: [(&str, &[&str]); 4] = [
        ("no_stats", &["--no_stats"]),
        ("stats_text", &["--stats_text"]),
        ("stats_file", &["--stats_file", named]),
        ("sf", &["--sf", bare]),
    ];
    for (label, flags) in &runs {
        let (text, stderr) = annotate(dir.path(), label, "default", flags);
        assert_eq!(
            without_volatile_headers(&text),
            without_volatile_headers(&reference),
            "{flags:?} changes the output"
        );
        assert!(!stderr.contains("WARNING"), "{flags:?}: {stderr}");
    }
    let written: BTreeSet<String> = fs::read_dir(dir.path())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    let expected: BTreeSet<String> = ["home", "in.vcf", "plain.default"]
        .into_iter()
        .map(str::to_string)
        .chain(runs.iter().map(|(label, _)| format!("{label}.default")))
        .collect();
    assert_eq!(written, expected, "no statistics file beside the outputs");
}

/// A spelling of a flag this version lacks is refused as a syntax error naming
/// the full flag, exit 2.
///
/// Documented: "An Ensembl VEP option this version lacks is refused by its full
/// name in every spelling (`-stats_html` and `--stats_h` are refused as
/// `--stats_html`, exit 2)" (Compatibility with Ensembl VEP, Option spellings).
fn assert_unknown_argument(spellings: &[&str], flag: &str) {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let cache = corpus_cache();
    for spelling in spellings {
        let output = run(
            dir.path(),
            "refused",
            "default",
            &["--json_cache", cache.to_str().unwrap()],
            &[spelling],
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(2), "{spelling}: {stderr}");
        assert_eq!(
            stderr.lines().next().unwrap_or_default(),
            format!("error: unexpected argument '{flag}' found"),
            "{spelling}: {stderr}"
        );
        assert!(
            !dir.path().join("refused.default").exists(),
            "{spelling}: no output is written"
        );
    }
}

/// `--stats_html` is refused as an unknown argument.
///
/// Documented: "`--stats_html` is refused as an unknown argument
/// (`error: unexpected argument '--stats_html' found`, exit 2)" (What vep-rs does
/// not do, The statistics report). Ensembl VEP declares it (`Config.pm:171`) and
/// writes the HTML report under it (`BaseRunner.pm:311-315`).
#[test]
fn stats_html_is_refused_as_an_unknown_argument() {
    assert_unknown_argument(
        &["--stats_html", "-stats_html", "--stats_h"],
        "--stats_html",
    );
}

/// `--all_refseq` is refused as an unknown argument.
///
/// Documented: "`--all_refseq` and `--bam` are refused as unknown arguments (exit
/// 2)" (What vep-rs does not do, RefSeq and merged caches; the `--bam` half is
/// [`bam_is_refused_as_an_unknown_argument`]). Ensembl VEP declares it
/// (`Config.pm:98`) and under it keeps the transcripts of a RefSeq or merged cache
/// whose stable id lacks the `^[A-Z]{2}_\d+` accession shape
/// (`AnnotationType/Transcript.pm:236-247`).
#[test]
fn all_refseq_is_refused_as_an_unknown_argument() {
    assert_unknown_argument(&["--all_refseq", "-all_refseq"], "--all_refseq");
}

/// `--bam` is refused as an unknown argument (the sentence of
/// [`all_refseq_is_refused_as_an_unknown_argument`]). Ensembl VEP declares it
/// (`Config.pm:253`) and edits the cache's transcript models from the aligned
/// reads of the BAM file it names (`AnnotationType/Transcript.pm:391-399, :421`).
#[test]
fn bam_is_refused_as_an_unknown_argument() {
    assert_unknown_argument(&["--bam", "-bam"], "--bam");
}

/// A cache root holding the corpus cache as `<root>/<species_dir>/116_GRCh38`
/// alone.
fn cache_root(dir: &Path, species_dir: &str) -> PathBuf {
    let root = dir.join("cache_root");
    copy_dir(&corpus_cache(), &root.join(species_dir).join("116_GRCh38"));
    root
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// Asserts `flag` resolves the cache under `<root>/<species_dir>/116_GRCh38`,
/// that the run annotates as the corpus cache does, and that the plain species
/// directory, which the root lacks, is the cache error without the flag.
fn assert_species_directory(flag: &str, species_dir: &str) {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let root = cache_root(dir.path(), species_dir);
    let layout = [
        "--cache",
        "--dir_cache",
        root.to_str().unwrap(),
        "--cache_version",
        "116",
        "--species",
        "homo_sapiens",
    ];
    let output = run(dir.path(), "suffixed", "default", &layout, &[flag]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{flag}: {stderr}");
    let text = fs::read_to_string(dir.path().join("suffixed.default")).unwrap();
    let leaf = root.join(species_dir).join("116_GRCh38");
    assert!(
        text.lines()
            .any(|line| line == format!("## Using cache in {}", leaf.display())),
        "{flag}: the header names the suffixed directory:\n{text}"
    );
    let (reference, _) = annotate(dir.path(), "plain", "default", &[]);
    assert_eq!(
        without_volatile_headers(&text),
        without_volatile_headers(&reference),
        "{flag}: the suffixed directory annotates as the corpus cache"
    );

    let output = run(dir.path(), "unsuffixed", "default", &layout, &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "without {flag}: {stderr}");
    assert!(
        stderr.lines().any(|line| {
            line == format!(
                "ERROR: [cache] Cache directory {} not found",
                root.join("homo_sapiens").display()
            )
        }),
        "without {flag}: {stderr}"
    );
}

/// `--refseq` selects the `_refseq` species directory of the cache root.
///
/// Documented: "`--refseq` and `--merged` select the cache, the species directory
/// suffixed `_refseq` or `_merged`, as Ensembl VEP does" (What vep-rs does not do,
/// RefSeq and merged caches). Ensembl VEP appends `_refseq` to the species
/// directory name under the flag (`CacheDir.pm:312`) and dies when the directory
/// is absent (`:318`).
#[test]
fn refseq_selects_the_refseq_species_directory() {
    assert_species_directory("--refseq", "homo_sapiens_refseq");
}

/// `--merged` selects the `_merged` species directory of the cache root
/// (`CacheDir.pm:312`, `:318`; the sentence of
/// [`refseq_selects_the_refseq_species_directory`]).
#[test]
fn merged_selects_the_merged_species_directory() {
    assert_species_directory("--merged", "homo_sapiens_merged");
}
