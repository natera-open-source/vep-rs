// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two runs of the documented difference
//! `refseq_edited_model_shifted_insertion_bases`: its record's own run under
//! `--refseq --hgvs --shift_3prime 1` and its further run under
//! `--refseq --hgvs`, each over the record's three insertions on RefSeq models
//! edited to the genome. One test runs the `vep` binary on the record's input
//! under the flags of each run and holds the rows the record keeps for vep-rs;
//! the other holds the commands and the table of the entry's page to the
//! record's input and to the rows it keeps for both engines. The input is read
//! from the record under `scripts/divergences/reproductions`, so each line is
//! spelt in one file.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

const RECORD: &str =
    "scripts/divergences/reproductions/refseq_edited_model_shifted_insertion_bases.toml";

const PAGE: &str = "docs/divergences/refseq-edited-model-shifted-insertion-bases.md";

/// The cache and the reference the record names, both of one committed corpus.
const JSON_CACHE: &str = "tests/golden/116/GRCh37-refseq/cache_root/homo_sapiens_refseq/116_GRCh37";
const FASTA: &str = "tests/golden/116/GRCh37-refseq/reference.fa.gz";

/// The two edited models whose rows the record compares: `NM_006933.7` on the
/// forward strand and `NM_001395858.1` on the reverse strand.
const FEATURES: [&str; 2] = ["NM_001395858.1", "NM_006933.7"];

/// The record's runs, its own first: the table that holds a run's flags and
/// the table that holds the rows of both engines for it.
const RUNS: [(&str, &str); 2] = [("[run]", "[expected]"), ("[also.run]", "[also.expected]")];

/// A file or directory of the repository, by its path from the root.
fn committed(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

/// The text of the entry's record.
fn record() -> String {
    std::fs::read_to_string(committed(RECORD)).unwrap()
}

/// The lines of the table of `record` that `header` opens, up to the next line
/// that opens a table.
fn table<'a>(record: &'a str, header: &str) -> Vec<&'a str> {
    let lines: Vec<&str> = record
        .lines()
        .skip_while(|line| *line != header)
        .skip(1)
        .take_while(|line| !line.starts_with('['))
        .collect();
    assert!(!lines.is_empty(), "{RECORD} holds no table {header}");
    lines
}

/// The strings of the list `key` of the table `header`, which closes on the
/// first line that ends with a bracket.
fn list(record: &str, header: &str, key: &str) -> Vec<String> {
    let opening = format!("{key} = [");
    let lines = table(record, header);
    let start = lines
        .iter()
        .position(|line| line.starts_with(&opening))
        .unwrap_or_else(|| panic!("{RECORD}: {header} holds no list {key}"));
    let mut text = String::new();
    for line in &lines[start..] {
        text.push_str(line);
        if line.trim_end().ends_with(']') {
            break;
        }
    }
    text.split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect()
}

/// The input lines of the record: the `text` of its `[input]` table, a VCF
/// without its header, which every run of the record reads.
fn input_lines(record: &str) -> Vec<&str> {
    table(record, "[input]")
        .into_iter()
        .skip_while(|line| *line != "text = \"\"\"")
        .skip(1)
        .take_while(|line| *line != "\"\"\"")
        .collect()
}

/// The rows of the two models that the binary writes for the record's input
/// under `flags`, each as the record's `[compare]` table renders it (Location,
/// Allele, Feature, `BAM_EDIT`, HGVSc and `HGVS_OFFSET`), sorted as the
/// registry tool sorts them. The run ends with exit code 0.
fn compared_rows(dir: &Path, name: &str, record: &str, flags: &[String]) -> Vec<String> {
    let input = dir.join("records.vcf");
    let records: String = input_lines(record)
        .iter()
        .map(|line| format!("{line}\n"))
        .collect();
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join(format!("{name}.txt"));
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(committed(JSON_CACHE))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh37"])
        .args(["--format", "vcf"])
        .args(["--force_overwrite", "--quiet", "--no_stats"])
        .arg("--fasta")
        .arg(committed(FASTA))
        .args(flags);
    let output = cmd.stdin(Stdio::null()).output().unwrap();
    assert!(
        output.status.success(),
        "{flags:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    let mut rows: Vec<String> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<&str>>())
        .filter(|cells| FEATURES.contains(&cells[4]))
        .map(|cells| {
            let extra = |key: &str| {
                cells[13]
                    .split(';')
                    .find_map(|item| item.strip_prefix(key)?.strip_prefix('='))
                    .unwrap_or("-")
            };
            [
                cells[1],
                cells[2],
                cells[4],
                extra("BAM_EDIT"),
                extra("HGVSc"),
                extra("HGVS_OFFSET"),
            ]
            .join(" | ")
        })
        .collect();
    rows.sort();
    rows
}

/// vep-rs spells the bases of a shifted insertion as the shift rotates them on
/// every row of an edited model, under `--hgvs` alone and under
/// `--shift_3prime 1`: the coding row of `NM_006933.7`, its 3' UTR row and the
/// 5' UTR row of the reverse-strand `NM_001395858.1`. The `Allele` cell of the
/// coding row follows the shift under the flag and is the written one without
/// it. Ensembl VEP 116.3 overwrites the notation's alleles with the alleles as
/// the record writes them on a model that carries an `_rna_edit` attribute
/// (`TranscriptVariationAllele.pm:1508-1513`) and rotates them only where
/// `codon` runs `shift_feature_seqs` (line 820), on a coding row under the
/// flag; the record keeps its rows of both runs. The rows of each run are the
/// ones the record keeps for vep-rs, and the record's two runs name these
/// flags, this cache and this reference.
#[test]
fn both_runs_of_the_record_spell_the_rotated_bases_on_every_row() {
    let record = record();
    let cache = table(&record, "[cache]");
    assert!(cache.contains(&format!("json_cache = \"{JSON_CACHE}\"").as_str()));
    assert!(cache.contains(&format!("fasta = \"{FASTA}\"").as_str()));
    let expected: [(&[&str], [&str; 3]); 2] = [
        (
            &["--refseq", "--hgvs", "--shift_3prime", "1"],
            [
                "21:16031362-16031363 | TC | NM_001395858.1 | OK | NM_001395858.1:c.-57_-56insAG | -1",
                "21:35467700-35467701 | AG | NM_006933.7 | OK | NM_006933.7:c.204_205insAG | 1",
                "21:35471002-35471003 | GA | NM_006933.7 | OK | NM_006933.7:c.*1349_*1350insAG | 1",
            ],
        ),
        (
            &["--refseq", "--hgvs"],
            [
                "21:16031362-16031363 | TC | NM_001395858.1 | OK | NM_001395858.1:c.-57_-56insAG | -1",
                "21:35467700-35467701 | GA | NM_006933.7 | OK | NM_006933.7:c.204_205insAG | 1",
                "21:35471002-35471003 | GA | NM_006933.7 | OK | NM_006933.7:c.*1349_*1350insAG | 1",
            ],
        ),
    ];
    let dir = tempfile::tempdir().unwrap();
    for (place, ((run, kept), (flags, rows))) in RUNS.iter().zip(expected).enumerate() {
        let recorded = list(&record, run, "flags");
        assert_eq!(recorded, flags, "{run}");
        assert_eq!(list(&record, kept, "vep_rs"), rows, "{kept}");
        assert_eq!(
            compared_rows(dir.path(), &format!("run_{place}"), &record, &recorded),
            rows,
            "{flags:?}"
        );
    }
}

/// The notation of the row at `location` among the rows `kept` holds for
/// `engine`: its HGVSc cell without the transcript's name.
fn notation(record: &str, kept: &str, engine: &str, location: &str) -> String {
    let rows = list(record, kept, engine);
    let row = rows
        .iter()
        .find(|row| row.split(" | ").next() == Some(location))
        .unwrap_or_else(|| panic!("{RECORD}: {kept} holds no {engine} row at {location}"));
    let hgvsc = row.split(" | ").nth(4).unwrap();
    hgvsc.split_once(':').unwrap().1.to_string()
}

/// The entry's page gives its reader every record the reproduction runs: each
/// input line of the record is a `printf` line of the commands of `Try it
/// yourself`, whose `grep` names both models, and the table under them holds
/// the notations the record keeps for both engines over its two runs. A record
/// whose notation differs between the runs has one row a run, named by the
/// run's flags; every other record has one row, named by its line, that stands
/// for either run.
#[test]
fn the_page_table_names_every_record_and_run_of_the_reproduction() {
    let record = record();
    let page = std::fs::read_to_string(committed(PAGE)).unwrap();
    let has = |wanted: &dyn Fn(&str) -> bool| page.lines().any(wanted);
    let models = page
        .lines()
        .find(|line| line.starts_with("grep ") && line.contains(" out.txt "))
        .unwrap_or_else(|| panic!("{PAGE} holds no grep over out.txt"));
    for feature in FEATURES {
        assert!(
            models.contains(feature),
            "{PAGE}: {models} leaves out {feature}"
        );
    }
    let lines = input_lines(&record);
    assert_eq!(lines.len(), 3, "{RECORD}");
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        let (chrom, pos, reference, alternate) = (cells[0], cells[1], cells[3], cells[4]);
        let written = format!("printf '{}\\n' >> ins.vcf", line.replace('\t', "\\t"));
        assert!(
            has(&|row| row == written),
            "{PAGE} lacks the line {written}"
        );
        let location = format!("{chrom}:{pos}-{}", pos.parse::<u64>().unwrap() + 1);
        let of_run = |kept: &str| {
            ["ensembl_vep", "vep_rs"].map(|engine| notation(&record, kept, engine, &location))
        };
        let (own, further) = (of_run(RUNS[0].1), of_run(RUNS[1].1));
        if own == further {
            let opening = format!("| `{chrom} {pos} {reference} {alternate}` (");
            let closing = format!("), either run | `{}` | `{}` |", own[0], own[1]);
            assert!(
                has(&|row| row.starts_with(&opening) && row.ends_with(&closing)),
                "{PAGE} lacks a row {opening}...{closing}"
            );
        } else {
            for ((run, _), notations) in RUNS.iter().zip([own, further]) {
                let flags = list(&record, run, "flags").join(" ");
                let wanted = format!("| `{flags}` | `{}` | `{}` |", notations[0], notations[1]);
                assert!(has(&|row| row == wanted), "{PAGE} lacks the row {wanted}");
            }
        }
    }
}
