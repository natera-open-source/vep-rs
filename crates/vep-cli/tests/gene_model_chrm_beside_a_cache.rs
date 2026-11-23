// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row the `vep` binary writes for a transcript of a gene-model file
//! whose lines write the mitochondrial sequence `chrM`, read beside a GRCh37
//! cache: the run of the reproduction record of the entry
//! `file_beside_cache_chrm_transcripts_dropped`.
//!
//! Ensembl VEP hands each transcript record of such a file on under the
//! file's name less a leading `chr`, `M` (`Parser/GXF.pm:88-94`,
//! `BaseGXF.pm:282`), and keeps a record only where that name is a key of the
//! buffer's table of positions or resolves to one, through a synonym or with
//! `chr` added or removed (`AnnotationSource.pm:288-311`,
//! `BaseVEP.pm:697-734`). Beside a cache the table's key is the cache's `MT`,
//! and the cache's synonyms name `chrM` and not `M`, so no transcript of the
//! file is annotated. vep-rs names the file's `chrM` as the cache's `MT` and
//! writes the row.
//!
//! Every input is made up: the gene of
//! `tests/fixtures/gene-model-chrm-beside-cache`, the reference of
//! `tests/fixtures/mitochondrial-frame-check` and the record's line, which is
//! read from the reproduction record so that it is spelt in one file. The
//! cache is the committed GRCh37 slice of
//! `tests/fixtures/bayesdel-mitochondrial-contig`, which holds `MT` and a
//! synonyms table that links `MT` and `chrM`.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The entry's reproduction record, by its path from the repository's root.
const RECORD: &str =
    "scripts/divergences/reproductions/file_beside_cache_chrm_transcripts_dropped.toml";

/// The made-up gene-model file's name, which is its `SOURCE` and the name of
/// the column the tab format adds for it.
const MODEL_FILE: &str = "chrm_gene_synthetic.gtf.gz";

/// A file or directory of the repository, by its path from the root.
fn committed(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

fn fixture(path: &str) -> String {
    committed(&format!("crates/vep-cli/tests/fixtures/{path}"))
        .display()
        .to_string()
}

/// The lines of the record's `[input]` text, a VCF without its header.
fn record_input() -> String {
    let record = std::fs::read_to_string(committed(RECORD)).unwrap();
    let lines: Vec<&str> = record
        .lines()
        .skip_while(|line| *line != "[input]")
        .skip_while(|line| *line != "text = \"\"\"")
        .skip(1)
        .take_while(|line| *line != "\"\"\"")
        .collect();
    assert_eq!(lines.len(), 1, "{RECORD}: one input line");
    format!("{}\n", lines[0])
}

/// The rows the record states for vep-rs, each as its compared columns
/// joined by ` | `: the strings of the `vep_rs` list of its `[expected]`
/// table.
fn record_rows_of_vep_rs() -> Vec<String> {
    let record = std::fs::read_to_string(committed(RECORD)).unwrap();
    let line = record
        .lines()
        .skip_while(|line| *line != "[expected]")
        .find(|line| line.starts_with("vep_rs = ["))
        .unwrap_or_else(|| panic!("{RECORD}: no vep_rs list"));
    line.trim_start_matches("vep_rs = [")
        .trim_end_matches(']')
        .split("\", \"")
        .map(|row| row.trim_matches('"').to_string())
        .collect()
}

/// The items the record names under `code` of its `[implementation]` table.
fn record_code_items() -> Vec<String> {
    let record = std::fs::read_to_string(committed(RECORD)).unwrap();
    record
        .lines()
        .skip_while(|line| *line != "[implementation]")
        .skip_while(|line| *line != "code = [")
        .skip(1)
        .take_while(|line| *line != "]")
        .map(|line| {
            line.trim()
                .trim_end_matches(',')
                .trim_matches('"')
                .to_string()
        })
        .collect()
}

/// One run's column names, its rows as their tab-separated cells and the
/// lines of the warnings file it left.
struct Run {
    columns: Vec<String>,
    rows: Vec<Vec<String>>,
    warnings: Vec<String>,
}

/// Runs the binary on the record's input beside the GRCh37 cache slice and
/// over the made-up reference with `flags`, writing the tab format under
/// `dir`.
fn run(dir: &Path, flags: &[&str]) -> Run {
    let input = dir.join("record.vcf");
    std::fs::write(&input, format!("{VCF_HEADER}{}", record_input())).unwrap();
    let out = dir.join("out.txt");
    let warnings_file = format!("{}_warnings.txt", out.display());
    for stale in [out.display().to_string(), warnings_file.clone()] {
        let _ = std::fs::remove_file(stale);
    }
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args([
            "--json_cache",
            &fixture("bayesdel-mitochondrial-contig/json_cache"),
            "--fasta",
            &fixture("mitochondrial-frame-check/reference_synthetic.fa"),
            "--assembly",
            "GRCh37",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
            "--tab",
        ])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{flags:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    let cells = |line: &str| line.split('\t').map(str::to_string).collect::<Vec<_>>();
    let mut columns = Vec::new();
    let mut rows = Vec::new();
    for line in text.lines().filter(|line| !line.starts_with("##")) {
        match line.strip_prefix('#') {
            Some(names) => columns = cells(names),
            None => rows.push(cells(line)),
        }
    }
    let warnings = std::fs::read_to_string(warnings_file)
        .unwrap_or_default()
        .lines()
        .map(str::to_string)
        .collect();
    Run {
        columns,
        rows,
        warnings,
    }
}

/// The record is on `MT` at position 150, inside the one exon, 101 to 200,
/// of the made-up transcript `MITOT1`, which the file writes on `chrM`: base
/// 150 - 101 + 1 = 50 of a forward-strand transcript of a non-coding biotype.
/// The run writes that transcript's row under the file's name, every row the
/// run without the file writes, and no warning; the row is the one the
/// reproduction record states for vep-rs.
#[test]
fn a_chrm_gene_beside_a_cache_keeps_its_row_on_mt() {
    let dir = tempfile::tempdir().unwrap();
    let without_the_file = run(dir.path(), &[]);
    assert!(!without_the_file.rows.is_empty());

    let model = fixture(&format!("gene-model-chrm-beside-cache/{MODEL_FILE}"));
    let beside = run(dir.path(), &["--gtf", &model]);
    assert_eq!(beside.warnings, [] as [&str; 0]);
    assert_eq!(
        beside.columns,
        [
            "Uploaded_variation",
            "Location",
            "Allele",
            "Gene",
            "Feature",
            "Feature_type",
            "Consequence",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Amino_acids",
            "Codons",
            "Existing_variation",
            "IMPACT",
            "DISTANCE",
            "STRAND",
            "FLAGS",
            "SOURCE",
            MODEL_FILE,
        ]
    );

    let source = beside.columns.len() - 2;
    let (of_the_file, of_the_cache): (Vec<_>, Vec<_>) = beside
        .rows
        .iter()
        .partition(|row| row[source] == MODEL_FILE);
    let expected: Vec<String> = [
        "chrm_gene_exon_snv",
        "MT:150",
        "A",
        "MITOG1",
        "MITOT1",
        "Transcript",
        "non_coding_transcript_exon_variant",
        "50",
        "-",
        "-",
        "-",
        "-",
        "-",
        "MODIFIER",
        "-",
        "1",
        "-",
        MODEL_FILE,
        "-",
    ]
    .map(str::to_string)
    .to_vec();
    assert_eq!(of_the_file, [&expected]);

    // The cache's rows carry `-` in the two columns the file adds.
    let cache_rows: Vec<Vec<String>> = of_the_cache
        .iter()
        .map(|row| {
            assert_eq!(row[source..], ["-", "-"], "{row:?}");
            row[..source].to_vec()
        })
        .collect();
    assert_eq!(cache_rows, without_the_file.rows);

    let compared = [
        "Location",
        "Allele",
        "Gene",
        "Feature",
        "Consequence",
        "cDNA_position",
        "SOURCE",
    ]
    .map(|name| {
        let column = beside.columns.iter().position(|c| c == name).unwrap();
        expected[column].as_str()
    })
    .join(" | ");
    assert_eq!(record_rows_of_vep_rs(), [compared]);
}

/// Three functions place a file's `chrM` transcript on the cache's `MT`.
/// `SequenceNames::name` names the region `MT` as the file is read. Behind
/// it, `load` hands each region's name to
/// `ChromosomeSynonyms::source_chr_name`, which gives the name the cache
/// holds: a `chrM` through the cache's synonym, a bare `M` as `MT`. The row
/// stands while the first does, and falls with the first and either of the
/// other two taken back, so the record names all three.
#[test]
fn the_record_names_the_functions_that_place_a_chrm_transcript_on_mt() {
    assert_eq!(
        record_code_items(),
        [
            "crates/vep-cache-builder/src/gff3.rs::SequenceNames::name",
            "crates/vep-cli/src/gene_model.rs::load",
            "crates/vep-core/src/synonyms.rs::ChromosomeSynonyms::source_chr_name",
        ]
    );
}

/// The text of a page below `docs/`, its whitespace collapsed so that a
/// sentence wrapped over several lines is matched as one.
fn page(path: &str) -> String {
    std::fs::read_to_string(committed(&format!("docs/{path}")))
        .unwrap_or_else(|error| panic!("docs/{path}: {error}"))
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

/// The entry's page states the difference, the lines of Ensembl VEP's source
/// behind it, the row this file pins and the runs that restore the row in
/// Ensembl VEP; the index row states the entry, its tier and its counts.
#[test]
fn the_page_and_the_index_row_state_the_dropped_chrm_transcripts() {
    let path = "divergences/file-beside-cache-chrm-transcripts-dropped.md";
    let text = page(path);
    for phrase in [
        "Beside a cache, Ensembl VEP writes no row for a transcript of a GTF or GFF3 whose lines \
         name the mitochondrial sequence `chrM`, and vep-rs writes the transcript's rows, with \
         the file's name in `SOURCE`.",
        "It then looks each transcript up under the file's name less its `chr`, `M`, a name the \
         synonyms table does not hold.",
        "The transcript is dropped without a warning.",
        "vep-rs reads a file's `chrM` as the cache's `MT` and annotates the file's transcripts \
         beside the cache's.",
        "`Bio/EnsEMBL/IO/Parser/GXF.pm`, lines 88-94, `get_seqname`",
        "`Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, line 282",
        "`InputBuffer.pm`, lines 463-498, `min_max`: the buffer keeps the first table it is \
         handed.",
        "`AnnotationSource.pm`, lines 288-311, `_check_overlap`",
        "`BaseVEP.pm`, lines 689-737, `get_source_chr_name`",
        "The class is documented only. No comparison sets its rows aside",
        "hold 1,808 lines outside their comments, each on chromosome 21, so no committed corpus \
         carries the shape",
        "chrm_gene_exon_snv MT:150 A MITOG1 MITOT1 Transcript non_coding_transcript_exon_variant \
         50 - - - - - IMPACT=MODIFIER;STRAND=1;SOURCE=chrm_gene_synthetic.gtf.gz",
        "| `chrM` | the cache's table, named again with `--synonyms` | none |",
        "| `chrM` | the cache's table and a line `MT M`, with `--synonyms` | written |",
        "| `MT` | the cache's table | written |",
        "so it is base 150 - 101 + 1 = 50 of the transcript",
        "The same input record and file read beside the GRCh38 cache give the same outputs: the \
         row of `MITOT1` in vep-rs, none in Ensembl VEP 116.3, and the row in both where the \
         file's lines name `MT`.",
    ] {
        assert!(text.contains(phrase), "docs/{path} lacks {phrase:?}");
    }

    let index = page("intended-divergences.md");
    for phrase in [
        "| [`file_beside_cache_chrm_transcripts_dropped`]\
         (divergences/file-beside-cache-chrm-transcripts-dropped.md) | Beside a cache, Ensembl \
         VEP reads the lines of a GTF or GFF3 that write the mitochondrial sequence `chrM`,",
        "vep-rs reads the file's `chrM` as the cache's `MT` and writes those rows, the ones \
         Ensembl VEP writes when the file names the sequence `MT` | 3 · 0 of 1,808 table rows · \
         5 flag sets, 0 corpora | 116.3 |",
    ] {
        assert!(
            index.contains(phrase),
            "docs/intended-divergences.md lacks {phrase:?}"
        );
    }
}
