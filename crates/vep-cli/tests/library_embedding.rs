// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The proof behind `docs/developers/embedding.md`: the re-exported surface of
//! `vep_cli` writes the bytes the `vep` binary writes. One `Annotator` built
//! from the binary's own command line annotates each record of a corpus, and
//! the `JsonOutputFormatter` renders it as `vep --json` does.

mod common;

use vep_cli::vcf_parser::parse_vcf_line;
use vep_cli::{Annotator, Config, InputVariant, JsonOutputFormatter};
use vep_io::output::fields::FieldOptions;

use common::{corpora, default_run, input_records, output_path, run_vep};

/// The JSON rows the library writes for the GRCh37 corpus equal, byte for byte,
/// the file `vep --json` writes for it from the same command line: the binary's
/// reader parses each data line, so the comparison proves the annotator and the
/// formatter. The binary writes no JSON object for a record whose span exceeds
/// `--max_sv_size`, so the caller leaves those records out, as the page says;
/// their number is the number of `--max_sv_size` warnings the binary wrote, and
/// every other record of the corpus has its row.
#[test]
fn annotator_json_rows_equal_the_vep_binarys_json_output_on_the_corpus() {
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh37")
        .expect("tests/golden/116/GRCh37");
    let run = default_run(corpus);
    let out_dir = tempfile::tempdir().unwrap();
    let binary_rows = run_vep(&run, "json", out_dir.path(), &[]);
    let warnings = std::fs::read_to_string(format!(
        "{}_warnings.txt",
        output_path(&run, "json", out_dir.path()).display()
    ))
    .unwrap();
    let oversize_warnings = warnings
        .lines()
        .filter(|l| l.contains("is bigger than --max_sv_size"))
        .count();

    // The binary's command line without its input and output files: the cache,
    // the species and the corpus flags, which name the assembly.
    let mut argv = vec!["vep".to_string(), "--offline".to_string()];
    argv.extend(corpus.cache_args());
    argv.extend(["--species".to_string(), "homo_sapiens".to_string()]);
    argv.extend(run.flags.iter().cloned());
    let config = Config::from_argv(argv).unwrap();
    let max_sv_size = config.max_sv_size.expect("the default limit");
    let annotator = Annotator::from_config(&config).unwrap();
    let formatter = JsonOutputFormatter::new(
        Default::default(),
        annotator.plugin_header_fields(),
        Some(corpus.assembly.clone()),
    );

    let records = input_records(corpus);
    let mut library_rows = String::new();
    let mut oversize = 0;
    for line in &records {
        let mut record = parse_vcf_line(line, true).unwrap();
        if record
            .iter()
            .any(|v| v.end.saturating_sub(v.start) > max_sv_size)
        {
            oversize += 1;
            continue;
        }
        annotator.annotate_batch(&mut record);
        let alleles: Vec<&InputVariant> = record.iter().collect();
        let value = formatter.format_record_value(&alleles);
        library_rows.push_str(&serde_json::to_string(&value).unwrap());
        library_rows.push('\n');
    }
    assert_eq!(oversize, oversize_warnings);
    assert!(oversize > 0, "the corpus carries records past the limit");
    assert_eq!(library_rows.lines().count(), records.len() - oversize);

    let differing = library_rows
        .lines()
        .zip(binary_rows.lines())
        .enumerate()
        .filter(|(_, (library_row, binary_row))| library_row != binary_row)
        .map(|(i, (library_row, binary_row))| {
            format!("row {i}\n  library: {library_row}\n  binary:  {binary_row}")
        })
        .take(3)
        .collect::<Vec<_>>()
        .join("\n");
    assert!(differing.is_empty(), "{differing}");
    assert_eq!(library_rows, binary_rows);
}

/// The same proof over the corpus whose run loads the LoFTEE built-in: the
/// library's `Annotator::from_config` loads the plugin from the command line
/// and attaches the exon and intron context LoFTEE reads to every transcript
/// before the plugin runs, so the `lof` and `lof_info` fields of its rows equal
/// the binary's, byte for byte, with the confidence spread the corpus carries.
#[test]
fn annotator_json_rows_equal_the_vep_binarys_json_output_with_the_loftee_plugin() {
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-plugin-loftee")
        .expect("tests/golden/116/GRCh38-plugin-loftee");
    let run = default_run(corpus);
    let out_dir = tempfile::tempdir().unwrap();
    let binary_rows = run_vep(&run, "json", out_dir.path(), &[]);

    let fasta = corpus
        .dir
        .join(corpus.manifest["fasta"].as_str().expect("a corpus FASTA"));
    let mut argv = vec!["vep".to_string(), "--offline".to_string()];
    argv.extend(corpus.cache_args());
    argv.extend([
        "--species".to_string(),
        "homo_sapiens".to_string(),
        "--assembly".to_string(),
        corpus.assembly.clone(),
    ]);
    argv.extend(run.flags.iter().cloned());
    argv.extend(["--fasta".to_string(), fasta.display().to_string()]);
    let config = Config::from_argv(argv).unwrap();
    let annotator = Annotator::from_config(&config).unwrap();
    let formatter = JsonOutputFormatter::new(
        Default::default(),
        annotator.plugin_header_fields(),
        Some(corpus.assembly.clone()),
    );

    let records = input_records(corpus);
    let mut library_rows = String::new();
    for line in &records {
        let mut record = parse_vcf_line(line, true).unwrap();
        annotator.annotate_batch(&mut record);
        let alleles: Vec<&InputVariant> = record.iter().collect();
        let value = formatter.format_record_value(&alleles);
        library_rows.push_str(&serde_json::to_string(&value).unwrap());
        library_rows.push('\n');
    }
    assert_eq!(library_rows.lines().count(), records.len());
    for verdict in ["\"lof\":\"HC\"", "\"lof\":\"LC\""] {
        assert!(
            library_rows.contains(verdict),
            "the library wrote no {verdict} field"
        );
    }

    let differing = library_rows
        .lines()
        .zip(binary_rows.lines())
        .enumerate()
        .filter(|(_, (library_row, binary_row))| library_row != binary_row)
        .map(|(i, (library_row, binary_row))| {
            format!("row {i}\n  library: {library_row}\n  binary:  {binary_row}")
        })
        .take(3)
        .collect::<Vec<_>>()
        .join("\n");
    assert!(differing.is_empty(), "{differing}");
    assert_eq!(library_rows, binary_rows);
}

/// The corpus whose cache names the alignment file its models were edited
/// against: resolving it switches `--use_transcript_ref` and the `BAM_EDIT`
/// field on, on the annotator's copy of the configuration. The caller's
/// configuration keeps both off and the annotator returns both on. The field
/// options that carry the annotator's two values write the rows of
/// `vep --json` byte for byte, where the options that carry the caller's write
/// no `given_ref` and no `used_ref`.
#[test]
fn annotator_returns_the_flags_a_bam_edited_cache_settles() {
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-refseq")
        .expect("tests/golden/116/GRCh38-refseq");
    let run = default_run(corpus);
    let out_dir = tempfile::tempdir().unwrap();
    let binary_rows = run_vep(&run, "json", out_dir.path(), &[]);

    let fasta = corpus
        .dir
        .join(corpus.manifest["fasta"].as_str().expect("a corpus FASTA"));
    let mut argv = vec!["vep".to_string(), "--offline".to_string()];
    argv.extend(corpus.cache_args());
    argv.extend([
        "--species".to_string(),
        "homo_sapiens".to_string(),
        "--assembly".to_string(),
        corpus.assembly.clone(),
    ]);
    argv.extend(run.flags.iter().cloned());
    argv.extend(["--fasta".to_string(), fasta.display().to_string()]);
    let config = Config::from_argv(argv).unwrap();
    let max_sv_size = config.max_sv_size.expect("the default limit");
    assert!(!config.use_transcript_ref, "the command line names no flag");
    assert!(!config.bam_edited);

    let annotator = Annotator::from_config(&config).unwrap();
    assert!(annotator.use_transcript_ref());
    assert!(annotator.bam_edited());
    assert_eq!(
        annotator.use_transcript_ref(),
        annotator.resources().effects_config.use_transcript_ref
    );

    let records = input_records(corpus);
    let rows_under = |use_transcript_ref: bool, bam_edited: bool| {
        let formatter = JsonOutputFormatter::new(
            FieldOptions {
                refseq: config.refseq,
                use_transcript_ref,
                bam_edited,
                ..Default::default()
            },
            annotator.plugin_header_fields(),
            Some(corpus.assembly.clone()),
        );
        let mut rows = String::new();
        for line in &records {
            let mut record = parse_vcf_line(line, true).unwrap();
            if record
                .iter()
                .any(|v| v.end.saturating_sub(v.start) > max_sv_size)
            {
                continue;
            }
            annotator.annotate_batch(&mut record);
            let alleles: Vec<&InputVariant> = record.iter().collect();
            let value = formatter.format_record_value(&alleles);
            rows.push_str(&serde_json::to_string(&value).unwrap());
            rows.push('\n');
        }
        rows
    };
    let settled = rows_under(annotator.use_transcript_ref(), annotator.bam_edited());
    let differing = settled
        .lines()
        .zip(binary_rows.lines())
        .enumerate()
        .filter(|(_, (library_row, binary_row))| library_row != binary_row)
        .map(|(i, (library_row, binary_row))| {
            format!("row {i}\n  library: {library_row}\n  binary:  {binary_row}")
        })
        .take(3)
        .collect::<Vec<_>>()
        .join("\n");
    assert!(differing.is_empty(), "{differing}");
    assert_eq!(settled, binary_rows);

    let unsettled = rows_under(config.use_transcript_ref, config.bam_edited);
    for key in ["\"given_ref\":", "\"used_ref\":"] {
        assert!(settled.contains(key), "the settled options write no {key}");
        assert!(
            !unsettled.contains(key),
            "the caller's own flags write {key}"
        );
    }
}
