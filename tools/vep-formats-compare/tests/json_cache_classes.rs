// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `data/json_cache_classes.tsv` names the classes for which the command line opens the
//! candidate's JSON cache, each under the table it builds for the class.

mod support;

use std::path::{Path, PathBuf};

use support::{golden, run, scratch};
use vep_formats_compare::cache_classes::{classes, classes_reading, CacheTable};
use vep_formats_compare::classes::tables;

/// The formats of an Ensembl VEP output a comparison reads as text, each with the committed
/// reference output of `GRCh38-release116`.
const FORMATS: [(&str, &str); 4] = [
    ("default", "default.txt.gz"),
    ("tab", "tab.txt.gz"),
    ("vcf", "vcf.vcf.gz"),
    ("json", "json.jsonl.gz"),
];

/// The class whose table is built beside `--reference-fasta` alone.
const SEQUENCE_EDGE_CLASS: &str = "hgvs_shift_sequence_edge_window";

fn corpus() -> PathBuf {
    golden().join("GRCh38-release116")
}

/// The names of `data/classes.tsv`, in the file's order.
fn names_of_the_class_file() -> Vec<String> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("data/classes.tsv");
    let text = std::fs::read_to_string(path).expect("the class file is read");
    text.lines()
        .skip(1)
        .filter_map(|line| line.split('\t').next())
        .map(str::to_string)
        .collect()
}

/// A directory that is no cache, and an indexed FASTA of one made-up sequence.
struct Planted {
    directory: PathBuf,
    no_cache: PathBuf,
    fasta: PathBuf,
}

fn planted(name: &str) -> Planted {
    let directory = scratch(name);
    let no_cache = directory.join("no-cache");
    std::fs::create_dir(&no_cache).expect("the directory is made");
    let fasta = directory.join("reference.fa");
    let bases = "N".repeat(60);
    std::fs::write(&fasta, format!(">21\n{bases}\n")).expect("the FASTA is written");
    std::fs::write(directory.join("reference.fa.fai"), "21\t60\t4\t60\t61\n")
        .expect("the index is written");
    Planted {
        directory,
        no_cache,
        fasta,
    }
}

/// Whether a comparison of the committed output of `format` with itself under
/// `--masks class` opens the cache it is given. The directory named holds no `transcripts/`,
/// which each reader of the cache refuses with the same line before it reads a shard; a
/// comparison that opens no cache ends with exit code 0, the two outputs being one file.
fn opens_the_cache(planted: &Planted, format: (&str, &str), class: &str, fasta: bool) -> bool {
    let output = corpus().join("expected").join(format.1);
    let mut args: Vec<String> = [
        "--format",
        format.0,
        "--expected",
        &output.to_string_lossy(),
        "--actual",
        &output.to_string_lossy(),
        "--input",
        &corpus().join("variants.vcf").to_string_lossy(),
        "--masks",
        class,
        "--json-cache",
        &planted.no_cache.to_string_lossy(),
    ]
    .iter()
    .map(|word| word.to_string())
    .collect();
    if fasta {
        args.push("--reference-fasta".to_string());
        args.push(planted.fasta.to_string_lossy().into_owned());
    }
    let ran = run(&planted.directory, &args);
    let refusal = format!(
        "--json-cache: {} has no transcripts/ directory",
        planted.no_cache.display()
    );
    let refused = ran.err().contains(&refusal);
    assert_eq!(
        ran.code,
        if refused { 2 } else { 0 },
        "{class} in the {} format: {}",
        format.0,
        ran.err()
    );
    refused
}

#[test]
fn the_file_names_each_class_under_the_table_built_for_it() {
    let held: Vec<(&str, &str)> = classes()
        .iter()
        .map(|(class, table)| (*class, table.name()))
        .collect();
    assert_eq!(
        held,
        [
            ("breakend_mate_local_read", "structural_variants"),
            ("bnd_own_chromosome_rows_lost", "structural_variants"),
            ("hgvs_shift_sequence_edge_window", "sequence_edges"),
            (
                "hgvsp_stop_named_on_incomplete_coding_sequence",
                "partial_codon_models"
            ),
            (
                "structural_deletion_stop_lost_without_stop_codon",
                "partial_codon_models"
            ),
            (
                "stop_codon_completed_by_polyadenylation",
                "partial_codon_models"
            ),
            (
                "stop_codon_overlap_without_stop_codon",
                "partial_codon_models"
            ),
            (
                "stop_codon_behind_last_residue_without_stop_codon",
                "partial_codon_models"
            ),
            (
                "start_codon_overlap_without_start_codon",
                "partial_codon_models"
            ),
            (
                "hgvsp_frameshift_count_at_the_first_residue",
                "partial_codon_models"
            ),
            ("start_cooccurrence_swap", "partial_codon_models"),
            ("stop_codon_unshifted_span", "partial_codon_models"),
            (
                "stop_retained_without_stop_codon_test",
                "partial_codon_models"
            ),
            ("stop_codon_overlap_genomic_span", "partial_codon_models"),
            (
                "stop_retained_three_base_stop_codon_test",
                "partial_codon_models"
            ),
            (
                "stop_codon_removed_without_stop_codon_test",
                "partial_codon_models"
            ),
            ("stop_codon_overlap_insertion_range", "partial_codon_models"),
            (
                "plugin_value_follows_splice_term_past_the_record_span",
                "partial_codon_models"
            ),
            (
                "plugin_value_follows_stop_codon_term_pair",
                "partial_codon_models"
            ),
            ("refseq_failed_edit_model_hgvsc_offset", "length_edits"),
            ("organelle_peptide_paths_translation_table", "codon_tables"),
            (
                "downstream_utr_frame_past_coding_end",
                "downstream_sequences"
            ),
            ("downstream_start_past_own_codon", "downstream_sequences"),
            (
                "plugin_hgvs_shifted_cds_coordinates",
                "downstream_sequences"
            ),
            ("hgvs_shift_shared_across_alleles", "downstream_sequences"),
            (
                "refseq_edited_model_inserted_span_codons",
                "downstream_sequences"
            ),
            (
                "hgvsp_stop_codon_protein_comparison",
                "downstream_sequences"
            ),
            (
                "hgvsp_frameshift_count_past_a_stop_read_as_a_residue",
                "downstream_sequences"
            ),
            (
                "stop_codon_settled_before_unknown_bases",
                "downstream_sequences"
            ),
            (
                "hgvs_clipped_reference_against_whole_alternate",
                "transcript_edges"
            ),
            (
                "hgvs_insertion_left_where_the_shift_stops",
                "transcript_edges"
            ),
            ("hgvsp_shifted_indel_written_span_gate", "transcript_edges"),
            ("hgvsp_intron_end_insertion", "intron_ends"),
            ("hgvsp_intron_start_insertion", "intron_ends"),
            ("maxentscan_overlapping_windows", "splice_windows"),
        ]
    );
    assert_eq!(
        classes_reading(CacheTable::SequenceEdges),
        [SEQUENCE_EDGE_CLASS]
    );
}

/// The structural-variant comparator opens a cache it is given for each class of its table,
/// the one that asks nothing of the cache among them, so the comparison of such a class
/// opens the cache though the file has no line for it.
#[test]
fn a_comparison_opens_the_cache_for_the_classes_of_the_file_and_for_no_other_class() {
    let planted = planted("json-cache-classes");
    let listed: Vec<&str> = classes().iter().map(|(class, _)| *class).collect();
    let structural = tables().structural_variant_classes();
    let names = names_of_the_class_file();
    assert!(names.len() > listed.len());
    for class in &names {
        let class = class.as_str();
        let expected = listed.contains(&class) || structural.contains(&class);
        assert_eq!(
            opens_the_cache(&planted, FORMATS[0], class, true),
            expected,
            "{class}"
        );
    }
    for class in &listed {
        for format in &FORMATS[1..] {
            assert!(
                opens_the_cache(&planted, *format, class, true),
                "{class} in the {} format",
                format.0
            );
        }
    }
    for format in FORMATS {
        assert!(
            !opens_the_cache(&planted, format, SEQUENCE_EDGE_CLASS, false),
            "the sequence-edge table is built beside --reference-fasta alone"
        );
    }
}
