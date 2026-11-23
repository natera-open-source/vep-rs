// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary against reference bytes.
//!
//! `tests/fixtures/caches/` holds small Storable caches written by perl 5.34.0 with
//! Storable 3.23, and `tests/fixtures/expected/<case>/` the reference conversion of one
//! of them, made on the same perl with JSON::XS 4.04: the arguments (`args`, one a
//! line), standard output, standard error, the exit code and every file of the output
//! directory (`out/`), with the files an earlier run had left there (`pre/`). Standard
//! error is held as the tool prints it, under the first two rules of the README's
//! "Where the tool does not write what Perl writes" (no ` at <file> line <n>.` after a
//! message, no line of `gzip -dc` itself), with the program of a usage line as
//! `<program>`. Each test runs the binary with the same arguments in a directory that
//! holds the same cache and asks for the same bytes.

mod common;

use common::{copy_tree, fixtures, run_in, tree, Scratch, BINARY};

fn lines(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn case(id: &str) {
    let expected = fixtures().join("expected").join(id);
    let arguments: Vec<String> = std::fs::read_to_string(expected.join("args"))
        .expect("args")
        .lines()
        .map(str::to_string)
        .collect();
    let scratch = Scratch::new(id);
    if let Some(cache) = arguments.first() {
        // The cache is the first component of the first argument.
        let top = cache.split('/').next().expect("a path");
        let source = fixtures().join("caches").join(top);
        if source.is_dir() {
            copy_tree(&source, &scratch.path().join(top));
        }
        // A directory without a file is not kept with the fixtures: the list names them.
        let empty = std::fs::read_to_string(fixtures().join("caches/EMPTY_DIRECTORIES"))
            .expect("EMPTY_DIRECTORIES");
        for directory in empty
            .lines()
            .filter(|line| line.split('/').next() == Some(top))
        {
            std::fs::create_dir_all(scratch.path().join(directory)).expect("empty directory");
        }
    }
    if let Some(output) = arguments.get(1) {
        let output = scratch.path().join(output);
        std::fs::create_dir_all(&output).expect("output directory");
        if expected.join("pre").is_dir() {
            copy_tree(&expected.join("pre"), &output);
        }
    }
    let argument_refs: Vec<&str> = arguments.iter().map(String::as_str).collect();
    let run = run_in(scratch.path(), &argument_refs);

    let exit: i32 = std::fs::read_to_string(expected.join("exit"))
        .expect("exit")
        .trim()
        .parse()
        .expect("an exit code");
    assert_eq!(
        run.status.code(),
        Some(exit),
        "{id}: exit code; stderr: {}",
        lines(&run.stderr)
    );
    let stdout = std::fs::read(expected.join("stdout")).expect("stdout");
    assert_eq!(lines(&run.stdout), lines(&stdout), "{id}: standard output");
    assert_eq!(run.stdout, stdout, "{id}: standard output bytes");
    // The usage line names the program as it was run.
    let recorded = std::fs::read(expected.join("stderr")).expect("stderr");
    let marker = b"<program>";
    let stderr = match recorded
        .windows(marker.len())
        .position(|window| window == marker)
    {
        Some(at) => [
            &recorded[..at],
            BINARY.as_bytes(),
            &recorded[at + marker.len()..],
        ]
        .concat(),
        None => recorded,
    };
    assert_eq!(lines(&run.stderr), lines(&stderr), "{id}: standard error");
    assert_eq!(run.stderr, stderr, "{id}: standard error bytes");

    if let Some(output) = arguments.get(1) {
        let written = tree(&scratch.path().join(output));
        let wanted = tree(&expected.join("out"));
        let names = |files: &std::collections::BTreeMap<String, Vec<u8>>| {
            files.keys().cloned().collect::<Vec<_>>()
        };
        assert_eq!(
            names(&written),
            names(&wanted),
            "{id}: files of the output directory"
        );
        for (name, bytes) in &wanted {
            assert!(
                written[name] == *bytes,
                "{id}: {name} differs ({} bytes written, {} wanted)",
                written[name].len(),
                bytes.len()
            );
        }
    }
}

macro_rules! cases {
    ($($name:ident => $id:literal,)*) => {
        $(
            #[test]
            fn $name() {
                case($id);
            }
        )*

        #[test]
        fn every_expected_directory_has_its_test() {
            let mut listed: Vec<&str> = vec![$($id),*];
            listed.sort_unstable();
            let mut present: Vec<String> = std::fs::read_dir(fixtures().join("expected"))
                .expect("expected directory")
                .map(|entry| entry.expect("entry").file_name().to_string_lossy().into_owned())
                .collect();
            present.sort();
            assert_eq!(present, listed);
        }
    };
}

cases! {
    // Prediction markers, regulatory shards, mapper pairs and transcript sets.
    marker_without_a_source_line_is_dropped => "polyphen-only",
    markers_with_their_source_lines_stay => "both-versions",
    markers_without_a_written_matrix_are_dropped => "no-matrices",
    regulatory_objects_activity_codes_motifs_and_the_matrices_table => "regulatory",
    empty_regulatory_shard_enables_the_written_tree => "empty-regulatory",
    input_regulatory_marker_without_shards_is_dropped => "missing-regulatory",
    activity_outside_the_five_values_is_rejected_by_name => "badactivity",
    matrix_that_differs_from_its_first_copy_is_rejected_by_name => "matrixdiffers",
    corrupt_regulatory_shard_leaves_no_metadata_of_an_earlier_run => "corrupt-regulatory",
    motif_without_a_binding_matrix_is_rejected_by_name => "no-matrix",
    edited_mapper_pairs_are_read_from_the_storable_mapper => "folded-mapper",
    mapper_object_without_pairs_falls_back_to_the_exons => "empty-mapper",
    refseq_directory_converts_with_its_source_type_and_bam_line => "sets-refseq",
    merged_directory_converts_with_its_source_type => "sets-merged",
    ensembl_directory_converts_with_its_source_type => "sets-ensembl",
    // A biotype held as an object: its name, in every kind a cache can store one as.
    biotype_stored_as_an_object_is_written_as_its_name => "biotype-object",
    biotype_object_is_written_by_the_kind_of_its_name => "biotype-shapes",
    // The sequence side file: one chromosome directory a shape of the slice, its length
    // and its circular flag (Storable writes a whole double as an integer), and a length
    // under the UTF-8 flag that is not UTF-8, before and after the chromosome states one.
    sequence_regions_list_a_linear_a_circular_and_an_unstated_sequence => "regions",
    failed_run_leaves_no_sequence_regions_that_stood_before_it => "regions-corrupt",
    sequence_regions_by_the_kind_of_slice_length_and_circular_flag => "regions-shapes",
    damaged_length_before_any_other_fails_its_shard => "regions-damaged-first",
    damaged_length_after_the_chromosome_states_one_is_not_read => "regions-damaged-later",
    // Every field of a transcript in every scalar kind a cache can store it as.
    transcript_fields_in_every_stored_kind => "scalars",
    native_order_doubles_and_integers => "native",
    roots_that_are_an_array_a_scalar_and_a_blessed_hash => "roots",
    shards_that_cannot_be_read_or_written_fail_the_run_by_name => "bad-shards",
    // Chromosome directories, the contig order, variation files and synonyms.
    every_directory_is_a_chromosome_in_contig_order => "contigs",
    third_argument_names_the_chromosomes_to_convert => "contigs-filter",
    chromosome_that_is_not_a_directory_is_named => "contigs-absent",
    empty_first_field_of_the_chromosome_list_is_a_name => "contigs-empty-field",
    third_argument_zero_converts_every_chromosome => "contigs-zero",
    all_vars_without_its_index_stops_the_run => "no-index",
    sereal_serialiser_line_is_refused => "sereal-info",
    sereal_shards_are_refused => "sereal-shards",
    cache_without_a_transcript_shard_is_refused => "no-shards",
    cache_without_info_txt_converts_with_a_warning => "no-info",
    info_txt_lines_as_ensembl_vep_reads_them => "info-lines",
    cache_directory_that_is_absent_is_named => "absent",
    cache_directory_that_is_a_file_is_named => "not-a-directory",
    // Regulatory shards that are refused, each by its message.
    regulatory_root_that_is_no_hash => "reg-root-array",
    regulatory_chromosome_entry_that_is_no_hash => "reg-chromosome-entry",
    regulatory_feature_list_that_is_no_array => "reg-feature-list",
    regulatory_feature_that_is_no_object => "reg-object",
    regulatory_feature_that_lacks_a_key => "reg-lacks-key",
    regulatory_feature_of_an_unknown_class => "reg-unknown-type",
    motif_that_lacks_its_strand => "reg-lacks-strand",
    motif_that_lacks_its_score => "reg-lacks-score",
    motif_that_lacks_its_sequence => "reg-lacks-sequence",
    motif_whose_sequence_does_not_fill_its_span => "reg-span",
    motif_whose_binding_matrix_is_no_object => "reg-invalid-matrix",
    binding_matrix_that_lacks_its_stable_id => "reg-matrix-id",
    binding_matrix_in_another_unit => "reg-unit",
    binding_matrix_whose_complexes_are_no_list => "reg-complexes",
    binding_matrix_whose_elements_are_no_hash => "reg-elements",
    binding_matrix_whose_length_is_not_its_positions => "reg-length",
    binding_matrix_that_lacks_a_position => "reg-position",
    binding_matrix_that_lacks_a_base => "reg-base",
    binding_matrix_with_neither_elements_nor_a_length => "reg-neither",
    regulatory_values_in_every_stored_kind => "regulatory-shapes",
    // The command line.
    no_argument_prints_the_usage_line => "usage-none",
    one_argument_prints_the_usage_line => "usage-one",
    argument_that_is_a_false_string_prints_the_usage_line => "usage-false",
}
