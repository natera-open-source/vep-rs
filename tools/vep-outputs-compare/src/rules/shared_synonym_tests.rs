// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The names a synonyms file gives to two or more held sequences, the sequences an Ensembl
//! VEP cache directory holds, and the class's pass over fabricated discordant rows with and
//! without vep-rs's output directory.

use super::*;
use crate::TempDir;

fn held(names: &[&str]) -> HashSet<String> {
    names.iter().map(|name| name.to_string()).collect()
}

fn names(listed: &[&str]) -> HashSet<String> {
    held(listed)
}

/// One discordant row of the output `variants.txt`.
fn row(source: &str, location: &str, feature: &str, consequence: &str) -> Row {
    let feature_type = if feature == "-" { "-" } else { "Transcript" };
    let cells = [
        "variants.txt",
        source,
        location,
        "A",
        feature,
        feature_type,
        consequence,
    ];
    Row {
        file_name: cells[0].to_string(),
        source: source.to_string(),
        location: location.to_string(),
        allele: cells[3].to_string(),
        feature: feature.to_string(),
        feature_type: feature_type.to_string(),
        consequence_set: consequence.to_string(),
        fields: cells.iter().map(|cell| cell.to_string()).collect(),
    }
}

fn scope_of(shared: Option<&[&str]>) -> ReferenceRecordScope {
    let mut scope = ReferenceRecordScope::default();
    if let Some(shared) = shared {
        scope.record_shared_synonyms(names(shared));
    }
    scope
}

/// The header of a reference output whose command line names no synonyms file.
const HEADER: &str = "## VEP command-line: vep --cache --offline --no_stats\n\
                      #Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";

/// A reference and a vep-rs output directory, each holding `variants.txt`: the reference's
/// is `reference_header` alone, vep-rs's the column line and `rust_rows`.
fn output_dirs(reference_header: &str, rust_rows: &str) -> TempDir {
    let dirs = TempDir::create(&std::env::temp_dir(), "shared_synonym_test_").unwrap();
    for side in ["perl", "rust"] {
        std::fs::create_dir(dirs.path.join(side)).unwrap();
    }
    std::fs::write(dirs.path.join("perl/variants.txt"), reference_header).unwrap();
    let columns = HEADER.lines().nth(1).unwrap();
    let rust = format!("{columns}\n{rust_rows}");
    std::fs::write(dirs.path.join("rust/variants.txt"), rust).unwrap();
    dirs
}

/// The statistics of a pass over `rows` beside the two output directories of `dirs`, and
/// the rows it set aside.
fn pass_beside(
    rows: &[Row],
    shared: Option<&[&str]>,
    dirs: &TempDir,
) -> (Option<OneSidedStats>, Vec<usize>) {
    let mut masked = HashSet::new();
    let outputs = (dirs.path.join("perl"), dirs.path.join("rust"));
    let scope = scope_of(shared);
    let stats =
        apply_shared_synonym_rule(rows, &mut masked, &scope, Some((&outputs.0, &outputs.1)))
            .unwrap();
    let mut set_aside: Vec<usize> = masked.into_iter().collect();
    set_aside.sort_unstable();
    (stats, set_aside)
}

/// [`pass_beside`] a reference output whose header names no synonyms file and a vep-rs
/// output with no row.
fn pass(rows: &[Row], shared: Option<&[&str]>) -> (Option<OneSidedStats>, Vec<usize>) {
    pass_beside(rows, shared, &output_dirs(HEADER, ""))
}

fn counts(stats: &OneSidedStats) -> [u64; 5] {
    [
        stats.count,
        stats.excluded,
        stats.records,
        stats.perl,
        stats.rust,
    ]
}

/// A name the cache holds no sequence by is of the class where two or more of the names the
/// file pairs it with are held, whichever side of a line it stands on, and each held name
/// counts once. A name with one held partner, a held name and a name whose partners the
/// cache lacks are not.
#[test]
fn a_name_is_shared_where_two_or_more_held_sequences_carry_it() {
    let text = "A\tshared\nshared\tB\nA\tshared\nA\tonly_a\nonly_a\tunheld\nA\tB\nC\tunheld\n\
                B\tunheld\nnowhere\talso_nowhere\n";
    assert_eq!(
        names_of_several_sequences(text, &held(&["A", "B", "C"])),
        names(&["shared", "unheld"])
    );
    assert_eq!(
        names_of_several_sequences(text, &held(&["A"])),
        HashSet::new()
    );
    assert_eq!(
        names_of_several_sequences("", &held(&["A"])),
        HashSet::new()
    );
}

/// The file is read as the engines read it (`BaseVEP.pm:646-656`): fields split on runs of
/// whitespace, a carriage return among them, every later field paired with the first and
/// not with each other, and a line that opens with whitespace pairing its fields with the
/// empty name.
#[test]
fn the_file_is_split_on_whitespace_and_pairs_each_field_with_the_first() {
    let sequences = held(&["A", "B"]);
    assert_eq!(
        names_of_several_sequences("A  shared\r\nB\x0Bshared\r\n", &sequences),
        names(&["shared"])
    );
    assert_eq!(
        names_of_several_sequences("first A B\n", &sequences),
        names(&["first"])
    );
    assert_eq!(
        names_of_several_sequences("A first second\nB first\n", &sequences),
        names(&["first"])
    );
    assert_eq!(
        names_of_several_sequences(" A B\n", &sequences),
        names(&[""])
    );
}

/// The held sequences are the cache directory's own subdirectories that open with no dot
/// (`CacheDir.pm:431`), whatever a converted cache lists: a file and a dotted directory are
/// no sequence. A cache without the synonyms file gives no name, and a directory that is
/// not there is an error that names it.
#[test]
fn the_reference_cache_gives_its_held_sequences_and_its_synonyms() {
    let cache = TempDir::create(&std::env::temp_dir(), "shared_synonym_test_").unwrap();
    for sequence in ["A", "B", ".hidden"] {
        std::fs::create_dir(cache.path.join(sequence)).unwrap();
    }
    std::fs::write(cache.path.join("C"), "a file").unwrap();
    assert_eq!(shared_synonym_names(&cache.path).unwrap(), HashSet::new());
    std::fs::write(
        cache.path.join(SYNONYMS_FILE),
        "A\tshared\nB\tshared\nA\tdotted\n.hidden\tdotted\nA\tfiled\nC\tfiled\nA\tB\n",
    )
    .unwrap();
    assert_eq!(
        shared_synonym_names(&cache.path).unwrap(),
        names(&["shared"])
    );
    let missing = cache.path.join("absent");
    let refusal = shared_synonym_names(&missing).unwrap_err().to_string();
    assert!(
        refusal.starts_with(&missing.display().to_string()),
        "{refusal}"
    );
}

/// The reference's rows of a record under a shared name, where vep-rs wrote none, are set
/// aside and counted on the reference's side, one record a Location.
#[test]
fn the_reference_rows_of_a_record_vep_rs_skipped_are_set_aside() {
    let rows = [
        row("missing_in_rust", "shared:399", "T1", "intron_variant"),
        row(
            "missing_in_rust",
            "shared:399",
            "T2",
            "upstream_gene_variant",
        ),
        row(
            "missing_in_rust",
            "shared:800-802",
            "-",
            "intergenic_variant",
        ),
    ];
    let (stats, set_aside) = pass(&rows, Some(&["shared"]));
    assert_eq!(counts(&stats.unwrap()), [3, 3, 2, 3, 0]);
    assert_eq!(set_aside, [0, 1, 2]);
}

/// A sequence name that holds a colon is read whole, up to the Location's last colon: the
/// rows under a shared name of that shape are set aside, and the rows of another name that
/// only opens with a shared one are not.
#[test]
fn a_name_that_holds_a_colon_is_read_up_to_the_last_colon() {
    let rows = [
        row("missing_in_rust", "scaffold:7:399", "T1", "intron_variant"),
        row(
            "missing_in_rust",
            "shared:other:399",
            "T2",
            "intron_variant",
        ),
    ];
    let (stats, set_aside) = pass(&rows, Some(&["scaffold:7", "shared"]));
    assert_eq!(counts(&stats.unwrap()), [1, 1, 1, 1, 0]);
    assert_eq!(set_aside, [0]);
}

/// A record vep-rs wrote a row for is one it did not skip: its rows stay, the reference's
/// among them. So do the rows of a name the file gives to one sequence, and vep-rs's rows of
/// a record the reference wrote nothing for.
#[test]
fn a_record_vep_rs_wrote_a_row_for_stays_a_difference() {
    let rows = [
        row("missing_in_rust", "shared:399", "-", "intergenic_variant"),
        row("extra_in_rust", "shared:399", "T1", "intron_variant"),
        row("missing_in_rust", "other:399", "T2", "intron_variant"),
        row("extra_in_rust", "shared:500", "T3", "intron_variant"),
    ];
    let (stats, set_aside) = pass(&rows, Some(&["shared"]));
    assert_eq!(counts(&stats.unwrap()), [0, 0, 0, 0, 0]);
    assert_eq!(set_aside, [] as [usize; 0]);
}

/// A collector without the names applies nothing and states no class; a row an earlier rule
/// set aside is not counted again.
#[test]
fn no_names_no_class_and_a_row_set_aside_before_is_not_counted() {
    let rows = [
        row("missing_in_rust", "shared:399", "T1", "intron_variant"),
        row("missing_in_rust", "shared:399", "T2", "intron_variant"),
    ];
    let (stats, set_aside) = pass(&rows, None);
    assert!(stats.is_none());
    assert_eq!(set_aside, [] as [usize; 0]);

    let mut masked = HashSet::from([0]);
    let scope = scope_of(Some(&["shared"]));
    let dirs = output_dirs(HEADER, "");
    let outputs = (dirs.path.join("perl"), dirs.path.join("rust"));
    let stats =
        apply_shared_synonym_rule(&rows, &mut masked, &scope, Some((&outputs.0, &outputs.1)))
            .unwrap();
    assert_eq!(counts(&stats.unwrap()), [1, 1, 1, 1, 0]);
    assert_eq!(masked, HashSet::from([0, 1]));
}

/// A row vep-rs wrote at the Location that the reference wrote too, and so is no discordant
/// row, keeps the reference's other rows open; a Location vep-rs's output holds no row at is
/// set aside.
#[test]
fn a_row_both_engines_wrote_at_the_location_keeps_the_rows_open() {
    let dirs = output_dirs(
        HEADER,
        "r1\tshared:399\tA\tG1\tT9\tTranscript\tintron_variant\n\
         r3\tother:5\tA\tG1\tT9\tTranscript\tintron_variant\n",
    );
    let rows = [
        row("missing_in_rust", "shared:399", "T1", "intron_variant"),
        row("missing_in_rust", "shared:700", "T2", "intron_variant"),
    ];
    let (stats, set_aside) = pass_beside(&rows, Some(&["shared"]), &dirs);
    assert_eq!(counts(&stats.unwrap()), [1, 1, 1, 1, 0]);
    assert_eq!(set_aside, [1]);
}

/// Only a reference run that read the cache's own synonyms file is of the class: one whose
/// header states a command line that names `--synonyms`, in either spelling, and one whose
/// header states no command line, keep their rows open. Without the two output directories
/// the rule applies nothing.
#[test]
fn a_reference_run_given_a_synonyms_file_or_stating_no_command_keeps_its_rows_open() {
    let rows = [row("missing_in_rust", "shared:399", "T1", "intron_variant")];
    let columns = HEADER.lines().nth(1).unwrap();
    for header in [
        format!("## VEP command-line: vep --cache --synonyms [PATH]/names.txt\n{columns}\n"),
        format!("## VEP command-line: vep --synonyms=names.txt --cache\n{columns}\n"),
        format!("{columns}\n"),
        String::new(),
    ] {
        let dirs = output_dirs(&header, "");
        let (stats, set_aside) = pass_beside(&rows, Some(&["shared"]), &dirs);
        assert_eq!(counts(&stats.unwrap()), [0, 0, 0, 0, 0], "{header:?}");
        assert_eq!(set_aside, [] as [usize; 0], "{header:?}");
    }
    let named_elsewhere = "## VEP command-line: vep --cache --var_synonyms\n";
    let dirs = output_dirs(&format!("{named_elsewhere}{columns}\n"), "");
    let (stats, set_aside) = pass_beside(&rows, Some(&["shared"]), &dirs);
    assert_eq!(counts(&stats.unwrap()), [1, 1, 1, 1, 0]);
    assert_eq!(set_aside, [0]);

    let mut masked = HashSet::new();
    let scope = scope_of(Some(&["shared"]));
    let stats = apply_shared_synonym_rule(&rows, &mut masked, &scope, None).unwrap();
    assert_eq!(counts(&stats.unwrap()), [0, 0, 0, 0, 0]);
    assert!(masked.is_empty());
}
