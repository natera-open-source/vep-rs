// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `check`, `row`, `keys` and `put` over table directories built for each test.

mod common;

use std::fs;
use std::path::Path;

use common::{
    bacteria_rows, digest_of, failure, genome_rows, gunzip, gzip, outcome_rows, scratch, set_rows,
    table, text, tool, valid_table, ARCHAEON_KEY, BACTERIUM_KEY, CEREVISIAE_KEY, HEADER, HUMAN_KEY,
    HUMAN_MERGED_KEY, POMBE_KEY, UNCLASSIFIED,
};
use vep_genome_table::{genome_key, put};

fn dir_arg(dir: &Path) -> &str {
    dir.to_str().unwrap()
}

/// `check` of `dir` fails with exactly `expected` on standard error.
fn refused(dir: &Path, flags: &[&str], expected: &str) {
    let mut args = vec!["check", dir_arg(dir)];
    args.extend_from_slice(flags);
    let run = tool(&args);
    assert_eq!(run.stderr, expected);
    assert_eq!(run.stdout, "");
    assert_eq!(run.code, 1);
}

#[test]
fn a_valid_table_of_six_rows_across_three_sets_passes() {
    let dir = valid_table("valid");
    let run = tool(&["check", dir_arg(&dir)]);
    assert_eq!(run.stdout, "6 rows in 3 sets: ok\n");
    assert_eq!(run.stderr, "");
    assert_eq!(run.code, 0);
}

#[test]
fn a_table_of_headers_alone_passes_with_no_rows() {
    let dir = table("headers_alone", &[], &[], &[UNCLASSIFIED.to_owned()], &[]);
    for flags in [&[][..], &["--complete"][..]] {
        let mut args = vec!["check", dir_arg(&dir)];
        args.extend_from_slice(flags);
        let run = tool(&args);
        assert_eq!(run.stdout, "0 rows in 0 sets: ok\n");
        assert_eq!(run.code, 0);
    }
}

#[test]
fn a_row_with_a_field_too_many_is_refused() {
    let mut genomes = genome_rows();
    genomes[2].push_str("\textra");
    let dir = table(
        "field_too_many",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "genomes.tsv",
            4,
            "10 fields where the table has 9 columns",
        ),
    );
}

#[test]
fn a_row_with_a_field_too_few_is_refused() {
    let mut bacteria = bacteria_rows();
    let short = bacteria[0].rsplit_once('\t').unwrap().0.to_owned();
    bacteria[0] = short;
    let dir = table(
        "field_too_few",
        &genome_rows(),
        &bacteria,
        &outcome_rows(),
        &set_rows(),
    );
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "bacteria.tsv.gz",
            2,
            "8 fields where the table has 9 columns",
        ),
    );
}

#[test]
fn a_bad_transcripts_value_is_refused() {
    let mut genomes = genome_rows();
    genomes[3] = genomes[3].replace("\tmerged\t", "\tboth\t");
    let dir = table(
        "bad_transcripts",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "genomes.tsv",
            5,
            "transcripts \"both\" is not ensembl, refseq or merged",
        ),
    );
}

#[test]
fn an_undeclared_outcome_is_refused() {
    let dir = table(
        "undeclared_outcome",
        &genome_rows(),
        &bacteria_rows(),
        &[UNCLASSIFIED.to_owned()],
        &set_rows(),
    );
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "bacteria.tsv.gz",
            3,
            "outcome excluded:no-transcripts is not declared in outcomes.tsv",
        ),
    );
}

#[test]
fn a_key_held_twice_is_refused() {
    let mut genomes = genome_rows();
    let again = genomes[2].replace("\t31\t", "\t32\t");
    genomes.insert(3, again);
    let dir = table(
        "key_twice",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    let first = dir.join("genomes.tsv");
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "genomes.tsv",
            5,
            &format!("the key {HUMAN_KEY} is also at {}:4", first.display()),
        ),
    );
}

#[test]
fn a_row_out_of_key_order_is_refused() {
    let mut genomes = genome_rows();
    genomes.swap(0, 1);
    let dir = table(
        "out_of_order",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "genomes.tsv",
            3,
            &format!("out of key order: {CEREVISIAE_KEY} follows {POMBE_KEY}"),
        ),
    );
}

#[test]
fn a_bacterial_row_in_the_plain_file_is_refused() {
    let mut bacteria = bacteria_rows();
    let mut genomes = genome_rows();
    genomes.insert(0, bacteria.remove(0));
    let dir = table(
        "bacterial_row_in_plain_file",
        &genomes,
        &bacteria,
        &outcome_rows(),
        &set_rows(),
    );
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "genomes.tsv",
            2,
            "a row of the set bacteria_0_collection belongs in bacteria.tsv.gz",
        ),
    );
}

#[test]
fn a_plain_row_in_the_bacterial_file_is_refused() {
    let mut bacteria = bacteria_rows();
    let mut genomes = genome_rows();
    bacteria.push(genomes.remove(0));
    let dir = table(
        "plain_row_in_bacterial_file",
        &genomes,
        &bacteria,
        &outcome_rows(),
        &set_rows(),
    );
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "bacteria.tsv.gz",
            4,
            "a row of the set fungi belongs in genomes.tsv",
        ),
    );
}

#[test]
fn an_exact_row_without_a_digest_is_refused() {
    let mut genomes = genome_rows();
    genomes[2] = genomes[2].replace(&digest_of('4'), "-");
    let dir = table(
        "exact_without_digest",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "genomes.tsv",
            4,
            "outcome exact needs a sha256, not -",
        ),
    );
}

#[test]
fn a_differs_row_may_leave_out_its_digest() {
    let mut genomes = genome_rows();
    genomes[1] = genomes[1].replace(&digest_of('3'), "-");
    let dir = table(
        "differs_without_digest",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    let run = tool(&["check", dir_arg(&dir)]);
    assert_eq!(run.stdout, "6 rows in 3 sets: ok\n");
    assert_eq!(run.code, 0);
}

#[test]
fn an_exact_row_without_records_is_refused() {
    let mut genomes = genome_rows();
    genomes[0] = genomes[0].replace("\t30\t", "\t0\t");
    let dir = table(
        "exact_without_records",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "genomes.tsv",
            2,
            "outcome exact needs records above zero",
        ),
    );
}

#[test]
fn every_failure_is_one_line_in_file_then_line_order() {
    let mut genomes = genome_rows();
    genomes[3] = genomes[3].replace("\tmerged\t", "\tboth\t");
    genomes[0] = genomes[0].replace("\t63\t", "\t063\t");
    let dir = table(
        "two_files_fail",
        &genomes,
        &bacteria_rows(),
        &[UNCLASSIFIED.to_owned()],
        &set_rows(),
    );
    let expected = [
        failure(
            &dir,
            "genomes.tsv",
            2,
            "cache_version \"063\" is not an integer above zero",
        ),
        failure(
            &dir,
            "genomes.tsv",
            5,
            "transcripts \"both\" is not ensembl, refseq or merged",
        ),
        failure(
            &dir,
            "bacteria.tsv.gz",
            3,
            "outcome excluded:no-transcripts is not declared in outcomes.tsv",
        ),
    ]
    .concat();
    refused(&dir, &[], &expected);
}

#[test]
fn a_file_that_opens_with_another_header_is_refused() {
    let dir = valid_table("other_header");
    let renamed = text(&HEADER.replace("sha256", "checksum"), &genome_rows());
    fs::write(dir.join("genomes.tsv"), renamed).unwrap();
    refused(
        &dir,
        &[],
        &failure(
            &dir,
            "genomes.tsv",
            1,
            "the first line is not the tab-separated header genome assembly transcripts set \
             cache_version seed records outcome sha256",
        ),
    );
}

#[test]
fn a_missing_file_is_named_at_line_zero() {
    let dir = valid_table("missing_file");
    fs::remove_file(dir.join("sets.tsv")).unwrap();
    let run = tool(&["check", dir_arg(&dir)]);
    let opening = format!(
        "ERROR: [vep-genome-table] {}:0: cannot be read: ",
        dir.join("sets.tsv").display()
    );
    assert!(run.stderr.starts_with(&opening), "{}", run.stderr);
    assert_eq!(run.stderr.lines().count(), 1);
    assert_eq!(run.code, 1);
}

#[test]
fn the_outcomes_and_sets_files_are_held_to_their_own_grammar() {
    let dir = table(
        "side_files",
        &genome_rows(),
        &bacteria_rows(),
        &[
            UNCLASSIFIED.to_owned(),
            "exact\tThe outputs are equal.".to_owned(),
            "excluded:no-transcripts\t".to_owned(),
            UNCLASSIFIED.to_owned(),
        ],
        &[
            "bacteria_0_collection\t63\t63\t2".to_owned(),
            "fungi\t63\t63".to_owned(),
            "vertebrates\t116\t116\ttwo".to_owned(),
            "bacteria_0_collection\t63\t63\t2".to_owned(),
        ],
    );
    let expected = [
        failure(
            &dir,
            "bacteria.tsv.gz",
            3,
            "outcome excluded:no-transcripts is not declared in outcomes.tsv",
        ),
        failure(
            &dir,
            "outcomes.tsv",
            3,
            "outcome \"exact\" is not differs:<name> or excluded:<name>",
        ),
        failure(
            &dir,
            "outcomes.tsv",
            4,
            "outcome excluded:no-transcripts has no meaning",
        ),
        failure(
            &dir,
            "outcomes.tsv",
            5,
            "outcome differs:unclassified is declared twice",
        ),
        failure(&dir, "sets.tsv", 3, "3 fields where the file has 4 columns"),
        failure(
            &dir,
            "sets.tsv",
            4,
            "caches \"two\" is not an unsigned integer",
        ),
        failure(
            &dir,
            "sets.tsv",
            5,
            "set bacteria_0_collection is listed twice",
        ),
    ]
    .concat();
    refused(&dir, &[], &expected);
}

#[test]
fn complete_passes_when_every_set_has_its_caches() {
    let dir = valid_table("complete");
    let run = tool(&["check", dir_arg(&dir), "--complete"]);
    assert_eq!(run.stdout, "6 rows in 3 sets: ok\n");
    assert_eq!(run.code, 0);
}

#[test]
fn complete_fails_on_a_missing_row() {
    let mut genomes = genome_rows();
    genomes.pop();
    let dir = table(
        "complete_missing_row",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    let plain = tool(&["check", dir_arg(&dir)]);
    assert_eq!(plain.stdout, "5 rows in 3 sets: ok\n");
    refused(
        &dir,
        &["--complete"],
        &failure(
            &dir,
            "sets.tsv",
            4,
            "set vertebrates has 1 rows where caches is 2",
        ),
    );
}

#[test]
fn complete_fails_on_a_set_the_sets_file_lacks() {
    let mut sets = set_rows();
    sets.remove(1);
    let dir = table(
        "complete_unlisted_set",
        &genome_rows(),
        &bacteria_rows(),
        &outcome_rows(),
        &sets,
    );
    let plain = tool(&["check", dir_arg(&dir)]);
    assert_eq!(plain.stdout, "6 rows in 3 sets: ok\n");
    refused(
        &dir,
        &["--complete"],
        &failure(
            &dir,
            "genomes.tsv",
            2,
            "set fungi is not in sets.tsv (2 rows name it)",
        ),
    );
}

#[test]
fn put_of_the_same_rows_twice_leaves_both_files_byte_identical() {
    let dir = table("put_twice", &[], &[], &outcome_rows(), &set_rows());
    let mut shuffled = genome_rows();
    shuffled.extend(bacteria_rows());
    shuffled.reverse();
    let rows = scratch("put_twice_rows").join("rows.tsv");
    fs::write(&rows, text(HEADER, &shuffled)).unwrap();

    let first = tool(&["put", dir_arg(&dir), rows.to_str().unwrap()]);
    assert_eq!(first.stdout, "6 rows in 3 sets: 6 added, 0 replaced\n");
    assert_eq!(first.code, 0);
    let plain = fs::read(dir.join("genomes.tsv")).unwrap();
    let compressed = fs::read(dir.join("bacteria.tsv.gz")).unwrap();
    assert_eq!(
        String::from_utf8(plain.clone()).unwrap(),
        text(HEADER, &genome_rows())
    );
    assert_eq!(gunzip(&compressed), text(HEADER, &bacteria_rows()));
    // RFC 1952: byte 3 is FLG, whose FNAME bit is 0x08, and bytes 4 to 7 are MTIME.
    assert_eq!(compressed[3] & 0x08, 0);
    assert_eq!(compressed[4..8], [0, 0, 0, 0]);

    let second = tool(&["put", dir_arg(&dir), rows.to_str().unwrap()]);
    assert_eq!(second.stdout, "6 rows in 3 sets: 0 added, 6 replaced\n");
    assert_eq!(fs::read(dir.join("genomes.tsv")).unwrap(), plain);
    assert_eq!(fs::read(dir.join("bacteria.tsv.gz")).unwrap(), compressed);
    assert_eq!(
        tool(&["check", dir_arg(&dir), "--complete"]).stdout,
        "6 rows in 3 sets: ok\n"
    );
}

#[test]
fn put_of_a_changed_row_replaces_it() {
    let dir = valid_table("put_changed");
    let changed = genome_rows()[1]
        .replace("differs:unclassified", "exact")
        .replace(&digest_of('3'), &digest_of('9'));
    let rows = scratch("put_changed_rows").join("row.tsv");
    fs::write(&rows, text(HEADER, std::slice::from_ref(&changed))).unwrap();

    let run = tool(&["put", dir_arg(&dir), rows.to_str().unwrap()]);
    assert_eq!(run.stdout, "6 rows in 3 sets: 0 added, 1 replaced\n");
    let mut expected = genome_rows();
    expected[1] = changed;
    assert_eq!(
        fs::read_to_string(dir.join("genomes.tsv")).unwrap(),
        text(HEADER, &expected)
    );
    let outcome = tool(&["row", dir_arg(&dir), POMBE_KEY, "--field", "outcome"]);
    assert_eq!(outcome.stdout, "exact\n");
}

#[test]
fn put_takes_a_gzip_row_file_and_the_last_row_given_for_a_key() {
    let dir = valid_table("put_gzip");
    let earlier = genome_rows()[2].replace("\t31\t", "\t40\t");
    let later = genome_rows()[2].replace("\t31\t", "\t41\t");
    let rows_dir = scratch("put_gzip_rows");
    let first = rows_dir.join("first.tsv.gz");
    let second = rows_dir.join("second.tsv");
    fs::write(&first, gzip(text(HEADER, &[earlier]).as_bytes())).unwrap();
    fs::write(&second, text(HEADER, &[later])).unwrap();

    let run = tool(&[
        "put",
        dir_arg(&dir),
        first.to_str().unwrap(),
        second.to_str().unwrap(),
    ]);
    assert_eq!(run.stdout, "6 rows in 3 sets: 0 added, 1 replaced\n");
    let records = tool(&["row", dir_arg(&dir), HUMAN_KEY, "--field", "records"]);
    assert_eq!(records.stdout, "41\n");
}

#[test]
fn put_refuses_a_row_that_fails_the_field_rules_and_writes_nothing() {
    let dir = valid_table("put_refused");
    let plain = fs::read(dir.join("genomes.tsv")).unwrap();
    let compressed = fs::read(dir.join("bacteria.tsv.gz")).unwrap();
    let good = format!(
        "mus_musculus\tGRCm39\tensembl\tvertebrates\t116\t11\t30\texact\t{}",
        digest_of('6')
    );
    let bad = good.replace("\tensembl\t", "\tgencode\t");
    let rows = scratch("put_refused_rows").join("rows.tsv");
    fs::write(&rows, text(HEADER, &[good, bad])).unwrap();

    let run = tool(&["put", dir_arg(&dir), rows.to_str().unwrap()]);
    assert_eq!(
        run.stderr,
        format!(
            "ERROR: [vep-genome-table] {}:3: transcripts \"gencode\" is not ensembl, refseq or merged\n",
            rows.display()
        )
    );
    assert_eq!(run.stdout, "");
    assert_eq!(run.code, 1);
    assert_eq!(fs::read(dir.join("genomes.tsv")).unwrap(), plain);
    assert_eq!(fs::read(dir.join("bacteria.tsv.gz")).unwrap(), compressed);
    let names: Vec<String> = fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    assert_eq!(names.len(), 4, "{names:?}");
}

#[test]
fn put_that_cannot_write_one_file_leaves_the_table_and_no_partial_file() {
    let dir = valid_table("put_unwritable");
    let plain = fs::read(dir.join("genomes.tsv")).unwrap();
    let compressed = fs::read(dir.join("bacteria.tsv.gz")).unwrap();
    // `put` writes each file as `<file>.<process id>.tmp` before renaming it, and this test
    // calls it in its own process: a directory of that name stops the second file's write.
    let blocked = format!("bacteria.tsv.gz.{}.tmp", std::process::id());
    fs::create_dir(dir.join(&blocked)).unwrap();
    let changed = genome_rows()[0].replace("\t30\t", "\t33\t");
    let rows = scratch("put_unwritable_rows").join("row.tsv");
    fs::write(&rows, text(HEADER, &[changed])).unwrap();

    let failures = put(&dir, &[rows]).unwrap_err();
    assert_eq!(failures.len(), 1, "{failures:?}");
    assert_eq!(failures[0].file, dir.join("bacteria.tsv.gz"));
    assert_eq!(failures[0].line, 0);
    assert!(
        failures[0].what.starts_with("cannot be written: "),
        "{}",
        failures[0].what
    );
    assert_eq!(fs::read(dir.join("genomes.tsv")).unwrap(), plain);
    assert_eq!(fs::read(dir.join("bacteria.tsv.gz")).unwrap(), compressed);
    let mut names: Vec<String> = fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    assert_eq!(
        names,
        [
            "bacteria.tsv.gz",
            blocked.as_str(),
            "genomes.tsv",
            "outcomes.tsv",
            "sets.tsv"
        ]
    );
}

#[test]
fn a_key_of_a_genome_that_begins_with_an_underscore_round_trips_through_row() {
    let dir = valid_table("underscore_genome");
    let key = genome_key(
        "bacteria_0_collection",
        "_example_bacterium_gca_000000001",
        "ASM1v1",
        "ensembl",
    );
    assert_eq!(key, BACTERIUM_KEY);
    let row = tool(&["row", dir_arg(&dir), &key]);
    assert_eq!(row.stdout, format!("{HEADER}\n{}\n", bacteria_rows()[0]));
    assert_eq!(row.code, 0);
    let genome = tool(&["row", dir_arg(&dir), &key, "--field=genome"]);
    assert_eq!(genome.stdout, "_example_bacterium_gca_000000001\n");
}

#[test]
fn row_of_a_key_no_row_has_exits_1_with_one_error_line() {
    let dir = valid_table("row_unknown_key");
    let key = "vertebrates__mus_musculus__GRCm39__ensembl";
    let run = tool(&["row", dir_arg(&dir), key]);
    assert_eq!(
        run.stderr,
        format!(
            "ERROR: [vep-genome-table] {}: no row has the key {key}\n",
            dir.display()
        )
    );
    assert_eq!(run.stdout, "");
    assert_eq!(run.code, 1);
}

#[test]
fn keys_lists_every_key_in_order_and_filters_by_set_and_outcome() {
    let dir = valid_table("keys");
    let all = tool(&["keys", dir_arg(&dir)]);
    assert_eq!(
        all.stdout,
        [
            BACTERIUM_KEY,
            ARCHAEON_KEY,
            CEREVISIAE_KEY,
            POMBE_KEY,
            HUMAN_KEY,
            HUMAN_MERGED_KEY
        ]
        .map(|key| format!("{key}\n"))
        .concat()
    );
    let fungi = tool(&["keys", dir_arg(&dir), "--set", "fungi"]);
    assert_eq!(fungi.stdout, format!("{CEREVISIAE_KEY}\n{POMBE_KEY}\n"));
    let differs = tool(&["keys", dir_arg(&dir), "--outcome", "differs"]);
    assert_eq!(differs.stdout, format!("{POMBE_KEY}\n"));
    let excluded = tool(&["keys", dir_arg(&dir), "--outcome=excluded:"]);
    assert_eq!(excluded.stdout, format!("{ARCHAEON_KEY}\n"));
    let exact = tool(&[
        "keys",
        dir_arg(&dir),
        "--set",
        "vertebrates",
        "--outcome",
        "exact",
    ]);
    assert_eq!(exact.stdout, format!("{HUMAN_KEY}\n{HUMAN_MERGED_KEY}\n"));
}

#[test]
fn row_and_keys_refuse_a_table_that_fails_its_row_rules() {
    let mut genomes = genome_rows();
    genomes.swap(0, 1);
    let dir = table(
        "query_refused",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    let expected = failure(
        &dir,
        "genomes.tsv",
        3,
        &format!("out of key order: {CEREVISIAE_KEY} follows {POMBE_KEY}"),
    );
    let keys = tool(&["keys", dir_arg(&dir)]);
    assert_eq!(keys.stderr, expected);
    assert_eq!(keys.code, 1);
    let row = tool(&["row", dir_arg(&dir), HUMAN_KEY]);
    assert_eq!(row.stderr, expected);
    assert_eq!(row.code, 1);
}

#[test]
fn row_keys_and_put_print_the_failures_of_a_file_in_line_order() {
    let mut genomes = genome_rows();
    genomes[1] = genomes[0].replace("\t30\t", "\t33\t");
    genomes[3] = genomes[3].replace("\tmerged\t", "\tboth\t");
    let dir = table(
        "query_two_failures",
        &genomes,
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    );
    let first = dir.join("genomes.tsv");
    let expected = [
        failure(
            &dir,
            "genomes.tsv",
            3,
            &format!("the key {CEREVISIAE_KEY} is also at {}:2", first.display()),
        ),
        failure(
            &dir,
            "genomes.tsv",
            5,
            "transcripts \"both\" is not ensembl, refseq or merged",
        ),
    ]
    .concat();
    let rows = scratch("query_two_failures_rows").join("row.tsv");
    fs::write(&rows, text(HEADER, &bacteria_rows()[..1])).unwrap();
    for args in [
        &["keys", dir_arg(&dir)][..],
        &["row", dir_arg(&dir), HUMAN_KEY][..],
        &["put", dir_arg(&dir), rows.to_str().unwrap()][..],
    ] {
        let run = tool(args);
        assert_eq!(run.stderr, expected, "{args:?}");
        assert_eq!(run.stdout, "", "{args:?}");
        assert_eq!(run.code, 1, "{args:?}");
    }
}

#[test]
fn a_command_line_the_tool_cannot_read_exits_2() {
    let dir = valid_table("usage");
    for args in [
        &[][..],
        &["rows", dir_arg(&dir)][..],
        &["check"][..],
        &["check", dir_arg(&dir), "--all"][..],
        &["row", dir_arg(&dir), HUMAN_KEY, "--field", "key"][..],
        &["row", dir_arg(&dir), HUMAN_KEY, "--field"][..],
        &["keys", dir_arg(&dir), "--set", "fungi", "--set", "plants"][..],
        &["put", dir_arg(&dir)][..],
        &["digest"][..],
    ] {
        let run = tool(args);
        assert_eq!(run.code, 2, "{args:?}");
        assert!(
            run.stderr.starts_with("ERROR: [vep-genome-table] "),
            "{args:?}: {}",
            run.stderr
        );
        assert_eq!(run.stdout, "", "{args:?}");
    }
}

#[cfg(unix)]
#[test]
fn an_argument_that_is_not_utf8_exits_2() {
    use std::ffi::OsStr;
    use std::os::unix::ffi::OsStrExt;
    use std::process::Command;

    let output = Command::new(env!("CARGO_BIN_EXE_vep-genome-table"))
        .arg("digest")
        .arg(OsStr::from_bytes(b"output\xff.txt"))
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.starts_with("ERROR: [vep-genome-table] the argument output"),
        "{stderr}"
    );
    assert!(output.stdout.is_empty());
    assert_eq!(output.status.code(), Some(2));
}

/// `put` and `check` take an exact row of the collection `set`, which goes to the plain
/// file, and `check --complete` counts it once `sets.tsv` lists the set.
fn a_collection_row_is_taken(name: &str, set: &str, genome: &str) {
    let mut sets = set_rows();
    sets.push(format!("{set}\t63\t63\t1"));
    let dir = table(
        name,
        &genome_rows(),
        &bacteria_rows(),
        &outcome_rows(),
        &sets,
    );
    let row = format!(
        "{genome}\tASM1v1\tensembl\t{set}\t63\t11\t27\texact\t{}",
        digest_of('7')
    );
    let rows = scratch(&format!("{name}_rows")).join("row.tsv");
    fs::write(&rows, text(HEADER, std::slice::from_ref(&row))).unwrap();

    let run = tool(&["put", dir_arg(&dir), rows.to_str().unwrap()]);
    assert_eq!(run.stderr, "");
    assert_eq!(run.stdout, "7 rows in 4 sets: 1 added, 0 replaced\n");
    assert_eq!(run.code, 0);
    let plain = fs::read_to_string(dir.join("genomes.tsv")).unwrap();
    assert!(plain.contains(&format!("{row}\n")), "{plain}");
    let bacteria = fs::read(dir.join("bacteria.tsv.gz")).unwrap();
    assert_eq!(gunzip(&bacteria), text(HEADER, &bacteria_rows()));

    let check = tool(&["check", dir_arg(&dir), "--complete"]);
    assert_eq!(check.stderr, "");
    assert_eq!(check.stdout, "7 rows in 4 sets: ok\n");
    let key = genome_key(set, genome, "ASM1v1", "ensembl");
    let keys = tool(&["keys", dir_arg(&dir), "--set", set]);
    assert_eq!(keys.stdout, format!("{key}\n"));
}

#[test]
fn check_and_put_take_a_row_of_a_fungal_collection() {
    a_collection_row_is_taken(
        "fungal_collection",
        "fungi_ascomycota1_collection",
        "example_fungus_gca_000000003",
    );
}

#[test]
fn check_and_put_take_a_row_of_a_protist_collection() {
    a_collection_row_is_taken(
        "protist_collection",
        "protists_alveolata1_collection",
        "example_protist_gca_000000004",
    );
}

#[test]
fn check_and_put_refuse_a_set_of_no_family() {
    const FORMS: &str = "vertebrates, metazoa, plants, fungi, protists, pangenome, \
                         bacteria_<N>_collection, fungi_<name>_collection or \
                         protists_<name>_collection";
    for (index, set) in [
        "archaea_1_collection",
        "fungi_Ascomycota1_collection",
        "protists__collection",
    ]
    .into_iter()
    .enumerate()
    {
        let row = format!(
            "example_genome_gca_000000005\tASM1v1\tensembl\t{set}\t63\t11\t27\texact\t{}",
            digest_of('8')
        );
        let what = format!("set \"{set}\" is not {FORMS}");

        let dir = valid_table(&format!("no_family_put_{index}"));
        let plain = fs::read(dir.join("genomes.tsv")).unwrap();
        let rows = scratch(&format!("no_family_rows_{index}")).join("row.tsv");
        fs::write(&rows, text(HEADER, std::slice::from_ref(&row))).unwrap();
        let run = tool(&["put", dir_arg(&dir), rows.to_str().unwrap()]);
        assert_eq!(
            run.stderr,
            format!("ERROR: [vep-genome-table] {}:2: {what}\n", rows.display())
        );
        assert_eq!(run.code, 1);
        assert_eq!(fs::read(dir.join("genomes.tsv")).unwrap(), plain);

        let mut genomes = genome_rows();
        genomes.insert(0, row);
        let dir = table(
            &format!("no_family_check_{index}"),
            &genomes,
            &bacteria_rows(),
            &outcome_rows(),
            &set_rows(),
        );
        refused(&dir, &[], &failure(&dir, "genomes.tsv", 2, &what));
    }
}
