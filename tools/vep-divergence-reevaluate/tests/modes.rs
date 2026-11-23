// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The modes that run nothing: `--list`, `--check`, the help, and the command lines the
//! binary refuses.

mod common;

use common::{Scenario, FIXTURES};

const HEADER: &str = "| Entry | In one sentence | Tier and count | Verified against | Page |\n| --- | --- | --- | --- | --- |\n";
const MALFORMED_SECTION: &str = "Differences only malformed input reaches";
const UNOBSERVED_BASIS: &str = "well-formed under VCF 4.2; not seen in the measured datasets";
const SIGHTING_BASIS: &str =
    "seen in a structural-variant file of a population dataset: 4,158 pairs";

fn fixture_records() -> String {
    format!("{FIXTURES}/reproductions")
}

#[test]
fn check_holds_the_page_and_the_reproductions_to_one_reference() {
    let scenario = Scenario::new("check");
    let body: String = [("still_present_fixture", "116.3"), ("closed_fixture", "116.3")]
        .iter()
        .map(|(key, verified)| {
            format!("| [`{key}`](x.md) | s | 3 · none on the cuts · 1 flag set, 0 corpora | {verified} | [P](x.md) |\n")
        })
        .collect();
    let page = scenario.write(
        "page.md",
        format!("# Index\n\n## Section\n\n{HEADER}{body}\n## How to read the counts\n\nprose\n")
            .as_bytes(),
    );
    let records = fixture_records();
    let ran = scenario.run(&["--check", "--page", &page, "--reproductions", &records]);
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert!(ran
        .stdout
        .contains("2 entries have a reproduction verified against 116.3"));
    scenario.golden("check", &ran);
    let ran = scenario.run(&[
        "--check",
        "--page",
        &page,
        "--reproductions",
        &records,
        "--reference",
        "116.2",
    ]);
    assert_eq!(ran.code, 1);
    assert!(ran.stderr.contains("ERROR: [reevaluate]"));
    scenario.golden("check_other_reference", &ran);
    let missing = scenario.path("missing.md");
    let ran = scenario.run(&["--check", "--page", &missing, "--reproductions", &records]);
    assert_eq!(ran.code, 1);
    scenario.golden("check_missing_page", &ran);
    // A page that cannot be read as an index stops the check with the line it names.
    let table_first = scenario.write("table-first.md", format!("{HEADER}\n## Late\n").as_bytes());
    let ran = scenario.run(&[
        "--check",
        "--page",
        &table_first,
        "--reproductions",
        &records,
    ]);
    assert_eq!(ran.code, 1);
    scenario.golden("check_unreadable_page", &ran);
    let not_utf8 = scenario.write("latin.md", b"## Caf\xe9\n");
    let ran = scenario.run(&["--check", "--page", &not_utf8, "--reproductions", &records]);
    assert_eq!(ran.code, 1);
    scenario.golden("check_page_not_utf8", &ran);
}

/// `--check` passes cells that count cache transcripts against one total, whatever total a
/// cell of another unit names, and refuses two totals in one line that names each with its
/// rows and its first entry, in page order.
#[test]
fn check_refuses_an_index_whose_cache_transcripts_cells_name_two_totals() {
    let scenario = Scenario::new("totals");
    let page = |rows: &[(&str, &str)]| {
        let body: String = rows
            .iter()
            .map(|(key, figure)| {
                format!("| [`{key}`](x.md) | s | 3 · {figure} · 1 flag set, 0 corpora | 116.3 | [P](x.md) |\n")
            })
            .collect();
        scenario.write(
            "page.md",
            format!(
                "# Index\n\n## Section\n\n{HEADER}{body}\n## How to read the counts\n\nprose\n"
            )
            .as_bytes(),
        )
    };
    let records = fixture_records();
    let path = page(&[
        ("still_present_fixture", "0 of 1,200 cache transcripts"),
        ("closed_fixture", "0 of 1,200 cache transcripts"),
    ]);
    let check_args = ["--check", "--page", &path, "--reproductions", &records];
    let ran = scenario.run(&check_args);
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert!(ran
        .stdout
        .contains("2 entries have a reproduction verified against 116.3"));
    scenario.golden("totals_one", &ran);
    page(&[
        ("still_present_fixture", "0 of 1,200 cache transcripts"),
        ("closed_fixture", "0 of 90 regulatory features"),
    ]);
    let ran = scenario.run(&check_args);
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    scenario.golden("totals_other_unit", &ran);

    page(&[
        ("still_present_fixture", "0 of 1,200 cache transcripts"),
        ("closed_fixture", "0 of 1,250 cache transcripts"),
    ]);
    let ran = scenario.run(&check_args);
    assert_eq!(ran.code, 1);
    assert_eq!(
        ran.stderr,
        "ERROR: [reevaluate] the page's `Tier and count` cells name 2 totals of cache transcripts: 1,200 in 1 row, the first still_present_fixture; 1,250 in 1 row, the first closed_fixture\n"
    );
    scenario.golden("totals_two", &ran);

    let registry = scenario.mkdir("registry");
    for name in ["first_entry", "second_entry", "third_entry"] {
        write_record(&scenario, name, Some("unobserved"), Some(UNOBSERVED_BASIS));
    }
    page(&[
        ("first_entry", "0 of 1,200 cache transcripts"),
        ("second_entry", "0 of 1,250 cache transcripts"),
        ("third_entry", "0 of 1,200 cache transcripts"),
    ]);
    let ran = scenario.run(&["--check", "--page", &path, "--reproductions", &registry]);
    assert_eq!(ran.code, 1);
    scenario.golden("totals_two_of_three", &ran);
}

/// A copy of the still-present fixture under another name with the given input class and
/// basis; None leaves the field out.
fn write_record(scenario: &Scenario, name: &str, input_class: Option<&str>, basis: Option<&str>) {
    let fixture =
        std::fs::read_to_string(format!("{}/still_present_fixture.toml", fixture_records()))
            .expect("the fixture record");
    let fields = [
        ("entry", Some(name)),
        ("input_class", input_class),
        ("input_class_basis", basis),
    ];
    let mut lines = Vec::new();
    for line in fixture.lines() {
        match fields
            .iter()
            .find(|(key, _)| line.starts_with(&format!("{key} = ")))
        {
            None => lines.push(line.to_string()),
            Some((key, Some(value))) => lines.push(format!("{key} = \"{value}\"")),
            Some((_, None)) => {}
        }
    }
    scenario.write(
        &format!("registry/{name}.toml"),
        format!("{}\n", lines.join("\n")).as_bytes(),
    );
}

/// `--check` refuses a reproduction without its input class or its basis and one naming a
/// fourth class, and holds each class to the page.
#[test]
fn check_holds_each_input_class_to_the_pages_tier_and_section() {
    let scenario = Scenario::new("classes");
    let seen = "2 · 98.6 per million · 15 flag sets, 20 corpora";
    let unseen = "3 · none on the cuts · 1 flag set, 0 corpora";
    let row =
        |key: &str, cell: &str| format!("| [`{key}`](x.md) | s | {cell} | 116.3 | [P](x.md) |\n");
    let page = scenario.write(
        "page.md",
        format!(
            "# Index\n\n## Section\n\n{HEADER}{}{}\n## {MALFORMED_SECTION}\n\n### Section\n\n{HEADER}{}\n## How to read the counts\n\nprose\n",
            row("seen_entry", seen),
            row("unseen_entry", unseen),
            row("broken_entry", unseen)
        )
        .as_bytes(),
    );
    let registry = scenario.mkdir("registry");
    let check_args = ["--check", "--page", &page, "--reproductions", &registry];
    write_record(&scenario, "seen_entry", Some("observed"), Some(seen));
    write_record(
        &scenario,
        "unseen_entry",
        Some("unobserved"),
        Some(UNOBSERVED_BASIS),
    );
    write_record(
        &scenario,
        "broken_entry",
        Some("malformed"),
        Some("POS is an Integer (VCF 4.2 specification, section 1.4.1, POS): the POS is text that is no number"),
    );
    let ran = scenario.run(&check_args);
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert!(ran.stdout.contains(
        "3 entries have a reproduction verified against 116.3: 1 observed, 1 unobserved, 1 malformed"
    ));
    scenario.golden("classes", &ran);

    write_record(&scenario, "unseen_entry", None, None);
    let ran = scenario.run(&check_args);
    assert_eq!(ran.code, 1);
    assert!(ran
        .stderr
        .contains("ERROR: [reevaluate] unseen_entry: input_class is missing"));
    scenario.golden("classes_missing", &ran);
    write_record(
        &scenario,
        "unseen_entry",
        Some("synthetic"),
        Some(UNOBSERVED_BASIS),
    );
    let ran = scenario.run(&check_args);
    assert_eq!(ran.code, 1);
    assert!(ran
        .stderr
        .contains("is not one of observed, unobserved, malformed"));
    scenario.golden("classes_fourth", &ran);

    // A tier 3 entry whose page records its shape in real data is observed.
    write_record(
        &scenario,
        "unseen_entry",
        Some("observed"),
        Some(SIGHTING_BASIS),
    );
    let ran = scenario.run(&check_args);
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert!(ran.stdout.contains(
        "3 entries have a reproduction verified against 116.3: 2 observed, 0 unobserved, 1 malformed"
    ));
    scenario.golden("classes_sighting", &ran);
    // A class the page's tier and section contradict is one line each, the page's line
    // named.
    write_record(
        &scenario,
        "seen_entry",
        Some("malformed"),
        Some("the rule the input breaks"),
    );
    write_record(
        &scenario,
        "broken_entry",
        Some("unobserved"),
        Some(UNOBSERVED_BASIS),
    );
    let ran = scenario.run(&check_args);
    assert_eq!(ran.code, 1);
    scenario.golden("classes_contradicted", &ran);
}

/// `--check` reports an entry whose record has no `implementation` table and refuses a
/// table whose shape the schema refuses, a path the repository lacks and an item its file
/// does not define, one line each; the fixture records' own tables resolve.
#[test]
fn check_holds_each_implementation_table_to_the_repository() {
    let scenario = Scenario::new("implementation");
    let registry = scenario.mkdir("registry");
    let names = [
        "bare_entry",
        "lost_entry",
        "mapped_entry",
        "shapeless_entry",
    ];
    let rows: String = names
        .iter()
        .map(|key| {
            format!("| [`{key}`](x.md) | s | 3 · none on the cuts · 1 flag set, 0 corpora | 116.3 | [P](x.md) |\n")
        })
        .collect();
    let page = scenario.write(
        "page.md",
        format!("# Index\n\n## Section\n\n{HEADER}{rows}\n## How to read the counts\n\nprose\n")
            .as_bytes(),
    );
    let source = "tools/vep-divergence-reevaluate/tests/fixtures/implementation/recorded.rs";
    let tables = [
        ("bare_entry", String::new()),
        (
            "lost_entry",
            format!(
                "\n[implementation]\ncode = [\"tools/vep-divergence-reevaluate/tests/fixtures/implementation/absent.rs::downstream_row\", \"{source}::upstream_row\"]\ntests = [\"{source}::the_downstream_row_is_written\"]\n"
            ),
        ),
        (
            "mapped_entry",
            format!(
                "\n[implementation]\nnone = \"The row is the ordinary one.\"\ntests = [\"{source}::tests::the_downstream_row_is_written\"]\n"
            ),
        ),
        (
            "shapeless_entry",
            format!("\n[implementation]\ncode = [\"{source}\"]\nnone = \"The row is the ordinary one.\"\n"),
        ),
    ];
    for (name, table) in &tables {
        write_record(&scenario, name, Some("unobserved"), Some(UNOBSERVED_BASIS));
        let record = scenario.read(&format!("registry/{name}.toml"));
        let (record, _) = record
            .split_once("\n[implementation]\n")
            .expect("the fixture record's table");
        scenario.write(
            &format!("registry/{name}.toml"),
            format!("{record}\n{table}").as_bytes(),
        );
    }
    let ran = scenario.run(&["--check", "--page", &page, "--reproductions", &registry]);
    assert_eq!(ran.code, 1);
    assert_eq!(ran.stdout, "");
    assert_eq!(
        ran.stderr.lines().collect::<Vec<_>>(),
        [
            "ERROR: [reevaluate] bare_entry: the record has no implementation table".to_string(),
            "ERROR: [reevaluate] lost_entry: implementation.code names 'tools/vep-divergence-reevaluate/tests/fixtures/implementation/absent.rs::downstream_row': the repository holds no file tools/vep-divergence-reevaluate/tests/fixtures/implementation/absent.rs".to_string(),
            format!("ERROR: [reevaluate] lost_entry: implementation.code names '{source}::upstream_row': {source} defines no item upstream_row"),
            format!("ERROR: [reevaluate] lost_entry: implementation.tests names '{source}::the_downstream_row_is_written': {source} defines no item the_downstream_row_is_written"),
            format!("ERROR: [reevaluate] shapeless_entry: implementation.code names '{source}', which is not <path>::<item>"),
            "ERROR: [reevaluate] shapeless_entry: implementation holds code and none: none is for an entry with no code of its own".to_string(),
            "ERROR: [reevaluate] shapeless_entry: implementation needs tests".to_string(),
        ]
    );
    scenario.golden("check_implementation", &ran);
}

#[test]
fn list_prints_every_reproduction() {
    let scenario = Scenario::new("list");
    let ran = scenario.run(&["--list", "--reproductions", &fixture_records()]);
    assert_eq!(ran.code, 0);
    assert_eq!(
        ran.stdout,
        "closed_fixture\t116.3\tdocs/intended-divergences.md\nstill_present_fixture\t116.3\tdocs/intended-divergences.md\n"
    );
    scenario.golden("list", &ran);
    // A record that lacks a field lists `?` in its place, and --list comes before --check.
    scenario.write("records/bare.toml", b"notes = 'no entry, no page'\n");
    scenario.write(
        "records/typed.toml",
        b"verified_against = 116.2\npage = ['a', 1]\n",
    );
    let ran = scenario.run(&[
        "--check",
        "--list",
        "--reproductions",
        &scenario.path("records"),
    ]);
    assert_eq!(ran.code, 0);
    scenario.golden("list_bare", &ran);
}

#[test]
fn running_needs_an_engine_and_a_binary() {
    let scenario = Scenario::new("needs");
    let ran = scenario.run(&["--vep", "/bin/true"]);
    assert_eq!(ran.code, 2);
    assert!(ran.stderr.ends_with(
        "vep-divergence-reevaluate: error: --vep and --image (or --ensembl-vep) are required to run\n"
    ));
    assert_eq!(ran.stdout, "");
    let records = fixture_records();
    let base = [
        "--vep",
        "/bin/true",
        "--ensembl-vep",
        "/bin/true",
        "--reproductions",
        &records,
    ];
    let with = |extra: &[&str]| {
        let mut args = base.to_vec();
        args.extend(extra);
        scenario.run(&args)
    };
    let ran = with(&["--extra-dir", "loftee"]);
    assert_eq!(ran.code, 1);
    assert!(ran.stderr.contains("--extra-dir takes NAME=PATH"));
    scenario.golden("extra_dir_no_path", &ran);
    let ran = with(&["--extra-dir", &format!("Loftee={}", scenario.root)]);
    assert_eq!(ran.code, 1);
    scenario.golden("extra_dir_name", &ran);
    let ran = with(&[
        "--extra-dir",
        &format!("loftee={}", scenario.path("absent")),
    ]);
    assert_eq!(ran.code, 1);
    scenario.golden("extra_dir_absent", &ran);
    let ran = with(&["--tool-dir", &scenario.path("absent")]);
    assert_eq!(ran.code, 1);
    assert!(ran.stderr.contains("--tool-dir"));
    scenario.golden("tool_dir_absent", &ran);
}

/// Every reproduction of the repository validates, is verified against the release the
/// index's `Verified against` cells name, and matches a row of the index page: `--check`
/// finds the records and the page by its default paths.
#[test]
fn the_committed_reproductions_are_valid_and_cover_the_page() {
    let scenario = Scenario::new("committed");
    let ran = scenario.run(&["--check", "--reference", "116.3"]);
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert_eq!(ran.stderr, "");
    let (count, rest) = ran.stdout.split_once(' ').expect("a count");
    assert!(
        count.parse::<usize>().is_ok_and(|n| n > 0),
        "{}",
        ran.stdout
    );
    assert!(
        rest.starts_with("entries have a reproduction verified against 116.3: "),
        "{}",
        ran.stdout
    );
    let listed = scenario.run(&["--list"]);
    assert_eq!(listed.code, 0);
    assert_eq!(listed.stdout.lines().count().to_string(), count);
}

#[test]
fn help_states_every_option() {
    let scenario = Scenario::new("help");
    for flag in ["-h", "--help"] {
        let ran = scenario.run(&[flag]);
        scenario.golden("help", &ran);
    }
}

#[test]
fn a_command_line_that_is_refused_exits_2() {
    let scenario = Scenario::new("refused");
    for (args, message) in [
        (&["--vep"][..], "argument --vep: expected one argument"),
        (
            &["--jobs", "x"][..],
            "argument --jobs: invalid int value: 'x'",
        ),
        (
            &["--keep=1"][..],
            "argument --keep: ignored explicit argument '1'",
        ),
        (
            &["--re", "1"][..],
            "ambiguous option: --re could match --reference, --reproductions",
        ),
        (&["stray", "--list"][..], "unrecognized arguments: stray"),
    ] {
        let ran = scenario.run(args);
        assert_eq!(ran.code, 2, "{args:?}");
        assert_eq!(ran.stdout, "");
        assert!(ran
            .stderr
            .starts_with("usage: vep-divergence-reevaluate [-h] "));
        assert!(
            ran.stderr
                .ends_with(&format!("\nvep-divergence-reevaluate: error: {message}\n")),
            "{}",
            ran.stderr
        );
    }
}
