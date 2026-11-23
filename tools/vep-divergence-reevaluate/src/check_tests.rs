// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::page::MALFORMED_SECTION;
use crate::record::load_reproductions;
use crate::repo::repo_root;
use crate::testing::{Scratch, SIGHTING_BASIS, UNOBSERVED_BASIS};

const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/reproductions");
const HEADER: &str = "| Entry | In one sentence | Tier and count | Verified against | Page |\n| --- | --- | --- | --- | --- |\n";

fn errors(page: &str, reference: &str, directory: &str) -> Vec<String> {
    check(page, reference, directory, &repo_root(), &Schema::new())
        .expect("a page and records that read")
}

#[test]
fn check_holds_the_page_and_the_reproductions_to_one_reference() {
    let scratch = Scratch::new("reference");
    let page = |rows: &[(&str, &str)]| {
        let body: String = rows
            .iter()
            .map(|(key, verified)| {
                format!("| [`{key}`](x.md) | s | 3 · none on the cuts · 1 flag set, 0 corpora | {verified} | [P](x.md) |\n")
            })
            .collect();
        scratch.write(
            "page.md",
            format!(
                "# Index\n\n## Section\n\n{HEADER}{body}\n## How to read the counts\n\nprose\n"
            )
            .as_bytes(),
        )
    };
    let both = [
        ("still_present_fixture", "116.3"),
        ("closed_fixture", "116.3"),
    ];
    let path = page(&both);
    assert_eq!(errors(&path, "116.3", FIXTURES), Vec::<String>::new());
    page(&[
        ("still_present_fixture", "116.3"),
        ("closed_fixture", "116.1"),
        ("orphan_entry", "116.3"),
    ]);
    assert_eq!(
        errors(&path, "116.3", FIXTURES),
        [
            "closed_fixture: the page says verified against '116.1', the reference is 116.3"
                .to_string(),
            format!("orphan_entry: no reproduction under {FIXTURES}"),
        ]
    );
    page(&[("still_present_fixture", "116.3")]);
    assert_eq!(
        errors(&path, "116.3", FIXTURES),
        ["closed_fixture: a reproduction with no row on the page"]
    );
    page(&both);
    assert_eq!(
        errors(&path, "116.2", FIXTURES),
        [
            "still_present_fixture: verified against 116.3, the reference is 116.2",
            "still_present_fixture: the page says verified against '116.3', the reference is 116.2",
            "closed_fixture: verified against 116.3, the reference is 116.2",
            "closed_fixture: the page says verified against '116.3', the reference is 116.2",
        ]
    );
    let missing = scratch.path("missing.md");
    assert_eq!(
        check(&missing, "116.3", FIXTURES, &repo_root(), &Schema::new()),
        Err(format!("[Errno 2] No such file or directory: '{missing}'"))
    );
}

/// `--check` passes cells that count cache transcripts against one total, whatever total a
/// cell of another unit names, and refuses two totals in one line that names each with its
/// rows and its first entry, in page order.
#[test]
fn check_refuses_an_index_whose_cache_transcripts_cells_name_two_totals() {
    let scratch = Scratch::new("totals");
    let page = |rows: &[(&str, &str)]| {
        let body: String = rows
            .iter()
            .map(|(key, figure)| {
                format!("| [`{key}`](x.md) | s | 3 · {figure} · 1 flag set, 0 corpora | 116.3 | [P](x.md) |\n")
            })
            .collect();
        scratch.write(
            "page.md",
            format!(
                "# Index\n\n## Section\n\n{HEADER}{body}\n## How to read the counts\n\nprose\n"
            )
            .as_bytes(),
        )
    };
    let none = Vec::<String>::new();
    let path = page(&[
        ("still_present_fixture", "0 of 1,200 cache transcripts"),
        ("closed_fixture", "0 of 1,200 cache transcripts"),
    ]);
    assert_eq!(errors(&path, "116.3", FIXTURES), none);
    page(&[
        ("still_present_fixture", "0 of 1,200 cache transcripts"),
        ("closed_fixture", "0 of 90 regulatory features"),
    ]);
    assert_eq!(errors(&path, "116.3", FIXTURES), none);

    page(&[
        ("still_present_fixture", "0 of 1,200 cache transcripts"),
        ("closed_fixture", "0 of 1,250 cache transcripts"),
    ]);
    assert_eq!(
        errors(&path, "116.3", FIXTURES),
        ["the page's `Tier and count` cells name 2 totals of cache transcripts: 1,200 in 1 row, the first still_present_fixture; 1,250 in 1 row, the first closed_fixture"]
    );

    let registry = scratch.mkdir("registry");
    for name in ["first_entry", "second_entry", "third_entry"] {
        write_record(&scratch, name, Some("unobserved"), Some(UNOBSERVED_BASIS));
    }
    page(&[
        ("first_entry", "0 of 1,200 cache transcripts"),
        ("second_entry", "0 of 1,250 cache transcripts"),
        ("third_entry", "0 of 1,200 cache transcripts"),
    ]);
    assert_eq!(
        errors(&path, "116.3", &registry),
        ["the page's `Tier and count` cells name 2 totals of cache transcripts: 1,200 in 2 rows, the first first_entry; 1,250 in 1 row, the first second_entry"]
    );
}

/// The order checker's table of units pairs the basis of a count over the committed caches'
/// transcripts with the unit this check reads at the end of a figure.
#[test]
fn the_cache_transcripts_unit_is_the_one_the_order_checker_writes() {
    let checker = std::fs::read_to_string(format!(
        "{}/scripts/divergences/check_divergence_order.py",
        repo_root()
    ))
    .expect("the order checker");
    let pair = format!("(\"transcripts of the committed caches\", \"{CACHE_TRANSCRIPTS}\")");
    assert!(checker.contains(&pair), "the order checker lacks {pair}");
}

/// A copy of the still-present fixture under another name with the given input class and
/// basis; None leaves the field out.
fn write_record(scratch: &Scratch, name: &str, input_class: Option<&str>, basis: Option<&str>) {
    let fixture = std::fs::read_to_string(format!("{FIXTURES}/still_present_fixture.toml"))
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
    scratch.write(
        &format!("registry/{name}.toml"),
        format!("{}\n", lines.join("\n")).as_bytes(),
    );
}

#[test]
fn check_holds_each_input_class_to_the_pages_tier_and_section() {
    let scratch = Scratch::new("classes");
    let seen = "2 · 98.6 per million · 15 flag sets, 20 corpora";
    let unseen = "3 · none on the cuts · 1 flag set, 0 corpora";
    let page_path = scratch.path("page.md");
    let registry = scratch.mkdir("registry");
    let table = |rows: &[(&str, &str)]| -> String {
        let body: String = rows
            .iter()
            .map(|(key, cell)| format!("| [`{key}`](x.md) | s | {cell} | 116.3 | [P](x.md) |\n"))
            .collect();
        format!("{HEADER}{body}")
    };
    let write_page = |main: &[(&str, &str)], malformed: &[(&str, &str)]| {
        let mut text = format!("# Index\n\n## Section\n\n{}", table(main));
        if !malformed.is_empty() {
            text.push_str(&format!(
                "\n## {MALFORMED_SECTION}\n\n### Section\n\n{}",
                table(malformed)
            ));
        }
        text.push_str("\n## How to read the counts\n\nprose\n");
        scratch.write("page.md", text.as_bytes());
    };
    let line = |key: &str| -> usize {
        let text = std::fs::read_to_string(&page_path).unwrap();
        let wanted = format!("`{key}`");
        text.lines().position(|l| l.contains(&wanted)).unwrap() + 1
    };
    let found = || errors(&page_path, "116.3", &registry);
    let none = Vec::<String>::new();

    write_record(&scratch, "seen_entry", Some("observed"), Some(seen));
    write_record(
        &scratch,
        "unseen_entry",
        Some("unobserved"),
        Some(UNOBSERVED_BASIS),
    );
    write_record(
        &scratch,
        "broken_entry",
        Some("malformed"),
        Some("POS is an Integer (VCF 4.2 specification, section 1.4.1, POS): the POS is text that is no number"),
    );
    let main = [("seen_entry", seen), ("unseen_entry", unseen)];
    let malformed = [("broken_entry", unseen)];
    write_page(&main, &malformed);
    assert_eq!(found(), none);
    let classes = input_classes(&load_reproductions(None, &registry).unwrap());
    assert_eq!(
        crate::page::class_counts(classes.iter().map(|(_, class)| class.as_str())),
        "1 observed, 1 unobserved, 1 malformed"
    );

    // A cell that states no tier is not read as tier 3: the basis is held to it as written.
    write_page(
        &[("seen_entry", "tier 2"), ("unseen_entry", unseen)],
        &malformed,
    );
    assert_eq!(
        found(),
        [format!(
            "seen_entry: input_class_basis is '{seen}', the page's `Tier and count` cell is 'tier 2'"
        )]
    );
    write_page(&main, &malformed);

    write_record(&scratch, "unseen_entry", None, None);
    assert_eq!(
        found(),
        [
            "unseen_entry: input_class is missing",
            "unseen_entry: input_class_basis is missing"
        ]
    );
    write_record(&scratch, "unseen_entry", Some("unobserved"), None);
    assert_eq!(found(), ["unseen_entry: input_class_basis is missing"]);
    write_record(
        &scratch,
        "unseen_entry",
        Some("synthetic"),
        Some(UNOBSERVED_BASIS),
    );
    assert_eq!(
        found(),
        ["unseen_entry: input_class 'synthetic' is not one of observed, unobserved, malformed"]
    );
    write_record(
        &scratch,
        "unseen_entry",
        Some("unobserved"),
        Some(UNOBSERVED_BASIS),
    );

    let stale = "2 · 98.6 per million · 14 flag sets, 20 corpora";
    write_record(&scratch, "seen_entry", Some("observed"), Some(stale));
    assert_eq!(
        found(),
        [format!(
            "seen_entry: input_class_basis is '{stale}', the page's `Tier and count` cell is '{seen}'"
        )]
    );
    write_record(
        &scratch,
        "seen_entry",
        Some("unobserved"),
        Some(UNOBSERVED_BASIS),
    );
    assert_eq!(
        found(),
        [format!(
            "line {}: seen_entry is tier 2 with input_class 'unobserved': a tier 1 or tier 2 entry is observed",
            line("seen_entry")
        )]
    );
    write_record(&scratch, "seen_entry", Some("synthetic"), Some(seen));
    assert_eq!(
        found(),
        ["seen_entry: input_class 'synthetic' is not one of observed, unobserved, malformed"]
    );
    write_record(&scratch, "seen_entry", Some("observed"), Some(seen));

    // A tier 3 entry whose page records its shape in real data is observed, and its basis
    // names that data in place of its cell.
    write_record(
        &scratch,
        "unseen_entry",
        Some("observed"),
        Some(SIGHTING_BASIS),
    );
    assert_eq!(found(), none);
    let classes = input_classes(&load_reproductions(None, &registry).unwrap());
    assert_eq!(
        crate::page::class_counts(classes.iter().map(|(_, class)| class.as_str())),
        "2 observed, 0 unobserved, 1 malformed"
    );
    write_record(&scratch, "unseen_entry", Some("observed"), Some(seen));
    assert_eq!(
        found(),
        [format!(
            "unseen_entry: input_class_basis is '{seen}': an observed entry of tier 3 reads `seen in <the real data its page names>: <the figure its page measures there>`"
        )]
    );
    write_record(&scratch, "unseen_entry", Some("observed"), None);
    assert_eq!(found(), ["unseen_entry: input_class_basis is missing"]);
    write_record(&scratch, "unseen_entry", Some("observed"), Some(unseen));
    assert_eq!(
        found(),
        ["unseen_entry: input_class_basis of an observed entry is neither a tier 1 or tier 2 `Tier and count` cell nor `seen in <the real data its page names>: <the figure its page measures there>`"]
    );
    write_record(
        &scratch,
        "seen_entry",
        Some("observed"),
        Some(SIGHTING_BASIS),
    );
    write_record(
        &scratch,
        "unseen_entry",
        Some("observed"),
        Some(SIGHTING_BASIS),
    );
    assert_eq!(
        found(),
        [format!(
            "seen_entry: input_class_basis is '{SIGHTING_BASIS}', the page's `Tier and count` cell is '{seen}'"
        )]
    );
    write_record(&scratch, "seen_entry", Some("observed"), Some(seen));

    write_page(
        &[("seen_entry", seen)],
        &[("broken_entry", unseen), ("unseen_entry", unseen)],
    );
    assert_eq!(
        found(),
        [format!(
            "line {}: unseen_entry stands under '{MALFORMED_SECTION}' with input_class 'observed'",
            line("unseen_entry")
        )]
    );
    write_record(
        &scratch,
        "unseen_entry",
        Some("unobserved"),
        Some(UNOBSERVED_BASIS),
    );
    assert_eq!(
        found(),
        [format!(
            "line {}: unseen_entry stands under '{MALFORMED_SECTION}' with input_class 'unobserved'",
            line("unseen_entry")
        )]
    );
    write_page(
        &[
            ("seen_entry", seen),
            ("broken_entry", unseen),
            ("unseen_entry", unseen),
        ],
        &[],
    );
    assert_eq!(
        found(),
        [format!(
            "line {}: broken_entry is malformed and stands under 'Section', outside '{MALFORMED_SECTION}'",
            line("broken_entry")
        )]
    );
}

/// The registry record `name` with its `implementation` table written again from `table`,
/// the lines under the table's heading; `None` takes the table out.
fn write_table(scratch: &Scratch, name: &str, table: Option<&str>) {
    write_record(scratch, name, Some("unobserved"), Some(UNOBSERVED_BASIS));
    let path = scratch.path(&format!("registry/{name}.toml"));
    let text = std::fs::read_to_string(&path).expect("the record");
    let (record, _) = text
        .split_once("\n[implementation]\n")
        .expect("the fixture record's table");
    let table = table.map_or_else(String::new, |table| format!("\n[implementation]\n{table}"));
    std::fs::write(&path, format!("{record}\n{table}")).expect("the record");
}

/// `--check` reports an entry whose record has no `implementation` table and refuses a
/// row that names a path the repository lacks or an item its file does not define.
#[test]
fn check_holds_each_implementation_table_to_the_repository() {
    let scratch = Scratch::new("implementation");
    let registry = scratch.mkdir("registry");
    let tree = scratch.mkdir("tree");
    scratch.write(
        "tree/crates/engine/src/rows.rs",
        b"pub struct Writer;\n\nimpl Writer {\n    pub fn row(&self) -> u8 {\n        1\n    }\n}\n\n#[cfg(test)]\nmod tests {\n    #[test]\n    fn the_row_is_written() {}\n}\n",
    );
    scratch.write(
        "tree/scripts/test_rows.py",
        b"class TestRows:\n    def test_row(self):\n        pass\n\n\ndef test_plain():\n    pass\n",
    );
    scratch.write(
        "tree/scripts/test_rows.sh",
        b"check_rows() {\n    true\n}\nfunction check_more {\n    true\n}\n",
    );
    scratch.write("tree/docs/page.md", b"# A page\n");
    let names = ["mapped_entry", "other_entry"];
    let rows: String = names
        .iter()
        .map(|key| {
            format!("| [`{key}`](x.md) | s | 3 · none on the cuts · 1 flag set, 0 corpora | 116.3 | [P](x.md) |\n")
        })
        .collect();
    let page = scratch.write(
        "page.md",
        format!("# Index\n\n## Section\n\n{HEADER}{rows}\n## How to read the counts\n\nprose\n")
            .as_bytes(),
    );
    let against_the_tree = || {
        check(&page, "116.3", &registry, &tree, &Schema::new())
            .expect("a page and records that read")
    };
    let resolved = "code = [\"crates/engine/src/rows.rs::Writer::row\"]\ntests = [\"crates/engine/src/rows.rs::tests::the_row_is_written\", \"scripts/test_rows.py::test_plain\", \"scripts/test_rows.py::TestRows\", \"scripts/test_rows.sh::check_rows\", \"scripts/test_rows.sh::check_more\"]\nrules = [\"crates/engine/src/rows.rs::Writer\"]\n";
    write_table(&scratch, "mapped_entry", Some(resolved));
    write_table(
        &scratch,
        "other_entry",
        Some("none = \"The row is the ordinary one.\"\ntests = [\"crates/engine/src/rows.rs::tests::the_row_is_written\"]\n"),
    );
    assert_eq!(against_the_tree(), Vec::<String>::new());

    write_table(&scratch, "other_entry", None);
    assert_eq!(
        against_the_tree(),
        ["other_entry: the record has no implementation table"]
    );

    write_table(
        &scratch,
        "other_entry",
        Some("code = [\"crates/engine/src/gone.rs::Writer::row\"]\ntests = [\"crates/engine/src/rows.rs::tests::the_row_is_written\"]\n"),
    );
    assert_eq!(
        against_the_tree(),
        ["other_entry: implementation.code names 'crates/engine/src/gone.rs::Writer::row': the repository holds no file crates/engine/src/gone.rs"]
    );

    write_table(
        &scratch,
        "other_entry",
        Some("code = [\"crates/engine/src/rows.rs::Writer::column\", \"crates/engine/src/rows.rs::row\"]\ntests = [\"crates/engine/src/rows.rs::the_row_is_written\", \"scripts/test_rows.py::test_row_missing\", \"scripts/test_rows.sh::check\", \"docs/page.md::heading\", \"crates/engine/src::rows\"]\n"),
    );
    assert_eq!(
        against_the_tree(),
        [
            "other_entry: implementation.code names 'crates/engine/src/rows.rs::Writer::column': crates/engine/src/rows.rs defines no item Writer::column",
            "other_entry: implementation.code names 'crates/engine/src/rows.rs::row': crates/engine/src/rows.rs defines no item row",
            "other_entry: implementation.tests names 'crates/engine/src/rows.rs::the_row_is_written': crates/engine/src/rows.rs defines no item the_row_is_written",
            "other_entry: implementation.tests names 'scripts/test_rows.py::test_row_missing': scripts/test_rows.py defines no function or class test_row_missing",
            "other_entry: implementation.tests names 'scripts/test_rows.sh::check': scripts/test_rows.sh defines no function or class check",
            "other_entry: implementation.tests names 'docs/page.md::heading': docs/page.md is neither a Rust, a Python nor a shell file",
            "other_entry: implementation.tests names 'crates/engine/src::rows': crates/engine/src is neither a Rust, a Python nor a shell file",
        ]
    );

    // A row the schema refuses for its shape is named once, by the schema's line.
    write_table(
        &scratch,
        "other_entry",
        Some("code = [\"crates/engine/src/rows.rs\"]\ntests = [\"crates/engine/src/rows.rs::tests::the_row_is_written\"]\n"),
    );
    assert_eq!(
        against_the_tree(),
        ["other_entry: implementation.code names 'crates/engine/src/rows.rs', which is not <path>::<item>"]
    );
}

/// Every reproduction of the repository validates, is verified against the reference, and
/// matches a row of the index page.
#[test]
fn the_committed_reproductions_are_valid_and_cover_the_page() {
    let repo = repo_root();
    let directory = format!("{repo}/scripts/divergences/reproductions");
    let repros = load_reproductions(None, &directory).expect("the committed records");
    assert!(!repros.is_empty(), "no reproduction committed");
    let schema = Schema::new();
    for repro in &repros {
        assert_eq!(
            schema.validate(repro),
            Vec::<String>::new(),
            "{}",
            repro.name
        );
    }
    assert_eq!(
        errors(
            &format!("{repo}/docs/intended-divergences.md"),
            "116.3",
            &directory
        ),
        Vec::<String>::new()
    );
}

/// The repository paths a `prepare` step names: each `$REPO/<path>` as the step writes
/// it, then each `/<mount>/<file>` word under a directory the step mounts from the
/// repository (`-v "$REPO/<directory>:/<mount>"`), as that directory's file. A path
/// ends at white space, a quote, a colon or a shell operator. A path that holds a
/// shell variable or a wildcard names no one file and is left out.
fn repository_paths(step: &str) -> Vec<String> {
    let ends = |c: char| {
        c.is_whitespace()
            || matches!(
                c,
                '"' | '\'' | ':' | ';' | '(' | ')' | '|' | '&' | '<' | '>' | '`'
            )
    };
    let mut paths = Vec::new();
    let mut mounts = Vec::new();
    for (at, marker) in step.match_indices("$REPO/") {
        let rest = &step[at + marker.len()..];
        let end = rest.find(ends).unwrap_or(rest.len());
        let path = &rest[..end];
        if path.contains(['$', '*']) {
            continue;
        }
        paths.push(path.to_string());
        if let Some(mount) = rest[end..].strip_prefix(":/") {
            let name = &mount[..mount.find(ends).unwrap_or(mount.len())];
            mounts.push((format!("/{name}/"), path));
        }
    }
    for word in step.split(ends) {
        for (mount, directory) in &mounts {
            if let Some(file) = word
                .strip_prefix(mount.as_str())
                .filter(|file| !file.contains(['$', '*']))
            {
                paths.push(format!("{directory}/{file}"));
            }
        }
    }
    paths
}

#[test]
fn repository_paths_reads_a_step_and_what_it_mounts() {
    assert_eq!(
        repository_paths(concat!(
            r#"docker run --rm -v "$REPO/a/b:/in:ro" -v "$STAGE:/stage" "$IMAGE" perl /in/c.pl "#,
            r#"/stage/cache > "$STAGE/log" 2>&1; cp "$REPO/d/$name.txt" $REPO/e/*.gz .; "#,
            r#"cargo run --manifest-path "$REPO/tools/t/Cargo.toml" -- "$STAGE/cache""#,
        )),
        ["a/b", "tools/t/Cargo.toml", "a/b/c.pl"]
    );
    assert_eq!(repository_paths("cp $REPO/f .; cat /f"), ["f"]);
    assert_eq!(
        repository_paths(concat!(
            r#"docker run -v "$REPO/a:/in" "$IMAGE" sh -c 'cat /in/b; ls /in/$name /in/*.pl' "#,
            r#"<$REPO/c>"$STAGE/d""#,
        )),
        ["a", "c", "a/b"]
    );
}

/// A `prepare` step of a committed reproduction names only files and directories the
/// repository holds, also where it reads one through a directory it mounts into a
/// container. A step that names another fails when the entry is run, which `--check`
/// does not show, since it runs no step.
#[test]
fn the_committed_prepare_steps_name_paths_the_repository_holds() {
    let repo = repo_root();
    let directory = format!("{repo}/scripts/divergences/reproductions");
    let repros = load_reproductions(None, &directory).expect("the committed records");
    let mut named = 0;
    let mut absent = Vec::new();
    for repro in &repros {
        let Some(Value::List(steps)) = repro.data.get("prepare") else {
            continue;
        };
        for path in steps
            .iter()
            .filter_map(Value::as_str)
            .flat_map(repository_paths)
        {
            named += 1;
            if !std::path::Path::new(&format!("{repo}/{path}")).exists() {
                absent.push(format!("{}: {path}", repro.name));
            }
        }
    }
    assert!(named > 0, "no prepare step names a repository path");
    assert_eq!(absent, Vec::<String>::new());
}

/// Whether a command can leave before its input ends: `head`; an `awk` whose program
/// exits outside its `END` rule; a `grep` whose leading options make it quiet or stop it
/// at a count; a `sed` whose script opens a command with a line number and `q`.
fn leaves_before_its_input_ends(command: &str) -> bool {
    let assigns = |word: &str| {
        word.split_once('=').is_some_and(|(name, _)| {
            !name.is_empty() && name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')
        })
    };
    let mut words = command
        .split_whitespace()
        .skip_while(|word| *word == "command" || assigns(word));
    let Some(name) = words.next() else {
        return false;
    };
    match name.rsplit('/').next().unwrap_or(name) {
        "head" => true,
        "awk" | "gawk" | "mawk" | "nawk" => awk_exits_outside_its_end_rule(command),
        "grep" | "egrep" | "fgrep" | "zgrep" => {
            words
                .take_while(|word| word.starts_with('-'))
                .any(|option| match option.strip_prefix("--") {
                    Some(long) => {
                        matches!(long, "quiet" | "silent") || long.starts_with("max-count")
                    }
                    None => option.contains(['q', 'm']),
                })
        }
        "sed" => {
            let bytes = command.as_bytes();
            (1..bytes.len()).any(|at| {
                let digits = bytes[..at]
                    .iter()
                    .rev()
                    .take_while(|c| c.is_ascii_digit())
                    .count();
                bytes[at] == b'q'
                    && digits > 0
                    && at.checked_sub(digits + 1).is_some_and(|before| {
                        matches!(bytes[before], b' ' | b'\t' | b'\'' | b'"' | b';' | b'{')
                    })
                    && !bytes.get(at + 1).is_some_and(u8::is_ascii_alphanumeric)
            })
        }
        _ => false,
    }
}

/// Whether an `awk` command holds an `exit` outside an `END` rule: the outermost brace
/// left open before each `exit` opens its rule's action, and the pattern before that
/// brace is read.
fn awk_exits_outside_its_end_rule(command: &str) -> bool {
    let bytes = command.as_bytes();
    let in_a_name = |c: &u8| c.is_ascii_alphanumeric() || *c == b'_';
    command.match_indices("exit").any(|(at, word)| {
        if at.checked_sub(1).is_some_and(|i| in_a_name(&bytes[i]))
            || bytes.get(at + word.len()).is_some_and(in_a_name)
        {
            return false;
        }
        let mut closed = 0;
        let mut action = None;
        for i in (0..at).rev() {
            match bytes[i] {
                b'}' => closed += 1,
                b'{' if closed == 0 => action = Some(i),
                b'{' => closed -= 1,
                _ => {}
            }
        }
        !action.is_some_and(|i| command[..i].trim_end().ends_with("END"))
    })
}

/// A statement's text past the words that can stand before its first command (`{`, `!`,
/// `then`, `do`, `else` and the four that follow), and whether one of them is `if`,
/// `elif`, `while` or `until`, whose statement is a condition. The text may open with the
/// operator that ends the statement before it.
fn past_the_opening_words(statement: &str) -> (&str, bool) {
    let apart = |c: char| c.is_whitespace() || matches!(c, ';' | '&' | '|' | ')');
    let mut rest = statement.trim_start_matches(apart);
    let mut condition = false;
    loop {
        let (word, after) = rest.split_at(rest.find(char::is_whitespace).unwrap_or(rest.len()));
        match word {
            "if" | "elif" | "while" | "until" => condition = true,
            "{" | "!" | "then" | "do" | "else" => {}
            _ => return (rest, condition),
        }
        rest = after.trim_start();
    }
}

/// Whether a statement opens with an assignment (`name=`, `name+=`), so that a command
/// substitution anywhere in its values gives it its status: an assignment alone has the
/// status of the last substitution it holds. A substitution that is a word of a command,
/// a declaration among them (`local name=$(`), gives its status to nothing.
fn opens_with_an_assignment(statement: &str) -> bool {
    let (rest, _) = past_the_opening_words(statement);
    let value = rest.trim_start_matches(|c: char| c.is_ascii_alphanumeric() || c == '_');
    value.len() < rest.len() && (value.starts_with('=') || value.starts_with("+="))
}

/// The readers of a `prepare` step that can end their writer before its last write
/// where the pipeline's status is read: each command after a pipe sign that can leave
/// before its input ends. The tool runs a step under `bash -euo pipefail`,
/// so a writer that has a write left when its reader leaves is ended by the closed
/// pipe and the pipeline's status is the writer's, 141, by the timing of the two.
///
/// The step is read as far as its quotes, comments, here-documents, subshells and
/// command substitutions go. Two pipelines are left out: one that follows `set +o
/// pipefail` in its own subshell or substitution, whose status is its last command's,
/// and one in a command substitution that is a word of a command, which gives its
/// status to nothing, unless `&&` or `||` follows the pipeline there or `if`, `elif`,
/// `while` or `until` opens its statement.
fn readers_that_leave_early(step: &str) -> Vec<String> {
    struct Frame {
        pipefail: bool,
        status_read: bool,
        double_quoted: bool,
        reader_from: Option<usize>,
        // Where the statement being read begins, and the early readers of its pipeline.
        statement_from: usize,
        early: Vec<String>,
    }
    let within = |frame: &Frame, status_read: bool, statement_from: usize| Frame {
        pipefail: frame.pipefail,
        status_read,
        double_quoted: false,
        reader_from: None,
        statement_from,
        early: Vec::new(),
    };
    let bytes = step.as_bytes();
    let mut found = Vec::new();
    // The reader a pipe sign or its pipeline's end closes. bash takes a backslash and
    // the line break after it out of a command.
    let reader = |frame: &mut Frame, end: usize| {
        if let Some(from) = frame.reader_from.take() {
            let text = step[from..end].replace("\\\n", "");
            let command = text.trim();
            if frame.pipefail && leaves_before_its_input_ends(command) {
                frame.early.push(command.to_string());
            }
        }
    };
    // The end of a pipeline: its early readers count where its status is read.
    let mut close = |frame: &mut Frame, end: usize, joined: bool| {
        reader(frame, end);
        let (_, condition) = past_the_opening_words(&step[frame.statement_from..end]);
        if frame.status_read || joined || condition {
            found.append(&mut frame.early);
        }
        frame.early.clear();
        frame.statement_from = end;
    };
    let mut frames = vec![Frame {
        pipefail: true,
        status_read: true,
        double_quoted: false,
        reader_from: None,
        statement_from: 0,
        early: Vec::new(),
    }];
    // The here-documents whose bodies begin after the line being read: each delimiter,
    // and whether leading tabs are taken off the body's lines.
    let mut here_documents: Vec<(String, bool)> = Vec::new();
    let mut at = 0;
    while at < bytes.len() {
        let next = bytes.get(at + 1).copied();
        let nested = frames.len() > 1;
        let frame = frames.last_mut().expect("the step's own frame");
        // The command a substitution stands in begins after the last pipe sign.
        let command_from = frame.reader_from.unwrap_or(frame.statement_from);
        if frame.double_quoted {
            match bytes[at] {
                b'\\' => at += 1,
                b'"' => frame.double_quoted = false,
                b'$' if next == Some(b'(') => {
                    let assigned = opens_with_an_assignment(&step[command_from..at]);
                    let inner = within(frame, frame.status_read && assigned, at + 2);
                    frames.push(inner);
                    at += 1;
                }
                _ => {}
            }
            at += 1;
            continue;
        }
        match bytes[at] {
            b'\\' => at += 1,
            b'\'' => {
                let escapes = at > 0 && bytes[at - 1] == b'$';
                at += 1;
                while at < bytes.len() && bytes[at] != b'\'' {
                    at += if escapes && bytes[at] == b'\\' { 2 } else { 1 };
                }
            }
            b'"' => frame.double_quoted = true,
            // A comment opens where `#` begins a word.
            b'#' if at == 0
                || bytes[at - 1].is_ascii_whitespace()
                || matches!(bytes[at - 1], b';' | b'&' | b'|' | b'(' | b')') =>
            {
                at = step[at..].find('\n').map_or(bytes.len(), |end| at + end);
                continue;
            }
            b'<' if next == Some(b'<') && bytes.get(at + 2) == Some(&b'<') => at += 2,
            b'<' if next == Some(b'<') => {
                let word = step[at + 2..].trim_start_matches([' ', '\t']);
                let tabs = word.starts_with('-');
                let word = word.trim_start_matches('-').trim_start_matches([' ', '\t']);
                let ends =
                    |c: char| c.is_whitespace() || matches!(c, ';' | '|' | '&' | ')' | '<' | '>');
                let end = word.find(ends).unwrap_or(word.len());
                let delimiter = word[..end].replace(['\'', '"', '\\'], "");
                here_documents.push((delimiter, tabs));
                at += 1;
            }
            b'$' if next == Some(b'(') => {
                let assigned = opens_with_an_assignment(&step[command_from..at]);
                let inner = within(frame, frame.status_read && assigned, at + 2);
                frames.push(inner);
                at += 1;
            }
            b'(' => {
                // The elements of an array assignment (`name=(`) stand in its statement.
                let array = at > 0 && bytes[at - 1] == b'=';
                let from = if array { command_from } else { at + 1 };
                let inner = within(frame, frame.status_read, from);
                frames.push(inner);
            }
            b')' => {
                close(frame, at, false);
                if frames.len() > 1 {
                    frames.pop();
                }
            }
            b'|' if next == Some(b'|') => {
                close(frame, at, true);
                at += 1;
            }
            b'|' => {
                reader(frame, at);
                at += usize::from(next == Some(b'&'));
                frame.reader_from = Some(at + 1);
            }
            b'&' => {
                close(frame, at, next == Some(b'&'));
                at += usize::from(next == Some(b'&'));
            }
            b';' => close(frame, at, false),
            // The body of a group or of a function opens with a statement of its own.
            b'{' if next.is_some_and(|c| c.is_ascii_whitespace()) => {
                frame.statement_from = at + 1;
            }
            b'\n' => {
                // A line that ends at a pipe sign, or at a comment after one, leaves the
                // pipeline's next command to a later line.
                let awaited = frame.reader_from.is_some_and(|from| {
                    let begun = step[from..at].trim_start();
                    begun.is_empty() || begun.starts_with('#')
                });
                if !awaited {
                    close(frame, at, false);
                }
                // A `<<` whose delimiter no later line holds opens no here-document, as
                // the shift of an arithmetic expansion.
                for (delimiter, tabs) in here_documents.drain(..) {
                    let mut end = at;
                    while end < bytes.len() {
                        let from = end + 1;
                        end = step[from..]
                            .find('\n')
                            .map_or(bytes.len(), |line| from + line);
                        let line = &step[from..end];
                        let line = if tabs {
                            line.trim_start_matches('\t')
                        } else {
                            line
                        };
                        if line == delimiter {
                            at = end;
                            break;
                        }
                    }
                }
                if awaited {
                    frame.reader_from = Some((at + 1).min(bytes.len()));
                } else {
                    frame.statement_from = at;
                }
            }
            // The step's own shell keeps `pipefail` past the body of a function it defines
            // and does not call, so the option is followed in a subshell or substitution.
            b's' if nested && !(at > 0 && bytes[at - 1].is_ascii_alphanumeric()) => {
                if step[at..].starts_with("set +o pipefail") {
                    frame.pipefail = false;
                } else if step[at..].starts_with("set -o pipefail") {
                    frame.pipefail = true;
                }
            }
            _ => {}
        }
        at += 1;
    }
    for frame in &mut frames {
        close(frame, bytes.len(), false);
    }
    found
}

#[test]
fn readers_that_leave_early_reads_where_a_pipelines_status_is_read() {
    let early = |step: &str| readers_that_leave_early(step);
    // A reader that leaves early, as a statement, in a condition and as an assignment's
    // value.
    assert_eq!(
        early(r#"gunzip -c "$CORPUS/a.gz" | awk '/^>MT/ { exit } { print }' > "$STAGE/a""#),
        [r#"awk '/^>MT/ { exit } { print }' > "$STAGE/a""#]
    );
    assert_eq!(
        early("if tabix a.gz 21:83-85 2> /dev/null | grep -q .; then exit 1; fi"),
        ["grep -q ."]
    );
    assert_eq!(
        early("first=$(zcat a.gz | head -n 1)\nlast=\"$(zcat a.gz | sed 3q)\""),
        ["head -n 1", "sed 3q"]
    );
    assert_eq!(
        early("if true; then first=$(zcat a.gz | head -n 1); fi"),
        ["head -n 1"]
    );
    assert_eq!(
        early("zcat a.gz | LC_ALL=C grep -m 1 x > b; zcat a.gz |& command fgrep --quiet x"),
        ["LC_ALL=C grep -m 1 x > b", "command fgrep --quiet x"]
    );
    // An assignment takes its status from a substitution anywhere in its values, in the
    // elements of an array, in the arm of a `case` and in the body of a function.
    assert_eq!(
        early(concat!(
            r#"line="a $(zcat a.gz | head -n 1)"; rows+=$(zcat a.gz | sed 3q); "#,
            "last=($(zcat a.gz | head -n 2))",
        )),
        ["head -n 1", "sed 3q", "head -n 2"]
    );
    assert_eq!(
        early(concat!(
            "case $x in a) first=$(zcat a.gz | head -n 1) ;; esac\n",
            "row() { first=$(zcat a.gz | sed 3q); }",
        )),
        ["head -n 1", "sed 3q"]
    );
    // A pipeline goes on past a line break after its pipe sign, a comment there and a
    // here-document's body, and past a backslash before a line break.
    assert_eq!(
        early(concat!(
            "zcat a.gz |\n  head -n 1\nzcat a.gz | # the first\n  sed \\\n3q\n",
            "cat <<GFF |\na | head\nGFF\ngrep -q a",
        )),
        ["head -n 1", "sed 3q", "grep -q a"]
    );
    // A reader that reads to the end of its input.
    for step in [
        r#"zcat a.gz | awk "NR <= 144395" > "$STAGE/cut.fa""#,
        "zcat a.gz | awk '/^>MT/ { past = 1 } !past' > b",
        "zcat a.gz | awk '$1 == id { found = 1 } END { exit !found }'",
        "zcat a.gz | awk '{ n++ } END { if (n < 2) { exit 1 } }'",
        "zcat a.gz | grep . > /dev/null && zcat a.gz | sed 's/1q/2q/' | bgzip > b",
        "printf '%s\\n' \"$cols\" | awk -F, '{ print $2 }' || exit 1",
    ] {
        assert_eq!(early(step), Vec::<String>::new(), "{step}");
    }
    // A pipeline whose status no statement reads: after `set +o pipefail` in its own
    // subshell, and in a command substitution that is a word of a command.
    for step in [
        r#"sequences() { (set +o pipefail; gzip -dc "$1" | head -c 8 | tail -c 4); }"#,
        r#"[ "$(sed -n 2p a.gff3 | head -c 1 | od -A n -t x1)" = 09 ] || exit 1"#,
        "local first=$(zcat a.gz | head -n 1)",
        "echo first=$(zcat a.gz | head -n 1)",
        "redo first=$(zcat a.gz | head -n 1)",
        "export rows=($(zcat a.gz | head -n 2))",
        "echo \"\u{e9}=$(zcat a.gz | head -n 1)\"",
    ] {
        assert_eq!(early(step), Vec::<String>::new(), "{step}");
    }
    // `set +o pipefail` in the body of a function is not the step's own shell's, and
    // `set -o pipefail` brings the option back.
    assert_eq!(
        early(concat!(
            "row() { set +o pipefail; zcat \"$1\" | head -n 1; }\n",
            "(set +o pipefail; zcat a.gz | head -n 2; set -o pipefail; zcat a.gz | head -n 3)\n",
            "zcat a.gz | head -n 4",
        )),
        ["head -n 1", "head -n 3", "head -n 4"]
    );
    // In a substitution that is a word of a command, a pipeline whose status `&&`, `||`
    // or a condition reads counts with every reader it holds.
    assert_eq!(
        early(r#"[ -n "$(zcat a.gz | grep -q x && echo found)" ]"#),
        ["grep -q x"]
    );
    assert_eq!(
        early(concat!(
            r#"echo "$(if zcat a.gz | grep -q x; then echo found; fi)" "#,
            r#""$(zcat a.gz | head -n 1 | cat || echo none)""#,
        )),
        ["grep -q x", "head -n 1"]
    );
    // Text that is no command of the step: a quoted pipe, a comment and the body of a
    // here-document.
    for step in [
        r#"awk -F'|' '{ print $1 }' a > b; echo "a | head" 'b | head' $'it\'s | head'"#,
        "sort a > b # then | head -n 1\ncat b",
        "python3 - <<'PY'\nrows = open('a').read() | head\n# it's | head\nPY\ncat b",
        "cat <<-GFF > a.gff3\n\t1\t.\tgene | head\n\tGFF\ncat a.gff3",
    ] {
        assert_eq!(early(step), Vec::<String>::new(), "{step}");
    }
    assert_eq!(
        early("cat <<'GFF' > a.gff3\nit's a line\nGFF\nzcat a.gz | head -n 1 <<< \"$x\""),
        ["head -n 1 <<< \"$x\""]
    );
    // The command after a here-document whose delimiter stands apart from `<<-` or ends
    // at a redirection, after the shift of an arithmetic expansion and after a comment
    // that follows `;`.
    for step in [
        "cat <<- GFF\n\tit's | head\n\tGFF\nzcat a.gz | head -n 1",
        "cat <<GFF>a.gff3\nit's | head\nGFF\nzcat a.gz | head -n 1",
        "n=$((1 << 4));# it's a shift\nzcat a.gz | head -n 1",
    ] {
        assert_eq!(early(step), ["head -n 1"], "{step}");
    }
}

/// No `prepare` step of a committed reproduction reads the status of a pipeline whose
/// reader can leave before its writer's last write. Such a step fails by timing, with
/// exit code 141 and no line of its own, and the entry then reads `error` in a run that
/// changed no row of either engine; `--check` does not show it, since it runs no step.
#[test]
fn no_committed_prepare_step_reads_the_status_of_a_writer_its_reader_can_end() {
    let directory = format!("{}/scripts/divergences/reproductions", repo_root());
    let repros = load_reproductions(None, &directory).expect("the committed records");
    let mut piped = 0;
    let mut early = Vec::new();
    for repro in &repros {
        let Some(Value::List(steps)) = repro.data.get("prepare") else {
            continue;
        };
        for step in steps.iter().filter_map(Value::as_str) {
            piped += usize::from(step.contains('|'));
            early.extend(
                readers_that_leave_early(step)
                    .into_iter()
                    .map(|reader| format!("{}: {reader}", repro.name)),
            );
        }
    }
    assert!(piped > 0, "no prepare step holds a pipe sign");
    assert_eq!(early, Vec::<String>::new());
}
