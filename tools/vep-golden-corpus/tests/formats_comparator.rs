// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field --spliceai-table` on the made-up run
//! `tests/fixtures/spliceai/record_of_two_entries/` with the format comparator named in each
//! way a run names one: the build the option names, the one the variable names, the release
//! build of the tree the executable stands in. A run that finds no build is refused with
//! nothing written, the command runs once a run over every field that reaches the class,
//! and a run whose command leaves no answers ends on one line. The stand-ins for the
//! format comparator are shell scripts the tests write.

mod classify_support;
mod common;
mod recorded_runs;
mod score_files;

use classify_support::counted;
use common::{
    binary_in_tree, clean, ended, formats_comparator, refused, started, Run, BIN,
    FORMATS_COMPARATOR, ROOT,
};
use recorded_runs::{Call, Recorded};
use score_files::{score_file, DELETION, SUBSTITUTION};
use std::fs;
use std::os::unix::fs::{symlink, PermissionsExt};
use std::path::Path;
use std::process::Command;
use vep_golden_corpus::TempDir;

const FIELDS: [&str; 2] = ["SpliceAI_pred", "SpliceAI_cutoff"];
const OPTION: &str = "--formats-comparator";

/// How a refusal of a run that finds no build ends.
const HOW_A_BUILD_IS_HAD: &str = ": build the tool with 'cargo build --release --locked \
     --manifest-path tools/vep-formats-compare/Cargo.toml' or name a built one with \
     --formats-comparator or in VEP_FORMATS_COMPARATOR";

/// The run of two entries with its score file of one record line.
fn scored_run() -> Recorded {
    let run = Recorded::new("spliceai/record_of_two_entries");
    score_file(&run.root().join("joined"), DELETION, true);
    run
}

/// The arguments of `classify` on the run's corpus over `fields` and `tables`, with no
/// build of the format comparator named.
fn arguments(run: &Recorded, fields: &[&str], tables: &[&str]) -> Vec<String> {
    let mut arguments = run.arguments(&Call {
        directory: "",
        outputs: ("vep.txt", "rs.txt"),
        fields,
        tables,
        summary: (0, "{}"),
        expected: "",
    });
    if !tables.is_empty() {
        arguments.truncate(arguments.len() - 2);
    }
    arguments
}

/// The arguments of the call the fixture's second manifest is written by.
fn asked(run: &Recorded) -> Vec<String> {
    arguments(run, &FIELDS, &["joined/indel.synthetic.vcf.gz"])
}

/// A command of `program` from the package root with no build named in its environment.
fn command(program: &str, arguments: &[String]) -> Command {
    let mut command = Command::new(program);
    command
        .current_dir(ROOT)
        .env_remove("COLUMNS")
        .env_remove(FORMATS_COMPARATOR)
        .args(arguments);
    command
}

/// The line a run ends on that writes no key divergence and `count` field divergences.
fn line(count: usize, summary: &str) -> String {
    let keys = counted((0, "{}"), (0, "{}"), "[]");
    let keys = keys.strip_suffix('\n').unwrap();
    format!("{keys}; field divergences {count} ({summary})\n")
}

/// Holds a run to the fixture's second call: the three fields under the class, each row's
/// value under the names of the other formats behind them, and the manifest's bytes.
fn holds(run: &Recorded, ran: Run, how: &str) {
    assert_eq!(
        ran,
        clean(0, &line(15, "{'spliceai_multi_entry_cell': 15}")),
        "{how}"
    );
    let expected = common::fixture("spliceai/record_of_two_entries/expected.2.json");
    assert_eq!(
        run.manifest(""),
        fs::read_to_string(expected).unwrap(),
        "{how}"
    );
}

/// Writes the shell script `body` as the executable `name` under `directory`; its path.
fn script(directory: &Path, name: &str, body: &str) -> String {
    let path = directory.join(name);
    fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
    fs::set_permissions(&path, fs::Permissions::from_mode(0o755)).unwrap();
    started(&path);
    path.to_str().unwrap().to_string()
}

/// A tree of its own under the package's test directory: the tool's binary inside it, and
/// the format tool's directory with, where `build` names one, a link to that build as its
/// release build. The scratch directory, the binary and the release build's path.
fn tree(build: Option<&str>) -> (TempDir, String, String) {
    let scratch = TempDir::create(Path::new(env!("CARGO_TARGET_TMPDIR")), "tree-").unwrap();
    let own = binary_in_tree(&scratch.path);
    let release = scratch
        .path
        .join("tools/vep-formats-compare/target/release/vep-formats-compare");
    fs::create_dir_all(release.parent().unwrap()).unwrap();
    if let Some(build) = build {
        symlink(build, &release).unwrap();
    }
    let release = release.to_str().unwrap().to_string();
    (scratch, own, release)
}

#[test]
fn the_build_asked_is_the_options_then_the_variables_then_the_trees_own() {
    let run = scored_run();
    let build = formats_comparator();
    let asked = asked(&run);
    let with_option = |name: &str| {
        let mut arguments = asked.clone();
        arguments.extend([OPTION.to_string(), name.to_string()]);
        arguments
    };
    holds(
        &run,
        ended(&mut command(BIN, &with_option(&build))),
        "the option",
    );
    holds(
        &run,
        ended(command(BIN, &asked).env(FORMATS_COMPARATOR, &build)),
        "the variable",
    );
    // The option is read before the variable, and an empty value names none.
    let absent = run.path("absent");
    holds(
        &run,
        ended(command(BIN, &with_option(&build)).env(FORMATS_COMPARATOR, &absent)),
        "the option before the variable",
    );
    holds(
        &run,
        ended(command(BIN, &with_option("")).env(FORMATS_COMPARATOR, &build)),
        "an empty option",
    );
    // A bare name is a command of PATH.
    let commands = run.root().join("commands");
    fs::create_dir(&commands).unwrap();
    symlink(&build, commands.join("formats-comparator-of-this-test")).unwrap();
    holds(
        &run,
        ended(command(BIN, &with_option("formats-comparator-of-this-test")).env("PATH", &commands)),
        "a command of PATH",
    );
    // With no name, the release build of the tree the executable stands in, and with a
    // name that build is not asked.
    let (_scratch, own, _) = tree(Some(&build));
    holds(&run, ended(&mut command(&own, &asked)), "the tree's build");
    holds(
        &run,
        ended(command(&own, &asked).env(FORMATS_COMPARATOR, "")),
        "an empty variable",
    );
    let (_scratch, unbuilt, _) = tree(None);
    holds(
        &run,
        ended(command(&unbuilt, &asked).env(FORMATS_COMPARATOR, &build)),
        "the variable before the tree",
    );
}

#[test]
fn a_run_that_finds_no_build_is_refused_with_nothing_written() {
    let run = scored_run();
    let before = run.manifest("");
    let asked = asked(&run);
    let own = |found: String| {
        refused(
            1,
            &format!("ERROR: [vep-golden-corpus] {found}{HOW_A_BUILD_IS_HAD}"),
        )
    };
    let not_executable = |name: &str, by: &str| {
        own(format!(
            "the format comparator {name} of {by} is not an executable"
        ))
    };
    let (absent, plain) = (run.path("absent"), run.path("vep.txt"));
    let directory = run.path("joined");
    for name in [&absent, &plain, &directory, "no-such-command-of-this-test"] {
        let mut arguments = asked.clone();
        arguments.extend([OPTION.to_string(), name.to_string()]);
        assert_eq!(
            ended(&mut command(BIN, &arguments)),
            not_executable(name, OPTION),
            "{name}"
        );
        assert_eq!(
            ended(command(BIN, &asked).env(FORMATS_COMPARATOR, name)),
            not_executable(name, FORMATS_COMPARATOR),
            "{name}"
        );
        assert_eq!(run.manifest(""), before, "{name}");
    }
    // A tree without the release build, whatever else stands in it.
    let (scratch, unbuilt, release) = tree(None);
    let unbuilt_tree = own(format!(
        "the format comparator {release} is not an executable"
    ));
    assert_eq!(ended(&mut command(&unbuilt, &asked)), unbuilt_tree);
    let debug = scratch
        .path
        .join("tools/vep-formats-compare/target/debug/vep-formats-compare");
    fs::create_dir_all(debug.parent().unwrap()).unwrap();
    symlink(formats_comparator(), debug).unwrap();
    assert_eq!(ended(&mut command(&unbuilt, &asked)), unbuilt_tree);
    assert_eq!(run.manifest(""), before);
    // A run that names no score file, and one that compares no field, look for no build.
    let unasked = |fields: &[&str], tables: &[&str]| {
        let mut arguments = arguments(&run, fields, tables);
        arguments.extend([OPTION.to_string(), absent.clone()]);
        ended(command(BIN, &arguments).env(FORMATS_COMPARATOR, &absent))
    };
    assert_eq!(
        unasked(&FIELDS, &[]),
        clean(0, &line(3, "{'unexplained_residual': 3}"))
    );
    assert_eq!(
        unasked(&[], &["joined/indel.synthetic.vcf.gz"]),
        clean(0, &counted((0, "{}"), (0, "{}"), "[]"))
    );
}

#[test]
fn the_command_runs_once_a_run_over_every_field_that_reaches_the_class() {
    let run = scored_run();
    score_file(&run.root().join("other"), SUBSTITUTION, true);
    run.write(
        "c/json_cache/transcripts/21/1-1000000.json",
        br#"[{"stable_id": "ENST1", "gene_symbol": "GENEA"}, {"stable_id": "ENST2", "gene_symbol": "GENEB"}]"#,
    );
    // A stand-in that keeps what it is started with and hands on to the build.
    let kept = run.root().join("kept");
    fs::create_dir(&kept).unwrap();
    let stand_in = script(
        run.root(),
        "logging-comparator",
        &format!(
            r#"[ "$1" = --help ] && exit 0
kept='{}'
echo started >> "$kept/runs.txt"
printf '%s\n' "$@" > "$kept/arguments.txt"
for argument in "$@"; do
    case "$argument" in
        --questions=*) cp "${{argument#--questions=}}" "$kept/questions.tsv" ;;
    esac
done
exec '{}' "$@""#,
            kept.display(),
            formats_comparator()
        ),
    );
    let tables = [
        "joined/indel.synthetic.vcf.gz",
        "other/snv.synthetic.vcf.gz",
    ];
    let with_stand_in = |fields: &[&str]| {
        let mut arguments = arguments(&run, fields, &tables);
        arguments.extend([OPTION.to_string(), stand_in.clone()]);
        ended(&mut command(BIN, &arguments))
    };
    let read = |name: &str| fs::read_to_string(kept.join(name)).unwrap();
    holds(&run, with_stand_in(&FIELDS), "a stand-in before the build");
    assert_eq!(read("runs.txt"), "started\n");
    let started_with = read("arguments.txt");
    let started_with: Vec<&str> = started_with.lines().collect();
    assert_eq!(
        started_with[..5],
        [
            "field-class".to_string(),
            "--class=spliceai_multi_entry_cell".to_string(),
            format!("--spliceai-table={}", run.path(tables[0])),
            format!("--spliceai-table={}", run.path(tables[1])),
            format!("--json-cache={}", run.path("c/json_cache")),
        ]
    );
    assert_eq!(started_with.len(), 7);
    let questions = Path::new(started_with[5].strip_prefix("--questions=").unwrap());
    let answers = Path::new(started_with[6].strip_prefix("--answers=").unwrap());
    assert_eq!(questions.parent(), answers.parent());
    assert!(
        !questions.parent().unwrap().exists(),
        "the questions leave with the run"
    );
    // One question a differing field, the keys in their order and the fields in the
    // order of `--field`, each with the twelve fields of both rows.
    let asked = read("questions.tsv");
    let asked: Vec<Vec<&str>> = asked
        .lines()
        .map(|line| line.split('\t').collect())
        .collect();
    let named: Vec<(&str, &str, &str, &str)> = asked
        .iter()
        .map(|cells| (cells[0], cells[1], cells[2], cells[3]))
        .collect();
    assert_eq!(
        named,
        [
            ("21:5001-5002", "-", "ENST1", "SpliceAI_pred"),
            ("21:5001-5002", "-", "ENST2", "SpliceAI_pred"),
            ("21:5001-5002", "-", "ENST2", "SpliceAI_cutoff"),
        ]
    );
    for cells in &asked {
        assert_eq!(cells.len(), 30);
    }
    let (first, last) = (&asked[0], &asked[2]);
    assert_eq!(first[4], first[6], "Ensembl VEP's value is its row's");
    assert_eq!(first[5], first[18], "vep-rs's value is its row's");
    assert_eq!((first[16], first[17]), ("FAIL", "GENEA"));
    assert_eq!((last[4], last[5]), ("-", "PASS"));
    assert_eq!((last[6], last[16], last[17]), ("-", "-", "GENEB"));
    assert_eq!((last[28], last[29]), ("PASS", "GENEB"));
    // A run under score files in which no compared field differs asks no question, and
    // the command runs once all the same: it reads the score files whatever it is asked.
    assert_eq!(with_stand_in(&["SYMBOL"]), clean(0, &line(0, "{}")));
    assert_eq!(read("runs.txt"), "started\nstarted\n");
    assert_eq!(read("questions.tsv"), "");
}

#[test]
fn a_run_whose_command_leaves_no_answers_is_refused_with_nothing_written() {
    let run = scored_run();
    let before = run.manifest("");
    let asked = asked(&run);
    let stand_ins = run.root().join("stand-ins");
    fs::create_dir(&stand_ins).unwrap();
    let ends_on = |name: &str, body: &str, detail: &dyn Fn(&str) -> String| {
        let stand_in = script(&stand_ins, name, body);
        let mut arguments = asked.clone();
        arguments.extend([OPTION.to_string(), stand_in.clone()]);
        let expected = format!("ERROR: [vep-golden-corpus] {}", detail(&stand_in));
        assert_eq!(
            ended(&mut command(BIN, &arguments)),
            refused(1, &expected),
            "{name}"
        );
        assert_eq!(run.manifest(""), before, "{name}");
    };
    let answers = r#"for argument in "$@"; do
    case "$argument" in
        --answers=*) answers="${argument#--answers=}" ;;
    esac
done"#;
    let unanswered = |stand_in: &str| {
        format!(
            "the format comparator {stand_in} did not answer each of the 3 questions of \
             spliceai_multi_entry_cell with 0 or 1"
        )
    };
    ends_on("silent", "exit 0", &|stand_in| {
        format!(
            "the format comparator {stand_in} ended with exit code 0 and left no answers \
             (No such file or directory (os error 2))"
        )
    });
    ends_on(
        "one-answer",
        &format!("{answers}\nprintf '1\\n' > \"$answers\""),
        &unanswered,
    );
    ends_on(
        "another-answer",
        &format!("{answers}\nprintf '1\\n2\\n0\\n' > \"$answers\""),
        &unanswered,
    );
    // A refusal of the command is the run's, in the command's words.
    ends_on(
        "refusing",
        "echo 'a line before' >&2\n\
         echo 'ERROR: [compare_vep_formats] --json-cache: c/json_cache has no transcripts/ directory' >&2\n\
         exit 2",
        &|_| "--json-cache: c/json_cache has no transcripts/ directory".to_string(),
    );
    // Any other ending names the build and how it ended.
    ends_on(
        "ending",
        "echo 'ERROR: [compare_vep_formats] embedded null byte' >&2\nexit 1",
        &|stand_in| {
            format!("the format comparator {stand_in} ended with exit code 1: embedded null byte")
        },
    );
    ends_on(
        "another-command",
        "echo 'usage: another [-h]' >&2\necho 'another: error: unrecognized arguments: field-class' >&2\nexit 2",
        &|stand_in| {
            format!(
                "the format comparator {stand_in} ended with exit code 2: another: error: \
                 unrecognized arguments: field-class"
            )
        },
    );
    ends_on(
        "killed",
        "[ \"$1\" = --help ] && exit 0\nkill -9 $$",
        &|stand_in| {
            format!(
                "the format comparator {stand_in} was ended by a signal and wrote no error line"
            )
        },
    );
    // A program the system cannot start.
    let unstarted = stand_ins.join("unstarted");
    fs::write(&unstarted, "#!/no/such/interpreter-of-this-test\n").unwrap();
    fs::set_permissions(&unstarted, fs::Permissions::from_mode(0o755)).unwrap();
    let unstarted = unstarted.to_str().unwrap();
    let mut arguments = asked.clone();
    arguments.extend([OPTION.to_string(), unstarted.to_string()]);
    assert_eq!(
        ended(&mut command(BIN, &arguments)),
        refused(
            1,
            &format!(
                "ERROR: [vep-golden-corpus] the format comparator {unstarted} did not start \
                 (No such file or directory (os error 2))"
            )
        )
    );
    assert_eq!(run.manifest(""), before);
}
