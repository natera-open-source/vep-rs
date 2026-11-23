// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The first word of the command line: `score` and `compare` each take the words after it,
//! and any other first word is refused with the usage line.

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-plugin-pairs");
const USAGE: &str = "usage: vep-plugin-pairs {score,compare} ...\n";

fn run(words: &[&str]) -> Output {
    Command::new(BIN).args(words).output().unwrap()
}

fn text(bytes: &[u8]) -> &str {
    std::str::from_utf8(bytes).unwrap()
}

#[test]
fn first_word_names_the_command_and_the_rest_is_its_command_line() {
    for (command, usage) in [
        ("score", "usage: vep-plugin-pairs score "),
        ("compare", "usage: vep-plugin-pairs compare "),
    ] {
        let help = run(&[command, "--help"]);
        assert_eq!(help.status.code(), Some(0), "{command}");
        assert!(text(&help.stdout).starts_with(usage), "{command}");
        let refused = run(&[command, "--no-such-option"]);
        assert_eq!(refused.status.code(), Some(2), "{command}");
        let error = text(&refused.stderr);
        assert!(error.starts_with(usage), "{error}");
        assert!(
            error.contains(&format!("vep-plugin-pairs {command}: error: ")),
            "{error}"
        );
    }
}

#[test]
fn help_lists_both_commands() {
    for word in ["-h", "--help"] {
        let help = run(&[word]);
        assert_eq!(help.status.code(), Some(0));
        let printed = text(&help.stdout);
        assert!(printed.starts_with(USAGE), "{printed}");
        assert!(printed.contains("\n  score  ") && printed.contains("\n  compare  "));
        assert!(help.stderr.is_empty());
    }
}

#[test]
fn command_line_without_a_command_is_refused_with_the_usage_line() {
    let none = run(&[]);
    assert_eq!(none.status.code(), Some(2));
    assert_eq!(
        text(&none.stderr),
        format!(
            "{USAGE}ERROR: [vep-plugin-pairs] the first word names the command: score or compare\n"
        )
    );
    assert!(none.stdout.is_empty());
    // An option where the command stands names no command, whatever follows it.
    let option = run(&["--reference", "reference.txt", "score"]);
    assert_eq!(option.status.code(), Some(2));
    assert_eq!(
        text(&option.stderr),
        format!("{USAGE}ERROR: [vep-plugin-pairs] no command --reference: the first word is score or compare\n")
    );
    assert!(option.stdout.is_empty());
}

/// The writing end of a pipe whose reading end is closed: a write to it fails.
fn pipe_without_a_reader() -> std::io::PipeWriter {
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    writer
}

#[test]
fn a_standard_stream_without_a_reader_leaves_the_exit_status_as_it_is() {
    for word in ["-h", "--help"] {
        let help = Command::new(BIN)
            .arg(word)
            .stdout(pipe_without_a_reader())
            .output()
            .unwrap();
        assert_eq!(help.status.code(), Some(0), "{word}");
        assert!(help.stderr.is_empty(), "{}", text(&help.stderr));
    }
    for words in [&[][..], &["annotate"]] {
        let refused = Command::new(BIN)
            .args(words)
            .stderr(pipe_without_a_reader())
            .output()
            .unwrap();
        assert_eq!(refused.status.code(), Some(2), "{words:?}");
        assert!(refused.stdout.is_empty());
    }
}
