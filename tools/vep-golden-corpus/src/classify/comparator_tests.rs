// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the format comparator's place: the build the option, the variable and the
//! executable's tree name, the refusal of a name that is no executable, and what a run of
//! the command that left no answers says.

use super::*;
use std::fs;

fn scratch(tag: &str) -> TempDir {
    TempDir::create(&std::env::temp_dir(), &format!("golden_comparator_{tag}_")).unwrap()
}

/// A file of the mode `mode` at `path`, with the directories above it.
fn file(path: &Path, mode: u32) {
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::write(path, "").unwrap();
    fs::set_permissions(path, fs::Permissions::from_mode(mode)).unwrap();
}

/// A tree under `root` with the tool's directory, an executable of this tool inside it and,
/// where `built`, the tool's release build: the two paths.
fn tree(root: &Path, built: bool) -> (PathBuf, PathBuf) {
    let executable = root.join("tools/vep-golden-corpus/target/debug/vep-golden-corpus");
    file(&executable, 0o755);
    let build = root.join(TOOL_DIRECTORY).join(RELEASE_BUILD);
    fs::create_dir_all(root.join(TOOL_DIRECTORY)).unwrap();
    if built {
        file(&build, 0o755);
    }
    (executable, build)
}

fn tail() -> String {
    format!(
        ": build the tool with '{BUILD_COMMAND}' or name a built one with {OPTION} or in {VARIABLE}"
    )
}

fn refused(named: &Named<'_>) -> String {
    let failure = locate(named).unwrap_err();
    assert_eq!(failure.code, 1);
    let line = failure.line.strip_prefix("ERROR: [vep-golden-corpus] ");
    let found = line.unwrap().strip_suffix(&tail());
    found
        .expect("the line ends on how a build is had")
        .to_string()
}

#[test]
fn the_build_is_the_options_then_the_variables_then_the_trees() {
    let scratch = scratch("order");
    let (executable, build) = tree(&scratch.path, true);
    let (by_option, by_variable) = (scratch.path.join("a/one"), scratch.path.join("b/two"));
    file(&by_option, 0o700);
    file(&by_variable, 0o755);
    let all = Named {
        option: Some(by_option.as_os_str()),
        variable: Some(by_variable.as_os_str()),
        executable: Some(&executable),
        search: None,
    };
    assert_eq!(locate(&all), Ok(by_option.clone()));
    let unnamed = Named {
        option: None,
        ..all
    };
    assert_eq!(locate(&unnamed), Ok(by_variable.clone()));
    // An empty value names none.
    let empty = Named {
        option: Some(OsStr::new("")),
        variable: Some(OsStr::new("")),
        ..all
    };
    assert_eq!(locate(&empty), Ok(build.clone()));
    let own_tree = Named {
        executable: Some(&executable),
        ..Named::default()
    };
    assert_eq!(locate(&own_tree), Ok(build));
    assert_eq!(
        BUILD_COMMAND.split(' ').next_back(),
        Some("tools/vep-formats-compare/Cargo.toml")
    );
}

#[test]
fn a_bare_name_is_a_command_of_the_search_path() {
    let scratch = scratch("search");
    let (first, second) = (scratch.path.join("first"), scratch.path.join("second"));
    file(&first.join("vep-formats-compare"), 0o644);
    file(&second.join("vep-formats-compare"), 0o755);
    fs::create_dir_all(first.join("directory")).unwrap();
    file(&second.join("directory"), 0o755);
    let search = std::env::join_paths([&first, &second]).unwrap();
    let named = |name: &'static str| Named {
        variable: Some(OsStr::new(name)),
        search: Some(&search),
        ..Named::default()
    };
    // The first executable file of the name, past a file that is none and a directory.
    assert_eq!(
        locate(&named("vep-formats-compare")),
        Ok(second.join("vep-formats-compare"))
    );
    assert_eq!(locate(&named("directory")), Ok(second.join("directory")));
    assert_eq!(
        refused(&named("absent")),
        "the format comparator absent of VEP_FORMATS_COMPARATOR is not an executable"
    );
    let unsearched = Named {
        search: None,
        ..named("vep-formats-compare")
    };
    assert_eq!(
        refused(&unsearched),
        "the format comparator vep-formats-compare of VEP_FORMATS_COMPARATOR is not an executable"
    );
}

#[test]
fn a_name_that_is_no_executable_is_refused_and_nothing_stands_in_for_it() {
    let scratch = scratch("refused");
    let (executable, build) = tree(&scratch.path, true);
    let plain = scratch.path.join("plain");
    file(&plain, 0o644);
    let absent = scratch.path.join("absent");
    for (path, by) in [
        (&plain, OPTION),
        (&absent, OPTION),
        (&scratch.path, OPTION),
        (&plain, VARIABLE),
    ] {
        let name = Some(path.as_os_str());
        let named = Named {
            option: name.filter(|_| by == OPTION),
            variable: if by == VARIABLE {
                name
            } else {
                Some(build.as_os_str())
            },
            executable: Some(&executable),
            search: None,
        };
        assert_eq!(
            refused(&named),
            format!(
                "the format comparator {} of {by} is not an executable",
                path.display()
            )
        );
    }
}

#[test]
fn the_trees_build_is_the_release_build_under_the_nearest_directory_that_holds_the_tool() {
    let scratch = scratch("tree");
    let (executable, build) = tree(&scratch.path.join("outer"), false);
    let own_tree = |executable: &'_ Path| {
        let named = Named {
            executable: Some(executable),
            ..Named::default()
        };
        locate(&named).map_err(|_| refused(&named))
    };
    assert_eq!(
        own_tree(&executable),
        Err(format!(
            "the format comparator {} is not an executable",
            build.display()
        ))
    );
    file(&build, 0o644);
    assert!(own_tree(&executable).is_err());
    file(&build, 0o755);
    assert_eq!(own_tree(&executable), Ok(build.clone()));
    // A tree inside a tree is the executable's own, built or not.
    let (inner, inner_build) = tree(&scratch.path.join("outer/nested/inner"), false);
    assert_eq!(
        own_tree(&inner),
        Err(format!(
            "the format comparator {} is not an executable",
            inner_build.display()
        ))
    );
    // An executable that stands in no tree, and a run that knows no executable.
    let outside = scratch.path.join("bin/vep-golden-corpus");
    file(&outside, 0o755);
    let no_tree = format!(
        "no format comparator is named and this executable stands in no tree that holds {TOOL_DIRECTORY}"
    );
    assert_eq!(own_tree(&outside), Err(no_tree.clone()));
    assert_eq!(refused(&Named::default()), no_tree);
}

#[test]
fn an_option_is_one_argument_and_an_ended_run_says_how_it_ended() {
    assert_eq!(
        option("--questions", Path::new("scratch/questions.tsv")),
        OsString::from("--questions=scratch/questions.tsv")
    );
    assert_eq!(option("--class", "-x"), OsString::from("--class=-x"));
    let ended = |code: Option<i32>, error: &str| Ended {
        code,
        error: error.to_string(),
    };
    let refused = ended(
        Some(2),
        "ERROR: [compare_vep_formats] cannot read q.tsv: gone\n\n",
    );
    assert_eq!(
        refused.last_line(),
        "ERROR: [compare_vep_formats] cannot read q.tsv: gone"
    );
    assert_eq!(refused.detail(), Some("cannot read q.tsv: gone"));
    assert_eq!(refused.how(), "ended with exit code 2");
    let usage = ended(Some(2), "usage: x\nx: error: unrecognized arguments: y\n");
    assert_eq!(usage.last_line(), "x: error: unrecognized arguments: y");
    assert_eq!(usage.detail(), None);
    let silent = ended(None, "");
    assert_eq!((silent.last_line(), silent.detail()), ("", None));
    assert_eq!(silent.how(), "was ended by a signal");
    // The failure of a run no caller words: the build, how it ended and its last line.
    let comparator = Comparator {
        program: PathBuf::from("bin/vep-formats-compare"),
    };
    let line = |ended: &Ended| {
        let failure = comparator.ended(ended);
        assert_eq!(failure.code, 1);
        failure.line
    };
    let named = "ERROR: [vep-golden-corpus] the format comparator bin/vep-formats-compare";
    assert_eq!(
        line(&refused),
        format!("{named} ended with exit code 2: cannot read q.tsv: gone")
    );
    assert_eq!(
        line(&usage),
        format!("{named} ended with exit code 2: x: error: unrecognized arguments: y")
    );
    assert_eq!(
        line(&silent),
        format!("{named} was ended by a signal and wrote no error line")
    );
}

#[test]
fn the_answers_of_a_run_are_one_zero_or_one_a_question() {
    let comparator = Comparator {
        program: PathBuf::from("bin/vep-formats-compare"),
    };
    let accepted = |questions: usize, answered: &str| {
        let accepted = comparator.accepted("a_class", questions, answered);
        accepted.map_err(|failure| failure.line)
    };
    assert_eq!(accepted(3, "1\n0\n1\n"), Ok(vec![true, false, true]));
    assert_eq!(accepted(0, ""), Ok(Vec::new()));
    let refused = |questions: usize| {
        Err(format!(
            "ERROR: [vep-golden-corpus] the format comparator bin/vep-formats-compare did not \
             answer each of the {questions} questions of a_class with 0 or 1"
        ))
    };
    assert_eq!(accepted(2, "1\n"), refused(2));
    assert_eq!(accepted(1, "1\n0\n"), refused(1));
    assert_eq!(accepted(2, "1\nyes\n"), refused(2));
}
