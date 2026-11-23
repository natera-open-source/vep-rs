// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pieces of `generated`: the tree a command runs from, the files a scratch directory
//! holds, the last line of a failed command, the command's ends and the lines of `--check`.

use super::*;
use crate::json::loads;
use std::fs;
use std::os::unix::fs::{symlink, PermissionsExt};

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-generated-test-").unwrap()
}

#[test]
fn a_command_runs_from_the_nearest_tree_above_the_executable() {
    let scratch = scratch();
    let outer = scratch.path.join("outer");
    let inner = outer.join("vendor/inner");
    for tree in [&outer, &inner] {
        fs::create_dir_all(tree.join(TOOL_DIRECTORY).join("target/release")).unwrap();
    }
    let build = |tree: &Path| {
        tree.join(TOOL_DIRECTORY)
            .join("target/release/vep-golden-corpus")
    };
    assert_eq!(tree_of(Some(&build(&outer))), Ok(outer.as_path()));
    assert_eq!(tree_of(Some(&build(&inner))), Ok(inner.as_path()));
    // A build that stands beside the tool's directory is in the tree all the same.
    assert_eq!(
        tree_of(Some(&outer.join("bin/vep-golden-corpus"))),
        Ok(outer.as_path())
    );
    let refusal = Failure::own(
        "this executable stands in no tree that holds tools/vep-golden-corpus, so there is no repository root to run a generated command from",
    );
    let elsewhere = scratch.path.join("elsewhere/bin/vep-golden-corpus");
    assert_eq!(tree_of(Some(&elsewhere)), Err(refusal.clone()));
    assert_eq!(tree_of(None), Err(refusal));
}

#[test]
fn the_files_of_a_directory_are_its_regular_files_and_the_links_to_one() {
    let scratch = scratch();
    let root = scratch.path.join("root");
    fs::create_dir_all(root.join("set/deep")).unwrap();
    fs::create_dir_all(root.join("empty")).unwrap();
    fs::write(root.join("top.txt"), "t").unwrap();
    fs::write(root.join(".hidden"), "h").unwrap();
    fs::write(root.join("set/b"), "b").unwrap();
    fs::write(root.join("set/deep/a"), "a").unwrap();
    // A link to a file is a file under its own name; a link to a directory is not entered
    // and a link to nothing is no file.
    symlink("set/b", root.join("link-to-file")).unwrap();
    symlink("set", root.join("link-to-directory")).unwrap();
    symlink("absent", root.join("link-to-nothing")).unwrap();
    let mut found = Vec::new();
    files_under(&root, "", &mut found).unwrap();
    let mut names: Vec<&str> = found.iter().map(|(below, _)| below.as_str()).collect();
    names.sort_unstable();
    assert_eq!(
        names,
        [".hidden", "link-to-file", "set/b", "set/deep/a", "top.txt"]
    );
    for (below, path) in &found {
        assert_eq!(path, &root.join(below));
    }
    // A directory that cannot be listed is passed over, as is one that is absent.
    let mut none = Vec::new();
    files_under(&root.join("absent"), "", &mut none).unwrap();
    files_under(&root.join("top.txt"), "", &mut none).unwrap();
    assert!(none.is_empty());
}

#[test]
fn the_files_of_a_directory_come_in_the_order_rglob_yields_them() {
    let scratch = scratch();
    // A host lists a directory in an order of its own, so the tree is made both ways.
    for (name, made) in [
        ("one_way", ["top", "set/b", "set/deep/a", "side/c"]),
        ("the_other", ["side/c", "set/deep/a", "set/b", "top"]),
    ] {
        let root = scratch.path.join(name);
        for file in made {
            let path = root.join(file);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, "").unwrap();
        }
        let mut found = Vec::new();
        files_under(&root, "", &mut found).unwrap();
        let mut names: Vec<&str> = found.iter().map(|(below, _)| below.as_str()).collect();
        // The file of the directory comes first, then those of the two directories it
        // holds, in the host's order, then the file of the directory one of them holds.
        names[1..3].sort_unstable();
        assert_eq!(names, ["top", "set/b", "side/c", "set/deep/a"], "{name}");
    }
}

/// Whether this user looks into a directory of any mode, as the superuser does.
fn looks_into_any_directory(scratch: &TempDir) -> bool {
    let sealed = scratch.path.join("sealed");
    fs::create_dir(&sealed).unwrap();
    fs::write(sealed.join("file"), "").unwrap();
    fs::set_permissions(&sealed, fs::Permissions::from_mode(0o000)).unwrap();
    let looked = fs::metadata(sealed.join("file")).is_ok();
    fs::set_permissions(&sealed, fs::Permissions::from_mode(0o700)).unwrap();
    looked
}

#[test]
fn what_a_command_seals_leaves_with_the_scratch_directory() {
    let scratch = scratch();
    // The script notes the directory it writes under in the file of its second argument
    // and gives the directory `set` the mode of its third.
    let body = "printf '%s' \"$1\" > \"$2\"\nmkdir \"$1/set\"\n: > \"$1/set/a\"\n: > \"$1/top\"\nchmod \"$3\" \"$1/set\"\n";
    let (corpus, tree, script) = stage(&scratch, "seal.sh", body);
    let noted = scratch.path.join("noted");
    let noted = noted.to_str().unwrap();
    let run = |mode: &str| {
        let command = words(&["/bin/sh", &script, "{generated}", noted, mode]);
        let outcome = generated_files(&command, &corpus, &tree);
        let written = fs::read_to_string(noted).unwrap();
        assert!(!Path::new(&written).exists(), "{written} is left");
        (outcome, written)
    };
    // The superuser looks into a directory of any mode, so the two ends are any other
    // user's.
    let superuser = looks_into_any_directory(&scratch);
    // A directory that is listed and not searched names its files and denies a look at them.
    let (outcome, written) = run("444");
    if !superuser {
        assert_eq!(
            outcome,
            Err(failure(format!(
                "`/bin/sh {script} {{generated}} {noted} 444` wrote a file that cannot be read ([Errno 13] Permission denied: '{written}/set/a')"
            )))
        );
    }
    // A directory that cannot be listed is passed over.
    let (outcome, _) = run("000");
    if !superuser {
        let digests = outcome.unwrap();
        let names: Vec<&str> = digests.keys().map(String::as_str).collect();
        assert_eq!(names, ["top"]);
    }
}

#[test]
fn the_last_line_of_a_failed_command_is_the_last_that_holds_anything() {
    let cases: [(&[u8], &str); 9] = [
        (b"", "it wrote nothing on standard error"),
        (b" \n\t\r\n", "it wrote nothing on standard error"),
        (b"one", "one"),
        (b"first line\nno set written\n", "no set written"),
        (b"first\r\nsecond\r\n\r\n", "second"),
        (b"first\rsecond", "second"),
        (b"  kept  inside  \n", "kept  inside"),
        (b"a\x0cb\x1cc", "c"),
        (
            b"bytes \xff\xfe outside\n",
            "bytes \u{fffd}\u{fffd} outside",
        ),
    ];
    for (error, line) in cases {
        assert_eq!(last_error_line(error), line, "{error:?}");
    }
}

/// A corpus directory and a tree under `scratch`, and the script `name` of `body`.
fn stage(scratch: &TempDir, name: &str, body: &str) -> (String, PathBuf, String) {
    let corpus = scratch.path.join("corpus");
    let tree = scratch.path.join("tree");
    fs::create_dir_all(&corpus).unwrap();
    fs::create_dir_all(&tree).unwrap();
    let script = scratch.path.join(name);
    fs::write(&script, body).unwrap();
    let text = |path: &Path| path.to_str().unwrap().to_string();
    (text(&corpus), tree, text(&script))
}

fn words(words: &[&str]) -> Vec<String> {
    words.iter().map(|word| word.to_string()).collect()
}

#[test]
fn a_command_runs_from_the_tree_with_its_two_directories_named() {
    let scratch = scratch();
    let body = "printf '%s' \"$(pwd -P)\" > \"$1/cwd\"\nprintf '%s' \"$2\" > \"$1/corpus\"\ncat > \"$1/input\"\nprintf '%s' \"$3\" > \"$1/both\"\n";
    let (corpus, tree, script) = stage(&scratch, "write.sh", body);
    let linked = scratch.path.join("linked");
    symlink(&corpus, &linked).unwrap();
    let command = words(&[
        "/bin/sh",
        &script,
        "{generated}",
        "{corpus}",
        "{corpus}:{generated}:{corpus}",
    ]);
    let digests = generated_files(&command, linked.to_str().unwrap(), &tree).unwrap();
    let resolved = fs::canonicalize(&corpus).unwrap();
    let resolved = resolved.to_str().unwrap();
    let names: Vec<&str> = digests.keys().map(String::as_str).collect();
    assert_eq!(names, ["both", "corpus", "cwd", "input"]);
    // The corpus is named with its links resolved, the command reads nothing on standard
    // input, and a word carries each directory wherever it names one.
    assert_eq!(digests["corpus"], sha256_hex(resolved.as_bytes()));
    let root = fs::canonicalize(&tree).unwrap();
    assert_eq!(
        digests["cwd"],
        sha256_hex(root.to_str().unwrap().as_bytes())
    );
    assert_eq!(digests["input"], sha256_hex(b""));
    assert_ne!(digests["both"], digests["corpus"]);
}

#[test]
fn a_command_that_ends_otherwise_is_refused_by_name() {
    let scratch = scratch();
    let (corpus, tree, script) = stage(&scratch, "end.sh", "eval \"$2\"\n");
    let run = |ending: &str| {
        let command = words(&["/bin/sh", &script, "{generated}", ending]);
        generated_files(&command, &corpus, &tree)
    };
    let named = |ending: &str, how: &str| {
        Err(failure(format!(
            "`/bin/sh {script} {{generated}} {ending}` {how}"
        )))
    };
    assert_eq!(
        run("exit 0"),
        named("exit 0", "wrote no file under {generated}")
    );
    assert_eq!(
        run("mkdir \"$1/empty\""),
        named("mkdir \"$1/empty\"", "wrote no file under {generated}")
    );
    assert_eq!(
        run("exit 7"),
        named(
            "exit 7",
            "failed with status 7: it wrote nothing on standard error"
        )
    );
    assert_eq!(
        run("echo kept >&2; echo last >&2; exit 1"),
        named(
            "echo kept >&2; echo last >&2; exit 1",
            "failed with status 1: last"
        )
    );
    // A run a signal ended has the signal's number below zero for its status.
    assert_eq!(
        run("kill -9 $$"),
        named(
            "kill -9 $$",
            "failed with status -9: it wrote nothing on standard error"
        )
    );
    let unstarted = |program: &str, reason: &str| {
        let command = words(&[program, "{generated}"]);
        assert_eq!(
            generated_files(&command, &corpus, &tree),
            Err(failure(format!(
                "`{program} {{generated}}` does not start ({reason})"
            ))),
            "{program}"
        );
    };
    let absent = scratch.path.join("absent");
    let absent = absent.to_str().unwrap();
    unstarted(
        absent,
        &format!("[Errno 2] No such file or directory: '{absent}'"),
    );
    unstarted(
        "no-such-program-on-the-path",
        "[Errno 2] No such file or directory: 'no-such-program-on-the-path'",
    );
    // A file without the mode of a program, and a directory, are denied.
    fs::set_permissions(&script, fs::Permissions::from_mode(0o644)).unwrap();
    unstarted(
        &script,
        &format!("[Errno 13] Permission denied: '{script}'"),
    );
    unstarted(
        &corpus,
        &format!("[Errno 13] Permission denied: '{corpus}'"),
    );
    unstarted("a\0b", "embedded null byte");
    assert_eq!(
        generated_files(&words(&["/bin/sh", "a\0b"]), &corpus, &tree),
        Err(failure(
            "`/bin/sh a\0b` does not start (embedded null byte)"
        ))
    );
}

#[test]
fn a_check_names_each_path_that_differs_in_path_order() {
    let recorded =
        loads(r#"{"set/b": "bb", "set/a": "aa", "gone/z": "zz", "set/moved": "mm"}"#).unwrap();
    let recorded = recorded.as_object().unwrap();
    let written: BTreeMap<String, String> = [
        ("set/a", "aa"),
        ("set/b", "other"),
        ("set/moved", "mm"),
        ("added", "dd"),
    ]
    .into_iter()
    .map(|(path, digest)| (path.to_string(), digest.to_string()))
    .collect();
    assert_eq!(
        differences(recorded, &written),
        [
            "generated file added is written and not recorded",
            "generated file gone/z is recorded and not written",
            "generated file set/b hashes as other, the manifest records bb",
        ]
    );
    let same: BTreeMap<String, String> = recorded
        .iter()
        .map(|(path, digest)| (path.to_string(), digest.as_str().unwrap().to_string()))
        .collect();
    assert_eq!(differences(recorded, &same), Vec::<String>::new());
    // A recorded path spelt with a `.` or an empty segment is the file a run writes without
    // them, named once.
    let spelt = loads(r#"{"./set/a": "aa", "set//b": "bb", "set/./c": "cc"}"#).unwrap();
    let spelt = spelt.as_object().unwrap();
    let written: BTreeMap<String, String> = [("set/a", "aa"), ("set/b", "other")]
        .into_iter()
        .map(|(path, digest)| (path.to_string(), digest.to_string()))
        .collect();
    assert_eq!(
        differences(spelt, &written),
        [
            "generated file set/b hashes as other, the manifest records bb",
            "generated file set/c is recorded and not written",
        ]
    );
}

#[test]
fn a_check_takes_a_directory_for_no_file_and_holds_a_file_to_each_of_its_spellings() {
    let written: BTreeMap<String, String> = [("set/a".to_string(), "aa".to_string())].into();
    let lines = |recorded: &str| {
        let recorded = loads(recorded).unwrap();
        differences(recorded.as_object().unwrap(), &written)
    };
    // A path that ends on a `.` or an empty segment names a directory, which is no file a
    // run writes: it is recorded and not written under its own spelling.
    let unrecorded = "generated file set/a is written and not recorded";
    assert_eq!(
        lines(r#"{"set/a/": "aa"}"#),
        [
            unrecorded,
            "generated file set/a/ is recorded and not written"
        ]
    );
    assert_eq!(
        lines(r#"{"./set//a/.": "aa"}"#),
        [
            "generated file ./set//a/. is recorded and not written",
            unrecorded
        ]
    );
    assert_eq!(
        lines(r#"{".": "aa", "set/a": "aa"}"#),
        ["generated file . is recorded and not written"]
    );
    // A file recorded under two spellings is held to the digest of each, whichever the
    // manifest names first.
    let differs = ["generated file set/a hashes as aa, the manifest records bb"];
    assert_eq!(lines(r#"{"set/a": "aa", "./set/a": "bb"}"#), differs);
    assert_eq!(lines(r#"{"./set/a": "bb", "set/a": "aa"}"#), differs);
    assert_eq!(
        lines(r#"{"set/a": "aa", "set//a": "aa"}"#),
        Vec::<String>::new()
    );
}
