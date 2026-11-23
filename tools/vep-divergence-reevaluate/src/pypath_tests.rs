// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::testing::{cases, Scratch};
use crate::value::Value;

#[test]
fn a_path_is_written_and_taken_apart_as_pathlib_does() {
    for case in cases("path") {
        let path = case.get_str("path").expect("a path");
        let text = normalise(path);
        let expect = |key: &str| case.get_str(key).expect("an expected text");
        assert_eq!(text, expect("str"), "str(Path({path:?}))");
        assert_eq!(name(&text), expect("name"), "Path({path:?}).name");
        assert_eq!(parent(&text), expect("parent"), "Path({path:?}).parent");
        assert_eq!(
            suffix(name(&text)),
            expect("suffix"),
            "Path({path:?}).suffix"
        );
        assert_eq!(stem(name(&text)), expect("stem"), "Path({path:?}).stem");
        let Some(Value::List(joins)) = case.get("joins") else {
            panic!("no joins");
        };
        for pair in joins {
            let Value::List(pair) = pair else {
                panic!("a pair")
            };
            let other = pair[0].as_str().unwrap();
            assert_eq!(
                join(&text, other),
                pair[1].as_str().unwrap(),
                "{path:?} / {other:?}"
            );
        }
        let Some(Value::List(renames)) = case.get("with_name") else {
            panic!("no renames");
        };
        for triple in renames {
            let Value::List(triple) = triple else {
                panic!("a triple")
            };
            let new = triple[0].as_str().unwrap();
            let expected = match (&triple[1], &triple[2]) {
                (Value::Str(renamed), _) => Ok(renamed.clone()),
                (_, Value::Str(error)) => Err(error.clone()),
                other => panic!("neither a path nor an error: {other:?}"),
            };
            assert_eq!(
                with_name(&text, new),
                expected,
                "Path({path:?}).with_name({new:?})"
            );
        }
    }
}

#[test]
fn a_path_is_resolved_through_its_links_and_past_what_does_not_exist() {
    let scratch = Scratch::new("realpath");
    let real = scratch.mkdir("real/inner");
    scratch.write("real/inner/file.txt", b"x");
    std::os::unix::fs::symlink(&real, scratch.path("link")).unwrap();
    std::os::unix::fs::symlink("real/inner/../inner", scratch.path("relative")).unwrap();
    std::os::unix::fs::symlink("loop-b", scratch.path("loop-a")).unwrap();
    std::os::unix::fs::symlink("loop-a", scratch.path("loop-b")).unwrap();
    let root = &scratch.0;
    assert_eq!(
        resolve(&scratch.path("link/file.txt")),
        format!("{real}/file.txt")
    );
    assert_eq!(
        resolve(&scratch.path("relative/file.txt")),
        format!("{real}/file.txt")
    );
    assert_eq!(
        resolve(&scratch.path("link/../x")),
        format!("{root}/real/x")
    );
    assert_eq!(
        resolve(&scratch.path("absent/../y/./z//")),
        format!("{root}/y/z")
    );
    assert_eq!(
        resolve(&scratch.path("absent/deeper")),
        format!("{root}/absent/deeper")
    );
    assert_eq!(
        resolve(&scratch.path("loop-a/x")),
        format!("{root}/loop-a/x")
    );
    assert_eq!(resolve("/"), "/");
    assert_eq!(resolve("/.."), "/");
    let here = std::env::current_dir().unwrap();
    assert_eq!(
        resolve("no-such-entry-of-the-directory"),
        format!(
            "{}/no-such-entry-of-the-directory",
            here.canonicalize().unwrap().display()
        )
    );
}

#[test]
fn a_relative_path_is_made_absolute_with_no_link_resolved() {
    let scratch = Scratch::new("absolute");
    let real = scratch.mkdir("real");
    std::os::unix::fs::symlink(&real, scratch.path("link")).unwrap();
    // An absolute path stands as written, its link and its `..` kept.
    let through_link = scratch.path("link/../link/file");
    assert_eq!(absolute(&through_link), through_link);
    let here = std::env::current_dir().unwrap().display().to_string();
    assert_eq!(absolute("some/file"), format!("{here}/some/file"));
    assert_eq!(absolute("../file"), format!("{here}/../file"));
    assert_eq!(absolute("."), here);
}

#[test]
fn a_failure_on_a_path_reads_as_an_oserror() {
    let scratch = Scratch::new("oserror");
    let missing = scratch.path("missing.txt");
    let error = std::fs::read(&missing).unwrap_err();
    assert_eq!(
        os_error(&error, &missing),
        format!("[Errno 2] No such file or directory: '{missing}'")
    );
    let error = std::fs::read(&scratch.0).unwrap_err();
    assert_eq!(
        os_error(&error, &scratch.0),
        format!("[Errno 21] Is a directory: '{}'", scratch.0)
    );
    let file = scratch.write("f", b"x");
    let under = format!("{file}/sub");
    let error = std::fs::create_dir_all(&under).unwrap_err();
    assert_eq!(
        os_error(&error, &under),
        format!("[Errno 20] Not a directory: '{under}'")
    );
}
