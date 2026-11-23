// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The files of a corpus: the order of the walk, the files that are read as text, the
//! reason of one that cannot be read, and the counts across the files of one directory.
//! Every identifier here is made up, and is written as its letters and its digits apart.

use super::super::gzip::compress;
use super::*;
use crate::scratch::TempDir;
use std::fs;

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-stored-").unwrap()
}

/// A made-up identifier: its letters, then its digits.
fn made_up(letters: &str, digits: &str) -> String {
    format!("{letters}{digits}")
}

fn place(directory: &TempDir, name: &str, data: &[u8]) -> String {
    let path = directory.path.join(name);
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::write(&path, data).unwrap();
    path.to_str().unwrap().to_string()
}

#[test]
fn the_files_of_a_directory_come_in_the_order_of_their_parts() {
    let scratch = scratch();
    let root = scratch.path.to_str().unwrap();
    for name in ["a/b/f.txt", "a.d/g", "a.txt", ".hidden", "z/y"] {
        place(&scratch, name, b"x");
    }
    fs::create_dir(scratch.path.join("empty")).unwrap();
    std::os::unix::fs::symlink("a", scratch.path.join("link_dir")).unwrap();
    std::os::unix::fs::symlink("a.txt", scratch.path.join("link_file")).unwrap();
    std::os::unix::fs::symlink("nowhere", scratch.path.join("broken")).unwrap();
    // `a/b/f.txt` comes before `a.d/g` and `a.txt`: a part is compared whole, where the
    // path as one text would put `.` before `/`. A link to a directory is not walked.
    let expected = [".hidden", "a/b/f.txt", "a.d/g", "a.txt", "link_file", "z/y"];
    let expected: Vec<String> = expected
        .iter()
        .map(|name| format!("{root}/{name}"))
        .collect();
    assert_eq!(corpus_files(root).unwrap(), expected);
    assert_eq!(
        corpus_files(&format!("{root}/empty")).unwrap(),
        Vec::<String>::new()
    );
}

#[test]
fn a_file_name_outside_utf8_is_refused_with_its_directory() {
    use std::os::unix::ffi::OsStrExt;
    let scratch = scratch();
    let held = scratch.path.join("expected");
    fs::create_dir(&held).unwrap();
    fs::write(held.join(std::ffi::OsStr::from_bytes(b"tab\xff.txt")), b"x").unwrap();
    assert_eq!(
        corpus_files(scratch.path.to_str().unwrap()).unwrap_err(),
        Failure::own(format!(
            "{} holds a file whose name is not UTF-8",
            held.display()
        ))
    );
}

#[test]
fn a_path_reads_below_the_directory_that_holds_it() {
    assert_eq!(
        below("/tmp/corpus/expected/tab.txt", "/tmp/corpus"),
        "expected/tab.txt"
    );
    assert_eq!(
        below("/tmp/corpus/expected/tab.txt", "/tmp"),
        "corpus/expected/tab.txt"
    );
    assert_eq!(below("corpus/tab.txt", "corpus"), "tab.txt");
    assert_eq!(below(".hidden/tab.txt", "."), ".hidden/tab.txt");
    assert_eq!(below("/tab.txt", "/"), "tab.txt");
    assert_eq!(below("corpus2/tab.txt", "corpus"), "corpus2/tab.txt");
}

#[test]
fn a_file_is_read_as_text_where_it_can_hold_an_identifier() {
    let scratch = scratch();
    let first = made_up("COSV", "12");
    let line = format!("x|{first}|y\n");
    let plain = place(&scratch, "plain.txt", line.as_bytes());
    assert_eq!(
        read_corpus_text(&plain).unwrap(),
        Some(CorpusText {
            raw: line.clone().into_bytes(),
            gzipped: false,
            blocked: false,
            text: line.clone(),
        })
    );
    // A gzip file is inflated whatever its name.
    let stored = compress(line.as_bytes());
    let named = place(&scratch, "stored.bin", &stored);
    assert_eq!(
        read_corpus_text(&named).unwrap(),
        Some(CorpusText {
            raw: stored,
            gzipped: true,
            blocked: false,
            text: line.clone(),
        })
    );
    // A member whose extra field opens with `BC` is a block of a block-compressed file.
    let mut block = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x00\x00".to_vec();
    block.extend(&compress(line.as_bytes())[10..]);
    let blocked = place(&scratch, "blocked.gz", &block);
    let read = read_corpus_text(&blocked).unwrap().unwrap();
    assert!(read.gzipped && read.blocked && read.text == line);
    // No identifier, a zero byte, and bytes outside UTF-8: nothing to read.
    for (name, data) in [
        ("none.txt", b"rs1\trs2\n".to_vec()),
        ("index.csi", [b"\0".as_slice(), line.as_bytes()].concat()),
        ("latin.txt", [b"\xe9 ".as_slice(), line.as_bytes()].concat()),
        ("empty.txt", Vec::new()),
    ] {
        let path = place(&scratch, name, &data);
        assert_eq!(read_corpus_text(&path).unwrap(), None, "{name}");
    }
}

#[test]
fn a_file_that_cannot_be_read_is_refused_with_its_path() {
    let scratch = scratch();
    let absent = scratch.path.join("absent.txt");
    let absent = absent.to_str().unwrap();
    assert_eq!(
        read_corpus_text(absent).unwrap_err(),
        failure(format!(
            "{absent} is not a readable file ([Errno 2] No such file or directory: '{absent}')"
        ))
    );
    let cut = place(&scratch, "cut.gz", &compress(b"rs1\n")[..12]);
    assert_eq!(
        read_corpus_text(&cut).unwrap_err(),
        failure(format!(
            "{cut} is not a readable file (Compressed file ended before the end-of-stream \
             marker was reached)"
        ))
    );
    let mut trailed = compress(b"rs1\n");
    trailed.extend(b"no");
    let trailed = place(&scratch, "trailed.gz", &trailed);
    assert_eq!(
        read_corpus_text(&trailed).unwrap_err(),
        failure(format!(
            "{trailed} is not a readable file (Not a gzipped file (b'no'))"
        ))
    );
}

#[test]
fn the_files_of_one_directory_are_counted_as_one_corpus() {
    let scratch = scratch();
    let root = scratch.path.to_str().unwrap().to_string();
    let first = made_up("COSV", "12");
    let second = made_up("CM", "654321");
    let other = made_up("CI", "654321");
    place(&scratch, "a.txt", format!("rs1,{second}\n").as_bytes());
    place(
        &scratch,
        "b/alone.txt",
        format!("{second}\t1\n{other}\t2\n").as_bytes(),
    );
    place(
        &scratch,
        "c.txt.gz",
        &compress(format!("x {first} {first}\n").as_bytes()),
    );
    place(
        &scratch,
        "d.txt",
        format!("rs1,{}\n", made_up("CM", "000000")).as_bytes(),
    );
    let files = corpus_files(&root).unwrap();
    let Counted { counts, seen } = restricted_identifier_counts(&files).unwrap();
    // The identifier `a.txt` lists is counted alone in a cell of `b/alone.txt`; the one no
    // file lists is not, and a file that holds the fabricated form alone has no count.
    assert_eq!(
        counts,
        [
            (format!("{root}/a.txt"), 1),
            (format!("{root}/b/alone.txt"), 1),
            (format!("{root}/c.txt.gz"), 2),
        ]
    );
    assert_eq!(seen, BTreeSet::from([first, second]));
    let refusal = "under restricted redistribution terms in another form than the fabricated \
                   one; the identifiers subcommand replaces them";
    assert_eq!(
        restricted_identifier_errors(&files, Some(&root)).unwrap(),
        [
            format!("a.txt holds 1 identifier {refusal}"),
            format!("b/alone.txt holds 1 identifier {refusal}"),
            format!("c.txt.gz holds 2 identifiers {refusal}"),
        ]
    );
    let lone = [format!("{root}/b/alone.txt")];
    assert_eq!(
        restricted_identifier_errors(&lone, None).unwrap(),
        Vec::<String>::new()
    );
}
