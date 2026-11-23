// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `generated` on a copy of the two-record corpus, from a copy of the binary inside a
//! tree of the test's own: `--command` records the words and the digests of the files they
//! write, `--check` holds the record against a second run, and a command that does not
//! start, fails, writes no file or writes one that cannot be read is refused by name. The
//! commands are shell scripts the test writes.

mod common;

use common::{
    binary_in_tree, clean, ended, failed, manifest, manifest_file, manifest_text, names, object,
    refused, scratch, two_record_corpus, write_manifest_line, Run, ROOT,
};
use sha2::{Digest, Sha256};
use std::fs;
use std::path::Path;
use std::process::Command;
use vep_golden_corpus::json::{Object, Value};
use vep_golden_corpus::TempDir;

/// The shell that runs the scripts.
const SHELL: &str = "/bin/sh";

/// Writes under the directory of its first argument, for the corpus of its second: two
/// files, a count of the corpus's record lines in a directory below, and the directory it
/// runs from. A third argument fails it with that argument on its last line of standard
/// error, and a directory named `nothing` is left unwritten.
const GENERATOR: &str = r#"out=$1
corpus=$2
if [ $# -gt 2 ]; then
    printf 'first line\n%s\n' "$3" >&2
    exit 3
fi
if [ "$(basename "$out")" = nothing ]; then
    exit 0
fi
mkdir -p "$out/deep"
printf 'b\n' > "$out/b.txt"
printf 'a\n' > "$out/a.txt"
printf '%s\n' "$(grep -c '' "$corpus/variants.vcf")" > "$out/deep/records.txt"
printf '%s' "$(pwd -P)" > "$out/cwd.txt"
"#;

/// Writes one file of mode 0 under the directory of its first argument.
const SEALED: &str = r#"mkdir "$1"
: > "$1/a.txt"
chmod 0 "$1/a.txt"
"#;

/// The sha256 of `data` in hex digits.
fn sha256(data: &[u8]) -> String {
    Sha256::digest(data)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// Runs the binary `binary` from the package root with `args`.
fn run_in_tree(binary: &str, args: &[&str]) -> Run {
    ended(
        Command::new(binary)
            .current_dir(ROOT)
            .env_remove("COLUMNS")
            .args(args),
    )
}

/// Whether this user reads a file of mode 0, as the superuser does.
fn reads_a_sealed_file(directory: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    let sealed = directory.join("sealed");
    fs::write(&sealed, "").unwrap();
    fs::set_permissions(&sealed, fs::Permissions::from_mode(0o000)).unwrap();
    fs::File::open(&sealed).is_ok()
}

#[test]
fn generated_writes_and_checks_the_block() {
    // `generated --command` runs the words from the repository root, `{generated}` a
    // scratch directory and `{corpus}` the corpus, and writes the manifest's `generated`
    // block: the words as given and the digest of every file under that directory, in path
    // order. `--check` runs the recorded command, writes nothing and fails, naming each
    // path, on a digest that differs, a recorded file not written and a written file not
    // recorded. A command that does not start, fails, writes no file or writes one that
    // cannot be read is refused by name, a recorded word that holds a null byte among those
    // that do not start. Every word after `--command` is the command's, a `--` and an
    // option of the subcommand among them.
    let scratch = scratch("generated");
    let out = two_record_corpus(&scratch.path);
    let path = manifest_file(&out);
    let named = path.to_str().unwrap();
    let untouched = manifest_text(&out);
    let tree = TempDir::create(Path::new(env!("CARGO_TARGET_TMPDIR")), "tree-").unwrap();
    let binary = binary_in_tree(&tree.path);
    let root = fs::canonicalize(&tree.path).unwrap();
    let generated = |args: &[&str]| {
        let mut all = vec!["generated", "--corpus", &out];
        all.extend(args);
        run_in_tree(&binary, &all)
    };
    let generator = common::write(&scratch.path, "write_set.sh", GENERATOR);
    let command = [SHELL, &generator, "{generated}/set", "{corpus}"];
    let shown = command.join(" ");
    let with_command = |words: &[&str]| {
        let mut args = vec!["--command"];
        args.extend(words);
        generated(&args)
    };
    assert_eq!(
        with_command(&command),
        clean(0, &format!("{out}: 4 generated files recorded\n"))
    );
    let vcf = fs::read_to_string(Path::new(&out).join("variants.vcf")).unwrap();
    let record_count = format!("{}\n", vcf.lines().count());
    let digests = [
        ("set/a.txt", sha256(b"a\n")),
        ("set/b.txt", sha256(b"b\n")),
        ("set/cwd.txt", sha256(root.to_str().unwrap().as_bytes())),
        ("set/deep/records.txt", sha256(record_count.as_bytes())),
    ];
    let mut written = manifest(&out);
    let block = written.remove("generated").expect("the block");
    assert_eq!(
        block.as_object().unwrap().keys().collect::<Vec<_>>(),
        ["command", "sha256"]
    );
    assert_eq!(names(block.as_object().unwrap(), "command"), command);
    let held = block.get("sha256").and_then(Value::as_object).unwrap();
    let held: Vec<(&str, &str)> = held
        .iter()
        .map(|(file, digest)| (file, digest.as_str().unwrap()))
        .collect();
    let expected: Vec<(&str, &str)> = digests
        .iter()
        .map(|(file, digest)| (*file, digest.as_str()))
        .collect();
    assert_eq!(held, expected);
    assert_eq!(written, object(&untouched));
    let recorded = manifest_text(&out);

    let check = |extra: &[&str]| {
        let mut args = vec!["--check"];
        args.extend(extra);
        generated(&args)
    };
    assert_eq!(
        check(&[]),
        clean(
            0,
            &format!("{out}: 4 generated files written, 0 differ from the manifest's record\n")
        )
    );
    assert_eq!(manifest_text(&out), recorded);

    let zeros = "0".repeat(64);
    let mut moved = object(&recorded);
    let moved_digests = moved
        .get_mut("generated")
        .and_then(|block| block.as_object_mut()?.get_mut("sha256")?.as_object_mut())
        .unwrap();
    moved_digests.insert("set/a.txt", zeros.as_str());
    moved_digests.remove("set/b.txt");
    moved_digests.insert("set/gone.txt", "1".repeat(64));
    write_manifest_line(&out, &moved);
    let differing = manifest_text(&out);
    assert_eq!(
        check(&[]),
        Run {
            code: 1,
            stdout: format!(
                "{out}: 4 generated files written, 3 differ from the manifest's record\n"
            ),
            stderr: format!(
                "ERROR: [build_golden_corpus] generated file set/a.txt hashes as {}, the manifest records {zeros}\n\
                 ERROR: [build_golden_corpus] generated file set/b.txt is written and not recorded\n\
                 ERROR: [build_golden_corpus] generated file set/gone.txt is recorded and not written\n",
                digests[0].1
            ),
        }
    );
    assert_eq!(manifest_text(&out), differing);
    // The command is the manifest's own under --check, and a block is owed to check one.
    let mut again = vec!["--command"];
    again.extend(command);
    assert_eq!(
        check(&again),
        failed("generated --check runs the manifest's own command: give no --command")
    );
    fs::write(&path, &untouched).unwrap();
    assert_eq!(
        check(&[]),
        failed(&format!("{named} has no generated block to check"))
    );
    let wanted = "generated wants --command <word>... with no empty word (or --check)";
    assert_eq!(generated(&[]), failed(wanted));
    assert_eq!(with_command(&[]), failed(wanted));
    assert_eq!(with_command(&[SHELL, ""]), failed(wanted));
    let mut failing = command.to_vec();
    failing.push("no set written");
    assert_eq!(
        with_command(&failing),
        failed(&format!(
            "`{shown} no set written` failed with status 3: no set written"
        ))
    );
    // The generator fails on a third argument and names it: the `--` reaches it as a word.
    let mut marked = command.to_vec();
    marked.extend(["--", "--check"]);
    assert_eq!(
        with_command(&marked),
        failed(&format!("`{shown} -- --check` failed with status 3: --"))
    );
    assert_eq!(
        with_command(&[SHELL, &generator, "{generated}/nothing", "{corpus}"]),
        failed(&format!(
            "`{SHELL} {generator} {{generated}}/nothing {{corpus}}` wrote no file under {{generated}}"
        ))
    );
    assert_eq!(manifest_text(&out), untouched);
    // The subcommand reads no abbreviated option: `--command` is spelled whole or refused.
    let mut abbreviated = vec!["--comm"];
    abbreviated.extend(command);
    let ran = generated(&abbreviated);
    assert_eq!((ran.code, ran.stdout.as_str()), (2, ""));
    assert!(
        ran.stderr.ends_with(&format!(
            "vep-golden-corpus: error: unrecognized arguments: --comm {shown}\n"
        )),
        "{}",
        ran.stderr
    );
    assert_eq!(manifest_text(&out), untouched);
    let absent = common::under(&scratch.path, "no_such_program");
    assert_eq!(
        with_command(&[&absent, "{generated}"]),
        failed(&format!(
            "`{absent} {{generated}}` does not start ([Errno 2] No such file or directory: '{absent}')"
        ))
    );
    assert_eq!(manifest_text(&out), untouched);
    // The superuser reads a file of any mode, so the refusal is any other user's.
    if !reads_a_sealed_file(&scratch.path) {
        let sealed = common::write(&scratch.path, "sealed.sh", SEALED);
        let ran = with_command(&[SHELL, &sealed, "{generated}/set"]);
        let opening = format!(
            "ERROR: [build_golden_corpus] `{SHELL} {sealed} {{generated}}/set` wrote a file that cannot be read ([Errno 13] Permission denied: '"
        );
        assert_eq!((ran.code, ran.stdout.as_str()), (1, ""));
        assert!(
            ran.stderr.starts_with(&opening) && ran.stderr.ends_with("/set/a.txt')\n"),
            "{}",
            ran.stderr
        );
        assert_eq!(manifest_text(&out), untouched);
    }
    let mut unstartable: Object = object(&recorded);
    let words = unstartable
        .get_mut("generated")
        .and_then(|block| block.as_object_mut())
        .unwrap();
    words.insert("command", Value::strings(&[SHELL, "a\0b", "{generated}"]));
    write_manifest_line(&out, &unstartable);
    assert_eq!(
        check(&[]),
        refused(
            1,
            &format!(
                "ERROR: [build_golden_corpus] `{SHELL} a\0b {{generated}}` does not start (embedded null byte)"
            )
        )
    );
}
