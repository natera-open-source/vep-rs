// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How `vep` settles a `--plugin` argument no built-in claims and no library
//! answers to, through the binary on the pruned golden corpus: an Ensembl VEP
//! plugin vep-rs does not provide ends the run with the plugin component's error
//! naming the reason, with or without `--safe`; a dylib of that name under
//! `--dir_plugins` loads first; a genuinely unknown name keeps the warning and
//! the run; a Perl module, by path or as `<dir>/<Name>.pm`, is refused by name;
//! a remote data file fails a built-in's instantiation; and a library of another
//! ABI version is dropped with a warning, an error under `--safe`.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use vep_builtins::UNPORTED_ENSEMBL_PLUGINS;

mod common;

const DYLIB_EXT: &str = if cfg!(target_os = "macos") {
    "dylib"
} else {
    "so"
};

/// The sentence every unported-plugin refusal ends with.
const SUPPORT_TABLE: &str =
    "The plugin support table in docs/plugins.md lists every Ensembl VEP plugin and its status.";

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37")
}

/// Runs `vep` on the corpus without `--quiet`, so warnings reach stderr, in the
/// default output format.
fn vep_on_corpus(out: &Path, extra: &[&str]) -> Output {
    vep_on_corpus_in(None, out, extra)
}

/// [`vep_on_corpus`] with the process started in `working_dir` when one is given.
fn vep_on_corpus_in(working_dir: Option<&Path>, out: &Path, extra: &[&str]) -> Output {
    let corpus = corpus_dir();
    let mut command = common::isolated_vep_anywhere();
    if let Some(dir) = working_dir {
        command.current_dir(dir);
    }
    command
        .arg("-i")
        .arg(corpus.join("variants.vcf"))
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh37"])
        .args(["--force_overwrite", "--no_stats"])
        .args(extra)
        .output()
        .expect("vep binary runs")
}

fn stderr_of(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

fn data_lines(out: &Path) -> Vec<String> {
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

/// Asserts that `output` is the plugin component's refusal carrying `message` as
/// one whole stderr line, exit 1, with no output file written.
fn assert_refused(output: &Output, out: &Path, message: &str) {
    let stderr = stderr_of(output);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    let line = format!("ERROR: [plugin] {message}");
    assert!(stderr.lines().any(|l| l == line), "{line}\n{stderr}");
    assert!(
        !out.exists(),
        "no output is written when a plugin is refused"
    );
}

/// Every entry of the unported table ends the run with its reason, the pointer
/// to the support table and exit 1, the plugin given bare or with parameters,
/// with or without `--safe`.
#[test]
fn every_unported_ensembl_plugin_name_refuses_with_its_reason() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    assert!(!UNPORTED_ENSEMBL_PLUGINS.is_empty());
    for (name, reason) in UNPORTED_ENSEMBL_PLUGINS {
        let message = format!(
            "{name} is an Ensembl VEP plugin vep-rs does not provide: {reason} {SUPPORT_TABLE}"
        );
        let bare = vep_on_corpus(&out, &["--plugin", name]);
        assert_refused(&bare, &out, &message);
        let with_params = format!("{name},file=x,1");
        let safe = vep_on_corpus(&out, &["--plugin", &with_params, "--safe"]);
        assert_refused(&safe, &out, &message);
    }
}

/// The working directory is never read for a bare name: a directory, an empty
/// file and a Perl-shaped file, each named like an unported plugin as the
/// plugin's own archive unpacks beside the input, leave the refusal in place
/// with its documented line, bare or with parameters, with or without `--safe`.
/// Ensembl VEP 116.2 finds a plugin module through `@INC`, onto which only the
/// plugins directories are pushed (Runner.pm:921, :931-933).
#[test]
fn a_working_directory_entry_named_like_an_unported_plugin_keeps_the_refusal() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let (name, reason) = UNPORTED_ENSEMBL_PLUGINS
        .iter()
        .find(|(name, _)| *name == "GeneSplicer")
        .expect("GeneSplicer is an unported plugin");
    let message = format!(
        "{name} is an Ensembl VEP plugin vep-rs does not provide: {reason} {SUPPORT_TABLE}"
    );
    type Entry = (&'static str, fn(&Path));
    let entries: [Entry; 3] = [
        ("directory", |entry| std::fs::create_dir(entry).unwrap()),
        ("empty-file", |entry| std::fs::write(entry, "").unwrap()),
        ("module-shaped-file", |entry| {
            std::fs::write(entry, "package GeneSplicer;\n1;\n").unwrap()
        }),
    ];
    let with_params = format!("{name},binary=/x/genesplicer,training=/x/human");
    for (kind, create) in entries {
        let working_dir = tmp.path().join(kind);
        std::fs::create_dir(&working_dir).unwrap();
        create(&working_dir.join(name));
        for plugin in [*name, with_params.as_str()] {
            for safe in [&[][..], &["--safe"][..]] {
                let mut args = vec!["--plugin", plugin];
                args.extend_from_slice(safe);
                let output = vep_on_corpus_in(Some(&working_dir), &out, &args);
                assert_refused(&output, &out, &message);
                assert!(
                    !stderr_of(&output).contains("Failed to compile plugin"),
                    "{kind}: {}",
                    stderr_of(&output)
                );
            }
        }
    }
}

/// `--safe` changes nothing about the refusal: the same line and the same exit.
#[test]
fn safe_does_not_change_the_unported_plugin_refusal() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let (name, reason) = UNPORTED_ENSEMBL_PLUGINS
        .iter()
        .find(|(name, _)| *name == "CAPICE")
        .expect("CAPICE is an unported plugin");
    let message = format!(
        "{name} is an Ensembl VEP plugin vep-rs does not provide: {reason} {SUPPORT_TABLE}"
    );
    let plain = vep_on_corpus(&out, &["--plugin", "CAPICE,file=x"]);
    let safe = vep_on_corpus(&out, &["--plugin", "CAPICE,file=x", "--safe"]);
    assert_refused(&plain, &out, &message);
    assert_refused(&safe, &out, &message);
    let refusal = |output: &Output| -> Vec<String> {
        stderr_of(output)
            .lines()
            .filter(|l| l.starts_with("ERROR: "))
            .map(str::to_owned)
            .collect()
    };
    assert_eq!(refusal(&plain), refusal(&safe));
}

/// The profile's `examples/` directory when the `echo_plugin` example dylib has
/// been built there (`cargo test` on a selection that includes `vep-plugin`).
fn echo_plugin() -> Option<PathBuf> {
    let exe = std::env::current_exe().unwrap();
    let examples = exe.parent().unwrap().parent().unwrap().join("examples");
    let dylib = examples.join(format!("libecho_plugin.{DYLIB_EXT}"));
    if dylib.is_file() {
        return Some(dylib);
    }
    // Under CI the whole workspace is under test, so an absent example is a build
    // defect rather than a narrower selection.
    assert!(
        std::env::var_os("CI").is_none(),
        "echo_plugin example dylib not found at {}",
        dylib.display()
    );
    eprintln!(
        "Skipping: echo_plugin example dylib not found at {}",
        dylib.display()
    );
    None
}

/// A library under `--dir_plugins` named like an unported plugin loads: the
/// table is consulted only once no built-in and no library answers to the
/// name. The process runs in a directory holding an entry of that name, which
/// is neither the plugin nor in its way.
#[test]
fn a_dylib_named_like_an_unported_plugin_loads() {
    let Some(echo) = echo_plugin() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let plugins = tmp.path().join("plugins");
    std::fs::create_dir(&plugins).unwrap();
    std::fs::copy(&echo, plugins.join(format!("libCAPICE.{DYLIB_EXT}"))).unwrap();
    let working_dir = tmp.path().join("work");
    std::fs::create_dir_all(working_dir.join("CAPICE")).unwrap();
    let out = tmp.path().join("default.txt");
    let output = vep_on_corpus_in(
        Some(&working_dir),
        &out,
        &[
            "--plugin",
            "CAPICE,hello",
            "--dir_plugins",
            plugins.to_str().unwrap(),
        ],
    );
    let stderr = stderr_of(&output);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(!stderr.contains("does not provide"), "{stderr}");

    let transcript_rows: Vec<Vec<String>> = data_lines(&out)
        .iter()
        .map(|line| line.split('\t').map(str::to_owned).collect::<Vec<_>>())
        .filter(|cols| cols[5] == "Transcript")
        .collect();
    assert!(!transcript_rows.is_empty(), "corpus has transcript rows");
    for cols in &transcript_rows {
        let (feature, extra) = (&cols[4], &cols[13]);
        let fields: Vec<&str> = extra.split(';').collect();
        assert!(
            fields.contains(&format!("ECHO_FEATURE={feature}").as_str()),
            "row {} lacks ECHO_FEATURE={feature}: {extra}",
            cols[0]
        );
        assert!(
            fields.contains(&"ECHO_PARAM=hello"),
            "row {} lacks ECHO_PARAM=hello: {extra}",
            cols[0]
        );
    }
}

/// A name that is neither a built-in, a library, an unported plugin nor a Perl
/// module keeps the warning and the run (Runner.pm:929-939).
#[test]
fn a_genuinely_unknown_name_still_warns_and_skips() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let output = vep_on_corpus(&out, &["--plugin", "NoSuchPlugin,param1"]);
    let stderr = stderr_of(&output);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(
        stderr.contains("WARNING: Failed to load plugin NoSuchPlugin: NoSuchPlugin is not a built-in plugin and no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under "),
        "{stderr}"
    );
    assert!(!data_lines(&out).is_empty(), "the run still annotates");
}

/// A Perl module is refused by name before any open: given by path, or as a bare
/// name that resolves to `<dir>/<Name>.pm` under `--dir_plugins`, with or
/// without `--safe`.
#[test]
fn a_perl_module_path_is_refused_as_not_runnable() {
    let tmp = tempfile::tempdir().unwrap();
    let plugins = tmp.path().join("plugins");
    std::fs::create_dir(&plugins).unwrap();
    let module = plugins.join("Foo.pm");
    std::fs::write(
        &module,
        "package Foo;\nuse base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);\n1;\n",
    )
    .unwrap();
    let out = tmp.path().join("default.txt");
    let message = format!(
        "{} is an Ensembl VEP Perl plugin; vep-rs runs no Perl plugin, see docs/plugins.md",
        module.display()
    );
    let by_path = format!("{},param", module.display());
    for safe in [&[][..], &["--safe"][..]] {
        let mut args = vec!["--plugin", by_path.as_str()];
        args.extend_from_slice(safe);
        assert_refused(&vep_on_corpus(&out, &args), &out, &message);

        let mut args = vec![
            "--plugin",
            "Foo",
            "--dir_plugins",
            plugins.to_str().unwrap(),
        ];
        args.extend_from_slice(safe);
        assert_refused(&vep_on_corpus(&out, &args), &out, &message);
    }
}

/// A built-in whose data file is a URL fails its instantiation with the reason
/// naming the file: a warnings-file line and the run goes on without the plugin,
/// the plugin component's error under `--safe` (Runner.pm:946-954).
#[test]
fn a_remote_data_file_is_refused_by_name() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    const REASON: &str = "Failed to instantiate plugin REVEL: remote data files are not supported; download https://example.invalid/x.tsv.gz and index it locally";
    let plugin = "REVEL,file=https://example.invalid/x.tsv.gz";
    let warned = vep_on_corpus(&out, &["--plugin", plugin]);
    let stderr = stderr_of(&warned);
    assert!(
        warned.status.success(),
        "exit {:?}\n{stderr}",
        warned.status.code()
    );
    let warning = format!("WARNING: {REASON}");
    assert!(stderr.lines().any(|l| l == warning), "{stderr}");
    let warnings_file = tmp.path().join("default.txt_warnings.txt");
    let warnings = std::fs::read_to_string(&warnings_file).unwrap();
    assert!(warnings.lines().any(|l| l == warning), "{warnings}");
    assert!(!data_lines(&out).is_empty(), "the run still annotates");
    std::fs::remove_file(&out).unwrap();

    let fatal = vep_on_corpus(&out, &["--plugin", plugin, "--safe"]);
    assert_refused(&fatal, &out, REASON);
}

/// A library reporting ABI version 1 against the version 2 host is dropped as a
/// module that fails to compile (Runner.pm:929-939), the run going on without
/// it; under `--safe` it is the plugin component's error. The fixture is the
/// `echo_plugin` example compiled with its `abi-one` feature, under which
/// `vep_plugin_abi_version` returns 1.
#[test]
fn safe_makes_an_abi_mismatch_fatal() {
    let source =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../vep-plugin/examples/echo_plugin.rs");
    let out_dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(&out_dir).unwrap();
    let fixture = out_dir.join(format!(
        "libecho_plugin_abi_one_cli_{}.{DYLIB_EXT}",
        std::process::id()
    ));
    let rustc = std::env::var_os("RUSTC").unwrap_or_else(|| "rustc".into());
    let compiled = Command::new(rustc)
        .args(["--edition", "2021", "--crate-type", "cdylib"])
        .args(["--crate-name", "echo_plugin_abi_one_cli"])
        .args(["--cfg", "feature=\"abi-one\""])
        .arg("-o")
        .arg(&fixture)
        .arg(&source)
        .output()
        .expect("rustc runs");
    assert!(
        compiled.status.success(),
        "compiling the abi-one fixture failed:\n{}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let plugin = fixture.to_str().unwrap();
    let message = format!(
        "Failed to compile plugin {plugin}: incompatible plugin ABI version: expected 2, got 1"
    );

    let warned = vep_on_corpus(&out, &["--plugin", plugin]);
    let stderr = stderr_of(&warned);
    assert!(
        warned.status.success(),
        "exit {:?}\n{stderr}",
        warned.status.code()
    );
    assert!(
        stderr.lines().any(|l| l == format!("WARNING: {message}")),
        "{stderr}"
    );
    assert!(!data_lines(&out).is_empty(), "the run still annotates");
    std::fs::remove_file(&out).unwrap();

    let fatal = vep_on_corpus(&out, &["--plugin", plugin, "--safe"]);
    assert_refused(&fatal, &out, &message);
    std::fs::remove_file(&fixture).unwrap();
}
