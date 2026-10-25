// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line as the `vep` binary reads it: `VEP_*` environment
//! variables behind the command line (`Config.pm:707-708, 1127-1168` at
//! release/116.2), long options abbreviated to a unique prefix, written after
//! one dash or in another case (Getopt::Long's defaults, `Config.pm:62`), and
//! the short flags, which clap reads.
//!
//! Every run sets `HOME` to an empty directory, so no `vep.ini` of the
//! machine's own reaches the binary.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use vep_cli::args::Args;
use vep_cli::argv::{parser_options, unresolved_perl_spellings, VEP_PARAMS};

fn corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116")
}

/// The binary with an empty home directory and no `VEP_*` variable of the
/// machine's own.
fn vep(home: &Path) -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    cmd.stdin(Stdio::null());
    cmd.env_clear();
    cmd.env("HOME", home);
    if let Some(path) = std::env::var_os("PATH") {
        cmd.env("PATH", path);
    }
    cmd
}

/// The golden corpus's own arguments, output path left to the caller.
fn corpus_args() -> Vec<String> {
    let dir = corpus_dir();
    vec![
        "-i".to_string(),
        dir.join("variants.vcf").display().to_string(),
        "--offline".to_string(),
        "--json_cache".to_string(),
        dir.join("json_cache").display().to_string(),
        "--assembly".to_string(),
        "GRCh38".to_string(),
        "--force_overwrite".to_string(),
        "--no_stats".to_string(),
        "--max_sv_size".to_string(),
        "-1".to_string(),
        "--quiet".to_string(),
    ]
}

fn json_cache() -> String {
    corpus_dir().join("json_cache").display().to_string()
}

/// A switch of `@VEP_PARAMS` the parser does not declare, by any name, for the
/// refusal tests; `None` once every one is declared.
fn undeclared_perl_switch() -> Option<String> {
    let mut command = Args::command_with_error_table();
    command.build();
    let parser = parser_options(&command);
    VEP_PARAMS
        .iter()
        .filter(|spec| !spec.contains(['=', ':']))
        .map(|spec| spec.split('|').next().unwrap_or(spec).to_string())
        .find(|name| !parser.iter().any(|option| option.names.contains(name)))
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The first line of the cache information `--show_cache_info` prints.
const CACHE_INFO_FIRST_LINE: &str = "1000genomes\tphase3\n";

fn assert_prints_cache_info(output: &Output, what: &str) {
    assert_eq!(
        output.status.code(),
        Some(0),
        "{what}: {}",
        text(&output.stderr)
    );
    assert!(
        text(&output.stdout).starts_with(CACHE_INFO_FIRST_LINE),
        "{what}: {}",
        text(&output.stdout)
    );
}

/// The header line of a tab output, which the default format lacks.
fn tab_header(output_file: &Path) -> bool {
    fs::read_to_string(output_file)
        .unwrap()
        .lines()
        .any(|l| l.starts_with("#Uploaded_variation\tLocation\tAllele"))
}

/// `VEP_<FLAG>=<value>` sets a flag the command line leaves unset
/// (`Config.pm:1158`): `VEP_DATABASE=1` reaches the refusal of `--database`;
/// `VEP_TAB=1` selects the tab output; `VEP_SHOW_CACHE_INFO=1` prints the cache
/// information and exits 0 as `--show_cache_info` does, the one flag the
/// environment sets here and not in Ensembl VEP, whose driver reads it from the
/// typed options alone (vep:39), as the CLI reference's compatibility page
/// states. A value that is not UTF-8 is read, with the replacement character for
/// the byte: `%ENV` holds bytes and Config.pm:1137 skips the empty value alone,
/// so `VEP_PROTEIN` set to the byte `0xff` switches `--protein` on.
#[test]
fn vep_env_sets_unset_flag() {
    let home = tempfile::tempdir().unwrap();
    let output = vep(home.path())
        .args(["--offline", "--json_cache", &json_cache()])
        .env("VEP_SHOW_CACHE_INFO", "1")
        .output()
        .unwrap();
    assert_prints_cache_info(&output, "VEP_SHOW_CACHE_INFO=1");

    let output = vep(home.path())
        .args(["-i", "/nonexistent/input.vcf"])
        .env("VEP_DATABASE", "1")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(
        text(&output.stderr).starts_with("ERROR: [cli] --database is not supported"),
        "{}",
        text(&output.stderr)
    );

    let out = home.path().join("out.txt");
    let mut args = corpus_args();
    args.extend(["-o".to_string(), out.display().to_string()]);
    let output = vep(home.path())
        .args(&args)
        .env("VEP_TAB", "1")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(tab_header(&out), "VEP_TAB=1 selects the tab output");

    #[cfg(unix)]
    {
        use std::ffi::OsStr;
        use std::os::unix::ffi::OsStrExt;
        let protein = |value: Option<&OsStr>| -> bool {
            let mut command = vep(home.path());
            command.args(&args);
            if let Some(value) = value {
                command.env("VEP_PROTEIN", value);
            }
            let output = command.output().unwrap();
            assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
            fs::read_to_string(&out).unwrap().contains("ENSP=")
        };
        assert!(!protein(None), "no ENSP without the flag");
        assert!(
            protein(Some(OsStr::from_bytes(b"\xff"))),
            "VEP_PROTEIN set to a byte outside UTF-8 switches --protein on"
        );
    }
}

/// A flag the command line gives keeps the command line's value whatever the
/// environment says (`Config.pm:693-697`): `--output_format tab` beside
/// `VEP_OUTPUT_FORMAT=vcf` writes tab output, and `--species homo_sapiens`
/// beside `VEP_SPECIES=mus_musculus` still finds the cache.
#[test]
fn vep_env_never_overrides_command_line() {
    let home = tempfile::tempdir().unwrap();
    let out = home.path().join("out.txt");
    let mut args = corpus_args();
    args.extend([
        "--output_format".to_string(),
        "tab".to_string(),
        "--species".to_string(),
        "homo_sapiens".to_string(),
        "-o".to_string(),
        out.display().to_string(),
    ]);
    let output = vep(home.path())
        .args(&args)
        .env("VEP_OUTPUT_FORMAT", "vcf")
        .env("VEP_SPECIES", "mus_musculus")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(
        tab_header(&out),
        "the command line's --output_format tab wins"
    );
    assert!(
        !fs::read_to_string(&out)
            .unwrap()
            .starts_with("##fileformat=VCF"),
        "VEP_OUTPUT_FORMAT=vcf did not win"
    );
}

/// A `VEP_*` variable naming no flag is ignored, as `is_valid_param` drops it
/// (`Config.pm:1143-1148`): an unknown name, an alias (`VEP_FA`, which Ensembl
/// VEP stores under a key nothing reads), a prefix, an empty value, and
/// `VEP_FASTA_STORAGE`, which the FASTA reader owns; the run proceeds as if
/// none were set.
#[test]
fn vep_env_unknown_key_is_ignored() {
    let home = tempfile::tempdir().unwrap();
    let output = vep(home.path())
        .args([
            "--offline",
            "--json_cache",
            &json_cache(),
            "--show_cache_info",
            "--quiet",
        ])
        .env("VEP_NO_SUCH_FLAG", "1")
        .env("VEP_FA", "/nonexistent/reference.fa")
        .env("VEP_INPUT_F", "/nonexistent/input.vcf")
        .env("VEP_DATABASE", "")
        .env("VEP_FASTA_STORAGE", "memory")
        .env("VEP_CONFIG", "/nonexistent/vep.ini")
        .env("vep_database", "1")
        .output()
        .unwrap();
    assert_prints_cache_info(&output, "unknown VEP_* variables");
    assert!(output.stderr.is_empty(), "{}", text(&output.stderr));
}

/// A unique prefix of a long option names it, as Getopt::Long's `auto_abbrev`
/// does (`Config.pm:62`): `--show_cache_i`, `--offl`, `--json_c` read as the
/// full flags, and `--force` stays the alias of `--force_overwrite`.
#[test]
fn unique_prefix_resolves() {
    let home = tempfile::tempdir().unwrap();
    let output = vep(home.path())
        .args(["--offl", "--json_c", &json_cache(), "--show_cache_i"])
        .output()
        .unwrap();
    assert_prints_cache_info(&output, "unique prefixes");

    let out = home.path().join("out.txt");
    let output = vep(home.path())
        .args([
            "-i",
            &corpus_dir().join("variants.vcf").display().to_string(),
            "--offl",
            "--json_c",
            &json_cache(),
            "--assem",
            "GRCh38",
            "--force",
            "--no_st",
            "--max_sv",
            "-1",
            "--qui",
            "--output_fo",
            "tab",
            "-o",
            &out.display().to_string(),
        ])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(tab_header(&out));
}

/// A prefix of several long options that Getopt::Long resolves no better is an
/// argument syntax error, exit 2 with clap's `error:` and no `ERROR: [` line:
/// `--fo` (fork, format, force_overwrite), `-fo` after one dash, `--hg` (the
/// hgvs flags), `-sh` and `-ou`, which are not `-s h` and `-o u`. A prefix only
/// vep-rs's own flags make ambiguous reads as Getopt::Long reads it over
/// `@VEP_PARAMS`: `--js` is `--json` and writes JSON, `--pa` is `--password`
/// and reaches its refusal.
#[test]
fn ambiguous_prefix_is_usage_error() {
    let home = tempfile::tempdir().unwrap();
    for argv in [["--fo"], ["-fo"], ["--hg"], ["-sh"], ["-ou"]] {
        let output = vep(home.path()).args(argv).output().unwrap();
        assert_eq!(output.status.code(), Some(2), "{argv:?}");
        let err = text(&output.stderr);
        assert!(err.starts_with("error: "), "{argv:?}: {err}");
        assert!(!err.contains("ERROR: ["), "{argv:?}: {err}");
    }

    let out = home.path().join("out.txt");
    let mut args = corpus_args();
    args.extend([
        "--js".to_string(),
        "-o".to_string(),
        out.display().to_string(),
    ]);
    let output = vep(home.path()).args(&args).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(
        fs::read_to_string(&out).unwrap().starts_with('{'),
        "--js is --json"
    );

    let output = vep(home.path()).args(["--pa", "secret"]).output().unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(
        text(&output.stderr).starts_with("ERROR: [cli] --password is not supported"),
        "{}",
        text(&output.stderr)
    );
}

/// An Ensembl VEP option this version lacks is refused by its full name in any
/// spelling (exit 2), never read as a cluster of shorts (`-stats_html` as
/// `-s tats_html`); the option is taken from `@VEP_PARAMS` at run time, so the
/// test holds as flags are declared. `-af_exac` reaches the error table's
/// refusal (exit 1) as `--af_exac` does.
#[test]
fn single_dash_unknown_option_is_refused_by_name() {
    let home = tempfile::tempdir().unwrap();
    if let Some(name) = undeclared_perl_switch() {
        let named = format!("'--{name}'");
        let glued_output = &name[1..];
        for spelling in [format!("-{name}"), format!("--{name}")] {
            let output = vep(home.path())
                .args([&spelling, "-i", "/nonexistent/x.vcf"])
                .current_dir(home.path())
                .output()
                .unwrap();
            let err = text(&output.stderr);
            assert_eq!(output.status.code(), Some(2), "{spelling}: {err}");
            assert!(
                err.starts_with("error: unexpected argument ") && err.contains(&named),
                "{spelling}: {err}"
            );
        }
        assert!(
            !home.path().join(glued_output).exists(),
            "-{name} did not write an output file named {glued_output}"
        );
    }

    let output = vep(home.path()).args(["-af_exac"]).output().unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(
        text(&output.stderr).starts_with("ERROR: [cli] --af_exac was removed"),
        "{}",
        text(&output.stderr)
    );
}

/// Every spelling of an Ensembl VEP option that Getopt::Long resolves against
/// `@VEP_PARAMS` (`Config.pm:73-276`), every name and alias and every
/// prefix of each, resolves to the same option here when the option is
/// declared: the pre-pass falls back to Getopt::Long's reading where a flag
/// vep-rs adds (`--json_cache`, the Parquet flags, `--version`) makes a prefix
/// ambiguous, and reads an alias the parser lacks (`id`). The assertion
/// message lists any spelling that does not, as the CLI reference does.
#[test]
fn every_perl_unique_prefix_resolves_the_same() {
    let mut command = Args::command_with_error_table();
    command.build();
    let unresolved = unresolved_perl_spellings(&command);
    let listed: Vec<String> = unresolved
        .iter()
        .map(|u| {
            format!(
                "{} -> {} (vep-rs: {})",
                u.spelling,
                u.perl,
                if u.candidates.is_empty() {
                    "no flag".to_string()
                } else {
                    u.candidates.join(", ")
                }
            )
        })
        .collect();
    assert!(
        unresolved.is_empty(),
        "Perl-unique spellings vep-rs does not resolve:\n  {}",
        listed.join("\n  ")
    );
}

/// A long option after one dash is the long option, as Getopt::Long reads it
/// with `bundling` off (`Config.pm:62`): `-offline -json_cache <dir>
/// -show_cache_info` prints the cache information, `-form vcf` is the
/// prefix, `-fa` the alias, `-ve` the one option the prefix names, and a run
/// spelled that way annotates as the double-dash one does.
#[test]
fn single_dash_long_option_is_rewritten() {
    let home = tempfile::tempdir().unwrap();
    let output = vep(home.path())
        .args(["-offline", "-json_cache", &json_cache(), "-show_cache_info"])
        .output()
        .unwrap();
    assert_prints_cache_info(&output, "single-dash long options");

    let out = home.path().join("out.txt");
    let output = vep(home.path())
        .args([
            "-input_file",
            &corpus_dir().join("variants.vcf").display().to_string(),
            "-offline",
            "-json_cache",
            &json_cache(),
            "-assembly",
            "GRCh38",
            "-force_overwrite",
            "-no_stats",
            "-max_sv_size",
            "-1",
            "-quiet",
            "-tab",
            "-output_file",
            &out.display().to_string(),
        ])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(tab_header(&out));

    let output = vep(home.path())
        .args([
            "-form",
            "vcf",
            "-fa",
            "/nonexistent/reference.fa",
            "-i",
            "/nonexistent/x.vcf",
        ])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
    assert!(
        text(&output.stderr).contains("ERROR: ["),
        "{}",
        text(&output.stderr)
    );
    assert!(
        !text(&output.stderr).starts_with("error:"),
        "-form and -fa parsed as flags: {}",
        text(&output.stderr)
    );

    // `-ve` is `verbose`, the one Ensembl VEP option the prefix names, not the
    // cluster `-v -e`: `--everything` would warn that it implies `--regulatory`.
    let output = vep(home.path())
        .args([
            "-ve",
            "--offline",
            "--json_cache",
            &json_cache(),
            "--show_cache_info",
        ])
        .output()
        .unwrap();
    assert_prints_cache_info(&output, "-ve");
    assert!(
        !text(&output.stderr).contains("implies --regulatory"),
        "-ve read as -v -e: {}",
        text(&output.stderr)
    );
}

/// Short flags and short clusters are clap's: `-i`, `-o`, `-q` alone; `-qv` is
/// `-q -v`, which the option rules refuse together (`Config.pm:620`); `-V` is
/// the version and `-v` the verbosity; a short with its value attached
/// (`-ob.txt`) is that short's value; a short taking a value takes the next
/// token whatever it looks like.
#[test]
fn short_flags_unchanged() {
    let home = tempfile::tempdir().unwrap();
    let out = home.path().join("out.txt");
    let output = vep(home.path())
        .args([
            "-i",
            &corpus_dir().join("variants.vcf").display().to_string(),
            "-o",
            &out.display().to_string(),
            "-q",
            "--offline",
            "--json_cache",
            &json_cache(),
            "--no_stats",
            "--max_sv_size",
            "-1",
        ])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(out.is_file());

    let output = vep(home.path()).args(["-qv"]).output().unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        text(&output.stderr),
        "ERROR: [cli] Can't use --quiet and --verbose together\n"
    );

    let output = vep(home.path()).arg("-V").output().unwrap();
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(
        text(&output.stdout),
        format!("vep {}\n", env!("CARGO_PKG_VERSION"))
    );

    let glued = home.path().join("glued.txt");
    let mut args = corpus_args();
    args.push(format!("-o{}", glued.display()));
    let output = vep(home.path()).args(&args).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(
        glued.is_file(),
        "-o with its value attached names the output file"
    );
    assert!(
        !text(&output.stderr).contains("More than one output file"),
        "{}",
        text(&output.stderr)
    );

    let output = vep(home.path())
        .args(["-i", "-cache", "--offline", "--json_cache", &json_cache()])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
    assert!(
        text(&output.stderr).contains("ERROR: [input] File \"-cache\" does not exist"),
        "{}",
        text(&output.stderr)
    );
}

/// A long option in another case names the flag, as Getopt::Long's
/// `ignore_case` does (`Config.pm:62`): `--OFFLINE --Json_Cache <dir>
/// --SHOW_CACHE_INFO` prints the cache information; a value keeps its case.
#[test]
fn long_option_case_is_folded() {
    let home = tempfile::tempdir().unwrap();
    let output = vep(home.path())
        .args([
            "--OFFLINE",
            "--Json_Cache",
            &json_cache(),
            "--SHOW_CACHE_INFO",
        ])
        .output()
        .unwrap();
    assert_prints_cache_info(&output, "upper-case long options");

    let output = vep(home.path())
        .args([
            "--Offline",
            "--Json_Cache",
            &json_cache(),
            "--Assembly",
            "grch38",
            "-i",
            "/nonexistent/x.vcf",
        ])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
    assert!(
        text(&output.stderr).contains("(grch38)"),
        "the value keeps its case: {}",
        text(&output.stderr)
    );
}

/// `--config <file>` on the command line reads the file through the binary:
/// its flags act (`offline 1`, `json_cache <dir>` and `show_cache_info 1` print
/// the cache information, the last a reading of this version's: Ensembl VEP's
/// driver tests the typed flag alone, vep:39, so a file's `show_cache_info`
/// prints nothing there, as the CLI reference's compatibility page states); a
/// key of an Ensembl VEP option this version lacks
/// (a switch of `@VEP_PARAMS` taken at run time) is refused as the flag is, exit
/// 2 naming it, where a key no option has (`no_such_flag`) is dropped; and a file
/// that cannot be opened is the one-line `config` error, exit 1 (`Config.pm:1075`).
#[test]
fn config_file_is_read_through_the_binary() {
    let home = tempfile::tempdir().unwrap();
    let file = home.path().join("flags.ini");
    fs::write(
        &file,
        format!(
            "# flags of the run\noffline 1\njson_cache {}\nshow_cache_info 1\n",
            json_cache()
        ),
    )
    .unwrap();
    let output = vep(home.path())
        .args(["--config", file.to_str().unwrap()])
        .output()
        .unwrap();
    assert_prints_cache_info(&output, "--config file");

    if let Some(name) = undeclared_perl_switch() {
        let unported = home.path().join("unported.ini");
        fs::write(&unported, format!("offline 1\n{name} 1\nno_such_flag 1\n")).unwrap();
        let output = vep(home.path())
            .args([
                "--config",
                unported.to_str().unwrap(),
                "--json_cache",
                &json_cache(),
                "--show_cache_info",
            ])
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2), "{}", text(&output.stderr));
        assert!(
            text(&output.stderr).contains(&format!("unexpected argument '--{name}'")),
            "{}",
            text(&output.stderr)
        );
    }

    let missing = home.path().join("absent.ini");
    let output = vep(home.path())
        .args([
            "--config",
            missing.to_str().unwrap(),
            "-i",
            "/nonexistent/x.vcf",
        ])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        text(&output.stderr),
        format!(
            "ERROR: [config] Could not open config file \"{}\"\n",
            missing.display()
        )
    );
    assert!(output.stdout.is_empty());
}

/// `<dir>/vep.ini` is read when it is a file, `--dir` naming the directory and
/// `~/.vep` by default (`Config.pm:704-705`).
#[test]
fn vep_ini_is_read_through_the_binary() {
    let home = tempfile::tempdir().unwrap();
    let dot_vep = home.path().join(".vep");
    fs::create_dir_all(&dot_vep).unwrap();
    fs::write(
        dot_vep.join("vep.ini"),
        format!(
            "offline 1\njson_cache {}\nshow_cache_info 1\n",
            json_cache()
        ),
    )
    .unwrap();
    let output = vep(home.path()).output().unwrap();
    assert_eq!(output.status.code(), Some(0));
    assert!(
        text(&output.stdout).starts_with("vep "),
        "an empty command line prints the help before any file is read"
    );
    let output = vep(home.path()).args(["--offline"]).output().unwrap();
    assert_prints_cache_info(&output, "~/.vep/vep.ini");

    let other = home.path().join("other");
    fs::create_dir_all(&other).unwrap();
    fs::write(other.join("vep.ini"), "show_cache_info 1\n").unwrap();
    let output = vep(home.path())
        .args([
            "--dir",
            other.to_str().unwrap(),
            "--offline",
            "--json_cache",
            &json_cache(),
        ])
        .output()
        .unwrap();
    assert_prints_cache_info(&output, "--dir vep.ini");
}

/// A repeatable key alone on its line of a `--config` file (`plugin`,
/// `custom`, `transcript_filter`) is dropped here (`read_config_file` in
/// `argv.rs`), where Ensembl VEP 116.2 keeps the empty string for it
/// (`Config.pm:1094-1097`, `my $value = join(' ', @split);` then
/// `push @{$config->{$key}}, $value;`): a bare `plugin` line reaches
/// `my ($module, @params) = split /,/, $plugin_config;` and the `use $module`
/// of `Runner.pm:928-935`, which warns `Failed to compile plugin`; a bare
/// `custom` line dies at `AnnotationSourceAdaptor.pm:376` (`No file was added
/// for custom annotation source`, exit 2) and a bare `transcript_filter` line
/// dies (exit 255). A scalar or switch key alone on its line is dropped by both
/// (`Config.pm:1101`, `$config->{$key} ||= $split[0];` with no value). Here
/// `protein 1` above a bare `plugin` line exits 0, annotates with `ENSP`,
/// raises no plugin warning, and prints `--protein` and no `--plugin` in the
/// command-line header.
#[test]
fn a_repeatable_key_alone_on_its_config_line_is_dropped() {
    let home = tempfile::tempdir().unwrap();
    let file = home.path().join("flags.ini");
    fs::write(&file, "protein 1\nplugin\n").unwrap();
    let out = home.path().join("out.txt");
    let mut args = corpus_args();
    args.extend([
        "--config".to_string(),
        file.display().to_string(),
        "-o".to_string(),
        out.display().to_string(),
    ]);
    let output = vep(home.path()).args(&args).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let warnings =
        fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default();
    assert!(
        !warnings.to_lowercase().contains("plugin"),
        "the bare plugin line raised a plugin warning: {warnings}"
    );
    let annotated = fs::read_to_string(&out).unwrap();
    let command_line = annotated
        .lines()
        .find(|l| l.starts_with("## VEP command-line:"))
        .unwrap();
    assert_eq!(
        command_line,
        "## VEP command-line: vep --assembly GRCh38 --config [PATH]/flags.ini --force_overwrite \
         --input_file [PATH]/variants.vcf --json_cache [PATH]/json_cache --max_sv_size -1 \
         --no_stats --offline --output_file [PATH]/out.txt --protein --quiet"
    );
    assert!(annotated.contains("ENSP="), "the file's protein line acts");
}

/// `VEP_AF_EXAC=1` reaches the refusal `--af_exac` reaches typed, the identical
/// line and exit 1, as `VEP_DATABASE=1` reaches `--database`'s: the key names a
/// flag the error table refuses before parsing, which `key_storage` in `argv.rs`
/// stores under its own name and `tokens_of` lays out as the switch `--af_exac`
/// for the scan in `config::parse_command_line`. `VEP_AF_EXAC=0` adds nothing,
/// a switch set to a false value. Ensembl VEP 116.2 ignores the variable:
/// `read_config_from_environment` keeps a key only when `is_valid_param`
/// matches it to an option of `@VEP_PARAMS` (Config.pm:1143-1147), which has no
/// `af_exac`.
#[test]
fn vep_env_af_exac_reaches_the_refusal() {
    let home = tempfile::tempdir().unwrap();
    let run = |value: &str| {
        vep(home.path())
            .args(["-i", "/nonexistent/input.vcf", "--quiet"])
            .env("VEP_AF_EXAC", value)
            .output()
            .unwrap()
    };
    let refused = run("1");
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr),
        "ERROR: [cli] --af_exac was removed in Ensembl VEP release 107: use --af_gnomade\n"
    );

    let off = run("0");
    let stderr = text(&off.stderr);
    assert_eq!(
        off.status.code(),
        Some(1),
        "the run fails later, on its cache: {stderr}"
    );
    assert!(!stderr.contains("af_exac"), "{stderr}");
}
