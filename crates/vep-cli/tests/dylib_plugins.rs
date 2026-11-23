// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The dylib plugin path through the `vep` binary on the pruned golden corpus: an
//! unresolvable `--plugin` name and a file that is not a plugin are each warned
//! about and skipped, the file an error under `--safe`; the `echo_plugin` example
//! resolved under `--dir_plugins` lands its fields in the Extra column of every
//! transcript row and nowhere else, its header info reaches every format's header
//! and columns, and the `echo_all_rows_plugin` example, which lists every kind
//! of row, runs on the regulatory, motif and intergenic rows too. A plugin that
//! fails at run time, the `faulty_plugin` example or a built-in whose data file
//! breaks under it, costs the row its own fields and one warning, never the run.

use std::path::{Path, PathBuf};
use std::process::Output;

mod common;

/// The warning an unresolvable name writes, up to the plugins directory it names.
const SKIPPED_WARNING: &str =
    "WARNING: Failed to load plugin NoSuchPlugin: NoSuchPlugin is not a built-in plugin and no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under ";

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37")
}

/// Runs `vep` on the corpus without `--quiet`, so warnings reach stderr; the
/// default output format unless `extra` names another.
fn vep_on_corpus(out: &Path, extra: &[&str]) -> Output {
    vep_on(&corpus_dir(), out, extra)
}

/// [`vep_on_corpus`] over any GRCh37 corpus directory.
fn vep_on(corpus: &Path, out: &Path, extra: &[&str]) -> Output {
    common::isolated_vep_anywhere()
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

fn data_lines(out: &Path) -> Vec<String> {
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

#[test]
fn unresolvable_plugin_name_is_warned_and_skipped() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let output = vep_on_corpus(&out, &["--plugin", "NoSuchPlugin,param1"]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(stderr.contains(SKIPPED_WARNING), "{stderr}");
    assert!(!data_lines(&out).is_empty(), "the run still annotates");
}

#[test]
fn unresolvable_plugin_name_is_skipped_with_a_plugin_dir_too() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let output = vep_on_corpus(
        &out,
        &["--plugin", "NoSuchPlugin", "--dir_plugins", "/no/such/dir"],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(stderr.contains(SKIPPED_WARNING), "{stderr}");
}

/// A file that is not a shared library is dropped with Ensembl VEP's warning for
/// a plugin module that fails to compile (Runner.pm:929-939) and the run goes on;
/// under `--safe` it is the plugin component's error and no output is written.
#[test]
fn a_file_that_is_not_a_plugin_is_dropped_and_fatal_under_safe() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let not_a_plugin = Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let not_a_plugin = not_a_plugin.to_str().unwrap();
    let plugins_dir = corpus_dir();
    let plugins_dir = plugins_dir.to_str().unwrap();
    let output = vep_on_corpus(
        &out,
        &["--plugin", not_a_plugin, "--dir_plugins", plugins_dir],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(
        stderr.contains(&format!(
            "WARNING: Failed to compile plugin {not_a_plugin}: {not_a_plugin}: "
        )),
        "{stderr}"
    );
    assert!(!data_lines(&out).is_empty(), "the run still annotates");
    std::fs::remove_file(&out).unwrap();

    let output = vep_on_corpus(
        &out,
        &[
            "--plugin",
            not_a_plugin,
            "--dir_plugins",
            plugins_dir,
            "--safe",
        ],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.contains(&format!(
            "ERROR: [plugin] Failed to compile plugin {not_a_plugin}: {not_a_plugin}: "
        )),
        "{stderr}"
    );
    assert!(
        !out.exists(),
        "no output is written when --safe refuses a plugin"
    );
}

/// The profile's `examples/` directory when the `echo_plugin`,
/// `echo_all_rows_plugin` and `faulty_plugin` example dylibs have been built
/// there. `cargo test` on a selection that includes `vep-plugin` (the workspace,
/// or `-p vep-plugin -p vep-cli`) builds them; `-p vep-cli` alone does not.
fn echo_plugin_dir() -> Option<PathBuf> {
    let exe = std::env::current_exe().unwrap();
    let examples = exe.parent().unwrap().parent().unwrap().join("examples");
    let ext = if cfg!(target_os = "macos") {
        "dylib"
    } else {
        "so"
    };
    let dylib = examples.join(format!("libecho_plugin.{ext}"));
    if [
        "libecho_plugin",
        "libecho_all_rows_plugin",
        "libfaulty_plugin",
    ]
    .iter()
    .all(|stem| examples.join(format!("{stem}.{ext}")).is_file())
    {
        return Some(examples);
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

#[test]
fn echo_plugin_under_dir_plugins_annotates_every_transcript_row() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let output = vep_on_corpus(
        &out,
        &[
            "--plugin",
            "echo_plugin,hello",
            "--dir_plugins",
            dir.to_str().unwrap(),
        ],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(!stderr.contains("not found"), "{stderr}");

    let lines = data_lines(&out);
    let transcript_rows: Vec<Vec<&str>> = lines
        .iter()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|cols| cols[5] == "Transcript")
        .collect();
    assert!(!transcript_rows.is_empty(), "corpus has transcript rows");
    for cols in &transcript_rows {
        let (feature, extra) = (cols[4], cols[13]);
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

/// A library named by its file name in the working directory loads: the token
/// carries a suffix, so it is read as a file and handed to the loader under `./`,
/// which opens it where a bare file name would be searched for on the loader's
/// own path. No Ensembl VEP line governs a file token: `use $module` compiles a
/// bare word alone (Runner.pm:931-933).
#[test]
fn a_library_named_by_its_file_name_in_the_working_directory_loads() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let ext = if cfg!(target_os = "macos") {
        "dylib"
    } else {
        "so"
    };
    let tmp = tempfile::tempdir().unwrap();
    let working_dir = tmp.path().join("work");
    std::fs::create_dir(&working_dir).unwrap();
    let file_name = format!("libecho_plugin.{ext}");
    std::fs::copy(dir.join(&file_name), working_dir.join(&file_name)).unwrap();
    let out = tmp.path().join("default.txt");
    let plugin = format!("{file_name},hello");
    let output = common::isolated_vep_anywhere()
        .current_dir(&working_dir)
        .arg("-i")
        .arg(corpus_dir().join("variants.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus_dir().join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh37"])
        .args(["--force_overwrite", "--no_stats", "--plugin", &plugin])
        .output()
        .expect("vep binary runs");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(!stderr.contains("Failed to"), "{stderr}");
    let annotated = data_lines(&out)
        .iter()
        .filter(|line| line.contains("ECHO_PARAM=hello"))
        .count();
    assert!(annotated > 0, "no row carries the plugin's field");
}

/// The header lines of `out`, `##` and `#` alike.
fn header_lines(out: &Path) -> Vec<String> {
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

/// The keys of the `## <KEY> : <description>` lines between `## Extra column
/// keys:` and the command line of a default-format header.
fn extra_column_keys(header: &[String]) -> Vec<String> {
    header
        .iter()
        .skip_while(|line| *line != "## Extra column keys:")
        .skip(1)
        .take_while(|line| !line.starts_with("## VEP command-line"))
        .map(|line| {
            line.trim_start_matches("## ")
                .split(" : ")
                .next()
                .unwrap()
                .to_string()
        })
        .collect()
}

/// Runs `vep` with `echo_plugin,hello` loaded under `--dir_plugins` plus `flags`,
/// and asserts the run succeeded with the plugin loaded.
fn run_echo(out: &Path, dir: &Path, flags: &[&str]) {
    let mut args = vec![
        "--plugin",
        "echo_plugin,hello",
        "--dir_plugins",
        dir.to_str().unwrap(),
    ];
    args.extend_from_slice(flags);
    let output = vep_on_corpus(out, &args);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(!stderr.contains("Failed to"), "{stderr}");
}

/// The dylib's `vep_plugin_header_info` reaches every format as the built-in
/// plugins' header info does (`OutputFactory/VEP_output.pm` and `Tab.pm`
/// `description_headers`, `Tab.pm` and `VCF.pm` `fields`, `VCF.pm` plugin
/// headers): the default header describes the two fields after the flag-field
/// lines, `--tab` describes them and adds the columns, `--vcf` adds a `##KEY=`
/// line per field and the CSQ Format entries, and the Parquet output has the
/// columns.
#[test]
fn echo_plugin_fields_are_described_in_every_header() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    const FEATURE_LINE: &str = "## ECHO_FEATURE : The transcript id the consequence names";
    const PARAM_LINE: &str = "## ECHO_PARAM : The first plugin parameter";

    let out = tmp.path().join("default.txt");
    run_echo(&out, &dir, &[]);
    let header = header_lines(&out);
    let keys = extra_column_keys(&header);
    assert!(keys.contains(&"FLAGS".to_string()), "{keys:?}");
    assert_eq!(
        &keys[keys.len() - 2..],
        ["ECHO_FEATURE", "ECHO_PARAM"],
        "{keys:?}"
    );
    assert!(header.contains(&FEATURE_LINE.to_string()), "{header:#?}");
    assert!(header.contains(&PARAM_LINE.to_string()), "{header:#?}");

    let out = tmp.path().join("tab.txt");
    run_echo(&out, &dir, &["--tab"]);
    let header = header_lines(&out);
    assert!(header.contains(&FEATURE_LINE.to_string()), "{header:#?}");
    assert!(header.contains(&PARAM_LINE.to_string()), "{header:#?}");
    let columns: Vec<&str> = header.last().unwrap().split('\t').collect();
    assert_eq!(columns[0], "#Uploaded_variation");
    assert_eq!(
        &columns[columns.len() - 2..],
        ["ECHO_FEATURE", "ECHO_PARAM"]
    );
    let feature_col = columns.iter().position(|c| *c == "Feature").unwrap();
    let echo_feature_col = columns.len() - 2;
    let echo_param_col = columns.len() - 1;
    let rows = data_lines(&out);
    assert!(!rows.is_empty());
    for row in &rows {
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols.len(), columns.len(), "{row}");
        if cols[5] == "Transcript" {
            assert_eq!(cols[echo_feature_col], cols[feature_col], "{row}");
            assert_eq!(cols[echo_param_col], "hello", "{row}");
        } else {
            assert_eq!(cols[echo_feature_col], "-", "{row}");
            assert_eq!(cols[echo_param_col], "-", "{row}");
        }
    }

    let out = tmp.path().join("out.vcf");
    run_echo(&out, &dir, &["--vcf"]);
    let header = header_lines(&out);
    assert!(
        header.contains(&"##ECHO_FEATURE=The transcript id the consequence names".to_string()),
        "{header:#?}"
    );
    assert!(
        header.contains(&"##ECHO_PARAM=The first plugin parameter".to_string()),
        "{header:#?}"
    );
    let info = header
        .iter()
        .find(|line| line.starts_with("##INFO=<ID=CSQ,"))
        .unwrap();
    let format = info
        .split("Format: ")
        .nth(1)
        .unwrap()
        .trim_end_matches("\">");
    let csq_fields: Vec<&str> = format.split('|').collect();
    assert_eq!(
        &csq_fields[csq_fields.len() - 2..],
        ["ECHO_FEATURE", "ECHO_PARAM"]
    );
    let feature_field = csq_fields.iter().position(|f| *f == "Feature").unwrap();
    let rows = data_lines(&out);
    let transcript_entries = rows
        .iter()
        .flat_map(|row| {
            let info = row.split('\t').nth(7).unwrap();
            let csq = info
                .split(';')
                .find_map(|kv| kv.strip_prefix("CSQ="))
                .unwrap_or("");
            csq.split(',')
                .map(|entry| entry.split('|').map(str::to_owned).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        })
        .filter(|entry| entry.len() == csq_fields.len() && entry[feature_field].starts_with("ENST"))
        .collect::<Vec<_>>();
    assert!(!transcript_entries.is_empty());
    for entry in &transcript_entries {
        assert_eq!(
            entry[csq_fields.len() - 2],
            entry[feature_field],
            "{entry:?}"
        );
        assert_eq!(entry[csq_fields.len() - 1], "hello", "{entry:?}");
    }

    // The Parquet output takes the tab columns, so the plugin's fields end its
    // schema.
    #[cfg(feature = "parquet")]
    {
        use parquet::file::reader::FileReader;
        let parquet_dir = tmp.path().join("out.parquet");
        run_echo(&parquet_dir, &dir, &["--output_format", "parquet"]);
        let parts = vep_io::output::parquet::part_files(&parquet_dir).unwrap();
        let reader = parquet::file::serialized_reader::SerializedFileReader::new(
            std::fs::File::open(&parts[0]).unwrap(),
        )
        .unwrap();
        let columns: Vec<String> = reader
            .metadata()
            .file_metadata()
            .schema()
            .get_fields()
            .iter()
            .map(|field| field.name().to_string())
            .collect();
        assert_eq!(
            &columns[columns.len() - 2..],
            ["ECHO_FEATURE", "ECHO_PARAM"]
        );
    }
}

/// `echo_plugin` lists `Transcript` alone, so the host never hands it the
/// intergenic row (OutputFactory.pm:2232, `check_feature_type(... || 'Intergenic')`).
#[test]
fn a_dylib_declaring_transcript_only_skips_intergenic_variants() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    run_echo(&out, &dir, &[]);
    let rows = data_lines(&out);
    let intergenic: Vec<&String> = rows
        .iter()
        .filter(|row| row.split('\t').nth(6).unwrap() == "intergenic_variant")
        .collect();
    assert!(!intergenic.is_empty(), "corpus has intergenic rows");
    for row in intergenic {
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols[5], "-", "{row}");
        assert!(!cols[13].contains("ECHO_"), "{row}");
    }
}

/// `echo_all_rows_plugin` lists every kind in `vep_plugin_feature_types`, so it
/// runs on every row: the intergenic row's consequence document is `null`, and
/// its result lands in the variant's fields, which that row's Extra column prints.
#[test]
fn a_dylib_listing_every_kind_runs_on_every_row_and_intergenic_variants() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let output = vep_on_corpus(
        &out,
        &[
            "--plugin",
            "echo_all_rows_plugin",
            "--dir_plugins",
            dir.to_str().unwrap(),
        ],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(!stderr.contains("Failed to"), "{stderr}");
    assert!(header_lines(&out)
        .contains(&"## ECHO_ROW : The kind of row: the feature type, or intergenic".to_string()));
    let rows = data_lines(&out);
    let (mut intergenic, mut transcript) = (0, 0);
    for row in &rows {
        let cols: Vec<&str> = row.split('\t').collect();
        let fields: Vec<&str> = cols[13].split(';').collect();
        match cols[5] {
            "Transcript" => {
                transcript += 1;
                assert!(fields.contains(&"ECHO_ROW=Transcript"), "{row}");
            }
            "-" => {
                intergenic += 1;
                assert_eq!(cols[6], "intergenic_variant", "{row}");
                assert!(fields.contains(&"ECHO_ROW=intergenic"), "{row}");
            }
            other => panic!("unexpected feature type {other} in {row}"),
        }
    }
    assert!(
        intergenic > 0 && transcript > 0,
        "{intergenic} intergenic, {transcript} transcript rows"
    );
}

/// Plugin header lines come plugin by plugin in load order, each plugin's keys
/// sorted (OutputFactory.pm:2416-2430 `get_plugin_headers`): a built-in's line
/// precedes a later dylib's, and the dylib's keys sort regardless of the order it
/// listed them in.
#[test]
fn header_keys_are_sorted_per_plugin_in_load_order() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let scores = tmp.path().join("pli.tsv");
    std::fs::write(&scores, "gene\tpLI\nSOD1\t0.5\n").unwrap();
    let out = tmp.path().join("default.txt");
    let pli = format!("pLI,{}", scores.display());
    let output = vep_on_corpus(
        &out,
        &[
            "--plugin",
            &pli,
            "--plugin",
            "echo_plugin,hello",
            "--dir_plugins",
            dir.to_str().unwrap(),
        ],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(!stderr.contains("Failed to"), "{stderr}");
    let keys = extra_column_keys(&header_lines(&out));
    let position = |key: &str| {
        keys.iter()
            .position(|k| k == key)
            .unwrap_or_else(|| panic!("{key} missing from {keys:?}"))
    };
    assert!(
        position("pLI_gene_value") < position("ECHO_FEATURE"),
        "{keys:?}"
    );
    assert!(
        position("ECHO_FEATURE") < position("ECHO_PARAM"),
        "{keys:?}"
    );
}

/// A plugin runs on the rows whose feature type it lists (`OutputFactory.pm:2232`
/// `check_feature_type`): the example's default `Transcript` list reaches no
/// regulatory or motif row, a list naming both kinds reaches both, and a plugin
/// that lists `MotifFeature` switches regulatory annotation on for the run
/// without `--regulatory` (`Runner.pm:998-1002`), where a `Transcript` list
/// leaves it off.
#[test]
fn plugins_run_on_the_rows_their_feature_types_declare() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let cache = common::regulatory_cache(tmp.path());
    let input = tmp.path().join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tinside\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    let run = |plugin: &str, regulatory: bool| -> Vec<(String, String, String)> {
        let out = tmp
            .path()
            .join(format!("out-{}-{regulatory}.txt", plugin.len()));
        let mut command = common::isolated_vep(tmp.path());
        command
            .arg("-i")
            .arg(&input)
            .arg("-o")
            .arg(&out)
            .args(["--json_cache", cache.to_str().unwrap()])
            .args(["--assembly", "GRCh38", "--no_fasta", "--no_stats"])
            .args(["--force_overwrite", "--plugin", plugin])
            .args(["--dir_plugins", dir.to_str().unwrap()]);
        if regulatory {
            command.arg("--regulatory");
        }
        let output = command.output().unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{stderr}");
        assert!(!stderr.contains("Failed to"), "{stderr}");
        data_lines(&out)
            .iter()
            .map(|line| {
                let cols: Vec<&str> = line.split('\t').collect();
                let echo = cols[13]
                    .split(';')
                    .find_map(|f| f.strip_prefix("ECHO_FEATURE="))
                    .unwrap_or("-")
                    .to_string();
                (cols[4].to_string(), cols[5].to_string(), echo)
            })
            .collect()
    };
    let rows = |echo_regulatory: &str, echo_motif: &str| {
        vec![
            (
                "ENSR_FIXTURE".to_string(),
                "RegulatoryFeature".to_string(),
                echo_regulatory.to_string(),
            ),
            (
                "ENSM_FORWARD".to_string(),
                "MotifFeature".to_string(),
                echo_motif.to_string(),
            ),
            ("-".to_string(), "-".to_string(), "-".to_string()),
        ]
    };
    assert_eq!(run("echo_plugin", true), rows("-", "-"));
    assert_eq!(
        run(
            "echo_plugin,feature_types=RegulatoryFeature;MotifFeature",
            true
        ),
        rows("ENSR_FIXTURE", "ENSM_FORWARD")
    );
    assert_eq!(
        run("echo_plugin,feature_types=MotifFeature", false),
        rows("-", "ENSM_FORWARD")
    );
    assert_eq!(
        run("echo_plugin", false),
        vec![("-".to_string(), "-".to_string(), "-".to_string())]
    );
}

/// The `WARNING:` lines of the warnings file beside `out`.
fn warning_lines(out: &Path) -> Vec<String> {
    std::fs::read_to_string(format!("{}_warnings.txt", out.display()))
        .unwrap_or_default()
        .lines()
        .map(str::to_owned)
        .collect()
}

/// `faulty_plugin` returns `[1]` on its first call and NULL on its second, then
/// `FAULTY_OK=1`: `run_plugins` warns `did not return a hashref, output ignored!`
/// for a result that is not a hash and `went wrong: <reason>` for a failed call
/// (OutputFactory.pm:2244, 2253), keeps the line and goes on, and the warning
/// sink writes each text once (BaseVEP.pm:810-813). The run exits 0 with every
/// row of a plugin-less run, the two warnings once each, and the field on every
/// transcript row but the two the failing calls served.
#[test]
fn a_dylib_that_returns_a_non_object_is_warned_and_the_run_goes_on() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let plain = tmp.path().join("plain.txt");
    assert!(vep_on_corpus(&plain, &[]).status.success());
    let out = tmp.path().join("faulty.txt");
    let output = vep_on_corpus(
        &out,
        &[
            "--plugin",
            "faulty_plugin",
            "--dir_plugins",
            dir.to_str().unwrap(),
        ],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(!stderr.contains("Failed to"), "{stderr}");

    let rows = data_lines(&out);
    assert_eq!(rows.len(), data_lines(&plain).len(), "every row is written");
    let warnings = warning_lines(&out);
    let count = |text: &str| warnings.iter().filter(|line| *line == text).count();
    assert_eq!(
        count("WARNING: Plugin 'faulty' did not return a hashref, output ignored!"),
        1,
        "{warnings:#?}"
    );
    assert_eq!(
        count("WARNING: Plugin 'faulty' went wrong: plugin 'faulty' run returned null"),
        1,
        "{warnings:#?}"
    );
    assert_eq!(
        warnings
            .iter()
            .filter(|line| line.contains("'faulty'"))
            .count(),
        2,
        "{warnings:#?}"
    );
    assert!(header_lines(&out).contains(
        &"## FAULTY_OK : 1 on every row after the plugin's two failing calls".to_string()
    ));

    let (mut with_field, mut without) = (0, 0);
    for row in &rows {
        let cols: Vec<&str> = row.split('\t').collect();
        let has_field = cols[13].split(';').any(|field| field == "FAULTY_OK=1");
        match cols[5] {
            "Transcript" if has_field => with_field += 1,
            "Transcript" => without += 1,
            _ => assert!(!has_field, "{row}"),
        }
    }
    assert_eq!(without, 2, "the two failing calls served two rows");
    assert!(with_field > 0, "the later calls annotate");
}

/// A built-in that fails on every variant of a batch the index maps to a
/// block costs its own fields and one `Plugin '<Name>' went wrong: <reason>`
/// warning (OutputFactory.pm:2253), and the built-in after it annotates the
/// batch: REVEL over a data file whose index names blocks the file lacks (the
/// file cut after its header block, so the plugin loads and every read of a
/// block fails; the record past the file's last row reads none), then pLI
/// over a values file naming a gene of the corpus.
#[test]
fn a_builtin_failure_is_one_warning_and_the_other_builtins_annotate() {
    let corpus =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37-plugin-revel");
    let tmp = tempfile::tempdir().unwrap();
    let revel = std::fs::read(corpus.join("plugin_data/revel.tsv.gz")).unwrap();
    // The first bgzf block holds the header line alone; its length is the
    // block's BSIZE field plus one.
    let block_len = usize::from(u16::from_le_bytes([revel[16], revel[17]])) + 1;
    let data = tmp.path().join("revel.tsv.gz");
    std::fs::write(&data, &revel[..block_len]).unwrap();
    std::fs::copy(
        corpus.join("plugin_data/revel.tsv.gz.tbi"),
        tmp.path().join("revel.tsv.gz.tbi"),
    )
    .unwrap();
    let values = tmp.path().join("pli.tsv");
    std::fs::write(&values, "gene\tpLI\nCBS\t0.5\n").unwrap();

    let out = tmp.path().join("default.txt");
    let revel_arg = format!("REVEL,file={}", data.display());
    let pli_arg = format!("pLI,{}", values.display());
    let output = vep_on(
        &corpus,
        &out,
        &["--symbol", "--plugin", &revel_arg, "--plugin", &pli_arg],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(!stderr.contains("Failed to"), "{stderr}");

    let warnings = warning_lines(&out);
    let went_wrong: Vec<&String> = warnings
        .iter()
        .filter(|line| line.starts_with("WARNING: Plugin 'REVEL' went wrong: "))
        .collect();
    assert_eq!(went_wrong.len(), 1, "{warnings:#?}");
    assert!(
        went_wrong[0].contains("does not describe"),
        "{}",
        went_wrong[0]
    );
    assert_eq!(
        warnings
            .iter()
            .filter(|line| line.contains("Plugin '"))
            .count(),
        1,
        "{warnings:#?}"
    );

    let rows = data_lines(&out);
    let cbs_rows: Vec<&String> = rows
        .iter()
        .filter(|row| {
            row.split('\t')
                .nth(13)
                .unwrap()
                .split(';')
                .any(|field| field == "SYMBOL=CBS")
        })
        .collect();
    assert!(!cbs_rows.is_empty(), "the corpus annotates CBS");
    for row in &cbs_rows {
        assert!(row.contains("pLI_gene_value=0.50"), "{row}");
    }
    for row in &rows {
        assert!(!row.contains("REVEL="), "{row}");
    }
}

/// The bgzf blocks of `file` as `(offset, length)`, each length the `BSIZE`
/// field of the block's header plus one.
fn bgzf_blocks(file: &[u8]) -> Vec<(usize, usize)> {
    let mut blocks = Vec::new();
    let mut offset = 0;
    while offset < file.len() {
        let length = usize::from(u16::from_le_bytes([file[offset + 16], file[offset + 17]])) + 1;
        blocks.push((offset, length));
        offset += length;
    }
    blocks
}

/// The start of a default-format row's `Location` column.
fn row_start(row: &str) -> u64 {
    let location = row.split('\t').nth(1).unwrap();
    let span = location.split_once(':').unwrap().1;
    span.split('-').next().unwrap().parse().unwrap()
}

/// `row` with its `REVEL` field removed from the Extra column.
fn without_revel(row: &str) -> String {
    let mut columns: Vec<&str> = row.split('\t').collect();
    let extra: Vec<&str> = columns[13]
        .split(';')
        .filter(|field| !field.starts_with("REVEL="))
        .collect();
    let extra = extra.join(";");
    columns[13] = &extra;
    columns.join("\t")
}

/// Ensembl VEP 116.3's REVEL reads each line's own span (REVEL.pm:168) inside
/// `run_plugins`'s per-line `eval` (OutputFactory.pm:2234-2235), and its tabix
/// iterator ends at a block it cannot inflate, so a corrupt block costs the
/// variants read from it their REVEL field and no other. A batch whose prefetch
/// fails is retried variant by variant: over the corpus's REVEL file with its
/// second data block corrupted, every buffer size writes the same REVEL fields,
/// the intact file's on every variant before the block and none from the block,
/// every other column as the intact file's, and one `went wrong` warning naming
/// the first variant the plugin failed on, the deletion ending on the last line
/// before the block, whose read runs into it, and the number of variants the
/// plugin failed on over the run. The warning reads the same at the default
/// buffer, at a buffer of two (the corpus split into batches that fail in part,
/// in whole and not at all) and at a buffer of one: the sink sums the batches
/// into the one line it writes for the plugin, as `warning_msg`'s seen set
/// leaves one line for a plugin dying the same way on every line
/// (BaseVEP.pm:810-813).
#[test]
fn a_corrupt_block_costs_the_variants_read_from_it_their_field_alone() {
    let corpus =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37-plugin-revel");
    let tmp = tempfile::tempdir().unwrap();
    let intact_file = corpus.join("plugin_data/revel.tsv.gz");
    let mut revel = std::fs::read(&intact_file).unwrap();
    let blocks = bgzf_blocks(&revel);
    // The header block, two data blocks and the empty end-of-file block.
    assert_eq!(blocks.len(), 4, "{blocks:?}");
    let (offset, length) = blocks[2];
    // The deflate payload, between the 18-byte header and the 8-byte trailer.
    revel[offset + 18..offset + length - 8].fill(0);
    let data = tmp.path().join("revel.tsv.gz");
    std::fs::write(&data, &revel).unwrap();
    std::fs::copy(
        corpus.join("plugin_data/revel.tsv.gz.tbi"),
        tmp.path().join("revel.tsv.gz.tbi"),
    )
    .unwrap();

    let run = |name: &str, extra: &[&str]| {
        let out = tmp.path().join(name);
        let output = vep_on(&corpus, &out, extra);
        assert!(
            output.status.success(),
            "exit {:?}\n{}",
            output.status.code(),
            String::from_utf8_lossy(&output.stderr)
        );
        let went_wrong: Vec<String> = warning_lines(&out)
            .into_iter()
            .filter(|line| line.starts_with("WARNING: Plugin 'REVEL' went wrong: "))
            .collect();
        (data_lines(&out), went_wrong)
    };
    let intact_arg = format!("REVEL,file={}", intact_file.display());
    let corrupt_arg = format!("REVEL,file={}", data.display());
    let (intact, intact_warnings) = run("intact.txt", &["--plugin", &intact_arg]);
    assert!(intact_warnings.is_empty(), "{intact_warnings:#?}");
    let (batched, batched_warnings) = run("batched.txt", &["--plugin", &corrupt_arg]);
    let (in_twos, twos_warnings) = run(
        "twos.txt",
        &["--buffer_size", "2", "--plugin", &corrupt_arg],
    );
    let (one_by_one, single_warnings) = run(
        "single.txt",
        &["--buffer_size", "1", "--plugin", &corrupt_arg],
    );

    let revel_rows = |rows: &[String]| -> Vec<String> {
        rows.iter()
            .filter(|row| row.contains("REVEL="))
            .cloned()
            .collect()
    };
    // The corrupt block holds the file's rows from GRCh37 position 44479076 on.
    let (before_block, in_block): (Vec<String>, Vec<String>) = revel_rows(&intact)
        .into_iter()
        .partition(|row| row_start(row) < 44_479_076);
    assert!(!before_block.is_empty() && !in_block.is_empty());
    let plain =
        |rows: &[String]| -> Vec<String> { rows.iter().map(|r| without_revel(r)).collect() };
    for rows in [&batched, &in_twos, &one_by_one] {
        assert_eq!(revel_rows(rows), before_block);
        assert_eq!(plain(rows), plain(&intact));
    }

    assert_eq!(batched_warnings.len(), 1, "{batched_warnings:#?}");
    let (note, count) = batched_warnings[0]
        .rsplit_once("; ")
        .expect("the note then the count");
    assert!(
        note.ends_with(
            "; the plugin annotated nothing on the variants it failed on, the first at \
             21:44479062-44479075"
        ),
        "{}",
        batched_warnings[0]
    );
    // The deletion whose read runs into the block and the variants whose spans
    // are read from it, six in all whatever the buffer size.
    assert_eq!(count, "6 variants affected", "{}", batched_warnings[0]);
    assert_eq!(twos_warnings, batched_warnings);
    assert_eq!(single_warnings, batched_warnings);
}

/// A built-in loaded twice is two plugins of one name (Runner.pm:925-954 loads
/// one instance per `--plugin` argument), each failing on its own file with
/// its own reason, so the run's warnings file carries one `went wrong` line per
/// instance, each naming its own file and count, where one line keyed on the
/// name alone would carry the first file's reason and the summed count.
#[test]
fn two_instances_of_a_failing_plugin_each_write_their_own_line() {
    let corpus =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37-plugin-revel");
    let tmp = tempfile::tempdir().unwrap();
    let intact = std::fs::read(corpus.join("plugin_data/revel.tsv.gz")).unwrap();
    let blocks = bgzf_blocks(&intact);
    let (offset, length) = blocks[2];
    let mut corrupt = intact.clone();
    corrupt[offset + 18..offset + length - 8].fill(0);
    let mut args = Vec::new();
    for name in ["first", "second"] {
        let data = tmp.path().join(format!("{name}.tsv.gz"));
        std::fs::write(&data, &corrupt).unwrap();
        std::fs::copy(
            corpus.join("plugin_data/revel.tsv.gz.tbi"),
            tmp.path().join(format!("{name}.tsv.gz.tbi")),
        )
        .unwrap();
        args.push(format!("REVEL,file={}", data.display()));
    }
    let out = tmp.path().join("two.txt");
    let output = vep_on(&corpus, &out, &["--plugin", &args[0], "--plugin", &args[1]]);
    assert!(
        output.status.success(),
        "exit {:?}\n{}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr)
    );
    let went_wrong: Vec<String> = warning_lines(&out)
        .into_iter()
        .filter(|line| line.starts_with("WARNING: Plugin 'REVEL' went wrong: "))
        .collect();
    assert_eq!(went_wrong.len(), 2, "{went_wrong:#?}");
    for (line, name) in went_wrong.iter().zip(["first", "second"]) {
        assert!(line.contains(&format!("{name}.tsv.gz")), "{line}");
        assert!(line.ends_with("6 variants affected"), "{line}");
    }
}

/// A plugin that lists `MotifFeature` switches regulatory annotation on for the
/// run and the status line names it (`Runner.pm:999-1001`, the shape of
/// `t/Runner.t:741-759`, where a loaded plugin's `feature_types` is what the
/// runner reads); it then sees the motif rows and no other. A plugin listing
/// `Intergenic` runs on the intergenic row of a record that overlaps no
/// transcript, its regulatory rows notwithstanding, and its field prints on that
/// row alone (`OutputFactory.pm:2232`); a plugin listing every kind prints on
/// every row.
#[test]
fn a_plugin_declaring_motif_features_switches_regulatory_on_and_sees_motif_rows() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let cache = common::regulatory_cache(tmp.path());
    let input = tmp.path().join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tinside\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    // (Feature_type, the plugin keys of the row's Extra column)
    let run = |plugins: &[&str], flags: &[&str]| -> (String, Vec<(String, Vec<String>)>) {
        let out = tmp
            .path()
            .join(format!("out-{}.txt", plugins.join("").len()));
        let mut command = common::isolated_vep(tmp.path());
        command
            .arg("-i")
            .arg(&input)
            .arg("-o")
            .arg(&out)
            .args(["--json_cache", cache.to_str().unwrap()])
            .args(["--assembly", "GRCh38", "--no_fasta", "--no_stats"])
            .args(["--force_overwrite", "--dir_plugins", dir.to_str().unwrap()]);
        for plugin in plugins {
            command.args(["--plugin", plugin]);
        }
        command.args(flags);
        let output = command.output().unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        assert!(output.status.success(), "{stderr}");
        assert!(!stderr.contains("Failed to"), "{stderr}");
        let rows = data_lines(&out)
            .iter()
            .map(|line| {
                let cols: Vec<&str> = line.split('\t').collect();
                let keys = cols[13]
                    .split(';')
                    .filter(|f| f.starts_with("ECHO_"))
                    .map(str::to_string)
                    .collect();
                (cols[5].to_string(), keys)
            })
            .collect();
        (stderr, rows)
    };
    let keys = |fields: &[&str]| fields.iter().map(|f| f.to_string()).collect::<Vec<_>>();

    let (stderr, rows) = run(&["echo_plugin,feature_types=MotifFeature"], &[]);
    assert!(
        stderr.contains("Fetching regulatory features for plugin: echo"),
        "{stderr}"
    );
    assert_eq!(
        rows,
        [
            ("RegulatoryFeature".to_string(), keys(&[])),
            (
                "MotifFeature".to_string(),
                keys(&["ECHO_FEATURE=ENSM_FORWARD"])
            ),
            ("-".to_string(), keys(&[])),
        ]
    );

    // The intergenic row's consequence document is `null`, so the plugin's
    // `ECHO_FEATURE` is empty there, which the default format writes as `KEY=`
    // (VEP_output.pm:107-119), beside its parameter echo.
    let (stderr, rows) = run(
        &["echo_plugin,hello,feature_types=Intergenic"],
        &["--regulatory"],
    );
    assert!(!stderr.contains("Fetching regulatory features"), "{stderr}");
    assert_eq!(
        rows,
        [
            ("RegulatoryFeature".to_string(), keys(&[])),
            ("MotifFeature".to_string(), keys(&[])),
            (
                "-".to_string(),
                keys(&["ECHO_FEATURE=", "ECHO_PARAM=hello"])
            ),
        ]
    );

    // A field the plugin's header info never declared stays on the intergenic
    // row it was written on as well: Ensembl VEP keeps a plugin's result in the
    // hash of the row it ran on (OutputFactory.pm:2229-2244).
    let (_, rows) = run(
        &["echo_plugin,hello,undeclared,feature_types=Intergenic"],
        &["--regulatory"],
    );
    assert_eq!(
        rows,
        [
            ("RegulatoryFeature".to_string(), keys(&[])),
            ("MotifFeature".to_string(), keys(&[])),
            (
                "-".to_string(),
                keys(&["ECHO_FEATURE=", "ECHO_PARAM=hello", "ECHO_UNDECLARED=u"])
            ),
        ]
    );

    let (stderr, rows) = run(&["echo_all_rows_plugin"], &[]);
    assert!(
        stderr.contains("Fetching regulatory features for plugin: echo_all_rows"),
        "{stderr}"
    );
    assert_eq!(
        rows,
        [
            (
                "RegulatoryFeature".to_string(),
                keys(&["ECHO_ROW=RegulatoryFeature"])
            ),
            ("MotifFeature".to_string(), keys(&["ECHO_ROW=MotifFeature"])),
            ("-".to_string(), keys(&["ECHO_ROW=intergenic"])),
        ]
    );
}
