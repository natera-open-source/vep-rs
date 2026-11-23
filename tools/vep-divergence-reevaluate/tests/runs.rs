// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs of the binary over fixture reproductions through stand-in engines: no annotation
//! runs, a stand-in writes the output a fixture holds.

mod common;

use common::{sorted_results, Scenario, FIXTURES};

fn fixture_records() -> String {
    format!("{FIXTURES}/reproductions")
}

fn fixture_outputs() -> String {
    format!("{FIXTURES}/outputs")
}

/// The lines of an arguments file a stand-in engine wrote.
fn arguments(scenario: &Scenario, relative: &str) -> Vec<String> {
    scenario
        .read(relative)
        .split('\n')
        .map(str::to_string)
        .collect()
}

fn after<'a>(args: &'a [String], flag: &str) -> &'a str {
    let index = args
        .iter()
        .position(|arg| arg == flag)
        .unwrap_or_else(|| panic!("no {flag} in {args:?}"));
    &args[index + 1]
}

fn has(args: &[String], flag: &str) -> bool {
    args.iter().any(|arg| arg == flag)
}

#[test]
fn the_two_fixture_reproductions_run_through_fake_engines() {
    let scenario = Scenario::new("fixtures");
    let ensembl = scenario.stand_in("ensembl-vep", &fixture_outputs());
    let vep_rs = scenario.stand_in("vep-rs", &fixture_outputs());
    let root = scenario.cache_root();
    let records = fixture_records();
    let work = scenario.path("work");
    let results = scenario.path("results.json");
    let args = [
        "--ensembl-vep",
        &ensembl,
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &records,
        "--work-dir",
        &work,
        "--json",
        &results,
    ];
    let ran = scenario.run(&args);
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert!(ran
        .stdout
        .contains("still_present_fixture\tstill-present\n"));
    assert!(ran.stdout.contains(
        "closed_fixture\tclosed\tEnsembl VEP ['ENSP2:p.Leu4%3D'], vep-rs ['ENSP2:p.Leu4%3D']\n"
    ));
    for row in ["| still-present | 1 |", "| closed | 1 |", "| total | 2 |"] {
        assert!(ran.stdout.contains(row), "{row}");
    }
    // The closed entry's result lists what implements it; the entry still present lists
    // nothing.
    let source = "tools/vep-divergence-reevaluate/tests/fixtures/implementation/recorded.rs";
    for row in [
        format!("| closed_fixture | closed | code | {source}::PLAIN_EQUALS_SIGN |"),
        format!("| closed_fixture | closed | tests | {source}::tests::the_equals_sign_is_written_plain |"),
        format!("| closed_fixture | closed | rules | {source}::is_escaped_equals_sign |"),
    ] {
        assert!(ran.stdout.contains(&row), "{row}");
    }
    assert!(!ran.stdout.contains("| still_present_fixture |"));
    scenario.golden("fixtures", &ran);
    scenario.golden_file("fixtures.json", "results.json");
    let recorded = scenario.read("results.json");
    assert!(recorded.contains("\"verdict\": \"closed\""));
    assert_eq!(recorded.matches("\"implementation\": {").count(), 1);
    assert!(recorded.contains(&format!(
        "\"implementation\": {{\n   \"code\": [\n    \"{source}::PLAIN_EQUALS_SIGN\"\n   ],"
    )));
    assert!(recorded.contains(
        "\"vep_rs\": [\n   \"21:100 | A | ENST1 | downstream_gene_variant | 5000\"\n  ]"
    ));
    assert!(recorded.contains("\"entry\": \"still_present_fixture\",\n  \"verdict\": \"still-present\",\n  \"ensembl_vep\": [],"));
    // The fixed arguments reached both engines, with vep-rs's cache and Ensembl VEP's
    // --dir_cache root.
    let out = "work/still_present_fixture/out";
    let ensembl_args = arguments(&scenario, &format!("{out}/ensembl-vep/default.txt.args"));
    let vep_rs_args = arguments(&scenario, &format!("{out}/vep-rs/default.txt.args"));
    assert_eq!(after(&ensembl_args, "--dir_cache"), root);
    for flag in ["--cache_version", "--fork", "--minimal", "--no_intergenic"] {
        assert!(has(&ensembl_args, flag), "{flag}");
    }
    assert!(has(&vep_rs_args, "--json_cache") && has(&vep_rs_args, "--no_stats"));
    assert!(!has(&vep_rs_args, "--fork"));
    assert_eq!(after(&vep_rs_args, "--format"), "vcf");
    scenario.golden_file(
        "fixtures.ensembl-vep.args",
        &format!("{out}/ensembl-vep/default.txt.args"),
    );
    scenario.golden_file(
        "fixtures.vep-rs.args",
        &format!("{out}/vep-rs/default.txt.args"),
    );
    scenario.golden_file(
        "fixtures.ensembl-vep.command",
        &format!("{out}/ensembl-vep/command.txt"),
    );
    scenario.golden_file(
        "fixtures.vep-rs.command",
        &format!("{out}/vep-rs/command.txt"),
    );
    let written = scenario.read("work/still_present_fixture/stage/input/record.vcf");
    assert!(written.starts_with(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t99\tfx1\t"
    ));
    scenario.golden_file(
        "fixtures.input",
        "work/still_present_fixture/stage/input/record.vcf",
    );

    let mut expecting = args.to_vec();
    expecting.push("--expect-all-present");
    let ran = scenario.run(&expecting);
    assert_eq!(ran.code, 1);
    scenario.golden("fixtures_expect_all_present", &ran);
    expecting.extend(["--entry", "still_present_fixture"]);
    let ran = scenario.run(&expecting);
    assert_eq!(ran.code, 0);
    scenario.golden("fixtures_one_entry", &ran);
    let mut unknown = args.to_vec();
    unknown.extend(["--entry", "unknown_fixture"]);
    let ran = scenario.run(&unknown);
    assert_eq!(ran.code, 1);
    assert!(ran.stderr.contains("no reproduction for unknown_fixture"));
    scenario.golden("fixtures_unknown_entry", &ran);
}

/// A reproduction of a `filter` run whose Ensembl VEP stand-in writes the PATH it runs
/// with.
const PATH_RECORD: &str = r#"entry = "path_fixture"
page = "docs/intended-divergences.md"
verified_against = "116.2"
image = "ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da"
input_class = "unobserved"
input_class_basis = "well-formed under VCF 4.2; not seen in the measured datasets"

[input]
text = "21\t100\trec\tG\tA\t.\t.\t.\n"

[run]
tool = "filter"
args = ["-i", "{input}", "-o", "{out}"]
tools = ["samtools"]

[compare]
kind = "lines"
pattern = "."
mask_paths = false

[expected]
ensembl_vep = ["a"]
vep_rs = ["b"]
"#;

/// A tool directory goes first on the local engine's PATH, alone when PATH is unset; a
/// local Ensembl VEP runs `filter_vep` beside its `vep`, and vep-rs runs `vep filter`.
#[test]
fn a_tool_directory_goes_first_on_the_local_engines_path() {
    let scenario = Scenario::new("path");
    scenario.write("records/path_fixture.toml", PATH_RECORD.as_bytes());
    scenario.script("tools/vep", "exit 9");
    let write_path = "while [ $# -gt 0 ]; do if [ \"$1\" = \"-o\" ]; then printf '%s' \"$PATH\" > \"$2\"; fi; shift; done";
    scenario.script("tools/filter_vep", write_path);
    let vep_rs = scenario.script("rs/vep", "printf '%s\\n' \"$1\" > \"$PWD/default.txt\"");
    scenario.write("bin/samtools", b"");
    let bin = scenario.path("bin");
    let ensembl = scenario.path("tools/vep");
    let records = scenario.path("records");
    let work = scenario.path("work");
    let args = [
        "--ensembl-vep",
        &ensembl,
        "--vep",
        &vep_rs,
        "--reproductions",
        &records,
        "--work-dir",
        &work,
        "--tool-dir",
        &bin,
    ];
    let ran = scenario.run_with(&args, |command| {
        command.env("PATH", "/usr/bin:/bin");
    });
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert!(
        ran.stdout.starts_with(&format!(
            "path_fixture\tchanged\tEnsembl VEP ['{bin}:/usr/bin:/bin'], vep-rs ['filter']\n"
        )),
        "{}",
        ran.stdout
    );
    scenario.golden("tool_dir_path", &ran);
    assert_eq!(
        scenario.read("work/path_fixture/out/ensembl-vep/command.txt"),
        format!(
            "{}/tools/filter_vep -i {work}/path_fixture/stage/input/record.vcf -o {work}/path_fixture/out/ensembl-vep/default.txt\nexit 0\n",
            scenario.root
        )
    );
    let ran = scenario.run_with(&args, |command| {
        command.env_remove("PATH");
    });
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert!(
        ran.stdout.starts_with(&format!(
            "path_fixture\tchanged\tEnsembl VEP ['{bin}'], vep-rs ['filter']\n"
        )),
        "{}",
        ran.stdout
    );
    scenario.golden("tool_dir_path_unset", &ran);
    // Without the directory the entry cannot run.
    let ran = scenario.run(&args[..8]);
    scenario.golden("tool_dir_missing", &ran);
}

#[test]
fn an_engine_that_dies_reads_as_an_error_in_a_run() {
    let scenario = Scenario::new("dying");
    let dying = scenario.script(
        "dying",
        "echo 'ERROR: [cache] no transcript shards' >&2\nexit 2",
    );
    let vep_rs = scenario.stand_in("vep-rs", &fixture_outputs());
    let root = scenario.cache_root();
    let records = fixture_records();
    let work = scenario.path("work");
    let ran = scenario.run(&[
        "--ensembl-vep",
        &dying,
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &records,
        "--work-dir",
        &work,
        "--entry",
        "still_present_fixture",
        "--expect-all-present",
    ]);
    assert_eq!(ran.code, 1);
    assert!(ran.stdout.contains(
        "still_present_fixture\terror\tensembl-vep exited 2: ERROR: [cache] no transcript shards\n"
    ));
    assert!(ran.stdout.contains("| error | 1 |"));
    scenario.golden("dying_engine", &ran);

    let empty = scenario.mkdir("empty");
    let on_empty = [
        "--ensembl-vep",
        &dying,
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &empty,
        "--expect-all-present",
    ];
    let ran = scenario.run(&on_empty);
    assert_eq!(ran.code, 1);
    assert!(ran.stderr.contains("no reproduction under"));
    scenario.golden("no_reproduction", &ran);

    // A file that is not TOML: the line names the file, then says what the reader refused.
    let bad = scenario.write("empty/bad.toml", b"entry = [\n");
    let ran = scenario.run(&["--list", "--reproductions", &empty]);
    assert_eq!(ran.code, 1);
    assert!(
        ran.stderr
            .starts_with(&format!("ERROR: [reevaluate] {bad}: ")),
        "{}",
        ran.stderr
    );
    assert_eq!(ran.stderr.lines().count(), 1);
    assert_eq!(ran.stdout, "");
    std::fs::remove_file(bad).unwrap();
    let latin = scenario.write("empty/latin.toml", b"entry = \"caf\xe9\"\n");
    let ran = scenario.run(&["--list", "--reproductions", &empty]);
    assert_eq!(ran.code, 1);
    assert!(ran.stderr.contains("latin.toml: "));
    scenario.golden("not_utf8", &ran);
    std::fs::remove_file(latin).unwrap();

    let text = std::fs::read_to_string(format!("{records}/still_present_fixture.toml"))
        .unwrap()
        .replace(
            "entry = \"still_present_fixture\"",
            "entry = \"number_fixture\"",
        )
        .replace("json_cache = \"{stage}/json_cache\"", "json_cache = 5");
    scenario.write("empty/number_fixture.toml", text.as_bytes());
    let ran = scenario.run(&on_empty);
    assert_eq!(ran.code, 1);
    assert!(ran
        .stdout
        .contains("number_fixture\terror\tcache.json_cache is not a string\n"));
    scenario.golden("schema_error", &ran);
}

#[test]
fn a_cache_the_roots_lack_and_a_failing_prepare_step_are_errors() {
    let scenario = Scenario::new("cache");
    let ensembl = scenario.stand_in("ensembl-vep", &fixture_outputs());
    let vep_rs = scenario.stand_in("vep-rs", &fixture_outputs());
    let empty_root = scenario.mkdir("empty");
    let records = fixture_records();
    let work = scenario.path("work");
    let ran = scenario.run(&[
        "--ensembl-vep",
        &ensembl,
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &empty_root,
        "--reproductions",
        &records,
        "--work-dir",
        &work,
        "--entry",
        "closed_fixture",
    ]);
    assert_eq!(ran.code, 0);
    assert!(ran.stdout.contains(
        "closed_fixture\terror\tno Ensembl VEP cache homo_sapiens/116_GRCh38 with chromosomes 21"
    ));
    scenario.golden("cache_missing", &ran);

    let text = std::fs::read_to_string(format!("{records}/closed_fixture.toml"))
        .unwrap()
        .replace(
            "entry = \"closed_fixture\"",
            "entry = \"broken_fixture\"\nprepare = [\"false\"]",
        );
    scenario.write("broken/broken_fixture.toml", text.as_bytes());
    let broken = scenario.path("broken");
    let root = scenario.cache_root();
    let work = scenario.path("work2");
    let ran = scenario.run(&[
        "--ensembl-vep",
        &ensembl,
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &broken,
        "--work-dir",
        &work,
        "--expect-all-present",
    ]);
    assert_eq!(ran.code, 1);
    assert!(ran
        .stdout
        .contains("broken_fixture\terror\tprepare step failed (1): false\n"));
    scenario.golden("prepare_fails", &ran);
}

/// `cache.perl_cache` names a root the prepare steps write, so it is resolved after them,
/// with the roots' cache in PERL_CACHE; IMAGE is empty without `--image`.
#[test]
fn a_prepared_cache_is_read_after_the_prepare_steps() {
    let scenario = Scenario::new("prepared");
    let ensembl = scenario.stand_in("ensembl-vep", &fixture_outputs());
    let vep_rs = scenario.stand_in("vep-rs", &fixture_outputs());
    let root = scenario.cache_root();
    let step = r#"mkdir -p "$STAGE/perl-cache/homo_sapiens" && cp -al "$PERL_CACHE" "$STAGE/perl-cache/homo_sapiens/116_GRCh38" && echo "$PERL_CACHE" > "$STAGE/root.txt" && printf "%s" "${IMAGE-unset}" > "$STAGE/image.txt""#;
    let text = std::fs::read_to_string(format!("{}/still_present_fixture.toml", fixture_records()))
        .unwrap()
        .replace(
            "entry = \"still_present_fixture\"",
            &format!("entry = \"prepared_fixture\"\nprepare = ['{step}']"),
        )
        .replace(
            "json_cache = \"{stage}/json_cache\"",
            "json_cache = \"{stage}/json_cache\"\nperl_cache = \"{stage}/perl-cache\"",
        );
    scenario.write("repros/prepared_fixture.toml", text.as_bytes());
    let repros = scenario.path("repros");
    let work = scenario.path("work");
    let ran = scenario.run(&[
        "--ensembl-vep",
        &ensembl,
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &repros,
        "--work-dir",
        &work,
        "--expect-all-present",
    ]);
    assert_eq!(ran.code, 0, "{}{}", ran.stdout, ran.stderr);
    assert!(ran.stdout.contains("prepared_fixture\tstill-present\n"));
    scenario.golden("prepared_cache", &ran);
    assert_eq!(
        scenario.read("work/prepared_fixture/stage/root.txt").trim(),
        format!("{root}/homo_sapiens/116_GRCh38")
    );
    assert_eq!(scenario.read("work/prepared_fixture/stage/image.txt"), "");
    let args = arguments(
        &scenario,
        "work/prepared_fixture/out/ensembl-vep/default.txt.args",
    );
    assert_eq!(
        after(&args, "--dir_cache"),
        format!("{work}/prepared_fixture/stage/perl-cache")
    );
}

/// The reproductions of `tests/fixtures/kinds/`, one a compare kind and its keys, through
/// stand-in engines that write the outputs beside them: every kind's reader, the
/// normalising of its values and each of the four verdicts.
#[test]
fn every_compare_kind_reads_its_value_from_a_run() {
    let scenario = Scenario::new("kinds");
    let outputs = format!("{FIXTURES}/kinds/outputs");
    let ensembl = scenario.entry_stand_in("ensembl-vep", &outputs);
    let vep_rs = scenario.entry_stand_in("vep-rs", &outputs);
    for sibling in ["filter_vep", "variant_recoder"] {
        std::fs::copy(
            &ensembl,
            scenario.path(&format!("stand-in/ensembl-vep/{sibling}")),
        )
        .unwrap();
    }
    let root = scenario.cache_root();
    let records = format!("{FIXTURES}/kinds/reproductions");
    let work = scenario.path("work");
    let results = scenario.path("results.json");
    let args = [
        "--ensembl-vep",
        &ensembl,
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &records,
        "--work-dir",
        &work,
        "--json",
        &results,
    ];
    let ran = scenario.run(&args);
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    scenario.golden("kinds", &ran);
    scenario.golden_file("kinds.json", "results.json");
    for verdict in ["still-present", "closed", "changed", "error"] {
        let none = format!("| {verdict} | 0 |");
        assert!(!ran.stdout.contains(&none), "no entry reads {verdict}");
    }
    for (name, relative) in [
        (
            "kinds.stdin",
            "work/stdout_kind/out/ensembl-vep/default.txt.stdin",
        ),
        (
            "kinds.stdin.args",
            "work/stdout_kind/out/vep-rs/default.txt.args",
        ),
        (
            "kinds.tab.args",
            "work/lines_unmasked_kind/out/ensembl-vep/tab.txt.args",
        ),
        (
            "kinds.vcf-line.args",
            "work/vcf_line_kind/out/ensembl-vep/vcf.vcf.args",
        ),
        (
            "kinds.json.args",
            "work/json_kind/out/vep-rs/json.jsonl.args",
        ),
        (
            "kinds.tool.ensembl-vep.command",
            "work/tool_kind/out/ensembl-vep/command.txt",
        ),
        (
            "kinds.tool.vep-rs.command",
            "work/tool_kind/out/vep-rs/command.txt",
        ),
        ("kinds.tool.input", "work/tool_kind/stage/input/ids.txt"),
        (
            "kinds.exit.command",
            "work/exit_kind/out/ensembl-vep/command.txt",
        ),
        (
            "kinds.stderr.stderr",
            "work/stderr_kind/out/ensembl-vep/stderr.txt",
        ),
        (
            "kinds.stdout.stdout",
            "work/stdout_kind/out/ensembl-vep/stdout.txt",
        ),
    ] {
        scenario.golden_file(name, relative);
    }
    // With several jobs the entries end in any order and the report is the same.
    let mut parallel = args.to_vec();
    parallel.extend(["--jobs", "4", "--expect-all-present"]);
    let again = scenario.run(&parallel);
    assert_eq!(again.code, 1);
    assert_eq!(sorted_results(&again.stdout), sorted_results(&ran.stdout));
    scenario.golden_file("kinds.json", "results.json");
}

/// A stand-in for the container runtime: writes its arguments beside the output and runs
/// the stand-in engine on what follows the image.
fn docker_stand_in(scenario: &Scenario, engine: &str) -> String {
    let body = format!(
        r#"printf '%s\n' "$@" >> "$PWD/docker.args"
[ "$1" = run ] || exit 0
shift
while [ $# -gt 0 ]; do
	case "$1" in
	--name | --user | -e | -v | -w) shift 2 ;;
	--rm | -i) shift ;;
	*) break ;;
	esac
done
shift 2
exec "{engine}" "$@""#
    );
    scenario.script("docker-bin/docker", &body);
    scenario.path("docker-bin")
}

/// A text with each container name, `reevaluate-<entry>-<process id>-<eight hex digits>`,
/// written without the two parts that differ from run to run.
fn mask_run(text: &str) -> String {
    let mut out = String::new();
    let mut rest = text;
    while let Some(start) = rest.find("reevaluate-") {
        let at_word_start =
            start == 0 || matches!(rest.as_bytes()[start - 1], b' ' | b'\'' | b'\n');
        let end = rest[start..]
            .find(['\'', ' ', '\n'])
            .map_or(rest.len(), |length| start + length);
        let parts: Vec<&str> = rest[start..end].rsplitn(3, '-').collect();
        let is_name = at_word_start
            && parts.len() == 3
            && parts[0].len() == 8
            && parts[1].bytes().all(|b| b.is_ascii_digit());
        if is_name {
            out.push_str(&rest[..start]);
            out.push_str(parts[2]);
            out.push_str("-PID-HEX");
        } else {
            out.push_str(&rest[..end]);
        }
        rest = &rest[end..];
    }
    out.push_str(rest);
    out
}

fn mask_user(text: &str, user: &str) -> String {
    text.replace(user, "UID:GID")
}

#[test]
fn the_image_route_runs_the_container_command() {
    use std::os::unix::fs::MetadataExt;
    let scenario = Scenario::new("image");
    let outputs = format!("{FIXTURES}/kinds/outputs");
    let ensembl = scenario.entry_stand_in("ensembl-vep", &outputs);
    let vep_rs = scenario.entry_stand_in("vep-rs", &outputs);
    let docker_dir = docker_stand_in(&scenario, &ensembl);
    let path = format!("{docker_dir}:/usr/bin:/bin");
    let root = scenario.cache_root();
    let records = format!("{FIXTURES}/kinds/reproductions");
    let work = scenario.path("work");
    let plugins = scenario.mkdir("plugins");
    let own = std::fs::metadata(&scenario.root).unwrap();
    let user = format!("{}:{}", own.uid(), own.gid());
    let args = [
        "--image",
        "ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da",
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &root,
        "--plugins-dir",
        &plugins,
        "--reproductions",
        &records,
        "--work-dir",
        &work,
        "--entry",
        "field_kind",
        "--entry",
        "stdout_kind",
        "--entry",
        "tool_kind",
        "--entry",
        "exit_kind",
    ];
    let ran = scenario.run_with(&args, |command| {
        command.env("PATH", &path);
    });
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    scenario.golden("image", &ran);
    for entry in ["field_kind", "stdout_kind", "tool_kind", "exit_kind"] {
        let command = scenario.read(&format!("work/{entry}/out/ensembl-vep/command.txt"));
        scenario.golden_text(
            &format!("image.{entry}.command"),
            &mask_user(&mask_run(&scenario.placeholders(&command)), &user),
        );
        let recorded = scenario.read(&format!("work/{entry}/out/ensembl-vep/docker.args"));
        scenario.golden_text(
            &format!("image.{entry}.docker-args"),
            &mask_user(&mask_run(&scenario.placeholders(&recorded)), &user),
        );
    }
    // The cache root, the plugin modules and the work directory written from the calling
    // directory: each mount names the same host path.
    let mut relative = args.to_vec();
    for (option, directory) in [
        ("--perl-cache-dir", "perl"),
        ("--plugins-dir", "plugins"),
        ("--work-dir", "work"),
    ] {
        let index = relative
            .iter()
            .position(|arg| *arg == option)
            .expect("the option");
        relative[index + 1] = directory;
    }
    let ran = scenario.run_with(&relative, |command| {
        command.env("PATH", &path);
    });
    scenario.golden("image", &ran);
    for entry in ["field_kind", "stdout_kind", "tool_kind", "exit_kind"] {
        let recorded = scenario.read(&format!("work/{entry}/out/ensembl-vep/docker.args"));
        scenario.golden_text(
            &format!("image.{entry}.docker-args"),
            &mask_user(&mask_run(&scenario.placeholders(&recorded)), &user),
        );
    }
    // Without a container runtime on PATH nothing runs.
    let empty = scenario.mkdir("no-docker");
    let ran = scenario.run_with(&args, |command| {
        command.env("PATH", &empty);
    });
    assert_eq!(ran.code, 1);
    scenario.golden("image_without_docker", &ran);
}

#[test]
fn an_engine_run_is_ended_at_the_time_allowed() {
    let scenario = Scenario::new("timeout");
    let sleeper = scenario.script("sleeper/vep", "sleep 30");
    let root = scenario.cache_root();
    let records = fixture_records();
    let work = scenario.path("work");
    let args = [
        "--ensembl-vep",
        &sleeper,
        "--vep",
        &sleeper,
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &records,
        "--work-dir",
        &work,
        "--timeout",
        "1",
        "--jobs",
        "2",
    ];
    let started = std::time::Instant::now();
    let ran = scenario.run(&args);
    assert!(started.elapsed().as_secs() < 20, "the runs were not ended");
    assert_eq!(ran.code, 0);
    scenario.golden_text("timeout.stdout", &sorted_results(&ran.stdout));
    assert_eq!(ran.stderr, "");
}

#[test]
fn outputs_are_kept_when_asked() {
    let scenario = Scenario::new("keep");
    let ensembl = scenario.stand_in("ensembl-vep", &fixture_outputs());
    let vep_rs = scenario.stand_in("vep-rs", &fixture_outputs());
    let root = scenario.cache_root();
    let records = fixture_records();
    let temporary = scenario.mkdir("tmp");
    let args = [
        "--ensembl-vep",
        &ensembl,
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &records,
    ];
    let entries = |directory: &str| -> Vec<String> {
        std::fs::read_dir(directory)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect()
    };
    // A temporary work directory is gone when the run ends.
    let ran = scenario.run_with(&args, |command| {
        command.env("TMPDIR", &temporary);
    });
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert!(!ran.stdout.contains("outputs under"));
    assert!(entries(&temporary).is_empty());
    // With --keep it stays and the last line names it.
    let mut keeping = args.to_vec();
    keeping.push("--keep");
    let ran = scenario.run_with(&keeping, |command| {
        command.env("TMPDIR", &temporary);
    });
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    let kept = entries(&temporary);
    assert_eq!(kept.len(), 1, "{kept:?}");
    assert!(kept[0].starts_with("reevaluate-") && kept[0].len() == "reevaluate-".len() + 8);
    assert!(ran
        .stdout
        .ends_with(&format!("outputs under {temporary}/{}\n", kept[0])));
    assert!(std::path::Path::new(&format!(
        "{temporary}/{}/closed_fixture/out/vep-rs/default.txt",
        kept[0]
    ))
    .is_file());
}

/// Each engine runs in its entry's output directory, so a path the caller wrote from the
/// directory the tool was called in must reach it whole: the run is the one the same
/// paths give when written from the root.
#[test]
fn paths_written_from_the_calling_directory_run_as_the_same_paths_from_the_root() {
    let scenario = Scenario::new("relative");
    scenario.stand_in("ensembl-vep", &fixture_outputs());
    scenario.stand_in("vep-rs", &fixture_outputs());
    scenario.cache_root();
    let records = fixture_records();
    let ran = scenario.run(&[
        "--ensembl-vep",
        "./fake-ensembl-vep",
        "--vep",
        "./fake-vep-rs",
        "--perl-cache-dir",
        "perl",
        "--reproductions",
        &records,
        "--work-dir",
        "work",
        "--json",
        "results.json",
    ]);
    scenario.golden("fixtures", &ran);
    scenario.golden_file("fixtures.json", "results.json");
    let out = "work/still_present_fixture/out";
    for engine in ["ensembl-vep", "vep-rs"] {
        scenario.golden_file(
            &format!("fixtures.{engine}.args"),
            &format!("{out}/{engine}/default.txt.args"),
        );
        scenario.golden_file(
            &format!("fixtures.{engine}.command"),
            &format!("{out}/{engine}/command.txt"),
        );
    }

    // A binary under a directory of the calling directory, and one above it.
    let nested = scenario.mkdir("nested/deeper");
    for binary in ["nested/fake-vep-rs", "../../fake-vep-rs"] {
        let directory = if binary.starts_with("..") {
            nested.as_str()
        } else {
            scenario.root.as_str()
        };
        if !binary.starts_with("..") {
            std::fs::copy(scenario.path("fake-vep-rs"), scenario.path(binary))
                .expect("a copy of the stand-in");
        }
        let work = scenario.path(&format!("work-{}", binary.len()));
        let ensembl = scenario.path("fake-ensembl-vep");
        let root = scenario.path("perl");
        let args = [
            "--ensembl-vep",
            &ensembl,
            "--vep",
            binary,
            "--perl-cache-dir",
            &root,
            "--reproductions",
            &records,
            "--work-dir",
            &work,
        ];
        let ran = scenario.run_with(&args, |command| {
            command.current_dir(directory);
        });
        assert_eq!(ran.code, 0, "{}", ran.stderr);
        assert!(
            ran.stdout.contains("| error | 0 |\n| total | 2 |"),
            "--vep {binary}: {}",
            ran.stdout
        );
    }

    // The plugin modules of a record that loads a plugin.
    let record = std::fs::read_to_string(format!("{records}/closed_fixture.toml"))
        .expect("the fixture record")
        .replace("[\"--hgvs\"]", "[\"--hgvs\", \"--plugin\", \"Made\"]");
    scenario.write("records/closed_fixture.toml", record.as_bytes());
    scenario.mkdir("modules");
    let ran = scenario.run(&[
        "--ensembl-vep",
        "./fake-ensembl-vep",
        "--vep",
        "./fake-vep-rs",
        "--perl-cache-dir",
        "perl",
        "--plugins-dir",
        "modules",
        "--reproductions",
        "records",
        "--work-dir",
        "work-plugin",
    ]);
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    let given = arguments(
        &scenario,
        "work-plugin/closed_fixture/out/ensembl-vep/default.txt.args",
    );
    assert_eq!(after(&given, "--dir_plugins"), scenario.path("modules"));

    // A script behind its interpreter is the last word of the command.
    let ran = scenario.run(&[
        "--ensembl-vep",
        "sh ./fake-ensembl-vep",
        "--vep",
        "./fake-vep-rs",
        "--perl-cache-dir",
        "perl",
        "--reproductions",
        &records,
        "--work-dir",
        "work",
    ]);
    scenario.golden("fixtures", &ran);
    let command = scenario.read(&format!("{out}/ensembl-vep/command.txt"));
    assert!(
        command.starts_with(&format!("sh {}/fake-ensembl-vep -i ", scenario.root)),
        "{command}"
    );
}

/// A binary named without a directory is a command name, found on `PATH`, also where the
/// calling directory holds a file of that name.
#[test]
fn a_binary_named_without_a_directory_is_looked_up_on_the_path() {
    let scenario = Scenario::new("path-lookup");
    let ensembl = scenario.stand_in("ensembl-vep", &fixture_outputs());
    let vep_rs = scenario.stand_in("vep-rs", &fixture_outputs());
    let root = scenario.cache_root();
    let records = fixture_records();
    let work = scenario.path("work");
    let on_path = scenario.mkdir("on-path");
    std::fs::rename(&vep_rs, format!("{on_path}/vep-stand-in")).expect("the stand-in moved");
    // A file of the same name in the calling directory, which writes no output.
    scenario.script("vep-stand-in", "exit 1");
    let args = [
        "--ensembl-vep",
        &ensembl,
        "--vep",
        "vep-stand-in",
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &records,
        "--work-dir",
        &work,
    ];
    let path = std::env::var("PATH").expect("a PATH");
    let ran = scenario.run_with(&args, |command| {
        command.env("PATH", format!("{on_path}:{path}"));
    });
    assert_eq!(ran.code, 0, "{}", ran.stderr);
    assert!(
        ran.stdout.contains("| error | 0 |\n| total | 2 |"),
        "{}",
        ran.stdout
    );
    let command = scenario.read("work/still_present_fixture/out/vep-rs/command.txt");
    assert!(command.starts_with("vep-stand-in -i "), "{command}");
}

/// A reproduction of two runs: the record's own under one flag and a further one, of
/// another record under another flag, in its `also` list.
const TWO_RUNS_RECORD: &str = r#"entry = "two_runs_fixture"
page = "docs/intended-divergences.md"
verified_against = "116.3"
image = "ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
input_class = "unobserved"
input_class_basis = "well-formed under VCF 4.2; not seen in the measured datasets"

[input]
text = "21\t99\tfx1\tGT\tGA,GTC\t.\t.\t.\n"

[cache]
assembly = "GRCh38"
json_cache = "{stage}/json_cache"

[run]
flags = ["--minimal"]

[compare]
kind = "rows"
where = { record = "fx1" }
columns = ["Location", "Feature", "DISTANCE"]

[expected]
ensembl_vep = []
vep_rs = ["21:100 | ENST1 | 5000"]

[[also]]

[also.input]
text = "21\t199\tfx2\tC\tT\t.\t.\t.\n"

[also.run]
flags = ["--hgvs"]

[also.compare]
kind = "field"
field = "HGVSp"
where = { record = "fx2" }

[also.expected]
ensembl_vep = ["ENSP2:p.Leu4%3D"]
vep_rs = ["ENSP2:p.Leu4="]
"#;

const DEFAULT_HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

/// Each run of a record is made in its own directory with its own input and flags, and the
/// entry reads `still-present` only while every run reads its recorded values. A further
/// run works under `also-1` of the entry's directory, by which the stand-in engines know
/// it.
#[test]
fn a_records_further_run_is_made_and_judged_beside_its_own() {
    let scenario = Scenario::new("also");
    let outputs = scenario.path("outputs");
    let ensembl = scenario.entry_stand_in("ensembl-vep", &outputs);
    let vep_rs = scenario.entry_stand_in("vep-rs", &outputs);
    let own_row = "fx1\t21:100\tA\tENSG1\tENST1\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=5000\n";
    let further_row = |hgvsp: &str| {
        format!("fx2\t21:200\tT\tENSG2\tENST2\tTranscript\tsynonymous_variant\t10\t10\t4\tL\tctA/ctG\t-\tIMPACT=LOW;HGVSp={hgvsp}\n")
    };
    let output = |relative: &str, rows: &str| {
        scenario.write(
            &format!("outputs/{relative}/default.txt"),
            format!("{DEFAULT_HEADER}{rows}").as_bytes(),
        );
    };
    output("two_runs_fixture/ensembl-vep", "");
    output("two_runs_fixture/vep-rs", own_row);
    output("also-1/ensembl-vep", &further_row("ENSP2:p.Leu4%3D"));
    output("also-1/vep-rs", &further_row("ENSP2:p.Leu4="));
    scenario.write("records/two_runs_fixture.toml", TWO_RUNS_RECORD.as_bytes());
    let root = scenario.cache_root();
    let records = scenario.path("records");
    let work = scenario.path("work");
    let results = scenario.path("results.json");
    let args = [
        "--ensembl-vep",
        &ensembl,
        "--vep",
        &vep_rs,
        "--perl-cache-dir",
        &root,
        "--reproductions",
        &records,
        "--work-dir",
        &work,
        "--json",
        &results,
        "--expect-all-present",
    ];
    let ran = scenario.run(&args);
    assert_eq!(ran.code, 0, "{}{}", ran.stdout, ran.stderr);
    assert!(ran.stdout.contains("two_runs_fixture\tstill-present\n"));
    assert!(ran.stdout.contains("| total | 1 |"));
    // The entry's values are its runs' values, the record's own run first.
    let recorded = scenario.read("results.json");
    assert!(recorded.contains("\"ensembl_vep\": [\n   \"ENSP2:p.Leu4%3D\"\n  ]"));
    assert!(recorded
        .contains("\"vep_rs\": [\n   \"21:100 | ENST1 | 5000\",\n   \"ENSP2:p.Leu4=\"\n  ]"));
    // Each run was given its own flags and its own input, in its own directory.
    for (directory, flag, other, record) in [
        (
            "work/two_runs_fixture",
            "--minimal",
            "--hgvs",
            "21\t99\tfx1\t",
        ),
        (
            "work/two_runs_fixture/also-1",
            "--hgvs",
            "--minimal",
            "21\t199\tfx2\t",
        ),
    ] {
        for engine in ["ensembl-vep", "vep-rs"] {
            let given = arguments(
                &scenario,
                &format!("{directory}/out/{engine}/default.txt.args"),
            );
            assert!(has(&given, flag), "{directory} {engine} {flag}");
            assert!(!has(&given, other), "{directory} {engine} {other}");
        }
        let written = scenario.read(&format!("{directory}/stage/input/record.vcf"));
        assert!(written.contains(record), "{directory}: {written}");
    }

    // The further run's engines agree: the entry differs in another way than recorded.
    output("also-1/vep-rs", &further_row("ENSP2:p.Leu4%3D"));
    let ran = scenario.run(&args);
    assert_eq!(ran.code, 1);
    assert!(
        ran.stdout.contains(
            "two_runs_fixture\tchanged\talso[1]: Ensembl VEP ['ENSP2:p.Leu4%3D'], vep-rs ['ENSP2:p.Leu4%3D']\n"
        ),
        "{}",
        ran.stdout
    );

    // A further run that cannot be made is the entry's error.
    std::fs::remove_file(scenario.path("outputs/also-1/vep-rs/default.txt")).unwrap();
    let ran = scenario.run(&args);
    assert_eq!(ran.code, 1);
    assert!(
        ran.stdout.contains("two_runs_fixture\terror\talso[1]: "),
        "{}",
        ran.stdout
    );
}
