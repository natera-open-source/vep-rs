// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::repo::repo_root;
use crate::testing::{repro, Scratch, DELETE};
use std::os::unix::fs::{MetadataExt, PermissionsExt};

fn container_paths<'a>(
    repro: &'a Reproduction,
    scratch: &Scratch,
    plugins: Option<&str>,
) -> Paths<'a> {
    Paths::new(
        repro,
        &scratch.path("stage"),
        &scratch.path("out"),
        plugins,
        true,
        &[],
        &repo_root(),
    )
}

fn message<T: std::fmt::Debug>(result: Run<T>) -> String {
    result.expect_err("a refusal").0
}

fn words(list: &[&str]) -> Vec<String> {
    list.iter().map(|w| w.to_string()).collect()
}

/// An executable shell script.
fn script(scratch: &Scratch, relative: &str, body: &str) -> String {
    let path = scratch.write(relative, format!("#!/bin/sh\n{body}\n").as_bytes());
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
    path
}

#[test]
fn write_input_keeps_a_header_the_text_brings_and_adds_samples() {
    let scratch = Scratch::new("input");
    let own_header = repro(&[(
        "input.text",
        r###""##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t1\n""###,
    )]);
    let written = write_input(&own_header, &container_paths(&own_header, &scratch, None)).unwrap();
    assert_eq!(
        std::fs::read_to_string(written).unwrap(),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t1\n"
    );
    let sampled = repro(&[
        ("input.text", r#""21\t1\tr\tA\tG\t.\t.\t.\tGT\t0/1\t1/1""#),
        ("input.samples", "['S1', 'S2']"),
        (
            "input.header",
            r#"['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']"#,
        ),
    ]);
    let sampled_paths = container_paths(&sampled, &scratch, None);
    let written = write_input(&sampled, &sampled_paths).unwrap();
    assert_eq!(written, scratch.path("stage/input/record.vcf"));
    assert_eq!(
        std::fs::read_to_string(&written).unwrap(),
        "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n21\t1\tr\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\n"
    );
    assert_eq!(input_chromosomes(&sampled, &written), ["21"]);
    let ensembl = repro(&[
        ("input.text", "'21 1 1 A/G +'"),
        ("input.format", "'ensembl'"),
    ]);
    let written = write_input(&ensembl, &sampled_paths).unwrap();
    assert_eq!(pypath::name(&written), "record.txt");
    assert!(input_chromosomes(&ensembl, &written).is_empty());
    assert_eq!(
        input_chromosomes(
            &repro(&[("cache.chromosomes", "['MT']")]),
            &scratch.path("none")
        ),
        ["MT"]
    );
    // A named file is resolved, not written.
    let named = repro(&[
        ("input.text", DELETE),
        ("input.path", "'{stage}/made/by/a/step.vcf'"),
    ]);
    assert_eq!(
        write_input(&named, &container_paths(&named, &scratch, None)).unwrap(),
        scratch.path("stage/made/by/a/step.vcf")
    );
    // A first column with no tab or blank after it keeps its line end.
    let bare = scratch.write("bare.vcf", b"#h\n\nMT\n21 5\n21\t7\n");
    assert_eq!(input_chromosomes(&repro(&[]), &bare), ["MT\n", "21"]);
}

#[test]
fn find_perl_cache_takes_the_first_root_holding_the_chromosomes() {
    let scratch = Scratch::new("cache");
    let first = scratch.path("first");
    let second = scratch.path("second");
    scratch.mkdir("first/homo_sapiens/116_GRCh38/21");
    scratch.mkdir("second/homo_sapiens/116_GRCh38/21");
    scratch.mkdir("second/homo_sapiens/116_GRCh38/MT");
    scratch.mkdir("second/homo_sapiens_merged/116_GRCh38/21");
    let roots = [first.clone(), second.clone()];
    let find = |edits: &[(&str, &str)], chroms: &[&str]| {
        find_perl_cache(&repro(edits), &roots, &words(chroms))
    };
    assert_eq!(
        find(&[], &["21"]).unwrap(),
        Some(format!("{first}/homo_sapiens/116_GRCh38"))
    );
    assert_eq!(
        find(&[], &["21", "MT"]).unwrap(),
        Some(format!("{second}/homo_sapiens/116_GRCh38"))
    );
    assert_eq!(
        find(&[("cache.source", "'merged'")], &["21"]).unwrap(),
        Some(format!("{second}/homo_sapiens_merged/116_GRCh38"))
    );
    assert_eq!(
        message(find(&[("cache.source", "'refseq'")], &["21"])),
        format!("no Ensembl VEP cache homo_sapiens_refseq/116_GRCh38 with chromosomes 21 under the --perl-cache-dir roots (tried {first}/homo_sapiens_refseq/116_GRCh38, {second}/homo_sapiens_refseq/116_GRCh38)")
    );
    assert_eq!(
        message(find_perl_cache(&repro(&[]), &[], &[])),
        "no Ensembl VEP cache homo_sapiens/116_GRCh38 under the --perl-cache-dir roots (tried none)"
    );
    // A reproduction that prepares its own cache may run with no root cache at all.
    let staged_repro = repro(&[
        ("cache.perl_cache", "'{stage}/perl-cache'"),
        ("cache.source", "'refseq'"),
    ]);
    assert_eq!(
        find_perl_cache(&staged_repro, &roots, &words(&["21"])).unwrap(),
        None
    );
    let staged_paths = container_paths(&staged_repro, &scratch, None);
    assert_eq!(
        message(prepared_perl_cache(&staged_repro, &staged_paths)),
        "cache.perl_cache {stage}/perl-cache holds no homo_sapiens_refseq/116_GRCh38"
    );
    let staged = scratch.mkdir("stage/perl-cache/homo_sapiens_refseq/116_GRCh38");
    assert_eq!(
        prepared_perl_cache(&staged_repro, &staged_paths).unwrap(),
        staged
    );
    assert_eq!(
        cache_relative(&repro(&[
            ("cache.species", "'mus_musculus'"),
            ("cache.assembly", "'GRCm39'"),
            ("cache.version", "'116'")
        ])),
        "mus_musculus/116_GRCm39"
    );
}

#[test]
fn substitute_fills_placeholders_per_engine_and_mounts_what_a_word_names() {
    let scratch = Scratch::new("substitute");
    let record = repro(&[(
        "cache.json_cache",
        "'tests/golden/116/GRCh38-custom/json_cache'",
    )]);
    let plugins = scratch.path("plugins");
    let mut paths = container_paths(&record, &scratch, Some(&plugins));
    let corpus = format!("{}/tests/golden/116/GRCh38-custom", repo_root());
    let input = scratch.path("in.vcf");
    let output = scratch.path("out/default.txt");
    let fill = |paths: &mut Paths, word: &str, engine: &str| {
        paths.substitute(word, engine, &input, &output)
    };
    assert_eq!(
        fill(&mut paths, "file={corpus}/custom/sites.vcf.gz", "vep-rs").unwrap(),
        format!("file={corpus}/custom/sites.vcf.gz")
    );
    assert_eq!(
        fill(
            &mut paths,
            "file={corpus}/custom/sites.vcf.gz",
            "ensembl-vep"
        )
        .unwrap(),
        "file=/work/corpus/custom/sites.vcf.gz"
    );
    assert_eq!(
        paths.mounts,
        [("/work/corpus".to_string(), corpus.clone(), false)]
    );
    assert_eq!(
        fill(&mut paths, "{plugins}", "vep-rs").unwrap(),
        format!("{corpus}/plugin_data")
    );
    assert_eq!(
        fill(&mut paths, "{plugins}", "ensembl-vep").unwrap(),
        "/work/plugins"
    );
    assert_eq!(
        fill(&mut paths, "{stage}/x", "ensembl-vep").unwrap(),
        "/work/stage/x"
    );
    assert_eq!(
        fill(&mut paths, "{input} {out}", "ensembl-vep").unwrap(),
        "/work/in/in.vcf /work/out/default.txt"
    );
    assert_eq!(
        fill(&mut paths, "--plain", "ensembl-vep").unwrap(),
        "--plain"
    );
    assert_eq!(
        paths.mounts,
        [
            ("/work/corpus".to_string(), corpus.clone(), false),
            ("/work/plugins".to_string(), plugins.clone(), false),
            ("/work/stage".to_string(), scratch.path("stage"), true),
        ]
    );
    assert_eq!(
        message(fill(&mut paths, "{loftee}/LoF.pm", "ensembl-vep")),
        "{loftee} in '{loftee}/LoF.pm' names no directory: pass --extra-dir loftee=<path>"
    );
    let loftee = scratch.path("loftee");
    let mut extras = Paths::new(
        &record,
        &scratch.path("stage"),
        &scratch.path("out"),
        Some(&plugins),
        true,
        &[("loftee".to_string(), loftee.clone())],
        &repo_root(),
    );
    assert_eq!(
        fill(&mut extras, "loftee_path:{loftee}", "vep-rs").unwrap(),
        format!("loftee_path:{loftee}")
    );
    assert_eq!(
        fill(&mut extras, "loftee_path:{loftee}", "ensembl-vep").unwrap(),
        "loftee_path:/work/extra/loftee"
    );
    assert_eq!(
        extras.mounts,
        [("/work/extra/loftee".to_string(), loftee, false)]
    );
    let bare = repro(&[
        ("run.tool", "'filter'"),
        ("run.args", "['x']"),
        ("cache.json_cache", DELETE),
        ("cache.assembly", DELETE),
    ]);
    let mut without_corpus = container_paths(&bare, &scratch, None);
    assert_eq!(
        message(fill(&mut without_corpus, "{corpus}/x", "vep-rs")),
        "{corpus} in '{corpus}/x' names no directory: the reproduction has no corpus"
    );
    assert_eq!(
        message(fill(&mut without_corpus, "{plugins}", "ensembl-vep")),
        "{plugins} in '{plugins}' needs --plugins-dir"
    );
    // A local Ensembl VEP reads the host paths and mounts nothing.
    let mut local = Paths::new(
        &record,
        &scratch.path("stage"),
        &scratch.path("out"),
        Some(&plugins),
        false,
        &[],
        &repo_root(),
    );
    assert_eq!(
        fill(&mut local, "{plugins}/x {repo}/y", "ensembl-vep").unwrap(),
        format!("{plugins}/x {}/y", repo_root())
    );
    assert!(local.mounts.is_empty());
}

#[test]
fn tool_directories_come_from_the_tool_dir_option() {
    let scratch = Scratch::new("tools");
    scratch.write("tools/samtools", b"");
    let tools = scratch.path("tools");
    let dirs = [scratch.mkdir("empty"), tools.clone()];
    assert_eq!(
        tool_directories(&repro(&[("run.tools", "['samtools']")]), &dirs).unwrap(),
        [tools]
    );
    assert!(tool_directories(&repro(&[]), &dirs).unwrap().is_empty());
    assert_eq!(
        message(tool_directories(
            &repro(&[("run.tools", "['bcftools']")]),
            &dirs
        )),
        "the reproduction needs the tool bcftools: pass --tool-dir <directory holding it>"
    );
}

#[test]
fn docker_command_mounts_read_only_unless_writable_and_passes_stdin() {
    let record = repro(&[]);
    let mut paths = Paths::new(&record, "/stage", "/out", None, true, &[], &repo_root());
    paths.mount("/host/cache", "/work/cache/homo_sapiens/116_GRCh38", false);
    paths.mount("/host/out", "/work/out", true);
    let seed = [("PERL_HASH_SEED", "0".to_string())];
    let command = docker_command(
        "img@sha256:0",
        &paths,
        "vep",
        &words(&["-i", "x"]),
        &seed,
        true,
        Some("reevaluate-x-1"),
    );
    let own = std::fs::metadata("/proc/self").unwrap();
    let user = format!("{}:{}", own.uid(), own.gid());
    assert_eq!(
        command,
        [
            "docker",
            "run",
            "--rm",
            "--name",
            "reevaluate-x-1",
            "--user",
            &user,
            "-e",
            "HOME=/tmp",
            "-i",
            "-e",
            "PERL_HASH_SEED=0",
            "-v",
            "/host/cache:/work/cache/homo_sapiens/116_GRCh38:ro",
            "-v",
            "/host/out:/work/out",
            "-w",
            "/work/out",
            "img@sha256:0",
            "vep",
            "-i",
            "x"
        ]
    );
    let quiet = docker_command("img@sha256:0", &paths, "vep", &[], &[], false, None);
    assert!(!quiet[..8].contains(&"-i".to_string()));
    assert_eq!(quiet[..5], ["docker", "run", "--rm", "--user", &user]);
    // A mount named twice keeps its last host and stays writable once it was.
    paths.mount("/other/out", "/work/out", false);
    assert!(paths
        .mounts
        .contains(&("/work/out".to_string(), "/other/out".to_string(), true)));
}

#[test]
fn engine_args_follow_the_fixed_argument_lists() {
    let scratch = Scratch::new("args");
    let repo = repo_root();
    let record = repro(&[
        (
            "cache.fasta",
            "'tests/golden/116/GRCh37-hgvs/reference.fa.gz'",
        ),
        ("cache.assembly", "'GRCh37'"),
        (
            "cache.json_cache",
            "'tests/golden/116/GRCh37-hgvs/json_cache'",
        ),
        ("run.flags", "['--hgvs', '--plugin', 'pLI']"),
        ("run.stats", "true"),
        ("run.buffer_size", "10"),
    ]);
    let plugins = scratch.path("plugins");
    let mut with_plugins = container_paths(&record, &scratch, Some(&plugins));
    let input = scratch.path("in.vcf");
    let output = scratch.path("out/tab.txt");
    let perl = scratch.path("perl/homo_sapiens/116_GRCh37");
    let vep_rs = engine_args(
        &record,
        "vep-rs",
        &mut with_plugins,
        &input,
        &output,
        Some(&perl),
        "tab",
    )
    .unwrap();
    assert_eq!(
        vep_rs,
        [
            "-i",
            &input,
            "-o",
            &output,
            "--offline",
            "--json_cache",
            &format!("{repo}/tests/golden/116/GRCh37-hgvs/json_cache"),
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh37",
            "--format",
            "vcf",
            "--buffer_size",
            "10",
            "--force_overwrite",
            "--quiet",
            "--stats_text",
            "--fasta",
            &format!("{repo}/tests/golden/116/GRCh37-hgvs/reference.fa.gz"),
            "--hgvs",
            "--plugin",
            "pLI",
            "--tab"
        ]
    );
    let ensembl = engine_args(
        &record,
        "ensembl-vep",
        &mut with_plugins,
        &input,
        &output,
        Some(&perl),
        "tab",
    )
    .unwrap();
    assert_eq!(
        ensembl,
        [
            "-i",
            "/work/in/in.vcf",
            "-o",
            "/work/out/tab.txt",
            "--offline",
            "--cache",
            "--dir_cache",
            "/work/cache",
            "--cache_version",
            "116",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh37",
            "--format",
            "vcf",
            "--buffer_size",
            "10",
            "--fork",
            "1",
            "--force_overwrite",
            "--quiet",
            "--stats_text",
            "--fasta",
            "/work/ref/reference.fa",
            "--hgvs",
            "--plugin",
            "pLI",
            "--dir_plugins",
            "/work/plugins",
            "--tab"
        ]
    );
    assert!(pypath::exists(&scratch.path("stage/ref/reference.fa")));
    assert!(pypath::exists(&scratch.path("stage/ref/reference.fa.fai")));
    assert!(with_plugins.mounts.contains(&(
        "/work/cache/homo_sapiens/116_GRCh37".to_string(),
        perl.clone(),
        false
    )));
    let args = |edits: &[(&str, &str)], engine: &str, cache: Option<&str>, format: &str| {
        let record = repro(edits);
        let mut paths = Paths::new(
            &record,
            &scratch.path("stage"),
            &scratch.path("out"),
            Some(&plugins),
            true,
            &[],
            &repo,
        );
        engine_args(&record, engine, &mut paths, &input, &output, cache, format).unwrap()
    };
    let has = |args: &[String], word: &str| args.iter().any(|a| a == word);
    let grch37 = args(
        &[("run.flags", "['--grch37']")],
        "vep-rs",
        Some(&perl),
        "default",
    );
    assert!(!has(&grch37, "--assembly"));
    let unnamed = args(
        &[("run.assembly_flag", "false")],
        "ensembl-vep",
        Some(&perl),
        "default",
    );
    assert!(!has(&unnamed, "--assembly"));
    let guess = args(
        &[
            ("input.text", "'21:100'"),
            ("input.format", "'guess'"),
            ("run.format_flag", "false"),
        ],
        "vep-rs",
        Some(&perl),
        "json",
    );
    assert!(!has(&guess, "--format") && !has(&guess, "--json"));
    let database = [("run.perl_mode", "'database'")];
    let db_args = args(&database, "ensembl-vep", None, "default");
    assert_eq!(db_args[4..7], ["--database", "--db_version", "116"]);
    assert!(!has(&db_args, "--offline") && !has(&db_args, "--dir_cache"));
    assert!(has(
        &args(&database, "vep-rs", None, "default"),
        "--json_cache"
    ));
    let plugin = repro(&[("run.flags", "['--plugin', 'pLI']")]);
    let mut bare = container_paths(&plugin, &scratch, None);
    assert_eq!(
        message(engine_args(
            &plugin,
            "ensembl-vep",
            &mut bare,
            &input,
            &output,
            Some(&perl),
            "default"
        )),
        "a reproduction that loads a plugin needs --plugins-dir, the Ensembl VEP plugin modules"
    );
    let own = repro(&[(
        "run.flags",
        "['--dir_plugins', '{stage}/mods', '--plugin', 'LoF']",
    )]);
    let mut bare = container_paths(&own, &scratch, None);
    let own_dir = engine_args(
        &own,
        "ensembl-vep",
        &mut bare,
        &input,
        &output,
        Some(&perl),
        "default",
    )
    .unwrap();
    assert_eq!(own_dir.iter().filter(|a| *a == "--dir_plugins").count(), 1);
    let gxf = [
        (
            "cache.json_cache",
            "'tests/golden/116/GRCh38-gtf-builder/json_cache'",
        ),
        ("cache.gxf", "'{corpus}/chr21_genes.gtf.gz'"),
    ];
    assert_eq!(
        args(&gxf, "ensembl-vep", None, "default")[4..6],
        ["--gtf", "/work/corpus/chr21_genes.gtf.gz"]
    );
    assert!(has(&args(&gxf, "vep-rs", None, "default"), "--json_cache"));
    // A file Ensembl VEP reads beside a committed cache is mounted by its directory.
    let beside = args(
        &[
            ("cache.gxf", "'tests/golden/116/GRCh38/genes.gff3.gz'"),
            ("cache.gxf_kind", "'gff'"),
        ],
        "ensembl-vep",
        None,
        "vcf",
    );
    assert_eq!(beside[4..6], ["--gff", "/work/gxf/genes.gff3.gz"]);
    assert_eq!(beside.last().map(String::as_str), Some("--vcf"));
}

#[test]
fn a_local_ensembl_vep_runs_its_sibling_tools() {
    let scratch = Scratch::new("siblings");
    let touch = "while [ $# -gt 0 ]; do if [ \"$1\" = \"-o\" ]; then : > \"$2\"; fi; shift; done";
    for name in ["vep", "filter_vep", "variant_recoder"] {
        script(&scratch, &format!("tools/{name}"), touch);
    }
    let tools = scratch.path("tools");
    let vep = format!("{tools}/vep");
    let options = |ensembl_vep: String| RunOptions {
        ensembl_vep: Some(ensembl_vep),
        vep: vep.clone(),
        timeout: 60,
        repo: repo_root(),
        ..RunOptions::default()
    };
    let tool_run = |tool: &'static str| {
        [
            ("run.tool", tool),
            ("run.args", "['-i', '{input}', '-o', '{out}']"),
            ("cache.assembly", DELETE),
            ("cache.json_cache", DELETE),
        ]
    };
    for (tool, program) in [("filter", "filter_vep"), ("recode", "variant_recoder")] {
        let quoted = format!("'{tool}'");
        let edits = [
            ("run.tool", quoted.as_str()),
            ("run.args", "['-i', '{input}', '-o', '{out}']"),
            ("cache.assembly", DELETE),
            ("cache.json_cache", DELETE),
        ];
        let record = repro(&edits);
        let stage = scratch.mkdir(&format!("stage-{tool}"));
        let out = scratch.path(&format!("out-{tool}"));
        let mut paths = Paths::new(&record, &stage, &out, None, false, &[], &repo_root());
        let input = write_input(&record, &paths).unwrap();
        let ensembl = run_engine(
            "ensembl-vep",
            &record,
            &mut paths,
            &options(vep.clone()),
            None,
            &input,
            "default",
        )
        .unwrap();
        assert_eq!(ensembl.command[0], format!("{tools}/{program}"));
        assert_eq!(ensembl.exit_code, 0);
        assert_eq!(
            std::fs::read_to_string(format!("{out}/ensembl-vep/command.txt")).unwrap(),
            format!("{tools}/{program} -i {input} -o {out}/ensembl-vep/default.txt\nexit 0\n")
        );
        let vep_rs = run_engine(
            "vep-rs",
            &record,
            &mut paths,
            &options(vep.clone()),
            None,
            &input,
            "default",
        )
        .unwrap();
        assert_eq!(vep_rs.command[..2], [vep.clone(), tool.to_string()]);
    }
    let record = repro(&tool_run("'filter'"));
    let stage = scratch.mkdir("stage-env");
    let mut paths = Paths::new(
        &record,
        &stage,
        &scratch.path("out-env"),
        None,
        false,
        &[],
        &repo_root(),
    );
    let input = write_input(&record, &paths).unwrap();
    let run = run_engine(
        "ensembl-vep",
        &record,
        &mut paths,
        &options(format!("/usr/bin/env {vep}")),
        None,
        &input,
        "default",
    )
    .unwrap();
    assert_eq!(
        run.command[..2],
        ["/usr/bin/env".to_string(), format!("{tools}/filter_vep")]
    );
    // A command that cannot be split, or that does not exist, is the entry's error.
    let refused = run_engine(
        "ensembl-vep",
        &record,
        &mut paths,
        &options("'unclosed".to_string()),
        None,
        &input,
        "default",
    );
    assert_eq!(message(refused), "No closing quotation");
    let absent = run_engine(
        "ensembl-vep",
        &record,
        &mut paths,
        &options(scratch.path("absent/vep")),
        None,
        &input,
        "default",
    );
    assert_eq!(
        message(absent),
        format!(
            "[Errno 2] No such file or directory: '{}'",
            scratch.path("absent/filter_vep")
        )
    );
}

#[test]
fn an_engine_that_outlasts_the_time_allowed_is_ended() {
    let scratch = Scratch::new("timeout");
    let sleeper = script(&scratch, "tools/vep", "sleep 30");
    let record = repro(&[]);
    let stage = scratch.mkdir("stage");
    let out = scratch.path("out");
    let mut paths = Paths::new(&record, &stage, &out, None, false, &[], &repo_root());
    let input = write_input(&record, &paths).unwrap();
    let options = RunOptions {
        ensembl_vep: Some(sleeper.clone()),
        vep: sleeper.clone(),
        timeout: 1,
        repo: repo_root(),
        ..RunOptions::default()
    };
    let started = std::time::Instant::now();
    let refused = message(run_engine(
        "vep-rs", &record, &mut paths, &options, None, &input, "default",
    ));
    assert!(started.elapsed().as_secs() < 20, "the run was not ended");
    assert!(
        refused.starts_with(&format!("Command '['{sleeper}', '-i', '{input}', '-o', ")),
        "{refused}"
    );
    assert!(
        refused.ends_with("]' timed out after 1 seconds"),
        "{refused}"
    );
    assert!(!pypath::exists(&format!("{out}/vep-rs/command.txt")));
}

#[test]
fn a_time_limit_past_what_the_clock_holds_ends_no_run() {
    let command = words(&["sh", "-c", "exit 3"]);
    let finished = run_command(&command, None, &environment(), ".", Some(i64::MAX));
    assert!(matches!(finished, Ok(Finished { code: 3, .. })));
}
