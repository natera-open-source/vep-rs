// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Refusals a cache's contents raise once a run has been configured, under the
//! `cache` component of the error contract.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn vep() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    cmd.stdin(Stdio::null());
    cmd
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// Copies the directory tree at `from` to `to`.
fn copy_tree(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.path().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), &target).unwrap();
        }
    }
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

fn assert_refused_without_gnomade(output: &Output) {
    let err = stderr(output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let error_lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(
        error_lines,
        ["ERROR: [cache] gnomAD data is not available in this cache"],
        "{err}"
    );
}

/// A cache at `dir/<name>` whose one variation file is
/// `variations/21/1-1000000_var.gz`, holding one record whose cells follow
/// `variation_cols` (the fixed columns, then `populations`), with an empty
/// transcript shard so the cache resolves.
fn variation_cache(dir: &Path, name: &str, populations: &[&str]) -> PathBuf {
    let cache = dir.join(name);
    let shard = cache.join("transcripts").join("21").join("1-1000000.json");
    std::fs::create_dir_all(shard.parent().unwrap()).unwrap();
    std::fs::write(&shard, b"[]").unwrap();
    let region_file = cache.join("variations").join("21").join("1-1000000_var.gz");
    std::fs::create_dir_all(region_file.parent().unwrap()).unwrap();
    let mut cells = vec!["rs1", "0", "0", "500", "500", "A/G", "1"];
    cells.extend(populations.iter().map(|_| "G:0.5"));
    let mut encoder = flate2::write::GzEncoder::new(
        std::fs::File::create(&region_file).unwrap(),
        flate2::Compression::default(),
    );
    writeln!(encoder, "{}", cells.join(" ")).unwrap();
    encoder.finish().unwrap();
    let mut columns = vec![
        "variation_name",
        "failed",
        "somatic",
        "start",
        "end",
        "allele_string",
        "strand",
    ];
    columns.extend_from_slice(populations);
    std::fs::write(
        cache.join("info.json"),
        format!(
            r#"{{"assembly":"GRCh38","source_versions":{{"dbSNP":"156"}},"variation_cols":{}}}"#,
            serde_json::to_string(&columns).unwrap()
        ),
    )
    .unwrap();
    cache
}

/// Runs `vep` over one variant on chromosome 21 against `cache` with `flags`.
fn run(dir: &Path, cache: &Path, flags: &[&str]) -> Output {
    let input = dir.join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t500\t.\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    run_on(dir, &input, cache, flags)
}

/// Runs `vep` on `input`, which need not exist, against `cache` with `flags`.
fn run_on(dir: &Path, input: &Path, cache: &Path, flags: &[&str]) -> Output {
    let out = dir.join("out.txt");
    let _ = std::fs::remove_file(&out);
    vep()
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .args([
            "--offline",
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .args(flags)
        .output()
        .unwrap()
}

/// `--af_gnomad` on a cache whose `variation_cols` lack `gnomADe` ends the run
/// with Ensembl VEP's message, `ERROR: gnomAD data is not available in this
/// cache` (CacheDir.pm:188-197 at release/116.2), under the `cache` component.
/// The check reads that spelling alone: `--af_gnomade` on the same cache runs,
/// `--everything` skips the check (reachable with `--safe`, since
/// `--everything` switches `--af_gnomade` on and the pair is otherwise
/// refused, Config.pm:366-370, 619, 961), and a cache carrying `gnomADe`
/// passes. `--show_cache_info` builds the annotation sources before it prints
/// (BaseRunner.pm:415), so the refusal fires there too.
#[test]
fn af_gnomad_without_gnomade_column_is_ensembl_error() {
    let dir = tempfile::tempdir().unwrap();
    let without = variation_cache(dir.path(), "without_gnomade", &["AF", "gnomADg"]);
    let refused = run(dir.path(), &without, &["--af_gnomad"]);
    assert_refused_without_gnomade(&refused);
    assert!(
        !dir.path().join("out.txt").exists(),
        "no output for a refused run"
    );

    let info_only = run(dir.path(), &without, &["--af_gnomad", "--show_cache_info"]);
    assert_refused_without_gnomade(&info_only);
    assert!(
        info_only.stdout.is_empty(),
        "the cache info is not printed before the refusal"
    );
    // The parser's arguments build the sources before it opens the file
    // (Runner.pm:827), so the refusal precedes a missing input's, as it does
    // on an annotation run, whose sources are built before the input opens
    // (Runner.pm:128, 136).
    for flags in [
        &["--af_gnomad", "--show_cache_info"][..],
        &["--af_gnomad"][..],
    ] {
        let missing_input = run_on(dir.path(), &dir.path().join("missing.vcf"), &without, flags);
        assert_refused_without_gnomade(&missing_input);
    }

    let paired = run(dir.path(), &without, &["--af_gnomad", "--everything"]);
    let err = stderr(&paired);
    assert_eq!(paired.status.code(), Some(1), "{err}");
    assert!(
        err.contains("ERROR: [cli] Can't use --af_gnomade and --af_gnomad together"),
        "{err}"
    );

    for flags in [
        &["--af_gnomade"][..],
        &["--af_gnomad", "--everything", "--safe"][..],
        &["--af_gnomadg", "--max_af"][..],
    ] {
        let output = run(dir.path(), &without, flags);
        let err = stderr(&output);
        assert_eq!(output.status.code(), Some(0), "{flags:?}: {err}");
        assert!(!err.contains("ERROR:"), "{flags:?}: {err}");
    }

    let with = variation_cache(
        dir.path(),
        "with_gnomade",
        &["AF", "gnomADe", "gnomADe_NFE"],
    );
    let info = run(dir.path(), &with, &["--af_gnomad", "--show_cache_info"]);
    assert_eq!(info.status.code(), Some(0), "{}", stderr(&info));
    assert_eq!(String::from_utf8_lossy(&info.stdout), "dbSNP\t156\n");
    let output = run(dir.path(), &with, &["--af_gnomad"]);
    let err = stderr(&output);
    assert_eq!(output.status.code(), Some(0), "{err}");
    let written = std::fs::read_to_string(dir.path().join("out.txt")).unwrap();
    let row = written.lines().find(|l| !l.starts_with('#')).unwrap();
    assert!(
        row.contains("gnomADe_AF=0.5") && row.contains("gnomADe_NFE_AF=0.5"),
        "the alias prints the gnomAD exome fields:\n{row}"
    );

    // The declared columns decide, not the files behind them (CacheDir.pm:179-197
    // reads `variation_cols`): a cache that declares columns without gnomADe and
    // holds no variation file is refused all the same.
    std::fs::remove_dir_all(without.join("variations")).unwrap();
    let declared_only = run(dir.path(), &without, &["--af_gnomad"]);
    assert_refused_without_gnomade(&declared_only);
}

/// A region file whose second line is one byte over the line cap ends the run
/// under the `cache` component as the variant's region is read, naming the
/// file and the line; nothing is written past the header. The same file with
/// the line at the cap annotates the variant with both records.
#[test]
fn a_variation_line_past_the_cap_ends_the_run_under_cache() {
    let cap = vep_io::input::MAX_LINE_BYTES;
    let dir = tempfile::tempdir().unwrap();
    let cache = variation_cache(dir.path(), "long_line", &["AF"]);
    let region_file = cache.join("variations").join("21").join("1-1000000_var.gz");
    let write_region = |second_line_len: usize| {
        let mut second = String::from("rs2 0 0 500 500 A/G 1 G:0.5 ");
        second.push_str(&"x".repeat(second_line_len - second.len()));
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(&region_file).unwrap(),
            flate2::Compression::default(),
        );
        writeln!(encoder, "rs1 0 0 500 500 A/G 1 G:0.5\n{second}").unwrap();
        encoder.finish().unwrap();
    };

    write_region(cap + 1);
    let output = run(dir.path(), &cache, &["--check_existing"]);
    let err = stderr(&output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let error_lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(
        error_lines,
        [format!(
            "ERROR: [cache] reading {}: line 2 is longer than the {cap}-byte limit for one line",
            region_file.display()
        )],
        "{err}"
    );
    let written = std::fs::read_to_string(dir.path().join("out.txt")).unwrap_or_default();
    assert!(
        written.lines().all(|l| l.starts_with('#')),
        "no row past the header:\n{written}"
    );

    write_region(cap);
    let output = run(dir.path(), &cache, &["--check_existing"]);
    assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
    let written = std::fs::read_to_string(dir.path().join("out.txt")).unwrap();
    let row = written.lines().find(|l| !l.starts_with('#')).unwrap();
    assert!(row.contains("rs1") && row.contains("rs2"), "{row}");
}

/// A transcript shard the cache cannot read ends a VCF run under the `cache`
/// component, the component of a cache whose content is unusable, with the
/// chromosome and the loader's refusal of the shard: the record is sound, so
/// neither `input` nor the run's own component applies. The corpus cache is
/// copied without its checksum manifest, so the record's chromosome load
/// reaches the shard with its gzip stream cut; every shard of the chromosome
/// is read for any record on it, so the record's position need not fall in the
/// cut one. Each pipeline mode loads the chromosome in its own coordinator and
/// both refuse alike; nothing past the header reaches the output. The line
/// ends in the inflater's own words for the cut, which are its to choose.
#[test]
fn a_transcript_shard_the_cache_cannot_read_ends_a_vcf_run_under_cache() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    copy_tree(&golden("GRCh38-input-hgvs").join("json_cache"), &cache);
    std::fs::remove_file(cache.join("CHECKSUMS.sha256")).unwrap();
    let shard = cache.join("transcripts/21/42000001-43000000.json.gz");
    let bytes = std::fs::read(&shard).unwrap();
    std::fs::write(&shard, &bytes[..bytes.len() / 2]).unwrap();
    let input = dir.path().join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30282169\t.\tC\tA\t.\t.\t.\n",
    )
    .unwrap();
    for mode in ["streaming", "batched"] {
        let output = vep()
            .arg("-i")
            .arg(&input)
            .args(["-o", "STDOUT", "--offline", "--json_cache"])
            .arg(&cache)
            .args([
                "--assembly",
                "GRCh38",
                "--no_stats",
                "--quiet",
                "--pipeline_mode",
                mode,
            ])
            .arg("--warning_file")
            .arg(dir.path().join("warnings.txt"))
            .output()
            .unwrap();
        let err = stderr(&output);
        assert_eq!(output.status.code(), Some(1), "{mode}: {err}");
        let prefix = format!(
            "ERROR: [cache] Cannot read the transcript cache for chr21: Failed to inflate {}: ",
            shard.display()
        );
        assert!(err.starts_with(&prefix), "{mode}: {err}");
        assert_eq!(err.lines().count(), 1, "{mode}: {err}");
        let stdout = String::from_utf8_lossy(&output.stdout);
        assert!(
            stdout.lines().all(|line| line.starts_with('#')),
            "{mode}: {stdout}"
        );
    }
}
