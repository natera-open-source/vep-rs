// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cache checksum manifest at run time: every cache file `vep` reads is
//! checked against `CHECKSUMS.sha256`, a changed or unlisted file ends the
//! run with one `ERROR: [cache] checksum mismatch for <relative path> in
//! <cache dir>` line, a listed file that is missing ends it as the cache is
//! resolved, and a cache without the manifest warns once and annotates as
//! the same cache with it does.

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use vep_core::cache_manifest::{CacheManifest, MANIFEST_FILE};

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// Copies the directory tree at `from` to `to`.
fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.path().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), &target).unwrap();
        }
    }
}

/// A copy of the corpus's `json_cache` under `dir/cache`.
fn cache_copy(corpus: &str, dir: &Path) -> PathBuf {
    let copy = dir.join("cache");
    copy_tree(&golden(corpus).join("json_cache"), &copy);
    copy
}

/// Runs `vep` quietly on `input` against `cache` with `flags`, writing
/// `dir/out.txt` and its warnings file beside it.
fn run(dir: &Path, input: &Path, cache: &Path, assembly: &str, flags: &[&str]) -> Output {
    let out = dir.join("out.txt");
    let _ = fs::remove_file(&out);
    let _ = fs::remove_file(warnings_file(dir));
    Command::new(env!("CARGO_BIN_EXE_vep"))
        .stdin(Stdio::null())
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .args([
            "--offline",
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            assembly,
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .args(flags)
        .output()
        .unwrap()
}

fn warnings_file(dir: &Path) -> PathBuf {
    dir.join("out.txt_warnings.txt")
}

fn output(dir: &Path) -> String {
    fs::read_to_string(dir.join("out.txt")).unwrap()
}

/// The output without the header line that names the run's wall-clock time,
/// so two runs compare whole whatever second each fell in.
fn output_without_timestamp(dir: &Path) -> String {
    output(dir)
        .lines()
        .filter(|line| !line.starts_with("## Output produced at "))
        .collect::<Vec<_>>()
        .join("\n")
}

/// The output without its `##` header lines, which name the run's paths.
fn rows(dir: &Path) -> String {
    output(dir)
        .lines()
        .filter(|line| !line.starts_with("##"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn warnings(dir: &Path) -> String {
    fs::read_to_string(warnings_file(dir)).unwrap_or_default()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The run succeeded and the warnings file carries no manifest warning.
fn assert_clean(output: &Output, dir: &Path) {
    assert_eq!(output.status.code(), Some(0), "{}", stderr(output));
    assert!(
        !warnings(dir).contains("checksum manifest"),
        "{}",
        warnings(dir)
    );
}

/// The run ended with the one error line of the contract for `relative`.
fn assert_mismatch(output: &Output, relative: &str, cache: &Path) {
    let err = stderr(output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let error_lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(
        error_lines,
        [format!(
            "ERROR: [cache] checksum mismatch for {relative} in {}",
            cache.display()
        )
        .as_str()],
        "{err}"
    );
}

/// Replaces the byte at `at` of the file at `path` with its complement.
fn flip_byte(path: &Path, at: usize) {
    let mut bytes = fs::read(path).unwrap();
    bytes[at] ^= 0xff;
    fs::write(path, bytes).unwrap();
}

/// The first shard of chromosome 21, by name.
fn first_shard(cache: &Path) -> PathBuf {
    let mut shards: Vec<PathBuf> = fs::read_dir(cache.join("transcripts/21"))
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .collect();
    shards.sort();
    shards.remove(0)
}

const COLOCATED_FLAGS: [&str; 3] = ["--check_existing", "--af", "--pubmed"];

/// Every committed corpus cache carries a manifest `sha256sum` accepts, and an
/// intact cache annotates with no manifest warning. Without the manifest the
/// same cache annotates byte-identically and the warnings file carries the
/// one line `cache <dir> carries no checksum manifest`.
#[test]
fn intact_cache_annotates_and_a_manifestless_cache_warns_once() {
    let corpora = fs::read_dir(golden("")).unwrap();
    let mut checked = 0;
    for entry in corpora {
        let json_cache = entry.unwrap().path().join("json_cache");
        if !json_cache.is_dir() {
            continue;
        }
        let manifest = CacheManifest::load(&json_cache).unwrap();
        assert!(manifest.is_present(), "{}", json_cache.display());
        checked += 1;
    }
    assert!(checked > 30, "{checked} corpus caches");

    let dir = tempfile::tempdir().unwrap();
    let cache = cache_copy("GRCh38-colocated", dir.path());
    let input = golden("GRCh38-colocated").join("variants.vcf");
    let with = run(dir.path(), &input, &cache, "GRCh38", &COLOCATED_FLAGS);
    assert_clean(&with, dir.path());
    let expected = output_without_timestamp(dir.path());
    assert!(expected.lines().any(|l| !l.starts_with('#')), "{expected}");

    fs::remove_file(cache.join(MANIFEST_FILE)).unwrap();
    let without = run(dir.path(), &input, &cache, "GRCh38", &COLOCATED_FLAGS);
    assert_eq!(without.status.code(), Some(0), "{}", stderr(&without));
    assert_eq!(output_without_timestamp(dir.path()), expected);
    let warning_lines: Vec<String> = warnings(dir.path())
        .lines()
        .filter(|line| line.contains("checksum manifest"))
        .map(str::to_string)
        .collect();
    assert_eq!(
        warning_lines,
        [format!(
            "WARNING: cache {} carries no checksum manifest",
            cache.display()
        )],
        "{}",
        warnings(dir.path())
    );
}

/// A changed byte in a transcript shard, in the variation index, in
/// `info.json`, a truncated `all_vars.gz` and a shard the manifest does not
/// list each end the run with the contract's line naming the file.
#[test]
fn changed_truncated_and_unlisted_cache_files_are_refused() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_copy("GRCh38-colocated", dir.path());
    let input = golden("GRCh38-colocated").join("variants.vcf");
    let manifest = fs::read(cache.join(MANIFEST_FILE)).unwrap();
    let pristine = dir.path().join("pristine");
    copy_tree(&cache, &pristine);
    let restore = |cache: &Path| {
        fs::remove_dir_all(cache).unwrap();
        copy_tree(&pristine, cache);
        assert_eq!(fs::read(cache.join(MANIFEST_FILE)).unwrap(), manifest);
    };

    let shard = first_shard(&cache);
    flip_byte(&shard, 40);
    let refused = run(dir.path(), &input, &cache, "GRCh38", &COLOCATED_FLAGS);
    assert_mismatch(
        &refused,
        &format!(
            "transcripts/21/{}",
            shard.file_name().unwrap().to_str().unwrap()
        ),
        &cache,
    );

    restore(&cache);
    flip_byte(&cache.join("variations/21/all_vars.gz.csi"), 20);
    let refused = run(dir.path(), &input, &cache, "GRCh38", &COLOCATED_FLAGS);
    assert_mismatch(&refused, "variations/21/all_vars.gz.csi", &cache);

    restore(&cache);
    let all_vars = cache.join("variations/21/all_vars.gz");
    let bytes = fs::read(&all_vars).unwrap();
    fs::write(&all_vars, &bytes[..bytes.len() / 2]).unwrap();
    let refused = run(dir.path(), &input, &cache, "GRCh38", &COLOCATED_FLAGS);
    assert_mismatch(&refused, "variations/21/all_vars.gz", &cache);

    restore(&cache);
    let mut info = fs::OpenOptions::new()
        .append(true)
        .open(cache.join("info.json"))
        .unwrap();
    writeln!(info).unwrap();
    drop(info);
    let refused = run(dir.path(), &input, &cache, "GRCh38", &COLOCATED_FLAGS);
    assert_mismatch(&refused, "info.json", &cache);
    let info_run = Command::new(env!("CARGO_BIN_EXE_vep"))
        .args(["--show_cache_info", "--json_cache", cache.to_str().unwrap()])
        .output()
        .unwrap();
    assert_mismatch(&info_run, "info.json", &cache);

    restore(&cache);
    let extra = cache.join("transcripts/21/47000001-48000000.json.gz");
    fs::copy(first_shard(&cache), &extra).unwrap();
    let refused = run(dir.path(), &input, &cache, "GRCh38", &COLOCATED_FLAGS);
    assert_mismatch(&refused, "transcripts/21/47000001-48000000.json.gz", &cache);

    restore(&cache);
    fs::write(cache.join(MANIFEST_FILE), "not a manifest\n").unwrap();
    let refused = run(dir.path(), &input, &cache, "GRCh38", &COLOCATED_FLAGS);
    let err = stderr(&refused);
    assert_eq!(refused.status.code(), Some(1), "{err}");
    assert!(
        err.lines().any(|line| line
            == format!(
                "ERROR: [cache] {}: line 1 is not a sha256sum line",
                cache.join(MANIFEST_FILE).display()
            )),
        "{err}"
    );
}

/// A file the manifest lists but the directory lacks, a shard or a variation
/// file of a copy that stopped short, ends the run before anything is read,
/// naming the file; an unlisted extra file is reported only when read.
#[test]
fn missing_listed_files_are_refused_at_resolution() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_copy("GRCh38-colocated", dir.path());
    let input = golden("GRCh38-colocated").join("variants.vcf");
    let shard = first_shard(&cache);
    let shard_name = shard.file_name().unwrap().to_str().unwrap().to_string();
    fs::remove_file(&shard).unwrap();
    let refused = run(dir.path(), &input, &cache, "GRCh38", &COLOCATED_FLAGS);
    let err = stderr(&refused);
    assert_eq!(refused.status.code(), Some(1), "{err}");
    assert_eq!(
        err.lines().filter(|l| l.starts_with("ERROR: [")).collect::<Vec<_>>(),
        [format!(
            "ERROR: [cache] transcripts/21/{shard_name} is listed in CHECKSUMS.sha256 but missing from {}",
            cache.display()
        )
        .as_str()],
        "{err}"
    );
    let info_run = Command::new(env!("CARGO_BIN_EXE_vep"))
        .current_dir(dir.path())
        .args(["--show_cache_info", "--json_cache", cache.to_str().unwrap()])
        .output()
        .unwrap();
    assert_eq!(info_run.status.code(), Some(1), "{}", stderr(&info_run));

    let cache = dir.path().join("variations-gone");
    copy_tree(&golden("GRCh38-colocated").join("json_cache"), &cache);
    fs::remove_dir_all(cache.join("variations")).unwrap();
    let refused = run(dir.path(), &input, &cache, "GRCh38", &[]);
    let err = stderr(&refused);
    assert_eq!(refused.status.code(), Some(1), "{err}");
    assert!(
        err.lines().any(|l| l
            == format!(
                "ERROR: [cache] variations/21/all_vars.gz is listed in CHECKSUMS.sha256 but missing from {}",
                cache.display()
            )),
        "{err}"
    );
}

/// The cache's `chr_synonyms.txt` is checked when it is read, and a
/// `--synonyms` file of the user's is not a cache file: an edited copy of the
/// same text passes through the flag and is refused in the cache.
#[test]
fn synonyms_file_of_the_cache_is_checked() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_copy("GRCh38-synonyms", dir.path());
    let input = golden("GRCh38-synonyms").join("variants.vcf");
    let synonyms = cache.join("chr_synonyms.txt");
    assert!(synonyms.is_file());
    let clean = run(dir.path(), &input, &cache, "GRCh38", &[]);
    assert_clean(&clean, dir.path());
    let expected = rows(dir.path());
    assert!(expected.contains("NC_000021.9:"), "{expected}");

    let edited = format!("{}\n", fs::read_to_string(&synonyms).unwrap());
    let mine = dir.path().join("mine.txt");
    fs::write(&mine, &edited).unwrap();
    let through_flag = run(
        dir.path(),
        &input,
        &cache,
        "GRCh38",
        &["--synonyms", mine.to_str().unwrap()],
    );
    assert_clean(&through_flag, dir.path());
    assert_eq!(rows(dir.path()), expected);

    fs::write(&synonyms, &edited).unwrap();
    let refused = run(dir.path(), &input, &cache, "GRCh38", &[]);
    assert_mismatch(&refused, "chr_synonyms.txt", &cache);
}

/// `--format id` reads each variation file start to end, so a changed
/// `all_vars.gz` is refused there by the digest of the whole file.
#[test]
fn id_input_checks_the_whole_variation_file() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_copy("GRCh38-input-id", dir.path());
    let input = golden("GRCh38-input-id").join("variants.id.txt");
    let flags = ["--format", "id"];
    let clean = run(dir.path(), &input, &cache, "GRCh38", &flags);
    assert_clean(&clean, dir.path());
    assert!(output(dir.path()).lines().any(|l| !l.starts_with('#')));

    let all_vars = cache.join("variations/21/all_vars.gz");
    let mut file = fs::OpenOptions::new().append(true).open(&all_vars).unwrap();
    // An empty gzip member appended: the stream still inflates to the same records.
    let empty = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default())
        .finish()
        .unwrap();
    file.write_all(&empty).unwrap();
    drop(file);
    let refused = run(dir.path(), &input, &cache, "GRCh38", &flags);
    assert_mismatch(&refused, "variations/21/all_vars.gz", &cache);
}

/// A FASTA found in the cache directory is read by byte range: its `.fai` and
/// `.gzi`, which are read whole, are hashed, and the FASTA must be listed. The
/// manifest the builder writes over the same files passes; a changed `.fai`
/// and a FASTA the manifest does not list are refused; the same files named
/// by `--fasta` are the user's and are read unchecked.
#[test]
fn fasta_inside_the_cache_is_checked_through_its_index() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_copy("GRCh37-hgvs", dir.path());
    let corpus = golden("GRCh37-hgvs");
    let input = corpus.join("variants.vcf");
    let outside = run(
        dir.path(),
        &input,
        &cache,
        "GRCh37",
        &[
            "--hgvs",
            "--fasta",
            corpus.join("reference.fa.gz").to_str().unwrap(),
        ],
    );
    assert_clean(&outside, dir.path());
    let expected = rows(dir.path());
    assert!(expected.contains("HGVSc="), "{expected}");

    for name in [
        "reference.fa.gz",
        "reference.fa.gz.fai",
        "reference.fa.gz.gzi",
    ] {
        fs::copy(corpus.join(name), cache.join(name)).unwrap();
    }
    let unlisted = run(dir.path(), &input, &cache, "GRCh37", &["--hgvs"]);
    assert_mismatch(&unlisted, "reference.fa.gz", &cache);

    CacheManifest::write(&cache).unwrap();
    let inside = run(dir.path(), &input, &cache, "GRCh37", &["--hgvs"]);
    assert_clean(&inside, dir.path());
    assert_eq!(rows(dir.path()), expected);

    let fai = cache.join("reference.fa.gz.fai");
    let mut file = fs::OpenOptions::new().append(true).open(&fai).unwrap();
    writeln!(file).unwrap();
    drop(file);
    let refused = run(dir.path(), &input, &cache, "GRCh37", &["--hgvs"]);
    assert_mismatch(&refused, "reference.fa.gz.fai", &cache);
    let through_flag = run(
        dir.path(),
        &input,
        &cache,
        "GRCh37",
        &[
            "--hgvs",
            "--fasta",
            cache.join("reference.fa.gz").to_str().unwrap(),
        ],
    );
    assert_clean(&through_flag, dir.path());
}

/// A plain FASTA in the cache directory is read by byte range through its
/// `.fai`, so a copy that lost its tail is refused by the length the index
/// requires, where a changed `.fai` is refused by its digest.
#[test]
fn plain_fasta_inside_the_cache_is_checked_for_its_length() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_copy("GRCh37-hgvs", dir.path());
    let corpus = golden("GRCh37-hgvs");
    let input = corpus.join("variants.vcf");
    let mut plain = Vec::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(corpus.join("reference.fa.gz")).unwrap())
        .read_to_end(&mut plain)
        .unwrap();
    fs::write(cache.join("reference.fa"), &plain).unwrap();
    fs::copy(
        corpus.join("reference.fa.fai"),
        cache.join("reference.fa.fai"),
    )
    .unwrap();
    CacheManifest::write(&cache).unwrap();
    let whole = run(dir.path(), &input, &cache, "GRCh37", &["--hgvs"]);
    assert_clean(&whole, dir.path());
    assert!(output(dir.path()).contains("HGVSc="));

    fs::write(cache.join("reference.fa"), &plain[..plain.len() / 2]).unwrap();
    let refused = run(dir.path(), &input, &cache, "GRCh37", &["--hgvs"]);
    assert_mismatch(&refused, "reference.fa", &cache);
}

/// `--show_cache_info` prints the same lines with and without the manifest.
#[test]
fn show_cache_info_is_unchanged_by_the_manifest() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_copy("GRCh38-colocated", dir.path());
    let show = |cache: &Path| {
        let output = Command::new(env!("CARGO_BIN_EXE_vep"))
            .current_dir(dir.path())
            .args(["--show_cache_info", "--json_cache", cache.to_str().unwrap()])
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        String::from_utf8(output.stdout).unwrap()
    };
    let with = show(&cache);
    assert!(with.contains("dbSNP\t"), "{with}");
    fs::remove_file(cache.join(MANIFEST_FILE)).unwrap();
    assert_eq!(show(&cache), with);
}
