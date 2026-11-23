// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs of the `vep` binary on caches whose version is not the Ensembl release
//! and whose species and assembly are not human's two builds, laid out under a
//! cache root as Ensembl VEP lays them out
//! (`<root>/<species>/<cache_version>_<assembly>`, `CacheDir.pm:303-376`): an
//! Ensembl Genomes cache, in a directory named for that division's release, and
//! a human pangenome cache, under its own species directory, version and
//! assembly name. Each cache is a copy of the synthetic species cache
//! (`fixtures/species/json_cache`, described by its README) with its metadata
//! renamed, so every run annotates the same two transcripts and its rows are
//! those of the fixture itself.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The species directory of the T2T-CHM13 pangenome cache.
const PANGENOME: &str = "homo_sapiens_gca009914755v4";

/// The assembly name that cache's directory carries.
const T2T: &str = "T2T-CHM13v2.0";

/// Five records on the fixture's one chromosome: upstream of both transcripts,
/// inside the coding transcript, intronic, in its 3' UTR and far from both.
const RECORDS: &str = "\
1\t400\tup\tG\tC\t.\t.\t.
1\t1014\tmis\tC\tT\t.\t.\t.
1\t1100\tintron\tA\tG\t.\t.\t.
1\t1250\tutr3\tG\tA\t.\t.\t.
1\t50000\tfar\tG\tA\t.\t.\t.
";

/// The rows the five records annotate to: two transcripts for each of the
/// first four, one intergenic row for the fifth.
const ROWS: usize = 9;

/// The synthetic cache of another species.
fn species_cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/species/json_cache")
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

/// A copy of the species cache at `<root>/<species>/<version>_<assembly>`
/// whose `info.json` names `species`, `assembly` and `version`; without
/// `variants` the variation file and its declared columns are dropped. The
/// checksum manifest is written over the copy, so a run reads it with no
/// manifest warning.
fn layout_cache(
    root: &Path,
    species: &str,
    version: u32,
    assembly: &str,
    variants: bool,
) -> PathBuf {
    let cache = root.join(species).join(format!("{version}_{assembly}"));
    copy_tree(&species_cache(), &cache);
    let info_path = cache.join("info.json");
    let mut info: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(&info_path).unwrap()).unwrap();
    info["species"] = species.into();
    info["assembly"] = assembly.into();
    info["cache_version"] = version.into();
    info["source_versions"]["assembly"] = assembly.into();
    if !variants {
        fs::remove_dir_all(cache.join("variations")).unwrap();
        info.as_object_mut().unwrap().remove("variation_cols");
    }
    fs::write(&info_path, info.to_string()).unwrap();
    vep_core::cache_manifest::CacheManifest::write(&cache).unwrap();
    cache
}

fn write_input(dir: &Path) {
    fs::write(
        dir.join("in.vcf"),
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{RECORDS}"),
    )
    .unwrap();
}

/// Runs the binary on `in.vcf` under `dir` with `flags`, which name the cache,
/// writing `<dir>/<label>.out` and `<dir>/<label>_warnings.txt`; nothing is
/// asserted about the exit status.
fn run(dir: &Path, label: &str, flags: &[&str]) -> Output {
    common::isolated_vep(dir)
        .arg("-i")
        .arg(dir.join("in.vcf"))
        .arg("-o")
        .arg(dir.join(format!("{label}.out")))
        .arg("--warning_file")
        .arg(dir.join(format!("{label}_warnings.txt")))
        .args(["--offline", "--no_stats", "--force_overwrite"])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

/// [`run`] asserting success; the output text and the warnings file's text
/// (empty when the run wrote none).
fn annotate(dir: &Path, label: &str, flags: &[&str]) -> (String, String) {
    let output = run(dir, label, flags);
    let err = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{flags:?}: {err}");
    let text = fs::read_to_string(dir.join(format!("{label}.out"))).unwrap();
    let warnings =
        fs::read_to_string(dir.join(format!("{label}_warnings.txt"))).unwrap_or_default();
    (text, warnings)
}

/// Asserts the run exits 1 with `message` as its one `ERROR: [cache]` line and
/// writes no output file.
fn assert_cache_refusal(dir: &Path, label: &str, flags: &[&str], message: &str) {
    let output = run(dir, label, flags);
    let err = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{flags:?}: {err}");
    let error_lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(error_lines, [format!("ERROR: [cache] {message}")], "{err}");
    assert!(
        !dir.join(format!("{label}.out")).exists(),
        "{flags:?}: no output for a refused run"
    );
}

fn data_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|line| !line.starts_with('#')).collect()
}

/// The header lines that begin with `prefix`.
fn header_lines<'a>(text: &'a str, prefix: &str) -> Vec<&'a str> {
    text.lines()
        .filter(|line| line.starts_with(prefix))
        .collect()
}

/// The `## <source> version <value>` lines the header prints from the cache's
/// metadata (`OutputFactory/BaseTab.pm:130`).
fn source_version_lines(text: &str) -> Vec<&str> {
    text.lines()
        .filter(|line| {
            line.starts_with("## ")
                && !line.starts_with("## Using ")
                && !line.contains(" : ")
                && line.contains(" version ")
        })
        .collect()
}

/// The rows of the five records on the fixture cache itself, named directly.
fn fixture_rows(dir: &Path) -> Vec<String> {
    let cache = species_cache();
    let (text, _) = annotate(
        dir,
        "fixture",
        &[
            "--json_cache",
            cache.to_str().unwrap(),
            "--species",
            "synthetic_species",
            "--assembly",
            "SYN1",
        ],
    );
    let rows: Vec<String> = data_lines(&text).iter().map(|r| r.to_string()).collect();
    assert_eq!(rows.len(), ROWS, "{text}");
    rows
}

/// An Ensembl Genomes cache annotates under a `--cache_version` naming its own
/// release, and the header names the cache where Ensembl VEP's does: the
/// `## Using cache in` line carries the directory, `63_<assembly>`
/// (`OutputFactory/BaseTab.pm:114`), and of its version lines the release
/// number stands in the API line alone (`:116`), which states the software's
/// version in both engines. Without the flag the run stops for want of a cache
/// of the release (`CacheDir.pm:323, 336`), and `--json_cache` on the directory
/// itself needs no version.
#[test]
fn a_cache_of_another_release_annotates_under_its_own_version() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let root = dir.path().join("root");
    let leaf = layout_cache(&root, "synthetic_species", 63, "SYN1", true);
    let root_s = root.to_str().unwrap();
    let species = ["--dir_cache", root_s, "--species", "synthetic_species"];

    let mut versioned = species.to_vec();
    versioned.extend_from_slice(&["--cache_version", "63"]);
    let (text, warnings) = annotate(dir.path(), "versioned", &versioned);
    assert_eq!(warnings, "", "no warning");
    assert_eq!(
        data_lines(&text),
        fixture_rows(dir.path()),
        "the rows are the fixture's"
    );
    assert_eq!(
        header_lines(&text, "## Using cache in "),
        [format!("## Using cache in {}", leaf.display())],
        "{text}"
    );
    assert_eq!(
        header_lines(&text, "## Using API version "),
        [format!(
            "## Using API version {}, DB version ?",
            vep_core::VEP_VERSION
        )],
        "{text}"
    );
    assert_eq!(
        source_version_lines(&text),
        [
            "## assembly version SYN1",
            "## genebuild version synthetic_1"
        ],
        "the cache's own source versions, none of them the release:\n{text}"
    );

    assert_cache_refusal(
        dir.path(),
        "default_version",
        &species,
        &format!(
            "No cache found for synthetic_species, version {}",
            vep_core::VEP_VERSION
        ),
    );

    let (direct, warnings) = annotate(
        dir.path(),
        "direct",
        &[
            "--json_cache",
            leaf.to_str().unwrap(),
            "--species",
            "synthetic_species",
        ],
    );
    assert_eq!(warnings, "", "no warning");
    assert_eq!(data_lines(&direct), data_lines(&text));
}

/// A human pangenome cache annotates under its species directory, its own
/// version and its assembly name, given or read from the directory name
/// (`CacheDir.pm:344-356`): the rows are those the same transcripts get under
/// any other species name, and the header names the cache directory and the
/// assembly. Another `--assembly` is refused (`:347-353`), and `--refseq`,
/// which no pangenome cache exists for, stops on the species directory the
/// root lacks (`:312-318`).
#[test]
fn a_pangenome_cache_annotates_under_its_species_directory_and_assembly() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let root = dir.path().join("root");
    let leaf = layout_cache(&root, PANGENOME, 107, T2T, false);
    let root_s = root.to_str().unwrap();
    let by_name = [
        "--dir_cache",
        root_s,
        "--species",
        PANGENOME,
        "--cache_version",
        "107",
    ];
    let rows = fixture_rows(dir.path());

    let (text, warnings) = annotate(dir.path(), "by_name", &by_name);
    assert_eq!(warnings, "", "no warning");
    assert_eq!(data_lines(&text), rows, "the rows are the fixture's");
    assert_eq!(
        header_lines(&text, "## Using cache in "),
        [format!("## Using cache in {}", leaf.display())],
        "{text}"
    );
    assert_eq!(
        source_version_lines(&text),
        [
            "## assembly version T2T-CHM13v2.0",
            "## genebuild version synthetic_1"
        ],
        "{text}"
    );

    let mut with_assembly = by_name.to_vec();
    with_assembly.extend_from_slice(&["--assembly", T2T]);
    let (text, warnings) = annotate(dir.path(), "with_assembly", &with_assembly);
    assert_eq!(warnings, "", "no warning");
    assert_eq!(data_lines(&text), rows);

    let mut other_assembly = by_name.to_vec();
    other_assembly.extend_from_slice(&["--assembly", "GRCh38"]);
    assert_cache_refusal(
        dir.path(),
        "other_assembly",
        &other_assembly,
        "Cache assembly version (T2T-CHM13v2.0) and database or selected assembly version (GRCh38) do not match",
    );

    let mut refseq = by_name.to_vec();
    refseq.push("--refseq");
    assert_cache_refusal(
        dir.path(),
        "refseq",
        &refseq,
        &format!(
            "Cache directory {} not found",
            root.join(format!("{PANGENOME}_refseq")).display()
        ),
    );
}

/// What the checks written for human do with a pangenome species name. The
/// species defaults to `homo_sapiens` (`Config.pm:290`), so a run naming no
/// species looks for that directory. The spellings refused as a species are
/// human's two alone (`CacheDir.pm:310-311`): `homo_sapiens_refseq` is refused
/// on this root as on any other, with advice naming `homo_sapiens`, while the
/// pangenome name with a set suffix is looked up as a directory. LoFTEE reads
/// neither species nor assembly to name its default GERP file (`LoF.pm` line
/// 130), so the file it looks for under `loftee_path` is the GRCh38 one and the
/// plugin is not loaded without it; a run that needs no GERP file loads the
/// plugin with the warning that the assembly gives it no hint, which names
/// the cache where the name was read from it and `--assembly` where the flag
/// gave it.
#[test]
fn the_checks_written_for_human_on_a_pangenome_species_name() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let root = dir.path().join("root");
    layout_cache(&root, PANGENOME, 107, T2T, false);
    let root_s = root.to_str().unwrap();
    let versioned = ["--dir_cache", root_s, "--cache_version", "107"];
    let with = |more: &[&'static str]| -> Vec<&str> {
        let mut flags = versioned.to_vec();
        flags.extend_from_slice(more);
        flags
    };

    assert_cache_refusal(
        dir.path(),
        "default_species",
        &versioned,
        &format!(
            "Cache directory {} not found",
            root.join("homo_sapiens").display()
        ),
    );
    assert_cache_refusal(
        dir.path(),
        "human_suffix",
        &with(&["--species", "homo_sapiens_refseq"]),
        "Should not use homo_sapiens_refseq as --species. Try using flags --refseq or --merged with --species homo_sapiens",
    );
    assert_cache_refusal(
        dir.path(),
        "pangenome_suffix",
        &with(&["--species", "homo_sapiens_gca009914755v4_merged"]),
        &format!(
            "Cache directory {} not found",
            root.join(format!("{PANGENOME}_merged")).display()
        ),
    );

    let loftee = dir.path().join("loftee");
    fs::create_dir_all(&loftee).unwrap();
    let derived = format!("LoF,loftee_path:{}/", loftee.display());
    let mut flags = with(&["--species", PANGENOME, "--plugin"]);
    flags.push(&derived);
    let (text, warnings) = annotate(dir.path(), "loftee_default", &flags);
    let not_loaded = format!(
        "WARNING: Failed to instantiate plugin LoF: LoFTEE: GERP data file not found: \
         '{}/gerp_conservation_scores.homo_sapiens.GRCh38.bw'",
        loftee.display()
    );
    let lines: Vec<&str> = warnings.lines().collect();
    assert_eq!(lines.len(), 1, "{warnings}");
    assert!(lines[0].starts_with(&not_loaded), "{warnings}");
    assert!(header_lines(&text, "## LoF").is_empty(), "{text}");
    assert_eq!(data_lines(&text).len(), ROWS, "{text}");

    let (text, warnings) = annotate(
        dir.path(),
        "loftee_no_gerp",
        &with(&[
            "--species",
            PANGENOME,
            "--plugin",
            "LoF,use_gerp_end_trunc:0",
        ]),
    );
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        ["WARNING: the cache's assembly 'T2T-CHM13v2.0' is not an assembly LoF recognises; it runs without an assembly hint"],
    );
    assert!(
        !header_lines(&text, "## LoF : ").is_empty(),
        "the plugin is loaded:\n{text}"
    );
    assert_eq!(data_lines(&text).len(), ROWS, "{text}");

    let (text, warnings) = annotate(
        dir.path(),
        "loftee_no_gerp_named_assembly",
        &with(&[
            "--species",
            PANGENOME,
            "--assembly",
            T2T,
            "--plugin",
            "LoF,use_gerp_end_trunc:0",
        ]),
    );
    assert_eq!(
        warnings.lines().collect::<Vec<_>>(),
        ["WARNING: --assembly 'T2T-CHM13v2.0' is not an assembly LoF recognises; it runs without an assembly hint"],
    );
    assert!(
        !header_lines(&text, "## LoF : ").is_empty(),
        "the plugin is loaded:\n{text}"
    );
    assert_eq!(data_lines(&text).len(), ROWS, "{text}");
}

/// The assembly names that give the built-in plugins their hint, as the help
/// of `--assembly` lists them: `GRCh37`, `hg19`, `37`, `GRCh38`, `hg38` and
/// `38`, in any letter case (`Assembly::parse` of `vep-core`, whose own tests
/// pin which of the two assemblies each name gives). A cache whose own name is
/// one of them loads a built-in plugin with no warning, and so does a cache
/// that names no assembly under a flag that gives one. A name outside the
/// list, a human name with a patch suffix or another species' name that
/// carries a human assembly's digits, draws the one warning that the plugin
/// runs without a hint, naming the cache or the flag the name was read from.
#[test]
fn the_assembly_names_the_help_lists_give_the_plugins_their_hint() {
    const HINTED: [&str; 10] = [
        "GRCh37", "hg19", "37", "GRCh38", "hg38", "38", "grch37", "HG19", "GRCH38", "Hg38",
    ];
    const UNHINTED: [&str; 2] = ["GRCh38.p14", "GRCm38"];
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    // The warnings of one run that loads a built-in plugin on the cache at
    // `leaf`, with `assembly` as the flag's value.
    let plugin_run = |label: &str, leaf: &Path, assembly: Option<&str>| -> Vec<String> {
        let mut flags = vec![
            "--json_cache",
            leaf.to_str().unwrap(),
            "--species",
            "synthetic_species",
            "--plugin",
            "TSSDistance",
        ];
        if let Some(name) = assembly {
            flags.extend(["--assembly", name]);
        }
        let (text, warnings) = annotate(dir.path(), label, &flags);
        assert!(
            !header_lines(&text, "## TSSDistance : ").is_empty(),
            "the plugin is loaded:\n{text}"
        );
        assert_eq!(data_lines(&text).len(), ROWS, "{text}");
        warnings.lines().map(str::to_string).collect()
    };
    // One root for each name, so that two names differing in letter case
    // alone never share a directory.
    let named = |n: usize, name: &str| -> PathBuf {
        let root = dir.path().join(format!("named_{n}"));
        layout_cache(&root, "synthetic_species", 116, name, false)
    };
    let unnamed = layout_cache(
        &dir.path().join("unnamed"),
        "synthetic_species",
        116,
        "SYN1",
        false,
    );
    let info_path = unnamed.join("info.json");
    let mut info: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(&info_path).unwrap()).unwrap();
    info.as_object_mut().unwrap().remove("assembly");
    fs::write(&info_path, info.to_string()).unwrap();
    vep_core::cache_manifest::CacheManifest::write(&unnamed).unwrap();

    for (n, name) in HINTED.iter().enumerate() {
        let leaf = named(n, name);
        let none: [&str; 0] = [];
        assert_eq!(
            plugin_run(&format!("cache_{n}"), &leaf, None),
            none,
            "the cache's own name {name}"
        );
        assert_eq!(
            plugin_run(&format!("flag_{n}"), &unnamed, Some(name)),
            none,
            "--assembly {name} on a cache that names none"
        );
    }
    for (n, name) in UNHINTED.iter().enumerate() {
        let leaf = named(HINTED.len() + n, name);
        assert_eq!(
            plugin_run(&format!("cache_unhinted_{n}"), &leaf, None),
            [format!(
                "WARNING: the cache's assembly '{name}' is not an assembly TSSDistance \
                 recognises; it runs without an assembly hint"
            )]
        );
        assert_eq!(
            plugin_run(&format!("flag_unhinted_{n}"), &unnamed, Some(name)),
            [format!(
                "WARNING: --assembly '{name}' is not an assembly TSSDistance recognises; it \
                 runs without an assembly hint"
            )]
        );
    }
}
