// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Resolution of the annotation cache directory before the input is opened.
//!
//! `--json_cache` (alias `--full_cache_dir`) names the cache directory itself;
//! otherwise the directory is found in Ensembl VEP's layout under the cache root,
//! `<dir_cache>/<species>[_refseq|_merged]/<cache_version>_<assembly>`, as
//! `CacheDir.pm` resolves it (`dir`, release/116.2 lines 298-379). In both cases
//! the directory must hold `transcripts/` and `info.json`, `--assembly` must
//! agree with the cache (`CacheDir.pm:347` against the directory name, `:404-408`
//! against the cache's metadata, which the directory name's own assembly must
//! also match), and an indexed `.fa` or `.fa.gz` inside the directory is picked
//! up unless `--fasta` or `--no_fasta` is given (`CacheDir.pm:235-251`; Ensembl
//! VEP builds a missing index itself, while the FASTA reader here needs the
//! `.fai`, and the `.gzi` of a compressed file, beside it). The chromosome
//! synonyms are read from `--synonyms`, else from the cache's own
//! `chr_synonyms.txt` (`CacheDir.pm:254`, `Runner.pm:131`), except under
//! `--show_cache_info`, which prints before Ensembl VEP reads them. Every
//! refusal is a [`CliError`]: the cache's under the `cache` component, an
//! unreadable `--synonyms` file under `config`.
//!
//! The cache's checksum manifest (`CHECKSUMS.sha256`, [`CacheManifest`]) is
//! read here, once, before `info.json`, and every cache file the run then reads
//! whole is checked against it; a cache without one warns once and loads.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use tracing::warn;
use vep_core::cache_info::CacheInfo;
use vep_core::cache_manifest::{CacheManifest, ManifestError};
use vep_core::synonyms::ChromosomeSynonyms;

use crate::config::{CliError, Config};

/// The cache a run reads, resolved and validated.
#[derive(Debug, Clone)]
pub struct ResolvedCache {
    /// The cache directory, as the output headers name it.
    pub dir: String,
    /// The cache's checksum manifest, which every cache file the run reads is
    /// checked against; one that is not present checks nothing.
    pub manifest: Arc<CacheManifest>,
    /// The cache's `info.json`.
    pub info: CacheInfo,
    /// The assembly of the run: the cache's own, else the one the directory name
    /// carries, else `--assembly` read as a canonical name.
    pub assembly: Option<String>,
    /// An indexed `.fa` or `.fa.gz` found in the cache directory, when neither
    /// `--fasta` nor `--no_fasta` was given.
    pub fasta: Option<String>,
    /// The chromosome synonyms of the run: the `--synonyms` file, else the
    /// cache's `chr_synonyms.txt`, else an empty table. A user file replaces the
    /// cache's; the two are never merged.
    pub synonyms: Arc<ChromosomeSynonyms>,
}

/// Perl's `species` parameter spellings that Ensembl VEP refuses as a species
/// directory (`CacheDir.pm:310-311`): the suffix comes from `--refseq` or `--merged`.
const SUFFIXED_SPECIES: [&str; 2] = ["homo_sapiens_refseq", "homo_sapiens_merged"];

fn cache_error(message: impl Into<String>) -> CliError {
    CliError::new("cache", message)
}

/// `hg19` and `hg38` are accepted for `--assembly` as the canonical names they
/// stand for; every other value is compared as written, as Ensembl VEP compares it.
pub(crate) fn canonical_assembly(value: &str) -> &str {
    match value {
        "hg19" => "GRCh37",
        "hg38" => "GRCh38",
        other => other,
    }
}

/// Resolves and validates the cache the configuration names.
pub fn resolve(config: &Config) -> Result<ResolvedCache, CliError> {
    let config_assembly = config.assembly.as_deref().map(canonical_assembly);
    let (dir, dir_assembly) = match &config.json_cache {
        Some(dir) => {
            // CacheDir.pm:109.
            if !Path::new(dir).is_dir() {
                return Err(cache_error(format!("Directory {dir} does not exist")));
            }
            (dir.clone(), None)
        }
        None => layout_dir(config, config_assembly)?,
    };
    validate_leaf(&dir)?;
    let manifest = Arc::new(load_manifest(&dir)?);
    let info = read_info(&dir, &manifest)?;
    // CacheDir.pm:407-408 compares the configuration's assembly after `dir` set it
    // from the directory name (:356) with the metadata's whenever it is present, an
    // empty value included (:466 keeps the empty field, :475 stores it); Ensembl
    // VEP's text names the Storable cache's info.txt, the file a JSON cache carries
    // as info.json.
    let wanted = config_assembly.or(dir_assembly.as_deref());
    if let (Some(wanted), Some(found)) = (wanted, info.assembly.as_deref()) {
        if wanted != found {
            return Err(cache_error(format!(
                "Mismatch in assembly versions from config ({wanted}) and cache info.json file ({found})"
            )));
        }
    }
    let assembly = info
        .assembly
        .clone()
        .filter(|a| !a.is_empty())
        .or(dir_assembly)
        .or_else(|| config_assembly.map(str::to_string));
    let fasta = if config.fasta.is_some() || config.no_fasta {
        None
    } else {
        detect_fasta(&dir)
    };
    // Ensembl VEP reads the synonyms in `Runner::init` (Runner.pm:131), which a
    // `--show_cache_info` run exits before (vep:39-43), so an unreadable file
    // never stops the cache information from printing.
    let synonyms = if config.show_cache_info {
        ChromosomeSynonyms::new()
    } else {
        load_synonyms(config.synonyms.as_deref(), &dir, &manifest)?
    };
    Ok(ResolvedCache {
        dir,
        manifest,
        info,
        assembly,
        fasta,
        synonyms: Arc::new(synonyms),
    })
}

/// The run's chromosome synonyms: the `--synonyms` file when given, else the
/// cache's `chr_synonyms.txt` when it exists (`CacheDir.pm:254`), checked
/// against the manifest, else an empty table. A file that cannot be read is
/// refused with the reason the operating system gives, a refusal of this
/// engine's own: Ensembl VEP's two-argument `open` (`BaseVEP.pm:644`) fails on a
/// missing file but opens a directory and loads an empty table from it.
fn load_synonyms(
    user_file: Option<&str>,
    cache_dir: &str,
    manifest: &CacheManifest,
) -> Result<ChromosomeSynonyms, CliError> {
    let cache_file = Path::new(cache_dir).join("chr_synonyms.txt");
    let (path, component) = match user_file {
        Some(file) => (PathBuf::from(file), "config"),
        None if cache_file.exists() => (cache_file, "cache"),
        None => return Ok(ChromosomeSynonyms::new()),
    };
    let text = std::fs::read_to_string(&path).map_err(|e| {
        CliError::new(
            component,
            format!(
                "Could not read synonyms file {}: {}",
                path.display(),
                os_error_text(&e)
            ),
        )
    })?;
    if user_file.is_none() {
        manifest
            .verify_bytes(&path, text.as_bytes())
            .map_err(manifest_error)?;
    }
    Ok(ChromosomeSynonyms::parse(&text))
}

/// The cache's checksum manifest, read once for the run, every file it lists
/// present; a cache without one warns once and is read unchecked.
fn load_manifest(dir: &str) -> Result<CacheManifest, CliError> {
    let manifest = CacheManifest::load(dir).map_err(manifest_error)?;
    if !manifest.is_present() {
        warn!("cache {dir} carries no checksum manifest");
    }
    manifest.require_complete().map_err(manifest_error)?;
    Ok(manifest)
}

/// A manifest failure as the run reports it: under the `cache` component.
pub(crate) fn manifest_error(error: ManifestError) -> CliError {
    cache_error(error.to_string())
}

/// An I/O error's text as the C library words it (`No such file or directory`),
/// the form Perl's `$!` prints, without the `(os error N)` suffix Rust appends.
fn os_error_text(error: &std::io::Error) -> String {
    let text = error.to_string();
    match text.rfind(" (os error ") {
        Some(at) if text.ends_with(')') => text[..at].to_string(),
        _ => text,
    }
}

/// The cache directory under the root in Ensembl VEP's layout, with the assembly
/// the directory name carries when a single entry matches the cache version
/// (`CacheDir.pm:303-376`).
fn layout_dir(
    config: &Config,
    config_assembly: Option<&str>,
) -> Result<(String, Option<String>), CliError> {
    let species = config.species.as_str();
    if SUFFIXED_SPECIES.contains(&species) {
        return Err(cache_error(format!(
            "Should not use {species} as --species. Try using flags --refseq or --merged with --species homo_sapiens"
        )));
    }
    let mut species_dir = species.to_string();
    if config.refseq {
        species_dir.push_str("_refseq");
    }
    if config.merged {
        species_dir.push_str("_merged");
    }
    let dir = collapse_slashes(&format!("{}/{species_dir}", config.dir_cache));
    if !Path::new(&dir).exists() {
        return Err(cache_error(format!("Cache directory {dir} not found")));
    }
    let version = config.cache_version.to_string();
    let mut matched: Vec<String> = std::fs::read_dir(&dir)
        .map_err(|e| cache_error(format!("{dir}: {e}")))?
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .filter(|name| !name.starts_with('.') && name.starts_with(&version))
        .collect();
    matched.sort();
    let (leaf, dir_assembly) = match matched.as_slice() {
        [] => {
            return Err(cache_error(format!(
                "No cache found for {species_dir}, version {version}"
            )))
        }
        [only] => {
            let dir_assembly = assembly_in_name(only);
            if let (Some(wanted), Some(found)) = (config_assembly, dir_assembly) {
                if wanted != found {
                    return Err(cache_error(format!(
                        "Cache assembly version ({found}) and database or selected assembly version ({wanted}) do not match"
                    )));
                }
            }
            (only.clone(), dir_assembly.map(str::to_string))
        }
        several => match config_assembly {
            None => {
                let prefix = format!("{version}_");
                let possibles: Vec<&str> = several
                    .iter()
                    .map(|name| name.strip_prefix(&prefix).unwrap_or(name))
                    .collect();
                return Err(cache_error(format!(
                    "Multiple assemblies found for cache version {version} ({}) - specify one using --assembly [assembly]",
                    possibles.join(", ")
                )));
            }
            Some(wanted) => {
                let leaf = format!("{version}_{wanted}");
                if !Path::new(&dir).join(&leaf).exists() {
                    return Err(cache_error(format!(
                        "No cache found for {species}, version {version}, assembly {wanted}"
                    )));
                }
                (leaf, None)
            }
        },
    };
    Ok((collapse_slashes(&format!("{dir}/{leaf}")), dir_assembly))
}

/// The assembly a cache directory name carries: what follows the first `_` that a
/// digit precedes (`CacheDir.pm:344`, `/\d+\_(.+)/`), or nothing.
fn assembly_in_name(name: &str) -> Option<&str> {
    let bytes = name.as_bytes();
    let underscore =
        (1..bytes.len()).find(|&i| bytes[i] == b'_' && bytes[i - 1].is_ascii_digit())?;
    let assembly = &name[underscore + 1..];
    (!assembly.is_empty()).then_some(assembly)
}

/// Consecutive `/` collapsed to one (`CacheDir.pm:373`).
fn collapse_slashes(path: &str) -> String {
    let mut out = String::with_capacity(path.len());
    for c in path.chars() {
        if c == '/' && out.ends_with('/') {
            continue;
        }
        out.push(c);
    }
    out
}

/// A cache directory holds `transcripts/` and `info.json`.
fn validate_leaf(dir: &str) -> Result<(), CliError> {
    let path = Path::new(dir);
    if !path.join("transcripts").is_dir() {
        return Err(cache_error(format!(
            "no cache found at {dir}: no transcripts directory"
        )));
    }
    if !path.join("info.json").is_file() {
        return Err(cache_error(format!(
            "no cache found at {dir}: no info.json"
        )));
    }
    Ok(())
}

/// The cache's `info.json`, checked against the manifest, a JSON object (a list
/// would deserialize into a metadata-less default).
fn read_info(dir: &str, manifest: &CacheManifest) -> Result<CacheInfo, CliError> {
    let path = Path::new(dir).join("info.json");
    let file_error =
        |cause: &dyn std::fmt::Display| cache_error(format!("{}: {cause}", path.display()));
    let text = std::fs::read_to_string(&path).map_err(|e| file_error(&e))?;
    manifest
        .verify_bytes(&path, text.as_bytes())
        .map_err(manifest_error)?;
    let value: serde_json::Value = serde_json::from_str(&text).map_err(|e| file_error(&e))?;
    if !value.is_object() {
        return Err(file_error(&"not a JSON object"));
    }
    serde_json::from_value(value).map_err(|e| file_error(&e))
}

/// The first `.fa` or `.fa.gz` of the cache directory, by name, that has its
/// indexes beside it (`CacheDir.pm:241`).
fn detect_fasta(dir: &str) -> Option<String> {
    vep_fasta::first_indexed_fasta(Path::new(dir)).map(|path| path.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;

    const INFO: &str = r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116,
        "source_versions":{"assembly":"GRCh38.p14","dbSNP":"156"}}"#;

    /// A cache leaf with `transcripts/` and the given `info.json` text.
    fn make_leaf(path: &Path, info: Option<&str>) {
        fs::create_dir_all(path.join("transcripts").join("21")).unwrap();
        if let Some(info) = info {
            fs::write(path.join("info.json"), info).unwrap();
        }
    }

    /// A cache root holding `homo_sapiens/<leaf>` for each leaf name.
    fn make_root(leaves: &[&str]) -> tempfile::TempDir {
        let root = tempfile::tempdir().unwrap();
        for leaf in leaves {
            let info = INFO.replace("GRCh38", assembly_in_name(leaf).unwrap_or("GRCh38"));
            make_leaf(&root.path().join("homo_sapiens").join(leaf), Some(&info));
        }
        root
    }

    fn config(flags: &[&str]) -> Config {
        let mut argv = vec!["vep"];
        argv.extend_from_slice(flags);
        Config::from_argv(argv).unwrap()
    }

    fn resolve_flags(flags: &[&str]) -> Result<ResolvedCache, CliError> {
        resolve(&config(flags))
    }

    fn error(flags: &[&str]) -> String {
        resolve_flags(flags).unwrap_err().to_string()
    }

    fn joined(root: &Path, rest: &str) -> String {
        format!("{}/{rest}", root.display())
    }

    /// `--dir_cache <root>` resolves `<root>/homo_sapiens/116_GRCh38`
    /// (`CacheDir.pm:315, 341`), the directory's assembly is read from its name
    /// (`:344-356`), consecutive `/` collapse (`:373`), and `--species`
    /// names the directory (`:309, 315`).
    #[test]
    fn perl_layout_resolves_leaf() {
        let root = make_root(&["116_GRCh38"]);
        let leaf = joined(root.path(), "homo_sapiens/116_GRCh38");
        let resolved = resolve_flags(&["--dir_cache", root.path().to_str().unwrap()]).unwrap();
        assert_eq!(resolved.dir, leaf);
        assert_eq!(resolved.assembly.as_deref(), Some("GRCh38"));
        assert_eq!(resolved.info.source_versions["dbSNP"], "156");
        assert!(resolved.fasta.is_none());
        for spelling in [
            format!("{}/", root.path().display()),
            format!("{}//", root.path().display()),
        ] {
            assert_eq!(
                resolve_flags(&["--dir_cache", &spelling]).unwrap().dir,
                leaf,
                "{spelling}"
            );
        }
        let explicit = resolve_flags(&[
            "--dir_cache",
            root.path().to_str().unwrap(),
            "--assembly",
            "GRCh38",
        ])
        .unwrap();
        assert_eq!(explicit.dir, leaf);
    }

    /// Without `--dir_cache` the root is `--dir` (`AnnotationSourceAdaptor.pm:127`).
    #[test]
    fn dir_is_the_default_cache_root() {
        let root = make_root(&["116_GRCh38"]);
        let resolved = resolve_flags(&["--dir", root.path().to_str().unwrap()]).unwrap();
        assert_eq!(resolved.dir, joined(root.path(), "homo_sapiens/116_GRCh38"));
    }

    /// `--dir_cache` wins over `--dir` (`AnnotationSourceAdaptor.pm:127`), and a
    /// directory value of `--cache` stands for `--dir_cache` when that flag is
    /// absent (`Config.pm:717`), never over it.
    #[test]
    fn dir_cache_overrides_dir() {
        let root = make_root(&["116_GRCh38"]);
        let other = tempfile::tempdir().unwrap();
        let leaf = joined(root.path(), "homo_sapiens/116_GRCh38");
        let root_s = root.path().to_str().unwrap();
        let other_s = other.path().to_str().unwrap();
        assert_eq!(
            resolve_flags(&["--dir", other_s, "--dir_cache", root_s])
                .unwrap()
                .dir,
            leaf
        );
        assert_eq!(resolve_flags(&["--cache", root_s]).unwrap().dir, leaf);
        assert_eq!(
            resolve_flags(&["--dir", other_s, "--cache", root_s])
                .unwrap()
                .dir,
            leaf
        );
        assert_eq!(
            resolve_flags(&["--cache", other_s, "--dir_cache", root_s])
                .unwrap()
                .dir,
            leaf
        );
        // A --cache value that is not a directory selects nothing (Config.pm:717 `-d`).
        assert_eq!(
            resolve_flags(&["--cache", "no-such-directory", "--dir", root_s])
                .unwrap()
                .dir,
            leaf
        );
    }

    /// `root_dir => dir_cache || dir` (`AnnotationSourceAdaptor.pm:127`): an
    /// empty or `0` `--dir_cache` is Perl-false, so the root is `--dir`, and a
    /// directory given to `--cache` is assigned over it (`Config.pm:717`, `||=`).
    #[test]
    fn empty_dir_cache_falls_back_to_dir() {
        let root = make_root(&["116_GRCh38"]);
        let root_s = root.path().to_str().unwrap();
        let leaf = joined(root.path(), "homo_sapiens/116_GRCh38");
        for value in ["", "0"] {
            assert_eq!(
                resolve_flags(&["--dir", root_s, "--dir_cache", value])
                    .unwrap()
                    .dir,
                leaf,
                "{value:?}"
            );
            assert_eq!(
                resolve_flags(&["--cache", root_s, "--dir_cache", value])
                    .unwrap()
                    .dir,
                leaf,
                "{value:?}"
            );
        }
    }

    /// `cache_version || db_version || software_version` (`CacheDir.pm:323`):
    /// `--cache_version 0` is Perl-false and resolves the release's directory,
    /// as no flag would.
    #[test]
    fn cache_version_zero_resolves_the_release() {
        let root = make_root(&["115_GRCh38", "116_GRCh38"]);
        let root_s = root.path().to_str().unwrap();
        let leaf = joined(root.path(), "homo_sapiens/116_GRCh38");
        assert_eq!(
            resolve_flags(&["--dir_cache", root_s, "--cache_version", "0"])
                .unwrap()
                .dir,
            leaf
        );
        assert_eq!(resolve_flags(&["--dir_cache", root_s]).unwrap().dir, leaf);
    }

    /// A root without the species directory is `Cache directory <dir> not found`
    /// (`CacheDir.pm:318`); a species directory without the version is
    /// `No cache found for <species>, version <v>` (`:336`).
    #[test]
    fn missing_species_dir_is_fatal() {
        let root = tempfile::tempdir().unwrap();
        let root_s = root.path().to_str().unwrap();
        assert_eq!(
            error(&["--dir_cache", root_s]),
            format!(
                "ERROR: [cache] Cache directory {} not found",
                joined(root.path(), "homo_sapiens")
            )
        );
        fs::create_dir_all(root.path().join("homo_sapiens")).unwrap();
        assert_eq!(
            error(&["--dir_cache", root_s]),
            "ERROR: [cache] No cache found for homo_sapiens, version 116"
        );
        assert_eq!(
            error(&["--dir_cache", root_s, "--merged"]),
            format!(
                "ERROR: [cache] Cache directory {} not found",
                joined(root.path(), "homo_sapiens_merged")
            )
        );
    }

    /// `--cache_version` is the directory prefix (`CacheDir.pm:323, 329`); two
    /// matching directories need `--assembly` (`:361-363`), which then names the
    /// leaf (`:368`) or fails with `:369`'s text; dotfiles are skipped (`:326`).
    #[test]
    fn cache_version_selects_directory() {
        let root = make_root(&["115_GRCh38", "116_GRCh37", "116_GRCh38"]);
        fs::write(root.path().join("homo_sapiens/.116_hidden"), "").unwrap();
        let root_s = root.path().to_str().unwrap();
        let resolved = resolve_flags(&["--dir_cache", root_s, "--cache_version", "115"]).unwrap();
        assert_eq!(resolved.dir, joined(root.path(), "homo_sapiens/115_GRCh38"));
        assert_eq!(resolved.assembly.as_deref(), Some("GRCh38"));
        assert_eq!(
            error(&["--dir_cache", root_s]),
            "ERROR: [cache] Multiple assemblies found for cache version 116 (GRCh37, GRCh38) - specify one using --assembly [assembly]"
        );
        assert_eq!(
            resolve_flags(&["--dir_cache", root_s, "--assembly", "GRCh37"])
                .unwrap()
                .dir,
            joined(root.path(), "homo_sapiens/116_GRCh37")
        );
        assert_eq!(
            error(&["--dir_cache", root_s, "--assembly", "NCBI36"]),
            "ERROR: [cache] No cache found for homo_sapiens, version 116, assembly NCBI36"
        );
        assert_eq!(
            error(&["--dir_cache", root_s, "--cache_version", "114"]),
            "ERROR: [cache] No cache found for homo_sapiens, version 114"
        );
    }

    /// `--species` names the directory under the root, suffixed `_refseq` or
    /// `_merged` under those flags (`CacheDir.pm:309-315`).
    #[test]
    fn species_names_the_directory() {
        let root = tempfile::tempdir().unwrap();
        for species in ["mus_musculus", "homo_sapiens_refseq", "homo_sapiens_merged"] {
            make_leaf(
                &root.path().join(species).join("116_GRCm39"),
                Some(&INFO.replace("GRCh38", "GRCm39")),
            );
        }
        let root_s = root.path().to_str().unwrap();
        let resolved =
            resolve_flags(&["--dir_cache", root_s, "--species", "mus_musculus"]).unwrap();
        assert_eq!(resolved.dir, joined(root.path(), "mus_musculus/116_GRCm39"));
        assert_eq!(resolved.assembly.as_deref(), Some("GRCm39"));
        assert_eq!(
            resolve_flags(&["--dir_cache", root_s, "--refseq"])
                .unwrap()
                .dir,
            joined(root.path(), "homo_sapiens_refseq/116_GRCm39")
        );
        assert_eq!(
            resolve_flags(&["--dir_cache", root_s, "--merged"])
                .unwrap()
                .dir,
            joined(root.path(), "homo_sapiens_merged/116_GRCm39")
        );
        assert_eq!(
            error(&["--dir_cache", root_s, "--species", "human"]),
            format!(
                "ERROR: [cache] Cache directory {} not found",
                joined(root.path(), "human")
            )
        );
    }

    /// A species already carrying the set suffix is refused with Ensembl VEP's text
    /// (`CacheDir.pm:310-311`), its two lines joined into one.
    #[test]
    fn species_refseq_suffix_is_refused() {
        let root = make_root(&["116_GRCh38"]);
        for species in SUFFIXED_SPECIES {
            assert_eq!(
                error(&["--dir_cache", root.path().to_str().unwrap(), "--species", species]),
                format!(
                    "ERROR: [cache] Should not use {species} as --species. Try using flags --refseq or --merged with --species homo_sapiens"
                )
            );
        }
    }

    /// `--assembly` is compared exactly and case-sensitively: against the directory
    /// name under the layout (`CacheDir.pm:347-353`) and against `info.json` under
    /// `--json_cache` (`:407-408`, Ensembl VEP's text with the file renamed);
    /// `grch38` is not `GRCh38`.
    #[test]
    fn assembly_mismatch_with_info_json_is_fatal() {
        let root = make_root(&["116_GRCh38"]);
        let root_s = root.path().to_str().unwrap();
        for wrong in ["GRCh37", "grch38"] {
            assert_eq!(
                error(&["--dir_cache", root_s, "--assembly", wrong]),
                format!(
                    "ERROR: [cache] Cache assembly version (GRCh38) and database or selected assembly version ({wrong}) do not match"
                )
            );
        }
        let leaf = joined(root.path(), "homo_sapiens/116_GRCh38");
        for wrong in ["GRCh37", "grch38"] {
            assert_eq!(
                error(&["--json_cache", &leaf, "--assembly", wrong]),
                format!(
                    "ERROR: [cache] Mismatch in assembly versions from config ({wrong}) and cache info.json file (GRCh38)"
                )
            );
        }
        assert!(resolve_flags(&["--json_cache", &leaf, "--assembly", "GRCh38"]).is_ok());
        // A cache whose info.json names no assembly is checked against nothing.
        let bare = tempfile::tempdir().unwrap();
        make_leaf(bare.path(), Some(r#"{"species":"homo_sapiens"}"#));
        let resolved = resolve_flags(&[
            "--json_cache",
            bare.path().to_str().unwrap(),
            "--assembly",
            "GRCh37",
        ])
        .unwrap();
        assert_eq!(resolved.assembly.as_deref(), Some("GRCh37"));
    }

    /// An `info.json` whose `assembly` is present and empty is compared, not read
    /// as absent: `CacheDir.pm:466` keeps the empty field of an `assembly<TAB>`
    /// line, `:475` stores it, and `:407-408` refuse `--assembly GRCh38` against
    /// it with an empty value in the text. Without `--assembly` and under a
    /// directory the layout gives no assembly, nothing is compared and the run's
    /// assembly is unset.
    #[test]
    fn empty_assembly_in_info_json_is_compared_and_refused() {
        let dir = tempfile::tempdir().unwrap();
        make_leaf(
            dir.path(),
            Some(r#"{"species":"homo_sapiens","assembly":""}"#),
        );
        let path = dir.path().to_str().unwrap();
        assert_eq!(
            error(&["--json_cache", path, "--assembly", "GRCh38"]),
            "ERROR: [cache] Mismatch in assembly versions from config (GRCh38) and cache info.json file ()"
        );
        let resolved = resolve_flags(&["--json_cache", path]).unwrap();
        assert_eq!(resolved.info.assembly.as_deref(), Some(""));
        assert_eq!(resolved.assembly, None);
    }

    /// The assembly the directory name carries is the run's (`CacheDir.pm:356`)
    /// and is checked against `info.json` as `--assembly` would be (`:404-408`):
    /// a `116_GRCh38` directory whose `info.json` says `GRCh37` is refused with
    /// no `--assembly` given.
    #[test]
    fn directory_assembly_is_checked_against_info_json() {
        let root = tempfile::tempdir().unwrap();
        make_leaf(
            &root.path().join("homo_sapiens").join("116_GRCh38"),
            Some(&INFO.replace("GRCh38", "GRCh37")),
        );
        assert_eq!(
            error(&["--dir_cache", root.path().to_str().unwrap()]),
            "ERROR: [cache] Mismatch in assembly versions from config (GRCh38) and cache info.json file (GRCh37)"
        );
    }

    /// The run's assembly is the cache's: inferred from the directory name under
    /// the layout (`CacheDir.pm:356`) and read from `info.json` under `--json_cache`,
    /// with no `--assembly`; `hg19` and `hg38` stand for the canonical names.
    #[test]
    fn assembly_inferred_from_cache_feeds_plugin_hint() {
        let root = make_root(&["116_GRCh37"]);
        let root_s = root.path().to_str().unwrap();
        let resolved = resolve_flags(&["--dir_cache", root_s]).unwrap();
        assert_eq!(resolved.assembly.as_deref(), Some("GRCh37"));
        assert_eq!(
            crate::annotator::resolve_assembly(resolved.assembly.as_deref()),
            Some(vep_core::assembly::Assembly::Grch37)
        );
        let leaf = joined(root.path(), "homo_sapiens/116_GRCh37");
        assert_eq!(
            resolve_flags(&["--json_cache", &leaf])
                .unwrap()
                .assembly
                .as_deref(),
            Some("GRCh37")
        );
        for flags in [
            &["--dir_cache", root_s, "--assembly", "hg19"][..],
            &["--json_cache", &leaf, "--assembly", "hg19"][..],
        ] {
            assert_eq!(
                resolve_flags(flags).unwrap().assembly.as_deref(),
                Some("GRCh37"),
                "{flags:?}"
            );
        }
        assert_eq!(
            error(&["--dir_cache", root_s, "--assembly", "hg38"]),
            "ERROR: [cache] Cache assembly version (GRCh37) and database or selected assembly version (GRCh38) do not match"
        );
    }

    /// `--full_cache_dir` is `--json_cache` (`AnnotationSourceAdaptor.pm:128`
    /// passes it as the directory itself): the same leaf, no layout resolution.
    #[test]
    fn full_cache_dir_is_json_cache_alias() {
        let root = make_root(&["116_GRCh38"]);
        let leaf = joined(root.path(), "homo_sapiens/116_GRCh38");
        let by_alias = resolve_flags(&["--full_cache_dir", &leaf]).unwrap();
        let by_name = resolve_flags(&["--json_cache", &leaf]).unwrap();
        assert_eq!(by_alias.dir, by_name.dir);
        assert_eq!(by_alias.dir, leaf);
        assert_eq!(by_alias.assembly, by_name.assembly);
        let elsewhere = tempfile::tempdir().unwrap();
        assert_eq!(
            resolve_flags(&[
                "--full_cache_dir",
                &leaf,
                "--dir_cache",
                elsewhere.path().to_str().unwrap()
            ])
            .unwrap()
            .dir,
            leaf,
            "the directory itself wins over the root"
        );
    }

    /// A resolved layout leaf without `transcripts/` is no cache.
    #[test]
    fn leaf_without_transcripts_dir_is_fatal() {
        let root = tempfile::tempdir().unwrap();
        let leaf = root.path().join("homo_sapiens").join("116_GRCh38");
        fs::create_dir_all(&leaf).unwrap();
        fs::write(leaf.join("info.json"), INFO).unwrap();
        assert_eq!(
            error(&["--dir_cache", root.path().to_str().unwrap()]),
            format!(
                "ERROR: [cache] no cache found at {}: no transcripts directory",
                leaf.display()
            )
        );
    }

    /// A `--json_cache` directory without `transcripts/` is no cache.
    #[test]
    fn json_cache_missing_transcripts_is_fatal() {
        let dir = tempfile::tempdir().unwrap();
        fs::write(dir.path().join("info.json"), INFO).unwrap();
        let path = dir.path().to_str().unwrap();
        assert_eq!(
            error(&["--json_cache", path]),
            format!("ERROR: [cache] no cache found at {path}: no transcripts directory")
        );
    }

    /// A `--json_cache` directory without `info.json` is no cache, and a malformed
    /// `info.json` names the file.
    #[test]
    fn json_cache_missing_info_json_is_fatal() {
        let dir = tempfile::tempdir().unwrap();
        make_leaf(dir.path(), None);
        let path = dir.path().to_str().unwrap();
        assert_eq!(
            error(&["--json_cache", path]),
            format!("ERROR: [cache] no cache found at {path}: no info.json")
        );
        fs::write(dir.path().join("info.json"), "{not json").unwrap();
        let malformed = error(&["--json_cache", path]);
        assert!(
            malformed.starts_with(&format!("ERROR: [cache] {path}/info.json: ")),
            "{malformed}"
        );
    }

    /// A `--json_cache` path that is not a directory is refused with Ensembl VEP's
    /// text (`CacheDir.pm:109`).
    #[test]
    fn json_cache_typo_path_is_fatal() {
        let dir = tempfile::tempdir().unwrap();
        make_leaf(dir.path(), Some(INFO));
        let typo = format!("{}x", dir.path().display());
        assert_eq!(
            error(&["--json_cache", &typo]),
            format!("ERROR: [cache] Directory {typo} does not exist")
        );
        let file: PathBuf = dir.path().join("info.json");
        assert_eq!(
            error(&["--json_cache", file.to_str().unwrap()]),
            format!("ERROR: [cache] Directory {} does not exist", file.display())
        );
    }

    /// The first indexed `.fa` or `.fa.gz` of the cache directory, by name, is the
    /// run's FASTA (`CacheDir.pm:235-249`), unless `--fasta` names one or
    /// `--no_fasta` is given; a `.fa` without its `.fai` and a bgzipped `.fa.gz`
    /// without its `.gzi` (which Ensembl VEP indexes itself) are left alone, since
    /// the reader here takes neither.
    #[test]
    fn no_fasta_suppresses_autodetect() {
        let dir = tempfile::tempdir().unwrap();
        make_leaf(dir.path(), Some(INFO));
        for name in [
            "0.fa.gz.fai",
            "a.fa",
            "a.fa.fai",
            "b.fa",
            "b.fa.fai",
            "notes.fasta",
        ] {
            fs::write(dir.path().join(name), "").unwrap();
        }
        // The scan reads the bytes, not the name, for the .gzi requirement, so
        // this file is a BGZF stream (its EOF block alone) where the others are
        // empty.
        let bgzf = noodles::bgzf::io::Writer::new(Vec::new()).finish().unwrap();
        fs::write(dir.path().join("0.fa.gz"), bgzf).unwrap();
        let path = dir.path().to_str().unwrap();
        assert_eq!(
            resolve_flags(&["--json_cache", path]).unwrap().fasta,
            Some(dir.path().join("a.fa").display().to_string()),
            "a bgzipped .fa.gz without its .gzi is skipped"
        );
        fs::write(dir.path().join("0.fa.gz.gzi"), "").unwrap();
        assert_eq!(
            resolve_flags(&["--json_cache", path]).unwrap().fasta,
            Some(dir.path().join("0.fa.gz").display().to_string()),
            "an indexed .fa.gz is taken by name"
        );
        fs::remove_file(dir.path().join("0.fa.gz.gzi")).unwrap();
        assert!(resolve_flags(&["--json_cache", path, "--no_fasta"])
            .unwrap()
            .fasta
            .is_none());
        assert!(
            resolve_flags(&["--json_cache", path, "--fasta", "/elsewhere/ref.fa"])
                .unwrap()
                .fasta
                .is_none()
        );
        fs::remove_file(dir.path().join("a.fa.fai")).unwrap();
        assert_eq!(
            resolve_flags(&["--json_cache", path]).unwrap().fasta,
            Some(dir.path().join("b.fa").display().to_string()),
            "an unindexed .fa is skipped"
        );
        fs::remove_file(dir.path().join("b.fa")).unwrap();
        assert!(
            resolve_flags(&["--json_cache", path])
                .unwrap()
                .fasta
                .is_none(),
            "nothing indexed is left"
        );
    }

    /// The cache's `chr_synonyms.txt` is read when `--synonyms` is absent
    /// (`CacheDir.pm:254`); a `--synonyms` file replaces it whole, so a name only
    /// the cache's file links resolves to nothing under the user's file; a cache
    /// without the file gives an empty table; and a file that cannot be read is
    /// refused with the system's reason, under the component of whoever named the
    /// file. A missing file is Ensembl VEP's refusal too (`BaseVEP.pm:644`); a
    /// directory is refused here alone, where that two-argument `open` succeeds
    /// and loads an empty table, and so is `-`, a file name here, which that
    /// `open` reads as standard input (`docs/cache-setup.md` states the
    /// difference).
    #[test]
    fn synonyms_come_from_the_flag_else_the_cache_file() {
        let dir = tempfile::tempdir().unwrap();
        make_leaf(dir.path(), Some(INFO));
        let path = dir.path().to_str().unwrap();
        let holds = |name: &str| name == "21";
        assert!(resolve_flags(&["--json_cache", path])
            .unwrap()
            .synonyms
            .is_empty());

        fs::write(dir.path().join("chr_synonyms.txt"), "21\tNC_000021.9\n").unwrap();
        let cache_table = resolve_flags(&["--json_cache", path]).unwrap().synonyms;
        assert_eq!(cache_table.source_chr_name("NC_000021.9", holds), "21");

        let user_file = dir.path().join("mine.txt");
        fs::write(&user_file, "21 CM000683.2\n").unwrap();
        let user_table = resolve_flags(&[
            "--json_cache",
            path,
            "--synonyms",
            user_file.to_str().unwrap(),
        ])
        .unwrap()
        .synonyms;
        assert_eq!(user_table.source_chr_name("CM000683.2", holds), "21");
        assert_eq!(
            user_table.source_chr_name("NC_000021.9", holds),
            "NC_000021.9",
            "the cache's file is not merged in"
        );

        let absent = dir.path().join("absent.txt");
        assert_eq!(
            error(&["--json_cache", path, "--synonyms", absent.to_str().unwrap()]),
            format!(
                "ERROR: [config] Could not read synonyms file {}: No such file or directory",
                absent.display()
            )
        );
        assert_eq!(
            error(&["--json_cache", path, "--synonyms", "-"]),
            "ERROR: [config] Could not read synonyms file -: No such file or directory"
        );
        fs::remove_file(dir.path().join("chr_synonyms.txt")).unwrap();
        fs::create_dir(dir.path().join("chr_synonyms.txt")).unwrap();
        assert_eq!(
            error(&["--json_cache", path]),
            format!(
                "ERROR: [cache] Could not read synonyms file {}/chr_synonyms.txt: Is a directory",
                path
            )
        );
    }

    /// An `info.json` that is not a JSON object is refused, not read as an empty
    /// metadata set.
    #[test]
    fn info_json_must_be_an_object() {
        let dir = tempfile::tempdir().unwrap();
        make_leaf(dir.path(), Some("[]"));
        let path = dir.path().to_str().unwrap();
        assert_eq!(
            error(&["--json_cache", path]),
            format!("ERROR: [cache] {path}/info.json: not a JSON object")
        );
    }

    /// The resolution reads `CHECKSUMS.sha256` once and checks `info.json` and
    /// the cache's `chr_synonyms.txt` against it as it reads them: a changed
    /// byte in either is `checksum mismatch for <relative path> in <cache dir>`
    /// under the `cache` component, a listed file that is missing is named, a
    /// `--synonyms` file of the user's is read unchecked, a manifest of another
    /// shape names its line, and a cache without the file resolves with a
    /// manifest that is not present.
    #[test]
    fn manifest_checks_info_json_and_the_cache_synonyms() {
        let dir = tempfile::tempdir().unwrap();
        make_leaf(dir.path(), Some(INFO));
        fs::write(dir.path().join("chr_synonyms.txt"), "21\tNC_000021.9\n").unwrap();
        let path = dir.path().to_str().unwrap();
        assert!(!resolve_flags(&["--json_cache", path])
            .unwrap()
            .manifest
            .is_present());

        CacheManifest::write(dir.path()).unwrap();
        let resolved = resolve_flags(&["--json_cache", path]).unwrap();
        assert!(resolved.manifest.is_present());
        assert_eq!(resolved.manifest.len(), 2);
        assert_eq!(resolved.info.assembly.as_deref(), Some("GRCh38"));

        fs::write(dir.path().join("chr_synonyms.txt"), "21\tNC_000021.9\n\n").unwrap();
        assert_eq!(
            error(&["--json_cache", path]),
            format!("ERROR: [cache] checksum mismatch for chr_synonyms.txt in {path}")
        );
        let mine = dir.path().join("mine.txt");
        fs::write(&mine, "21\tCM000683.2\n").unwrap();
        assert!(
            resolve_flags(&["--json_cache", path, "--synonyms", mine.to_str().unwrap()]).is_ok()
        );

        fs::write(dir.path().join("info.json"), format!("{INFO} ")).unwrap();
        assert_eq!(
            error(&["--json_cache", path]),
            format!("ERROR: [cache] checksum mismatch for info.json in {path}")
        );

        fs::write(dir.path().join("info.json"), INFO).unwrap();
        fs::remove_file(dir.path().join("chr_synonyms.txt")).unwrap();
        assert_eq!(
            error(&["--json_cache", path]),
            format!(
                "ERROR: [cache] chr_synonyms.txt is listed in CHECKSUMS.sha256 but missing from {path}"
            )
        );
        fs::write(dir.path().join("chr_synonyms.txt"), "21\tNC_000021.9\n").unwrap();

        fs::write(dir.path().join("CHECKSUMS.sha256"), "x\n").unwrap();
        assert_eq!(
            error(&["--json_cache", path]),
            format!("ERROR: [cache] {path}/CHECKSUMS.sha256: line 1 is not a sha256sum line")
        );
    }

    #[test]
    fn assembly_is_read_after_the_first_digit_run_and_underscore() {
        assert_eq!(assembly_in_name("116_GRCh38"), Some("GRCh38"));
        assert_eq!(assembly_in_name("116_GRCh37_extra"), Some("GRCh37_extra"));
        assert_eq!(assembly_in_name("116"), None);
        assert_eq!(assembly_in_name("116_"), None);
        assert_eq!(assembly_in_name("v116_GRCh38"), Some("GRCh38"));
        assert_eq!(assembly_in_name("x_116_GRCh38"), Some("GRCh38"));
    }
}
