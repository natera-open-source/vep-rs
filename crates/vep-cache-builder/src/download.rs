// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Download Ensembl GFF3, protein FASTA, variation VCF and MySQL table dump
//! files, or resolve local paths.

use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use tracing::{info, warn};

use crate::core_dumps::{table_path, CORE_TABLES, VARIATION_TABLES};
use crate::variation::GenomeWideInputs;

/// Resolve input file paths: download from Ensembl FTP if not provided locally.
///
/// Returns `(gff3_path, fasta_path)`.
pub fn resolve_inputs(
    local_gff3: &Option<String>,
    local_fasta: &Option<String>,
    release: u32,
    assembly: &str,
    species: &str,
) -> Result<(PathBuf, PathBuf)> {
    let gff3_path = match local_gff3 {
        Some(p) => {
            let path = PathBuf::from(p);
            if !path.exists() {
                bail!("Local GFF3 file not found: {}", path.display());
            }
            path
        }
        None => download_gff3(release, assembly, species)?,
    };

    let fasta_path = match local_fasta {
        Some(p) => {
            let path = PathBuf::from(p);
            if !path.exists() {
                bail!("Local protein FASTA not found: {}", path.display());
            }
            path
        }
        None => download_protein_fasta(release, assembly, species)?,
    };

    Ok((gff3_path, fasta_path))
}

/// Build the Ensembl FTP URL for the GFF3 file.
///
/// Ensembl froze the GRCh37 gene set at release 87: every `grch37/release-N/`
/// directory from 88 onward carries `<Species>.GRCh37.87.gff3.gz`, not a file
/// named for `N`, so the gene-set version in the file name is capped at 87.
fn gff3_url(release: u32, assembly: &str, species: &str) -> String {
    let species_cap = capitalize_species(species);
    if assembly == "GRCh37" {
        let gene_set = release.min(87);
        format!(
            "https://ftp.ensembl.org/pub/grch37/release-{}/gff3/{}/{}.GRCh37.{}.gff3.gz",
            release, species, species_cap, gene_set
        )
    } else {
        format!(
            "https://ftp.ensembl.org/pub/release-{}/gff3/{}/{}.{}.{}.gff3.gz",
            release, species, species_cap, assembly, release
        )
    }
}

/// Build the Ensembl FTP URL for the protein FASTA file.
fn protein_fasta_url(release: u32, assembly: &str, species: &str) -> String {
    let species_cap = capitalize_species(species);
    if assembly == "GRCh37" {
        format!(
            "https://ftp.ensembl.org/pub/grch37/release-{}/fasta/{}/pep/{}.GRCh37.pep.all.fa.gz",
            release, species, species_cap
        )
    } else {
        format!(
            "https://ftp.ensembl.org/pub/release-{}/fasta/{}/pep/{}.{}.pep.all.fa.gz",
            release, species, species_cap, assembly
        )
    }
}

/// The MySQL directory of one database on the Ensembl FTP site:
/// `pub/release-N/mysql/<species>_<kind>_N_<assembly number>/`, under
/// `pub/grch37/` for GRCh37, where the assembly number is the trailing digits
/// of the assembly name (`GRCh38` is 38, `GRCh37` 37).
fn mysql_dir_url(release: u32, assembly: &str, species: &str, kind: &str) -> Result<String> {
    let name_without_digits = assembly.trim_end_matches(|c: char| c.is_ascii_digit());
    let digits = &assembly[name_without_digits.len()..];
    if digits.is_empty() {
        bail!("Assembly {assembly} has no trailing number to name its MySQL database");
    }
    let tree = if assembly == "GRCh37" {
        "pub/grch37"
    } else {
        "pub"
    };
    Ok(format!(
        "https://ftp.ensembl.org/{tree}/release-{release}/mysql/{species}_{kind}_{release}_{digits}/"
    ))
}

/// The URL of one core table dump (`<table>.txt.gz`).
pub fn core_dump_url(release: u32, assembly: &str, species: &str, table: &str) -> Result<String> {
    Ok(format!(
        "{}{table}.txt.gz",
        mysql_dir_url(release, assembly, species, "core")?
    ))
}

/// The URL of one variation table dump (`<table>.txt.gz`); a release may
/// publish the GRCh37 tables uncompressed, so a caller falls back to `.txt`.
pub fn variation_dump_url(
    release: u32,
    assembly: &str,
    species: &str,
    table: &str,
) -> Result<String> {
    Ok(format!(
        "{}{table}.txt.gz",
        mysql_dir_url(release, assembly, species, "variation")?
    ))
}

/// Make `dir` hold every table dump the cache build reads: a core table
/// missing from it is downloaded (an error when that fails), and a missing
/// variation table is downloaded as `.txt.gz`, then as `.txt` when the first
/// is a 404, and left absent with a warning when both are 404, so a species
/// without a variation database still builds. Any other failure stops the
/// build.
pub fn resolve_core_dumps(dir: &Path, release: u32, assembly: &str, species: &str) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("Failed to create {}", dir.display()))?;
    for table in CORE_TABLES {
        if table_path(dir, table).is_some() {
            continue;
        }
        let url = core_dump_url(release, assembly, species, table)?;
        download_file(&url, &dir.join(format!("{table}.txt.gz")))
            .with_context(|| format!("Failed to download core table {table}"))?;
    }
    for table in VARIATION_TABLES {
        if table_path(dir, table).is_some() {
            continue;
        }
        let gz_url = variation_dump_url(release, assembly, species, table)?;
        match download_file(&gz_url, &dir.join(format!("{table}.txt.gz"))) {
            Ok(()) => continue,
            Err(e) if is_not_found(&e) => {}
            Err(e) => return Err(e.context(format!("Failed to download variation table {table}"))),
        }
        let txt_url = gz_url.trim_end_matches(".gz");
        match download_file(txt_url, &dir.join(format!("{table}.txt"))) {
            Ok(()) => {}
            Err(e) if is_not_found(&e) => warn!(
                "Variation table {table} is not published for this release and assembly; gene_phenotype is not written"
            ),
            Err(e) => return Err(e.context(format!("Failed to download variation table {table}"))),
        }
    }
    Ok(())
}

/// Capitalize species name for Ensembl file naming (e.g., "homo_sapiens" -> "Homo_sapiens").
fn capitalize_species(species: &str) -> String {
    let mut chars = species.chars();
    match chars.next() {
        None => String::new(),
        Some(first) => {
            let upper: String = first.to_uppercase().collect();
            format!("{}{}", upper, chars.as_str())
        }
    }
}

/// Download GFF3 from Ensembl FTP.
fn download_gff3(release: u32, assembly: &str, species: &str) -> Result<PathBuf> {
    let url = gff3_url(release, assembly, species);
    let dest = std::env::temp_dir().join(format!(
        "vep_cache_builder_{}_{}_{}.gff3.gz",
        species, assembly, release
    ));

    if dest.exists() {
        info!("Using cached GFF3: {}", dest.display());
        return Ok(dest);
    }

    download_file(&url, &dest)?;
    Ok(dest)
}

/// Download protein FASTA from Ensembl FTP.
fn download_protein_fasta(release: u32, assembly: &str, species: &str) -> Result<PathBuf> {
    let url = protein_fasta_url(release, assembly, species);
    let dest = std::env::temp_dir().join(format!(
        "vep_cache_builder_{}_{}_{}.pep.all.fa.gz",
        species, assembly, release
    ));

    if dest.exists() {
        info!("Using cached protein FASTA: {}", dest.display());
        return Ok(dest);
    }

    download_file(&url, &dest)?;
    Ok(dest)
}

/// Default chromosomes for homo_sapiens.
const HUMAN_CHROMOSOMES: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

/// The variation VCF inputs of a build: one per-chromosome germline VCF per
/// chromosome, and the genome-wide files read beside them.
pub struct VariationInputs {
    /// `(chromosome, path)` of each per-chromosome VCF, by chromosome name.
    pub per_chromosome: Vec<(String, PathBuf)>,
    pub genome_wide: GenomeWideInputs,
}

/// The genome-wide variation VCFs Ensembl publishes beside the per-chromosome
/// files, by the name each has in the FTP directory and in a local directory.
/// The 1000 Genomes file is human-only and exists for GRCh38 alone.
const THOUSAND_GENOMES_FILE: &str = "1000GENOMES-phase_3.vcf.gz";
const SOMATIC_SUFFIX: &str = "_somatic.vcf.gz";
const GERMLINE_SUFFIXES: [&str; 2] = [
    "_phenotype_associated.vcf.gz",
    "_clinically_associated.vcf.gz",
];

/// Resolve the variation VCF inputs: scan a local directory or download from
/// the Ensembl FTP site. A genome-wide file the directory or the site lacks is
/// left out with a warning naming what the cache will lack; a downloaded file's
/// index (`.csi`) is fetched with it, so each chromosome's records are read
/// through the index.
pub fn resolve_variation_inputs(
    local_vcf_dir: &Option<String>,
    release: u32,
    assembly: &str,
    species: &str,
    chromosomes: Option<&[String]>,
) -> Result<VariationInputs> {
    match local_vcf_dir {
        Some(dir) => {
            let dir_path = PathBuf::from(dir);
            if !dir_path.is_dir() {
                bail!(
                    "Local variation VCF directory not found: {}",
                    dir_path.display()
                );
            }
            Ok(VariationInputs {
                per_chromosome: scan_local_vcfs(&dir_path)?,
                genome_wide: scan_local_genome_wide(&dir_path, species)?,
            })
        }
        None => {
            let chroms: Vec<&str> = match chromosomes {
                Some(c) => c.iter().map(|s| s.as_str()).collect(),
                None => HUMAN_CHROMOSOMES.to_vec(),
            };
            Ok(VariationInputs {
                per_chromosome: download_variation_vcfs(release, assembly, species, &chroms)?,
                genome_wide: download_genome_wide(release, assembly, species)?,
            })
        }
    }
}

/// Scan a local directory for variation VCF files matching `*-chr*.vcf.gz`.
fn scan_local_vcfs(dir: &PathBuf) -> Result<Vec<(String, PathBuf)>> {
    let mut results = Vec::new();
    for entry in std::fs::read_dir(dir)
        .with_context(|| format!("Failed to read directory: {}", dir.display()))?
    {
        let entry = entry?;
        let path = entry.path();
        let name = match path.file_name().and_then(|n| n.to_str()) {
            Some(n) => n.to_string(),
            None => continue,
        };
        if name.ends_with(".vcf.gz") && name.contains("-chr") {
            if let Some(chrom) = extract_chromosome_from_filename(&name) {
                info!(
                    "Found local variation VCF for chr{}: {}",
                    chrom,
                    path.display()
                );
                results.push((chrom, path));
            }
        }
    }
    if results.is_empty() {
        bail!(
            "No variation VCF files matching *-chr*.vcf.gz found in {}",
            dir.display()
        );
    }
    results.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(results)
}

/// The genome-wide files present in a local directory under their Ensembl names.
fn scan_local_genome_wide(dir: &Path, species: &str) -> Result<GenomeWideInputs> {
    genome_wide_inputs(species, |name| {
        let path = dir.join(name);
        if !path.is_file() {
            return Ok(None);
        }
        info!("Found local variation VCF: {}", path.display());
        Ok(Some(path))
    })
}

/// The genome-wide inputs of `species`, each file looked up by its Ensembl name
/// through `locate`: the 1000 Genomes file, the somatic file, then the germline
/// files. A file `locate` does not find is left out with a warning naming what
/// the cache lacks without it.
fn genome_wide_inputs(
    species: &str,
    mut locate: impl FnMut(&str) -> Result<Option<PathBuf>>,
) -> Result<GenomeWideInputs> {
    let frequencies = locate(THOUSAND_GENOMES_FILE)?;
    let somatic = locate(&format!("{species}{SOMATIC_SUFFIX}"))?;
    let mut germline = Vec::new();
    for suffix in GERMLINE_SUFFIXES {
        let name = format!("{species}{suffix}");
        match locate(&name)? {
            Some(path) => germline.push(path),
            None => warn!("No {name}: a variant only it carries is absent from the cache"),
        }
    }
    if frequencies.is_none() {
        warn!("No {THOUSAND_GENOMES_FILE}: the cache carries no 1000 Genomes frequency column");
    }
    if somatic.is_none() {
        warn!("No somatic variation VCF: no record of the cache is flagged somatic");
    }
    Ok(GenomeWideInputs {
        frequencies,
        somatic,
        germline,
    })
}

/// Extract chromosome name from a filename like "homo_sapiens-chr21.vcf.gz".
fn extract_chromosome_from_filename(filename: &str) -> Option<String> {
    let stem = filename.strip_suffix(".vcf.gz")?;
    let idx = stem.rfind("-chr")?;
    let chr_part = &stem[idx + 4..];
    if chr_part.is_empty() {
        return None;
    }
    Some(chr_part.to_string())
}

/// The Ensembl FTP directory of a species' variation VCFs.
///
/// Ensembl names the variation VCFs with the lowercase species directory name
/// (`homo_sapiens-chr21.vcf.gz`), unlike the capitalised GFF3 and FASTA files.
fn variation_vcf_dir_url(release: u32, assembly: &str, species: &str) -> String {
    if assembly == "GRCh37" {
        format!("https://ftp.ensembl.org/pub/grch37/release-{release}/variation/vcf/{species}/")
    } else {
        format!("https://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/")
    }
}

/// Build the Ensembl FTP URL for a single chromosome's variation VCF.
pub fn variation_vcf_url(release: u32, assembly: &str, species: &str, chromosome: &str) -> String {
    format!(
        "{}{species}-chr{chromosome}.vcf.gz",
        variation_vcf_dir_url(release, assembly, species)
    )
}

/// Download per-chromosome variation VCFs from Ensembl FTP.
///
/// Downloads each chromosome's VCF and its `.csi` index in parallel using rayon.
/// Returns `Vec<(chromosome, local_path)>`.
pub fn download_variation_vcfs(
    release: u32,
    assembly: &str,
    species: &str,
    chromosomes: &[&str],
) -> Result<Vec<(String, PathBuf)>> {
    info!(
        "Downloading variation VCFs for {} chromosomes ({} release-{} {})",
        chromosomes.len(),
        species,
        release,
        assembly
    );

    let results: Result<Vec<(String, PathBuf)>> = chromosomes
        .par_iter()
        .map(|chrom| {
            let vcf_url = variation_vcf_url(release, assembly, species, chrom);
            let vcf_dest = std::env::temp_dir().join(format!(
                "vep_cache_builder_{}_{}_{}_chr{}.vcf.gz",
                species, assembly, release, chrom
            ));
            download_cached(&vcf_url, &vcf_dest)
                .with_context(|| format!("Failed to download variation VCF for chr{}", chrom))?;
            // The CSI index is optional: some chromosomes have none.
            if let Err(e) = download_index(&vcf_url, &vcf_dest) {
                info!("CSI index not available for chr{}: {}", chrom, e);
            }
            Ok((chrom.to_string(), vcf_dest))
        })
        .collect();

    let mut results = results?;
    results.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(results)
}

/// Downloads the genome-wide variation VCFs and their indexes; a file the site
/// does not publish for this release and assembly (the 1000 Genomes file on
/// GRCh37) is left out.
fn download_genome_wide(release: u32, assembly: &str, species: &str) -> Result<GenomeWideInputs> {
    let dir_url = variation_vcf_dir_url(release, assembly, species);
    genome_wide_inputs(species, |name| {
        let url = format!("{dir_url}{name}");
        let dest = std::env::temp_dir().join(format!(
            "vep_cache_builder_{species}_{assembly}_{release}_{name}"
        ));
        match download_cached(&url, &dest) {
            Ok(()) => {}
            Err(e) if is_not_found(&e) => return Ok(None),
            Err(e) => return Err(e.context(format!("Failed to download {name}"))),
        }
        download_index(&url, &dest)
            .with_context(|| format!("Failed to download the index of {name}"))?;
        Ok(Some(dest))
    })
}

/// Downloads `url` to `dest` unless `dest` already exists.
fn download_cached(url: &str, dest: &Path) -> Result<()> {
    if dest.exists() {
        info!("Using cached download: {}", dest.display());
        return Ok(());
    }
    download_file(url, dest)
}

/// Downloads the `.csi` index of the file at `url` beside `dest`.
fn download_index(url: &str, dest: &Path) -> Result<()> {
    let mut index_dest = dest.as_os_str().to_owned();
    index_dest.push(".csi");
    download_cached(&format!("{url}.csi"), Path::new(&index_dest))
}

/// An HTTP response status other than success, kept as a typed error so a
/// caller can tell a missing file (404) from a failed transfer.
#[derive(Debug)]
struct HttpStatus(reqwest::StatusCode);

impl std::fmt::Display for HttpStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HTTP status {}", self.0)
    }
}

impl std::error::Error for HttpStatus {}

/// Whether an error is an HTTP 404 from [`download_file`].
fn is_not_found(error: &anyhow::Error) -> bool {
    error
        .downcast_ref::<HttpStatus>()
        .is_some_and(|status| status.0 == reqwest::StatusCode::NOT_FOUND)
}

/// Download a file from a URL to a local path, streaming the body into
/// `<dest>.part` and renaming it at the end, so an interrupted transfer never
/// leaves a truncated file under the final name (the part file is removed).
/// The client has a connect timeout, TCP keepalive against a dead peer and a
/// four-hour bound on the whole transfer, about the time a 1.5 GB table dump
/// takes at 100 KB/s.
fn download_file(url: &str, dest: &Path) -> Result<()> {
    info!("Downloading: {}", url);

    let client = reqwest::blocking::Client::builder()
        .connect_timeout(std::time::Duration::from_secs(60))
        .tcp_keepalive(std::time::Duration::from_secs(60))
        .timeout(std::time::Duration::from_secs(4 * 3600))
        .build()
        .context("Failed to build the HTTP client")?;
    let mut response = client
        .get(url)
        .send()
        .with_context(|| format!("Failed to download {}", url))?;
    if !response.status().is_success() {
        return Err(anyhow::Error::new(HttpStatus(response.status()))
            .context(format!("Download failed: {}", url)));
    }

    let part = dest.with_file_name(format!(
        "{}.part",
        dest.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("download")
    ));
    let mut file = std::fs::File::create(&part)
        .with_context(|| format!("Failed to create {}", part.display()))?;
    let written = response
        .copy_to(&mut file)
        .with_context(|| format!("Failed to read response from {}", url))
        .and_then(|bytes| {
            file.flush()
                .with_context(|| format!("Failed to write {}", part.display()))?;
            Ok(bytes)
        });
    let bytes = match written {
        Ok(bytes) => bytes,
        Err(e) => {
            let _ = std::fs::remove_file(&part);
            return Err(e);
        }
    };
    std::fs::rename(&part, dest)
        .with_context(|| format!("Failed to move {} into place", part.display()))?;

    info!("Downloaded {} ({} bytes)", dest.display(), bytes);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_a_404_reads_as_not_found() {
        let not_found = anyhow::Error::new(HttpStatus(reqwest::StatusCode::NOT_FOUND))
            .context("Download failed: x");
        assert!(is_not_found(&not_found));
        let forbidden = anyhow::Error::new(HttpStatus(reqwest::StatusCode::FORBIDDEN))
            .context("Download failed: x");
        assert!(!is_not_found(&forbidden));
        assert!(!is_not_found(&anyhow::anyhow!("connection reset")));
    }

    #[test]
    fn test_gff3_url_grch38() {
        let url = gff3_url(115, "GRCh38", "homo_sapiens");
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz"
        );
    }

    #[test]
    fn test_gff3_url_grch37() {
        let url = gff3_url(115, "GRCh37", "homo_sapiens");
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/grch37/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz"
        );
        let url = gff3_url(87, "GRCh37", "homo_sapiens");
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz"
        );
    }

    #[test]
    fn test_protein_fasta_url_grch38() {
        let url = protein_fasta_url(115, "GRCh38", "homo_sapiens");
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
        );
    }

    #[test]
    fn test_protein_fasta_url_grch37() {
        let url = protein_fasta_url(115, "GRCh37", "homo_sapiens");
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/grch37/release-115/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz"
        );
    }

    #[test]
    fn core_and_variation_dump_urls() {
        assert_eq!(
            core_dump_url(116, "GRCh38", "homo_sapiens", "transcript_attrib").unwrap(),
            "https://ftp.ensembl.org/pub/release-116/mysql/homo_sapiens_core_116_38/transcript_attrib.txt.gz"
        );
        assert_eq!(
            core_dump_url(116, "GRCh37", "homo_sapiens", "gene").unwrap(),
            "https://ftp.ensembl.org/pub/grch37/release-116/mysql/homo_sapiens_core_116_37/gene.txt.gz"
        );
        assert_eq!(
            variation_dump_url(116, "GRCh38", "homo_sapiens", "phenotype_feature").unwrap(),
            "https://ftp.ensembl.org/pub/release-116/mysql/homo_sapiens_variation_116_38/phenotype_feature.txt.gz"
        );
        assert_eq!(
            core_dump_url(116, "GRCm39", "mus_musculus", "meta").unwrap(),
            "https://ftp.ensembl.org/pub/release-116/mysql/mus_musculus_core_116_39/meta.txt.gz"
        );
        assert!(core_dump_url(116, "hg", "homo_sapiens", "meta").is_err());
    }

    #[test]
    fn test_capitalize_species() {
        assert_eq!(capitalize_species("homo_sapiens"), "Homo_sapiens");
        assert_eq!(capitalize_species("mus_musculus"), "Mus_musculus");
        assert_eq!(capitalize_species(""), "");
    }

    #[test]
    fn test_variation_vcf_url_grch38() {
        let url = variation_vcf_url(115, "GRCh38", "homo_sapiens", "21");
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/release-115/variation/vcf/homo_sapiens/homo_sapiens-chr21.vcf.gz"
        );
    }

    #[test]
    fn test_variation_vcf_url_grch37() {
        let url = variation_vcf_url(115, "GRCh37", "homo_sapiens", "21");
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/grch37/release-115/variation/vcf/homo_sapiens/homo_sapiens-chr21.vcf.gz"
        );
    }

    #[test]
    fn test_variation_vcf_url_x_chromosome() {
        let url = variation_vcf_url(115, "GRCh38", "homo_sapiens", "X");
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/release-115/variation/vcf/homo_sapiens/homo_sapiens-chrX.vcf.gz"
        );
    }

    #[test]
    fn test_variation_vcf_url_mt() {
        let url = variation_vcf_url(115, "GRCh37", "homo_sapiens", "MT");
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/grch37/release-115/variation/vcf/homo_sapiens/homo_sapiens-chrMT.vcf.gz"
        );
    }

    #[test]
    fn test_extract_chromosome_from_filename() {
        assert_eq!(
            extract_chromosome_from_filename("Homo_sapiens-chr21.vcf.gz"),
            Some("21".to_string())
        );
        assert_eq!(
            extract_chromosome_from_filename("Homo_sapiens-chrX.vcf.gz"),
            Some("X".to_string())
        );
        assert_eq!(
            extract_chromosome_from_filename("Homo_sapiens-chrMT.vcf.gz"),
            Some("MT".to_string())
        );
        assert_eq!(extract_chromosome_from_filename("random_file.txt"), None);
        // A `.vcf.gz` name without a `-chr` marker yields None.
        assert_eq!(extract_chromosome_from_filename("nothing.vcf.gz"), None);
    }

    #[test]
    fn test_human_chromosomes_count() {
        assert_eq!(HUMAN_CHROMOSOMES.len(), 25);
        assert_eq!(HUMAN_CHROMOSOMES[0], "1");
        assert_eq!(HUMAN_CHROMOSOMES[21], "22");
        assert_eq!(HUMAN_CHROMOSOMES[22], "X");
        assert_eq!(HUMAN_CHROMOSOMES[23], "Y");
        assert_eq!(HUMAN_CHROMOSOMES[24], "MT");
    }
}
