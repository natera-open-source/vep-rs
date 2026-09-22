// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Download Ensembl GFF3, protein FASTA, and variation VCF files, or resolve local paths.

use std::io::Write;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use tracing::info;

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

/// Resolve variation VCF inputs: scan a local directory or download from Ensembl FTP.
///
/// Returns `Vec<(chromosome, local_path)>`.
pub fn resolve_variation_inputs(
    local_vcf_dir: &Option<String>,
    release: u32,
    assembly: &str,
    species: &str,
    chromosomes: Option<&[String]>,
) -> Result<Vec<(String, PathBuf)>> {
    match local_vcf_dir {
        Some(dir) => {
            let dir_path = PathBuf::from(dir);
            if !dir_path.is_dir() {
                bail!(
                    "Local variation VCF directory not found: {}",
                    dir_path.display()
                );
            }
            scan_local_vcfs(&dir_path)
        }
        None => {
            let chroms: Vec<&str> = match chromosomes {
                Some(c) => c.iter().map(|s| s.as_str()).collect(),
                None => HUMAN_CHROMOSOMES.to_vec(),
            };
            download_variation_vcfs(release, assembly, species, &chroms)
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

/// Build the Ensembl FTP URL for a single chromosome's variation VCF.
///
/// Ensembl names the variation VCFs with the lowercase species directory name
/// (`homo_sapiens-chr21.vcf.gz`), unlike the capitalised GFF3 and FASTA files.
pub fn variation_vcf_url(release: u32, assembly: &str, species: &str, chromosome: &str) -> String {
    if assembly == "GRCh37" {
        format!(
            "https://ftp.ensembl.org/pub/grch37/release-{}/variation/vcf/{}/{}-chr{}.vcf.gz",
            release, species, species, chromosome
        )
    } else {
        format!(
            "https://ftp.ensembl.org/pub/release-{}/variation/vcf/{}/{}-chr{}.vcf.gz",
            release, species, species, chromosome
        )
    }
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
            let csi_url = format!("{}.csi", vcf_url);

            let vcf_dest = std::env::temp_dir().join(format!(
                "vep_cache_builder_{}_{}_{}_chr{}.vcf.gz",
                species, assembly, release, chrom
            ));
            let csi_dest = std::env::temp_dir().join(format!(
                "vep_cache_builder_{}_{}_{}_chr{}.vcf.gz.csi",
                species, assembly, release, chrom
            ));

            if vcf_dest.exists() {
                info!(
                    "Using cached variation VCF for chr{}: {}",
                    chrom,
                    vcf_dest.display()
                );
            } else {
                download_file(&vcf_url, &vcf_dest).with_context(|| {
                    format!("Failed to download variation VCF for chr{}", chrom)
                })?;
            }

            if csi_dest.exists() {
                info!(
                    "Using cached CSI index for chr{}: {}",
                    chrom,
                    csi_dest.display()
                );
            } else {
                // The CSI index is optional: some chromosomes have none.
                if let Err(e) = download_file(&csi_url, &csi_dest) {
                    info!("CSI index not available for chr{}: {}", chrom, e);
                }
            }

            Ok((chrom.to_string(), vcf_dest))
        })
        .collect();

    let mut results = results?;
    results.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(results)
}

/// Download a file from a URL to a local path.
fn download_file(url: &str, dest: &PathBuf) -> Result<()> {
    info!("Downloading: {}", url);

    let response =
        reqwest::blocking::get(url).with_context(|| format!("Failed to download {}", url))?;

    if !response.status().is_success() {
        bail!("Download failed with status {}: {}", response.status(), url);
    }

    let bytes = response
        .bytes()
        .with_context(|| format!("Failed to read response from {}", url))?;

    let mut file = std::fs::File::create(dest)
        .with_context(|| format!("Failed to create {}", dest.display()))?;
    file.write_all(&bytes)
        .with_context(|| format!("Failed to write {}", dest.display()))?;

    info!("Downloaded {} ({} bytes)", dest.display(), bytes.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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
