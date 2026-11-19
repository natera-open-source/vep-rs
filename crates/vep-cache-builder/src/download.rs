// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Download Ensembl GFF3, protein FASTA, variation VCF and MySQL table dump
//! files, or resolve local paths.

use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use tracing::{info, warn};

use crate::core_dumps::{table_path, CORE_TABLES, VARIATION_TABLES};
use crate::variation::GenomeWideInputs;

/// The GFF3 to read: the local file when one is given, else the release's
/// Ensembl GFF3, downloaded from the Ensembl FTP site.
pub fn resolve_gff3(
    local: Option<&str>,
    release: u32,
    assembly: &str,
    species: &str,
) -> Result<PathBuf> {
    match local {
        Some(p) => {
            let path = PathBuf::from(p);
            if !path.exists() {
                bail!("Local GFF3 file not found: {}", path.display());
            }
            Ok(path)
        }
        None => download_gff3(release, assembly, species),
    }
}

/// The protein FASTA to read: the local file when one is given, else the
/// release's Ensembl protein FASTA, downloaded from the Ensembl FTP site.
pub fn resolve_protein_fasta(
    local: Option<&str>,
    release: u32,
    assembly: &str,
    species: &str,
) -> Result<PathBuf> {
    match local {
        Some(p) => {
            let path = PathBuf::from(p);
            if !path.exists() {
                bail!("Local protein FASTA not found: {}", path.display());
            }
            Ok(path)
        }
        None => download_protein_fasta(release, assembly, species),
    }
}

/// The Ensembl site: the default `--site`, and where a build named one input
/// at a time downloads from.
pub const ENSEMBL_SITE: &str = "https://ftp.ensembl.org";

/// The release directory of the site, relative to its root: `pub/release-N`,
/// or `pub/grch37/release-N` for GRCh37, which Ensembl publishes in a tree of
/// its own.
pub fn release_dir(release: u32, assembly: &str) -> String {
    if assembly == "GRCh37" {
        format!("pub/grch37/release-{release}")
    } else {
        format!("pub/release-{release}")
    }
}

/// The gene-set version in the GFF3 and GTF file names. Ensembl froze the
/// GRCh37 gene set at release 87: every `grch37/release-N/` directory from 88
/// onward carries `<Species>.GRCh37.87.gff3.gz`, not a file named for `N`, so
/// the version is capped at 87 there.
pub fn gene_set(release: u32, assembly: &str) -> u32 {
    if assembly == "GRCh37" {
        release.min(87)
    } else {
        release
    }
}

/// The GFF3 of a release: its directory relative to the site root and its
/// file name.
pub fn gff3_location(release: u32, assembly: &str, species: &str) -> (String, String) {
    (
        format!("{}/gff3/{species}", release_dir(release, assembly)),
        format!(
            "{}.{assembly}.{}.gff3.gz",
            capitalize_species(species),
            gene_set(release, assembly)
        ),
    )
}

/// The GTF of a release: its directory relative to the site root and its file
/// name.
pub fn gtf_location(release: u32, assembly: &str, species: &str) -> (String, String) {
    (
        format!("{}/gtf/{species}", release_dir(release, assembly)),
        format!(
            "{}.{assembly}.{}.gtf.gz",
            capitalize_species(species),
            gene_set(release, assembly)
        ),
    )
}

/// The protein FASTA of a release: its directory relative to the site root and
/// its file name.
pub fn protein_fasta_location(release: u32, assembly: &str, species: &str) -> (String, String) {
    (
        format!("{}/fasta/{species}/pep", release_dir(release, assembly)),
        format!("{}.{assembly}.pep.all.fa.gz", capitalize_species(species)),
    )
}

/// The genome FASTA directories of a release relative to the site root, in the
/// order the installer tries them: `dna_index/`, whose files come with their
/// `.fai` and `.gzi`, then `dna/`.
pub fn genome_fasta_dirs(release: u32, assembly: &str, species: &str) -> [String; 2] {
    let root = release_dir(release, assembly);
    [
        format!("{root}/fasta/{species}/dna_index"),
        format!("{root}/fasta/{species}/dna"),
    ]
}

/// The URL of the GFF3 file.
fn gff3_url(release: u32, assembly: &str, species: &str) -> String {
    let (dir, name) = gff3_location(release, assembly, species);
    format!("{ENSEMBL_SITE}/{dir}/{name}")
}

/// The URL of the protein FASTA file.
fn protein_fasta_url(release: u32, assembly: &str, species: &str) -> String {
    let (dir, name) = protein_fasta_location(release, assembly, species);
    format!("{ENSEMBL_SITE}/{dir}/{name}")
}

/// The MySQL directory of one database, relative to the site root:
/// `pub/release-N/mysql/<species>_<kind>_N_<assembly number>`, under
/// `pub/grch37/` for GRCh37, where the assembly number is the trailing digits
/// of the assembly name (`GRCh38` is 38, `GRCh37` 37).
pub fn mysql_dir(release: u32, assembly: &str, species: &str, kind: &str) -> Result<String> {
    let name_without_digits = assembly.trim_end_matches(|c: char| c.is_ascii_digit());
    let digits = &assembly[name_without_digits.len()..];
    if digits.is_empty() {
        bail!("Assembly {assembly} has no trailing number to name its MySQL database");
    }
    Ok(format!(
        "{}/mysql/{species}_{kind}_{release}_{digits}",
        release_dir(release, assembly)
    ))
}

/// The URL of one database's MySQL directory, with a trailing slash.
fn mysql_dir_url(release: u32, assembly: &str, species: &str, kind: &str) -> Result<String> {
    Ok(format!(
        "{ENSEMBL_SITE}/{}/",
        mysql_dir(release, assembly, species, kind)?
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
                "Variation table {table} is not published for this release and assembly; gene_phenotype is 0 on every transcript"
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

/// The chromosomes Ensembl publishes a per-chromosome variation VCF for under
/// `homo_sapiens`.
const HUMAN_CHROMOSOMES: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

/// The chromosomes a build without `--chromosomes` reads variation VCFs for:
/// [`HUMAN_CHROMOSOMES`] under `homo_sapiens` whatever the transcript model,
/// the GFF3's `##sequence-region` names for any other species
/// ([`sequence_region_names`]). A GTF declares no sequence regions, so another
/// species whose model is a GTF alone has no default.
pub fn default_chromosomes(species: &str, gff3: Option<&Path>) -> Result<Vec<String>> {
    if species == "homo_sapiens" {
        return Ok(HUMAN_CHROMOSOMES.iter().map(|c| c.to_string()).collect());
    }
    let Some(gff3) = gff3 else {
        bail!(
            "the transcript model is a GTF, which declares no ##sequence-region, so the chromosomes of {species} are unknown; name them with --chromosomes"
        );
    };
    let names = sequence_region_names(gff3)?;
    if names.is_empty() {
        bail!(
            "{} declares no ##sequence-region, so the chromosomes of {species} are unknown; name them with --chromosomes",
            gff3.display()
        );
    }
    Ok(names)
}

/// The names of the `##sequence-region <name> <start> <end>` pragmas in a
/// GFF3's header block, in file order: the lines before the first feature
/// line, where Ensembl's GFF3 files declare every sequence region they
/// annotate. The file is read plain or gzipped by its name.
pub fn sequence_region_names(gff3: &Path) -> Result<Vec<String>> {
    let file =
        std::fs::File::open(gff3).with_context(|| format!("Cannot open {}", gff3.display()))?;
    let reader: Box<dyn BufRead> = if gff3.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut names = Vec::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("Failed to read {}", gff3.display()))?;
        if line.is_empty() {
            continue;
        }
        if !line.starts_with('#') {
            break;
        }
        if let Some(pragma) = line.strip_prefix("##sequence-region") {
            if let Some(name) = pragma.split_whitespace().next() {
                names.push(name.to_string());
            }
        }
    }
    Ok(names)
}

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

/// The names of the genome-wide variation VCFs of `species`, the 1000 Genomes
/// file first.
pub fn genome_wide_names(species: &str) -> Vec<String> {
    let mut names = vec![THOUSAND_GENOMES_FILE.to_string()];
    names.push(format!("{species}{SOMATIC_SUFFIX}"));
    names.extend(GERMLINE_SUFFIXES.iter().map(|s| format!("{species}{s}")));
    names
}

/// Resolve the variation VCF inputs: scan a local directory or download from
/// the Ensembl FTP site. A genome-wide file the directory or the site lacks is
/// left out with a warning naming what the cache will lack; a downloaded file's
/// index (`.csi`) is fetched with it, so each chromosome's records are read
/// through the index. The chromosomes downloaded are `chromosomes`, else
/// [`default_chromosomes`] of the species and its GFF3 (`None` when the
/// transcript model is a GTF).
pub fn resolve_variation_inputs(
    local_vcf_dir: &Option<String>,
    release: u32,
    assembly: &str,
    species: &str,
    gff3: Option<&Path>,
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
            let chroms: Vec<String> = match chromosomes {
                Some(c) => c.to_vec(),
                None => default_chromosomes(species, gff3)?,
            };
            let chroms: Vec<&str> = chroms.iter().map(String::as_str).collect();
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

/// The directory of a species' variation VCFs, relative to the site root.
///
/// Ensembl names the variation VCFs with the lowercase species directory name
/// (`homo_sapiens-chr21.vcf.gz`), unlike the capitalised GFF3 and FASTA files.
pub fn variation_vcf_dir(release: u32, assembly: &str, species: &str) -> String {
    format!("{}/variation/vcf/{species}", release_dir(release, assembly))
}

/// The name of one chromosome's variation VCF in that directory.
pub fn variation_vcf_name(species: &str, chromosome: &str) -> String {
    format!("{species}-chr{chromosome}.vcf.gz")
}

/// The URL of a species' variation VCF directory, with a trailing slash.
fn variation_vcf_dir_url(release: u32, assembly: &str, species: &str) -> String {
    format!(
        "{ENSEMBL_SITE}/{}/",
        variation_vcf_dir(release, assembly, species)
    )
}

/// Build the Ensembl FTP URL for a single chromosome's variation VCF.
pub fn variation_vcf_url(release: u32, assembly: &str, species: &str, chromosome: &str) -> String {
    format!(
        "{}{}",
        variation_vcf_dir_url(release, assembly, species),
        variation_vcf_name(species, chromosome)
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
pub fn is_not_found(error: &anyhow::Error) -> bool {
    error
        .downcast_ref::<HttpStatus>()
        .is_some_and(|status| status.0 == reqwest::StatusCode::NOT_FOUND)
}

/// Sends a GET for `url` and returns the response once its status is a
/// success; any other status is an [`HttpStatus`] error. The client has a
/// connect timeout, TCP keepalive against a dead peer and a four-hour bound on
/// the whole transfer, about the time a 1.5 GB table dump takes at 100 KB/s.
fn http_get(url: &str) -> Result<reqwest::blocking::Response> {
    let client = reqwest::blocking::Client::builder()
        .connect_timeout(std::time::Duration::from_secs(60))
        .tcp_keepalive(std::time::Duration::from_secs(60))
        .timeout(std::time::Duration::from_secs(4 * 3600))
        .build()
        .context("Failed to build the HTTP client")?;
    let response = client
        .get(url)
        .send()
        .with_context(|| format!("Failed to download {}", url))?;
    if !response.status().is_success() {
        return Err(anyhow::Error::new(HttpStatus(response.status()))
            .context(format!("Download failed: {}", url)));
    }
    Ok(response)
}

/// The text of the file at `url`, or `None` when the server answers 404.
pub fn fetch_text(url: &str) -> Result<Option<String>> {
    let response = match http_get(url) {
        Ok(response) => response,
        Err(e) if is_not_found(&e) => return Ok(None),
        Err(e) => return Err(e),
    };
    let text = response
        .text()
        .with_context(|| format!("Failed to read response from {}", url))?;
    Ok(Some(text))
}

/// Download a file from a URL to a local path, streaming the body into
/// `<dest>.part` and renaming it at the end, so an interrupted transfer never
/// leaves a truncated file under the final name (the part file is removed).
pub fn download_file(url: &str, dest: &Path) -> Result<()> {
    info!("Downloading: {}", url);

    let mut response = http_get(url)?;

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

    /// The release tree the installer walks: `pub/release-N` or
    /// `pub/grch37/release-N`, the GRCh37 gene set capped at 87 in the GFF3 and
    /// GTF names, the species capitalised in every file name but the variation
    /// VCFs', and `dna_index/` before `dna/`.
    #[test]
    fn release_tree_locations() {
        assert_eq!(
            gff3_location(116, "GRCh38", "homo_sapiens"),
            (
                "pub/release-116/gff3/homo_sapiens".to_string(),
                "Homo_sapiens.GRCh38.116.gff3.gz".to_string()
            )
        );
        assert_eq!(
            gtf_location(116, "GRCh37", "homo_sapiens"),
            (
                "pub/grch37/release-116/gtf/homo_sapiens".to_string(),
                "Homo_sapiens.GRCh37.87.gtf.gz".to_string()
            )
        );
        assert_eq!(
            protein_fasta_location(116, "GRCm39", "mus_musculus"),
            (
                "pub/release-116/fasta/mus_musculus/pep".to_string(),
                "Mus_musculus.GRCm39.pep.all.fa.gz".to_string()
            )
        );
        assert_eq!(
            genome_fasta_dirs(116, "GRCh37", "homo_sapiens"),
            [
                "pub/grch37/release-116/fasta/homo_sapiens/dna_index".to_string(),
                "pub/grch37/release-116/fasta/homo_sapiens/dna".to_string()
            ]
        );
        assert_eq!(
            mysql_dir(116, "GRCh38", "homo_sapiens", "variation").unwrap(),
            "pub/release-116/mysql/homo_sapiens_variation_116_38"
        );
        assert_eq!(
            variation_vcf_dir(116, "GRCh38", "homo_sapiens"),
            "pub/release-116/variation/vcf/homo_sapiens"
        );
        assert_eq!(
            variation_vcf_name("homo_sapiens", "MT"),
            "homo_sapiens-chrMT.vcf.gz"
        );
        assert_eq!(
            genome_wide_names("homo_sapiens"),
            [
                "1000GENOMES-phase_3.vcf.gz",
                "homo_sapiens_somatic.vcf.gz",
                "homo_sapiens_phenotype_associated.vcf.gz",
                "homo_sapiens_clinically_associated.vcf.gz"
            ]
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

    /// A GFF3 header in Ensembl's layout: the version pragma, one
    /// `##sequence-region` per region, the `#!` metadata lines, then the
    /// features. The names come out in file order, a pragma after the first
    /// feature line is not read, and a gzipped file reads the same.
    #[test]
    fn sequence_region_names_follow_the_header_block() {
        let dir = tempfile::tempdir().unwrap();
        let text = "##gff-version 3\n\
                    ##sequence-region   2L 1 23513712\n\
                    ##sequence-region   X 1 23542271\n\
                    #!genome-build FlyBase BDGP6.54\n\
                    #!genome-version BDGP6.54\n\
                    \n\
                    2L\tFlyBase\tchromosome\t1\t23513712\t.\t.\t.\tID=chromosome:2L\n\
                    ##sequence-region   Y 1 3667352\n";
        let plain = dir.path().join("fly.gff3");
        std::fs::write(&plain, text).unwrap();
        assert_eq!(sequence_region_names(&plain).unwrap(), ["2L", "X"]);

        let gz = dir.path().join("fly.gff3.gz");
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(&gz).unwrap(),
            flate2::Compression::fast(),
        );
        encoder.write_all(text.as_bytes()).unwrap();
        encoder.finish().unwrap();
        assert_eq!(sequence_region_names(&gz).unwrap(), ["2L", "X"]);

        let featureless = dir.path().join("empty.gff3");
        std::fs::write(&featureless, "##gff-version 3\n").unwrap();
        assert!(sequence_region_names(&featureless).unwrap().is_empty());
        assert!(sequence_region_names(&dir.path().join("missing.gff3")).is_err());
    }

    /// Without `--chromosomes`, `homo_sapiens` reads the 25 chromosomes
    /// Ensembl publishes per-chromosome VCFs for, whatever its GFF3 declares;
    /// another species reads its GFF3's `##sequence-region` names, and a GFF3
    /// declaring none is an error naming the flag that supplies them.
    #[test]
    fn default_chromosomes_follow_the_species() {
        let dir = tempfile::tempdir().unwrap();
        let gff3 = dir.path().join("genes.gff3");
        std::fs::write(
            &gff3,
            "##gff-version 3\n##sequence-region   2L 1 23513712\n##sequence-region   mitochondrion_genome 1 19524\n",
        )
        .unwrap();
        let human = default_chromosomes("homo_sapiens", Some(&gff3)).unwrap();
        assert_eq!(human.len(), 25);
        assert_eq!(human[0], "1");
        assert_eq!(
            default_chromosomes("drosophila_melanogaster", Some(&gff3)).unwrap(),
            ["2L", "mitochondrion_genome"]
        );
        let bare = dir.path().join("bare.gff3");
        std::fs::write(&bare, "##gff-version 3\n").unwrap();
        let err = default_chromosomes("drosophila_melanogaster", Some(&bare))
            .unwrap_err()
            .to_string();
        assert!(err.contains("--chromosomes"), "{err}");
        assert!(
            err.contains("chromosomes of drosophila_melanogaster are unknown"),
            "{err}"
        );
    }

    /// With a GTF as the transcript model there is no GFF3 to read sequence
    /// regions from, so another species has no default chromosomes: the error
    /// says a GTF declares none and names the flag that supplies them.
    #[test]
    fn default_chromosomes_without_a_gff3_refuse_another_species() {
        let err = default_chromosomes("drosophila_melanogaster", None)
            .unwrap_err()
            .to_string();
        assert!(err.contains("--chromosomes"), "{err}");
        assert!(err.contains("GTF"), "{err}");
        assert!(err.contains("declares no ##sequence-region"), "{err}");
        assert!(
            err.contains("chromosomes of drosophila_melanogaster are unknown"),
            "{err}"
        );
    }

    /// `homo_sapiens` keeps its fixed 25 chromosomes without a GFF3, so a
    /// human cache built from a GTF reads the same variation VCFs as one
    /// built from the GFF3.
    #[test]
    fn default_chromosomes_without_a_gff3_are_human_for_homo_sapiens() {
        let human = default_chromosomes("homo_sapiens", None).unwrap();
        assert_eq!(human.len(), 25);
        assert_eq!(human[0], "1");
        assert_eq!(human[22], "X");
        assert_eq!(human[24], "MT");
    }
}
