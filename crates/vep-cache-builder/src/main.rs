// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Build VEP JSON transcript cache directly from Ensembl GFF3 + protein FASTA.
//!
//! This tool downloads (or reads local) GFF3 and protein FASTA files from
//! Ensembl and produces the JSON cache directory structure consumed by
//! `vep-cli`'s `json_cache.rs` loader.
//!
//! ## SIFT/PolyPhen predictions
//!
//! With `--include-predictions`, connects to Ensembl's public MySQL server
//! to fetch pre-computed SIFT/PolyPhen prediction matrices. These are embedded
//! in the JSON cache alongside each protein-coding transcript.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use tracing::info;

mod builder;
mod download;
mod fasta;
mod gff3;
mod gtf;
mod mapper;
pub mod predictions;
pub mod variation;

#[derive(Parser)]
#[command(name = "vep-cache-builder")]
#[command(about = "Build VEP JSON transcript cache from Ensembl GFF3 + protein FASTA")]
struct Args {
    /// Ensembl release number (e.g., 115)
    #[arg(long)]
    release: u32,

    /// Genome assembly (e.g., GRCh37, GRCh38)
    #[arg(long)]
    assembly: String,

    /// Species name (default: homo_sapiens)
    #[arg(long, default_value = "homo_sapiens")]
    species: String,

    /// Output directory for JSON cache
    #[arg(long)]
    output_dir: String,

    /// Skip download, use local GFF3 file
    #[arg(long)]
    gff3: Option<String>,

    /// Skip download, use local protein FASTA
    #[arg(long)]
    fasta: Option<String>,

    /// Ensembl GTF of the same release, the source of the `cds_start_NF` and
    /// `cds_end_NF` transcript attributes (FLAGS in the output) that the GFF3
    /// does not carry
    #[arg(long)]
    gtf: Option<String>,

    /// Include variation data in the cache
    #[arg(long)]
    include_variations: bool,

    /// Skip download, use local variation VCF directory
    #[arg(long)]
    variation_vcf_dir: Option<String>,

    /// Specific chromosomes to include (default: all)
    #[arg(long, value_delimiter = ',')]
    chromosomes: Option<Vec<String>>,

    /// Fetch SIFT/PolyPhen prediction matrices from Ensembl MySQL
    #[arg(long)]
    include_predictions: bool,

    /// MySQL host for prediction fetching (default: ensembldb.ensembl.org)
    #[arg(long, default_value = predictions::DEFAULT_MYSQL_HOST)]
    mysql_host: String,

    /// MySQL port for prediction fetching (default: 3306)
    #[arg(long, default_value_t = predictions::DEFAULT_MYSQL_PORT)]
    mysql_port: u16,

    /// Reference genome FASTA (indexed .fa + .fai) for computing CDS sequences.
    /// When provided, enables full coding consequence resolution (missense, synonymous, etc.)
    /// by pre-computing translateable_seq for each protein-coding transcript.
    #[arg(long)]
    genome_fasta: Option<String>,
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_max_level(tracing::Level::INFO)
        .with_writer(std::io::stderr)
        .init();

    let args = Args::parse();
    let output_dir = PathBuf::from(&args.output_dir);

    let (gff3_path, fasta_path) = download::resolve_inputs(
        &args.gff3,
        &args.fasta,
        args.release,
        &args.assembly,
        &args.species,
    )
    .context("Failed to resolve input files")?;

    info!("Parsing protein FASTA: {}", fasta_path.display());
    let protein_sequences =
        fasta::parse_protein_fasta(&fasta_path).context("Failed to parse protein FASTA")?;
    info!(
        "Loaded {} protein sequences from FASTA",
        protein_sequences.len()
    );

    info!("Parsing GFF3: {}", gff3_path.display());
    let genes = gff3::parse_gff3(&gff3_path).context("Failed to parse GFF3")?;
    let transcript_count: usize = genes.values().map(|g| g.transcripts.len()).sum();
    info!(
        "Parsed {} genes, {} transcripts from GFF3",
        genes.len(),
        transcript_count
    );

    let prediction_data = if args.include_predictions {
        info!("Fetching SIFT/PolyPhen predictions from Ensembl MySQL...");
        let data = predictions::fetch_predictions(
            &protein_sequences,
            args.release,
            &args.assembly,
            &args.mysql_host,
            args.mysql_port,
        )
        .context("Failed to fetch SIFT/PolyPhen predictions")?;
        info!("Fetched predictions for {} transcripts", data.len());
        Some(data)
    } else {
        None
    };

    let genome_fasta = if let Some(ref fasta_path) = args.genome_fasta {
        info!("Loading genome FASTA: {}", fasta_path);
        let fasta = vep_fasta::IndexedFasta::from_path(fasta_path)
            .context("Failed to load genome FASTA (needs .fa + .fai index)")?;
        Some(fasta)
    } else {
        None
    };

    let gtf_tags = if let Some(ref gtf_path) = args.gtf {
        info!("Parsing GTF transcript tags: {}", gtf_path);
        let tags = gtf::parse_gtf_tags(Path::new(gtf_path)).context("Failed to parse GTF")?;
        info!("Read tags for {} transcripts from GTF", tags.tags.len());
        Some(tags)
    } else {
        None
    };

    info!("Building JSON cache in: {}", output_dir.display());
    let stats = builder::build_cache(
        &genes,
        &protein_sequences,
        &output_dir,
        prediction_data.as_ref(),
        genome_fasta.as_ref(),
        gtf_tags.as_ref(),
    )
    .context("Failed to build cache")?;
    info!(
        "Cache built: {} transcripts across {} chromosomes, {} region files",
        stats.transcript_count, stats.chromosome_count, stats.region_file_count
    );

    // Source versions the output headers print, from the GTF header when given.
    let mut source_versions = std::collections::HashMap::new();
    if let Some(tags) = &gtf_tags {
        if let Some(v) = tags.header.get("genome-build") {
            source_versions.insert("assembly".to_string(), v.clone());
        }
        if let Some(v) = tags.header.get("genebuild-last-updated") {
            source_versions.insert("genebuild".to_string(), v.clone());
        }
    }
    builder::generate_info_json(
        &output_dir,
        &args.species,
        &args.assembly,
        args.release,
        args.include_predictions,
        source_versions,
    )
    .context("Failed to generate info.json")?;

    if args.include_variations {
        let variation_inputs = download::resolve_variation_inputs(
            &args.variation_vcf_dir,
            args.release,
            &args.assembly,
            &args.species,
            args.chromosomes.as_deref(),
        )
        .context("Failed to resolve variation inputs")?;
        info!("Resolved {} variation VCF files", variation_inputs.len());

        let mut total_variations = 0usize;
        let mut total_region_files = 0usize;

        for (chrom, path) in &variation_inputs {
            info!("Parsing variation VCF for chr{}: {}", chrom, path.display());
            let records = variation::parse_variation_vcf(path, false)
                .with_context(|| format!("Failed to parse variation VCF for chr{}", chrom))?;

            let stats = variation::build_variation_cache(&records, &output_dir, chrom)
                .with_context(|| format!("Failed to build variation cache for chr{}", chrom))?;

            total_variations += stats.variation_count;
            total_region_files += stats.region_file_count;
        }

        info!(
            "Variation cache built: {} variations across {} region files",
            total_variations, total_region_files
        );
    }

    Ok(())
}
