// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Build VEP JSON transcript cache directly from Ensembl GFF3 + protein FASTA.
//!
//! This tool downloads (or reads local) GFF3 and protein FASTA files from
//! Ensembl and produces the JSON cache directory structure consumed by
//! `vep-cli`'s `json_cache.rs` loader.
//!
//! ## Core table dumps
//!
//! With `--core-dumps <DIR>`, the Ensembl core and variation MySQL table
//! dumps of the release supply what the GFF3 lacks: the MANE accessions,
//! APPRIS, every kept transcript attribute, the CCDS, RefSeq and UniProt
//! xrefs, the gene symbol's source and HGNC id, the protein features and
//! the gene phenotype count. A table missing from the directory is
//! downloaded into it.
//!
//! ## SIFT/PolyPhen predictions
//!
//! With `--include-predictions`, connects to Ensembl's public MySQL server
//! to fetch pre-computed SIFT/PolyPhen prediction matrices. These are embedded
//! in the JSON cache alongside each protein-coding transcript.
//!
//! ## Variation files
//!
//! With `--include-variations`, Ensembl's variation VCFs (the per-chromosome
//! germline files, the genome-wide somatic, phenotype- and clinically-associated
//! files and the 1000 Genomes frequency file) are read into the variation files
//! Ensembl VEP's own cache carries: `variations/<chr>/all_vars.gz` with its
//! `.csi`, the columns listed by `info.json`'s `variation_cols`.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result};
use clap::Parser;
use tracing::{info, warn};

mod builder;
mod core_dumps;
mod download;
mod fasta;
mod gff3;
mod gtf;
mod info;
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

    /// Directory of the release's Ensembl core and variation MySQL table dumps
    /// (`<table>.txt.gz`); a table missing from it is downloaded there. Supplies
    /// the transcript attributes, xrefs, protein features and gene phenotype
    /// counts the GFF3 lacks
    #[arg(long)]
    core_dumps: Option<String>,

    /// Include variation data in the cache: Ensembl's per-chromosome variation
    /// VCFs with the genome-wide somatic, phenotype-associated,
    /// clinically-associated and 1000 Genomes frequency VCFs, downloaded or
    /// read from `--variation-vcf-dir`
    #[arg(long)]
    include_variations: bool,

    /// Skip download, use local variation VCF directory: `<species>-chr<N>.vcf.gz`
    /// per chromosome and, under their Ensembl names, the genome-wide files
    /// (`1000GENOMES-phase_3.vcf.gz`, `<species>_somatic.vcf.gz`,
    /// `<species>_phenotype_associated.vcf.gz`,
    /// `<species>_clinically_associated.vcf.gz`), each with its `.csi` index
    #[arg(long)]
    variation_vcf_dir: Option<String>,

    /// Directory of the release's Ensembl variation MySQL table dumps
    /// (`source`, `seq_region`, `attrib_type`, `phenotype`, `phenotype_feature`
    /// and `phenotype_feature_attrib`, each `.txt.gz` or `.txt`): the ClinVar
    /// allele and somatic classification columns of the variation records,
    /// which no variation VCF carries. Read with `--include-variations`
    #[arg(long)]
    variation_dumps: Option<String>,

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

/// A failed build prints one line, `ERROR: [cache-builder] <message>`, on
/// standard error and exits 1.
fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("ERROR: [cache-builder] {error:#}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<()> {
    tracing_subscriber::fmt()
        .with_max_level(tracing::Level::INFO)
        .with_writer(std::io::stderr)
        .init();

    let args = Args::parse();
    let output_dir = PathBuf::from(&args.output_dir);

    if !args.include_variations {
        let ignored: Vec<&str> = [
            ("--variation-vcf-dir", args.variation_vcf_dir.is_some()),
            ("--variation-dumps", args.variation_dumps.is_some()),
            ("--chromosomes", args.chromosomes.is_some()),
        ]
        .into_iter()
        .filter_map(|(flag, given)| given.then_some(flag))
        .collect();
        if !ignored.is_empty() {
            warn!(
                "{} read only with --include-variations; no variation record is built",
                ignored.join(", ")
            );
        }
    }

    let (gff3_path, fasta_path) = download::resolve_inputs(
        &args.gff3,
        &args.fasta,
        args.release,
        &args.assembly,
        &args.species,
    )
    .context("Failed to resolve input files")?;

    let core_dumps = if let Some(ref dir) = args.core_dumps {
        let dir = Path::new(dir);
        download::resolve_core_dumps(dir, args.release, &args.assembly, &args.species)
            .context("Failed to resolve the core table dumps")?;
        info!("Reading core table dumps: {}", dir.display());
        let dumps =
            core_dumps::CoreDumps::load(dir).context("Failed to read the core table dumps")?;
        Some(dumps)
    } else {
        None
    };

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
        core_dumps.as_ref(),
    )
    .context("Failed to build cache")?;
    info!(
        "Cache built: {} transcripts across {} chromosomes, {} region files",
        stats.transcript_count, stats.chromosome_count, stats.region_file_count
    );

    let variation_cols = if args.include_variations {
        let variation_inputs = download::resolve_variation_inputs(
            &args.variation_vcf_dir,
            args.release,
            &args.assembly,
            &args.species,
            args.chromosomes.as_deref(),
        )
        .context("Failed to resolve variation inputs")?;
        info!(
            "Resolved {} per-chromosome variation VCF files",
            variation_inputs.per_chromosome.len()
        );

        let clinvar = if let Some(ref dir) = args.variation_dumps {
            info!("Reading the ClinVar phenotype features: {}", dir);
            let chromosomes: Vec<String> = variation_inputs
                .per_chromosome
                .iter()
                .map(|(chrom, _)| chrom.clone())
                .collect();
            let attribs = variation::ClinVarAttribs::load(Path::new(dir), &chromosomes)
                .context("Failed to read the variation table dumps")?;
            Some(attribs)
        } else {
            None
        };

        let columns =
            variation::variation_columns(variation_inputs.genome_wide.frequencies.is_some());
        let mut total_variations = 0usize;
        for (chrom, path) in &variation_inputs.per_chromosome {
            info!(
                "Reading the variation records of chr{}: {}",
                chrom,
                path.display()
            );
            total_variations += variation::write_chromosome(
                chrom,
                path,
                &variation_inputs.genome_wide,
                clinvar.as_ref(),
                &columns,
                &output_dir,
            )
            .with_context(|| format!("Failed to write the variation file of chr{}", chrom))?;
        }
        info!(
            "Variation files written: {} records across {} chromosomes",
            total_variations,
            variation_inputs.per_chromosome.len()
        );
        Some(columns)
    } else {
        None
    };

    // Source versions the output headers print: the core `meta` table's, as
    // Ensembl VEP's own cache info carries them, else the GTF header's.
    let mut source_versions = core_dumps
        .as_ref()
        .map(|d| d.source_versions())
        .unwrap_or_default();
    if let Some(tags) = &gtf_tags {
        if let Some(v) = tags.header.get("genome-build") {
            source_versions
                .entry("assembly".to_string())
                .or_insert_with(|| v.clone());
        }
        if let Some(v) = tags.header.get("genebuild-last-updated") {
            source_versions
                .entry("genebuild".to_string())
                .or_insert_with(|| v.clone());
        }
    }
    if let Some(ref dir) = args.core_dumps {
        info::write_chr_synonyms(&output_dir, Path::new(dir))
            .context("Failed to write chr_synonyms.txt")?;
    }
    info::CacheInfoJson::new(
        &args.species,
        &args.assembly,
        args.release,
        &stats,
        variation_cols,
        source_versions,
    )
    .write(&output_dir)
    .context("Failed to write info.json")?;

    let listed = vep_core::cache_manifest::CacheManifest::write(&output_dir)
        .context("Failed to write CHECKSUMS.sha256")?;
    info!("Checksum manifest written for {} files", listed);

    Ok(())
}
