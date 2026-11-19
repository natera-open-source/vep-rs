// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Build VEP JSON transcript cache directly from a GFF3 or GTF transcript
//! model and a protein or genome FASTA.
//!
//! This tool downloads (or reads local) GFF3 and protein FASTA files from
//! Ensembl and produces the JSON cache directory structure consumed by
//! `vep-cli`'s `json_cache.rs` loader. A local `--gff3` or `--gtf` may be any
//! GFF3 or GTF, read as Ensembl VEP reads a `--gff` or `--gtf` annotation
//! source; a build reading the release's protein FASTA or core dumps reads an
//! Ensembl GFF3 for the cache's content instead ([`gff3::Reading`]). With
//! `--genome-fasta`, no `--fasta` and no `--include-predictions` the peptides are translated
//! from the genome, so a cache builds from a transcript model and a genome alone.
//!
//! ## One-command installation
//!
//! `--auto <species> <assembly> --release N` fetches every input of the release
//! from the Ensembl site (or `--site`), checks each against the `CHECKSUMS`
//! file of its directory, builds the cache as
//! `<dir_cache>/<species>/<release>_<assembly>` and places the indexed genome
//! FASTA inside it, the layout `vep --dir_cache` resolves. `--include-predictions`
//! and `--include-variations` add the matrices and the variation files to such
//! a build.
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
//!
//! ## Regulatory tree
//!
//! With `--regulatory-gff`, `--emars-gff` and `--motif-gff`, the Ensembl
//! regulatory annotation GFF files are baked into the cache's regulatory tree
//! (`regulatory/<chr>/<start>-<end>.json` shards and the matrix stubs of
//! `regulatory/binding_matrices.json`), the features a run with `--regulatory`
//! annotates against as a run with `--regulatory_gff` annotates against the
//! files themselves.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{bail, Context, Result};
use clap::Parser;
use tracing::{info, warn};

mod builder;
mod core_dumps;
mod download;
mod fasta;
mod fasta_index;
mod gff3;
mod gtf;
mod info;
mod install;
mod mapper;
pub mod predictions;
mod regulatory;
pub mod variation;

/// The flags of a `--auto` installation; every flag of a build named one input
/// at a time is refused beside any of them.
const INSTALL_FLAGS: [&str; 3] = ["auto", "dir_cache", "site"];

#[derive(Parser)]
#[command(name = "vep-cache-builder")]
#[command(
    about = "Build VEP JSON transcript cache from a GFF3 or GTF and a protein or genome FASTA, or install a release in one command with --auto"
)]
struct Args {
    /// Ensembl release number (e.g., 115)
    #[arg(long)]
    release: u32,

    /// Genome assembly (e.g., GRCh37, GRCh38)
    #[arg(long, required_unless_present = "auto", conflicts_with_all = INSTALL_FLAGS)]
    assembly: Option<String>,

    /// Species name (default: homo_sapiens)
    #[arg(long, default_value = "homo_sapiens", conflicts_with_all = INSTALL_FLAGS)]
    species: String,

    /// Output directory for JSON cache
    #[arg(long, required_unless_present = "auto", conflicts_with_all = INSTALL_FLAGS)]
    output_dir: Option<String>,

    /// Install the release in one command: fetch every input of the species and
    /// assembly from the Ensembl site (or `--site`), check each against the
    /// `CHECKSUMS` file of its directory, build the cache as
    /// `<dir_cache>/<species>/<release>_<assembly>` and place the indexed genome
    /// FASTA inside it, where `vep --dir_cache` finds both. The downloads land
    /// under `<dir_cache>/tmp/` and are removed once the cache is in place
    #[arg(long, num_args = 2, value_names = ["SPECIES", "ASSEMBLY"])]
    auto: Option<Vec<String>>,

    /// The cache root `--auto` installs under (default: `$HOME/.vep`), the
    /// directory a run names with `--dir_cache`
    #[arg(long = "dir_cache", visible_alias = "dir-cache", requires = "auto")]
    dir_cache: Option<PathBuf>,

    /// Where `--auto` fetches from: a URL or a local directory holding the
    /// Ensembl site's `pub/release-<N>/` tree (default: https://ftp.ensembl.org)
    #[arg(long, default_value = download::ENSEMBL_SITE, requires = "auto")]
    site: String,

    /// Local GFF3 as the transcript model, instead of the Ensembl download:
    /// Ensembl's own, or any GFF3, read as Ensembl VEP reads a `--gff` source
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    gff3: Option<String>,

    /// Local protein FASTA, the peptide of each transcript, instead of the
    /// Ensembl download. With `--genome-fasta`, without this flag and without
    /// `--include-predictions` the peptides are translated from the genome and
    /// no protein FASTA is read
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    fasta: Option<String>,

    /// Local GTF. Beside `--gff3`, the Ensembl GTF of the same release, the
    /// source of the `cds_start_NF` and `cds_end_NF` transcript attributes
    /// (FLAGS in the output) that the GFF3 does not carry. Without `--gff3`,
    /// the transcript model itself, read as Ensembl VEP reads a `--gtf` source
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    gtf: Option<String>,

    /// Directory of the release's Ensembl core and variation MySQL table dumps
    /// (`<table>.txt.gz`); a table missing from it is downloaded there. Supplies
    /// the transcript attributes, xrefs, protein features and gene phenotype
    /// counts the GFF3 lacks
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    core_dumps: Option<String>,

    /// The codon table of a sequence region, as `<region>=<table>`, where the
    /// core dumps state none: `MT=2` for a vertebrate's mitochondrion built
    /// without `--core-dumps`, `Mito=3` for yeast's, `mitochondrion_genome=5`
    /// for the fly's, `Pt=11,Mt=1` for the Arabidopsis plastid and
    /// mitochondrion. A region named like an organelle genome whose table
    /// neither the dumps nor this flag states refuses the build. Repeatable or
    /// comma-joined
    #[arg(
        long = "codon-table",
        value_name = "REGION=TABLE",
        value_delimiter = ','
    )]
    codon_table: Vec<String>,

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
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    variation_vcf_dir: Option<String>,

    /// Directory of the release's Ensembl variation MySQL table dumps
    /// (`source`, `seq_region`, `attrib_type`, `phenotype`, `phenotype_feature`
    /// and `phenotype_feature_attrib`, each `.txt.gz` or `.txt`): the ClinVar
    /// allele and somatic classification columns of the variation records,
    /// which no variation VCF carries. Read with `--include-variations`
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    variation_dumps: Option<String>,

    /// The chromosomes whose variation VCFs are read (default: the 25 human
    /// chromosomes under homo_sapiens, the GFF3's `##sequence-region` names
    /// for another species; a GTF declares none, so another species whose
    /// model is a GTF alone names them here)
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
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    genome_fasta: Option<String>,

    /// Ensembl regulatory feature GFF3 (`<species>.<assembly>.regulatory_features.v<N>.gff3.gz`
    /// under the site's `regulation/<species>/<assembly>/annotation/`), written as the
    /// cache's regulatory shards: promoters, enhancers, CTCF binding sites and open
    /// chromatin regions, each with its type as BIOTYPE, read by a run with --regulatory
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    regulatory_gff: Option<String>,

    /// Ensembl EMAR GFF (`EMARs.v<N>.gff.gz`): regulatory features of type EMAR, written
    /// beside those of --regulatory-gff
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    emars_gff: Option<String>,

    /// Ensembl motif feature GFF3 (`motif_features.v<N>.gff3.gz`): each motif names its
    /// binding matrix by identifier, and the cache carries a stub of the matrix (its span
    /// and transcription factors, no counts), so motif rows carry MOTIF_NAME, MOTIF_POS and
    /// TRANSCRIPTION_FACTORS and no HIGH_INF_POS or MOTIF_SCORE_CHANGE
    #[arg(long, conflicts_with_all = INSTALL_FLAGS)]
    motif_gff: Option<String>,

    /// Write each promoter of --regulatory-gff at its extended_start and extended_end, the
    /// wider region the Ensembl files record beside the core promoter
    #[arg(long, requires = "regulatory_gff")]
    extended_promoters: bool,
}

/// The MySQL server the SIFT and PolyPhen matrices are fetched from.
#[derive(Debug, Clone)]
pub struct MysqlServer {
    pub host: String,
    pub port: u16,
}

/// The variation inputs of a build: the VCFs, and the variation database's
/// ClinVar tables when read.
pub struct VariationBuild {
    pub inputs: download::VariationInputs,
    pub dumps: Option<PathBuf>,
}

/// Every input of one cache build, each resolved to a local path.
pub struct BuildInputs {
    pub species: String,
    pub assembly: String,
    pub release: u32,
    /// The transcript model's GFF3; `None` when the GTF is the model.
    pub gff3: Option<PathBuf>,
    /// The protein FASTA; `None` translates each peptide from the genome FASTA.
    pub protein_fasta: Option<PathBuf>,
    /// The GTF: the transcript model when `gff3` is `None`, else the source of
    /// the attributes and header the GFF3 lacks.
    pub gtf: Option<PathBuf>,
    /// A directory holding every core table dump the build reads.
    pub core_dumps: Option<PathBuf>,
    /// `--codon-table`: sequence region name to codon table, for regions the
    /// dumps state none for.
    pub codon_tables: HashMap<String, u8>,
    pub genome_fasta: Option<PathBuf>,
    pub predictions: Option<MysqlServer>,
    pub variations: Option<VariationBuild>,
    /// The regulatory GFF files baked into the cache's regulatory tree; `None`
    /// writes no tree.
    pub regulatory: Option<regulatory::RegulatoryInputs>,
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

    let predictions = args.include_predictions.then(|| MysqlServer {
        host: args.mysql_host.clone(),
        port: args.mysql_port,
    });
    let codon_tables = builder::codon_table_flag(&args.codon_table)?;

    if let Some(auto) = &args.auto {
        let [species, assembly] = auto.as_slice() else {
            bail!("--auto takes a species and an assembly");
        };
        let dir_cache = match &args.dir_cache {
            Some(dir) => dir.clone(),
            None => default_dir_cache()?,
        };
        let options = install::InstallOptions {
            species: species.clone(),
            assembly: assembly.clone(),
            release: args.release,
            dir_cache,
            site: install::Site::parse(&args.site)?,
            include_variations: args.include_variations,
            chromosomes: args.chromosomes.clone(),
            codon_tables,
            predictions,
        };
        install::install(&options).context("Failed to install the release")?;
        return Ok(());
    }

    let assembly = args
        .assembly
        .clone()
        .expect("clap requires --assembly without --auto");
    let output_dir = PathBuf::from(
        args.output_dir
            .as_deref()
            .expect("clap requires --output-dir without --auto"),
    );

    // The transcript model is the GFF3 when one is given or downloaded, else
    // the GTF; the GTF beside a GFF3 supplies tags alone.
    let gtf_is_the_model = args.gff3.is_none() && args.gtf.is_some();
    let gff3 = if gtf_is_the_model {
        None
    } else {
        Some(
            download::resolve_gff3(args.gff3.as_deref(), args.release, &assembly, &args.species)
                .context("Failed to resolve the GFF3")?,
        )
    };
    // The predictions are keyed by the peptide's MD5, so they need the protein
    // FASTA's residues; otherwise a genome stands in for it.
    let peptides_from_genome =
        args.fasta.is_none() && args.genome_fasta.is_some() && !args.include_predictions;
    let protein_fasta = if peptides_from_genome {
        None
    } else {
        Some(
            download::resolve_protein_fasta(
                args.fasta.as_deref(),
                args.release,
                &assembly,
                &args.species,
            )
            .context("Failed to resolve the protein FASTA")?,
        )
    };
    let gtf = args.gtf.as_deref().map(PathBuf::from);

    let core_dumps = match &args.core_dumps {
        Some(dir) => {
            let dir = PathBuf::from(dir);
            download::resolve_core_dumps(&dir, args.release, &assembly, &args.species)
                .context("Failed to resolve the core table dumps")?;
            Some(dir)
        }
        None => None,
    };

    let variations = if args.include_variations {
        let inputs = download::resolve_variation_inputs(
            &args.variation_vcf_dir,
            args.release,
            &assembly,
            &args.species,
            gff3.as_deref(),
            args.chromosomes.as_deref(),
        )
        .context("Failed to resolve variation inputs")?;
        Some(VariationBuild {
            inputs,
            dumps: args.variation_dumps.as_deref().map(PathBuf::from),
        })
    } else {
        None
    };

    let regulatory = regulatory::RegulatoryInputs {
        regulatory_gff: args.regulatory_gff.as_deref().map(PathBuf::from),
        emars_gff: args.emars_gff.as_deref().map(PathBuf::from),
        motif_gff: args.motif_gff.as_deref().map(PathBuf::from),
        extended_promoters: args.extended_promoters,
    };

    let inputs = BuildInputs {
        species: args.species.clone(),
        assembly,
        release: args.release,
        gff3,
        protein_fasta,
        gtf,
        core_dumps,
        codon_tables,
        genome_fasta: args.genome_fasta.as_deref().map(PathBuf::from),
        predictions,
        variations,
        regulatory: (!regulatory.is_empty()).then_some(regulatory),
    };
    build(&inputs, &output_dir)
}

/// `$HOME/.vep`, the cache root Ensembl VEP's installer defaults to.
fn default_dir_cache() -> Result<PathBuf> {
    std::env::var_os("HOME")
        .map(|home| PathBuf::from(home).join(".vep"))
        .context("--dir_cache is needed when HOME is not set")
}

/// Builds the cache from `inputs` into `output_dir`: the transcript shards, the
/// regulatory tree, the variation files, `chr_synonyms.txt`, `info.json` and,
/// last, the checksum manifest over every file the directory then holds. A
/// model from which no transcript is read, or whose every transcript a build
/// of the cache's content sets aside, fails the build before this function
/// writes anything, naming the model's file.
pub fn build(inputs: &BuildInputs, output_dir: &Path) -> Result<()> {
    let core_dumps = match &inputs.core_dumps {
        Some(dir) => {
            info!("Reading core table dumps: {}", dir.display());
            let dumps =
                core_dumps::CoreDumps::load(dir).context("Failed to read the core table dumps")?;
            Some(dumps)
        }
        None => None,
    };

    let protein_sequences = match &inputs.protein_fasta {
        Some(fasta_path) => {
            info!("Parsing protein FASTA: {}", fasta_path.display());
            let sequences =
                fasta::parse_protein_fasta(fasta_path).context("Failed to parse protein FASTA")?;
            if sequences.is_empty() {
                warn!(
                    "{} holds no record naming a transcript (`transcript:` token) with a \
                     sequence: the file supplies no peptide, and the build is still the cache's \
                     content",
                    fasta_path.display()
                );
            }
            info!("Loaded {} protein sequences from FASTA", sequences.len());
            Some(sequences)
        }
        None => {
            info!("No protein FASTA: each peptide is translated from the genome FASTA");
            None
        }
    };

    let reading = if builder::is_cache_content(protein_sequences.as_ref(), core_dumps.as_ref()) {
        gff3::Reading::EnsemblCache
    } else {
        gff3::Reading::GffSource
    };
    let (model_path, genes, gtf_model_header) = match &inputs.gff3 {
        Some(gff3_path) => {
            info!("Parsing GFF3: {}", gff3_path.display());
            let genes = gff3::parse_gff3(gff3_path, reading).context("Failed to parse GFF3")?;
            (gff3_path.as_path(), genes, None)
        }
        None => {
            let gtf_path = inputs.gtf.as_deref().expect("the GTF is the model");
            info!("Parsing GTF transcript model: {}", gtf_path.display());
            let model = gtf::parse_gtf_model(gtf_path, reading).context("Failed to parse GTF")?;
            (gtf_path, model.genes, Some(model.header))
        }
    };
    let transcript_count: usize = genes.values().map(|g| g.transcripts.len()).sum();
    // A cache is its transcripts: a model that yields none refuses the build
    // here, before this function writes a file or a directory.
    if transcript_count == 0 {
        bail!("No transcript was read from {}", model_path.display());
    }
    info!(
        "Parsed {} genes, {} transcripts",
        genes.len(),
        transcript_count
    );
    // A build of the cache's content that sets every transcript aside has no
    // shard to write, so it is refused here as a model with none is.
    if reading == gff3::Reading::EnsemblCache {
        let set_aside = builder::SetAside::of(&genes, core_dumps.as_ref());
        if set_aside.total() == transcript_count {
            bail!(
                "No transcript was kept from {}: {transcript_count} read, {} set aside as \
                 artifact and {} as readthrough, as Ensembl VEP's cache holds neither",
                model_path.display(),
                set_aside.artifact,
                set_aside.readthrough
            );
        }
    }
    // Settled before the predictions are fetched, so an unknown table refuses
    // the build before any MySQL round trip; no context line is added, so the
    // refusal prints as its own ERROR line.
    let codon_tables = builder::codon_tables(&genes, core_dumps.as_ref(), &inputs.codon_tables)?;

    let prediction_data = match &inputs.predictions {
        Some(server) => {
            info!("Fetching SIFT/PolyPhen predictions from Ensembl MySQL...");
            // Both callers resolve a protein FASTA for a build with predictions:
            // `run` through `peptides_from_genome`, `install::install` always.
            let peptides = protein_sequences
                .as_ref()
                .expect("a build with predictions reads the protein FASTA");
            let data = predictions::fetch_predictions(
                peptides,
                inputs.release,
                &inputs.assembly,
                &server.host,
                server.port,
            )
            .context("Failed to fetch SIFT/PolyPhen predictions")?;
            info!("Fetched predictions for {} transcripts", data.len());
            Some(data)
        }
        None => None,
    };

    let genome_fasta = match &inputs.genome_fasta {
        Some(path) => {
            info!("Loading genome FASTA: {}", path.display());
            let fasta = vep_fasta::IndexedFasta::from_path(path)
                .context("Failed to load genome FASTA (needs .fa + .fai index)")?;
            Some(fasta)
        }
        None => None,
    };

    let gtf_tags = match &inputs.gtf {
        Some(path) if inputs.gff3.is_some() => {
            info!("Parsing GTF transcript tags: {}", path.display());
            let tags = gtf::parse_gtf_tags(path).context("Failed to parse GTF")?;
            info!("Read tags for {} transcripts from GTF", tags.tags.len());
            Some(tags)
        }
        _ => None,
    };
    let gtf_header = gtf_model_header.or_else(|| gtf_tags.as_ref().map(|tags| tags.header.clone()));

    info!("Building JSON cache in: {}", output_dir.display());
    let stats = builder::build_cache(
        &genes,
        model_path,
        output_dir,
        &builder::BuildSources {
            protein_sequences: protein_sequences.as_ref(),
            prediction_data: prediction_data.as_ref(),
            genome_fasta: genome_fasta.as_ref(),
            gtf_tags: gtf_tags.as_ref(),
            core_dumps: core_dumps.as_ref(),
        },
        &codon_tables,
    )
    .context("Failed to build cache")?;
    info!(
        "Cache built: {} transcripts across {} chromosomes, {} region files",
        stats.transcript_count, stats.chromosome_count, stats.region_file_count
    );

    if let Some(regulatory) = &inputs.regulatory {
        regulatory::build_regulatory(regulatory, output_dir)
            .context("Failed to build the regulatory tree")?;
    }

    let variation_cols = match &inputs.variations {
        Some(variations) => Some(write_variation_files(variations, output_dir)?),
        None => None,
    };

    // Source versions the output headers print: the core `meta` table's, as
    // Ensembl VEP's own cache info carries them, else the GTF header's.
    let mut source_versions = core_dumps
        .as_ref()
        .map(|d| d.source_versions())
        .unwrap_or_default();
    if let Some(header) = &gtf_header {
        if let Some(v) = header.get("genome-build") {
            source_versions
                .entry("assembly".to_string())
                .or_insert_with(|| v.clone());
        }
        if let Some(v) = header.get("genebuild-last-updated") {
            source_versions
                .entry("genebuild".to_string())
                .or_insert_with(|| v.clone());
        }
    }
    if let Some(dir) = &inputs.core_dumps {
        info::write_chr_synonyms(output_dir, dir).context("Failed to write chr_synonyms.txt")?;
    }
    info::CacheInfoJson::new(
        &inputs.species,
        &inputs.assembly,
        inputs.release,
        &stats,
        variation_cols,
        source_versions,
        core_dumps.as_ref().map(|dumps| dumps.karyotype_ranks()),
    )
    .with_regulatory(inputs.regulatory.is_some())
    .write(output_dir)
    .context("Failed to write info.json")?;

    let listed = vep_core::cache_manifest::CacheManifest::write(output_dir)
        .context("Failed to write CHECKSUMS.sha256")?;
    info!("Checksum manifest written for {} files", listed);

    Ok(())
}

/// Writes the variation file of each chromosome in `variations` into
/// `output_dir` (`variations/<chr>/all_vars.gz` with its `.csi`), the ClinVar
/// columns filled from the variation table dumps when given, and returns the
/// columns `info.json` lists as `variation_cols`.
fn write_variation_files(variations: &VariationBuild, output_dir: &Path) -> Result<Vec<String>> {
    let variation_inputs = &variations.inputs;
    info!(
        "Resolved {} per-chromosome variation VCF files",
        variation_inputs.per_chromosome.len()
    );

    let clinvar = match &variations.dumps {
        Some(dir) => {
            info!("Reading the ClinVar phenotype features: {}", dir.display());
            let chromosomes: Vec<String> = variation_inputs
                .per_chromosome
                .iter()
                .map(|(chrom, _)| chrom.clone())
                .collect();
            let attribs = variation::ClinVarAttribs::load(dir, &chromosomes)
                .context("Failed to read the variation table dumps")?;
            Some(attribs)
        }
        None => None,
    };

    let columns = variation::variation_columns(variation_inputs.genome_wide.frequencies.is_some());
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
            output_dir,
        )
        .with_context(|| format!("Failed to write the variation file of chr{}", chrom))?;
    }
    info!(
        "Variation files written: {} records across {} chromosomes",
        total_variations,
        variation_inputs.per_chromosome.len()
    );
    Ok(columns)
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    fn parse(argv: &[&str]) -> Result<Args, clap::Error> {
        Args::try_parse_from(std::iter::once("vep-cache-builder").chain(argv.iter().copied()))
    }

    /// `--auto` takes the species and the assembly as its two values and
    /// releases `--assembly` and `--output-dir` from being required; without it
    /// both stay required. `--dir_cache` (or `--dir-cache`) and `--site` go
    /// with `--auto` alone, and the local-file flags, `--species` and
    /// `--output-dir` are refused beside it; the prediction and variation
    /// switches combine with it.
    #[test]
    fn auto_flags_parse_as_the_installer_options() {
        Args::command().debug_assert();
        let args = parse(&[
            "--auto",
            "homo_sapiens",
            "GRCh38",
            "--release",
            "116",
            "--dir-cache",
            "/tmp/vep",
            "--site",
            "/site",
            "--include-predictions",
            "--include-variations",
            "--chromosomes",
            "21,MT",
            "--codon-table",
            "Pt=11,Mt=1",
            "--codon-table",
            "MT=2",
        ])
        .unwrap();
        assert_eq!(
            args.auto.as_deref(),
            Some(&["homo_sapiens".to_string(), "GRCh38".to_string()][..])
        );
        assert_eq!(args.dir_cache.as_deref(), Some(Path::new("/tmp/vep")));
        assert_eq!(args.site, "/site");
        assert!(args.include_predictions && args.include_variations);
        assert_eq!(
            args.chromosomes,
            Some(vec!["21".to_string(), "MT".to_string()])
        );
        assert_eq!(args.codon_table, ["Pt=11", "Mt=1", "MT=2"]);
        let plain = parse(&["--auto", "homo_sapiens", "GRCh38", "--release", "116"]).unwrap();
        assert_eq!(plain.site, download::ENSEMBL_SITE);
        assert!(plain.dir_cache.is_none() && plain.assembly.is_none());

        for refused in [
            vec!["--auto", "homo_sapiens", "--release", "116"],
            vec!["--release", "116", "--output-dir", "out"],
            vec!["--release", "116", "--assembly", "GRCh38"],
            vec![
                "--dir_cache",
                "/tmp/vep",
                "--release",
                "116",
                "--assembly",
                "GRCh38",
                "--output-dir",
                "out",
            ],
            vec![
                "--site",
                "/m",
                "--release",
                "116",
                "--assembly",
                "GRCh38",
                "--output-dir",
                "out",
            ],
            vec![
                "--auto",
                "homo_sapiens",
                "GRCh38",
                "--release",
                "116",
                "--assembly",
                "GRCh38",
            ],
            vec![
                "--auto",
                "homo_sapiens",
                "GRCh38",
                "--release",
                "116",
                "--species",
                "homo_sapiens",
            ],
            vec![
                "--auto",
                "homo_sapiens",
                "GRCh38",
                "--release",
                "116",
                "--output-dir",
                "out",
            ],
            vec![
                "--auto",
                "homo_sapiens",
                "GRCh38",
                "--release",
                "116",
                "--gff3",
                "a.gff3",
            ],
            vec![
                "--auto",
                "homo_sapiens",
                "GRCh38",
                "--release",
                "116",
                "--genome-fasta",
                "a.fa",
            ],
            vec![
                "--auto",
                "homo_sapiens",
                "GRCh38",
                "--release",
                "116",
                "--core-dumps",
                "d",
            ],
        ] {
            assert!(parse(&refused).is_err(), "{refused:?} parses");
        }
        let explicit = parse(&[
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--output-dir",
            "out",
        ])
        .unwrap();
        assert_eq!(explicit.assembly.as_deref(), Some("GRCh38"));
        assert_eq!(explicit.species, "homo_sapiens");
    }
}
