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
//! FASTA inside it, the layout `vep --dir_cache` resolves. The installation
//! brings the variation files, the SIFT and PolyPhen matrices and the
//! regulatory tree wherever the release publishes a file of them for the
//! genome, the matrices from the variation database's table dumps, and ends
//! on one line a part it left out. `--include-predictions` fetches the
//! matrices from Ensembl's MySQL server instead. A species of an Ensembl
//! Genomes division is found on the Ensembl Genomes site under that site's
//! own release number.
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
//! With `--include-variations`, Ensembl's variation VCFs (for `homo_sapiens`
//! the per-chromosome germline files, the genome-wide somatic, phenotype- and
//! clinically-associated files and the 1000 Genomes frequency file; for any
//! other species its one genome-wide file, `<species>.vcf.gz`) are read into
//! the variation files Ensembl VEP's own cache carries:
//! `variations/<chr>/all_vars.gz` with its `.csi`, the columns listed by
//! `info.json`'s `variation_cols`. An installation also reads, where the
//! release publishes them, the tables of the variation database that hold a
//! variant's citations, its synonyms and its failed mark, and the GVF file
//! of the failed variants: it fills the `pubmed`, `var_synonyms` and
//! `failed` cells of the records and adds a record for each failed variant
//! the VCFs leave out. The tables are read whole whatever `--chromosomes`
//! names.
//!
//! ## Regulatory tree
//!
//! With `--regulatory-gff`, `--emars-gff` and `--motif-gff`, the Ensembl
//! regulatory annotation GFF files are baked into the cache's regulatory tree
//! (`regulatory/<chr>/<start>-<end>.json` shards and the matrix stubs of
//! `regulatory/binding_matrices.json`), the features a run with `--regulatory`
//! annotates against as a run with `--regulatory_gff` annotates against the
//! files themselves. With a genome FASTA each motif carries its sequence. An
//! installation also reads, where the release publishes them, the activity
//! file beside the GFF files and the matrix tables of the funcgen database,
//! and bakes each regulatory feature's activity by cell type and each
//! matrix's counts.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{bail, Context, Result};
use clap::Parser;
use tracing::{info, warn};
use vep_core::stderr;

use vep_cache_builder::{
    builder, core_dumps, count, dump_rows, fasta, fasta_index, gff3, gtf, info, regulatory,
    variant_cells, variation,
};

mod download;
mod install;
mod left_out;
pub mod predictions;

use count::counted;

/// The flags of a `--auto` installation; every flag of a build named one input
/// at a time is refused beside any of them.
const INSTALL_FLAGS: [&str; 3] = ["auto", "dir_cache", "site"];

#[derive(Parser)]
#[command(name = "vep-cache-builder", version = env!("CARGO_PKG_VERSION"))]
#[command(
    about = "Build VEP JSON transcript cache from a GFF3 or GTF and a protein or genome FASTA, or install a release in one command with --auto"
)]
struct Args {
    /// Ensembl release number (e.g., 115); for a species of Ensembl Genomes,
    /// the Ensembl Genomes release (63 beside Ensembl 116)
    #[arg(long)]
    release: u32,

    /// Genome assembly (e.g., GRCh37, GRCh38)
    #[arg(long, required_unless_present = "auto", conflicts_with_all = INSTALL_FLAGS)]
    assembly: Option<String>,

    /// Species name, the Ensembl directory name: names the downloads and becomes
    /// `species` in `info.json` (default: homo_sapiens, as `vep --species` has)
    #[arg(long, default_value = "homo_sapiens", conflicts_with_all = INSTALL_FLAGS)]
    species: String,

    /// Output directory for JSON cache. A file the build downloads lands under
    /// its `tmp/` and is removed once the cache is built
    #[arg(long, required_unless_present = "auto", conflicts_with_all = INSTALL_FLAGS)]
    output_dir: Option<String>,

    /// Install the release in one command: fetch every input of the species and
    /// assembly from the Ensembl site (or `--site`), check each against the
    /// `CHECKSUMS` file of its directory, build the cache as
    /// `<dir_cache>/<species>/<release>_<assembly>` and place the indexed genome
    /// FASTA inside it, where `vep --dir_cache` finds both. The downloads land
    /// under `<dir_cache>/tmp/` and are removed once the cache is in place. The
    /// cache carries the known variants with the 1000 Genomes frequencies, the
    /// SIFT and PolyPhen matrices and the regulatory and motif features
    /// wherever the release publishes a file of them for the genome, and the
    /// installation ends on one line for each part it left out, with the
    /// reason. A
    /// species of the plants, metazoa, fungi or protists division is fetched
    /// from the Ensembl Genomes site; a genome of a collection (every bacterium)
    /// and an assembly of the rapid release are refused before any download
    #[arg(long, num_args = 2, value_names = ["SPECIES", "ASSEMBLY"])]
    auto: Option<Vec<String>>,

    /// The cache root `--auto` installs under (default: `$HOME/.vep`), the
    /// directory a run names with `--dir_cache`
    #[arg(long = "dir_cache", visible_alias = "dir-cache", requires = "auto")]
    dir_cache: Option<PathBuf>,

    /// Where `--auto` fetches from: a URL or a local directory holding the
    /// `pub/release-<N>/` tree of the Ensembl site or of the Ensembl Genomes
    /// site (default: <https://ftp.ensembl.org>, and
    /// <https://ftp.ensemblgenomes.ebi.ac.uk> for a species the first does not
    /// list)
    #[arg(long, requires = "auto")]
    site: Option<String>,

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

    /// Include variation data in the cache, from Ensembl's variation VCFs,
    /// downloaded or read from `--variation-vcf-dir`: for homo_sapiens the
    /// per-chromosome files with the genome-wide somatic, phenotype-associated,
    /// clinically-associated and 1000 Genomes frequency VCFs; for any other
    /// species its one genome-wide file, `<species>.vcf.gz`. An `--auto`
    /// installation includes it without the flag
    #[arg(long)]
    include_variations: bool,

    /// Skip download, use local variation VCF directory: `<species>-chr<N>.vcf.gz`
    /// per chromosome, or, in a directory holding no per-chromosome file, the
    /// one genome-wide `<species>.vcf.gz`; and, under their Ensembl names, the
    /// further genome-wide files (`1000GENOMES-phase_3.vcf.gz`,
    /// `<species>_somatic.vcf.gz`, `<species>_phenotype_associated.vcf.gz`,
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

    /// The chromosomes whose variation records are read (default: the 25 human
    /// chromosomes under homo_sapiens; for another species every sequence its
    /// genome-wide variation VCF holds records on). The flag names the
    /// per-chromosome files of homo_sapiens and narrows the read of another
    /// species' one genome-wide file, with `--auto` as in a build from named
    /// inputs
    #[arg(long, value_delimiter = ',')]
    chromosomes: Option<Vec<String>>,

    /// Fetch SIFT/PolyPhen prediction matrices from Ensembl MySQL, which serves
    /// homo_sapiens on GRCh37 and GRCh38. An `--auto` installation without
    /// the flag reads the matrices of any species from the release's variation
    /// table dumps, and with it from the server instead: the route to the
    /// GRCh37 matrices, whose key table the GRCh37 tree does not publish
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

/// Where a build reads the SIFT and PolyPhen matrices from.
#[derive(Debug, Clone)]
pub enum PredictionSource {
    /// Ensembl's MySQL server, under `--include-predictions`.
    Server(MysqlServer),
    /// A directory holding the dumps of [`predictions::PREDICTION_TABLES`],
    /// the files an `--auto` installation fetches from the release.
    Dumps(PathBuf),
}

/// The variation inputs of a build: the VCFs, the variation database's
/// ClinVar tables when read and, in an installation, the files the `pubmed`,
/// `var_synonyms` and `failed` cells are read from.
pub struct VariationBuild {
    pub inputs: download::VariationInputs,
    pub dumps: Option<PathBuf>,
    /// The tables and the failed-variants GVF an installation found for the
    /// three cells; `None` in a build named one input at a time, which reads
    /// none.
    pub cells: Option<variant_cells::VariantCellInputs>,
    /// Where an installation reads the known variants of a sequence its
    /// variation VCFs hold no record on; `None` in a build named one input
    /// at a time and where the release publishes no such table.
    pub features: Option<SequenceFeatures>,
}

/// The `variation_feature` table of an installation's variation database,
/// read for the sequences of the cache that no variation VCF holds a record
/// on. The rule is one of the layout of the variation files: it applies
/// where every germline record is read from one genome-wide file
/// ([`variation::GermlineVcf::Genome`]), which can hold no record on a
/// sequence the database holds known variants on. Where the files are one a
/// chromosome, a sequence without a file is one the release publishes no
/// known variant of, and nothing is read for it.
pub struct SequenceFeatures {
    /// The sequences `--chromosomes` names: with them, the only ones in
    /// question.
    pub named: Option<Vec<String>>,
    /// Places the tables and names them. Called once, and only when a
    /// sequence is in question.
    pub fetch: Box<dyn Fn() -> Result<variation::FeatureTables>>,
}

/// What a build read of the `variation_feature` table.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SequenceFeatureFacts {
    /// The sequences of the cache that no variation VCF holds a record on,
    /// among those `--chromosomes` names when it names any. With none the
    /// table is neither fetched nor read.
    pub in_question: Vec<String>,
    /// Each of those sequences whose variation file holds a record, with
    /// the number of its records.
    pub brought: Vec<(String, usize)>,
}

/// What a build wrote that its caller reads back.
pub struct Built {
    /// What the build read for the `pubmed`, `var_synonyms` and `failed`
    /// cells of its known variants and what it wrote of them; `None` for a
    /// build given no file of them.
    pub variant_cells: Option<variant_cells::VariantCellFacts>,
    /// What the build read of the `variation_feature` table; `None` for a
    /// build given none and for one whose variation files are one a
    /// chromosome.
    pub sequence_features: Option<SequenceFeatureFacts>,
    /// What the build wrote of the regulatory tree; `None` for a build given
    /// no regulatory file.
    pub regulatory: Option<regulatory::RegulatoryStats>,
    /// What the cache leaves out (`left_out::of`): the names `info.json`
    /// lists and the lines the caller logs last.
    pub left_out: left_out::LeftOut,
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
    pub predictions: Option<PredictionSource>,
    pub variations: Option<VariationBuild>,
    /// The regulatory GFF files baked into the cache's regulatory tree, with
    /// the activity file and the funcgen table dumps an installation reads
    /// beside them; `None` writes no tree.
    pub regulatory: Option<regulatory::RegulatoryInputs>,
    /// The work directory under the output directory, which holds the
    /// downloads of a build named one input at a time: [`build`] removes it
    /// before the checksum manifest lists the directory's files. `None` when
    /// the caller's work directory lies elsewhere.
    pub work: Option<download::WorkDir>,
}

/// A failed build prints one line, `ERROR: [cache-builder] <message>`, on
/// standard error and exits 1. A line standard error refuses, a log line or
/// that one, is dropped and changes neither the build nor its exit status.
fn main() -> ExitCode {
    stderr::drop_refused_writes();
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            stderr::write(format_args!("ERROR: [cache-builder] {error:#}\n"));
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<()> {
    // The subscriber reports a write its writer fails with `eprintln!`, which
    // panics on a stream that refuses it; `Stream` fails no write.
    tracing_subscriber::fmt()
        .with_ansi(log_lines_carry_colour(
            std::io::IsTerminal::is_terminal(&std::io::stderr()),
            std::env::var_os("NO_COLOR").as_deref(),
        ))
        .with_max_level(tracing::Level::INFO)
        .with_writer(|| stderr::Stream)
        .init();

    let args = Args::parse();

    // An installation reads the variation data by itself, and `--chromosomes`
    // with it.
    if !args.include_variations && args.auto.is_none() {
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

    if args.include_predictions {
        // Refused here, before any input is resolved or downloaded.
        let (species, assembly) = match args.auto.as_deref() {
            Some([species, assembly]) => (species.as_str(), assembly.as_str()),
            _ => (
                args.species.as_str(),
                args.assembly.as_deref().unwrap_or_default(),
            ),
        };
        predictions::prediction_database(species, args.release, assembly)?;
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
            sites: match &args.site {
                Some(site) => vec![install::Site::parse(site)?],
                None => install::Site::defaults(),
            },
            chromosomes: args.chromosomes.clone(),
            codon_tables,
            prediction_server: predictions,
        };
        let work = install::work_dir(&options);
        return install::install(&options)
            .map(|_| ())
            .context("Failed to install the release")
            .map_err(|error| work.kept(error));
    }

    build_named_inputs(&args, predictions, codon_tables, &install::Site::defaults())
}

/// A build named one input at a time: the cache of the inputs the flags name
/// in `--output-dir`. A file the flags leave out is fetched into the work
/// directory under the output directory, which [`build`] removes once the
/// cache is built; a build that fails with a download in it keeps it and
/// names it in its error. The GFF3 and the protein FASTA come from the first
/// of `sites` that publishes the genome, a variation VCF from the first of
/// `sites`.
fn build_named_inputs(
    args: &Args,
    predictions: Option<MysqlServer>,
    codon_tables: HashMap<String, u8>,
    sites: &[install::Site],
) -> Result<()> {
    let assembly = args
        .assembly
        .clone()
        .expect("clap requires --assembly without --auto");
    let output_dir = PathBuf::from(
        args.output_dir
            .as_deref()
            .expect("clap requires --output-dir without --auto"),
    );
    let site = sites.first().context("no site to read the release from")?;
    let published = install::Published::new(sites.to_vec(), &args.species, &assembly, args.release);
    let work =
        download::WorkDir::of_output_dir(&output_dir, &args.species, args.release, &assembly);
    let from = download::Downloads {
        published: &published,
        site,
        work: &work,
    };
    let left_out = named_inputs(args, &assembly, predictions, codon_tables, &from)
        .and_then(|inputs| build(&inputs, &output_dir))
        .map(|built| built.left_out)
        .map_err(|error| work.kept(error))?;
    // The order of an installation's closing lines: the parts, the sequences
    // without a known variant, then the cells.
    let sequences = install::sequences_warning_of_a_build(&args.species, &assembly, &output_dir)?;
    for line in left_out
        .lines
        .iter()
        .chain(&sequences)
        .chain(&left_out.cell_lines)
    {
        warn!("{line}");
    }
    Ok(())
}

/// The inputs the flags name one at a time, each resolved to a local path: a
/// file the flags leave out is fetched through `from`.
fn named_inputs(
    args: &Args,
    assembly: &str,
    predictions: Option<MysqlServer>,
    codon_tables: HashMap<String, u8>,
    from: &download::Downloads<'_>,
) -> Result<BuildInputs> {
    // The transcript model is the GFF3 when one is given or downloaded, else
    // the GTF; the GTF beside a GFF3 supplies tags alone.
    let gtf_is_the_model = args.gff3.is_none() && args.gtf.is_some();
    let gff3 = if gtf_is_the_model {
        None
    } else {
        Some(
            download::resolve_gff3(
                args.gff3.as_deref(),
                args.release,
                assembly,
                &args.species,
                from,
            )
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
                assembly,
                &args.species,
                from,
            )
            .context("Failed to resolve the protein FASTA")?,
        )
    };
    let gtf = args.gtf.as_deref().map(PathBuf::from);

    let core_dumps = match &args.core_dumps {
        Some(dir) => {
            let dir = PathBuf::from(dir);
            download::resolve_core_dumps(&dir, &|| from.published.database_urls())
                .context("Failed to resolve the core table dumps")?;
            Some(dir)
        }
        None => None,
    };

    let variations = if args.include_variations {
        let inputs = download::resolve_variation_inputs(
            &args.variation_vcf_dir,
            args.release,
            assembly,
            &args.species,
            args.chromosomes.as_deref(),
            from,
        )
        .context("Failed to resolve variation inputs")?;
        Some(VariationBuild {
            inputs,
            dumps: args.variation_dumps.as_deref().map(PathBuf::from),
            cells: None,
            features: None,
        })
    } else {
        None
    };

    let regulatory = regulatory::RegulatoryInputs {
        regulatory_gff: args.regulatory_gff.as_deref().map(PathBuf::from),
        emars_gff: args.emars_gff.as_deref().map(PathBuf::from),
        motif_gff: args.motif_gff.as_deref().map(PathBuf::from),
        extended_promoters: args.extended_promoters,
        activity: None,
        funcgen_dumps: None,
    };

    Ok(BuildInputs {
        species: args.species.clone(),
        assembly: assembly.to_string(),
        release: args.release,
        gff3,
        protein_fasta,
        gtf,
        core_dumps,
        codon_tables,
        genome_fasta: args.genome_fasta.as_deref().map(PathBuf::from),
        predictions: predictions.map(PredictionSource::Server),
        variations,
        regulatory: (!regulatory.is_empty()).then_some(regulatory),
        work: Some(from.work.clone()),
    })
}

/// `$HOME/.vep`, the cache root Ensembl VEP's installer defaults to.
fn default_dir_cache() -> Result<PathBuf> {
    std::env::var_os("HOME")
        .map(|home| PathBuf::from(home).join(".vep"))
        .context("--dir_cache is needed when HOME is not set")
}

/// Builds the cache from `inputs` into `output_dir`: the transcript shards, the
/// regulatory tree, the variation files, `chr_synonyms.txt`, `info.json` and,
/// last, once the work directory of `inputs` is removed, the checksum manifest
/// over every file the directory then holds. A
/// model from which no transcript is read, or whose every transcript a build
/// of the cache's content sets aside, fails the build before this function
/// writes anything, naming the model's file; so does a GRCh37 model whose
/// mitochondrial coding transcripts do not read on the genome FASTA record
/// they are cut from
/// (`builder::grch37_mitochondrial_regions_off_their_record`). The build writes into
/// `output_dir` itself, beside any file the directory holds, and the
/// manifest lists those files too. Before its first file the manifest is one
/// line that lists no file (`CacheManifest::mark_unfinished`), so a build
/// that fails or is killed after it leaves a directory `vep` refuses by name
/// and a build of the same directory completes. Returns what the build read
/// and wrote for the cells of its known variants that the variation
/// database's tables carry, what it wrote of the regulatory tree, and what
/// the cache leaves out (`left_out::of`), whose lines the caller writes as
/// the last warnings of its log; `info.json` lists the same parts by name
/// under `left_out`, and none of a cell whose source the build read.
pub fn build(inputs: &BuildInputs, output_dir: &Path) -> Result<Built> {
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
            info!(
                "Loaded {} from FASTA",
                counted(sequences.len(), "protein sequence")
            );
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
    // Opened before the model is read, whose sequences are named by what the
    // FASTA holds; a FASTA that does not open is reported where the build
    // first reads it, after the model's own refusals.
    let genome_fasta = inputs.genome_fasta.as_ref().map(|path| {
        info!("Loading genome FASTA: {}", path.display());
        vep_fasta::IndexedFasta::from_path(path)
    });
    let sequence_names = match &genome_fasta {
        Some(Ok(fasta)) => gff3::SequenceNames::of_fasta(fasta, Some(&inputs.assembly)),
        _ => gff3::SequenceNames::default(),
    };
    let (model_path, genes, gtf_model_header) = match &inputs.gff3 {
        Some(gff3_path) => {
            info!("Parsing GFF3: {}", gff3_path.display());
            let genes = gff3::parse_gff3(gff3_path, reading, &sequence_names)
                .context("Failed to parse GFF3")?;
            (gff3_path.as_path(), genes, None)
        }
        None => {
            let gtf_path = inputs.gtf.as_deref().expect("the GTF is the model");
            info!("Parsing GTF transcript model: {}", gtf_path.display());
            let model = gtf::parse_gtf_model(gtf_path, reading, &sequence_names)
                .context("Failed to parse GTF")?;
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
        "Parsed {}, {}",
        counted(genes.len(), "gene"),
        counted(transcript_count, "transcript")
    );
    // A region kept apart from GRCh37's `MT` is built under the name the
    // model writes; the build says so once a region, in name order. A
    // region is the sequence a transcript's own line names.
    let kept_apart: std::collections::BTreeSet<&str> = genes
        .values()
        .flat_map(|gene| &gene.transcripts)
        .map(|transcript| transcript.chr.as_str())
        .filter(|region| sequence_names.keeps_apart_from_mt(region))
        .collect();
    for region in kept_apart {
        warn!(
            "Keeping {region} of {} under its written name: {}, and the genome FASTA does \
             neither",
            model_path.display(),
            vep_cache_builder::gene_model::GRCH37_MITOCHONDRIAL_NAMING
        );
    }
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
    let codon_tables = builder::codon_tables(
        &genes,
        core_dumps.as_ref(),
        &inputs.codon_tables,
        &sequence_names,
    )?;

    // Both callers resolve a protein FASTA for a build with predictions:
    // `run` through `peptides_from_genome`, `install::install` always.
    let peptides_of_the_fasta = || {
        protein_sequences
            .as_ref()
            .expect("a build with predictions reads the protein FASTA")
    };
    let prediction_data = match &inputs.predictions {
        Some(PredictionSource::Server(server)) => {
            info!("Fetching SIFT/PolyPhen predictions from Ensembl MySQL...");
            let data = predictions::fetch_predictions(
                peptides_of_the_fasta(),
                &inputs.species,
                inputs.release,
                &inputs.assembly,
                &server.host,
                server.port,
            )
            .context("Failed to fetch SIFT/PolyPhen predictions")?;
            info!(
                "Fetched predictions for {}",
                counted(data.len(), "transcript")
            );
            Some(data)
        }
        Some(PredictionSource::Dumps(dir)) => {
            info!(
                "Reading SIFT/PolyPhen predictions from the variation table dumps: {}",
                dir.display()
            );
            let data = predictions::read_prediction_dumps(peptides_of_the_fasta(), dir)
                .context("Failed to read the SIFT/PolyPhen predictions")?;
            Some(data)
        }
        None => None,
    };

    let genome_fasta = genome_fasta
        .transpose()
        .context("Failed to load genome FASTA (needs .fa + .fai index)")?;

    let gtf_tags = match &inputs.gtf {
        Some(path) if inputs.gff3.is_some() => {
            info!("Parsing GTF transcript tags: {}", path.display());
            let tags = gtf::parse_gtf_tags(path).context("Failed to parse GTF")?;
            info!(
                "Read tags for {} from GTF",
                counted(tags.names(), "transcript")
            );
            Some(tags)
        }
        _ => None,
    };
    let gtf_header = gtf_model_header.or_else(|| gtf_tags.as_ref().map(|tags| tags.header.clone()));

    let sources = builder::BuildSources {
        protein_sequences: protein_sequences.as_ref(),
        prediction_data: prediction_data.as_ref(),
        genome_fasta: genome_fasta.as_ref(),
        gtf_tags: gtf_tags.as_ref(),
        core_dumps: core_dumps.as_ref(),
    };
    // A GRCh37 mitochondrial region whose coding transcripts do not read on
    // the FASTA record they are cut from refuses the build here, before a
    // file is written: its sequences would be another sequence's than its
    // models count on. The test reads the sources the records are built
    // from, the GTF's tags among them.
    let off_record = builder::grch37_mitochondrial_regions_off_their_record(
        &genes,
        &sources,
        &codon_tables,
        Some(&inputs.assembly),
    );
    if let Some((region, found)) = off_record.iter().next() {
        bail!(
            "sequence region {region} of {} does not fit the genome FASTA: {}. On GRCh37 two \
             mitochondrial sequences are in use, the assembly's of 16,569 bases and the \
             16,571-base one of the UCSC hg19 naming: build from a genome FASTA that holds the \
             sequence this model counts on, or from a model without the region's lines",
            model_path.display(),
            found.finding()
        );
    }

    info!("Building JSON cache in: {}", output_dir.display());
    let stats = builder::build_cache(&genes, model_path, output_dir, &sources, &codon_tables)
        .context("Failed to build cache")?;
    info!(
        "Cache built: {} across {}, {}",
        counted(stats.transcript_count, "transcript"),
        counted(stats.chromosome_count, "chromosome"),
        counted(stats.region_file_count, "region file")
    );
    // The transcript phase's inputs are released before the regulatory and
    // variation phases; of the core dumps only what `info.json` prints is kept.
    let core_source_versions = core_dumps
        .as_ref()
        .map(|d| d.source_versions())
        .unwrap_or_default();
    let karyotype_ranks = core_dumps
        .as_ref()
        .map(|dumps| dumps.karyotype_ranks().clone());
    let core_dumps_read = core_dumps.is_some();
    drop(core_dumps);
    drop(genes);
    drop(protein_sequences);
    drop(prediction_data);
    drop(gtf_tags);

    let regulatory_tree = match &inputs.regulatory {
        Some(regulatory) => Some(
            regulatory::build_regulatory(
                regulatory,
                output_dir,
                &sequence_names,
                genome_fasta.as_ref(),
            )
            .context("Failed to build the regulatory tree")?,
        ),
        None => None,
    };

    let (variation_cols, variant_cells, sequence_features) = match &inputs.variations {
        Some(variations) => {
            let files = write_variation_files(variations, output_dir, &stats.chromosomes)?;
            (Some(files.columns), files.cells, files.features)
        }
        None => (None, None, None),
    };

    // Source versions the output headers print: the core `meta` table's, as
    // Ensembl VEP's own cache info carries them, else the GTF header's.
    let mut source_versions = core_source_versions;
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
    // The version Ensembl VEP's cache info carries as `source_regbuild`
    // (`Pipeline/DumpVEP/Dumper/Regulation.pm:90-91`).
    let regbuild = regulatory_tree
        .as_ref()
        .and_then(|tree| tree.regbuild.clone());
    if let Some(version) = regbuild {
        source_versions.insert("regbuild".to_string(), version);
    }
    if let Some(dir) = &inputs.core_dumps {
        info::write_chr_synonyms(output_dir, dir).context("Failed to write chr_synonyms.txt")?;
    }
    let left = left_out::of(&left_out::Build {
        species: &inputs.species,
        assembly: &inputs.assembly,
        core_dumps: core_dumps_read,
        reading,
        unnamed_translations: stats.unnamed_translations,
        variation: variation_cols.is_some(),
        panel_frequency: variation_cols
            .as_ref()
            .is_some_and(|columns| columns.iter().any(|column| column == "AF")),
        regulatory_features: inputs
            .regulatory
            .as_ref()
            .is_some_and(|files| files.regulatory_gff.is_some() || files.emars_gff.is_some()),
        motif_features: inputs
            .regulatory
            .as_ref()
            .is_some_and(|files| files.motif_gff.is_some()),
        pubmed: variant_cells
            .as_ref()
            .is_some_and(|cells| cells.pubmed_tables_read),
        var_synonyms: variant_cells
            .as_ref()
            .is_some_and(|cells| cells.synonym_tables_read),
        failed_variants: variant_cells
            .as_ref()
            .is_some_and(|cells| cells.failed_gvf_read || cells.failed_variation_joined),
        cell_type_activity: regulatory_tree
            .as_ref()
            .is_some_and(|tree| tree.activity_on_every_feature),
        motif_scores: regulatory_tree
            .as_ref()
            .is_some_and(|tree| tree.counts_on_every_matrix && tree.sequence_on_every_motif),
    });
    info::CacheInfoJson::new(
        &inputs.species,
        &inputs.assembly,
        inputs.release,
        &stats,
        variation_cols,
        source_versions,
        karyotype_ranks.as_ref(),
    )
    .with_regulatory(inputs.regulatory.is_some())
    .with_cell_types(
        regulatory_tree
            .as_ref()
            .and_then(|tree| tree.cell_types.clone()),
    )
    .with_left_out(&left.names)
    .write(output_dir)
    .context("Failed to write info.json")?;

    if let Some(work) = &inputs.work {
        work.remove()?;
    }
    let listed = vep_core::cache_manifest::CacheManifest::write(output_dir)
        .context("Failed to write CHECKSUMS.sha256")?;
    info!("Checksum manifest written for {}", counted(listed, "file"));

    Ok(Built {
        variant_cells,
        sequence_features,
        regulatory: regulatory_tree,
        left_out: left,
    })
}

/// What [`write_variation_files`] wrote.
struct VariationFiles {
    /// The columns `info.json` lists as `variation_cols`.
    columns: Vec<String>,
    /// What was read and written for the `pubmed`, `var_synonyms` and
    /// `failed` cells.
    cells: Option<variant_cells::VariantCellFacts>,
    /// What was read of the `variation_feature` table.
    features: Option<SequenceFeatureFacts>,
}

/// Writes the variation file of each chromosome in `variations` into
/// `output_dir` (`variations/<chr>/all_vars.gz` with its `.csi`), the ClinVar
/// columns filled from the variation table dumps when given and the `pubmed`,
/// `var_synonyms` and `failed` cells from the files named for them.
///
/// Where the germline records come from one genome-wide file and the
/// variation database's `variation_feature` table was found, each of
/// `cache_sequences`, the sequences the cache holds transcripts on, that no
/// variation VCF holds a record on is written from that table instead; the
/// table is fetched only when there is such a sequence.
fn write_variation_files(
    variations: &VariationBuild,
    output_dir: &Path,
    cache_sequences: &[String],
) -> Result<VariationFiles> {
    let variation_inputs = &variations.inputs;
    info!(
        "Resolved the variation VCFs of {}",
        counted(variation_inputs.per_chromosome.len(), "chromosome")
    );
    let in_question = match &variations.features {
        Some(features) => Some(sequences_without_a_vcf_record(
            variation_inputs,
            cache_sequences,
            features.named.as_deref(),
        )?),
        None => None,
    };
    let from_tables: &[String] = in_question.as_deref().unwrap_or_default();
    let mut chromosomes: Vec<String> = variation_inputs
        .per_chromosome
        .iter()
        .map(|(chrom, _)| chrom.clone())
        .collect();
    for sequence in from_tables {
        if !chromosomes.contains(sequence) {
            chromosomes.push(sequence.clone());
        }
    }

    let clinvar = match &variations.dumps {
        Some(dir) => {
            info!("Reading the ClinVar phenotype features: {}", dir.display());
            let attribs = variation::ClinVarAttribs::load(dir, &chromosomes)
                .context("Failed to read the variation table dumps")?;
            Some(attribs)
        }
        None => None,
    };
    let cells = match &variations.cells {
        Some(inputs) => Some(
            variant_cells::VariantCells::read(inputs, &chromosomes).context(
                "Failed to read the PubMed identifiers, variant synonyms and failed variants",
            )?,
        ),
        None => None,
    };

    let columns = variation::variation_columns(variation_inputs.genome_wide.frequencies.is_some());
    let mut total_variations = 0usize;
    let mut files = 0usize;
    for (chrom, vcf) in &variation_inputs.per_chromosome {
        if from_tables.contains(chrom) {
            continue;
        }
        info!(
            "Reading the variation records of chromosome {}: {}",
            chrom,
            vcf.path().display()
        );
        total_variations += variation::write_chromosome(
            chrom,
            vcf,
            &variation_inputs.genome_wide,
            clinvar.as_ref(),
            cells.as_ref(),
            &columns,
            output_dir,
        )
        .with_context(|| format!("Failed to write the variation file of chromosome {}", chrom))?;
        files += 1;
    }

    let mut brought = Vec::new();
    if let (Some(features), false) = (&variations.features, from_tables.is_empty()) {
        let tables = (features.fetch)().context("Failed to fetch the variation_feature table")?;
        let written = variation::write_sequences_from_tables(
            &tables,
            from_tables,
            &variation_inputs.genome_wide,
            clinvar.as_ref(),
            cells.as_ref(),
            &columns,
            output_dir,
        )
        .context("Failed to write the variation files of the variation_feature table")?;
        for (sequence, records) in written {
            if records > 0 {
                total_variations += records;
                files += 1;
                brought.push((sequence, records));
            } else if variation_inputs
                .per_chromosome
                .iter()
                .any(|(chrom, _)| *chrom == sequence)
            {
                // A sequence `--chromosomes` names has a variation file
                // whatever it holds.
                variation::AllVarsWriter::create(output_dir, &sequence, &columns)?.finish()?;
                files += 1;
            }
        }
    }
    info!(
        "Variation files written: {} across {}",
        counted(total_variations, "record"),
        counted(files, "chromosome")
    );
    Ok(VariationFiles {
        columns,
        cells: cells.map(|cells| cells.facts()),
        features: in_question.map(|in_question| SequenceFeatureFacts {
            in_question,
            brought,
        }),
    })
}

/// The sequences of `cache_sequences` that no variation VCF of `inputs`
/// holds a record on, among `named` when it names any; none where the
/// germline records come from one file a chromosome.
fn sequences_without_a_vcf_record(
    inputs: &download::VariationInputs,
    cache_sequences: &[String],
    named: Option<&[String]>,
) -> Result<Vec<String>> {
    let mut vcfs: Vec<&Path> = Vec::new();
    for (_, vcf) in &inputs.per_chromosome {
        match vcf {
            variation::GermlineVcf::Chromosome(_) => return Ok(Vec::new()),
            variation::GermlineVcf::Genome(path) if !vcfs.contains(&path.as_path()) => {
                vcfs.push(path)
            }
            variation::GermlineVcf::Genome(_) => {}
        }
    }
    vcfs.extend(inputs.genome_wide.somatic.as_deref());
    vcfs.extend(inputs.genome_wide.germline.iter().map(PathBuf::as_path));
    let mut covered: Vec<String> = Vec::new();
    for vcf in vcfs {
        covered.extend(variation::sequence_names(vcf)?);
    }
    Ok(cache_sequences
        .iter()
        .filter(|sequence| !covered.contains(sequence))
        .filter(|sequence| named.is_none_or(|named| named.contains(sequence)))
        .cloned()
        .collect())
}

/// Whether the log lines carry colour codes: only where standard error is a
/// terminal and `NO_COLOR` is unset or empty, so a log sent to a file or a
/// pipe is plain text. `with_ansi` replaces the subscriber's own reading of
/// `NO_COLOR`, so the variable is read here.
fn log_lines_carry_colour(terminal: bool, no_color: Option<&std::ffi::OsStr>) -> bool {
    terminal && no_color.is_none_or(std::ffi::OsStr::is_empty)
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    /// Colour codes go to a terminal alone, and not there once `NO_COLOR` holds
    /// a value; a file or a pipe takes none whatever the variable holds.
    #[test]
    fn log_lines_carry_colour_on_a_terminal_alone_and_not_under_no_color() {
        let set = std::ffi::OsStr::new("1");
        let empty = std::ffi::OsStr::new("");
        assert!(log_lines_carry_colour(true, None));
        assert!(log_lines_carry_colour(true, Some(empty)));
        assert!(!log_lines_carry_colour(true, Some(set)));
        for no_color in [None, Some(empty), Some(set)] {
            assert!(!log_lines_carry_colour(false, no_color));
        }
    }

    /// The sequences a build reads the `variation_feature` table for are
    /// those of the cache that no variation VCF holds a record on, the
    /// genome-wide files read beside the germline one among them, and under
    /// `--chromosomes` those of them it names. Where the germline files are
    /// one a chromosome there is none, whatever those files hold.
    #[test]
    fn the_sequences_read_from_the_table_are_those_one_genome_wide_vcf_leaves_without_a_record() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = |name: &str, sequences: &[&str]| -> PathBuf {
            let records: String = sequences
                .iter()
                .map(|sequence| format!("{sequence}\t100\tmv100\tA\tG\t.\t.\t.\n"))
                .collect();
            let path = dir.path().join(name);
            std::fs::write(
                &path,
                format!(
                    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     {records}"
                ),
            )
            .unwrap();
            path
        };
        let names =
            |names: &[&str]| -> Vec<String> { names.iter().map(|n| n.to_string()).collect() };
        let cache = names(&["1", "2", "3", "MT"]);
        let genome = vcf("made_up.vcf", &["1", "2", "9"]);
        let one_file =
            |built: &[&str], genome_wide: variation::GenomeWideInputs| download::VariationInputs {
                per_chromosome: built
                    .iter()
                    .map(|chrom| {
                        (
                            chrom.to_string(),
                            variation::GermlineVcf::Genome(genome.clone()),
                        )
                    })
                    .collect(),
                genome_wide,
            };
        let none = variation::GenomeWideInputs::default;

        let inputs = one_file(&["1", "2", "9"], none());
        assert_eq!(
            sequences_without_a_vcf_record(&inputs, &cache, None).unwrap(),
            names(&["3", "MT"])
        );
        // A sequence `--chromosomes` names and the VCF holds no record on.
        let inputs = one_file(&["2", "MT"], none());
        assert_eq!(
            sequences_without_a_vcf_record(&inputs, &cache, Some(&names(&["2", "MT"]))).unwrap(),
            names(&["MT"])
        );
        assert_eq!(
            sequences_without_a_vcf_record(&inputs, &cache, Some(&names(&["1", "2"]))).unwrap(),
            names(&[])
        );
        // A record of a genome-wide file read beside the germline one.
        let beside = variation::GenomeWideInputs {
            somatic: Some(vcf("made_up_somatic.vcf", &["3"])),
            germline: vec![vcf("made_up_associated.vcf", &["MT"])],
            ..Default::default()
        };
        let inputs = one_file(&["1", "2", "9"], beside);
        assert_eq!(
            sequences_without_a_vcf_record(&inputs, &cache, None).unwrap(),
            names(&[])
        );

        let a_file_a_chromosome = download::VariationInputs {
            per_chromosome: vec![(
                "1".to_string(),
                variation::GermlineVcf::Chromosome(vcf("made_up-chr1.vcf", &["1"])),
            )],
            genome_wide: none(),
        };
        assert_eq!(
            sequences_without_a_vcf_record(&a_file_a_chromosome, &cache, None).unwrap(),
            names(&[])
        );
    }

    fn parse(argv: &[&str]) -> Result<Args, clap::Error> {
        Args::try_parse_from(std::iter::once("vep-cache-builder").chain(argv.iter().copied()))
    }

    /// A local copy of the site's tree holding the committed chromosome 21
    /// slice as release 116 of `homo_sapiens` GRCh38: the GFF3 and, with
    /// `peptides`, the protein FASTA, each gzipped under its published name,
    /// with the directory of the core database, by which the genome is found
    /// in a tree that publishes no species list.
    fn slice_site(root: &Path, peptides: bool) -> install::Site {
        let genome = Genome {
            layout: download::Layout::Release,
            release: 116,
            assembly: "GRCh38",
            species: "homo_sapiens",
            core_db: "homo_sapiens_core_116_38",
        };
        gene_set_site(root, &genome, peptides)
    }

    /// A genome as a release tree publishes it.
    struct Genome {
        layout: download::Layout,
        release: u32,
        assembly: &'static str,
        species: &'static str,
        core_db: &'static str,
    }

    /// [`slice_site`] for any genome: the same slice under the names `genome`
    /// is published by.
    fn gene_set_site(root: &Path, genome: &Genome, peptides: bool) -> install::Site {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21");
        let Genome {
            layout,
            release,
            assembly,
            species,
            core_db,
        } = genome;
        let core = layout.mysql_dir(*release, core_db);
        std::fs::create_dir_all(root.join(core.unwrap())).unwrap();
        let mut files = vec![(
            layout.gff3_location(*release, assembly, species),
            "slice.gff3",
        )];
        if peptides {
            files.push((
                layout.protein_fasta_location(*release, assembly, species),
                "peptides.fa",
            ));
        }
        for ((dir, name), source) in files {
            let dir = root.join(dir);
            std::fs::create_dir_all(&dir).unwrap();
            let mut encoder = flate2::write::GzEncoder::new(
                std::fs::File::create(dir.join(name)).unwrap(),
                flate2::Compression::fast(),
            );
            std::io::Write::write_all(&mut encoder, &std::fs::read(fixture.join(source)).unwrap())
                .unwrap();
            encoder.finish().unwrap();
        }
        install::Site::Dir(root.to_path_buf())
    }

    /// A build that names no GFF3 and no protein FASTA fetches both into
    /// `tmp/vep_cache_builder/<species>_<release>_<assembly>` under its output
    /// directory. A
    /// build that fails keeps that directory and names it at the end of its
    /// error; the next build reuses the file it holds, and a build that
    /// succeeds removes the directory and the `tmp` it leaves empty before the
    /// checksum manifest is written, so the manifest lists no download. An
    /// output directory itself named `tmp` builds the same way.
    #[test]
    fn a_build_removes_its_work_directory_and_a_failed_build_keeps_it() {
        for name in ["GRCh38", "tmp"] {
            let root = tempfile::tempdir().unwrap();
            let site_dir = root.path().join("site");
            let output_dir = root.path().join("caches").join(name);
            let work_dir = output_dir.join("tmp/vep_cache_builder/homo_sapiens_116_GRCh38");
            let args = parse(&[
                "--release",
                "116",
                "--assembly",
                "GRCh38",
                "--output-dir",
                output_dir.to_str().unwrap(),
            ])
            .unwrap();
            let build = |site: &install::Site| {
                build_named_inputs(&args, None, HashMap::new(), std::slice::from_ref(site))
            };

            let error = build(&slice_site(&site_dir, false)).unwrap_err();
            let error = format!("{error:#}");
            assert!(
                error.starts_with("Failed to resolve the protein FASTA: "),
                "{error}"
            );
            assert!(
                error.ends_with(&format!(
                    "is not on the site; the work directory {} is kept for the next run",
                    work_dir.display()
                )),
                "{error}"
            );
            let gff3 = work_dir.join("vep_cache_builder_homo_sapiens_GRCh38_116.gff3.gz");
            assert!(gff3.is_file());
            assert!(!output_dir.join("info.json").exists());

            let site = slice_site(&site_dir, true);
            let (dir, file) =
                download::Layout::Release.gff3_location(116, "GRCh38", "homo_sapiens");
            std::fs::remove_file(site_dir.join(dir).join(file)).unwrap();
            build(&site).unwrap();
            assert!(output_dir.join("info.json").is_file());
            assert!(!output_dir.join("tmp").exists());
            let manifest = std::fs::read_to_string(output_dir.join("CHECKSUMS.sha256")).unwrap();
            assert!(manifest.contains("  info.json"), "{manifest}");
            assert!(!manifest.contains("tmp/"), "{manifest}");
            assert!(!manifest.contains("vep_cache_builder_"), "{manifest}");
            // The check `vep` makes of a cache directory before it reads it.
            let listed = vep_core::cache_manifest::CacheManifest::load(&output_dir).unwrap();
            assert!(listed.is_present());
            listed.require_complete().unwrap();
        }
    }

    /// A build from local files alone downloads nothing and creates no work
    /// directory; one a failed build kept in the output directory is removed
    /// before the manifest is written, and a file of the output directory's
    /// own `tmp` that is not the work directory stays.
    #[test]
    fn a_build_from_local_files_clears_a_kept_work_directory() {
        let root = tempfile::tempdir().unwrap();
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21");
        let output_dir = root.path().join("cache");
        let work_dir = output_dir.join("tmp/vep_cache_builder/homo_sapiens_116_GRCh38");
        std::fs::create_dir_all(&work_dir).unwrap();
        std::fs::write(work_dir.join("vep_cache_builder_kept.vcf.gz"), b"kept").unwrap();
        std::fs::write(output_dir.join("tmp/notes.txt"), b"notes").unwrap();
        let args = parse(&[
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--gff3",
            fixture.join("slice.gff3").to_str().unwrap(),
            "--fasta",
            fixture.join("peptides.fa").to_str().unwrap(),
            "--output-dir",
            output_dir.to_str().unwrap(),
        ])
        .unwrap();
        // No file is fetched, so the site is never read.
        let site = install::Site::Dir(root.path().join("no-site"));
        build_named_inputs(&args, None, HashMap::new(), std::slice::from_ref(&site)).unwrap();
        assert!(!work_dir.exists());
        assert!(output_dir.join("tmp/notes.txt").is_file());
        let manifest = std::fs::read_to_string(output_dir.join("CHECKSUMS.sha256")).unwrap();
        assert!(!manifest.contains("vep_cache_builder_"), "{manifest}");
        assert!(manifest.contains("  tmp/notes.txt"), "{manifest}");
    }

    /// The paths a checksum manifest lists, in its order.
    fn listed_paths(output_dir: &Path) -> Vec<String> {
        std::fs::read_to_string(output_dir.join("CHECKSUMS.sha256"))
            .unwrap()
            .lines()
            .map(|line| line.split_once("  ").unwrap().1.to_string())
            .collect()
    }

    /// The files under `dir`, at any depth, whose name is `name`.
    fn files_named(dir: &Path, name: &str) -> usize {
        std::fs::read_dir(dir)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .map(|path| {
                if path.is_dir() {
                    files_named(&path, name)
                } else {
                    usize::from(path.file_name().is_some_and(|file| file == name))
                }
            })
            .sum()
    }

    /// A build that succeeds removes the downloads a failed build of another
    /// release, assembly or species kept in its output directory: its
    /// checksum manifest lists the files of a build into an empty directory
    /// and no other, no `tmp` is left, and the directory passes the check
    /// `vep` makes of it.
    #[test]
    fn a_build_removes_the_downloads_a_failed_build_of_another_genome_kept() {
        let root = tempfile::tempdir().unwrap();
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21");
        let from_local_files = |output_dir: &Path| {
            let args = parse(&[
                "--release",
                "116",
                "--assembly",
                "GRCh38",
                "--gff3",
                fixture.join("slice.gff3").to_str().unwrap(),
                "--fasta",
                fixture.join("peptides.fa").to_str().unwrap(),
                "--output-dir",
                output_dir.to_str().unwrap(),
            ])
            .unwrap();
            // No file is fetched, so the site is never read.
            let site = install::Site::Dir(root.path().join("no-site"));
            build_named_inputs(&args, None, HashMap::new(), std::slice::from_ref(&site)).unwrap();
        };
        let clean = root.path().join("clean");
        from_local_files(&clean);
        let expected = listed_paths(&clean);
        assert_eq!(
            expected,
            [
                "info.json",
                "transcripts/21/10000001-11000000.json",
                "transcripts/21/13000001-14000000.json",
                "transcripts/21/17000001-18000000.json",
                "transcripts/21/18000001-19000000.json",
                "transcripts/21/30000001-31000000.json",
                "transcripts/21/40000001-41000000.json",
                "transcripts/21/5000001-6000000.json",
            ]
        );

        let other_genomes = [
            Genome {
                layout: download::Layout::Release,
                release: 115,
                assembly: "GRCh38",
                species: "homo_sapiens",
                core_db: "homo_sapiens_core_115_38",
            },
            Genome {
                layout: download::Layout::Grch37,
                release: 116,
                assembly: "GRCh37",
                species: "homo_sapiens",
                core_db: "homo_sapiens_core_116_37",
            },
            Genome {
                layout: download::Layout::Release,
                release: 116,
                assembly: "GRCm39",
                species: "mus_musculus",
                core_db: "mus_musculus_core_116_39",
            },
        ];
        for genome in &other_genomes {
            let Genome {
                release,
                assembly,
                species,
                ..
            } = genome;
            let case = root.path().join(format!("{species}_{release}_{assembly}"));
            let output_dir = case.join("cache");
            let args = parse(&[
                "--release",
                &release.to_string(),
                "--assembly",
                assembly,
                "--species",
                species,
                "--output-dir",
                output_dir.to_str().unwrap(),
            ])
            .unwrap();
            // The site publishes the genome's GFF3 and no protein FASTA: the
            // build fetches the one and fails on the other.
            let site = gene_set_site(&case.join("site"), genome, false);
            let error =
                build_named_inputs(&args, None, HashMap::new(), std::slice::from_ref(&site))
                    .unwrap_err();
            let error = format!("{error:#}");
            assert!(error.ends_with(" is kept for the next run"), "{error}");
            let download = format!("vep_cache_builder_{species}_{assembly}_{release}.gff3.gz");
            assert_eq!(
                files_named(&output_dir.join("tmp"), &download),
                1,
                "{species} {release} {assembly}"
            );

            from_local_files(&output_dir);
            assert_eq!(
                listed_paths(&output_dir),
                expected,
                "{species} {release} {assembly}"
            );
            assert!(!output_dir.join("tmp").exists());
            // The check `vep` makes of a cache directory before it reads it.
            let listed = vep_core::cache_manifest::CacheManifest::load(&output_dir).unwrap();
            assert!(listed.is_present());
            listed.require_complete().unwrap();
        }
    }

    /// The line `vep` ends on when it is run against `cache`.
    fn vep_refusal(dir: &Path, cache: &Path) -> String {
        let input = dir.join("input.vcf");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             21\t1500\trecord\tA\tG\t.\t.\t.\n",
        )
        .unwrap();
        let out = dir.join("out.txt");
        let argv = [
            "vep".to_string(),
            "--offline".to_string(),
            "--json_cache".to_string(),
            cache.display().to_string(),
            "--no_stats".to_string(),
            "--quiet".to_string(),
            "--force_overwrite".to_string(),
            "-i".to_string(),
            input.display().to_string(),
            "-o".to_string(),
            out.display().to_string(),
        ];
        let config = vep_cli::Config::from_argv(argv).expect("the command line is accepted");
        let error = vep_cli::runner::Runner::new(config)
            .and_then(|runner| runner.run())
            .expect_err("vep refuses the directory");
        vep_cli::config::CliError::from_run(&error)
    }

    /// A build that fetched its inputs and then fails after its first shard
    /// keeps both what it downloaded and the mark of a build that has not
    /// ended: the work directory holds the two files and is named at the end
    /// of the error, the manifest is the one line that lists no file, loading
    /// the output directory is refused by name, and `vep` ends on that
    /// refusal.
    #[test]
    fn a_build_that_fails_after_its_downloads_keeps_them_and_the_unfinished_mark() {
        let root = tempfile::tempdir().unwrap();
        let site = slice_site(&root.path().join("site"), true);
        let output_dir = root.path().join("cache");
        let work_dir = output_dir.join("tmp/vep_cache_builder/homo_sapiens_116_GRCh38");
        // The regulatory bake, which follows the transcript shards, refuses an
        // `ID` that decodes to bytes that are not UTF-8.
        let regulatory = root.path().join("refused.gff3");
        std::fs::write(
            &regulatory,
            "21\tEnsembl\tenhancer\t1000\t1100\t.\t.\t.\tID=caf%FF\n",
        )
        .unwrap();
        let args = parse(&[
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--regulatory-gff",
            regulatory.to_str().unwrap(),
            "--output-dir",
            output_dir.to_str().unwrap(),
        ])
        .unwrap();

        let error = build_named_inputs(&args, None, HashMap::new(), std::slice::from_ref(&site))
            .unwrap_err();
        let error = format!("{error:#}");
        assert!(
            error.starts_with("Failed to build the regulatory tree: "),
            "{error}"
        );
        assert!(
            error.ends_with(&format!(
                "; the work directory {} is kept for the next run",
                work_dir.display()
            )),
            "{error}"
        );
        for name in [
            "vep_cache_builder_homo_sapiens_GRCh38_116.gff3.gz",
            "vep_cache_builder_homo_sapiens_GRCh38_116.pep.all.fa.gz",
        ] {
            assert!(work_dir.join(name).is_file(), "{name}");
        }
        assert!(output_dir.join("transcripts").is_dir());

        let refusal = format!(
            "{} holds a cache whose build did not finish: build it again",
            output_dir.display()
        );
        let unfinished = vep_core::cache_manifest::CacheManifest::load(&output_dir).unwrap_err();
        assert!(
            matches!(
                unfinished,
                vep_core::cache_manifest::ManifestError::Unfinished { .. }
            ),
            "{unfinished}"
        );
        assert_eq!(unfinished.to_string(), refusal);
        assert_eq!(
            vep_refusal(root.path(), &output_dir),
            format!("ERROR: [cache] {refusal}")
        );
    }

    /// A build from named inputs takes `--species` from the default of
    /// `vep --species`, Ensembl VEP's own (`Config.pm:290`), so a cache built
    /// and a run made without the flag name the same species. `--auto` has no
    /// default: like `INSTALL.pl --AUTO c`, which dies with `No species
    /// specified` (`INSTALL.pl:299-300`), it takes the species as a value and
    /// is refused without one.
    #[test]
    fn species_defaults_to_the_runs_default_and_auto_has_none() {
        let explicit = parse(&[
            "--release",
            "116",
            "--assembly",
            "GRCm39",
            "--output-dir",
            "out",
        ])
        .unwrap();
        assert_eq!(explicit.species, "homo_sapiens");
        let named = parse(&[
            "--release",
            "116",
            "--assembly",
            "GRCm39",
            "--species",
            "mus_musculus",
            "--output-dir",
            "out",
        ])
        .unwrap();
        assert_eq!(named.species, "mus_musculus");
        for refused in [
            vec!["--auto", "--release", "116"],
            vec!["--auto", "GRCm39", "--release", "116"],
        ] {
            assert!(parse(&refused).is_err(), "{refused:?} parses");
        }
        let auto = parse(&["--auto", "mus_musculus", "GRCm39", "--release", "116"]).unwrap();
        assert_eq!(
            auto.auto.as_deref(),
            Some(&["mus_musculus".to_string(), "GRCm39".to_string()][..])
        );
    }

    /// INSTALL.pl's `--HELP` (`:170`) is `--help`, which prints the usage and
    /// builds nothing; its plugin and BioPerl options (`--NO_PLUGINS`,
    /// `--SKIP_PLUGINS`, `--BIOPERL`; `:173`, `:175`, `:166`) are no flag of
    /// the builder in either letter case.
    #[test]
    fn installer_help_is_help_and_its_plugin_and_bioperl_options_are_no_flag() {
        use clap::error::ErrorKind;
        let kind = |argv: &[&str]| parse(argv).err().map(|error| error.kind());
        assert_eq!(kind(&["--help"]), Some(ErrorKind::DisplayHelp));
        for option in [
            &["--NO_PLUGINS"][..],
            &["--SKIP_PLUGINS", "CADD"],
            &["--BIOPERL", "bioperl.tar.gz"],
        ] {
            for name in [option[0].to_string(), option[0].to_lowercase()] {
                let mut argv = vec!["--auto", "homo_sapiens", "GRCh38", "--release", "116"];
                argv.push(&name);
                argv.extend(&option[1..]);
                assert_eq!(kind(&argv), Some(ErrorKind::UnknownArgument), "{name}");
            }
        }
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
        assert_eq!(args.site.as_deref(), Some("/site"));
        assert!(args.include_predictions && args.include_variations);
        assert_eq!(
            args.chromosomes,
            Some(vec!["21".to_string(), "MT".to_string()])
        );
        assert_eq!(args.codon_table, ["Pt=11", "Mt=1", "MT=2"]);
        let plain = parse(&["--auto", "homo_sapiens", "GRCh38", "--release", "116"]).unwrap();
        assert!(plain.site.is_none());
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
