// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One-command installation of a release (`--auto`): every input is fetched
//! from the Ensembl site, or a copy of its tree, and checked against the
//! `CHECKSUMS` file of its directory; the cache is built in Ensembl VEP's
//! layout, `<dir_cache>/<species>/<release>_<assembly>`, with the indexed
//! genome FASTA inside it, where `vep` finds both through `--dir_cache`.
//!
//! The downloads land under `<dir_cache>/tmp/<species>_<release>_<assembly>/`
//! and are removed once the cache is in place; a downloaded file whose checksum
//! still matches is reused by the next run, so an interrupted installation
//! resumes where it stopped. A release directory that already exists is never
//! written into.

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use tracing::{info, warn};

use crate::core_dumps::{table_path, CORE_TABLES, VARIATION_TABLES};
use crate::download;
use crate::{fasta_index, BuildInputs, MysqlServer, VariationBuild};

/// Where the release's files are read from: the Ensembl site or another server
/// publishing the same tree, or a local directory holding that tree.
#[derive(Debug, Clone)]
pub enum Site {
    Url(String),
    Dir(PathBuf),
}

impl Site {
    /// A URL (`http://` or `https://`, a trailing slash dropped) or an existing
    /// directory.
    pub fn parse(value: &str) -> Result<Self> {
        if value.starts_with("http://") || value.starts_with("https://") {
            return Ok(Self::Url(value.trim_end_matches('/').to_string()));
        }
        let dir = PathBuf::from(value);
        if !dir.is_dir() {
            bail!("--site {value} is neither a URL nor a directory");
        }
        Ok(Self::Dir(dir))
    }

    /// The text of the file at `relative`, or `None` when the site has no such
    /// file.
    fn read_text(&self, relative: &str) -> Result<Option<String>> {
        match self {
            Self::Url(root) => download::fetch_text(&format!("{root}/{relative}")),
            Self::Dir(root) => {
                let path = root.join(relative);
                if !path.is_file() {
                    return Ok(None);
                }
                Ok(Some(std::fs::read_to_string(&path).with_context(|| {
                    format!("Failed to read {}", path.display())
                })?))
            }
        }
    }

    /// Places the file at `relative` at `dest`: downloaded from a URL, copied
    /// from a directory.
    fn place(&self, relative: &str, dest: &Path) -> Result<()> {
        match self {
            Self::Url(root) => download::download_file(&format!("{root}/{relative}"), dest),
            Self::Dir(root) => {
                let source = root.join(relative);
                if !source.is_file() {
                    bail!("{} is not on the site", source.display());
                }
                info!("Copying {}", source.display());
                std::fs::copy(&source, dest).with_context(|| {
                    format!("Failed to copy {} to {}", source.display(), dest.display())
                })?;
                Ok(())
            }
        }
    }
}

/// What `--auto` installs.
pub struct InstallOptions {
    pub species: String,
    pub assembly: String,
    pub release: u32,
    /// The cache root the release directory is created under.
    pub dir_cache: PathBuf,
    pub site: Site,
    pub include_variations: bool,
    /// The chromosomes whose variation VCFs are fetched;
    /// [`download::default_chromosomes`] of the species and its GFF3 when
    /// `None`.
    pub chromosomes: Option<Vec<String>>,
    /// `--codon-table`: the codon table of each sequence region the core
    /// dumps state none for.
    pub codon_tables: HashMap<String, u8>,
    /// The server the SIFT and PolyPhen matrices come from, under
    /// `--include-predictions`.
    pub predictions: Option<MysqlServer>,
}

/// One `CHECKSUMS` line: the BSD `sum` of the file (a 16-bit checksum and its
/// size in 1 KiB blocks, rounded up).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Checksum {
    pub sum: u16,
    pub blocks: u64,
}

impl std::fmt::Display for Checksum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:05} {}", self.sum, self.blocks)
    }
}

/// The `CHECKSUMS` file of one directory of the site: one line per file,
/// `<sum> <blocks> <name>`, the directory's own listing.
#[derive(Debug, Default)]
pub struct Listing {
    entries: HashMap<String, Checksum>,
}

impl Listing {
    /// Parses the file's text; a line without the three fields, or whose first
    /// two are not numbers, is passed over.
    pub fn parse(text: &str) -> Self {
        let entries = text
            .lines()
            .filter_map(|line| {
                let mut fields = line.split_whitespace();
                let sum = fields.next()?.parse::<u16>().ok()?;
                let blocks = fields.next()?.parse::<u64>().ok()?;
                let name = fields.next()?;
                Some((name.to_string(), Checksum { sum, blocks }))
            })
            .collect();
        Self { entries }
    }

    pub fn get(&self, name: &str) -> Option<Checksum> {
        self.entries.get(name).copied()
    }

    pub fn contains(&self, name: &str) -> bool {
        self.entries.contains_key(name)
    }

    /// The listed names, sorted.
    pub fn names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.entries.keys().map(String::as_str).collect();
        names.sort_unstable();
        names
    }
}

/// The BSD `sum` of a file: the 16-bit checksum (each byte added after the
/// running value rotates right by one bit) and the size in 1 KiB blocks, rounded
/// up, the two numbers a `CHECKSUMS` line carries.
pub fn bsd_sum(path: &Path) -> Result<Checksum> {
    let mut file =
        File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
    let mut buffer = vec![0u8; 1 << 20];
    let mut sum: u16 = 0;
    let mut bytes: u64 = 0;
    loop {
        let read = file
            .read(&mut buffer)
            .with_context(|| format!("Failed to read {}", path.display()))?;
        if read == 0 {
            break;
        }
        for &byte in &buffer[..read] {
            sum = sum.rotate_right(1).wrapping_add(u16::from(byte));
        }
        bytes += read as u64;
    }
    Ok(Checksum {
        sum,
        blocks: bytes.div_ceil(1024),
    })
}

/// Refuses `path` when its BSD sum differs from the `CHECKSUMS` line for
/// `name`.
fn verify(path: &Path, name: &str, expected: Checksum) -> Result<()> {
    let actual = bsd_sum(path)?;
    if actual != expected {
        bail!(
            "checksum for {} ({actual}) does not match the CHECKSUMS line for {name} ({expected})",
            path.display()
        );
    }
    Ok(())
}

/// The files of one directory of the site, placed and verified on demand.
#[derive(Debug)]
struct Directory<'a> {
    site: &'a Site,
    relative: String,
    listing: Listing,
}

impl Directory<'_> {
    /// Places the listed file `name` in `dest_dir` and verifies it against the
    /// listing; a file already there that verifies is kept, one that does not
    /// is fetched again. Refuses a name the listing lacks.
    fn fetch(&self, name: &str, dest_dir: &Path) -> Result<PathBuf> {
        let Some(expected) = self.listing.get(name) else {
            bail!("{name} is not listed in {}/CHECKSUMS", self.relative);
        };
        std::fs::create_dir_all(dest_dir)
            .with_context(|| format!("Failed to create {}", dest_dir.display()))?;
        let dest = dest_dir.join(name);
        if dest.is_file() {
            if verify(&dest, name, expected).is_ok() {
                info!("Verified {}", dest.display());
                return Ok(dest);
            }
            warn!(
                "{} does not match its CHECKSUMS line; fetching it again",
                dest.display()
            );
            std::fs::remove_file(&dest)?;
        }
        self.site
            .place(&format!("{}/{name}", self.relative), &dest)?;
        if let Err(error) = verify(&dest, name, expected) {
            let _ = std::fs::remove_file(&dest);
            return Err(error);
        }
        info!("Verified {}", dest.display());
        Ok(dest)
    }
}

/// Reads the `CHECKSUMS` of the directory at `relative`; `None` when the site
/// has no such directory.
fn directory<'a>(site: &'a Site, relative: &str) -> Result<Option<Directory<'a>>> {
    let text = site
        .read_text(&format!("{relative}/CHECKSUMS"))
        .with_context(|| format!("Failed to read {relative}/CHECKSUMS"))?;
    Ok(text.map(|text| Directory {
        site,
        relative: relative.to_string(),
        listing: Listing::parse(&text),
    }))
}

/// As [`directory`], refusing a directory the site lacks.
fn required_directory<'a>(site: &'a Site, relative: &str) -> Result<Directory<'a>> {
    directory(site, relative)?.ok_or_else(|| {
        anyhow::anyhow!("{relative}/CHECKSUMS is not on the site: is the release published?")
    })
}

/// The release directory `vep` resolves: `<dir_cache>/<species>/<release>_<assembly>`.
pub fn release_dir(options: &InstallOptions) -> PathBuf {
    options
        .dir_cache
        .join(&options.species)
        .join(format!("{}_{}", options.release, options.assembly))
}

/// The genome FASTA of the release in `listing`: the primary assembly, else the
/// toplevel file, neither soft- nor hard-masked, as Ensembl VEP's installer
/// chooses it.
pub fn genome_fasta_name(listing: &Listing) -> Option<String> {
    let candidates: Vec<&str> = listing
        .names()
        .into_iter()
        .filter(|name| !name.contains("_sm.") && !name.contains("_rm."))
        .collect();
    ["primary_assembly.fa.gz", "toplevel.fa.gz"]
        .iter()
        .find_map(|suffix| {
            candidates
                .iter()
                .find(|name| name.ends_with(suffix))
                .map(|name| name.to_string())
        })
}

/// Installs the release: fetches and verifies every input, places the indexed
/// genome FASTA in the cache directory under construction, builds the cache
/// there and moves it to [`release_dir`]. Returns that directory.
pub fn install(options: &InstallOptions) -> Result<PathBuf> {
    let species = options.species.as_str();
    let assembly = options.assembly.as_str();
    let release = options.release;
    let site = &options.site;
    let dir_cache = &options.dir_cache;
    let target = release_dir(options);
    if target.exists() {
        bail!(
            "{} exists: remove it to install the release again",
            target.display()
        );
    }
    let work_dir = dir_cache
        .join("tmp")
        .join(format!("{species}_{release}_{assembly}"));
    let inputs_dir = work_dir.join("inputs");
    let core_dumps_dir = work_dir.join("core_dumps");
    let cache_dir = work_dir.join("cache");
    if cache_dir.exists() {
        std::fs::remove_dir_all(&cache_dir)
            .with_context(|| format!("Failed to remove {}", cache_dir.display()))?;
    }
    std::fs::create_dir_all(&cache_dir)
        .with_context(|| format!("Failed to create {}", cache_dir.display()))?;
    info!(
        "Installing {species} {assembly} release {release} into {}",
        target.display()
    );

    let (gff3_dir, gff3_name) = download::gff3_location(release, assembly, species);
    let gff3 = required_directory(site, &gff3_dir)?.fetch(&gff3_name, &inputs_dir)?;
    let (gtf_dir, gtf_name) = download::gtf_location(release, assembly, species);
    let gtf = required_directory(site, &gtf_dir)?.fetch(&gtf_name, &inputs_dir)?;
    let (pep_dir, pep_name) = download::protein_fasta_location(release, assembly, species);
    let protein_fasta = required_directory(site, &pep_dir)?.fetch(&pep_name, &inputs_dir)?;

    let core = required_directory(
        site,
        &download::mysql_dir(release, assembly, species, "core")?,
    )?;
    for table in CORE_TABLES {
        fetch_table(&core, table, &core_dumps_dir)?.ok_or_else(|| {
            anyhow::anyhow!("core table {table} is not published for this release")
        })?;
    }
    let variation_db = directory(
        site,
        &download::mysql_dir(release, assembly, species, "variation")?,
    )?;
    match &variation_db {
        Some(variation) => {
            for table in VARIATION_TABLES {
                if fetch_table(variation, table, &core_dumps_dir)?.is_none() {
                    warn!(
                        "Variation table {table} is not published for this release and assembly; gene_phenotype is 0 on every transcript"
                    );
                }
            }
        }
        None => warn!(
            "No variation database is published for this release and assembly; gene_phenotype is 0 on every transcript"
        ),
    }

    let genome_fasta =
        install_genome_fasta(site, species, assembly, release, &inputs_dir, &cache_dir)?;

    let variations = if options.include_variations {
        Some(fetch_variations(
            options,
            variation_db.as_ref(),
            &work_dir,
            &core_dumps_dir,
            &gff3,
        )?)
    } else {
        None
    };

    let inputs = BuildInputs {
        species: species.to_string(),
        assembly: assembly.to_string(),
        release,
        gff3: Some(gff3),
        protein_fasta: Some(protein_fasta),
        gtf: Some(gtf),
        core_dumps: Some(core_dumps_dir),
        codon_tables: options.codon_tables.clone(),
        genome_fasta: Some(genome_fasta),
        predictions: options.predictions.clone(),
        variations,
        regulatory: None,
    };
    crate::build(&inputs, &cache_dir)?;

    let species_dir = target.parent().expect("the release directory has a parent");
    std::fs::create_dir_all(species_dir)
        .with_context(|| format!("Failed to create {}", species_dir.display()))?;
    std::fs::rename(&cache_dir, &target).with_context(|| {
        format!(
            "Failed to move {} to {}",
            cache_dir.display(),
            target.display()
        )
    })?;
    std::fs::remove_dir_all(&work_dir)
        .with_context(|| format!("Failed to remove {}", work_dir.display()))?;
    let _ = std::fs::remove_dir(dir_cache.join("tmp"));
    info!(
        "Installed {}: run vep with --dir_cache {} --species {species} --cache_version {release} --assembly {assembly}",
        target.display(),
        dir_cache.display()
    );
    Ok(target)
}

/// Places `<table>.txt.gz`, else `<table>.txt`, from a MySQL directory in
/// `dest_dir`; `None` when the directory lists neither.
fn fetch_table(directory: &Directory<'_>, table: &str, dest_dir: &Path) -> Result<Option<PathBuf>> {
    for name in [format!("{table}.txt.gz"), format!("{table}.txt")] {
        if directory.listing.contains(&name) {
            return directory.fetch(&name, dest_dir).map(Some);
        }
    }
    Ok(None)
}

/// Places the release's genome FASTA in `cache_dir` as an indexed BGZF file:
/// the `dna_index/` file with the `.fai` and `.gzi` published beside it, else
/// the `dna/` file recompressed and indexed here. Returns the FASTA's path.
fn install_genome_fasta(
    site: &Site,
    species: &str,
    assembly: &str,
    release: u32,
    inputs_dir: &Path,
    cache_dir: &Path,
) -> Result<PathBuf> {
    let [indexed_dir, plain_dir] = download::genome_fasta_dirs(release, assembly, species);
    if let Some(indexed) = directory(site, &indexed_dir)? {
        if let Some(name) = genome_fasta_name(&indexed.listing) {
            let fai = format!("{name}.fai");
            let gzi = format!("{name}.gzi");
            if indexed.listing.contains(&fai) && indexed.listing.contains(&gzi) {
                // Kept with the other inputs, so a build that fails after the
                // download does not fetch the genome again.
                for file in [&name, &fai, &gzi] {
                    let staged = indexed.fetch(file, inputs_dir)?;
                    link_or_copy(&staged, &cache_dir.join(file))?;
                }
                let fasta = cache_dir.join(&name);
                vep_fasta::IndexedFasta::from_path(&fasta).with_context(|| {
                    format!("The indexed FASTA {} does not open", fasta.display())
                })?;
                return Ok(fasta);
            }
        }
    }
    let plain = required_directory(site, &plain_dir)?;
    let Some(name) = genome_fasta_name(&plain.listing) else {
        bail!("{plain_dir}/CHECKSUMS lists no primary assembly or toplevel genome FASTA");
    };
    let source = plain.fetch(&name, inputs_dir)?;
    let target = cache_dir.join(&name);
    fasta_index::install_fasta(&source, &target)?;
    Ok(target)
}

/// Places `source` at `dest` as a hard link, or a copy where the file system
/// allows no link, replacing anything at `dest`.
fn link_or_copy(source: &Path, dest: &Path) -> Result<()> {
    if dest.exists() {
        std::fs::remove_file(dest)
            .with_context(|| format!("Failed to remove {}", dest.display()))?;
    }
    if std::fs::hard_link(source, dest).is_err() {
        std::fs::copy(source, dest).with_context(|| {
            format!("Failed to copy {} to {}", source.display(), dest.display())
        })?;
    }
    Ok(())
}

/// Places the variation VCFs of the chromosomes (`--chromosomes`, else
/// [`download::default_chromosomes`] of the species and `gff3`) and the
/// genome-wide files the site publishes, each with its `.csi`, under
/// `<work_dir>/variation_vcf/`, and the variation database's ClinVar tables
/// under `<work_dir>/variation_dumps/`; a ClinVar table already fetched into
/// `core_dumps_dir` is linked from there, not fetched again.
fn fetch_variations(
    options: &InstallOptions,
    variation_db: Option<&Directory<'_>>,
    work_dir: &Path,
    core_dumps_dir: &Path,
    gff3: &Path,
) -> Result<VariationBuild> {
    let species = options.species.as_str();
    let assembly = options.assembly.as_str();
    let release = options.release;
    let vcf_dir = work_dir.join("variation_vcf");
    let vcfs = required_directory(
        &options.site,
        &download::variation_vcf_dir(release, assembly, species),
    )?;
    let chromosomes = match &options.chromosomes {
        Some(list) => list.clone(),
        None => download::default_chromosomes(species, Some(gff3))?,
    };
    for chromosome in &chromosomes {
        let name = download::variation_vcf_name(species, chromosome);
        if !vcfs.listing.contains(&name) {
            bail!(
                "{name} is not published for this release and assembly; name the chromosomes with --chromosomes"
            );
        }
        vcfs.fetch(&name, &vcf_dir)?;
        fetch_index(&vcfs, &name, &vcf_dir)?;
    }
    for name in download::genome_wide_names(species) {
        if vcfs.listing.contains(&name) {
            vcfs.fetch(&name, &vcf_dir)?;
            fetch_index(&vcfs, &name, &vcf_dir)?;
        }
    }
    let inputs = download::resolve_variation_inputs(
        &Some(vcf_dir.display().to_string()),
        release,
        assembly,
        species,
        Some(gff3),
        None,
    )?;

    let dumps = match variation_db {
        Some(variation) => {
            let dumps_dir = work_dir.join("variation_dumps");
            std::fs::create_dir_all(&dumps_dir)
                .with_context(|| format!("Failed to create {}", dumps_dir.display()))?;
            for table in CLINVAR_TABLES {
                // Only a table fetched from the variation database is shared:
                // the core database's `seq_region` and `attrib_type` in
                // `core_dumps_dir` are other tables of the same names.
                let fetched = if VARIATION_TABLES.contains(&table) {
                    table_path(core_dumps_dir, table)
                } else {
                    None
                };
                match fetched {
                    Some(fetched) => {
                        let name = fetched.file_name().expect("a table dump has a file name");
                        let linked = dumps_dir.join(name);
                        info!("Linking {} to {}", fetched.display(), linked.display());
                        link_or_copy(&fetched, &linked)?;
                    }
                    None => {
                        fetch_table(variation, table, &dumps_dir)?.ok_or_else(|| {
                            anyhow::anyhow!(
                                "variation table {table} is not published for this release"
                            )
                        })?;
                    }
                }
            }
            Some(dumps_dir)
        }
        None => {
            warn!("No variation database is published: the ClinVar allele and somatic classification columns are empty");
            None
        }
    };
    Ok(VariationBuild { inputs, dumps })
}

/// The variation tables the ClinVar columns are read from.
const CLINVAR_TABLES: [&str; 6] = [
    "source",
    "seq_region",
    "attrib_type",
    "phenotype",
    "phenotype_feature",
    "phenotype_feature_attrib",
];

/// Places the `.csi` of a VCF when the directory lists one.
fn fetch_index(directory: &Directory<'_>, name: &str, dest_dir: &Path) -> Result<()> {
    let index = format!("{name}.csi");
    if directory.listing.contains(&index) {
        directory.fetch(&index, dest_dir)?;
    } else {
        info!("No index is published for {name}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(dir: &Path, name: &str, bytes: &[u8]) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }

    /// The checksum and block count `sum` prints for the same bytes: the empty
    /// file, one byte, two bytes (the rotate before the add), a block boundary
    /// and one byte past it.
    #[test]
    fn bsd_sum_matches_the_sum_command() {
        let dir = tempfile::tempdir().unwrap();
        let cases: [(&[u8], u16, u64); 6] = [
            (b"", 0, 0),
            (b"a", 97, 1),
            (b"ab", 32914, 1),
            (b"The quick brown fox jumps over the lazy dog\n", 25281, 1),
            (&[0u8; 1024], 0, 1),
            (&[0u8; 1025], 0, 2),
        ];
        for (i, (bytes, sum, blocks)) in cases.iter().enumerate() {
            let path = write(dir.path(), &format!("case{i}"), bytes);
            assert_eq!(
                bsd_sum(&path).unwrap(),
                Checksum {
                    sum: *sum,
                    blocks: *blocks
                },
                "case {i}"
            );
        }
        assert_eq!(
            Checksum {
                sum: 405,
                blocks: 7
            }
            .to_string(),
            "00405 7"
        );
    }

    /// A listing keeps the leading-zero checksum as a number, splits on tabs
    /// as on spaces, and passes over lines that are not three fields.
    #[test]
    fn listing_parses_checksums_lines() {
        let listing = Listing::parse(
            "04440 6188 Homo_sapiens.GRCh38.116.abinitio.gff3.gz\n\
             02333 105838 Homo_sapiens.GRCh38.116.gff3.gz\n\
             42797 849535\tHomo_sapiens.GRCh37.dna.primary_assembly.fa.gz\n\
             a line that is not a checksum\n\
             12345\n\n",
        );
        assert_eq!(
            listing.get("Homo_sapiens.GRCh38.116.gff3.gz"),
            Some(Checksum {
                sum: 2333,
                blocks: 105838
            })
        );
        assert_eq!(
            listing.get("Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"),
            Some(Checksum {
                sum: 42797,
                blocks: 849535
            })
        );
        assert!(!listing.contains("12345"));
        assert_eq!(listing.names().len(), 3);
    }

    /// The primary assembly is preferred to the toplevel file and the masked
    /// files are never chosen; a listing without a genome file names none.
    #[test]
    fn genome_fasta_prefers_the_primary_assembly() {
        let listing = Listing::parse(
            "1 1 Homo_sapiens.GRCh38.dna.toplevel.fa.gz\n\
             1 1 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz\n\
             1 1 Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz\n\
             1 1 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz\n",
        );
        assert_eq!(
            genome_fasta_name(&listing).as_deref(),
            Some("Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
        );
        let toplevel_only = Listing::parse(
            "1 1 Homo_sapiens.GRCh38.dna.toplevel.fa.gz\n\
             1 1 Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai\n\
             1 1 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz\n",
        );
        assert_eq!(
            genome_fasta_name(&toplevel_only).as_deref(),
            Some("Homo_sapiens.GRCh38.dna.toplevel.fa.gz")
        );
        assert_eq!(genome_fasta_name(&Listing::parse("1 1 README\n")), None);
    }

    /// A directory site places a listed file by copy and verifies it, keeps a
    /// placed copy that verifies, replaces one that does not, and refuses a
    /// mismatch naming the file and both checksums.
    #[test]
    fn directory_fetch_verifies_against_the_listing() {
        let root = tempfile::tempdir().unwrap();
        let data_dir = root.path().join("pub/release-116/gff3/homo_sapiens");
        std::fs::create_dir_all(&data_dir).unwrap();
        write(&data_dir, "good.gz", b"ab");
        write(&data_dir, "bad.gz", b"ab");
        write(&data_dir, "CHECKSUMS", b"32914 1 good.gz\n00001 1 bad.gz\n");
        let site = Site::parse(root.path().to_str().unwrap()).unwrap();
        let dir = required_directory(&site, "pub/release-116/gff3/homo_sapiens").unwrap();
        let dest = tempfile::tempdir().unwrap();

        let placed = dir.fetch("good.gz", dest.path()).unwrap();
        assert_eq!(std::fs::read(&placed).unwrap(), b"ab");
        std::fs::write(&placed, b"ab").unwrap();
        assert_eq!(dir.fetch("good.gz", dest.path()).unwrap(), placed);
        std::fs::write(&placed, b"stale").unwrap();
        dir.fetch("good.gz", dest.path()).unwrap();
        assert_eq!(std::fs::read(&placed).unwrap(), b"ab");

        let error = dir.fetch("bad.gz", dest.path()).unwrap_err().to_string();
        assert!(
            error.contains("bad.gz") && error.contains("32914 1") && error.contains("00001 1"),
            "{error}"
        );
        assert!(!dest.path().join("bad.gz").exists());
        let error = dir.fetch("absent.gz", dest.path()).unwrap_err().to_string();
        assert!(error.contains("absent.gz is not listed"), "{error}");

        assert!(directory(&site, "pub/release-116/gtf/homo_sapiens")
            .unwrap()
            .is_none());
        let error = required_directory(&site, "pub/release-116/gtf/homo_sapiens")
            .unwrap_err()
            .to_string();
        assert!(error.contains("gtf/homo_sapiens/CHECKSUMS"), "{error}");
    }

    /// A URL keeps its scheme and loses a trailing slash; anything else must be
    /// a directory.
    #[test]
    fn site_is_a_url_or_a_directory() {
        assert!(matches!(
            Site::parse("https://ftp.ensembl.org/").unwrap(),
            Site::Url(url) if url == "https://ftp.ensembl.org"
        ));
        let dir = tempfile::tempdir().unwrap();
        assert!(matches!(
            Site::parse(dir.path().to_str().unwrap()).unwrap(),
            Site::Dir(_)
        ));
        let error = Site::parse("/no/such/directory").unwrap_err().to_string();
        assert!(error.contains("neither a URL nor a directory"), "{error}");
    }
}
