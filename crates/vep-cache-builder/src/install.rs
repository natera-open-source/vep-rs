// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One-command installation of a release (`--auto`): every input is fetched
//! from the Ensembl site, or a copy of its tree, and checked against the
//! `CHECKSUMS` file of its directory; the cache is built in Ensembl VEP's
//! layout, `<dir_cache>/<species>/<release>_<assembly>`, with the indexed
//! genome FASTA inside it, where `vep` finds both through `--dir_cache`.
//!
//! The genome is found first ([`locate`]): the release's species list says
//! which tree a species is published in and names its core database, on the
//! Ensembl site, the Ensembl Genomes site or the one site `--site` names. A
//! release tree with no species list that lists the genome's GFF3 is the
//! genome's, and no later site is asked; a site of several that cannot be
//! read is passed over with a warning. Every input is then looked up in its
//! directory's listing, so a genome whose tree lacks one is refused before
//! anything is fetched, as a genome of a collection (one core database for
//! all its genomes) and an assembly of the rapid release (no table dumps)
//! are. A file a tree may lack (the release's own `CHECKSUMS`, an index
//! page, the variation database, `dna_index/`) reads as absent on a server
//! that answers 403 for it as on one that answers 404.
//!
//! Beside the transcripts the installation brings the parts of Ensembl VEP's
//! own cache that a file of the release carries: the known variants with the
//! 1000 Genomes frequencies (the variation VCFs), the SIFT and PolyPhen
//! matrices (the variation database's table dumps) and the regulatory and
//! motif features (the regulatory annotation GFF3 files). A part the site
//! publishes no file of for the genome is left out, and the installation
//! ends on one line a part left out, naming it and the reason ([`close`]),
//! then on one that names the sequences of the cache it read no variation
//! record on, then on one for each reason a cell of a part it brought holds
//! no value.
//! The known variants take their PubMed identifiers, their synonyms and
//! their failed flag from tables of the variation database, and the failed
//! variants the VCFs leave out from the failed-variants GVF, wherever the
//! release publishes each (`VariantCellSources`); the tables are fetched
//! whole whatever `--chromosomes` names. Where the germline variants come
//! in one genome-wide VCF, a sequence of the cache that no variation VCF
//! holds a record on takes its known variants from the variation
//! database's `variation_feature` table, which is fetched only when there
//! is such a sequence (`SequenceFeatureTables`).
//!
//! The downloads land under `<dir_cache>/tmp/<species>_<release>_<assembly>/`
//! and are removed once the cache is in place; a downloaded file whose checksum
//! still matches is reused by the next run, so an interrupted installation
//! resumes where it stopped. A release directory that already exists is never
//! written into.

use std::cell::OnceCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Component, Path, PathBuf};

use anyhow::{bail, Context, Result};
use tracing::{info, warn};

use crate::core_dumps::{table_path, CORE_TABLES, VARIATION_TABLES};
use crate::count::counted;
use crate::download::{self, Layout};
use crate::left_out::any_of;
use crate::{
    fasta_index, predictions, regulatory, variant_cells, variation, BuildInputs, MysqlServer,
    PredictionSource, SequenceFeatures, VariationBuild,
};

/// Where the release's files are read from: the Ensembl site, the Ensembl
/// Genomes site or another server publishing the tree of either, or a local
/// directory holding such a tree.
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

    /// The sites an installation without `--site` looks a species up on, in
    /// order: Ensembl's, then Ensembl Genomes'.
    pub fn defaults() -> Vec<Self> {
        vec![
            Self::Url(download::ENSEMBL_SITE.to_string()),
            Self::Url(download::ENSEMBL_GENOMES_SITE.to_string()),
        ]
    }

    /// The address of `relative` on the site: its URL, or its path.
    pub(crate) fn address(&self, relative: &str) -> String {
        match self {
            Self::Url(root) => format!("{root}/{relative}"),
            Self::Dir(root) => root.join(relative).display().to_string(),
        }
    }

    /// The bytes of the file at `relative`, or `None` when the site has no
    /// such file.
    fn read_bytes(&self, relative: &str) -> Result<Option<Vec<u8>>> {
        match self {
            Self::Url(root) => download::fetch_bytes(&format!("{root}/{relative}")),
            Self::Dir(root) => {
                let path = root.join(relative);
                if !path.is_file() {
                    return Ok(None);
                }
                Ok(Some(std::fs::read(&path).with_context(|| {
                    format!("Failed to read {}", path.display())
                })?))
            }
        }
    }

    /// The names of the entries of the directory at `relative`, read from a
    /// server's index page of it; `None` when the site has no such directory
    /// or, on a server, serves no index page of it.
    fn list(&self, relative: &str) -> Result<Option<Vec<String>>> {
        match self {
            Self::Url(root) => download::fetch_index(&format!("{root}/{relative}/")),
            Self::Dir(root) => {
                let path = root.join(relative);
                if !path.is_dir() {
                    return Ok(None);
                }
                let mut entries = Vec::new();
                for entry in std::fs::read_dir(&path)
                    .with_context(|| format!("Failed to read {}", path.display()))?
                {
                    let name = entry
                        .with_context(|| format!("Failed to read {}", path.display()))?
                        .file_name();
                    if let Some(name) = name.to_str() {
                        entries.push(name.to_string());
                    }
                }
                Ok(Some(entries))
            }
        }
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
    pub(crate) fn place(&self, relative: &str, dest: &Path) -> Result<()> {
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
    /// The sites the species is looked up on, in order: [`Site::defaults`],
    /// or the one `--site` names.
    pub sites: Vec<Site>,
    /// The chromosomes whose variation records are read: under
    /// `homo_sapiens` the per-chromosome files fetched,
    /// [`download::default_chromosomes`] when `None`; under any other
    /// species the sequences read from its one genome-wide file, every
    /// sequence it holds records on when `None`.
    pub chromosomes: Option<Vec<String>>,
    /// `--codon-table`: the codon table of each sequence region the core
    /// dumps state none for.
    pub codon_tables: HashMap<String, u8>,
    /// The server the SIFT and PolyPhen matrices come from under
    /// `--include-predictions`; `None` reads them from the variation
    /// database's table dumps, where the release publishes them.
    pub prediction_server: Option<MysqlServer>,
}

/// A part of Ensembl VEP's own cache of a genome that the installed cache
/// lacks, with the reason: the site publishes no file of it for the genome,
/// or the installation found no file that carries it, and then the route
/// that does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeftOut {
    pub part: String,
    pub reason: String,
}

impl LeftOut {
    fn new(part: &str, reason: String) -> Self {
        Self {
            part: part.to_string(),
            reason,
        }
    }
}

/// A part of an installation looked up on the site: its files, or the reason
/// the site has none for the genome.
enum Part<T> {
    Published(T),
    Unpublished(String),
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
    /// Parses the file's text; a line without the three fields, whose first
    /// two are not numbers, or whose name is not one bare file name, is passed
    /// over.
    pub fn parse(text: &str) -> Self {
        let entries = text
            .lines()
            .filter_map(|line| {
                let mut fields = line.split_whitespace();
                let sum = fields.next()?.parse::<u16>().ok()?;
                let blocks = fields.next()?.parse::<u64>().ok()?;
                let name = fields.next()?;
                // A name is joined onto the destination directory, so an
                // absolute path, `..` or a separator would place a file outside it.
                let mut parts = Path::new(name).components();
                match (parts.next(), parts.next()) {
                    (Some(Component::Normal(part)), None) if part == name => {}
                    _ => return None,
                }
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
        sum = add_to_sum(sum, &buffer[..read]);
        bytes += read as u64;
    }
    Ok(Checksum {
        sum,
        blocks: bytes.div_ceil(1024),
    })
}

/// The BSD checksum `sum` carried over `bytes`.
fn add_to_sum(sum: u16, bytes: &[u8]) -> u16 {
    bytes.iter().fold(sum, |sum, &byte| {
        sum.rotate_right(1).wrapping_add(u16::from(byte))
    })
}

/// The BSD `sum` of a file held in memory.
fn bsd_sum_of(bytes: &[u8]) -> Checksum {
    Checksum {
        sum: add_to_sum(0, bytes),
        blocks: (bytes.len() as u64).div_ceil(1024),
    }
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

/// As [`directory`], for a directory a tree may lack: a server's 403 for its
/// `CHECKSUMS` reads as a 404 does, since a server may answer so for a file
/// it does not hold.
fn optional_directory<'a>(site: &'a Site, relative: &str) -> Result<Option<Directory<'a>>> {
    match directory(site, relative) {
        Err(error) if download::is_forbidden(&error) => {
            info!(
                "{} answers 403 for {relative}/CHECKSUMS, read as a directory it does not publish",
                site.address("")
            );
            Ok(None)
        }
        read => read,
    }
}

/// The work directory of the installation:
/// `<dir_cache>/tmp/<species>_<release>_<assembly>`.
pub fn work_dir(options: &InstallOptions) -> download::WorkDir {
    download::WorkDir::of_release(
        &options.dir_cache,
        &options.species,
        options.release,
        &options.assembly,
    )
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

/// A genome found on a site: the tree its files sit in and its core database,
/// `None` in a tree with no table dumps.
#[derive(Debug, Clone)]
pub struct Located {
    pub site: Site,
    pub layout: Layout,
    pub core_db: Option<String>,
}

/// The refusal of a genome whose tree lacks an input the build reads: one
/// line naming the genome, what its tree lacks and the routes that build its
/// cache.
fn lacking(site: &Site, species: &str, assembly: &str, what: &str) -> anyhow::Error {
    anyhow::Error::new(download::NotPublished(format!(
        "{species} {assembly} cannot be installed from {}: {what}; convert Ensembl VEP's cache of it instead (docs/cache-conversion.md), or build from local copies of its files with --gff3, --fasta, --genome-fasta and --codon-table",
        site.address("")
    )))
}

/// A lookup's answer that the genome is not on the site, as apart from a
/// failed read.
fn not_published(message: String) -> anyhow::Error {
    anyhow::Error::new(download::NotPublished(message))
}

impl Located {
    /// Refuses a genome whose tree has no table dumps of its own to read: a
    /// genome of a collection, whose core database holds every genome of the
    /// collection in the same tables, which the table reader does not split
    /// by genome, and an assembly of the rapid release.
    pub fn buildable(&self, species: &str, assembly: &str, release: u32) -> Result<()> {
        let what = match &self.layout {
            Layout::Collection { collection, .. } => format!(
                "it has no core database of its own ({} holds every genome of {collection}, and a shared database is not split by genome)",
                self.core_db.as_deref().unwrap_or(collection)
            ),
            Layout::RapidRelease { .. } => format!(
                "{}/ publishes no core table dumps",
                self.layout.root(release)
            ),
            Layout::Release | Layout::Grch37 | Layout::Division { .. } => return Ok(()),
        };
        Err(lacking(&self.site, species, assembly, &what))
    }
}

/// The text of the species list `name` of a release, when the release
/// directory's `CHECKSUMS` lists it; the file is refused when its BSD `sum`
/// differs from that line.
fn species_list(site: &Site, release: u32, name: &str) -> Result<Option<String>> {
    let root = Layout::Release.root(release);
    let Some(directory) = optional_directory(site, &root)? else {
        return Ok(None);
    };
    let Some(expected) = directory.listing.get(name) else {
        return Ok(None);
    };
    let relative = format!("{root}/{name}");
    let Some(bytes) = site.read_bytes(&relative)? else {
        bail!("{relative} is listed in {root}/CHECKSUMS and is not on the site");
    };
    let actual = bsd_sum_of(&bytes);
    if actual != expected {
        bail!(
            "checksum for {} ({actual}) does not match the CHECKSUMS line for {name} ({expected})",
            site.address(&relative)
        );
    }
    Ok(Some(String::from_utf8_lossy(&bytes).into_owned()))
}

/// The core database of `species` in the `mysql/` directory of a release
/// tree that publishes no species list, found by listing the directory.
fn listed_core_database(
    site: &Site,
    layout: &Layout,
    species: &str,
    assembly: &str,
    release: u32,
) -> Result<String> {
    let mysql = format!("{}/mysql", layout.root(release));
    let entries = site.list(&mysql)?;
    core_database_among(
        &mysql,
        entries.as_deref(),
        species,
        assembly,
        release,
        &|database| Ok(optional_directory(site, &format!("{mysql}/{database}"))?.is_some()),
    )
}

/// The core database of `species` among the `entries` of the directory
/// `mysql`: the one named for the species and the release. A server that
/// serves no index page lists nothing (`None`); there the database named for
/// the assembly's number ([`download::numbered_core_database`]) is taken when
/// `holds` finds its directory on the site.
fn core_database_among(
    mysql: &str,
    entries: Option<&[String]>,
    species: &str,
    assembly: &str,
    release: u32,
    holds: &dyn Fn(&str) -> Result<bool>,
) -> Result<String> {
    let Some(entries) = entries else {
        if let Some(database) = download::numbered_core_database(species, assembly, release) {
            if holds(&database)? {
                return Ok(database);
            }
        }
        return Err(not_published(format!(
            "{mysql}/ is not on the site: is the release published?"
        )));
    };
    match download::core_databases(entries, species, release).as_slice() {
        [database] => Ok(database.to_string()),
        [] => Err(not_published(format!(
            "{mysql}/ holds no core database of {species} for release {release}"
        ))),
        several => bail!(
            "{mysql}/ holds more than one core database of {species}: {}",
            several.join(", ")
        ),
    }
}

/// What the species lists a site publishes for a release say of a species.
enum Listed {
    /// A list names the species: its tree and its core database.
    Named(Located),
    /// The addresses of the lists read, none of which names the species.
    Unnamed(Vec<String>),
    /// The site publishes no species list for the release.
    NoList,
}

/// Reads the species lists of `site` for the release
/// ([`download::SPECIES_LISTS`]) for the row of `species`.
fn listed(site: &Site, species: &str, release: u32) -> Result<Listed> {
    let mut read = Vec::new();
    for name in download::SPECIES_LISTS {
        let Some(text) = species_list(site, release, name)? else {
            continue;
        };
        let list = site.address(&format!("{}/{name}", Layout::Release.root(release)));
        let entry = download::species_entry(&text, species)
            .with_context(|| format!("Failed to read {list}"))?;
        match entry {
            Some(entry) => {
                return Ok(Listed::Named(Located {
                    site: site.clone(),
                    layout: entry.layout(species),
                    core_db: Some(entry.core_db),
                }))
            }
            None => read.push(list),
        }
    }
    Ok(if read.is_empty() {
        Listed::NoList
    } else {
        Listed::Unnamed(read)
    })
}

/// The genome in the release tree of a site that publishes no species list:
/// found when the tree's GFF3 directory lists the genome's GFF3, whose name
/// carries the species, the assembly and the release, and its `mysql/`
/// directory holds a core database of the species. A release tree holding
/// another assembly of the species under the same release number is not the
/// genome's.
fn held_in_release_tree(
    site: &Site,
    species: &str,
    assembly: &str,
    release: u32,
) -> Result<Option<Located>> {
    let layout = Layout::Release;
    let (dir, name) = layout.gff3_location(release, assembly, species);
    let lists_the_gff3 =
        optional_directory(site, &dir)?.is_some_and(|directory| directory.listing.contains(&name));
    if !lists_the_gff3 {
        return Ok(None);
    }
    match listed_core_database(site, &layout, species, assembly, release) {
        Ok(core_db) => Ok(Some(Located {
            site: site.clone(),
            layout,
            core_db: Some(core_db),
        })),
        Err(error) if download::is_not_published(&error) => Ok(None),
        Err(error) => Err(error),
    }
}

/// The sites of a lookup that could not be read, each with the error of the
/// read that failed.
#[derive(Default)]
struct Unread {
    sites: Vec<(usize, String)>,
}

impl Unread {
    /// Records a failed read of the site at `nth`, with a warning; the
    /// lookup goes on to the sites left.
    fn pass_over(&mut self, nth: usize, site: &Site, error: &anyhow::Error) {
        let address = site.address("");
        warn!("{address} could not be read and is passed over: {error:#}");
        self.sites
            .push((nth, format!("{address} could not be read ({error:#})")));
    }

    fn holds(&self, nth: usize) -> bool {
        self.sites.iter().any(|(site, _)| *site == nth)
    }

    /// `error` as it stands when every site was read; else a failed lookup,
    /// no [`download::NotPublished`] answer, that ends naming each site that
    /// could not be read.
    fn named_in(&self, error: anyhow::Error) -> anyhow::Error {
        if self.sites.is_empty() {
            return error;
        }
        let unread: Vec<&str> = self.sites.iter().map(|(_, note)| note.as_str()).collect();
        anyhow::anyhow!("{error:#}; {}", unread.join("; "))
    }
}

/// Finds the genome on the first of `sites` that publishes it.
///
/// GRCh37 is the first site's tree of that name. Any other genome is looked
/// up in each site's species lists of the release
/// ([`download::SPECIES_LISTS`]), whose row gives its layout and its core
/// database. Where the first of several sites publishes no list, its release
/// tree is the genome's when it lists the genome's GFF3 and holds its core
/// database, and no later site is asked. A species no list names, called for
/// an assembly accession, is looked for in the rapid release, where the
/// latest gene-set directory is taken. A species that a list was read for and
/// none names is refused; where every site is read and none publishes a list
/// for the release, the genome is taken to sit in the release tree of the
/// first site and its core database is found in that tree's `mysql/`
/// directory, or by the assembly's number on a server that lists no
/// directory.
///
/// A site of several that cannot be read is passed over with a warning, and
/// a genome found on no site is then refused naming it; one site alone stops
/// the lookup on a failed read.
pub fn locate(sites: &[Site], species: &str, assembly: &str, release: u32) -> Result<Located> {
    let first = sites.first().context("no site to read the release from")?;
    if assembly == "GRCh37" {
        let layout = Layout::Grch37;
        let core_db = listed_core_database(first, &layout, species, assembly, release)?;
        return Ok(Located {
            site: first.clone(),
            layout,
            core_db: Some(core_db),
        });
    }
    let several = sites.len() > 1;
    let mut lists_read = Vec::new();
    let mut unread = Unread::default();
    for (nth, site) in sites.iter().enumerate() {
        let found = listed(site, species, release).and_then(|listed| match listed {
            Listed::Named(located) => Ok(Some(located)),
            Listed::Unnamed(lists) => {
                lists_read.extend(lists);
                Ok(None)
            }
            Listed::NoList if nth == 0 && several => {
                held_in_release_tree(site, species, assembly, release)
            }
            Listed::NoList => Ok(None),
        });
        match found {
            Ok(Some(located)) => return Ok(located),
            Ok(None) => {}
            Err(error) if several => unread.pass_over(nth, site, &error),
            Err(error) => return Err(error),
        }
    }
    if let Some((organism, accession)) = download::rapid_release_name(species) {
        for (nth, site) in sites.iter().enumerate() {
            if unread.holds(nth) {
                continue;
            }
            let gene_sets =
                format!("pub/rapid-release/species/{organism}/{accession}/ensembl/geneset");
            let dated = match site.list(&gene_sets) {
                Ok(dated) => dated,
                Err(error) if several => {
                    unread.pass_over(nth, site, &error);
                    continue;
                }
                Err(error) => return Err(error),
            };
            let Some(latest) = dated.and_then(|dated| dated.into_iter().max()) else {
                continue;
            };
            return Ok(Located {
                site: site.clone(),
                layout: Layout::RapidRelease {
                    organism,
                    accession,
                    gene_set: latest,
                },
                core_db: None,
            });
        }
    }
    if !lists_read.is_empty() {
        return Err(unread.named_in(not_published(format!(
            "{species} is not in the species list of release {release} ({}); a species of Ensembl Genomes takes the Ensembl Genomes release number in --release",
            lists_read.join(", ")
        ))));
    }
    if unread.holds(0) {
        return Err(unread.named_in(anyhow::anyhow!(
            "{species} {assembly} was not found for release {release}"
        )));
    }
    let layout = Layout::Release;
    let core_db = listed_core_database(first, &layout, species, assembly, release)
        .map_err(|error| unread.named_in(error))?;
    // The species list of a site that could not be read may name the species,
    // and a list outranks a release tree that does not list the genome's GFF3.
    if !unread.sites.is_empty() {
        return Err(unread.named_in(anyhow::anyhow!(
            "{species} {assembly} was not found for release {release}"
        )));
    }
    Ok(Located {
        site: first.clone(),
        layout,
        core_db: Some(core_db),
    })
}

/// The genome of a build named one input at a time, located on `sites` the
/// first time a file has to be downloaded, so a build from local files alone
/// reads no site.
pub struct Published<'a> {
    sites: Vec<Site>,
    species: &'a str,
    assembly: &'a str,
    release: u32,
    located: OnceCell<Located>,
}

impl<'a> Published<'a> {
    pub fn new(sites: Vec<Site>, species: &'a str, assembly: &'a str, release: u32) -> Self {
        Self {
            sites,
            species,
            assembly,
            release,
            located: OnceCell::new(),
        }
    }

    fn located(&self) -> Result<&Located> {
        if let Some(located) = self.located.get() {
            return Ok(located);
        }
        let located = locate(&self.sites, self.species, self.assembly, self.release)?;
        located.buildable(self.species, self.assembly, self.release)?;
        Ok(self.located.get_or_init(|| located))
    }

    /// The genome's GFF3: the site that publishes it and the file's path
    /// under that site's root.
    pub fn gff3(&self) -> Result<(&Site, String)> {
        let located = self.located()?;
        let (dir, name) = located
            .layout
            .gff3_location(self.release, self.assembly, self.species);
        Ok((&located.site, format!("{dir}/{name}")))
    }

    /// The genome's protein FASTA: the site that publishes it and the file's
    /// path under that site's root.
    pub fn protein_fasta(&self) -> Result<(&Site, String)> {
        let located = self.located()?;
        let (dir, name) =
            located
                .layout
                .protein_fasta_location(self.release, self.assembly, self.species);
        Ok((&located.site, format!("{dir}/{name}")))
    }

    /// The URLs of the genome's core and variation MySQL directories.
    pub fn database_urls(&self) -> Result<download::DatabaseUrls> {
        let located = self.located()?;
        let (core, variation) = located.databases(self.species, self.release)?;
        let url = |dir: String| format!("{}/", located.site.address(&dir));
        Ok(download::DatabaseUrls {
            core: url(core),
            variation: variation.map(url),
        })
    }
}

impl Located {
    /// The MySQL directories of the genome's core database and of the
    /// variation database named beside it, which a site may not publish.
    fn databases(&self, species: &str, release: u32) -> Result<(String, Option<String>)> {
        let directories = self.core_db.as_deref().and_then(|core_db| {
            let core = self.layout.mysql_dir(release, core_db)?;
            let variation = download::variation_database(core_db, species)
                .and_then(|database| self.layout.mysql_dir(release, &database));
            Some((core, variation))
        });
        directories.with_context(|| format!("no core database is published for {species}"))
    }

    /// The MySQL directory of the funcgen database named beside the genome's
    /// core database, which a site may not publish. `None` where the core
    /// database is not named for the species and in a tree with no table
    /// dumps.
    fn funcgen_dir(&self, species: &str, release: u32) -> Option<String> {
        let database = download::funcgen_database(self.core_db.as_deref()?, species)?;
        self.layout.mysql_dir(release, &database)
    }
}

/// One input of an installation: a file its directory's listing names.
struct Input<'a> {
    directory: Directory<'a>,
    name: String,
}

impl Input<'_> {
    fn fetch(&self, dest_dir: &Path) -> Result<PathBuf> {
        self.directory.fetch(&self.name, dest_dir)
    }
}

/// The genome FASTA of an installation: the `dna_index/` file with the `.fai`
/// and `.gzi` published beside it, else the `dna/` file.
struct GenomeFasta<'a> {
    file: Input<'a>,
    indexed: bool,
}

/// Every input of an installation, each found in its directory's listing
/// before anything is fetched.
struct Inputs<'a> {
    gff3: Input<'a>,
    gtf: Input<'a>,
    protein_fasta: Input<'a>,
    core: Directory<'a>,
    variation: Option<Directory<'a>>,
    genome: GenomeFasta<'a>,
}

impl<'a> Inputs<'a> {
    /// Opens the directory of each input of the located genome and looks the
    /// input up in its listing; a directory or a file the tree lacks refuses
    /// the genome, naming it.
    fn find(located: &'a Located, species: &str, assembly: &str, release: u32) -> Result<Self> {
        let site = &located.site;
        let layout = &located.layout;
        let lacking = |what: String| lacking(site, species, assembly, &what);
        let open = |relative: &str| -> Result<Directory<'a>> {
            directory(site, relative)?.ok_or_else(|| {
                lacking(format!(
                    "{relative}/CHECKSUMS is not on the site (is the release published?)"
                ))
            })
        };
        let input = |(relative, name): (String, String)| -> Result<Input<'a>> {
            let directory = open(&relative)?;
            if !directory.listing.contains(&name) {
                return Err(lacking(format!(
                    "{name} is not listed in {relative}/CHECKSUMS"
                )));
            }
            Ok(Input { directory, name })
        };

        let gff3 = input(layout.gff3_location(release, assembly, species))?;
        let gtf = input(layout.gtf_location(release, assembly, species))?;
        let protein_fasta = input(layout.protein_fasta_location(release, assembly, species))?;

        let (core_dir, variation_dir) = located.databases(species, release)?;
        let core = open(&core_dir)?;
        for table in CORE_TABLES {
            if table_name(&core, table).is_none() {
                return Err(lacking(format!(
                    "core table {table} is not listed in {core_dir}/CHECKSUMS"
                )));
            }
        }
        let variation = match variation_dir {
            Some(relative) => optional_directory(site, &relative)?,
            None => None,
        };

        let dirs = layout.genome_fasta_dirs(release, species);
        let (plain_dir, indexed_dirs) = dirs.split_last().context("no genome FASTA directory")?;
        let mut genome = None;
        for relative in indexed_dirs {
            let Some(directory) = optional_directory(site, relative)? else {
                continue;
            };
            let Some(name) = genome_fasta_name(&directory.listing) else {
                continue;
            };
            let listing = &directory.listing;
            if listing.contains(&format!("{name}.fai")) && listing.contains(&format!("{name}.gzi"))
            {
                genome = Some(GenomeFasta {
                    file: Input { directory, name },
                    indexed: true,
                });
                break;
            }
        }
        let genome = match genome {
            Some(genome) => genome,
            None => {
                let directory = open(plain_dir)?;
                let Some(name) = genome_fasta_name(&directory.listing) else {
                    return Err(lacking(format!(
                        "{plain_dir}/CHECKSUMS lists no primary assembly or toplevel genome FASTA"
                    )));
                };
                GenomeFasta {
                    file: Input { directory, name },
                    indexed: false,
                }
            }
        };

        Ok(Self {
            gff3,
            gtf,
            protein_fasta,
            core,
            variation,
            genome,
        })
    }
}

/// Installs the release: finds the genome and every input in its tree, fetches
/// and verifies each, places the indexed genome FASTA in the cache directory
/// under construction, builds the cache there and moves it to [`release_dir`].
/// Returns that directory.
pub fn install(options: &InstallOptions) -> Result<PathBuf> {
    let species = options.species.as_str();
    let assembly = options.assembly.as_str();
    let release = options.release;
    let dir_cache = &options.dir_cache;
    let target = release_dir(options);
    if target.exists() {
        bail!(
            "{} exists: remove it to install the release again",
            target.display()
        );
    }
    let located = locate(&options.sites, species, assembly, release)?;
    located.buildable(species, assembly, release)?;
    let inputs = Inputs::find(&located, species, assembly, release)?;
    let mut unpublished = Unpublished::default();
    let vcfs = match VariationVcfs::find(&located, options)? {
        Part::Published(vcfs) => Some(vcfs),
        Part::Unpublished(reason) => {
            unpublished.variants = Some(reason);
            None
        }
    };
    let variant_cells = match &vcfs {
        Some(_) => Some(VariantCellSources::find(
            &located,
            options,
            inputs.variation.as_ref(),
        )?),
        None => None,
    };
    let sequence_features = SequenceFeatureTables::find(inputs.variation.as_ref());
    let matrices_in_the_dumps = match &options.prediction_server {
        Some(_) => false,
        None => match prediction_tables(inputs.variation.as_ref(), options) {
            Part::Published(()) => true,
            Part::Unpublished(reason) => {
                unpublished.matrices = Some(reason);
                false
            }
        },
    };
    let regulatory_files = match RegulatoryFiles::find(&located, options)? {
        Part::Published(files) => Some(files),
        Part::Unpublished(reason) => {
            unpublished.regulatory = Some(reason);
            None
        }
    };
    let regulatory_cells = match &regulatory_files {
        Some(files) => Some(RegulatoryCells::find(&located, options, files)?),
        None => None,
    };

    let work = work_dir(options);
    let work_dir = work.path();
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
    info!(
        "Reading {} with the tables of {}",
        located.site.address(&located.layout.root(release)),
        inputs.core.relative
    );

    let gff3 = inputs.gff3.fetch(&inputs_dir)?;
    let gtf = inputs.gtf.fetch(&inputs_dir)?;
    let protein_fasta = inputs.protein_fasta.fetch(&inputs_dir)?;

    for table in CORE_TABLES {
        fetch_table(&inputs.core, table, &core_dumps_dir)?.ok_or_else(|| {
            anyhow::anyhow!("core table {table} is not published for this release")
        })?;
    }
    match &inputs.variation {
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

    let genome_fasta = install_genome_fasta(&inputs.genome, &inputs_dir, &cache_dir)?;

    let variations = match &vcfs {
        Some(vcfs) => {
            let mut variations = fetch_variations(
                options,
                vcfs,
                inputs.variation.as_ref(),
                work_dir,
                &core_dumps_dir,
            )?;
            if let Some(sources) = &variant_cells {
                variations.cells = sources.fetch(inputs.variation.as_ref(), work_dir)?;
            }
            // One genome-wide file can hold no record on a sequence the
            // database holds known variants on; a file a chromosome is
            // published for each chromosome that has any.
            let one_genome_wide_file = variations
                .inputs
                .per_chromosome
                .iter()
                .all(|(_, vcf)| matches!(vcf, variation::GermlineVcf::Genome(_)));
            if one_genome_wide_file {
                match &sequence_features {
                    Ok(tables) => {
                        variations.features = Some(tables.fetched_when_asked(
                            &located.site,
                            work_dir,
                            options.chromosomes.clone(),
                        ));
                    }
                    Err(reason) => info!(
                        "No variation_feature table is read for a sequence the variation VCF \
                         holds no record on: {reason}"
                    ),
                }
            }
            Some(variations)
        }
        None => None,
    };

    let predictions = match &options.prediction_server {
        Some(server) => Some(PredictionSource::Server(server.clone())),
        None if matrices_in_the_dumps => {
            let variation = inputs
                .variation
                .as_ref()
                .expect("the tables were found in the variation database's listing");
            let dumps_dir = work_dir.join("variation_dumps");
            for table in predictions::PREDICTION_TABLES {
                fetch_table(variation, table, &dumps_dir)?.ok_or_else(|| {
                    anyhow::anyhow!("variation table {table} is not published for this release")
                })?;
            }
            Some(PredictionSource::Dumps(dumps_dir))
        }
        None => None,
    };
    let regulatory = match &regulatory_files {
        Some(files) => {
            let regulatory_dir = work_dir.join("regulatory");
            let mut inputs = files.fetch(&regulatory_dir)?;
            if let Some(cells) = &regulatory_cells {
                let dumps_dir = work_dir.join("funcgen_dumps");
                cells.fetch(files, &regulatory_dir, &dumps_dir, &mut inputs)?;
            }
            Some(inputs)
        }
        None => None,
    };

    let build_inputs = BuildInputs {
        species: species.to_string(),
        assembly: assembly.to_string(),
        release,
        gff3: Some(gff3),
        protein_fasta: Some(protein_fasta),
        gtf: Some(gtf),
        core_dumps: Some(core_dumps_dir),
        codon_tables: options.codon_tables.clone(),
        genome_fasta: Some(genome_fasta),
        predictions,
        variations,
        regulatory,
        work: None,
    };
    let built = crate::build(&build_inputs, &cache_dir)?;
    if let Some(cells) = &built.variant_cells {
        info!(
            "The known variants carry PubMed identifiers on {}, variant synonyms on {} and the \
             failed flag on {}, {} of them added from the failed-variants GVF",
            counted(cells.records_with_pubmed, "record"),
            counted(cells.records_with_synonyms, "record"),
            counted(cells.records_flagged_failed, "record"),
            cells.records_added_from_gvf
        );
    }
    match &built.sequence_features {
        Some(features) if features.in_question.is_empty() => info!(
            "The variation_feature table is not fetched: the variation VCF holds a record on \
             every sequence of the cache{}",
            if options.chromosomes.is_some() {
                " that --chromosomes names"
            } else {
                ""
            }
        ),
        Some(features) => info!(
            "The variation_feature table gave {} on {}, which the variation VCF holds no record \
             on",
            counted(
                features
                    .brought
                    .iter()
                    .map(|(_, records)| records)
                    .sum::<usize>(),
                "known variant"
            ),
            named_sequences(&features.in_question)
        ),
        None => {}
    }
    let left_out_of_the_build = &built.left_out;
    let record = built_record(&cache_dir)?;
    let without_variants = match &vcfs {
        Some(_) => sequences_without_known_variants(&cache_dir, &record)?,
        None => Vec::new(),
    };
    let left_out = left_out(
        options,
        &record,
        &unpublished,
        &without_variants,
        &Brought {
            vcfs: vcfs.as_ref(),
            variation_db: inputs.variation.as_ref(),
            variant_cells: variant_cells.as_ref(),
            regulatory_files: regulatory_files.as_ref(),
            regulatory_cells: regulatory_cells.as_ref(),
            regulatory_tree: built.regulatory.as_ref(),
        },
    );

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
    work.remove()?;
    info!(
        "Installed {}: run vep with --dir_cache {} --species {species} --cache_version {release} --assembly {assembly}",
        target.display(),
        dir_cache.display()
    );
    close(species, assembly, &left_out, &left_out_of_the_build.lines);
    Ok(target)
}

/// What an installed cache lacks of Ensembl VEP's own cache of the genome.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Missing {
    /// The parts left out: the site publishes no file of each for the genome.
    pub parts: Vec<LeftOut>,
    /// The sequences of the cache that hold no known variant: those the
    /// installation read no variation record on, in a variation VCF or in
    /// the `variation_feature` table, which it reads for a sequence the one
    /// genome-wide VCF of a species holds no record on.
    pub sequences: Vec<LeftOut>,
    /// The cells of a part that was brought which the cache holds no value
    /// for: the gnomAD frequencies of `homo_sapiens`, which no file of the
    /// release carries, and each cell whose file or table the installation
    /// did not find, one entry a reason.
    pub cells: Vec<LeftOut>,
}

/// The last lines of an installation, after the line that names the installed
/// directory. With no part left out, one line says the cache holds each part
/// that a file of the release carries; else one warning for each part left
/// out, "No `<part>` in the cache of `<species>` `<assembly>`: `<reason>`".
/// The same warning follows for each entry of `missing.sequences`, then of
/// `missing.cells`, and then one warning for each line of
/// `left_out_of_the_build`, the lines every build ends on
/// (`crate::left_out::LeftOut::lines`).
fn close(species: &str, assembly: &str, missing: &Missing, left_out_of_the_build: &[String]) {
    if missing.parts.is_empty() {
        info!("{}", every_part_line(species, assembly));
    }
    let lacks = missing
        .parts
        .iter()
        .chain(&missing.sequences)
        .chain(&missing.cells);
    for lack in lacks {
        warn!("{}", left_out_line(species, assembly, lack));
    }
    for line in left_out_of_the_build {
        warn!("{line}");
    }
}

/// The line of an installation that left no part out. It names the parts in
/// the order [`left_out`] lists them; the 1000 Genomes frequencies are
/// `homo_sapiens`' alone, the one species the release publishes them for.
fn every_part_line(species: &str, assembly: &str) -> String {
    let variants = if species == download::HUMAN_SPECIES {
        "known variants with their 1000 Genomes frequencies"
    } else {
        "known variants"
    };
    format!(
        "The cache of {species} {assembly} holds each part of Ensembl VEP's own cache that a \
         file of the release carries: {variants}, SIFT and PolyPhen predictions, regulatory \
         features and motif features"
    )
}

/// The route to a cell the installation did not bring.
const CONVERSION_CARRIES: &str =
    "a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

/// The reasons the site publishes no file of a part for the genome, each
/// `None` for a part whose files were found.
#[derive(Debug, Default)]
struct Unpublished {
    variants: Option<String>,
    matrices: Option<String>,
    regulatory: Option<String>,
}

/// The record of the cache built in `cache_dir`, its `info.json`: the keys
/// `variation_cols`, `sift`, `polyphen` and `regulatory` say which parts the
/// cache carries.
fn built_record(cache_dir: &Path) -> Result<vep_core::cache_info::CacheInfo> {
    let path = cache_dir.join("info.json");
    let file = File::open(&path).with_context(|| format!("Cannot open {}", path.display()))?;
    serde_json::from_reader(std::io::BufReader::new(file))
        .with_context(|| format!("Failed to read {}", path.display()))
}

/// The sequences of the cache built in `cache_dir` that hold no known
/// variant, in the order of its record's `contig_order`: each sequence the
/// cache holds transcripts on whose `variations/<sequence>/all_vars.gz` is
/// absent or holds no record.
fn sequences_without_known_variants(
    cache_dir: &Path,
    record: &vep_core::cache_info::CacheInfo,
) -> Result<Vec<String>> {
    let mut without = Vec::new();
    for sequence in &record.contig_order {
        let file = cache_dir
            .join("variations")
            .join(sequence)
            .join("all_vars.gz");
        if !file.is_file() || !holds_a_record(&file)? {
            without.push(sequence.clone());
        }
    }
    Ok(without)
}

/// Whether the gzip file at `path` inflates to a byte at least: a variation
/// file of no record is the end-of-file block alone.
fn holds_a_record(path: &Path) -> Result<bool> {
    let file = File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut first = [0u8; 1];
    let read = flate2::read::MultiGzDecoder::new(file)
        .read(&mut first)
        .with_context(|| format!("Failed to read {}", path.display()))?;
    Ok(read > 0)
}

/// The most sequences a closing line names; the rest are counted.
const NAMED_SEQUENCES: usize = 10;

/// The entry for the sequences of a cache on which `reader`, the installation
/// or the build, read no variation record (`without_variants`, from
/// [`sequences_without_known_variants`]); `None` where every sequence holds
/// one.
fn sequences_left_out(without_variants: &[String], reader: &str) -> Option<LeftOut> {
    (!without_variants.is_empty()).then(|| {
        LeftOut::new(
            &format!("known variants on {}", named_sequences(without_variants)),
            format!(
                "the {reader} read no variation record there, in a variation VCF or in the \
                 variation_feature table; a conversion of Ensembl VEP's cache carries any known \
                 variant that cache holds there (docs/cache-conversion.md)"
            ),
        )
    })
}

/// The warning of a build named one input at a time, whose cache stands in
/// `cache_dir`, for the sequences it read no variation record on: the
/// warning an installation ends on ([`close`]), naming the build. `None`
/// for a cache without variation files and where every sequence holds a
/// record.
pub(crate) fn sequences_warning_of_a_build(
    species: &str,
    assembly: &str,
    cache_dir: &Path,
) -> Result<Option<String>> {
    let record = built_record(cache_dir)?;
    if record.variation_cols.is_empty() {
        return Ok(None);
    }
    let without_variants = sequences_without_known_variants(cache_dir, &record)?;
    Ok(sequences_left_out(&without_variants, "build")
        .map(|lack| left_out_line(species, assembly, &lack)))
}

/// One closing warning: "No `<part>` in the cache of `<species>`
/// `<assembly>`: `<reason>`".
fn left_out_line(species: &str, assembly: &str, lack: &LeftOut) -> String {
    let LeftOut { part, reason } = lack;
    format!("No {part} in the cache of {species} {assembly}: {reason}")
}

/// `sequences` as a closing line names them: one by its name, more by their
/// count with the first [`NAMED_SEQUENCES`] of their names.
fn named_sequences(sequences: &[String]) -> String {
    let named = sequences
        .iter()
        .take(NAMED_SEQUENCES)
        .map(String::as_str)
        .collect::<Vec<_>>()
        .join(", ");
    match sequences.len() {
        1 => format!("sequence {named}"),
        n if n <= NAMED_SEQUENCES => format!("{n} sequences ({named})"),
        n => format!("{n} sequences ({named} and {} more)", n - NAMED_SEQUENCES),
    }
}

/// What the installed cache lacks of Ensembl VEP's own cache of the genome,
/// each list in the order known variants, matrices, regulatory tree. The
/// parts: each one the site
/// publishes no file of (`unpublished`) and each one its files hold nothing
/// of for this genome (read from the built cache's `record`). The
/// sequences: those of the cache the built variation files hold no record
/// on (`without_variants`, from [`sequences_without_known_variants`]), in
/// one entry; a variation file holds the records of the variation VCFs or,
/// for a sequence those hold none on, of the `variation_feature` table. The
/// cells: with variation files, the gnomAD frequencies of `homo_sapiens`,
/// which no file of the release carries; then each cell the record names
/// under `left_out` (`vep_core::cache_info::left_out`), with the reason the
/// lookup of its file or table kept (`brought`), the cells of one reason in
/// one entry: the PubMed identifiers, the variant synonyms and the failed
/// variants, then the cell-type activity and the matrix counts of the
/// regulatory tree. A cell the record does not name was brought and draws
/// no entry.
fn left_out(
    options: &InstallOptions,
    record: &vep_core::cache_info::CacheInfo,
    unpublished: &Unpublished,
    without_variants: &[String],
    brought: &Brought<'_, '_>,
) -> Missing {
    let human = options.species == download::HUMAN_SPECIES;
    let mut left_out = Vec::new();
    let mut sequences = Vec::new();
    let mut cells = Vec::new();

    match (&unpublished.variants, brought.vcfs) {
        (Some(reason), _) => left_out.push(LeftOut::new(
            if human {
                "known variants or their frequencies"
            } else {
                "known variants"
            },
            reason.clone(),
        )),
        (None, Some(vcfs)) => {
            let no_frequency = !record.variation_cols.iter().any(|column| column == "AF");
            if human && no_frequency {
                left_out.push(LeftOut::new(
                    "1000 Genomes frequencies",
                    format!(
                        "{}/CHECKSUMS lists no {}; {CONVERSION_CARRIES}",
                        vcfs.directory.relative,
                        download::THOUSAND_GENOMES_FILE
                    ),
                ));
            }
            sequences.extend(sequences_left_out(without_variants, "installation"));
            if human {
                cells.push(LeftOut::new(
                    "gnomAD frequencies",
                    format!("no file of the release carries them; {CONVERSION_CARRIES}"),
                ));
            }
            let lacking = brought
                .variant_cells
                .map(|sources| sources.left_out(record))
                .unwrap_or_default();
            cells.extend(cells_by_reason(&lacking, any_of));
        }
        (None, None) => {}
    }

    match &unpublished.matrices {
        Some(reason) => left_out.push(LeftOut::new("SIFT or PolyPhen predictions", reason.clone())),
        None => {
            let source = match (&options.prediction_server, brought.variation_db) {
                (None, Some(variation)) => variation.relative.clone(),
                _ => "Ensembl's MySQL server".to_string(),
            };
            for (tool, held) in [("SIFT", &record.sift), ("PolyPhen", &record.polyphen)] {
                if held.is_none() {
                    left_out.push(LeftOut::new(
                        &format!("{tool} predictions"),
                        format!("{source} holds no {tool} matrix of a protein of this genome"),
                    ));
                }
            }
        }
    }

    match (&unpublished.regulatory, brought.regulatory_files) {
        (Some(reason), _) => left_out.push(LeftOut::new(
            "regulatory features or motif features",
            reason.clone(),
        )),
        (None, Some(files)) => {
            for (part, kind, name) in [
                ("regulatory features", REGULATORY_FEATURES, &files.features),
                ("motif features", MOTIF_FEATURES, &files.motifs),
            ] {
                if name.is_none() {
                    left_out.push(LeftOut::new(
                        part,
                        format!(
                            "{}/CHECKSUMS lists no {kind} GFF3",
                            files.directory.relative
                        ),
                    ));
                }
            }
            let lacking = brought
                .regulatory_cells
                .map(|lookups| lookups.left_out(record, brought.regulatory_tree))
                .unwrap_or_default();
            cells.extend(cells_by_reason(&lacking, |cells| cells.join(", ")));
        }
        (None, None) => {}
    }
    Missing {
        parts: left_out,
        sequences,
        cells,
    }
}

/// What [`left_out`] reads of an installation beside the built cache's
/// record: what was found on the site for each part and for the cells of
/// each, and what the build wrote of the regulatory tree. Each is `None`
/// for a part the installation did not bring.
#[derive(Default)]
struct Brought<'r, 's> {
    /// The variation VCFs.
    vcfs: Option<&'r VariationVcfs<'s>>,
    /// The variation database.
    variation_db: Option<&'r Directory<'s>>,
    /// The lookups of the cells of the known variants.
    variant_cells: Option<&'r VariantCellSources<'s>>,
    /// The regulatory annotation files.
    regulatory_files: Option<&'r RegulatoryFiles<'s>>,
    /// The lookups of the cells of the regulatory tree.
    regulatory_cells: Option<&'r RegulatoryCells<'s>>,
    /// What the build wrote of the regulatory tree.
    regulatory_tree: Option<&'r regulatory::RegulatoryStats>,
}

/// One entry for each distinct reason of `cells`, each a cell with the
/// reason the cache holds no value for it, in the order the reasons first
/// stand: the entry names through `joined` the cells that share the reason
/// and ends on the route to them.
fn cells_by_reason(cells: &[(&str, String)], joined: impl Fn(&[&str]) -> String) -> Vec<LeftOut> {
    let mut reasons: Vec<&str> = Vec::new();
    for (_, reason) in cells {
        if !reasons.contains(&reason.as_str()) {
            reasons.push(reason);
        }
    }
    reasons
        .into_iter()
        .map(|reason| {
            let sharing: Vec<&str> = cells
                .iter()
                .filter(|(_, of_cell)| of_cell == reason)
                .map(|(cell, _)| *cell)
                .collect();
            LeftOut::new(&joined(&sharing), format!("{reason}; {CONVERSION_CARRIES}"))
        })
        .collect()
}

/// Whether the variation database of the release lists the dump of each table
/// the matrices are read from ([`predictions::PREDICTION_TABLES`]); else the
/// reason they cannot be read from it. For a genome `--include-predictions`
/// serves, the reason names that switch.
fn prediction_tables(variation: Option<&Directory<'_>>, options: &InstallOptions) -> Part<()> {
    let lacking = match variation {
        None => {
            Some("no variation database is published for this release and assembly".to_string())
        }
        Some(variation) => predictions::PREDICTION_TABLES
            .into_iter()
            .find(|table| table_name(variation, table).is_none())
            .map(|table| format!("{}/CHECKSUMS lists no {table} table", variation.relative)),
    };
    let Some(mut reason) = lacking else {
        return Part::Published(());
    };
    let served =
        predictions::prediction_database(&options.species, options.release, &options.assembly);
    if served.is_ok() {
        reason.push_str(
            "; --include-predictions fetches them from Ensembl's MySQL server (docs/cache-predictions.md)",
        );
    }
    Part::Unpublished(reason)
}

/// The kinds of regulatory annotation file an installation reads, as their
/// names carry them: `<Species>.<assembly>.<kind>.v<N>.gff3.gz`.
const REGULATORY_FEATURES: &str = "regulatory_features";
const MOTIF_FEATURES: &str = "motif_features";

/// The regulatory annotation files of an installation: the `annotation/`
/// directory of the genome's assembly and the regulatory-feature and
/// motif-feature GFF3 its listing names, each found before anything is
/// fetched. The EMAR file beside them is not read: Ensembl VEP's own cache
/// holds no EMAR.
struct RegulatoryFiles<'a> {
    directory: Directory<'a>,
    features: Option<String>,
    motifs: Option<String>,
}

impl<'a> RegulatoryFiles<'a> {
    /// Opens the `annotation/` directory of the assembly under the species'
    /// regulation directory, the first of [`assembly_directories`] the site
    /// holds.
    fn find(located: &'a Located, options: &InstallOptions) -> Result<Part<Self>> {
        let species = options.species.as_str();
        let assembly = options.assembly.as_str();
        let layout = &located.layout;
        let Some(species_dir) = layout.regulation_dir(options.release, species) else {
            return Ok(Part::Unpublished(format!(
                "{} publishes no regulatory annotation",
                layout.root(options.release)
            )));
        };
        let listed = located
            .site
            .list(&species_dir)
            .with_context(|| format!("Failed to list {species_dir}"))?;
        let mut found = None;
        let mut first = None;
        for name in assembly_directories(assembly, listed.as_deref()) {
            let relative = format!("{species_dir}/{name}/annotation");
            if let Some(directory) = optional_directory(&located.site, &relative)? {
                found = Some((directory, relative));
                break;
            }
            first.get_or_insert(relative);
        }
        let Some((directory, relative)) = found else {
            return Ok(Part::Unpublished(format!(
                "{}/CHECKSUMS is not on the site",
                first.unwrap_or(species_dir)
            )));
        };
        let features = regulatory_file(&directory.listing, REGULATORY_FEATURES);
        let motifs = regulatory_file(&directory.listing, MOTIF_FEATURES);
        if features.is_none() && motifs.is_none() {
            return Ok(Part::Unpublished(format!(
                "{relative}/CHECKSUMS lists no {REGULATORY_FEATURES} or {MOTIF_FEATURES} GFF3"
            )));
        }
        Ok(Part::Published(Self {
            directory,
            features,
            motifs,
        }))
    }

    /// Places the files in `dest_dir` and returns them as the regulatory
    /// inputs of the build, each promoter at its core span, as Ensembl VEP's
    /// own cache holds it.
    fn fetch(&self, dest_dir: &Path) -> Result<regulatory::RegulatoryInputs> {
        let fetch = |name: &Option<String>| -> Result<Option<PathBuf>> {
            name.as_deref()
                .map(|name| self.directory.fetch(name, dest_dir))
                .transpose()
        };
        Ok(regulatory::RegulatoryInputs {
            regulatory_gff: fetch(&self.features)?,
            emars_gff: None,
            motif_gff: fetch(&self.motifs)?,
            extended_promoters: false,
            activity: None,
            funcgen_dumps: None,
        })
    }
}

/// The kind of the activity file beside the regulatory GFF3 files, as its
/// name carries it: `<Species>.<assembly>.regulatory_activity.v<N>.tsv.gz`.
const REGULATORY_ACTIVITY: &str = "regulatory_activity";

/// A file or a table looked up on the site: what was found, or the reason
/// nothing was, worded as the reason of a part left out is.
type Lookup<T> = std::result::Result<T, String>;

/// What an installation reads, beside the GFF3 files, for the cells of the
/// regulatory tree those files do not carry: the activity file of the
/// annotation directory, and tables of the funcgen database named beside
/// the core database ([`download::funcgen_database`]). Each is looked up
/// before anything is fetched. One the site does not publish is held as the
/// reason, logged in one line, and the tree is baked without it.
struct RegulatoryCells<'a> {
    /// The activity file the annotation directory's listing names: the
    /// cell-type activity of the regulatory features.
    activity: Lookup<String>,
    /// The directory of the funcgen database.
    funcgen: Lookup<FuncgenTables<'a>>,
}

/// The directory of a funcgen database and whether its listing names the
/// tables an installation reads from it.
struct FuncgenTables<'a> {
    directory: Directory<'a>,
    /// Both of [`regulatory::MATRIX_TABLES`]: the counts of the matrices the
    /// motif features name.
    matrices: Lookup<()>,
    /// Both of [`regulatory::COMPLEX_TABLES`]: the factors of those matrices.
    complexes: Lookup<()>,
    /// [`regulatory::REGULATORY_BUILD_TABLE`]: the regulatory build's version.
    regulatory_build: Lookup<()>,
}

impl<'a> RegulatoryCells<'a> {
    /// Looks the activity file up in the listing of the annotation directory
    /// of `files`, the one of the highest version, and opens the funcgen
    /// database's directory, reading in its listing which tables it names.
    fn find(
        located: &'a Located,
        options: &InstallOptions,
        files: &RegulatoryFiles<'_>,
    ) -> Result<Self> {
        let annotation = &files.directory;
        let activity = versioned_file(&annotation.listing, REGULATORY_ACTIVITY, ".tsv.gz")
            .ok_or_else(|| {
                format!(
                    "{}/CHECKSUMS lists no {REGULATORY_ACTIVITY} file",
                    annotation.relative
                )
            });
        let funcgen = match located.funcgen_dir(&options.species, options.release) {
            Some(relative) => match optional_directory(&located.site, &relative)? {
                Some(directory) => Ok(FuncgenTables::of(directory)),
                None => Err(format!("{relative}/CHECKSUMS is not on the site")),
            },
            None => Err("no funcgen database is named beside the core database".to_string()),
        };
        let cells = Self { activity, funcgen };
        for (what, reason) in cells.not_found() {
            info!("No {what} is read for the regulatory tree: {reason}");
        }
        Ok(cells)
    }

    /// What the site does not publish of the files and tables looked up,
    /// each with the reason: the activity file, the funcgen database or,
    /// in its directory, a group of tables. The factors' tables are named
    /// only where the matrices' are listed, which they are read beside.
    fn not_found(&self) -> Vec<(&'static str, &str)> {
        let mut not_found = Vec::new();
        if let Err(reason) = &self.activity {
            not_found.push(("cell-type activity", reason.as_str()));
        }
        match &self.funcgen {
            Err(reason) => not_found.push(("table of the funcgen database", reason.as_str())),
            Ok(tables) => {
                let complexes = match &tables.matrices {
                    Ok(()) => &tables.complexes,
                    Err(_) => &Ok(()),
                };
                for (what, listed) in [
                    ("binding matrix table", &tables.matrices),
                    ("transcription factor complex table", complexes),
                    ("regulatory build table", &tables.regulatory_build),
                ] {
                    if let Err(reason) = listed {
                        not_found.push((what, reason.as_str()));
                    }
                }
            }
        }
        not_found
    }

    /// The cells of the regulatory tree that `record`, the built cache's,
    /// names as left out, each with the reason, the activity first. The
    /// activity's reason is its file's lookup. The matrix counts' is the
    /// funcgen database's lookup, else that of the matrices' tables, else
    /// what `tree`, the build's own account, found short: a matrix a motif
    /// names without counts in those tables, a motif whose span the genome
    /// FASTA does not hold.
    fn left_out(
        &self,
        record: &vep_core::cache_info::CacheInfo,
        tree: Option<&regulatory::RegulatoryStats>,
    ) -> Vec<(&'static str, String)> {
        use vep_core::cache_info::left_out as name;
        let scores = match &self.funcgen {
            Err(reason) => Some(reason.clone()),
            Ok(tables) => match &tables.matrices {
                Err(reason) => Some(reason.clone()),
                Ok(()) => tree.and_then(|tree| {
                    let mut short = Vec::new();
                    if !tree.counts_on_every_matrix {
                        short.push(format!(
                            "{} holds no counts of a matrix a motif feature names",
                            tables.directory.relative
                        ));
                    }
                    if !tree.sequence_on_every_motif {
                        short.push(
                            "the genome FASTA does not hold the span of a motif feature"
                                .to_string(),
                        );
                    }
                    (!short.is_empty()).then(|| short.join(" and "))
                }),
            },
        };
        [
            (
                name::CELL_TYPE_ACTIVITY,
                "cell-type activity of regulatory features",
                self.activity.as_ref().err().cloned(),
            ),
            (
                name::MOTIF_SCORES,
                "HIGH_INF_POS or MOTIF_SCORE_CHANGE of motif features",
                scores,
            ),
        ]
        .into_iter()
        .filter(|(token, _, _)| record.leaves_out(token))
        .filter_map(|(_, cell, reason)| Some((cell, reason?)))
        .collect()
    }

    /// Places what the tree of `files` reads beside its GFF3 files and names
    /// it in `inputs`: the activity file in `regulatory_dir`, where the
    /// regulatory-feature GFF3 was found; in `dumps_dir` the matrices'
    /// tables, where the motif-feature GFF3 was found, with the factors'
    /// tables where the listing names them, and the regulatory build's
    /// table.
    fn fetch(
        &self,
        files: &RegulatoryFiles<'_>,
        regulatory_dir: &Path,
        dumps_dir: &Path,
        inputs: &mut regulatory::RegulatoryInputs,
    ) -> Result<()> {
        if let (Some(_), Ok(name)) = (&files.features, &self.activity) {
            inputs.activity = Some(files.directory.fetch(name, regulatory_dir)?);
        }
        let Ok(funcgen) = &self.funcgen else {
            return Ok(());
        };
        let mut tables: Vec<&str> = Vec::new();
        if files.motifs.is_some() && funcgen.matrices.is_ok() {
            tables.extend(regulatory::MATRIX_TABLES);
            if funcgen.complexes.is_ok() {
                tables.extend(regulatory::COMPLEX_TABLES);
            }
        }
        if funcgen.regulatory_build.is_ok() {
            tables.push(regulatory::REGULATORY_BUILD_TABLE);
        }
        for table in &tables {
            fetch_table(&funcgen.directory, table, dumps_dir)?.ok_or_else(|| {
                anyhow::anyhow!("funcgen table {table} is not published for this release")
            })?;
        }
        if !tables.is_empty() {
            inputs.funcgen_dumps = Some(dumps_dir.to_path_buf());
        }
        Ok(())
    }
}

impl<'a> FuncgenTables<'a> {
    /// The tables the listing of `directory` names, each group held as
    /// found or as the reason naming the first of its tables the listing
    /// lacks.
    fn of(directory: Directory<'a>) -> Self {
        let listed = |tables: &[&str]| -> Lookup<()> {
            match tables
                .iter()
                .find(|table| table_name(&directory, table).is_none())
            {
                Some(table) => Err(format!(
                    "{}/CHECKSUMS lists no {table} table",
                    directory.relative
                )),
                None => Ok(()),
            }
        };
        let matrices = listed(&regulatory::MATRIX_TABLES);
        let complexes = listed(&regulatory::COMPLEX_TABLES);
        let regulatory_build = listed(&[regulatory::REGULATORY_BUILD_TABLE]);
        Self {
            directory,
            matrices,
            complexes,
            regulatory_build,
        }
    }
}

/// The directories an assembly's regulatory annotation may stand under, in
/// the order they are tried. The directory is named for the assembly, or for
/// the last dot-separated part of its name (chicken's
/// `bGalGal1.mat.broiler.GRCg7b` stands under `GRCg7b`). With the species
/// directory's `listed` entries it is the one of them so named, else the
/// assembly's name; a site that gives no listing is asked for both names.
fn assembly_directories(assembly: &str, listed: Option<&[String]>) -> Vec<String> {
    let last_part = assembly.rsplit_once('.').map(|(_, last)| last);
    match listed {
        Some(listed) => {
            let named = listed
                .iter()
                .find(|name| name.as_str() == assembly)
                .or_else(|| listed.iter().find(|name| Some(name.as_str()) == last_part));
            vec![named.map_or(assembly, String::as_str).to_string()]
        }
        None => std::iter::once(assembly)
            .chain(last_part)
            .map(str::to_string)
            .collect(),
    }
}

/// The file of `kind` a regulatory annotation directory lists,
/// `<anything>.<kind>.v<N>.gff3.gz`: the one of the highest `N`.
fn regulatory_file(listing: &Listing, kind: &str) -> Option<String> {
    versioned_file(listing, kind, ".gff3.gz")
}

/// The file of `kind` and `extension` a regulatory annotation directory
/// lists, `<anything>.<kind>.v<N><extension>`: the one of the highest `N`.
fn versioned_file(listing: &Listing, kind: &str, extension: &str) -> Option<String> {
    listing
        .names()
        .into_iter()
        .filter_map(|name| {
            let (head, version) = name.strip_suffix(extension)?.rsplit_once(".v")?;
            if !head.ends_with(&format!(".{kind}")) {
                return None;
            }
            Some((version.parse::<u32>().ok()?, name.to_string()))
        })
        .max()
        .map(|(_, name)| name)
}

/// The name a MySQL directory lists `table` under: `<table>.txt.gz`, else
/// `<table>.txt`; `None` when it lists neither.
fn table_name(directory: &Directory<'_>, table: &str) -> Option<String> {
    [format!("{table}.txt.gz"), format!("{table}.txt")]
        .into_iter()
        .find(|name| directory.listing.contains(name))
}

/// Places `<table>.txt.gz`, else `<table>.txt`, from a MySQL directory in
/// `dest_dir`; `None` when the directory lists neither.
fn fetch_table(directory: &Directory<'_>, table: &str, dest_dir: &Path) -> Result<Option<PathBuf>> {
    match table_name(directory, table) {
        Some(name) => directory.fetch(&name, dest_dir).map(Some),
        None => Ok(None),
    }
}

/// Places the release's genome FASTA in `cache_dir` as an indexed BGZF file:
/// the `dna_index/` file with the `.fai` and `.gzi` published beside it, else
/// the `dna/` file recompressed and indexed here. A `dna_index/` file that
/// does not open through its published indexes is installed as a `dna/` file
/// is, with a warning that carries the refusal. Returns the FASTA's path.
fn install_genome_fasta(
    genome: &GenomeFasta<'_>,
    inputs_dir: &Path,
    cache_dir: &Path,
) -> Result<PathBuf> {
    let name = &genome.file.name;
    let target = cache_dir.join(name);
    // Kept with the other inputs, so a build that fails after the download
    // does not fetch the genome again.
    let source = genome.file.fetch(inputs_dir)?;
    if genome.indexed {
        link_or_copy(&source, &target)?;
        let mut placed = vec![target.clone()];
        for suffix in ["fai", "gzi"] {
            let index = format!("{name}.{suffix}");
            let staged = genome.file.directory.fetch(&index, inputs_dir)?;
            let dest = cache_dir.join(&index);
            link_or_copy(&staged, &dest)?;
            placed.push(dest);
        }
        let Err(refusal) = vep_fasta::IndexedFasta::from_path(&target) else {
            return Ok(target);
        };
        // The storage mode comes from the environment, so indexes written
        // here would not open the file either.
        if matches!(refusal, vep_fasta::FastaError::InvalidStorageMode(_)) {
            return Err(refusal)
                .with_context(|| format!("The indexed FASTA {} does not open", target.display()));
        }
        warn!(
            "The .fai and .gzi published beside {name} do not open it, so both are written \
             here from the file itself: {refusal}"
        );
        // A file placed as a link shares its bytes with the staged download,
        // which a write through the link would change.
        for file in &placed {
            std::fs::remove_file(file)
                .with_context(|| format!("Failed to remove {}", file.display()))?;
        }
    }
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

/// The variation VCFs of an installation: the directory of the genome's
/// variation VCFs and the germline files fetched from it, each found in its
/// listing before anything is fetched.
struct VariationVcfs<'a> {
    directory: Directory<'a>,
    germline: Vec<String>,
}

impl<'a> VariationVcfs<'a> {
    /// Opens the genome's variation directory and looks up the germline
    /// files of [`download::germline_vcf_names`] in its listing. A tree that
    /// publishes no such directory, or a listing that names none of the
    /// files, publishes no known variant of the genome. A listing that names
    /// some of `homo_sapiens`' per-chromosome files and lacks another, or
    /// lacks one `--chromosomes` names, refuses the installation.
    fn find(located: &'a Located, options: &InstallOptions) -> Result<Part<Self>> {
        let species = options.species.as_str();
        let Some(published) = located.layout.variation_vcf_dir(options.release, species) else {
            return Ok(Part::Unpublished(format!(
                "{} publishes no variation VCF",
                located.layout.root(options.release)
            )));
        };
        let Some(directory) = optional_directory(&located.site, &published)? else {
            return Ok(Part::Unpublished(format!(
                "{published}/CHECKSUMS is not on the site"
            )));
        };
        let germline = download::germline_vcf_names(species, options.chromosomes.as_deref());
        let per_chromosome = species == download::HUMAN_SPECIES;
        let named = per_chromosome && options.chromosomes.is_some();
        if !named && !germline.iter().any(|name| directory.listing.contains(name)) {
            let files = if per_chromosome {
                format!("none of the per-chromosome files {species}-chr<N>.vcf.gz")
            } else {
                format!("no {}", download::genome_vcf_name(species))
            };
            return Ok(Part::Unpublished(format!(
                "{published}/CHECKSUMS lists {files}"
            )));
        }
        if let Some(name) = germline
            .iter()
            .find(|name| !directory.listing.contains(name))
        {
            bail!(
                "{name} is not published for this release and assembly; name the chromosomes with --chromosomes"
            );
        }
        Ok(Part::Published(Self {
            directory,
            germline,
        }))
    }
}

/// The tables an installation reads the known variants of a sequence from
/// where the one genome-wide variation VCF of a species holds no record on
/// it: `variation_feature` and `seq_region` of the variation database, with
/// `attrib` where the listing names it, each under the name the listing
/// gives it. Looked up before anything is fetched; fetched when the build
/// finds such a sequence.
struct SequenceFeatureTables {
    /// The directory of the variation database.
    relative: String,
    /// The listing's lines of the tables.
    listed: HashMap<String, Checksum>,
}

impl SequenceFeatureTables {
    /// Reads in the listing of the variation database whether it names
    /// `variation_feature` and `seq_region`; else the reason, the first of
    /// the two it lacks.
    fn find(variation_db: Option<&Directory<'_>>) -> Lookup<Self> {
        let Some(variation) = variation_db else {
            return Err(
                "no variation database is published for this release and assembly".to_string(),
            );
        };
        let mut listed = HashMap::new();
        for (table, needed) in [
            ("variation_feature", true),
            ("seq_region", true),
            ("attrib", false),
        ] {
            match table_name(variation, table) {
                Some(name) => {
                    let line = variation.listing.get(&name).expect("a listed name");
                    listed.insert(name, line);
                }
                None if needed => {
                    return Err(format!(
                        "{}/CHECKSUMS lists no {table} table",
                        variation.relative
                    ))
                }
                None => {}
            }
        }
        Ok(Self {
            relative: variation.relative.clone(),
            listed,
        })
    }

    /// The tables as an input of the build: placed under
    /// `<work_dir>/variation_dumps/` from `site` when the build asks, with
    /// `<work_dir>/variation_features/` as their scratch directory.
    /// `named` is what `--chromosomes` names.
    fn fetched_when_asked(
        &self,
        site: &Site,
        work_dir: &Path,
        named: Option<Vec<String>>,
    ) -> SequenceFeatures {
        let (site, relative, listed) = (site.clone(), self.relative.clone(), self.listed.clone());
        let (dumps_dir, scratch) = (
            work_dir.join("variation_dumps"),
            work_dir.join("variation_features"),
        );
        SequenceFeatures {
            named,
            fetch: Box::new(move || {
                let directory = Directory {
                    site: &site,
                    relative: relative.clone(),
                    listing: Listing {
                        entries: listed.clone(),
                    },
                };
                let table = |table: &str| -> Result<PathBuf> {
                    fetch_table(&directory, table, &dumps_dir)?.ok_or_else(|| {
                        anyhow::anyhow!("variation table {table} is not published for this release")
                    })
                };
                Ok(variation::FeatureTables {
                    variation_feature: table("variation_feature")?,
                    seq_region: table("seq_region")?,
                    attrib: fetch_table(&directory, "attrib", &dumps_dir)?,
                    scratch: scratch.clone(),
                })
            }),
        }
    }
}

/// The files an installation reads, beside the variation VCFs, for the
/// cells of the known variants no VCF carries: tables of the variation
/// database for the PubMed identifiers, the synonyms and the failed flag
/// ([`variant_cells::PUBMED_TABLES`], [`variant_cells::SYNONYM_TABLES`],
/// [`variant_cells::FAILED_TABLES`]), and the failed-variants GVF of the
/// genome's GVF directory. Each is looked up before anything is fetched. One
/// the site does not publish is held as the reason, logged in one line, and
/// the variants are written without what it carries; none ends the
/// installation.
struct VariantCellSources<'a> {
    /// Every table of the `pubmed` cell is listed.
    pubmed: Lookup<()>,
    /// Every table of the `var_synonyms` cell is listed.
    synonyms: Lookup<()>,
    /// Both tables of the `failed` cell are listed.
    failed_variation: Lookup<()>,
    /// The failed-variants GVF, in the listing of its directory.
    failed_gvf: Lookup<Input<'a>>,
}

impl<'a> VariantCellSources<'a> {
    /// Reads in the listing of the variation database which of the tables it
    /// names, and looks the failed-variants GVF up in the listing of the
    /// genome's GVF directory.
    fn find(
        located: &'a Located,
        options: &InstallOptions,
        variation_db: Option<&Directory<'_>>,
    ) -> Result<Self> {
        let listed = |tables: &[&str]| -> Lookup<()> {
            let Some(variation) = variation_db else {
                return Err(
                    "no variation database is published for this release and assembly".to_string(),
                );
            };
            match tables
                .iter()
                .find(|table| table_name(variation, table).is_none())
            {
                Some(table) => Err(format!(
                    "{}/CHECKSUMS lists no {table} table",
                    variation.relative
                )),
                None => Ok(()),
            }
        };
        let species = options.species.as_str();
        let name = download::failed_gvf_name(species);
        let failed_gvf = match located.layout.variation_gvf_dir(options.release, species) {
            Some(relative) => match optional_directory(&located.site, &relative)? {
                Some(directory) if directory.listing.contains(&name) => {
                    Ok(Input { directory, name })
                }
                Some(_) => Err(format!("{relative}/CHECKSUMS lists no {name}")),
                None => Err(format!("{relative}/CHECKSUMS is not on the site")),
            },
            None => Err(format!(
                "{} publishes no variation GVF",
                located.layout.root(options.release)
            )),
        };
        let sources = Self {
            pubmed: listed(&variant_cells::PUBMED_TABLES),
            synonyms: listed(&variant_cells::SYNONYM_TABLES),
            failed_variation: listed(&variant_cells::FAILED_TABLES),
            failed_gvf,
        };
        for (what, reason) in sources.not_found() {
            info!("No {what} is read for the known variants: {reason}");
        }
        Ok(sources)
    }

    /// What the site does not publish of the tables and the file looked up,
    /// each with the reason.
    fn not_found(&self) -> Vec<(&'static str, &str)> {
        let mut not_found = Vec::new();
        for (what, listed) in [
            ("PubMed identifier", &self.pubmed),
            ("variant synonym", &self.synonyms),
            ("failed_variation table", &self.failed_variation),
        ] {
            if let Err(reason) = listed {
                not_found.push((what, reason.as_str()));
            }
        }
        if let Err(reason) = &self.failed_gvf {
            not_found.push(("failed-variants GVF", reason.as_str()));
        }
        not_found
    }

    /// The cells of the known variants that `record`, the built cache's,
    /// names as left out, each with the reason its lookup kept, in the order
    /// PubMed identifiers, variant synonyms, failed variants. The failed
    /// variants are read from the GVF, else flagged from `failed_variation`:
    /// the record names them where neither was found, and their reason is
    /// the GVF's.
    fn left_out(&self, record: &vep_core::cache_info::CacheInfo) -> Vec<(&'static str, String)> {
        use vep_core::cache_info::left_out as name;
        [
            (
                name::PUBMED,
                "PubMed identifiers",
                self.pubmed.as_ref().err(),
            ),
            (
                name::VAR_SYNONYMS,
                "variant synonyms",
                self.synonyms.as_ref().err(),
            ),
            (
                name::FAILED_VARIANTS,
                "failed variants",
                self.failed_gvf.as_ref().err(),
            ),
        ]
        .into_iter()
        .filter(|(token, _, _)| record.leaves_out(token))
        .filter_map(|(_, cell, reason)| Some((cell, reason?.clone())))
        .collect()
    }

    /// Places the tables of each cell whose every table is listed under
    /// `<work_dir>/variation_dumps/`, each table once, and the
    /// failed-variants GVF under `<work_dir>/variation_gvf/`, and returns
    /// them as the inputs of the cells, with `<work_dir>/variant_cells/` as
    /// their scratch directory; `None` where nothing was found.
    fn fetch(
        &self,
        variation_db: Option<&Directory<'_>>,
        work_dir: &Path,
    ) -> Result<Option<variant_cells::VariantCellInputs>> {
        let dumps_dir = work_dir.join("variation_dumps");
        let mut tables: Vec<&str> = Vec::new();
        for (listed, group) in [
            (&self.pubmed, &variant_cells::PUBMED_TABLES[..]),
            (&self.synonyms, &variant_cells::SYNONYM_TABLES[..]),
            (&self.failed_variation, &variant_cells::FAILED_TABLES[..]),
        ] {
            for table in group.iter().filter(|_| listed.is_ok()) {
                if !tables.contains(table) {
                    tables.push(table);
                }
            }
        }
        let mut fetched: HashMap<&str, PathBuf> = HashMap::new();
        if let Some(variation) = variation_db {
            for table in tables {
                let path = fetch_table(variation, table, &dumps_dir)?.ok_or_else(|| {
                    anyhow::anyhow!("variation table {table} is not published for this release")
                })?;
                fetched.insert(table, path);
            }
        }
        let failed_gvf = match &self.failed_gvf {
            Ok(gvf) => Some(gvf.fetch(&work_dir.join("variation_gvf"))?),
            Err(_) => None,
        };
        if fetched.is_empty() && failed_gvf.is_none() {
            return Ok(None);
        }
        let mut table = |name: &str| fetched.remove(name);
        Ok(Some(variant_cells::VariantCellInputs {
            variation: table("variation"),
            variation_citation: table("variation_citation"),
            publication: table("publication"),
            variation_synonym: table("variation_synonym"),
            source: table("source"),
            failed_variation: table("failed_variation"),
            failed_gvf,
            scratch: work_dir.join("variant_cells"),
        }))
    }
}

/// Places the germline files of `vcfs` and the further genome-wide files its
/// directory lists, each with its `.csi`, under `<work_dir>/variation_vcf/`,
/// and the variation database's ClinVar tables under
/// `<work_dir>/variation_dumps/`; a ClinVar table already fetched into
/// `core_dumps_dir` is linked from there, not fetched again. `--chromosomes`
/// names the per-chromosome files of `homo_sapiens` and narrows the read of
/// any other species' one genome-wide file; it narrows no table.
fn fetch_variations(
    options: &InstallOptions,
    vcfs: &VariationVcfs<'_>,
    variation_db: Option<&Directory<'_>>,
    work_dir: &Path,
    core_dumps_dir: &Path,
) -> Result<VariationBuild> {
    let species = options.species.as_str();
    let vcf_dir = work_dir.join("variation_vcf");
    let directory = &vcfs.directory;
    for name in &vcfs.germline {
        directory.fetch(name, &vcf_dir)?;
        fetch_index(directory, name, &vcf_dir)?;
    }
    for name in download::genome_wide_names(species) {
        if directory.listing.contains(&name) && !vcfs.germline.contains(&name) {
            directory.fetch(&name, &vcf_dir)?;
            fetch_index(directory, &name, &vcf_dir)?;
        }
    }
    let inputs =
        download::local_variation_inputs(&vcf_dir, species, options.chromosomes.as_deref())?;

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
    Ok(VariationBuild {
        inputs,
        dumps,
        cells: None,
        features: None,
    })
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

    /// A listed name that is not one bare file name, absolute, climbing out
    /// or naming a subdirectory, is passed over, so no fetch places a file
    /// outside its destination directory.
    #[test]
    fn listing_passes_over_names_that_are_not_a_file_name() {
        let listing = Listing::parse(
            "12345 1 /home/user/notes.primary_assembly.fa.gz\n\
             12345 1 ../../../x.primary_assembly.fa.gz\n\
             12345 1 sub/x.primary_assembly.fa.gz\n\
             12345 1 ..\n\
             12345 1 .\n\
             12345 1 Homo_sapiens.GRCh38.116.gff3.gz\n",
        );
        assert_eq!(listing.names(), ["Homo_sapiens.GRCh38.116.gff3.gz"]);
        assert_eq!(genome_fasta_name(&listing), None);
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
        let dir = directory(&site, "pub/release-116/gff3/homo_sapiens")
            .unwrap()
            .unwrap();
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

    /// Writes `relative` under `root` and adds its line to the `CHECKSUMS` of
    /// its directory.
    fn publish(root: &Path, relative: &str, bytes: &[u8]) {
        let path = root.join(relative);
        let dir = path.parent().unwrap();
        std::fs::create_dir_all(dir).unwrap();
        std::fs::write(&path, bytes).unwrap();
        let name = path.file_name().unwrap().to_str().unwrap();
        let mut checksums = std::fs::read_to_string(dir.join("CHECKSUMS")).unwrap_or_default();
        checksums.push_str(&format!("{} {name}\n", bsd_sum_of(bytes)));
        std::fs::write(dir.join("CHECKSUMS"), checksums).unwrap();
    }

    fn site_at(root: &Path) -> Site {
        Site::parse(root.to_str().unwrap()).unwrap()
    }

    const VERTEBRATES: &str = "#name\tspecies\tdivision\ttaxonomy_id\tassembly\tassembly_accession\tgenebuild\tvariation\tmicroarray\tpan_compara\tpeptide_compara\tgenome_alignments\tother_alignments\tcore_db\tspecies_id\n\
        Rat\trattus_norvegicus\tEnsemblVertebrates\t10116\tGRCr8\tGCA_036323735.1\t2024-02\tY\tN\tN\tY\tY\tY\trattus_norvegicus_core_116_1\t1\t\n";

    const GENOMES: &str = "#name\tspecies\tdivision\ttaxonomy_id\tassembly\tassembly_accession\tgenebuild\tvariation\tmicroarray\tpan_compara\tpeptide_compara\tgenome_alignments\tother_alignments\tcore_db\tspecies_id\n\
        A plant\tplanta_ficta\tEnsemblPlants\t1\tPf1.0\tGCA_000000002.1\t2020-01\tN\tN\tN\tN\tN\tN\tplanta_ficta_core_63_116_4\t1\n\
        A bacterium\tbacterium_fictum_gca_000000001\tEnsemblBacteria\t2\tASM1v1\tGCA_000000001.1\t2020-01\tN\tN\tN\tN\tN\tN\tbacteria_7_collection_core_63_116_1\t12\n";

    /// A species list names the core database whatever the assembly is
    /// called: GRCr8's is numbered 1. The list is believed only when it sums
    /// to its `CHECKSUMS` line, and a species it does not name is refused,
    /// naming the list.
    #[test]
    fn species_list_names_the_database_the_assembly_does_not() {
        let root = tempfile::tempdir().unwrap();
        publish(
            root.path(),
            "pub/release-116/species_EnsemblVertebrates.txt",
            VERTEBRATES.as_bytes(),
        );
        let site = site_at(root.path());

        let located = locate(
            std::slice::from_ref(&site),
            "rattus_norvegicus",
            "GRCr8",
            116,
        )
        .unwrap();
        assert_eq!(located.layout, Layout::Release);
        assert_eq!(
            located.core_db.as_deref(),
            Some("rattus_norvegicus_core_116_1")
        );
        assert_eq!(
            located.databases("rattus_norvegicus", 116).unwrap(),
            (
                "pub/release-116/mysql/rattus_norvegicus_core_116_1".to_string(),
                Some("pub/release-116/mysql/rattus_norvegicus_variation_116_1".to_string())
            )
        );
        located
            .buildable("rattus_norvegicus", "GRCr8", 116)
            .unwrap();

        let error = locate(std::slice::from_ref(&site), "rattus_rattus", "Rrat1", 116)
            .unwrap_err()
            .to_string();
        assert!(
            error.starts_with("rattus_rattus is not in the species list of release 116 (")
                && error.contains("pub/release-116/species_EnsemblVertebrates.txt")
                && error.contains("--release"),
            "{error}"
        );

        std::fs::write(
            root.path()
                .join("pub/release-116/species_EnsemblVertebrates.txt"),
            VERTEBRATES.replace("_116_1", "_116_2"),
        )
        .unwrap();
        let error = locate(&[site], "rattus_norvegicus", "GRCr8", 116)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("does not match the CHECKSUMS line for species_EnsemblVertebrates.txt"),
            "{error}"
        );
    }

    /// A tree that publishes no species list has its `mysql/` directory read
    /// for the core database: the one named for the species and the release,
    /// whatever number ends it, in the release tree and in GRCh37's. No such
    /// database, two of them and a missing `mysql/` directory are each
    /// refused.
    #[test]
    fn a_tree_without_a_species_list_is_searched_for_the_database() {
        let root = tempfile::tempdir().unwrap();
        for database in [
            "pub/release-116/mysql/bos_taurus_core_116_2",
            "pub/release-116/mysql/bos_taurus_variation_116_2",
            "pub/release-116/mysql/sus_scrofa_core_116_111",
            "pub/release-116/mysql/sus_scrofa_core_116_112",
            "pub/grch37/release-116/mysql/homo_sapiens_core_116_37",
        ] {
            std::fs::create_dir_all(root.path().join(database)).unwrap();
        }
        let sites = [site_at(root.path())];

        let cow = locate(&sites, "bos_taurus", "ARS-UCD2.0", 116).unwrap();
        assert_eq!(cow.layout, Layout::Release);
        assert_eq!(cow.core_db.as_deref(), Some("bos_taurus_core_116_2"));
        let human = locate(&sites, "homo_sapiens", "GRCh37", 116).unwrap();
        assert_eq!(human.layout, Layout::Grch37);
        assert_eq!(human.core_db.as_deref(), Some("homo_sapiens_core_116_37"));

        let error = locate(&sites, "ovis_aries", "Oar1", 116)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("pub/release-116/mysql/ holds no core database of ovis_aries"),
            "{error}"
        );
        let error = locate(&sites, "sus_scrofa", "Sscrofa11.1", 116)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("more than one core database of sus_scrofa")
                && error.contains("sus_scrofa_core_116_111")
                && error.contains("sus_scrofa_core_116_112"),
            "{error}"
        );
        let error = locate(&sites, "bos_taurus", "ARS-UCD2.0", 115)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("pub/release-115/mysql/ is not on the site"),
            "{error}"
        );
    }

    /// On a server that serves no index page the `mysql/` directory lists
    /// nothing: the database named for the assembly's number is taken when
    /// the site holds its directory, and refused when it does not or when the
    /// assembly's name ends in no number. A directory that lists its entries
    /// is never asked for a name it does not list.
    #[test]
    fn a_server_without_index_pages_is_asked_for_the_numbered_database() {
        let mysql = "pub/release-116/mysql";
        let held = |database: &str| Ok(database == "homo_sapiens_core_116_38");
        let among = |entries: Option<&[String]>,
                     species: &str,
                     assembly: &str,
                     holds: &dyn Fn(&str) -> Result<bool>| {
            core_database_among(mysql, entries, species, assembly, 116, holds)
        };
        assert_eq!(
            among(None, "homo_sapiens", "GRCh38", &held).unwrap(),
            "homo_sapiens_core_116_38"
        );
        for (species, assembly) in [("rattus_norvegicus", "GRCr8"), ("bos_taurus", "ARS-UCD")] {
            let error = among(None, species, assembly, &held).unwrap_err();
            assert!(download::is_not_published(&error), "{error}");
            assert!(
                error
                    .to_string()
                    .contains("pub/release-116/mysql/ is not on the site"),
                "{error}"
            );
        }
        let listed = ["rattus_norvegicus_core_116_1".to_string()];
        assert_eq!(
            among(Some(&listed), "rattus_norvegicus", "GRCr8", &|_| {
                panic!("the listing answers")
            })
            .unwrap(),
            "rattus_norvegicus_core_116_1"
        );
    }

    /// A species the first site's list does not name is found on the second,
    /// in the tree of its division when its core database is its own and in
    /// its collection's when the database is shared; the collection's genome
    /// is refused, naming the database and both routes that build its cache.
    #[test]
    fn an_ensembl_genomes_species_is_placed_by_the_species_list() {
        let ensembl = tempfile::tempdir().unwrap();
        publish(
            ensembl.path(),
            "pub/release-63/species_EnsemblVertebrates.txt",
            VERTEBRATES.as_bytes(),
        );
        let genomes = tempfile::tempdir().unwrap();
        publish(
            genomes.path(),
            "pub/release-63/species.txt",
            GENOMES.as_bytes(),
        );
        let sites = [site_at(ensembl.path()), site_at(genomes.path())];

        let plant = locate(&sites, "planta_ficta", "Pf1.0", 63).unwrap();
        assert!(matches!(&plant.site, Site::Dir(dir) if dir == genomes.path()));
        assert_eq!(
            plant.layout,
            Layout::Division {
                division: "plants".to_string()
            }
        );
        assert_eq!(
            plant.databases("planta_ficta", 63).unwrap(),
            (
                "pub/release-63/plants/mysql/planta_ficta_core_63_116_4".to_string(),
                Some("pub/release-63/plants/mysql/planta_ficta_variation_63_116_4".to_string())
            )
        );
        plant.buildable("planta_ficta", "Pf1.0", 63).unwrap();

        let species = "bacterium_fictum_gca_000000001";
        let bacterium = locate(&sites, species, "ASM1v1", 63).unwrap();
        assert_eq!(
            bacterium.layout,
            Layout::Collection {
                division: "bacteria".to_string(),
                collection: "bacteria_7_collection".to_string()
            }
        );
        let error = bacterium
            .buildable(species, "ASM1v1", 63)
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            format!(
                "bacterium_fictum_gca_000000001 ASM1v1 cannot be installed from {}/: it has no core database of its own (bacteria_7_collection_core_63_116_1 holds every genome of bacteria_7_collection, and a shared database is not split by genome); convert Ensembl VEP's cache of it instead (docs/cache-conversion.md), or build from local copies of its files with --gff3, --fasta, --genome-fasta and --codon-table",
                genomes.path().display()
            )
        );
    }

    /// A build named one input at a time takes its download addresses from
    /// the located genome: the GFF3, the protein FASTA and the two MySQL
    /// directories of a species of a division, each under the site's root; a
    /// genome of a collection is refused there as under `--auto`.
    #[test]
    fn a_build_from_flags_downloads_from_the_located_tree() {
        let root = tempfile::tempdir().unwrap();
        publish(
            root.path(),
            "pub/release-63/species.txt",
            GENOMES.as_bytes(),
        );
        let at = |relative: &str| root.path().join(relative).display().to_string();

        let address = |(site, relative): (&Site, String)| site.address(&relative);

        let plant = Published::new(vec![site_at(root.path())], "planta_ficta", "Pf1.0", 63);
        assert_eq!(
            address(plant.gff3().unwrap()),
            at("pub/release-63/plants/gff3/planta_ficta/Planta_ficta.Pf1.0.63.gff3.gz")
        );
        assert_eq!(
            address(plant.protein_fasta().unwrap()),
            at("pub/release-63/plants/fasta/planta_ficta/pep/Planta_ficta.Pf1.0.pep.all.fa.gz")
        );
        let databases = plant.database_urls().unwrap();
        assert_eq!(
            databases.core,
            format!(
                "{}/",
                at("pub/release-63/plants/mysql/planta_ficta_core_63_116_4")
            )
        );
        assert_eq!(
            databases.variation,
            Some(format!(
                "{}/",
                at("pub/release-63/plants/mysql/planta_ficta_variation_63_116_4")
            ))
        );

        let bacterium = Published::new(
            vec![site_at(root.path())],
            "bacterium_fictum_gca_000000001",
            "ASM1v1",
            63,
        );
        let error = bacterium.database_urls().err().unwrap();
        assert!(download::is_not_published(&error), "{error}");
        assert!(
            error
                .to_string()
                .contains("it has no core database of its own"),
            "{error}"
        );
        let unlisted = Published::new(vec![site_at(root.path())], "planta_ignota", "Pi1", 63);
        assert!(download::is_not_published(
            &unlisted.database_urls().err().unwrap()
        ));
    }

    /// A species no list names, called for an assembly accession, is found in
    /// the rapid release under that accession, with its latest gene set, and
    /// refused for the table dumps that tree does not publish.
    #[test]
    fn a_rapid_release_assembly_is_found_by_its_accession_and_refused() {
        let root = tempfile::tempdir().unwrap();
        let tree = "pub/rapid-release/species/Homo_sapiens/GCA_000000001.1/ensembl";
        for gene_set in ["2022_06", "2022_07"] {
            std::fs::create_dir_all(root.path().join(tree).join("geneset").join(gene_set)).unwrap();
        }
        let sites = [site_at(root.path())];
        let species = "homo_sapiens_gca000000001v1";

        let located = locate(&sites, species, "Made-up1.0", 107).unwrap();
        assert_eq!(
            located.layout,
            Layout::RapidRelease {
                organism: "Homo_sapiens".to_string(),
                accession: "GCA_000000001.1".to_string(),
                gene_set: "2022_07".to_string()
            }
        );
        assert_eq!(located.core_db, None);
        let error = located
            .buildable(species, "Made-up1.0", 107)
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            format!(
                "homo_sapiens_gca000000001v1 Made-up1.0 cannot be installed from {}/: {tree}/ publishes no core table dumps; convert Ensembl VEP's cache of it instead (docs/cache-conversion.md), or build from local copies of its files with --gff3, --fasta, --genome-fasta and --codon-table",
                root.path().display()
            )
        );
    }

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, bytes).unwrap();
        encoder.finish().unwrap()
    }

    /// Lays the committed chromosome 21 slice out under `root` as release
    /// 116 of `species` on `assembly` in the tree `layout`, each directory
    /// with its `CHECKSUMS` and none with a species list: the GFF3, the GTF,
    /// the peptides, the core tables as `database` and, under `dna/`, the
    /// committed chromosome 21 reference as the genome.
    fn publish_genome(root: &Path, layout: &Layout, species: &str, assembly: &str, database: &str) {
        let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
        let fixture = manifest.join("tests/fixtures/grch38_chr21");
        let gzipped = |name: &str| gzip(&std::fs::read(fixture.join(name)).unwrap());
        let at = |(dir, name): (String, String)| format!("{dir}/{name}");
        publish(
            root,
            &at(layout.gff3_location(116, assembly, species)),
            &gzipped("slice.gff3"),
        );
        publish(
            root,
            &at(layout.gtf_location(116, assembly, species)),
            &gzipped("slice.gtf"),
        );
        publish(
            root,
            &at(layout.protein_fasta_location(116, assembly, species)),
            &gzipped("peptides.fa"),
        );
        let mysql = layout.mysql_dir(116, database).unwrap();
        for table in CORE_TABLES {
            publish(
                root,
                &format!("{mysql}/{table}.txt.gz"),
                &gzipped(&format!("{table}.txt")),
            );
        }
        let dna = layout.genome_fasta_dirs(116, species).pop().unwrap();
        publish(
            root,
            &format!(
                "{dna}/{}.{assembly}.dna.toplevel.fa.gz",
                download::capitalize_species(species)
            ),
            &std::fs::read(manifest.join("../../tests/golden/116/GRCh37-hgvs/reference.fa.gz"))
                .unwrap(),
        );
    }

    /// A site that cannot be read, laid out under `root`: a directory whose
    /// `CHECKSUMS` of release 116 is no text, so the read of it fails where
    /// an absent file reads as a site that publishes no species list.
    fn unreadable_site(root: &Path) -> Site {
        let release = root.join("pub/release-116");
        std::fs::create_dir_all(&release).unwrap();
        std::fs::write(release.join("CHECKSUMS"), [0xff, 0xfe, 0xfd]).unwrap();
        site_at(root)
    }

    fn options(
        dir_cache: &Path,
        sites: Vec<Site>,
        species: &str,
        assembly: &str,
    ) -> InstallOptions {
        InstallOptions {
            species: species.to_string(),
            assembly: assembly.to_string(),
            release: 116,
            dir_cache: dir_cache.to_path_buf(),
            sites,
            chromosomes: None,
            codon_tables: HashMap::new(),
            prediction_server: None,
        }
    }

    /// GRCh37, and a vertebrate of a release whose tree publishes no species
    /// list, each install from the first site while the second cannot be
    /// read: the GRCh37 tree is the first site's alone, and a release tree
    /// that lists the genome's GFF3 and holds its core database is not
    /// looked for on a later site.
    #[test]
    fn a_tree_without_a_species_list_installs_while_the_second_site_cannot_be_read() {
        for (layout, species, assembly, database) in [
            (
                Layout::Grch37,
                "homo_sapiens",
                "GRCh37",
                "homo_sapiens_core_116_37",
            ),
            (
                Layout::Release,
                "danio_rerio",
                "GRCz11",
                "danio_rerio_core_116_11",
            ),
        ] {
            let root = tempfile::tempdir().unwrap();
            let tree = root.path().join("site");
            publish_genome(&tree, &layout, species, assembly, database);
            let dir_cache = root.path().join("vep");
            let sites = vec![
                site_at(&tree),
                unreadable_site(&root.path().join("unreadable")),
            ];

            let installed = install(&options(&dir_cache, sites, species, assembly))
                .unwrap_or_else(|error| panic!("{species} {assembly}: {error:#}"));
            assert_eq!(
                installed,
                dir_cache.join(species).join(format!("116_{assembly}"))
            );
            assert!(installed.join("info.json").is_file(), "{assembly}");
            assert!(installed.join("transcripts/21").is_dir(), "{assembly}");
            assert!(!dir_cache.join("tmp").exists(), "{assembly}");
        }
    }

    /// A species neither site holds, with one of them unreadable, is refused
    /// before anything is fetched: the error says what the readable tree
    /// lacks and names the site that could not be read, it is no
    /// not-published answer, and the cache root is never created. The lookup
    /// of a species the readable site's list does not name, and one with the
    /// unreadable site first, fail likewise: each error names the unreadable
    /// site and is no not-published answer.
    #[test]
    fn a_species_on_no_site_is_refused_naming_the_site_that_could_not_be_read() {
        let root = tempfile::tempdir().unwrap();
        let tree = root.path().join("site");
        publish_genome(
            &tree,
            &Layout::Release,
            "danio_rerio",
            "GRCz11",
            "danio_rerio_core_116_11",
        );
        let dir_cache = root.path().join("vep");
        let unreadable = unreadable_site(&root.path().join("unreadable"));
        let sites = vec![site_at(&tree), unreadable.clone()];

        let error = install(&options(&dir_cache, sites, "ovis_aries", "Oar1")).unwrap_err();
        assert!(!download::is_not_published(&error));
        let text = format!("{error:#}");
        assert!(
            text.starts_with(&format!(
                "pub/release-116/mysql/ holds no core database of ovis_aries for release 116; \
                 {} could not be read (Failed to read pub/release-116/CHECKSUMS: ",
                unreadable.address("")
            )),
            "{text}"
        );
        assert!(!dir_cache.exists());

        let listing = root.path().join("listing");
        let list = "pub/release-116/species_EnsemblVertebrates.txt";
        publish(&listing, list, VERTEBRATES.as_bytes());
        for (sites, refusal) in [
            (
                [site_at(&listing), unreadable.clone()],
                format!(
                    "ovis_aries is not in the species list of release 116 ({}); ",
                    listing.join(list).display()
                ),
            ),
            (
                [unreadable.clone(), site_at(&tree)],
                "ovis_aries Oar1 was not found for release 116; ".to_string(),
            ),
        ] {
            let error = locate(&sites, "ovis_aries", "Oar1", 116).unwrap_err();
            assert!(!download::is_not_published(&error), "{error:#}");
            let text = format!("{error:#}");
            assert!(
                text.starts_with(&refusal)
                    && text.contains(&format!(
                        "; {} could not be read (Failed to read pub/release-116/CHECKSUMS: ",
                        unreadable.address("")
                    )),
                "{text}"
            );
        }
    }

    /// With several sites the first tree is taken for the genome only when it
    /// lists the genome's GFF3: a release tree holding a core database of
    /// the species and no GFF3 of the assembly, as an old release of the
    /// Ensembl site does beside an Ensembl Genomes release of the same
    /// number, leaves the genome to the site whose species list names it,
    /// and the genome is refused, naming that site, when it cannot be read.
    /// A first site that cannot be read is passed over for the site whose
    /// list names the species. One site alone still stops on a failed read.
    #[test]
    fn a_later_site_is_asked_only_for_a_genome_the_first_tree_does_not_hold() {
        let genomes = tempfile::tempdir().unwrap();
        publish(
            genomes.path(),
            "pub/release-116/species.txt",
            GENOMES
                .replace("planta_ficta_core_63_116_4", "danio_rerio_core_116_169_11")
                .replace("planta_ficta", "danio_rerio")
                .as_bytes(),
        );
        let second = site_at(genomes.path());
        let division = Layout::Division {
            division: "plants".to_string(),
        };

        let holding = tempfile::tempdir().unwrap();
        publish_genome(
            holding.path(),
            &Layout::Release,
            "danio_rerio",
            "GRCz11",
            "danio_rerio_core_116_11",
        );
        let sites = [site_at(holding.path()), second.clone()];
        let located = locate(&sites, "danio_rerio", "GRCz11", 116).unwrap();
        assert!(matches!(&located.site, Site::Dir(dir) if dir == holding.path()));
        assert_eq!(located.layout, Layout::Release);
        assert_eq!(located.core_db.as_deref(), Some("danio_rerio_core_116_11"));

        let older = tempfile::tempdir().unwrap();
        std::fs::create_dir_all(
            older
                .path()
                .join("pub/release-116/mysql/danio_rerio_core_116_9"),
        )
        .unwrap();
        let sites = [site_at(older.path()), second.clone()];
        let located = locate(&sites, "danio_rerio", "GRCz11", 116).unwrap();
        assert!(matches!(&located.site, Site::Dir(dir) if dir == genomes.path()));
        assert_eq!(located.layout, division);

        let unread = tempfile::tempdir().unwrap();
        let unreadable = unreadable_site(unread.path());
        let located = locate(&[unreadable.clone(), second], "danio_rerio", "GRCz11", 116).unwrap();
        assert!(matches!(&located.site, Site::Dir(dir) if dir == genomes.path()));
        assert_eq!(located.layout, division);

        let sites = [site_at(older.path()), unreadable.clone()];
        let error = locate(&sites, "danio_rerio", "GRCz11", 116).unwrap_err();
        assert!(!download::is_not_published(&error));
        let text = format!("{error:#}");
        assert!(
            text.starts_with(&format!(
                "danio_rerio GRCz11 was not found for release 116; \
                 {} could not be read (Failed to read pub/release-116/CHECKSUMS: ",
                unreadable.address("")
            )),
            "{text}"
        );

        let error = locate(&[unreadable], "danio_rerio", "GRCz11", 116).unwrap_err();
        assert!(
            format!("{error:#}").starts_with("Failed to read pub/release-116/CHECKSUMS: "),
            "{error:#}"
        );
    }

    /// A listing settles the assembly's directory: the entry named for the
    /// assembly, else the one named for the last dot-separated part of its
    /// name, else the assembly's own name. Without a listing both names are
    /// tried, the assembly's first, and an assembly without a dot has one.
    #[test]
    fn the_assembly_directory_is_named_for_the_assembly_or_its_last_part() {
        let listed = |names: &[&str]| {
            names
                .iter()
                .map(|name| name.to_string())
                .collect::<Vec<_>>()
        };
        let chicken = "bGalGal1.mat.broiler.GRCg7b";
        assert_eq!(
            assembly_directories(chicken, Some(&listed(&["GRCg6a", "GRCg7b"]))),
            ["GRCg7b"]
        );
        assert_eq!(
            assembly_directories(chicken, Some(&listed(&["GRCg7b", chicken]))),
            [chicken]
        );
        assert_eq!(
            assembly_directories(chicken, Some(&listed(&["GRCg6a"]))),
            [chicken]
        );
        assert_eq!(
            assembly_directories("GRCh38", Some(&listed(&["GRCh37", "GRCh38"]))),
            ["GRCh38"]
        );
        assert_eq!(assembly_directories(chicken, None), [chicken, "GRCg7b"]);
        assert_eq!(assembly_directories("GRCh38", None), ["GRCh38"]);
    }

    /// The regulatory file of a kind is the one of the highest version, the
    /// versions compared as numbers. A file of another kind, one whose kind
    /// only ends as the kind asked for, one of another extension and one
    /// whose version is no number are passed over.
    #[test]
    fn regulatory_file_is_the_highest_version_of_its_kind() {
        let listing = Listing::parse(
            "1 1 Homo_sapiens.GRCh38.regulatory_features.v115.gff3.gz\n\
             1 1 Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz\n\
             1 1 Homo_sapiens.GRCh38.regulatory_features.v9.gff3.gz\n\
             1 1 Homo_sapiens.GRCh38.motif_features.v117.gff3.gz\n\
             1 1 Homo_sapiens.GRCh38.other_regulatory_features.v118.gff3.gz\n\
             1 1 Homo_sapiens.GRCh38.regulatory_features.v119.gff.gz\n\
             1 1 Homo_sapiens.GRCh38.regulatory_features.v120.gff3.gz.tbi\n\
             1 1 Homo_sapiens.GRCh38.regulatory_features.vNext.gff3.gz\n\
             1 1 Homo_sapiens.GRCh38.EMARs.v116.gff.gz\n\
             1 1 README\n",
        );
        assert_eq!(
            regulatory_file(&listing, REGULATORY_FEATURES).as_deref(),
            Some("Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz")
        );
        assert_eq!(
            regulatory_file(&listing, MOTIF_FEATURES).as_deref(),
            Some("Homo_sapiens.GRCh38.motif_features.v117.gff3.gz")
        );
        assert_eq!(regulatory_file(&listing, "EMARs"), None);
        assert_eq!(
            regulatory_file(&Listing::parse("1 1 README\n"), MOTIF_FEATURES),
            None
        );
    }

    /// The activity file of an annotation directory is the one of the highest
    /// version among the names `<anything>.regulatory_activity.v<N>.tsv.gz`.
    /// A file of that kind and another extension, and a file of another
    /// kind, are passed over.
    #[test]
    fn the_activity_file_is_the_highest_version_of_its_kind() {
        let listing = Listing::parse(
            "1 1 Homo_sapiens.GRCh38.regulatory_activity.v114.tsv.gz\n\
             1 1 Homo_sapiens.GRCh38.regulatory_activity.v116.tsv.gz\n\
             1 1 Homo_sapiens.GRCh38.regulatory_activity.v9.tsv.gz\n\
             1 1 Homo_sapiens.GRCh38.regulatory_activity.v117.tsv\n\
             1 1 Homo_sapiens.GRCh38.regulatory_activity.v118.tsv.gz.tbi\n\
             1 1 Homo_sapiens.GRCh38.other_regulatory_activity.v119.tsv.gz\n\
             1 1 Homo_sapiens.GRCh38.regulatory_features.v120.gff3.gz\n",
        );
        assert_eq!(
            versioned_file(&listing, REGULATORY_ACTIVITY, ".tsv.gz").as_deref(),
            Some("Homo_sapiens.GRCh38.regulatory_activity.v116.tsv.gz")
        );
        assert_eq!(
            versioned_file(
                &Listing::parse("1 1 Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz\n"),
                REGULATORY_ACTIVITY,
                ".tsv.gz"
            ),
            None
        );
    }

    /// The funcgen database stands beside the core database, in the tree
    /// that publishes the core database: GRCh37's in the GRCh37 tree. A
    /// genome of a collection, whose core database is not named for it, and
    /// an assembly of the rapid release have none.
    #[test]
    fn the_funcgen_directory_stands_beside_the_core_database() {
        let located = |layout: Layout, core_db: Option<&str>| Located {
            site: Site::Url("https://a.site".to_string()),
            layout,
            core_db: core_db.map(str::to_string),
        };
        assert_eq!(
            located(Layout::Release, Some("homo_sapiens_core_116_38"))
                .funcgen_dir("homo_sapiens", 116)
                .as_deref(),
            Some("pub/release-116/mysql/homo_sapiens_funcgen_116_38")
        );
        assert_eq!(
            located(Layout::Grch37, Some("homo_sapiens_core_116_37"))
                .funcgen_dir("homo_sapiens", 116)
                .as_deref(),
            Some("pub/grch37/release-116/mysql/homo_sapiens_funcgen_116_37")
        );
        let collection = Layout::Collection {
            division: "bacteria".to_string(),
            collection: "bacteria_7_collection".to_string(),
        };
        assert_eq!(
            located(collection, Some("bacteria_7_collection_core_63_116_1"))
                .funcgen_dir("bacterium_fictum_gca_000000001", 63),
            None
        );
        assert_eq!(
            located(Layout::Release, None).funcgen_dir("homo_sapiens", 116),
            None
        );
    }

    /// What an installation reads beside the regulatory GFF3 files is looked
    /// up before anything is fetched, and each file or table the tree does
    /// not publish is held as its reason: the activity file the annotation
    /// listing lacks, the funcgen directory the site lacks, and in that
    /// directory the first table of a group its listing lacks, the factors'
    /// tables being named only once the matrices' are listed. What is found
    /// is fetched for the files it belongs to: the activity file with the
    /// regulatory-feature GFF3, the matrices' and the factors' tables with
    /// the motif-feature GFF3, the regulatory build's table with either.
    #[test]
    fn the_files_of_the_regulatory_cells_are_found_or_held_as_their_reasons() {
        let root = tempfile::tempdir().unwrap();
        let annotation = "pub/release-116/regulation/homo_sapiens/GRCh38/annotation";
        let funcgen = "pub/release-116/mysql/homo_sapiens_funcgen_116_38";
        let named = "Homo_sapiens.GRCh38";
        for kind in [REGULATORY_FEATURES, MOTIF_FEATURES] {
            publish(
                root.path(),
                &format!("{annotation}/{named}.{kind}.v116.gff3.gz"),
                kind.as_bytes(),
            );
        }
        let site = site_at(root.path());
        let located = Located {
            site: site.clone(),
            layout: Layout::Release,
            core_db: Some("homo_sapiens_core_116_38".to_string()),
        };
        let options = options(root.path(), vec![site.clone()], "homo_sapiens", "GRCh38");
        let files = || match RegulatoryFiles::find(&located, &options).unwrap() {
            Part::Published(files) => files,
            Part::Unpublished(reason) => panic!("{reason}"),
        };
        let not_found = || -> Vec<String> {
            RegulatoryCells::find(&located, &options, &files())
                .unwrap()
                .not_found()
                .into_iter()
                .map(|(what, reason)| format!("{what}: {reason}"))
                .collect()
        };

        assert_eq!(
            not_found(),
            [
                format!(
                    "cell-type activity: {annotation}/CHECKSUMS lists no regulatory_activity file"
                ),
                format!("table of the funcgen database: {funcgen}/CHECKSUMS is not on the site"),
            ]
        );

        let activity = format!("{named}.regulatory_activity.v116.tsv.gz");
        publish(
            root.path(),
            &format!("{annotation}/{activity}"),
            b"activity",
        );
        publish(
            root.path(),
            &format!("{funcgen}/binding_matrix.txt.gz"),
            b"m",
        );
        assert_eq!(
            not_found(),
            [
                format!("binding matrix table: {funcgen}/CHECKSUMS lists no binding_matrix_frequencies table"),
                format!("regulatory build table: {funcgen}/CHECKSUMS lists no regulatory_build table"),
            ]
        );

        publish(
            root.path(),
            &format!("{funcgen}/binding_matrix_frequencies.txt"),
            b"f",
        );
        publish(
            root.path(),
            &format!("{funcgen}/regulatory_build.txt.gz"),
            b"b",
        );
        publish(
            root.path(),
            &format!("{funcgen}/binding_matrix_transcription_factor_complex.txt.gz"),
            b"l",
        );
        assert_eq!(
            not_found(),
            [format!(
                "transcription factor complex table: {funcgen}/CHECKSUMS lists no transcription_factor_complex table"
            )]
        );

        let work = tempfile::tempdir().unwrap();
        let fetched = |files: &RegulatoryFiles<'_>, to: &str| {
            let cells = RegulatoryCells::find(&located, &options, files).unwrap();
            let (regulatory_dir, dumps_dir) = (
                work.path().join(to).join("regulatory"),
                work.path().join(to).join("funcgen_dumps"),
            );
            let mut inputs = regulatory::RegulatoryInputs::default();
            cells
                .fetch(files, &regulatory_dir, &dumps_dir, &mut inputs)
                .unwrap();
            let mut tables: Vec<String> = match std::fs::read_dir(&dumps_dir) {
                Ok(entries) => entries
                    .map(|entry| entry.unwrap().file_name().into_string().unwrap())
                    .collect(),
                Err(_) => Vec::new(),
            };
            tables.sort();
            assert_eq!(inputs.funcgen_dumps.is_some(), !tables.is_empty());
            (
                inputs
                    .activity
                    .map(|path| path.strip_prefix(work.path()).unwrap().to_path_buf()),
                tables,
            )
        };
        let without_complexes = [
            "binding_matrix.txt.gz",
            "binding_matrix_frequencies.txt",
            "regulatory_build.txt.gz",
        ];
        assert_eq!(
            fetched(&files(), "one_complex_table"),
            (
                Some(Path::new("one_complex_table/regulatory").join(&activity)),
                without_complexes.map(String::from).to_vec()
            )
        );

        publish(
            root.path(),
            &format!("{funcgen}/transcription_factor_complex.txt.gz"),
            b"c",
        );
        assert_eq!(not_found(), Vec::<String>::new());
        assert_eq!(
            fetched(&files(), "every_table"),
            (
                Some(Path::new("every_table/regulatory").join(&activity)),
                [
                    "binding_matrix.txt.gz",
                    "binding_matrix_frequencies.txt",
                    "binding_matrix_transcription_factor_complex.txt.gz",
                    "regulatory_build.txt.gz",
                    "transcription_factor_complex.txt.gz",
                ]
                .map(String::from)
                .to_vec()
            )
        );
        let no_motif_file = RegulatoryFiles {
            motifs: None,
            ..files()
        };
        assert_eq!(
            fetched(&no_motif_file, "no_motif_file"),
            (
                Some(Path::new("no_motif_file/regulatory").join(&activity)),
                vec!["regulatory_build.txt.gz".to_string()]
            )
        );
        let no_feature_file = RegulatoryFiles {
            features: None,
            ..files()
        };
        assert_eq!(fetched(&no_feature_file, "no_feature_file").0, None);
    }

    /// The directory `relative` of `site`, with a listing that names no file.
    fn unlisted<'a>(site: &'a Site, relative: &str) -> Directory<'a> {
        Directory {
            site,
            relative: relative.to_string(),
            listing: Listing::default(),
        }
    }

    /// The line of an installation that left no part out names the
    /// 1000 Genomes frequencies for `homo_sapiens` and for no other species,
    /// whose variation VCF carries none.
    #[test]
    fn the_line_of_every_part_names_the_frequencies_for_human_alone() {
        assert_eq!(
            every_part_line("homo_sapiens", "GRCh38"),
            "The cache of homo_sapiens GRCh38 holds each part of Ensembl VEP's own cache that a \
             file of the release carries: known variants with their 1000 Genomes frequencies, \
             SIFT and PolyPhen predictions, regulatory features and motif features"
        );
        assert_eq!(
            every_part_line("danio_rerio", "GRCz11"),
            "The cache of danio_rerio GRCz11 holds each part of Ensembl VEP's own cache that a \
             file of the release carries: known variants, SIFT and PolyPhen predictions, \
             regulatory features and motif features"
        );
    }

    /// A built cache's record that names `names` under `left_out`, with the
    /// SIFT marker and variation columns that hold no `AF`.
    fn record_naming(names: &[&str]) -> vep_core::cache_info::CacheInfo {
        let mut info = vep_core::cache_info::CacheInfo::default();
        info.sift = Some("b".to_string());
        info.variation_cols = vec!["chr".to_string(), "variation_name".to_string()];
        info.left_out = names.iter().map(|name| name.to_string()).collect();
        info
    }

    /// The lookups of the cells of the known variants, each a reason, `None`
    /// for one that was found: the tables of the PubMed identifiers, of the
    /// synonyms and of the failed flag, and the failed-variants GVF.
    fn variant_cell_lookups<'a>(
        site: &'a Site,
        reasons: [Option<&str>; 4],
    ) -> VariantCellSources<'a> {
        let [pubmed, synonyms, failed_variation, failed_gvf] = reasons;
        let table = |reason: Option<&str>| reason.map_or(Ok(()), |reason| Err(reason.to_string()));
        VariantCellSources {
            pubmed: table(pubmed),
            synonyms: table(synonyms),
            failed_variation: table(failed_variation),
            failed_gvf: match failed_gvf {
                Some(reason) => Err(reason.to_string()),
                None => Ok(Input {
                    directory: unlisted(site, "pub/release-116/variation/gvf/homo_sapiens"),
                    name: "homo_sapiens_failed.gvf.gz".to_string(),
                }),
            },
        }
    }

    /// The route every entry of a cell ends on.
    const ROUTE: &str =
        "a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";
    const NO_VARIATION_TABLE: &str =
        "pub/release-116/mysql/homo_sapiens_variation_116_38/CHECKSUMS lists no variation table";
    const NO_GVF_LISTING: &str =
        "pub/release-116/variation/gvf/homo_sapiens/CHECKSUMS is not on the site";
    const NO_ACTIVITY_FILE: &str = "pub/release-116/regulation/homo_sapiens/GRCh38/annotation/CHECKSUMS lists no regulatory_activity file";
    const NO_FUNCGEN_LISTING: &str =
        "pub/release-116/mysql/homo_sapiens_funcgen_116_38/CHECKSUMS is not on the site";

    /// What a built cache lacks of a part whose files were found is read
    /// from its record. A `homo_sapiens` cache whose variation columns hold
    /// no `AF` lacks the 1000 Genomes frequencies, naming the file the
    /// variation directory does not list; under `--include-predictions` a
    /// tool the cache holds no matrix of is left out naming Ensembl's MySQL
    /// server; a regulatory part read without a motif file lacks the motif
    /// features. The cells follow the names of the record's `left_out`:
    /// the gnomAD frequencies, which no file of the release carries, then
    /// the PubMed identifiers and the synonyms in one entry with the reason
    /// they share, the failed variants with the reason of the GVF, and the
    /// activity of the regulatory features with the reason of its file; a
    /// tree read without a regulatory-feature file lacks the regulatory
    /// features and the cells of the motif features alone. With no file of
    /// any part, the three reasons the site gave stand in the order of the
    /// parts and no cell follows.
    #[test]
    fn what_a_built_cache_lacks_is_read_from_its_record() {
        let root = tempfile::tempdir().unwrap();
        let site = site_at(root.path());
        let vcfs = VariationVcfs {
            directory: unlisted(&site, "pub/release-116/variation/vcf/homo_sapiens"),
            germline: Vec::new(),
        };
        let variant_cells = variant_cell_lookups(
            &site,
            [
                Some(NO_VARIATION_TABLE),
                Some(NO_VARIATION_TABLE),
                Some(NO_VARIATION_TABLE),
                Some(NO_GVF_LISTING),
            ],
        );
        let regulatory_cells = RegulatoryCells {
            activity: Err(NO_ACTIVITY_FILE.to_string()),
            funcgen: Err(NO_FUNCGEN_LISTING.to_string()),
        };
        let regulatory = RegulatoryFiles {
            directory: unlisted(
                &site,
                "pub/release-116/regulation/homo_sapiens/GRCh38/annotation",
            ),
            features: Some("a regulatory-feature file".to_string()),
            motifs: None,
        };
        let options = InstallOptions {
            prediction_server: Some(MysqlServer {
                host: "a server".to_string(),
                port: 3306,
            }),
            ..options(root.path(), vec![site.clone()], "homo_sapiens", "GRCh38")
        };
        let record = record_naming(&[
            "thousand_genomes_frequencies",
            "gnomad_frequencies",
            "pubmed",
            "var_synonyms",
            "failed_variants",
            "cell_type_activity",
        ]);
        let left = |part: &str, reason: &str| LeftOut::new(part, reason.to_string());

        let missing = left_out(
            &options,
            &record,
            &Unpublished::default(),
            &[],
            &Brought {
                vcfs: Some(&vcfs),
                variant_cells: Some(&variant_cells),
                regulatory_files: Some(&regulatory),
                regulatory_cells: Some(&regulatory_cells),
                ..Default::default()
            },
        );
        assert_eq!(
            missing.parts,
            [
                left(
                    "1000 Genomes frequencies",
                    "pub/release-116/variation/vcf/homo_sapiens/CHECKSUMS lists no 1000GENOMES-phase_3.vcf.gz; a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)"
                ),
                left(
                    "PolyPhen predictions",
                    "Ensembl's MySQL server holds no PolyPhen matrix of a protein of this genome"
                ),
                left(
                    "motif features",
                    "pub/release-116/regulation/homo_sapiens/GRCh38/annotation/CHECKSUMS lists no motif_features GFF3"
                ),
            ]
        );
        assert_eq!(
            missing.cells,
            [
                left(
                    "gnomAD frequencies",
                    &format!("no file of the release carries them; {ROUTE}")
                ),
                left(
                    "PubMed identifiers or variant synonyms",
                    &format!("{NO_VARIATION_TABLE}; {ROUTE}")
                ),
                left("failed variants", &format!("{NO_GVF_LISTING}; {ROUTE}")),
                left(
                    "cell-type activity of regulatory features",
                    &format!("{NO_ACTIVITY_FILE}; {ROUTE}")
                ),
            ]
        );

        let motifs_alone = RegulatoryFiles {
            features: None,
            motifs: Some("a motif-feature file".to_string()),
            ..regulatory
        };
        let missing = left_out(
            &options,
            &record_naming(&["motif_scores"]),
            &Unpublished::default(),
            &[],
            &Brought {
                regulatory_files: Some(&motifs_alone),
                regulatory_cells: Some(&regulatory_cells),
                ..Default::default()
            },
        );
        assert_eq!(
            missing.parts,
            [
                left(
                    "PolyPhen predictions",
                    "Ensembl's MySQL server holds no PolyPhen matrix of a protein of this genome"
                ),
                left(
                    "regulatory features",
                    "pub/release-116/regulation/homo_sapiens/GRCh38/annotation/CHECKSUMS lists no regulatory_features GFF3"
                ),
            ]
        );
        assert_eq!(
            missing.cells,
            [left(
                "HIGH_INF_POS or MOTIF_SCORE_CHANGE of motif features",
                &format!("{NO_FUNCGEN_LISTING}; {ROUTE}")
            )]
        );

        let unpublished = Unpublished {
            variants: Some("no variation file".to_string()),
            matrices: Some("no table".to_string()),
            regulatory: Some("no annotation file".to_string()),
        };
        let missing = left_out(&options, &record, &unpublished, &[], &Brought::default());
        assert_eq!(
            missing,
            Missing {
                parts: vec![
                    left("known variants or their frequencies", "no variation file"),
                    left("SIFT or PolyPhen predictions", "no table"),
                    left(
                        "regulatory features or motif features",
                        "no annotation file"
                    ),
                ],
                sequences: Vec::new(),
                cells: Vec::new(),
            }
        );
    }

    /// The sequences of a built cache that hold no known variant are those
    /// of its record's `contig_order`, in that order, with no variation file
    /// or with one of no record; a variation file of a sequence the record
    /// does not list counts for nothing.
    #[test]
    fn sequences_without_known_variants_are_read_from_the_built_cache() {
        use std::io::Write;

        let cache = tempfile::tempdir().unwrap();
        let variation_file = |sequence: &str, records: &str| {
            let dir = cache.path().join("variations").join(sequence);
            std::fs::create_dir_all(&dir).unwrap();
            let mut file =
                noodles::bgzf::io::Writer::new(File::create(dir.join("all_vars.gz")).unwrap());
            file.write_all(records.as_bytes()).unwrap();
            file.finish().unwrap();
        };
        variation_file("1", "1\trs1\t.\t.\t100\t.\tA/G\n");
        variation_file("MT", "");
        variation_file("off_the_model", "off_the_model\trs2\t.\t.\t5\t.\tA/G\n");
        let mut record = vep_core::cache_info::CacheInfo::default();
        record.contig_order = ["1", "X", "MT"].map(String::from).to_vec();

        assert_eq!(
            sequences_without_known_variants(cache.path(), &record).unwrap(),
            ["X", "MT"]
        );
        let no_sequence = vep_core::cache_info::CacheInfo::default();
        assert!(sequences_without_known_variants(cache.path(), &no_sequence)
            .unwrap()
            .is_empty());
    }

    /// The sequences of the cache that hold no known variant stand in one
    /// entry between the parts and the cells: one sequence by its name, up to
    /// ten by their count and names, more by their count, the first ten
    /// names and the count of the rest. The entry stands only for an
    /// installation that brings the variation files, and for no cache that
    /// holds a known variant on each of its sequences.
    #[test]
    fn sequences_without_a_known_variant_are_named_ten_at_most() {
        let root = tempfile::tempdir().unwrap();
        let site = site_at(root.path());
        let vcfs = VariationVcfs {
            directory: unlisted(&site, "pub/release-116/variation/vcf/danio_rerio"),
            germline: Vec::new(),
        };
        let options = options(root.path(), vec![site.clone()], "danio_rerio", "GRCz11");
        let mut record = vep_core::cache_info::CacheInfo::default();
        record.left_out = ["pubmed", "var_synonyms", "failed_variants"]
            .map(String::from)
            .to_vec();
        const NO_DATABASE: &str =
            "no variation database is published for this release and assembly";
        let variant_cells = variant_cell_lookups(&site, [Some(NO_DATABASE); 4]);
        let unpublished = Unpublished {
            matrices: Some("no table".to_string()),
            regulatory: Some("no annotation file".to_string()),
            ..Default::default()
        };
        let named = |count: usize| {
            let without: Vec<String> = (1..=count).map(|n| format!("s{n}")).collect();
            left_out(
                &options,
                &record,
                &unpublished,
                &without,
                &Brought {
                    vcfs: Some(&vcfs),
                    variant_cells: Some(&variant_cells),
                    ..Default::default()
                },
            )
        };
        const REASON: &str = "the installation read no variation record there, in a variation VCF or in the variation_feature table; a conversion of Ensembl VEP's cache carries any known variant that cache holds there (docs/cache-conversion.md)";
        let entry = |part: &str| vec![LeftOut::new(part, REASON.to_string())];

        assert_eq!(named(1).sequences, entry("known variants on sequence s1"));
        assert_eq!(
            named(2).sequences,
            entry("known variants on 2 sequences (s1, s2)")
        );
        assert_eq!(
            named(10).sequences,
            entry("known variants on 10 sequences (s1, s2, s3, s4, s5, s6, s7, s8, s9, s10)")
        );
        assert_eq!(
            named(12).sequences,
            entry(
                "known variants on 12 sequences (s1, s2, s3, s4, s5, s6, s7, s8, s9, s10 and 2 more)"
            )
        );
        let brought = named(12);
        assert_eq!(brought.parts.len(), 2);
        assert_eq!(
            brought.cells,
            [LeftOut::new(
                "PubMed identifiers, variant synonyms or failed variants",
                format!("{NO_DATABASE}; {ROUTE}")
            )]
        );

        assert!(named(0).sequences.is_empty());
        let no_file = Unpublished {
            variants: Some("no variation file".to_string()),
            ..Default::default()
        };
        let without = ["s1".to_string()];
        let missing = left_out(&options, &record, &no_file, &without, &Brought::default());
        assert!(missing.sequences.is_empty());
    }

    /// The cells of the known variants an installation names are those the
    /// built cache's record names, each with the reason its lookup kept and
    /// the cells of one reason in one entry, in the order PubMed
    /// identifiers, variant synonyms, failed variants. Three reasons draw
    /// three entries; a cell the record does not name draws none, whatever
    /// its lookup kept, as the failed variants a `failed_variation` table
    /// flagged where the site has no GVF. The gnomAD frequencies stand
    /// first, for `homo_sapiens` alone.
    #[test]
    fn the_cells_of_the_known_variants_are_named_by_the_reason_each_lookup_kept() {
        let root = tempfile::tempdir().unwrap();
        let site = site_at(root.path());
        let vcfs = VariationVcfs {
            directory: unlisted(&site, "pub/release-116/variation/vcf/homo_sapiens"),
            germline: Vec::new(),
        };
        let unpublished = Unpublished {
            matrices: Some("no table".to_string()),
            regulatory: Some("no annotation file".to_string()),
            ..Default::default()
        };
        let cells = |species: &str, names: &[&str], reasons: [Option<&str>; 4]| {
            let lookups = variant_cell_lookups(&site, reasons);
            left_out(
                &options(root.path(), vec![site.clone()], species, "an assembly"),
                &record_naming(names),
                &unpublished,
                &[],
                &Brought {
                    vcfs: Some(&vcfs),
                    variant_cells: Some(&lookups),
                    ..Default::default()
                },
            )
            .cells
        };
        let entry = |part: &str, reason: &str| LeftOut::new(part, format!("{reason}; {ROUTE}"));
        let gnomad = LeftOut::new(
            "gnomAD frequencies",
            format!("no file of the release carries them; {ROUTE}"),
        );
        let three = ["pubmed", "var_synonyms", "failed_variants"];
        const NO_CITATION_TABLE: &str = "a database/CHECKSUMS lists no variation_citation table";
        const NO_SYNONYM_TABLE: &str = "a database/CHECKSUMS lists no variation_synonym table";

        assert_eq!(
            cells(
                "homo_sapiens",
                &three,
                [
                    Some(NO_CITATION_TABLE),
                    Some(NO_SYNONYM_TABLE),
                    Some(NO_SYNONYM_TABLE),
                    Some(NO_GVF_LISTING)
                ]
            ),
            [
                gnomad.clone(),
                entry("PubMed identifiers", NO_CITATION_TABLE),
                entry("variant synonyms", NO_SYNONYM_TABLE),
                entry("failed variants", NO_GVF_LISTING),
            ]
        );
        assert_eq!(
            cells(
                "homo_sapiens",
                &three,
                [
                    Some(NO_GVF_LISTING),
                    Some(NO_VARIATION_TABLE),
                    Some(NO_VARIATION_TABLE),
                    Some(NO_GVF_LISTING)
                ]
            ),
            [
                gnomad.clone(),
                entry("PubMed identifiers or failed variants", NO_GVF_LISTING),
                entry("variant synonyms", NO_VARIATION_TABLE),
            ]
        );
        assert_eq!(
            cells(
                "danio_rerio",
                &["var_synonyms"],
                [None, Some(NO_SYNONYM_TABLE), None, Some(NO_GVF_LISTING)]
            ),
            [entry("variant synonyms", NO_SYNONYM_TABLE)]
        );
        assert_eq!(
            cells("homo_sapiens", &[], [None, None, None, None]),
            [gnomad]
        );
        assert!(cells("danio_rerio", &[], [None, None, None, None]).is_empty());
    }

    /// The cells of the regulatory tree an installation names are those the
    /// built cache's record names: the activity with the reason of its
    /// file's lookup and the matrix counts with that of the funcgen
    /// database's, the two in one entry where the reason is one. Where the
    /// matrices' tables were read, the matrix counts take their reason from
    /// what the build found short: a matrix without counts, a motif without
    /// its sequence, or both. A cell the record does not name draws no
    /// entry.
    #[test]
    fn the_cells_of_the_regulatory_tree_are_named_by_the_reason_of_each() {
        let root = tempfile::tempdir().unwrap();
        let site = site_at(root.path());
        let files = RegulatoryFiles {
            directory: unlisted(
                &site,
                "pub/release-116/regulation/homo_sapiens/GRCh38/annotation",
            ),
            features: Some("a regulatory-feature file".to_string()),
            motifs: Some("a motif-feature file".to_string()),
        };
        let options = options(root.path(), vec![site.clone()], "homo_sapiens", "GRCh38");
        let unpublished = Unpublished {
            variants: Some("no variation file".to_string()),
            matrices: Some("no table".to_string()),
            ..Default::default()
        };
        let cells = |names: &[&str],
                     lookups: &RegulatoryCells<'_>,
                     tree: Option<&regulatory::RegulatoryStats>| {
            left_out(
                &options,
                &record_naming(names),
                &unpublished,
                &[],
                &Brought {
                    regulatory_files: Some(&files),
                    regulatory_cells: Some(lookups),
                    regulatory_tree: tree,
                    ..Default::default()
                },
            )
            .cells
        };
        let entry = |part: &str, reason: &str| LeftOut::new(part, format!("{reason}; {ROUTE}"));
        let both = ["cell_type_activity", "motif_scores"];
        const ACTIVITY: &str = "cell-type activity of regulatory features";
        const SCORES: &str = "HIGH_INF_POS or MOTIF_SCORE_CHANGE of motif features";

        let neither = RegulatoryCells {
            activity: Err(NO_ACTIVITY_FILE.to_string()),
            funcgen: Err(NO_FUNCGEN_LISTING.to_string()),
        };
        assert_eq!(
            cells(&both, &neither, None),
            [
                entry(ACTIVITY, NO_ACTIVITY_FILE),
                entry(SCORES, NO_FUNCGEN_LISTING)
            ]
        );
        assert_eq!(
            cells(&["motif_scores"], &neither, None),
            [entry(SCORES, NO_FUNCGEN_LISTING)]
        );
        assert!(cells(&[], &neither, None).is_empty());

        let one_reason = RegulatoryCells {
            activity: Err("one reason".to_string()),
            funcgen: Err("one reason".to_string()),
        };
        assert_eq!(
            cells(&both, &one_reason, None),
            [entry(
                "cell-type activity of regulatory features, HIGH_INF_POS or MOTIF_SCORE_CHANGE \
                 of motif features",
                "one reason"
            )]
        );

        let funcgen = "pub/release-116/mysql/homo_sapiens_funcgen_116_38";
        let tables = |matrices: Lookup<()>| RegulatoryCells {
            activity: Ok("an activity file".to_string()),
            funcgen: Ok(FuncgenTables {
                directory: unlisted(&site, funcgen),
                matrices,
                complexes: Ok(()),
                regulatory_build: Ok(()),
            }),
        };
        let no_matrix_table = format!("{funcgen}/CHECKSUMS lists no binding_matrix table");
        assert_eq!(
            cells(&both, &tables(Err(no_matrix_table.clone())), None),
            [entry(SCORES, &no_matrix_table)]
        );
        let tree = |counts: bool, sequences: bool| regulatory::RegulatoryStats {
            activity_on_every_feature: true,
            counts_on_every_matrix: counts,
            sequence_on_every_motif: sequences,
            ..Default::default()
        };
        let no_counts = format!("{funcgen} holds no counts of a matrix a motif feature names");
        let no_span = "the genome FASTA does not hold the span of a motif feature";
        assert_eq!(
            cells(&["motif_scores"], &tables(Ok(())), Some(&tree(false, true))),
            [entry(SCORES, &no_counts)]
        );
        assert_eq!(
            cells(&["motif_scores"], &tables(Ok(())), Some(&tree(true, false))),
            [entry(SCORES, no_span)]
        );
        assert_eq!(
            cells(
                &["motif_scores"],
                &tables(Ok(())),
                Some(&tree(false, false))
            ),
            [entry(SCORES, &format!("{no_counts} and {no_span}"))]
        );
        assert!(cells(&[], &tables(Ok(())), Some(&tree(true, true))).is_empty());
    }

    /// Reads the Ensembl site: the release's species list names a core
    /// database the assembly name does not give, the GRCh37 tree's `mysql/`
    /// index page names its own, and every input of the GRCh38 installation is
    /// in its directory's listing.
    #[test]
    #[ignore = "reads https://ftp.ensembl.org"]
    fn live_ensembl_release_names_its_databases() {
        let sites = Site::defaults();
        let rat = locate(&sites, "rattus_norvegicus", "GRCr8", 116).unwrap();
        assert!(matches!(&rat.site, Site::Url(url) if url == download::ENSEMBL_SITE));
        assert_eq!(rat.layout, Layout::Release);
        assert_eq!(rat.core_db.as_deref(), Some("rattus_norvegicus_core_116_1"));

        let grch37 = locate(&sites, "homo_sapiens", "GRCh37", 116).unwrap();
        assert_eq!(grch37.layout, Layout::Grch37);
        assert_eq!(grch37.core_db.as_deref(), Some("homo_sapiens_core_116_37"));

        let human = locate(&sites, "homo_sapiens", "GRCh38", 116).unwrap();
        assert_eq!(human.core_db.as_deref(), Some("homo_sapiens_core_116_38"));
        let inputs = Inputs::find(&human, "homo_sapiens", "GRCh38", 116).unwrap();
        assert_eq!(inputs.gff3.name, "Homo_sapiens.GRCh38.116.gff3.gz");
        assert_eq!(
            inputs.core.relative,
            "pub/release-116/mysql/homo_sapiens_core_116_38"
        );
        assert!(inputs.variation.is_some());
        assert!(inputs.genome.indexed);
    }

    /// Reads the Ensembl site's listings, and no file they name: for
    /// `homo_sapiens` GRCh38 the variation database lists the four tables
    /// the matrices are read from, the regulatory annotation directory lists
    /// a regulatory-feature and a motif-feature GFF3 of the release and the
    /// activity file beside them, the funcgen database lists every table an
    /// installation reads from it, and the variation directory lists the
    /// chromosome 21 file and the 1000 Genomes file.
    #[test]
    #[ignore = "reads https://ftp.ensembl.org"]
    fn live_ensembl_release_lists_the_files_of_each_part_an_installation_brings() {
        let (species, assembly) = ("homo_sapiens", "GRCh38");
        let located = locate(&Site::defaults(), species, assembly, 116).unwrap();
        let options = InstallOptions {
            chromosomes: Some(vec!["21".to_string()]),
            ..options(Path::new("unused"), Site::defaults(), species, assembly)
        };

        let inputs = Inputs::find(&located, species, assembly, 116).unwrap();
        let variation = inputs.variation.as_ref().expect("a variation database");
        assert_eq!(
            variation.relative,
            "pub/release-116/mysql/homo_sapiens_variation_116_38"
        );
        for table in predictions::PREDICTION_TABLES {
            assert!(table_name(variation, table).is_some(), "{table}");
        }
        assert!(matches!(
            prediction_tables(Some(variation), &options),
            Part::Published(())
        ));

        let files = match RegulatoryFiles::find(&located, &options).unwrap() {
            Part::Published(files) => files,
            Part::Unpublished(reason) => panic!("{reason}"),
        };
        assert_eq!(
            files.directory.relative,
            "pub/release-116/regulation/homo_sapiens/GRCh38/annotation"
        );
        assert_eq!(
            files.features.as_deref(),
            Some("Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz")
        );
        assert_eq!(
            files.motifs.as_deref(),
            Some("Homo_sapiens.GRCh38.motif_features.v116.gff3.gz")
        );
        let cells = RegulatoryCells::find(&located, &options, &files).unwrap();
        assert_eq!(
            cells.activity.as_deref(),
            Ok("Homo_sapiens.GRCh38.regulatory_activity.v116.tsv.gz")
        );
        let funcgen = cells.funcgen.as_ref().expect("a funcgen database");
        assert_eq!(
            funcgen.directory.relative,
            "pub/release-116/mysql/homo_sapiens_funcgen_116_38"
        );
        assert_eq!(cells.not_found(), []);

        let vcfs = match VariationVcfs::find(&located, &options).unwrap() {
            Part::Published(vcfs) => vcfs,
            Part::Unpublished(reason) => panic!("{reason}"),
        };
        assert_eq!(
            vcfs.directory.relative,
            "pub/release-116/variation/vcf/homo_sapiens"
        );
        assert_eq!(vcfs.germline, ["homo_sapiens-chr21.vcf.gz"]);
        assert!(vcfs
            .directory
            .listing
            .contains(download::THOUSAND_GENOMES_FILE));
    }

    /// What an installation reads for the `pubmed`, `var_synonyms` and
    /// `failed` cells is looked up before anything is fetched, and each
    /// source the tree does not publish is held as its reason: no variation
    /// database, the first table of a cell its listing lacks, the GVF
    /// directory the site lacks, the failed-variants GVF its listing lacks.
    /// What is found is fetched: the tables of each cell whose every table
    /// is listed, each once and under the name the listing gives it, and the
    /// GVF; where nothing is found there are no inputs.
    #[test]
    fn the_files_of_the_variant_cells_are_found_or_held_as_their_reasons() {
        let root = tempfile::tempdir().unwrap();
        let database = "pub/release-116/mysql/homo_sapiens_variation_116_38";
        let gvf = "pub/release-116/variation/gvf/homo_sapiens";
        let site = site_at(root.path());
        let located = Located {
            site: site.clone(),
            layout: Layout::Release,
            core_db: Some("homo_sapiens_core_116_38".to_string()),
        };
        let options = options(root.path(), vec![site.clone()], "homo_sapiens", "GRCh38");
        let work = tempfile::tempdir().unwrap();
        // The reasons, and the inputs fetched into a directory of `work`:
        // the tables placed, and which inputs are named.
        let found = |to: &str| -> (Vec<String>, Vec<String>, Vec<&'static str>) {
            let variation = optional_directory(&site, database).unwrap();
            let sources = VariantCellSources::find(&located, &options, variation.as_ref()).unwrap();
            let not_found = sources
                .not_found()
                .into_iter()
                .map(|(what, reason)| format!("{what}: {reason}"))
                .collect();
            let work_dir = work.path().join(to);
            let inputs = sources.fetch(variation.as_ref(), &work_dir).unwrap();
            let mut tables: Vec<String> = match std::fs::read_dir(work_dir.join("variation_dumps"))
            {
                Ok(entries) => entries
                    .map(|entry| entry.unwrap().file_name().into_string().unwrap())
                    .collect(),
                Err(_) => Vec::new(),
            };
            tables.sort();
            let named = match &inputs {
                Some(inputs) => {
                    assert_eq!(inputs.scratch, work_dir.join("variant_cells"));
                    assert_eq!(
                        inputs.failed_gvf.as_deref(),
                        inputs
                            .failed_gvf
                            .as_ref()
                            .map(|_| work_dir.join("variation_gvf/homo_sapiens_failed.gvf.gz"))
                            .as_deref()
                    );
                    [
                        ("variation", &inputs.variation),
                        ("variation_citation", &inputs.variation_citation),
                        ("publication", &inputs.publication),
                        ("variation_synonym", &inputs.variation_synonym),
                        ("source", &inputs.source),
                        ("failed_variation", &inputs.failed_variation),
                        ("failed_gvf", &inputs.failed_gvf),
                    ]
                    .into_iter()
                    .filter(|(_, path)| path.as_ref().is_some_and(|path| path.is_file()))
                    .map(|(name, _)| name)
                    .collect()
                }
                None => Vec::new(),
            };
            (not_found, tables, named)
        };
        let none: Vec<String> = Vec::new();

        let no_database = "no variation database is published for this release and assembly";
        assert_eq!(
            found("nothing"),
            (
                vec![
                    format!("PubMed identifier: {no_database}"),
                    format!("variant synonym: {no_database}"),
                    format!("failed_variation table: {no_database}"),
                    format!("failed-variants GVF: {gvf}/CHECKSUMS is not on the site"),
                ],
                none.clone(),
                Vec::new()
            )
        );

        // A database without `variation`, two of its tables uncompressed, and
        // a GVF directory without the failed variants' file.
        for table in ["source.txt", "publication.txt", "failed_variation.txt"] {
            publish(
                root.path(),
                &format!("{database}/{table}"),
                table.as_bytes(),
            );
        }
        publish(root.path(), &format!("{gvf}/homo_sapiens.gvf.gz"), b"g");
        let no_table = format!("{database}/CHECKSUMS lists no variation table");
        assert_eq!(
            found("no_variation_table"),
            (
                vec![
                    format!("PubMed identifier: {no_table}"),
                    format!("variant synonym: {no_table}"),
                    format!("failed_variation table: {no_table}"),
                    format!(
                        "failed-variants GVF: {gvf}/CHECKSUMS lists no homo_sapiens_failed.gvf.gz"
                    ),
                ],
                none.clone(),
                Vec::new()
            )
        );

        publish(root.path(), &format!("{database}/variation.txt.gz"), b"v");
        publish(
            root.path(),
            &format!("{gvf}/homo_sapiens_failed.gvf.gz"),
            b"f",
        );
        assert_eq!(
            found("failed_alone"),
            (
                vec![
                    format!(
                        "PubMed identifier: {database}/CHECKSUMS lists no variation_citation table"
                    ),
                    format!(
                        "variant synonym: {database}/CHECKSUMS lists no variation_synonym table"
                    ),
                ],
                vec![
                    "failed_variation.txt".to_string(),
                    "variation.txt.gz".to_string()
                ],
                vec!["variation", "failed_variation", "failed_gvf"]
            )
        );

        for table in ["variation_citation", "variation_synonym"] {
            publish(
                root.path(),
                &format!("{database}/{table}.txt.gz"),
                table.as_bytes(),
            );
        }
        assert_eq!(
            found("every_source"),
            (
                none,
                [
                    "failed_variation.txt",
                    "publication.txt",
                    "source.txt",
                    "variation.txt.gz",
                    "variation_citation.txt.gz",
                    "variation_synonym.txt.gz",
                ]
                .map(String::from)
                .to_vec(),
                vec![
                    "variation",
                    "variation_citation",
                    "publication",
                    "variation_synonym",
                    "source",
                    "failed_variation",
                    "failed_gvf"
                ]
            )
        );
    }

    /// The tables read for a sequence the one genome-wide variation VCF
    /// holds no record on are looked up in the variation database's listing
    /// before anything is fetched, and held as the reason where the site
    /// publishes no database or the listing lacks `variation_feature` or
    /// `seq_region`. Found, they are placed only when the build asks: each
    /// under the name the listing gives it, `attrib` with them where the
    /// listing names it.
    #[test]
    fn the_feature_tables_are_found_or_held_as_their_reason_and_fetched_when_asked() {
        let root = tempfile::tempdir().unwrap();
        let database = "pub/release-116/mysql/danio_rerio_variation_116_11";
        let site = site_at(root.path());
        let found = || {
            let variation = optional_directory(&site, database).unwrap();
            SequenceFeatureTables::find(variation.as_ref())
        };
        let reason = || found().err();

        assert_eq!(
            reason().as_deref(),
            Some("no variation database is published for this release and assembly")
        );
        publish(root.path(), &format!("{database}/seq_region.txt.gz"), b"s");
        assert_eq!(
            reason(),
            Some(format!(
                "{database}/CHECKSUMS lists no variation_feature table"
            ))
        );
        std::fs::remove_dir_all(root.path().join(database)).unwrap();
        publish(
            root.path(),
            &format!("{database}/variation_feature.txt"),
            b"v",
        );
        assert_eq!(
            reason(),
            Some(format!("{database}/CHECKSUMS lists no seq_region table"))
        );

        publish(root.path(), &format!("{database}/seq_region.txt.gz"), b"s");
        let work = tempfile::tempdir().unwrap();
        let placed = |to: &str| {
            let work_dir = work.path().join(to);
            let features =
                found()
                    .unwrap()
                    .fetched_when_asked(&site, &work_dir, Some(vec!["21".to_string()]));
            assert_eq!(features.named.as_deref(), Some(&["21".to_string()][..]));
            // Nothing is placed before the build asks.
            assert!(!work_dir.exists());
            let tables = (features.fetch)().unwrap();
            assert_eq!(tables.scratch, work_dir.join("variation_features"));
            let name = |path: &Path| {
                assert!(path.is_file());
                path.strip_prefix(work_dir.join("variation_dumps"))
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .to_string()
            };
            (
                name(&tables.variation_feature),
                name(&tables.seq_region),
                tables.attrib.as_deref().map(name),
            )
        };
        assert_eq!(
            placed("no_attrib"),
            (
                "variation_feature.txt".to_string(),
                "seq_region.txt.gz".to_string(),
                None
            )
        );
        publish(root.path(), &format!("{database}/attrib.txt.gz"), b"a");
        assert_eq!(
            placed("every_table"),
            (
                "variation_feature.txt".to_string(),
                "seq_region.txt.gz".to_string(),
                Some("attrib.txt.gz".to_string())
            )
        );
    }

    /// Reads the listings of the Ensembl sites, and no file they name. The
    /// GRCh38 variation database lists every table of the three cells and
    /// its GVF directory the failed variants' file. The GRCh37 tree lists
    /// that file in its own GVF directory and no `variation` table. A
    /// species of an Ensembl Genomes division lists both in its division's
    /// tree.
    #[test]
    #[ignore = "reads https://ftp.ensembl.org and https://ftp.ensemblgenomes.ebi.ac.uk"]
    fn live_trees_list_the_files_of_the_variant_cells() {
        let found = |species: &str, assembly: &str, release: u32| {
            let located = locate(&Site::defaults(), species, assembly, release).unwrap();
            let options = InstallOptions {
                release,
                ..options(Path::new("unused"), Site::defaults(), species, assembly)
            };
            let inputs = Inputs::find(&located, species, assembly, release).unwrap();
            let sources =
                VariantCellSources::find(&located, &options, inputs.variation.as_ref()).unwrap();
            let gvf = sources
                .failed_gvf
                .as_ref()
                .map(|gvf| format!("{}/{}", gvf.directory.relative, gvf.name))
                .map_err(String::clone);
            let not_found: Vec<String> = sources
                .not_found()
                .into_iter()
                .map(|(what, reason)| format!("{what}: {reason}"))
                .collect();
            (gvf, not_found)
        };

        assert_eq!(
            found("homo_sapiens", "GRCh38", 116),
            (
                Ok(
                    "pub/release-116/variation/gvf/homo_sapiens/homo_sapiens_failed.gvf.gz"
                        .to_string()
                ),
                Vec::new()
            )
        );
        let no_table =
            "pub/grch37/release-116/mysql/homo_sapiens_variation_116_37/CHECKSUMS lists no variation table";
        assert_eq!(
            found("homo_sapiens", "GRCh37", 116),
            (
                Ok(
                    "pub/grch37/release-116/variation/gvf/homo_sapiens/homo_sapiens_failed.gvf.gz"
                        .to_string()
                ),
                vec![
                    format!("PubMed identifier: {no_table}"),
                    format!("variant synonym: {no_table}"),
                    format!("failed_variation table: {no_table}"),
                ]
            )
        );
        assert_eq!(
            found("saccharomyces_cerevisiae", "R64-1-1", 63),
            (
                Ok("pub/release-63/fungi/variation/gvf/saccharomyces_cerevisiae/saccharomyces_cerevisiae_failed.gvf.gz"
                    .to_string()),
                Vec::new()
            )
        );
    }

    /// Reads the Ensembl site's listings, and no file they name: for
    /// `homo_sapiens` GRCh37 the regulatory annotation directory of the
    /// release tree lists the three files, and the funcgen database of the
    /// GRCh37 tree lists every table an installation reads from it.
    #[test]
    #[ignore = "reads https://ftp.ensembl.org"]
    fn live_grch37_lists_its_regulatory_files_in_two_trees() {
        let (species, assembly) = ("homo_sapiens", "GRCh37");
        let located = locate(&Site::defaults(), species, assembly, 116).unwrap();
        let options = options(Path::new("unused"), Site::defaults(), species, assembly);
        let files = match RegulatoryFiles::find(&located, &options).unwrap() {
            Part::Published(files) => files,
            Part::Unpublished(reason) => panic!("{reason}"),
        };
        assert_eq!(
            files.directory.relative,
            "pub/release-116/regulation/homo_sapiens/GRCh37/annotation"
        );
        assert!(files.features.is_some() && files.motifs.is_some());
        let cells = RegulatoryCells::find(&located, &options, &files).unwrap();
        assert_eq!(
            cells.activity.as_deref(),
            Ok("Homo_sapiens.GRCh37.regulatory_activity.v116.tsv.gz")
        );
        let funcgen = cells.funcgen.as_ref().expect("a funcgen database");
        assert_eq!(
            funcgen.directory.relative,
            "pub/grch37/release-116/mysql/homo_sapiens_funcgen_116_37"
        );
        assert_eq!(cells.not_found(), []);
    }

    /// Reads the Ensembl Genomes site: a plant is found there with no
    /// `--site`, under the Ensembl Genomes release, with the core database the
    /// species list names, and every input of its installation is in its
    /// directory's listing.
    #[test]
    #[ignore = "reads https://ftp.ensemblgenomes.ebi.ac.uk"]
    fn live_ensembl_genomes_division_lists_every_input() {
        let species = "arabidopsis_thaliana";
        let located = locate(&Site::defaults(), species, "TAIR10", 63).unwrap();
        assert!(matches!(&located.site, Site::Url(url) if url == download::ENSEMBL_GENOMES_SITE));
        assert_eq!(
            located.layout,
            Layout::Division {
                division: "plants".to_string()
            }
        );
        assert_eq!(
            located.core_db.as_deref(),
            Some("arabidopsis_thaliana_core_63_116_11")
        );
        let inputs = Inputs::find(&located, species, "TAIR10", 63).unwrap();
        assert_eq!(inputs.gff3.name, "Arabidopsis_thaliana.TAIR10.63.gff3.gz");
        assert_eq!(inputs.gtf.name, "Arabidopsis_thaliana.TAIR10.63.gtf.gz");
        assert_eq!(
            inputs.protein_fasta.name,
            "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
        );
        assert_eq!(
            inputs.genome.file.name,
            "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
        );
        assert!(inputs.genome.indexed);
        assert!(inputs.variation.is_some());
    }

    /// Reads the Ensembl Genomes site: a bacterial genome is placed in its
    /// collection and refused for the database it shares.
    #[test]
    #[ignore = "reads https://ftp.ensemblgenomes.ebi.ac.uk"]
    fn live_collection_genome_is_refused() {
        let species = "escherichia_coli_str_k_12_substr_mg1655_gca_000005845";
        let located = locate(&Site::defaults(), species, "ASM584v2", 63).unwrap();
        assert_eq!(
            located.layout,
            Layout::Collection {
                division: "bacteria".to_string(),
                collection: "bacteria_0_collection".to_string()
            }
        );
        let error = located
            .buildable(species, "ASM584v2", 63)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("bacteria_0_collection_core_63_116_1 holds every genome"),
            "{error}"
        );
    }

    /// Reads the Ensembl site: a human pangenome assembly, in no species
    /// list, is found in the rapid release by its accession and refused for
    /// the table dumps that tree lacks.
    #[test]
    #[ignore = "reads https://ftp.ensembl.org"]
    fn live_rapid_release_assembly_is_refused() {
        let species = "homo_sapiens_gca009914755v4";
        let located = locate(&Site::defaults(), species, "T2T-CHM13v2.0", 107).unwrap();
        assert_eq!(
            located.layout,
            Layout::RapidRelease {
                organism: "Homo_sapiens".to_string(),
                accession: "GCA_009914755.4".to_string(),
                gene_set: "2022_07".to_string()
            }
        );
        let error = located
            .buildable(species, "T2T-CHM13v2.0", 107)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains(
                "pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/ publishes no core table dumps"
            ),
            "{error}"
        );
    }
}
