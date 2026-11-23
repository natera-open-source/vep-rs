// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Download Ensembl GFF3, protein FASTA, variation VCF and MySQL table dump
//! files, or resolve local paths. A GFF3, a protein FASTA or a variation VCF
//! that a build downloads lands in the build's [`WorkDir`], never in the
//! system's temp directory; a table dump lands in the directory
//! `--core-dumps` names.
//!
//! [`Layout`] is the tree a genome's files are published in: the directory
//! and the file name of each input under the root of its site. A release's
//! species list ([`species_entry`]) names the division and the core database
//! of every genome of the release, which places a species in its layout.

use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use tracing::{info, warn};

use crate::core_dumps::{table_path, CORE_TABLES, VARIATION_TABLES};
use crate::count::counted;
use crate::install::{Published, Site};
use crate::variation::{self, GenomeWideInputs, GermlineVcf};

/// The GFF3 to read: the local file when one is given, else the release's
/// Ensembl GFF3, fetched through `from`, which looks the genome up only when
/// the file is not already in the work directory.
pub fn resolve_gff3(
    local: Option<&str>,
    release: u32,
    assembly: &str,
    species: &str,
    from: &Downloads<'_>,
) -> Result<PathBuf> {
    match local {
        Some(p) => {
            let path = PathBuf::from(p);
            if !path.exists() {
                bail!("Local GFF3 file not found: {}", path.display());
            }
            Ok(path)
        }
        None => download_gff3(from, release, assembly, species),
    }
}

/// The protein FASTA to read: the local file when one is given, else the
/// release's Ensembl protein FASTA, fetched through `from`, which looks the
/// genome up only when the file is not already in the work directory.
pub fn resolve_protein_fasta(
    local: Option<&str>,
    release: u32,
    assembly: &str,
    species: &str,
    from: &Downloads<'_>,
) -> Result<PathBuf> {
    match local {
        Some(p) => {
            let path = PathBuf::from(p);
            if !path.exists() {
                bail!("Local protein FASTA not found: {}", path.display());
            }
            Ok(path)
        }
        None => download_protein_fasta(from, release, assembly, species),
    }
}

/// The Ensembl site: the vertebrate release and, in a tree of its own, GRCh37.
pub const ENSEMBL_SITE: &str = "https://ftp.ensembl.org";

/// The Ensembl Genomes site: the plants, metazoa, fungi, protists and bacteria
/// divisions.
pub const ENSEMBL_GENOMES_SITE: &str = "https://ftp.ensemblgenomes.ebi.ac.uk";

/// The species lists a release directory publishes, each named in that
/// directory's `CHECKSUMS`: the Ensembl site's list of its vertebrates and the
/// Ensembl Genomes site's list over all its divisions.
pub const SPECIES_LISTS: [&str; 2] = ["species_EnsemblVertebrates.txt", "species.txt"];

/// The tree a genome's files are published in, each path relative to the root
/// of its site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Layout {
    /// `pub/release-<N>/`, the Ensembl release.
    Release,
    /// `pub/grch37/release-<N>/`, the tree Ensembl keeps GRCh37 in; its gene
    /// set is the release-87 one.
    Grch37,
    /// `pub/release-<N>/<division>/`, a species of an Ensembl Genomes division
    /// with a core database of its own; `N` is the Ensembl Genomes release,
    /// the number its file names carry.
    Division { division: String },
    /// A genome of a collection: the division's tree with the collection's
    /// directory above the species', and one core database for every genome
    /// of the collection.
    Collection {
        division: String,
        collection: String,
    },
    /// `pub/rapid-release/species/<organism>/<accession>/ensembl/`, an
    /// assembly of the rapid release: a dated gene-set directory and a genome
    /// directory, each listing its files in `md5sum.txt`, and no table dump.
    RapidRelease {
        organism: String,
        accession: String,
        gene_set: String,
    },
}

impl Layout {
    /// The tree of `assembly` on the Ensembl site: GRCh37's own, else the
    /// release.
    pub fn ensembl(assembly: &str) -> Self {
        if assembly == "GRCh37" {
            Self::Grch37
        } else {
            Self::Release
        }
    }

    /// The directory the tree hangs from.
    pub fn root(&self, release: u32) -> String {
        match self {
            Self::Release => format!("pub/release-{release}"),
            Self::Grch37 => format!("pub/grch37/release-{release}"),
            Self::Division { division } | Self::Collection { division, .. } => {
                format!("pub/release-{release}/{division}")
            }
            Self::RapidRelease {
                organism,
                accession,
                ..
            } => format!("pub/rapid-release/species/{organism}/{accession}/ensembl"),
        }
    }

    /// The directory of one kind of file (`gff3`, `gtf`, `fasta`) of a
    /// species, under the collection's directory for a genome of one.
    fn species_dir(&self, release: u32, kind: &str, species: &str) -> String {
        let root = self.root(release);
        match self {
            Self::Collection { collection, .. } => format!("{root}/{kind}/{collection}/{species}"),
            _ => format!("{root}/{kind}/{species}"),
        }
    }

    /// The gene-set version in the GFF3 and GTF file names: the release, the
    /// rapid release's dated directory, and 87 at most on GRCh37. Ensembl
    /// froze the GRCh37 gene set at release 87: every `grch37/release-N/`
    /// directory from 88 onward carries `<Species>.GRCh37.87.gff3.gz`, not a
    /// file named for `N`.
    pub fn gene_set(&self, release: u32) -> String {
        match self {
            Self::Grch37 => release.min(87).to_string(),
            Self::RapidRelease { gene_set, .. } => gene_set.clone(),
            _ => release.to_string(),
        }
    }

    /// A file of the gene set, by its directory and its name: `kind` is the
    /// directory and extension of a release tree (`gff3`, `gtf`), `dated` the
    /// rapid release's name for the same file.
    fn gene_set_file(
        &self,
        release: u32,
        assembly: &str,
        species: &str,
        kind: &str,
        dated: &str,
    ) -> (String, String) {
        match self {
            Self::RapidRelease {
                organism,
                accession,
                gene_set,
            } => (
                format!("{}/geneset/{gene_set}", self.root(release)),
                format!("{organism}-{accession}-{gene_set}-{dated}"),
            ),
            _ => (
                self.species_dir(release, kind, species),
                format!(
                    "{}.{assembly}.{}.{kind}.gz",
                    capitalize_species(species),
                    self.gene_set(release)
                ),
            ),
        }
    }

    /// The GFF3: its directory and its file name.
    pub fn gff3_location(&self, release: u32, assembly: &str, species: &str) -> (String, String) {
        self.gene_set_file(release, assembly, species, "gff3", "genes.gff3.gz")
    }

    /// The GTF: its directory and its file name.
    pub fn gtf_location(&self, release: u32, assembly: &str, species: &str) -> (String, String) {
        self.gene_set_file(release, assembly, species, "gtf", "genes.gtf.gz")
    }

    /// The protein FASTA: its directory and its file name.
    pub fn protein_fasta_location(
        &self,
        release: u32,
        assembly: &str,
        species: &str,
    ) -> (String, String) {
        match self {
            Self::RapidRelease { .. } => {
                self.gene_set_file(release, assembly, species, "pep", "pep.fa.gz")
            }
            _ => (
                format!("{}/pep", self.species_dir(release, "fasta", species)),
                format!("{}.{assembly}.pep.all.fa.gz", capitalize_species(species)),
            ),
        }
    }

    /// The genome FASTA directories, in the order the installer tries them:
    /// `dna_index/`, whose files come with their `.fai` and `.gzi`, then
    /// `dna/`; the rapid release has one, `genome/`.
    pub fn genome_fasta_dirs(&self, release: u32, species: &str) -> Vec<String> {
        match self {
            Self::RapidRelease { .. } => vec![format!("{}/genome", self.root(release))],
            _ => {
                let fasta = self.species_dir(release, "fasta", species);
                vec![format!("{fasta}/dna_index"), format!("{fasta}/dna")]
            }
        }
    }

    /// The directory of the MySQL table dumps of `database`; `None` in the
    /// rapid release, which publishes none.
    pub fn mysql_dir(&self, release: u32, database: &str) -> Option<String> {
        match self {
            Self::RapidRelease { .. } => None,
            _ => Some(format!("{}/mysql/{database}", self.root(release))),
        }
    }

    /// The directory of a species' variation VCFs. `None` for a genome of a
    /// collection, whose division publishes none, and in the rapid release,
    /// which keeps them under dated directories this layout does not name.
    pub fn variation_vcf_dir(&self, release: u32, species: &str) -> Option<String> {
        match self {
            Self::Collection { .. } | Self::RapidRelease { .. } => None,
            _ => Some(format!("{}/variation/vcf/{species}", self.root(release))),
        }
    }

    /// The directory of a species' variation GVF files, beside that of its
    /// variation VCFs and `None` where that one is.
    pub fn variation_gvf_dir(&self, release: u32, species: &str) -> Option<String> {
        match self {
            Self::Collection { .. } | Self::RapidRelease { .. } => None,
            _ => Some(format!("{}/variation/gvf/{species}", self.root(release))),
        }
    }

    /// The directory of a species' regulatory annotation, which holds one
    /// directory an assembly, each with the GFF files under `annotation/`.
    /// GRCh37's stands in the release tree beside GRCh38's: the GRCh37 tree
    /// keeps a copy whose `CHECKSUMS` holds MD5 lines, not the BSD sums the
    /// installer checks a file against. `None` for a genome of a collection
    /// and in the rapid release.
    pub fn regulation_dir(&self, release: u32, species: &str) -> Option<String> {
        match self {
            Self::Release | Self::Grch37 => {
                Some(format!("pub/release-{release}/regulation/{species}"))
            }
            Self::Division { .. } => Some(format!("{}/regulation/{species}", self.root(release))),
            Self::Collection { .. } | Self::RapidRelease { .. } => None,
        }
    }
}

/// The release directory of the Ensembl site, relative to its root:
/// `pub/release-N`, or `pub/grch37/release-N` for GRCh37, which Ensembl
/// publishes in a tree of its own.
pub fn release_dir(release: u32, assembly: &str) -> String {
    Layout::ensembl(assembly).root(release)
}

/// A genome's row in a release's species list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpeciesEntry {
    /// `EnsemblVertebrates`, `EnsemblPlants`, `EnsemblBacteria` and so on.
    pub division: String,
    /// The name of the genome's core database, its directory under `mysql/`.
    pub core_db: String,
}

impl SpeciesEntry {
    /// The layout the row places `species` in: the release tree for a
    /// vertebrate, else the tree of its division, the directory named for the
    /// division without its `Ensembl` prefix, in lower case. A core database
    /// not named for the species is a collection's, named
    /// `<collection>_core_<...>`.
    pub fn layout(&self, species: &str) -> Layout {
        let division = self
            .division
            .strip_prefix("Ensembl")
            .unwrap_or(&self.division)
            .to_lowercase();
        if division == "vertebrates" {
            return Layout::Release;
        }
        if self.core_db.starts_with(&format!("{species}_core_")) {
            return Layout::Division { division };
        }
        let collection = self.core_db.split("_core_").next().unwrap_or_default();
        Layout::Collection {
            division,
            collection: collection.to_string(),
        }
    }
}

/// The row of `species` in the text of a species list, a tab-separated table
/// whose first line names its columns behind a `#`; the columns are read by
/// name. `None` when no row is the species'; an error when the header lacks a
/// column read here or the row stops before it.
pub fn species_entry(text: &str, species: &str) -> Result<Option<SpeciesEntry>> {
    let mut lines = text.lines();
    let header: Vec<&str> = lines
        .next()
        .unwrap_or_default()
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|title| *title == name)
            .with_context(|| format!("the species list has no {name} column"))
    };
    let (species_at, division_at, core_db_at) =
        (column("species")?, column("division")?, column("core_db")?);
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.get(species_at) != Some(&species) {
            continue;
        }
        let field = |at: usize| {
            fields
                .get(at)
                .map(|value| value.to_string())
                .with_context(|| format!("the species list's row of {species} is cut short"))
        };
        return Ok(Some(SpeciesEntry {
            division: field(division_at)?,
            core_db: field(core_db_at)?,
        }));
    }
    Ok(None)
}

/// The variation database beside the core database `core_db` of `species`:
/// the same name with `variation` in place of `core`. `None` when the core
/// database is not named for the species.
pub fn variation_database(core_db: &str, species: &str) -> Option<String> {
    database_beside(core_db, species, "variation")
}

/// The funcgen database beside the core database `core_db` of `species`: the
/// same name with `funcgen` in place of `core`. `None` when the core
/// database is not named for the species.
pub fn funcgen_database(core_db: &str, species: &str) -> Option<String> {
    database_beside(core_db, species, "funcgen")
}

/// The database of `kind` named as the core database `core_db` of `species`
/// is, `kind` standing in place of `core`.
fn database_beside(core_db: &str, species: &str, kind: &str) -> Option<String> {
    let rest = core_db.strip_prefix(species)?.strip_prefix("_core_")?;
    Some(format!("{species}_{kind}_{rest}"))
}

/// The core databases of `species` among the entries of a release tree's
/// `mysql/` directory: `<species>_core_<release>_<n>`, `n` a number the
/// assembly name does not give (`rattus_norvegicus_core_116_1` is GRCr8).
pub fn core_databases<'a>(entries: &'a [String], species: &str, release: u32) -> Vec<&'a str> {
    let prefix = format!("{species}_core_{release}_");
    entries
        .iter()
        .map(String::as_str)
        .filter(|entry| {
            entry
                .strip_prefix(&prefix)
                .is_some_and(|n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()))
        })
        .collect()
}

/// The core database the release trees name for an assembly whose name ends
/// in a number, `<species>_core_<release>_<digits>` (`GRCh38` is 38, `GRCh37`
/// 37); `None` for a name with no trailing digits. Most species' databases
/// end in another number, so the name is tried only on a server that lists
/// neither its species nor its directories.
pub fn numbered_core_database(species: &str, assembly: &str, release: u32) -> Option<String> {
    let digits = &assembly[assembly
        .trim_end_matches(|c: char| c.is_ascii_digit())
        .len()..];
    (!digits.is_empty()).then(|| format!("{species}_core_{release}_{digits}"))
}

/// The rapid-release directory of a species named for its assembly's
/// accession, `<organism>_gca<digits>v<version>`: the organism capitalised and
/// the accession as `GCA_<digits>.<version>`. `None` for any other name.
pub fn rapid_release_name(species: &str) -> Option<(String, String)> {
    let (organism, tail) = species.rsplit_once("_gca")?;
    let (digits, version) = tail.split_once('v')?;
    let is_number = |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
    if organism.is_empty() || !is_number(digits) || !is_number(version) {
        return None;
    }
    Some((
        capitalize_species(organism),
        format!("GCA_{digits}.{version}"),
    ))
}

/// The entries of a web server's directory index page: the target of each
/// link to a child of the directory, without the slash that ends a
/// directory's. The links to the parent, to a sort order and to another site
/// are passed over.
pub fn index_entries(html: &str) -> Vec<String> {
    html.split("href=\"")
        .skip(1)
        .filter_map(|rest| rest.split('"').next())
        .filter(|target| !target.starts_with(['?', '/', '.']) && !target.contains("://"))
        .map(|target| target.trim_end_matches('/').to_string())
        .filter(|entry| !entry.is_empty())
        .collect()
}

/// A genome that no site publishes, or that its site publishes without an
/// input the build reads: what a lookup answers when the genome is not there,
/// kept as a typed error so a caller can tell it from a failed read.
#[derive(Debug)]
pub struct NotPublished(pub String);

impl std::fmt::Display for NotPublished {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for NotPublished {}

/// Whether an error is a [`NotPublished`] answer.
pub fn is_not_published(error: &anyhow::Error) -> bool {
    error.downcast_ref::<NotPublished>().is_some()
}

/// The URLs of a genome's core and variation MySQL directories, each ending
/// in a slash; no variation URL where the core database's name gives none.
pub struct DatabaseUrls {
    pub core: String,
    pub variation: Option<String>,
}

/// Make `dir` hold every table dump the cache build reads: a core table
/// missing from it is downloaded (an error when that fails), and a missing
/// variation table is downloaded as `.txt.gz`, then as `.txt` when the first
/// is a 404, and left absent with a warning when both are 404, so a species
/// without a variation database still builds. Any other failure stops the
/// build. `databases` is asked for the genome's MySQL directories only when a
/// table is missing; when the core tables are all there and it answers that
/// no site publishes a database of the genome's own ([`NotPublished`]), the
/// variation tables are left absent with the same warning.
pub fn resolve_core_dumps(dir: &Path, databases: &dyn Fn() -> Result<DatabaseUrls>) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("Failed to create {}", dir.display()))?;
    let missing = |tables: &[&str]| tables.iter().any(|table| table_path(dir, table).is_none());
    let core_missing = missing(&CORE_TABLES);
    if !core_missing && !missing(&VARIATION_TABLES) {
        return Ok(());
    }
    let urls = match databases() {
        Ok(urls) => urls,
        Err(e) if !core_missing && is_not_published(&e) => {
            warn!("No variation database is published for this genome; gene_phenotype is 0 on every transcript");
            return Ok(());
        }
        Err(e) => return Err(e),
    };
    for table in CORE_TABLES {
        if table_path(dir, table).is_some() {
            continue;
        }
        let url = format!("{}{table}.txt.gz", urls.core);
        download_file(&url, &dir.join(format!("{table}.txt.gz")))
            .with_context(|| format!("Failed to download core table {table}"))?;
    }
    let Some(variation) = &urls.variation else {
        warn!("No variation database is named for this genome; gene_phenotype is 0 on every transcript");
        return Ok(());
    };
    for table in VARIATION_TABLES {
        if table_path(dir, table).is_some() {
            continue;
        }
        let gz_url = format!("{variation}{table}.txt.gz");
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
pub fn capitalize_species(species: &str) -> String {
    let mut chars = species.chars();
    match chars.next() {
        None => String::new(),
        Some(first) => {
            let upper: String = first.to_uppercase().collect();
            format!("{}{}", upper, chars.as_str())
        }
    }
}

/// The work directory of one build: every file the build downloads lands
/// under it, except a table dump of a build named one input at a time, which
/// lands in the directory `--core-dumps` names. An `--auto` installation has
/// `<cache root>/tmp/<species>_<release>_<assembly>` ([`Self::of_release`])
/// and a build named one input at a time
/// `<output directory>/tmp/vep_cache_builder/<species>_<release>_<assembly>`
/// ([`Self::of_output_dir`]). The directory is removed once the cache is built
/// and kept when the build fails with a file in it, so the next run reuses
/// what it holds.
#[derive(Debug, Clone)]
pub struct WorkDir {
    path: PathBuf,
    /// What [`Self::remove`] removes: `path`, or the directory above it that
    /// holds the work directory of every build into one output directory.
    cleared: PathBuf,
}

impl WorkDir {
    pub fn of_release(root: &Path, species: &str, release: u32, assembly: &str) -> Self {
        let path = root
            .join("tmp")
            .join(format!("{species}_{release}_{assembly}"));
        Self {
            cleared: path.clone(),
            path,
        }
    }

    /// The work directory of a build into `output_dir`:
    /// `tmp/vep_cache_builder/<species>_<release>_<assembly>` under it.
    /// [`Self::remove`] removes `tmp/vep_cache_builder` whole. The checksum
    /// manifest of a built cache lists every file under its directory, so a
    /// build that succeeds has to remove the downloads a failed build of
    /// another release, assembly or species kept there as well as its own.
    pub fn of_output_dir(output_dir: &Path, species: &str, release: u32, assembly: &str) -> Self {
        let cleared = output_dir.join("tmp").join("vep_cache_builder");
        Self {
            path: cleared.join(format!("{species}_{release}_{assembly}")),
            cleared,
        }
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Removes the directory with everything under it, then `tmp` itself when
    /// that leaves it empty; for a build into an output directory, the
    /// directory removed is the `tmp/vep_cache_builder` that holds the work
    /// directory of every build into it. A directory that does not exist is
    /// left alone, and so is its `tmp`.
    pub fn remove(&self) -> Result<()> {
        if !self.cleared.exists() {
            return Ok(());
        }
        std::fs::remove_dir_all(&self.cleared)
            .with_context(|| format!("Failed to remove {}", self.cleared.display()))?;
        if let Some(tmp) = self.cleared.parent() {
            let _ = std::fs::remove_dir(tmp);
        }
        Ok(())
    }

    /// The error of a failed build: `error` followed by the work directory
    /// the build leaves behind when it holds a file the next run reuses. A
    /// work directory that holds no file is removed and `error` returned
    /// alone, as it is when the build created none.
    pub fn kept(&self, error: anyhow::Error) -> anyhow::Error {
        if !self.path.exists() {
            return error;
        }
        if !holds_a_file(&self.path) {
            self.remove_own();
            return error;
        }
        anyhow::anyhow!(
            "{error:#}; the work directory {} is kept for the next run",
            self.path.display()
        )
    }

    /// Removes this build's own directory, then each directory above it, up
    /// to `tmp`, that this leaves empty. The work directory another build
    /// into the same output directory kept stays, where [`Self::remove`]
    /// would remove it.
    fn remove_own(&self) {
        let _ = std::fs::remove_dir_all(&self.path);
        if self.cleared != self.path {
            let _ = std::fs::remove_dir(&self.cleared);
        }
        if let Some(tmp) = self.cleared.parent() {
            let _ = std::fs::remove_dir(tmp);
        }
    }

    /// Places the file at `relative` on `site` at `dest`, a path of this
    /// directory. The directory is created once the site is known to hold
    /// that file, here for a local copy of the site's tree and by
    /// [`download_file`] for a URL. A file a local copy lacks is an error
    /// [`is_not_found`] reads as a 404.
    fn fetch(&self, site: &Site, relative: &str, dest: &Path) -> Result<()> {
        if let Site::Dir(root) = site {
            let source = root.join(relative);
            if !source.is_file() {
                return Err(anyhow::Error::new(Absent(source)));
            }
            std::fs::create_dir_all(&self.path)
                .with_context(|| format!("Failed to create {}", self.path.display()))?;
        }
        site.place(relative, dest)
    }

    /// As [`Self::fetch`], keeping a `dest` that already exists.
    fn fetch_cached(&self, site: &Site, relative: &str, dest: &Path) -> Result<()> {
        if dest.exists() {
            info!("Using cached download: {}", dest.display());
            return Ok(());
        }
        self.fetch(site, relative, dest)
    }

    /// Fetches the `.csi` index of the file at `relative` beside `dest`.
    fn fetch_index(&self, site: &Site, relative: &str, dest: &Path) -> Result<()> {
        let mut index_dest = dest.as_os_str().to_owned();
        index_dest.push(".csi");
        self.fetch_cached(site, &format!("{relative}.csi"), Path::new(&index_dest))
    }
}

/// Whether `dir` holds a file at any depth. A directory that cannot be read
/// counts as holding one, so it is never taken for empty.
fn holds_a_file(dir: &Path) -> bool {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return true;
    };
    entries.into_iter().any(|entry| match entry {
        Ok(entry) => {
            let path = entry.path();
            !path.is_dir() || holds_a_file(&path)
        }
        Err(_) => true,
    })
}

/// Where a build named one input at a time fetches a file it was not given,
/// and the work directory the file lands in.
#[derive(Clone, Copy)]
pub struct Downloads<'a> {
    /// The genome on the site that publishes it, which names its GFF3 and
    /// its protein FASTA.
    pub published: &'a Published<'a>,
    /// The site of the variation VCFs, whose release tree
    /// [`variation_vcf_dir`] names.
    pub site: &'a Site,
    pub work: &'a WorkDir,
}

/// A file a local copy of the site's tree does not hold.
#[derive(Debug)]
struct Absent(PathBuf);

impl std::fmt::Display for Absent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} is not on the site", self.0.display())
    }
}

impl std::error::Error for Absent {}

/// A downloaded file's path in the work directory: a name of the species,
/// assembly and release it belongs to, then `suffix`.
fn download_dest(
    work: &WorkDir,
    species: &str,
    assembly: &str,
    release: u32,
    suffix: &str,
) -> PathBuf {
    work.path().join(format!(
        "vep_cache_builder_{species}_{assembly}_{release}{suffix}"
    ))
}

/// Where the downloaded GFF3 lands.
fn gff3_dest(work: &WorkDir, species: &str, assembly: &str, release: u32) -> PathBuf {
    download_dest(work, species, assembly, release, ".gff3.gz")
}

/// Where the downloaded protein FASTA lands.
fn protein_fasta_dest(work: &WorkDir, species: &str, assembly: &str, release: u32) -> PathBuf {
    download_dest(work, species, assembly, release, ".pep.all.fa.gz")
}

/// Where the downloaded variation VCF of `chromosome` lands.
fn variation_vcf_dest(
    work: &WorkDir,
    species: &str,
    assembly: &str,
    release: u32,
    chromosome: &str,
) -> PathBuf {
    download_dest(
        work,
        species,
        assembly,
        release,
        &format!("_chr{chromosome}.vcf.gz"),
    )
}

/// Where the downloaded genome-wide variation VCF `name` lands.
fn genome_wide_dest(
    work: &WorkDir,
    species: &str,
    assembly: &str,
    release: u32,
    name: &str,
) -> PathBuf {
    download_dest(work, species, assembly, release, &format!("_{name}"))
}

/// Fetches the genome's GFF3 into the work directory unless the directory
/// holds it.
fn download_gff3(
    from: &Downloads<'_>,
    release: u32,
    assembly: &str,
    species: &str,
) -> Result<PathBuf> {
    let dest = gff3_dest(from.work, species, assembly, release);

    if dest.exists() {
        info!("Using cached GFF3: {}", dest.display());
        return Ok(dest);
    }

    let (site, relative) = from.published.gff3()?;
    from.work.fetch(site, &relative, &dest)?;
    Ok(dest)
}

/// Fetches the genome's protein FASTA into the work directory unless the
/// directory holds it.
fn download_protein_fasta(
    from: &Downloads<'_>,
    release: u32,
    assembly: &str,
    species: &str,
) -> Result<PathBuf> {
    let dest = protein_fasta_dest(from.work, species, assembly, release);

    if dest.exists() {
        info!("Using cached protein FASTA: {}", dest.display());
        return Ok(dest);
    }

    let (site, relative) = from.published.protein_fasta()?;
    from.work.fetch(site, &relative, &dest)?;
    Ok(dest)
}

/// The species Ensembl publishes variation VCFs for in the human layout: one
/// germline file a chromosome ([`variation_vcf_name`]) and the genome-wide
/// somatic, phenotype-associated, clinically-associated and 1000 Genomes
/// files. Every other species with variation data has one genome-wide germline
/// file ([`genome_vcf_name`]) and none of the others.
pub const HUMAN_SPECIES: &str = "homo_sapiens";

/// The chromosomes Ensembl publishes a per-chromosome variation VCF for under
/// `homo_sapiens`.
const HUMAN_CHROMOSOMES: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

/// The chromosomes a build without `--chromosomes` reads a per-chromosome
/// variation VCF for: [`HUMAN_CHROMOSOMES`] under `homo_sapiens`, and none
/// under any other species, whose germline records come from its genome-wide
/// file for every sequence that file holds records on. The transcript model
/// has no say: its sequence names say which sequences the model annotates,
/// not which files Ensembl publishes.
pub fn default_chromosomes(species: &str) -> Vec<String> {
    if species == HUMAN_SPECIES {
        HUMAN_CHROMOSOMES.iter().map(|c| c.to_string()).collect()
    } else {
        Vec::new()
    }
}

/// The germline variation VCFs a download of `species` fetches from the
/// directory of its variation VCFs: under `homo_sapiens` the file of each of
/// `chromosomes`, else of [`default_chromosomes`]; under any other species
/// its one genome-wide file, whose read `chromosomes` narrows.
pub fn germline_vcf_names(species: &str, chromosomes: Option<&[String]>) -> Vec<String> {
    if species != HUMAN_SPECIES {
        return vec![genome_vcf_name(species)];
    }
    let chromosomes = match chromosomes {
        Some(list) => list.to_vec(),
        None => default_chromosomes(species),
    };
    chromosomes
        .iter()
        .map(|chromosome| variation_vcf_name(species, chromosome))
        .collect()
}

/// The variation VCF inputs of a build: the germline VCF each chromosome's
/// records are read from, and the genome-wide files read beside them.
pub struct VariationInputs {
    /// `(chromosome, germline VCF)` of each chromosome built, by chromosome
    /// name: the chromosome's own file, or the species' one genome-wide file
    /// under every chromosome.
    pub per_chromosome: Vec<(String, GermlineVcf)>,
    pub genome_wide: GenomeWideInputs,
}

/// The genome-wide variation VCFs Ensembl publishes beside the per-chromosome
/// files of `homo_sapiens`, by the name each has in the FTP directory and in a
/// local directory. The 1000 Genomes file exists for GRCh38 alone.
pub const THOUSAND_GENOMES_FILE: &str = "1000GENOMES-phase_3.vcf.gz";
const SOMATIC_SUFFIX: &str = "_somatic.vcf.gz";
const GERMLINE_SUFFIXES: [&str; 2] = [
    "_phenotype_associated.vcf.gz",
    "_clinically_associated.vcf.gz",
];

/// The name of the GVF file of the failed variants of `species` in the
/// directory of its variation GVF files: every feature of a variant the
/// variation database marks as failed.
pub fn failed_gvf_name(species: &str) -> String {
    format!("{species}_failed.gvf.gz")
}

/// The names of the genome-wide variation VCFs of `species`, the 1000 Genomes
/// file first; for a species other than `homo_sapiens` its one germline file
/// ([`genome_vcf_name`]) is the last.
pub fn genome_wide_names(species: &str) -> Vec<String> {
    let mut names = vec![THOUSAND_GENOMES_FILE.to_string()];
    names.push(format!("{species}{SOMATIC_SUFFIX}"));
    names.extend(GERMLINE_SUFFIXES.iter().map(|s| format!("{species}{s}")));
    if species != HUMAN_SPECIES {
        names.push(genome_vcf_name(species));
    }
    names
}

/// Resolve the variation VCF inputs: scan a local directory
/// ([`scan_local_vcfs`]) or fetch the files through `from`. A genome-wide
/// file the directory or the site lacks is left out, a warning under
/// `homo_sapiens` naming what the cache will lack; a downloaded file's index
/// (`.csi`) is fetched with it, so each chromosome's records are read through
/// the index. The download is the per-chromosome files of `chromosomes`, else
/// of [`default_chromosomes`], under `homo_sapiens`, and the one genome-wide
/// germline file under any other species, read for `chromosomes`, else for
/// every sequence it holds records on.
pub fn resolve_variation_inputs(
    local_vcf_dir: &Option<String>,
    release: u32,
    assembly: &str,
    species: &str,
    chromosomes: Option<&[String]>,
    from: &Downloads<'_>,
) -> Result<VariationInputs> {
    match local_vcf_dir {
        Some(dir) => local_variation_inputs(Path::new(dir), species, chromosomes),
        None => {
            let per_chromosome = if species == HUMAN_SPECIES {
                let chroms: Vec<String> = match chromosomes {
                    Some(c) => c.to_vec(),
                    None => default_chromosomes(species),
                };
                let chroms: Vec<&str> = chroms.iter().map(String::as_str).collect();
                download_variation_vcfs(from, release, assembly, species, &chroms)?
                    .into_iter()
                    .map(|(chrom, path)| (chrom, GermlineVcf::Chromosome(path)))
                    .collect()
            } else {
                genome_vcf_inputs(
                    download_genome_vcf(from, release, assembly, species)?,
                    chromosomes,
                )?
            };
            Ok(VariationInputs {
                per_chromosome,
                genome_wide: download_genome_wide(from, release, assembly, species)?,
            })
        }
    }
}

/// The variation VCF inputs a local directory holds: its germline files
/// ([`scan_local_vcfs`]) and the genome-wide files present under their
/// Ensembl names.
pub fn local_variation_inputs(
    dir: &Path,
    species: &str,
    chromosomes: Option<&[String]>,
) -> Result<VariationInputs> {
    if !dir.is_dir() {
        bail!(
            "Local variation VCF directory not found: {}: {}",
            dir.display(),
            looked_for(species)
        );
    }
    Ok(VariationInputs {
        per_chromosome: scan_local_vcfs(dir, species, chromosomes)?,
        genome_wide: scan_local_genome_wide(dir, species)?,
    })
}

/// The genome-wide germline file `vcf` under each chromosome it is read for:
/// `chromosomes`, else every sequence the file holds records on, by name.
fn genome_vcf_inputs(
    vcf: PathBuf,
    chromosomes: Option<&[String]>,
) -> Result<Vec<(String, GermlineVcf)>> {
    let mut names = match chromosomes {
        Some(list) => list.to_vec(),
        None => variation::sequence_names(&vcf)?,
    };
    if names.is_empty() {
        bail!("{} holds no variation record", vcf.display());
    }
    names.sort();
    Ok(names
        .into_iter()
        .map(|chrom| (chrom, GermlineVcf::Genome(vcf.clone())))
        .collect())
}

/// The two forms of germline file a directory is searched for, as the refusal
/// of a directory that is missing or holds neither says them.
fn looked_for(species: &str) -> String {
    format!(
        "looked for *-chr*.vcf.gz, one file a chromosome, and {}, the one genome-wide file \
         Ensembl publishes for a species other than homo_sapiens",
        genome_vcf_name(species)
    )
}

/// The germline VCFs of a local directory: its `*-chr*.vcf.gz` files, one a
/// chromosome, or, where it holds none, the species' genome-wide file
/// ([`genome_vcf_name`]) read for `chromosomes`, else for every sequence it
/// holds records on. A directory with neither is refused naming both.
fn scan_local_vcfs(
    dir: &Path,
    species: &str,
    chromosomes: Option<&[String]>,
) -> Result<Vec<(String, GermlineVcf)>> {
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
                results.push((chrom, GermlineVcf::Chromosome(path)));
            }
        }
    }
    if results.is_empty() {
        let genome_vcf = dir.join(genome_vcf_name(species));
        if !genome_vcf.is_file() {
            bail!(
                "No variation VCF of {species} in {}: {}",
                dir.display(),
                looked_for(species)
            );
        }
        info!(
            "Found local genome-wide variation VCF: {}",
            genome_vcf.display()
        );
        return genome_vcf_inputs(genome_vcf, chromosomes);
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
/// files. A file `locate` does not find is left out, saying what the cache
/// lacks without it: a warning under `homo_sapiens`, the one species Ensembl
/// publishes the files for, and an information line under any other.
fn genome_wide_inputs(
    species: &str,
    mut locate: impl FnMut(&str) -> Result<Option<PathBuf>>,
) -> Result<GenomeWideInputs> {
    let absent = |message: String| {
        if species == HUMAN_SPECIES {
            warn!("{message}");
        } else {
            info!("{message}");
        }
    };
    let frequencies = locate(THOUSAND_GENOMES_FILE)?;
    let somatic = locate(&format!("{species}{SOMATIC_SUFFIX}"))?;
    let mut germline = Vec::new();
    for suffix in GERMLINE_SUFFIXES {
        let name = format!("{species}{suffix}");
        match locate(&name)? {
            Some(path) => germline.push(path),
            None => absent(format!(
                "No {name}: a variant only it carries is absent from the cache"
            )),
        }
    }
    if frequencies.is_none() {
        absent(format!(
            "No {THOUSAND_GENOMES_FILE}: the cache carries no 1000 Genomes frequency column"
        ));
    }
    if somatic.is_none() {
        absent("No somatic variation VCF: no record of the cache is flagged somatic".to_string());
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

/// The name of one chromosome's variation VCF in that directory, the layout
/// Ensembl publishes for `homo_sapiens` alone.
pub fn variation_vcf_name(species: &str, chromosome: &str) -> String {
    format!("{species}-chr{chromosome}.vcf.gz")
}

/// The name of the one genome-wide germline variation VCF Ensembl publishes
/// for a species other than `homo_sapiens` (`mus_musculus.vcf.gz`).
pub fn genome_vcf_name(species: &str) -> String {
    format!("{species}.vcf.gz")
}

/// One chromosome's variation VCF, relative to the site root.
fn variation_vcf_path(release: u32, assembly: &str, species: &str, chromosome: &str) -> String {
    format!(
        "{}/{}",
        variation_vcf_dir(release, assembly, species),
        variation_vcf_name(species, chromosome)
    )
}

/// Fetches the per-chromosome variation VCFs into the work directory.
///
/// Fetches each chromosome's VCF and its `.csi` index in parallel using rayon.
/// Returns `Vec<(chromosome, local_path)>`.
fn download_variation_vcfs(
    from: &Downloads<'_>,
    release: u32,
    assembly: &str,
    species: &str,
    chromosomes: &[&str],
) -> Result<Vec<(String, PathBuf)>> {
    info!(
        "Downloading variation VCFs for {} ({} release-{} {})",
        counted(chromosomes.len(), "chromosome"),
        species,
        release,
        assembly
    );

    let (site, work) = (from.site, from.work);
    let results: Result<Vec<(String, PathBuf)>> = chromosomes
        .par_iter()
        .map(|chrom| {
            let vcf_path = variation_vcf_path(release, assembly, species, chrom);
            let vcf_dest = variation_vcf_dest(work, species, assembly, release, chrom);
            work.fetch_cached(site, &vcf_path, &vcf_dest)
                .with_context(|| format!("Failed to download variation VCF for chr{}", chrom))?;
            // The CSI index is optional: some chromosomes have none.
            if let Err(e) = work.fetch_index(site, &vcf_path, &vcf_dest) {
                info!("CSI index not available for chr{}: {}", chrom, e);
            }
            Ok((chrom.to_string(), vcf_dest))
        })
        .collect();

    let mut results = results?;
    results.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(results)
}

/// Fetches the genome-wide germline variation VCF of `species` and its index
/// into the work directory. A species the site publishes no such file for is
/// refused naming the file looked for.
fn download_genome_vcf(
    from: &Downloads<'_>,
    release: u32,
    assembly: &str,
    species: &str,
) -> Result<PathBuf> {
    let name = genome_vcf_name(species);
    let path = format!("{}/{name}", variation_vcf_dir(release, assembly, species));
    let dest = genome_wide_dest(from.work, species, assembly, release, &name);
    match from.work.fetch_cached(from.site, &path, &dest) {
        Ok(()) => {}
        Err(e) if is_not_found(&e) => bail!(
            "Ensembl publishes no variation VCF for {species} on {assembly} at release \
             {release}: looked for {}, the one genome-wide file of a species other than \
             homo_sapiens",
            from.site.address(&path)
        ),
        Err(e) => return Err(e.context(format!("Failed to download {name}"))),
    }
    from.work
        .fetch_index(from.site, &path, &dest)
        .with_context(|| format!("Failed to download the index of {name}"))?;
    Ok(dest)
}

/// Fetches the genome-wide variation VCFs and their indexes into the work
/// directory; a file the site does not publish for this release and assembly
/// (the 1000 Genomes file on GRCh37) is left out.
fn download_genome_wide(
    from: &Downloads<'_>,
    release: u32,
    assembly: &str,
    species: &str,
) -> Result<GenomeWideInputs> {
    let dir = variation_vcf_dir(release, assembly, species);
    genome_wide_inputs(species, |name| {
        let path = format!("{dir}/{name}");
        let dest = genome_wide_dest(from.work, species, assembly, release, name);
        match from.work.fetch_cached(from.site, &path, &dest) {
            Ok(()) => {}
            Err(e) if is_not_found(&e) => return Ok(None),
            Err(e) => return Err(e.context(format!("Failed to download {name}"))),
        }
        from.work
            .fetch_index(from.site, &path, &dest)
            .with_context(|| format!("Failed to download the index of {name}"))?;
        Ok(Some(dest))
    })
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

/// Whether an error is an HTTP 404 from [`download_file`], or a file missing
/// from a local copy of the site's tree.
pub fn is_not_found(error: &anyhow::Error) -> bool {
    error
        .downcast_ref::<HttpStatus>()
        .is_some_and(|status| status.0 == reqwest::StatusCode::NOT_FOUND)
        || error.downcast_ref::<Absent>().is_some()
}

/// Whether an error is an HTTP 403, the answer of a server that withholds a
/// file and of one that answers so, not 404, for a file it does not hold.
pub fn is_forbidden(error: &anyhow::Error) -> bool {
    error
        .downcast_ref::<HttpStatus>()
        .is_some_and(|status| status.0 == reqwest::StatusCode::FORBIDDEN)
}

/// Whether an error is the answer of a server that serves no index page of a
/// directory: HTTP 404, or 403 where listing is switched off.
fn serves_no_index(error: &anyhow::Error) -> bool {
    error.downcast_ref::<HttpStatus>().is_some_and(|status| {
        matches!(
            status.0,
            reqwest::StatusCode::NOT_FOUND | reqwest::StatusCode::FORBIDDEN
        )
    })
}

/// The entries of the directory at `url`, which ends in a slash, read from
/// the server's index page of it; `None` when the server serves none.
pub fn fetch_index(url: &str) -> Result<Option<Vec<String>>> {
    let response = match http_get(url) {
        Ok(response) => response,
        Err(e) if serves_no_index(&e) => return Ok(None),
        Err(e) => return Err(e),
    };
    let html = response
        .text()
        .with_context(|| format!("Failed to read response from {}", url))?;
    Ok(Some(index_entries(&html)))
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

/// The bytes of the file at `url`, or `None` when the server answers 404.
pub fn fetch_bytes(url: &str) -> Result<Option<Vec<u8>>> {
    let response = match http_get(url) {
        Ok(response) => response,
        Err(e) if is_not_found(&e) => return Ok(None),
        Err(e) => return Err(e),
    };
    let bytes = response
        .bytes()
        .with_context(|| format!("Failed to read response from {}", url))?;
    Ok(Some(bytes.to_vec()))
}

/// Download a file from a URL to a local path, streaming the body into
/// `<dest>.part` and renaming it at the end, so an interrupted transfer never
/// leaves a truncated file under the final name (the part file is removed).
/// The directory of `dest` is created once the server answers with the file,
/// so a refused request leaves no directory behind.
pub fn download_file(url: &str, dest: &Path) -> Result<()> {
    info!("Downloading: {}", url);

    let mut response = http_get(url)?;
    if let Some(dir) = dest.parent() {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("Failed to create {}", dir.display()))?;
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

    info!("Downloaded {} ({})", dest.display(), counted(bytes, "byte"));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The URL the Ensembl site serves `relative` at, as `Site::Url` forms it.
    fn on_the_site(relative: &str) -> String {
        format!("{ENSEMBL_SITE}/{relative}")
    }

    #[test]
    fn only_a_404_reads_as_not_found() {
        let not_found = anyhow::Error::new(HttpStatus(reqwest::StatusCode::NOT_FOUND))
            .context("Download failed: x");
        assert!(is_not_found(&not_found));
        let forbidden = anyhow::Error::new(HttpStatus(reqwest::StatusCode::FORBIDDEN))
            .context("Download failed: x");
        assert!(!is_not_found(&forbidden));
        assert!(is_forbidden(&forbidden));
        assert!(!is_forbidden(&not_found));
        assert!(!is_forbidden(&anyhow::anyhow!("connection reset")));
        assert!(!is_not_found(&anyhow::anyhow!("connection reset")));
        let absent = anyhow::Error::new(Absent(PathBuf::from("pub/release-116/x.vcf.gz")));
        assert!(is_not_found(&absent.context("Failed to download x")));
    }

    /// The work directory of a build into an output directory is
    /// `tmp/vep_cache_builder/<species>_<release>_<assembly>` under it.
    /// Removing it removes `tmp/vep_cache_builder` with the work directory a
    /// build of another release, assembly and species kept there, and leaves
    /// a file that stands beside it in `tmp`. A failed build that fetched
    /// nothing names no directory while another build's is kept.
    #[test]
    fn a_build_into_an_output_directory_clears_the_work_directory_of_every_build() {
        let root = tempfile::tempdir().unwrap();
        let output_dir = root.path().join("cache");
        let work = WorkDir::of_output_dir(&output_dir, "homo_sapiens", 116, "GRCh38");
        assert_eq!(
            work.path(),
            output_dir.join("tmp/vep_cache_builder/homo_sapiens_116_GRCh38")
        );
        let other = WorkDir::of_output_dir(&output_dir, "mus_musculus", 115, "GRCm39");
        assert_eq!(
            other.path(),
            output_dir.join("tmp/vep_cache_builder/mus_musculus_115_GRCm39")
        );
        let keep = |work: &WorkDir| {
            std::fs::create_dir_all(work.path()).unwrap();
            std::fs::write(gff3_dest(work, "mus_musculus", "GRCm39", 115), b"kept").unwrap();
        };

        keep(&other);
        let error = work.kept(anyhow::anyhow!("no site"));
        assert_eq!(format!("{error:#}"), "no site");
        work.remove().unwrap();
        assert!(!output_dir.join("tmp").exists());
        assert!(output_dir.is_dir());

        keep(&other);
        std::fs::write(output_dir.join("tmp/notes.txt"), b"notes").unwrap();
        work.remove().unwrap();
        assert!(!output_dir.join("tmp/vep_cache_builder").exists());
        assert!(output_dir.join("tmp/notes.txt").is_file());
    }

    /// The work directory is `tmp/<species>_<release>_<assembly>` under the
    /// root it is given. Each of the four downloads of a `homo_sapiens` build
    /// lands in it, under a name that carries its species, assembly and
    /// release, and none under `std::env::temp_dir()`.
    #[test]
    fn every_download_lands_in_the_work_directory() {
        let work = WorkDir::of_release(
            Path::new("/caches/work-directory-test"),
            "homo_sapiens",
            116,
            "GRCh38",
        );
        assert_eq!(
            work.path(),
            Path::new("/caches/work-directory-test/tmp/homo_sapiens_116_GRCh38")
        );
        let dests = [
            gff3_dest(&work, "homo_sapiens", "GRCh38", 116),
            protein_fasta_dest(&work, "homo_sapiens", "GRCh38", 116),
            variation_vcf_dest(&work, "homo_sapiens", "GRCh38", 116, "21"),
            genome_wide_dest(&work, "homo_sapiens", "GRCh38", 116, THOUSAND_GENOMES_FILE),
        ];
        let names = [
            "vep_cache_builder_homo_sapiens_GRCh38_116.gff3.gz",
            "vep_cache_builder_homo_sapiens_GRCh38_116.pep.all.fa.gz",
            "vep_cache_builder_homo_sapiens_GRCh38_116_chr21.vcf.gz",
            "vep_cache_builder_homo_sapiens_GRCh38_116_1000GENOMES-phase_3.vcf.gz",
        ];
        for (dest, name) in dests.iter().zip(names) {
            assert_eq!(dest, &work.path().join(name));
            assert!(
                !dest.starts_with(std::env::temp_dir()),
                "{}",
                dest.display()
            );
        }
    }

    /// The GFF3 of the Ensembl release tree, relative to the site root.
    fn gff3_path(release: u32, assembly: &str, species: &str) -> String {
        let (dir, name) = Layout::ensembl(assembly).gff3_location(release, assembly, species);
        format!("{dir}/{name}")
    }

    /// The protein FASTA of the Ensembl release tree, relative to the site
    /// root.
    fn protein_fasta_path(release: u32, assembly: &str, species: &str) -> String {
        let (dir, name) =
            Layout::ensembl(assembly).protein_fasta_location(release, assembly, species);
        format!("{dir}/{name}")
    }

    /// A local copy of the site's tree holding the GFF3, the protein FASTA,
    /// two per-chromosome VCFs (one with an index) and the somatic file of
    /// `homo_sapiens` GRCh38 release 116, with the directory of its core
    /// database, by which the genome is found in a tree that publishes no
    /// species list. Each file's bytes are its own path.
    fn local_site(root: &Path) -> Site {
        let core = Layout::Release.mysql_dir(116, "homo_sapiens_core_116_38");
        std::fs::create_dir_all(root.join(core.unwrap())).unwrap();
        let vcf_dir = variation_vcf_dir(116, "GRCh38", "homo_sapiens");
        let somatic = format!("{vcf_dir}/homo_sapiens{SOMATIC_SUFFIX}");
        let chr21 = variation_vcf_path(116, "GRCh38", "homo_sapiens", "21");
        for relative in [
            gff3_path(116, "GRCh38", "homo_sapiens"),
            protein_fasta_path(116, "GRCh38", "homo_sapiens"),
            variation_vcf_path(116, "GRCh38", "homo_sapiens", "22"),
            format!("{chr21}.csi"),
            chr21,
            format!("{somatic}.csi"),
            somatic,
        ] {
            let path = root.join(&relative);
            std::fs::create_dir_all(path.parent().unwrap()).unwrap();
            std::fs::write(&path, &relative).unwrap();
        }
        Site::Dir(root.to_path_buf())
    }

    /// A build that names no local file fetches each input into its work
    /// directory, which the first fetch creates: the GFF3, the protein FASTA,
    /// each chromosome's VCF with the index the site holds, and the
    /// genome-wide files the site holds. A second resolution reuses the
    /// files, and removing the work directory removes its `tmp` with it.
    #[test]
    fn a_missing_input_is_fetched_into_the_work_directory() {
        let root = tempfile::tempdir().unwrap();
        let site = local_site(&root.path().join("site"));
        let work = WorkDir::of_release(&root.path().join("caches"), "homo_sapiens", 116, "GRCh38");
        let published = Published::new(vec![site.clone()], "homo_sapiens", "GRCh38", 116);
        let from = Downloads {
            published: &published,
            site: &site,
            work: &work,
        };
        assert!(!work.path().exists());

        let gff3 = resolve_gff3(None, 116, "GRCh38", "homo_sapiens", &from).unwrap();
        assert_eq!(gff3, gff3_dest(&work, "homo_sapiens", "GRCh38", 116));
        assert_eq!(
            std::fs::read_to_string(&gff3).unwrap(),
            gff3_path(116, "GRCh38", "homo_sapiens")
        );
        let fasta = resolve_protein_fasta(None, 116, "GRCh38", "homo_sapiens", &from).unwrap();
        assert_eq!(
            fasta,
            protein_fasta_dest(&work, "homo_sapiens", "GRCh38", 116)
        );
        assert_eq!(
            std::fs::read_to_string(&fasta).unwrap(),
            protein_fasta_path(116, "GRCh38", "homo_sapiens")
        );

        let chromosomes = ["21".to_string(), "22".to_string()];
        let resolve_variations = || {
            resolve_variation_inputs(
                &None,
                116,
                "GRCh38",
                "homo_sapiens",
                Some(&chromosomes),
                &from,
            )
            .unwrap()
        };
        let variations = resolve_variations();
        let chr21 = variation_vcf_dest(&work, "homo_sapiens", "GRCh38", 116, "21");
        let chr22 = variation_vcf_dest(&work, "homo_sapiens", "GRCh38", 116, "22");
        assert_eq!(
            variations.per_chromosome,
            [
                ("21".to_string(), GermlineVcf::Chromosome(chr21.clone())),
                ("22".to_string(), GermlineVcf::Chromosome(chr22.clone()))
            ]
        );
        assert_eq!(
            std::fs::read_to_string(&chr21).unwrap(),
            variation_vcf_path(116, "GRCh38", "homo_sapiens", "21")
        );
        assert!(chr21.with_extension("gz.csi").is_file());
        assert!(!chr22.with_extension("gz.csi").exists());
        let somatic = genome_wide_dest(
            &work,
            "homo_sapiens",
            "GRCh38",
            116,
            "homo_sapiens_somatic.vcf.gz",
        );
        assert_eq!(variations.genome_wide.somatic.as_deref(), Some(&*somatic));
        assert!(somatic.with_extension("gz.csi").is_file());
        assert!(variations.genome_wide.frequencies.is_none());
        assert!(variations.genome_wide.germline.is_empty());
        for entry in std::fs::read_dir(work.path()).unwrap() {
            let name = entry.unwrap().file_name();
            assert!(
                name.to_string_lossy()
                    .starts_with("vep_cache_builder_homo_sapiens_GRCh38_116"),
                "{name:?}"
            );
        }

        // The site is gone, so a second resolution finds every file in the
        // work directory.
        std::fs::remove_dir_all(root.path().join("site")).unwrap();
        std::fs::create_dir(root.path().join("site")).unwrap();
        assert_eq!(
            resolve_gff3(None, 116, "GRCh38", "homo_sapiens", &from).unwrap(),
            gff3
        );
        assert_eq!(resolve_variations().per_chromosome.len(), 2);

        work.remove().unwrap();
        assert!(!root.path().join("caches/tmp").exists());
        assert!(root.path().join("caches").is_dir());
        work.remove().unwrap();
    }

    /// A required file the site lacks fails the resolution naming it; the
    /// error of a failed build names the work directory once a download has
    /// placed a file in it and is left as it is before that. A work
    /// directory holding no file, at any depth, is named by no error and is
    /// removed with the `tmp` it leaves empty.
    #[test]
    fn a_failed_build_names_the_work_directory_it_keeps() {
        let root = tempfile::tempdir().unwrap();
        let site_dir = root.path().join("site");
        let site = local_site(&site_dir);
        std::fs::remove_file(site_dir.join(protein_fasta_path(116, "GRCh38", "homo_sapiens")))
            .unwrap();
        let work = WorkDir::of_release(&root.path().join("caches"), "homo_sapiens", 116, "GRCh38");
        let published = Published::new(vec![site.clone()], "homo_sapiens", "GRCh38", 116);
        let from = Downloads {
            published: &published,
            site: &site,
            work: &work,
        };

        let error = resolve_protein_fasta(None, 116, "GRCh38", "homo_sapiens", &from).unwrap_err();
        assert!(
            format!("{error:#}").ends_with("Homo_sapiens.GRCh38.pep.all.fa.gz is not on the site"),
            "{error:#}"
        );
        assert!(!work.path().exists());
        assert_eq!(format!("{:#}", work.kept(error)), {
            let absent = site_dir.join(protein_fasta_path(116, "GRCh38", "homo_sapiens"));
            format!("{} is not on the site", absent.display())
        });

        std::fs::create_dir_all(work.path().join("inputs")).unwrap();
        std::fs::create_dir_all(work.path().join("cache")).unwrap();
        let error = anyhow::anyhow!("No transcript was read").context("Failed to build");
        assert_eq!(
            format!("{:#}", work.kept(error)),
            "Failed to build: No transcript was read"
        );
        assert!(!root.path().join("caches/tmp").exists());

        std::fs::create_dir_all(work.path().join("inputs")).unwrap();
        std::fs::write(work.path().join("inputs/kept.gz"), b"kept").unwrap();
        let error = anyhow::anyhow!("No transcript was read").context("Failed to build");
        assert!(format!("{:#}", work.kept(error)).ends_with("is kept for the next run"));
        work.remove().unwrap();

        resolve_gff3(None, 116, "GRCh38", "homo_sapiens", &from).unwrap();
        let error = anyhow::anyhow!("No transcript was read").context("Failed to build");
        assert_eq!(
            format!("{:#}", work.kept(error)),
            format!(
                "Failed to build: No transcript was read; the work directory {} is kept for the next run",
                work.path().display()
            )
        );
        // Another build's files keep the shared `tmp` in place.
        std::fs::create_dir(root.path().join("caches/tmp/other")).unwrap();
        work.remove().unwrap();
        assert!(!work.path().exists());
        assert!(root.path().join("caches/tmp/other").is_dir());
    }

    /// A failed build into an output directory whose own work directory
    /// holds no file removes that directory alone: the download a failed
    /// build of another genome kept beside it stays for that genome's next
    /// run, and that build's own failure still names its directory. Alone in
    /// the output directory, the empty work directory goes with the
    /// `tmp/vep_cache_builder` and the `tmp` it leaves empty.
    #[test]
    fn a_failed_build_that_kept_no_file_leaves_another_builds_download() {
        let root = tempfile::tempdir().unwrap();
        let output_dir = root.path().join("cache");
        let work = WorkDir::of_output_dir(&output_dir, "homo_sapiens", 116, "GRCh38");
        let other = WorkDir::of_output_dir(&output_dir, "mus_musculus", 115, "GRCm39");
        let kept_download = gff3_dest(&other, "mus_musculus", "GRCm39", 115);
        std::fs::create_dir_all(other.path()).unwrap();
        std::fs::write(&kept_download, b"kept").unwrap();

        std::fs::create_dir_all(work.path().join("inputs")).unwrap();
        let error = work.kept(anyhow::anyhow!("no site"));
        assert_eq!(format!("{error:#}"), "no site");
        assert!(!work.path().exists());
        assert_eq!(std::fs::read(&kept_download).unwrap(), b"kept");
        assert_eq!(
            format!("{:#}", other.kept(anyhow::anyhow!("no site"))),
            format!(
                "no site; the work directory {} is kept for the next run",
                other.path().display()
            )
        );

        std::fs::remove_dir_all(other.path()).unwrap();
        std::fs::create_dir_all(work.path().join("inputs")).unwrap();
        let error = work.kept(anyhow::anyhow!("no site"));
        assert_eq!(format!("{error:#}"), "no site");
        assert!(!output_dir.join("tmp").exists());
        assert!(output_dir.is_dir());
    }

    /// A directory's index page is not served when the server answers 404 or
    /// 403; any other failure is a failed read.
    #[test]
    fn a_404_or_a_403_reads_as_no_index_page() {
        for status in [
            reqwest::StatusCode::NOT_FOUND,
            reqwest::StatusCode::FORBIDDEN,
        ] {
            let error = anyhow::Error::new(HttpStatus(status)).context("Download failed: x");
            assert!(serves_no_index(&error), "{status}");
        }
        let unavailable = anyhow::Error::new(HttpStatus(reqwest::StatusCode::SERVICE_UNAVAILABLE))
            .context("Download failed: x");
        assert!(!serves_no_index(&unavailable));
        assert!(!serves_no_index(&anyhow::anyhow!("connection reset")));
    }

    /// The URL of a file on the Ensembl site, from its directory and name.
    fn ensembl_url((dir, name): (String, String)) -> String {
        format!("{ENSEMBL_SITE}/{dir}/{name}")
    }

    #[test]
    fn test_gff3_url_grch38() {
        let url =
            ensembl_url(Layout::ensembl("GRCh38").gff3_location(115, "GRCh38", "homo_sapiens"));
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz"
        );
    }

    #[test]
    fn test_gff3_url_grch37() {
        let layout = Layout::ensembl("GRCh37");
        let url = ensembl_url(layout.gff3_location(115, "GRCh37", "homo_sapiens"));
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/grch37/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz"
        );
        let url = ensembl_url(layout.gff3_location(87, "GRCh37", "homo_sapiens"));
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz"
        );
    }

    #[test]
    fn test_protein_fasta_url_grch38() {
        let url = ensembl_url(Layout::ensembl("GRCh38").protein_fasta_location(
            115,
            "GRCh38",
            "homo_sapiens",
        ));
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
        );
    }

    #[test]
    fn test_protein_fasta_url_grch37() {
        let url = ensembl_url(Layout::ensembl("GRCh37").protein_fasta_location(
            115,
            "GRCh37",
            "homo_sapiens",
        ));
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
        assert_eq!(Layout::ensembl("GRCh38"), Layout::Release);
        assert_eq!(Layout::ensembl("GRCm39"), Layout::Release);
        assert_eq!(Layout::ensembl("GRCh37"), Layout::Grch37);
        assert_eq!(release_dir(116, "GRCh38"), "pub/release-116");
        assert_eq!(release_dir(116, "GRCh37"), "pub/grch37/release-116");
        assert_eq!(
            Layout::Release.gff3_location(116, "GRCh38", "homo_sapiens"),
            (
                "pub/release-116/gff3/homo_sapiens".to_string(),
                "Homo_sapiens.GRCh38.116.gff3.gz".to_string()
            )
        );
        assert_eq!(
            Layout::Grch37.gtf_location(116, "GRCh37", "homo_sapiens"),
            (
                "pub/grch37/release-116/gtf/homo_sapiens".to_string(),
                "Homo_sapiens.GRCh37.87.gtf.gz".to_string()
            )
        );
        assert_eq!(
            Layout::Release.protein_fasta_location(116, "GRCm39", "mus_musculus"),
            (
                "pub/release-116/fasta/mus_musculus/pep".to_string(),
                "Mus_musculus.GRCm39.pep.all.fa.gz".to_string()
            )
        );
        assert_eq!(
            Layout::Grch37.genome_fasta_dirs(116, "homo_sapiens"),
            [
                "pub/grch37/release-116/fasta/homo_sapiens/dna_index".to_string(),
                "pub/grch37/release-116/fasta/homo_sapiens/dna".to_string()
            ]
        );
        assert_eq!(
            Layout::Release
                .mysql_dir(116, "homo_sapiens_variation_116_38")
                .unwrap(),
            "pub/release-116/mysql/homo_sapiens_variation_116_38"
        );
        assert_eq!(
            Layout::Release.variation_vcf_dir(116, "homo_sapiens"),
            Some("pub/release-116/variation/vcf/homo_sapiens".to_string())
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

    /// The GVF directory of a species stands beside the directory of its
    /// variation VCFs in each tree that has one, the GRCh37 tree's in that
    /// tree, and the failed variants' file is named for the species.
    #[test]
    fn variation_gvf_directory_of_each_tree() {
        assert_eq!(
            Layout::Release
                .variation_gvf_dir(116, "homo_sapiens")
                .as_deref(),
            Some("pub/release-116/variation/gvf/homo_sapiens")
        );
        assert_eq!(
            Layout::Grch37
                .variation_gvf_dir(116, "homo_sapiens")
                .as_deref(),
            Some("pub/grch37/release-116/variation/gvf/homo_sapiens")
        );
        let division = Layout::Division {
            division: "fungi".to_string(),
        };
        assert_eq!(
            division
                .variation_gvf_dir(63, "saccharomyces_cerevisiae")
                .as_deref(),
            Some("pub/release-63/fungi/variation/gvf/saccharomyces_cerevisiae")
        );
        let collection = Layout::Collection {
            division: "bacteria".to_string(),
            collection: "bacteria_7_collection".to_string(),
        };
        let rapid = Layout::RapidRelease {
            organism: "Organismus_fictus".to_string(),
            accession: "GCA_000000003.1".to_string(),
            gene_set: "2024_01".to_string(),
        };
        for layout in [collection, rapid] {
            assert_eq!(layout.variation_gvf_dir(63, "species_ficta"), None);
            assert_eq!(layout.variation_vcf_dir(63, "species_ficta"), None);
        }
        assert_eq!(
            failed_gvf_name("saccharomyces_cerevisiae"),
            "saccharomyces_cerevisiae_failed.gvf.gz"
        );
    }

    /// The regulatory annotation of a species stands under the release
    /// directory, for GRCh37 as for any other assembly of the release tree,
    /// and under the division's directory for a species of a division; a
    /// genome of a collection and an assembly of the rapid release have
    /// none.
    #[test]
    fn regulation_directory_of_each_tree() {
        for layout in [Layout::Release, Layout::Grch37] {
            assert_eq!(
                layout.regulation_dir(116, "homo_sapiens").as_deref(),
                Some("pub/release-116/regulation/homo_sapiens"),
                "{layout:?}"
            );
        }
        let division = Layout::Division {
            division: "plants".to_string(),
        };
        assert_eq!(
            division.regulation_dir(63, "planta_ficta").as_deref(),
            Some("pub/release-63/plants/regulation/planta_ficta")
        );
        let collection = Layout::Collection {
            division: "bacteria".to_string(),
            collection: "bacteria_7_collection".to_string(),
        };
        assert_eq!(
            collection.regulation_dir(63, "bacterium_fictum_gca_000000001"),
            None
        );
        let rapid_release = Layout::RapidRelease {
            organism: "Homo_sapiens".to_string(),
            accession: "GCA_000000001.1".to_string(),
            gene_set: "2022_07".to_string(),
        };
        assert_eq!(
            rapid_release.regulation_dir(107, "homo_sapiens_gca000000001v1"),
            None
        );
    }

    /// A directory holding every table reads no site. With a core table
    /// missing, the lookup's error stops the build. With the core tables there
    /// and the variation tables alone missing, a genome no site publishes
    /// builds without them, and any other error still stops the build.
    #[test]
    fn core_dumps_ask_for_the_databases_only_when_a_table_is_missing() {
        let dir = tempfile::tempdir().unwrap();
        for table in CORE_TABLES.iter().chain(&VARIATION_TABLES) {
            std::fs::write(dir.path().join(format!("{table}.txt")), "").unwrap();
        }
        let unpublished = || -> Result<DatabaseUrls> {
            Err(anyhow::Error::new(NotPublished("not there".to_string())))
        };
        let unreachable = || -> Result<DatabaseUrls> { bail!("connection reset") };
        resolve_core_dumps(dir.path(), &|| panic!("no table is missing")).unwrap();

        for table in VARIATION_TABLES {
            std::fs::remove_file(dir.path().join(format!("{table}.txt"))).unwrap();
        }
        resolve_core_dumps(dir.path(), &unpublished).unwrap();
        let error = resolve_core_dumps(dir.path(), &unreachable).unwrap_err();
        assert_eq!(error.to_string(), "connection reset");

        std::fs::remove_file(dir.path().join("gene.txt")).unwrap();
        let error = resolve_core_dumps(dir.path(), &unpublished).unwrap_err();
        assert!(is_not_published(&error.context("Failed to resolve")));
    }

    /// A species of an Ensembl Genomes division: the release tree's
    /// directories under `pub/release-<N>/<division>/`, every file named for
    /// that release.
    #[test]
    fn division_tree_locations() {
        let layout = Layout::Division {
            division: "plants".to_string(),
        };
        assert_eq!(layout.root(63), "pub/release-63/plants");
        assert_eq!(
            layout.gff3_location(63, "Pf1.0", "planta_ficta"),
            (
                "pub/release-63/plants/gff3/planta_ficta".to_string(),
                "Planta_ficta.Pf1.0.63.gff3.gz".to_string()
            )
        );
        assert_eq!(
            layout.gtf_location(63, "Pf1.0", "planta_ficta"),
            (
                "pub/release-63/plants/gtf/planta_ficta".to_string(),
                "Planta_ficta.Pf1.0.63.gtf.gz".to_string()
            )
        );
        assert_eq!(
            layout.protein_fasta_location(63, "Pf1.0", "planta_ficta"),
            (
                "pub/release-63/plants/fasta/planta_ficta/pep".to_string(),
                "Planta_ficta.Pf1.0.pep.all.fa.gz".to_string()
            )
        );
        assert_eq!(
            layout.genome_fasta_dirs(63, "planta_ficta"),
            [
                "pub/release-63/plants/fasta/planta_ficta/dna_index".to_string(),
                "pub/release-63/plants/fasta/planta_ficta/dna".to_string()
            ]
        );
        assert_eq!(
            layout
                .mysql_dir(63, "planta_ficta_core_63_116_4")
                .as_deref(),
            Some("pub/release-63/plants/mysql/planta_ficta_core_63_116_4")
        );
        assert_eq!(
            layout.variation_vcf_dir(63, "planta_ficta").as_deref(),
            Some("pub/release-63/plants/variation/vcf/planta_ficta")
        );
    }

    /// A genome of a collection: the division's tree with the collection's
    /// directory above the species', one database for the collection and no
    /// variation VCF directory.
    #[test]
    fn collection_tree_locations() {
        let layout = Layout::Collection {
            division: "bacteria".to_string(),
            collection: "bacteria_7_collection".to_string(),
        };
        let species = "bacterium_fictum_gca_000000001";
        assert_eq!(
            layout.gff3_location(63, "ASM1v1", species),
            (
                "pub/release-63/bacteria/gff3/bacteria_7_collection/bacterium_fictum_gca_000000001"
                    .to_string(),
                "Bacterium_fictum_gca_000000001.ASM1v1.63.gff3.gz".to_string()
            )
        );
        assert_eq!(
            layout.gtf_location(63, "ASM1v1", species).0,
            "pub/release-63/bacteria/gtf/bacteria_7_collection/bacterium_fictum_gca_000000001"
        );
        assert_eq!(
            layout.protein_fasta_location(63, "ASM1v1", species),
            (
                "pub/release-63/bacteria/fasta/bacteria_7_collection/bacterium_fictum_gca_000000001/pep"
                    .to_string(),
                "Bacterium_fictum_gca_000000001.ASM1v1.pep.all.fa.gz".to_string()
            )
        );
        assert_eq!(
            layout.genome_fasta_dirs(63, species)[1],
            "pub/release-63/bacteria/fasta/bacteria_7_collection/bacterium_fictum_gca_000000001/dna"
        );
        assert_eq!(
            layout
                .mysql_dir(63, "bacteria_7_collection_core_63_116_1")
                .as_deref(),
            Some("pub/release-63/bacteria/mysql/bacteria_7_collection_core_63_116_1")
        );
        assert_eq!(layout.variation_vcf_dir(63, species), None);
    }

    /// An assembly of the rapid release: every gene-set file named for the
    /// organism, the accession and the dated gene set, one genome directory,
    /// and no table dump.
    #[test]
    fn rapid_release_tree_locations() {
        let layout = Layout::RapidRelease {
            organism: "Homo_sapiens".to_string(),
            accession: "GCA_000000001.1".to_string(),
            gene_set: "2022_07".to_string(),
        };
        let species = "homo_sapiens_gca000000001v1";
        let root = "pub/rapid-release/species/Homo_sapiens/GCA_000000001.1/ensembl";
        assert_eq!(layout.root(107), root);
        assert_eq!(
            layout.gff3_location(107, "Made-up1.0", species),
            (
                format!("{root}/geneset/2022_07"),
                "Homo_sapiens-GCA_000000001.1-2022_07-genes.gff3.gz".to_string()
            )
        );
        assert_eq!(
            layout.gtf_location(107, "Made-up1.0", species),
            (
                format!("{root}/geneset/2022_07"),
                "Homo_sapiens-GCA_000000001.1-2022_07-genes.gtf.gz".to_string()
            )
        );
        assert_eq!(
            layout.protein_fasta_location(107, "Made-up1.0", species),
            (
                format!("{root}/geneset/2022_07"),
                "Homo_sapiens-GCA_000000001.1-2022_07-pep.fa.gz".to_string()
            )
        );
        assert_eq!(
            layout.genome_fasta_dirs(107, species),
            [format!("{root}/genome")]
        );
        assert_eq!(layout.mysql_dir(107, "any"), None);
        assert_eq!(layout.variation_vcf_dir(107, species), None);
        assert_eq!(
            rapid_release_name(species),
            Some(("Homo_sapiens".to_string(), "GCA_000000001.1".to_string()))
        );
        for other in [
            "homo_sapiens",
            "bacterium_fictum_gca_000000001",
            "formica_ficta_gca000000001v1rs",
            "_gca1v1",
        ] {
            assert_eq!(rapid_release_name(other), None, "{other}");
        }
    }

    /// A species list is read by its header's column names, whatever their
    /// order and whether or not a row ends in a tab: the row gives the
    /// division and the core database, which place a vertebrate in the release
    /// tree, a species with a database of its own in its division's tree and
    /// any other in its collection.
    #[test]
    fn species_list_places_a_species() {
        let list = "#name\tspecies\tdivision\ttaxonomy_id\tassembly\tcore_db\tspecies_id\n\
                    Rat\trattus_norvegicus\tEnsemblVertebrates\t10116\tGRCr8\trattus_norvegicus_core_116_1\t1\t\n\
                    A plant\tplanta_ficta\tEnsemblPlants\t1\tPf1.0\tplanta_ficta_core_63_116_4\t1\n\
                    A bacterium\tbacterium_fictum_gca_000000001\tEnsemblBacteria\t2\tASM1v1\tbacteria_7_collection_core_63_116_1\t12\n\
                    Cut short\tformica_ficta\tEnsemblMetazoa\n";
        let rat = species_entry(list, "rattus_norvegicus").unwrap().unwrap();
        assert_eq!(rat.core_db, "rattus_norvegicus_core_116_1");
        assert_eq!(rat.layout("rattus_norvegicus"), Layout::Release);
        let plant = species_entry(list, "planta_ficta").unwrap().unwrap();
        assert_eq!(
            plant,
            SpeciesEntry {
                division: "EnsemblPlants".to_string(),
                core_db: "planta_ficta_core_63_116_4".to_string()
            }
        );
        assert_eq!(
            plant.layout("planta_ficta"),
            Layout::Division {
                division: "plants".to_string()
            }
        );
        let bacterium = species_entry(list, "bacterium_fictum_gca_000000001")
            .unwrap()
            .unwrap();
        assert_eq!(
            bacterium.layout("bacterium_fictum_gca_000000001"),
            Layout::Collection {
                division: "bacteria".to_string(),
                collection: "bacteria_7_collection".to_string()
            }
        );
        assert_eq!(species_entry(list, "planta").unwrap(), None);
        assert_eq!(species_entry(list, "A plant").unwrap(), None);
        let error = species_entry(list, "formica_ficta")
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("row of formica_ficta is cut short"),
            "{error}"
        );
        let error = species_entry("#name\tspecies\tdivision\n", "planta_ficta")
            .unwrap_err()
            .to_string();
        assert!(error.contains("no core_db column"), "{error}");
    }

    /// The variation database carries the core database's name with
    /// `variation` for `core`, whatever numbers end it; a collection's
    /// database, not named for the species, has none beside it.
    #[test]
    fn variation_database_follows_the_core_database() {
        assert_eq!(
            variation_database("homo_sapiens_core_116_38", "homo_sapiens").as_deref(),
            Some("homo_sapiens_variation_116_38")
        );
        assert_eq!(
            variation_database("planta_ficta_core_63_116_4", "planta_ficta").as_deref(),
            Some("planta_ficta_variation_63_116_4")
        );
        assert_eq!(
            variation_database(
                "bacteria_7_collection_core_63_116_1",
                "bacterium_fictum_gca_000000001"
            ),
            None
        );
        assert_eq!(
            Layout::Grch37
                .mysql_dir(116, "homo_sapiens_core_116_37")
                .as_deref(),
            Some("pub/grch37/release-116/mysql/homo_sapiens_core_116_37")
        );
    }

    /// The funcgen database carries the core database's name with `funcgen`
    /// for `core`, whatever numbers end it, and stands in the tree the core
    /// database stands in, GRCh37's in the GRCh37 tree; a collection's
    /// database, not named for the species, has none beside it.
    #[test]
    fn funcgen_database_follows_the_core_database() {
        assert_eq!(
            funcgen_database("homo_sapiens_core_116_38", "homo_sapiens").as_deref(),
            Some("homo_sapiens_funcgen_116_38")
        );
        assert_eq!(
            funcgen_database("planta_ficta_core_63_116_4", "planta_ficta").as_deref(),
            Some("planta_ficta_funcgen_63_116_4")
        );
        assert_eq!(
            funcgen_database(
                "bacteria_7_collection_core_63_116_1",
                "bacterium_fictum_gca_000000001"
            ),
            None
        );
        let grch37 = funcgen_database("homo_sapiens_core_116_37", "homo_sapiens").unwrap();
        assert_eq!(
            Layout::Grch37.mysql_dir(116, &grch37).as_deref(),
            Some("pub/grch37/release-116/mysql/homo_sapiens_funcgen_116_37")
        );
    }

    /// The core database is picked out of a `mysql/` listing by the species
    /// and the release, whatever number ends it: a database whose number is
    /// not the assembly's trailing digits is found, and one of a species whose
    /// name the species' own begins, of another kind or of another release is
    /// not.
    #[test]
    fn core_database_is_found_in_a_mysql_listing() {
        let entries: Vec<String> = [
            "rattus_norvegicus_core_116_1",
            "rattus_norvegicus_variation_116_1",
            "rattus_norvegicus_otherfeatures_116_1",
            "mus_musculus_core_116_39",
            "mus_musculus_129s1svimj_core_116_1",
            "bos_taurus_core_115_2",
            "bos_taurus_core_116_x",
        ]
        .iter()
        .map(|entry| entry.to_string())
        .collect();
        assert_eq!(
            core_databases(&entries, "rattus_norvegicus", 116),
            ["rattus_norvegicus_core_116_1"]
        );
        assert_eq!(
            core_databases(&entries, "mus_musculus", 116),
            ["mus_musculus_core_116_39"]
        );
        assert_eq!(
            core_databases(&entries, "mus_musculus_129s1svimj", 116),
            ["mus_musculus_129s1svimj_core_116_1"]
        );
        assert!(core_databases(&entries, "bos_taurus", 116).is_empty());

        assert_eq!(
            numbered_core_database("homo_sapiens", "GRCh38", 116).as_deref(),
            Some("homo_sapiens_core_116_38")
        );
        assert_eq!(
            numbered_core_database("homo_sapiens", "GRCh37", 116).as_deref(),
            Some("homo_sapiens_core_116_37")
        );
        assert_eq!(numbered_core_database("bos_taurus", "ARS-UCD", 116), None);
    }

    /// A directory index page yields the names its links point at inside the
    /// directory, a directory's without its slash; the parent, the sort
    /// orders and a link elsewhere are not entries.
    #[test]
    fn index_entries_are_the_links_to_children() {
        let html = r#"<h1>Index of /pub/release-116/mysql</h1>
<tr><th><a href="?C=N;O=D">Name</a></th><th><a href="?C=M;O=A">Last modified</a></th></tr>
<tr><td><a href="/pub/release-116/">Parent Directory</a></td></tr>
<tr><td><a href="../">Up</a></td></tr>
<tr><td><a href="rattus_norvegicus_core_116_1/">rattus_norvegicus_core_116_1/</a></td></tr>
<tr><td><a href="CHECKSUMS">CHECKSUMS</a></td></tr>
<tr><td><a href="https://example.org/">elsewhere</a></td></tr>"#;
        assert_eq!(
            index_entries(html),
            ["rattus_norvegicus_core_116_1", "CHECKSUMS"]
        );
        assert!(index_entries("no links").is_empty());
    }

    #[test]
    fn test_capitalize_species() {
        assert_eq!(capitalize_species("homo_sapiens"), "Homo_sapiens");
        assert_eq!(capitalize_species("mus_musculus"), "Mus_musculus");
        assert_eq!(capitalize_species(""), "");
    }

    #[test]
    fn test_variation_vcf_url_grch38() {
        let url = on_the_site(&variation_vcf_path(115, "GRCh38", "homo_sapiens", "21"));
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/release-115/variation/vcf/homo_sapiens/homo_sapiens-chr21.vcf.gz"
        );
    }

    #[test]
    fn test_variation_vcf_url_grch37() {
        let url = on_the_site(&variation_vcf_path(115, "GRCh37", "homo_sapiens", "21"));
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/grch37/release-115/variation/vcf/homo_sapiens/homo_sapiens-chr21.vcf.gz"
        );
    }

    #[test]
    fn test_variation_vcf_url_x_chromosome() {
        let url = on_the_site(&variation_vcf_path(115, "GRCh38", "homo_sapiens", "X"));
        assert_eq!(
            url,
            "https://ftp.ensembl.org/pub/release-115/variation/vcf/homo_sapiens/homo_sapiens-chrX.vcf.gz"
        );
    }

    #[test]
    fn test_variation_vcf_url_mt() {
        let url = on_the_site(&variation_vcf_path(115, "GRCh37", "homo_sapiens", "MT"));
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

    /// Without `--chromosomes`, `homo_sapiens` reads the 25 chromosomes
    /// Ensembl publishes per-chromosome VCFs for. Another species has no
    /// per-chromosome default, since Ensembl publishes one genome-wide file
    /// for it: that file is the last of its genome-wide names, and is not
    /// among human's. A download fetches human's file of each chromosome
    /// named, else of each default one, and another species' one file
    /// whatever chromosomes are named.
    #[test]
    fn per_chromosome_files_are_the_default_for_human_alone() {
        let human = default_chromosomes("homo_sapiens");
        assert_eq!(human.len(), 25);
        assert_eq!(human[0], "1");
        assert_eq!(human[22], "X");
        assert_eq!(human[24], "MT");
        assert!(default_chromosomes("drosophila_melanogaster").is_empty());
        let named = ["21".to_string(), "X".to_string()];
        assert_eq!(
            germline_vcf_names("homo_sapiens", Some(&named)),
            ["homo_sapiens-chr21.vcf.gz", "homo_sapiens-chrX.vcf.gz"]
        );
        assert_eq!(germline_vcf_names("homo_sapiens", None).len(), 25);
        for chromosomes in [Some(&named[..]), None] {
            assert_eq!(
                germline_vcf_names("mus_musculus", chromosomes),
                ["mus_musculus.vcf.gz"]
            );
        }
        assert_eq!(genome_vcf_name("mus_musculus"), "mus_musculus.vcf.gz");
        assert_eq!(
            genome_wide_names("mus_musculus").last().map(String::as_str),
            Some("mus_musculus.vcf.gz")
        );
        assert!(!genome_wide_names("homo_sapiens").contains(&genome_vcf_name("homo_sapiens")));
    }

    /// Writes `records` under `name` in `dir` as a BGZF VCF with no index.
    fn write_vcf(dir: &Path, name: &str, records: &str) -> PathBuf {
        let path = dir.join(name);
        let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(&path).unwrap());
        write!(
            writer,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"
        )
        .unwrap();
        writer.finish().unwrap();
        path
    }

    /// A local directory is read by what it holds. Per-chromosome files
    /// (`*-chr*.vcf.gz`) are each a chromosome's own, for any species, and are
    /// read before a genome-wide file beside them. A directory with the
    /// species' one genome-wide file alone gives that file under every
    /// sequence it holds records on, by name, or under `--chromosomes` alone.
    /// A directory with neither, and a directory that does not exist, are
    /// each refused naming both names looked for.
    #[test]
    fn a_local_directory_is_read_per_chromosome_or_from_the_one_genome_wide_file() {
        const RECORDS: &str = "2\t100\trs1\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
                               19\t400\trs3\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n";
        let dir = tempfile::tempdir().unwrap();
        let empty = scan_local_vcfs(dir.path(), "mus_musculus", None)
            .unwrap_err()
            .to_string();
        assert_eq!(
            empty,
            format!(
                "No variation VCF of mus_musculus in {}: looked for *-chr*.vcf.gz, one file a \
                 chromosome, and mus_musculus.vcf.gz, the one genome-wide file Ensembl publishes \
                 for a species other than homo_sapiens",
                dir.path().display()
            )
        );

        let missing = dir.path().join("absent");
        let no_directory = local_variation_inputs(&missing, "mus_musculus", None)
            .err()
            .unwrap()
            .to_string();
        assert_eq!(
            no_directory,
            format!(
                "Local variation VCF directory not found: {}: {}",
                missing.display(),
                looked_for("mus_musculus")
            )
        );

        let genome = write_vcf(dir.path(), "mus_musculus.vcf.gz", RECORDS);
        assert_eq!(
            scan_local_vcfs(dir.path(), "mus_musculus", None).unwrap(),
            [
                ("19".to_string(), GermlineVcf::Genome(genome.clone())),
                ("2".to_string(), GermlineVcf::Genome(genome.clone())),
            ]
        );
        assert_eq!(
            scan_local_vcfs(dir.path(), "mus_musculus", Some(&["19".to_string()])).unwrap(),
            [("19".to_string(), GermlineVcf::Genome(genome))]
        );
        assert!(scan_local_vcfs(dir.path(), "rattus_norvegicus", None).is_err());

        let own = write_vcf(dir.path(), "mus_musculus-chr2.vcf.gz", RECORDS);
        assert_eq!(
            scan_local_vcfs(dir.path(), "mus_musculus", None).unwrap(),
            [("2".to_string(), GermlineVcf::Chromosome(own))]
        );
    }

    /// A build of a species other than human that names no variation
    /// directory fetches the species' one genome-wide germline file, with its
    /// index, into its work directory. A site without the file is refused
    /// naming the address looked for, and no work directory is created.
    #[test]
    fn the_one_genome_wide_file_is_fetched_into_the_work_directory() {
        const RECORDS: &str = "19\t400\trs3\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n";
        let root = tempfile::tempdir().unwrap();
        let site_dir = root.path().join("site");
        let vcf_dir = site_dir.join(variation_vcf_dir(116, "GRCm39", "mus_musculus"));
        std::fs::create_dir_all(&vcf_dir).unwrap();
        let site = Site::Dir(site_dir);
        let work = WorkDir::of_release(&root.path().join("caches"), "mus_musculus", 116, "GRCm39");
        let published = Published::new(vec![site.clone()], "mus_musculus", "GRCm39", 116);
        let from = Downloads {
            published: &published,
            site: &site,
            work: &work,
        };
        let chromosomes = ["19".to_string()];
        let resolve = || {
            resolve_variation_inputs(
                &None,
                116,
                "GRCm39",
                "mus_musculus",
                Some(&chromosomes),
                &from,
            )
        };

        let error = resolve().err().unwrap().to_string();
        assert_eq!(
            error,
            format!(
                "Ensembl publishes no variation VCF for mus_musculus on GRCm39 at release 116: \
                 looked for {}, the one genome-wide file of a species other than homo_sapiens",
                vcf_dir.join("mus_musculus.vcf.gz").display()
            )
        );
        assert!(!work.path().exists());

        write_vcf(&vcf_dir, "mus_musculus.vcf.gz", RECORDS);
        std::fs::write(vcf_dir.join("mus_musculus.vcf.gz.csi"), b"index").unwrap();
        let variations = resolve().unwrap();
        let dest = work
            .path()
            .join("vep_cache_builder_mus_musculus_GRCm39_116_mus_musculus.vcf.gz");
        assert_eq!(
            variations.per_chromosome,
            [("19".to_string(), GermlineVcf::Genome(dest.clone()))]
        );
        assert!(dest.is_file());
        assert!(dest.with_extension("gz.csi").is_file());
    }
}
