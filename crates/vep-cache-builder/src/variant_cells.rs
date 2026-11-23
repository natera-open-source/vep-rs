// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `pubmed`, `var_synonyms` and `failed` cells of the known variants,
//! which no variation VCF carries: read from the table dumps of the
//! variation database and from the GVF file a release publishes of its
//! failed variants.
//!
//! Ensembl VEP's dumper writes the three cells from the database
//! (`Pipeline/DumpVEP/InitDump.pm:194-201, 257-269`,
//! `AnnotationSource/Database/Variation.pm:167-179`):
//!
//! - `pubmed`: the PubMed identifiers of the publications that cite the
//!   variation, from `variation_citation` and `publication`, joined by `,`.
//!   The query states no order; the identifiers stand here in ascending
//!   order of the publication's identifier.
//! - `var_synonyms`: the variation's names in `variation_synonym`, under
//!   every source whose name does not hold `dbSNP` in any letter case, as
//!   `<source>::<name>,<name>` groups joined by `--`. The groups stand in the
//!   order of the source names compared without regard to the case of an
//!   ASCII letter. The published table does not record the order the dump's
//!   query returned the names of one source in: they stand here in the
//!   reverse of the table's row order, which is the order Ensembl VEP's
//!   cache of release 116 holds for most GRCh38 groups.
//! - `failed`: 1 for a variation `failed_variation` holds a row of.
//!
//! The dumper keys the first by the variation's name and the other two by
//! its identifier; a variation VCF carries the name alone, so each table is
//! joined to `variation`, whose names are unique, and the cells are looked
//! up by name.
//!
//! The dump ends each list at 1,024 bytes
//! ([`vep_core::variant::STORED_LIST_CAP`]), inside an item as readily as
//! after one, and `vep` prints such a value without its cut item. No cut
//! piece is written here: a list of more than 1,024 bytes is stored as
//! `vep` prints its first 1,024 bytes
//! ([`vep_core::variant::written_pubmed`],
//! [`vep_core::variant::written_synonyms`]), and a shorter one as it stands.
//! The dump cuts the names of one source at that length before it puts the
//! source's name before them, and then the whole text; the first 1,024
//! bytes of the whole text are the same either way.
//!
//! ## Memory
//!
//! `variation` and `variation_synonym` hold a row for every variant and for
//! every name of one, so neither is held in memory. Each is read once as a
//! stream of bytes; the rows that take part in a cell are written to sorted
//! runs under the scratch directory and merged by the variation's
//! identifier, and no table is taken to be sorted. The cells are stored in
//! one file under the scratch directory, with an index in memory of 16
//! bytes for each variant that has a cell.
//!
//! ## Failed variants the VCFs lack
//!
//! A variation VCF holds the variants kept on display, and most failed
//! variants are not. The failed-variants GVF holds every feature of every
//! failed variant: a feature whose name a record of its chromosome's VCFs
//! carries flags that record, and any other becomes a record of its own
//! (`crate::variation::chromosome_records`).

use std::collections::HashMap;
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, PoisonError};

use anyhow::{bail, Context, Result};
use tracing::info;
use vep_core::variant::{written_pubmed, written_synonyms, STORED_LIST_CAP};

use crate::count::counted;
use crate::dump_rows::{read_rows, unsigned, DumpRow};
use crate::gvf::{FailedVariants, GvfFeature};
use crate::sorted_runs::{RunWriter, Runs};
use crate::variation::VariationRecord;

/// The tables of a variation database the `pubmed` cell is read from.
pub const PUBMED_TABLES: [&str; 3] = ["variation", "variation_citation", "publication"];

/// The tables of a variation database the `var_synonyms` cell is read from.
pub const SYNONYM_TABLES: [&str; 3] = ["variation", "variation_synonym", "source"];

/// The tables of a variation database the `failed` cell is read from, beside
/// the failed-variants GVF.
pub const FAILED_TABLES: [&str; 2] = ["variation", "failed_variation"];

/// The bytes of rows buffered before a sorted run is written.
const RUN_BUFFER_BYTES: usize = 64 << 20;

/// The files the three cells are read from, each `None` where the release
/// publishes none. A cell is read where every file it needs is named:
/// `pubmed` from `variation`, `variation_citation` and `publication`;
/// `var_synonyms` from `variation`, `variation_synonym` and `source`;
/// `failed` from the failed-variants GVF and from `variation` with
/// `failed_variation`.
#[derive(Debug, Clone, Default)]
pub struct VariantCellInputs {
    /// The dump of `variation`: the name of each variation.
    pub variation: Option<PathBuf>,
    /// The dump of `variation_citation`: the publications citing a variation.
    pub variation_citation: Option<PathBuf>,
    /// The dump of `publication`: the PubMed identifier of a publication.
    pub publication: Option<PathBuf>,
    /// The dump of `variation_synonym`: the other names of a variation.
    pub variation_synonym: Option<PathBuf>,
    /// The dump of `source`: the name of each source of a synonym.
    pub source: Option<PathBuf>,
    /// The dump of `failed_variation`: the variations marked as failed.
    pub failed_variation: Option<PathBuf>,
    /// The failed-variants GVF, `<species>_failed.gvf.gz`.
    pub failed_gvf: Option<PathBuf>,
    /// The directory the sorted runs and the cell store are written under,
    /// created when it does not exist; the files are removed when the cells
    /// are dropped.
    pub scratch: PathBuf,
}

/// What a build read for the three cells and what it wrote of them.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VariantCellFacts {
    /// `variation`, `variation_citation` and `publication` were read.
    pub pubmed_tables_read: bool,
    /// `variation`, `variation_synonym` and `source` were read.
    pub synonym_tables_read: bool,
    /// The failed-variants GVF was read.
    pub failed_gvf_read: bool,
    /// `failed_variation` was joined to `variation`.
    pub failed_variation_joined: bool,
    /// The records written with a `pubmed` cell.
    pub records_with_pubmed: u64,
    /// The records written with a `var_synonyms` cell.
    pub records_with_synonyms: u64,
    /// The records written with `failed` 1, those added from the GVF among
    /// them.
    pub records_flagged_failed: u64,
    /// The records added from the failed-variants GVF: features whose name no
    /// record of their chromosome's VCFs carries.
    pub records_added_from_gvf: u64,
    /// The records whose `pubmed` list is longer than the stored length and
    /// was written as its whole items within it.
    pub pubmed_lists_at_whole_items: u64,
    /// The records whose `var_synonyms` list is longer than the stored
    /// length and was written as its whole items within it.
    pub synonym_lists_at_whole_items: u64,
}

/// The three cells of a build's known variants, looked up by name.
#[derive(Debug)]
pub struct VariantCells {
    store: Option<Store>,
    failed: Option<FailedVariants>,
    read: VariantCellFacts,
    with_pubmed: AtomicU64,
    with_synonyms: AtomicU64,
    flagged_failed: AtomicU64,
    added_from_gvf: AtomicU64,
    pubmed_at_whole_items: AtomicU64,
    synonyms_at_whole_items: AtomicU64,
}

impl VariantCells {
    /// Reads the files of `inputs`: joins the tables into the cell store and
    /// reads the failed-variants GVF, keeping its features on the sequences
    /// `chromosomes` names.
    pub fn read(inputs: &VariantCellInputs, chromosomes: &[String]) -> Result<Self> {
        Self::read_with(inputs, chromosomes, RUN_BUFFER_BYTES, name_hash)
    }

    fn read_with(
        inputs: &VariantCellInputs,
        chromosomes: &[String],
        run_buffer: usize,
        hash: fn(&[u8]) -> u64,
    ) -> Result<Self> {
        let tables = Tables::of(inputs);
        let store = match &tables {
            Some(tables) => match Joined::read(tables, &inputs.scratch, run_buffer)? {
                Some(joined) => Some(joined.store(&inputs.scratch, hash)?),
                None => None,
            },
            None => None,
        };
        let failed = match &inputs.failed_gvf {
            Some(path) => {
                info!("Reading the failed variants: {}", path.display());
                let failed = FailedVariants::load(path, chromosomes)
                    .context("Failed to read the failed-variants GVF")?;
                info!(
                    "Failed variants read: {}, {} on a sequence no variation file is written for",
                    counted(failed.features, "feature"),
                    failed.elsewhere
                );
                Some(failed)
            }
            None => None,
        };
        Ok(Self {
            store,
            failed,
            read: VariantCellFacts {
                pubmed_tables_read: tables.as_ref().is_some_and(|t| t.pubmed.is_some()),
                synonym_tables_read: tables.as_ref().is_some_and(|t| t.synonyms.is_some()),
                failed_gvf_read: inputs.failed_gvf.is_some(),
                failed_variation_joined: tables.as_ref().is_some_and(|t| t.failed.is_some()),
                ..Default::default()
            },
            with_pubmed: AtomicU64::new(0),
            with_synonyms: AtomicU64::new(0),
            flagged_failed: AtomicU64::new(0),
            added_from_gvf: AtomicU64::new(0),
            pubmed_at_whole_items: AtomicU64::new(0),
            synonyms_at_whole_items: AtomicU64::new(0),
        })
    }

    /// What was read, and what the records given their cells so far were
    /// written with.
    pub fn facts(&self) -> VariantCellFacts {
        let count = |counter: &AtomicU64| counter.load(Ordering::Relaxed);
        VariantCellFacts {
            records_with_pubmed: count(&self.with_pubmed),
            records_with_synonyms: count(&self.with_synonyms),
            records_flagged_failed: count(&self.flagged_failed),
            records_added_from_gvf: count(&self.added_from_gvf),
            pubmed_lists_at_whole_items: count(&self.pubmed_at_whole_items),
            synonym_lists_at_whole_items: count(&self.synonyms_at_whole_items),
            ..self.read.clone()
        }
    }

    /// Fills the three cells of `record` by its name, as
    /// `Dumper/Variation.pm:239, 252-253` prints them: `failed` 1 when the
    /// failed-variants GVF or the joined `failed_variation` names the
    /// variant, and its `pubmed` and `var_synonyms` texts. A record without
    /// a name takes none.
    pub(crate) fn apply(&self, record: &mut VariationRecord) -> Result<()> {
        let name = record.variation_name.as_str();
        if name.is_empty() {
            return Ok(());
        }
        let bump = |counter: &AtomicU64| {
            counter.fetch_add(1, Ordering::Relaxed);
        };
        let mut failed = record.failed != 0
            || self
                .failed
                .as_ref()
                .is_some_and(|failed| failed.names(name));
        if let Some(stored) = self.store.as_ref().map(|store| store.lookup(name)) {
            if let Some(cells) = stored? {
                failed |= cells.failed;
                if let Some(pubmed) = cells.pubmed {
                    bump(&self.with_pubmed);
                    if cells.pubmed_at_whole_items {
                        bump(&self.pubmed_at_whole_items);
                    }
                    record.pubmed = Some(pubmed);
                }
                if let Some(synonyms) = cells.synonyms {
                    bump(&self.with_synonyms);
                    if cells.synonyms_at_whole_items {
                        bump(&self.synonyms_at_whole_items);
                    }
                    record.var_synonyms = Some(synonyms);
                }
            }
        }
        if failed {
            record.failed = 1;
            bump(&self.flagged_failed);
        }
        Ok(())
    }

    /// The features of the failed-variants GVF on `chromosome`, by start
    /// then name.
    pub(crate) fn failed_features(&self, chromosome: &str) -> &[GvfFeature] {
        self.failed
            .as_ref()
            .map_or(&[], |failed| failed.on(chromosome))
    }

    /// Counts `records` more written from a feature of the failed-variants
    /// GVF.
    pub(crate) fn count_added(&self, records: u64) {
        self.added_from_gvf.fetch_add(records, Ordering::Relaxed);
    }
}

/// The tables a join reads: `variation` and, for each cell, the tables read
/// beside it.
struct Tables<'a> {
    variation: &'a Path,
    /// `variation_citation` and `publication`.
    pubmed: Option<(&'a Path, &'a Path)>,
    /// `variation_synonym` and `source`.
    synonyms: Option<(&'a Path, &'a Path)>,
    /// `failed_variation`.
    failed: Option<&'a Path>,
}

impl<'a> Tables<'a> {
    /// The tables `inputs` names for at least one cell; `None` without
    /// `variation` or without every table of any cell.
    fn of(inputs: &'a VariantCellInputs) -> Option<Self> {
        let path = |table: &'a Option<PathBuf>| table.as_deref();
        let tables = Self {
            variation: path(&inputs.variation)?,
            pubmed: path(&inputs.variation_citation).zip(path(&inputs.publication)),
            synonyms: path(&inputs.variation_synonym).zip(path(&inputs.source)),
            failed: path(&inputs.failed_variation),
        };
        (tables.pubmed.is_some() || tables.synonyms.is_some() || tables.failed.is_some())
            .then_some(tables)
    }
}

/// The identifier in the `nth` kept column of `row`, or an error naming the
/// file, the row and the column.
fn identifier(row: &DumpRow<'_>, nth: usize, path: &Path, column: &str) -> Result<u32> {
    row.value(nth).and_then(unsigned).with_context(|| {
        format!(
            "{} row {} holds no identifier in {column}: {}",
            path.display(),
            row.number,
            match row.value(nth) {
                Some(value) => format!("{:?}", String::from_utf8_lossy(value)),
                None => "NULL".to_string(),
            }
        )
    })
}

/// The sources whose synonyms a cell lists: every source of the `source`
/// dump at `path` whose name does not hold `dbSNP`, the dump's query
/// comparing under a collation that tells no letter case
/// (`InitDump.pm:262-265`: `name not like '%dbSNP%'`). Each is ranked by its
/// name, compared without regard to the case of an ASCII letter, the order
/// the dump's `group by` leaves the sources of one variation in.
struct KeptSources {
    /// The rank of each kept source, by its identifier.
    rank: HashMap<u32, u16>,
    /// The names of the kept sources, by rank.
    names: Vec<String>,
}

impl KeptSources {
    fn load(path: &Path) -> Result<Self> {
        let mut kept: Vec<(Vec<u8>, String, u32)> = Vec::new();
        read_rows(path, 8, &[0, 1], |row| {
            let id = identifier(row, 0, path, "source_id")?;
            let name = String::from_utf8_lossy(row.value(1).unwrap_or_default()).into_owned();
            let folded = name.to_ascii_uppercase();
            if !folded.contains("DBSNP") {
                kept.push((folded.into_bytes(), name, id));
            }
            Ok(())
        })?;
        kept.sort();
        if kept.len() > usize::from(u16::MAX) {
            bail!("{} names {} sources", path.display(), kept.len());
        }
        Ok(Self {
            rank: kept
                .iter()
                .enumerate()
                .map(|(rank, (_, _, id))| (*id, rank as u16))
                .collect(),
            names: kept.into_iter().map(|(_, name, _)| name).collect(),
        })
    }
}

/// The identifiers a cell is read for, sorted. A table mostly in the order
/// of its identifiers is tested against the neighbours of the last test;
/// any other row falls back to a search.
struct Wanted {
    ids: Vec<u32>,
    cursor: usize,
}

impl Wanted {
    fn new(mut ids: Vec<u32>) -> Self {
        ids.sort_unstable();
        ids.dedup();
        Self { ids, cursor: 0 }
    }

    fn contains(&mut self, id: u32) -> bool {
        let ids = &self.ids;
        // The first identifier not below `id` lies at `at` when the one
        // before `at` is below `id` and the one at `at` is not.
        let first_not_below =
            |at: usize| (at == 0 || ids[at - 1] < id) && ids.get(at).is_none_or(|next| *next >= id);
        let at = if first_not_below(self.cursor) {
            self.cursor
        } else if self.cursor < ids.len() && first_not_below(self.cursor + 1) {
            self.cursor + 1
        } else {
            ids.partition_point(|held| *held < id)
        };
        self.cursor = at;
        ids.get(at) == Some(&id)
    }
}

/// What the tables hold for the cells, each part in the order of the
/// variations' identifiers: the small tables in memory, the rows kept of the
/// two large ones as sorted runs.
struct Joined {
    /// (variation, publication, PubMed identifier) of each citation of a
    /// publication that has an identifier, by variation then publication.
    citations: Vec<(u32, u32, u32)>,
    /// The failed variations.
    failed: Vec<u32>,
    sources: KeptSources,
    /// The synonyms of the kept sources, by variation, then source, then
    /// against the table's row order.
    synonyms: Runs,
    /// The names of the variations a cell is read for.
    names: Runs,
}

impl Joined {
    /// Reads `tables`, spilling the kept rows of `variation_synonym` and of
    /// `variation` under `scratch` in runs of `run_buffer` bytes: `None` when
    /// no variation takes part in a cell.
    fn read(tables: &Tables<'_>, scratch: &Path, run_buffer: usize) -> Result<Option<Self>> {
        std::fs::create_dir_all(scratch)
            .with_context(|| format!("Failed to create {}", scratch.display()))?;

        let mut citations: Vec<(u32, u32, u32)> = Vec::new();
        if let Some((variation_citation, publication)) = tables.pubmed {
            info!("Reading the PubMed identifiers: {}", publication.display());
            let mut pmids: HashMap<u32, u32> = HashMap::new();
            read_rows(publication, 8, &[0, 3], |row| {
                let id = identifier(row, 0, publication, "publication_id")?;
                if row.value(1).is_some() {
                    pmids.insert(id, identifier(row, 1, publication, "pmid")?);
                }
                Ok(())
            })?;
            let rows = read_rows(variation_citation, 3, &[0, 1], |row| {
                let variation = identifier(row, 0, variation_citation, "variation_id")?;
                let publication = identifier(row, 1, variation_citation, "publication_id")?;
                if let Some(pmid) = pmids.get(&publication) {
                    citations.push((variation, publication, *pmid));
                }
                Ok(())
            })?;
            citations.sort_unstable();
            info!(
                "PubMed identifiers read: {} with one, {} of {} to one",
                counted(pmids.len(), "publication"),
                citations.len(),
                counted(rows, "citation")
            );
        }

        let mut failed: Vec<u32> = Vec::new();
        if let Some(failed_variation) = tables.failed {
            read_rows(failed_variation, 3, &[1], |row| {
                failed.push(identifier(row, 0, failed_variation, "variation_id")?);
                Ok(())
            })?;
            failed.sort_unstable();
            failed.dedup();
            info!(
                "Failed variations read: {} from {}",
                failed.len(),
                failed_variation.display()
            );
        }

        let mut sources = KeptSources {
            rank: HashMap::new(),
            names: Vec::new(),
        };
        let mut synonym_runs = Runs::default();
        let mut with_synonyms: Vec<u32> = Vec::new();
        if let Some((variation_synonym, source)) = tables.synonyms {
            sources = KeptSources::load(source)?;
            info!(
                "Reading the variant synonyms of {}: {}",
                counted(sources.names.len(), "source"),
                variation_synonym.display()
            );
            let mut runs = RunWriter::new(scratch, "synonyms", run_buffer);
            let rows = read_rows(variation_synonym, 5, &[1, 3, 4], |row| {
                let source = identifier(row, 1, variation_synonym, "source_id")?;
                let (Some(rank), Some(name)) = (sources.rank.get(&source), row.value(2)) else {
                    return Ok(());
                };
                let variation = identifier(row, 0, variation_synonym, "variation_id")?;
                if with_synonyms.last() != Some(&variation) {
                    with_synonyms.push(variation);
                }
                // By variation, then source, then against the row order.
                let key = ((u64::from(variation) << 16) | u64::from(*rank), !row.number);
                runs.push(key, name)
            })?;
            synonym_runs = runs.finish()?;
            info!(
                "Variant synonyms read: {} of {}, sorted in {}",
                synonym_runs.records,
                counted(rows, "row"),
                spilled(&synonym_runs)
            );
        }

        with_synonyms.extend(citations.iter().map(|(variation, _, _)| *variation));
        with_synonyms.extend(&failed);
        let mut wanted = Wanted::new(with_synonyms);
        if wanted.ids.is_empty() {
            info!("No variation of the tables takes part in a cell");
            return Ok(None);
        }

        let variation = tables.variation;
        info!(
            "Reading the names of {}: {}",
            counted(wanted.ids.len(), "variation"),
            variation.display()
        );
        let mut runs = RunWriter::new(scratch, "names", run_buffer);
        let rows = read_rows(variation, 12, &[0, 2], |row| {
            let id = identifier(row, 0, variation, "variation_id")?;
            match row.value(1) {
                Some(name) if wanted.contains(id) => runs.push((u64::from(id), 0), name),
                _ => Ok(()),
            }
        })?;
        let name_runs = runs.finish()?;
        info!(
            "Names read: {} of {}, sorted in {}",
            name_runs.records,
            counted(rows, "row"),
            spilled(&name_runs)
        );
        Ok(Some(Self {
            citations,
            failed,
            sources,
            synonyms: synonym_runs,
            names: name_runs,
        }))
    }

    /// Merges the parts by the variation's identifier into the cell store
    /// under `scratch`, its names indexed under `hash`: one record for each
    /// variation that has a name and a cell.
    fn store(self, scratch: &Path, hash: fn(&[u8]) -> u64) -> Result<Store> {
        let Self {
            citations,
            failed,
            sources,
            synonyms,
            names,
        } = self;
        let (name_runs, synonym_runs) = (names, synonyms);
        let mut writer =
            StoreWriter::create(&scratch.join("cells.store"), hash, name_runs.records)?;
        let mut names = name_runs.merge()?;
        let mut synonyms = synonym_runs.merge()?;
        let (mut citation, mut failure) = (0usize, 0usize);
        let mut name: Vec<u8> = Vec::new();
        let mut pubmed = String::new();
        let mut synonym_text = String::new();
        let mut last = None;
        while let Some(((id, _), read)) = names.head() {
            // The first row of an identifier stands.
            let repeated = last == Some(id);
            name.clear();
            name.extend_from_slice(read);
            names.advance()?;
            if repeated {
                continue;
            }
            last = Some(id);
            let id = id as u32;

            pubmed.clear();
            while let Some((variation, _, pmid)) = citations.get(citation) {
                if *variation > id {
                    break;
                }
                // Past the stored length nothing more of the list is read.
                if *variation == id && pubmed.len() <= STORED_LIST_CAP {
                    if !pubmed.is_empty() {
                        pubmed.push(',');
                    }
                    pubmed.push_str(&pmid.to_string());
                }
                citation += 1;
            }

            synonym_text.clear();
            let mut group = None;
            while let Some(((key, _), synonym)) = synonyms.head() {
                let (variation, rank) = ((key >> 16) as u32, (key & 0xffff) as usize);
                if variation > id {
                    break;
                }
                if variation == id && synonym_text.len() <= STORED_LIST_CAP {
                    if group == Some(rank) {
                        synonym_text.push(',');
                    } else {
                        if group.is_some() {
                            synonym_text.push_str("--");
                        }
                        synonym_text.push_str(&sources.names[rank]);
                        synonym_text.push_str("::");
                        group = Some(rank);
                    }
                    synonym_text.push_str(&String::from_utf8_lossy(synonym));
                }
                synonyms.advance()?;
            }

            while failed.get(failure).is_some_and(|variation| *variation < id) {
                failure += 1;
            }
            let cells = StoredCells {
                failed: failed.get(failure) == Some(&id),
                pubmed: stored_list(&pubmed, written_pubmed),
                synonyms: stored_list(&synonym_text, written_synonyms),
            };
            if cells.failed || cells.pubmed.is_some() || cells.synonyms.is_some() {
                writer.add(&name, &cells)?;
            }
        }
        let store = writer.finish()?;
        info!(
            "Cells stored for {}: {} in {}, {} of index in memory",
            counted(store.index.len(), "variation"),
            counted(store.bytes, "byte"),
            store.path.display(),
            counted(store.index_bytes(), "byte")
        );
        Ok(store)
    }
}

/// The run files of `runs`, for a log line.
fn spilled(runs: &Runs) -> String {
    format!(
        "{} of {} on disk",
        counted(runs.count(), "run"),
        counted(runs.bytes, "byte")
    )
}

/// A list as it is stored, and whether it stands at its whole items within
/// the stored length; `None` for a text Perl reads as false, empty or `0`,
/// which the dumper prints as an empty cell (`Dumper/Variation.pm:252-253`).
///
/// A text of at most [`STORED_LIST_CAP`] bytes is stored as it stands. A
/// longer one is stored as `written` reads its first [`STORED_LIST_CAP`]
/// bytes: without the item the length cuts. Where the length falls inside a
/// character, the item that character belongs to is the cut one. A text
/// ends before a tab or a line end inside it, where the dumper's read of its
/// own file of the lists ends it (`:482-484`).
fn stored_list(text: &str, written: fn(Option<&str>) -> Option<&str>) -> Option<(&str, bool)> {
    let mut cut = text.len().min(STORED_LIST_CAP);
    while !text.is_char_boundary(cut) {
        cut -= 1;
    }
    let head = &text[..cut];
    let (stored, bounded) = match head.find(['\t', '\n']) {
        Some(end) => (&head[..end], false),
        None if text.len() <= STORED_LIST_CAP => (text, false),
        None if cut == STORED_LIST_CAP => (written(Some(head)).unwrap_or(""), true),
        None => {
            // `written` tells a cut list by its length: the bytes of the cut
            // character are made up with bytes of no separator.
            let padded = format!("{head}{}", "x".repeat(STORED_LIST_CAP - cut));
            let whole = written(Some(&padded)).map_or(0, str::len);
            (&head[..whole.min(cut)], true)
        }
    };
    (!stored.is_empty() && stored != "0").then_some((stored, bounded))
}

/// The cells of one variation as the store holds them: each list with
/// whether it stands at its whole items within the stored length.
struct StoredCells<'a> {
    failed: bool,
    pubmed: Option<(&'a str, bool)>,
    synonyms: Option<(&'a str, bool)>,
}

/// The cells a lookup returns.
struct Cells {
    failed: bool,
    pubmed: Option<String>,
    synonyms: Option<String>,
    pubmed_at_whole_items: bool,
    synonyms_at_whole_items: bool,
}

/// The bits of a stored record's flags.
const FAILED: u8 = 1;
const PUBMED_AT_WHOLE_ITEMS: u8 = 2;
const SYNONYMS_AT_WHOLE_ITEMS: u8 = 4;

/// The bytes before a stored record's texts: the lengths of its name and of
/// its two lists, each 16 bits, and its flags.
const RECORD_HEADER_BYTES: usize = 7;

/// Writes the cell store: one record a variation, and the index of where
/// each lies.
struct StoreWriter {
    path: PathBuf,
    writer: BufWriter<File>,
    hash: fn(&[u8]) -> u64,
    /// (hash of the name, place of the record) of each record written.
    index: Vec<(u64, u64)>,
    bytes: u64,
}

impl StoreWriter {
    /// A store at `path` of at most `records` records, whose names are
    /// indexed under `hash`.
    fn create(path: &Path, hash: fn(&[u8]) -> u64, records: u64) -> Result<Self> {
        let file =
            File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
        Ok(Self {
            path: path.to_path_buf(),
            writer: BufWriter::with_capacity(1 << 20, file),
            hash,
            // Sized once: an index that grew by doubling would hold up to
            // twice its entries.
            index: Vec::with_capacity(usize::try_from(records).context("record count")?),
            bytes: 0,
        })
    }

    fn add(&mut self, name: &[u8], cells: &StoredCells<'_>) -> Result<()> {
        let pubmed = cells.pubmed.map_or("", |(text, _)| text);
        let synonyms = cells.synonyms.map_or("", |(text, _)| text);
        let length = |bytes: usize| {
            u16::try_from(bytes).with_context(|| {
                format!(
                    "a value of {bytes} bytes for the variation {}",
                    String::from_utf8_lossy(name)
                )
            })
        };
        let at_whole_items = |list: Option<(&str, bool)>| list.is_some_and(|(_, bounded)| bounded);
        let mut flags = 0;
        for (flag, set) in [
            (FAILED, cells.failed),
            (PUBMED_AT_WHOLE_ITEMS, at_whole_items(cells.pubmed)),
            (SYNONYMS_AT_WHOLE_ITEMS, at_whole_items(cells.synonyms)),
        ] {
            if set {
                flags |= flag;
            }
        }
        let mut header = [0u8; RECORD_HEADER_BYTES];
        header[0..2].copy_from_slice(&length(name.len())?.to_le_bytes());
        header[2..4].copy_from_slice(&length(pubmed.len())?.to_le_bytes());
        header[4..6].copy_from_slice(&length(synonyms.len())?.to_le_bytes());
        header[6] = flags;
        self.writer
            .write_all(&header)
            .and_then(|()| self.writer.write_all(name))
            .and_then(|()| self.writer.write_all(pubmed.as_bytes()))
            .and_then(|()| self.writer.write_all(synonyms.as_bytes()))
            .with_context(|| format!("Failed to write {}", self.path.display()))?;
        self.index.push(((self.hash)(name), self.bytes));
        self.bytes += (RECORD_HEADER_BYTES + name.len() + pubmed.len() + synonyms.len()) as u64;
        Ok(())
    }

    /// Ends the file and sorts the index over it.
    fn finish(mut self) -> Result<Store> {
        self.writer
            .flush()
            .with_context(|| format!("Failed to write {}", self.path.display()))?;
        drop(self.writer);
        let file = File::open(&self.path)
            .with_context(|| format!("Failed to open {}", self.path.display()))?;
        let mut index = self.index;
        index.sort_unstable();
        if index.len() > u32::MAX as usize {
            bail!("{} variations have a cell", index.len());
        }

        // About eight records a bucket, by the leading bits of the hash.
        let bits = (index.len() / 8)
            .next_power_of_two()
            .trailing_zeros()
            .min(24);
        let mut buckets = vec![0u32; (1usize << bits) + 1];
        for (hash, _) in &index {
            buckets[bucket_of(*hash, bits) + 1] += 1;
        }
        for bucket in 1..buckets.len() {
            buckets[bucket] += buckets[bucket - 1];
        }
        Ok(Store {
            path: self.path,
            file: Mutex::new(file),
            index,
            buckets,
            bits,
            hash: self.hash,
            bytes: self.bytes,
        })
    }
}

/// The bucket of `hash` among `2^bits`: its leading bits.
fn bucket_of(hash: u64, bits: u32) -> usize {
    if bits == 0 {
        0
    } else {
        (hash >> (64 - bits)) as usize
    }
}

/// The hash a name is indexed under.
fn name_hash(name: &[u8]) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    hasher.write(name);
    hasher.finish()
}

/// The cells of every variation that has one, in a file, with the index of
/// the file in memory: 16 bytes a variation, the hash of its name and the
/// place of its record, and a table of the index's buckets. A lookup reads
/// the file only for a name whose hash the index holds, and compares the
/// stored name, so two names of one hash are told apart. The file is removed
/// when the store is dropped.
#[derive(Debug)]
struct Store {
    path: PathBuf,
    file: Mutex<File>,
    /// (hash of the name, place of the record), sorted.
    index: Vec<(u64, u64)>,
    /// Where each bucket's entries begin in `index`, and the end of the last.
    buckets: Vec<u32>,
    bits: u32,
    hash: fn(&[u8]) -> u64,
    /// The bytes of the file.
    bytes: u64,
}

impl Store {
    /// The bytes the index takes in memory.
    fn index_bytes(&self) -> usize {
        self.index.len() * std::mem::size_of::<(u64, u64)>()
            + self.buckets.len() * std::mem::size_of::<u32>()
    }

    /// The cells stored under `name`.
    fn lookup(&self, name: &str) -> Result<Option<Cells>> {
        let hash = (self.hash)(name.as_bytes());
        let bucket = bucket_of(hash, self.bits);
        let entries = &self.index[self.buckets[bucket] as usize..self.buckets[bucket + 1] as usize];
        let first = entries.partition_point(|(held, _)| *held < hash);
        for (_, place) in entries[first..]
            .iter()
            .take_while(|(held, _)| *held == hash)
        {
            if let Some(cells) = self.read(*place, name)? {
                return Ok(Some(cells));
            }
        }
        Ok(None)
    }

    /// The cells of the record at `place`, when its name is `name`.
    fn read(&self, place: u64, name: &str) -> Result<Option<Cells>> {
        let failed = || format!("Failed to read {}", self.path.display());
        let mut file = self.file.lock().unwrap_or_else(PoisonError::into_inner);
        file.seek(SeekFrom::Start(place)).with_context(failed)?;
        let mut header = [0u8; RECORD_HEADER_BYTES];
        file.read_exact(&mut header).with_context(failed)?;
        let length = |at: usize| usize::from(u16::from_le_bytes([header[at], header[at + 1]]));
        let (name_bytes, pubmed_bytes, synonym_bytes) = (length(0), length(2), length(4));
        if name_bytes != name.len() {
            return Ok(None);
        }
        let mut texts = vec![0u8; name_bytes + pubmed_bytes + synonym_bytes];
        file.read_exact(&mut texts).with_context(failed)?;
        drop(file);
        if &texts[..name_bytes] != name.as_bytes() {
            return Ok(None);
        }
        let text =
            |bytes: &[u8]| (!bytes.is_empty()).then(|| String::from_utf8_lossy(bytes).into_owned());
        let flags = header[6];
        Ok(Some(Cells {
            failed: flags & FAILED != 0,
            pubmed: text(&texts[name_bytes..name_bytes + pubmed_bytes]),
            synonyms: text(&texts[name_bytes + pubmed_bytes..]),
            pubmed_at_whole_items: flags & PUBMED_AT_WHOLE_ITEMS != 0,
            synonyms_at_whole_items: flags & SYNONYMS_AT_WHOLE_ITEMS != 0,
        }))
    }
}

impl Drop for Store {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A row of `variation`: the identifier and the name among the schema's
    /// twelve columns.
    fn variation_row(id: u32, name: &str) -> String {
        format!("{id}\t1\t{name}\t\\N\t2\t0\t\\N\t\\N\t\\N\t\\N\t\\N\t1\n")
    }

    /// A row of `publication`: the identifier, a title that holds a tab and a
    /// line end behind their escapes, and the PubMed identifier or NULL.
    fn publication_row(id: u32, pmid: Option<u32>) -> String {
        let pmid = pmid.map_or("\\N".to_string(), |pmid| pmid.to_string());
        format!("{id}\ta title\\\tof two\\\nlines\tAn Author\t{pmid}\t\\N\t2020\t\\N\t\\N\n")
    }

    fn source_row(id: u32, name: &str) -> String {
        format!("{id}\t{name}\t1\t\\N\t\\N\t\\N\tgermline\tvariation_synonym\n")
    }

    fn synonym_row(id: u32, variation: u32, source: u32, name: &str) -> String {
        format!("{id}\t{variation}\t\\N\t{source}\t{name}\n")
    }

    /// Made-up tables written as plain dumps into a directory, and the
    /// inputs that name them.
    struct Dumps {
        dir: tempfile::TempDir,
        inputs: VariantCellInputs,
    }

    impl Dumps {
        fn new() -> Self {
            let dir = tempfile::tempdir().unwrap();
            let inputs = VariantCellInputs {
                scratch: dir.path().join("scratch"),
                ..Default::default()
            };
            Self { dir, inputs }
        }

        /// Writes `rows` as the dump of `table` and names it in the inputs.
        fn table(mut self, table: &str, rows: &str) -> Self {
            let path = self.dir.path().join(format!("{table}.txt"));
            std::fs::write(&path, rows).unwrap();
            let named = match table {
                "variation" => &mut self.inputs.variation,
                "variation_citation" => &mut self.inputs.variation_citation,
                "publication" => &mut self.inputs.publication,
                "variation_synonym" => &mut self.inputs.variation_synonym,
                "source" => &mut self.inputs.source,
                "failed_variation" => &mut self.inputs.failed_variation,
                other => panic!("no table {other}"),
            };
            *named = Some(path);
            self
        }

        fn read(&self) -> VariantCells {
            VariantCells::read(&self.inputs, &[]).unwrap()
        }
    }

    /// The `failed`, `pubmed` and `var_synonyms` cells a record named `name`
    /// takes.
    fn cells_of(cells: &VariantCells, name: &str) -> (u8, Option<String>, Option<String>) {
        let mut record = VariationRecord {
            variation_name: name.to_string(),
            ..Default::default()
        };
        cells.apply(&mut record).unwrap();
        (record.failed, record.pubmed, record.var_synonyms)
    }

    fn text(value: &str) -> Option<String> {
        Some(value.to_string())
    }

    /// Tables in no order of their identifiers. Variation 10 is cited by
    /// three publications, one without a PubMed identifier, and has synonyms
    /// under two kept sources and one named for dbSNP; variation 30 has the
    /// same two sources with their rows the other way round; variation 20 is
    /// cited by the publication without an identifier alone and is failed;
    /// variation 40 is failed for two reasons; variation 50 has a synonym of
    /// a source named for dbSNP in another letter case alone; variation 60 is
    /// cited by a publication whose identifier is 0; variation 70 is in no
    /// other table; and each table but `variation` holds a row of a
    /// variation `variation` lacks.
    fn made_up_tables() -> Dumps {
        let variation: String = [
            (30, "mv30"),
            (10, "mv10"),
            (50, "mv50"),
            (20, "mv20"),
            (40, "mv40"),
            (60, "mv60"),
            (70, "mv70"),
        ]
        .iter()
        .map(|(id, name)| variation_row(*id, name))
        .collect();
        let publication = [
            publication_row(7, Some(300)),
            publication_row(3, Some(900)),
            publication_row(5, None),
            publication_row(9, Some(0)),
        ]
        .concat();
        let citation = "10\t7\t\\N\n30\t3\t\\N\n10\t5\t\\N\n10\t3\t\\N\n20\t5\t\\N\n60\t9\t\\N\n\
                        99\t3\t\\N\n";
        let source = [
            source_row(1, "dbSNP"),
            source_row(2, "Archive DBSNP"),
            source_row(3, "Zeta"),
            source_row(4, "alpha"),
        ]
        .concat();
        let synonym = [
            synonym_row(1, 10, 3, "z1"),
            synonym_row(2, 30, 4, "a30"),
            synonym_row(3, 10, 1, "old10"),
            synonym_row(4, 10, 4, "a1"),
            synonym_row(5, 30, 3, "z30"),
            synonym_row(6, 10, 4, "a2"),
            synonym_row(7, 50, 2, "old50"),
            synonym_row(8, 10, 3, "z2"),
            synonym_row(9, 10, 4, "a3"),
            synonym_row(10, 98, 3, "zz"),
        ]
        .concat();
        Dumps::new()
            .table("variation", &variation)
            .table("publication", &publication)
            .table("variation_citation", citation)
            .table("source", &source)
            .table("variation_synonym", &synonym)
            .table(
                "failed_variation",
                "1\t40\t3\n2\t20\t1\n3\t40\t7\n4\t97\t1\n",
            )
    }

    /// The cells of the made-up tables. The PubMed identifiers stand in
    /// ascending order of the publication's identifier, whatever the order
    /// of the rows, and a publication without one adds nothing. The sources
    /// named for dbSNP, in either letter case, are left out; the kept ones
    /// stand in the order of their names compared without regard to case
    /// (`alpha` before `Zeta`, which byte order puts the other way), whichever
    /// source's rows come first; the names of one source stand in the
    /// reverse of the table's row order. A failed variation is flagged once
    /// however many rows name it. A list that is `0` and a variation of no
    /// cell give nothing, as does a name no table holds.
    #[test]
    fn the_tables_join_into_the_cells_of_each_variant() {
        let dumps = made_up_tables();
        let cells = dumps.read();
        assert_eq!(
            cells_of(&cells, "mv10"),
            (0, text("900,300"), text("alpha::a3,a2,a1--Zeta::z2,z1"))
        );
        assert_eq!(
            cells_of(&cells, "mv30"),
            (0, text("900"), text("alpha::a30--Zeta::z30"))
        );
        assert_eq!(cells_of(&cells, "mv20"), (1, None, None));
        assert_eq!(cells_of(&cells, "mv40"), (1, None, None));
        for none in ["mv50", "mv60", "mv70", "mv99", ""] {
            assert_eq!(cells_of(&cells, none), (0, None, None), "{none}");
        }
        assert_eq!(
            cells.facts(),
            VariantCellFacts {
                pubmed_tables_read: true,
                synonym_tables_read: true,
                failed_gvf_read: false,
                failed_variation_joined: true,
                records_with_pubmed: 2,
                records_with_synonyms: 2,
                records_flagged_failed: 2,
                records_added_from_gvf: 0,
                pubmed_lists_at_whole_items: 0,
                synonym_lists_at_whole_items: 0,
            }
        );
        // The store is removed with the cells.
        let store = dumps.inputs.scratch.join("cells.store");
        assert!(store.is_file());
        drop(cells);
        assert_eq!(std::fs::read_dir(&dumps.inputs.scratch).unwrap().count(), 0);
    }

    /// A run buffer that holds two rows spills the eight kept synonyms as
    /// four runs and the five names read as three, and the cells merged from
    /// the runs are those of one run each.
    #[test]
    fn rows_spilled_as_several_runs_merge_into_the_same_cells() {
        let dumps = made_up_tables();
        let tables = Tables::of(&dumps.inputs).unwrap();
        let joined = Joined::read(&tables, &dumps.inputs.scratch, 64)
            .unwrap()
            .unwrap();
        assert_eq!((joined.synonyms.count(), joined.synonyms.records), (4, 8));
        assert_eq!((joined.names.count(), joined.names.records), (3, 5));
        drop(joined);

        let spilled = VariantCells::read_with(&dumps.inputs, &[], 64, name_hash).unwrap();
        let whole = dumps.read();
        for name in ["mv10", "mv20", "mv30", "mv40", "mv50", "mv60", "mv70"] {
            assert_eq!(cells_of(&spilled, name), cells_of(&whole, name), "{name}");
        }
        assert_eq!(
            cells_of(&spilled, "mv10").2,
            text("alpha::a3,a2,a1--Zeta::z2,z1")
        );
    }

    /// Names whose hashes are equal are told apart by the stored name: under
    /// a hash that is the same for every name, and under one that is the
    /// name's length, each name still finds its own cells and a name of the
    /// same length that no table holds finds none.
    #[test]
    fn names_of_one_hash_find_their_own_cells() {
        let dumps = made_up_tables();
        let hashes: [fn(&[u8]) -> u64; 2] = [|_| 42, |name| name.len() as u64];
        for hash in hashes {
            let cells = VariantCells::read_with(&dumps.inputs, &[], 1 << 20, hash).unwrap();
            assert_eq!(cells_of(&cells, "mv10").1, text("900,300"));
            assert_eq!(cells_of(&cells, "mv30").1, text("900"));
            assert_eq!(cells_of(&cells, "mv40"), (1, None, None));
            assert_eq!(cells_of(&cells, "mv11"), (0, None, None));
            assert_eq!(cells_of(&cells, "mv100"), (0, None, None));
        }
    }

    /// `count` seven-digit identifiers joined by `,`, the first one `first`.
    fn identifiers(first: &str, count: usize) -> String {
        std::iter::once(first.to_string())
            .chain((1..count).map(|n| format!("{}", 1_000_000 + n)))
            .collect::<Vec<String>>()
            .join(",")
    }

    /// A PubMed list of at most 1,024 bytes is stored as it stands, one of
    /// exactly that length too. A longer one is stored as its first 1,024
    /// bytes read without their cut item: the text before the last comma
    /// where the length falls inside an identifier or at the end of one,
    /// and the 1,024 bytes themselves where it falls after a comma. A list
    /// that is `0`, or empty, is no cell.
    #[test]
    fn a_pubmed_list_over_the_stored_length_is_stored_at_its_whole_items() {
        let pubmed = |text: &str| stored_list(text, written_pubmed).map(|(s, b)| (s.len(), b));
        let short = identifiers("1000000", 100);
        assert_eq!(pubmed(&short), Some((799, false)));

        // 128 identifiers: 127 of seven digits and one of eight.
        let exact = identifiers("10000000", 128);
        assert_eq!(exact.len(), 1024);
        assert_eq!(
            stored_list(&exact, written_pubmed),
            Some((exact.as_str(), false))
        );

        // The length falls two digits into the 129th identifier.
        let inside = identifiers("10000", 300);
        let whole = identifiers("10000", 128);
        assert_eq!(whole.len(), 1021);
        assert_eq!(
            stored_list(&inside, written_pubmed),
            Some((whole.as_str(), true))
        );

        // The length falls at the end of the 128th identifier, which cannot
        // be told from a cut one.
        let at_an_end = format!("{exact},1000999");
        let before = identifiers("10000000", 127);
        assert_eq!(
            stored_list(&at_an_end, written_pubmed),
            Some((before.as_str(), true))
        );

        // The length falls after the comma that follows the 128th.
        let after_a_comma = identifiers("1000000", 300);
        let stored = stored_list(&after_a_comma, written_pubmed).unwrap();
        assert_eq!(stored, (&after_a_comma[..1024], true));
        assert!(stored.0.ends_with(','));

        assert_eq!(pubmed("0"), None);
        assert_eq!(pubmed(""), None);
        assert_eq!(pubmed("0,0"), Some((3, false)));
    }

    /// `count` names of eight bytes joined by `,`.
    fn names(count: usize) -> String {
        (0..count)
            .map(|n| format!("name{n:04}"))
            .collect::<Vec<String>>()
            .join(",")
    }

    /// A synonym text over the stored length is stored as `vep` prints its
    /// first 1,024 bytes. Where the length falls inside a name that follows
    /// another of its source, the text ends before that name's comma. Where
    /// it falls inside the only name of its source, or inside a source's
    /// name, that source goes with the `--` before it. Where it falls after
    /// a separator, the 1,024 bytes stand. Where it falls inside a character
    /// of two bytes, the name that character belongs to goes. A text whose
    /// first name alone passes the length is no cell, and a text of exactly
    /// 1,024 bytes is stored as it stands. A text ends before a tab or a
    /// line end inside a name.
    #[test]
    fn a_synonym_text_over_the_stored_length_is_stored_as_it_prints() {
        fn synonyms(text: &str) -> Option<(&str, bool)> {
            stored_list(text, written_synonyms)
        }
        // 6 bytes of source, then nine a name: the length falls one byte
        // into the 114th name.
        let inside_a_name = format!("alph::{}", names(200));
        let whole = format!("alph::{}", names(113));
        assert_eq!(whole.len(), 1022);
        assert_eq!(synonyms(&inside_a_name), Some((whole.as_str(), true)));

        // 1,013 bytes of the first source, then `--Zeta::` and three bytes
        // of its only name.
        let first = format!("alph::{}", names(112));
        assert_eq!(first.len(), 1013);
        let inside_an_only_name = format!("{first}--Zeta::{}", names(4));
        assert_eq!(synonyms(&inside_an_only_name), Some((first.as_str(), true)));

        let inside_a_source = format!("{first}--ALongSourceName::{}", names(4));
        assert_eq!(synonyms(&inside_a_source), Some((first.as_str(), true)));

        let after_a_separator = format!("alph::{}--Zeta::{}", names(113), names(4));
        let stored = synonyms(&after_a_separator).unwrap();
        assert_eq!(stored, (&after_a_separator[..1024], true));
        assert!(stored.0.ends_with("name0112--"));

        // The 1,024th byte is the second of a two-byte character.
        let inside_a_character = format!("{first},abcdefghi\u{e9}tude,name");
        assert!(!inside_a_character.is_char_boundary(1024));
        assert_eq!(synonyms(&inside_a_character), Some((first.as_str(), true)));

        let one_long_name = format!("alpha::{}", "x".repeat(1100));
        assert_eq!(synonyms(&one_long_name), None);

        let exact = format!("alph::{},n", names(113));
        assert_eq!(exact.len(), 1024);
        assert_eq!(synonyms(&exact), Some((exact.as_str(), false)));

        assert_eq!(synonyms("alpha::a\tb,c"), Some(("alpha::a", false)));
        assert_eq!(
            synonyms("alpha::a,b\nc--Zeta::z"),
            Some(("alpha::a,b", false))
        );
    }

    /// Lists over the stored length read from tables: a variation cited by
    /// 300 publications and named 200 times under one source is written
    /// with the whole items of each list's first 1,024 bytes, and the facts
    /// count one record of each kind; the synonyms stand from the last row
    /// back.
    #[test]
    fn long_lists_of_the_tables_are_stored_within_the_length() {
        let publication: String = (1..=300)
            .map(|id| publication_row(id, Some(if id == 1 { 10_000 } else { 1_000_000 + id })))
            .collect();
        let citation: String = (1..=300)
            .rev()
            .map(|id| format!("7\t{id}\t\\N\n"))
            .collect();
        let synonym: String = (0..200)
            .map(|n| synonym_row(n + 1, 7, 3, &format!("name{:04}", 199 - n)))
            .collect();
        let cells = Dumps::new()
            .table("variation", &variation_row(7, "mv7"))
            .table("publication", &publication)
            .table("variation_citation", &citation)
            .table("source", &source_row(3, "alph"))
            .table("variation_synonym", &synonym)
            .read();
        let (failed, pubmed, synonyms) = cells_of(&cells, "mv7");
        assert_eq!(failed, 0);
        let expected: Vec<String> = std::iter::once("10000".to_string())
            .chain((2..=128).map(|id| (1_000_000 + id).to_string()))
            .collect();
        assert_eq!(pubmed, Some(expected.join(",")));
        assert_eq!(synonyms, Some(format!("alph::{}", names(113))));
        let facts = cells.facts();
        assert_eq!(
            (
                facts.records_with_pubmed,
                facts.pubmed_lists_at_whole_items,
                facts.records_with_synonyms,
                facts.synonym_lists_at_whole_items
            ),
            (1, 1, 1, 1)
        );
    }

    /// A cell is read where every table of it is named. Without `variation`
    /// no table is read; a cell whose own table is not named stays empty
    /// while the others are read; and inputs that name no file read nothing.
    #[test]
    fn a_table_that_is_not_named_switches_its_cell_off() {
        let all = made_up_tables();
        let read = |inputs: &VariantCellInputs| {
            let cells = VariantCells::read(inputs, &[]).unwrap();
            let facts = cells.facts();
            (
                [
                    facts.pubmed_tables_read,
                    facts.synonym_tables_read,
                    facts.failed_variation_joined,
                    facts.failed_gvf_read,
                ],
                cells_of(&cells, "mv10"),
                cells_of(&cells, "mv40"),
            )
        };
        let nothing = (0, None, None);

        let no_variation = VariantCellInputs {
            variation: None,
            ..all.inputs.clone()
        };
        assert_eq!(
            read(&no_variation),
            ([false; 4], nothing.clone(), nothing.clone())
        );
        assert!(!no_variation.scratch.exists());

        let no_publication = VariantCellInputs {
            publication: None,
            ..all.inputs.clone()
        };
        assert_eq!(
            read(&no_publication),
            (
                [false, true, true, false],
                (0, None, text("alpha::a3,a2,a1--Zeta::z2,z1")),
                (1, None, None)
            )
        );

        let no_source = VariantCellInputs {
            source: None,
            failed_variation: None,
            ..all.inputs.clone()
        };
        assert_eq!(
            read(&no_source),
            (
                [true, false, false, false],
                (0, text("900,300"), None),
                nothing.clone()
            )
        );

        let failed_alone = VariantCellInputs {
            variation: all.inputs.variation.clone(),
            failed_variation: all.inputs.failed_variation.clone(),
            scratch: all.inputs.scratch.clone(),
            ..Default::default()
        };
        assert_eq!(
            read(&failed_alone),
            (
                [false, false, true, false],
                nothing.clone(),
                (1, None, None)
            )
        );

        let none = VariantCellInputs {
            scratch: all.inputs.scratch.clone(),
            ..Default::default()
        };
        assert_eq!(read(&none), ([false; 4], nothing.clone(), nothing));
    }

    /// A value that is no identifier where the schema has one refuses the
    /// read, naming the file, the row, the column and the value.
    #[test]
    fn a_value_that_is_no_identifier_refuses_the_read() {
        let dumps = made_up_tables().table("variation_citation", "10\t7\t\\N\nten\t3\t\\N\n");
        let path = dumps.inputs.variation_citation.clone().unwrap();
        let error = VariantCells::read(&dumps.inputs, &[]).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            format!(
                "{} row 2 holds no identifier in variation_id: \"ten\"",
                path.display()
            )
        );

        let dumps = made_up_tables().table("failed_variation", "1\t\\N\t3\n");
        let path = dumps.inputs.failed_variation.clone().unwrap();
        let error = VariantCells::read(&dumps.inputs, &[]).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            format!(
                "{} row 1 holds no identifier in variation_id: NULL",
                path.display()
            )
        );
    }

    /// The identifiers wanted are found whatever the order they are asked
    /// in: ascending, as a sorted table asks, and in no order.
    #[test]
    fn wanted_identifiers_are_found_in_any_order_of_asking() {
        let ids = vec![40, 7, 7, 900, 12, 41];
        let held = |id: u32| ids.contains(&id);
        let mut wanted = Wanted::new(ids.clone());
        assert_eq!(wanted.ids, [7, 12, 40, 41, 900]);
        for id in 0..1000 {
            assert_eq!(wanted.contains(id), held(id), "{id}");
        }
        for id in (0..1000).rev().chain([41, 3, 900, 0, 40, 40, 7, 1000, 12]) {
            assert_eq!(wanted.contains(id), held(id), "{id}");
        }
        assert!(!Wanted::new(Vec::new()).contains(5));
    }
}
