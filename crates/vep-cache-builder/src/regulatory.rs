// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The regulatory tree of a native cache, baked from the regulatory GFF files
//! the Ensembl site distributes (`--regulatory-gff`, `--emars-gff`,
//! `--motif-gff`): `regulatory/<chr>/<start>-<end>.json` shards in the shape
//! `vep-cli` reads, one feature per record, and `regulatory/binding_matrices.json`
//! holding an entry per motif matrix. The records read as Ensembl VEP's GFF
//! regulatory source reads them (`AnnotationSource/File/RegFeat.pm` at
//! release/116.3), through the vocabulary of [`vep_core::regulatory_gff`].
//!
//! From the GFF files alone a regulatory feature carries no activity and a
//! matrix's entry is a stub (its identifier, the motifs' span and their
//! factors, no counts), so a run with `--regulatory` over the cache gives the
//! rows a run with `--regulatory_gff` over the same files gives. Three further
//! inputs bring what a conversion of Ensembl VEP's own cache carries, each
//! where the build is given it: the activity file published beside the GFF
//! files gives every regulatory feature one activity code per cell type
//! ([`RegulatoryInputs::activity`]), the table dumps of the funcgen database
//! give every matrix its counts ([`RegulatoryInputs::funcgen_dumps`]), and the
//! genome FASTA gives every motif its sequence. A run reads `CELL_TYPE`,
//! `HIGH_INF_POS` and `MOTIF_SCORE_CHANGE` from those three.

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use serde::Serialize;
use tracing::{info, warn};
use vep_core::binding_matrix::NUCLEOTIDES;
use vep_core::codon::reverse_complement;
use vep_core::coordinate::is_cache_path_component;
use vep_core::regulatory_gff::{
    decode_attribute, decode_type, duplicate_identifier_warning, feature_id, ignored_bound_warning,
    transcription_factors, unsupported_type_warning, EXTENDED_BOUNDS_TYPES, MOTIF_FEATURE_TYPES,
    REGULATORY_FEATURE_TYPES,
};

use crate::builder::{region_start, REGION_SIZE};
use crate::core_dumps::{read_table, table_path};
use crate::count::{counted, counted_as, one_or};
use crate::gff3::SequenceNames;
use crate::line::{check_indexed_coordinates, index_build_places, lines, tab_columns};

/// The two funcgen tables a matrix's counts are read from: the matrices, and
/// the count of each nucleotide at each position of each.
pub const MATRIX_TABLES: [&str; 2] = ["binding_matrix", "binding_matrix_frequencies"];

/// The two funcgen tables a matrix's factors are read from: the complexes
/// that bind each matrix, and each complex's names.
pub const COMPLEX_TABLES: [&str; 2] = [
    "binding_matrix_transcription_factor_complex",
    "transcription_factor_complex",
];

/// The funcgen table the regulatory build's version is read from.
pub const REGULATORY_BUILD_TABLE: &str = "regulatory_build";

/// The regulatory inputs of a build.
#[derive(Debug, Clone, Default)]
pub struct RegulatoryInputs {
    /// The regulatory features (`regulatory_features.gff3.gz`).
    pub regulatory_gff: Option<PathBuf>,
    /// The EMAR features (`EMARs.gff.gz`), regulatory features of their own
    /// identifier space.
    pub emars_gff: Option<PathBuf>,
    /// The motif features (`motif_features.gff3.gz`).
    pub motif_gff: Option<PathBuf>,
    /// Write each promoter at its `extended_start` and `extended_end`.
    pub extended_promoters: bool,
    /// The cell-type activity of the regulatory features
    /// (`regulatory_activity.tsv.gz`): one row a feature, one column a cell
    /// type. `None` writes every feature with an empty activity string.
    pub activity: Option<PathBuf>,
    /// A directory of table dumps of the funcgen database. Where it holds
    /// both of [`MATRIX_TABLES`], each matrix a motif names is written with
    /// its counts; where it also holds both of [`COMPLEX_TABLES`], with the
    /// names of the complexes that bind it; where it holds
    /// [`REGULATORY_BUILD_TABLE`], the build's version is read. `None` writes
    /// each matrix as its stub.
    pub funcgen_dumps: Option<PathBuf>,
}

impl RegulatoryInputs {
    pub fn is_empty(&self) -> bool {
        self.regulatory_gff.is_none() && self.emars_gff.is_none() && self.motif_gff.is_none()
    }
}

/// What the regulatory build wrote.
#[derive(Debug, Default)]
pub struct RegulatoryStats {
    pub regulatory_features: usize,
    pub motif_features: usize,
    pub matrices: usize,
    pub shards: usize,
    pub chromosomes: BTreeSet<String>,
    /// The cell types of the activity file, in the order every activity
    /// string follows, the list `info.json` carries as `cell_types`; `None`
    /// when no activity file was read.
    pub cell_types: Option<Vec<String>>,
    /// Whether no regulatory feature was written without the activity an
    /// activity file gives it: one code per name of `cell_types`, `-` for
    /// each where the file holds no row of the feature. True of a tree with
    /// no regulatory feature.
    pub activity_on_every_feature: bool,
    /// Whether no matrix a motif names was written without its counts. True
    /// of a tree with no motif feature.
    pub counts_on_every_matrix: bool,
    /// Whether no motif feature was written without its sequence. True of a
    /// tree with no motif feature.
    pub sequence_on_every_motif: bool,
    /// The version of the regulatory build, read from
    /// [`REGULATORY_BUILD_TABLE`]; `None` when the build was given no such
    /// table or the table states no version.
    pub regbuild: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Regulatory,
    Motif,
}

/// A shard element in the shape `vep-cli`'s regulatory loader reads
/// (`regulatory/source.rs`), the shape `info.json` declares as
/// [`vep_core::cache_info::REGULATORY_SCHEMA`]: a regulatory feature with one
/// activity code per cell type the cache lists, an empty string where it
/// lists none; a motif naming its matrix by identifier, with its sequence on
/// its own strand where the build's genome FASTA holds its span, and with
/// neither a database id nor a score.
/// Both carry the file they came from (`source`, the input's slot: the
/// regulatory file 0, the EMAR file 1, the motif file 2) and the record's
/// place in that file (`ordinal`), the two keys Ensembl VEP's GFF source
/// gives a record (one identifier set per file, `RegFeat.pm:485`; records
/// met in file order, `:312-328`), which the shard's own order cannot carry
/// once a widened promoter is binned by its widened start.
#[derive(Serialize)]
#[serde(tag = "feature_type")]
enum ShardFeature {
    RegulatoryFeature {
        stable_id: String,
        start: u64,
        end: u64,
        biotype: String,
        activity: String,
        source: u8,
        ordinal: u32,
    },
    MotifFeature {
        stable_id: String,
        start: u64,
        end: u64,
        strand: i8,
        matrix: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        seq: Option<String>,
        source: u8,
        ordinal: u32,
    },
}

impl ShardFeature {
    fn start(&self) -> u64 {
        match self {
            ShardFeature::RegulatoryFeature { start, .. }
            | ShardFeature::MotifFeature { start, .. } => *start,
        }
    }
}

/// What the motif GFF gives of a matrix: the span of its motifs and their
/// factors, the stub Ensembl VEP builds for a GFF motif
/// (`RegFeat.pm:440-444`).
#[derive(Debug, PartialEq, Eq)]
struct MatrixStub {
    length: u64,
    transcription_factors: Vec<String>,
}

/// One entry of `binding_matrices.json`, in the form a conversion of Ensembl
/// VEP's cache writes it. A matrix read from the funcgen tables carries its
/// counts (`frequencies`, one `[A, C, G, T]` row a position, position 1
/// first, the integers of the table), its `name`, `source` and `threshold`,
/// each absent where the table holds NULL, and `length`, the number of its
/// positions. A matrix known from the motif GFF alone is its stub: `length`
/// and `transcription_factors`.
#[derive(Debug, Serialize, PartialEq)]
struct MatrixEntry {
    #[serde(skip_serializing_if = "Option::is_none")]
    frequencies: Option<Vec<[u64; 4]>>,
    length: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    source: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    threshold: Option<f64>,
    transcription_factors: Vec<String>,
}

impl MatrixEntry {
    fn of_stub(stub: MatrixStub) -> Self {
        Self {
            frequencies: None,
            length: stub.length,
            name: None,
            source: None,
            threshold: None,
            transcription_factors: stub.transcription_factors,
        }
    }
}

/// The features read so far, by chromosome, and the matrices they name.
#[derive(Default)]
struct Features {
    by_chromosome: BTreeMap<String, Vec<ShardFeature>>,
    matrices: BTreeMap<String, MatrixStub>,
    regulatory_features: usize,
    motif_features: usize,
}

/// Reads the files `inputs` names and writes the cache's regulatory tree under
/// `output_dir`: the shards, by chromosome and 1 Mb region of each feature's
/// start, and the matrices table. A motif whose matrix identifier an earlier
/// motif gave another span or factor list refuses the build, since the table
/// holds one entry per identifier. A feature's chromosome is named as `names`
/// names its seqid, the rule the transcript model's sequences are named by.
///
/// With an activity file every regulatory feature carries its activity
/// (`read_activity`). With the funcgen tables every matrix a motif names
/// that the tables hold carries its counts (`matrix_entries`), and the
/// matrices they lack, written as their stubs, are counted in one warning.
/// With `genome`, the build's genome FASTA, every motif whose span the FASTA
/// holds carries its sequence (`motif_sequence`), and the motifs whose span
/// it does not hold are counted in one warning. The returned facts say what
/// of the three the tree carries.
pub fn build_regulatory(
    inputs: &RegulatoryInputs,
    output_dir: &Path,
    names: &SequenceNames,
    genome: Option<&vep_fasta::IndexedFasta>,
) -> Result<RegulatoryStats> {
    let mut features = Features::default();
    for (source, (kind, path)) in [
        (Kind::Regulatory, &inputs.regulatory_gff),
        (Kind::Regulatory, &inputs.emars_gff),
        (Kind::Motif, &inputs.motif_gff),
    ]
    .into_iter()
    .enumerate()
    {
        if let Some(path) = path {
            info!("Reading regulatory GFF: {}", path.display());
            read_file(
                path,
                kind,
                source as u8,
                inputs.extended_promoters,
                names,
                &mut features,
            )
            .with_context(|| format!("Failed to read {}", path.display()))?;
        }
    }

    let cell_types = match &inputs.activity {
        Some(path) => {
            info!("Reading regulatory activity: {}", path.display());
            let activity = read_activity(path)
                .with_context(|| format!("Failed to read {}", path.display()))?;
            write_activity(&mut features, &activity, path);
            Some(activity.cell_types)
        }
        None => None,
    };
    let without_sequence = genome.map(|genome| write_sequences(&mut features, genome));
    if let Some(without_sequence) = without_sequence.filter(|count| *count > 0) {
        warn!(
            "{} written without a sequence: the genome FASTA does not hold the span",
            counted(without_sequence, "motif feature")
        );
    }
    let dumps = inputs.funcgen_dumps.as_deref();
    let holds = |tables: &[&str]| {
        dumps.filter(|dir| tables.iter().all(|table| table_path(dir, table).is_some()))
    };
    let matrices: BTreeMap<String, MatrixEntry> = match holds(&MATRIX_TABLES) {
        Some(dir) => {
            info!("Reading the binding matrices: {}", dir.display());
            matrix_entries(dir, features.matrices, holds(&COMPLEX_TABLES).is_some())
                .context("Failed to read the binding matrices")?
        }
        None => features
            .matrices
            .into_iter()
            .map(|(id, stub)| (id, MatrixEntry::of_stub(stub)))
            .collect(),
    };
    let regbuild = match holds(&[REGULATORY_BUILD_TABLE]) {
        Some(dir) => regulatory_build_version(dir)
            .context("Failed to read the regulatory build's version")?,
        None => None,
    };

    let root = output_dir.join("regulatory");
    let mut stats = RegulatoryStats {
        regulatory_features: features.regulatory_features,
        motif_features: features.motif_features,
        matrices: matrices.len(),
        activity_on_every_feature: cell_types.is_some() || features.regulatory_features == 0,
        cell_types,
        counts_on_every_matrix: matrices.values().all(|entry| entry.frequencies.is_some()),
        sequence_on_every_motif: features.motif_features == 0 || without_sequence == Some(0),
        regbuild,
        ..Default::default()
    };
    for (chr, chromosome_features) in &features.by_chromosome {
        let chr_dir = root.join(chr);
        if !is_cache_path_component(chr) || chr_dir.parent() != Some(root.as_path()) {
            bail!("invalid chromosome name '{chr}'");
        }
        std::fs::create_dir_all(&chr_dir)
            .with_context(|| format!("Failed to create {}", chr_dir.display()))?;
        let mut by_region: BTreeMap<u64, Vec<&ShardFeature>> = BTreeMap::new();
        for feature in chromosome_features {
            by_region
                .entry(region_start(feature.start()))
                .or_default()
                .push(feature);
        }
        for (region_start, shard) in &by_region {
            let path = chr_dir.join(format!(
                "{}-{}.json",
                region_start,
                region_start + REGION_SIZE - 1
            ));
            let json = serde_json::to_string_pretty(shard)
                .context("Failed to serialize regulatory features")?;
            std::fs::write(&path, json)
                .with_context(|| format!("Failed to write {}", path.display()))?;
            stats.shards += 1;
        }
        stats.chromosomes.insert(chr.clone());
        info!(
            "chr {}: {} in {}",
            chr,
            counted_as(
                chromosome_features.len(),
                "regulatory or motif feature",
                "regulatory and motif features"
            ),
            counted(by_region.len(), "region")
        );
    }
    std::fs::create_dir_all(&root)
        .with_context(|| format!("Failed to create {}", root.display()))?;
    let table = root.join("binding_matrices.json");
    let json = serde_json::to_string_pretty(&matrices)
        .context("Failed to serialize the binding matrices table")?;
    std::fs::write(&table, json).with_context(|| format!("Failed to write {}", table.display()))?;
    info!(
        "Regulatory tree written: {}, {} naming {}, {}",
        counted(stats.regulatory_features, "regulatory feature"),
        counted(stats.motif_features, "motif feature"),
        counted_as(stats.matrices, "matrix", "matrices"),
        counted(stats.shards, "shard")
    );
    Ok(stats)
}

/// The text of the file at `path`: gzip or bgzip compressed where its name
/// ends in `.gz`, else plain.
fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    Ok(if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

/// What an activity file holds: its cell types and each feature's codes.
struct Activity {
    /// The names of the header, every space written `_`, in byte order: the
    /// keys Ensembl VEP's cache holds a feature's activity under
    /// (`AnnotationSource/Database/RegFeat.pm:196-203`), in the order its
    /// `cell_types` line lists them (`:112-127`).
    cell_types: Vec<String>,
    /// The codes of each row by its feature identifier, one code per name of
    /// `cell_types`.
    codes: HashMap<String, String>,
}

/// The code a shard holds for a value of the activity file, the values of
/// the funcgen `regulatory_activity` table; `None` for any other text.
fn activity_code(value: &str) -> Option<char> {
    match value {
        "ACTIVE" => Some('A'),
        "INACTIVE" => Some('I'),
        "POISED" => Some('P'),
        "REPRESSED" => Some('R'),
        "NA" => Some('N'),
        _ => None,
    }
}

/// Reads an activity file: tab-separated, a header line of `feature_id`,
/// `feature_type` and one name a cell type, then one row a regulatory
/// feature holding its identifier, its type, which is not read, and its
/// activity in each cell type. The file is gzip or bgzip compressed, or
/// plain, and an empty line is passed over. Each of these refuses the read,
/// naming the line: a header that opens with other columns, a name the
/// header holds twice once its spaces are written `_`, a row of another
/// number of columns than the header, an activity that is none of the five
/// values (`activity_code`), and a second row of one feature.
fn read_activity(path: &Path) -> Result<Activity> {
    let file_name = path.display().to_string();
    let mut cell_types: Vec<String> = Vec::new();
    // The header column of each name of `cell_types`, counted from the
    // first name.
    let mut columns: Vec<usize> = Vec::new();
    let mut header_columns = None;
    let mut codes: HashMap<String, String> = HashMap::new();
    for (index, line) in lines(open_text(path)?, path).enumerate() {
        let line = line?;
        let line_number = index + 1;
        let text = line.strip_suffix(b"\r").unwrap_or(&line);
        if text.is_empty() {
            continue;
        }
        let text = std::str::from_utf8(text)
            .with_context(|| format!("line {line_number} of {file_name} is not UTF-8"))?;
        let cells: Vec<&str> = text.split('\t').collect();
        let Some(expected) = header_columns else {
            if cells.len() < 2 || cells[..2] != ["feature_id", "feature_type"] {
                bail!(
                    "line {line_number} of {file_name} opens with the columns {:?} where the \
                     header of an activity file opens with feature_id and feature_type",
                    &cells[..cells.len().min(2)]
                );
            }
            let names: Vec<String> = cells[2..]
                .iter()
                .map(|name| name.replace(' ', "_"))
                .collect();
            columns = (0..names.len()).collect();
            columns.sort_by(|a, b| names[*a].cmp(&names[*b]));
            let twice = columns
                .windows(2)
                .find(|pair| names[pair[0]] == names[pair[1]]);
            if let Some(pair) = twice {
                bail!(
                    "line {line_number} of {file_name} names the cell type {:?} twice, each \
                     space read as an underscore",
                    names[pair[0]]
                );
            }
            cell_types = columns
                .iter()
                .map(|column| names[*column].clone())
                .collect();
            header_columns = Some(cells.len());
            continue;
        };
        if cells.len() != expected {
            bail!(
                "line {line_number} of {file_name} has {} columns where the header has {expected}",
                cells.len()
            );
        }
        let feature = cells[0];
        let mut row = String::with_capacity(columns.len());
        for (name, column) in cell_types.iter().zip(&columns) {
            let value = cells[2 + column];
            let Some(code) = activity_code(value) else {
                bail!(
                    "line {line_number} of {file_name} holds the activity {value:?} of {feature} \
                     in {name}, which is not ACTIVE, INACTIVE, POISED, REPRESSED or NA"
                );
            };
            row.push(code);
        }
        if codes.insert(feature.to_string(), row).is_some() {
            bail!("line {line_number} of {file_name} is a second row of the feature {feature}");
        }
    }
    if header_columns.is_none() {
        bail!("{file_name} holds no header line");
    }
    Ok(Activity { cell_types, codes })
}

/// Gives every regulatory feature read so far its activity: the codes of its
/// row of the activity file, two records under one identifier taking the one
/// row, and `-` for every cell type where the file holds no row of the
/// feature, as Ensembl VEP's cache holds an empty map for a feature the
/// funcgen database states no activity of. A row whose identifier no
/// feature carries is passed over, and such rows are counted in one line.
fn write_activity(features: &mut Features, activity: &Activity, path: &Path) {
    let absent = "-".repeat(activity.cell_types.len());
    let mut read: HashSet<&str> = HashSet::new();
    let mut without_row = 0usize;
    for feature in features.by_chromosome.values_mut().flatten() {
        let ShardFeature::RegulatoryFeature {
            stable_id,
            activity: codes,
            ..
        } = feature
        else {
            continue;
        };
        match activity.codes.get_key_value(stable_id.as_str()) {
            Some((id, row)) => {
                read.insert(id.as_str());
                codes.clone_from(row);
            }
            None => {
                without_row += 1;
                codes.clone_from(&absent);
            }
        }
    }
    info!(
        "Regulatory activity read: {} in {}, {} without a row",
        counted(activity.codes.len(), "row"),
        counted(activity.cell_types.len(), "cell type"),
        counted(without_row, "regulatory feature")
    );
    let passed_over = activity.codes.len() - read.len();
    if passed_over > 0 {
        info!(
            "{} of {} passed over: no regulatory feature read carries {}",
            counted(passed_over, "row"),
            path.display(),
            one_or(passed_over, "its identifier", "their identifiers")
        );
    }
}

/// Gives every motif read so far the bases of its span in `genome`, on the
/// motif's own strand (`motif_sequence`), and returns the number of motifs
/// whose span `genome` does not hold, each left without a sequence.
fn write_sequences(features: &mut Features, genome: &vep_fasta::IndexedFasta) -> usize {
    let mut without = 0;
    for (chr, chromosome_features) in &mut features.by_chromosome {
        for feature in chromosome_features {
            if let ShardFeature::MotifFeature {
                start,
                end,
                strand,
                seq,
                ..
            } = feature
            {
                *seq = motif_sequence(genome, chr, *start, *end, *strand);
                without += usize::from(seq.is_none());
            }
        }
    }
    without
}

/// The sequence Ensembl VEP's cache holds for a motif feature, the feature's
/// own `seq` (`Pipeline/DumpVEP/Dumper/Regulation.pm:137`): the bases of
/// `start..=end` of `chr` in upper case, reverse-complemented for a motif on
/// the reverse strand. `None` where `genome` does not hold the whole span.
fn motif_sequence(
    genome: &vep_fasta::IndexedFasta,
    chr: &str,
    start: u64,
    end: u64,
    strand: i8,
) -> Option<String> {
    let bases = genome.sequence(chr, start, end)?;
    let bases = if strand < 0 {
        reverse_complement(&bases)
    } else {
        bases
    };
    String::from_utf8(bases).ok()
}

/// What the funcgen tables hold of one matrix a motif names.
struct TableMatrix {
    stable_id: String,
    name: Option<String>,
    threshold: Option<f64>,
    source: Option<String>,
    /// The count of each nucleotide, in the order of [`NUCLEOTIDES`], by
    /// position.
    counts: BTreeMap<u64, [Option<u64>; 4]>,
}

/// The entry of each matrix of `stubs`, the matrices the motifs name, from
/// the funcgen tables under `dir`. A matrix `binding_matrix` holds under its
/// stable identifier is written with its name, source and threshold and with
/// the counts `binding_matrix_frequencies` holds for it, as Ensembl VEP's
/// cache holds a matrix (funcgen `BindingMatrix.pm:240-261`, `_elements`,
/// and `:352-360`, `length`). A matrix the table lacks keeps its stub, and
/// such matrices are counted in one warning. With `complexes`, the factors
/// of a matrix read from the tables are those of `complex_names`; else they
/// are the motif GFF's, as a stub's are.
///
/// Each of these refuses the read, naming the matrix: a stable identifier
/// or a `binding_matrix_id` two rows of `binding_matrix` carry, a threshold
/// that is no number, a nucleotide other than A, C, G and T, a count stated
/// twice, positions that do not run from 1 without a gap, a position that
/// lacks a nucleotide, and a number of positions other than the span of the
/// matrix's motifs.
fn matrix_entries(
    dir: &Path,
    stubs: BTreeMap<String, MatrixStub>,
    complexes: bool,
) -> Result<BTreeMap<String, MatrixEntry>> {
    let [matrix_table, frequency_table] = MATRIX_TABLES;
    // The matrices by `binding_matrix_id`, the key the other tables name a
    // matrix by, and that key by stable identifier.
    let mut by_id: HashMap<String, TableMatrix> = HashMap::new();
    let mut ids: HashMap<String, String> = HashMap::new();
    for row in read_table(dir, matrix_table, 5)? {
        let row = row?;
        let (id, stable_id) = (row.field(0), row.field(4));
        if !stubs.contains_key(stable_id) {
            continue;
        }
        let threshold = match row.optional(2) {
            Some(text) => Some(
                text.parse::<f64>()
                    .ok()
                    .filter(|threshold| threshold.is_finite())
                    .with_context(|| {
                        format!("the threshold of binding matrix {stable_id} is no number: {text}")
                    })?,
            ),
            None => None,
        };
        let matrix = TableMatrix {
            stable_id: stable_id.to_string(),
            name: row.optional(1).map(str::to_string),
            threshold,
            source: row.optional(3).map(str::to_string),
            counts: BTreeMap::new(),
        };
        let twice = ids.insert(stable_id.to_string(), id.to_string()).is_some()
            | by_id.insert(id.to_string(), matrix).is_some();
        if twice {
            bail!(
                "{matrix_table} holds two rows of the binding matrix {stable_id}, or two under \
                 its identifier {id}"
            );
        }
    }
    for row in read_table(dir, frequency_table, 5)? {
        let row = row?;
        let Some(matrix) = by_id.get_mut(row.field(1)) else {
            continue;
        };
        let stable_id = matrix.stable_id.as_str();
        let of_matrix = || {
            format!(
                "{frequency_table} row {} of binding matrix {stable_id}",
                row.field(0)
            )
        };
        let position = row.number(2).with_context(of_matrix)?;
        let count = row.number(4).with_context(of_matrix)?;
        let nucleotide = row.field(3);
        let Some(at) = NUCLEOTIDES
            .iter()
            .position(|letter| nucleotide.as_bytes() == [*letter])
        else {
            bail!(
                "binding matrix {stable_id} holds a count of the nucleotide {nucleotide:?} at \
                 position {position}; the nucleotides are A, C, G and T"
            );
        };
        if matrix.counts.entry(position).or_default()[at]
            .replace(count)
            .is_some()
        {
            bail!(
                "binding matrix {stable_id} holds the count of {nucleotide} at position \
                 {position} twice"
            );
        }
    }
    let factors = complexes
        .then(|| complex_names(dir, |id| by_id.contains_key(id)))
        .transpose()?;

    let mut entries = BTreeMap::new();
    let mut lacking = 0usize;
    for (stable_id, stub) in stubs {
        let held = ids
            .get(&stable_id)
            .and_then(|id| Some((id, by_id.remove(id)?)));
        let Some((id, matrix)) = held else {
            lacking += 1;
            entries.insert(stable_id, MatrixEntry::of_stub(stub));
            continue;
        };
        let mut frequencies = Vec::with_capacity(matrix.counts.len());
        for (expected, (position, counts)) in (1u64..).zip(&matrix.counts) {
            if *position != expected {
                bail!(
                    "the positions of binding matrix {stable_id} do not run from 1 without a \
                     gap: position {position} stands where position {expected} would"
                );
            }
            let mut row = [0u64; 4];
            for (at, count) in counts.iter().enumerate() {
                row[at] = count.with_context(|| {
                    format!(
                        "binding matrix {stable_id} holds no count of {} at position {position}",
                        char::from(NUCLEOTIDES[at])
                    )
                })?;
            }
            frequencies.push(row);
        }
        let length = frequencies.len() as u64;
        if length != stub.length {
            bail!(
                "binding matrix {stable_id} holds {} where its motif features span {}",
                counted(length, "position"),
                counted(stub.length, "base")
            );
        }
        let transcription_factors = match &factors {
            Some(names) => names.get(id).cloned().unwrap_or_default(),
            None => stub.transcription_factors,
        };
        entries.insert(
            stable_id,
            MatrixEntry {
                frequencies: Some(frequencies),
                length,
                name: matrix.name,
                source: matrix.source,
                threshold: matrix.threshold,
                transcription_factors,
            },
        );
    }
    info!(
        "Binding matrices read: {} with counts",
        counted_as(entries.len() - lacking, "matrix", "matrices")
    );
    if lacking > 0 {
        warn!(
            "{} named by a motif feature {} not in the {matrix_table} table of {}, and {} \
             written without counts",
            counted_as(lacking, "binding matrix", "binding matrices"),
            one_or(lacking, "is", "are"),
            dir.display(),
            one_or(lacking, "is", "are")
        );
    }
    Ok(entries)
}

/// The display names of the complexes that bind each matrix whose
/// `binding_matrix_id` is `wanted`, by that identifier, in the order of the
/// complexes' identifiers, the factors Ensembl VEP's cache lists for a
/// matrix. Funcgen fetches the complexes a matrix's links name by a list of
/// their identifiers (`DBSQL/TranscriptionFactorComplexAdaptor.pm:125-145`),
/// so a link to a complex the table lacks names nothing.
fn complex_names(
    dir: &Path,
    wanted: impl Fn(&str) -> bool,
) -> Result<HashMap<String, Vec<String>>> {
    let [link_table, complex_table] = COMPLEX_TABLES;
    let mut links: HashMap<String, BTreeSet<u64>> = HashMap::new();
    for row in read_table(dir, link_table, 3)? {
        let row = row?;
        let matrix = row.field(1);
        if wanted(matrix) {
            links
                .entry(matrix.to_string())
                .or_default()
                .insert(row.number(2)?);
        }
    }
    let mut names: HashMap<u64, String> = HashMap::new();
    for row in read_table(dir, complex_table, 3)? {
        let row = row?;
        names.insert(row.number(0)?, row.field(2).to_string());
    }
    Ok(links
        .into_iter()
        .map(|(matrix, complexes)| {
            let listed = complexes
                .iter()
                .filter_map(|complex| names.get(complex).cloned())
                .collect();
            (matrix, listed)
        })
        .collect())
}

/// The version the funcgen `regulatory_build` table under `dir` states: that
/// of the row marked current (`is_current`), else of the table's first row.
/// `None` for a table with no row and for a NULL version.
fn regulatory_build_version(dir: &Path) -> Result<Option<String>> {
    let mut first = None;
    for row in read_table(dir, REGULATORY_BUILD_TABLE, 11)? {
        let row = row?;
        let version = row.optional(4).map(str::to_string);
        if row.field(9) == "1" {
            return Ok(version);
        }
        first.get_or_insert(version);
    }
    Ok(first.flatten())
}

/// Reads one file's records into `features` as `_record_to_regfeat` builds
/// them (`RegFeat.pm:349-385`): a type the kind does not accept is warned about
/// once and skipped (`:314-323`), its start and end read only as
/// [`check_indexed_coordinates`] reads a skipped line's, the reading a kept
/// line that holds a NUL byte is held to as well, and no other column read
/// at all, so a byte that is not UTF-8 refuses nothing on such a line, as it
/// refuses nothing there at run time; a record of an accepted type is UTF-8
/// or refused by its line; a line of fewer
/// than nine columns is refused, unless its sequence name holds a NUL byte
/// and the index build places it ([`index_build_places`]): that build files
/// the line under the name before the NUL byte, where a run finds no record
/// of the chromosome in it, so the line is skipped and its type, the third
/// column or none, warned about; a line of an accepted type and nine
/// columns that the build so files is skipped too, as a run passes over it;
/// a promoter is widened to its extended bounds
/// under `extended_promoters` (`:360-365`), a bound that is not a whole number
/// ignored with a warning; a motif's strand is the column's, `0` when it
/// holds none (`:369-370`), and a regulatory feature carries none, as the
/// shard shape has none; the identifier is `ID`, else `Name`, else the
/// chromosome as the file spells it with the coordinates and type
/// (`:575-584`); a record the file holds twice is kept once (`:628-640`), and
/// a distinct record under a used identifier is kept with the collision warned
/// about once per identifier (`:478-497`). Each feature carries `source` and
/// its `ordinal` in the file. The file is gzip or bgzip compressed, or plain.
/// A byte-order mark at the start of a line is taken off before the line is
/// read ([`crate::line::without_byte_order_marks`]).
fn read_file(
    path: &Path,
    kind: Kind,
    source: u8,
    extended_promoters: bool,
    names: &SequenceNames,
    features: &mut Features,
) -> Result<()> {
    let reader = open_text(path)?;
    let accepted: &[&str] = match kind {
        Kind::Regulatory => &REGULATORY_FEATURE_TYPES,
        Kind::Motif => &MOTIF_FEATURE_TYPES,
    };
    let file_name = path.display().to_string();
    let mut warned_types: HashSet<String> = HashSet::new();
    let mut seen = SeenRecords::default();
    // The identifiers this file has yielded, and the ones warned about: a
    // distinct record under a used identifier is a collision (`RegFeat.pm:487`),
    // warned about once per identifier and file.
    let mut seen_ids: HashSet<String> = HashSet::new();
    let mut warned_ids: HashSet<String> = HashSet::new();
    let mut ordinal: u32 = 0;
    for (index, line) in lines(reader, path).enumerate() {
        let line = line?;
        let text = line.strip_suffix(b"\r").unwrap_or(&line);
        if text.is_empty() || text[0] == b'#' {
            continue;
        }
        let line_number = index + 1;
        let tab_ended = || text.split(|byte| *byte == b'\t');
        let written = tab_ended().count();
        // The index build files a line it places under the name before the
        // NUL byte, and a run finds no record of that chromosome in it.
        let names_no_chromosome =
            || tab_ended().next().is_some_and(|name| name.contains(&0)) && index_build_places(text);
        if written < 9 {
            // `RegFeat.pm:314-322` warns of the line's type and skips it.
            if names_no_chromosome() {
                let feature_type = String::from_utf8_lossy(tab_ended().nth(2).unwrap_or_default());
                if warned_types.insert(feature_type.to_string()) {
                    warn!("{}", unsupported_type_warning(&feature_type, &file_name));
                }
                continue;
            }
            bail!(
                "line {line_number} of {file_name} has {written} columns where a GFF record has 9"
            );
        }
        let coordinate = |name: &str, value: &str| -> Result<u64> {
            value.parse().with_context(|| {
                format!(
                    "the {name} of line {line_number} of {file_name} is not a whole number: {value}"
                )
            })
        };
        let feature_type = decode_type(tab_ended().nth(2).unwrap_or_default());
        let feature_type = feature_type.as_ref();
        if !accepted.contains(&feature_type) {
            // A skipped line's start and end are held to the one rule the
            // transcript readers hold them to. That rule reads no byte that
            // is not UTF-8, so the line is handed to it with each replaced.
            check_indexed_coordinates(&String::from_utf8_lossy(text), coordinate)?;
            if warned_types.insert(feature_type.to_string()) {
                warn!("{}", unsupported_type_warning(feature_type, &file_name));
            }
            continue;
        }
        if names_no_chromosome() {
            continue;
        }
        let text = std::str::from_utf8(text)
            .with_context(|| format!("line {line_number} of {file_name} is not UTF-8"))?;
        // A NUL byte ends a column to the index build alone, so a kept line
        // that holds one is held to that build's reading of its start and
        // its end, as a skipped line is.
        if text.contains('\0') {
            check_indexed_coordinates(text, coordinate)?;
        }
        let columns: Vec<&str> = tab_columns(text).collect();
        let chr = names.name(columns[0]);
        if !is_cache_path_component(chr) {
            bail!(
                "invalid chromosome name '{}' at line {line_number}",
                columns[0]
            );
        }
        let mut start = coordinate("start", columns[3])?;
        let mut end = coordinate("end", columns[4])?;
        if start == 0 || end < start {
            bail!(
                "GFF coordinates {start}-{end} at line {line_number} are not 1-based and ordered"
            );
        }
        // Checked before any widening, on the start column the file is
        // sorted by.
        if !seen.insert(chr, start, text) {
            continue;
        }
        // The hash holds each value as written, the text Ensembl VEP's tests
        // read; a value is decoded where it is read (`decode_attribute`).
        let attributes = crate::gff3::parse_attributes(columns[8]);
        crate::gff3::check_attributes(columns[8], &attributes)
            .with_context(|| format!("line {line_number} of {file_name}"))?;
        let attribute = |key: &str| attributes.get(key).map(String::as_str);
        let decoded = |key: &str| attribute(key).map(|value| decode_attribute(value.as_bytes()));
        let mut ignored_bounds = Vec::new();
        if extended_promoters && EXTENDED_BOUNDS_TYPES.contains(&feature_type) {
            let bound = |key: &str| decoded(key).map(|value| (value.parse::<u64>(), value));
            match bound("extended_start") {
                Some((Ok(extended), _)) => start = start.min(extended),
                Some((Err(_), value)) => ignored_bounds.push(("extended_start", value)),
                None => {}
            }
            match bound("extended_end") {
                Some((Ok(extended), _)) => end = end.max(extended),
                Some((Err(_), value)) => ignored_bounds.push(("extended_end", value)),
                None => {}
            }
        }
        // A coordinate-derived identifier carries the file's own spelling of
        // the chromosome, the name the run-time reader gives the record
        // (`AnnotationSource.pm:179,198`, `RegFeat.pm:372,583`); the shard
        // directory alone takes the cache's. The identifier is chosen on the
        // values as written, and the attribute value chosen is decoded after;
        // the coordinate form is the file's text as the run-time reader
        // writes it.
        let stable_id = match feature_id(
            attribute("ID"),
            attribute("Name"),
            columns[0],
            start,
            end,
            feature_type,
        ) {
            Cow::Borrowed(value) => decode_attribute(value.as_bytes()).into_owned(),
            Cow::Owned(coordinates) => coordinates,
        };
        for (key, value) in ignored_bounds {
            warn!(
                "{}",
                ignored_bound_warning(key, &value, feature_type, &stable_id, &file_name)
            );
        }
        if !seen_ids.insert(stable_id.clone()) && warned_ids.insert(stable_id.clone()) {
            warn!("{}", duplicate_identifier_warning(&stable_id, &file_name));
        }
        let feature = match kind {
            Kind::Regulatory => {
                features.regulatory_features += 1;
                ShardFeature::RegulatoryFeature {
                    stable_id,
                    start,
                    end,
                    biotype: feature_type.to_string(),
                    activity: String::new(),
                    source,
                    ordinal,
                }
            }
            Kind::Motif => {
                let matrix = decoded("binding_matrix_id")
                    .unwrap_or_default()
                    .into_owned();
                let stub = MatrixStub {
                    length: end - start + 1,
                    transcription_factors: transcription_factors(
                        attribute("transcription_factor").map(str::as_bytes),
                    )
                    .into_iter()
                    .map(Cow::into_owned)
                    .collect(),
                };
                match features.matrices.get(&matrix) {
                    Some(earlier) if *earlier != stub => bail!(
                        "binding matrix {matrix} spans {} bases with factors {:?} at line {line_number} of {file_name} where an earlier motif gave it {} bases with factors {:?}; the table holds one entry per matrix",
                        stub.length,
                        stub.transcription_factors,
                        earlier.length,
                        earlier.transcription_factors
                    ),
                    Some(_) => {}
                    None => {
                        features.matrices.insert(matrix.clone(), stub);
                    }
                }
                features.motif_features += 1;
                ShardFeature::MotifFeature {
                    stable_id,
                    start,
                    end,
                    strand: match columns[6] {
                        "+" => 1,
                        "-" => -1,
                        _ => 0,
                    },
                    matrix,
                    seq: None,
                    source,
                    ordinal,
                }
            }
        };
        features
            .by_chromosome
            .entry(chr.to_string())
            .or_default()
            .push(feature);
        ordinal += 1;
    }
    Ok(())
}

/// The raw records of one chromosome and start position read so far: a
/// tabix-indexed file is sorted by the start column within a chromosome, so a
/// record the file holds twice meets its earlier copy here.
#[derive(Default)]
struct SeenRecords {
    chr: String,
    start: u64,
    records: HashSet<String>,
}

impl SeenRecords {
    /// Whether `record` at `start` on `chr` is new.
    fn insert(&mut self, chr: &str, start: u64, record: &str) -> bool {
        if (self.chr.as_str(), self.start) != (chr, start) {
            self.chr = chr.to_string();
            self.start = start;
            self.records.clear();
        }
        self.records.insert(record.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const FIXTURES: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../vep-cli/tests/fixtures/regulatory-gff"
    );

    fn write_gff(dir: &Path, name: &str, lines: &[impl AsRef<[u8]>]) -> PathBuf {
        let path = dir.join(name);
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(&path).unwrap(),
            flate2::Compression::default(),
        );
        for line in lines {
            encoder.write_all(line.as_ref()).unwrap();
            encoder.write_all(b"\n").unwrap();
        }
        encoder.finish().unwrap();
        path
    }

    fn shard(dir: &Path, chr: &str, name: &str) -> Vec<serde_json::Value> {
        serde_json::from_slice(&std::fs::read(dir.join("regulatory").join(chr).join(name)).unwrap())
            .unwrap()
    }

    /// The three fixture files write the shard layout the runtime reads: one
    /// shard per 1 Mb region of chromosome 21 holding the features starting
    /// there, each regulatory feature with its type as `biotype` and an empty
    /// activity string, each motif naming its matrix, every feature with its
    /// input's slot as `source` and its place in the file as `ordinal`, and a
    /// table of 50 stubs without counts; `--extended-promoters` widens the
    /// promoters' spans.
    #[test]
    fn gff_inputs_write_the_shard_layout() {
        let fixtures = Path::new(FIXTURES);
        let inputs = RegulatoryInputs {
            regulatory_gff: Some(fixtures.join("regulatory_features.gff3.gz")),
            emars_gff: Some(fixtures.join("EMARs.gff.gz")),
            motif_gff: Some(fixtures.join("motif_features.gff3.gz")),
            extended_promoters: false,
            ..Default::default()
        };
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(&inputs, out.path(), &SequenceNames::default(), None).unwrap();
        assert_eq!(
            (
                stats.regulatory_features,
                stats.motif_features,
                stats.matrices
            ),
            (516 + 279, 2200, 50)
        );
        assert_eq!(stats.chromosomes.iter().collect::<Vec<_>>(), ["21"]);
        let mut shards: Vec<String> = std::fs::read_dir(out.path().join("regulatory/21"))
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        shards.sort();
        assert_eq!(
            shards,
            [
                "10000001-11000000.json",
                "5000001-6000000.json",
                "6000001-7000000.json",
                "7000001-8000000.json",
                "8000001-9000000.json",
                "9000001-10000000.json"
            ]
        );
        assert_eq!(stats.shards, 6);
        let first = shard(out.path(), "21", "5000001-6000000.json");
        assert!(first.iter().all(|feature| {
            let start = feature["start"].as_u64().unwrap();
            (5_000_001..=6_000_000).contains(&start)
        }));
        let enhancer = first
            .iter()
            .find(|feature| feature["stable_id"] == "ENSR21_85ZF5")
            .unwrap();
        assert_eq!(
            enhancer,
            &serde_json::json!({
                "feature_type": "RegulatoryFeature",
                "stable_id": "ENSR21_85ZF5",
                "start": 5011320,
                "end": 5011550,
                "biotype": "enhancer",
                "activity": "",
                "source": 0,
                "ordinal": 0
            })
        );
        let emar = first
            .iter()
            .find(|feature| feature["stable_id"] == "ENSR00001793149")
            .unwrap();
        assert_eq!(
            (emar["biotype"].as_str(), emar["source"].as_u64()),
            (Some("EMAR"), Some(1))
        );
        let motif = first
            .iter()
            .find(|feature| feature["stable_id"] == "ENSM00000574652")
            .unwrap();
        assert_eq!(
            motif,
            &serde_json::json!({
                "feature_type": "MotifFeature",
                "stable_id": "ENSM00000574652",
                "start": 5010250,
                "end": 5010263,
                "strand": -1,
                "matrix": "ENSPFM0085",
                "source": 2,
                "ordinal": 0
            })
        );
        let table: BTreeMap<String, serde_json::Value> = serde_json::from_slice(
            &std::fs::read(out.path().join("regulatory/binding_matrices.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(table.len(), 50);
        assert_eq!(
            table["ENSPFM0085"],
            serde_json::json!({"length": 14, "transcription_factors": ["EBF1"]})
        );
        assert_eq!(
            table["ENSPFM0290"]["transcription_factors"],
            serde_json::json!(["HNF4A", "NR2F1", "RXRA"])
        );
        assert!(table
            .values()
            .all(|entry| entry.get("frequencies").is_none()));

        let widened = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                extended_promoters: true,
                ..inputs
            },
            widened.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap();
        assert_eq!(stats.regulatory_features, 516 + 279);
        let promoter = |dir: &Path| {
            shard(dir, "21", "5000001-6000000.json")
                .into_iter()
                .find(|feature| feature["stable_id"] == "ENSR21_B34QG")
                .map(|feature| {
                    (
                        feature["start"].as_u64().unwrap(),
                        feature["end"].as_u64().unwrap(),
                    )
                })
                .unwrap()
        };
        assert_eq!(promoter(out.path()), (5_553_169, 5_553_746));
        assert_eq!(promoter(widened.path()), (5_553_154, 5_553_746));
    }

    /// A record of an unsupported type is skipped, a `chr`-prefixed name is
    /// written under the cache's spelling, a record without `ID` or `Name` is
    /// named by the file's own spelling and its coordinates, as the run-time
    /// reader names it, a duplicate record is written once while two distinct
    /// records sharing an `ID` are both written in file order (the runtime
    /// collapses them per variant), an `ID` holding a second `=` is cut
    /// there, and a motif naming a matrix with another span refuses the
    /// build.
    #[test]
    fn records_read_as_the_runtime_reads_them() {
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_gff(
            dir.path(),
            "reg.gff3.gz",
            &[
                "##gff-version 3",
                "chr21\tEnsembl\texon\t500\t600\t.\t.\t.\tID=skipped",
                "chr21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tcolor=#faca00",
                "chr21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tcolor=#faca00",
                "chr21\tEnsembl\tCTCF_binding_site\t2000\t2020\t.\t-\t.\tName=CTCF_1",
                "chr21\tEnsembl\tenhancer\t4000\t4200\t.\t.\t.\tID=TWICE",
                "chr21\tEnsembl\tpromoter\t4150\t4300\t.\t.\t.\tID=TWICE",
                "chr21\tEnsembl\tenhancer\t5000\t5200\t.\t.\t.\tID=CUT=HERE",
            ],
        );
        let motifs = write_gff(
            dir.path(),
            "motifs.gff3.gz",
            &[
                "chr21\tEnsembl\tTF_binding_site\t3000\t3013\t5.1\t+\t.\tID=M1;binding_matrix_id=PFM;transcription_factor=A,B",
                "chr21\tEnsembl\tTF_binding_site\t3100\t3113\t5.1\t-\t.\tID=M2;binding_matrix_id=PFM;transcription_factor=A,B",
            ],
        );
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(regulatory),
                emars_gff: None,
                motif_gff: Some(motifs),
                extended_promoters: false,
                ..Default::default()
            },
            out.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap();
        assert_eq!(
            (
                stats.regulatory_features,
                stats.motif_features,
                stats.shards
            ),
            (5, 2, 1)
        );
        let ids: Vec<String> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| feature["stable_id"].as_str().unwrap().to_string())
            .collect();
        // `ID=CUT=HERE` names the feature `CUT`, the text between the first
        // and the second `=`, as the run-time reader and GFF3.pm:194 cut it.
        assert_eq!(
            ids,
            [
                "chr21:1000-1200:enhancer",
                "CTCF_1",
                "TWICE",
                "TWICE",
                "CUT",
                "M1",
                "M2"
            ]
        );

        // The duplicate check keys on the start column the file is sorted by,
        // so a promoter widened under --extended-promoters still meets its copy.
        let widened = write_gff(
            dir.path(),
            "widened.gff3.gz",
            &[
                "21\tEnsembl\tpromoter\t1000\t1200\t.\t.\t.\tID=P;extended_start=900",
                "21\tEnsembl\tenhancer\t1000\t1100\t.\t.\t.\tID=E",
                "21\tEnsembl\tpromoter\t1000\t1200\t.\t.\t.\tID=P;extended_start=900",
            ],
        );
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(widened),
                extended_promoters: true,
                ..Default::default()
            },
            out.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap();
        assert_eq!(stats.regulatory_features, 2);
        let spans: Vec<(String, u64)> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| {
                (
                    feature["stable_id"].as_str().unwrap().to_string(),
                    feature["start"].as_u64().unwrap(),
                )
            })
            .collect();
        assert_eq!(spans, [("P".to_string(), 900), ("E".to_string(), 1000)]);

        let conflicting = write_gff(
            dir.path(),
            "conflicting.gff3.gz",
            &[
                "21\tEnsembl\tTF_binding_site\t3000\t3013\t5.1\t+\t.\tID=M1;binding_matrix_id=PFM;transcription_factor=A",
                "21\tEnsembl\tTF_binding_site\t3100\t3114\t5.1\t-\t.\tID=M2;binding_matrix_id=PFM;transcription_factor=A",
            ],
        );
        let refused = tempfile::tempdir().unwrap();
        let error = build_regulatory(
            &RegulatoryInputs {
                motif_gff: Some(conflicting),
                ..Default::default()
            },
            refused.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap_err();
        assert!(
            format!("{error:#}")
                .contains("binding matrix PFM spans 15 bases with factors [\"A\"] at line 2"),
            "{error:#}"
        );
    }

    /// The type column is read as `get_type` reads it (`GFF3.pm:161-164`,
    /// `uri_unescape`), so a `%XX` escape in it is decoded before the
    /// include-list check (`RegFeat.pm:312-314`): `CTCF%5Fbinding%5Fsite` is
    /// a `CTCF_binding_site` feature with that biotype, `TF%5Fbinding%5Fsite`
    /// a motif, a `promot%65r` is widened and named by the decoded type, and
    /// a type unsupported once decoded is skipped.
    #[test]
    fn an_escaped_type_column_is_decoded_before_the_include_list_check() {
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_gff(
            dir.path(),
            "escaped.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=PLAIN_TYPE",
                "21\tEnsembl\tex%6Fn\t2000\t2100\t.\t.\t.\tID=exon1",
                "21\tEnsembl\tCTCF%5Fbinding%5Fsite\t3000\t3020\t.\t+\t.\tID=ESCAPED_TYPE",
                "21\tEnsembl\tpromot%65r\t4000\t4200\t.\t.\t.\textended_start=3900",
                "21\tEnsembl\tpromoter%5Fflanking%5Fregion\t5000\t5100\t.\t.\t.\tcolor=#faca00",
            ],
        );
        let motifs = write_gff(
            dir.path(),
            "escaped_motifs.gff3.gz",
            &["21\tEnsembl\tTF%5Fbinding%5Fsite\t6000\t6013\t.\t-\t.\tID=ESCAPED_MOTIF;binding_matrix_id=M;transcription_factor=TF"],
        );
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(regulatory),
                emars_gff: None,
                motif_gff: Some(motifs),
                extended_promoters: true,
                ..Default::default()
            },
            out.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap();
        assert_eq!((stats.regulatory_features, stats.motif_features), (4, 1));
        let decoded: Vec<(String, u64, Option<String>)> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| {
                (
                    feature["stable_id"].as_str().unwrap().to_string(),
                    feature["start"].as_u64().unwrap(),
                    feature["biotype"].as_str().map(str::to_string),
                )
            })
            .collect();
        assert_eq!(
            decoded,
            [
                ("PLAIN_TYPE".to_string(), 1000, Some("enhancer".to_string())),
                (
                    "ESCAPED_TYPE".to_string(),
                    3000,
                    Some("CTCF_binding_site".to_string())
                ),
                (
                    "21:3900-4200:promoter".to_string(),
                    3900,
                    Some("promoter".to_string())
                ),
                (
                    "21:5000-5100:promoter_flanking_region".to_string(),
                    5000,
                    Some("promoter_flanking_region".to_string())
                ),
                ("ESCAPED_MOTIF".to_string(), 6000, None),
            ]
        );
    }

    /// A type column that is not UTF-8, whether an escape spells the byte
    /// (`caf%FF`) or it stands in the column, is a type like any other once
    /// decoded: it matches no accepted type, so the record is skipped and the
    /// build goes on, as Ensembl VEP skips the record with its warning
    /// (`RegFeat.pm:314-322`) and as the run-time reader does; an escape
    /// decoding to text that is UTF-8 is skipped the same way, an escape that
    /// completes a raw first byte included. A record of an accepted type
    /// holding such a byte raw refuses the build, naming the line.
    #[test]
    fn a_type_column_that_is_not_utf8_is_skipped() {
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_gff(
            dir.path(),
            "escaped.gff3.gz",
            &[
                b"21\tEnsembl\tcaf%FF\t1000\t1100\t.\t.\t.\tID=ESCAPED".as_slice(),
                b"21\tEnsembl\tcaf\xFF\t1500\t1600\t.\t.\t.\tID=RAW",
                b"21\tEnsembl\tcaf%CA%80\t2000\t2100\t.\t.\t.\tID=UNSUPPORTED",
                b"21\tEnsembl\tmix\xCA%80\t2500\t2600\t.\t.\t.\tID=MIXED",
                b"21\tEnsembl\tenhancer\t3000\t3200\t.\t.\t.\tID=KEPT",
            ],
        );
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(regulatory),
                ..Default::default()
            },
            out.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap();
        assert_eq!(stats.regulatory_features, 1);
        let ids: Vec<String> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| feature["stable_id"].as_str().unwrap().to_string())
            .collect();
        assert_eq!(ids, ["KEPT"]);

        let attribute = write_gff(
            dir.path(),
            "attribute.gff3.gz",
            &[b"21\tEnsembl\tenhancer\t1000\t1100\t.\t.\t.\tID=caf\xCA%80".as_slice()],
        );
        let out = tempfile::tempdir().unwrap();
        let error = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(attribute),
                ..Default::default()
            },
            out.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap_err();
        let message = format!("{error:#}");
        assert!(
            message.contains(&format!(
                "line 1 of {} is not UTF-8",
                dir.path().join("attribute.gff3.gz").display()
            )),
            "{message}"
        );
    }

    /// A byte that is not UTF-8 in a line of a type the kind skips refuses
    /// nothing, whichever column holds it: the line is skipped for its type
    /// and the bake writes the features the run-time reader builds from the
    /// file, where Ensembl VEP skips the line with its warning and carries
    /// the byte (`RegFeat.pm:314-322`). The line's start and end are still
    /// held to the index build's rule, a start a number leads passing with
    /// such a byte behind its digits and one no number leads refusing the
    /// build with the replacement character standing for the byte.
    #[test]
    fn a_byte_that_is_not_utf8_in_a_skipped_line_refuses_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let bake = |name: &str, lines: &[&[u8]]| {
            let out = tempfile::tempdir().unwrap();
            build_regulatory(
                &RegulatoryInputs {
                    regulatory_gff: Some(write_gff(dir.path(), name, lines)),
                    ..Default::default()
                },
                out.path(),
                &SequenceNames::default(),
                None,
            )
            .map(|stats| stats.regulatory_features)
            .map_err(|error| format!("{error:#}"))
        };
        let kept = b"21\tEnsembl\tenhancer\t3000\t3200\t.\t.\t.\tID=KEPT".as_slice();
        assert_eq!(
            bake(
                "skipped.gff3.gz",
                &[
                    b"21\tEnsembl\texon\t1000\t1100\t.\t.\t.\tID=SKIPPED;note=caf\xFF",
                    b"21\tEns\xFFmbl\texon\t1200\t1300\t.\t.\t.\tID=SKIPPED",
                    b"2\xFF1\tEnsembl\texon\t1400\t1500\t.\t.\t.\tID=SKIPPED",
                    b"21\tEnsembl\texon\t1600\xFF\t1700\xFF\t.\xFF\t.\t.\tID=SKIPPED",
                    kept,
                ],
            ),
            Ok(1)
        );
        let message = bake(
            "refused.gff3.gz",
            &[
                kept,
                b"21\tEnsembl\texon\tab\xFFc\t3400\t.\t.\t.\tID=SKIPPED",
            ],
        )
        .unwrap_err();
        let refusal = format!(
            "the start of line 2 of {} is not a whole number: ab\u{FFFD}c",
            dir.path().join("refused.gff3.gz").display()
        );
        assert!(message.contains(&refusal), "{message}");
    }

    /// A line that opens with a tab, or with white space and a tab, is read
    /// at the columns the file writes, where the index build reads them: its
    /// type stands where a start is read, so the line refuses the bake by
    /// that type, the record's own type one the kind keeps or not, as it
    /// refuses the run-time reader's read of the file. A line so led whose
    /// type a number leads is skipped for the type its third column spells.
    #[test]
    fn a_line_that_opens_with_a_tab_is_refused_for_the_type_where_a_start_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let bake = |led: &str| {
            let regulatory = write_gff(
                dir.path(),
                "t.gff3.gz",
                &["21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=A21", led],
            );
            let out = tempfile::tempdir().unwrap();
            build_regulatory(
                &RegulatoryInputs {
                    regulatory_gff: Some(regulatory),
                    ..Default::default()
                },
                out.path(),
                &SequenceNames::default(),
                None,
            )
            .map(|stats| stats.regulatory_features)
            .map_err(|error| format!("{error:#}"))
        };
        let file = dir.path().join("t.gff3.gz").display().to_string();
        for (lead, feature_type) in [("\t", "enhancer"), ("\t", "exon"), (" \t", "enhancer")] {
            let led = format!("{lead}21\tEnsembl\t{feature_type}\t2000\t2100\t.\t.\t.\tID=LED");
            let message = bake(&led).unwrap_err();
            let refusal =
                format!("the start of line 2 of {file} is not a whole number: {feature_type}");
            assert!(message.contains(&refusal), "{message}");
        }
        assert_eq!(
            bake("\t21\tEnsembl\t5_region\t2000\t2100\t.\t.\t.\tID=LED"),
            Ok(1)
        );
    }

    /// Every `%XX` escape of a read attribute value is decoded, in either
    /// letter case and whether or not GFF3 reserves the character, after the
    /// tests on the value as written: `ID=ESC%3bAPE%2c1%3E2` names the feature
    /// `ESC;APE,1>2`, an empty `ID` falls to `Name=a%2Gb%2`, which stays
    /// written, a motif's `binding_matrix_id=ENSPFM%30` is the matrix
    /// `ENSPFM0` and its `transcription_factor=A%2cB,C` is split on the comma
    /// written, so the table lists `A,B` and `C`, the shape the run-time reader
    /// builds. A record without `ID` or `Name` is named by the file's own
    /// spelling of the chromosome as written, `chr%41:3000-3200:enhancer`,
    /// since the run-time reader decodes no coordinate-derived identifier, and
    /// is filed under that spelling whole: no digit or capital follows its
    /// `chr`, so it is no UCSC prefix. A
    /// value whose escapes decode to bytes that are not UTF-8 refuses the
    /// build, naming the line and the file, on a record of an accepted type
    /// alone: the skipped `exon` line's `ID=caf%FF` refuses nothing.
    #[test]
    fn every_escape_of_a_read_attribute_value_is_decoded_after_the_written_tests() {
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_gff(
            dir.path(),
            "escaped.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ESC%3bAPE%2c1%3E2",
                "21\tEnsembl\tpromoter\t2000\t2200\t.\t.\t.\tID=;Name=a%2Gb%2",
                "chr%41\tEnsembl\tenhancer\t3000\t3200\t.\t.\t.\tcolor=#faca00",
                "21\tEnsembl\texon\t4000\t4100\t.\t.\t.\tID=caf%FF",
            ],
        );
        let motifs = write_gff(
            dir.path(),
            "motifs.gff3.gz",
            &["21\tEnsembl\tTF_binding_site\t5000\t5013\t.\t+\t.\tID=m;binding_matrix_id=ENSPFM%30;transcription_factor=A%2cB,C"],
        );
        let out = tempfile::tempdir().unwrap();
        build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(regulatory),
                motif_gff: Some(motifs),
                ..Default::default()
            },
            out.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap();
        let ids: Vec<String> = shard(out.path(), "21", "1-1000000.json")
            .iter()
            .map(|feature| feature["stable_id"].as_str().unwrap().to_string())
            .collect();
        assert_eq!(ids, ["ESC;APE,1>2", "a%2Gb%2", "m"]);
        assert_eq!(
            shard(out.path(), "chr%41", "1-1000000.json")[0]["stable_id"],
            "chr%41:3000-3200:enhancer"
        );
        let table: serde_json::Value = serde_json::from_slice(
            &std::fs::read(out.path().join("regulatory/binding_matrices.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(
            table["ENSPFM0"],
            serde_json::json!({"length": 14, "transcription_factors": ["A,B", "C"]})
        );

        let refused = write_gff(
            dir.path(),
            "refused.gff3.gz",
            &["21\tEnsembl\tenhancer\t1000\t1100\t.\t.\t.\tID=caf%FF"],
        );
        let out = tempfile::tempdir().unwrap();
        let error = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(refused),
                ..Default::default()
            },
            out.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap_err();
        let message = format!("{error:#}");
        assert!(
            message.ends_with(&format!(
                "line 1 of {}: GFF3 attribute ID \"caf%FF\" decodes to text that is not UTF-8",
                dir.path().join("refused.gff3.gz").display()
            )),
            "{message}"
        );
    }

    /// A byte-order mark at the start of a line changes nothing to the bake
    /// ([`crate::gff3::without_byte_order_marks`]): a file with a mark before
    /// each of its lines, its comment among them, writes the files the same
    /// lines write without the marks.
    #[test]
    fn a_marked_file_bakes_as_the_file_without_its_marks() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "##gff-version 3",
            "20\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=A20",
            "21\tEnsembl\tpromoter\t1000\t1200\t.\t.\t.\tID=A21",
        ];
        let marked = lines.map(|line| format!("\u{feff}{line}"));
        let bake = |name: &str, lines: &[String]| {
            let out = dir.path().join(name);
            let inputs = RegulatoryInputs {
                regulatory_gff: Some(write_gff(dir.path(), &format!("{name}.gff3.gz"), lines)),
                ..Default::default()
            };
            build_regulatory(&inputs, &out, &SequenceNames::default(), None).unwrap();
            [
                "20/1-1000000.json",
                "21/1-1000000.json",
                "binding_matrices.json",
            ]
            .map(|file| std::fs::read(out.join("regulatory").join(file)).unwrap())
        };
        assert_eq!(
            bake("marked", &marked),
            bake("plain", &lines.map(String::from))
        );
    }

    /// A start or end no number leads, or an end below zero, refuses the bake
    /// on a line of a type the kind skips, on a file's first sequence and on a
    /// later one, with the text a record of an accepted type draws, the rule
    /// the transcript readers follow
    /// ([`vep_core::regulatory_gff::unindexed_coordinate`]). Nothing else of a
    /// skipped line is read: a start of 0, an end before its start that is
    /// not below zero, a seqid that could name no cache directory and a start
    /// a number leads without being one, an empty start before its end among
    /// them, refuse nothing there, where each refuses an `enhancer` line.
    #[test]
    fn a_skipped_line_is_refused_for_a_coordinate_the_index_build_fails_alone() {
        let dir = tempfile::tempdir().unwrap();
        // Two sequences, each an enhancer and one more line: line 2 of the
        // file is on the first sequence and line 4 on the later one.
        let bake = |second: &str, fourth: &str| {
            let regulatory = write_gff(
                dir.path(),
                "t.gff3.gz",
                &[
                    "20\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=A20",
                    second,
                    "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=A21",
                    fourth,
                ],
            );
            let out = tempfile::tempdir().unwrap();
            build_regulatory(
                &RegulatoryInputs {
                    regulatory_gff: Some(regulatory),
                    ..Default::default()
                },
                out.path(),
                &SequenceNames::default(),
                None,
            )
            .map(|stats| stats.regulatory_features)
            .map_err(|error| format!("{error:#}"))
        };
        let line = |chr: &str, feature_type: &str, start: &str, end: &str| {
            format!("{chr}\tEnsembl\t{feature_type}\t{start}\t{end}\t.\t.\t.\tID=S{chr}")
        };
        let sound_first = line("20", "exon", "2000", "2100");
        let sound_later = line("21", "exon", "2000", "2100");
        let file = dir.path().join("t.gff3.gz").display().to_string();
        for (start, end, column, value) in [
            ("abc", "2100", "start", "abc"),
            ("2000", "xyz", "end", "xyz"),
            ("", "xyz", "start", ""),
            ("2000", "", "end", ""),
            ("2000", "-9", "end", "-9"),
        ] {
            let message = bake(&line("20", "exon", start, end), &sound_later).unwrap_err();
            let refusal =
                format!("the {column} of line 2 of {file} is not a whole number: {value}");
            assert!(message.contains(&refusal), "{message}");
            let message = bake(&sound_first, &line("21", "exon", start, end)).unwrap_err();
            let refusal =
                format!("the {column} of line 4 of {file} is not a whole number: {value}");
            assert!(message.contains(&refusal), "{message}");
        }
        for (chr, start, end) in [
            ("20", "0", "2100"),
            ("20", "2100", "2000"),
            ("20", "2000x", "2100"),
            ("20", "-5", "2100"),
            ("20", " 2000", "2100"),
            ("20", "0x7D0", "2100"),
            ("20", "", "2100"),
            ("20", "2000", "-0"),
            ("20/x", "2000", "2100"),
        ] {
            let on_later = chr.replace("20", "21");
            for (second, fourth) in [
                (line(chr, "exon", start, end), sound_later.clone()),
                (sound_first.clone(), line(&on_later, "exon", start, end)),
            ] {
                assert_eq!(bake(&second, &fourth), Ok(2), "{chr:?} {start:?} {end:?}");
            }
        }
        for (chr, start, end, refusal) in [
            (
                "21",
                "2000x",
                "2100",
                format!("the start of line 4 of {file} is not a whole number: 2000x"),
            ),
            (
                "21",
                "",
                "2100",
                format!(
                    "the start of line 4 of {file} is not a whole number: : cannot parse integer \
                     from empty string"
                ),
            ),
            (
                "21",
                "0",
                "2100",
                "GFF coordinates 0-2100 at line 4 are not 1-based and ordered".to_string(),
            ),
            (
                "21",
                "2100",
                "2000",
                "GFF coordinates 2100-2000 at line 4 are not 1-based and ordered".to_string(),
            ),
            (
                "21/x",
                "2000",
                "2100",
                "invalid chromosome name '21/x' at line 4".to_string(),
            ),
        ] {
            let message = bake(&sound_first, &line(chr, "enhancer", start, end)).unwrap_err();
            assert!(message.contains(&refusal), "{message}");
        }
    }

    /// A line whose sequence name a NUL byte ends is filed by the index build
    /// under the name before the NUL byte, and a run with `--regulatory_gff`
    /// finds no record of the chromosome in it. The bake skips such a line
    /// the index build places, of an accepted type or a skipped one, its
    /// columns NUL bytes end, and writes the features of the file without
    /// it. A line that keeps its sequence name with a NUL byte in place of
    /// one later tab is a record of fewer than nine columns, refused naming
    /// the count, as is a line of NUL-ended columns whose start no number
    /// leads, which the index build fails.
    #[test]
    fn a_line_whose_sequence_name_a_nul_byte_ends_is_skipped() {
        let dir = tempfile::tempdir().unwrap();
        let bake = |lines: &[&str]| {
            let regulatory = write_gff(dir.path(), "t.gff3.gz", lines);
            let out = tempfile::tempdir().unwrap();
            build_regulatory(
                &RegulatoryInputs {
                    regulatory_gff: Some(regulatory),
                    ..Default::default()
                },
                out.path(),
                &SequenceNames::default(),
                None,
            )
            .map(|_| std::fs::read(out.path().join("regulatory/21/1-1000000.json")).unwrap())
            .map_err(|error| format!("{error:#}"))
        };
        let first = "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=A";
        let last = "21\tEnsembl\tenhancer\t3000\t3200\t.\t.\t.\tID=B";
        let without = bake(&[first, last]).unwrap();
        let ids: Vec<String> = serde_json::from_slice::<Vec<serde_json::Value>>(&without)
            .unwrap()
            .iter()
            .map(|feature| feature["stable_id"].as_str().unwrap().to_string())
            .collect();
        assert_eq!(ids, ["A", "B"]);
        let enhancer = "21\0Ensembl\0enhancer\x002000\x002100\0.\0.\0.\0ID=LOST";
        let exon = "21\0Ensembl\0exon\x002500\x002600\0.\0.\0.\0ID=SKIPPED";
        assert!(bake(&[first, enhancer, exon, last]).unwrap() == without);

        let file = dir.path().join("t.gff3.gz").display().to_string();
        // A kept line with one NUL byte in its sequence name or its source:
        // the index build reads the type as the start.
        for (line, start) in [
            (
                "21\0x\tEnsembl\tenhancer\t2000\t2100\t.\t.\t.\tID=C",
                "enhancer",
            ),
            (
                "21\tEns\0embl\tenhancer\t2000\t2100\t.\t.\t.\tID=C",
                "enhancer",
            ),
            (
                "21\tEnsembl\0\tenhancer\t2000\t2100\t.\t.\t.\tID=C",
                "enhancer",
            ),
        ] {
            let message = bake(&[first, line, last]).unwrap_err();
            let refusal = format!("the start of line 2 of {file} is not a whole number: {start}");
            assert!(message.contains(&refusal), "{message}");
        }
        // Behind the end the columns the build reads are the record's own.
        let color = "21\tEnsembl\tenhancer\t2000\t2100\t.\t.\t.\tID=C;color=#fa\0ca00";
        assert!(bake(&[first, color, last]).is_ok());
        for (line, columns) in [
            ("21\tEnsembl\tenhancer\x002000\t2100\t.\t.\t.\tID=C", 8),
            ("21\tEnsembl\tenhancer\t2000\t2100\0.\t.\t.\tID=C", 8),
            ("21\0Ensembl\0enhancer\0abc\x002100\0.\0.\0.\0ID=C", 1),
            ("21\0Ensembl\0enhancer", 1),
        ] {
            let message = bake(&[first, line, last]).unwrap_err();
            let refusal =
                format!("line 2 of {file} has {columns} columns where a GFF record has 9");
            assert!(message.ends_with(&refusal), "{message}");
        }
    }

    /// A line the kind skips that holds a NUL byte and a byte that is not
    /// UTF-8 is held to the index build's reading of its start and its end
    /// and to nothing else. It is skipped where that build parses both, as a
    /// line of nine columns and as one whose sequence name the NUL byte
    /// ends; where the NUL byte leaves the type's tail for a start, the line
    /// refuses the bake by that text, the replacement character standing for
    /// the other byte.
    #[test]
    fn a_skipped_line_holding_a_nul_byte_and_a_byte_that_is_not_utf8_is_held_to_the_index_build() {
        let dir = tempfile::tempdir().unwrap();
        let bake = |lines: &[&[u8]]| {
            let out = tempfile::tempdir().unwrap();
            build_regulatory(
                &RegulatoryInputs {
                    regulatory_gff: Some(write_gff(dir.path(), "t.gff3.gz", lines)),
                    ..Default::default()
                },
                out.path(),
                &SequenceNames::default(),
                None,
            )
            .map(|stats| stats.regulatory_features)
            .map_err(|error| format!("{error:#}"))
        };
        let kept = b"21\tEnsembl\tenhancer\t3000\t3200\t.\t.\t.\tID=KEPT".as_slice();
        for skipped in [
            b"21\tEns\xFFmbl\texon\t1000\t1100\t.\t.\t.\tID=SKIPPED;note=a\0b".as_slice(),
            b"21\0Ensembl\0exon\x001000\x001100\0.\0.\0.\0ID=caf\xFF",
        ] {
            assert_eq!(bake(&[skipped, kept]), Ok(1));
        }
        let message = bake(&[
            b"21\tEnsembl\tex\0\xFFn\t1000\t1100\t.\t.\t.\tID=SKIPPED".as_slice(),
            kept,
        ])
        .unwrap_err();
        let refusal = format!(
            "the start of line 1 of {} is not a whole number: \u{FFFD}n",
            dir.path().join("t.gff3.gz").display()
        );
        assert!(message.contains(&refusal), "{message}");
    }

    fn write_text(dir: &Path, name: &str, text: &str) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, text).unwrap();
        path
    }

    /// The `activity` of each regulatory feature of the tree under `dir`, by
    /// its identifier, in shard order.
    fn activities(dir: &Path) -> Vec<(String, String)> {
        shard(dir, "21", "1-1000000.json")
            .iter()
            .filter(|feature| feature["feature_type"] == "RegulatoryFeature")
            .map(|feature| {
                (
                    feature["stable_id"].as_str().unwrap().to_string(),
                    feature["activity"].as_str().unwrap().to_string(),
                )
            })
            .collect()
    }

    fn matrices_table(dir: &Path) -> serde_json::Value {
        serde_json::from_slice(
            &std::fs::read(dir.join("regulatory/binding_matrices.json")).unwrap(),
        )
        .unwrap()
    }

    /// Three regulatory records and one motif on chromosome 21: two features
    /// the made-up activity files hold a row of, and a CTCF binding site they
    /// hold none of.
    fn activity_gffs(dir: &Path) -> RegulatoryInputs {
        RegulatoryInputs {
            regulatory_gff: Some(write_gff(
                dir,
                "regulatory.gff3.gz",
                &[
                    "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=R1",
                    "21\tEnsembl\tCTCF_binding_site\t2000\t2020\t.\t.\t.\tID=NO_ROW",
                    "21\tEnsembl\topen_chromatin_region\t3000\t3100\t.\t.\t.\tID=R2",
                ],
            )),
            motif_gff: Some(write_gff(
                dir,
                "motifs.gff3.gz",
                &["21\tEnsembl\tTF_binding_site\t5000\t5002\t.\t+\t.\tID=M1;binding_matrix_id=PFM_A;transcription_factor=GFF_A"],
            )),
            ..Default::default()
        }
    }

    /// An activity file gives every regulatory feature one code per cell type.
    /// The record's list is the header's names with each space written `_`,
    /// in byte order, whatever the header's order: a capital sorts before a
    /// small letter, and a name keeps its comma and its brackets. Each of the
    /// five activities has its code, a row's type cell is not read (empty on
    /// the second row), a feature the file holds no row of carries `-` for
    /// every name, and a row whose feature the GFF does not hold is passed
    /// over. A gzip-compressed file reads as the plain one.
    #[test]
    fn an_activity_file_gives_each_feature_one_code_per_cell_type() {
        let dir = tempfile::tempdir().unwrap();
        let text = "feature_id\tfeature_type\ttissue two (m, 4 y)\tcell line one\t\
                    tissue one (f, 30 y)\tTissue three\torgan four\n\
                    R1\tenhancer\tACTIVE\tINACTIVE\tPOISED\tREPRESSED\tNA\n\
                    R2\t\tINACTIVE\tINACTIVE\tACTIVE\tINACTIVE\tINACTIVE\n\
                    NOT_IN_THE_GFF\tpromoter\tACTIVE\tACTIVE\tACTIVE\tACTIVE\tACTIVE\n";
        let plain = write_text(dir.path(), "activity.tsv", text);
        let compressed = write_gff(dir.path(), "activity.tsv.gz", &[text.trim_end()]);
        for activity in [plain, compressed] {
            let out = tempfile::tempdir().unwrap();
            let stats = build_regulatory(
                &RegulatoryInputs {
                    activity: Some(activity),
                    ..activity_gffs(dir.path())
                },
                out.path(),
                &SequenceNames::default(),
                None,
            )
            .unwrap();
            assert_eq!(
                stats.cell_types.as_deref(),
                Some(
                    [
                        "Tissue_three",
                        "cell_line_one",
                        "organ_four",
                        "tissue_one_(f,_30_y)",
                        "tissue_two_(m,_4_y)"
                    ]
                    .map(String::from)
                    .as_slice()
                )
            );
            assert!(stats.activity_on_every_feature);
            assert_eq!(
                activities(out.path()),
                [
                    ("R1".to_string(), "RINPA".to_string()),
                    ("NO_ROW".to_string(), "-----".to_string()),
                    ("R2".to_string(), "IIIAI".to_string()),
                ]
            );
            assert_eq!((stats.regulatory_features, stats.motif_features), (3, 1));
        }
    }

    /// Two records of the GFF under one identifier take the one row of the
    /// activity file, and a tree baked with an activity file and no
    /// regulatory feature still lists the file's cell types.
    #[test]
    fn records_sharing_an_identifier_share_its_activity_row() {
        let dir = tempfile::tempdir().unwrap();
        let activity = write_text(
            dir.path(),
            "activity.tsv",
            "feature_id\tfeature_type\ttissue one\nTWICE\tenhancer\tPOISED\n",
        );
        let regulatory = write_gff(
            dir.path(),
            "regulatory.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t4000\t4200\t.\t.\t.\tID=TWICE",
                "21\tEnsembl\tpromoter\t4150\t4300\t.\t.\t.\tID=TWICE",
            ],
        );
        let out = tempfile::tempdir().unwrap();
        build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(regulatory),
                activity: Some(activity.clone()),
                ..Default::default()
            },
            out.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap();
        assert_eq!(
            activities(out.path()),
            [
                ("TWICE".to_string(), "P".to_string()),
                ("TWICE".to_string(), "P".to_string())
            ]
        );

        let motifs_alone = RegulatoryInputs {
            regulatory_gff: None,
            activity: Some(activity),
            ..activity_gffs(dir.path())
        };
        let out = tempfile::tempdir().unwrap();
        let stats =
            build_regulatory(&motifs_alone, out.path(), &SequenceNames::default(), None).unwrap();
        assert_eq!(
            stats.cell_types.as_deref(),
            Some(["tissue_one".to_string()].as_slice())
        );
        assert!(stats.activity_on_every_feature);
    }

    /// Each fault of an activity file refuses the build in one error that
    /// names the file, the line and what was read: an activity that is none
    /// of the five values, in another letter case or empty; a row of another
    /// number of columns than the header; a name the header holds twice once
    /// its spaces are underscores; a second row of one feature; a header that
    /// does not open with `feature_id` and `feature_type`; a file with no
    /// line.
    #[test]
    fn an_activity_file_that_does_not_read_as_one_refuses_the_build() {
        let dir = tempfile::tempdir().unwrap();
        let header = "feature_id\tfeature_type\ttissue two\ttissue one (f, 30 y)\n";
        let refusal = |text: &str| {
            let activity = write_text(dir.path(), "activity.tsv", text);
            let out = tempfile::tempdir().unwrap();
            let error = build_regulatory(
                &RegulatoryInputs {
                    activity: Some(activity),
                    ..activity_gffs(dir.path())
                },
                out.path(),
                &SequenceNames::default(),
                None,
            )
            .unwrap_err();
            format!("{error:#}")
        };
        let file = dir.path().join("activity.tsv").display().to_string();
        for (text, fault) in [
            (
                format!("{header}R1\tenhancer\tACTIVE\tINACTIVE\nR2\t\tINACTIVE\tActive\n"),
                format!(
                    "line 3 of {file} holds the activity \"Active\" of R2 in \
                     tissue_one_(f,_30_y), which is not ACTIVE, INACTIVE, POISED, REPRESSED or NA"
                ),
            ),
            (
                format!("{header}R1\tenhancer\t\tINACTIVE\n"),
                format!(
                    "line 2 of {file} holds the activity \"\" of R1 in tissue_two, which is not \
                     ACTIVE, INACTIVE, POISED, REPRESSED or NA"
                ),
            ),
            (
                format!("{header}R1\tenhancer\tACTIVE\tINACTIVE\nR2\t\tINACTIVE\n"),
                format!("line 3 of {file} has 3 columns where the header has 4"),
            ),
            (
                format!("{header}R1\tenhancer\tACTIVE\tINACTIVE\tACTIVE\n"),
                format!("line 2 of {file} has 5 columns where the header has 4"),
            ),
            (
                "feature_id\tfeature_type\ttissue one\ttissue two\ttissue_one\n".to_string(),
                format!(
                    "line 1 of {file} names the cell type \"tissue_one\" twice, each space read \
                     as an underscore"
                ),
            ),
            (
                format!("{header}R1\tenhancer\tACTIVE\tINACTIVE\n\nR1\tenhancer\tACTIVE\tACTIVE\n"),
                format!("line 4 of {file} is a second row of the feature R1"),
            ),
            (
                "id\tfeature_type\ttissue one\nR1\tenhancer\tACTIVE\n".to_string(),
                format!(
                    "line 1 of {file} opens with the columns [\"id\", \"feature_type\"] where \
                     the header of an activity file opens with feature_id and feature_type"
                ),
            ),
            (String::new(), format!("{file} holds no header line")),
        ] {
            let message = refusal(&text);
            assert_eq!(
                message,
                format!("Failed to read {file}: {fault}"),
                "{text:?}"
            );
        }
    }

    /// The made-up funcgen tables of the matrix tests: `PFM_A` of three
    /// positions and `PFM_B` of two, whose threshold is NULL, both named by a
    /// motif, and a matrix no motif names. `frequencies` holds the rows of
    /// `binding_matrix_frequencies`, the counts of `PFM_A` when `None`.
    fn matrix_tables(dir: &Path, frequencies: Option<&str>) {
        write_text(
            dir,
            "binding_matrix.txt",
            "10\tname of A\t4.4\tSELEX\tPFM_A\n\
             11\tname of B\t\\N\tSELEX\tPFM_B\n\
             12\tname of no motif\t1.5\tSELEX\tPFM_UNNAMED\n",
        );
        // As the dumps order a matrix's rows: every position of A, then of C,
        // G and T.
        let whole = "1\t10\t1\tA\t100\n2\t10\t2\tA\t0\n3\t10\t3\tA\t7\n\
                     4\t10\t1\tC\t0\n5\t10\t2\tC\t50\n6\t10\t3\tC\t7\n\
                     7\t10\t1\tG\t0\n8\t10\t2\tG\t50\n9\t10\t3\tG\t7\n\
                     10\t10\t1\tT\t3\n11\t10\t2\tT\t0\n12\t10\t3\tT\t115596\n";
        let others = "13\t11\t1\tA\t1\n14\t11\t2\tA\t2\n15\t11\t1\tC\t3\n16\t11\t2\tC\t4\n\
                      17\t11\t1\tG\t5\n18\t11\t2\tG\t6\n19\t11\t1\tT\t7\n20\t11\t2\tT\t8\n\
                      21\t12\t5\tA\t9\n";
        write_text(
            dir,
            "binding_matrix_frequencies.txt",
            &format!("{}{others}", frequencies.unwrap_or(whole)),
        );
    }

    /// Motifs naming `PFM_A` (three bases), `PFM_B` (two) and, with
    /// `unlisted`, `PFM_C` (four), a matrix the made-up tables lack.
    fn matrix_motifs(dir: &Path, unlisted: bool) -> RegulatoryInputs {
        let mut lines = vec![
            "21\tEnsembl\tTF_binding_site\t5000\t5002\t.\t+\t.\tID=M1;binding_matrix_id=PFM_A;transcription_factor=GFF_A",
            "21\tEnsembl\tTF_binding_site\t5010\t5012\t.\t-\t.\tID=M2;binding_matrix_id=PFM_A;transcription_factor=GFF_A",
            "21\tEnsembl\tTF_binding_site\t6000\t6001\t.\t+\t.\tID=M3;binding_matrix_id=PFM_B;transcription_factor=GFF_B,GFF_B2",
        ];
        if unlisted {
            lines.push("21\tEnsembl\tTF_binding_site\t7000\t7003\t.\t+\t.\tID=M4;binding_matrix_id=PFM_C;transcription_factor=GFF_C");
        }
        RegulatoryInputs {
            motif_gff: Some(write_gff(dir, "motifs.gff3.gz", &lines)),
            ..Default::default()
        }
    }

    /// The funcgen tables give each matrix a motif names the entry a
    /// converted cache holds: one `[A, C, G, T]` row a position, position 1
    /// first, the table's integers; `length`, the number of positions; the
    /// name and source; the threshold, absent where the table holds NULL. A
    /// matrix no motif names is not written. Without the complex tables the
    /// factors are the motif GFF's; with them they are the display names of
    /// the matrix's complexes in the order of the complexes' identifiers,
    /// whatever the links' order, a link to a complex the table lacks naming
    /// nothing. A matrix the table lacks keeps its stub, and the tree is then
    /// one whose matrices do not all carry counts.
    #[test]
    fn the_funcgen_tables_give_each_matrix_its_counts() {
        let dir = tempfile::tempdir().unwrap();
        let dumps = dir.path().join("funcgen");
        std::fs::create_dir(&dumps).unwrap();
        matrix_tables(&dumps, None);
        let bake = |unlisted: bool| {
            let out = tempfile::tempdir().unwrap();
            let stats = build_regulatory(
                &RegulatoryInputs {
                    funcgen_dumps: Some(dumps.clone()),
                    ..matrix_motifs(dir.path(), unlisted)
                },
                out.path(),
                &SequenceNames::default(),
                None,
            )
            .unwrap();
            (matrices_table(out.path()), stats)
        };

        let (table, stats) = bake(false);
        assert_eq!(
            table,
            serde_json::json!({
                "PFM_A": {
                    "frequencies": [[100, 0, 0, 3], [0, 50, 50, 0], [7, 7, 7, 115596]],
                    "length": 3,
                    "name": "name of A",
                    "source": "SELEX",
                    "threshold": 4.4,
                    "transcription_factors": ["GFF_A"]
                },
                "PFM_B": {
                    "frequencies": [[1, 3, 5, 7], [2, 4, 6, 8]],
                    "length": 2,
                    "name": "name of B",
                    "source": "SELEX",
                    "transcription_factors": ["GFF_B", "GFF_B2"]
                }
            })
        );
        assert_eq!(stats.matrices, 2);
        assert!(stats.counts_on_every_matrix);

        let (table, stats) = bake(true);
        assert_eq!(
            table["PFM_C"],
            serde_json::json!({"length": 4, "transcription_factors": ["GFF_C"]})
        );
        assert_eq!(table["PFM_A"]["length"], 3);
        assert_eq!(stats.matrices, 3);
        assert!(!stats.counts_on_every_matrix);

        write_text(
            &dumps,
            "binding_matrix_transcription_factor_complex.txt",
            "1\t10\t7\n2\t10\t3\n3\t11\t99\n4\t12\t3\n",
        );
        write_text(
            &dumps,
            "transcription_factor_complex.txt",
            "3\tthird_factor\tThird::Factor\n7\tseventh\tSeventh\n8\tof_no_matrix\tEighth\n",
        );
        let (table, _) = bake(true);
        assert_eq!(
            table["PFM_A"]["transcription_factors"],
            serde_json::json!(["Third::Factor", "Seventh"])
        );
        assert_eq!(
            table["PFM_B"]["transcription_factors"],
            serde_json::json!([])
        );
        assert_eq!(
            table["PFM_C"]["transcription_factors"],
            serde_json::json!(["GFF_C"])
        );

        // One table of the pair is no pair: the factors are the GFF's.
        std::fs::remove_file(dumps.join("transcription_factor_complex.txt")).unwrap();
        let (table, _) = bake(false);
        assert_eq!(
            table["PFM_A"]["transcription_factors"],
            serde_json::json!(["GFF_A"])
        );
        // And without the counts' table every matrix is its stub.
        std::fs::remove_file(dumps.join("binding_matrix_frequencies.txt")).unwrap();
        let (table, stats) = bake(false);
        assert_eq!(
            table["PFM_A"],
            serde_json::json!({"length": 3, "transcription_factors": ["GFF_A"]})
        );
        assert!(!stats.counts_on_every_matrix);
    }

    /// Tables that do not hold a matrix whole refuse the build, naming the
    /// matrix: positions with a gap, a count stated twice, a position without
    /// one of the four nucleotides, a nucleotide that is none of the four, a
    /// matrix of another number of positions than its motifs span, one with
    /// no count at all, a count that is no number, a threshold that is no
    /// number, and a stable identifier two rows carry.
    #[test]
    fn funcgen_tables_that_do_not_hold_a_matrix_whole_refuse_the_build() {
        let dir = tempfile::tempdir().unwrap();
        let dumps = dir.path().join("funcgen");
        std::fs::create_dir(&dumps).unwrap();
        let refusal = |dumps: &Path| {
            let out = tempfile::tempdir().unwrap();
            let error = build_regulatory(
                &RegulatoryInputs {
                    funcgen_dumps: Some(dumps.to_path_buf()),
                    ..matrix_motifs(dir.path(), false)
                },
                out.path(),
                &SequenceNames::default(),
                None,
            )
            .unwrap_err();
            format!("{error:#}")
        };
        let position = |at: u32, counts: [&str; 4]| -> String {
            ["A", "C", "G", "T"]
                .iter()
                .zip(counts)
                .map(|(nucleotide, count)| format!("1\t10\t{at}\t{nucleotide}\t{count}\n"))
                .collect()
        };
        let sound = [1, 2, 3].map(|at| position(at, ["1", "2", "3", "4"]));
        for (frequencies, fault) in [
            (
                format!(
                    "{}{}{}",
                    sound[0],
                    sound[1],
                    position(4, ["1", "2", "3", "4"])
                ),
                "the positions of binding matrix PFM_A do not run from 1 without a gap: \
                 position 4 stands where position 3 would",
            ),
            (
                format!(
                    "{}{}{}",
                    position(0, ["1", "2", "3", "4"]),
                    sound[0],
                    sound[1]
                ),
                "the positions of binding matrix PFM_A do not run from 1 without a gap: \
                 position 0 stands where position 1 would",
            ),
            (
                format!("{}{}{}1\t10\t2\tG\t9\n", sound[0], sound[1], sound[2]),
                "binding matrix PFM_A holds the count of G at position 2 twice",
            ),
            (
                format!(
                    "{}{}1\t10\t3\tA\t1\n1\t10\t3\tC\t1\n1\t10\t3\tG\t1\n",
                    sound[0], sound[1]
                ),
                "binding matrix PFM_A holds no count of T at position 3",
            ),
            (
                format!("{}{}{}1\t10\t3\tN\t1\n", sound[0], sound[1], sound[2]),
                "binding matrix PFM_A holds a count of the nucleotide \"N\" at position 3; the \
                 nucleotides are A, C, G and T",
            ),
            (
                format!("{}{}", sound[0], sound[1]),
                "binding matrix PFM_A holds 2 positions where its motif features span 3 bases",
            ),
            (
                String::new(),
                "binding matrix PFM_A holds 0 positions where its motif features span 3 bases",
            ),
            (
                format!(
                    "{}{}{}",
                    sound[0],
                    sound[1],
                    position(3, ["1", "2", "many", "4"])
                ),
                "binding_matrix_frequencies row 1 of binding matrix PFM_A: \
                 binding_matrix_frequencies column 4 holds \"many\", not a number: invalid digit \
                 found in string",
            ),
        ] {
            matrix_tables(&dumps, Some(&frequencies));
            assert_eq!(
                refusal(&dumps),
                format!("Failed to read the binding matrices: {fault}"),
                "{frequencies:?}"
            );
        }

        matrix_tables(&dumps, None);
        for (rows, fault) in [
            (
                "10\tname of A\thigh\tSELEX\tPFM_A\n",
                "the threshold of binding matrix PFM_A is no number: high",
            ),
            (
                "10\tname of A\t4.4\tSELEX\tPFM_A\n13\tanother\t4.4\tSELEX\tPFM_A\n",
                "binding_matrix holds two rows of the binding matrix PFM_A, or two under its \
                 identifier 13",
            ),
        ] {
            write_text(&dumps, "binding_matrix.txt", rows);
            assert_eq!(
                refusal(&dumps),
                format!("Failed to read the binding matrices: {fault}"),
                "{rows:?}"
            );
        }
    }

    /// A genome of one record, `21`, of 24 bases on one line, the last eight
    /// in lower case.
    fn small_genome(dir: &Path) -> vep_fasta::IndexedFasta {
        let bases = "AACCGGTTACGTGGCAttgacgta";
        let fasta = write_text(dir, "genome.fa", &format!(">21\n{bases}\n"));
        write_text(dir, "genome.fa.fai", "21\t24\t4\t24\t25\n");
        vep_fasta::IndexedFasta::from_path(fasta).unwrap()
    }

    /// With a genome FASTA each motif carries the bases of its span in upper
    /// case: as the FASTA reads them on the forward strand and for a motif
    /// with no strand, their reverse complement on the reverse strand. A
    /// motif whose span runs past its record's end, or whose sequence the
    /// FASTA does not hold, is written without a sequence, and the tree is
    /// then one whose motifs do not all carry one. Without a genome FASTA no
    /// motif carries a sequence.
    #[test]
    fn a_motif_carries_the_bases_of_its_span_on_its_own_strand() {
        let dir = tempfile::tempdir().unwrap();
        let genome = small_genome(dir.path());
        let motif = |chr: &str, start: u32, end: u32, strand: &str, id: &str| {
            format!(
                "{chr}\tEnsembl\tTF_binding_site\t{start}\t{end}\t.\t{strand}\t.\t\
                 ID={id};binding_matrix_id=PFM{};transcription_factor=TF",
                end - start + 1
            )
        };
        let held = [
            motif("21", 3, 8, "+", "FORWARD"),
            motif("21", 9, 14, "-", "REVERSE"),
            motif("21", 17, 22, "+", "SMALL_LETTERS"),
            motif("chr21", 19, 24, "-", "SMALL_LETTERS_REVERSE"),
            motif("21", 1, 4, ".", "NO_STRAND"),
        ];
        let not_held = [
            motif("21", 22, 27, "+", "PAST_THE_END"),
            motif("22", 3, 8, "+", "NO_RECORD"),
        ];
        let bake = |lines: &[String], genome: Option<&vep_fasta::IndexedFasta>| {
            let out = tempfile::tempdir().unwrap();
            let stats = build_regulatory(
                &RegulatoryInputs {
                    motif_gff: Some(write_gff(dir.path(), "motifs.gff3.gz", lines)),
                    ..Default::default()
                },
                out.path(),
                &SequenceNames::default(),
                genome,
            )
            .unwrap();
            let sequences: Vec<(String, Option<String>)> = ["21", "22"]
                .iter()
                .filter(|chr| out.path().join("regulatory").join(chr).is_dir())
                .flat_map(|chr| shard(out.path(), chr, "1-1000000.json"))
                .map(|feature| {
                    (
                        feature["stable_id"].as_str().unwrap().to_string(),
                        feature
                            .get("seq")
                            .map(|seq| seq.as_str().unwrap().to_string()),
                    )
                })
                .collect();
            (sequences, stats.sequence_on_every_motif)
        };
        let named = |pairs: &[(&str, Option<&str>)]| -> Vec<(String, Option<String>)> {
            pairs
                .iter()
                .map(|(id, seq)| (id.to_string(), seq.map(str::to_string)))
                .collect()
        };

        let (sequences, every_motif) = bake(&held, Some(&genome));
        assert_eq!(
            sequences,
            named(&[
                ("FORWARD", Some("CCGGTT")),
                ("REVERSE", Some("CCACGT")),
                ("SMALL_LETTERS", Some("TTGACG")),
                ("SMALL_LETTERS_REVERSE", Some("TACGTC")),
                ("NO_STRAND", Some("AACC")),
            ])
        );
        assert!(every_motif);

        let mut lines = held[..1].to_vec();
        lines.extend(not_held);
        let (sequences, every_motif) = bake(&lines, Some(&genome));
        assert_eq!(
            sequences,
            named(&[
                ("FORWARD", Some("CCGGTT")),
                ("PAST_THE_END", None),
                ("NO_RECORD", None),
            ])
        );
        assert!(!every_motif);

        let (sequences, every_motif) = bake(&held[..2], None);
        assert_eq!(sequences, named(&[("FORWARD", None), ("REVERSE", None)]));
        assert!(!every_motif);
    }

    /// The cut of a motif's sequence equals the sequence Ensembl VEP's cache
    /// holds for it: every motif of the converted release 116 GRCh38 cache
    /// slice of the committed regulatory corpus whose span the corpus's
    /// reference holds as bases (the reference masks the rest of the
    /// chromosome) is cut from that reference as the cache has it, on both
    /// strands.
    #[test]
    fn motif_sequences_equal_those_of_ensembl_veps_cache() {
        let corpus =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-regulatory");
        let genome = vep_fasta::IndexedFasta::from_path(corpus.join("reference.fa.gz")).unwrap();
        let (mut forward, mut reverse, mut masked) = (0, 0, 0);
        let mut shards: Vec<PathBuf> = std::fs::read_dir(corpus.join("json_cache/regulatory/21"))
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .collect();
        shards.sort();
        for path in shards {
            let features: Vec<serde_json::Value> =
                serde_json::from_reader(MultiGzDecoder::new(std::fs::File::open(&path).unwrap()))
                    .unwrap();
            for feature in features {
                if feature["feature_type"] != "MotifFeature" {
                    continue;
                }
                let (start, end) = (
                    feature["start"].as_u64().unwrap(),
                    feature["end"].as_u64().unwrap(),
                );
                let strand = feature["strand"].as_i64().unwrap() as i8;
                let cut = motif_sequence(&genome, "21", start, end, strand).unwrap();
                if cut.contains('N') {
                    masked += 1;
                    continue;
                }
                assert_eq!(
                    Some(cut.as_str()),
                    feature["seq"].as_str(),
                    "{} of {}",
                    feature["stable_id"],
                    path.display()
                );
                match strand {
                    1 => forward += 1,
                    -1 => reverse += 1,
                    other => panic!("strand {other}"),
                }
            }
        }
        assert_eq!((forward, reverse, masked), (836, 780, 1237));
    }

    /// The regulatory build's version is the `version` of the row of
    /// `regulatory_build` marked current, else of the table's first row; a
    /// NULL version, a table with no row and a directory without the table
    /// give none.
    #[test]
    fn the_regulatory_build_version_is_the_current_rows() {
        let dir = tempfile::tempdir().unwrap();
        let dumps = dir.path().join("funcgen");
        std::fs::create_dir(&dumps).unwrap();
        let inputs = RegulatoryInputs {
            funcgen_dumps: Some(dumps.clone()),
            ..activity_gffs(dir.path())
        };
        let version = || {
            let out = tempfile::tempdir().unwrap();
            build_regulatory(&inputs, out.path(), &SequenceNames::default(), None)
                .unwrap()
                .regbuild
        };
        assert_eq!(version(), None);
        let row = |id: u32, version: &str, current: u32| {
            format!(
                "{id}\ta build\t116\ta build\t{version}\t2024-7\t2024-7\t19\t16\t{current}\t\\N\n"
            )
        };
        for (rows, expected) in [
            (row(1, "1.0", 0) + &row(2, "2.0", 1), Some("2.0")),
            (row(1, "1.0", 0) + &row(2, "2.0", 0), Some("1.0")),
            (row(1, "\\N", 1) + &row(2, "2.0", 0), None),
            (String::new(), None),
        ] {
            write_text(&dumps, "regulatory_build.txt", &rows);
            assert_eq!(version().as_deref(), expected, "{rows:?}");
        }
    }

    /// A build given neither an activity file, the funcgen tables nor a
    /// genome FASTA writes the tree of the GFF files alone, byte for byte: a
    /// regulatory feature with an empty activity string, a motif with no
    /// sequence, a matrix as its stub. Its facts say so: no cell type, and
    /// neither activity, counts nor sequences on what the tree holds.
    #[test]
    fn a_build_given_the_gff_files_alone_writes_their_tree_alone() {
        let dir = tempfile::tempdir().unwrap();
        let out = tempfile::tempdir().unwrap();
        let stats = build_regulatory(
            &RegulatoryInputs {
                regulatory_gff: Some(write_gff(
                    dir.path(),
                    "regulatory.gff3.gz",
                    &["21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=R1"],
                )),
                motif_gff: Some(write_gff(
                    dir.path(),
                    "motifs.gff3.gz",
                    &["21\tEnsembl\tTF_binding_site\t3000\t3002\t.\t-\t.\tID=M1;binding_matrix_id=PFM_A;transcription_factor=GFF_A"],
                )),
                ..Default::default()
            },
            out.path(),
            &SequenceNames::default(),
            None,
        )
        .unwrap();
        let written =
            |file: &str| std::fs::read_to_string(out.path().join("regulatory").join(file)).unwrap();
        assert_eq!(
            written("21/1-1000000.json"),
            r#"[
  {
    "feature_type": "RegulatoryFeature",
    "stable_id": "R1",
    "start": 1000,
    "end": 1200,
    "biotype": "enhancer",
    "activity": "",
    "source": 0,
    "ordinal": 0
  },
  {
    "feature_type": "MotifFeature",
    "stable_id": "M1",
    "start": 3000,
    "end": 3002,
    "strand": -1,
    "matrix": "PFM_A",
    "source": 2,
    "ordinal": 0
  }
]"#
        );
        assert_eq!(
            written("binding_matrices.json"),
            r#"{
  "PFM_A": {
    "length": 3,
    "transcription_factors": [
      "GFF_A"
    ]
  }
}"#
        );
        assert_eq!(stats.cell_types, None);
        assert_eq!(stats.regbuild, None);
        assert!(!stats.activity_on_every_feature);
        assert!(!stats.counts_on_every_matrix);
        assert!(!stats.sequence_on_every_motif);
    }
}
