// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Regulatory and motif features read at run time from the GFF files the
//! Ensembl site distributes (`--regulatory_gff`), as
//! `AnnotationSource/File/RegFeat.pm` reads them at release/116.3: each file
//! bgzip-compressed with its tabix index, one source per key of the flag
//! (`file`, `emars`, `motifs`), a regulatory file yielding regulatory features
//! of the eight accepted types and a motif file motif features whose binding
//! matrix is a stub known by its id, span and factors alone, so the motif rows
//! carry no HIGH_INF_POS or MOTIF_SCORE_CHANGE (`OutputFactory.pm:1944`).
//!
//! Ensembl VEP reads a file region by region, padding each query by 10 kb
//! under `--extended_promoters` so a promoter whose extended bounds reach into
//! a region is found (`RegFeat.pm:291-329`), and collapses the records two
//! regions both return by the digest of the raw record (`:600-640`). This
//! source reads a chromosome whole into one index, so no padding is needed and
//! only a record the file holds twice collapses.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use noodles::bgzf;
use noodles::core::region::Interval;
use noodles::csi::BinningIndex;
use rustc_hash::{FxHashMap, FxHashSet};
use vep_builtins::tabix::{find_index_path, MAX_LINE_BYTES};
use vep_builtins::PluginError;
use vep_core::binding_matrix::BindingMatrix;
use vep_core::regulatory_gff::{
    decode_attribute, decode_type, duplicate_identifier_warning, feature_id, ignored_bound_warning,
    transcription_factors, unindexed_coordinate, unsupported_type_warning, EXTENDED_BOUNDS_TYPES,
    MOTIF_FEATURE_TYPES, REGULATORY_FEATURE_TYPES,
};
use vep_core::synonyms::ChromosomeSynonyms;

use super::source::{Feature, FeatureKind};
use crate::config::{CliError, RegulatoryGff};
use crate::input_shape::{load_index, TabixIndex};

/// The warning `--cell_type` draws from a regulatory GFF source
/// (`RegFeat.pm:179-184`), written once however many regulatory files the flag
/// names (`AnnotationSourceAdaptor.pm:281-300`); a motif file alone draws none.
pub const CELL_TYPE_WARNING: &str = "--cell_type is ignored for regulatory features from --regulatory_gff; a GFF carries no epigenome activity data";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Regulatory,
    Motif,
}

impl Kind {
    fn accepted_types(self) -> &'static [&'static str] {
        match self {
            Kind::Regulatory => &REGULATORY_FEATURE_TYPES,
            Kind::Motif => &MOTIF_FEATURE_TYPES,
        }
    }
}

/// One file of the flag, open with its index.
struct GffFile {
    path: PathBuf,
    index: TabixIndex,
    kind: Kind,
    /// The file's position in the flag's source order, which its features
    /// carry ([`Feature::file`]).
    position: u8,
    /// The feature types this file was warned about, each once
    /// (`RegFeat.pm:319-322`).
    warned_types: Mutex<FxHashSet<String>>,
    /// The identifiers read from this file so far, across its chromosomes, as
    /// `_seen_ids` spans one source's run (`RegFeat.pm:485`); a chromosome
    /// read merges its own set here once.
    seen_ids: Mutex<FileIds>,
}

/// The identifiers one file has yielded and the chromosomes, under the file's
/// own spellings, whose reads were merged: a chromosome two files spell
/// differently is entered under each spelling and so read once per entry, and
/// the second read of its records is no collision (Ensembl VEP keys the
/// check on the record's digest, `RegFeat.pm:487`).
#[derive(Default)]
struct FileIds {
    ids: FxHashSet<Arc<str>>,
    merged: FxHashSet<String>,
}

/// The files `--regulatory_gff` names, in Ensembl VEP's source order.
pub(super) struct GffRegulatorySource {
    files: Vec<GffFile>,
    extended_promoters: bool,
    /// The cache's chromosome synonyms, through which a chromosome's name is
    /// matched to a file's own spelling of it, as `get_source_chr_name` matches
    /// it (`RegFeat.pm:296`, `BaseVEP.pm:689-738`).
    synonyms: Arc<ChromosomeSynonyms>,
}

impl GffRegulatorySource {
    /// Opens every file the flag names, under the `cli` component: a file that
    /// does not exist (`AnnotationSourceAdaptor.pm:266`), a remote URL, which
    /// vep-rs never reads, a file without its `.tbi` or `.csi` index and a
    /// file whose index cannot be read or carries no tabix header are refused
    /// by key and name, written `<key>=<path>`, a bare file name under the key
    /// `file`. `cell_type_given` draws [`CELL_TYPE_WARNING`] when a
    /// regulatory file is among them.
    pub(super) fn open(
        spec: &RegulatoryGff,
        extended_promoters: bool,
        cell_type_given: bool,
        synonyms: Arc<ChromosomeSynonyms>,
    ) -> Result<Self, CliError> {
        let mut files = Vec::new();
        for (key, file) in spec.sources() {
            if is_remote(file) {
                return Err(CliError::cli(format!(
                    "--regulatory_gff {key}={file} is a remote URL; vep-rs reads a local file alone"
                )));
            }
            let path = Path::new(file);
            if !path.is_file() {
                return Err(CliError::cli(format!(
                    "--regulatory_gff {key}={file} not found"
                )));
            }
            let index_path = find_index_path(path).map_err(|_| {
                CliError::cli(format!("--regulatory_gff {key}={file} has no tabix index"))
            })?;
            let index = load_index(&index_path).map_err(|error| {
                CliError::cli(format!(
                    "--regulatory_gff {key}={file}: {}",
                    engine_error_text(error)
                ))
            })?;
            if index.header().is_none() {
                return Err(CliError::cli(format!(
                    "--regulatory_gff {key}={file}: index {} carries no tabix header, so its records cannot be mapped to chromosomes; index the file with tabix",
                    index_path.display()
                )));
            }
            files.push(GffFile {
                path: path.to_path_buf(),
                index,
                kind: if key == "motifs" {
                    Kind::Motif
                } else {
                    Kind::Regulatory
                },
                position: u8::try_from(files.len()).expect("at most three files"),
                warned_types: Mutex::default(),
                seen_ids: Mutex::default(),
            });
        }
        if cell_type_given && files.iter().any(|file| file.kind == Kind::Regulatory) {
            crate::warnings::global().warn(CELL_TYPE_WARNING);
        }
        Ok(Self {
            files,
            extended_promoters,
            synonyms,
        })
    }

    /// Every chromosome name the files' indexes list, each once, in byte order.
    pub(super) fn chromosomes(&self) -> BTreeSet<String> {
        self.files
            .iter()
            .flat_map(|file| file.reference_names())
            .collect()
    }

    /// The features of one chromosome from every file, in file order: `name`
    /// is a spelling one of the files lists ([`Self::chromosomes`]), matched
    /// to each other file's own spelling through the synonyms, the spelling a
    /// coordinate-derived identifier carries whatever the cache or the input
    /// calls the chromosome (`AnnotationSource.pm:179,198` maps the input's
    /// name onto the file's before the region reaches `RegFeat.pm:372,583`); a
    /// file that lacks the chromosome under any spelling contributes none.
    pub(super) fn read_chromosome(&self, name: &str) -> anyhow::Result<Vec<Feature>> {
        let mut read = ChromosomeRead {
            extended_promoters: self.extended_promoters,
            features: Vec::new(),
            matrices: FxHashMap::default(),
            biotypes: FxHashSet::default(),
        };
        for file in &self.files {
            let names = file.reference_names();
            let file_chr = self
                .synonyms
                .source_chr_name(name, |candidate| names.contains(candidate));
            file.read_chromosome(&file_chr, Some(&mut read))?;
        }
        Ok(read.features)
    }

    /// Puts every record of every file through the refusals of a chromosome's
    /// read ([`GffFile::read_chromosome`]), file by file in the flag's source
    /// order and, within a file, chromosome by chromosome in its index's
    /// order, and returns the first refusal. No feature is kept and no warning
    /// drawn: a chromosome's warnings come with its read.
    pub(super) fn check_records(&self) -> anyhow::Result<()> {
        for file in &self.files {
            for name in file.header().reference_sequence_names() {
                file.read_chromosome(&String::from_utf8_lossy(name.as_ref()), None)?;
            }
        }
        Ok(())
    }
}

/// The two schemes of `^(ht|f)tp://` (`AnnotationSourceAdaptor.pm:267`),
/// `http://` and `ftp://`, and `https://` beside them.
fn is_remote(file: &str) -> bool {
    ["http://", "https://", "ftp://"]
        .iter()
        .any(|scheme| file.starts_with(scheme))
}

/// The motif stubs built so far, one per matrix id, span and factor list as
/// written, shared by every motif naming the same three.
type MatrixStubs = FxHashMap<(String, u64, Box<[u8]>), Arc<BindingMatrix>>;

/// The state of one chromosome's read across the files.
struct ChromosomeRead {
    extended_promoters: bool,
    features: Vec<Feature>,
    matrices: MatrixStubs,
    /// The regulatory feature types seen, each one allocation.
    biotypes: FxHashSet<Arc<str>>,
}

impl GffFile {
    fn header(&self) -> &noodles::csi::binning_index::index::Header {
        self.index.header().expect("checked when the file opened")
    }

    /// The reference names the index lists, as strings.
    fn reference_names(&self) -> FxHashSet<String> {
        self.header()
            .reference_sequence_names()
            .iter()
            .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
            .collect()
    }

    /// Reads every record of `file_chr`, the file's own spelling of the
    /// chromosome, through the index's chunks for the whole sequence; a file
    /// whose index lacks the name contributes none. Without a `read` to fill,
    /// each line meets the length limit and its record's refusals
    /// ([`Self::record`]) and nothing else: the file's identifiers stay as
    /// they are.
    fn read_chromosome(
        &self,
        file_chr: &str,
        mut read: Option<&mut ChromosomeRead>,
    ) -> anyhow::Result<()> {
        let names = self.header().reference_sequence_names();
        let Some(reference_sequence_id) = names.get_index_of(file_chr.as_bytes()) else {
            return Ok(());
        };
        let at_path = |error: std::io::Error| anyhow::anyhow!("{}: {error}", self.path.display());
        let chunks = self
            .index
            .query(reference_sequence_id, Interval::from(..))
            .map_err(at_path)?;
        let file = std::fs::File::open(&self.path).map_err(at_path)?;
        let mut reader = bgzf::io::Reader::new(file);
        let mut line = Vec::new();
        let mut seen = SeenRecords::default();
        for chunk in chunks {
            reader.seek(chunk.start()).map_err(at_path)?;
            while reader.virtual_position() < chunk.end() {
                line.clear();
                let read_bytes = (&mut reader)
                    .take(MAX_LINE_BYTES as u64 + 1)
                    .read_until(b'\n', &mut line)
                    .map_err(at_path)?;
                if read_bytes == 0 {
                    break;
                }
                if line.len() > MAX_LINE_BYTES && line.last() != Some(&b'\n') {
                    anyhow::bail!(
                        "{}: a line of chromosome {file_chr} is longer than the {MAX_LINE_BYTES}-byte limit for one line",
                        self.path.display()
                    );
                }
                let text = trim_line_end(&line);
                if text.is_empty() || text[0] == b'#' {
                    continue;
                }
                let kept = read.as_deref_mut().map(|read| (read, &mut seen));
                self.record(text, file_chr, kept)?;
            }
        }
        if read.is_none() {
            return Ok(());
        }
        // An identifier another chromosome of the file already used is a
        // collision too (`RegFeat.pm:485`, one set per source over the run).
        let mut file_ids = self
            .seen_ids
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if !file_ids.merged.insert(file_chr.to_string()) {
            return Ok(());
        }
        for id in seen.ids {
            if file_ids.ids.contains(&id) {
                self.warn_duplicate_identifier(&id);
            } else {
                file_ids.ids.insert(id);
            }
        }
        Ok(())
    }

    /// The warning two distinct records of this file sharing an identifier
    /// draw ([`duplicate_identifier_warning`]), once per identifier and file
    /// through the sink's seen set.
    fn warn_duplicate_identifier(&self, id: &str) {
        crate::warnings::global().warn(&duplicate_identifier_warning(
            id,
            &self.path.display().to_string(),
        ));
    }

    /// The warning a promoter's `extended_start` or `extended_end` draws under
    /// `--extended_promoters` when its value is not a whole number
    /// ([`ignored_bound_warning`]): the bound is ignored and the promoter
    /// keeps its own.
    fn warn_ignored_bound(
        &self,
        attribute: &str,
        value: &str,
        feature_type: &str,
        stable_id: &str,
    ) {
        crate::warnings::global().warn(&ignored_bound_warning(
            attribute,
            value,
            feature_type,
            stable_id,
            &self.path.display().to_string(),
        ));
    }

    /// The refusal of a line of `file_chr`'s read whose `name` column, the
    /// start or the end, holds `value`, no whole number; `record` names the
    /// line ([`record_name`]).
    fn not_a_whole_number(
        &self,
        name: &str,
        value: &[u8],
        record: &str,
        file_chr: &str,
    ) -> anyhow::Error {
        anyhow::anyhow!(
            "{}: the {name} of {record} of chromosome {file_chr} is not a whole number: {}",
            self.path.display(),
            String::from_utf8_lossy(value)
        )
    }

    /// One record as `_record_to_regfeat` builds it (`RegFeat.pm:349-385`): a
    /// line whose first column is not the chromosome's name, which a chunk's
    /// bounds alone do not keep out, is skipped where `tabix` indexes its
    /// start and its end and refused where it does not, as a line of any
    /// type is, so a line that opens with a tab refuses the read for the
    /// type that stands where a start is read; the type is the column with its
    /// `%XX` escapes decoded
    /// (`:312`, `get_type` of `GFF3.pm:161-164`), and a type the file does not
    /// accept is warned about once and skipped (`:314-323`), its start and end
    /// read no further than whether `tabix` indexes them
    /// ([`unindexed_coordinate`]): one it does not refuses the line as a kept
    /// line's does, where Ensembl VEP's query ends at a line whose end fails
    /// the parse and loses the records after it (`:311`; htslib 1.9
    /// `tbx.c:174-175`, `hts.c:2643`); `tabix` ends a column at a NUL byte
    /// as at a tab, so a line of any type that holds one is held to the
    /// start and the end `tabix` reads in it, as the cache builder's bake
    /// of the same file holds it, and a NUL byte in a sequence name, in a
    /// source or inside a type refuses the read for the type, or its tail,
    /// that stands where a start is read; a promoter is widened to its
    /// extended bounds
    /// under `--extended_promoters` (`:360-365`), a bound that is not a whole
    /// number ignored with a warning; a motif's strand is the
    /// column's, `0` when it holds none (`:369-370`), and a regulatory
    /// feature's is `0` as the cache's are, since no output reads it
    /// (`OutputFactory.pm:1876-1893` sets no STRAND); the identifier is `ID`,
    /// else `Name`, else the file's own spelling of the chromosome with the
    /// coordinates and type (`:575-584`, `$chr` the region name
    /// `AnnotationSource.pm:179,198` mapped onto the file's spelling); a record the
    /// file holds twice is kept once (`:628-640`), and a distinct record under
    /// an identifier the file already used is kept with the collision warned
    /// about (`:478-497`). With nothing `kept`, the read state and the records
    /// seen, the record meets its refusals alone: it draws no warning and
    /// builds no feature.
    fn record(
        &self,
        text: &[u8],
        file_chr: &str,
        kept: Option<(&mut ChromosomeRead, &mut SeenRecords)>,
    ) -> anyhow::Result<()> {
        let mut columns = text.split(|byte| *byte == b'\t');
        let mut column = |name: &str| {
            columns.next().ok_or_else(|| {
                anyhow::anyhow!(
                    "{}: a record of chromosome {file_chr} has no {name} column: {}",
                    self.path.display(),
                    String::from_utf8_lossy(text)
                )
            })
        };
        if column("sequence")? != file_chr.as_bytes() {
            // The index build files a line whose start it cannot parse among
            // the lines of the sequence before it (htslib 1.9 `tbx.c:235-236`).
            let Some((name, value)) = unindexed_coordinate(text) else {
                return Ok(());
            };
            let mut columns = text.split(|byte| *byte == b'\t').skip(2);
            let feature_type = decode_type(columns.next().unwrap_or_default());
            let attributes = Attributes::parse(columns.nth(5).unwrap_or_default());
            let record = record_name(&feature_type, attributes.get("ID").as_deref());
            return Err(self.not_a_whole_number(name, value, &record, file_chr));
        }
        column("source")?;
        let feature_type = decode_type(column("type")?).into_owned();
        let start_text = column("start")?;
        let end_text = column("end")?;
        column("score")?;
        let strand = match (self.kind, column("strand")?) {
            (Kind::Motif, b"+") => 1,
            (Kind::Motif, b"-") => -1,
            _ => 0,
        };
        column("phase")?;
        let attributes = Attributes::parse(column("attributes")?);
        // The record's name is built only on a refusal, never for a record the
        // reader keeps.
        let record = || record_name(&feature_type, attributes.get("ID").as_deref());
        let not_a_whole_number =
            |name: &str, value: &[u8]| self.not_a_whole_number(name, value, &record(), file_chr);
        if !self.kind.accepted_types().contains(&feature_type.as_str()) {
            if let Some((name, value)) = unindexed_coordinate(text) {
                return Err(not_a_whole_number(name, value));
            }
            if kept.is_none() {
                return Ok(());
            }
            let mut warned = self
                .warned_types
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            if warned.insert(feature_type.clone()) {
                crate::warnings::global().warn(&unsupported_type_warning(
                    &feature_type,
                    &self.path.display().to_string(),
                ));
            }
            return Ok(());
        }
        // A NUL byte ends a column to `tabix` alone, so only a line that
        // holds one can have a start and an end there that the columns
        // read below do not hold.
        if text.contains(&0) {
            if let Some((name, value)) = unindexed_coordinate(text) {
                return Err(not_a_whole_number(name, value));
            }
        }
        let id = attributes.get("ID");
        let coordinate = |name: &str, text: &[u8]| -> anyhow::Result<u64> {
            std::str::from_utf8(text)
                .ok()
                .and_then(|text| text.parse().ok())
                .ok_or_else(|| not_a_whole_number(name, text))
        };
        let mut start = coordinate("start", start_text)?;
        let mut end = coordinate("end", end_text)?;
        // A record spans a one-based ordered pair, the rule the cache builder's
        // bake of the same file holds a line to; a start of 0 or an end below
        // its start names no span.
        if start == 0 || end < start {
            anyhow::bail!(
                "{}: the coordinates {start}-{end} of {} of chromosome {file_chr} are not 1-based and ordered",
                self.path.display(),
                record()
            );
        }
        let Some((read, seen)) = kept else {
            return Ok(());
        };
        // The file is sorted by the start column, the key of the duplicate
        // check, which runs before any widening.
        if !seen.insert(start, text) {
            return Ok(());
        }
        let mut ignored_bounds = Vec::new();
        if read.extended_promoters && EXTENDED_BOUNDS_TYPES.contains(&feature_type.as_str()) {
            match attributes.number("extended_start") {
                Ok(Some(extended)) => start = start.min(extended),
                Ok(None) => {}
                Err(value) => ignored_bounds.push(("extended_start", value)),
            }
            match attributes.number("extended_end") {
                Ok(Some(extended)) => end = end.max(extended),
                Ok(None) => {}
                Err(value) => ignored_bounds.push(("extended_end", value)),
            }
        }
        // `Name` is read only when `ID` is absent or empty, the one case in
        // which `_record_get_id` (`RegFeat.pm:575-581`) reaches it.
        let name = match id.as_deref() {
            None | Some("") => attributes.get("Name"),
            Some(_) => None,
        };
        let stable_id: Arc<str> = Arc::from(feature_id(
            id.as_deref(),
            name.as_deref(),
            file_chr,
            start,
            end,
            &feature_type,
        ));
        for (attribute, value) in ignored_bounds {
            self.warn_ignored_bound(attribute, &value, &feature_type, &stable_id);
        }
        if !seen.ids.insert(Arc::clone(&stable_id)) {
            self.warn_duplicate_identifier(&stable_id);
        }
        let kind = match self.kind {
            Kind::Regulatory => FeatureKind::Regulatory {
                biotype: Some(shared_biotype(&mut read.biotypes, &feature_type)),
                activity: Box::from(""),
            },
            Kind::Motif => {
                let matrix_id = attributes.get("binding_matrix_id");
                let factors = attributes.written("transcription_factor");
                FeatureKind::Motif {
                    score: None,
                    seq: None,
                    matrix: stub_matrix(
                        &mut read.matrices,
                        matrix_id.as_deref().unwrap_or(""),
                        end.saturating_sub(start) + 1,
                        factors.unwrap_or(b""),
                    ),
                }
            }
        };
        read.features.push(Feature {
            stable_id,
            start,
            end,
            strand,
            dbid: None,
            file: self.position,
            ordinal: u32::try_from(read.features.len()).unwrap_or(u32::MAX),
            kind,
        });
        Ok(())
    }
}

/// The name a refusal gives a line: its type and its `ID`, its type alone
/// when it has none or an empty one.
fn record_name(feature_type: &str, id: Option<&str>) -> String {
    match id {
        None | Some("") => format!("a record of type {feature_type}"),
        Some(id) => format!("{feature_type} record {id}"),
    }
}

/// The one allocation of `feature_type` every regulatory feature of that type
/// on the chromosome shares as its biotype.
fn shared_biotype(biotypes: &mut FxHashSet<Arc<str>>, feature_type: &str) -> Arc<str> {
    if let Some(shared) = biotypes.get(feature_type) {
        return Arc::clone(shared);
    }
    let biotype: Arc<str> = Arc::from(feature_type);
    biotypes.insert(Arc::clone(&biotype));
    biotype
}

/// The binding matrix stub of a motif record (`RegFeat.pm:431-444`): the
/// `binding_matrix_id` as its stable id, the motif's span as its length, the
/// `transcription_factor` list as written split on its commas, each name
/// then decoded, as its complexes ([`transcription_factors`]), and no
/// counts, so the motif's rows carry MOTIF_NAME, MOTIF_POS and
/// TRANSCRIPTION_FACTORS and nothing scored. Motifs naming one id, span and
/// list share one stub.
fn stub_matrix(
    matrices: &mut MatrixStubs,
    id: &str,
    length: u64,
    factors: &[u8],
) -> Arc<BindingMatrix> {
    let key = (id.to_string(), length, Box::from(factors));
    if let Some(matrix) = matrices.get(&key) {
        return Arc::clone(matrix);
    }
    let factors: Arc<[Arc<str>]> = transcription_factors(Some(factors))
        .into_iter()
        .map(Arc::from)
        .collect();
    let matrix = Arc::new(
        BindingMatrix::new(id, length, Vec::new(), factors).expect("no counts to disagree"),
    );
    matrices.insert(key, Arc::clone(&matrix));
    matrix
}

/// The ninth column as `GFF3::get_attributes` reads it (ensembl-io
/// `GFF3.pm:190-198`): the `;`-separated pairs, each split on `=` into key
/// and value, the key the text before the first `=` and the value the text
/// between the first and the second (`:194`, a two-variable `split`), a
/// later pair overwriting an earlier one with its key, a pair without `=`
/// with the undefined value (`:195`), so the key reads as absent; a value's
/// escapes are decoded over its bytes, then read as text
/// ([`decode_attribute`]). A byte that is no UTF-8 in a value once decoded
/// reads as the replacement character, as every input's does, where Ensembl
/// VEP carries the byte into its rows.
struct Attributes<'a> {
    pairs: Vec<(&'a [u8], Option<&'a [u8]>)>,
}

impl<'a> Attributes<'a> {
    fn parse(column: &'a [u8]) -> Self {
        let pairs = column
            .split(|byte| *byte == b';')
            .filter(|pair| !pair.is_empty())
            .map(|pair| {
                let mut fields = pair.splitn(3, |byte| *byte == b'=');
                let key = fields.next().unwrap_or_default();
                let value = fields.next();
                (key, value)
            })
            .collect();
        Self { pairs }
    }

    /// The value the last pair written under `key` holds, as written; `None`
    /// for a key no pair names or whose last pair has no `=`.
    fn written(&self, key: &str) -> Option<&'a [u8]> {
        self.pairs
            .iter()
            .rev()
            .find(|(pair_key, _)| *pair_key == key.as_bytes())
            .and_then(|(_, value)| *value)
    }

    /// [`Self::written`] decoded ([`decode_attribute`]).
    fn get(&self, key: &str) -> Option<Cow<'a, str>> {
        self.written(key).map(decode_attribute)
    }

    /// The value under `key` as a coordinate ([`Bound`]): `None` for an
    /// absent key, and the value itself as the error when it is not a whole
    /// number.
    fn number(&self, key: &str) -> Bound<'a> {
        match self.get(key) {
            None => Ok(None),
            Some(value) => value.parse().map(Some).map_err(|_| value),
        }
    }
}

/// A coordinate attribute as read: absent, a whole number, or the text that
/// is not one.
type Bound<'a> = Result<Option<u64>, Cow<'a, str>>;

/// What one chromosome read of a file has seen: the raw records of the start
/// position in hand (a tabix-indexed file is sorted by start, so a record the
/// file holds twice meets its earlier copy here) and the identifier of every
/// feature built, so a distinct record under a used identifier meets it too.
#[derive(Default)]
struct SeenRecords {
    start: u64,
    records: FxHashSet<Box<[u8]>>,
    ids: FxHashSet<Arc<str>>,
}

impl SeenRecords {
    /// Whether `record` at `start` is new.
    fn insert(&mut self, start: u64, record: &[u8]) -> bool {
        if start != self.start {
            self.start = start;
            self.records.clear();
        }
        self.records.insert(Box::from(record))
    }
}

/// The engine's error reduced to its text: the engine names a plugin, and the
/// file is a regulatory file.
fn engine_error_text(error: PluginError) -> String {
    match error {
        PluginError::Init(message) | PluginError::Run(message) => message,
        PluginError::Io(error) => error.to_string(),
        other => other.detail(),
    }
}

/// `line` without its line end, `\n` or `\r\n`.
fn trim_line_end(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use std::fs;
    use std::io::Write;

    use noodles::core::Position;
    use noodles::csi::binning_index::index::header;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles::tabix;
    use vep_core::cache_info::CacheInfo;
    use vep_core::consequence::{Consequence, FeatureType};
    use vep_core::variant::{InputVariant, VariantClass};

    use super::super::fixtures::{resolved_cache, variant, write_fixture_cache};
    use super::super::source::{RegulatorySource, GFF_REPLACES_CACHE_WARNING};
    use super::*;
    use crate::cache_resolve::ResolvedCache;
    use crate::warnings::{global, sink_test_guard, Target};

    const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/regulatory-gff");

    fn fixture(name: &str) -> String {
        format!("{FIXTURES}/{name}")
    }

    /// The three fixture files under their keys.
    fn fixture_spec() -> RegulatoryGff {
        RegulatoryGff {
            file: Some(fixture("regulatory_features.gff3.gz")),
            emars: Some(fixture("EMARs.gff.gz")),
            motifs: Some(fixture("motif_features.gff3.gz")),
        }
    }

    /// The fixtures' cache at `dir` with metadata declaring no regulatory data.
    fn cache_without_regulatory(dir: &Path) -> ResolvedCache {
        ResolvedCache {
            info: CacheInfo::default(),
            ..resolved_cache(dir)
        }
    }

    fn open(spec: &RegulatoryGff, extended_promoters: bool) -> Arc<RegulatorySource> {
        let dir = tempfile::tempdir().unwrap();
        let source = RegulatorySource::open_gff(
            spec,
            extended_promoters,
            false,
            &cache_without_regulatory(dir.path()),
        )
        .unwrap();
        source.prewarm("21").unwrap();
        source
    }

    /// Writes `lines`, GFF records sorted by chromosome and start, as a
    /// bgzip-compressed file under `dir` with its `.tbi` beside it.
    fn write_indexed_gff(dir: &Path, name: &str, lines: &[&str]) -> String {
        let lines: Vec<&[u8]> = lines.iter().map(|line| line.as_bytes()).collect();
        write_indexed_gff_bytes(dir, name, &lines)
    }

    /// [`write_indexed_gff`] over lines given as bytes, so a line can carry a
    /// value that is not UTF-8; the first, fourth and fifth columns are UTF-8.
    /// A line is filed by the columns `tabix` reads, which end at a tab or a
    /// NUL byte, and one whose first column is empty or blank, or whose
    /// start or end `tabix` fails, under the sequence of the line before
    /// it, where `tabix` files a line it cannot parse (htslib 1.9
    /// `tbx.c:235-236`); a file's first line is filed under its own name.
    fn write_indexed_gff_bytes(dir: &Path, name: &str, lines: &[&[u8]]) -> String {
        let path = dir.join(name);
        let mut writer = bgzf::io::Writer::new(fs::File::create(&path).unwrap());
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(header::Builder::gff().build());
        let mut sequence = "";
        for line in lines {
            let start = writer.virtual_position();
            writer.write_all(line).unwrap();
            writer.write_all(b"\n").unwrap();
            let end = writer.virtual_position();
            if line.starts_with(b"#") {
                continue;
            }
            let columns: Vec<&str> = vep_core::regulatory_gff::index_columns(line)
                .map(|column| std::str::from_utf8(column).unwrap_or(""))
                .collect();
            // A column that is no one-based position is indexed as the first
            // base, so a query's chunk still holds the record the reader must
            // refuse; `tabix` files a start below one there too (`tbx.c:99`).
            let position = |column: &str| match column.parse::<usize>() {
                Ok(base) if base >= 1 => Position::try_from(base).unwrap(),
                _ => Position::MIN,
            };
            let (first, second) = (position(columns[3]), position(columns[4]));
            let parsed = sequence.is_empty() || unindexed_coordinate(line).is_none();
            if !columns[0].trim().is_empty() && parsed {
                sequence = columns[0];
            }
            indexer
                .add_record(
                    sequence,
                    first.min(second),
                    first.max(second),
                    Chunk::new(start, end),
                )
                .unwrap();
        }
        writer.finish().unwrap();
        tabix::fs::write(format!("{}.tbi", path.display()), &indexer.build()).unwrap();
        path.to_string_lossy().into_owned()
    }

    fn file_spec(path: String) -> RegulatoryGff {
        RegulatoryGff {
            file: Some(path),
            ..Default::default()
        }
    }

    /// The rows a one-base variant at `position` gains, as text.
    fn rows_at(source: &RegulatorySource, position: u64) -> String {
        let mut v = variant(position, position, b"A", b"C");
        source.annotate(&mut v);
        format!("{:?}", v.transcript_consequences)
    }

    /// One row: id, feature type, biotype, strand and terms.
    type RowSummary = (String, FeatureType, Option<String>, i8, Vec<Consequence>);

    /// The regulatory and motif rows a one-base variant at `position` gains.
    fn row_summary(source: &RegulatorySource, position: u64) -> Vec<RowSummary> {
        let mut v = variant(position, position, b"A", b"C");
        source.annotate(&mut v);
        let mut rows: Vec<_> = v
            .transcript_consequences
            .iter()
            .map(|row| {
                (
                    row.transcript_id().to_string(),
                    row.feature_type,
                    row.biotype().map(str::to_string),
                    row.strand,
                    row.consequences.to_vec(),
                )
            })
            .collect();
        rows.sort();
        rows
    }

    /// The three fixture files read as one source give, at every feature's
    /// first and last base, the rows a converted cache gives from the same
    /// records: the cache route writes each regulatory feature with an empty
    /// activity string under an empty cell-type list and each motif naming a
    /// stub matrix of its span and factors, the shape `vep-cache-builder`
    /// writes from the same files.
    #[test]
    fn gff_sources_build_the_same_rows_as_the_cache() {
        let gff = open(&fixture_spec(), false);
        let features = gff.features("21").unwrap();
        assert_eq!(features.len(), 516 + 279 + 2200);
        assert!(gff.matrices().is_none() && gff.cell_types().is_empty());

        // The cache route: one shard of every feature plus the stub table.
        let mut shard = Vec::new();
        let mut table: BTreeMap<String, serde_json::Value> = BTreeMap::new();
        for feature in features {
            let element = match &feature.kind {
                FeatureKind::Regulatory { biotype, .. } => serde_json::json!({
                    "feature_type": "RegulatoryFeature",
                    "stable_id": feature.stable_id.as_ref(),
                    "start": feature.start,
                    "end": feature.end,
                    "biotype": biotype.as_deref(),
                    "activity": "",
                }),
                FeatureKind::Motif { matrix, .. } => {
                    let factors: Vec<&str> = matrix
                        .transcription_factors
                        .iter()
                        .map(AsRef::as_ref)
                        .collect();
                    let entry = serde_json::json!({"length": matrix.length, "transcription_factors": factors});
                    let previous = table.insert(matrix.stable_id.to_string(), entry.clone());
                    assert!(
                        previous.is_none_or(|previous| previous == entry),
                        "{}",
                        matrix.stable_id
                    );
                    serde_json::json!({
                        "feature_type": "MotifFeature",
                        "stable_id": feature.stable_id.as_ref(),
                        "start": feature.start,
                        "end": feature.end,
                        "strand": feature.strand,
                        "matrix": matrix.stable_id.as_ref(),
                    })
                }
            };
            shard.push(element);
        }
        assert_eq!(table.len(), 50);
        let dir = tempfile::tempdir().unwrap();
        let chr = dir.path().join("regulatory/21");
        fs::create_dir_all(&chr).unwrap();
        fs::write(
            dir.path().join("regulatory/binding_matrices.json"),
            serde_json::to_vec(&table).unwrap(),
        )
        .unwrap();
        fs::write(
            chr.join("5000001-11000000.json"),
            serde_json::to_vec(&shard).unwrap(),
        )
        .unwrap();
        let mut cache = resolved_cache(dir.path());
        cache.info.cell_types = Vec::new();
        let from_cache = RegulatorySource::open(&cache, true).unwrap().unwrap();
        from_cache.prewarm("21").unwrap();
        assert_eq!(from_cache.features("21").unwrap().len(), features.len());

        let mut positions: BTreeSet<u64> = BTreeSet::new();
        for feature in features {
            positions.insert(feature.start);
            positions.insert(feature.end);
            positions.insert(feature.start.saturating_sub(1));
        }
        let mut compared = 0;
        for position in positions {
            let rows = rows_at(&gff, position);
            assert_eq!(rows, rows_at(&from_cache, position), "at 21:{position}");
            compared += usize::from(rows != "[]");
        }
        assert!(compared > 2_000, "{compared} positions with rows");
    }

    /// The identifier is the record's `ID`, else its `Name`, else the file's
    /// chromosome spelling with the feature's coordinates and type
    /// (`RegFeat.pm:575-584`), the coordinates as `--extended_promoters`
    /// widened them (`:357-372`); an empty `ID` falls through to `Name`.
    #[test]
    fn regulatory_feature_id_falls_back_to_coordinates() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "##gff-version 3",
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ENSR_WITH_ID;Name=named",
                "21\tEnsembl\tenhancer\t2000\t2200\t.\t.\t.\tName=ENSR_NAMED",
                "21\tEnsembl\tCTCF_binding_site\t3000\t3020\t.\t+\t.\tcolor=#40e0d0",
                "21\tEnsembl\tenhancer\t4000\t4200\t.\t.\t.\tID=;Name=ENSR_EMPTY_ID",
                "21\tEnsembl\tpromoter\t5000\t5500\t.\t.\t.\textended_start=4900;extended_end=5600",
            ],
        );
        let ids = |extended: bool| -> Vec<String> {
            open(&file_spec(path.clone()), extended)
                .features("21")
                .unwrap()
                .iter()
                .map(|feature| feature.stable_id.to_string())
                .collect()
        };
        assert_eq!(
            ids(false),
            [
                "ENSR_WITH_ID",
                "ENSR_NAMED",
                "21:3000-3020:CTCF_binding_site",
                "ENSR_EMPTY_ID",
                "21:5000-5500:promoter"
            ]
        );
        assert_eq!(ids(true)[4], "21:4900-5600:promoter");
        let source = open(&file_spec(path), false);
        let ctcf = &source.features("21").unwrap()[2];
        assert_eq!(
            ctcf.strand, 0,
            "a regulatory feature's strand reaches no output"
        );
        assert_eq!(
            row_summary(&source, 3010),
            [(
                "21:3000-3020:CTCF_binding_site".to_string(),
                FeatureType::RegulatoryFeature,
                Some("CTCF_binding_site".to_string()),
                0,
                vec![Consequence::RegulatoryRegionVariant]
            )]
        );
    }

    /// The ninth column split as `GFF3.pm:194` splits it: a value is the text
    /// between the first and the second `=`, so `ID=a=b` names the feature
    /// `a` and `binding_matrix_id=M=1` the matrix `M`; a later pair with the
    /// key wins, and a pair without `=` leaves its key absent, an earlier
    /// value included (`:195` writes the undefined value over it); the
    /// escapes of the cut value are decoded. The record with `ID=a=b`
    /// reaches the rows under the identifier `a`.
    #[test]
    fn attribute_values_are_cut_at_a_second_equals_sign() {
        let attributes = Attributes::parse(
            b"ID=a=b;Name=x=y=z;ID=c=d;flag;binding_matrix_id=M=1;k=%3D=v;extended_start=5;extended_start",
        );
        let get = |key| attributes.get(key);
        assert_eq!(get("ID").as_deref(), Some("c"));
        assert_eq!(get("Name").as_deref(), Some("x"));
        assert_eq!(get("flag"), None);
        assert_eq!(get("extended_start"), None);
        assert_eq!(attributes.number("extended_start"), Ok(None));
        assert_eq!(get("binding_matrix_id").as_deref(), Some("M"));
        assert_eq!(get("k").as_deref(), Some("="));
        assert_eq!(get("missing"), None);
        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "##gff-version 3",
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=a=b;Name=named",
            ],
        );
        let source = open(&file_spec(path), false);
        assert_eq!(source.features("21").unwrap()[0].stable_id.as_ref(), "a");
    }

    /// A record's coordinates are a one-based ordered pair, the rule the cache
    /// builder's bake of the same file holds a line to: a start of 0 or an end
    /// below its start refuses the chromosome's read, naming the file, the
    /// coordinates, the record by its type and `ID` (its type alone without
    /// one) and the chromosome, as a start or end that is not a whole number
    /// does, where Ensembl VEP builds the feature as written (`RegFeat.pm:291-325`
    /// reads what the index returns). A pair that is in order is read as it
    /// stands.
    #[test]
    fn a_record_whose_coordinates_are_not_a_one_based_ordered_pair_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let refused = |lines: &[&str], text: &str| {
            let path = write_indexed_gff(dir.path(), "features.gff3.gz", lines);
            let source = RegulatorySource::open_gff(
                &file_spec(path.clone()),
                false,
                false,
                &cache_without_regulatory(dir.path()),
            )
            .unwrap();
            let error = source.prewarm("21").unwrap_err();
            assert_eq!(format!("{error:#}"), format!("{path}: {text}"));
        };
        refused(
            &[
                "##gff-version 3",
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ENSR_A",
                "21\tEnsembl\tenhancer\t1400\t1300\t.\t.\t.\tID=ENSR_SWAPPED",
            ],
            "the coordinates 1400-1300 of enhancer record ENSR_SWAPPED of chromosome 21 are not 1-based and ordered",
        );
        refused(
            &[
                "##gff-version 3",
                "21\tEnsembl\tpromoter\t0\t1200\t.\t.\t.\tID=ENSR_ZERO",
            ],
            "the coordinates 0-1200 of promoter record ENSR_ZERO of chromosome 21 are not 1-based and ordered",
        );
        refused(
            &[
                "##gff-version 3",
                "21\tEnsembl\tenhancer\t1000\txyz\t.\t.\t.\tID=ENSR_END",
            ],
            "the end of enhancer record ENSR_END of chromosome 21 is not a whole number: xyz",
        );
        refused(
            &[
                "##gff-version 3",
                "21\tEnsembl\topen_chromatin_region\t1000\t900\t.\t.\t.\tName=ENSR_NAMED",
            ],
            "the coordinates 1000-900 of a record of type open_chromatin_region of chromosome 21 are not 1-based and ordered",
        );
        refused(
            &[
                "##gff-version 3",
                "21\tEnsembl\tTF_binding_site\t0\t900\t.\t.\t.\tID=;Name=ENSR_EMPTY_ID",
            ],
            "the coordinates 0-900 of a record of type TF_binding_site of chromosome 21 are not 1-based and ordered",
        );
        let path = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "##gff-version 3",
                "21\tEnsembl\tenhancer\t1000\t1000\t.\t.\t.\tID=ENSR_ONE_BASE",
            ],
        );
        let source = open(&file_spec(path), false);
        assert_eq!(source.features("21").unwrap().len(), 1);
    }

    /// The start and end of each enhancer [`enhancers_with_a_skipped_line`]
    /// writes on a sequence.
    const ENHANCER_SPANS: [(u64, u64); 3] =
        [(5011320, 5011550), (5013566, 5014200), (5015556, 5015987)];

    /// A file of two sequences, three enhancers on each, with `skipped`, a
    /// line of a type no regulatory file yields, after the first enhancer of
    /// `chr`.
    fn enhancers_with_a_skipped_line(chr: &str, skipped: &str) -> Vec<String> {
        let mut lines = vec!["##gff-version 3".to_string()];
        for sequence in ["20", "21"] {
            for (place, (start, end)) in ENHANCER_SPANS.into_iter().enumerate() {
                lines.push(format!(
                    "{sequence}\tEnsembl\tenhancer\t{start}\t{end}\t.\t.\t.\tID=ENSR{sequence}_{start}"
                ));
                if sequence == chr && place == 0 {
                    lines.push(skipped.to_string());
                }
            }
        }
        lines
    }

    /// The source over [`enhancers_with_a_skipped_line`], written under `dir`,
    /// and the file's path.
    fn open_with_a_skipped_line(
        dir: &Path,
        chr: &str,
        skipped: &str,
    ) -> (Arc<RegulatorySource>, String) {
        let lines = enhancers_with_a_skipped_line(chr, skipped);
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        let path = write_indexed_gff(dir, "features.gff3.gz", &lines);
        let source = RegulatorySource::open_gff(
            &file_spec(path.clone()),
            false,
            false,
            &cache_without_regulatory(dir),
        )
        .unwrap();
        (source, path)
    }

    /// The identifier and span of every feature read for `chr`.
    fn spans(source: &RegulatorySource, chr: &str) -> Vec<(String, u64, u64)> {
        source.prewarm(chr).unwrap();
        source
            .features(chr)
            .unwrap()
            .iter()
            .map(|feature| (feature.stable_id.to_string(), feature.start, feature.end))
            .collect()
    }

    /// A line of a type the file does not yield refuses the read of its
    /// sequence, a file's first or a later one, when `tabix` does not index
    /// its start or its end: a column no number leads, an empty one read on
    /// into the next, an end that reads below zero. The refusal is a kept
    /// line's for a coordinate that is not a whole number, naming the file,
    /// the column, the record by its type and `ID` (its type alone without
    /// one), the chromosome and the text. Past such an end (`xyz`, `-9`)
    /// Ensembl VEP 116.3 writes the rows of the records before the line and
    /// none for those after it in the region it reads, with exit 0: its query
    /// ends at the line (`RegFeat.pm:311`; htslib 1.9 `tbx.c:106` and `:142`
    /// fail the parse, `:174-175` hand the iterator no sequence id). Over
    /// such a start inside a sequence the release's `tabix` builds no index
    /// (`tbx.c:96`, `hts.c:1540`). The other sequence is read whole, and the
    /// refused line's type draws no warning.
    #[test]
    fn a_skipped_line_tabix_does_not_index_is_refused_on_either_sequence() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        for (start, end, column, value) in [
            ("5012000", "xyz", "end", "xyz"),
            ("5012000", "-9", "end", "-9"),
            ("abc", "5012100", "start", "abc"),
            ("", "xyz", "start", ""),
            ("5012000", "", "end", ""),
            ("5012000", "-", "end", "-"),
        ] {
            for (chr, other) in [("20", "21"), ("21", "20")] {
                let skipped =
                    format!("{chr}\tEnsembl\texon\t{start}\t{end}\t.\t.\t.\tID=SKIP{chr}");
                let (source, path) = open_with_a_skipped_line(dir.path(), chr, &skipped);
                let error = source.prewarm(chr).unwrap_err();
                assert_eq!(
                    format!("{error:#}"),
                    format!(
                        "{path}: the {column} of exon record SKIP{chr} of chromosome {chr} is \
                         not a whole number: {value}"
                    ),
                    "{start:?} {end:?} on {chr}"
                );
                assert_eq!(spans(&source, other).len(), 3, "{start:?} {end:?} on {chr}");
            }
        }
        let unnamed = "21\tEnsembl\texon\t5012000\txyz\t.\t.\t.\tName=SKIP";
        let (source, path) = open_with_a_skipped_line(dir.path(), "21", unnamed);
        assert_eq!(
            format!("{:#}", source.prewarm("21").unwrap_err()),
            format!(
                "{path}: the end of a record of type exon of chromosome 21 is not a whole \
                 number: xyz"
            )
        );
        assert!(global().take_collected().is_empty());
    }

    /// Nothing else of a skipped line's start and end is read: a start of 0,
    /// an end below its start, trailing text, a start below zero, a blank, a
    /// hexadecimal prefix, a `+`, zero padding, an empty start before its end
    /// and an end that reads zero after its minus sign are columns `tabix`
    /// indexes, and Ensembl VEP 116.3 reads a file holding one with the rows
    /// of every record. Each gives, on either sequence, the features of the
    /// file without the line and the one warning about its type. The same
    /// start on a line the file yields refuses the read.
    #[test]
    fn a_skipped_line_tabix_indexes_is_not_read() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let sound = |sequence: &str| -> Vec<(String, u64, u64)> {
            ENHANCER_SPANS
                .into_iter()
                .map(|(start, end)| (format!("ENSR{sequence}_{start}"), start, end))
                .collect()
        };
        for (start, end) in [
            ("5012000", "5012100"),
            ("0", "5012100"),
            ("5012100", "5012000"),
            ("5012000x", "5012100"),
            ("-5", "5012100"),
            (" 5012000", "5012100"),
            ("0x4C7B40", "5012100"),
            ("+5012000", "5012100"),
            ("005012000", "005012100"),
            ("5012000", "5012100x"),
            ("", "5012100"),
            ("5012000", "-0"),
        ] {
            for chr in ["20", "21"] {
                // A file of its own, so that its warning is no repeat.
                let dir = tempfile::tempdir().unwrap();
                let skipped =
                    format!("{chr}\tEnsembl\texon\t{start}\t{end}\t.\t.\t.\tID=SKIP{chr}");
                let (source, path) = open_with_a_skipped_line(dir.path(), chr, &skipped);
                for sequence in ["20", "21"] {
                    assert_eq!(
                        spans(&source, sequence),
                        sound(sequence),
                        "{start:?} {end:?} on {chr}, read of {sequence}"
                    );
                }
                let warnings: Vec<String> = global()
                    .take_collected()
                    .into_iter()
                    .map(|(_, text)| text)
                    .collect();
                assert_eq!(
                    warnings,
                    [format!(
                        "WARNING: Ignoring 'exon' feature type from {path}; not a supported \
                         regulatory feature type\n"
                    )],
                    "{start:?} {end:?} on {chr}"
                );
            }
        }
        let dir = tempfile::tempdir().unwrap();
        let kept = "21\tEnsembl\tenhancer\t-5\t5012100\t.\t.\t.\tID=ENSR_BELOW_ZERO";
        let (source, path) = open_with_a_skipped_line(dir.path(), "21", kept);
        assert_eq!(
            format!("{:#}", source.prewarm("21").unwrap_err()),
            format!(
                "{path}: the start of enhancer record ENSR_BELOW_ZERO of chromosome 21 is not a \
                 whole number: -5"
            )
        );
    }

    /// A line that opens with a tab, or with white space and a tab, holds
    /// its type where `tabix` reads a start: the index build fails it and
    /// files it among the lines of the sequence before it (htslib 1.9
    /// `tbx.c:96`, `:235-236`). The read of that sequence refuses the line
    /// as the cache builder's bake of the file does, naming the type as the
    /// start and the record by its third column as written, whether the
    /// record's own type is one the file yields or not; Ensembl VEP passes
    /// over the line with exit 0 and writes no row for its record
    /// (`RegFeat.pm:311`). The other sequence is read whole. A line so led
    /// whose type a number leads is one `tabix` indexes, as a sequence of its
    /// own, and is passed over.
    #[test]
    fn a_line_that_opens_with_a_tab_is_refused_for_the_type_where_a_start_is_read() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        for (lead, feature_type) in [("\t", "enhancer"), ("\t", "exon"), (" \t", "enhancer")] {
            for (chr, other) in [("20", "21"), ("21", "20")] {
                let led = format!(
                    "{lead}{chr}\tEnsembl\t{feature_type}\t5012000\t5012100\t.\t.\t.\tID=LED{chr}"
                );
                let (source, path) = open_with_a_skipped_line(dir.path(), chr, &led);
                assert_eq!(
                    format!("{:#}", source.prewarm(chr).unwrap_err()),
                    format!(
                        "{path}: the start of a record of type Ensembl of chromosome {chr} is \
                         not a whole number: {feature_type}"
                    ),
                    "{lead:?} {feature_type} on {chr}"
                );
                assert_eq!(
                    spans(&source, other).len(),
                    3,
                    "{lead:?} {feature_type} on {chr}"
                );
            }
        }
        let led = "\t21\tEnsembl\t5_region\t5012000\t5012100\t.\t.\t.\tID=LED21";
        let (source, _) = open_with_a_skipped_line(dir.path(), "21", led);
        assert_eq!(spans(&source, "21").len(), 3);
        assert!(global().take_collected().is_empty());
    }

    /// A NUL byte behind the end column, or one that ends the type, leaves
    /// `tabix` the record's own start and end and is a byte of the column
    /// it stands in, as it is to Ensembl VEP, whose rows carry it: an
    /// attribute value holds it, a type that ends in one is no accepted
    /// type and is skipped with its warning (`RegFeat.pm:314-322`), and a
    /// score, a strand or a phase that ends in one changes nothing, a
    /// motif's strand so written being neither `+` nor `-` (`:369-370`).
    #[test]
    fn a_nul_byte_in_a_record_is_read_as_ensembl_vep_reads_it() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_indexed_gff_bytes(
            dir.path(),
            "features.gff3.gz",
            &[
                b"21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ENSR\0_INSIDE;color=#faca00",
                b"21\tEnsembl\tenhancer\t2000\t2200\t.\t.\t.\tID=ENSR_VALUE_END\0;color=#faca00",
                b"21\tEnsembl\tenhancer\t3000\t3200\t.\t.\t.\tcolor=#faca00\0;ID=ENSR_OTHER_VALUE",
                b"21\tEnsembl\tenhancer\t4000\t4200\t.\t.\t.\tID=ENSR_COLUMN_END;color=#faca00\0",
                b"21\tEnsembl\tenhancer\0\t5000\t5200\t.\t.\t.\tID=ENSR_TYPE_END",
                b"21\tEnsembl\tenhancer\t6000\t6200\t.\0\t.\t.\tID=ENSR_SCORE_END",
                b"21\tEnsembl\tCTCF_binding_site\t7000\t7200\t.\t-\0\t.\tID=ENSR_STRAND_END",
                b"21\tEnsembl\tenhancer\t8000\t8200\t.\t.\t.\0\tID=ENSR_PHASE_END",
            ],
        );
        let motifs = write_indexed_gff_bytes(
            dir.path(),
            "motifs.gff3.gz",
            &[b"21\tEnsembl\tTF_binding_site\t9000\t9013\t.\t-\0\t.\tID=ENSM_STRAND_END;binding_matrix_id=M;transcription_factor=TF"],
        );
        let source = open(
            &RegulatoryGff {
                file: Some(regulatory.clone()),
                emars: None,
                motifs: Some(motifs),
            },
            false,
        );
        let read: Vec<(&str, u64, i8)> = source
            .features("21")
            .unwrap()
            .iter()
            .map(|feature| (feature.stable_id.as_ref(), feature.start, feature.strand))
            .collect();
        assert_eq!(
            read,
            [
                ("ENSR\0_INSIDE", 1000, 0),
                ("ENSR_VALUE_END\0", 2000, 0),
                ("ENSR_OTHER_VALUE", 3000, 0),
                ("ENSR_COLUMN_END", 4000, 0),
                ("ENSR_SCORE_END", 6000, 0),
                ("ENSR_STRAND_END", 7000, 0),
                ("ENSR_PHASE_END", 8000, 0),
                ("ENSM_STRAND_END", 9000, 0),
            ]
        );
        let warnings: Vec<String> = global()
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .collect();
        assert_eq!(
            warnings,
            [format!("WARNING: Ignoring 'enhancer\0' feature type from {regulatory}; not a supported regulatory feature type\n")]
        );
    }

    /// A start or an end that holds a NUL byte, behind its digits or among
    /// them, is not a whole number, and its record refuses the read as one
    /// with any other text behind its digits does. Ensembl VEP reads the
    /// digits before a byte that ends the cell, with Perl's `isn't numeric`
    /// warning (`RegFeat.pm:311`, `AnnotationSource.pm:306-309`), and writes
    /// no row for a record with the byte among the digits: its index, which
    /// ends a column at the byte (`tbx.c` `tbx_parse1`), returns the record
    /// to no query. A NUL byte that leads the cell leaves `tabix` an empty
    /// column whose reading stops at the byte: the record is refused for
    /// that column, as the cache builder's bake refuses it.
    #[test]
    fn a_nul_byte_in_a_start_or_an_end_is_not_a_whole_number() {
        let dir = tempfile::tempdir().unwrap();
        for (start, end, column, value) in [
            ("5012000\0", "5012100", "start", "5012000\0"),
            ("5012000\x0077", "5012100", "start", "5012000\x0077"),
            ("5012000", "5012100\0", "end", "5012100\0"),
            ("5012000", "5012\x00100", "end", "5012\x00100"),
            ("\x005012000", "5012100", "start", ""),
            ("5012000", "\x005012100", "end", ""),
        ] {
            let kept = format!("21\tEnsembl\tenhancer\t{start}\t{end}\t.\t.\t.\tID=ENSR_NUL");
            let (source, path) = open_with_a_skipped_line(dir.path(), "21", &kept);
            assert_eq!(
                format!("{:#}", source.prewarm("21").unwrap_err()),
                format!(
                    "{path}: the {column} of enhancer record ENSR_NUL of chromosome 21 is not a \
                     whole number: {value}"
                ),
                "{start:?} {end:?}"
            );
        }
    }

    /// `tabix` ends a column at a NUL byte as at a tab (`tbx.c:89`), so one
    /// in a sequence name or a source, or inside a type, puts the type or
    /// its tail where `tabix` reads a start, and the index build fails the
    /// line (`tbx.c:96`). The read of the sequence the line is filed under
    /// refuses it as the cache builder's bake of the file does, naming that
    /// text as the start, whether the record's type is one the file yields
    /// or not; Ensembl VEP reads no such line from its index and writes no
    /// row for its record. The other sequence is read whole. A line whose
    /// every column a NUL byte ends is one `tabix` indexes where its start
    /// and its end are numbers: it names no sequence of the index and is
    /// passed over, and one whose start no number leads is refused.
    #[test]
    fn a_nul_byte_before_the_start_column_refuses_the_read_for_the_start_tabix_reads() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let rest = "5012000\t5012100\t.\t.\t.\tID=ENSR_NUL";
        for (columns, record, start) in [
            (
                "21\0\tEnsembl\tenhancer",
                "enhancer record ENSR_NUL",
                "enhancer",
            ),
            (
                "2\x001\tEnsembl\tenhancer",
                "enhancer record ENSR_NUL",
                "enhancer",
            ),
            (
                "21\tEnsembl\0\tenhancer",
                "enhancer record ENSR_NUL",
                "enhancer",
            ),
            (
                "21\tEns\0embl\tenhancer",
                "enhancer record ENSR_NUL",
                "enhancer",
            ),
            (
                "21\tEnsembl\tenh\0ancer",
                "enh\0ancer record ENSR_NUL",
                "ancer",
            ),
            ("21\0\tEnsembl\texon", "exon record ENSR_NUL", "exon"),
            ("21\tEnsembl\0\texon", "exon record ENSR_NUL", "exon"),
            ("21\tEns\0embl\texon", "exon record ENSR_NUL", "exon"),
            ("21\tEnsembl\tex\0on", "ex\0on record ENSR_NUL", "on"),
        ] {
            let line = format!("{columns}\t{rest}");
            let (source, path) = open_with_a_skipped_line(dir.path(), "21", &line);
            assert_eq!(
                format!("{:#}", source.prewarm("21").unwrap_err()),
                format!(
                    "{path}: the start of {record} of chromosome 21 is not a whole number: {start}"
                ),
                "{columns:?}"
            );
            assert_eq!(spans(&source, "20").len(), 3, "{columns:?}");
        }
        let sound: Vec<(String, u64, u64)> = ENHANCER_SPANS
            .into_iter()
            .map(|(start, end)| (format!("ENSR21_{start}"), start, end))
            .collect();
        for feature_type in ["enhancer", "exon"] {
            let nul_ended = format!("21\tEnsembl\t{feature_type}\t{rest}").replace('\t', "\0");
            let (source, _) = open_with_a_skipped_line(dir.path(), "21", &nul_ended);
            assert_eq!(spans(&source, "21"), sound, "{feature_type}");
            let no_number = nul_ended.replace("5012000", "abc");
            let (source, path) = open_with_a_skipped_line(dir.path(), "21", &no_number);
            // The line has no third column for a type to be read from.
            assert_eq!(
                format!("{:#}", source.prewarm("21").unwrap_err()),
                format!(
                    "{path}: the start of a record of type  of chromosome 21 is not a whole \
                     number: abc"
                ),
                "{feature_type}"
            );
        }
        assert!(global().take_collected().is_empty());
    }

    /// The check of every record refuses what the read of the record's
    /// chromosome refuses, in the read's own words, and loads nothing: a
    /// record without its attributes column, an end that is not a whole
    /// number, a pair that is not ordered, a skipped line `tabix` does not
    /// index. Over a file it refuses nothing of, it draws no warning, and the
    /// read that follows draws the skipped type's warning once.
    #[test]
    fn the_check_of_every_record_refuses_as_a_read_does_and_draws_no_warning() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        for refused in [
            "21\tEnsembl\tenhancer\t5012000\t5012100\t.\t.\t.",
            "21\tEnsembl\tenhancer\t5012000\txyz\t.\t.\t.\tID=ENSR_END",
            "21\tEnsembl\tenhancer\t5012100\t5012000\t.\t.\t.\tID=ENSR_SWAPPED",
            "21\tEnsembl\texon\t5012000\txyz\t.\t.\t.\tID=SKIP",
        ] {
            let (source, path) = open_with_a_skipped_line(dir.path(), "21", refused);
            let checked = format!("{:#}", source.check_gff_records().unwrap_err());
            assert!(checked.starts_with(&format!("{path}: ")), "{checked}");
            assert!(source.features("20").is_none() && source.features("21").is_none());
            assert_eq!(
                checked,
                format!("{:#}", source.prewarm("21").unwrap_err()),
                "{refused}"
            );
        }
        assert!(global().take_collected().is_empty());

        let dir = tempfile::tempdir().unwrap();
        let skipped = "21\tEnsembl\texon\t5012000\t5012100\t.\t.\t.\tID=SKIP";
        let (source, path) = open_with_a_skipped_line(dir.path(), "21", skipped);
        source.check_gff_records().unwrap();
        assert!(global().take_collected().is_empty());
        assert!(source.features("20").is_none() && source.features("21").is_none());
        assert_eq!(spans(&source, "21").len(), 3);
        let warnings: Vec<String> = global()
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .collect();
        assert_eq!(
            warnings,
            [format!(
                "WARNING: Ignoring 'exon' feature type from {path}; not a supported regulatory \
                 feature type\n"
            )]
        );
    }

    /// A value the reader takes (`ID`, `Name`, `binding_matrix_id`,
    /// `transcription_factor`, `extended_start`, `extended_end`) that holds
    /// a byte that is no UTF-8 is read with the replacement character for the
    /// byte, as every input is, and the feature is built: Ensembl VEP carries
    /// the byte into its rows. A bound so read is not a whole number, so it is
    /// ignored with its warning. An escape beside such a byte completes the
    /// character with it ([`decode_attribute`] decodes over the bytes), where
    /// Ensembl VEP, decoding no attribute value (`GFF3.pm:190-198`), carries
    /// the byte and the escape as written.
    #[test]
    fn a_read_attribute_that_is_not_utf8_reads_as_the_replacement_character() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_indexed_gff_bytes(
            dir.path(),
            "features.gff3.gz",
            &[
                b"##gff-version 3",
                b"21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ENSR1;description=caf\xFF",
                b"21\tEnsembl\tpromoter\t2000\t2200\t.\t.\t.\tID=caf\xFF",
                b"21\tEnsembl\tpromoter\t3000\t3200\t.\t.\t.\tID=;Name=n\xFF",
                b"21\tEnsembl\tpromoter\t4000\t4200\t.\t.\t.\tID=p;extended_start=5\xFF",
            ],
        );
        let motifs = write_indexed_gff_bytes(
            dir.path(),
            "motifs.gff3.gz",
            &[b"21\tEnsembl\tTF_binding_site\t5000\t5013\t.\t+\t.\tID=m;binding_matrix_id=caf\xFF;transcription_factor=t\xFF,u"],
        );
        let source = open(
            &RegulatoryGff {
                file: Some(regulatory.clone()),
                emars: None,
                motifs: Some(motifs),
            },
            true,
        );
        let features = source.features("21").unwrap();
        let ids: Vec<(&str, u64)> = features
            .iter()
            .map(|feature| (feature.stable_id.as_ref(), feature.start))
            .collect();
        assert_eq!(
            ids,
            [
                ("ENSR1", 1000),
                ("caf\u{FFFD}", 2000),
                ("n\u{FFFD}", 3000),
                ("p", 4000),
                ("m", 5000)
            ]
        );
        let FeatureKind::Motif { matrix, .. } = &features[4].kind else {
            panic!("the motif file yields a motif");
        };
        assert_eq!(matrix.stable_id.as_ref(), "caf\u{FFFD}");
        assert_eq!(
            matrix.transcription_factors.as_ref(),
            [Arc::from("t\u{FFFD}"), Arc::from("u")]
        );
        let warnings: Vec<String> = global()
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .collect();
        assert_eq!(
            warnings,
            [format!("WARNING: Ignoring extended_start '5\u{FFFD}' of promoter p in {regulatory}; not a whole number\n")]
        );
        let attributes = Attributes::parse(b"ID=caf\xFF;Name=ok;a=caf\xCA%80;b=caf%CA\x80");
        assert_eq!(attributes.get("ID").as_deref(), Some("caf\u{FFFD}"));
        assert_eq!(attributes.get("Name").as_deref(), Some("ok"));
        assert_eq!(attributes.get("a").as_deref(), Some("caf\u{280}"));
        assert_eq!(attributes.get("b").as_deref(), Some("caf\u{280}"));
    }

    /// Every `%XX` escape of a read value is decoded, in either letter case
    /// and whether or not GFF3 reserves the character (`uri_unescape`, the
    /// decoding the format names for its ninth column): `ESC%3bAPE%2c1%3E2`
    /// names the feature `ESC;APE,1>2`, a `transcription_factor` list is split
    /// on its written commas before its members are decoded, so `A%2cB,C` lists
    /// `A,B` and `C`, a `%` not followed by two hexadecimal digits stays
    /// written, and `%252C` reads `%2C`. Ensembl VEP 116.3 prints each value as
    /// written (`GFF3.pm:190-198`), the registered class.
    #[test]
    fn every_escape_of_a_read_value_is_decoded_in_either_case() {
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "##gff-version 3",
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ESC%3bAPE%2c1%3E2",
                "21\tEnsembl\tpromoter\t2000\t2200\t.\t.\t.\tID=;Name=a%2Gb%2;description=x",
                "21\tEnsembl\tpromoter\t3000\t3200\t.\t.\t.\tID=%252C",
            ],
        );
        let motifs = write_indexed_gff(
            dir.path(),
            "motifs.gff3.gz",
            &["21\tEnsembl\tTF_binding_site\t5000\t5013\t.\t+\t.\tID=m;binding_matrix_id=ENSPFM%30;transcription_factor=A%2cB,C"],
        );
        let source = open(
            &RegulatoryGff {
                file: Some(regulatory),
                emars: None,
                motifs: Some(motifs),
            },
            false,
        );
        let features = source.features("21").unwrap();
        let ids: Vec<&str> = features
            .iter()
            .map(|feature| feature.stable_id.as_ref())
            .collect();
        assert_eq!(ids, ["ESC;APE,1>2", "a%2Gb%2", "%2C", "m"]);
        let FeatureKind::Motif { matrix, .. } = &features[3].kind else {
            panic!("the motif file yields a motif");
        };
        assert_eq!(matrix.stable_id.as_ref(), "ENSPFM0");
        assert_eq!(
            matrix.transcription_factors.as_ref(),
            [Arc::from("A,B"), Arc::from("C")]
        );
    }

    /// A coordinate-derived identifier carries the chromosome as the file
    /// spells it, whatever name the cache or the input gives the chromosome:
    /// Ensembl VEP maps the input's name onto the file's own spelling before
    /// the region is built (`AnnotationSource.pm:179,198`) and that region
    /// name is the `$chr` of `_record_get_id` (`RegFeat.pm:372,583`). A
    /// `chr21` file read under `21` or `chr21` names the record `chr21:...`,
    /// a `21` file read under either names it `21:...`.
    #[test]
    fn coordinate_identifier_carries_the_files_chromosome_spelling() {
        let dir = tempfile::tempdir().unwrap();
        let record = |chr: &str| {
            format!("{chr}\tEnsembl\tCTCF_binding_site\t3000\t3020\t.\t-\t.\tcolor=#40e0d0")
        };
        let prefixed = write_indexed_gff(dir.path(), "chr21.gff3.gz", &[&record("chr21")]);
        let bare = write_indexed_gff(dir.path(), "21.gff3.gz", &[&record("21")]);
        let ids_under = |path: &str, read_as: &str| -> Vec<String> {
            let cache_dir = tempfile::tempdir().unwrap();
            let source = RegulatorySource::open_gff(
                &file_spec(path.to_string()),
                false,
                false,
                &cache_without_regulatory(cache_dir.path()),
            )
            .unwrap();
            source.prewarm(read_as).unwrap();
            source
                .features(read_as)
                .unwrap()
                .iter()
                .map(|feature| feature.stable_id.to_string())
                .collect()
        };
        for read_as in ["21", "chr21"] {
            assert_eq!(
                ids_under(&prefixed, read_as),
                ["chr21:3000-3020:CTCF_binding_site"],
                "the chr21 file read under {read_as}"
            );
            assert_eq!(
                ids_under(&bare, read_as),
                ["21:3000-3020:CTCF_binding_site"],
                "the 21 file read under {read_as}"
            );
        }
        let source = open(&file_spec(prefixed), false);
        assert_eq!(
            row_summary(&source, 3010),
            [(
                "chr21:3000-3020:CTCF_binding_site".to_string(),
                FeatureType::RegulatoryFeature,
                Some("CTCF_binding_site".to_string()),
                0,
                vec![Consequence::RegulatoryRegionVariant]
            )]
        );
    }

    /// The type column is read as `get_type` reads it (`GFF3.pm:161-164`,
    /// `uri_unescape`), so a `%XX` escape in it is decoded before the
    /// include-list check (`RegFeat.pm:312-314`): `CTCF%5Fbinding%5Fsite` is
    /// a `CTCF_binding_site` feature with that BIOTYPE, `TF%5Fbinding%5Fsite`
    /// a motif, a `promot%65r` is widened and named by the decoded type, and
    /// a type unsupported once decoded is warned about by its decoded text.
    #[test]
    fn an_escaped_type_column_is_decoded_before_the_include_list_check() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=PLAIN_TYPE",
                "21\tEnsembl\tex%6Fn\t2000\t2100\t.\t.\t.\tID=exon1",
                "21\tEnsembl\tCTCF%5Fbinding%5Fsite\t3000\t3020\t.\t+\t.\tID=ESCAPED_TYPE",
                "21\tEnsembl\tpromot%65r\t4000\t4200\t.\t.\t.\textended_start=3900",
                "21\tEnsembl\tpromoter%5Fflanking%5Fregion\t5000\t5100\t.\t.\t.\tcolor=#faca00",
            ],
        );
        let motifs = write_indexed_gff(
            dir.path(),
            "motifs.gff3.gz",
            &["21\tEnsembl\tTF%5Fbinding%5Fsite\t6000\t6013\t.\t-\t.\tID=ESCAPED_MOTIF;binding_matrix_id=M;transcription_factor=TF"],
        );
        let source = open(
            &RegulatoryGff {
                file: Some(regulatory.clone()),
                emars: None,
                motifs: Some(motifs),
            },
            true,
        );
        let summary: Vec<(&str, u64, Option<&str>)> = source
            .features("21")
            .unwrap()
            .iter()
            .map(|feature| {
                let biotype = match &feature.kind {
                    FeatureKind::Regulatory { biotype, .. } => biotype.as_deref(),
                    FeatureKind::Motif { .. } => None,
                };
                (feature.stable_id.as_ref(), feature.start, biotype)
            })
            .collect();
        assert_eq!(
            summary,
            [
                ("PLAIN_TYPE", 1000, Some("enhancer")),
                ("ESCAPED_TYPE", 3000, Some("CTCF_binding_site")),
                ("21:3900-4200:promoter", 3900, Some("promoter")),
                (
                    "21:5000-5100:promoter_flanking_region",
                    5000,
                    Some("promoter_flanking_region")
                ),
                ("ESCAPED_MOTIF", 6000, None),
            ]
        );
        assert_eq!(
            row_summary(&source, 3010),
            [(
                "ESCAPED_TYPE".to_string(),
                FeatureType::RegulatoryFeature,
                Some("CTCF_binding_site".to_string()),
                0,
                vec![Consequence::RegulatoryRegionVariant]
            )]
        );
        let warnings: Vec<String> = global()
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .collect();
        assert_eq!(
            warnings,
            [format!("WARNING: Ignoring 'exon' feature type from {regulatory}; not a supported regulatory feature type\n")]
        );
    }

    /// A type column that is not UTF-8, whether an escape spells the byte
    /// (`caf%FF`) or it stands in the column, is a type like any other once
    /// decoded: it matches no accepted type, so the record is skipped with the
    /// warning naming the type as decoded, the replacement character standing
    /// for the byte, and the run goes on, as Ensembl VEP skips the record with
    /// its warning carrying the byte (`RegFeat.pm:314-322`); an escape decoding
    /// to text that is UTF-8 is warned about by that text, an escape that
    /// completes a raw first byte included.
    #[test]
    fn a_type_column_that_is_not_utf8_is_skipped_with_its_warning() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff_bytes(
            dir.path(),
            "features.gff3.gz",
            &[
                b"##gff-version 3",
                b"21\tEnsembl\tcaf%FF\t1000\t1100\t.\t.\t.\tID=ESCAPED",
                b"21\tEnsembl\tcaf\xFF\t2000\t2100\t.\t.\t.\tID=RAW",
                b"21\tEnsembl\tcaf%CA%80\t3000\t3100\t.\t.\t.\tID=UNSUPPORTED",
                b"21\tEnsembl\tmix\xCA%80\t3500\t3600\t.\t.\t.\tID=COMPLETED",
                b"21\tEnsembl\tmix%CA\x80\t3700\t3800\t.\t.\t.\tID=COMPLETED_TOO",
                b"21\tEnsembl\tenhancer\t4000\t4200\t.\t.\t.\tID=KEPT",
            ],
        );
        let source = open(&file_spec(path.clone()), false);
        let ids: Vec<&str> = source
            .features("21")
            .unwrap()
            .iter()
            .map(|feature| feature.stable_id.as_ref())
            .collect();
        assert_eq!(ids, ["KEPT"]);
        let warnings: Vec<String> = global()
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .collect();
        assert_eq!(
            warnings,
            [
                format!("WARNING: Ignoring 'caf\u{FFFD}' feature type from {path}; not a supported regulatory feature type\n"),
                format!("WARNING: Ignoring 'caf\u{280}' feature type from {path}; not a supported regulatory feature type\n"),
                format!("WARNING: Ignoring 'mix\u{280}' feature type from {path}; not a supported regulatory feature type\n"),
            ]
        );
    }

    /// A feature type the file does not yield is warned about once per type
    /// and file with Ensembl VEP's text (`RegFeat.pm:318-322`) and skipped; a
    /// motif file yields `TF_binding_site` alone (`:118-120`), a regulatory
    /// file the eight types of `:94-103`, EMAR among them.
    #[test]
    fn unsupported_type_is_warned_once() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let regulatory = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "21\tEnsembl\texon\t1000\t1200\t.\t.\t.\tID=exon1",
                "21\tEnsembl\tenhancer\t2000\t2200\t.\t.\t.\tID=ENSR_KEPT",
                "21\tEnsembl\texon\t3000\t3020\t.\t+\t.\tID=exon2",
                "21\tEnsembl\tEMAR\t4000\t4200\t.\t.\t.\tID=EMAR_KEPT",
            ],
        );
        let motifs = write_indexed_gff(
            dir.path(),
            "motifs.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ENSR_NOT_A_MOTIF",
                "21\tEnsembl\tTF_binding_site\t2000\t2010\t.\t+\t.\tID=ENSM_KEPT;binding_matrix_id=M;transcription_factor=TF",
            ],
        );
        let source = open(
            &RegulatoryGff {
                file: Some(regulatory.clone()),
                emars: None,
                motifs: Some(motifs.clone()),
            },
            false,
        );
        let ids: Vec<&str> = source
            .features("21")
            .unwrap()
            .iter()
            .map(|feature| feature.stable_id.as_ref())
            .collect();
        assert_eq!(ids, ["ENSR_KEPT", "ENSM_KEPT", "EMAR_KEPT"]);
        let warnings: Vec<String> = global()
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .collect();
        assert_eq!(
            warnings,
            [
                format!("WARNING: Ignoring 'exon' feature type from {regulatory}; not a supported regulatory feature type\n"),
                format!("WARNING: Ignoring 'enhancer' feature type from {motifs}; not a supported regulatory feature type\n"),
            ]
        );
    }

    /// An EMAR record is a regulatory feature whose type, `EMAR`, prints as
    /// its BIOTYPE, under its own identifier beside the regulatory feature of
    /// the same span; `emars=` alone is a valid source.
    #[test]
    fn emar_rows_carry_biotype_emar() {
        let both = open(&fixture_spec(), false);
        assert_eq!(
            row_summary(&both, 5011400),
            [
                (
                    "ENSR00001793149".to_string(),
                    FeatureType::RegulatoryFeature,
                    Some("EMAR".to_string()),
                    0,
                    vec![Consequence::RegulatoryRegionVariant]
                ),
                (
                    "ENSR21_85ZF5".to_string(),
                    FeatureType::RegulatoryFeature,
                    Some("enhancer".to_string()),
                    0,
                    vec![Consequence::RegulatoryRegionVariant]
                ),
            ]
        );
        let emars_alone = open(
            &RegulatoryGff {
                emars: Some(fixture("EMARs.gff.gz")),
                ..Default::default()
            },
            false,
        );
        assert_eq!(emars_alone.features("21").unwrap().len(), 279);
        assert_eq!(
            row_summary(&emars_alone, 5011400)
                .into_iter()
                .map(|row| row.2)
                .collect::<Vec<_>>(),
            [Some("EMAR".to_string())]
        );
    }

    /// A structural variant's regulatory and motif rows come in the order the
    /// sources' features were met, `file=` records then `emars=` records then
    /// `motifs=` records, each file in its own order: Ensembl VEP builds one
    /// source per key in that order (`AnnotationSourceAdaptor.pm:284-310`),
    /// pushes a structural overlap per feature met
    /// (`AnnotationType/RegFeat.pm:118-135`) and reads the lists as pushed
    /// (`StructuralVariationFeature.pm:495-509`). Each such row carries its
    /// feature's met rank as `id_rank`; a sequence variant's regulatory rows
    /// carry none, Ensembl VEP sorting those by stable id
    /// (`VariationFeature.pm:710`).
    #[test]
    fn structural_rows_come_in_the_sources_met_order() {
        let source = open(&fixture_spec(), false);
        let structural = |start: u64, end: u64| {
            let mut v = variant(start, end, b"N", b"<DEL>");
            v.is_structural = true;
            v.variant_class = VariantClass::StructuralDeletion;
            v.sv_end = Some(end);
            source.annotate(&mut v);
            v
        };
        let written = |v: &InputVariant| -> Vec<(String, FeatureType, u32)> {
            vep_io::output::fields::record_rows(&[v])
                .iter()
                .filter_map(|(_, row)| {
                    row.map(|row| {
                        (
                            row.transcript_id().to_string(),
                            row.feature_type,
                            row.id_rank,
                        )
                    })
                })
                .collect()
        };
        // The enhancer and the EMAR share one span, where a start-then-identifier
        // order would list the EMAR first (`ENSR0` sorts before `ENSR2`).
        let over_enhancer = structural(5_010_240, 5_011_600);
        let rows = written(&over_enhancer);
        let ids: Vec<&str> = rows.iter().map(|(id, _, _)| id.as_str()).collect();
        assert_eq!(
            ids,
            [
                "ENSR21_85ZF5",
                "ENSR00001793149",
                "ENSM00000574652",
                "ENSM00000574653",
                "ENSM00000574654",
                "ENSM00000574655"
            ]
        );
        assert!(rows.iter().all(|(_, _, rank)| *rank != 0), "{rows:?}");
        assert!(
            rows.windows(2).all(|pair| pair[0].2 < pair[1].2),
            "{rows:?}"
        );
        assert_eq!(rows[0].2, 1, "the regulatory file's first record");
        // The EMAR file's records follow every record of the regulatory file,
        // the CTCF site inside the enhancer's EMAR among them.
        let over_emar = written(&structural(5_013_560, 5_014_350));
        let regulatory: Vec<&str> = over_emar
            .iter()
            .filter(|(_, kind, _)| *kind == FeatureType::RegulatoryFeature)
            .map(|(id, _, _)| id.as_str())
            .collect();
        assert_eq!(
            regulatory,
            ["ENSR21_BZTD", "ENSR21_536SRH", "ENSR00001793150"]
        );
        let mut sequence = variant(5_011_400, 5_011_400, b"A", b"C");
        source.annotate(&mut sequence);
        assert_eq!(sequence.transcript_consequences.len(), 2);
        assert!(sequence
            .transcript_consequences
            .iter()
            .all(|row| row.id_rank == 0 && row.feature_type == FeatureType::RegulatoryFeature));
    }

    /// `--extended_promoters` widens a promoter to its `extended_start` and
    /// `extended_end` (`RegFeat.pm:360-365`): the promoter `ENSR21_B34QG`
    /// (5553169-5553746, extended from 5553154) gains a base in the extension
    /// under the flag and not without it; an enhancer carrying the attributes
    /// is never widened (`:106-108`); the attributes never narrow a feature.
    #[test]
    fn extended_promoters_widen_promoters_only() {
        let spec = RegulatoryGff {
            file: Some(fixture("regulatory_features.gff3.gz")),
            ..Default::default()
        };
        let core = open(&spec, false);
        let extended = open(&spec, true);
        let promoter = |source: &RegulatorySource| {
            source
                .features("21")
                .unwrap()
                .iter()
                .find(|feature| feature.stable_id.as_ref() == "ENSR21_B34QG")
                .map(|feature| (feature.start, feature.end))
                .unwrap()
        };
        assert_eq!(promoter(&core), (5_553_169, 5_553_746));
        assert_eq!(promoter(&extended), (5_553_154, 5_553_746));
        assert!(row_summary(&core, 5_553_160).is_empty());
        assert_eq!(row_summary(&extended, 5_553_160)[0].0, "ENSR21_B34QG");
        assert_eq!(
            row_summary(&extended, 5_553_153),
            row_summary(&core, 5_553_153)
        );
        let widened: Vec<&str> = extended
            .features("21")
            .unwrap()
            .iter()
            .zip(core.features("21").unwrap())
            .filter(|(wide, narrow)| (wide.start, wide.end) != (narrow.start, narrow.end))
            .map(|(wide, _)| wide.stable_id.as_ref())
            .collect();
        assert_eq!(widened.len(), 25, "{widened:?}");

        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ENH;extended_start=900;extended_end=1300",
                "21\tEnsembl\tpromoter\t2000\t2200\t.\t.\t.\tID=NARROWER;extended_start=2050;extended_end=2150",
            ],
        );
        let spans: Vec<(u64, u64)> = open(&file_spec(path), true)
            .features("21")
            .unwrap()
            .iter()
            .map(|feature| (feature.start, feature.end))
            .collect();
        assert_eq!(spans, [(1000, 1200), (2000, 2200)]);
    }

    /// A promoter's `extended_start` or `extended_end` whose value is not a
    /// whole number is ignored, the promoter keeping its own bound while the
    /// record's other bound still widens it, and is warned about once, naming
    /// the record and the attribute, however many times the file is read.
    /// Ensembl VEP 116.3 compares the value as a number (`RegFeat.pm:363-364`),
    /// so `x` reads as 0 there and the promoter starts at the chromosome's
    /// first base. Without the flag no bound is read and nothing is warned,
    /// and a type other than promoter never reads the attributes (`:106-108`).
    #[test]
    fn non_integer_extended_bound_is_ignored_with_one_warning() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "21\tEnsembl\tpromoter\t3000\t3200\t.\t.\t.\tID=NOT_A_NUMBER;extended_start=x;extended_end=3300",
                "21\tEnsembl\tpromoter\t4000\t4200\t.\t.\t.\textended_start=3900;extended_end=1.5",
                "21\tEnsembl\tenhancer\t5000\t5200\t.\t.\t.\tID=ENH;extended_start=x",
            ],
        );
        let spans = |extended: bool| -> Vec<(String, u64, u64)> {
            open(&file_spec(path.clone()), extended)
                .features("21")
                .unwrap()
                .iter()
                .map(|feature| (feature.stable_id.to_string(), feature.start, feature.end))
                .collect()
        };
        let warnings = || -> Vec<String> {
            global()
                .take_collected()
                .into_iter()
                .map(|(_, text)| text)
                .collect()
        };
        let widened = [
            ("NOT_A_NUMBER".to_string(), 3000, 3300),
            ("21:3900-4200:promoter".to_string(), 3900, 4200),
            ("ENH".to_string(), 5000, 5200),
        ];
        assert_eq!(spans(true), widened);
        assert_eq!(spans(true), widened, "a second read of the file");
        assert_eq!(
            warnings(),
            [
                format!("WARNING: Ignoring extended_start 'x' of promoter NOT_A_NUMBER in {path}; not a whole number\n"),
                format!("WARNING: Ignoring extended_end '1.5' of promoter 21:3900-4200:promoter in {path}; not a whole number\n"),
            ]
        );
        assert_eq!(
            spans(false),
            [
                ("NOT_A_NUMBER".to_string(), 3000, 3200),
                ("21:4000-4200:promoter".to_string(), 4000, 4200),
                ("ENH".to_string(), 5000, 5200),
            ]
        );
        assert!(warnings().is_empty(), "no bound is read without the flag");
    }

    /// A motif's matrix is a stub of its `binding_matrix_id`, span and
    /// `transcription_factor` list (`RegFeat.pm:431-444`): the row carries
    /// MOTIF_NAME, MOTIF_POS along the motif's strand and the factors, and no
    /// HIGH_INF_POS or MOTIF_SCORE_CHANGE (`OutputFactory.pm:1944-1949`); motifs
    /// naming one matrix share one stub.
    #[test]
    fn motif_stub_omits_high_inf_pos_and_score_change() {
        let source = open(
            &RegulatoryGff {
                motifs: Some(fixture("motif_features.gff3.gz")),
                ..Default::default()
            },
            false,
        );
        let mut v = variant(5_010_250, 5_010_250, b"A", b"C");
        source.annotate(&mut v);
        let mut rows: Vec<_> = v
            .transcript_consequences
            .iter()
            .map(|row| {
                let motif = row.motif.as_deref().unwrap();
                (
                    row.transcript_id().to_string(),
                    row.strand,
                    motif.name.to_string(),
                    motif.pos,
                    motif.high_inf_pos,
                    motif.score_change,
                    motif
                        .factors
                        .iter()
                        .map(|f| f.to_string())
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        rows.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(
            rows,
            [
                (
                    "ENSM00000574652".to_string(),
                    -1,
                    "ENSPFM0085".to_string(),
                    Some(14),
                    None,
                    None,
                    vec!["EBF1".to_string()]
                ),
                (
                    "ENSM00000574653".to_string(),
                    1,
                    "ENSPFM0085".to_string(),
                    Some(1),
                    None,
                    None,
                    vec!["EBF1".to_string()]
                ),
            ]
        );
        assert_eq!(
            v.most_severe_consequence,
            Some(Consequence::TfBindingSiteVariant)
        );
        let features = source.features("21").unwrap();
        let stubs: Vec<&Arc<BindingMatrix>> = features
            .iter()
            .filter_map(|feature| match &feature.kind {
                FeatureKind::Motif { matrix, .. } if matrix.stable_id.as_ref() == "ENSPFM0085" => {
                    Some(matrix)
                }
                _ => None,
            })
            .collect();
        assert!(stubs.len() > 1 && stubs.iter().all(|stub| Arc::ptr_eq(stub, stubs[0])));
        assert!(!stubs[0].has_frequencies() && stubs[0].length == 14);
        let three = features
            .iter()
            .find(|feature| feature.stable_id.as_ref() == "ENSM00000574814")
            .map(|feature| match &feature.kind {
                FeatureKind::Motif { matrix, .. } => matrix
                    .transcription_factors
                    .iter()
                    .map(|f| f.to_string())
                    .collect::<Vec<_>>(),
                _ => unreachable!(),
            })
            .unwrap();
        assert_eq!(three, ["HNF4A", "NR2F1", "RXRA"]);
        assert!(features.iter().all(|feature| feature.dbid.is_none()));
    }

    /// The warnings about `--regulatory_gff` and `--cell_type` the sink holds:
    /// a test that opens a cache without the guard writes its own warning to
    /// the process's sink meanwhile.
    fn gff_warnings() -> Vec<String> {
        global()
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .filter(|text| text.contains("regulatory_gff") || text.contains("cell_type"))
            .collect()
    }

    /// Under `--regulatory_gff` the GFF files replace the cache's regulatory
    /// tree, with the warning of `CacheDir.pm:161-164` when the cache holds
    /// one and none when it does not; `--cell_type` draws the source's
    /// warning once for any number of regulatory files and not for a motif
    /// file alone, and every row then carries the empty `CELL_TYPE` list; a
    /// missing file, a file without an index, a remote URL, a file whose
    /// index cannot be read and one whose index carries no tabix header are
    /// refused by key and name, an unreadable index with no word of a plugin.
    #[test]
    fn gff_wins_over_cache_shards_with_the_warning() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        write_fixture_cache(dir.path(), &[]);
        let cache = resolved_cache(dir.path());
        let spec = fixture_spec();
        let source = RegulatorySource::open_gff(&spec, false, false, &cache).unwrap();
        source.prewarm("21").unwrap();
        let ids: FxHashSet<&str> = source
            .features("21")
            .unwrap()
            .iter()
            .map(|feature| feature.stable_id.as_ref())
            .collect();
        assert!(!ids.contains("ENSR_FIXTURE") && ids.contains("ENSR21_85ZF5"));
        assert!(
            rows_at(&source, 106) == "[]",
            "the cache's enhancer at 100-120 is ignored"
        );
        assert_eq!(
            gff_warnings(),
            [format!("WARNING: {GFF_REPLACES_CACHE_WARNING}\n")]
        );

        let bare = tempfile::tempdir().unwrap();
        let without = cache_without_regulatory(bare.path());
        let with_cell_type = RegulatorySource::open_gff(&spec, false, true, &without).unwrap();
        // No cache source is built under the flag, so no name is checked
        // against a list (`CacheDir.pm:156-176`), and every regulatory and
        // motif row carries CELL_TYPE as the empty list (OutputFactory.pm:1888,
        // 1936, 1968-1978: a GFF feature has no activity map).
        assert!(with_cell_type.reads_gff());
        let cell_types = super::super::CellTypes::resolve(Some("HUVEC"), Some(&with_cell_type))
            .unwrap()
            .expect("the flag is on");
        with_cell_type.prewarm("21").unwrap();
        let mut v = variant(5011400, 5011400, b"A", b"C");
        with_cell_type.annotate_rows(&mut v, Some(&cell_types));
        assert!(!v.transcript_consequences.is_empty());
        for row in &v.transcript_consequences {
            assert!(
                matches!(row.cell_type.as_deref(), Some([])),
                "{}: {:?}",
                row.transcript_id(),
                row.cell_type
            );
        }
        RegulatorySource::open_gff(
            &RegulatoryGff {
                motifs: spec.motifs.clone(),
                ..Default::default()
            },
            false,
            true,
            &without,
        )
        .unwrap();
        assert_eq!(gff_warnings(), [format!("WARNING: {CELL_TYPE_WARNING}\n")]);

        let refused = |spec: RegulatoryGff| {
            RegulatorySource::open_gff(&spec, false, false, &without)
                .err()
                .unwrap()
                .to_string()
        };
        let missing = bare
            .path()
            .join("missing.gff3.gz")
            .to_string_lossy()
            .into_owned();
        assert_eq!(
            refused(RegulatoryGff {
                emars: Some(missing.clone()),
                ..Default::default()
            }),
            format!("ERROR: [cli] --regulatory_gff emars={missing} not found")
        );
        let unindexed = bare.path().join("unindexed.gff3.gz");
        fs::copy(fixture("EMARs.gff.gz"), &unindexed).unwrap();
        assert_eq!(
            refused(file_spec(unindexed.to_string_lossy().into_owned())),
            format!(
                "ERROR: [cli] --regulatory_gff file={} has no tabix index",
                unindexed.display()
            )
        );
        assert_eq!(
            refused(RegulatoryGff {
                motifs: Some("https://example.org/motifs.gff3.gz".into()),
                ..Default::default()
            }),
            "ERROR: [cli] --regulatory_gff motifs=https://example.org/motifs.gff3.gz is a remote URL; vep-rs reads a local file alone"
        );

        // The index reader's errors display with a sentence on a plugin before
        // their text (`PluginError`); the refusal carries the text alone.
        let unreadable = bare.path().join("unreadable.gff3.gz");
        fs::copy(fixture("EMARs.gff.gz"), &unreadable).unwrap();
        let unreadable_index = format!("{}.tbi", unreadable.display());
        fs::write(&unreadable_index, b"").unwrap();
        assert_eq!(
            refused(file_spec(unreadable.to_string_lossy().into_owned())),
            format!(
                "ERROR: [cli] --regulatory_gff file={}: tabix index {unreadable_index} is malformed: the index ends inside its magic",
                unreadable.display()
            )
        );
        fs::write(&unreadable_index, b"not a compressed index").unwrap();
        let message = refused(file_spec(unreadable.to_string_lossy().into_owned()));
        assert!(
            message.starts_with(&format!(
                "ERROR: [cli] --regulatory_gff file={}: failed to read tabix index {unreadable_index}: ",
                unreadable.display()
            )) && !message.contains("plugin"),
            "{message}"
        );

        let headerless = bare.path().join("headerless.gff3.gz");
        fs::copy(fixture("EMARs.gff.gz"), &headerless).unwrap();
        let headerless_index = format!("{}.csi", headerless.display());
        noodles::csi::fs::write(&headerless_index, &noodles::csi::Index::default()).unwrap();
        assert_eq!(
            refused(RegulatoryGff {
                emars: Some(headerless.to_string_lossy().into_owned()),
                ..Default::default()
            }),
            format!(
                "ERROR: [cli] --regulatory_gff emars={}: index {headerless_index} carries no tabix header, so its records cannot be mapped to chromosomes; index the file with tabix",
                headerless.display()
            )
        );
    }

    /// A chromosome is read under the file's own spelling whatever name the
    /// cache gives it, through the `chr` prefix and the cache's synonyms
    /// (`RegFeat.pm:296`, `BaseVEP.pm:689-738`); a record of another
    /// chromosome inside a shared block is not read; a record the file holds
    /// twice is indexed once (`:628-640`), and a distinct record under the
    /// same identifier is indexed too, with Ensembl VEP 116.3's collision
    /// warning once for the identifier (`RegFeat.pm:485-494`; the byte-identical
    /// copy shares its digest and draws none, `:487`).
    #[test]
    fn chromosome_names_resolve_and_duplicate_records_collapse() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "chr21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ONCE",
                "chr21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ONCE",
                "chr21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ONCE;color=#000000",
                "chr22\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=OTHER",
            ],
        );
        let source = open(&file_spec(path.clone()), false);
        let ids: Vec<&str> = source
            .features("21")
            .unwrap()
            .iter()
            .map(|feature| feature.stable_id.as_ref())
            .collect();
        assert_eq!(
            ids,
            ["ONCE", "ONCE"],
            "two distinct records sharing an id are both indexed; a variant inside one alone keeps that record's row"
        );
        assert_eq!(
            duplicate_identifier_warnings(),
            [format!(
                "WARNING: Duplicate regulatory feature identifier 'ONCE' in {path}; identifiers must be unique or annotations will be lost\n"
            )]
        );
        source.prewarm("chr22").unwrap();
        assert!(duplicate_identifier_warnings().is_empty());
        assert_eq!(source.features("chr22").unwrap().len(), 1);
        assert_eq!(source.features("22").unwrap().len(), 1);
        assert!(source.features("23").is_none());

        // A chromosome two files spell differently is entered under each
        // spelling and read once per entry: the second read meets the first's
        // identifiers again and is no collision.
        let emars = write_indexed_gff(
            dir.path(),
            "EMARs.gff.gz",
            &["22\tEnsembl\tEMAR\t1000\t1200\t.\t.\t.\tID=EMAR_OTHER"],
        );
        global().configure(Target::Memory, true, false);
        let two_spellings = open(
            &RegulatoryGff {
                file: Some(path.clone()),
                emars: Some(emars),
                motifs: None,
            },
            false,
        );
        two_spellings.prewarm("chr22").unwrap();
        two_spellings.prewarm("22").unwrap();
        assert_eq!(two_spellings.features("chr22").unwrap().len(), 2);
        assert_eq!(two_spellings.features("22").unwrap().len(), 2);
        assert_eq!(
            duplicate_identifier_warnings(),
            [format!(
                "WARNING: Duplicate regulatory feature identifier 'ONCE' in {path}; identifiers must be unique or annotations will be lost\n"
            )],
            "the chr21 collision alone, once; a re-read of chromosome 22 adds none"
        );

        // The duplicate check runs on the start column the file is sorted by: a
        // promoter widened under --extended_promoters keeps its place among the
        // records of its own start, so its second copy still meets the first.
        let widened = write_indexed_gff(
            dir.path(),
            "widened.gff3.gz",
            &[
                "21\tEnsembl\tpromoter\t1000\t1200\t.\t.\t.\tID=P;extended_start=900",
                "21\tEnsembl\tenhancer\t1000\t1100\t.\t.\t.\tID=E",
                "21\tEnsembl\tpromoter\t1000\t1200\t.\t.\t.\tID=P;extended_start=900",
            ],
        );
        for extended in [false, true] {
            let source = open(&file_spec(widened.clone()), extended);
            let spans: Vec<(&str, u64)> = source
                .features("21")
                .unwrap()
                .iter()
                .map(|feature| (feature.stable_id.as_ref(), feature.start))
                .collect();
            let promoter_start = if extended { 900 } else { 1000 };
            assert_eq!(spans, [("P", promoter_start), ("E", 1000)], "{extended}");
        }
        assert!(
            duplicate_identifier_warnings().is_empty(),
            "a byte-identical copy of a record is no collision"
        );
    }

    /// The collision warnings the sink holds.
    fn duplicate_identifier_warnings() -> Vec<String> {
        global()
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .filter(|text| text.contains("Duplicate regulatory feature identifier"))
            .collect()
    }

    /// A cache whose chromosome 21 holds the regulatory features
    /// `(source, biotype, start, end)`, every one under the identifier `X`,
    /// in the given order, each in the 1 Mb shard of its start with the
    /// source and its place in the list baked in, as `vep-cache-builder`
    /// bakes a file's records; loaded.
    fn baked_cache(
        features: &[(u8, &str, u64, u64)],
    ) -> (tempfile::TempDir, Arc<RegulatorySource>) {
        let dir = tempfile::tempdir().unwrap();
        let chr = dir.path().join("regulatory/21");
        fs::create_dir_all(&chr).unwrap();
        fs::write(dir.path().join("regulatory/binding_matrices.json"), b"{}").unwrap();
        let mut shards: BTreeMap<u64, Vec<serde_json::Value>> = BTreeMap::new();
        for (ordinal, (source, biotype, start, end)) in features.iter().enumerate() {
            let region = (start - 1) / 1_000_000 * 1_000_000 + 1;
            shards.entry(region).or_default().push(serde_json::json!({
                "feature_type": "RegulatoryFeature",
                "stable_id": "X",
                "start": start,
                "end": end,
                "biotype": biotype,
                "activity": "",
                "source": source,
                "ordinal": ordinal,
            }));
        }
        for (region, shard) in shards {
            fs::write(
                chr.join(format!("{region}-{}.json", region + 999_999)),
                serde_json::to_vec(&shard).unwrap(),
            )
            .unwrap();
        }
        let mut cache = resolved_cache(dir.path());
        cache.info.cell_types = Vec::new();
        let source = RegulatorySource::open(&cache, true).unwrap().unwrap();
        source.prewarm("21").unwrap();
        (dir, source)
    }

    /// Two distinct records of one file sharing an identifier are one feature
    /// to a sequence variant: inside both it gets one `RegulatoryFeature` row,
    /// the record later in the file (Ensembl VEP 116.3 keys the variant's
    /// regulatory overlaps on the stable id and the later feature's write
    /// wins, `VariationFeature.pm:913-919`, the features arriving in file
    /// order), inside one alone that record's row; a structural variant over
    /// both keeps both rows (`AnnotationType/RegFeat.pm:128-136` lists its
    /// overlaps). The same records baked into a cache's shard give the same
    /// rows. Across two files the identifiers are separate (`RegFeat.pm:485`
    /// keeps one set per source): a `file=` enhancer and an `emars=` EMAR under
    /// one identifier are two rows and no warning.
    #[test]
    fn distinct_records_sharing_an_identifier_give_a_variant_the_later_row() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=X",
                "21\tEnsembl\tpromoter\t1150\t1300\t.\t.\t.\tID=X",
            ],
        );
        let source = open(&file_spec(path.clone()), false);
        let biotypes = |source: &RegulatorySource, position: u64| -> Vec<Option<String>> {
            row_summary(source, position)
                .into_iter()
                .map(|(id, feature_type, biotype, _, _)| {
                    assert_eq!(
                        (id.as_str(), feature_type),
                        ("X", FeatureType::RegulatoryFeature)
                    );
                    biotype
                })
                .collect()
        };
        assert_eq!(biotypes(&source, 1100), [Some("enhancer".to_string())]);
        assert_eq!(biotypes(&source, 1180), [Some("promoter".to_string())]);
        assert_eq!(biotypes(&source, 1250), [Some("promoter".to_string())]);
        assert_eq!(
            duplicate_identifier_warnings(),
            [format!(
                "WARNING: Duplicate regulatory feature identifier 'X' in {path}; identifiers must be unique or annotations will be lost\n"
            )]
        );
        let mut deletion = variant(1000, 1300, b"A", b"<DEL>");
        deletion.is_structural = true;
        deletion.variant_class = VariantClass::StructuralDeletion;
        source.annotate(&mut deletion);
        let structural_rows = deletion
            .transcript_consequences
            .iter()
            .filter(|row| row.transcript_id() == "X")
            .count();
        assert_eq!(structural_rows, 2, "a structural variant keeps both rows");
        // A breakend with its own position in the earlier record and its mate
        // in the later one: two features, two rows
        // (`AnnotationType/RegFeat.pm:122-135` builds one overlap per feature a
        // breakend position lies in).
        let mut breakend = variant(1100, 1100, b"N", b"N[21:1250[");
        breakend.is_structural = true;
        breakend.variant_class = VariantClass::Translocation;
        breakend.sv_end = Some(1100);
        breakend.mate_id = Some("mate".into());
        breakend.mate_chr = Some("21".into());
        breakend.mate_pos = Some(1250);
        source.annotate(&mut breakend);
        let mut breakend_rows: Vec<(u64, u64)> = breakend
            .transcript_consequences
            .iter()
            .filter(|row| row.transcript_id() == "X")
            .map(|row| (row.feature_start, row.feature_end))
            .collect();
        breakend_rows.sort();
        assert_eq!(breakend_rows, [(1000, 1200), (1150, 1300)]);

        // The two records as `vep-cache-builder` bakes them into a shard, in
        // file order.
        let (_baked, from_cache) =
            baked_cache(&[(0, "enhancer", 1000, 1200), (0, "promoter", 1150, 1300)]);
        assert_eq!(from_cache.features("21").unwrap().len(), 2);
        for position in [1100, 1180, 1250] {
            assert_eq!(
                rows_at(&from_cache, position),
                rows_at(&source, position),
                "at 21:{position}"
            );
        }
        // Two records of one span differing in type are two features on both
        // routes, and the later record's row.
        let same_span = write_indexed_gff(
            dir.path(),
            "same_span.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=X",
                "21\tEnsembl\tCTCF_binding_site\t1000\t1200\t.\t.\t.\tID=X",
            ],
        );
        let same_span = open(&file_spec(same_span), false);
        let (_baked, same_span_cache) = baked_cache(&[
            (0, "enhancer", 1000, 1200),
            (0, "CTCF_binding_site", 1000, 1200),
        ]);
        assert_eq!(same_span_cache.features("21").unwrap().len(), 2);
        assert_eq!(
            biotypes(&same_span, 1100),
            [Some("CTCF_binding_site".to_string())]
        );
        assert_eq!(rows_at(&same_span_cache, 1100), rows_at(&same_span, 1100));
        assert_eq!(duplicate_identifier_warnings().len(), 1);

        let emars = write_indexed_gff(
            dir.path(),
            "EMARs.gff.gz",
            &["21\tEnsembl\tEMAR\t1000\t1200\t.\t.\t.\tID=X"],
        );
        let enhancer = write_indexed_gff(
            dir.path(),
            "enhancer.gff3.gz",
            &["21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=X"],
        );
        let two_files = open(
            &RegulatoryGff {
                file: Some(enhancer),
                emars: Some(emars),
                motifs: None,
            },
            false,
        );
        assert_eq!(
            biotypes(&two_files, 1100),
            [Some("EMAR".to_string()), Some("enhancer".to_string())]
        );
        assert!(
            duplicate_identifier_warnings().is_empty(),
            "one set of identifiers per file"
        );
    }

    /// The record later in the file keeps the row even when the index's sort
    /// by start puts it first: a promoter widened under `--extended_promoters`
    /// to a start below an earlier record's sorts before that record, while
    /// Ensembl VEP 116.3 reads the records in file order whatever the widening
    /// (`RegFeat.pm:312-328`) and the later write wins
    /// (`VariationFeature.pm:918`), so the promoter's row is the one written,
    /// on the file route and from a baked cache alike: each baked element
    /// carries its place in the file ([`Feature::ordinal`]), whichever shard
    /// its widened start put it in and however the shards' names sort.
    #[test]
    fn widened_promoter_sharing_an_identifier_keeps_the_file_order() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t5000\t5500\t.\t.\t.\tID=X",
                "21\tEnsembl\tpromoter\t5200\t5800\t.\t.\t.\tID=X;extended_start=4900;extended_end=5900",
            ],
        );
        for extended in [false, true] {
            let source = open(&file_spec(path.clone()), extended);
            let rows = row_summary(&source, 5300);
            assert_eq!(rows.len(), 1, "{extended}");
            assert_eq!(rows[0].2, Some("promoter".to_string()), "{extended}");
        }
        let (_baked, from_cache) =
            baked_cache(&[(0, "enhancer", 5000, 5500), (0, "promoter", 4900, 5900)]);
        assert_eq!(
            row_summary(&from_cache, 5300)[0].2,
            Some("promoter".to_string())
        );
        assert_eq!(duplicate_identifier_warnings().len(), 1);
        // Two shards whose byte order (`10000001-` before `9000001-`) is not
        // their numeric order still read in the file's order.
        let (_baked, two_shards) = baked_cache(&[
            (0, "enhancer", 9_999_500, 10_000_500),
            (0, "promoter", 10_000_100, 10_001_000),
        ]);
        assert_eq!(two_shards.features("21").unwrap().len(), 2);
        assert_eq!(
            row_summary(&two_shards, 10_000_200)[0].2,
            Some("promoter".to_string())
        );
        // A later promoter widened across a region boundary is baked into the
        // earlier shard; its baked ordinal still places it after the record it
        // follows in the file.
        let (_baked, earlier_shard) = baked_cache(&[
            (0, "enhancer", 20_000_001, 20_000_300),
            (0, "promoter", 19_999_900, 20_000_500),
        ]);
        assert_eq!(earlier_shard.features("21").unwrap().len(), 2);
        assert_eq!(
            row_summary(&earlier_shard, 20_000_100)[0].2,
            Some("promoter".to_string())
        );
    }

    /// A baked cache keeps what the file route keeps: two records of two input
    /// files sharing an identifier are two features and two rows, since an
    /// identifier is unique within its file alone (`RegFeat.pm:485`, one set
    /// per source); two distinct records of one file equal in every baked
    /// field but their place in the file are two features, so a structural
    /// variant over both keeps both rows (`AnnotationType/RegFeat.pm:128-136`
    /// lists its overlaps) and a sequence variant gets the later record's one
    /// row (`VariationFeature.pm:913-919`).
    #[test]
    fn baked_records_keep_their_file_and_place() {
        let (_baked, two_sources) =
            baked_cache(&[(0, "enhancer", 1000, 1200), (1, "EMAR", 1000, 1200)]);
        assert_eq!(two_sources.features("21").unwrap().len(), 2);
        let mut biotypes: Vec<Option<String>> = row_summary(&two_sources, 1100)
            .into_iter()
            .map(|row| row.2)
            .collect();
        biotypes.sort();
        assert_eq!(
            biotypes,
            [Some("EMAR".to_string()), Some("enhancer".to_string())]
        );

        let (_baked, twins) =
            baked_cache(&[(0, "enhancer", 1000, 1200), (0, "enhancer", 1000, 1200)]);
        assert_eq!(twins.features("21").unwrap().len(), 2);
        assert_eq!(row_summary(&twins, 1100).len(), 1);
        let mut deletion = variant(1000, 1200, b"A", b"<DEL>");
        deletion.is_structural = true;
        deletion.variant_class = VariantClass::StructuralDeletion;
        twins.annotate(&mut deletion);
        assert_eq!(
            deletion
                .transcript_consequences
                .iter()
                .filter(|row| row.transcript_id() == "X")
                .count(),
            2
        );
    }

    /// A structural variant's rows from a baked tree come in the order the
    /// files' records were met, the `source` slot then the `ordinal` of each
    /// element, as the run over the files writes them; a converted tree's
    /// elements carry neither and take their shard order, here the same order.
    /// A sequence variant's regulatory rows carry no rank from either tree.
    #[test]
    fn baked_tree_orders_a_structural_variants_rows_as_met() {
        // A regulatory file's enhancer (slot 0) and the EMAR file's region
        // (slot 1) that starts before it.
        let bake = |numbered: bool| -> (tempfile::TempDir, Arc<RegulatorySource>) {
            let dir = tempfile::tempdir().unwrap();
            let chr = dir.path().join("regulatory/21");
            fs::create_dir_all(&chr).unwrap();
            fs::write(dir.path().join("regulatory/binding_matrices.json"), b"{}").unwrap();
            let mut shard = Vec::new();
            for (source, id, start, end) in
                [(0, "ENSR21_B", 1000, 1200), (1, "ENSR00000001", 900, 1300)]
            {
                let mut element = serde_json::json!({
                    "feature_type": "RegulatoryFeature",
                    "stable_id": id,
                    "start": start,
                    "end": end,
                    "biotype": "enhancer",
                    "activity": "",
                });
                if numbered {
                    element["source"] = serde_json::json!(source);
                    element["ordinal"] = serde_json::json!(0);
                }
                shard.push(element);
            }
            fs::write(
                chr.join("1-1000000.json"),
                serde_json::to_vec(&shard).unwrap(),
            )
            .unwrap();
            let mut cache = resolved_cache(dir.path());
            cache.info.cell_types = Vec::new();
            let source = RegulatorySource::open(&cache, true).unwrap().unwrap();
            source.prewarm("21").unwrap();
            (dir, source)
        };
        let written = |source: &RegulatorySource| -> Vec<(String, u32)> {
            let mut deletion = variant(800, 1400, b"N", b"<DEL>");
            deletion.is_structural = true;
            deletion.variant_class = VariantClass::StructuralDeletion;
            deletion.sv_end = Some(1400);
            source.annotate(&mut deletion);
            vep_io::output::fields::record_rows(&[&deletion])
                .iter()
                .filter_map(|(_, row)| {
                    row.map(|row| (row.transcript_id().to_string(), row.id_rank))
                })
                .collect()
        };
        let (_baked, baked) = bake(true);
        assert_eq!(
            written(&baked),
            [("ENSR21_B".to_string(), 1), ("ENSR00000001".to_string(), 2)]
        );
        let (_converted, converted) = bake(false);
        assert_eq!(
            written(&converted),
            [("ENSR21_B".to_string(), 1), ("ENSR00000001".to_string(), 2)]
        );
        for source in [&baked, &converted] {
            let mut sequence = variant(1100, 1100, b"A", b"C");
            source.annotate(&mut sequence);
            assert_eq!(sequence.transcript_consequences.len(), 2);
            assert!(sequence
                .transcript_consequences
                .iter()
                .all(|row| row.id_rank == 0));
        }
    }

    /// An identifier repeated on two chromosomes of one file is a collision
    /// too, warned about once whichever chromosome is read first: Ensembl VEP
    /// 116.3 keeps one identifier set per source over the run
    /// (`RegFeat.pm:485`), and both features stay indexed.
    #[test]
    fn identifier_repeated_across_chromosomes_warns_once() {
        let _guard = sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff(
            dir.path(),
            "features.gff3.gz",
            &[
                "21\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ACROSS",
                "22\tEnsembl\tenhancer\t1000\t1200\t.\t.\t.\tID=ACROSS",
            ],
        );
        let cache_dir = tempfile::tempdir().unwrap();
        let cache = cache_without_regulatory(cache_dir.path());
        for order in [["21", "22"], ["22", "21"]] {
            global().configure(Target::Memory, true, false);
            // Opened without the helper's prewarm of chromosome 21, so the
            // chromosomes are read in `order`.
            let source =
                RegulatorySource::open_gff(&file_spec(path.clone()), false, false, &cache).unwrap();
            for chr in order {
                source.prewarm(chr).unwrap();
                assert_eq!(source.features(chr).unwrap().len(), 1, "{chr}");
            }
            assert_eq!(
                duplicate_identifier_warnings(),
                [format!(
                    "WARNING: Duplicate regulatory feature identifier 'ACROSS' in {path}; identifiers must be unique or annotations will be lost\n"
                )],
                "{order:?}"
            );
        }
    }

    /// Two distinct motif records sharing an identifier draw the collision
    /// warning once and keep both rows: Ensembl VEP 116.3 keys a variant's
    /// motif overlaps on the per-record database id, not the stable id
    /// (`MotifFeatureVariation.pm:70-76`, `RegFeat.pm:452`), while
    /// `_check_id_collision` runs for a motif source too (`RegFeat.pm:382`).
    #[test]
    fn distinct_motif_records_sharing_an_identifier_keep_both_rows_and_warn_once() {
        let _guard = sink_test_guard();
        global().configure(Target::Memory, true, false);
        let dir = tempfile::tempdir().unwrap();
        let path = write_indexed_gff(
            dir.path(),
            "motifs.gff3.gz",
            &[
                "21\tEnsembl\tTF_binding_site\t2000\t2010\t.\t+\t.\tID=M;binding_matrix_id=PFM;transcription_factor=TF",
                "21\tEnsembl\tTF_binding_site\t2005\t2015\t.\t-\t.\tID=M;binding_matrix_id=PFM;transcription_factor=TF",
            ],
        );
        let source = open(
            &RegulatoryGff {
                motifs: Some(path.clone()),
                ..Default::default()
            },
            false,
        );
        let rows: Vec<(String, FeatureType, i8)> = row_summary(&source, 2008)
            .into_iter()
            .map(|(id, feature_type, _, strand, _)| (id, feature_type, strand))
            .collect();
        assert_eq!(
            rows,
            [
                ("M".to_string(), FeatureType::MotifFeature, -1),
                ("M".to_string(), FeatureType::MotifFeature, 1),
            ]
        );
        assert_eq!(
            duplicate_identifier_warnings(),
            [format!(
                "WARNING: Duplicate regulatory feature identifier 'M' in {path}; identifiers must be unique or annotations will be lost\n"
            )]
        );
    }
}
