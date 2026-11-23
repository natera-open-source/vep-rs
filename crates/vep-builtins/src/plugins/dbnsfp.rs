// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! dbNSFP plugin: the columns the user selects from a tabix-indexed dbNSFP file,
//! on the transcript rows of a single-base substitution, as Ensembl VEP's
//! `dbNSFP.pm` writes them.
//!
//! **Usage:** `--plugin dbNSFP,/path/to/dbNSFP4.9a_grch38.gz[,/path/to/replacement_logic],col1,col2,...`
//! with any of `consequence=`, `output_prefix=`, `pep_match=` and
//! `transcript_match=` anywhere in the list (dbNSFP.pm:192, 586).
//!
//! | Parameter | Meaning |
//! | --- | --- |
//! | first positional | the bgzip-compressed, tabix-indexed dbNSFP file |
//! | second positional, when it names an existing file | a replacement-logic file (`colname from to` lines, `default` rows required) in place of the built-in `;` to `,` and `\|` to `&` (dbNSFP.pm:220-226, 545-566) |
//! | remaining positionals | the columns to write; `ALL` selects every column of the header (dbNSFP.pm:233-236); a name the header lacks prints `invalid_field` (dbNSFP.pm:239-249, 403-405) |
//! | `consequence=` | `ALL`, or SO terms joined by `&`; by default a row is annotated when its consequences include `missense_variant`, `stop_lost`, `stop_gained` or `start_lost` (dbNSFP.pm:157, 290-292, 591-597) |
//! | `output_prefix=` | truthy prefixes every output name with `dbNSFP_` (dbNSFP.pm:159, 257-261, 582) |
//! | `pep_match=` | on by default: a record's `aaref`/`aaalt` (`X` read as `*`) must equal the row's amino acid change (dbNSFP.pm:340-344, 583, 605-611) |
//! | `transcript_match=` | truthy keeps, in a `;`-joined transcript-specific column, the entry of the row's transcript alone (dbNSFP.pm:350-390, 584) |
//!
//! A record matches on position and the alternate allele, on the forward
//! strand (dbNSFP.pm:299-303, 335-338); only a single-base substitution whose
//! allele is `A`, `C`, `G` or `T` is annotated (dbNSFP.pm:297, 303), and a
//! record on `MT` or `chrMT` is looked up under `M` (dbNSFP.pm:313): on GRCh37,
//! where `M` names another sequence, only in a table whose rows state the
//! reference bases, and under `MT` otherwise. The first
//! matching record of a position is the row's record (dbNSFP.pm:391). A `.`
//! value is not written (dbNSFP.pm:408).
//!
//! The data file's name must carry a dbNSFP version, read from the path as given
//! exactly as dbNSFP.pm:197-211 reads it (`2.9` or `4.0b1` anywhere, else the
//! first `[3-9].N` or `[3-9].N.N` run, an `a` or `c` suffix preferred): a name
//! with none fails the plugin's initialisation, as the constructor's `die` does,
//! and a bare version without the suffix is accepted with the constructor's
//! warning on standard error (dbNSFP.pm:208).
//!
//! The file is read as the plugin's documentation prepares it: for GRCh38
//! sorted and indexed on `#chr` and `pos(1-based)`, for GRCh37 on `hg19_chr` and
//! `hg19_pos(1-based)` (dbNSFP.pm:69-74, 316-334). The position column is
//! resolved by name for the run's assembly and checked against the column the
//! tabix index was built on: the two preparations of one release carry the
//! same header and differ in their index, so a file indexed for the other
//! assembly would otherwise match nothing and exit 0. A file named as dbNSFP
//! 2.9 serves a GRCh37 run by its `pos(1-coor)` column (dbNSFP.pm:317-319) and
//! a file named as 4.0b1 a GRCh38 run by the same column beside
//! `hg19_pos(1-based)` (dbNSFP.pm:329-330); the Ensembl VEP plugin tells the
//! two apart by the version it read
//! from the file name (dbNSFP.pm:317, 329), this module by the columns the
//! header carries, so the version decides acceptance alone, never the layout.
//!
//! Transcript-specific columns are found as the Ensembl VEP plugin finds them: the
//! first record in `1:1008170-1082927` with a `;` in a value names every
//! `;`-bearing column outside the fixed list of columns whose lists follow no
//! transcript (dbNSFP.pm:161-181, 444-483). The region is read under the bare
//! chromosome name `1` alone, as the plugin's `tabix` query reads it
//! (dbNSFP.pm:449), so a file naming the chromosome `chr1` yields no record
//! there and `transcript_match=1` is refused, while its rows are still looked
//! up under either spelling. A README beside the data file, any entry whose
//! name carries `dbnsfp` and then `readme.txt` in any letter case
//! (dbNSFP.pm:499), adds the columns its descriptions call "separated by" or
//! tie to an "Ensembl gene id" (dbNSFP.pm:469-476) and supplies the header
//! descriptions, prefixed `(from <file base name>)` (dbNSFP.pm:268-280,
//! 485-543). Of several such entries the plugin reads the first its directory
//! listing returns; this module reads the first in name order.
//!
//! Two-phase plugin: [`BuiltinPlugin::prefetch`] batches the tabix query for a
//! buffer and [`BuiltinPlugin::annotate`] writes the per-transcript results.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::path::{Path, PathBuf};

use bstr::ByteSlice;
use indexmap::IndexMap;
use tracing::{debug, warn};

use crate::annotation_store::{check_table_or_store, open_best_store_with, AnnotationStore};
use crate::params::{Grammar, PluginParams};
use crate::tabix::{
    allele_match, indexed_position, parse_header, split_field_count, table_contig,
    BatchQueryResult, IndexedPosition, TabixAnnotator, TabixConfig, TabixRecord,
};
use crate::traits::{print_on_stderr, warn_on_stderr, BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::variant::InputVariant;

const NAME: &str = "dbNSFP";
const OUTPUT_PREFIX: &str = "dbNSFP_";
const NAMED_PARAMS: &[&str] = &[
    "consequence",
    "output_prefix",
    "pep_match",
    "transcript_match",
];
/// The consequences annotated without `consequence=` (dbNSFP.pm:157).
const DEFAULT_CONSEQUENCES: &[Consequence] = &[
    Consequence::MissenseVariant,
    Consequence::StopLost,
    Consequence::StopGained,
    Consequence::StartLost,
];
/// The region whose first `;`-bearing record names the transcript-specific
/// columns (dbNSFP.pm:162).
const DETECTION_REGION: (&str, u64, u64) = ("1", 1_008_170, 1_082_927);
/// `;`-separated columns whose entries follow no transcript, or not the
/// transcripts of `Ensembl_transcriptid` (dbNSFP.pm:166-181).
const NON_TRANSCRIPT_SPECIFIC: &[&str] = &[
    "MutationTaster_score",
    "MutationTaster_converted_rankscore",
    "MutationTaster_pred",
    "MutationTaster_model",
    "Interpro_domain",
    "MutPred_Top5features",
    "Transcript_id_VEST3",
    "Transcript_var_VEST3",
    "VEST3_score",
    "dbNSFP_POPMAX_AC",
    "dbNSFP_POPMAX_POP",
    "AllofUs_POPMAX_AC",
    "AllofUs_POPMAX_AN",
    "AllofUs_POPMAX_POP",
];
/// The value written for a selected column the header lacks (dbNSFP.pm:404).
const INVALID_FIELD: &str = "invalid_field";

/// A selected column: its header name, its output name and its position in a
/// record, `None` for a name the header lacks.
#[derive(Debug, Clone)]
struct SelectedColumn {
    name: String,
    output_name: String,
    index: Option<usize>,
    transcript_specific: bool,
}

/// Why the detection region yielded no transcript-specific columns, the
/// reason `transcript_match=1` is refused (dbNSFP.pm:613-615).
#[derive(Debug)]
enum DetectionFailure {
    /// No record of the region under the bare chromosome name carries a
    /// `;`-separated value.
    NoListRecord,
    /// The file spells the region's chromosome otherwise, `chr1` for `1`, so
    /// the query under the bare name reads none of its records.
    ChromosomeSpelled(String),
}

impl std::fmt::Display for DetectionFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (chr, start, end) = DETECTION_REGION;
        match self {
            Self::NoListRecord => write!(
                f,
                "no record of {chr}:{start}-{end} holds a `;`-separated value; run without \
                 transcript_match"
            ),
            Self::ChromosomeSpelled(spelling) => write!(
                f,
                "the file names chromosome {chr} `{spelling}`, and the detection query \
                 {chr}:{start}-{end} reads the bare name alone, as the Ensembl VEP plugin's \
                 does; run without transcript_match, or name the file's chromosomes without \
                 the prefix, `{chr}` for `{spelling}`"
            ),
        }
    }
}

/// Which rows the plugin annotates: every row, or the rows whose consequences
/// include one of a set of terms.
#[derive(Debug, Clone)]
enum ConsequenceFilter {
    All,
    Terms(Vec<Consequence>),
}

impl ConsequenceFilter {
    fn passes(&self, consequence: &TranscriptConsequence) -> bool {
        match self {
            Self::All => true,
            Self::Terms(terms) => consequence.consequences.iter().any(|c| terms.contains(c)),
        }
    }
}

/// dbNSFP multi-column annotation plugin.
pub struct DbNsfpPlugin {
    /// Annotation store, initialized in `init()`.
    store: Option<Box<dyn AnnotationStore>>,
    /// The file's column names.
    header: Vec<String>,
    /// The columns to write, in the order given (`ALL` in header order).
    selected: Vec<SelectedColumn>,
    /// Column the index names the chromosome by; a record carries its value in
    /// `chr`, not among its columns.
    chr_col: usize,
    /// Column of the assembly's position (dbNSFP.pm:316-334), the one the
    /// index is built on.
    pos_col: usize,
    /// Column of the alternate allele.
    alt_col: usize,
    /// Column of `Ensembl_transcriptid`.
    transcript_id_col: usize,
    /// Columns of the reference and alternate amino acids, read under `pep_match`.
    aaref_col: Option<usize>,
    aaalt_col: Option<usize>,
    filter: ConsequenceFilter,
    output_prefix: bool,
    pep_match: bool,
    transcript_match: bool,
    /// The `from, to` pairs applied to every value in order, and the pairs that
    /// replace them for one column.
    default_replacements: Vec<(String, String)>,
    column_replacements: HashMap<String, Vec<(String, String)>>,
    /// Header descriptions from the README beside the file, by column name.
    descriptions: HashMap<String, String>,
    /// The base name of the data file, `.gz` dropped, named in a README description.
    basename: String,
    assembly: Option<Assembly>,
    refseq_or_merged: bool,
}

impl DbNsfpPlugin {
    /// Create a new uninitialized dbNSFP plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            header: Vec::new(),
            selected: Vec::new(),
            chr_col: 0,
            pos_col: 0,
            alt_col: 0,
            transcript_id_col: 0,
            aaref_col: None,
            aaalt_col: None,
            filter: ConsequenceFilter::Terms(DEFAULT_CONSEQUENCES.to_vec()),
            output_prefix: false,
            pep_match: true,
            transcript_match: false,
            default_replacements: Vec::new(),
            column_replacements: HashMap::new(),
            descriptions: HashMap::new(),
            basename: String::new(),
            assembly: None,
            refseq_or_merged: false,
        }
    }

    fn column(&self, name: &str) -> Option<usize> {
        self.header.iter().position(|h| h == name)
    }

    /// The output name of a column: prefixed `dbNSFP_` under `output_prefix`
    /// unless it already is (dbNSFP.pm:257-261).
    fn output_name(&self, column: &str) -> String {
        if !self.output_prefix || column.starts_with(OUTPUT_PREFIX) {
            column.to_string()
        } else {
            format!("{OUTPUT_PREFIX}{column}")
        }
    }

    /// The `from, to` pairs a column's value goes through (dbNSFP.pm:410-413).
    fn replacements(&self, column: &str) -> &[(String, String)] {
        self.column_replacements
            .get(column)
            .map_or(&self.default_replacements[..], Vec::as_slice)
    }

    /// Whether `record` is the row's record: at the row's position, with the
    /// row's alternate allele and, under `pep_match`, its amino acid change
    /// (dbNSFP.pm:335-344). The position is the record's cell as
    /// `$pos == $vf->{start}` (`:336`) numifies it
    /// ([`allele_match::perl_position`]: `30396091.0` and `30396091x` are
    /// 30396091, `30396091.5` is no position).
    fn record_matches(
        &self,
        record: &TabixRecord,
        position: u64,
        alt: &[u8],
        amino_acids: Option<&str>,
    ) -> bool {
        if allele_match::perl_position(record, self.pos_col) != Some(position)
            || record.columns.get(self.alt_col).map(String::as_bytes) != Some(alt)
        {
            return false;
        }
        if !self.pep_match {
            return true;
        }
        let Some(row_change) = amino_acids else {
            return false;
        };
        let cell = |col: Option<usize>| {
            col.and_then(|c| record.columns.get(c))
                .map_or("", String::as_str)
        };
        let star_for_x = |byte: u8| if byte == b'X' { b'*' } else { byte };
        cell(self.aaref_col)
            .bytes()
            .chain(std::iter::once(b'/'))
            .chain(cell(self.aaalt_col).bytes())
            .map(star_for_x)
            .eq(row_change.bytes())
    }

    /// The selected columns of `record` for a row of transcript `transcript_id`.
    fn fields(&self, record: &TabixRecord, transcript_id: &str) -> IndexMap<String, String> {
        // The record's `Ensembl_transcriptid` entries and the row's place among
        // them are read under `transcript_match` alone (dbNSFP.pm:349-361).
        let transcript_ids: Vec<&str> = if self.transcript_match {
            record
                .columns
                .get(self.transcript_id_col)
                .map_or_else(Vec::new, |ids| split_entries(ids))
        } else {
            Vec::new()
        };
        let transcript_index = transcript_ids.iter().position(|id| *id == transcript_id);
        let mut fields = IndexMap::with_capacity(self.selected.len());
        // A cell past the row's last filled one is undefined and skipped
        // (dbNSFP.pm:428-431, :407); an empty cell before it is the empty
        // value, written as every value but `.` is (:408).
        let present = split_field_count(&record.columns);
        for column in &self.selected {
            let Some(column_index) = column.index else {
                fields.insert(column.output_name.clone(), INVALID_FIELD.to_string());
                continue;
            };
            let raw = if column_index == self.chr_col {
                &record.chr
            } else {
                match record.columns.get(column_index) {
                    Some(raw) if column_index < present => raw,
                    _ => continue,
                }
            };
            if raw == "." {
                continue;
            }
            let mut value = raw.as_str();
            if self.transcript_match && column.transcript_specific {
                let Some(transcript_index) = transcript_index else {
                    continue;
                };
                // A malformed list leaves the whole value and is reported on
                // every row, a bare `warn` whose text begins `ERROR:`
                // (dbNSFP.pm:378-385).
                let entries = split_entries(raw);
                if transcript_index >= entries.len() {
                    print_on_stderr(&format!(
                        "ERROR: Transcript index out of range for field {}",
                        column.name
                    ));
                } else if entries.len() != transcript_ids.len() {
                    print_on_stderr(&format!(
                        "ERROR: Number of transcript IDs does not match number of data entries for field {}",
                        column.name
                    ));
                } else {
                    value = entries[transcript_index];
                }
            }
            // The transcript's own entry is tested against `.` as the whole
            // value was (:388, then :408).
            if value == "." {
                continue;
            }
            let mut value = Cow::Borrowed(value);
            for (from, to) in self.replacements(&column.name) {
                if value.contains(from.as_str()) {
                    value = Cow::Owned(value.replace(from.as_str(), to));
                }
            }
            fields.insert(column.output_name.clone(), value.into_owned());
        }
        fields
    }

    /// The fields of one transcript row from `records`, the file's records at
    /// the row's position: nothing for a row the filter leaves out, a record
    /// outside the plugin's shape, or no matching record.
    fn annotate_row(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
        records: &[&TabixRecord],
    ) -> IndexMap<String, String> {
        if !self.filter.passes(consequence) || variant.start != variant.end {
            return IndexMap::new();
        }
        let (_, alt) = variant.forward_strand_alleles();
        if !matches!(alt.as_ref(), b"A" | b"C" | b"G" | b"T") {
            return IndexMap::new();
        }
        let Some(record) = records.iter().find(|record| {
            self.record_matches(
                record,
                variant.start,
                &alt,
                consequence.amino_acids.as_deref(),
            )
        }) else {
            return IndexMap::new();
        };
        self.fields(record, consequence.transcript_id())
    }

    /// Resolve the coordinate columns for the run's assembly and check the
    /// index was built on them.
    fn resolve_coordinates(&self, path: &Path) -> Result<(usize, usize), PluginError> {
        let hg19_pos_col = self.column("hg19_pos(1-based)");
        let (chr_col, pos_col) = match self.assembly {
            Some(Assembly::Grch37) => {
                let pos_col = hg19_pos_col
                    .or_else(|| self.column("pos(1-coor)"))
                    .ok_or_else(|| {
                        PluginError::Init(format!(
                            "{NAME}: {} has no hg19_pos(1-based) column (nor a pos(1-coor) column) to use with GRCh37",
                            path.display()
                        ))
                    })?;
                let chr_col = if hg19_pos_col.is_some() {
                    self.column("hg19_chr").or_else(|| self.column("chr"))
                } else {
                    self.column("chr")
                };
                (chr_col.unwrap_or(0), pos_col)
            }
            _ => {
                let mut pos_col = self.column("pos(1-based)");
                if pos_col.is_none() && hg19_pos_col.is_some() {
                    pos_col = self.column("pos(1-coor)");
                }
                let pos_col = pos_col.ok_or_else(|| {
                    PluginError::Init(format!(
                        "{NAME}: {} has no pos(1-based) column to use with GRCh38",
                        path.display()
                    ))
                })?;
                (self.column("chr").unwrap_or(0), pos_col)
            }
        };
        let pos_name = self.header[pos_col].as_str();
        let assembly = self
            .assembly
            .map_or_else(|| "GRCh38".to_string(), |a| a.to_string());
        if self.assembly.is_none() {
            warn!("{NAME}: no assembly given; reading the GRCh38 position column '{pos_name}'");
        }
        match indexed_position(path, pos_col)? {
            IndexedPosition::Unchecked => Ok((chr_col, pos_col)),
            IndexedPosition::Agrees { sequence_column } => Ok((sequence_column, pos_col)),
            IndexedPosition::Disagrees {
                indexed_column: indexed,
            } => {
                let indexed_name = self.header.get(indexed).map_or_else(
                    || format!("column {}", indexed + 1),
                    |name| format!("'{name}'"),
                );
                Err(PluginError::Init(format!(
                    "{NAME}: a {assembly} run reads '{pos_name}' (column {}), but {} is indexed on {indexed_name} (column {}); \
                     prepare and index the file for {assembly} as the plugin's documentation describes, or point --plugin {NAME} \
                     at the copy indexed for it",
                    pos_col + 1,
                    path.display(),
                    indexed + 1,
                )))
            }
        }
    }

    /// The transcript-specific columns: every `;`-bearing column of the first
    /// record in the detection region with a `;` in a value, outside the fixed
    /// list, plus the columns the README calls per-transcript; a
    /// [`DetectionFailure`] when the region holds no such record, whatever the
    /// README names, as the Ensembl VEP plugin reads the README inside the
    /// record loop (dbNSFP.pm:460-480). The region is read under the bare name
    /// `1` alone, as the plugin's `tabix` query reads it (dbNSFP.pm:449): the
    /// store resolves the name under the file's own spelling, so a record found
    /// under `chr1` is one that query does not see.
    fn detect_transcript_specific(
        &self,
        store: &dyn AnnotationStore,
    ) -> Result<Result<HashSet<String>, DetectionFailure>, PluginError> {
        let (chr, start, end) = DETECTION_REGION;
        let records = store.query(chr, start, end)?;
        let Some(record) = records
            .iter()
            .filter(|record| record.chr == chr)
            .find(|record| record.columns.iter().any(|cell| cell.contains(';')))
        else {
            let failure = match records.iter().find(|record| record.chr != chr) {
                Some(other) => DetectionFailure::ChromosomeSpelled(other.chr.clone()),
                None => DetectionFailure::NoListRecord,
            };
            return Ok(Err(failure));
        };
        let mut columns: HashSet<String> = self
            .header
            .iter()
            .enumerate()
            .filter(|(i, name)| {
                record
                    .columns
                    .get(*i)
                    .is_some_and(|cell| cell.contains(';'))
                    && !NON_TRANSCRIPT_SPECIFIC.contains(&name.as_str())
            })
            .map(|(_, name)| name.clone())
            .collect();
        for (name, description) in &self.descriptions {
            if !NON_TRANSCRIPT_SPECIFIC.contains(&name.as_str())
                && (contains_ignore_case(description, "separated by")
                    || contains_ignore_case(description, "Ensembl gene id"))
            {
                columns.insert(name.clone());
            }
        }
        Ok(Ok(columns))
    }

    /// The replacement pairs of a logic file (dbNSFP.pm:552-566): `colname from
    /// to` lines, a `colname` header line skipped, `default` rows required. A line
    /// is split as Perl's `split /\s+/` splits it: on runs of ASCII blanks, a run at
    /// the start opening an empty first value, so a line that opens with a blank
    /// fails the three-value test (dbNSFP.pm:557-559) as a blank line, a short line
    /// and a `0` value do. The file is read as bytes, as Ensembl VEP's plugin reads it,
    /// a non-UTF-8 value standing as the replacement character.
    fn read_replacement_logic(&mut self, path: &Path) -> Result<(), PluginError> {
        let bytes = std::fs::read(path).map_err(|e| {
            PluginError::Init(format!(
                "{NAME}: failed to read replacement logic file {}: {e}",
                path.display()
            ))
        })?;
        let text = String::from_utf8_lossy(&bytes);
        let mut by_column: HashMap<String, Vec<(String, String)>> = HashMap::new();
        for (index, line) in text.lines().enumerate() {
            if line.starts_with("colname") {
                continue;
            }
            let empty_first_field = line.starts_with(is_perl_blank_char).then_some("");
            let mut fields = empty_first_field.into_iter().chain(
                line.split(is_perl_blank_char)
                    .filter(|field| !field.is_empty()),
            );
            match (fields.next(), fields.next(), fields.next()) {
                (Some(column), Some(from), Some(to))
                    if perl_true(column) && perl_true(from) && perl_true(to) =>
                {
                    by_column
                        .entry(column.to_string())
                        .or_default()
                        .push((from.to_string(), to.to_string()));
                }
                _ => {
                    return Err(PluginError::Init(format!(
                        "{NAME}: replacement logic file {} line {} requires 3 values separated by whitespace in this format: colname from to",
                        path.display(),
                        index + 1
                    )));
                }
            }
        }
        let Some(default) = by_column.remove("default") else {
            return Err(PluginError::Init(format!(
                "{NAME}: no default replacement logic has been specified in {}",
                path.display()
            )));
        };
        self.default_replacements = default;
        self.column_replacements = by_column;
        Ok(())
    }

    /// The named parameters (dbNSFP.pm:578-620); every value is read as the
    /// Ensembl VEP plugin reads a Perl truth value: `0` and the empty value are false.
    fn read_named_params(&mut self, params: &PluginParams) -> Result<(), PluginError> {
        if let Some(value) = params.get("consequence") {
            self.filter = if value.eq_ignore_ascii_case("ALL") {
                ConsequenceFilter::All
            } else {
                ConsequenceFilter::Terms(
                    value
                        .split('&')
                        .filter_map(Consequence::from_so_term)
                        .collect(),
                )
            };
        }
        if let Some(value) = params.get("output_prefix") {
            self.output_prefix = perl_true(value);
        }
        if let Some(value) = params.get("pep_match") {
            self.pep_match = perl_true(value);
        }
        if let Some(value) = params.get("transcript_match") {
            self.transcript_match = perl_true(value);
        }
        if self.pep_match {
            for name in ["aaalt", "aaref"] {
                if self.column(name).is_none() {
                    return Err(PluginError::Init(format!(
                        "{NAME}: could not find the required column {name} for pep_match option in the dbNSFP file"
                    )));
                }
            }
            self.aaref_col = self.column("aaref");
            self.aaalt_col = self.column("aaalt");
        }
        Ok(())
    }

    /// The columns to write from the positional parameters after the files
    /// (dbNSFP.pm:231-249): `ALL` takes every header column; a name the header
    /// lacks is kept and warned about; no name at all is an error.
    fn select_columns(
        &mut self,
        names: &[String],
        transcript_specific: Option<&HashSet<String>>,
    ) -> Result<(), PluginError> {
        let mut selected: Vec<String> = Vec::new();
        let mut invalid: Vec<String> = Vec::new();
        for name in names {
            if name == "ALL" {
                selected = self.header.clone();
                break;
            }
            if !selected.contains(name) {
                selected.push(name.clone());
            }
            if self.column(name).is_none() {
                invalid.push(name.clone());
            }
        }
        if selected.is_empty() {
            return Err(PluginError::Init(format!(
                "{NAME}: no input columns were found in file header. Available columns are: {}",
                self.header.join(",")
            )));
        }
        if !invalid.is_empty() {
            warn_on_stderr(&format!(
                "the following columns were not found in file header: {}",
                invalid.join(",")
            ));
        }
        self.selected = selected
            .into_iter()
            .map(|name| SelectedColumn {
                output_name: self.output_name(&name),
                index: self.column(&name),
                transcript_specific: transcript_specific.is_some_and(|set| set.contains(&name)),
                name,
            })
            .collect();
        Ok(())
    }
}

/// A Perl truth value: false for `0` and the empty string, true otherwise.
fn perl_true(value: &str) -> bool {
    !(value.is_empty() || value == "0")
}

/// The dbNSFP version a data file's name carries, as dbNSFP.pm:197-211 sets it,
/// and whether it came from the fallback form without an `a` or `c` suffix,
/// which the constructor accepts with a warning (dbNSFP.pm:205-208).
#[derive(Debug, Clone, PartialEq, Eq)]
struct FileNameVersion {
    version: String,
    suffix_missing: bool,
}

/// Reads the version out of `file`, the data file's path as the `--plugin` value
/// gives it, with the constructor's four tests in its order (dbNSFP.pm:197-211):
/// `2.9` anywhere, `4.0b1` anywhere (version `4.0.1`), the leftmost
/// `[3-9]\.\d+(?:\.\d+)?` run followed by `a` or `c`, then the leftmost such run
/// without the suffix; `None` for a name meeting none, which the constructor
/// refuses.
fn version_from_file_name(file: &str) -> Option<FileNameVersion> {
    let (version, suffix_missing) = if file.contains("2.9") {
        ("2.9", false)
    } else if file.contains("4.0b1") {
        ("4.0.1", false)
    } else if let Some(run) = leftmost_version_run(file, true) {
        (run, false)
    } else {
        (leftmost_version_run(file, false)?, true)
    };
    Some(FileNameVersion {
        version: version.to_string(),
        suffix_missing,
    })
}

/// The leftmost `[3-9]\.\d+(?:\.\d+)?` run of `file`, with the `[ac]` that must
/// follow it when `suffixed`, as Perl's capture reads it. The quantifiers are
/// greedy and every shorter run ends on a digit or a `.`, so at a given start the
/// suffix can only follow the greedy run.
fn leftmost_version_run(file: &str, suffixed: bool) -> Option<&str> {
    let bytes = file.as_bytes();
    for (start, &byte) in bytes.iter().enumerate() {
        if !(b'3'..=b'9').contains(&byte) {
            continue;
        }
        let Some(end) = version_run_end(bytes, start) else {
            continue;
        };
        if !suffixed {
            return Some(&file[start..end]);
        }
        if matches!(bytes.get(end), Some(b'a' | b'c')) {
            return Some(&file[start..=end]);
        }
    }
    None
}

/// The end of the greedy `\.\d+(?:\.\d+)?` tail after the digit at `start`, or
/// `None` when no such tail follows it.
fn version_run_end(bytes: &[u8], start: usize) -> Option<usize> {
    let minor_end = dot_digits_end(bytes, start + 1)?;
    Some(dot_digits_end(bytes, minor_end).unwrap_or(minor_end))
}

/// The end of the greedy `\.\d+` group whose `.` is at `dot`, or `None` when no
/// such group begins there.
fn dot_digits_end(bytes: &[u8], dot: usize) -> Option<usize> {
    if bytes.get(dot) != Some(&b'.') {
        return None;
    }
    let digits = bytes[dot + 1..]
        .iter()
        .take_while(|byte| byte.is_ascii_digit())
        .count();
    if digits == 0 {
        return None;
    }
    Some(dot + 1 + digits)
}

/// A `;`-joined value's entries as Perl's `split(';', ...)` gives them: trailing
/// empty entries dropped.
fn split_entries(value: &str) -> Vec<&str> {
    let mut entries: Vec<&str> = value.split(';').collect();
    while entries.last() == Some(&"") {
        entries.pop();
    }
    entries
}

fn contains_ignore_case(text: &str, needle: &str) -> bool {
    text.to_ascii_lowercase()
        .contains(&needle.to_ascii_lowercase())
}

/// A byte the Ensembl VEP plugin's `\s` matches: it reads the README as bytes
/// with no Unicode layer (dbNSFP.pm:502), so `\s`, `\w` and `\d` are the ASCII
/// classes, with the vertical tab, and a multi-byte blank is text.
fn is_perl_blank(byte: u8) -> bool {
    byte == 0x0B || byte.is_ascii_whitespace()
}

/// `is_perl_blank` over a character of text: a no-break space or a next-line
/// character is text to the plugin's byte-wise `\s`.
fn is_perl_blank_char(c: char) -> bool {
    c.is_ascii() && is_perl_blank(c as u8)
}

/// `bytes` behind its leading blanks (`s/^\s+//`).
fn after_blanks(bytes: &[u8]) -> &[u8] {
    let blanks = bytes
        .iter()
        .take_while(|&&byte| is_perl_blank(byte))
        .count();
    &bytes[blanks..]
}

/// Whether `name` is a dbNSFP README's, as the Ensembl VEP plugin matches the
/// entry's bytes (`/dbnsfp.*readme\.txt/i`, dbNSFP.pm:499): `dbnsfp` and then
/// `readme.txt` on one line of the name, in any letter case, whatever follows;
/// the regex's `.` stops at a newline, so the two never straddle one.
fn is_readme_name(name: &OsStr) -> bool {
    let lower = name.as_encoded_bytes().to_ascii_lowercase();
    lower.split(|&byte| byte == b'\n').any(|line| {
        line.split_once_str(b"dbnsfp")
            .is_some_and(|(_, after_dbnsfp)| after_dbnsfp.contains_str(b"readme.txt"))
    })
}

/// The bytes of the dbNSFP README beside `data_path`: the first, in name order,
/// of the directory entries [`is_readme_name`] admits, where the Ensembl VEP
/// plugin reads the first its directory listing returns (dbNSFP.pm:499);
/// `None` without one or for a remote file.
fn readme_bytes(data_path: &Path) -> Option<Vec<u8>> {
    let readme = std::fs::read_dir(data_path.parent()?)
        .ok()?
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| path.file_name().is_some_and(is_readme_name))
        .min()?;
    std::fs::read(readme).ok()
}

/// The column descriptions of the dbNSFP README beside `data_path`
/// (dbNSFP.pm:485-543, [`readme_bytes`]), read as bytes: its numbered `N  column:
/// description` lines with their continuation lines, ASCII whitespace collapsed,
/// each text decoded from UTF-8 with the replacement character for a byte
/// sequence that is not. A numbered line's name ends at the first colon that a
/// blank and one more byte follow (`(.+?)\:\s+(.+)`, dbNSFP.pm:520), or, when
/// none does, at its first colon, which opens an entry when text follows it, a
/// line the Ensembl VEP plugin reads as no entry; a `0` name or text is kept
/// where the plugin's truth test drops it (dbNSFP.pm:523). The continuation
/// lines of a numbered line that opens no entry are dropped (dbNSFP.pm:521, 527).
fn readme_descriptions(data_path: &Path) -> HashMap<String, String> {
    let Some(bytes) = readme_bytes(data_path) else {
        return HashMap::new();
    };
    let mut descriptions: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    let mut column: Option<Vec<u8>> = None;
    let mut reading = false;
    for line in bytes.split(|&byte| byte == b'\n') {
        // `s/\r$//` (dbNSFP.pm:515) removes one carriage return.
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        let digits = line.iter().take_while(|byte| byte.is_ascii_digit()).count();
        let numbered = digits > 0 && line.get(digits).is_some_and(|&byte| is_perl_blank(byte));
        if numbered {
            reading = true;
            column = None;
            let body = after_blanks(&line[digits..]);
            let colons = || (0..body.len()).filter(|&index| body[index] == b':');
            let plugin_colon = colons().find(|&colon| {
                colon >= 1
                    && body.get(colon + 1).is_some_and(|&byte| is_perl_blank(byte))
                    && colon + 2 < body.len()
            });
            if let Some(colon) = plugin_colon.or_else(|| colons().next()) {
                let name = &body[..colon];
                let text = after_blanks(&body[colon + 1..]);
                if !name.is_empty() && (plugin_colon.is_some() || !text.is_empty()) {
                    column = Some(name.to_vec());
                    descriptions.insert(name.to_vec(), text.to_vec());
                }
            }
        } else if reading
            && line
                .iter()
                .any(|byte| byte.is_ascii_alphanumeric() || *byte == b'_')
        {
            if let Some(description) = column.as_ref().and_then(|name| descriptions.get_mut(name)) {
                description.push(b' ');
                description.extend_from_slice(after_blanks(line));
            }
        } else {
            reading = false;
        }
    }
    descriptions
        .into_iter()
        .map(|(name, description)| {
            (
                String::from_utf8_lossy(&name).into_owned(),
                String::from_utf8_lossy(&collapse_blanks(&description)).into_owned(),
            )
        })
        .collect()
}

/// `text` with each run of blanks one space and none at either end
/// (`s/\s+/ /g`, dbNSFP.pm:537, less the ends), in one pass: splitting into
/// words and joining them allocates a slice per word.
fn collapse_blanks(text: &[u8]) -> Vec<u8> {
    let mut collapsed = Vec::with_capacity(text.len());
    let mut blank_pending = false;
    for &byte in text {
        if is_perl_blank(byte) {
            blank_pending = !collapsed.is_empty();
        } else {
            if blank_pending {
                collapsed.push(b' ');
                blank_pending = false;
            }
            collapsed.push(byte);
        }
    }
    collapsed
}

impl Default for DbNsfpPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for DbNsfpPlugin {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DbNsfpPlugin")
            .field("selected", &self.selected)
            .field("filter", &self.filter)
            .field("output_prefix", &self.output_prefix)
            .field("pep_match", &self.pep_match)
            .field("transcript_match", &self.transcript_match)
            .finish_non_exhaustive()
    }
}

impl crate::traits::sealed::Sealed for DbNsfpPlugin {}

impl BuiltinPlugin for DbNsfpPlugin {
    fn name(&self) -> &str {
        NAME
    }

    fn header_info(&self) -> Vec<(String, String)> {
        self.selected
            .iter()
            .map(|column| {
                let description = match self.descriptions.get(&column.name) {
                    None => format!("{} from dbNSFP file", column.name),
                    Some(text) if text.is_empty() => format!("(from {})", self.basename),
                    Some(text) => format!("(from {}) {text}", self.basename),
                };
                (column.output_name.clone(), description)
            })
            .collect()
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse_with(
            NAME,
            Grammar::PositionalAndKeyValue { flags: &[] },
            NAMED_PARAMS,
            params,
        )?;
        let mut positional = params.positional().iter();
        let Some(file) = positional.next() else {
            return Err(PluginError::Init(format!(
                "{NAME}: no file given; the first parameter is the tabix-indexed dbNSFP file"
            )));
        };
        let path = PathBuf::from(file);
        self.basename = path
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| name.strip_suffix(".gz").unwrap_or(name).to_string())
            .unwrap_or_default();
        // Before the file is opened, as the constructor dies before `add_file`
        // (dbNSFP.pm:210, 215).
        let version = version_from_file_name(file).ok_or_else(|| {
            PluginError::Init(format!(
                "{NAME}: Could not retrieve dbNSFP version from filename {file}"
            ))
        })?;
        if version.suffix_missing {
            warn_on_stderr(&format!(
                "Specific suffix 'a' or 'c' not found for version {}. Proceeding with general version.",
                version.version
            ));
        }

        // `add_file` runs before the header is read (dbNSFP.pm:215).
        check_table_or_store(&path)?;
        self.header = parse_header(&path, 1)?.ok_or_else(|| {
            PluginError::Init(format!("{NAME}: could not read headers from {file}"))
        })?;
        for required in ["alt", "Ensembl_transcriptid"] {
            if self.column(required).is_none() {
                return Err(PluginError::Init(format!(
                    "{NAME}: could not find required column {required} in {file}"
                )));
            }
        }
        let (chr_col, start_col) = self.resolve_coordinates(&path)?;
        let config = TabixConfig {
            file_path: path.clone(),
            chr_col,
            start_col,
            end_col: None,
            ref_col: self.column("ref"),
            alt_col: self.column("alt"),
            zero_based: false,
        };
        // A row's last cell is written as `split /\t/` reads it off the line
        // with its `\r` alone removed (dbNSFP.pm:426-428), its trailing
        // whitespace included; `. ` is not the `.` that is skipped (:408).
        let mut store = open_best_store_with(
            &path,
            config,
            self.assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?;
        // The module names every record on `MT` or `chrMT` `M` to its table
        // (dbNSFP.pm:313), so on GRCh37 `M` is the first name tried; a contig of
        // the other pair of names than the record's is read where the table's
        // rows state the reference bases (`MitochondrialFrame`).
        if let Some(frame) = store.mitochondrial_frame_mut() {
            frame.set_first_name("M");
        }
        self.chr_col = chr_col;
        self.pos_col = start_col;
        self.alt_col = self.column("alt").unwrap_or(0);
        self.transcript_id_col = self.column("Ensembl_transcriptid").unwrap_or(0);
        self.descriptions = readme_descriptions(&path);
        let transcript_specific = self.detect_transcript_specific(store.as_ref())?;

        let mut columns: Vec<String> = positional.cloned().collect();
        match columns.first() {
            Some(logic) if Path::new(logic).exists() => {
                let logic = PathBuf::from(columns.remove(0));
                self.read_replacement_logic(&logic)?;
            }
            _ => {
                self.default_replacements = vec![
                    (";".to_string(), ",".to_string()),
                    ("|".to_string(), "&".to_string()),
                ];
            }
        }
        self.read_named_params(&params)?;
        if self.transcript_match {
            if let Err(failure) = &transcript_specific {
                return Err(PluginError::Init(format!(
                    "{NAME}: transcript_match parameter specified but transcript-specific \
                     field detection failed: {failure}"
                )));
            }
        }
        // After that error and before the column warning, as dbNSFP.pm:613-618
        // and :248 order them.
        if self.transcript_match && self.refseq_or_merged {
            warn_on_stderr(
                "dbNSFP transcript_match may fail with --refseq/--merged; consider \
                 transcript_match=0",
            );
        }
        self.select_columns(&columns, transcript_specific.as_ref().ok())?;

        debug!(
            columns = self.selected.len(),
            header_columns = self.header.len(),
            transcript_specific = transcript_specific.as_ref().map_or(0, HashSet::len),
            "{NAME} plugin initialized with {file}"
        );
        self.store = Some(store);
        Ok(())
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_reference_fasta(&mut self, fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>) {
        if let Some(store) = self.store.as_mut() {
            store.set_reference_fasta(fasta);
        }
    }

    fn set_transcript_set(&mut self, refseq: bool, merged: bool) {
        self.refseq_or_merged = refseq || merged;
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        // Only a single-base substitution is annotated, so only those are read.
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| self.annotates_variant(v) && v.start == v.end)
            .map(|v| (table_contig(v, self.assembly), v.start, v.end))
            .collect();
        if regions.is_empty() {
            return Ok(Some(PrefetchData::new(BatchQueryResult::empty())));
        }
        let batch = store.query_batch(&regions)?;
        Ok(Some(PrefetchData::new(batch)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch: BatchQueryResult = data.downcast()?;
        for variant in variants.iter_mut() {
            if !self.annotates_variant(variant) || variant.start != variant.end {
                continue;
            }
            let records = batch.lookup(
                table_contig(variant, self.assembly),
                variant.start,
                variant.end,
            );
            if records.is_empty() {
                continue;
            }
            for tc_idx in 0..variant.transcript_consequences.len() {
                let consequence = &variant.transcript_consequences[tc_idx];
                if !self.annotates_row(consequence) {
                    continue;
                }
                let fields = self.annotate_row(consequence, variant, &records);
                if !fields.is_empty() {
                    // An empty value here is a cell the row has (:407 skips
                    // the undefined ones), Ensembl VEP's empty string: an
                    // empty tab cell and `""` in JSON (Tab.pm:128,
                    // OutputFactory/JSON.pm:327), which the row's list names.
                    // A later plugin's value replaces an earlier one under
                    // the same key (OutputFactory.pm:2239-2240), its mark too.
                    let row = &mut variant.transcript_consequences[tc_idx];
                    row.plugin_empty_strings
                        .retain(|key| !fields.contains_key(key));
                    row.plugin_empty_strings.extend(
                        fields
                            .iter()
                            .filter(|(_, value)| value.is_empty())
                            .map(|(key, _)| key.clone()),
                    );
                    row.plugin_data.extend(fields);
                }
            }
        }
        Ok(())
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        if let Some(data) = self.prefetch(variants)? {
            self.annotate(variants, data)?;
        }
        Ok(())
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let store = self
            .store
            .as_ref()
            .ok_or_else(|| PluginError::Run(format!("{NAME} plugin not initialized")))?;
        if variant.start != variant.end {
            return Ok(IndexMap::new());
        }
        let records = store.query(
            table_contig(variant, self.assembly),
            variant.start,
            variant.end,
        )?;
        let records: Vec<&TabixRecord> = records.iter().collect();
        Ok(self.annotate_row(consequence, variant, &records))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use std::os::unix::ffi::OsStrExt;
    use std::sync::Arc;
    use vep_core::consequence::FeatureType;
    use vep_core::transcript::Transcript;

    /// The fixture's header: the coordinate, allele and amino acid columns the
    /// plugin reads by name (dbNSFP.pm:117-121, 320-328), then the columns the
    /// tests read.
    const HEADER: &str = "#chr\tpos(1-based)\tref\talt\taaref\taaalt\trs_dbSNP\thg19_chr\thg19_pos(1-based)\tEnsembl_transcriptid\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tMutationTaster_score\tInterpro_domain\tGERP++_RS";

    /// A made-up row of `HEADER`: `pos` on GRCh38 and `hg19` on GRCh37, the lists
    /// `;`-joined.
    fn row(
        site: (&str, u64, u64),
        alleles: (&str, &str),
        aa: (&str, &str),
        transcripts: &str,
        sift: &str,
        interpro: &str,
    ) -> String {
        let (chr, pos, hg19) = site;
        let (reference, alt) = alleles;
        format!(
            "{chr}\t{pos}\t{reference}\t{alt}\t{}\t{}\trs1\t{chr}\t{hg19}\t{transcripts}\t{sift}\tD;T\t0.9;0.1\t1;0.5\t{interpro}\t4.5",
            aa.0, aa.1
        )
    }

    /// The made-up row Ensembl VEP's plugin reads to find the `;`-joined columns
    /// (dbNSFP.pm:444-483), planted in `1:1008170-1082927` on both assemblies.
    fn detection_row() -> String {
        row(
            ("1", 1_045_000, 1_045_000),
            ("A", "G"),
            ("A", "G"),
            "ENST1;ENST2",
            "0.1;0.2",
            "Dom (1)|Other (2);Third (1)",
        )
    }

    fn grch38_file(dir: &Path, rows: &[String]) -> String {
        grch38_file_named(dir, "dbNSFP4.9a_grch38_synthetic.gz", rows)
    }

    fn grch38_file_named(dir: &Path, name: &str, rows: &[String]) -> String {
        let mut lines: Vec<&str> = vec![HEADER];
        lines.extend(rows.iter().map(String::as_str));
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
            .to_string_lossy()
            .to_string()
    }

    fn grch37_file(dir: &Path, rows: &[String]) -> String {
        let mut lines: Vec<&str> = vec![HEADER];
        lines.extend(rows.iter().map(String::as_str));
        write_tabix_fixture(
            dir,
            "dbNSFP4.9a_grch37_synthetic.gz",
            &lines,
            IndexSpec::standard().with_span_cols(8, 9, 9),
        )
        .to_string_lossy()
        .to_string()
    }

    /// The default fixture: the detection row, then a chromosome 21 site with two
    /// rows of different amino acid changes for `C`, a row for `T` and a row for
    /// `G` whose list lengths disagree.
    fn fixture(dir: &Path) -> String {
        grch38_file(
            dir,
            &[
                detection_row(),
                row(
                    ("21", 100, 200),
                    ("A", "C"),
                    ("R", "G"),
                    "ENST1;ENST2",
                    "0.01;0.02",
                    "Kinase (1)|P-loop (2);EGF (1)",
                ),
                row(
                    ("21", 100, 200),
                    ("A", "C"),
                    ("K", "E"),
                    "ENST3",
                    "0.5",
                    ".",
                ),
                row(
                    ("21", 100, 200),
                    ("A", "G"),
                    ("R", "X"),
                    "ENST1;ENST2",
                    "0.3;0.4;0.5",
                    ".",
                ),
                row(
                    ("21", 100, 200),
                    ("A", "T"),
                    ("R", "W"),
                    "ENST9",
                    "0.7",
                    ".",
                ),
            ],
        )
    }

    fn load(assembly: Option<Assembly>, params: &[&str]) -> Result<DbNsfpPlugin, PluginError> {
        let mut plugin = DbNsfpPlugin::new();
        plugin.set_assembly(assembly);
        let params: Vec<String> = params.iter().map(|p| p.to_string()).collect();
        plugin.init(&params)?;
        Ok(plugin)
    }

    fn variant(alt: &str) -> InputVariant {
        InputVariant::new(
            "21".into(),
            100,
            100,
            b"A".to_vec(),
            alt.as_bytes().to_vec(),
        )
    }

    fn consequence(
        transcript_id: &str,
        term: Consequence,
        amino_acids: Option<&str>,
    ) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.transcript = Some(Arc::new({
            let mut transcript = Transcript::default();
            transcript.stable_id = transcript_id.into();
            transcript
        }));
        consequence.consequences = [term].into_iter().collect();
        consequence.amino_acids = amino_acids.map(str::to_string);
        consequence
    }

    fn annotated(
        plugin: &DbNsfpPlugin,
        variant: &mut InputVariant,
        rows: Vec<TranscriptConsequence>,
    ) -> Vec<IndexMap<String, String>> {
        variant.transcript_consequences = rows;
        plugin.run_batch(std::slice::from_mut(variant)).unwrap();
        variant
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.clone())
            .collect()
    }

    fn pairs(fields: &IndexMap<String, String>) -> Vec<(&str, &str)> {
        fields
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect()
    }

    /// `add_file` checks for `<file>.tbi` before the header is read
    /// (dbNSFP.pm:215): a table lacking both is refused naming the index.
    #[test]
    fn table_without_its_tbi_is_refused_before_its_header_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "dbNSFP4.9a_grch38_synthetic.gz",
            &["#chr\tpos(1-based)", "21\t100"],
            IndexSpec::standard(),
        );
        std::fs::remove_file(dir.path().join("dbNSFP4.9a_grch38_synthetic.gz.tbi")).unwrap();
        let message = match load(None, &[&path.display().to_string()]) {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected an init error, got {other:?}"),
        };
        assert!(
            message.contains("dbNSFP4.9a_grch38_synthetic.gz.tbi not found"),
            "{message}"
        );
    }

    /// The output names are the column names, and `output_prefix=1` prefixes them
    /// with `dbNSFP_` (dbNSFP.pm:159, 257-261, 582); a name already prefixed stays.
    #[test]
    fn default_keys_have_no_prefix_and_output_prefix_adds_it() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let plain = load(Some(Assembly::Grch38), &[&file, "SIFT_score", "GERP++_RS"]).unwrap();
        assert_eq!(
            plain.header_info(),
            vec![
                (
                    "SIFT_score".to_string(),
                    "SIFT_score from dbNSFP file".to_string()
                ),
                (
                    "GERP++_RS".to_string(),
                    "GERP++_RS from dbNSFP file".to_string()
                ),
            ]
        );
        let prefixed = load(
            Some(Assembly::Grch38),
            &[&file, "output_prefix=1", "SIFT_score", "dbNSFP_POPMAX_AC"],
        )
        .unwrap();
        let names: Vec<String> = prefixed
            .header_info()
            .into_iter()
            .map(|(name, _)| name)
            .collect();
        assert_eq!(names, ["dbNSFP_SIFT_score", "dbNSFP_POPMAX_AC"]);
        let mut variant = variant("C");
        let rows = annotated(
            &prefixed,
            &mut variant,
            vec![consequence(
                "ENST1",
                Consequence::MissenseVariant,
                Some("R/G"),
            )],
        );
        assert_eq!(
            pairs(&rows[0]),
            [
                ("dbNSFP_SIFT_score", "0.01,0.02"),
                ("dbNSFP_POPMAX_AC", "invalid_field")
            ]
        );
    }

    /// Without `consequence=` only missense, stop-lost, stop-gained and start-lost
    /// rows are annotated (dbNSFP.pm:157, 290-292); `consequence=ALL` annotates every
    /// row and a `&`-joined list replaces the set (dbNSFP.pm:591-597).
    #[test]
    fn default_consequence_filter_and_consequence_all() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let rows = || {
            vec![
                consequence("ENST1", Consequence::MissenseVariant, Some("R/G")),
                consequence("ENST2", Consequence::SynonymousVariant, Some("R")),
                consequence("ENST3", Consequence::IntronVariant, None),
                consequence("ENST9", Consequence::StopGained, Some("R/*")),
            ]
        };
        let filtered = load(
            Some(Assembly::Grch38),
            &[&file, "pep_match=0", "SIFT_score"],
        )
        .unwrap();
        let mut variant = variant("C");
        let got = annotated(&filtered, &mut variant, rows());
        assert_eq!(
            got.iter()
                .map(|f| f.get("SIFT_score").map(String::as_str))
                .collect::<Vec<_>>(),
            [Some("0.01,0.02"), None, None, Some("0.01,0.02")]
        );
        let all = load(
            Some(Assembly::Grch38),
            &[&file, "consequence=all", "pep_match=0", "SIFT_score"],
        )
        .unwrap();
        let got = annotated(&all, &mut variant, rows());
        assert!(got
            .iter()
            .all(|f| f.get("SIFT_score").map(String::as_str) == Some("0.01,0.02")));
        let listed = load(
            Some(Assembly::Grch38),
            &[
                &file,
                "consequence=synonymous_variant&intron_variant",
                "pep_match=0",
                "SIFT_score",
            ],
        )
        .unwrap();
        let got = annotated(&listed, &mut variant, rows());
        assert_eq!(
            got.iter()
                .map(|f| f.get("SIFT_score").map(String::as_str))
                .collect::<Vec<_>>(),
            [None, Some("0.01,0.02"), Some("0.01,0.02"), None]
        );
    }

    /// Under the default `pep_match=1` a row takes the first record whose `aaref/aaalt`
    /// equals its amino acid change, `X` read as `*` (dbNSFP.pm:340-344); a row with
    /// another change or none is left alone; `pep_match=0` takes the first record of
    /// the allele (dbNSFP.pm:583).
    #[test]
    fn pep_match_default_requires_matching_peptides() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let matching = load(Some(Assembly::Grch38), &[&file, "SIFT_score"]).unwrap();
        let mut c = variant("C");
        let got = annotated(
            &matching,
            &mut c,
            vec![
                consequence("ENST1", Consequence::MissenseVariant, Some("R/G")),
                consequence("ENST3", Consequence::MissenseVariant, Some("K/E")),
                consequence("ENST2", Consequence::MissenseVariant, Some("R/W")),
                consequence("ENST1", Consequence::MissenseVariant, None),
            ],
        );
        assert_eq!(
            got.iter()
                .map(|f| f.get("SIFT_score").map(String::as_str))
                .collect::<Vec<_>>(),
            [Some("0.01,0.02"), Some("0.5"), None, None]
        );
        let mut g = variant("G");
        let got = annotated(
            &matching,
            &mut g,
            vec![consequence("ENST1", Consequence::StopGained, Some("R/*"))],
        );
        assert_eq!(
            got[0].get("SIFT_score").map(String::as_str),
            Some("0.3,0.4,0.5")
        );
        let unmatched = load(
            Some(Assembly::Grch38),
            &[&file, "pep_match=0", "SIFT_score"],
        )
        .unwrap();
        let got = annotated(
            &unmatched,
            &mut c,
            vec![consequence(
                "ENST2",
                Consequence::MissenseVariant,
                Some("R/W"),
            )],
        );
        assert_eq!(
            got[0].get("SIFT_score").map(String::as_str),
            Some("0.01,0.02")
        );
    }

    /// Under `transcript_match=1` a transcript-specific column keeps the entry of the
    /// row's transcript; a transcript the record does not list loses every
    /// transcript-specific column and keeps the others; a list whose length differs
    /// from the transcript list, or too short for the index, keeps its whole value
    /// (dbNSFP.pm:350-390).
    #[test]
    fn transcript_match_deletes_transcript_fields_when_unmatched() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let plugin = load(
            Some(Assembly::Grch38),
            &[
                &file,
                "transcript_match=1",
                "pep_match=0",
                "SIFT_score",
                "MutationTaster_score",
                "GERP++_RS",
            ],
        )
        .unwrap();
        let mut c = variant("C");
        let got = annotated(
            &plugin,
            &mut c,
            vec![
                consequence("ENST2", Consequence::MissenseVariant, Some("R/G")),
                consequence("ENST3", Consequence::MissenseVariant, Some("K/E")),
            ],
        );
        assert_eq!(
            pairs(&got[0]),
            [
                ("SIFT_score", "0.02"),
                ("MutationTaster_score", "1,0.5"),
                ("GERP++_RS", "4.5")
            ]
        );
        assert_eq!(
            pairs(&got[1]),
            [("MutationTaster_score", "1,0.5"), ("GERP++_RS", "4.5")]
        );
        let mut g = variant("G");
        let got = annotated(
            &plugin,
            &mut g,
            vec![consequence("ENST1", Consequence::StopGained, Some("R/*"))],
        );
        assert_eq!(
            got[0].get("SIFT_score").map(String::as_str),
            Some("0.3,0.4,0.5")
        );
    }

    /// `;` becomes `,` and `|` becomes `&` in every value by default (dbNSFP.pm:550-551,
    /// 410-417).
    #[test]
    fn semicolon_and_pipe_are_replaced_by_default() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let plugin = load(
            Some(Assembly::Grch38),
            &[&file, "SIFT_score", "Interpro_domain"],
        )
        .unwrap();
        let mut c = variant("C");
        let got = annotated(
            &plugin,
            &mut c,
            vec![consequence(
                "ENST1",
                Consequence::MissenseVariant,
                Some("R/G"),
            )],
        );
        assert_eq!(
            pairs(&got[0]),
            [
                ("SIFT_score", "0.01,0.02"),
                ("Interpro_domain", "Kinase (1)&P-loop (2),EGF (1)")
            ]
        );
    }

    /// A replacement-logic file given as the second positional parameter replaces the
    /// default pairs and adds per-column pairs (dbNSFP.pm:220-226, 552-566); a file
    /// without `default` rows or with a short line is refused.
    #[test]
    fn replacement_logic_file_overrides_per_column() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        // The file Ensembl VEP ships beside the plugin: the default pairs and
        // `?` and `+` for Interpro_domain.
        let shipped = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/fixtures/plugins/dbNSFP_replacement_logic"
        );
        let plugin = load(
            Some(Assembly::Grch38),
            &[&file, shipped, "SIFT_score", "Interpro_domain"],
        )
        .unwrap();
        let mut c = variant("C");
        let got = annotated(
            &plugin,
            &mut c,
            vec![consequence(
                "ENST1",
                Consequence::MissenseVariant,
                Some("R/G"),
            )],
        );
        assert_eq!(
            pairs(&got[0]),
            [
                ("SIFT_score", "0.01,0.02"),
                ("Interpro_domain", "Kinase (1)+P-loop (2)?EGF (1)")
            ]
        );

        let logic = dir.path().join("dbNSFP_replacement_logic");
        let logic_path = logic.to_string_lossy().to_string();
        std::fs::write(&logic, "colname from to\nInterpro_domain ; ?\n").unwrap();
        let err = load(Some(Assembly::Grch38), &[&file, &logic_path, "SIFT_score"]).unwrap_err();
        assert!(
            err.to_string().contains("no default replacement logic"),
            "{err}"
        );
        std::fs::write(&logic, "colname from to\ndefault ;\n").unwrap();
        let err = load(Some(Assembly::Grch38), &[&file, &logic_path, "SIFT_score"]).unwrap_err();
        assert!(err.to_string().contains("requires 3 values"), "{err}");
    }

    /// A logic line that opens with a blank is refused at load as Ensembl VEP's
    /// plugin refuses it: its `split /\s+/` opens an empty first field, so the
    /// three-value test fails (dbNSFP.pm:557-559). The blanks are the ASCII ones,
    /// so a line whose values a no-break space or a next-line character separates
    /// is one value; a blank line, a `0` value and a short line fail the same
    /// test. The refusal names the file and the line.
    #[test]
    fn replacement_logic_line_opening_with_a_blank_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let logic = dir.path().join("dbNSFP_replacement_logic");
        let logic_path = logic.to_string_lossy().to_string();
        let refused: [(&[u8], usize); 10] = [
            (b"colname from to\n  default ; ,\ndefault | &\n", 2),
            (b"colname from to\n\tdefault ; ,\ndefault | &\n", 2),
            (b" colname from to\ndefault ; ,\ndefault | &\n", 1),
            (
                b"colname from to\ndefault ; ,\ndefault | &\n  Interpro_domain ; ?\n",
                4,
            ),
            (
                "colname from to\ndefault\u{a0};\u{a0},\ndefault | &\n".as_bytes(),
                2,
            ),
            (
                "colname from to\ndefault\u{85};\u{85},\ndefault | &\n".as_bytes(),
                2,
            ),
            (b"colname from to\ndefault ; ,\n\ndefault | &\n", 3),
            (b"colname from to\ndefault ; ,\n   \ndefault | &\n", 3),
            (b"colname from to\ndefault 0 ,\ndefault | &\n", 2),
            (b"colname from to\ndefault ;\ndefault | &\n", 2),
        ];
        for (body, line) in refused {
            std::fs::write(&logic, body).unwrap();
            let err =
                load(Some(Assembly::Grch38), &[&file, &logic_path, "SIFT_score"]).unwrap_err();
            let expected =
                format!("replacement logic file {logic_path} line {line} requires 3 values");
            assert!(
                err.to_string().contains(&expected),
                "{}: {err}",
                String::from_utf8_lossy(body)
            );
        }
    }

    /// The shapes Ensembl VEP's plugin reads are read (dbNSFP.pm:553-561): a
    /// trailing blank, CRLF line endings, a fourth value (dropped), a row opening
    /// with `colname` (skipped), a vertical tab or form feed between values, a
    /// `0.0` value, and a non-UTF-8 value, which stands as the replacement
    /// character where Ensembl VEP's plugin prints the byte; a column's own `;` rule
    /// leaves its `|` alone.
    #[test]
    fn replacement_logic_file_shapes_ensembl_vep_reads_are_read() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let logic = dir.path().join("dbNSFP_replacement_logic");
        let logic_path = logic.to_string_lossy().to_string();
        let accepted: [(&[u8], &str, &str); 8] = [
            (
                b"colname from to\ndefault ; , \ndefault | &\n",
                "0.01,0.02",
                "Kinase (1)&P-loop (2),EGF (1)",
            ),
            (
                b"colname from to\r\ndefault ; ,\r\ndefault | &\r\n",
                "0.01,0.02",
                "Kinase (1)&P-loop (2),EGF (1)",
            ),
            (
                b"colname from to\ndefault ; ? extra\ndefault | &\n",
                "0.01?0.02",
                "Kinase (1)&P-loop (2)?EGF (1)",
            ),
            (
                b"colname from to\ncolnameX a b\ndefault ; ,\ndefault | &\n",
                "0.01,0.02",
                "Kinase (1)&P-loop (2),EGF (1)",
            ),
            (
                b"colname from to\ndefault\x0b;\x0b,\ndefault\x0c|\x0c&\n",
                "0.01,0.02",
                "Kinase (1)&P-loop (2),EGF (1)",
            ),
            (
                b"colname from to\ndefault ; 0.0\ndefault | &\n",
                "0.010.00.02",
                "Kinase (1)&P-loop (2)0.0EGF (1)",
            ),
            (
                b"colname from to\ndefault ; \xff\ndefault | &\n",
                "0.01\u{fffd}0.02",
                "Kinase (1)&P-loop (2)\u{fffd}EGF (1)",
            ),
            (
                b"colname from to\ndefault ; ,\ndefault | &\nInterpro_domain ; ?\n",
                "0.01,0.02",
                "Kinase (1)|P-loop (2)?EGF (1)",
            ),
        ];
        for (body, sift, interpro) in accepted {
            std::fs::write(&logic, body).unwrap();
            let plugin = load(
                Some(Assembly::Grch38),
                &[&file, &logic_path, "SIFT_score", "Interpro_domain"],
            )
            .unwrap_or_else(|err| panic!("{}: {err}", String::from_utf8_lossy(body)));
            let mut c = variant("C");
            let got = annotated(
                &plugin,
                &mut c,
                vec![consequence(
                    "ENST1",
                    Consequence::MissenseVariant,
                    Some("R/G"),
                )],
            );
            assert_eq!(
                pairs(&got[0]),
                [("SIFT_score", sift), ("Interpro_domain", interpro)],
                "{}",
                String::from_utf8_lossy(body)
            );
        }
    }

    /// `ALL` selects every column of the header in header order, the coordinate
    /// columns included (dbNSFP.pm:233-236); a `.` value is not written (dbNSFP.pm:408).
    #[test]
    fn all_selects_every_column() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let plugin = load(Some(Assembly::Grch38), &[&file, "ALL"]).unwrap();
        let names: Vec<String> = plugin
            .header_info()
            .into_iter()
            .map(|(name, _)| name)
            .collect();
        let header: Vec<&str> = HEADER.trim_start_matches('#').split('\t').collect();
        assert_eq!(names, header);
        let mut c = variant("C");
        let got = annotated(
            &plugin,
            &mut c,
            vec![consequence(
                "ENST3",
                Consequence::MissenseVariant,
                Some("K/E"),
            )],
        );
        assert_eq!(got[0].get("pos(1-based)").map(String::as_str), Some("100"));
        assert_eq!(
            got[0].get("Ensembl_transcriptid").map(String::as_str),
            Some("ENST3")
        );
        assert!(!got[0].contains_key("Interpro_domain"));
        assert_eq!(got[0].len(), header.len() - 1);
    }

    /// A GRCh37 run reads `hg19_chr` and `hg19_pos(1-based)`, the columns the GRCh37
    /// preparation indexes (dbNSFP.pm:73-74, 316-325); the same file indexed for the
    /// other assembly is refused, naming both columns.
    #[test]
    fn grch37_file_resolves_hg19_pos_column() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            detection_row(),
            row(
                ("21", 100, 200),
                ("A", "C"),
                ("R", "G"),
                "ENST1;ENST2",
                "0.01;0.02",
                ".",
            ),
        ];
        let file37 = grch37_file(dir.path(), &rows);
        let plugin37 = load(Some(Assembly::Grch37), &[&file37, "SIFT_score"]).unwrap();
        let mut c = InputVariant::new("21".into(), 200, 200, b"A".to_vec(), b"C".to_vec());
        let got = annotated(
            &plugin37,
            &mut c,
            vec![consequence(
                "ENST1",
                Consequence::MissenseVariant,
                Some("R/G"),
            )],
        );
        assert_eq!(
            got[0].get("SIFT_score").map(String::as_str),
            Some("0.01,0.02")
        );

        let err = load(Some(Assembly::Grch38), &[&file37, "SIFT_score"]).unwrap_err();
        let message = err.to_string();
        assert!(
            message.contains("'pos(1-based)' (column 2)")
                && message.contains("'hg19_pos(1-based)' (column 9)"),
            "{message}"
        );
        let file38 = grch38_file(dir.path(), &rows);
        let err = load(Some(Assembly::Grch37), &[&file38, "SIFT_score"]).unwrap_err();
        assert!(
            err.to_string()
                .contains("indexed on 'pos(1-based)' (column 2)"),
            "{err}"
        );
    }

    /// A named parameter the Ensembl VEP plugin does not know is an error
    /// (dbNSFP.pm:601).
    #[test]
    fn unknown_named_param_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let err = load(Some(Assembly::Grch38), &[&file, "pepmatch=0", "SIFT_score"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "plugin initialization failed: dbNSFP: unknown parameter 'pepmatch'; dbNSFP accepts consequence, output_prefix, pep_match, transcript_match"
        );
    }

    /// The version is read from the data file's name as dbNSFP.pm:197-211 reads
    /// it from the whole path: the corpus name carries `4.9a`; a name carrying
    /// none is refused with the constructor's text before the file is opened; a
    /// bare `4.9` is accepted with the suffix warning and annotates the same rows
    /// as the suffixed name. `2.9` and `4.0b1` are read anywhere in the path.
    #[test]
    fn version_is_read_from_the_file_name_as_ensembl_vep_reads_it() {
        let exact = |version: &str| {
            Some(FileNameVersion {
                version: version.to_string(),
                suffix_missing: false,
            })
        };
        let bare = |version: &str| {
            Some(FileNameVersion {
                version: version.to_string(),
                suffix_missing: true,
            })
        };
        assert_eq!(
            version_from_file_name("dbNSFP4.9a_grch38_synthetic.gz"),
            exact("4.9a")
        );
        assert_eq!(version_from_file_name("dbnsfp_synthetic.gz"), None);
        assert_eq!(
            version_from_file_name("dbNSFP4.9_grch38_synthetic.gz"),
            bare("4.9")
        );
        assert_eq!(version_from_file_name("/path/to/dbNSFP_synthetic.gz"), None);
        assert_eq!(
            version_from_file_name("/data/grch38/dbNSFP_v3.synthetic.gz"),
            None
        );
        assert_eq!(
            version_from_file_name("/data/dbNSFP2.9.1_synthetic.gz"),
            exact("2.9")
        );
        assert_eq!(
            version_from_file_name("/data/dbNSFP4.0b1_grch38_synthetic.gz"),
            exact("4.0.1")
        );
        assert_eq!(
            version_from_file_name("/data/dbNSFP4.9.1c_grch37_synthetic.gz"),
            exact("4.9.1c")
        );
        assert_eq!(
            version_from_file_name("/data/dbNSFP4.9.1_grch37_synthetic.gz"),
            bare("4.9.1")
        );
        assert_eq!(
            version_from_file_name("/releases/5.1a/dbNSFP_synthetic.gz"),
            exact("5.1a")
        );
        assert_eq!(
            version_from_file_name("dbNSFP43.9a_synthetic.gz"),
            exact("3.9a")
        );

        let dir = tempfile::tempdir().unwrap();
        let rows = [
            detection_row(),
            row(
                ("21", 100, 200),
                ("A", "C"),
                ("R", "G"),
                "ENST1;ENST2",
                "0.01;0.02",
                ".",
            ),
        ];
        let unversioned = grch38_file_named(dir.path(), "dbnsfp_synthetic.gz", &rows);
        let err = load(Some(Assembly::Grch38), &[&unversioned, "SIFT_score"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            format!(
                "plugin initialization failed: dbNSFP: Could not retrieve dbNSFP version from filename {unversioned}"
            )
        );

        let sift_of = |file: &str| {
            let plugin = load(Some(Assembly::Grch38), &[file, "SIFT_score"]).unwrap();
            let mut c = variant("C");
            let got = annotated(
                &plugin,
                &mut c,
                vec![consequence(
                    "ENST1",
                    Consequence::MissenseVariant,
                    Some("R/G"),
                )],
            );
            got[0].get("SIFT_score").cloned()
        };
        let suffixed = grch38_file(dir.path(), &rows);
        let bare_version = grch38_file_named(dir.path(), "dbNSFP4.9_grch38_synthetic.gz", &rows);
        assert_eq!(sift_of(&suffixed), Some("0.01,0.02".to_string()));
        assert_eq!(sift_of(&bare_version), sift_of(&suffixed));
    }

    /// The transcript-specific columns are read from the first `;`-bearing record in
    /// `1:1008170-1082927` (dbNSFP.pm:444-483); without one `transcript_match=1` is
    /// refused (dbNSFP.pm:613-615), and the columns of the fixed list stay whole.
    #[test]
    fn transcript_specific_columns_come_from_the_detection_region() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let plugin = load(
            Some(Assembly::Grch38),
            &[&file, "transcript_match=1", "SIFT_score"],
        )
        .unwrap();
        let specific: Vec<&str> = plugin
            .selected
            .iter()
            .filter(|c| c.transcript_specific)
            .map(|c| c.name.as_str())
            .collect();
        assert_eq!(specific, ["SIFT_score"]);
        let without = grch38_file(
            dir.path(),
            &[row(
                ("21", 100, 200),
                ("A", "C"),
                ("R", "G"),
                "ENST1;ENST2",
                "0.01;0.02",
                ".",
            )],
        );
        let err = load(
            Some(Assembly::Grch38),
            &[&without, "transcript_match=1", "SIFT_score"],
        )
        .unwrap_err()
        .to_string();
        assert!(
            err.contains("transcript-specific field detection failed")
                && err.contains("no record of 1:1008170-1082927 holds a `;`-separated value"),
            "{err}"
        );
        load(Some(Assembly::Grch38), &[&without, "SIFT_score"]).unwrap();
    }

    /// The detection region is read under the bare name `1` alone, as the
    /// Ensembl VEP plugin's `tabix` query reads it (dbNSFP.pm:162, 449): a file
    /// naming its contigs `chr1` and `chr21` yields no `;`-bearing record there,
    /// so `transcript_match=1` is refused (dbNSFP.pm:613-615) naming the file's
    /// spelling, while a run without it annotates its rows, looked up under
    /// either spelling as the base class looks them up
    /// (BaseVepTabixPlugin.pm:679-715).
    #[test]
    fn chr_prefixed_file_is_refused_under_transcript_match_as_the_plugin_refuses_it() {
        let dir = tempfile::tempdir().unwrap();
        let file = grch38_file(dir.path(), &chr_prefixed_rows("0.01;0.02"));
        assert_chr_prefixed_file_refused_then_annotated(&file, "0.01,0.02");
    }

    /// A `.vpd` store built from the `chr`-named file answers as the file does:
    /// its records carry the store's spelling, so the detection under the bare
    /// name finds none there either and `transcript_match=1` is refused, while a
    /// run without it annotates. The store holds the file's own rows, and the
    /// file's `.tbi` is removed once the store is written, so the store alone
    /// can answer.
    #[test]
    fn a_vpd_store_of_a_chr_prefixed_file_is_refused_under_transcript_match() {
        use crate::binary_store::{BinaryStoreWriter, ContentIdentity};

        let dir = tempfile::tempdir().unwrap();
        let rows = chr_prefixed_rows("0.01;0.02");
        let file = grch38_file(dir.path(), &rows);
        let header: Arc<Vec<String>> = Arc::new(
            HEADER
                .trim_start_matches('#')
                .split('\t')
                .map(str::to_string)
                .collect(),
        );
        let record = |line: &String| {
            let columns: Vec<String> = line.split('\t').map(str::to_string).collect();
            let start: u64 = columns[1].parse().unwrap();
            TabixRecord {
                chr: columns[0].clone(),
                start,
                end: start,
                columns,
                header: Some(Arc::clone(&header)),
            }
        };
        BinaryStoreWriter::new(Path::new(&file).with_extension("vpd"))
            .write(
                &header,
                &[
                    ("chr1".to_string(), vec![record(&rows[0])]),
                    ("chr21".to_string(), vec![record(&rows[1])]),
                ],
                &TabixConfig::default(),
                &ContentIdentity::of_file(Path::new(&file)).unwrap(),
            )
            .unwrap();
        std::fs::remove_file(format!("{file}.tbi")).unwrap();
        assert_chr_prefixed_file_refused_then_annotated(&file, "0.01,0.02");
    }

    /// The detection row and one chromosome 21 record with `sift` as its
    /// `SIFT_score`, both under `chr` names.
    fn chr_prefixed_rows(sift: &str) -> [String; 2] {
        [
            row(
                ("chr1", 1_045_000, 1_045_000),
                ("A", "G"),
                ("A", "G"),
                "ENST1;ENST2",
                "0.1;0.2",
                ".",
            ),
            row(
                ("chr21", 100, 200),
                ("A", "C"),
                ("R", "G"),
                "ENST1;ENST2",
                sift,
                ".",
            ),
        ]
    }

    /// `transcript_match=1` on `file` is refused naming the `chr1` spelling and
    /// the remedy; a load without it annotates the `21:100 A>C` record on
    /// `ENST1` with `sift`, its `chr` column the file's own spelling.
    fn assert_chr_prefixed_file_refused_then_annotated(file: &str, sift: &str) {
        let err = load(
            Some(Assembly::Grch38),
            &[file, "transcript_match=1", "SIFT_score"],
        )
        .unwrap_err()
        .to_string();
        assert!(
            err.contains("transcript-specific field detection failed")
                && err.contains("names chromosome 1 `chr1`")
                && err.contains("run without transcript_match"),
            "{err}"
        );
        let plugin = load(Some(Assembly::Grch38), &[file, "chr", "SIFT_score"]).unwrap();
        let mut c = variant("C");
        let got = annotated(
            &plugin,
            &mut c,
            vec![consequence(
                "ENST1",
                Consequence::MissenseVariant,
                Some("R/G"),
            )],
        );
        assert_eq!(pairs(&got[0]), [("chr", "chr21"), ("SIFT_score", sift)]);
    }

    /// The README is any directory entry whose name carries `dbnsfp` and then
    /// `readme.txt`, in any letter case, whatever follows (dbNSFP.pm:499,
    /// `/dbnsfp.*readme\.txt/i`): the made-up `dbNSFP4.9a.readme.txt.bak` is read, its text
    /// describes the column and, under `transcript_match=1`, its "separated by"
    /// marks a column the detection row holds whole, dropped on a transcript the
    /// record does not list (dbNSFP.pm:364-368). Of several such entries the
    /// Ensembl VEP plugin reads the first its directory listing returns; vep-rs
    /// reads the first in name order, whatever the order of creation.
    #[test]
    fn readme_name_is_matched_as_the_ensembl_vep_plugin_matches_it() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let stale = "1\tchr: A made-up text for the first column.\n16\tGERP++_RS: GERP++ RS score, values separated by \";\" for each transcript\n";
        std::fs::write(dir.path().join("dbNSFP4.9a.readme.txt.bak"), stale).unwrap();
        let params = [
            file.as_str(),
            "pep_match=0",
            "transcript_match=1",
            "GERP++_RS",
            "SIFT_score",
        ];
        let plugin = load(Some(Assembly::Grch38), &params).unwrap();
        assert_eq!(
            plugin.header_info()[0],
            (
                "GERP++_RS".to_string(),
                "(from dbNSFP4.9a_grch38_synthetic) GERP++ RS score, values separated by \";\" for each transcript".to_string()
            )
        );
        let unlisted = || {
            vec![consequence(
                "ENST1",
                Consequence::MissenseVariant,
                Some("R/W"),
            )]
        };
        let mut t = variant("T");
        assert_eq!(pairs(&annotated(&plugin, &mut t, unlisted())[0]), []);

        // The made-up entry first in name order is read, created after the other.
        std::fs::write(
            dir.path().join("dbNSFP4.9a.readme.txt"),
            "1\tchr: A made-up text for the first column.\n16\tGERP++_RS: GERP++ RS score\n",
        )
        .unwrap();
        let plugin = load(Some(Assembly::Grch38), &params).unwrap();
        assert_eq!(
            plugin.header_info()[0],
            (
                "GERP++_RS".to_string(),
                "(from dbNSFP4.9a_grch38_synthetic) GERP++ RS score".to_string()
            )
        );
        let mut t = variant("T");
        assert_eq!(
            pairs(&annotated(&plugin, &mut t, unlisted())[0]),
            [("GERP++_RS", "4.5")]
        );

        // The plugin matches the entry's bytes: a name that is not UTF-8 is
        // read, and the regex's `.` stops at a newline. A Linux file system
        // stores such a name's bytes as they are; other platforms refuse it.
        let mut names = vec![
            (b"DBNSFP4.9A.README.TXT.old".as_slice(), true),
            (b"my-dbnsfp-readme.txt", true),
            (b"readme.txt.dbNSFP", false),
            (b"dbNSFP_readme.text", false),
            (b"dbNSFP\n.readme.txt", false),
            (b"dbNSFP\ndbNSFP.readme.txt\n", true),
        ];
        if cfg!(target_os = "linux") {
            names.push((b"dbNSFP\xff.readme.txt", true));
        }
        for (name, read) in names {
            let dir = tempfile::tempdir().unwrap();
            let file = fixture(dir.path());
            std::fs::write(dir.path().join(OsStr::from_bytes(name)), stale).unwrap();
            assert_eq!(
                readme_descriptions(Path::new(&file)).contains_key("GERP++_RS"),
                read,
                "{}",
                name.escape_ascii()
            );
        }
    }

    /// A `dbNSFP*readme.txt` beside the data file supplies the header descriptions,
    /// prefixed with the file's base name, and marks the columns it calls
    /// "separated by" as transcript-specific (dbNSFP.pm:268-280, 469-476, 485-543).
    #[test]
    fn readme_beside_the_file_describes_columns() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        std::fs::write(
            dir.path().join("dbNSFP4.9a.readme.txt"),
            "A made-up list of columns:\n1\tchr: A made-up text for the first column.\n11\tSIFT_score: A made-up number for each transcript. The numbers are separated\n\t\tby \";\", in the order of the transcript list.\n\n16\tGERP++_RS: A made-up number for the site,  with a doubled blank in its text.\n",
        )
        .unwrap();
        let plugin = load(
            Some(Assembly::Grch38),
            &[&file, "transcript_match=1", "GERP++_RS", "SIFT_score"],
        )
        .unwrap();
        assert_eq!(
            plugin.header_info(),
            vec![
                ("GERP++_RS".to_string(), "(from dbNSFP4.9a_grch38_synthetic) A made-up number for the site, with a doubled blank in its text.".to_string()),
                ("SIFT_score".to_string(), "(from dbNSFP4.9a_grch38_synthetic) A made-up number for each transcript. The numbers are separated by \";\", in the order of the transcript list.".to_string()),
            ]
        );
    }

    /// A numbered README line opens a column's entry when a non-empty name, a
    /// colon, a blank and at least one more character follow its number
    /// (dbNSFP.pm:520-523): a line without a colon, with nothing or one blank after
    /// the colon, or with an empty name opens none, and its continuation lines are
    /// dropped rather than appended to the column before it (dbNSFP.pm:521, 527).
    /// Two blanks after the colon open an entry whose text is its continuation;
    /// the Ensembl VEP plugin keeps the blank the entry opened with and prints two
    /// spaces after the prefix where the `chr` row here expects the collapsed one.
    #[test]
    fn readme_line_without_a_column_opens_no_entry() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let readme = concat!(
            "A made-up list of columns:\n",
            "1\tchr:  \n",
            "\t\tA made-up text for the first column.\n",
            "14\tMutationTaster_score: MutationTaster p-value\n",
            "16\tGERP++_RS\n",
            "\t\tContinuation without a colon on the numbered line.\n",
            "12\tSIFT_pred:\n",
            "\t\tContinuation after a colon with nothing behind it.\n",
            "11\tSIFT_score: A made-up number for each transcript. The numbers are separated\n",
            "\t\tby \";\", in the order of the transcript list.\n",
            "13\t: Polyphen2 description with an empty name.\n",
            "\t\tContinuation after an empty name.\n",
            "13\tPolyphen2_HDIV_score: \n",
            "\t\tContinuation after a colon and one blank.\n",
        );
        std::fs::write(dir.path().join("dbNSFP4.9a.readme.txt"), readme).unwrap();
        let plugin = load(
            Some(Assembly::Grch38),
            &[
                &file,
                "transcript_match=1",
                "chr",
                "MutationTaster_score",
                "GERP++_RS",
                "SIFT_pred",
                "SIFT_score",
                "Polyphen2_HDIV_score",
            ],
        )
        .unwrap();
        let from = "(from dbNSFP4.9a_grch38_synthetic)";
        assert_eq!(
            plugin.header_info(),
            vec![
                ("chr".to_string(), format!("{from} A made-up text for the first column.")),
                (
                    "MutationTaster_score".to_string(),
                    format!("{from} MutationTaster p-value")
                ),
                ("GERP++_RS".to_string(), "GERP++_RS from dbNSFP file".to_string()),
                ("SIFT_pred".to_string(), "SIFT_pred from dbNSFP file".to_string()),
                (
                    "SIFT_score".to_string(),
                    format!("{from} A made-up number for each transcript. The numbers are separated by \";\", in the order of the transcript list.")
                ),
                (
                    "Polyphen2_HDIV_score".to_string(),
                    "Polyphen2_HDIV_score from dbNSFP file".to_string()
                ),
            ]
        );
    }

    /// The columns the README names and the header carries, in `names` order,
    /// with their header descriptions; the README is written beside `file`.
    fn described(dir: &Path, file: &str, readme: &[u8], names: &[&str]) -> Vec<(String, String)> {
        std::fs::write(dir.join("dbNSFP4.9a.readme.txt"), readme).unwrap();
        let mut params = vec![file];
        params.extend_from_slice(names);
        load(Some(Assembly::Grch38), &params).unwrap().header_info()
    }

    /// A `0` name and a `0` text are kept where the Ensembl VEP plugin's truth test
    /// (`$descs{$col} = $2 if $col && $2`, dbNSFP.pm:523) drops them: it prints
    /// `SIFT_score from dbNSFP file` and `SIFT_pred from dbNSFP file` for the two
    /// `0` texts here and, for `chr`, its continuation alone behind two blanks. `0.0`,
    /// `00` and `0 ` are true in Perl and print in both engines, the last with its
    /// trailing blank in Ensembl VEP.
    #[test]
    fn readme_zero_name_or_text_is_kept() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let readme = concat!(
            "A made-up list of columns:\n",
            "1\tchr: 0\n",
            "\t\tContinuation after a zero description.\n",
            "11\tSIFT_score: 0\n",
            "12\t0: Description of a column named zero.\n",
            "\t\tContinuation after a zero name.\n",
            "14\tMutationTaster_score: 0.0\n",
            "16\tGERP++_RS: 00\n",
            "13\tPolyphen2_HDIV_score: 0 \n",
            "12\tSIFT_pred: 0\n",
        );
        let names = [
            "chr",
            "SIFT_score",
            "SIFT_pred",
            "MutationTaster_score",
            "GERP++_RS",
            "Polyphen2_HDIV_score",
        ];
        let from = "(from dbNSFP4.9a_grch38_synthetic)";
        assert_eq!(
            described(dir.path(), &file, readme.as_bytes(), &names),
            vec![
                (
                    "chr".to_string(),
                    format!("{from} 0 Continuation after a zero description.")
                ),
                ("SIFT_score".to_string(), format!("{from} 0")),
                ("SIFT_pred".to_string(), format!("{from} 0")),
                ("MutationTaster_score".to_string(), format!("{from} 0.0")),
                ("GERP++_RS".to_string(), format!("{from} 00")),
                ("Polyphen2_HDIV_score".to_string(), format!("{from} 0")),
            ]
        );
        assert_eq!(
            readme_descriptions(Path::new(&file))
                .get("0")
                .map(String::as_str),
            Some("Description of a column named zero. Continuation after a zero name.")
        );
    }

    /// A numbered line's name ends at the first colon that a blank and one more
    /// byte follow, where the lazy `(.+?)\:\s+(.+)` (dbNSFP.pm:520) ends it: the
    /// `SIFT_score:x: ...` line names `SIFT_score:x`, so both engines print
    /// `SIFT_score from dbNSFP file`; a colon inside the text (`a:b pair`) and a text
    /// opening with a colon behind a tab are text.
    #[test]
    fn readme_name_ends_at_the_first_colon_a_blank_follows() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let readme = concat!(
            "A made-up list of columns:\n",
            "1\tchr: A made-up text for the first column.\n",
            "11\tSIFT_score:x: SIFT score with a colon inside its name.\n",
            "\t\tContinuation of the lazy name.\n",
            "12\tSIFT_pred: a:b pair\n",
            "13\tPolyphen2_HDIV_score:\t: description behind a tab-colon pair\n",
        );
        let names = ["chr", "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score"];
        let from = "(from dbNSFP4.9a_grch38_synthetic)";
        assert_eq!(
            described(dir.path(), &file, readme.as_bytes(), &names),
            vec![
                (
                    "chr".to_string(),
                    format!("{from} A made-up text for the first column.")
                ),
                (
                    "SIFT_score".to_string(),
                    "SIFT_score from dbNSFP file".to_string()
                ),
                ("SIFT_pred".to_string(), format!("{from} a:b pair")),
                (
                    "Polyphen2_HDIV_score".to_string(),
                    format!("{from} : description behind a tab-colon pair")
                ),
            ]
        );
        assert_eq!(
            readme_descriptions(Path::new(&file))
                .get("SIFT_score:x")
                .map(String::as_str),
            Some("SIFT score with a colon inside its name. Continuation of the lazy name.")
        );
    }

    /// The README is read as bytes, so the plugin's `\s`, `\w` and `\d` are the
    /// ASCII classes (dbNSFP.pm:517-531): a no-break space (U+00A0) or a next line
    /// (U+0085) after a number opens no numbered line, so the line is a continuation;
    /// a no-break space inside a text or opening a continuation is kept; a tab and a
    /// vertical tab collapse to one space; a continuation line without an ASCII word
    /// character ends the entry's continuations. Ensembl VEP prints every row here as
    /// expected except `SIFT_score` and `GERP++_RS`, whose colon a no-break space
    /// follows: it prints `SIFT_score from dbNSFP file` and `GERP++_RS from dbNSFP
    /// file` where vep-rs reads the text behind the colon, no-break spaces included
    /// (the documented divergence `dbnsfp_readme_description_text`).
    #[test]
    fn readme_blanks_are_the_ascii_ones() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let readme: &[u8] = b"A made-up list of columns:\n\
            1\tchr: A made-up text for the first column.\n\
            \xc2\xa0\xc2\xa0Continuation opening with no-break spaces.\n\
            12\xc2\x85SIFT_pred: a next-line character after the number.\n\
            11\tSIFT_score:\xc2\xa0SIFT score behind a no-break space.\n\
            14\tMutationTaster_score: tab\tinside, vt\x0binside, no-break\xc2\xa0inside\n\
            16\xc2\xa0GERP++_RS: a no-break space after the number.\n\
            13\tPolyphen2_HDIV_score: Polyphen2 score\n\
            \t\t\xc3\xa9\n\
            \t\tContinuation after a line of one non-ASCII letter.\n\
            16\tGERP++_RS:\xc2\xa0\xc2\xa0two no-break spaces then text\n";
        let names = [
            "chr",
            "SIFT_pred",
            "SIFT_score",
            "MutationTaster_score",
            "Polyphen2_HDIV_score",
            "GERP++_RS",
        ];
        let from = "(from dbNSFP4.9a_grch38_synthetic)";
        assert_eq!(
            described(dir.path(), &file, readme, &names),
            vec![
                (
                    "chr".to_string(),
                    format!("{from} A made-up text for the first column. \u{a0}\u{a0}Continuation opening with no-break spaces. 12\u{85}SIFT_pred: a next-line character after the number.")
                ),
                ("SIFT_pred".to_string(), "SIFT_pred from dbNSFP file".to_string()),
                (
                    "SIFT_score".to_string(),
                    format!("{from} \u{a0}SIFT score behind a no-break space.")
                ),
                (
                    "MutationTaster_score".to_string(),
                    format!("{from} tab inside, vt inside, no-break\u{a0}inside 16\u{a0}GERP++_RS: a no-break space after the number.")
                ),
                (
                    "Polyphen2_HDIV_score".to_string(),
                    format!("{from} Polyphen2 score")
                ),
                (
                    "GERP++_RS".to_string(),
                    format!("{from} \u{a0}\u{a0}two no-break spaces then text")
                ),
            ]
        );
    }

    /// A README that is not UTF-8 is read: the Ensembl VEP plugin reads and prints
    /// bytes (`Zyxw\xe4rt`, `\xfc`), where vep-rs prints the replacement character
    /// for each byte sequence that is not UTF-8 and the rest of the text as the
    /// plugin prints it.
    #[test]
    fn readme_that_is_not_utf8_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let readme: &[u8] = b"A made-up list of columns:\n\
            1\tchr: A made-up text for the first column.\n\
            11\tSIFT_score: SIFT score from Zyxw\xe4rt et al.\n\
            \t\tContinuation with a Latin-1 \xfc.\n\
            12\tSIFT_pred: prediction\n";
        let from = "(from dbNSFP4.9a_grch38_synthetic)";
        assert_eq!(
            described(dir.path(), &file, readme, &["chr", "SIFT_score", "SIFT_pred"]),
            vec![
                ("chr".to_string(), format!("{from} A made-up text for the first column.")),
                (
                    "SIFT_score".to_string(),
                    format!("{from} SIFT score from Zyxw\u{fffd}rt et al. Continuation with a Latin-1 \u{fffd}.")
                ),
                ("SIFT_pred".to_string(), format!("{from} prediction")),
            ]
        );
    }

    /// Text with no blank between the colon and it describes the column (the
    /// documented divergence `dbnsfp_readme_description_text`): Ensembl VEP's
    /// `\:\s+(.+)` (dbNSFP.pm:520) opens no entry for `MutationTaster_score:text` or
    /// `GERP++_RS:x:y:z` and prints `<column> from dbNSFP file` for both. A colon
    /// with nothing or one blank behind it opens none in either engine. The
    /// fixture README's lines end in `\r\n`, the `12` and `GERP++_RS:x:y:z`
    /// lines in `\r\r\n`, of which `s/\r$//` (dbNSFP.pm:515) strips one carriage
    /// return and leaves the other as a blank. Behind `x:y:z` that blank is a
    /// trailing one `collapse_blanks` drops; behind the bare `12` it is the blank
    /// that makes the line a numbered one, which opens no entry, so its
    /// continuation is dropped. Two blanks and nothing more open an entry with no
    /// text, printed as the prefix alone where Ensembl VEP prints the prefix and
    /// two blanks.
    #[test]
    fn readme_text_with_no_blank_after_the_colon_describes_the_column() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let readme = concat!(
            "A made-up list of columns:\r\n",
            "1\tchr: \r\n",
            "14\tMutationTaster_score:text\r\n",
            "\t\tContinuation behind the lenient line.\r\n",
            "12\r\r\n",
            "\t\tContinuation after a bare number.\r\n",
            "16\tGERP++_RS:x:y:z\r\r\n",
            "12\tSIFT_pred:\r\n",
            "13\tPolyphen2_HDIV_score:  \r\n",
        );
        let names = [
            "chr",
            "MutationTaster_score",
            "GERP++_RS",
            "SIFT_pred",
            "Polyphen2_HDIV_score",
        ];
        let from = "(from dbNSFP4.9a_grch38_synthetic)";
        assert_eq!(
            described(dir.path(), &file, readme.as_bytes(), &names),
            vec![
                ("chr".to_string(), "chr from dbNSFP file".to_string()),
                (
                    "MutationTaster_score".to_string(),
                    format!("{from} text Continuation behind the lenient line.")
                ),
                ("GERP++_RS".to_string(), format!("{from} x:y:z")),
                (
                    "SIFT_pred".to_string(),
                    "SIFT_pred from dbNSFP file".to_string()
                ),
                ("Polyphen2_HDIV_score".to_string(), from.to_string()),
            ]
        );
    }

    /// A dropped continuation that says "separated by" marks no column
    /// transcript-specific, as in Ensembl VEP, where the dropped line is appended
    /// under an undefined key (dbNSFP.pm:527) the transcript-specific scan
    /// (dbNSFP.pm:471-476) matches to no column: `GERP++_RS`, whole in the detection
    /// row, stays whole and prints on a row whose transcript the record does not
    /// list, where a transcript-specific column is dropped (dbNSFP.pm:364-366).
    #[test]
    fn dropped_continuation_marks_no_column_transcript_specific() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let readme = concat!(
            "A made-up list of columns:\n",
            "1\tchr: A made-up text for the first column.\n",
            "16\tGERP++_RS: GERP++ RS score\n",
            "17\tGERP++_NR\n",
            "\t\tA made-up continuation whose values are separated by \";\".\n",
            "11\tSIFT_score: SIFT score\n",
        );
        std::fs::write(dir.path().join("dbNSFP4.9a.readme.txt"), readme).unwrap();
        let plugin = load(
            Some(Assembly::Grch38),
            &[&file, "transcript_match=1", "GERP++_RS", "SIFT_score"],
        )
        .unwrap();
        let specific: Vec<&str> = plugin
            .selected
            .iter()
            .filter(|c| c.transcript_specific)
            .map(|c| c.name.as_str())
            .collect();
        assert_eq!(specific, ["SIFT_score"]);
        assert_eq!(
            plugin.header_info()[0],
            (
                "GERP++_RS".to_string(),
                "(from dbNSFP4.9a_grch38_synthetic) GERP++ RS score".to_string()
            )
        );
        let mut t = variant("T");
        let got = annotated(
            &plugin,
            &mut t,
            vec![consequence(
                "ENST1",
                Consequence::MissenseVariant,
                Some("R/W"),
            )],
        );
        assert_eq!(pairs(&got[0]), [("GERP++_RS", "4.5")]);
    }

    /// A column the header lacks is kept, warned about and written as `invalid_field`
    /// on every annotated row (dbNSFP.pm:239-249, 403-405); a parameter list that names
    /// no column is refused (dbNSFP.pm:244-246).
    #[test]
    fn unknown_column_prints_invalid_field_and_no_column_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let plugin = load(
            Some(Assembly::Grch38),
            &[&file, "SIFT_score", "NoSuchColumn"],
        )
        .unwrap();
        let mut c = variant("C");
        let got = annotated(
            &plugin,
            &mut c,
            vec![
                consequence("ENST1", Consequence::MissenseVariant, Some("R/G")),
                consequence("ENST2", Consequence::IntronVariant, None),
            ],
        );
        assert_eq!(
            pairs(&got[0]),
            [
                ("SIFT_score", "0.01,0.02"),
                ("NoSuchColumn", "invalid_field")
            ]
        );
        assert!(got[1].is_empty());
        let err = load(Some(Assembly::Grch38), &[&file]).unwrap_err();
        assert!(
            err.to_string()
                .contains("no input columns were found in file header"),
            "{err}"
        );
    }

    /// Only a single-base substitution whose allele is a base is annotated
    /// (dbNSFP.pm:297, 303), on the forward strand (dbNSFP.pm:301), and only on
    /// transcript rows (`feature_types` of `['Transcript']`).
    #[test]
    fn single_base_substitutions_on_transcript_rows_alone() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let plugin = load(
            Some(Assembly::Grch38),
            &[&file, "pep_match=0", "consequence=ALL", "SIFT_score"],
        )
        .unwrap();
        let mut deletion = InputVariant::new("21".into(), 100, 101, b"AC".to_vec(), b"A".to_vec());
        let got = annotated(
            &plugin,
            &mut deletion,
            vec![consequence("ENST1", Consequence::FrameshiftVariant, None)],
        );
        assert!(got[0].is_empty());
        let mut reverse = InputVariant::new("21".into(), 100, 100, b"T".to_vec(), b"G".to_vec());
        reverse.strand = vep_core::coordinate::Strand::Reverse;
        let got = annotated(
            &plugin,
            &mut reverse,
            vec![consequence("ENST1", Consequence::MissenseVariant, None)],
        );
        assert_eq!(
            got[0].get("SIFT_score").map(String::as_str),
            Some("0.01,0.02")
        );
        let mut c = variant("C");
        let mut regulatory = consequence("ENSR1", Consequence::RegulatoryRegionVariant, None);
        regulatory.feature_type = FeatureType::RegulatoryFeature;
        let got = annotated(
            &plugin,
            &mut c,
            vec![
                regulatory,
                consequence("ENST1", Consequence::IntronVariant, None),
            ],
        );
        assert!(got[0].is_empty());
        assert_eq!(
            got[1].get("SIFT_score").map(String::as_str),
            Some("0.01,0.02")
        );
        assert!(c.plugin_data.is_empty());
    }

    /// Ensembl VEP 116.3 never hands dbNSFP a structural variant: the module
    /// keeps the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2229). A one-base `<DEL>` row, whose
    /// start equals its end as a substitution's does, gets no field, and nor
    /// does a matching substitution carrying the structural mark, which pins
    /// the gate rather than the span test; the plain substitution is annotated.
    #[test]
    fn a_structural_variant_row_gets_no_field() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(dir.path());
        let plugin = load(
            Some(Assembly::Grch38),
            &[&file, "pep_match=0", "consequence=ALL", "SIFT_score"],
        )
        .unwrap();
        assert!(!plugin.annotates_structural_variants());
        let missense = || vec![consequence("ENST1", Consequence::MissenseVariant, None)];
        let mut deletion = variant("<DEL>");
        deletion.ref_allele = b"N".to_vec();
        deletion.is_structural = true;
        assert!(annotated(&plugin, &mut deletion, missense())[0].is_empty());
        let mut flagged = variant("C");
        flagged.is_structural = true;
        assert!(annotated(&plugin, &mut flagged, missense())[0].is_empty());
        let mut c = variant("C");
        assert_eq!(
            annotated(&plugin, &mut c, missense())[0]
                .get("SIFT_score")
                .map(String::as_str),
            Some("0.01,0.02")
        );
    }

    /// dbNSFP.pm:426-428 split the line htslib hands over, its `\r` alone
    /// removed, on `\t` and key the cells by header name; :408 skips a value
    /// equal to `.` and :417 copies the rest as they stand: a last cell
    /// `4.56 ` is written with its trailing space, `. ` is written since it is
    /// not `.`, and a cell holding a single space is written as that space.
    #[test]
    fn a_trailing_space_on_the_last_cell_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let ending = |alt: &str, aa: (&str, &str), gerp: &str| {
            let full = row(("21", 100, 200), ("A", alt), aa, "ENST1", "0.5", ".");
            let (head, _) = full
                .rsplit_once('\t')
                .expect("GERP++_RS is the last column");
            format!("{head}\t{gerp}")
        };
        let file = grch38_file(
            dir.path(),
            &[
                detection_row(),
                ending("C", ("R", "G"), "4.56 "),
                ending("G", ("R", "W"), ". "),
                ending("T", ("R", "C"), " "),
            ],
        );
        let plugin = load(Some(Assembly::Grch38), &[&file, "GERP++_RS"]).unwrap();
        let gerp = |alt: &str, change: &str| {
            let rows = annotated(
                &plugin,
                &mut variant(alt),
                vec![consequence(
                    "ENST1",
                    Consequence::MissenseVariant,
                    Some(change),
                )],
            );
            rows[0].get("GERP++_RS").cloned()
        };
        assert_eq!(gerp("C", "R/G"), Some("4.56 ".to_string()));
        assert_eq!(gerp("G", "R/W"), Some(". ".to_string()));
        assert_eq!(gerp("T", "R/C"), Some(" ".to_string()));
    }

    /// dbNSFP.pm:407-408 skip an undefined value and `.` alone, so an empty
    /// cell before the row's last filled one is written as the empty value,
    /// Ensembl VEP's defined empty string (an empty tab cell and `""` in JSON,
    /// Tab.pm:128 and OutputFactory/JSON.pm:327), which the row names in
    /// `plugin_empty_strings`; under `transcript_match` the transcript's own
    /// entry is what `:388` writes back and `:408` tests, so a `.` entry of a
    /// `;`-joined list is not written while a sibling transcript's `0.5` is.
    #[test]
    fn an_empty_cell_is_written_and_a_dot_list_entry_is_not() {
        let dir = tempfile::tempdir().unwrap();
        let sparse_row = row(
            ("21", 100, 200),
            ("A", "C"),
            ("R", "G"),
            "ENST1;ENST2",
            ".;0.5",
            ".",
        )
        .replace("\tD;T\t", "\t\t");
        let file = grch38_file(dir.path(), &[detection_row(), sparse_row]);
        let plugin = load(
            Some(Assembly::Grch38),
            &[&file, "transcript_match=1", "SIFT_score", "SIFT_pred"],
        )
        .unwrap();
        let mut c = variant("C");
        let rows = annotated(
            &plugin,
            &mut c,
            vec![
                consequence("ENST1", Consequence::MissenseVariant, Some("R/G")),
                consequence("ENST2", Consequence::MissenseVariant, Some("R/G")),
            ],
        );
        assert_eq!(pairs(&rows[0]), [("SIFT_pred", "")]);
        assert_eq!(pairs(&rows[1]), [("SIFT_score", "0.5"), ("SIFT_pred", "")]);
        let marks: Vec<&[String]> = c
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_empty_strings.as_slice())
            .collect();
        assert_eq!(
            marks,
            [["SIFT_pred".to_string()], ["SIFT_pred".to_string()]]
        );

        // An earlier plugin's mark under a key this plugin writes a value for
        // goes with the value it replaces (OutputFactory.pm:2239-2240); a mark
        // under a key it leaves alone stays.
        let mut marked = consequence("ENST2", Consequence::MissenseVariant, Some("R/G"));
        marked.plugin_empty_strings = strings(&["SIFT_score", "Other"]);
        let mut c = variant("C");
        let rows = annotated(&plugin, &mut c, vec![marked]);
        assert_eq!(pairs(&rows[0]), [("SIFT_score", "0.5"), ("SIFT_pred", "")]);
        assert_eq!(
            c.transcript_consequences[0].plugin_empty_strings,
            strings(&["Other", "SIFT_pred"])
        );
    }

    fn strings(texts: &[&str]) -> Vec<String> {
        texts.iter().map(|text| text.to_string()).collect()
    }

    /// A row's position cell is read as `$pos == $vf->{start}` (dbNSFP.pm:336)
    /// numifies it: `100.0`, `100x` and `1e2` are 100 and the row annotates
    /// the variant there; `100.5`, `abc` and `0` are no position. The tabix
    /// index places a cell at its leading integer, so the first file shapes
    /// reach the match and the exponent and the rest are read from records
    /// built here.
    #[test]
    fn a_position_cell_is_matched_as_perl_numifies_it() {
        let dir = tempfile::tempdir().unwrap();
        // A made-up row of `HEADER`, its position cell as given.
        let cell_row = |pos: &str, alt: &str, sift: &str| {
            format!("21\t{pos}\tA\t{alt}\tA\tG\trs1\t21\t100\tENST1\t{sift}\tD\t0.9\t1\t.\t4.5")
        };
        let rows = [
            detection_row(),
            cell_row("100.0", "G", "0.1"),
            cell_row("100.5", "T", "0.2"),
            cell_row("100x", "C", "0.3"),
        ];
        let file = grch38_file(dir.path(), &rows);
        let plugin = load(Some(Assembly::Grch38), &[&file, "SIFT_score"]).unwrap();
        let sift = |alt: &str| {
            let fields = annotated(
                &plugin,
                &mut variant(alt),
                vec![consequence(
                    "ENST1",
                    Consequence::MissenseVariant,
                    Some("A/G"),
                )],
            );
            fields[0].get("SIFT_score").cloned()
        };
        assert_eq!(sift("G"), Some("0.1".to_string()));
        assert_eq!(sift("T"), None);
        assert_eq!(sift("C"), Some("0.3".to_string()));

        let record = |pos: &str| TabixRecord {
            chr: "21".to_string(),
            start: 100,
            end: 100,
            columns: cell_row(pos, "G", "0.9")
                .split('\t')
                .map(str::to_string)
                .collect(),
            header: None,
        };
        let matched = |pos: &str| plugin.record_matches(&record(pos), 100, b"G", Some("A/G"));
        assert!(matched("1e2"));
        assert!(matched(" 100"));
        assert!(!matched("abc"));
        assert!(!matched("0"));
        assert!(!matched("100.5"));
    }

    /// The module names a record on `MT` or `chrMT` `M` to its table
    /// (dbNSFP.pm:313). On GRCh37, with the run's reference sequence, a table
    /// naming its mitochondrial contig `M` is read for such a record when its
    /// rows state the reference bases at their own positions, and not when
    /// they state the bases one position on: there the row at the record's
    /// site carries the record's alternate allele and amino acids, and its
    /// value is the neighbouring base's. Without a reference sequence neither
    /// table is read for such a record; a record on `chrM` or `M` reads both,
    /// and a table naming the contig `MT` is read for a record on `MT`
    /// whatever its rows state.
    #[test]
    fn on_grch37_a_table_named_m_is_read_for_a_record_on_mt_only_when_its_rows_state_the_reference()
    {
        use crate::test_fixtures::{
            made_up_alternate, made_up_reference, made_up_rows, repeated_base,
        };
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference(dir.path(), 1600);
        let site = repeated_base(&bases, 100);
        let base = bases[site - 1];
        let sift = |plugin: &DbNsfpPlugin, spelling: &str| -> Option<String> {
            let mut record = InputVariant::new(
                spelling.into(),
                site as u64,
                site as u64,
                vec![base],
                vec![made_up_alternate(base)],
            );
            annotated(
                plugin,
                &mut record,
                vec![consequence(
                    "ENST1",
                    Consequence::MissenseVariant,
                    Some("R/G"),
                )],
            )[0]
            .get("SIFT_score")
            .cloned()
        };
        let table = |name: &str, contig: &str, offset: usize| -> String {
            let mut rows = vec![detection_row()];
            rows.extend(made_up_rows(
                contig,
                &bases,
                10..1590,
                offset,
                |contig, position, reference, alt| {
                    row(
                        (contig, position as u64, position as u64),
                        (&reference.to_string(), &alt.to_string()),
                        ("R", "G"),
                        "ENST1",
                        &format!("0.{position}"),
                        ".",
                    )
                },
            ));
            let mut lines: Vec<&str> = vec![HEADER];
            lines.extend(rows.iter().map(String::as_str));
            write_tabix_fixture(
                dir.path(),
                name,
                &lines,
                IndexSpec::standard().with_span_cols(8, 9, 9),
            )
            .to_string_lossy()
            .to_string()
        };
        let loaded = |file: &str, fasta: Option<Arc<vep_fasta::IndexedFasta>>| {
            let mut plugin = load(Some(Assembly::Grch37), &[file, "SIFT_score"]).unwrap();
            plugin.set_reference_fasta(fasta);
            plugin
        };
        let own = Some(format!("0.{site}"));

        let stating = table("dbNSFP4.9a_stating_the_reference_synthetic.gz", "M", 0);
        let one_on = table("dbNSFP4.9a_one_position_on_synthetic.gz", "M", 1);
        for (file, with_reference) in [(&stating, own.clone()), (&one_on, None)] {
            let plugin = loaded(file, Some(Arc::clone(&reference)));
            for record in ["MT", "chrMT"] {
                assert_eq!(sift(&plugin, record), with_reference, "{file}, {record}");
            }
            for record in ["chrM", "M"] {
                assert_eq!(sift(&plugin, record), own, "{file}, {record}");
            }
            let plugin = loaded(file, None);
            assert_eq!(sift(&plugin, "MT"), None, "{file} without a reference");
            assert_eq!(sift(&plugin, "chrM"), own, "{file} without a reference");
        }

        let named_mt = table("dbNSFP4.9a_named_mt_one_position_on_synthetic.gz", "MT", 1);
        let plugin = loaded(&named_mt, None);
        assert_eq!(sift(&plugin, "MT"), own);
        assert_eq!(sift(&plugin, "chrM"), None);
    }
}
