// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Batch-oriented tabix query engine.
//!
//! [`TabixAnnotator`] wraps `noodles::tabix` to provide region-batched queries:
//! a batch's regions are sorted by start within each chromosome and merged
//! (bridging gaps of up to 1000 bp), and each merged region is read through one
//! index lookup, so the per-variant lookups that follow are in-memory range
//! lookups rather than I/O.
//!
//! Ensembl VEP source citations name modules of Ensembl/VEP_plugins release/116, or
//! of ensembl-variation release/116 (`Utils/Sequence.pm`, `BaseVepTabixPlugin.pm`).

use std::borrow::Cow;
use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use noodles::bgzf;
use noodles::core::Position;
use noodles::csi::binning_index::index::header::format::CoordinateSystem;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::binning_index::BinningIndex;
use noodles::tabix;
use rayon::prelude::*;
use tracing::debug;
use vep_core::assembly::Assembly;
use vep_core::variant::InputVariant;

use crate::PluginError;

/// The longest line read from a tabix-indexed file or a plugin's catalogue,
/// 64 MiB. A line of a population VCF with every INFO key runs to a few
/// hundred kilobytes, so the cap sits far above any record; a file with no
/// line end in it is refused at the cap instead of being read whole into one
/// line buffer.
#[doc(hidden)]
pub const MAX_LINE_BYTES: usize = 64 << 20;

/// The longest index (`.tbi` or `.csi`) read, 256 MiB decompressed. A
/// whole-genome index decompresses to a few megabytes; the cap bounds what an
/// index whose counts are then checked ([`check_index_counts`]) can make the
/// reader hold.
const MAX_INDEX_BYTES: usize = 256 << 20;

/// `BufRead::read_until(b'\n', line)` with the line's length capped at
/// [`MAX_LINE_BYTES`]: the bytes appended, 0 at the end of the input, or an
/// error built with `kind` naming `path`, for a read failure or, with the
/// line's place in the file, for a line still open past the cap.
fn read_line_bounded<R: BufRead>(
    reader: &mut R,
    line: &mut Vec<u8>,
    path: &Path,
    kind: fn(String) -> PluginError,
    place: impl FnOnce() -> String,
) -> Result<usize, PluginError> {
    // The cap plus the line end: a line of exactly the cap is read whole.
    let limit = MAX_LINE_BYTES as u64 + 1;
    let read = reader
        .take(limit)
        .read_until(b'\n', line)
        .map_err(|e| kind(format!("failed to read {}: {e}", path.display())))?;
    if read as u64 == limit && line.last() != Some(&b'\n') {
        return Err(kind(format!(
            "{}: {} is longer than the {MAX_LINE_BYTES}-byte limit for one line",
            path.display(),
            place()
        )));
    }
    Ok(read)
}

/// The lines of a text file as bytes, each without its line end (`\n` or
/// `\r\n`; the `\n` alone through `Self::chomped`), read one at a time with
/// the length of a line capped at [`MAX_LINE_BYTES`]: through the iterator,
/// each line its own `Vec`, or through [`Self::next_line`], one buffer reused.
/// A line is read as Perl's `<$fh>` reads one from a handle opened without an
/// encoding layer, which is how every Ensembl VEP plugin opens its flat file:
/// a byte sequence that is not UTF-8 is a sequence of bytes like any other,
/// split and matched as such, and a caller makes text of the cells it keeps
/// with [`lossy_text`] or `String::from_utf8_lossy`. Every error is
/// [`PluginError::Init`]: a plugin reads such a file while it opens, whether
/// the file's head alone or a whole catalogue.
#[doc(hidden)]
pub struct BoundedLines<R> {
    reader: R,
    path: PathBuf,
    line_number: u64,
    /// Strip the `\n` alone, as Perl's `chomp` does, so a CRLF file's line
    /// keeps its `\r`.
    chomp: bool,
    /// The line [`Self::next_line`] hands out, reused across lines.
    line: Vec<u8>,
}

impl<R: BufRead> BoundedLines<R> {
    /// The lines of `reader`; `path` is not opened, it names the file in
    /// every error.
    pub fn new(reader: R, path: &Path) -> Self {
        Self {
            reader,
            path: path.to_path_buf(),
            line_number: 0,
            chomp: false,
            line: Vec::new(),
        }
    }

    /// The lines of `reader` as an Ensembl VEP module's `chomp` leaves them: the `\n`
    /// stripped and a `\r` before it kept, for a caller whose `split` then
    /// treats that `\r` as a blank.
    pub(crate) fn chomped(reader: R, path: &Path) -> Self {
        Self {
            chomp: true,
            ..Self::new(reader, path)
        }
    }

    /// The next line, in a buffer of this reader's that the next call
    /// overwrites: `None` at the end of the input. A loader that keeps a cell
    /// or two of every line copies them out and allocates nothing per line,
    /// where the iterator hands each line over as its own `Vec`.
    pub fn next_line(&mut self) -> Option<Result<&[u8], PluginError>> {
        self.line.clear();
        self.line_number += 1;
        let number = self.line_number;
        match read_line_bounded(
            &mut self.reader,
            &mut self.line,
            &self.path,
            PluginError::Init,
            || format!("line {number}"),
        ) {
            Ok(0) => None,
            Ok(_) => {
                strip_line_end(&mut self.line, self.chomp);
                Some(Ok(&self.line))
            }
            Err(e) => Some(Err(e)),
        }
    }
}

/// Strips the line end `BufRead::lines` strips, `\n` or `\r\n`, or under
/// `chomp` the `\n` alone.
fn strip_line_end(line: &mut Vec<u8>, chomp: bool) {
    if line.last() == Some(&b'\n') {
        line.pop();
        if !chomp && line.last() == Some(&b'\r') {
            line.pop();
        }
    }
}

impl<R: BufRead> Iterator for BoundedLines<R> {
    type Item = Result<Vec<u8>, PluginError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = Vec::new();
        self.line_number += 1;
        let number = self.line_number;
        match read_line_bounded(
            &mut self.reader,
            &mut line,
            &self.path,
            PluginError::Init,
            || format!("line {number}"),
        ) {
            Ok(0) => None,
            Ok(_) => {
                strip_line_end(&mut line, self.chomp);
                Some(Ok(line))
            }
            Err(e) => Some(Err(e)),
        }
    }
}

/// The text of a line, or a cell, read from a flat file: valid UTF-8 as it
/// stands, moved rather than copied, and each byte sequence that is not UTF-8
/// as one U+FFFD, the form the query path gives a record's cells. An Ensembl
/// VEP plugin prints such a byte as it stands; vep-rs writes UTF-8, so the
/// replacement character marks the byte's place in the output.
#[doc(hidden)]
pub fn lossy_text(bytes: Vec<u8>) -> String {
    match String::from_utf8(bytes) {
        Ok(text) => text,
        Err(error) => String::from_utf8_lossy(error.as_bytes()).into_owned(),
    }
}

/// A parsed record from a tabix-indexed file.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct TabixRecord {
    /// Chromosome / sequence name.
    pub chr: String,
    /// 1-based start position.
    pub start: u64,
    /// 1-based end position (same as start for SNVs in most formats).
    pub end: u64,
    /// Raw column values from the tab-separated line.
    pub columns: Vec<String>,
    /// Column names (from header), if available.
    /// Wrapped in `Arc` to avoid cloning the header for every record.
    pub header: Option<Arc<Vec<String>>>,
}

impl TabixRecord {
    /// A record, from each of its fields.
    pub fn new(
        chr: String,
        start: u64,
        end: u64,
        columns: Vec<String>,
        header: Option<Arc<Vec<String>>>,
    ) -> Self {
        Self {
            chr,
            start,
            end,
            columns,
            header,
        }
    }

    /// Get a column value by header name.
    pub fn get(&self, name: &str) -> Option<&str> {
        let header = self.header.as_ref()?;
        let idx = header.iter().position(|h| h == name)?;
        self.columns.get(idx).map(|s| s.as_str())
    }

    /// Extract multiple fields into an IndexMap.
    pub fn extract_fields(&self, names: &[&str]) -> indexmap::IndexMap<String, String> {
        let mut result = indexmap::IndexMap::new();
        for name in names {
            if let Some(val) = self.get(name) {
                if !val.is_empty() && val != "." {
                    result.insert(name.to_string(), val.to_string());
                }
            }
        }
        result
    }
}

/// The number of fields Ensembl VEP's `split /\t/` reads from `cells`, a line
/// split on tabs: the trailing empty cells are dropped and a whitespace-only
/// cell is not empty (TranscriptAnnotator.pm:245-246; :109-110 for a header),
/// so a cell past the last filled one is undefined where an empty cell before
/// a filled one is a field holding the empty string.
pub(crate) fn split_field_count(cells: &[String]) -> usize {
    cells
        .iter()
        .rposition(|cell| !cell.is_empty())
        .map_or(0, |last| last + 1)
}

/// Result of a batch tabix query, indexed for fast per-variant lookup.
///
/// Records are stored in a `HashMap<chr, BTreeMap<start, Vec<TabixRecord>>>`
/// so that per-variant lookups are O(log N) BTreeMap range queries instead
/// of linear scans. A record is found by every range its `start..=end` span
/// overlaps, not only by the ranges containing its start.
pub struct BatchQueryResult {
    /// Records indexed by chromosome and start position.
    records_by_pos: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>>,
    /// Per chromosome, the longest `end - start` of any record: how far below a
    /// range's start a record keyed by its start can lie and still reach into
    /// the range.
    longest_span: HashMap<String, u64>,
}

impl BatchQueryResult {
    /// Create an empty batch result.
    pub(crate) fn empty() -> Self {
        Self::from_records(HashMap::new())
    }

    /// Create a batch result from pre-built records indexed by chromosome and position.
    pub fn from_records(records_by_pos: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>>) -> Self {
        let longest_span = records_by_pos
            .iter()
            .map(|(chr, by_start)| {
                let longest = by_start
                    .values()
                    .flatten()
                    .map(|record| record.end.saturating_sub(record.start))
                    .max()
                    .unwrap_or(0);
                (chr.clone(), longest)
            })
            .collect();
        Self {
            records_by_pos,
            longest_span,
        }
    }

    /// Look up all records whose span overlaps a given range.
    pub fn lookup(&self, chr: &str, start: u64, end: u64) -> Vec<&TabixRecord> {
        let Some(chr_map) = self.records_by_pos.get(chr) else {
            return Vec::new();
        };
        let (start, end) = normalize_range_bounds(start, end);
        let reach = self.longest_span.get(chr).copied().unwrap_or(0);

        chr_map
            .range(start.saturating_sub(reach)..=end)
            .flat_map(|(_start, records)| records.iter())
            .filter(|record| record.end >= start)
            .collect()
    }

    /// Check if the result is empty.
    pub(crate) fn is_empty(&self) -> bool {
        self.records_by_pos.is_empty()
    }
}

/// Configuration for how to parse a tabix file's columns.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct TabixConfig {
    /// Path to the tabix-indexed file.
    pub file_path: PathBuf,
    /// 0-based column index for chromosome (default: 0).
    pub chr_col: usize,
    /// 0-based column index for start position (default: 1).
    pub start_col: usize,
    /// 0-based column index for end position. If None, end = start.
    ///
    /// A record is matched by its `start..=end` span, as `tabix` matches it by
    /// the columns its index was built on (`-b` and `-e`), so this names the
    /// index's `-e` column and is `None` for an index built with `-b` and `-e`
    /// on the same column. With a `ref_col` the column is a VCF `INFO` column,
    /// read for its `END` key as `tabix -p vcf` reads it; without one it is a
    /// generic end column, read to its leading integer ([`Self::record_end`]).
    pub end_col: Option<usize>,
    /// 0-based column index for reference allele, if applicable.
    pub ref_col: Option<usize>,
    /// 0-based column index for alt allele, if applicable.
    pub alt_col: Option<usize>,
    /// Whether positions in the file are 0-based (like BED). Default: false (1-based).
    pub zero_based: bool,
}

impl Default for TabixConfig {
    /// The documented defaults: the chromosome in column 0, a 1-based start in
    /// column 1 that is the record's whole span, no allele columns and no path.
    fn default() -> Self {
        Self {
            file_path: PathBuf::new(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        }
    }
}

impl TabixConfig {
    /// The 1-based start of a row from its start cell, read as htslib's
    /// `tbx_parse1` reads the cell with `strtol` (`strtoll` from htslib 1.10):
    /// the leading decimal integer after any blanks and a `+`, whatever follows
    /// it (`100x` and `100.5` are 100), so a row carrying such a cell is
    /// indexed and returned at that position, where a module's own reading of
    /// the cell (`$_->{start} == $start`, BayesDel.pm:141) is the module's.
    /// `strtol` also takes a `-` and a hex or octal prefix (`0x1A`, `0103`);
    /// here a `-` leads no integer (a reader
    /// places such a row only under
    /// `TabixAnnotator::reading_rows_as_tabix_indexes_them`) and a prefix
    /// reads as decimal digits. `None` for a cell no integer leads, which htslib
    /// refuses to index, and for a zero-based start with no one-based form.
    pub fn record_start(&self, cell: &str) -> Option<u64> {
        let start = leading_integer(cell)?;
        if self.zero_based {
            start.checked_add(1)
        } else {
            Some(start)
        }
    }

    /// The 1-based end of a row starting at `start` from its `columns`, the
    /// span `tabix` matches the row by. Without an end column, or with the
    /// start column as the end column, the row ends at its start. A
    /// configuration naming a `ref_col` is the VCF preset, whose end column is
    /// the `INFO` column: the row ends at the `END` key read to its first
    /// non-digit when it lies at or past the start, which htslib through 1.21,
    /// Ensembl VEP's 1.9 among them, takes over the reference allele's span
    /// wherever it lies inside it (1.22 ends the row at the larger of the two);
    /// without one the row spans the reference allele from the start, one base
    /// where that cell is empty or missing (the span the start column alone
    /// gives a row, which `tbx_parse1` replaces only for a cell of one or more
    /// bases). The `INFO` column is never read as a number: its one digit-led
    /// key, `1000G`, is no end. A configuration without a `ref_col` is the
    /// generic preset, whose end cell is read to its leading integer as the
    /// start cell is (a BED end, exclusive and zero-based, reads as the
    /// inclusive one-based end as it stands); a cell no integer leads ends the
    /// row at its start.
    pub fn record_end(&self, start: u64, columns: &[String]) -> u64 {
        self.record_end_of(start, columns)
    }

    /// [`Self::record_end`] over cells of any text type, so a reader that
    /// keeps no row ends one from cells borrowed out of its line.
    pub fn record_end_of<S: AsRef<str>>(&self, start: u64, cells: &[S]) -> u64 {
        self.end_or(start, cells, |_| Some(start)).unwrap_or(start)
    }

    /// [`Self::record_end`], with `unread` deciding the end of a generic end
    /// cell no integer leads, which it is handed only once that read fails:
    /// `None` from it is no row at all.
    fn end_or<S: AsRef<str>>(
        &self,
        start: u64,
        columns: &[S],
        unread: impl FnOnce(&str) -> Option<u64>,
    ) -> Option<u64> {
        // An end column that is the start column names a single base, as
        // `tabix -b N -e N` does, so it takes the start's zero-based shift too.
        let Some(end_col) = self.end_col.filter(|&end_col| end_col != self.start_col) else {
            return Some(start);
        };
        let Some(ref_col) = self.ref_col else {
            return match columns.get(end_col).map(AsRef::as_ref) {
                Some(cell) => leading_integer(cell).or_else(|| unread(cell)),
                None => Some(start),
            };
        };
        let ref_allele_end = || {
            let ref_len = columns.get(ref_col).map_or(0, |cell| cell.as_ref().len()) as u64;
            start.saturating_add(ref_len.max(1) - 1)
        };
        // A row short of its INFO cell spans its reference allele, as the VCF
        // preset of tabix ends every row, so a query the allele reaches into
        // returns it.
        Some(
            columns
                .get(end_col)
                .and_then(|cell| info_end(cell.as_ref()))
                .filter(|&end| end >= start)
                .unwrap_or_else(ref_allele_end),
        )
    }
}

/// The leading decimal integer of `cell`: blanks of C's `isspace` skipped, an
/// optional `+`, then at least one digit, the rest of the cell left unread, as
/// C's `strtoll` reads a cell with neither a `-` nor a hex or octal prefix.
/// `None` when no digit leads the cell or the digits overflow.
#[doc(hidden)]
pub fn leading_integer(cell: &str) -> Option<u64> {
    let cell = cell.trim_start_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r']);
    let digits = cell.strip_prefix('+').unwrap_or(cell);
    let len = digits.bytes().take_while(u8::is_ascii_digit).count();
    digits[..len].parse().ok()
}

/// Whether a minus sign leads a number in `cell`: the blanks
/// [`leading_integer`] skips, a `-`, then a digit, which C's `strtol` reads
/// as a number below zero however many digits follow.
pub(crate) fn minus_leads_a_number(cell: &str) -> bool {
    cell.trim_start_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r'])
        .strip_prefix('-')
        .is_some_and(|digits| digits.starts_with(|digit: char| digit.is_ascii_digit()))
}

/// Batch-oriented tabix query engine.
///
/// Opens a tabix-indexed file and provides four query modes:
/// - `query_batch`: accepts multiple regions, merges them, and issues minimal
///   seeks. This is the primary hot path for plugin execution.
/// - `query_windows`: sweeps a batch's sorted windows in one pass, keeping the
///   records inside them alone.
/// - `query_windows_parallel`: the same windows, their span's blocks inflated
///   side by side on the current rayon pool.
/// - `query`: single-region convenience method.
pub struct TabixAnnotator {
    config: TabixConfig,
    /// Column header names parsed from the file header.
    header: Option<Arc<Vec<String>>>,
    /// Whether `header` is the file's first line, taken for a bare header
    /// because no integer leads its start cell ([`parse_header`]).
    bare_header: bool,
    /// Loaded once and reused across queries.
    index: Index,
    /// Where `index` was read from, for error messages.
    index_path: PathBuf,
    /// Whether a record line keeps the whitespace before its line end
    /// ([`Self::keeping_trailing_whitespace`]); by default it is trimmed.
    keep_trailing_whitespace: bool,
    /// The cells of a record built from the line's bytes, each with its
    /// decoder ([`Self::decoding_column_from_bytes`]); by default every cell
    /// is text.
    byte_columns: Vec<(usize, CellDecoder)>,
    /// Whether a row is a record exactly when the index build of `tabix`
    /// parses its start and end cells
    /// ([`Self::reading_rows_as_tabix_indexes_them`]).
    rows_as_tabix_indexes_them: bool,
    /// The run's assembly, which decides the spellings a mitochondrial name
    /// is tried under ([`plugin_table_candidates`]); none by default.
    assembly: Option<Assembly>,
    /// What decides, on GRCh37, whether the rows of a mitochondrial contig
    /// are read for a record of the other pair of names.
    frame: MitochondrialFrame,
}

/// A cell's text read from its bytes, as a format that spells bytes in
/// escapes reads a cell: every other cell is the bytes read as UTF-8 with
/// U+FFFD for a sequence that is not.
pub(crate) type CellDecoder = fn(&[u8]) -> Cow<'_, str>;

/// A loaded `.tbi` or `.csi` index; the query path is the same for both.
pub(crate) type Index = Box<dyn BinningIndex + Send + Sync>;

/// Where a line's start falls against the windows of one cluster of
/// [`TabixAnnotator::query_windows`].
enum LinePlacement {
    InWindow,
    BetweenWindows,
    PastCluster,
}

impl TabixAnnotator {
    /// Open a tabix-indexed file with the given configuration.
    pub fn open(config: TabixConfig) -> Result<Self, PluginError> {
        let file_path = config.file_path.clone();

        let index_path = find_index_path(&file_path)?;
        let index = load_index(&index_path)?;
        if index.header().is_none() {
            return Err(PluginError::Init(format!(
                "index {} carries no tabix header, so its records cannot be \
                 mapped to chromosomes; index the file with tabix",
                index_path.display()
            )));
        }

        let header_line = header_line(&file_path, config.start_col)?;
        let bare_header = header_line
            .as_ref()
            .is_some_and(|line| !line.starts_with('#'));
        let header = header_line.map(|line| Arc::new(header_names(&line)));
        std::fs::File::open(&file_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to open tabix data file {}: {e}",
                file_path.display()
            ))
        })?;

        debug!(
            file = %file_path.display(),
            index = %index_path.display(),
            header_cols = header.as_ref().map(|h| h.len()).unwrap_or(0),
            "opened tabix file"
        );

        let frame = MitochondrialFrame::comparing_column(config.ref_col);
        Ok(Self {
            config,
            header,
            bare_header,
            index,
            index_path,
            keep_trailing_whitespace: false,
            byte_columns: Vec::new(),
            rows_as_tabix_indexes_them: false,
            assembly: None,
            frame,
        })
    }

    /// Open a plugin's data table as Ensembl VEP's tabix plugins open theirs:
    /// [`check_plugin_table`] first, then the index htslib reads
    /// ([`find_index_path`]). A `--custom` file and the cache's variation
    /// files have no such check and open through [`Self::open`].
    pub(crate) fn open_plugin_table(config: TabixConfig) -> Result<Self, PluginError> {
        check_plugin_table(&config.file_path)?;
        Self::open(config)
    }

    /// `self` naming a query's chromosome to the table as a plugin's table is
    /// named it under `assembly` (`plugin_table_candidates`).
    #[doc(hidden)]
    pub fn on_assembly(mut self, assembly: Option<Assembly>) -> Self {
        self.assembly = assembly;
        self
    }

    /// The frame check's inputs, for the caller that holds the run's reference
    /// sequence or names the table's reference-allele column
    /// ([`MitochondrialFrame`]).
    #[doc(hidden)]
    pub fn frame_mut(&mut self) -> &mut MitochondrialFrame {
        &mut self.frame
    }

    /// The index's reference-sequence id for `chr` and the file's own spelling
    /// of the name, which its records carry; `None` for a chromosome the
    /// table is not read under.
    ///
    /// A mitochondrial name on GRCh37 is tried under each spelling of
    /// [`grch37_mitochondrial_candidates`] in turn. A contig of the name's own
    /// pair is read as it stands; a contig of the other pair is read only
    /// when its rows state the run's reference bases
    /// ([`MitochondrialFrame::admits`]). Every other name, and every name on
    /// another assembly, is resolved by [`plugin_table_candidates`] alone.
    fn resolve(&self, chr: &str) -> Result<Option<(usize, &bstr::BString)>, PluginError> {
        let Some(idx_header) = self.index.header() else {
            return Ok(None);
        };
        let names = idx_header.reference_sequence_names();
        if self.assembly != Some(Assembly::Grch37) || !is_mitochondrial_name(chr.trim()) {
            return Ok(resolve_reference_sequence(names, chr, self.assembly));
        }
        for (candidate, crosses) in grch37_mitochondrial_candidates(chr, self.frame.first_name) {
            let Some((id, file_chr)) = names.get_full(&bstr::BString::from(candidate.as_bytes()))
            else {
                continue;
            };
            let read = !crosses
                || self.frame.admits(
                    &self.config.file_path,
                    &candidate,
                    &|start, end, limit| self.region_rows(id, file_chr, start, end, limit),
                )?;
            if read {
                return Ok(Some((id, file_chr)));
            }
        }
        Ok(None)
    }

    /// `self` reading each record line as Ensembl VEP's tabix iterator hands it
    /// to a plugin: the line end removed and the whitespace before it kept, so a
    /// trailing cell that is empty or whitespace stays a cell of the record for
    /// a plugin that counts a row's cells as Ensembl VEP's `split /\t/` does.
    #[doc(hidden)]
    pub fn keeping_trailing_whitespace(mut self) -> Self {
        self.keep_trailing_whitespace = true;
        self
    }

    /// `self` building the cell at `column` of each record from the line's
    /// bytes through `decode` rather than from its text, for a format whose
    /// cell spells bytes in escapes: an escape there can complete a character
    /// whose first byte stands raw in the cell, as it does when a byte-wise
    /// parser decodes it, where the cell's text already holds U+FFFD for that
    /// byte. Each call adds a column, a later call on the same column replacing
    /// its decoder; the other cells, and every record of a reader without this,
    /// are text.
    pub(crate) fn decoding_column_from_bytes(mut self, column: usize, decode: CellDecoder) -> Self {
        self.byte_columns.retain(|(at, _)| *at != column);
        self.byte_columns.push((column, decode));
        self
    }

    /// `self` taking a row for a record exactly when the index build of
    /// `tabix` parses its start and end cells, for a format whose reader takes
    /// each cell's own number from the record. A row whose start cell a minus
    /// sign leads (`-5`) starts where a cell of `0` starts one, the first base
    /// `tabix` indexes it at (htslib 1.9 `tbx.c:99`), so every range from
    /// there to its end returns it; without this the row is dropped with every
    /// other row no integer leads. A row whose generic end cell no number
    /// leads, or a minus sign leads, is no record: the build fails the line
    /// (`tbx.c:105-106`, `:142`) or, for an end that reads zero, files it
    /// where no range returns it; without this it ends at its start
    /// ([`TabixConfig::record_end`]). End digits past the range end the row at
    /// the last position.
    ///
    /// The first line of a file with no `#` header is a row for this reader
    /// whenever its start cell starts one (`Self::cell_start`), so a line
    /// [`parse_header`] took for a bare header because a minus sign leads
    /// that cell is no header here. A batch lists a row that starts below the
    /// first base among the rows written at the first base, in the file's
    /// order ([`Self::query_batch`]).
    pub(crate) fn reading_rows_as_tabix_indexes_them(mut self) -> Self {
        self.rows_as_tabix_indexes_them = true;
        let first_line_is_a_row = self.bare_header
            && self
                .header
                .as_ref()
                .and_then(|cells| cells.get(self.config.start_col))
                .is_some_and(|cell| self.cell_start(cell).is_some());
        if first_line_is_a_row {
            self.header = None;
            self.bare_header = false;
        }
        self
    }

    /// `text` without its line end, and without the whitespace before it unless
    /// [`Self::keeping_trailing_whitespace`] was set.
    fn line_text<'a>(&self, text: &'a str) -> &'a str {
        if self.keep_trailing_whitespace {
            text.trim_end_matches(['\r', '\n'])
        } else {
            text.trim_end()
        }
    }

    /// A record line as [`Self::line_text`] cuts it: its text, read with U+FFFD
    /// for a sequence that is not UTF-8, and the bytes of the line it was cut
    /// from. The characters cut from the end are whitespace, which the lossy
    /// read copies byte for byte, so the text's cut is the bytes' cut.
    fn line_parts<'a>(&self, line: &'a [u8]) -> (Cow<'a, str>, &'a [u8]) {
        let text = String::from_utf8_lossy(line);
        let kept = self.line_text(&text).len();
        let cut = text.len() - kept;
        let bytes = &line[..line.len() - cut];
        debug_assert_eq!(String::from_utf8_lossy(bytes), text[..kept]);
        let text = match text {
            Cow::Borrowed(text) => Cow::Borrowed(&text[..kept]),
            Cow::Owned(mut text) => {
                text.truncate(kept);
                Cow::Owned(text)
            }
        };
        (text, bytes)
    }

    /// Batch query: accepts a slice of (chr, start, end) regions, merges
    /// overlapping regions, issues minimal tabix seeks, and returns all
    /// records indexed by position.
    pub fn query_batch(
        &self,
        regions: &[(&str, u64, u64)],
    ) -> Result<BatchQueryResult, PluginError> {
        if regions.is_empty() {
            return Ok(BatchQueryResult::empty());
        }

        let mut by_chr: HashMap<&str, Vec<(u64, u64)>> = HashMap::new();
        for &(chr, start, end) in regions {
            let (start, end) = normalize_range_bounds(start, end);
            by_chr.entry(chr).or_default().push((start, end));
        }

        let mut all_records: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>> = HashMap::new();

        for (chr, mut intervals) in by_chr {
            intervals.sort_by_key(|&(s, _)| s);
            let merged = merge_intervals(&intervals);

            let mut previous_end: Option<u64> = None;
            for &(region_start, region_end) in &merged {
                let records = self.query_tabix_region(chr, region_start, region_end)?;
                for record in records {
                    // A record starting at or before the previous region's end
                    // and reaching into this one overlaps both, and the previous
                    // region's read already collected it.
                    if previous_end.is_some_and(|end| record.start <= end) {
                        continue;
                    }
                    // `tabix` files a row that starts below the first base at
                    // the first base (htslib 1.9 `tbx.c:99`) and hands it out
                    // among the rows written there, in the file's order.
                    let filed_at = if self.rows_as_tabix_indexes_them {
                        record.start.max(1)
                    } else {
                        record.start
                    };
                    all_records
                        .entry(chr.to_string())
                        .or_default()
                        .entry(filed_at)
                        .or_default()
                        .push(record);
                }
                previous_end = Some(region_end);
            }
        }

        Ok(BatchQueryResult::from_records(all_records))
    }

    /// Query a single region. Convenience wrapper. The file is read as tabix
    /// indexed it, sorted by start within a chromosome, so the scan of a chunk
    /// ends at the first record past the interval.
    pub fn query(&self, chr: &str, start: u64, end: u64) -> Result<Vec<TabixRecord>, PluginError> {
        let (start, end) = normalize_range_bounds(start, end);
        self.query_tabix_region(chr, start, end)
    }

    /// The records whose start lies in any of `windows`, sorted and disjoint
    /// `(start, end)` pairs of one chromosome, read in one pass over the file
    /// for an index built on the start column alone (`end_col` none or the
    /// start column; any other index answers window by window).
    ///
    /// Windows within `GAP_BRIDGE` bases of each other share one index
    /// lookup and one scan of its chunks, and the reader stays where the
    /// previous scan ended when the next chunk holds that position, so a block
    /// is inflated once however many windows fall in it. A line between two
    /// windows is skipped by its start alone and nothing between them is kept,
    /// so a sparse batch holds its own records only.
    pub fn query_windows(
        &self,
        chr: &str,
        windows: &[(u64, u64)],
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let mut records = Vec::new();
        let start_is_span = self
            .config
            .end_col
            .is_none_or(|end_col| end_col == self.config.start_col);
        if !start_is_span {
            for &(start, end) in windows {
                records.extend(self.query(chr, start, end)?);
            }
            return Ok(records);
        }
        let Some((ref_seq_id, file_chr)) = self.resolve(chr)? else {
            return Ok(records);
        };
        let file = std::fs::File::open(&self.config.file_path).map_err(|e| self.read_error(e))?;
        let mut bgzf_reader = bgzf::io::Reader::new(file);
        let last_held = self.last_held();
        let mut line = Vec::new();
        let mut last_chunk_end: Option<bgzf::VirtualPosition> = None;
        // The start and bytes of the line read past the previous cluster's last
        // window, carried to the next cluster, which may hold it.
        let mut carried: Option<(u64, Vec<u8>)> = None;
        let mut first = 0;
        while first < windows.len() {
            let past = cluster_end(windows, first);
            let cluster = &windows[first..past];
            first = past;
            let (start, end) = (cluster[0].0, cluster[cluster.len() - 1].1);
            let start_pos = Position::try_from(start.max(1) as usize)
                .map_err(|e| PluginError::Run(format!("invalid start position: {e}")))?;
            let end_pos = Position::try_from(end.max(1) as usize)
                .map_err(|e| PluginError::Run(format!("invalid end position: {e}")))?;
            let chunks = match self.index.query(ref_seq_id, (start_pos..=end_pos).into()) {
                Ok(chunks) => chunks,
                // The index refuses a span that ends past its last addressable
                // position, though the windows short of it may hold records.
                Err(_) => {
                    for &(start, end) in cluster {
                        records.extend(self.query(chr, start, end)?);
                    }
                    continue;
                }
            };
            // The window the next line can fall in; the lines come sorted.
            let mut window = 0;
            let mut place = |line_start: u64| -> LinePlacement {
                while window < cluster.len() && cluster[window].1 < line_start {
                    window += 1;
                }
                if window == cluster.len() {
                    LinePlacement::PastCluster
                } else if line_start >= cluster[window].0 {
                    LinePlacement::InWindow
                } else {
                    LinePlacement::BetweenWindows
                }
            };
            'cluster: for chunk in chunks {
                let position = bgzf_reader.virtual_position();
                let continues = last_chunk_end.is_some_and(|last_end| {
                    chunk.start() <= position && position <= last_end && last_end <= chunk.end()
                });
                if continues {
                    if let Some((line_start, bytes)) = carried.take() {
                        match place(line_start) {
                            LinePlacement::PastCluster => {
                                carried = Some((line_start, bytes));
                                break 'cluster;
                            }
                            LinePlacement::BetweenWindows => {}
                            LinePlacement::InWindow => {
                                let text = String::from_utf8_lossy(&bytes);
                                if let Some(record) = self.record_on(&text, &bytes, file_chr) {
                                    records.push(record);
                                }
                            }
                        }
                    }
                } else {
                    carried = None;
                    self.seek_to_chunk_start(&mut bgzf_reader, chunk.start())?;
                }
                last_chunk_end = Some(chunk.end());
                while bgzf_reader.virtual_position() < chunk.end() {
                    line.clear();
                    let bytes_read = self.read_line(&mut bgzf_reader, &mut line)?;
                    if bytes_read == 0 {
                        break;
                    }
                    let (text, bytes) = self.line_parts(&line);
                    if text.is_empty() || text.starts_with('#') {
                        continue;
                    }
                    let Some(line_start) = self
                        .line_start(&text)
                        .filter(|&line_start| line_start <= last_held)
                    else {
                        continue;
                    };
                    match place(line_start) {
                        LinePlacement::PastCluster => {
                            carried = Some((line_start, bytes.to_vec()));
                            break 'cluster;
                        }
                        LinePlacement::BetweenWindows => {}
                        LinePlacement::InWindow => {
                            if let Some(record) = self.record_on(&text, bytes, file_chr) {
                                records.push(record);
                            }
                        }
                    }
                }
            }
        }
        Ok(records)
    }

    /// The records [`query_windows`](Self::query_windows) returns for `windows`
    /// (sorted and disjoint `(start, end)` pairs of one chromosome), read block
    /// by block on the current rayon pool for an index built on the start
    /// column alone (any other index answers through `query_windows`).
    ///
    /// The windows form the clusters `query_windows` reads, one index lookup
    /// per cluster; the chunks of every cluster are merged where they share a
    /// block, so an index tile several clusters fall in is read once, as the
    /// line-by-line read keeps its place in it, and the blocks between two
    /// clusters apart in the file are not read at all. Each merged span's
    /// blocks are read once, in slices of a bounded byte count, inflated side by
    /// side, and every line is kept when its start lies in a window. A line
    /// that crosses a block boundary belongs to the block it starts in, and one
    /// more block past each span's end is read so the span's last line is
    /// whole. The records come back in file order.
    ///
    /// `query_windows` reads a cluster of windows on one thread; a table with a
    /// row at every base hands a batch's span to that thread as tens of blocks,
    /// and the threads the batch's annotation would run on wait for it. A
    /// cluster the index refuses, a line left open by the blocks read, or a
    /// block that does not inflate sends the windows to `query_windows`.
    pub fn query_windows_parallel(
        &self,
        chr: &str,
        windows: &[(u64, u64)],
    ) -> Result<Vec<TabixRecord>, PluginError> {
        self.query_windows_in_slices(chr, windows, SLICE_BYTES)
    }

    /// [`query_windows_parallel`](Self::query_windows_parallel) with the slice
    /// size it reads a span in, so that a test can drive a span through many
    /// slices.
    fn query_windows_in_slices(
        &self,
        chr: &str,
        windows: &[(u64, u64)],
        slice_bytes: usize,
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let start_is_span = self
            .config
            .end_col
            .is_none_or(|end_col| end_col == self.config.start_col);
        if !start_is_span {
            return self.query_windows(chr, windows);
        }
        let Some((ref_seq_id, file_chr)) = self.resolve(chr)? else {
            return Ok(Vec::new());
        };
        let mut chunks = Vec::new();
        let mut first = 0;
        while first < windows.len() {
            let past = cluster_end(windows, first);
            let (start, end) = (windows[first].0, windows[past - 1].1);
            first = past;
            let start_pos = Position::try_from(start.max(1) as usize)
                .map_err(|e| PluginError::Run(format!("invalid start position: {e}")))?;
            let end_pos = Position::try_from(end.max(1) as usize)
                .map_err(|e| PluginError::Run(format!("invalid end position: {e}")))?;
            match self.index.query(ref_seq_id, (start_pos..=end_pos).into()) {
                Ok(found) => chunks.extend(found),
                // The index refuses a span that ends past its last addressable
                // position, though the windows short of it may hold records.
                Err(_) => return self.query_windows(chr, windows),
            }
        }
        let mut records = Vec::new();
        for (span_start, span_end) in merge_chunks(&chunks) {
            match self.span_records_by_block(
                span_start,
                span_end,
                slice_bytes,
                windows,
                file_chr,
            )? {
                Some(found) => records.extend(found),
                None => return self.query_windows(chr, windows),
            }
        }
        Ok(records)
    }

    /// The records of one merged span of
    /// [`query_windows_parallel`](Self::query_windows_parallel) whose start
    /// lies in a window, or `None` when the span is to be read line by line
    /// instead: the blocks read leave a line open, or a block does not inflate
    /// (the line-by-line read reaches a damaged block only when a window needs
    /// it, and raises its own error there).
    ///
    /// The span's blocks, the one past its end included, are read in order in
    /// slices of about `slice_bytes` compressed bytes, each inflated and
    /// scanned side by side before the next is read, so a span of any length
    /// holds one slice at a time. The last block of a slice is read for the
    /// line that ends in it alone and opens the next slice, which begins after
    /// that line when the block before it ends inside one.
    fn span_records_by_block(
        &self,
        span_start: bgzf::VirtualPosition,
        span_end: bgzf::VirtualPosition,
        slice_bytes: usize,
        windows: &[(u64, u64)],
        file_chr: &bstr::BString,
    ) -> Result<Option<Vec<TabixRecord>>, PluginError> {
        let not_a_block = |offset: u64| {
            PluginError::Run(format!(
                "index {} does not describe {}: no block at offset {offset}",
                self.index_path.display(),
                self.config.file_path.display()
            ))
        };
        let start_block = span_start.compressed();
        let end_block = span_end.compressed();
        let file = std::fs::File::open(&self.config.file_path).map_err(|e| self.read_error(e))?;
        let file_len = file.metadata().map_err(|e| self.read_error(e))?.len();
        if end_block > file_len || start_block > end_block {
            return Err(not_a_block(end_block));
        }
        let mut blocks = BlockReader::new(file, start_block).map_err(|e| self.read_error(e))?;
        // A span whose last line ends with a block ends at the byte past it, so
        // the end block then holds none of the span's bytes and may be missing
        // (a file with no EOF block) or cut short.
        let end_block_holds_span_bytes = end_block == start_block || span_end.uncompressed() != 0;

        let mut records = Vec::new();
        let mut slice = Slice {
            blocks: Vec::new(),
            start_block,
            first_offset: usize::from(span_start.uncompressed()),
            open_line_before: false,
        };
        let mut slice_len = 0;
        let (last_is_tail, file_ended) = loop {
            let block_offset = blocks.offset();
            let needed = block_offset < end_block
                || (block_offset == end_block && end_block_holds_span_bytes);
            match blocks.next_block() {
                Ok(Some(block)) => {
                    slice_len += block.len();
                    slice.blocks.push(block);
                    if block_offset > end_block {
                        break (true, false);
                    }
                    if slice_len >= slice_bytes {
                        let Some(found) =
                            self.scan_slice(&mut slice, true, false, windows, file_chr)?
                        else {
                            return Ok(None);
                        };
                        records.extend(found);
                        slice_len = slice.blocks.iter().map(Vec::len).sum();
                    }
                }
                Ok(None) if needed => return Err(not_a_block(block_offset)),
                Ok(None) => break (false, true),
                Err(e) if needed => {
                    return Err(match e.kind() {
                        std::io::ErrorKind::InvalidData => not_a_block(block_offset),
                        _ => self.read_error(e),
                    })
                }
                // The block past the span's end is cut short or is not a
                // block: the span's own last line is judged without it.
                Err(_) => break (false, false),
            }
        };
        let Some(found) =
            self.scan_slice(&mut slice, last_is_tail, file_ended, windows, file_chr)?
        else {
            return Ok(None);
        };
        records.extend(found);
        Ok(Some(records))
    }

    /// The records of one slice whose start lies in a window, or `None` when a
    /// line is left open or a block does not inflate. The blocks are inflated
    /// side by side, cut into the lines that start in each
    /// ([`block_line_runs`]), and each block's lines are placed and parsed side
    /// by side. A slice that is not the span's last keeps its last block, read
    /// for the line that ends in it alone, as the first block of the next slice.
    fn scan_slice(
        &self,
        slice: &mut Slice,
        last_is_tail: bool,
        file_ended: bool,
        windows: &[(u64, u64)],
        file_chr: &bstr::BString,
    ) -> Result<Option<Vec<TabixRecord>>, PluginError> {
        let Ok(inflated) = slice
            .blocks
            .par_iter()
            .map(|block| self.inflate_block(block))
            .collect::<Result<Vec<Vec<u8>>, PluginError>>()
        else {
            return Ok(None);
        };
        if !slice.open_line_before
            && inflated
                .first()
                .is_some_and(|text| slice.first_offset > text.len())
        {
            return Err(PluginError::Run(format!(
                "index {} does not describe {}: no block at offset {} holds byte {}",
                self.index_path.display(),
                self.config.file_path.display(),
                slice.start_block,
                slice.first_offset
            )));
        }
        let (runs, incomplete, open_line_before_tail) = block_line_runs(
            &inflated,
            slice.first_offset,
            slice.open_line_before,
            last_is_tail,
            file_ended,
        );
        if incomplete {
            return Ok(None);
        }
        let last_held = self.last_held();
        let per_run: Vec<Vec<TabixRecord>> = runs
            .par_iter()
            .map(|run| {
                let mut found = Vec::new();
                // The lines of a run come sorted by start, so one cursor over
                // the sorted windows, set by a binary search for the first line,
                // places every line.
                let mut window = None;
                let lines = run
                    .body
                    .split(|&b| b == b'\n')
                    .chain(run.straddle.as_deref());
                for line in lines {
                    // Most lines lie between the windows, so a line is placed
                    // from its start column alone before its text is decoded
                    // and its record built.
                    if line.is_empty() || line.starts_with(b"#") {
                        continue;
                    }
                    let Some(line_start) = self
                        .line_start_bytes(line)
                        .filter(|&line_start| line_start <= last_held)
                    else {
                        continue;
                    };
                    let at = window.get_or_insert_with(|| {
                        windows.partition_point(|&(_, end)| end < line_start)
                    });
                    while *at < windows.len() && windows[*at].1 < line_start {
                        *at += 1;
                    }
                    if *at < windows.len() && windows[*at].0 <= line_start {
                        let (text, bytes) = self.line_parts(line);
                        if let Some(record) = self.record_on(&text, bytes, file_chr) {
                            found.push(record);
                        }
                    }
                }
                found
            })
            .collect();
        if last_is_tail && slice.blocks.len() > 1 {
            // The tail opens the next slice; the next slice's first line starts
            // after the line the blocks before the tail end inside, if any.
            slice.open_line_before = open_line_before_tail;
            slice.first_offset = 0;
            slice.blocks.drain(..slice.blocks.len() - 1);
        }
        Ok(Some(per_run.into_iter().flatten().collect()))
    }

    /// The text of one bgzf block, its buffer sized from the block's `ISIZE`
    /// trailer (RFC 1952 § 2.3.1: the last four bytes, the inflated size),
    /// capped at the 65,536 bytes a bgzf block may inflate to so that a
    /// corrupt trailer reserves no more.
    fn inflate_block(&self, block: &[u8]) -> Result<Vec<u8>, PluginError> {
        let inflated_size = block
            .last_chunk::<4>()
            .map_or(0, |trailer| u32::from_le_bytes(*trailer));
        let mut text = Vec::with_capacity((inflated_size as usize).min(1 << 16));
        bgzf::io::Reader::new(std::io::Cursor::new(block))
            .read_to_end(&mut text)
            .map_err(|e| self.read_error(e))?;
        Ok(text)
    }

    /// Query a region using the tabix index and bgzf reader.
    ///
    /// Uses the binning index to find relevant chunks in the bgzf file, seeks
    /// to each chunk, and reads the records whose span overlaps the query
    /// region, as `tabix <file> <region>` does.
    fn query_tabix_region(
        &self,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let (start, end) = normalize_range_bounds(start, end);
        let Some((ref_seq_id, file_chr)) = self.resolve(chr)? else {
            return Ok(Vec::new()); // Chromosome not in index
        };
        self.region_rows(ref_seq_id, file_chr, start, end, usize::MAX)
    }

    /// The records of the index's sequence `ref_seq_id`, which the file names
    /// `file_chr`, whose span overlaps `start..=end`, in the file's order, the
    /// read ending once `limit` of them are held.
    fn region_rows(
        &self,
        ref_seq_id: usize,
        file_chr: &bstr::BString,
        start: u64,
        end: u64,
        limit: usize,
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let mut records = Vec::new();
        if limit == 0 {
            return Ok(records);
        }

        // `tabix` reads an interval ending past the last position the index
        // holds to it and finds nothing in one starting past it (htslib
        // `reg2bins` clamps the end). Refusing the interval would drop every
        // record of a merged region that one wide span reaches past the limit.
        let last_held = self.last_held();
        if start > last_held {
            return Ok(records);
        }
        let end = end.min(last_held);
        // The index lookup refuses a bound past `last_held - 1`, which lies in
        // the smallest bin holding `last_held`, so that position names the
        // same bins and stands in for it; the record filter keeps the real end.
        let last_lookup = last_held - 1;
        let start_pos = Position::try_from(start.min(last_lookup).max(1) as usize)
            .map_err(|e| PluginError::Run(format!("invalid start position: {e}")))?;
        let end_pos = Position::try_from(end.min(last_lookup).max(1) as usize)
            .map_err(|e| PluginError::Run(format!("invalid end position: {e}")))?;
        let interval = start_pos..=end_pos;

        let chunks = match self.index.query(ref_seq_id, interval.into()) {
            Ok(chunks) => chunks,
            Err(_) => return Ok(records), // Region not found in index
        };

        if chunks.is_empty() {
            return Ok(records);
        }

        // Each query opens its own handle: a duplicated handle shares the file
        // offset, so two threads seeking the same file would read each other's
        // blocks.
        let file = std::fs::File::open(&self.config.file_path).map_err(|e| self.read_error(e))?;
        let mut bgzf_reader = bgzf::io::Reader::new(file);

        // No `BufReader` wraps the reader: a bulk fill runs past a chunk's end
        // into bytes the next chunk's seek reads again, so records are returned
        // twice and a line cut at the fill's end is glued to the next chunk's
        // first line.

        // A record is found by its span, so a start-only index (`end_col` the
        // start column) can drop a line by its start alone.
        let start_is_span = self
            .config
            .end_col
            .is_none_or(|end_col| end_col == self.config.start_col);
        let mut line = Vec::new();
        'chunks: for chunk in chunks {
            self.seek_to_chunk_start(&mut bgzf_reader, chunk.start())?;

            while bgzf_reader.virtual_position() < chunk.end() {
                line.clear();
                let bytes_read = self.read_line(&mut bgzf_reader, &mut line)?;
                if bytes_read == 0 {
                    break;
                }

                let (text, bytes) = self.line_parts(&line);
                if text.is_empty() || text.starts_with('#') {
                    continue;
                }

                // A bin of a dense file holds thousands of records, so the
                // start is read before the record is built: the file is sorted
                // by start, so the first line past the interval ends the read,
                // as htslib's iterator ends there (the chunks come in file
                // order, so every later chunk holds later starts), and on a
                // start-only index the lines before it need no record either.
                let Some(line_start) = self
                    .line_start(&text)
                    .filter(|&line_start| line_start <= last_held)
                else {
                    continue;
                };
                if line_start > end {
                    break 'chunks;
                }
                if start_is_span && line_start < start {
                    continue;
                }

                let Some(record) = self.parse_line(&text, bytes) else {
                    continue;
                };
                // The index is bin-granular, so a chunk legitimately carries
                // records outside the requested interval; its bounds are also
                // all that keeps another chromosome's records out of the read,
                // so the chromosome is checked as well.
                if record.chr.as_bytes() == file_chr.as_slice()
                    && record.start <= end
                    && record.end >= start
                {
                    records.push(record);
                    if records.len() >= limit {
                        break 'chunks;
                    }
                }
            }
        }

        Ok(records)
    }

    /// Position `reader` at `chunk_start`, refusing an offset the block does
    /// not hold.
    ///
    /// The chunk's block is entered at offset 0 so that its length is known
    /// before the chunk's offset into it is applied: an index built for another
    /// file names offsets no block of this one holds. A seek past the end of
    /// the file reads no block and leaves the one read before in place, so the
    /// block reached must report the offset asked for; a chunk starts at a
    /// record, so a block with no bytes is refused as well.
    ///
    /// One chunk start names a block without bytes and is read: `tabix` gives
    /// the record that opens the block after an empty one (the end marker
    /// inside a file joined from two compressed streams) the empty block's
    /// offset, and reads on from the next block that holds bytes (htslib
    /// `bgzf_read_block` skips empty blocks). The seek skips them likewise, so
    /// the block reached is taken when only empty blocks lie between the offset
    /// asked for and it.
    fn seek_to_chunk_start(
        &self,
        reader: &mut bgzf::io::Reader<std::fs::File>,
        chunk_start: bgzf::VirtualPosition,
    ) -> Result<(), PluginError> {
        let (block_offset, offset_in_block): (u64, u16) = chunk_start.into();
        let offset_in_block = usize::from(offset_in_block);
        let block_start = bgzf::VirtualPosition::try_from((block_offset, 0))
            .map_err(|e| PluginError::Run(format!("invalid index offset: {e}")))?;
        reader.seek(block_start).map_err(|e| self.read_error(e))?;
        let block_len = reader.fill_buf().map_err(|e| self.read_error(e))?.len();
        let reached = reader.virtual_position();
        let at_chunk_start = reached == block_start
            || (offset_in_block == 0
                && reached.uncompressed() == 0
                && self.only_empty_blocks_between(block_offset, reached.compressed()));
        if block_len == 0 || !at_chunk_start || offset_in_block > block_len {
            return Err(PluginError::Run(format!(
                "index {} does not describe {}: no block at offset {block_offset} holds \
                 byte {offset_in_block}",
                self.index_path.display(),
                self.config.file_path.display()
            )));
        }
        reader.consume(offset_in_block);
        Ok(())
    }

    /// Whether the file holds whole bgzf blocks from compressed offset `from`
    /// up to `to`, at least one, each inflating to no bytes (its `ISIZE`
    /// trailer, RFC 1952 § 2.3.1, is zero). Bytes that are no block, a block
    /// that holds text, and a file that ends before `to` are all `false`.
    fn only_empty_blocks_between(&self, from: u64, to: u64) -> bool {
        if from >= to {
            return false;
        }
        let Ok(file) = std::fs::File::open(&self.config.file_path) else {
            return false;
        };
        let Ok(mut blocks) = BlockReader::new(file, from) else {
            return false;
        };
        while blocks.offset() < to {
            match blocks.next_block() {
                Ok(Some(block)) if block.ends_with(&[0; 4]) => {}
                _ => return false,
            }
        }
        blocks.offset() == to
    }

    /// The last position the index holds, 2^(min_shift + 3 depth): 2^29 for a
    /// `.tbi`. A row that starts past it is not one `tabix` filed in the
    /// file's order, so every read passes over it as over a line that is no
    /// row: it ends no read and moves no window on, and the rows after it are
    /// read as they stand.
    #[doc(hidden)]
    pub fn last_held(&self) -> u64 {
        let address_bits = u32::from(self.index.min_shift()) + 3 * u32::from(self.index.depth());
        1u64.checked_shl(address_bits).unwrap_or(u64::MAX)
    }

    /// The error for an I/O failure while reading the data file.
    fn read_error(&self, e: std::io::Error) -> PluginError {
        PluginError::Run(format!(
            "failed to read {}: {e}",
            self.config.file_path.display()
        ))
    }

    /// Reads the next line of the data file into `line` with its length
    /// capped ([`read_line_bounded`]); a refused line is named by the virtual
    /// position it starts at.
    fn read_line(
        &self,
        reader: &mut bgzf::io::Reader<std::fs::File>,
        line: &mut Vec<u8>,
    ) -> Result<usize, PluginError> {
        let start = reader.virtual_position();
        read_line_bounded(
            reader,
            line,
            &self.config.file_path,
            PluginError::Run,
            || {
                let (block, offset): (u64, u16) = start.into();
                format!("the line at block offset {block} byte {offset}")
            },
        )
    }

    /// [`line_start`](Self::line_start) on the line's bytes: the start column
    /// is cut at the tabs, which no multi-byte character contains, so the two
    /// agree on every line.
    fn line_start_bytes(&self, line: &[u8]) -> Option<u64> {
        let field = line.split(|&b| b == b'\t').nth(self.config.start_col)?;
        // A line that crosses a block boundary keeps its line end, which the
        // leading-integer read leaves unread with anything else after the
        // digits, a byte that is not UTF-8 among them.
        self.cell_start(&String::from_utf8_lossy(field))
    }

    /// The 1-based start of a line, read from its start column alone
    /// ([`Self::cell_start`]); `None` for a line `parse_line` would reject on
    /// that column, among them a zero-based start with no one-based form.
    fn line_start(&self, line: &str) -> Option<u64> {
        let field = line.split('\t').nth(self.config.start_col)?;
        self.cell_start(field)
    }

    /// The 1-based start of a row from its start cell
    /// ([`TabixConfig::record_start`]) or, under
    /// [`Self::reading_rows_as_tabix_indexes_them`], the start a cell of `0`
    /// gives for a cell a minus sign leads, whatever its digits. That cell is
    /// read only once no integer leads it, so a row an integer leads costs
    /// what it costs without the option.
    fn cell_start(&self, cell: &str) -> Option<u64> {
        self.config.record_start(cell).or_else(|| {
            (self.rows_as_tabix_indexes_them && minus_leads_a_number(cell))
                .then_some(u64::from(self.config.zero_based))
        })
    }

    /// The end of a generic end cell no integer leads: the row's start, or
    /// under [`Self::reading_rows_as_tabix_indexes_them`] the last position
    /// for digits past the range and no row for any other cell.
    fn unread_end(&self, start: u64, cell: &str) -> Option<u64> {
        if !self.rows_as_tabix_indexes_them {
            return Some(start);
        }
        let cell = cell.trim_start_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r']);
        let digits = cell.strip_prefix('+').unwrap_or(cell);
        digits
            .starts_with(|digit: char| digit.is_ascii_digit())
            .then_some(u64::MAX)
    }

    /// Parse a tab-separated line into a TabixRecord; `bytes` is the line
    /// `line` was read from ([`Self::line_parts`]), which the cell of
    /// [`Self::decoding_column_from_bytes`] is built from.
    ///
    /// A bare (un-prefixed) header row is rejected naturally: no integer leads
    /// its start-column field, so `line_start` below returns `None`. `\r` is
    /// stripped per field because a line ending in CRLF (`dbscSNV.pm:191`
    /// strips the `\r`) would otherwise weld a carriage return onto the last
    /// column's value, corrupting the score that plugins read by name.
    fn parse_line(&self, line: &str, bytes: &[u8]) -> Option<TabixRecord> {
        if line.starts_with('#') {
            return None;
        }

        let cells = line.split('\t').map(|s| s.trim_end_matches('\r'));
        let mut columns: Vec<String> = if self.byte_columns.is_empty() {
            cells.map(str::to_string).collect()
        } else {
            // A tab is ASCII, so the bytes split into the text's cells.
            cells
                .zip(bytes.split(|&byte| byte == b'\t'))
                .enumerate()
                .map(|(index, (cell, raw))| {
                    match self.byte_columns.iter().find(|(at, _)| *at == index) {
                        Some((_, decode)) => decode(without_carriage_returns(raw)).into_owned(),
                        None => cell.to_string(),
                    }
                })
                .collect()
        };
        if columns.len() <= self.config.chr_col || columns.len() <= self.config.start_col {
            return None;
        }

        let chr = std::mem::take(&mut columns[self.config.chr_col]);
        let start = self.line_start(line)?;
        let end = self
            .config
            .end_or(start, &columns, |cell| self.unread_end(start, cell))?;

        Some(TabixRecord {
            chr,
            start,
            end,
            columns,
            header: self.header.clone(),
        })
    }

    /// The record of `line` when it lies on `file_chr`: a chunk's bounds alone
    /// do not keep another chromosome's lines out of a read.
    fn record_on(&self, line: &str, bytes: &[u8], file_chr: &bstr::BString) -> Option<TabixRecord> {
        self.parse_line(line, bytes)
            .filter(|record| record.chr.as_bytes() == file_chr.as_slice())
    }

    /// Get the column header names, if available.
    pub fn header(&self) -> Option<&[String]> {
        self.header.as_ref().map(|h| h.as_slice())
    }

    /// Whether the index names `chr`, under any spelling a query resolves
    /// ([`plugin_table_candidates`]): a query on a chromosome it does not
    /// name reads nothing, as `tabix` prints nothing for one, so a plugin that
    /// tells a chromosome the file lacks from a span with no records asks here.
    pub(crate) fn holds_chromosome(&self, chr: &str) -> bool {
        self.resolve(chr).is_ok_and(|found| found.is_some())
    }

    /// Get the config.
    pub(crate) fn config(&self) -> &TabixConfig {
        &self.config
    }
}

/// `cell` without the `\r` bytes at its end, the bytes' form of the per-field
/// `trim_end_matches('\r')` of [`TabixAnnotator::parse_line`].
fn without_carriage_returns(mut cell: &[u8]) -> &[u8] {
    while let Some(stripped) = cell.strip_suffix(b"\r") {
        cell = stripped;
    }
    cell
}

/// The `END` key of a `;`-separated `key=value` column (a VCF `INFO` column):
/// the value's leading integer as [`leading_integer`] reads it, blanks and a
/// `+` skipped as htslib's VCF preset skips them (`strtoll` from htslib 1.10),
/// `None` for a column without the key or a value no digit leads.
fn info_end(column: &str) -> Option<u64> {
    column
        .split(';')
        .find_map(|pair| pair.strip_prefix("END="))
        .and_then(leading_integer)
}

/// A gap of up to this many bases between two regions is read with them: one
/// longer read beats a second seek, and one index lookup serves them both.
const GAP_BRIDGE: u64 = 1000;

/// About how many compressed bytes of a span the block-parallel read holds at
/// once with their text: a few hundred blocks, enough to keep a pool busy, and
/// about half a megabase of a table with three rows at every base.
const SLICE_BYTES: usize = 16 << 20;

/// The index past the cluster of `windows` (sorted, disjoint) that starts at
/// `first`: the windows from `first` on, each within [`GAP_BRIDGE`] of the one
/// before it.
fn cluster_end(windows: &[(u64, u64)], first: usize) -> usize {
    let mut past = first + 1;
    while past < windows.len() && windows[past].0 <= windows[past - 1].1.saturating_add(GAP_BRIDGE)
    {
        past += 1;
    }
    past
}

/// Reads into `buf` until it is full or the reader ends: the bytes read.
fn read_up_to<R: Read>(reader: &mut R, buf: &mut [u8]) -> std::io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(filled)
}

/// The bytes of a bgzf block header: the gzip header (10) and its extra
/// field's length (2); the `BC` subfield that follows them carries the block
/// size (RFC 1952 § 2.3.1.1; the BGZF format of the SAM specification § 4.1).
const BGZF_HEADER_BYTES: usize = 12;

/// The size in bytes of the bgzf block whose first bytes are `bytes`, read
/// from the `BC` subfield of its header (the size less one); `None` when the
/// bytes are not a bgzf block header.
fn bgzf_block_size(bytes: &[u8]) -> Option<usize> {
    let header = bytes.get(..BGZF_HEADER_BYTES)?;
    if header[0] != 0x1f || header[1] != 0x8b || header[2] != 8 || header[3] & 4 == 0 {
        return None;
    }
    let xlen = usize::from(u16::from_le_bytes([header[10], header[11]]));
    let extra = bytes.get(BGZF_HEADER_BYTES..BGZF_HEADER_BYTES + xlen)?;
    let mut at = 0;
    while at + 4 <= extra.len() {
        let slen = usize::from(u16::from_le_bytes([extra[at + 2], extra[at + 3]]));
        if extra[at] == b'B' && extra[at + 1] == b'C' && slen == 2 {
            let bsize = u16::from_le_bytes([*extra.get(at + 4)?, *extra.get(at + 5)?]);
            return Some(usize::from(bsize) + 1);
        }
        at += 4 + slen;
    }
    None
}

/// The chunks' spans as sorted `(start, end)` virtual positions, two spans
/// merged when they share a block: a block is read and its lines kept once.
fn merge_chunks(chunks: &[Chunk]) -> Vec<(bgzf::VirtualPosition, bgzf::VirtualPosition)> {
    let mut spans: Vec<_> = chunks.iter().map(|c| (c.start(), c.end())).collect();
    spans.sort_unstable();
    let mut merged: Vec<(bgzf::VirtualPosition, bgzf::VirtualPosition)> = Vec::new();
    for (start, end) in spans {
        match merged.last_mut() {
            Some((_, last_end)) if start.compressed() <= last_end.compressed() => {
                if end > *last_end {
                    *last_end = end;
                }
            }
            _ => merged.push((start, end)),
        }
    }
    merged
}

/// Whole bgzf blocks read one after another from a compressed offset of a
/// file, each as its compressed bytes.
struct BlockReader {
    reader: BufReader<std::fs::File>,
    /// The compressed offset of the next block.
    offset: u64,
}

impl BlockReader {
    fn new(file: std::fs::File, offset: u64) -> std::io::Result<Self> {
        use std::io::{Seek, SeekFrom};
        let mut reader = BufReader::with_capacity(1 << 20, file);
        reader.seek(SeekFrom::Start(offset))?;
        Ok(Self { reader, offset })
    }

    fn offset(&self) -> u64 {
        self.offset
    }

    /// The next block, or `None` where the file ends cleanly at a block start.
    /// A file that ends inside a block is an `UnexpectedEof` error, and bytes
    /// that are not a bgzf block header an `InvalidData` one.
    fn next_block(&mut self) -> std::io::Result<Option<Vec<u8>>> {
        let mut block = vec![0u8; BGZF_HEADER_BYTES];
        let read = read_up_to(&mut self.reader, &mut block)?;
        if read == 0 {
            return Ok(None);
        }
        let short = || std::io::Error::from(std::io::ErrorKind::UnexpectedEof);
        if read < BGZF_HEADER_BYTES {
            return Err(short());
        }
        let extra_len = usize::from(u16::from_le_bytes([block[10], block[11]]));
        block.resize(BGZF_HEADER_BYTES + extra_len, 0);
        if read_up_to(&mut self.reader, &mut block[BGZF_HEADER_BYTES..])? < extra_len {
            return Err(short());
        }
        let header_len = block.len();
        let size = bgzf_block_size(&block)
            .filter(|&size| size >= header_len)
            .ok_or_else(|| std::io::Error::from(std::io::ErrorKind::InvalidData))?;
        block.resize(size, 0);
        if read_up_to(&mut self.reader, &mut block[header_len..])? < size - header_len {
            return Err(short());
        }
        self.offset += size as u64;
        Ok(Some(block))
    }
}

/// The compressed blocks of one slice of a span and where its first line
/// starts: `first_offset` into the first block's text, or after the line the
/// block before the slice ends inside when `open_line_before`. `start_block`
/// is the compressed offset of the span's first block, named when the span's
/// offset into it lies past its text.
struct Slice {
    blocks: Vec<Vec<u8>>,
    start_block: u64,
    first_offset: usize,
    open_line_before: bool,
}

/// The lines a block holds that start in it: `body` is the block's text from
/// its first line start to its last line end, and `straddle` the line that
/// starts in the block and ends in a later one, made whole.
struct LineRun<'a> {
    body: &'a [u8],
    straddle: Option<Vec<u8>>,
}

/// One [`LineRun`] per block of consecutive block texts that starts a line,
/// whether a line was left incomplete, and whether the blocks scanned end
/// inside a line (an empty block ends no line and leaves that as it was). A
/// block's first line start is `first_offset` in the first block; in every
/// other it is the block's first byte, or the byte after its first line end
/// when the block before ends inside a line (that line is the earlier block's
/// straddle; `open_line_before` says the same of the block before the first).
/// A block that only continues a line starts no run. A final line with no line
/// end is the file's last line when `file_ended` and is kept; otherwise it runs
/// past the blocks read and is reported incomplete. Under `last_is_tail` the
/// last block was read for the line that ends in it alone (the block past a
/// span's end) and starts no run; the blocks scanned end before it.
fn block_line_runs(
    blocks: &[Vec<u8>],
    first_offset: usize,
    open_line_before: bool,
    last_is_tail: bool,
    file_ended: bool,
) -> (Vec<LineRun<'_>>, bool, bool) {
    let mut runs = Vec::new();
    let mut incomplete = false;
    // Whether the bytes before the current block end inside a line.
    let mut open_line = open_line_before;
    for (i, block) in blocks.iter().enumerate() {
        if last_is_tail && i + 1 == blocks.len() {
            break;
        }
        let begin = if open_line {
            match block.iter().position(|&b| b == b'\n') {
                Some(at) => at + 1,
                None => continue,
            }
        } else if i == 0 {
            first_offset.min(block.len())
        } else {
            0
        };
        if let Some(&last) = block.last() {
            open_line = last != b'\n';
        }
        let rest = &block[begin..];
        let body_end = rest
            .iter()
            .rposition(|&b| b == b'\n')
            .map_or(0, |at| at + 1);
        let mut straddle = None;
        if body_end < rest.len() {
            let mut line = rest[body_end..].to_vec();
            let mut ended = false;
            for next in blocks.iter().skip(i + 1) {
                match next.iter().position(|&b| b == b'\n') {
                    Some(at) => {
                        line.extend_from_slice(&next[..=at]);
                        ended = true;
                        break;
                    }
                    None => line.extend_from_slice(next),
                }
            }
            if ended || file_ended {
                straddle = Some(line);
            } else {
                incomplete = true;
            }
        }
        if body_end > 0 || straddle.is_some() {
            runs.push(LineRun {
                body: &rest[..body_end],
                straddle,
            });
        }
    }
    (runs, incomplete, open_line)
}

/// Merge overlapping/adjacent intervals. Input must be sorted by start.
fn merge_intervals(intervals: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if intervals.is_empty() {
        return Vec::new();
    }

    let mut merged: Vec<(u64, u64)> = Vec::new();
    let (mut cur_start, mut cur_end) = intervals[0];

    for &(s, e) in &intervals[1..] {
        // A small gap is bridged too: one longer read beats a second seek.
        if s <= cur_end.saturating_add(GAP_BRIDGE) {
            cur_end = cur_end.max(e);
        } else {
            merged.push((cur_start, cur_end));
            cur_start = s;
            cur_end = e;
        }
    }
    merged.push((cur_start, cur_end));
    merged
}

fn normalize_range_bounds(start: u64, end: u64) -> (u64, u64) {
    if start <= end {
        (start, end)
    } else {
        (end, start)
    }
}

/// Contig-name spellings to try, in order, when resolving a query chromosome
/// against a data file's own naming.
///
/// Annotation files disagree: a table may use bare Ensembl names (`21`), the
/// spelling under which `CADD.pm` queries chromosome 1 (CADD.pm:147), while
/// AlphaMissense uses UCSC names (`chr21`). A query that matches
/// neither returns zero records **silently**, indistinguishable from a file with
/// no record at that locus, so both spellings (plus the `M`/`MT` mitochondrial
/// variants) are tried before giving up.
///
/// Public so the binary (`.vpd`) store applies the identical resolution.
pub(crate) fn chromosome_query_candidates(chr: &str) -> Vec<String> {
    let mut candidates: Vec<String> = Vec::new();
    let mut push_unique = |value: String| {
        if !value.is_empty() && !candidates.iter().any(|existing| existing == &value) {
            candidates.push(value);
        }
    };

    let trimmed = chr.trim().to_string();
    push_unique(trimmed.clone());

    let lower = trimmed.to_ascii_lowercase();
    if lower.starts_with("chr") {
        let stripped = trimmed[3..].to_string();
        push_unique(stripped.clone());
        if stripped.eq_ignore_ascii_case("m") {
            push_unique("MT".to_string());
        }
        if stripped.eq_ignore_ascii_case("mt") {
            push_unique("M".to_string());
        }
    } else {
        push_unique(format!("chr{trimmed}"));
        if trimmed.eq_ignore_ascii_case("m") {
            push_unique("MT".to_string());
            push_unique("chrM".to_string());
            push_unique("chrMT".to_string());
        } else if trimmed.eq_ignore_ascii_case("mt") {
            push_unique("M".to_string());
            push_unique("chrM".to_string());
            push_unique("chrMT".to_string());
        }
    }
    candidates
}

/// Whether `name` spells the mitochondrial contig: `M` or `MT`, bare or after
/// `chr`, in any letter case.
pub(crate) fn is_mitochondrial_name(name: &str) -> bool {
    let bare = match name.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &name[3..],
        _ => name,
    };
    bare.eq_ignore_ascii_case("m") || bare.eq_ignore_ascii_case("mt")
}

/// The spellings a plugin's table is asked for `chr` under with no check of
/// its rows, in order.
///
/// [`chromosome_query_candidates`], except for a mitochondrial name on GRCh37.
/// There `MT` and `chrMT` name the assembly's mitochondrial sequence and
/// `chrM` and `M` the one of the UCSC hg19 naming, another sequence whose
/// positions do not stand at the assembly's, so a row under one pair is a row
/// of a record under the other only in a table that counts on the record's
/// sequence. The name is tried as Ensembl VEP's tabix plugins try it
/// (`BaseVepTabixPlugin.pm:693-708`): as written, then without its `chr`
/// prefix or with `chr` added, which stays inside the pair; a prefix in
/// another letter case is tried in lower case as well. A reader tries the
/// other pair after these, behind its frame check
/// ([`grch37_mitochondrial_candidates`], [`MitochondrialFrame`]).
pub(crate) fn plugin_table_candidates(chr: &str, assembly: Option<Assembly>) -> Vec<String> {
    let name = chr.trim();
    if assembly != Some(Assembly::Grch37) || !is_mitochondrial_name(name) {
        return chromosome_query_candidates(chr);
    }
    let bare = match name.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &name[3..],
        _ => name,
    };
    let mut candidates = vec![name.to_string()];
    for other in [bare.to_string(), format!("chr{bare}")] {
        if !candidates.contains(&other) {
            candidates.push(other);
        }
    }
    candidates
}

/// The contig name a plugin asks its table for a record's rows under.
///
/// The cache's name of the record's contig, which the table's reader spells
/// every way ([`chromosome_query_candidates`]), except for a mitochondrial
/// record on GRCh37: there the name is the record's own spelling, the one an
/// Ensembl VEP module hands its base class (`$vf->{chr}`), which tells a
/// reader opened on that assembly the pair of names the record is of: it
/// tries the name one `chr` step ([`plugin_table_candidates`]) and reads a
/// contig of the other pair only where the table's rows state the reference
/// bases ([`MitochondrialFrame`]).
pub(crate) fn table_contig(variant: &InputVariant, assembly: Option<Assembly>) -> &str {
    if assembly == Some(Assembly::Grch37) && is_mitochondrial_name(&variant.chr) {
        &variant.original_chr
    } else {
        &variant.chr
    }
}

/// Whether `name` is of the pair `MT` and `chrMT`, the names GRCh37 gives its
/// own mitochondrial sequence; the other pair is `M` and `chrM`.
fn names_the_assembly_sequence(name: &str) -> bool {
    let bare = match name.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &name[3..],
        _ => name,
    };
    bare.eq_ignore_ascii_case("mt")
}

/// The spellings a table is asked for the mitochondrial name `chr` under on
/// GRCh37, in order, each with whether it is of the other pair of names than
/// `chr`.
///
/// `first`, the name a plugin's Ensembl VEP module gives every mitochondrial
/// record itself, leads with its `chr` variant. The record's own pair follows
/// ([`plugin_table_candidates`]), then the two names of the other pair, the
/// one with `chr`'s own prefix form first.
pub(crate) fn grch37_mitochondrial_candidates(
    chr: &str,
    first: Option<&str>,
) -> Vec<(String, bool)> {
    let name = chr.trim();
    let own_pair = names_the_assembly_sequence(name);
    let prefixed = name
        .get(..3)
        .is_some_and(|prefix| prefix.eq_ignore_ascii_case("chr"));
    let other_bare = if own_pair { "M" } else { "MT" };
    let other_pair = if prefixed {
        [format!("chr{other_bare}"), other_bare.to_string()]
    } else {
        [other_bare.to_string(), format!("chr{other_bare}")]
    };
    let mut candidates: Vec<(String, bool)> = Vec::new();
    let named_first = first
        .map(|first| plugin_table_candidates(first, Some(Assembly::Grch37)))
        .unwrap_or_default();
    for candidate in named_first
        .into_iter()
        .chain(plugin_table_candidates(name, Some(Assembly::Grch37)))
        .chain(other_pair)
    {
        if !candidates.iter().any(|(held, _)| *held == candidate) {
            let crosses = names_the_assembly_sequence(&candidate) != own_pair;
            candidates.push((candidate, crosses));
        }
    }
    candidates
}

/// The name the frame check reads the run's mitochondrial sequence under: the
/// one GRCh37 gives it, which the reference reader spells as its file does.
const ASSEMBLY_MITOCHONDRIAL_NAME: &str = "MT";

/// The equal stretches of the reference sequence the frame check reads a
/// table's rows from.
pub(crate) const FRAME_STRETCHES: u64 = 16;

/// The rows the frame check reads of one stretch, at most: its first ones.
pub(crate) const FRAME_ROWS_A_STRETCH: usize = 8;

/// The fewest positions at which the frame check compares rows before a
/// table is read across the names.
pub(crate) const FRAME_POSITIONS: usize = 8;

/// The rows of one contig of a table whose span overlaps `start..=end`, in the
/// table's order, at most `limit` of them: how the frame check reads a table
/// or the store built from it.
pub(crate) type FrameRows<'a> =
    dyn Fn(u64, u64, usize) -> Result<Vec<TabixRecord>, PluginError> + 'a;

/// What the frame check found of one mitochondrial contig of a table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FrameVerdict {
    /// Every compared row states the reference bases at its own position, over
    /// at least [`FRAME_POSITIONS`] positions: the contig is read.
    StatesTheReference,
    /// A compared row states another base than the reference holds at its
    /// position, or stands past the end of the sequence: the contig is not
    /// read.
    StatesAnotherSequence,
    /// The comparison could not be made, for the reason given: the contig is
    /// not read.
    NotCompared(&'static str),
}

/// The frame check of a plugin's table on GRCh37.
///
/// There `MT` and `chrMT` name the assembly's mitochondrial sequence and
/// `chrM` and `M` another sequence, of another length, whose positions do not
/// stand at the assembly's. A contig a table names under one pair is read for
/// a record of the other pair only where the table's rows state the bases the
/// run's reference sequence holds at the rows' own positions; the reference
/// allele of a row is the base it states.
///
/// The check reads the first `FRAME_ROWS_A_STRETCH` rows of each of
/// `FRAME_STRETCHES` equal stretches of the reference sequence, and one row
/// past its end, through the table's index. A row is compared when its
/// reference allele is made of `A`, `C`, `G` and `T` alone and the reference
/// holds such bases under it. The contig is read when every compared row
/// states the reference and the compared rows stand at
/// `FRAME_POSITIONS` positions or more, so one base that happens to agree
/// admits no table. The verdict is made once for each contig and kept; it is
/// logged with the table's path and the rows compared.
///
/// Without a reference sequence, or for a table with no reference-allele
/// column, the comparison cannot be made: the contig is not read, as Ensembl
/// VEP's base class reads none of the other pair
/// (`BaseVepTabixPlugin.pm:693-708`), and one warning names the table.
#[derive(Debug, Default)]
#[doc(hidden)]
pub struct MitochondrialFrame {
    reference: Option<Arc<vep_fasta::IndexedFasta>>,
    allele_column: Option<usize>,
    first_name: Option<&'static str>,
    own_pair_alone: bool,
    verdicts: Mutex<Vec<(String, FrameVerdict)>>,
}

impl MitochondrialFrame {
    /// A check that reads a row's reference allele from `allele_column`, a
    /// 0-based column, with no reference sequence yet.
    pub(crate) fn comparing_column(allele_column: Option<usize>) -> Self {
        Self {
            allele_column,
            ..Self::default()
        }
    }

    /// Hands the check the run's reference sequence.
    pub(crate) fn set_reference(&mut self, reference: Option<Arc<vep_fasta::IndexedFasta>>) {
        self.reference = reference;
    }

    /// Names the 0-based column that holds a row's reference allele.
    pub(crate) fn set_allele_column(&mut self, allele_column: Option<usize>) {
        self.allele_column = allele_column;
    }

    /// Names the spelling the plugin's Ensembl VEP module gives every
    /// mitochondrial record, which is then the first one tried
    /// ([`grch37_mitochondrial_candidates`]).
    pub(crate) fn set_first_name(&mut self, name: &'static str) {
        self.first_name = Some(name);
    }

    /// Has the reader keep every mitochondrial name to its own pair on GRCh37,
    /// with no comparison made and no warning: the rule of a source that
    /// Ensembl VEP reads through its own resolver (`BaseVEP.pm:689-738`),
    /// which reaches no name of the other pair but through a synonym.
    pub fn keep_to_own_pair(&mut self) {
        self.own_pair_alone = true;
    }

    /// The spelling [`Self::set_first_name`] named, if any.
    pub(crate) fn first_name(&self) -> Option<&'static str> {
        self.first_name
    }

    /// Whether the rows `table` holds under `contig` are read for a
    /// mitochondrial record of the other pair of names: the verdict kept for
    /// the contig, made from `rows` at the first call.
    pub(crate) fn admits(
        &self,
        table: &Path,
        contig: &str,
        rows: &FrameRows<'_>,
    ) -> Result<bool, PluginError> {
        if self.own_pair_alone {
            return Ok(false);
        }
        let mut verdicts = self
            .verdicts
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if let Some((_, verdict)) = verdicts.iter().find(|(held, _)| held == contig) {
            return Ok(*verdict == FrameVerdict::StatesTheReference);
        }
        let comparison = self.compare(rows)?;
        match comparison.verdict {
            FrameVerdict::NotCompared(reason) => tracing::warn!(
                "{}: its contig {contig} is not read for a mitochondrial record written {}: on \
                 GRCh37 the names MT and chrMT are one sequence and chrM and M another, and the \
                 table's rows could not be compared with the reference sequence ({reason})",
                table.display(),
                if names_the_assembly_sequence(contig) {
                    "chrM or M"
                } else {
                    "MT or chrMT"
                }
            ),
            FrameVerdict::StatesTheReference | FrameVerdict::StatesAnotherSequence => {}
        }
        tracing::info!(
            table = %table.display(),
            contig,
            rows_compared = comparison.compared,
            rows_stating_the_reference = comparison.stating,
            positions = comparison.positions,
            outcome = match comparison.verdict {
                FrameVerdict::StatesTheReference => "read across the mitochondrial names",
                FrameVerdict::StatesAnotherSequence | FrameVerdict::NotCompared(_) =>
                    "read under its own pair of names alone",
            },
            "mitochondrial rows compared with the reference sequence"
        );
        verdicts.push((contig.to_string(), comparison.verdict));
        Ok(comparison.verdict == FrameVerdict::StatesTheReference)
    }

    /// Reads the rows and compares them ([`Self`]).
    fn compare(&self, rows: &FrameRows<'_>) -> Result<FrameComparison, PluginError> {
        let not_compared = |reason| FrameComparison {
            verdict: FrameVerdict::NotCompared(reason),
            ..FrameComparison::default()
        };
        let Some(column) = self.allele_column else {
            return Ok(not_compared("the table has no reference-allele column"));
        };
        let Some(reference) = &self.reference else {
            return Ok(not_compared(
                "no reference sequence was given; pass --fasta",
            ));
        };
        let (Some(name), Some(length)) = (
            reference.record_name(ASSEMBLY_MITOCHONDRIAL_NAME),
            reference.sequence_length(ASSEMBLY_MITOCHONDRIAL_NAME),
        ) else {
            return Ok(not_compared(
                "the reference sequence holds no mitochondrial sequence",
            ));
        };
        // A sequence the reference names `chrM` or `M` at another length than
        // the assembly's is the other naming's own sequence, which a run's
        // synonym table can hand a read of `MT`.
        if !names_the_assembly_sequence(name) && length != vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH {
            return Ok(not_compared(
                "the reference sequence holds its mitochondrial sequence under the other pair \
                 of names, at another length than the assembly's",
            ));
        }
        let mut comparison = FrameComparison::default();
        // A span that reaches past the end of the sequence is read again
        // there, so a position is counted by its start, once.
        let mut positions: Vec<u64> = Vec::new();
        let mut compare = |row: &TabixRecord| {
            let Some(allele) = row.columns.get(column).map(|cell| cell.trim()) else {
                return;
            };
            if allele.is_empty() || !allele.bytes().all(is_base) {
                return;
            }
            let end = row.start.saturating_add(allele.len() as u64 - 1);
            let stated = if row.start == 0 || end > length {
                false
            } else {
                match reference.sequence(ASSEMBLY_MITOCHONDRIAL_NAME, row.start, end) {
                    Some(bases) if bases.iter().copied().all(is_base) => {
                        bases.eq_ignore_ascii_case(allele.as_bytes())
                    }
                    _ => return,
                }
            };
            comparison.compared += 1;
            comparison.stating += usize::from(stated);
            if !positions.contains(&row.start) {
                positions.push(row.start);
            }
        };
        let stretch = length.div_ceil(FRAME_STRETCHES).max(1);
        let mut start = 1;
        while start <= length {
            let end = (start + stretch - 1).min(length);
            for row in rows(start, end, FRAME_ROWS_A_STRETCH)?
                .iter()
                .filter(|row| row.start >= start)
            {
                compare(row);
            }
            start = end + 1;
        }
        for row in &rows(length + 1, u64::MAX, 1)? {
            compare(row);
        }
        comparison.positions = positions.len();
        comparison.verdict = if comparison.stating < comparison.compared {
            FrameVerdict::StatesAnotherSequence
        } else if comparison.positions < FRAME_POSITIONS {
            FrameVerdict::NotCompared("too few of its rows state a reference allele to tell")
        } else {
            FrameVerdict::StatesTheReference
        };
        Ok(comparison)
    }
}

/// Whether `byte` is one of the four bases, in either letter case.
fn is_base(byte: u8) -> bool {
    matches!(byte.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// What [`MitochondrialFrame::compare`] counted, with its verdict.
#[derive(Debug)]
struct FrameComparison {
    compared: usize,
    stating: usize,
    positions: usize,
    verdict: FrameVerdict,
}

impl Default for FrameComparison {
    fn default() -> Self {
        Self {
            compared: 0,
            stating: 0,
            positions: 0,
            verdict: FrameVerdict::NotCompared("no row was read"),
        }
    }
}

/// Resolve a query chromosome to the index's reference-sequence id and the
/// file's own spelling of the name, which its records carry.
fn resolve_reference_sequence<'a>(
    ref_seq_names: &'a indexmap::IndexSet<bstr::BString>,
    chr: &str,
    assembly: Option<Assembly>,
) -> Option<(usize, &'a bstr::BString)> {
    for candidate in plugin_table_candidates(chr, assembly) {
        let chr_bytes = bstr::BString::from(candidate.as_bytes());
        if let Some(found) = ref_seq_names.get_full(&chr_bytes) {
            return Some(found);
        }
    }
    None
}

/// Ensembl VEP's check on a plugin's data table before anything of it is read
/// (`BaseVepTabixPlugin.pm`, `check_file`): the table exists and `<file>.tbi`
/// sits beside it. A table indexed by `tabix -C` alone, or whose index is named
/// on the file's stem, is refused naming the missing `.tbi`, though htslib
/// alone would read it. A plugin that reads the table's header or index before
/// opening it runs this check first, as Ensembl VEP's plugins do.
pub(crate) fn check_plugin_table(data_path: &Path) -> Result<(), PluginError> {
    if !data_path.exists() {
        return Err(PluginError::Init(format!(
            "data file {} not found",
            data_path.display()
        )));
    }
    let tbi = PathBuf::from(format!("{}.tbi", data_path.display()));
    if !tbi.exists() {
        return Err(PluginError::Init(format!(
            "tabix index {} not found; a plugin data file needs its .tbi index beside \
             it, so index it with tabix without -C",
            tbi.display()
        )));
    }
    Ok(())
}

/// The paths the index of a data file is looked for at, in the order htslib
/// searches: `<file>.csi` (the format `tabix -C` writes), `<stem>.csi`,
/// `<file>.tbi`, `<stem>.tbi`, where `<stem>` drops the data file's last
/// extension; a data file without an extension has the two appended names.
#[doc(hidden)]
pub fn index_path_candidates(data_path: &Path) -> Vec<PathBuf> {
    let mut candidates = Vec::with_capacity(4);
    for extension in ["csi", "tbi"] {
        candidates.push(PathBuf::from(format!(
            "{}.{extension}",
            data_path.display()
        )));
        if data_path.extension().is_some() {
            candidates.push(data_path.with_extension(extension));
        }
    }
    candidates
}

/// Find the index file for a given data file: the first of
/// [`index_path_candidates`] that exists.
#[doc(hidden)]
pub fn find_index_path(data_path: &Path) -> Result<PathBuf, PluginError> {
    index_path_candidates(data_path)
        .into_iter()
        .find(|candidate| candidate.exists())
        .ok_or_else(|| {
            PluginError::Init(format!(
                "tabix index not found for {path}. Expected {path}.csi or {path}.tbi",
                path = data_path.display(),
            ))
        })
}

/// Load a tabix index from a `.tbi` or `.csi` file, by its extension, through
/// `read_tabix_index` or `read_csi_index`.
pub(crate) fn load_index(index_path: &Path) -> Result<Index, PluginError> {
    if index_path.extension().is_some_and(|ext| ext == "csi") {
        Ok(Box::new(read_csi_index(index_path)?))
    } else {
        Ok(Box::new(read_tabix_index(index_path)?))
    }
}

/// The query configuration of a table as its tabix index describes it: the
/// index's sequence, start and end columns (`tabix -s`, `-b`, `-e`) and a
/// zero-based start where the index is a BED one, so a table is queried as
/// `tabix` queries it, whatever its own header says. The configuration names
/// no allele column.
///
/// `no_header` words the refusal of an index that carries no tabix header: it
/// is handed the index's path and returns the caller's whole message.
///
/// A table with no index takes the configuration a `.vpd` store beside it
/// records (`config_from_store_beside` of the annotation store), the columns
/// the store was built on, so a plugin that reads its columns here before it
/// opens its store loads from an accepted store as the other store-backed
/// plugins do. A caller that reads the table through tabix alone checks the
/// table's index before it calls.
pub(crate) fn indexed_table_config(
    data_path: &Path,
    no_header: impl FnOnce(&Path) -> String,
) -> Result<TabixConfig, PluginError> {
    let index_path = match find_index_path(data_path) {
        Ok(index_path) => index_path,
        Err(no_index) => {
            return crate::annotation_store::config_from_store_beside(data_path, no_index)
        }
    };
    let index = load_index(&index_path)?;
    let header = index
        .header()
        .ok_or_else(|| PluginError::Init(no_header(&index_path)))?;
    Ok(TabixConfig {
        file_path: data_path.to_path_buf(),
        chr_col: header.reference_sequence_name_index(),
        start_col: header.start_position_index(),
        end_col: header.end_position_index(),
        ref_col: None,
        alt_col: None,
        zero_based: header.format().coordinate_system() == CoordinateSystem::Bed,
    })
}

/// Reads a `.tbi` index with its counts checked (`checked_index_stream`).
pub(crate) fn read_tabix_index(index_path: &Path) -> Result<tabix::Index, PluginError> {
    read_checked_index(index_path, false, |stream| {
        tabix::io::Reader::new(stream).read_index()
    })
}

/// Reads a `.csi` index with its counts checked (`checked_index_stream`).
pub(crate) fn read_csi_index(index_path: &Path) -> Result<noodles::csi::Index, PluginError> {
    read_checked_index(index_path, true, |stream| {
        noodles::csi::io::Reader::new(stream).read_index()
    })
}

/// What `parse` reads from the index at `index_path`, a `.csi` index under
/// `csi` and a `.tbi` index otherwise. `parse` is handed the index as the
/// BGZF stream an index reader takes, held in memory, and only once the
/// index's counts are checked (`checked_index_stream`), so a reader that
/// allocates for each count it meets reads counts that fit the index.
///
/// An index that cannot be opened or read is [`PluginError::Run`], naming the
/// index and the read error, and so is an error `parse` returns. An index
/// that decompresses past the size limit, or that the count check refuses, is
/// [`PluginError::Init`], naming the index and what is malformed in it;
/// `parse` is not called for it.
#[doc(hidden)]
pub fn read_checked_index<T>(
    index_path: &Path,
    csi: bool,
    parse: impl FnOnce(std::io::Cursor<Vec<u8>>) -> std::io::Result<T>,
) -> Result<T, PluginError> {
    let stream = checked_index_stream(index_path, csi)?;
    parse(stream).map_err(|e| index_read_error(index_path, e))
}

fn index_read_error(index_path: &Path, e: std::io::Error) -> PluginError {
    PluginError::Run(format!(
        "failed to read tabix index {}: {e}",
        index_path.display()
    ))
}

/// The bytes of the index at `index_path` as the BGZF stream the index readers
/// take, decompressed first, to at most [`MAX_INDEX_BYTES`], with every count
/// in them checked against the bytes that follow it ([`check_index_counts`])
/// before a reader, which allocates for each count as it meets it, reads the
/// same bytes.
fn checked_index_stream(
    index_path: &Path,
    csi: bool,
) -> Result<std::io::Cursor<Vec<u8>>, PluginError> {
    let read_error = |e: std::io::Error| index_read_error(index_path, e);
    let malformed = |why: String| {
        PluginError::Init(format!(
            "tabix index {} is malformed: {why}",
            index_path.display()
        ))
    };
    let file = std::fs::File::open(index_path).map_err(read_error)?;
    let mut bytes = Vec::new();
    bgzf::io::Reader::new(file)
        .take(MAX_INDEX_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(read_error)?;
    if bytes.len() > MAX_INDEX_BYTES {
        return Err(malformed(format!(
            "it decompresses to more than {MAX_INDEX_BYTES} bytes"
        )));
    }
    check_index_counts(&bytes, csi).map_err(malformed)?;
    let mut encoder = bgzf::io::Writer::new(Vec::new());
    encoder.write_all(&bytes).map_err(read_error)?;
    Ok(std::io::Cursor::new(encoder.finish().map_err(read_error)?))
}

/// Checks every count of a decompressed `.tbi` (or, under `csi`, `.csi`)
/// index against the bytes that follow it, in the layout the tabix and CSI
/// index format specifications give: the reference count, a name block's
/// length, each reference's bin count, each bin's chunk count and, in a
/// `.tbi`, each reference's interval count.
/// A count whose entries, at their smallest size, would run past the end of
/// the index is the error, naming the count; the bytes themselves are left to
/// the index reader.
fn check_index_counts(bytes: &[u8], csi: bool) -> Result<(), String> {
    const CHUNK_BYTES: usize = 16;
    const INTERVAL_BYTES: usize = 8;
    let mut cursor = IndexCursor { bytes, at: 0 };
    cursor.skip(4, "magic")?;
    if csi {
        cursor.skip(8, "min_shift and depth")?;
        let aux = cursor.count("l_aux")?;
        cursor.skip(aux, "aux")?;
    }
    let references = cursor.count("n_ref")?;
    if !csi {
        // The tabix header: six integers, then the name block.
        cursor.skip(24, "header")?;
        let names = cursor.count("l_nm")?;
        cursor.skip(names, "names")?;
    }
    // A reference is at least its bin count and, in a `.tbi`, its interval
    // count; a bin at least its id, its chunk count and, in a `.csi`, its
    // offset.
    let (reference_bytes, bin_bytes) = if csi { (4, 16) } else { (8, 8) };
    cursor.fits(references, reference_bytes, "n_ref")?;
    for _ in 0..references {
        let bins = cursor.count("n_bin")?;
        cursor.fits(bins, bin_bytes, "n_bin")?;
        for _ in 0..bins {
            cursor.skip(if csi { 12 } else { 4 }, "bin")?;
            let chunks = cursor.count("n_chunk")?;
            cursor.skip_entries(chunks, CHUNK_BYTES, "n_chunk")?;
        }
        if !csi {
            let intervals = cursor.count("n_intv")?;
            cursor.skip_entries(intervals, INTERVAL_BYTES, "n_intv")?;
        }
    }
    Ok(())
}

/// A position in the decompressed bytes of an index, for
/// [`check_index_counts`].
struct IndexCursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl IndexCursor<'_> {
    /// The bytes after the position.
    fn left(&self) -> usize {
        self.bytes.len() - self.at
    }

    /// Moves past `n` bytes of `what`, or reports that the index ends inside
    /// them.
    fn skip(&mut self, n: usize, what: &str) -> Result<(), String> {
        if n > self.left() {
            return Err(format!("the index ends inside its {what}"));
        }
        self.at += n;
        Ok(())
    }

    /// Reads the little-endian `int32_t` count `name`, refusing a negative one.
    fn count(&mut self, name: &str) -> Result<usize, String> {
        let Some(raw) = self.bytes.get(self.at..self.at + 4) else {
            return Err(format!("the index ends inside its {name}"));
        };
        self.at += 4;
        let value = i32::from_le_bytes(raw.try_into().expect("four bytes"));
        usize::try_from(value).map_err(|_| format!("{name} is {value}"))
    }

    /// Checks that `count` entries of at least `entry_bytes` each fit in the
    /// bytes left.
    fn fits(&self, count: usize, entry_bytes: usize, name: &str) -> Result<(), String> {
        let needed = count
            .checked_mul(entry_bytes)
            .ok_or_else(|| format!("{name} is {count}, more entries than the index can hold"))?;
        if needed > self.left() {
            return Err(format!(
                "{name} is {count}, more entries than the {} bytes left in the index can hold",
                self.left()
            ));
        }
        Ok(())
    }

    /// Moves past `count` entries of exactly `entry_bytes` each, once they
    /// [`fit`](Self::fits).
    fn skip_entries(&mut self, count: usize, entry_bytes: usize, name: &str) -> Result<(), String> {
        self.fits(count, entry_bytes, name)?;
        self.at += count * entry_bytes;
        Ok(())
    }
}

/// Read the 0-based column index that a tabix index was built on.
///
/// This is the ground truth for which column a query resolves against, and it is
/// independent of anything in the data file's header. A plugin that parses a
/// different column than the index was built on stamps every record with a
/// position from the wrong column, so the engine's own position filter discards
/// all of them: a silent zero, exit 0.
///
/// That is exactly what happens when a run is pointed at the other assembly's
/// copy of a dual-coordinate file: REVEL ships one table carrying both `hg19_pos`
/// (column 2) and `grch38_pos` (column 3), so the two copies carry the same
/// header and differ in this index field. No header inspection can tell them apart.
///
/// Returns `None` when no index is present, or when the index carries no header
/// block. Callers treat that as "cannot verify" rather than "verified compatible".
#[cfg(test)]
pub(crate) fn indexed_start_column(data_path: &Path) -> Result<Option<usize>, PluginError> {
    // Resolved through the same helper as the query path, so the appended
    // `<file>.tbi` form is tried and an extensionless data file (REVEL unzips to
    // a bare `revel_with_transcript_ids`) resolves.
    let Ok(index_path) = find_index_path(data_path) else {
        return Ok(None);
    };
    let index = load_index(&index_path)?;
    Ok(index.header().map(|h| h.start_position_index()))
}

/// The 0-based sequence-name and start columns a tabix index was built on
/// (`tabix -s` and `-b`), the pair `indexed_start_column` reads the second
/// of: `None` when no index is present or the index carries no header block.
pub(crate) fn indexed_columns(data_path: &Path) -> Result<Option<(usize, usize)>, PluginError> {
    let Ok(index_path) = find_index_path(data_path) else {
        return Ok(None);
    };
    let index = load_index(&index_path)?;
    Ok(index
        .header()
        .map(|h| (h.reference_sequence_name_index(), h.start_position_index())))
}

/// How a table's tabix index stands to the column a plugin reads a record's
/// position from ([`indexed_position`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum IndexedPosition {
    /// No index is present, or the index carries no header block, so the
    /// column cannot be checked.
    Unchecked,
    /// The index's start column is the plugin's; `sequence_column` is the
    /// index's 0-based sequence-name column (`tabix -s`).
    Agrees {
        /// The 0-based sequence-name column of the index.
        sequence_column: usize,
    },
    /// The index's start column is another one.
    Disagrees {
        /// The 0-based start column of the index (`tabix -b`).
        indexed_column: usize,
    },
}

/// The file's tabix index checked against `position_column`, the 0-based column
/// a plugin reads a record's position from ([`indexed_columns`]). The answer
/// is data: a plugin that resolved the column by its own rule words its own
/// refusal of a disagreement and takes the sequence-name column of an index
/// that agrees.
pub(crate) fn indexed_position(
    path: &Path,
    position_column: usize,
) -> Result<IndexedPosition, PluginError> {
    Ok(match indexed_columns(path)? {
        None => IndexedPosition::Unchecked,
        Some((sequence_column, start_column)) if start_column == position_column => {
            IndexedPosition::Agrees { sequence_column }
        }
        Some((_, indexed_column)) => IndexedPosition::Disagrees { indexed_column },
    })
}

/// The 0-based index of the column named `wanted_column_name` in `header`, the
/// column the plugin reads a record's position from, checked against the start
/// column of the file's tabix index ([`indexed_position`]).
///
/// A file carrying coordinates on more than one assembly (REVEL's `hg19_pos` and `grch38_pos`,
/// dbNSFP's `pos(1-based)` and `hg19_pos(1-based)`, dbscSNV's `pos` and `hg38_pos`) is indexed
/// on one of them, and Ensembl VEP's plugins choose the column by the run's assembly
/// (dbNSFP.pm:316-334, dbscSNV.pm:209-231) and die when the header lacks it. An index built
/// on the other column passes every header check and annotates nothing at exit 0,
/// since the engine filters records on the position it re-parses from the chosen
/// column, so that disagreement is refused here; a file with no index, or an
/// index with no header block, cannot be checked and is accepted. The message
/// names no plugin: the caller prefixes its own name.
pub(crate) fn resolve_indexed_position_column(
    path: &Path,
    header: &[String],
    wanted_column_name: &str,
) -> Result<usize, PluginError> {
    let column = header
        .iter()
        .position(|name| name.trim() == wanted_column_name)
        .ok_or_else(|| {
            PluginError::Init(format!(
                "the data file {} has no '{wanted_column_name}' column (header: {})",
                path.display(),
                header.join(",")
            ))
        })?;
    if let IndexedPosition::Disagrees {
        indexed_column: indexed,
    } = indexed_position(path, column)?
    {
        let indexed_name = header
            .get(indexed)
            .map(|name| name.trim().to_string())
            .unwrap_or_else(|| format!("column {}", indexed + 1));
        return Err(PluginError::Init(format!(
            "{} is indexed on '{indexed_name}' (column {}), not on '{wanted_column_name}' \
             (column {}); index the file on the '{wanted_column_name}' column",
            path.display(),
            indexed + 1,
            column + 1,
        )));
    }
    Ok(column)
}

/// The last of the file's leading `#` lines, as `tabix -H` prints them: the
/// lines are read from the top of the file and the first line that is empty or
/// does not start with `#` ends the preamble. `None` when no `#` line opens
/// the file.
///
/// The line is returned whole, every `#` included, for a plugin that applies
/// its own rule to the raw header text where [`parse_header`] strips the marks.
pub(crate) fn last_preamble_line(path: &Path) -> Result<Option<String>, PluginError> {
    let file = std::fs::File::open(path)
        .map_err(|e| PluginError::Init(format!("failed to open {}: {e}", path.display())))?;
    let reader = bgzf::io::Reader::new(BufReader::new(file));

    let mut last = None;
    for line_result in BoundedLines::new(reader, path) {
        let line = line_result?;
        if line.is_empty() || !line.starts_with(b"#") {
            break;
        }
        last = Some(line);
    }
    Ok(last.map(lossy_text))
}

/// Read one field from the file's first data row, by 0-based column index.
///
/// Some databases declare their build in a data column rather than in a comment
/// preamble or in their column names: AlphaMissense carries a `genome` column
/// holding the literal UCSC spelling (`hg19` / `hg38`) on every row.
/// [`parse_header`] reads the column names alone and cannot see it, so a
/// wrong-build file passes the header check and then annotates nothing, because
/// its positions do not exist in the requested build.
///
/// Returns `None` when the file has no data row or the row is too short. Reads only
/// as far as the first data line, so it costs one bgzf block.
pub(crate) fn first_data_row_field(path: &Path, col: usize) -> Result<Option<String>, PluginError> {
    let file = std::fs::File::open(path)
        .map_err(|e| PluginError::Init(format!("failed to open {}: {e}", path.display())))?;
    let reader = bgzf::io::Reader::new(BufReader::new(file));

    for line_result in BoundedLines::new(reader, path) {
        let line = lossy_text(line_result?);
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        return Ok(line
            .split('\t')
            .nth(col)
            .map(|f| f.trim().trim_end_matches('\r').to_string()));
    }
    Ok(None)
}

/// The header lines of a bgzipped file, the lines `tabix -h` prints: every
/// leading line opening `#`, each without its line end.
pub(crate) fn header_lines(path: &Path) -> Result<Vec<String>, PluginError> {
    let file = std::fs::File::open(path)
        .map_err(|e| PluginError::Init(format!("failed to open {}: {e}", path.display())))?;
    let reader = bgzf::io::Reader::new(BufReader::new(file));
    let mut lines = Vec::new();
    for line in BoundedLines::new(reader, path) {
        let line = line?;
        if !line.starts_with(b"#") {
            break;
        }
        lines.push(lossy_text(line));
    }
    Ok(lines)
}

/// Parse column names from a bgzf-compressed tabix file's header.
///
/// Handles two conventions, because plugins resolve their score columns **by
/// name** and a file whose header fails to parse annotates nothing at all
/// (every [`TabixRecord::get`] returns `None`) without raising an error:
///
/// 1. **`#`-prefixed** (the common case, e.g. the `#Chrom` line of CADD.pm:161).
///    All leading `#` lines are consumed and the last one is taken as the column
///    names, so a multi-line provenance preamble is handled.
/// 2. **Bare, un-prefixed** (the dbscSNV plugin documents a header line that
///    opens with `c`, indexed with `tabix -c c`: dbscSNV.pm:55-62, with no
///    leading `#`). Detected by checking the column tabix indexed as the start
///    position: a data row's cell leads with an integer, as `tbx_parse1` read
///    it to index the row ([`leading_integer`], the reading the record reader
///    takes), a header row's does not. That test is exact rather than
///    heuristic, which is why `start_col` is a parameter.
///
/// `\r` is stripped so a header line ending in CRLF does not leave a carriage
/// return welded to the final column name.
pub(crate) fn parse_header(
    path: &Path,
    start_col: usize,
) -> Result<Option<Vec<String>>, PluginError> {
    Ok(header_line(path, start_col)?.map(|line| header_names(&line)))
}

/// The line [`parse_header`] reads the column names from, as the file writes
/// it: the last of the leading `#` lines or, in a file with none, a first
/// line whose start cell no integer leads.
fn header_line(path: &Path, start_col: usize) -> Result<Option<String>, PluginError> {
    let file = std::fs::File::open(path)
        .map_err(|e| PluginError::Init(format!("failed to open {}: {e}", path.display())))?;
    let reader = bgzf::io::Reader::new(BufReader::new(file));

    let mut last_header_line = None;
    for line_result in BoundedLines::new(reader, path) {
        let line = lossy_text(line_result?);
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            last_header_line = Some(line.to_string());
            continue;
        }
        if last_header_line.is_none() && !field_leads_with_integer(line, start_col) {
            last_header_line = Some(line.to_string());
        }
        break;
    }

    Ok(last_header_line)
}

/// The column names of a header line: its cells less the leading `#` marks
/// and each cell's trailing `\r`.
fn header_names(line: &str) -> Vec<String> {
    line.trim_start_matches('#')
        .split('\t')
        .map(|s| s.trim_end_matches('\r').to_string())
        .collect()
}

/// True when the tab-separated `line`'s field at `col` leads with an integer,
/// the reading `tbx_parse1` indexed the row by and [`TabixAnnotator`] reads a
/// record's start with.
///
/// Distinguishes a bare header row from a data row: every row tabix indexed
/// opens its start cell with an integer, whatever follows it (`100x`, `100.5`),
/// and a row whose cell does not was never indexed as a data row.
fn field_leads_with_integer(line: &str, col: usize) -> bool {
    line.split('\t')
        .nth(col)
        .is_some_and(|cell| leading_integer(cell).is_some())
}

/// Allele matching utilities for tabix plugins.
pub(crate) mod allele_match {
    use super::{split_field_count, TabixRecord};
    use crate::perl_format::{numify_arithmetic, perl_number};
    use crate::PluginError;
    use std::borrow::Cow;
    use vep_core::allele_trim::{self, strip_shared_first_base, CaseRule};
    use vep_core::sv_routing::perl_true;
    use vep_core::variant::InputVariant;

    pub use vep_core::allele_trim::trim_directions;

    /// The position a module's numeric comparison reads from the record's start
    /// cell (`$_->{start} == $start`, BayesDel.pm:141, PrimateAI.pm:175-176):
    /// the cell's number as Perl's `==` numifies it (`numify_arithmetic`), 0
    /// where none leads the cell, as a `u64` when it is a whole non-negative
    /// number and `None` otherwise, since a fractional or negative number
    /// equals no variant position. The record's `start` is the cell's leading
    /// integer, where tabix indexed the row; the two differ on a cell with a
    /// fraction or an exponent (`100.5`, `1e2`).
    pub fn perl_position(record: &TabixRecord, start_col: usize) -> Option<u64> {
        let number = record
            .columns
            .get(start_col)
            .and_then(|cell| numify_arithmetic(cell))
            .unwrap_or(0.0);
        (number >= 0.0 && number.fract() == 0.0 && number <= u64::MAX as f64)
            .then_some(number as u64)
    }

    /// Normalize a VCF-style allele for comparison.
    /// Trims shared prefix and suffix, returning the prefix offset.
    pub fn normalize_alleles(file_ref: &str, file_alt: &str) -> (String, String, u64) {
        let ref_bytes = file_ref.as_bytes().to_vec();
        let alt_bytes = file_alt.as_bytes().to_vec();
        let mut ref_slice = ref_bytes.as_slice();
        let mut alt_slice = alt_bytes.as_slice();

        let mut prefix_len = 0usize;
        let min_len = ref_slice.len().min(alt_slice.len());
        while prefix_len < min_len
            && ref_slice[prefix_len].eq_ignore_ascii_case(&alt_slice[prefix_len])
        {
            prefix_len += 1;
        }

        ref_slice = &ref_slice[prefix_len..];
        alt_slice = &alt_slice[prefix_len..];

        let mut suffix_len = 0usize;
        let min_post_prefix = ref_slice.len().min(alt_slice.len());
        while suffix_len < min_post_prefix
            && ref_slice[ref_slice.len() - 1 - suffix_len]
                .eq_ignore_ascii_case(&alt_slice[alt_slice.len() - 1 - suffix_len])
        {
            suffix_len += 1;
        }

        if suffix_len > 0 {
            ref_slice = &ref_slice[..ref_slice.len() - suffix_len];
            alt_slice = &alt_slice[..alt_slice.len() - suffix_len];
        }

        let ref_str = if ref_slice.is_empty() {
            "-".to_string()
        } else {
            String::from_utf8_lossy(ref_slice).to_string()
        };
        let alt_str = if alt_slice.is_empty() {
            "-".to_string()
        } else {
            String::from_utf8_lossy(alt_slice).to_string()
        };

        (ref_str, alt_str, prefix_len as u64)
    }

    /// Whether a tabix record matches a variant by allele and position: the
    /// record's alleles, trimmed of their shared prefix and suffix
    /// ([`normalize_alleles`]), equal the variant's in either case, and its
    /// position cell as a module's `==` numifies it ([`perl_position`]), moved
    /// past that prefix, equals the variant's start. The file's alleles are on
    /// the forward strand, so a record given on the reverse strand is compared
    /// through its forward-strand alleles.
    pub fn matches_allele_with_position(
        record: &TabixRecord,
        variant: &InputVariant,
        pos_col: usize,
        ref_col: usize,
        alt_col: usize,
    ) -> bool {
        matched_allele_offset(record, variant, ref_col, alt_col).is_some_and(|offset| {
            perl_position(record, pos_col).and_then(|position| position.checked_add(offset))
                == Some(variant.start)
        })
    }

    /// The position of one side of `get_matched_variant_alleles` as the key
    /// `"$ref\_$alt\_$pos"` holds it (ensembl-variation
    /// Utils/Sequence.pm:1207, 1232).
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum KeyPosition<'a> {
        /// The variant's start, a number, plus the bases its trimming moved.
        Number(u64),
        /// A row's position cell as `split` read it, and the bases moved off
        /// its alleles' start, each a `$start++` (`trim_sequences`, `:998`,
        /// `:1013`).
        Cell { text: &'a str, moved: u64 },
    }

    impl KeyPosition<'_> {
        /// The position after one more shared base left the alleles' start.
        #[cfg(test)]
        pub(crate) fn moved(self) -> Self {
            match self {
                Self::Number(start) => Self::Number(start + 1),
                Self::Cell { text, moved } => Self::Cell {
                    text,
                    moved: moved + 1,
                },
            }
        }

        /// The position after `bases` more shared bases left the alleles' start.
        pub fn moved_by(self, bases: u64) -> Self {
            match self {
                Self::Number(start) => Self::Number(start + bases),
                Self::Cell { text, moved } => Self::Cell {
                    text,
                    moved: moved + bases,
                },
            }
        }

        /// Whether two keys hold the same position text (`Self::key_text`);
        /// `==` compares the representation instead. A cell with no base moved
        /// is in the key as its text, equal to a number only as that number's
        /// plain decimal (`30396092`; not `30396092.0`, `+30396092` or
        /// `030396092`). A cell with a base moved was incremented after
        /// `$start + (length($ref) - 1)` (`:971`) numified it, so it is in the
        /// key as the number Perl's `+` reads from it (`numify_arithmetic`,
        /// `0` where nothing numeric leads the cell) plus the bases moved,
        /// equal to a number when that sum is the number exactly: `30396092.0`
        /// moved once is `30396093`, `30396092.5` is no position.
        pub fn equals(self, other: Self) -> bool {
            match (self, other) {
                (Self::Number(start), Self::Number(other_start)) => start == other_start,
                (Self::Cell { text, moved: 0 }, Self::Number(start))
                | (Self::Number(start), Self::Cell { text, moved: 0 }) => {
                    is_plain_decimal(text, start)
                }
                _ => self.key_text() == other.key_text(),
            }
        }

        /// The position as the key's text holds it: a number's decimal, an
        /// unmoved cell's text, or a moved cell's number plus its moves written
        /// as Perl writes a number ([`perl_number`]).
        fn key_text(self) -> String {
            match self {
                Self::Number(start) => start.to_string(),
                Self::Cell { text, moved: 0 } => text.to_string(),
                Self::Cell { text, moved } => {
                    perl_number(numify_arithmetic(text).unwrap_or(0.0) + moved as f64)
                }
            }
        }
    }

    /// `$start++` on a cell no arithmetic has read yet (gnomADMt.pm:310): the
    /// empty string is `1`; a cell of letters then digits, either alone
    /// (`abc`, `21`, `09`), takes Perl's string increment, each character
    /// stepping within its range and carrying leftwards (`abd`, `22`, `10`,
    /// `Az` to `Ba`, `zz` to `aaa`); any other cell is its number as Perl's
    /// `+` reads it (`numify_arithmetic`, `0` where nothing numeric leads
    /// the cell) plus one, written as Perl writes a number ([`perl_number`]:
    /// `21.0` is `22`, `42.5` is `43.5`, `6x` is `7`).
    pub fn perl_increment(text: &str) -> String {
        if text.is_empty() {
            return "1".to_string();
        }
        let letters = text.bytes().take_while(u8::is_ascii_alphabetic).count();
        let magic = text.bytes().skip(letters).all(|byte| byte.is_ascii_digit());
        if !magic {
            return perl_number(numify_arithmetic(text).unwrap_or(0.0) + 1.0);
        }
        let mut bytes = text.as_bytes().to_vec();
        let mut index = bytes.len();
        loop {
            if index == 0 {
                let carried = match bytes[0] {
                    b'0' => b'1',
                    b'a' => b'a',
                    _ => b'A',
                };
                bytes.insert(0, carried);
                break;
            }
            index -= 1;
            match bytes[index] {
                b'z' => bytes[index] = b'a',
                b'Z' => bytes[index] = b'A',
                b'9' => bytes[index] = b'0',
                byte => {
                    bytes[index] = byte + 1;
                    break;
                }
            }
        }
        String::from_utf8(bytes).expect("ASCII letters and digits")
    }

    /// Whether `text` is `value` written as Perl stringifies an integer: its
    /// decimal digits alone, with no sign, point, exponent, blank or leading
    /// zero.
    fn is_plain_decimal(text: &str, value: u64) -> bool {
        let canonical = text == "0" || !text.starts_with('0');
        canonical
            && !text.is_empty()
            && text.bytes().all(|byte| byte.is_ascii_digit())
            && text.parse::<u64>() == Ok(value)
    }

    /// The prefix offset of the record's normalized alleles when they equal the
    /// variant's forward-strand alleles, case-insensitively; `None` when a
    /// column is absent or the alleles differ.
    fn matched_allele_offset(
        record: &TabixRecord,
        variant: &InputVariant,
        ref_col: usize,
        alt_col: usize,
    ) -> Option<u64> {
        let file_ref = record.columns.get(ref_col)?.as_str();
        let file_alt = record.columns.get(alt_col)?.as_str();
        let (norm_ref, norm_alt, offset) = normalize_alleles(file_ref, file_alt);
        let (var_ref, var_alt) = variant.forward_strand_alleles();
        let alleles_match = norm_ref.as_bytes().eq_ignore_ascii_case(&var_ref)
            && norm_alt.as_bytes().eq_ignore_ascii_case(&var_alt);
        alleles_match.then_some(offset)
    }

    /// An empty allele as Ensembl writes it.
    pub fn dashed(allele: &[u8]) -> Vec<u8> {
        allele_trim::dashed(allele).to_vec()
    }

    /// `trim_sequences` (`Utils/Sequence.pm:965-1038`) with empty alleles written
    /// `-`, compared as written ([`allele_trim::trim_sequences`]): the trimmed
    /// pair and its start, moved up by each base that left the left end.
    pub fn trim_sequences(
        reference: &[u8],
        alt: &[u8],
        start: i64,
        end_first: bool,
    ) -> (Vec<u8>, Vec<u8>, i64) {
        let trimmed = allele_trim::trim_sequences(reference, alt, end_first, CaseRule::AsWritten);
        (
            trimmed.reference.to_vec(),
            trimmed.alt.to_vec(),
            start + trimmed.left as i64,
        )
    }

    /// The minimal forms of a reference and alternate pair at `start`, as
    /// `get_matched_variant_alleles` (`Utils/Sequence.pm:1144-1247`) keys them:
    /// the pair trimmed from the left then the right and, when either allele is
    /// longer than one base, from the right then the left too. Two variants
    /// match when a form of one equals a form of the other.
    pub fn minimal_forms(reference: &[u8], alt: &[u8], start: i64) -> Vec<(Vec<u8>, Vec<u8>, i64)> {
        trim_directions(reference, alt)
            .iter()
            .map(|end_first| trim_sequences(reference, alt, start, *end_first))
            .collect()
    }

    /// Whether `record` and `variant` are one change as `get_matched_variant_alleles`
    /// matches them (ensembl-variation Utils/Sequence.pm:1144-1247): both sides'
    /// alleles are trimmed of the bases they share, from the left then the right
    /// and, when either allele is longer than one base, from the right then the
    /// left as well (`trim_sequences`, `:965-1038`), and the two match when a
    /// trimming of each gives the same reference, alternate and position
    /// ([`KeyPosition`]: the row's cell as text, or as a number once a shared
    /// base moved it; a row whose cell is Perl-false, `0` or empty, matches
    /// nothing, `:1173`). The variant's own alleles are trimmed too, so a
    /// substitution written with an unchanged leading base (`GG/GA`) matches
    /// the file's row for the one base that changes, as
    /// [`matches_allele_with_position`], which trims the file's alleles alone,
    /// does not. The variant's alleles are compared on the forward strand, as
    /// the file writes them.
    pub fn matched_variant_alleles(
        record: &TabixRecord,
        variant: &InputVariant,
        pos_col: usize,
        ref_col: usize,
        alt_col: usize,
    ) -> bool {
        let (Some(cell), Some(file_ref), Some(file_alt)) = (
            record.columns.get(pos_col),
            record.columns.get(ref_col),
            record.columns.get(alt_col),
        ) else {
            return false;
        };
        if !perl_true(cell) {
            return false;
        }
        let (var_ref, var_alt) = variant.forward_strand_alleles();
        let file_keys = TrimmedKeys::of(
            file_ref.as_bytes(),
            file_alt.as_bytes(),
            KeyPosition::Cell {
                text: cell,
                moved: 0,
            },
            CaseRule::EitherCase,
        );
        let variant_keys = TrimmedKeys::of(
            &var_ref,
            &var_alt,
            KeyPosition::Number(variant.start),
            CaseRule::EitherCase,
        );
        file_keys.any_equal(&variant_keys, CaseRule::EitherCase)
    }

    /// The row's reference and alternate cells as a module whose `parse_data`
    /// splits the line into a cells array (`split /\t/, $line`) hands them to
    /// `get_matched_variant_alleles`: `None` when either lies past the
    /// row's last filled cell, the undefined value that split leaves
    /// ([`split_field_count`]), which `trim_sequences` dies on before any
    /// comparison (ensembl-variation Utils/Sequence.pm:968, reached from
    /// `:1230` for every row the module's loop puts to the match). An empty cell
    /// before a filled one is the defined empty string and compares as `-`. A
    /// module splitting into a list of scalars keeps its trailing empty cells
    /// (Perl supplies a limit one above the list's length), so for it a cell is
    /// undefined only past the row's last physical cell, which `columns.get`
    /// tells.
    pub fn defined_alleles(
        record: &TabixRecord,
        ref_col: usize,
        alt_col: usize,
    ) -> Option<(&str, &str)> {
        let present = split_field_count(&record.columns);
        let cell = |column: usize| {
            record
                .columns
                .get(column)
                .filter(|_| column < present)
                .map(String::as_str)
        };
        Some((cell(ref_col)?, cell(alt_col)?))
    }

    /// The row as the `parse_data` of gnomADMt (gnomADMt.pm:306-314) and
    /// AlphaMissense (AlphaMissense.pm:293-300) hands it to the match: when its
    /// alleles differ in length and share their first base under `case`, that
    /// base goes, an allele left empty or as the one-character zero is written
    /// `-`, the record's start moves up one and the
    /// position cell steps as `$start++` steps a cell nothing has numified
    /// ([`perl_increment`]: `21.0` is `22`, `42.5` is `43.5`, `0100` is `0101`),
    /// the text the match then compares. The trim inside the match then reads
    /// the shortened pair, so a repeat indel written with its anchor base never
    /// keys the anchor's own start, which the right-then-left trim of the
    /// unstripped pair would reach. Any other row is handed on as it is.
    pub fn vcf_adjusted_row(
        record: &TabixRecord,
        pos_col: usize,
        ref_col: usize,
        alt_col: usize,
        case: CaseRule,
    ) -> Cow<'_, TabixRecord> {
        let (Some(reference), Some(alt)) =
            (record.columns.get(ref_col), record.columns.get(alt_col))
        else {
            return Cow::Borrowed(record);
        };
        let Some((stripped_reference, stripped_alt)) =
            strip_shared_first_base(reference.as_bytes(), alt.as_bytes(), case)
        else {
            return Cow::Borrowed(record);
        };
        // A first base that is one byte of a longer character leaves no text.
        let as_cell = |allele: &[u8]| {
            std::str::from_utf8(allele).map_or_else(|_| "-".to_string(), str::to_string)
        };
        let (stripped_reference, stripped_alt) =
            (as_cell(stripped_reference), as_cell(stripped_alt));
        let mut adjusted = record.clone();
        adjusted.start = adjusted.start.saturating_add(1);
        if let Some(cell) = adjusted.columns.get_mut(pos_col) {
            *cell = perl_increment(cell);
        }
        adjusted.columns[ref_col] = stripped_reference;
        adjusted.columns[alt_col] = stripped_alt;
        Cow::Owned(adjusted)
    }

    /// Whether the row's position cell is Perl-false, `0`, empty or missing:
    /// `!$b->{pos}` (Utils/Sequence.pm:1173) returns the match on such a cell
    /// before `trim_sequences` dies on an undefined allele (`:1230`), so the row
    /// matches nothing and fails nothing, whatever its allele cells hold.
    pub fn perl_false_position(record: &TabixRecord, pos_col: usize) -> bool {
        record
            .columns
            .get(pos_col)
            .is_none_or(|cell| !perl_true(cell))
    }

    /// [`matched_variant_alleles`] for a row whose alleles [`defined_alleles`]
    /// read, or `None` for a row Ensembl VEP's match dies on: a module whose
    /// loop reaches such a row before its match loses its fields for the line,
    /// so the caller fails the variant instead of skipping the row. A row whose
    /// position cell is Perl-false matches nothing before either
    /// ([`perl_false_position`]).
    pub fn matched_or_undefined(
        record: &TabixRecord,
        variant: &InputVariant,
        pos_col: usize,
        ref_col: usize,
        alt_col: usize,
    ) -> Option<bool> {
        if perl_false_position(record, pos_col) {
            return Some(false);
        }
        defined_alleles(record, ref_col, alt_col)
            .map(|_| matched_variant_alleles(record, variant, pos_col, ref_col, alt_col))
    }

    /// The failure of `plugin` on a row of `file` with no reference or alternate
    /// allele cell, the row Ensembl VEP's `get_matched_variant_alleles` dies on
    /// (Utils/Sequence.pm:968): the plugin's `went wrong` warning names the row.
    pub fn undefined_allele_row(plugin: &str, file: &str, record: &TabixRecord) -> PluginError {
        PluginError::Run(format!(
            "{plugin}: the row of {file} at {}:{} has no reference or alternate allele cell",
            record.chr, record.start
        ))
    }

    /// An allele pair trimmed of its shared bases, with the position the trim
    /// leaves: borrowed from the alleles it was cut from, `-` when nothing is
    /// left, as `trim_sequences` writes it.
    type TrimmedKey<'a> = (&'a [u8], &'a [u8], KeyPosition<'a>);

    /// The `(reference, alternate, start)` of an allele pair after each trimming
    /// order `get_matched_variant_alleles` tries (`_get_trim_directions`,
    /// Utils/Sequence.pm:1250-1258): the left-then-right trim, and the
    /// right-then-left trim when either allele is longer than one base.
    struct TrimmedKeys<'a> {
        left_first: TrimmedKey<'a>,
        right_first: Option<TrimmedKey<'a>>,
    }

    impl<'a> TrimmedKeys<'a> {
        /// The keys of a pair whose shared bases are found under `case`
        /// ([`allele_trim::trim_sequences`]), the position moved past each base
        /// that left the alleles' start.
        fn of(
            reference: &'a [u8],
            alternate: &'a [u8],
            position: KeyPosition<'a>,
            case: CaseRule,
        ) -> Self {
            let key = |end_first: bool| -> TrimmedKey<'a> {
                let trimmed = allele_trim::trim_sequences(reference, alternate, end_first, case);
                (
                    trimmed.reference,
                    trimmed.alt,
                    position.moved_by(trimmed.left as u64),
                )
            };
            let mut orders = trim_directions(reference, alternate).iter();
            Self {
                left_first: key(*orders.next().expect("the left-first order")),
                right_first: orders.next().map(|end_first| key(*end_first)),
            }
        }

        fn keys(&self) -> impl Iterator<Item = &TrimmedKey<'a>> {
            std::iter::once(&self.left_first).chain(self.right_first.iter())
        }

        /// Whether a trimming of each side gives the same reference, alternate
        /// and position, the alleles compared under `case`.
        fn any_equal(&self, other: &TrimmedKeys<'_>, case: CaseRule) -> bool {
            self.keys().any(|(reference, alternate, position)| {
                other
                    .keys()
                    .any(|(other_reference, other_alternate, other_position)| {
                        case.same_allele(reference, other_reference)
                            && case.same_allele(alternate, other_alternate)
                            && position.equals(*other_position)
                    })
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::variant::InputVariant;

    /// A committed fixture under `tests/data`, each written by `bgzip --binary`
    /// and indexed twice, by `tabix` (`.tbi`) and `tabix -C` (`.csi`):
    /// `spans.tsv.gz` (`-s 1 -b 2 -e 3`; records that start before a query and
    /// reach into it) and `blocks.tsv.gz` (`-s 1 -b 2 -e 2`; 3,400 rows on
    /// chromosome 21 at positions `5010000 + 7i`, enough that the row at
    /// 5031749 straddles the boundary between the file's two bgzf blocks).
    /// `limit.tsv.gz` (`-s 1 -b 2 -e 3`) carries a `.tbi` alone.
    /// `limit_csi.tsv.gz` (`-s 1 -b 2 -e 3`) and `limit_csi_start.tsv.gz`
    /// (`-s 1 -b 2 -e 2`) are indexed by `tabix -C` and carry a `.csi` alone.
    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    /// Copy a fixture and one of its indexes into `dir`, so that only
    /// that index sits beside the data file.
    fn fixture_with_index(dir: &Path, name: &str, index_extension: &str) -> PathBuf {
        let data = dir.join(name);
        std::fs::copy(fixture(name), &data).unwrap();
        let index = format!("{name}.{index_extension}");
        std::fs::copy(fixture(&index), dir.join(&index)).unwrap();
        data
    }

    fn config(file_path: PathBuf, end_col: Option<usize>) -> TabixConfig {
        TabixConfig {
            file_path,
            chr_col: 0,
            start_col: 1,
            end_col,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        }
    }

    fn open(file_path: PathBuf, end_col: Option<usize>) -> TabixAnnotator {
        TabixAnnotator::open(config(file_path, end_col)).unwrap()
    }

    /// Ensembl VEP's `split /\t/` keeps an interior empty field and drops the
    /// trailing empty fields alone, a whitespace-only field being a field.
    #[test]
    fn split_field_count_drops_the_trailing_empty_cells_alone() {
        let strings =
            |texts: &[&str]| -> Vec<String> { texts.iter().map(|text| text.to_string()).collect() };
        assert_eq!(split_field_count(&strings(&["a", "", "b", "", ""])), 3);
        assert_eq!(split_field_count(&strings(&["a", "", " "])), 3);
        assert_eq!(split_field_count(&strings(&["", ""])), 0);
        assert_eq!(split_field_count(&[]), 0);
    }

    /// The records as `tabix` prints them, with the chromosome column restored
    /// in place (`parse_line` moves it out of `columns`).
    fn lines<'a>(records: impl IntoIterator<Item = &'a TabixRecord>) -> Vec<String> {
        records
            .into_iter()
            .map(|record| {
                let mut columns = record.columns.clone();
                columns[0] = record.chr.clone();
                columns.join("\t")
            })
            .collect()
    }

    /// The lines `query_windows` returns for `windows` of chromosome 21, which
    /// `query_windows_parallel` returns too.
    fn window_lines(ann: &TabixAnnotator, windows: &[(u64, u64)]) -> Vec<String> {
        let by_sweep = lines(&ann.query_windows("21", windows).unwrap());
        assert_eq!(
            lines(&ann.query_windows_parallel("21", windows).unwrap()),
            by_sweep,
            "{windows:?}"
        );
        by_sweep
    }

    /// Expected lines are the output of `tabix spans.tsv.gz <region>` for each
    /// region: a record is printed when its start..=end overlaps the region,
    /// so a record starting before the region start is included. Read here
    /// through the `.tbi`.
    #[test]
    fn record_spanning_query_start_is_returned() {
        let dir = tempfile::tempdir().unwrap();
        let ann = open(
            fixture_with_index(dir.path(), "spans.tsv.gz", "tbi"),
            Some(2),
        );

        // tabix spans.tsv.gz 21:200-210
        assert_eq!(
            lines(&ann.query("21", 200, 210).unwrap()),
            [
                "21\t150\t260\tspan150_260",
                "21\t180\t320\tspan180_320",
                "21\t200\t200\tp200",
            ]
        );
        // tabix spans.tsv.gz 21:250-255
        assert_eq!(
            lines(&ann.query("21", 250, 255).unwrap()),
            ["21\t150\t260\tspan150_260", "21\t180\t320\tspan180_320"]
        );
        // tabix spans.tsv.gz 21:300-300
        assert_eq!(
            lines(&ann.query("21", 300, 300).unwrap()),
            ["21\t180\t320\tspan180_320", "21\t300\t300\tp300"]
        );
        // tabix spans.tsv.gz 21:321-399 prints nothing.
        assert!(ann.query("21", 321, 399).unwrap().is_empty());
        // tabix spans.tsv.gz 22:200-200
        assert_eq!(
            lines(&ann.query("22", 200, 200).unwrap()),
            ["22\t150\t260\tchr22_span150_260"]
        );

        // The batch path reads the merged region once and answers each
        // variant's lookup by span.
        let batch = ann
            .query_batch(&[("21", 200, 210), ("21", 250, 255), ("22", 200, 200)])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 200, 210)),
            [
                "21\t150\t260\tspan150_260",
                "21\t180\t320\tspan180_320",
                "21\t200\t200\tp200",
            ]
        );
        assert_eq!(
            lines(batch.lookup("21", 250, 255)),
            ["21\t150\t260\tspan150_260", "21\t180\t320\tspan180_320"]
        );
        assert_eq!(
            lines(batch.lookup("22", 200, 200)),
            ["22\t150\t260\tchr22_span150_260"]
        );
    }

    /// Regions farther apart than `merge_intervals` bridges are read separately,
    /// and a record longer than the gap between them overlaps both reads. The
    /// batch holds it once, so a lookup returns it once, as `tabix` prints each
    /// line once.
    #[test]
    fn batch_holds_a_record_spanning_two_merged_regions_once() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = [
            "#chrom\tstart\tend\tname",
            "21\t100\t3000\tspan100_3000",
            "21\t150\t150\tp150",
            "21\t2500\t2500\tp2500",
        ];
        let path = write_tabix_fixture(dir.path(), "long.tsv.gz", &lines_in, IndexSpec::standard());
        let ann = open(path, Some(2));

        let batch = ann
            .query_batch(&[("21", 2500, 2500), ("21", 150, 150)])
            .unwrap();
        assert_eq!(batch.records_by_pos["21"][&100].len(), 1);
        assert_eq!(
            lines(batch.lookup("21", 150, 150)),
            ["21\t100\t3000\tspan100_3000", "21\t150\t150\tp150"]
        );
        assert_eq!(
            lines(batch.lookup("21", 2500, 2500)),
            ["21\t100\t3000\tspan100_3000", "21\t2500\t2500\tp2500"]
        );
        assert!(batch.lookup("21", 3001, 3001).is_empty());
    }

    /// Both chromosomes of `spans.tsv.gz` share one bgzf block, and chromosome
    /// 22's record follows chromosome 21's last one at the same positions, so a
    /// read that runs past chromosome 21's chunk takes it in. `tabix
    /// spans.tsv.gz 21:150-260` prints the four chromosome 21 records and
    /// `tabix spans.tsv.gz 22:150-260` the one chromosome 22 record.
    #[test]
    fn records_of_another_chromosome_in_the_same_block_are_not_returned() {
        let ann = open(fixture("spans.tsv.gz"), Some(2));
        assert_eq!(
            lines(&ann.query("21", 150, 260).unwrap()),
            [
                "21\t150\t260\tspan150_260",
                "21\t180\t320\tspan180_320",
                "21\t190\t195\tspan190_195",
                "21\t200\t200\tp200",
            ]
        );
        assert_eq!(
            lines(&ann.query("22", 150, 260).unwrap()),
            ["22\t150\t260\tchr22_span150_260"]
        );
        // tabix spans.tsv.gz 22:100-100 prints nothing: p100 is chromosome 21's.
        assert!(ann.query("22", 100, 100).unwrap().is_empty());
    }

    /// `tabix spans.tsv.gz 1:100-200` prints nothing, as `22:100-100` does, and
    /// only the index tells the two apart: the file names 21 and 22 (under the
    /// spellings a query resolves) and no other sequence.
    #[test]
    fn holds_chromosome_reads_the_index_names() {
        let ann = open(fixture("spans.tsv.gz"), Some(2));
        assert!(ann.query("1", 100, 200).unwrap().is_empty());
        assert!(ann.holds_chromosome("21"));
        assert!(ann.holds_chromosome("chr21"));
        assert!(ann.holds_chromosome("22"));
        assert!(!ann.holds_chromosome("1"));
        assert!(!ann.holds_chromosome("chr1"));
    }

    /// With only `<file>.csi` beside the data file, the index is found and
    /// answers as `tabix` does with that index: `tabix spans.tsv.gz 21:200-210`
    /// and `tabix blocks.tsv.gz 21:5033779-5033793` print the lines below.
    #[test]
    fn csi_index_is_found_beside_tbi() {
        let dir = tempfile::tempdir().unwrap();

        let spans = fixture_with_index(dir.path(), "spans.tsv.gz", "csi");
        assert_eq!(
            find_index_path(&spans).unwrap(),
            dir.path().join("spans.tsv.gz.csi")
        );
        assert_eq!(indexed_start_column(&spans).unwrap(), Some(1));
        let ann = open(spans, Some(2));
        assert_eq!(
            lines(&ann.query("21", 200, 210).unwrap()),
            [
                "21\t150\t260\tspan150_260",
                "21\t180\t320\tspan180_320",
                "21\t200\t200\tp200",
            ]
        );

        let blocks = fixture_with_index(dir.path(), "blocks.tsv.gz", "csi");
        let ann = open(blocks, None);
        assert_eq!(
            lines(&ann.query("21", 5033779, 5033793).unwrap()),
            [
                "21\t5033779\tC\tG\t0.843",
                "21\t5033786\tG\tT\t0.762",
                "21\t5033793\tT\tA\t0.681",
            ]
        );
        assert_eq!(ann.query("21", 5010000, 5033793).unwrap().len(), 3400);
    }

    /// The committed fixture carries both indexes; `.csi` is the one opened, as
    /// htslib opens it (`tabix` reads `<file>.csi` before `<file>.tbi`), so a
    /// `.tbi` left beside a later `tabix -C` index is not the one read.
    #[test]
    fn csi_preferred_when_both_exist() {
        let data = fixture("spans.tsv.gz");
        assert!(data.with_file_name("spans.tsv.gz.tbi").exists());
        assert_eq!(
            find_index_path(&data).unwrap(),
            data.with_file_name("spans.tsv.gz.csi")
        );

        let ann = open(data, Some(2));
        assert_eq!(
            lines(&ann.query("21", 100, 100).unwrap()),
            ["21\t100\t100\tp100"]
        );

        let orphan = tempfile::tempdir().unwrap();
        let no_index = orphan.path().join("spans.tsv.gz");
        std::fs::copy(fixture("spans.tsv.gz"), &no_index).unwrap();
        let path = no_index.display();
        assert_eq!(
            find_index_path(&no_index).unwrap_err().to_string(),
            format!(
                "plugin initialization failed: tabix index not found for {path}. \
                 Expected {path}.csi or {path}.tbi"
            )
        );
    }

    /// htslib also finds an index named after the data file's stem
    /// (`stem.tsv.tbi` beside `stem.tsv.gz`): `tabix stem.tsv.gz 21:100-100`
    /// prints `21 99 A` from such a pair, and only the `.tbi` is read.
    #[test]
    fn index_named_after_the_data_file_stem_is_found() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = ["#chrom\tpos\tname", "21\t99\tA", "21\t100\tB", "21\t199\tC"];
        let data = write_tabix_fixture(dir.path(), "stem.tsv.gz", &lines_in, IndexSpec::standard());
        std::fs::rename(
            dir.path().join("stem.tsv.gz.tbi"),
            dir.path().join("stem.tsv.tbi"),
        )
        .unwrap();

        assert_eq!(
            find_index_path(&data).unwrap(),
            dir.path().join("stem.tsv.tbi")
        );
        let ann = TabixAnnotator::open(TabixConfig {
            zero_based: true,
            ..config(data, None)
        })
        .unwrap();
        assert_eq!(lines(&ann.query("21", 100, 100).unwrap()), ["21\t99\tA"]);
    }

    /// Ensembl VEP's tabix plugins read a table only when `<file>.tbi` sits
    /// beside it (`BaseVepTabixPlugin.pm`, `check_file`): a table indexed by
    /// `tabix -C` alone, or whose index is named on the file's stem, is
    /// refused at initialization naming the missing `.tbi`, where htslib, and
    /// so a `--custom` file, reads the index it finds.
    #[test]
    fn plugin_table_without_its_tbi_is_refused_naming_it() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let csi_only = fixture_with_index(dir.path(), "spans.tsv.gz", "csi");
        let refusal = |data: &Path| {
            TabixAnnotator::open_plugin_table(config(data.to_path_buf(), Some(2)))
                .err()
                .map(|e| e.to_string())
        };
        assert_eq!(
            refusal(&csi_only),
            Some(format!(
                "plugin initialization failed: tabix index {path}.tbi not found; a plugin \
                 data file needs its .tbi index beside it, so index it with tabix without -C",
                path = csi_only.display()
            ))
        );
        assert_eq!(
            lines(&open(csi_only, Some(2)).query("21", 100, 100).unwrap()),
            ["21\t100\t100\tp100"]
        );

        let lines_in = ["#chrom\tpos\tname", "21\t99\tA", "21\t100\tB"];
        let stem = write_tabix_fixture(dir.path(), "stem.tsv.gz", &lines_in, IndexSpec::standard());
        std::fs::rename(
            dir.path().join("stem.tsv.gz.tbi"),
            dir.path().join("stem.tsv.tbi"),
        )
        .unwrap();
        let message = refusal(&stem).unwrap();
        assert!(message.contains("stem.tsv.gz.tbi not found"), "{message}");
        assert!(TabixAnnotator::open(config(stem, Some(2))).is_ok());

        let unindexed = dir.path().join("unindexed.tsv.gz");
        std::fs::copy(fixture("spans.tsv.gz"), &unindexed).unwrap();
        let message = refusal(&unindexed).unwrap();
        assert!(
            message.contains("unindexed.tsv.gz.tbi not found"),
            "{message}"
        );

        let missing = dir.path().join("missing.tsv.gz");
        assert_eq!(
            refusal(&missing),
            Some(format!(
                "plugin initialization failed: data file {} not found",
                missing.display()
            ))
        );
    }

    /// With `<file>.tbi` beside the table the plugin opener reads the index
    /// htslib reads: `<file>.csi` when both are present, else the `.tbi`.
    #[test]
    fn plugin_table_with_its_tbi_opens_through_the_index_htslib_reads() {
        let dir = tempfile::tempdir().unwrap();
        let tbi_only = fixture_with_index(dir.path(), "spans.tsv.gz", "tbi");
        let ann = TabixAnnotator::open_plugin_table(config(tbi_only, Some(2))).unwrap();
        assert_eq!(ann.index_path, dir.path().join("spans.tsv.gz.tbi"));
        assert_eq!(
            lines(&ann.query("21", 100, 100).unwrap()),
            ["21\t100\t100\tp100"]
        );

        let both = fixture("spans.tsv.gz");
        let ann = TabixAnnotator::open_plugin_table(config(both.clone(), Some(2))).unwrap();
        assert_eq!(ann.index_path, both.with_file_name("spans.tsv.gz.csi"));
        assert_eq!(
            lines(&ann.query("21", 100, 100).unwrap()),
            ["21\t100\t100\tp100"]
        );
    }

    /// A record line is trimmed of its trailing whitespace, so a trailing cell
    /// that is empty or whitespace is not a column of the record, unless the
    /// reader keeps the whitespace before the line end, as Ensembl VEP's tabix
    /// iterator hands the line to a plugin: then every cell of the line is a
    /// column, through `query` and `query_batch` alike.
    #[test]
    fn trailing_whitespace_is_trimmed_unless_the_reader_keeps_it() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = ["#chrom\tpos\tname\tnote", "21\t100\tA\t", "21\t101\tB\t  "];
        let data =
            write_tabix_fixture(dir.path(), "cells.tsv.gz", &lines_in, IndexSpec::standard());

        let trimmed = open(data.clone(), None);
        assert_eq!(
            lines(&trimmed.query("21", 100, 101).unwrap()),
            ["21\t100\tA", "21\t101\tB"]
        );
        let kept = open(data, None).keeping_trailing_whitespace();
        assert_eq!(
            lines(&kept.query("21", 100, 101).unwrap()),
            ["21\t100\tA\t", "21\t101\tB\t  "]
        );
        let batch = kept.query_batch(&[("21", 100, 101)]).unwrap();
        assert_eq!(
            lines(batch.lookup("21", 100, 101)),
            ["21\t100\tA\t", "21\t101\tB\t  "]
        );
    }

    /// A reader told to build a column from the line's bytes hands that cell
    /// to its decoder as the file holds it, `\r` trimmed, so an escape there can
    /// complete a raw first byte (`decode_type`: every `%XX` spelled as its
    /// byte, then the bytes read as text); every other cell, and every cell of
    /// a reader without the decoder, is the line's text read with U+FFFD for a
    /// byte that is not UTF-8, the line's trailing whitespace of any script and
    /// a cell's `\r` trimmed, through `query`, `query_windows`,
    /// `query_windows_parallel` and `query_batch` alike. A reader told two
    /// columns builds each from the bytes through its own decoder (the second
    /// here spelling a byte that is not UTF-8 as its escape) and leaves the
    /// rest text; a column told twice takes the later decoder.
    #[test]
    fn a_column_decoded_from_bytes_leaves_the_other_cells_as_text() {
        use crate::custom::gff_overlay::attributes_text;
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        use vep_core::regulatory_gff::decode_type;

        let lines_in: &[&[u8]] = &[
            b"#chrom\tpos\tkind\tnote",
            b"21\t100\tmix\xCA%80\tplain",
            b"21\t101\tmix%CA\x80\r\tcaf\xCA%80\r",
            b"21\t102\tcaf%FF\t\xFF\xC2\xA0 ",
            b"21\t103\t\xFFa%41\tz\xE3\x80\x80",
        ];
        let dir = tempfile::tempdir().unwrap();
        let data = write_tabix_fixture(dir.path(), "bytes.tsv.gz", lines_in, IndexSpec::standard());
        let plain = open(data.clone(), None);
        let from_bytes = open(data.clone(), None).decoding_column_from_bytes(2, decode_type);
        let expected_plain = [
            "21\t100\tmix\u{FFFD}%80\tplain",
            "21\t101\tmix%CA\u{FFFD}\tcaf\u{FFFD}%80",
            "21\t102\tcaf%FF\t\u{FFFD}",
            "21\t103\t\u{FFFD}a%41\tz",
        ];
        let expected_decoded = [
            "21\t100\tmix\u{280}\tplain",
            "21\t101\tmix\u{280}\tcaf\u{FFFD}%80",
            "21\t102\tcaf\u{FFFD}\t\u{FFFD}",
            "21\t103\t\u{FFFD}aA\tz",
        ];
        assert_eq!(lines(&plain.query("21", 100, 103).unwrap()), expected_plain);
        assert_eq!(
            lines(&from_bytes.query("21", 100, 103).unwrap()),
            expected_decoded
        );
        let two_columns = open(data, None)
            .decoding_column_from_bytes(3, decode_type)
            .decoding_column_from_bytes(2, decode_type)
            .decoding_column_from_bytes(3, attributes_text);
        assert_eq!(
            lines(&two_columns.query("21", 100, 103).unwrap()),
            [
                "21\t100\tmix\u{280}\tplain",
                "21\t101\tmix\u{280}\tcaf%CA%80",
                "21\t102\tcaf\u{FFFD}\t%FF",
                "21\t103\t\u{FFFD}aA\tz",
            ]
        );
        assert_eq!(
            lines(
                &from_bytes
                    .query_windows("21", &[(100, 101), (103, 103)])
                    .unwrap()
            ),
            [
                expected_decoded[0],
                expected_decoded[1],
                expected_decoded[3]
            ]
        );
        assert_eq!(
            lines(
                &from_bytes
                    .query_windows_parallel("21", &[(100, 101), (103, 103)])
                    .unwrap()
            ),
            [
                expected_decoded[0],
                expected_decoded[1],
                expected_decoded[3]
            ]
        );
        let batch = from_bytes.query_batch(&[("21", 100, 103)]).unwrap();
        assert_eq!(lines(batch.lookup("21", 100, 103)), expected_decoded);
    }

    /// A zero-based file whose index names one position column (`tabix -0 -s 1
    /// -b 2 -e 2`) matches on the shifted position whether the end column is
    /// given as that same column or left out: `tabix bed0.tsv.gz 21:100-100`
    /// prints `21 99 A`, 21:101-101 `21 100 B`, 21:200-200 `21 199 C`, 21:99-99
    /// nothing, and 21:1-1000 all three rows.
    #[test]
    fn zero_based_point_column_matches_on_the_shifted_position() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = ["#chrom\tpos\tname", "21\t99\tA", "21\t100\tB", "21\t199\tC"];
        let data = write_tabix_fixture(dir.path(), "bed0.tsv.gz", &lines_in, IndexSpec::standard());

        for end_col in [Some(1), None] {
            let ann = TabixAnnotator::open(TabixConfig {
                zero_based: true,
                ..config(data.clone(), end_col)
            })
            .unwrap();
            assert_eq!(lines(&ann.query("21", 100, 100).unwrap()), ["21\t99\tA"]);
            assert_eq!(lines(&ann.query("21", 101, 101).unwrap()), ["21\t100\tB"]);
            assert_eq!(lines(&ann.query("21", 200, 200).unwrap()), ["21\t199\tC"]);
            assert!(ann.query("21", 99, 99).unwrap().is_empty());
            assert_eq!(ann.query("21", 1, 1000).unwrap().len(), 3);
            let batch = ann.query_batch(&[("21", 100, 100)]).unwrap();
            assert_eq!(lines(batch.lookup("21", 100, 100)), ["21\t99\tA"]);
        }
    }

    /// A start cell with text after its digits is read to the digits, as
    /// htslib's `strtol` reads it, so `tabix` indexes and prints the row at
    /// that position (`tabix -s1 -b2 -e2` on this file prints `21 100x B` and
    /// `21 100.5 C` for `21:100-100`); a cell no digit leads is no record. The
    /// record keeps the cell as written for a module's own reading of it.
    #[test]
    fn a_start_cell_is_read_to_its_leading_integer() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = [
            "#chrom\tpos\tname",
            "21\t99\tA",
            "21\t100x\tB",
            "21\t100.5\tC",
            "21\t 101\tD",
            "21\t102\tE",
        ];
        let data = write_tabix_fixture(
            dir.path(),
            "prefix.tsv.gz",
            &lines_in,
            IndexSpec::standard(),
        );
        let ann = open(data, None);
        let at_100 = ann.query("21", 100, 100).unwrap();
        assert_eq!(lines(&at_100), ["21\t100x\tB", "21\t100.5\tC"]);
        let record = &at_100[0];
        assert_eq!(
            (record.start, record.end, record.columns[1].as_str()),
            (100, 100, "100x")
        );
        assert_eq!(lines(&ann.query("21", 101, 101).unwrap()), ["21\t 101\tD"]);
        assert_eq!(ann.line_start("21\tx100\tF"), None);
        assert_eq!(ann.line_start("21\t+100\tF"), Some(100));
        assert_eq!(ann.line_start("21\t\tF"), None);
        assert_eq!(leading_integer("99999999999999999999"), None);
    }

    /// A row whose start cell a minus sign leads is no record by default, the
    /// reading every plugin's table gets. Under
    /// [`TabixAnnotator::reading_rows_as_tabix_indexes_them`] it is one that
    /// starts where a cell of `0` starts a row, at 0 or, zero-based, at 1, as
    /// `tabix` indexes it at the first base (htslib 1.9 `tbx.c:99`): every
    /// range from there to its end returns it, in the file's order, and the
    /// record keeps the cell as written for its number, digits past the
    /// integer range included. A `-` no digit follows leads no number for
    /// either reader.
    #[test]
    fn a_start_below_zero_is_a_record_only_for_a_reader_that_keeps_it() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let data = write_tabix_fixture(
            dir.path(),
            "below.tsv.gz",
            &[
                "21\t-5\t300\tbelow",
                "21\t-3x\t120\tbelow_too",
                "21\t100\t100\tp100",
                "21\t200\t260\tspan",
            ],
            IndexSpec::standard().with_span_cols(1, 2, 3),
        );
        let dropping = open(data.clone(), Some(2));
        assert!(dropping.query("21", 1, 1).unwrap().is_empty());
        assert_eq!(
            lines(&dropping.query("21", 100, 100).unwrap()),
            ["21\t100\t100\tp100"]
        );
        assert_eq!(dropping.line_start("21\t-5\t300"), None);

        let keeping = open(data.clone(), Some(2)).reading_rows_as_tabix_indexes_them();
        let at_first_base = keeping.query("21", 1, 1).unwrap();
        assert_eq!(
            lines(&at_first_base),
            ["21\t-5\t300\tbelow", "21\t-3x\t120\tbelow_too"]
        );
        let record = &at_first_base[0];
        assert_eq!(
            (record.start, record.end, record.columns[1].as_str()),
            (0, 300, "-5")
        );
        assert_eq!(
            lines(&keeping.query("21", 100, 100).unwrap()),
            [
                "21\t-5\t300\tbelow",
                "21\t-3x\t120\tbelow_too",
                "21\t100\t100\tp100"
            ]
        );
        assert_eq!(
            lines(&keeping.query("21", 121, 150).unwrap()),
            ["21\t-5\t300\tbelow"]
        );
        assert!(keeping.query("21", 301, 310).unwrap().is_empty());
        let batch = keeping
            .query_batch(&[("21", 250, 250), ("21", 110, 110)])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 110, 110)),
            ["21\t-5\t300\tbelow", "21\t-3x\t120\tbelow_too"]
        );
        assert_eq!(
            lines(batch.lookup("21", 250, 250)),
            ["21\t-5\t300\tbelow", "21\t200\t260\tspan"]
        );
        for cell in ["-5", " -5", "-0", "-5x", "-1e2", "-9223372036854775809"] {
            assert_eq!(keeping.line_start(&format!("21\t{cell}\t300")), Some(0));
            assert_eq!(
                keeping.line_start_bytes(format!("21\t{cell}\t300\n").as_bytes()),
                Some(0)
            );
        }
        for cell in ["-", "-x", "- 5", "--5", "+-5", "-+5"] {
            assert_eq!(keeping.line_start(&format!("21\t{cell}\t300")), None);
        }

        let zero_based = TabixAnnotator::open(TabixConfig {
            zero_based: true,
            ..config(data, Some(2))
        })
        .unwrap()
        .reading_rows_as_tabix_indexes_them();
        let record = &zero_based.query("21", 1, 1).unwrap()[0];
        assert_eq!((record.start, record.columns[1].as_str()), (1, "-5"));
    }

    /// A row whose generic end cell no number leads, or a minus sign leads,
    /// ends at its start by default, the reading every plugin's table gets.
    /// Under [`TabixAnnotator::reading_rows_as_tabix_indexes_them`] it is no
    /// record, as the index build of `tabix` fails the line (htslib 1.9
    /// `tbx.c:105-106`, `:142`), the rows around it are read as they are by
    /// default, and end digits past the range end a row at the last position.
    #[test]
    fn a_row_whose_end_tabix_fails_is_a_record_only_for_a_reader_that_falls_back() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let data = write_tabix_fixture(
            dir.path(),
            "ends.tsv.gz",
            &[
                "21\t50\t60\tbefore",
                "21\t100\txyz\tno_number",
                "21\t100\t-9\tminus",
                "21\t100\t\tempty",
                "21\t100\t+300\tplus",
                "21\t100\t99999999999999999999\tpast",
                "21\t120\t130\tafter",
            ],
            IndexSpec::standard().with_span_cols(1, 2, 3),
        );
        let spans = |records: &[TabixRecord]| -> Vec<(String, u64, u64)> {
            records
                .iter()
                .map(|record| (record.columns[3].clone(), record.start, record.end))
                .collect()
        };
        let named = |rows: &[(&str, u64, u64)]| -> Vec<(String, u64, u64)> {
            rows.iter()
                .map(|(name, start, end)| (name.to_string(), *start, *end))
                .collect()
        };

        let falling_back = open(data.clone(), Some(2));
        assert_eq!(
            spans(&falling_back.query("21", 100, 100).unwrap()),
            named(&[
                ("no_number", 100, 100),
                ("minus", 100, 100),
                ("empty", 100, 100),
                ("plus", 100, 300),
                ("past", 100, 100),
            ])
        );

        let as_tabix = open(data, Some(2)).reading_rows_as_tabix_indexes_them();
        assert_eq!(
            spans(&as_tabix.query("21", 100, 100).unwrap()),
            named(&[("plus", 100, 300), ("past", 100, u64::MAX)])
        );
        assert_eq!(
            spans(&as_tabix.query("21", 55, 125).unwrap()),
            named(&[
                ("before", 50, 60),
                ("plus", 100, 300),
                ("past", 100, u64::MAX),
                ("after", 120, 130),
            ])
        );
    }

    /// A minus sign leads a number where a digit follows it, after any blank
    /// and whatever follows the digit, however many digits there are.
    #[test]
    fn a_minus_sign_leads_a_number_where_a_digit_follows_it() {
        for cell in [
            "-5",
            " -5",
            "\x0b-5",
            "-5x",
            "-5.9",
            "-1e2",
            "-0",
            "-0x1A",
            "-9223372036854775809",
        ] {
            assert!(minus_leads_a_number(cell), "{cell:?}");
        }
        for cell in [
            "", "5", "+5", "-", "-x", "-.5", "- 5", "--5", "+-5", "-+5", "x-5",
        ] {
            assert!(!minus_leads_a_number(cell), "{cell:?}");
        }
    }

    /// `record_end` reads a generic end cell to its leading integer, as the
    /// generic preset of `tabix` reads it (`300x` is 300); under the VCF preset,
    /// the configuration naming a `ref_col`, it takes the INFO column's `END` at
    /// or past the start, blanks before the value skipped (`END= 101`), and
    /// never reads the column as a number (an INFO opening with the digit-led
    /// key `1000G`, or that key alone, is no end), falling back on the
    /// reference allele's span from the start, one base for an empty or missing
    /// reference cell; without an end column the row ends at its start.
    #[test]
    fn record_end_reads_the_end_cell_as_tabix_does() {
        let cells =
            |texts: &[&str]| -> Vec<String> { texts.iter().map(|t| t.to_string()).collect() };
        let generic = config(PathBuf::from("rows.tsv.gz"), Some(2));
        assert_eq!(generic.record_end(100, &cells(&["21", "100", "300"])), 300);
        assert_eq!(generic.record_end(100, &cells(&["21", "100", "300x"])), 300);
        assert_eq!(generic.record_end(100, &cells(&["21", "100", " 300"])), 300);
        assert_eq!(generic.record_end(100, &cells(&["21", "100", "x"])), 100);
        assert_eq!(generic.record_end(100, &cells(&["21", "100"])), 100);
        let vcf = TabixConfig {
            end_col: Some(7),
            ref_col: Some(3),
            ..generic.clone()
        };
        let row = |info: &str| cells(&["21", "100", ".", "ACGT", "A", ".", ".", info]);
        assert_eq!(vcf.record_end(100, &row("END=101")), 101);
        assert_eq!(vcf.record_end(100, &row("END= 101")), 101);
        assert_eq!(vcf.record_end(100, &row("END=500;SVTYPE=DEL")), 500);
        assert_eq!(vcf.record_end(100, &row("1000G;END=300")), 300);
        assert_eq!(vcf.record_end(100, &row("1000G")), 103);
        assert_eq!(vcf.record_end(100, &row("END=50")), 103);
        assert_eq!(vcf.record_end(100, &row("DP=3")), 103);
        assert_eq!(
            vcf.record_end(100, &cells(&["21", "100", ".", "ACGT"])),
            103
        );
        assert_eq!(vcf.record_end(100, &cells(&["21", "100", "."])), 100);
        assert_eq!(
            vcf.record_end(100, &cells(&["21", "100", ".", "", "A", ".", ".", "DP=3"])),
            100
        );
        let point = TabixConfig {
            end_col: None,
            ..generic
        };
        assert_eq!(point.record_end(100, &cells(&["21", "100", "300"])), 100);
    }

    /// A row whose start column is the `u64` maximum. Zero-based, it has no
    /// one-based start, so a query over the file prints the rows before it and
    /// never a row at position 0, the start's sum modulo 2^64; one-based, with
    /// a four-base REF and an end column naming no end, its end saturates at
    /// the maximum.
    #[test]
    fn start_at_the_u64_maximum_is_skipped_zero_based_and_its_end_saturates_one_based() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let max_row = format!("21\t{}\tACGT\tDP=3", u64::MAX);
        let lines_in = [
            "#chrom\tpos\tref\tinfo",
            "21\t99\tA\tDP=1",
            "21\t100\tACGT\tDP=2",
            max_row.as_str(),
        ];
        let data = write_tabix_fixture(dir.path(), "max.tsv.gz", &lines_in, IndexSpec::standard());

        let zero_based = TabixAnnotator::open(TabixConfig {
            zero_based: true,
            ..config(data.clone(), None)
        })
        .unwrap();
        let before_it = ["21\t99\tA\tDP=1", "21\t100\tACGT\tDP=2"];
        assert_eq!(lines(&zero_based.query("21", 0, 1000).unwrap()), before_it);
        assert_eq!(
            lines(&zero_based.query_windows("21", &[(0, 1000)]).unwrap()),
            before_it
        );
        assert_eq!(zero_based.line_start(&max_row), None);
        assert!(zero_based
            .parse_line(&max_row, max_row.as_bytes())
            .is_none());

        let one_based = TabixAnnotator::open(TabixConfig {
            ref_col: Some(2),
            ..config(data, Some(3))
        })
        .unwrap();
        let record = one_based.parse_line(&max_row, max_row.as_bytes()).unwrap();
        assert_eq!((record.start, record.end), (u64::MAX, u64::MAX));
    }

    /// A `.csi` written without the tabix header (the shape an index of a BCF
    /// has) names no chromosome, so every query through it would return
    /// nothing; opening the file is refused instead.
    #[test]
    fn csi_without_a_tabix_header_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let data = fixture_with_index(dir.path(), "spans.tsv.gz", "csi");
        let csi = dir.path().join("spans.tsv.gz.csi");
        noodles::csi::fs::write(&csi, &noodles::csi::Index::default()).unwrap();

        let err = match TabixAnnotator::open(config(data, Some(2))) {
            Ok(_) => panic!("a header-less index must be refused"),
            Err(e) => e.to_string(),
        };
        assert_eq!(
            err,
            format!(
                "plugin initialization failed: index {} carries no tabix header, so its \
                 records cannot be mapped to chromosomes; index the file with tabix",
                csi.display()
            )
        );
    }

    /// Write `lines_in` as one bgzf block with a `.tbi` whose one bin holds one
    /// chunk per entry of `chunks_in` (a range of line indexes and the
    /// compressed offset the chunk names as its block; the data's one block
    /// sits at 0), the shape a span-indexed file's bins give a query: several
    /// chunks apart in the file.
    fn write_multi_chunk_fixture(
        dir: &Path,
        lines_in: &[&str],
        chunks_in: &[(std::ops::Range<usize>, u64)],
    ) -> PathBuf {
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use std::io::Write;

        let data_path = dir.join("chunks.tsv.gz");
        let body: String = lines_in.iter().map(|l| format!("{l}\n")).collect();
        let mut writer = bgzf::io::Writer::new(std::fs::File::create(&data_path).unwrap());
        writer.write_all(body.as_bytes()).unwrap();
        writer.try_finish().unwrap();

        let mut line_offsets = vec![0usize];
        for line in lines_in {
            line_offsets.push(line_offsets.last().unwrap() + line.len() + 1);
        }
        let position = |block_offset: u64, offset: usize| {
            bgzf::VirtualPosition::try_from((block_offset, u16::try_from(offset).unwrap())).unwrap()
        };
        let chunks: Vec<Chunk> = chunks_in
            .iter()
            .map(|(range, block_offset)| {
                Chunk::new(
                    position(*block_offset, line_offsets[range.start]),
                    position(*block_offset, line_offsets[range.end]),
                )
            })
            .collect();

        write_tbi(&data_path, chunks);
        data_path
    }

    /// Write `<data_path>.tbi` for a span-indexed file of chromosome 21 whose
    /// one bin holds `chunks`.
    fn write_tbi(data_path: &Path, chunks: Vec<Chunk>) {
        let first = chunks[0].start();
        write_tbi_bins(data_path, vec![(0, chunks)], vec![first]);
    }

    /// Write `<data_path>.tbi` for a span-indexed file of chromosome 21 with
    /// `bins` (bin number, its chunks) and `linear` the least offset of each
    /// 16 kb tile from the first.
    fn write_tbi_bins(
        data_path: &Path,
        bins: Vec<(usize, Vec<Chunk>)>,
        linear: Vec<bgzf::VirtualPosition>,
    ) {
        use noodles::csi::binning_index::index::header;
        use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
        use noodles::csi::binning_index::index::reference_sequence::Bin;
        use noodles::csi::binning_index::index::ReferenceSequence;

        let bins = bins
            .into_iter()
            .map(|(number, chunks)| (number, Bin::new(chunks)))
            .collect();
        let reference = ReferenceSequence::<LinearIndex>::new(bins, linear, None);
        let header = header::Builder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(Some(2))
            .set_reference_sequence_names(["21".as_bytes().into()].into_iter().collect())
            .build();
        let index = tabix::Index::builder()
            .set_header(header)
            .set_reference_sequences(vec![reference])
            .build();
        tabix::fs::write(format!("{}.tbi", data_path.display()), &index).unwrap();
    }

    /// A chunk that ends 7 bytes short of a buffered read's 8,192-byte fill
    /// (817 lines of 10 bytes and one of 15 after the header) is followed by
    /// lines at 20005: a read that fills past the chunk's end cuts the first of
    /// them to `21\t2000`, a fragment that parses as a record inside a wide
    /// query, so `tabix` would print 818 lines and the engine 819.
    #[test]
    fn fragment_of_a_line_cut_at_a_buffer_boundary_is_not_a_record() {
        let dir = tempfile::tempdir().unwrap();
        let mut rows: Vec<String> = vec!["#chrom\tpos\tname".to_string()];
        rows.extend((1000..1817).map(|pos| format!("21\t{pos}\tX")));
        rows.push("21\t1817\tpadded".to_string());
        rows.extend((20005..20045).map(|pos| format!("21\t{pos}\tX")));
        let lines_in: Vec<&str> = rows.iter().map(String::as_str).collect();
        assert_eq!(
            lines_in[1..819].iter().map(|l| l.len() + 1).sum::<usize>(),
            8185
        );
        let data = write_multi_chunk_fixture(dir.path(), &lines_in, &[(1..819, 0)]);
        let ann = open(data, Some(1));

        let found = ann.query("21", 1, 16384).unwrap();
        assert_eq!(found.len(), 818);
        assert!(found
            .iter()
            .all(|r| r.columns.len() == 3 && (1000..=1817).contains(&r.start)));
    }

    /// An index built for another file names offsets its blocks do not hold;
    /// the query is refused with both files named instead of reading past the
    /// block. A chunk whose offset the block does hold is read as usual.
    #[test]
    fn index_of_another_file_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("blocks.tsv.gz");
        std::fs::copy(fixture("spans.tsv.gz"), &data).unwrap();
        std::fs::copy(
            fixture("blocks.tsv.gz.tbi"),
            dir.path().join("blocks.tsv.gz.tbi"),
        )
        .unwrap();
        let ann = open(data.clone(), None);

        assert!(ann.query("21", 5010000, 5010007).unwrap().is_empty());
        let err = match ann.query("21", 5033779, 5033793) {
            Ok(_) => panic!("an offset past the block must be refused"),
            Err(e) => e.to_string(),
        };
        let path = data.display();
        assert_eq!(
            err,
            format!(
                "plugin execution failed: index {path}.tbi does not describe {path}: no block \
                 at offset 0 holds byte 59707"
            )
        );
    }

    /// An index whose chunk begins past the end of the data file reaches no
    /// block; a chunk begins at a record, so the pair is refused rather than
    /// answered with nothing. The seek past the end leaves the block read for
    /// an earlier chunk in place, so a second chunk past the end is likewise
    /// refused, not read again from that block.
    #[test]
    fn index_reaching_past_the_end_of_the_data_file_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let lines_in = ["#chrom\tstart\tend\tname", "21\t100\t100\tp100"];
        let past_the_end = 1_000_000;
        let refused = |data: &Path| {
            let ann = open(data.to_path_buf(), Some(2));
            let err = match ann.query("21", 100, 100) {
                Ok(_) => panic!("a chunk past the end of the file must be refused"),
                Err(e) => e.to_string(),
            };
            let path = data.display();
            assert_eq!(
                err,
                format!(
                    "plugin execution failed: index {path}.tbi does not describe {path}: no \
                     block at offset {past_the_end} holds byte 0"
                )
            );
        };

        let first_chunk_past_the_end =
            write_multi_chunk_fixture(dir.path(), &lines_in, &[(0..2, past_the_end)]);
        refused(&first_chunk_past_the_end);

        let other_dir = tempfile::tempdir().unwrap();
        let second_chunk_past_the_end = write_multi_chunk_fixture(
            other_dir.path(),
            &lines_in,
            &[(0..2, 0), (0..2, past_the_end)],
        );
        refused(&second_chunk_past_the_end);
    }

    /// Writes `decompressed` as a BGZF stream at `path`.
    fn write_bgzf(path: &Path, decompressed: &[u8]) {
        let mut writer = bgzf::io::Writer::new(std::fs::File::create(path).unwrap());
        writer.write_all(decompressed).unwrap();
        writer.try_finish().unwrap();
    }

    /// The decompressed bytes of a `.tbi` up to its first reference's bin
    /// count: a header with the format and columns of `spans.tsv.gz.tbi` and
    /// one reference, `21`, then `n_bin`.
    fn tbi_head(n_bin: i32) -> Vec<u8> {
        let mut bytes = b"TBI\x01".to_vec();
        for value in [1i32, 0, 1, 2, 3, i32::from(b'#'), 0, 3] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        bytes.extend_from_slice(b"21\0");
        bytes.extend_from_slice(&n_bin.to_le_bytes());
        bytes
    }

    /// An index whose bin, interval, chunk or reference count claims more
    /// entries than the bytes after it can hold is refused when the file is
    /// opened, naming the count, before the index reader allocates for it; so
    /// is a negative count, a name block running past the end, a `.csi` with
    /// the same shape, and an index decompressing past the size cap. An index
    /// whose counts fit opens and answers queries.
    #[test]
    fn index_counts_are_checked_against_its_size() {
        let dir = tempfile::tempdir().unwrap();
        let data = fixture_with_index(dir.path(), "spans.tsv.gz", "tbi");
        let tbi_path = dir.path().join("spans.tsv.gz.tbi");
        let refused = |index_path: &Path, decompressed: &[u8], why: &str| {
            write_bgzf(index_path, decompressed);
            let error = TabixAnnotator::open(config(data.clone(), Some(2)))
                .err()
                .map(|e| e.to_string())
                .unwrap_or_else(|| panic!("an index with {why} is refused"));
            assert_eq!(
                error,
                format!(
                    "plugin initialization failed: tabix index {} is malformed: {why}",
                    index_path.display()
                )
            );
        };
        let mut many_bins = tbi_head(i32::MAX);
        many_bins.extend_from_slice(&[0; 4]);
        refused(
            &tbi_path,
            &many_bins,
            "n_bin is 2147483647, more entries than the 4 bytes left in the index can hold",
        );
        let mut negative_bins = tbi_head(-1);
        negative_bins.extend_from_slice(&[0; 4]);
        refused(&tbi_path, &negative_bins, "n_bin is -1");
        let mut many_chunks = tbi_head(1);
        many_chunks.extend_from_slice(&4681u32.to_le_bytes());
        many_chunks.extend_from_slice(&i32::MAX.to_le_bytes());
        refused(
            &tbi_path,
            &many_chunks,
            "n_chunk is 2147483647, more entries than the 0 bytes left in the index can hold",
        );
        let mut many_intervals = tbi_head(0);
        many_intervals.extend_from_slice(&1_000_000i32.to_le_bytes());
        many_intervals.extend_from_slice(&[0; 8]);
        refused(
            &tbi_path,
            &many_intervals,
            "n_intv is 1000000, more entries than the 8 bytes left in the index can hold",
        );
        let mut many_references = tbi_head(0);
        many_references[4..8].copy_from_slice(&i32::MAX.to_le_bytes());
        refused(
            &tbi_path,
            &many_references,
            "n_ref is 2147483647, more entries than the 4 bytes left in the index can hold",
        );
        let mut long_names = tbi_head(0);
        long_names[32..36].copy_from_slice(&1_000i32.to_le_bytes());
        refused(&tbi_path, &long_names, "the index ends inside its names");
        refused(
            &tbi_path,
            &tbi_head(0)[..39],
            "n_ref is 1, more entries than the 0 bytes left in the index can hold",
        );
        let csi_path = dir.path().join("spans.tsv.gz.csi");
        let mut csi = b"CSI\x01".to_vec();
        for value in [14i32, 5, 0, 1, i32::MAX] {
            csi.extend_from_slice(&value.to_le_bytes());
        }
        refused(
            &csi_path,
            &csi,
            "n_bin is 2147483647, more entries than the 0 bytes left in the index can hold",
        );
        std::fs::remove_file(&csi_path).unwrap();
        {
            let mut writer = bgzf::io::Writer::new(std::fs::File::create(&tbi_path).unwrap());
            let block = vec![0u8; 1 << 16];
            for _ in 0..(MAX_INDEX_BYTES >> 16) + 1 {
                writer.write_all(&block).unwrap();
            }
            writer.try_finish().unwrap();
        }
        let error = TabixAnnotator::open(config(data.clone(), Some(2)))
            .err()
            .map(|e| e.to_string())
            .expect("an oversize index is refused");
        assert_eq!(
            error,
            format!(
                "plugin initialization failed: tabix index {} is malformed: it decompresses to \
                 more than {MAX_INDEX_BYTES} bytes",
                tbi_path.display()
            )
        );
        std::fs::copy(fixture("spans.tsv.gz.tbi"), &tbi_path).unwrap();
        let records = open(data, Some(2)).query("21", 100, 100).unwrap();
        assert!(!records.is_empty());
    }

    /// A data line longer than the cap is refused where it is met: by the
    /// header read when the file is opened, and by a query naming the line's
    /// block offset when it sits past the header.
    #[test]
    fn data_line_past_the_cap_is_refused() {
        use noodles::csi::binning_index::index::header;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
        use noodles::csi::binning_index::index::reference_sequence::Bin;
        use noodles::csi::binning_index::index::ReferenceSequence;
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("long.tsv.gz");
        // The header read when the file is opened stops at the first data line,
        // so the long line is the second one, 35 bytes in.
        let mut body = b"#chrom\tstart\tend\tname\n21\t50\t50\tp50\n21\t100\t100\tp100\t".to_vec();
        body.extend(std::iter::repeat_n(b'x', MAX_LINE_BYTES));
        body.push(b'\n');
        write_bgzf(&data, &body);
        let position = |offset: u64| {
            bgzf::VirtualPosition::try_from((0u64, u16::try_from(offset).unwrap())).unwrap()
        };
        let chunk = Chunk::new(position(35), position(60));
        let mut bins = indexmap::IndexMap::new();
        bins.insert(4681, Bin::new(vec![chunk]));
        let reference = ReferenceSequence::<LinearIndex>::new(bins, vec![chunk.start()], None);
        let header = header::Builder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(Some(2))
            .set_reference_sequence_names(["21".as_bytes().into()].into_iter().collect())
            .build();
        let index = tabix::Index::builder()
            .set_header(header)
            .set_reference_sequences(vec![reference])
            .build();
        tabix::fs::write(format!("{}.tbi", data.display()), &index).unwrap();
        let error = open(data.clone(), Some(2))
            .query("21", 100, 100)
            .err()
            .map(|e| e.to_string())
            .expect("the long line is refused");
        assert_eq!(
            error,
            format!(
                "plugin execution failed: {}: the line at block offset 0 byte 35 is longer than \
                 the {MAX_LINE_BYTES}-byte limit for one line",
                data.display()
            )
        );

        let mut header_only = Vec::new();
        header_only.extend(std::iter::repeat_n(b'#', MAX_LINE_BYTES + 1));
        header_only.push(b'\n');
        write_bgzf(&data, &header_only);
        let error = TabixAnnotator::open(config(data.clone(), Some(2)))
            .err()
            .map(|e| e.to_string())
            .expect("the long header line is refused");
        assert_eq!(
            error,
            format!(
                "plugin initialization failed: {}: line 1 is longer than the \
                 {MAX_LINE_BYTES}-byte limit for one line",
                data.display()
            )
        );
    }

    /// `new` strips `\n` or `\r\n` as `BufRead::lines` does; `chomped` strips
    /// the `\n` alone, as Perl's `chomp` does, and keeps the `\r`.
    #[test]
    fn bounded_lines_strip_the_line_end_or_the_newline_alone() {
        let path = Path::new("scores.txt");
        let text = b"a\tb\r\nc\r\nd\n\r\ne";
        let lines: Vec<Vec<u8>> = BoundedLines::new(&text[..], path)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(lines, [&b"a\tb"[..], b"c", b"d", b"", b"e"]);
        let lines: Vec<Vec<u8>> = BoundedLines::chomped(&text[..], path)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(lines, [&b"a\tb\r"[..], b"c\r", b"d", b"\r", b"e"]);
    }

    /// A line holding a byte that is not UTF-8 (Latin-1 `é`, 0xE9) is handed
    /// over as its bytes, as `<$fh>` hands it to an Ensembl VEP plugin, and the text
    /// made of it keeps every other byte, writing the one byte as U+FFFD;
    /// valid UTF-8 is kept as it stands.
    #[test]
    fn bounded_lines_hand_a_line_that_is_not_utf8_as_its_bytes() {
        let path = Path::new("scores.txt");
        let text = b"# caf\xE9\nSOD1\t0.44\xE9\tcaf\xC3\xA9\n";
        let lines: Vec<Vec<u8>> = BoundedLines::new(&text[..], path)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(lines, [&b"# caf\xE9"[..], b"SOD1\t0.44\xE9\tcaf\xC3\xA9"]);
        let mut reused = BoundedLines::new(&text[..], path);
        assert_eq!(reused.next_line().unwrap().unwrap(), lines[0]);
        assert_eq!(reused.next_line().unwrap().unwrap(), lines[1]);
        assert!(reused.next_line().is_none());
        assert_eq!(lossy_text(lines[0].clone()), "# caf\u{FFFD}");
        assert_eq!(lossy_text(lines[1].clone()), "SOD1\t0.44\u{FFFD}\tcafé");
        assert_eq!(lossy_text(b"caf\xC3\xA9".to_vec()), "café");
    }

    /// The header readers read a bgzipped file whose preamble holds a byte that
    /// is not UTF-8 as `tabix -H` prints it: the last preamble line, the column
    /// names and the first data row's field are all read, the byte written as
    /// U+FFFD in the text they return.
    #[test]
    fn header_readers_read_a_preamble_that_is_not_utf8() {
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("synthetic_scores.tsv.gz");
        write_bgzf(
            &data,
            b"##CADD GRCh38-v1.7 made-up fixtur\xE9\n#Chrom\tPos\tRef\tAlt\tcaf\xE9\n21\t5\tA\tC\t0.5\xE9\n",
        );
        assert_eq!(
            last_preamble_line(&data).unwrap().as_deref(),
            Some("#Chrom\tPos\tRef\tAlt\tcaf\u{FFFD}")
        );
        assert_eq!(
            parse_header(&data, 1).unwrap().unwrap(),
            ["Chrom", "Pos", "Ref", "Alt", "caf\u{FFFD}"]
        );
        assert_eq!(
            first_data_row_field(&data, 4).unwrap().as_deref(),
            Some("0.5\u{FFFD}")
        );
    }

    /// The bare-header test reads the start cell as the record reader does, to
    /// its leading integer: a table with no header whose first row's start cell
    /// carries a fraction or a suffix has no header, and the reader returns the
    /// row at the integer `tabix` indexed it by, where a bare header row, whose
    /// start cell opens with no integer, is the header.
    #[test]
    fn a_bare_header_is_told_from_a_data_row_by_its_leading_integer() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let headerless = write_tabix_fixture(
            dir.path(),
            "headerless.tsv.gz",
            &["21\t100.5\tA\tC\t0.5", "21\t200x\tA\tG\t0.7"],
            IndexSpec::standard(),
        );
        assert_eq!(parse_header(&headerless, 1).unwrap(), None);
        let ann = open(headerless, None);
        assert!(ann.header().is_none());
        assert_eq!(
            lines(&ann.query("21", 100, 200).unwrap()),
            ["21\t100.5\tA\tC\t0.5", "21\t200x\tA\tG\t0.7"]
        );

        let bare = write_tabix_fixture(
            dir.path(),
            "bare.tsv.gz",
            &["chr\tpos\tref\talt\tscore", "21\t100.5\tA\tC\t0.5"],
            IndexSpec::standard().bare_header(),
        );
        assert_eq!(
            parse_header(&bare, 1).unwrap().unwrap(),
            ["chr", "pos", "ref", "alt", "score"]
        );
        let ann = open(bare, None);
        assert_eq!(
            ann.query("21", 100, 100).unwrap()[0].get("score"),
            Some("0.5")
        );
    }

    /// A headerless file whose first row's start cell a minus sign leads has
    /// no header under [`TabixAnnotator::reading_rows_as_tabix_indexes_them`],
    /// the reader that takes the row for a record, where a plugin's table
    /// reads the line as a bare header and drops the row. A bare header whose
    /// start cell leads no number, and a `#` header whatever that cell holds,
    /// stay the header for that reader too.
    #[test]
    fn a_headerless_files_minus_led_first_row_is_a_row() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let first_row = ["21", "-5", "300", "below"];
        let headerless = write_tabix_fixture(
            dir.path(),
            "headerless.tsv.gz",
            &["21\t-5\t300\tbelow", "21\t100\t100\tp100"],
            IndexSpec::standard().with_span_cols(1, 2, 3),
        );
        assert_eq!(parse_header(&headerless, 1).unwrap().unwrap(), first_row);
        let table = open(headerless.clone(), Some(2));
        assert_eq!(table.header().unwrap(), first_row);
        assert!(table.query("21", 1, 1).unwrap().is_empty());

        let as_tabix = open(headerless, Some(2)).reading_rows_as_tabix_indexes_them();
        assert!(as_tabix.header().is_none());
        let at_first_base = as_tabix.query("21", 1, 1).unwrap();
        assert_eq!(lines(&at_first_base), ["21\t-5\t300\tbelow"]);
        assert!(at_first_base[0].header.is_none());

        let bare = write_tabix_fixture(
            dir.path(),
            "bare.tsv.gz",
            &["chr\tstart\tend\tname", "21\t-5\t300\tbelow"],
            IndexSpec::standard().with_span_cols(1, 2, 3).bare_header(),
        );
        let as_tabix = open(bare, Some(2)).reading_rows_as_tabix_indexes_them();
        assert_eq!(as_tabix.header().unwrap(), ["chr", "start", "end", "name"]);
        assert_eq!(
            as_tabix.query("21", 1, 1).unwrap()[0].get("name"),
            Some("below")
        );

        let hashed = write_tabix_fixture(
            dir.path(),
            "hashed.tsv.gz",
            &["#chr\t-1\tend\tname", "21\t-5\t300\tbelow"],
            IndexSpec::standard().with_span_cols(1, 2, 3),
        );
        let as_tabix = open(hashed, Some(2)).reading_rows_as_tabix_indexes_them();
        assert_eq!(as_tabix.header().unwrap(), ["chr", "-1", "end", "name"]);
        assert_eq!(
            as_tabix.query("21", 1, 1).unwrap()[0].get("name"),
            Some("below")
        );
    }

    /// A line of exactly the cap is read whole, whether a line end or the end
    /// of the input closes it, and the line after it follows; one byte more is
    /// refused by its line number.
    #[test]
    fn bounded_lines_admit_a_line_of_exactly_the_cap() {
        let path = Path::new("catalogue.tsv");
        let mut text: Vec<u8> = std::iter::repeat_n(b'x', MAX_LINE_BYTES).collect();
        text.extend(b"\ny\n");
        let lines: Vec<Vec<u8>> = BoundedLines::new(&text[..], path)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0].len(), MAX_LINE_BYTES);
        assert_eq!(lines[1], b"y");

        let lines: Vec<Vec<u8>> = BoundedLines::new(&text[..MAX_LINE_BYTES], path)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0].len(), MAX_LINE_BYTES);

        let mut text: Vec<u8> = std::iter::repeat_n(b'x', MAX_LINE_BYTES + 1).collect();
        text.push(b'\n');
        let error = BoundedLines::new(&text[..], path)
            .next()
            .unwrap()
            .unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "plugin initialization failed: catalogue.tsv: line 1 is longer than the \
                 {MAX_LINE_BYTES}-byte limit for one line"
            )
        );
    }

    /// The one-pass sweep over a batch's windows answers exactly as the
    /// per-window query does: windows of one bin, a window straddling the
    /// 16 kb bin boundary and the bgzf block boundary, windows holding no
    /// record, and a window past the last record.
    #[test]
    fn query_windows_equals_the_per_window_query() {
        let dir = tempfile::tempdir().unwrap();
        let blocks = fixture_with_index(dir.path(), "blocks.tsv.gz", "csi");
        let ann = open(blocks, None);
        // Records sit at 5010000 + 7i; the 16 kb bin boundary at 5029888 falls
        // between the rows at 5029887 and 5029894, and the bgzf block boundary
        // falls at the row at 5031749.
        let mut windows: Vec<(u64, u64)> = Vec::new();
        for i in [
            0u64, 1, 2, 5, 6, 7, 100, 2840, 2841, 2842, 3105, 3106, 3107, 3399,
        ] {
            let start = 5_010_000 + 7 * i;
            windows.push((start - 1, start + 1));
        }
        windows.push((5_010_002, 5_010_005)); // between two records
        windows.push((5_029_885, 5_029_895)); // across the bin boundary
        windows.push((5_031_745, 5_031_752)); // across the block boundary
        windows.push((5_033_800, 5_033_900)); // past the last record
        windows.sort_unstable();
        // The sweep takes disjoint windows, as a caller merges them.
        let mut disjoint: Vec<(u64, u64)> = Vec::new();
        for (start, end) in windows {
            match disjoint.last_mut() {
                Some(last) if start <= last.1 => last.1 = last.1.max(end),
                _ => disjoint.push((start, end)),
            }
        }
        let windows = disjoint;
        let mut by_query: Vec<String> = Vec::new();
        for &(start, end) in &windows {
            by_query.extend(lines(&ann.query("21", start, end).unwrap()));
        }
        let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
        assert_eq!(by_sweep, by_query);
        assert!(by_sweep.iter().any(|l| l.starts_with("21\t5029887\t")));
        assert!(by_sweep.iter().any(|l| l.starts_with("21\t5029894\t")));
        assert!(by_sweep.iter().any(|l| l.starts_with("21\t5031749\t")));
        assert_eq!(
            by_sweep
                .iter()
                .filter(|l| l.starts_with("21\t5010000\t"))
                .count(),
            1
        );
        assert!(ann.query_windows("21", &[]).unwrap().is_empty());
        assert!(ann.query_windows("22", &windows).unwrap().is_empty());
    }

    /// The block-parallel read returns the rows `query_windows` returns, in its
    /// order: windows of one bin, a window across the 16 kb bin boundary, one
    /// across the bgzf block boundary (the row at 5031749 straddles it), windows
    /// holding no record, a window past the last record, single-base windows at
    /// every row, no window, and a chromosome the file lacks.
    #[test]
    fn query_windows_parallel_equals_query_windows() {
        let dir = tempfile::tempdir().unwrap();
        let blocks = fixture_with_index(dir.path(), "blocks.tsv.gz", "csi");
        let ann = open(blocks, None);
        let mut windows: Vec<(u64, u64)> = Vec::new();
        for i in [
            0u64, 1, 2, 5, 6, 7, 100, 2840, 2841, 2842, 3105, 3106, 3107, 3399,
        ] {
            let start = 5_010_000 + 7 * i;
            windows.push((start - 1, start + 1));
        }
        windows.push((5_010_002, 5_010_005));
        windows.push((5_029_885, 5_029_895));
        windows.push((5_031_745, 5_031_752));
        windows.push((5_033_800, 5_033_900));
        windows.sort_unstable();
        let mut disjoint: Vec<(u64, u64)> = Vec::new();
        for (start, end) in windows {
            match disjoint.last_mut() {
                Some(last) if start <= last.1 => last.1 = last.1.max(end),
                _ => disjoint.push((start, end)),
            }
        }
        let by_sweep = lines(&ann.query_windows("21", &disjoint).unwrap());
        let by_blocks = lines(&ann.query_windows_parallel("21", &disjoint).unwrap());
        assert_eq!(by_blocks, by_sweep);
        assert!(by_blocks.iter().any(|l| l.starts_with("21\t5031749\t")));

        let every_row: Vec<(u64, u64)> = (0..3400u64)
            .map(|i| (5_010_000 + 7 * i, 5_010_000 + 7 * i))
            .collect();
        let all = lines(&ann.query_windows_parallel("21", &every_row).unwrap());
        assert_eq!(all.len(), 3400);
        assert_eq!(all, lines(&ann.query_windows("21", &every_row).unwrap()));
        assert_eq!(all[0], "21\t5010000\tA\tC\t0.000");
        assert_eq!(all[3399], "21\t5033793\tT\tA\t0.681");

        let one = lines(
            &ann.query_windows_parallel("21", &[(5_031_749, 5_031_749)])
                .unwrap(),
        );
        assert_eq!(one.len(), 1);
        assert!(one[0].starts_with("21\t5031749\t"));
        assert!(ann.query_windows_parallel("21", &[]).unwrap().is_empty());
        assert!(ann
            .query_windows_parallel("22", &disjoint)
            .unwrap()
            .is_empty());
    }

    /// A line read whole across a block boundary keeps its line end, so a
    /// start column that is the file's last column is read to its digits with
    /// the line end left unread, as `line_start` reads it; a byte that is not
    /// UTF-8 after the digits is left unread the same way, so the two readers
    /// place such a line alike.
    #[test]
    fn line_start_bytes_reads_a_last_column_with_its_line_end() {
        let ann = open(fixture("blocks.tsv.gz"), None);
        assert_eq!(
            ann.line_start_bytes(b"21\t5031749\tA\tC\t0.5"),
            Some(5_031_749)
        );
        assert_eq!(ann.line_start_bytes(b"21\t5031749\n"), Some(5_031_749));
        assert_eq!(ann.line_start_bytes(b"21\t5031749\r\n"), Some(5_031_749));
        assert_eq!(ann.line_start_bytes(b"21\t5031749 \n"), Some(5_031_749));
        assert_eq!(ann.line_start_bytes(b"#chr\tpos"), None);
        assert_eq!(ann.line_start_bytes(b"21"), None);
        let latin1 = b"21\t5031749\xE9\tA\n";
        assert_eq!(ann.line_start_bytes(latin1), Some(5_031_749));
        assert_eq!(
            ann.line_start(&String::from_utf8_lossy(latin1)),
            ann.line_start_bytes(latin1)
        );
        assert_eq!(ann.line_start_bytes(b"21\t\xE95031749\tA\n"), None);
    }

    /// An index that does not describe the data file is refused by the
    /// block-parallel read as by the line-by-line one: a chunk whose block lies
    /// past the end of the file; a longer file's index over a shorter file,
    /// whose chunk names a block at or past the shorter file's end; a chunk
    /// whose block exists but holds fewer bytes than the chunk's offset into it;
    /// and a chunk that reaches past the file's end, refused from the file's
    /// length before any block of it is read.
    #[test]
    fn query_windows_parallel_refuses_an_index_the_file_does_not_describe() {
        let dir = tempfile::tempdir().unwrap();
        let past_the_end = 1_000_000;
        let data = write_multi_chunk_fixture(
            dir.path(),
            &["#chrom\tstart\tend\tname", "21\t100\t100\tp100"],
            &[(0..2, past_the_end)],
        );
        let ann = open(data.clone(), None);
        let err = ann
            .query_windows_parallel("21", &[(100, 100)])
            .expect_err("a chunk past the end of the file must be refused")
            .to_string();
        let path = data.display();
        assert_eq!(
            err,
            format!(
                "plugin execution failed: index {path}.tbi does not describe {path}: no block \
                 at offset {past_the_end}"
            )
        );

        let other = tempfile::tempdir().unwrap();
        let data = other.path().join("blocks.tsv.gz");
        std::fs::copy(fixture("spans.tsv.gz"), &data).unwrap();
        std::fs::copy(
            fixture("blocks.tsv.gz.tbi"),
            other.path().join("blocks.tsv.gz.tbi"),
        )
        .unwrap();
        let ann = open(data.clone(), None);
        let err = ann
            .query_windows_parallel("21", &[(5_033_779, 5_033_793)])
            .expect_err("a block the shorter file lacks must be refused")
            .to_string();
        let path = data.display();
        let prefix = format!(
            "plugin execution failed: index {path}.tbi does not describe {path}: no block \
             at offset "
        );
        let named_offset: u64 = err
            .strip_prefix(&prefix)
            .unwrap_or_else(|| panic!("{err}"))
            .parse()
            .unwrap();
        assert!(
            named_offset >= std::fs::metadata(&data).unwrap().len(),
            "{err}"
        );

        let third = tempfile::tempdir().unwrap();
        let data = third.path().join("chunks.tsv.gz");
        write_bgzf(&data, b"#chrom\tstart\tend\tname\n21\t100\t100\tp100\n");
        let position = |offset: u16| bgzf::VirtualPosition::try_from((0u64, offset)).unwrap();
        write_tbi(&data, vec![Chunk::new(position(60_000), position(60_001))]);
        let ann = open(data.clone(), None);
        let err = ann
            .query_windows_parallel("21", &[(100, 100)])
            .expect_err("an offset past the block's text must be refused")
            .to_string();
        let path = data.display();
        assert_eq!(
            err,
            format!(
                "plugin execution failed: index {path}.tbi does not describe {path}: no block \
                 at offset 0 holds byte 60000"
            )
        );

        let span_to = |block: u64| bgzf::VirtualPosition::try_from((block, 0u16)).unwrap();
        let past = 1u64 << 20;
        write_tbi(&data, vec![Chunk::new(position(0), span_to(past))]);
        let ann = open(data.clone(), None);
        let err = ann
            .query_windows_parallel("21", &[(100, 100)])
            .expect_err("a span past the end of the file must be refused")
            .to_string();
        assert_eq!(
            err,
            format!(
                "plugin execution failed: index {path}.tbi does not describe {path}: no block \
                 at offset {past}"
            )
        );
    }

    /// A block past the span's end whose header is sound but whose body does
    /// not inflate sends the span to the line-by-line read, which never enters
    /// that block and answers the span's own rows.
    #[test]
    fn query_windows_parallel_defers_to_the_sweep_when_a_block_does_not_inflate() {
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("damaged.tsv.gz");
        let first = b"#chrom\tpos\tref\talt\tscore\n21\t100\tA\tC\t0.5\n21\t100\tA\tG\t0.6\n";
        let second = b"21\t40000\tA\tC\t0.7\n21\t40001\tA\tC\t0.8\n";
        let offsets = write_blocks(&data, &[first, second]);
        let mut bytes = std::fs::read(&data).unwrap();
        // The deflate body follows the header's fixed bytes and its six-byte
        // `BC` subfield.
        let body = usize::try_from(offsets[1]).unwrap() + BGZF_HEADER_BYTES + 6 + 3;
        bytes[body] ^= 0xff;
        bytes[body + 1] ^= 0xff;
        std::fs::write(&data, &bytes).unwrap();
        let header_len = first.iter().position(|&b| b == b'\n').unwrap() + 1;
        let position = |offset: usize| {
            bgzf::VirtualPosition::try_from((0u64, u16::try_from(offset).unwrap())).unwrap()
        };
        write_tbi(
            &data,
            vec![Chunk::new(position(header_len), position(first.len()))],
        );
        let ann = open(data, None);
        let damaged =
            &bytes[usize::try_from(offsets[1]).unwrap()..usize::try_from(offsets[2]).unwrap()];
        assert!(ann.inflate_block(damaged).is_err());
        let windows = [(100, 100)];
        let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
        assert_eq!(by_sweep, ["21\t100\tA\tC\t0.5", "21\t100\tA\tG\t0.6"]);
        assert_eq!(
            lines(&ann.query_windows_parallel("21", &windows).unwrap()),
            by_sweep
        );
    }

    /// A span read in slices of one block answers as the one-slice read and the
    /// line-by-line read do: every line that crosses a slice boundary is made
    /// whole from the block that opens the next slice and scanned once.
    #[test]
    fn query_windows_parallel_reads_a_span_in_slices_of_one_block() {
        let dir = tempfile::tempdir().unwrap();
        let ann = open(fixture_with_index(dir.path(), "blocks.tsv.gz", "tbi"), None);
        let every_row: Vec<(u64, u64)> = (0..3400)
            .map(|i| (5_010_000 + 7 * i, 5_010_000 + 7 * i))
            .collect();
        for windows in [
            vec![
                (5_010_000, 5_010_007),
                (5_014_000, 5_014_010),
                (5_033_790, 5_033_793),
            ],
            every_row,
        ] {
            let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
            assert!(by_sweep.len() >= 4, "{by_sweep:?}");
            assert_eq!(
                lines(&ann.query_windows_in_slices("21", &windows, 1).unwrap()),
                by_sweep
            );
            assert_eq!(
                lines(&ann.query_windows_parallel("21", &windows).unwrap()),
                by_sweep
            );
        }
    }

    /// An empty block between two blocks (the end marker of a concatenated
    /// file) ends no line: read in slices of one block, the slice that opens
    /// with the block after it starts at that block's first byte and keeps the
    /// line there.
    #[test]
    fn query_windows_parallel_keeps_the_line_after_an_empty_block() {
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("concatenated.tsv.gz");
        let offsets = write_blocks(&data, &[b"21\t100\tX\n", b"21\t200\tX\n21\t300\tX\n"]);
        let bytes = std::fs::read(&data).unwrap();
        let eof_block = &bytes[usize::try_from(offsets[2]).unwrap()..];
        let second = usize::try_from(offsets[1]).unwrap();
        let mut concatenated = bytes[..second].to_vec();
        concatenated.extend_from_slice(eof_block);
        concatenated.extend_from_slice(&bytes[second..]);
        std::fs::write(&data, &concatenated).unwrap();
        let position = |block: usize| {
            bgzf::VirtualPosition::try_from((u64::try_from(block).unwrap(), 0u16)).unwrap()
        };
        let eof_at = concatenated.len() - eof_block.len();
        write_tbi(&data, vec![Chunk::new(position(0), position(eof_at))]);
        let ann = open(data, None);
        let windows = [(100, 300)];
        let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
        assert_eq!(by_sweep, ["21\t100\tX", "21\t200\tX", "21\t300\tX"]);
        assert_eq!(
            lines(&ann.query_windows_in_slices("21", &windows, 1).unwrap()),
            by_sweep
        );
        assert_eq!(
            lines(&ann.query_windows_parallel("21", &windows).unwrap()),
            by_sweep
        );
    }

    /// Write `first` and `second` as one bgzf block each with `empty_blocks`
    /// empty blocks between them, the shape of a file joined from compressed
    /// streams, then the EOF block; the compressed offsets of the first empty
    /// block, of `second`'s block and of the EOF block.
    fn write_joined_streams(
        path: &Path,
        first: &[u8],
        second: &[u8],
        empty_blocks: usize,
    ) -> (u64, u64, u64) {
        let offsets = write_blocks(path, &[first, second]);
        let bytes = std::fs::read(path).unwrap();
        let eof_block = &bytes[usize::try_from(offsets[2]).unwrap()..];
        let second_at = usize::try_from(offsets[1]).unwrap();
        let mut joined = bytes[..second_at].to_vec();
        for _ in 0..empty_blocks {
            joined.extend_from_slice(eof_block);
        }
        joined.extend_from_slice(&bytes[second_at..]);
        std::fs::write(path, &joined).unwrap();
        let moved = u64::try_from(empty_blocks * eof_block.len()).unwrap();
        (offsets[1], offsets[1] + moved, offsets[2] + moved)
    }

    /// `tabix` gives the record that opens the block after an empty one the
    /// empty block's offset as its chunk start, and reads on from the next
    /// block that holds bytes (htslib `bgzf_read_block`). A chunk that starts
    /// at one empty block, or at two, is read from the block after them, by
    /// the line-by-line reads and the block-by-block one alike. The two
    /// blocks' rows lie in two 16 kb tiles, each with its own bin and chunk,
    /// so a query inside the second tile enters the file at the empty block.
    #[test]
    fn a_chunk_that_starts_at_an_empty_block_is_read_from_the_next_block() {
        let position = |block: u64| bgzf::VirtualPosition::try_from((block, 0u16)).unwrap();
        for empty_blocks in [1, 2] {
            let dir = tempfile::tempdir().unwrap();
            let data = dir.path().join("joined.tsv.gz");
            let (empty_at, _, eof_at) = write_joined_streams(
                &data,
                b"21\t100\tX\n",
                b"21\t20000\tX\n21\t20100\tX\n",
                empty_blocks,
            );
            write_tbi_bins(
                &data,
                vec![
                    (
                        tile_bin(100),
                        vec![Chunk::new(position(0), position(empty_at))],
                    ),
                    (
                        tile_bin(20_000),
                        vec![Chunk::new(position(empty_at), position(eof_at))],
                    ),
                ],
                vec![position(0), position(empty_at)],
            );
            let ann = open(data, None);
            let after = ["21\t20000\tX", "21\t20100\tX"];
            assert_eq!(lines(&ann.query("21", 20_000, 20_100).unwrap()), after);
            assert_eq!(
                lines(&ann.query("21", 100, 20_100).unwrap()),
                ["21\t100\tX", "21\t20000\tX", "21\t20100\tX"]
            );
            let windows = [(20_000, 20_100)];
            assert_eq!(lines(&ann.query_windows("21", &windows).unwrap()), after);
            assert_eq!(
                lines(&ann.query_windows_parallel("21", &windows).unwrap()),
                after
            );
            let batch = ann.query_batch(&[("21", 20_100, 20_100)]).unwrap();
            assert_eq!(lines(batch.lookup("21", 20_100, 20_100)), ["21\t20100\tX"]);
        }
    }

    /// An empty block is taken as a chunk's start only for the block after it:
    /// a chunk at the file's last block, empty with nothing after it, and an
    /// offset past the first byte of an empty block, which holds none, name
    /// no record and are refused.
    #[test]
    fn a_chunk_at_an_empty_block_that_opens_no_record_is_an_error() {
        let position =
            |block: u64, offset: u16| bgzf::VirtualPosition::try_from((block, offset)).unwrap();
        let refused = |chunk_start: bgzf::VirtualPosition| {
            let dir = tempfile::tempdir().unwrap();
            let data = dir.path().join("joined.tsv.gz");
            write_joined_streams(&data, b"21\t100\tX\n", b"21\t200\tX\n21\t300\tX\n", 1);
            let file_end = position(std::fs::metadata(&data).unwrap().len(), 0);
            write_tbi(&data, vec![Chunk::new(chunk_start, file_end)]);
            let err = match open(data.clone(), None).query("21", 200, 300) {
                Ok(_) => panic!("a chunk that opens no record must be refused"),
                Err(e) => e.to_string(),
            };
            let path = data.display();
            let (block, offset): (u64, u16) = chunk_start.into();
            assert_eq!(
                err,
                format!(
                    "plugin execution failed: index {path}.tbi does not describe {path}: no \
                     block at offset {block} holds byte {offset}"
                )
            );
        };
        let dir = tempfile::tempdir().unwrap();
        let (empty_at, _, eof_at) = write_joined_streams(
            &dir.path().join("offsets.tsv.gz"),
            b"21\t100\tX\n",
            b"21\t200\tX\n21\t300\tX\n",
            1,
        );
        refused(position(eof_at, 0));
        refused(position(empty_at, 1));
    }

    /// Two clusters more than `GAP_BRIDGE` apart inside one index tile name the
    /// tile's chunk twice; the chunks merge into one span, read once, and the
    /// records are the line-by-line read's.
    #[test]
    fn query_windows_parallel_equals_query_windows_on_two_clusters_in_one_tile() {
        let dir = tempfile::tempdir().unwrap();
        let ann = open(fixture_with_index(dir.path(), "blocks.tsv.gz", "tbi"), None);
        let windows = [
            (5_010_000, 5_010_007),
            (5_014_000, 5_014_010),
            (5_020_000, 5_020_010),
            (5_033_790, 5_033_793),
        ];
        let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
        assert!(by_sweep.len() >= 5, "{by_sweep:?}");
        assert_eq!(
            lines(&ann.query_windows_parallel("21", &windows).unwrap()),
            by_sweep
        );
    }

    /// Chunks that share a block merge into one span (two clusters in one tile
    /// name the same chunk; a bin's chunk may end in the block the next bin's
    /// begins in); chunks in different blocks stay apart.
    #[test]
    fn merge_chunks_joins_chunks_that_share_a_block() {
        let position =
            |block: u64, offset: u16| bgzf::VirtualPosition::try_from((block, offset)).unwrap();
        let chunks = [
            Chunk::new(position(700, 0), position(900, 40)),
            Chunk::new(position(0, 10), position(300, 0)),
            Chunk::new(position(0, 10), position(300, 0)),
            Chunk::new(position(300, 0), position(500, 20)),
        ];
        assert_eq!(
            merge_chunks(&chunks),
            [
                (position(0, 10), position(500, 20)),
                (position(700, 0), position(900, 40)),
            ]
        );
    }

    /// An index built on a span (`-e` a column of its own) answers through
    /// `query_windows`, so the two agree on `spans.tsv.gz` too.
    #[test]
    fn query_windows_parallel_defers_to_query_windows_on_a_span_index() {
        let dir = tempfile::tempdir().unwrap();
        let ann = open(
            fixture_with_index(dir.path(), "spans.tsv.gz", "tbi"),
            Some(2),
        );
        let windows = [(150, 160), (200, 210), (300, 300)];
        assert_eq!(
            lines(&ann.query_windows_parallel("21", &windows).unwrap()),
            lines(&ann.query_windows("21", &windows).unwrap())
        );
    }

    /// `bgzf_block_size` reads the `BC` subfield of a block header and refuses
    /// bytes that are not one; `block_line_runs` gives each block the lines
    /// that start in it, whole across the block boundary, drops a block with no
    /// line end into its predecessor's run, keeps a final line with no line
    /// end only when the file ends with it, reporting it incomplete otherwise,
    /// and says whether the blocks scanned end inside a line, an empty block
    /// leaving that as it was.
    #[test]
    fn block_helpers_cut_at_boundaries_and_keep_lines_whole() {
        let mut header = vec![0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0];
        header.extend_from_slice(&[b'B', b'C', 2, 0, 0x0f, 0x01]);
        assert_eq!(bgzf_block_size(&header), Some(0x010f + 1));
        assert_eq!(bgzf_block_size(&header[..11]), None);
        assert_eq!(bgzf_block_size(b"not a block header at all"), None);
        let mut other = header.clone();
        other[12] = b'X';
        assert_eq!(bgzf_block_size(&other), None);

        let blocks: Vec<Vec<u8>> = vec![
            b"skip\nfirst\nsec".to_vec(),
            b"ond\nthird-starts-here-and-goes-on".to_vec(),
            b"-and-on".to_vec(),
            b"-done\nfourth\nfif".to_vec(),
            b"th".to_vec(),
        ];
        let run_texts = |runs: &[LineRun<'_>]| -> Vec<Vec<u8>> {
            runs.iter()
                .map(|run| {
                    let mut all = run.body.to_vec();
                    all.extend(run.straddle.iter().flatten());
                    all
                })
                .collect()
        };
        assert_eq!(
            run_texts(&block_line_runs(&blocks, 5, false, false, false).0),
            vec![
                b"first\nsecond\n".to_vec(),
                b"third-starts-here-and-goes-on-and-on-done\n".to_vec(),
                b"fourth\n".to_vec(),
            ]
        );
        let two = [b"a\nb\n".to_vec(), b"c\n".to_vec()];
        assert_eq!(
            run_texts(&block_line_runs(&two, 0, false, false, false).0),
            vec![b"a\nb\n".to_vec(), b"c\n".to_vec()]
        );
        assert_eq!(
            run_texts(&block_line_runs(&two, 0, false, true, false).0),
            vec![b"a\nb\n".to_vec()]
        );
        assert_eq!(
            run_texts(
                &block_line_runs(
                    &[b"a\nb".to_vec(), b"c\nd\n".to_vec()],
                    0,
                    false,
                    true,
                    false
                )
                .0
            ),
            vec![b"a\nbc\n".to_vec()]
        );
        // A final line with no line end is the file's last line when the file
        // ends there, and an incomplete line otherwise.
        let tail = [b"a\ntail".to_vec()];
        let (runs, incomplete, open_line) = block_line_runs(&tail, 0, false, false, false);
        assert_eq!(run_texts(&runs), vec![b"a\n".to_vec()]);
        assert!(incomplete);
        assert!(open_line);
        let (runs, incomplete, open_line) = block_line_runs(&tail, 0, false, false, true);
        assert_eq!(run_texts(&runs), vec![b"a\ntail".to_vec()]);
        assert!(!incomplete);
        assert!(open_line);
        // The blocks scanned end inside a line when the last byte of the last
        // non-empty one is not a line end; the tail is not scanned.
        let open_after = |blocks: &[Vec<u8>], open_before: bool, last_is_tail: bool| {
            block_line_runs(blocks, 0, open_before, last_is_tail, false).2
        };
        assert!(!open_after(
            &[b"a\n".to_vec(), b"".to_vec(), b"c".to_vec()],
            false,
            true
        ));
        assert!(open_after(
            &[b"a\nb".to_vec(), b"".to_vec(), b"c\n".to_vec()],
            false,
            true
        ));
        assert!(open_after(&[b"".to_vec(), b"c\n".to_vec()], true, true));
        assert!(!open_after(&[b"".to_vec(), b"c\n".to_vec()], false, true));
        assert!(!open_after(
            &[b"a\nb".to_vec(), b"c\n".to_vec()],
            false,
            false
        ));
        assert!(block_line_runs(&[b"tail".to_vec()], 0, false, false, false)
            .0
            .is_empty());
        // A first block that continues a line of the block before it starts
        // its first run after that line's end, or starts none without one.
        assert_eq!(
            run_texts(&block_line_runs(&two, 0, true, false, false).0),
            vec![b"b\n".to_vec(), b"c\n".to_vec()]
        );
        assert_eq!(
            run_texts(
                &block_line_runs(
                    &[b"no line end".to_vec(), b"x\ny\n".to_vec()],
                    0,
                    true,
                    false,
                    false
                )
                .0
            ),
            vec![b"y\n".to_vec()]
        );
    }

    /// Write each of `blocks` as one bgzf block of `path`, then the EOF block;
    /// the compressed offset of each block, the EOF block's last.
    fn write_blocks(path: &Path, blocks: &[&[u8]]) -> Vec<u64> {
        let mut writer = bgzf::io::Writer::new(std::fs::File::create(path).unwrap());
        let mut offsets = Vec::new();
        for block in blocks {
            offsets.push(writer.virtual_position().compressed());
            writer.write_all(block).unwrap();
            writer.flush().unwrap();
        }
        offsets.push(writer.virtual_position().compressed());
        writer.try_finish().unwrap();
        offsets
    }

    /// The bin of a single-base record at `position`: the 16 kb tile's bin on
    /// the deepest level of a `.tbi`.
    fn tile_bin(position: u64) -> usize {
        4681 + usize::try_from((position - 1) >> 14).unwrap()
    }

    /// A file cut inside the block past a span's end block: the span's own
    /// bytes are whole, so its records are the line-by-line read's, which never
    /// enters that block.
    #[test]
    fn query_windows_parallel_reads_a_span_whose_tail_block_is_cut_short() {
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("cut.tsv.gz");
        let first = b"#chrom\tpos\tref\talt\tscore\n21\t100\tA\tC\t0.5\n21\t100\tA\tG\t0.6\n";
        let second = b"21\t40000\tA\tC\t0.7\n";
        let offsets = write_blocks(&data, &[first, second]);
        let bytes = std::fs::read(&data).unwrap();
        let cut_inside_second = usize::try_from(offsets[1]).unwrap() + 20;
        std::fs::write(&data, &bytes[..cut_inside_second]).unwrap();
        let header_len = first.iter().position(|&b| b == b'\n').unwrap() + 1;
        let position = |offset: usize| {
            bgzf::VirtualPosition::try_from((0u64, u16::try_from(offset).unwrap())).unwrap()
        };
        write_tbi(
            &data,
            vec![Chunk::new(position(header_len), position(first.len()))],
        );
        let ann = open(data, None);
        let windows = [(100, 100)];
        let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
        assert_eq!(by_sweep, ["21\t100\tA\tC\t0.5", "21\t100\tA\tG\t0.6"]);
        assert_eq!(
            lines(&ann.query_windows_parallel("21", &windows).unwrap()),
            by_sweep
        );
    }

    /// The block past a span's end block holds no whole bgzf header: the file
    /// is cut inside the header's fixed bytes, inside its extra field, or the
    /// block's first byte is not the gzip magic. The span's own bytes are
    /// whole, so its records are the line-by-line read's, which never enters
    /// that block.
    #[test]
    fn query_windows_parallel_reads_a_span_whose_tail_block_has_no_whole_header() {
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("cut.tsv.gz");
        let first = b"#chrom\tpos\tref\talt\tscore\n21\t100\tA\tC\t0.5\n21\t100\tA\tG\t0.6\n";
        let second = b"21\t40000\tA\tC\t0.7\n";
        let offsets = write_blocks(&data, &[first, second]);
        let bytes = std::fs::read(&data).unwrap();
        let tail = usize::try_from(offsets[1]).unwrap();
        let header_len = first.iter().position(|&b| b == b'\n').unwrap() + 1;
        let position = |offset: usize| {
            bgzf::VirtualPosition::try_from((0u64, u16::try_from(offset).unwrap())).unwrap()
        };
        write_tbi(
            &data,
            vec![Chunk::new(position(header_len), position(first.len()))],
        );
        let mut not_bgzf = bytes.clone();
        not_bgzf[tail] = b'~';
        for damaged in [
            &bytes[..tail + BGZF_HEADER_BYTES - 7],
            &bytes[..tail + BGZF_HEADER_BYTES + 2],
            &not_bgzf[..],
        ] {
            std::fs::write(&data, damaged).unwrap();
            let ann = open(data.clone(), None);
            let windows = [(100, 100)];
            let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
            assert_eq!(by_sweep, ["21\t100\tA\tC\t0.5", "21\t100\tA\tG\t0.6"]);
            assert_eq!(
                lines(&ann.query_windows_parallel("21", &windows).unwrap()),
                by_sweep,
                "{} bytes",
                damaged.len()
            );
        }
    }

    /// `blocks.tsv.gz` with its EOF block cut off, under its own `.tbi`, whose
    /// last chunk ends at the offset the EOF block had: the end of the file.
    /// Both reads answer the file's first and last rows.
    #[test]
    fn query_windows_parallel_reads_a_file_whose_eof_block_is_missing() {
        let dir = tempfile::tempdir().unwrap();
        let data = fixture_with_index(dir.path(), "blocks.tsv.gz", "tbi");
        let bytes = std::fs::read(&data).unwrap();
        std::fs::write(&data, &bytes[..bytes.len() - 28]).unwrap();
        let ann = open(data, None);
        let windows = [(5_010_000, 5_010_007), (5_033_779, 5_033_793)];
        let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
        assert_eq!(by_sweep.len(), 5);
        assert_eq!(by_sweep[4], "21\t5033793\tT\tA\t0.681");
        assert_eq!(
            lines(&ann.query_windows_parallel("21", &windows).unwrap()),
            by_sweep
        );
    }

    /// Two clusters more than `GAP_BRIDGE` apart read their own blocks alone:
    /// the block between them (the fourth, overwritten so that it is not bgzf)
    /// is never reached, where a read of the whole span from the first window
    /// to the last would refuse it. Each bin's chunk ends where the next block
    /// starts, as `tabix` records a chunk that ends with its block.
    #[test]
    fn query_windows_parallel_reads_the_blocks_of_each_cluster_alone() {
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("tiles.tsv.gz");
        let rows = |first: u64| -> Vec<u8> {
            (first..first + 5)
                .map(|position| format!("21\t{position}\tX\n"))
                .collect::<String>()
                .into_bytes()
        };
        let tiles: Vec<Vec<u8>> = [1, 2, 3, 4, 5]
            .iter()
            .map(|&tile| rows(16_384 * tile + 1))
            .collect();
        let offsets = write_blocks(&data, &tiles.iter().map(Vec::as_slice).collect::<Vec<_>>());
        {
            use std::io::{Seek, SeekFrom};
            let mut file = std::fs::OpenOptions::new().write(true).open(&data).unwrap();
            file.seek(SeekFrom::Start(offsets[3])).unwrap();
            let garbage = vec![b'~'; usize::try_from(offsets[4] - offsets[3]).unwrap()];
            file.write_all(&garbage).unwrap();
        }
        let position = |block: u64| bgzf::VirtualPosition::try_from((block, 0)).unwrap();
        let bins: Vec<(usize, Vec<Chunk>)> = [0usize, 1, 2, 4]
            .iter()
            .map(|&i| {
                let chunk = Chunk::new(position(offsets[i]), position(offsets[i + 1]));
                (tile_bin(16_384 * (i as u64 + 1) + 1), vec![chunk])
            })
            .collect();
        // Tile 0 holds no row and takes the first row's offset; tile 4's rows
        // were overwritten and it repeats tile 3's, as `tabix` fills a tile
        // without a row.
        let linear = [0usize, 0, 1, 2, 2, 4]
            .iter()
            .map(|&block| position(offsets[block]))
            .collect();
        write_tbi_bins(&data, bins, linear);
        let ann = open(data, None);

        let windows = [(16_387, 16_387), (16_389, 16_389), (81_922, 81_923)];
        let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
        assert_eq!(
            by_sweep,
            [
                "21\t16387\tX",
                "21\t16389\tX",
                "21\t81922\tX",
                "21\t81923\tX"
            ]
        );
        assert_eq!(
            lines(&ann.query_windows_parallel("21", &windows).unwrap()),
            by_sweep
        );
    }

    /// A line that starts in the block past a chunk's end block and runs
    /// through the one after it is left open by the blocks read for the chunk,
    /// so the cluster answers through `query_windows`; read as its own chunk,
    /// the line comes whole from the three blocks it spans.
    #[test]
    fn query_windows_parallel_answers_a_cluster_whose_blocks_leave_a_line_open() {
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("long.tsv.gz");
        let mut first = b"21\t1000\tX\n21\t1002\tX\n21\t1004\tX\n".to_vec();
        let long_line_at = first.len();
        first.extend_from_slice(b"21\t20000\t");
        first.extend(std::iter::repeat_n(b'X', 40));
        let middle = vec![b'X'; 60_000];
        let mut last = vec![b'X'; 100];
        last.push(b'\n');
        let after_long_line = last.len();
        last.extend_from_slice(b"21\t40000\tX\n21\t40002\tX\n");
        let offsets = write_blocks(&data, &[&first, &middle, &last]);
        let position = |block: u64, offset: usize| {
            bgzf::VirtualPosition::try_from((block, u16::try_from(offset).unwrap())).unwrap()
        };
        let bins = vec![
            (
                tile_bin(1000),
                vec![Chunk::new(
                    position(offsets[0], 0),
                    position(offsets[0], long_line_at),
                )],
            ),
            (
                tile_bin(20000),
                vec![Chunk::new(
                    position(offsets[0], long_line_at),
                    position(offsets[2], after_long_line),
                )],
            ),
            (
                tile_bin(40000),
                vec![Chunk::new(
                    position(offsets[2], after_long_line),
                    position(offsets[3], 0),
                )],
            ),
        ];
        let linear = vec![
            position(offsets[0], 0),
            position(offsets[0], long_line_at),
            position(offsets[2], after_long_line),
        ];
        write_tbi_bins(&data, bins, linear);
        let ann = open(data, None);

        for windows in [
            vec![(1000, 1000), (1004, 1004)],
            vec![(20000, 20000)],
            vec![(1002, 1002), (20000, 20000), (40002, 40002)],
        ] {
            let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
            assert_eq!(by_sweep.len(), windows.len(), "{windows:?}");
            assert_eq!(
                lines(&ann.query_windows_parallel("21", &windows).unwrap()),
                by_sweep,
                "{windows:?}"
            );
        }
        let long = ann.query_windows_parallel("21", &[(20000, 20000)]).unwrap();
        assert_eq!(long[0].columns[2].len(), 60_140);
    }

    /// A cluster whose span ends past the index's last addressable position
    /// (2^29 for a `.tbi`) still answers as the per-window query does: the
    /// index refuses the span, the window short of that position holds two
    /// records, and the window reaching past it is read to that position, as
    /// `tabix blocks.tsv.gz 21:5010500-600000000` prints the file's 3328
    /// records from 5010504 on.
    #[test]
    fn query_windows_answers_a_cluster_the_index_refuses_window_by_window() {
        let dir = tempfile::tempdir().unwrap();
        let blocks = fixture_with_index(dir.path(), "blocks.tsv.gz", "tbi");
        let ann = open(blocks, None);
        let windows = [(5_010_000, 5_010_010), (5_010_500, 600_000_000)];
        let mut by_query: Vec<String> = Vec::new();
        for &(start, end) in &windows {
            by_query.extend(lines(&ann.query("21", start, end).unwrap()));
        }
        assert_eq!(by_query.len(), 2 + 3328);
        assert_eq!(
            by_query[..3],
            [
                "21\t5010000\tA\tC\t0.000",
                "21\t5010007\tC\tG\t0.919",
                "21\t5010504\tA\tC\t0.168"
            ]
        );
        assert_eq!(by_query[by_query.len() - 1], "21\t5033793\tT\tA\t0.681");
        assert_eq!(
            lines(&ann.query("21", 5_010_500, 600_000_000).unwrap()),
            lines(&ann.query("21", 5_010_500, 1 << 29).unwrap())
        );
        assert_eq!(lines(&ann.query_windows("21", &windows).unwrap()), by_query);
    }

    /// A batch whose one wide span ends past the index's last addressable
    /// position merges with its neighbours' spans into one region; the region
    /// is read to that position, so the neighbours keep their records.
    #[test]
    fn query_batch_keeps_a_merged_region_one_span_reaches_past_the_index_limit() {
        let dir = tempfile::tempdir().unwrap();
        let blocks = fixture_with_index(dir.path(), "blocks.tsv.gz", "tbi");
        let ann = open(blocks, None);
        let batch = ann
            .query_batch(&[
                ("21", 5_010_000, 600_000_000),
                ("21", 5_010_007, 5_010_007),
                ("21", 5_033_793, 5_033_793),
            ])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 5_010_007, 5_010_007)),
            ["21\t5010007\tC\tG\t0.919"]
        );
        assert_eq!(
            lines(batch.lookup("21", 5_033_793, 5_033_793)),
            ["21\t5033793\tT\tA\t0.681"]
        );
    }

    /// `limit.tsv.gz` holds the rows `tabix` indexes at 2^29, the last
    /// position a `.tbi` holds: one ending there, one starting there and one
    /// ending a base short of it. Expected lines are the output of
    /// `tabix limit.tsv.gz <region>` for each region, which prints nothing
    /// for a region starting past that position.
    #[test]
    fn query_returns_the_rows_at_the_last_position_a_tbi_holds() {
        let dir = tempfile::tempdir().unwrap();
        let limit = fixture_with_index(dir.path(), "limit.tsv.gz", "tbi");
        let ann = open(limit, Some(2));
        let last: u64 = 1 << 29;
        let one_short = "21\t536869000\t536870911\tends_one_short";
        let at_last = [
            "21\t536870000\t536870912\tends_at_limit",
            "21\t536870912\t536870912\tstarts_at_limit",
        ];
        assert_eq!(lines(&ann.query("21", last, last).unwrap()), at_last);
        assert_eq!(
            lines(&ann.query("21", last - 1, last - 1).unwrap()),
            [one_short, at_last[0]]
        );
        assert_eq!(lines(&ann.query("21", last, 600_000_000).unwrap()), at_last);
        assert_eq!(
            lines(&ann.query("21", 1_500, 600_000_000).unwrap()),
            [
                "21\t1000\t2000\tordinary",
                one_short,
                at_last[0],
                at_last[1]
            ]
        );
        assert!(ann.query("21", last + 1, last + 1).unwrap().is_empty());
        let batch = ann
            .query_batch(&[("21", 1_500, 1_500), ("21", last, last)])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 1_500, 1_500)),
            ["21\t1000\t2000\tordinary"]
        );
        assert_eq!(lines(batch.lookup("21", last, last)), at_last);
        // Each region as a window: a span index answers a window as `query`
        // answers its region.
        assert_eq!(window_lines(&ann, &[(last, last)]), at_last);
        assert_eq!(
            window_lines(&ann, &[(last - 1, last - 1)]),
            [one_short, at_last[0]]
        );
        assert_eq!(window_lines(&ann, &[(last, 600_000_000)]), at_last);
        assert_eq!(
            window_lines(&ann, &[(1_500, 600_000_000)]),
            [
                "21\t1000\t2000\tordinary",
                one_short,
                at_last[0],
                at_last[1]
            ]
        );
        assert!(window_lines(&ann, &[(last + 1, last + 1)]).is_empty());
        assert_eq!(
            window_lines(&ann, &[(1_500, 1_500), (last, last)]),
            ["21\t1000\t2000\tordinary", at_last[0], at_last[1]]
        );
    }

    /// `limit_csi.tsv.gz`, written by `bgzip --binary` and indexed by
    /// `tabix -C -s 1 -b 2 -e 3`, holds the rows at 2^38, the last position
    /// that `.csi` holds (`min_shift` 14, depth 8): one ending there, one
    /// starting there and one ending a base short of it. Expected lines are
    /// the output of `tabix limit_csi.tsv.gz <region>` for each region, which
    /// prints nothing for a region starting past that position.
    #[test]
    fn query_returns_the_rows_at_the_last_position_a_csi_holds() {
        let dir = tempfile::tempdir().unwrap();
        let limit = fixture_with_index(dir.path(), "limit_csi.tsv.gz", "csi");
        let ann = open(limit, Some(2));
        assert_eq!((ann.index.min_shift(), ann.index.depth()), (14, 8));
        let last: u64 = 1 << 38;
        let one_short = "21\t274877906000\t274877906943\tends_one_short";
        let at_last = [
            "21\t274877906500\t274877906944\tends_at_limit",
            "21\t274877906944\t274877906944\tstarts_at_limit",
        ];
        assert_eq!(lines(&ann.query("21", last, last).unwrap()), at_last);
        assert_eq!(
            lines(&ann.query("21", last - 1, last - 1).unwrap()),
            [one_short, at_last[0]]
        );
        assert_eq!(
            lines(&ann.query("21", last, 300_000_000_000).unwrap()),
            at_last
        );
        assert_eq!(
            lines(&ann.query("21", 1_500, 300_000_000_000).unwrap()),
            [
                "21\t1000\t2000\tordinary",
                one_short,
                at_last[0],
                at_last[1]
            ]
        );
        assert!(ann.query("21", last + 1, last + 1).unwrap().is_empty());
        let batch = ann
            .query_batch(&[("21", 1_500, 1_500), ("21", last, last)])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 1_500, 1_500)),
            ["21\t1000\t2000\tordinary"]
        );
        assert_eq!(lines(batch.lookup("21", last, last)), at_last);
        // Each region as a window: a span index answers a window as `query`
        // answers its region.
        assert_eq!(window_lines(&ann, &[(last, last)]), at_last);
        assert_eq!(
            window_lines(&ann, &[(last - 1, last - 1)]),
            [one_short, at_last[0]]
        );
        assert_eq!(window_lines(&ann, &[(last, 300_000_000_000)]), at_last);
        assert_eq!(
            window_lines(&ann, &[(1_500, 300_000_000_000)]),
            [
                "21\t1000\t2000\tordinary",
                one_short,
                at_last[0],
                at_last[1]
            ]
        );
        assert!(window_lines(&ann, &[(last + 1, last + 1)]).is_empty());
        assert_eq!(
            window_lines(&ann, &[(1_500, 1_500), (last, last)]),
            ["21\t1000\t2000\tordinary", at_last[0], at_last[1]]
        );
    }

    /// `limit_csi_start.tsv.gz`, written by `bgzip --binary` and indexed by
    /// `tabix -C -s 1 -b 2 -e 2`, holds two rows at 2^38, the last position
    /// that `.csi` holds (`min_shift` 14, depth 8), one a base short of it and
    /// one at 1000. An index on the start column alone takes `query_windows`
    /// past the window-by-window answer of a span index, to the read of a
    /// batch in one pass. Expected lines are the output of
    /// `tabix limit_csi_start.tsv.gz <region>` for the region of each window,
    /// in window order. The index lookup takes no bound past 2^38 - 1, so
    /// `query` answers a cluster of windows that reaches 2^38 and the sweep a
    /// cluster short of it.
    #[test]
    fn query_windows_reads_the_last_position_a_csi_holds_in_one_pass() {
        let dir = tempfile::tempdir().unwrap();
        let limit = fixture_with_index(dir.path(), "limit_csi_start.tsv.gz", "csi");
        let ann = open(limit, None);
        assert_eq!((ann.index.min_shift(), ann.index.depth()), (14, 8));
        let last: u64 = 1 << 38;
        let one_short = "21\t274877906943\tone_short";
        let at_last = [
            "21\t274877906944\tat_limit",
            "21\t274877906944\talso_at_limit",
        ];
        let lookup = |position: u64| {
            let position = Position::try_from(position as usize).unwrap();
            ann.index.query(0, (position..=position).into())
        };
        assert!(lookup(last - 1).is_ok());
        assert!(lookup(last).is_err());

        assert_eq!(window_lines(&ann, &[(last, last)]), at_last);
        assert_eq!(window_lines(&ann, &[(last - 1, last - 1)]), [one_short]);
        assert_eq!(window_lines(&ann, &[(last, 300_000_000_000)]), at_last);
        assert_eq!(
            window_lines(&ann, &[(1_500, 300_000_000_000)]),
            [one_short, at_last[0], at_last[1]]
        );
        assert!(window_lines(&ann, &[(last + 1, last + 1)]).is_empty());
        // One batch of two clusters, the second reaching 2^38 and past it.
        assert_eq!(
            window_lines(
                &ann,
                &[
                    (1_000, 1_000),
                    (last - 1, last - 1),
                    (last, last),
                    (last + 1, last + 10)
                ]
            ),
            ["21\t1000\tordinary", one_short, at_last[0], at_last[1]]
        );
    }

    /// A query answered from two chunks returns every record of both, each
    /// once. A span record sits in a bin above the point records', so its chunk
    /// lies apart from theirs in the file, and a read buffered in bulk past its
    /// chunk's end takes in the second chunk's records before the seek to that
    /// chunk reads them again.
    #[test]
    fn records_from_a_second_index_chunk_are_returned() {
        let dir = tempfile::tempdir().unwrap();
        // Points every 100 bp from 100 to 47900, with one record spanning
        // 16000..40000 inserted in start order: 480 records, about 10 KB, so
        // the bytes from the span record to the end fit one 8 KB fill.
        let mut rows: Vec<String> = vec!["#chrom\tstart\tend\tname".to_string()];
        for start in (100..48000).step_by(100) {
            if start == 16000 {
                rows.push("21\t16000\t40000\tlong16000_40000".to_string());
            }
            rows.push(format!("21\t{start}\t{start}\tp{start}"));
        }
        let lines_in: Vec<&str> = rows.iter().map(String::as_str).collect();
        let long = lines_in.iter().position(|l| l.contains("long")).unwrap();
        let first_in_window = lines_in
            .iter()
            .position(|l| l.starts_with("21\t32800\t"))
            .unwrap();
        let data = write_multi_chunk_fixture(
            dir.path(),
            &lines_in,
            &[(long..long + 1, 0), (first_in_window..lines_in.len(), 0)],
        );
        let ann = open(data, Some(2));

        // The precondition the test rests on: the index answers with two chunks.
        let chunks = ann
            .index
            .query(
                0,
                (Position::new(39000).unwrap()..=Position::new(39000).unwrap()).into(),
            )
            .unwrap();
        assert_eq!(chunks.len(), 2);

        assert_eq!(
            lines(&ann.query("21", 39000, 39000).unwrap()),
            [
                "21\t16000\t40000\tlong16000_40000",
                "21\t39000\t39000\tp39000"
            ]
        );
        // The long record and the 152 points from 32800 to 47900, each once.
        let window = lines(&ann.query("21", 32800, 47900).unwrap());
        assert_eq!(window.len(), 153);
        assert_eq!(window[0], "21\t16000\t40000\tlong16000_40000");
        assert_eq!(window[1], "21\t32800\t32800\tp32800");
        let batch = ann
            .query_batch(&[("21", 39000, 39000), ("21", 40000, 40000)])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 39000, 39000)),
            [
                "21\t16000\t40000\tlong16000_40000",
                "21\t39000\t39000\tp39000"
            ]
        );
        assert_eq!(
            lines(batch.lookup("21", 40000, 40000)),
            [
                "21\t16000\t40000\tlong16000_40000",
                "21\t40000\t40000\tp40000"
            ]
        );
    }

    /// The row at 5031749 straddles the boundary between the fixture's two
    /// bgzf blocks, read here through the `.tbi`. Expected lines are the output
    /// of `tabix blocks.tsv.gz <region>` for each region, and
    /// `tabix blocks.tsv.gz 21 | wc -l` is 3400.
    #[test]
    fn bgzf_block_boundary_record_is_returned() {
        let dir = tempfile::tempdir().unwrap();
        let data = fixture_with_index(dir.path(), "blocks.tsv.gz", "tbi");
        assert_eq!(
            find_index_path(&data).unwrap(),
            dir.path().join("blocks.tsv.gz.tbi")
        );
        let ann = open(data, None);

        // tabix blocks.tsv.gz 21:5031749-5031749
        assert_eq!(
            lines(&ann.query("21", 5031749, 5031749).unwrap()),
            ["21\t5031749\tT\tG\t0.333"]
        );
        // tabix blocks.tsv.gz 21:5031735-5031763
        assert_eq!(
            lines(&ann.query("21", 5031735, 5031763).unwrap()),
            [
                "21\t5031735\tC\tA\t0.495",
                "21\t5031742\tG\tC\t0.414",
                "21\t5031749\tT\tG\t0.333",
                "21\t5031756\tA\tC\t0.252",
                "21\t5031763\tC\tG\t0.171",
            ]
        );
        // tabix blocks.tsv.gz 21:5031748-5031748 and 21:5031750-5031750 print nothing.
        assert!(ann.query("21", 5031748, 5031748).unwrap().is_empty());
        assert!(ann.query("21", 5031750, 5031750).unwrap().is_empty());
        // tabix blocks.tsv.gz 21:5010000-5010007
        assert_eq!(
            lines(&ann.query("21", 5010000, 5010007).unwrap()),
            ["21\t5010000\tA\tC\t0.000", "21\t5010007\tC\tG\t0.919"]
        );
        // tabix blocks.tsv.gz 21:5033779-5033793, entirely in the second block.
        assert_eq!(
            lines(&ann.query("21", 5033779, 5033793).unwrap()),
            [
                "21\t5033779\tC\tG\t0.843",
                "21\t5033786\tG\tT\t0.762",
                "21\t5033793\tT\tA\t0.681",
            ]
        );
        assert_eq!(ann.query("21", 5010000, 5033793).unwrap().len(), 3400);
        // `tabix blocks.tsv.gz 21` prints 3400 lines: a wide query returns each
        // record once and no fragment of a line cut by the reader.
        let every = ann.query("21", 1, 10_000_000).unwrap();
        assert_eq!(every.len(), 3400);
        assert!(every
            .iter()
            .all(|r| r.columns.len() == 5 && r.start >= 5010000));

        let batch = ann
            .query_batch(&[("21", 5031749, 5031749), ("21", 5033779, 5033793)])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 5031749, 5031749)),
            ["21\t5031749\tT\tG\t0.333"]
        );
        assert_eq!(batch.lookup("21", 5033779, 5033793).len(), 3);
    }

    /// A region query must return every record in range, not only the first.
    ///
    /// Records in a small file share one bgzf block, so a chunk's end is that
    /// block's boundary. A loop that compares the reader's `virtual_position()`
    /// against `chunk.end()` after each `read_line` exits holding one record,
    /// because `BufReader` fills in bulk and the position reaches the boundary as
    /// soon as the first line is pulled. `tabix f.gz 21:100-102` on an
    /// htslib-built copy of this fixture returns 3 rows.
    #[test]
    fn region_query_returns_every_record_in_range() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = ["#chrom\tpos\tval", "21\t100\t1", "21\t101\t2", "21\t102\t3"];
        let path = write_tabix_fixture(dir.path(), "multi.txt.gz", &lines, IndexSpec::standard());
        let ann = TabixAnnotator::open(TabixConfig {
            file_path: path,
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        })
        .unwrap();

        let recs = ann.query("21", 100, 102).unwrap();
        assert_eq!(
            recs.len(),
            3,
            "all three records share one bgzf block and must all be returned"
        );
        assert_eq!(
            recs.iter().map(|r| r.start).collect::<Vec<_>>(),
            vec![100, 101, 102]
        );

        // The region filter must still exclude out-of-range records that the
        // coarse, bin-granular index legitimately includes in the chunk.
        let narrowed = ann.query("21", 101, 101).unwrap();
        assert_eq!(narrowed.len(), 1);
        assert_eq!(narrowed[0].start, 101);

        let batch = ann.query_batch(&[("21", 100, 102)]).unwrap();
        assert_eq!(batch.lookup("21", 100, 102).len(), 3);
    }

    /// A VCF's records end at their INFO `END` when `end_col` names the INFO
    /// column and `ref_col` the REF column, as `tabix -p vcf` (htslib 1.11 or
    /// later) reads them (`tabix -p vcf` on this file prints the deletion for
    /// `21:250-250` and `21:300-300`, nothing for `21:301-301`); an `END` is
    /// read to its first non-digit (`END=165abc` ends at 165) and one before
    /// the start is ignored, as htslib 1.11 and later ignore it with one
    /// warning (`tabix` 1.22 indexes the file and prints the record for
    /// `21:153-153` by its reference allele's span), where 1.10 refuses to
    /// index the file (`hts_idx_push`: `end < begin`) and 1.9 indexes the
    /// inverted span silently, so `tabix` 1.9 prints the record for no
    /// position, `21:153-153` included; a record without a usable `END` ends
    /// at its start plus its REF allele's length less one (`tabix -p vcf`
    /// prints the four-base `del_ref` for `21:143-143` and not for
    /// `21:144-144`, and the `END=50` record for `21:153-153`), one base for a
    /// row with no REF cell. Without a `ref_col` the column is a generic end
    /// column, which no INFO text leads with an integer, so every record ends
    /// at its start.
    #[test]
    fn vcf_info_end_key_is_the_record_end() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "21\t100\tsv_del\tN\t<deletion>\t.\t.\tdbVar;SVTYPE=deletion;END=300;IMPRECISE",
            "21\t120\tsv_no_end\tN\t<insertion>\t.\t.\tSVTYPE=insertion",
            "21\t130\tsv_bad_end\tN\t<insertion>\t.\t.\tSVTYPE=insertion;END=x",
            "21\t140\tdel_ref\tACGT\tA\t.\t.\tDP=3",
            "21\t145\tend_in_ref\tACGT\tA\t.\t.\tEND=146",
            "21\t150\tend_before_start\tACGT\t<deletion>\t.\t.\tSVTYPE=deletion;END=50",
            "21\t160\tend_with_tail\tN\t<deletion>\t.\t.\tSVTYPE=deletion;END=165abc",
            "21\t170\tshort_row\tTAC",
            "21\t180\tno_ref",
        ];
        let path = write_tabix_fixture(dir.path(), "svs.vcf.gz", &lines_in, IndexSpec::standard());
        let spans_at = |ann: &TabixAnnotator, pos: u64| -> Vec<String> {
            ann.query("21", pos, pos)
                .unwrap()
                .iter()
                .map(|r| format!("{}:{}-{}", r.columns[2], r.start, r.end))
                .collect()
        };

        let without_ref = open(path.clone(), Some(7));
        assert_eq!(spans_at(&without_ref, 100), ["sv_del:100-100"]);
        assert!(spans_at(&without_ref, 250).is_empty());
        assert_eq!(spans_at(&without_ref, 120), ["sv_no_end:120-120"]);
        assert_eq!(spans_at(&without_ref, 130), ["sv_bad_end:130-130"]);
        assert!(spans_at(&without_ref, 141).is_empty());

        let with_ref = TabixAnnotator::open(TabixConfig {
            ref_col: Some(3),
            ..config(path, Some(7))
        })
        .unwrap();
        assert_eq!(
            spans_at(&with_ref, 143),
            ["sv_del:100-300", "del_ref:140-143"]
        );
        assert_eq!(spans_at(&with_ref, 144), ["sv_del:100-300"]);
        // An `END` inside the reference allele's span is the end, as htslib
        // through 1.21, Ensembl VEP's 1.9 among them, takes the key over the
        // allele (`tabix -p vcf` prints the record for `21:146-146` and nothing
        // for `21:147-147`); 1.22 ends the record at the larger of the two.
        assert_eq!(
            spans_at(&with_ref, 146),
            ["sv_del:100-300", "end_in_ref:145-146"]
        );
        assert_eq!(spans_at(&with_ref, 147), ["sv_del:100-300"]);
        assert_eq!(spans_at(&with_ref, 131), ["sv_del:100-300"]);
        assert_eq!(
            spans_at(&with_ref, 153),
            ["sv_del:100-300", "end_before_start:150-153"]
        );
        assert_eq!(spans_at(&with_ref, 154), ["sv_del:100-300"]);
        assert_eq!(
            spans_at(&with_ref, 165),
            ["sv_del:100-300", "end_with_tail:160-165"]
        );
        assert_eq!(spans_at(&without_ref, 150), ["end_before_start:150-150"]);
        // A row short of its end cell spans its reference allele, as the VCF
        // preset ends it, so the query two bases past its start still returns it.
        assert_eq!(
            spans_at(&with_ref, 172),
            ["sv_del:100-300", "short_row:170-172"]
        );
        assert_eq!(spans_at(&with_ref, 173), ["sv_del:100-300"]);
        assert!(spans_at(&without_ref, 171).is_empty());
        // A row with no reference cell keeps the one base its start column
        // gives it, the span htslib's `tbx_parse1` replaces only for a REF cell
        // of one or more bases, so a query from its start returns it and one
        // ending the base before it does not; without a reference column it is
        // one base as well.
        let spans_in = |ann: &TabixAnnotator, start: u64, end: u64| -> Vec<String> {
            ann.query("21", start, end)
                .unwrap()
                .iter()
                .map(|r| format!("{}:{}-{}", r.columns[2], r.start, r.end))
                .collect()
        };
        assert_eq!(
            spans_in(&with_ref, 180, 182),
            ["sv_del:100-300", "no_ref:180-180"]
        );
        assert_eq!(spans_in(&with_ref, 179, 179), ["sv_del:100-300"]);
        assert_eq!(spans_in(&without_ref, 180, 182), ["no_ref:180-180"]);

        assert_eq!(info_end("SVTYPE=deletion;END=300;IMPRECISE"), Some(300));
        assert_eq!(info_end("END=300"), Some(300));
        assert_eq!(info_end("END=300abc"), Some(300));
        assert_eq!(info_end("SVTYPE=deletion;END=x"), None);
        assert_eq!(info_end("SVTYPE=deletion"), None);
        assert_eq!(info_end("CIEND=1;SEND=5"), None);
    }

    #[test]
    fn test_merge_intervals_empty() {
        assert!(merge_intervals(&[]).is_empty());
    }

    #[test]
    fn test_merge_intervals_single() {
        let merged = merge_intervals(&[(100, 200)]);
        assert_eq!(merged, vec![(100, 200)]);
    }

    #[test]
    fn test_merge_intervals_overlapping() {
        // All three intervals are within 1000bp gap, so they merge into one.
        let merged = merge_intervals(&[(100, 200), (150, 300), (500, 600)]);
        assert_eq!(merged, vec![(100, 600)]);
    }

    #[test]
    fn test_merge_intervals_adjacent_within_gap() {
        let merged = merge_intervals(&[(100, 200), (800, 900)]);
        assert_eq!(merged, vec![(100, 900)]);
    }

    #[test]
    fn test_merge_intervals_far_apart() {
        let merged = merge_intervals(&[(100, 200), (5000, 6000)]);
        assert_eq!(merged, vec![(100, 200), (5000, 6000)]);
    }

    #[test]
    fn test_allele_normalize_snv() {
        let (r, a, offset) = allele_match::normalize_alleles("A", "G");
        assert_eq!(r, "A");
        assert_eq!(a, "G");
        assert_eq!(offset, 0);
    }

    #[test]
    fn test_allele_normalize_deletion() {
        let (r, a, offset) = allele_match::normalize_alleles("ACGT", "A");
        assert_eq!(r, "CGT");
        assert_eq!(a, "-");
        assert_eq!(offset, 1);
    }

    #[test]
    fn test_allele_normalize_insertion() {
        let (r, a, offset) = allele_match::normalize_alleles("A", "ACGT");
        assert_eq!(r, "-");
        assert_eq!(a, "CGT");
        assert_eq!(offset, 1);
    }

    #[test]
    fn test_allele_normalize_suffix_trim() {
        let (r, a, offset) = allele_match::normalize_alleles("TAC", "TC");
        assert_eq!(r, "A");
        assert_eq!(a, "-");
        assert_eq!(offset, 1);
    }

    #[test]
    fn test_chr_candidate_resolution() {
        let cands = chromosome_query_candidates("1");
        assert!(cands.contains(&"1".to_string()));
        assert!(cands.contains(&"chr1".to_string()));

        let cands_mt = chromosome_query_candidates("chrM");
        assert!(cands_mt.contains(&"chrM".to_string()));
        assert!(cands_mt.contains(&"M".to_string()));
        assert!(cands_mt.contains(&"MT".to_string()));
    }

    /// On GRCh37 `MT` and `chrMT` name the assembly's mitochondrial sequence
    /// and `chrM` and `M` another one, so the names tried with no check of
    /// the table's rows are those `BaseVepTabixPlugin.pm:693-708` tries,
    /// inside the name's pair; every other name, GRCh38 and a run without an
    /// assembly keep every spelling.
    #[test]
    fn a_mitochondrial_name_stays_in_its_pair_on_grch37() {
        let grch37 = Some(Assembly::Grch37);
        for (name, tried) in [
            ("MT", ["MT", "chrMT"]),
            ("chrMT", ["chrMT", "MT"]),
            ("chrM", ["chrM", "M"]),
            ("M", ["M", "chrM"]),
            ("mt", ["mt", "chrmt"]),
        ] {
            assert_eq!(plugin_table_candidates(name, grch37), tried, "{name}");
        }
        assert_eq!(
            plugin_table_candidates("ChrM", grch37),
            ["ChrM", "M", "chrM"]
        );
        assert_eq!(
            plugin_table_candidates("21", grch37),
            chromosome_query_candidates("21")
        );
        for assembly in [Some(Assembly::Grch38), None] {
            for name in ["MT", "chrM"] {
                assert_eq!(
                    plugin_table_candidates(name, assembly),
                    chromosome_query_candidates(name),
                    "{name} under {assembly:?}"
                );
            }
        }
        for name in ["M", "MT", "chrM", "chrMT", "ChrMt"] {
            assert!(is_mitochondrial_name(name), "{name}");
        }
        for name in ["1", "chr1", "MT2", "chrMito", "Mito", ""] {
            assert!(!is_mitochondrial_name(name), "{name}");
        }
    }

    /// A mitochondrial record on GRCh37 is named to its table by its own
    /// spelling (`$vf->{chr}`), which says which pair of names it is of; every
    /// other record, and every record on GRCh38 or without an assembly, by the
    /// cache's name of its contig.
    #[test]
    fn table_contig_is_the_records_own_spelling_for_a_mitochondrial_record_on_grch37() {
        let record = |spelling: &str| {
            InputVariant::new(spelling.into(), 3460, 3460, b"G".to_vec(), b"A".to_vec())
        };
        for spelling in ["MT", "chrMT", "chrM", "M"] {
            let record = record(spelling);
            assert_eq!(table_contig(&record, Some(Assembly::Grch37)), spelling);
            for assembly in [Some(Assembly::Grch38), None] {
                assert_eq!(
                    table_contig(&record, assembly),
                    "MT",
                    "{spelling} under {assembly:?}"
                );
            }
        }
        assert_eq!(table_contig(&record("chr21"), Some(Assembly::Grch37)), "21");
    }

    use crate::test_fixtures::{made_up_alternate, made_up_reference, made_up_rows, repeated_base};

    /// A made-up five-column table naming its one contig `contig`, with one
    /// row at each of `positions` ([`made_up_rows`]).
    fn made_up_table(
        dir: &Path,
        name: &str,
        contig: &str,
        bases: &[u8],
        positions: impl Iterator<Item = usize>,
        offset: usize,
    ) -> PathBuf {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        let mut lines = vec!["#chr\tpos\tref\talt\tscore".to_string()];
        lines.extend(made_up_rows(
            contig,
            bases,
            positions,
            offset,
            |contig, position, reference, alt| {
                format!("{contig}\t{position}\t{reference}\t{alt}\t0.{position}")
            },
        ));
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
    }

    fn allele_config(file_path: PathBuf) -> TabixConfig {
        TabixConfig {
            file_path,
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: Some(2),
            alt_col: Some(3),
            zero_based: false,
        }
    }

    /// The table at `path` opened on GRCh37 with the made-up `reference`.
    fn on_grch37(path: &Path, reference: Option<&Arc<vep_fasta::IndexedFasta>>) -> TabixAnnotator {
        let mut table = TabixAnnotator::open(allele_config(path.to_path_buf()))
            .unwrap()
            .on_assembly(Some(Assembly::Grch37));
        table.frame_mut().set_reference(reference.cloned());
        table
    }

    /// On GRCh37 a contig the table names under the other pair of names than
    /// the record's is read once the table's rows are found to state the
    /// reference bases at their own positions: a record on `MT` or `chrMT`
    /// reads a table naming `chrM` or `M`, and a record on `chrM` or `M` one
    /// naming `MT` or `chrMT`.
    #[test]
    fn on_grch37_a_table_under_the_other_pair_is_read_when_its_rows_state_the_reference() {
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference(dir.path(), 1600);
        for (contig, records) in [
            ("chrM", ["MT", "chrMT"]),
            ("M", ["MT", "chrMT"]),
            ("MT", ["chrM", "M"]),
            ("chrMT", ["chrM", "M"]),
        ] {
            let path = made_up_table(
                dir.path(),
                &format!("stating_{contig}_synthetic.tsv.gz"),
                contig,
                &bases,
                (13..1590).step_by(13),
                0,
            );
            let table = on_grch37(&path, Some(&reference));
            for record in records {
                let rows = table.query(record, 130, 130).unwrap();
                assert_eq!(rows.len(), 1, "record {record}, table {contig}");
                assert_eq!(rows[0].chr, contig);
                assert_eq!(rows[0].columns[4], "0.130");
                assert!(table.holds_chromosome(record), "{record} against {contig}");
                assert_eq!(table.query_windows(record, &[(130, 130)]).unwrap().len(), 1);
                assert_eq!(
                    table
                        .query_batch(&[(record, 130, 130)])
                        .unwrap()
                        .lookup(record, 130, 130)
                        .len(),
                    1
                );
            }
        }
    }

    /// A table whose rows state the bases one position on counts on another
    /// sequence: on GRCh37 no row of it is read for a record of the other
    /// pair of names, the row included whose alleles are the record's because
    /// the reference repeats its base there. A record of the table's own pair
    /// reads it, as `BaseVepTabixPlugin.pm:693-708` does.
    #[test]
    fn on_grch37_a_table_one_position_on_is_not_read_across_the_pairs() {
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference(dir.path(), 1600);
        let repeated = repeated_base(&bases, 100);
        let path = made_up_table(
            dir.path(),
            "one_position_on_synthetic.tsv.gz",
            "chrM",
            &bases,
            10..1590,
            1,
        );
        let table = on_grch37(&path, Some(&reference));
        for record in ["MT", "chrMT"] {
            for position in [130, repeated as u64] {
                assert!(
                    table.query(record, position, position).unwrap().is_empty(),
                    "record {record} at {position}"
                );
            }
            assert!(!table.holds_chromosome(record));
        }
        for record in ["chrM", "M"] {
            let rows = table
                .query(record, repeated as u64, repeated as u64)
                .unwrap();
            assert_eq!(rows.len(), 1, "record {record}");
            assert_eq!(rows[0].columns[2].as_bytes(), [bases[repeated - 1]]);
        }
    }

    /// The check takes more than rows that agree: a table of rows at fewer
    /// than eight positions is not read across the pairs though each of its
    /// rows states the reference, a row past the end of the reference
    /// sequence refuses a table whose other rows state it, and one row that
    /// states another base refuses a table.
    #[test]
    fn on_grch37_a_few_agreeing_rows_do_not_pass_and_one_row_of_another_sequence_refuses() {
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference(dir.path(), 1600);
        let read = |name: &str, lines: Vec<String>| {
            use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
            let path = write_tabix_fixture(dir.path(), name, &lines, IndexSpec::standard());
            let table = on_grch37(&path, Some(&reference));
            table.query("MT", 130, 130).unwrap().len()
        };
        let row = |position: usize, offset: usize| {
            let reference = bases[position - 1 + offset];
            format!(
                "chrM\t{position}\t{}\t{}\t0.{position}",
                reference as char,
                made_up_alternate(reference) as char
            )
        };
        let header = "#chr\tpos\tref\talt\tscore".to_string();

        let seven: Vec<String> = std::iter::once(header.clone())
            .chain((0..7).map(|i| row(130 + 100 * i, 0)))
            .collect();
        assert_eq!(read("seven_positions_synthetic.tsv.gz", seven), 0);

        let eight: Vec<String> = std::iter::once(header.clone())
            .chain((0..8).map(|i| row(130 + 100 * i, 0)))
            .collect();
        assert_eq!(read("eight_positions_synthetic.tsv.gz", eight), 1);

        let mut past_the_end: Vec<String> = std::iter::once(header.clone())
            .chain((13..1590).step_by(13).map(|position| row(position, 0)))
            .collect();
        past_the_end.push("chrM\t1601\tA\tC\t0.1601".to_string());
        assert_eq!(read("past_the_end_synthetic.tsv.gz", past_the_end), 0);

        let one_other_base: Vec<String> = std::iter::once(header)
            .chain((13..1590).step_by(13).map(|position| {
                if position == 650 {
                    let stated = made_up_alternate(bases[position - 1]);
                    format!("chrM\t{position}\t{}\tN\t0.{position}", stated as char)
                } else {
                    row(position, 0)
                }
            }))
            .collect();
        assert_eq!(read("one_other_base_synthetic.tsv.gz", one_other_base), 0);
    }

    /// A row is one position however often the check reads it: a table of
    /// rows at seven positions, one of them a span that reaches past the end
    /// of the reference sequence and is read again there, is not read across
    /// the pairs though each of its rows states the reference.
    #[test]
    fn on_grch37_a_span_read_again_past_the_end_is_one_position() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference(dir.path(), 1600);
        let mut lines = vec!["#chr\tpos\tend\tref\talt\tscore".to_string()];
        lines.extend((0..7).map(|i| {
            let position = 130 + 100 * i;
            let end = if i == 5 { 1700 } else { position };
            let reference = bases[position - 1];
            format!(
                "chrM\t{position}\t{end}\t{}\t{}\t0.{position}",
                reference as char,
                made_up_alternate(reference) as char
            )
        }));
        let path = write_tabix_fixture(
            dir.path(),
            "seven_positions_one_span_synthetic.tsv.gz",
            &lines,
            IndexSpec::standard().with_span_cols(1, 2, 3),
        );
        let mut table = TabixAnnotator::open(config(path, Some(2)))
            .unwrap()
            .on_assembly(Some(Assembly::Grch37));
        table.frame_mut().set_reference(Some(reference));
        table.frame_mut().set_allele_column(Some(3));
        assert!(table.query("MT", 130, 130).unwrap().is_empty());
        assert_eq!(table.query("chrM", 130, 130).unwrap().len(), 1);
    }

    /// Without a reference sequence the comparison cannot be made: on GRCh37
    /// the table is read under the record's own pair alone, as
    /// `BaseVepTabixPlugin.pm:693-708` reads it, and one warning names the
    /// table, however many records ask.
    #[test]
    fn without_a_reference_sequence_the_pairs_stay_apart_and_one_warning_names_the_table() {
        use crate::test_fixtures::Warnings;
        let dir = tempfile::tempdir().unwrap();
        let (_, bases) = made_up_reference(dir.path(), 1600);
        let path = made_up_table(
            dir.path(),
            "no_reference_synthetic.tsv.gz",
            "chrM",
            &bases,
            (13..1590).step_by(13),
            0,
        );
        let table = on_grch37(&path, None);
        let warnings = Arc::new(Warnings::default());
        tracing::subscriber::with_default(Arc::clone(&warnings), || {
            for record in ["MT", "chrMT", "MT"] {
                assert!(
                    table.query(record, 130, 130).unwrap().is_empty(),
                    "{record}"
                );
            }
            assert_eq!(table.query("chrM", 130, 130).unwrap().len(), 1);
            assert_eq!(table.query("M", 130, 130).unwrap().len(), 1);
        });
        let warnings = warnings.0.lock().unwrap();
        assert_eq!(warnings.len(), 1, "{warnings:?}");
        assert_eq!(
            warnings[0],
            format!(
                "{}: its contig chrM is not read for a mitochondrial record written MT or \
                 chrMT: on GRCh37 the names MT and chrMT are one sequence and chrM and M \
                 another, and the table's rows could not be compared with the reference \
                 sequence (no reference sequence was given; pass --fasta)",
                path.display()
            )
        );
    }

    /// A reference that holds its mitochondrial sequence as `chrM` at another
    /// length than the assembly's holds the other naming's sequence: a table
    /// whose rows state its bases is not read for a record on `MT`, and the
    /// warning says why.
    #[test]
    fn a_reference_of_the_other_naming_and_length_is_not_compared_with() {
        use crate::test_fixtures::{made_up_reference_named, Warnings};
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference_named(dir.path(), "chrM", 1600);
        assert_eq!(reference.record_name("MT"), Some("chrM"));
        let path = made_up_table(
            dir.path(),
            "other_naming_synthetic.tsv.gz",
            "chrM",
            &bases,
            (13..1590).step_by(13),
            0,
        );
        let table = on_grch37(&path, Some(&reference));
        let warnings = Arc::new(Warnings::default());
        tracing::subscriber::with_default(Arc::clone(&warnings), || {
            assert!(table.query("MT", 130, 130).unwrap().is_empty());
            assert_eq!(table.query("chrM", 130, 130).unwrap().len(), 1);
        });
        let warnings = warnings.0.lock().unwrap();
        assert_eq!(warnings.len(), 1, "{warnings:?}");
        assert!(
            warnings[0].ends_with(
                "(the reference sequence holds its mitochondrial sequence under the other pair \
                 of names, at another length than the assembly's)"
            ),
            "{}",
            warnings[0]
        );
    }

    /// A reader told to keep to the record's own pair reads no contig of the
    /// other pair on GRCh37 though its rows state the reference, and warns of
    /// nothing; on GRCh38 it reads every spelling.
    #[test]
    fn a_reader_kept_to_the_own_pair_makes_no_comparison() {
        use crate::test_fixtures::Warnings;
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference(dir.path(), 1600);
        let path = made_up_table(
            dir.path(),
            "own_pair_alone_synthetic.tsv.gz",
            "M",
            &bases,
            (13..1590).step_by(13),
            0,
        );
        let warnings = Arc::new(Warnings::default());
        tracing::subscriber::with_default(Arc::clone(&warnings), || {
            let mut table = on_grch37(&path, Some(&reference));
            table.frame_mut().keep_to_own_pair();
            assert!(table.query("MT", 130, 130).unwrap().is_empty());
            assert!(table.query("chrMT", 130, 130).unwrap().is_empty());
            assert_eq!(table.query("chrM", 130, 130).unwrap().len(), 1);
            let mut table = TabixAnnotator::open(allele_config(path.clone()))
                .unwrap()
                .on_assembly(Some(Assembly::Grch38));
            table.frame_mut().keep_to_own_pair();
            assert_eq!(table.query("MT", 130, 130).unwrap().len(), 1);
        });
        assert!(warnings.0.lock().unwrap().is_empty());
    }

    /// A table with no reference-allele column states no base to compare: on
    /// GRCh37 it is read under the record's own pair alone, with or without a
    /// reference sequence, and one warning names the table and the reason.
    #[test]
    fn a_table_with_no_reference_column_stays_in_the_records_pair() {
        use crate::test_fixtures::Warnings;
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference(dir.path(), 1600);
        let path = made_up_table(
            dir.path(),
            "no_allele_column_synthetic.tsv.gz",
            "chrM",
            &bases,
            (13..1590).step_by(13),
            0,
        );
        let mut table = TabixAnnotator::open(config(path.clone(), None))
            .unwrap()
            .on_assembly(Some(Assembly::Grch37));
        table.frame_mut().set_reference(Some(reference));
        let warnings = Arc::new(Warnings::default());
        tracing::subscriber::with_default(Arc::clone(&warnings), || {
            assert!(table.query("MT", 130, 130).unwrap().is_empty());
            assert!(table.query("MT", 143, 143).unwrap().is_empty());
            assert_eq!(table.query("chrM", 130, 130).unwrap().len(), 1);
        });
        let warnings = warnings.0.lock().unwrap();
        assert_eq!(warnings.len(), 1, "{warnings:?}");
        assert!(
            warnings[0].starts_with(&format!("{}: its contig chrM is not read", path.display()))
                && warnings[0].ends_with("(the table has no reference-allele column)"),
            "{}",
            warnings[0]
        );
    }

    /// On GRCh38, and in a run without an assembly, the four names are one
    /// sequence: a record, named there by the cache's `MT` whatever its own
    /// spelling ([`table_contig`]), reads every spelling of the table with no
    /// comparison made, a table one position on included, and no reference
    /// sequence is asked for.
    #[test]
    fn off_grch37_every_spelling_is_read_as_before() {
        use crate::test_fixtures::Warnings;
        let dir = tempfile::tempdir().unwrap();
        let (_, bases) = made_up_reference(dir.path(), 1600);
        let warnings = Arc::new(Warnings::default());
        tracing::subscriber::with_default(Arc::clone(&warnings), || {
            for contig in ["MT", "chrMT", "M", "chrM"] {
                let path = made_up_table(
                    dir.path(),
                    &format!("one_position_on_{contig}_synthetic.tsv.gz"),
                    contig,
                    &bases,
                    10..1590,
                    1,
                );
                for assembly in [Some(Assembly::Grch38), None] {
                    let table = TabixAnnotator::open(allele_config(path.clone()))
                        .unwrap()
                        .on_assembly(assembly);
                    for spelling in ["MT", "chrMT", "M", "chrM"] {
                        let record = InputVariant::new(
                            spelling.into(),
                            130,
                            130,
                            vec![bases[130]],
                            b"N".to_vec(),
                        );
                        assert_eq!(
                            table
                                .query(table_contig(&record, assembly), 130, 130)
                                .unwrap()
                                .len(),
                            1,
                            "record {spelling}, table {contig}, {assembly:?}"
                        );
                    }
                }
            }
        });
        assert!(warnings.0.lock().unwrap().is_empty());
    }

    /// The store built from a table gets the table's verdict, since it holds
    /// the table's rows: on GRCh37 it is read across the pairs when its rows
    /// state the reference and not when they state the bases one position on.
    #[test]
    fn on_grch37_a_store_is_read_across_the_pairs_as_its_table_is() {
        use crate::annotation_store::{open_backend, AnnotationStore, Backend};
        use crate::binary_store::{BinaryStoreWriter, ContentIdentity};
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference(dir.path(), 1600);
        for (name, offset, read) in [("stating", 0, 1), ("one_on", 1, 0)] {
            let path = made_up_table(
                dir.path(),
                &format!("{name}_store_synthetic.tsv.gz"),
                "chrM",
                &bases,
                10..1590,
                offset,
            );
            let table = on_grch37(&path, Some(&reference));
            let rows = table.query("chrM", 1, 1600).unwrap();
            let header: Vec<String> = table.header().unwrap().to_vec();
            let rows: Vec<TabixRecord> = rows
                .into_iter()
                .map(|mut row| {
                    row.columns[0] = "chrM".to_string();
                    row
                })
                .collect();
            BinaryStoreWriter::new(dir.path().join(format!("{name}_store_synthetic.tsv.vpd")))
                .write(
                    &header,
                    &[("chrM".to_string(), rows)],
                    &allele_config(path.clone()),
                    &ContentIdentity::of_file(&path).unwrap(),
                )
                .unwrap();
            let backend = open_backend(
                &path,
                allele_config(path.clone()),
                Some(Assembly::Grch37),
                |table| table,
            )
            .unwrap();
            let Backend::Store(mut store) = backend else {
                panic!("{name}: the store beside the table was not opened");
            };
            store.set_reference_fasta(Some(Arc::clone(&reference)));
            assert_eq!(store.query("MT", 130, 130).unwrap().len(), read, "{name}");
            assert_eq!(
                store
                    .query_batch(&[("chrMT", 130, 130)])
                    .unwrap()
                    .lookup("chrMT", 130, 130)
                    .len(),
                read,
                "{name}"
            );
            assert_eq!(store.query("chrM", 130, 130).unwrap().len(), 1, "{name}");
            assert_eq!(table.query("MT", 130, 130).unwrap().len(), read, "{name}");
        }
    }

    /// The names a mitochondrial record's rows are asked under on GRCh37: its
    /// own pair first, then the other pair, marked as read only after the
    /// check; a name a plugin's module gives every such record leads.
    #[test]
    fn on_grch37_the_other_pair_follows_the_records_own_and_a_modules_name_leads() {
        let names = |chr: &str, first: Option<&str>| -> Vec<(String, bool)> {
            grch37_mitochondrial_candidates(chr, first)
        };
        let listed = |pairs: &[(&str, bool)]| -> Vec<(String, bool)> {
            pairs
                .iter()
                .map(|(name, crosses)| (name.to_string(), *crosses))
                .collect()
        };
        assert_eq!(
            names("MT", None),
            listed(&[("MT", false), ("chrMT", false), ("M", true), ("chrM", true)])
        );
        assert_eq!(
            names("chrMT", None),
            listed(&[("chrMT", false), ("MT", false), ("chrM", true), ("M", true)])
        );
        assert_eq!(
            names("chrM", None),
            listed(&[("chrM", false), ("M", false), ("chrMT", true), ("MT", true)])
        );
        assert_eq!(
            names("M", None),
            listed(&[("M", false), ("chrM", false), ("MT", true), ("chrMT", true)])
        );
        assert_eq!(
            names("MT", Some("M")),
            listed(&[("M", true), ("chrM", true), ("MT", false), ("chrMT", false)])
        );
        assert_eq!(
            names("chrM", Some("M")),
            listed(&[("M", false), ("chrM", false), ("chrMT", true), ("MT", true)])
        );
    }

    #[test]
    fn test_matches_allele_with_position_requires_normalized_start() {
        let record = TabixRecord {
            chr: "1".to_string(),
            start: 100,
            end: 102,
            columns: vec!["1".into(), "100".into(), "TAC".into(), "TC".into()],
            header: None,
        };
        let variant = InputVariant::new("1".into(), 101, 101, b"A".to_vec(), b"-".to_vec());
        assert!(allele_match::matches_allele_with_position(
            &record, &variant, 1, 2, 3
        ));

        let variant_wrong_pos =
            InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"-".to_vec());
        assert!(!allele_match::matches_allele_with_position(
            &record,
            &variant_wrong_pos,
            1,
            2,
            3
        ));
    }

    /// A record given on the reverse strand is matched through its forward-strand
    /// alleles, as `get_matched_variant_alleles` reverse-complements them
    /// (Utils/Sequence.pm:1182, 1196): `T/C` on strand -1 matches the file's `A G`
    /// row and not its `T C` row.
    #[test]
    fn reverse_strand_record_matches_the_complemented_file_row() {
        use vep_core::coordinate::Strand;
        let row = |reference: &str, alt: &str| TabixRecord {
            chr: "21".to_string(),
            start: 100,
            end: 100,
            columns: vec!["21".into(), "100".into(), reference.into(), alt.into()],
            header: None,
        };
        let mut reverse = InputVariant::new("21".into(), 100, 100, b"T".to_vec(), b"C".to_vec());
        reverse.strand = Strand::Reverse;
        assert!(allele_match::matches_allele_with_position(
            &row("A", "G"),
            &reverse,
            1,
            2,
            3
        ));
        assert!(!allele_match::matches_allele_with_position(
            &row("T", "C"),
            &reverse,
            1,
            2,
            3
        ));
    }

    /// A row whose allele cells are written in lower case matches the record it
    /// describes, whose alleles the parser upper-cases: the shared bases are
    /// found and the trimmed alleles compared in either letter case, where
    /// Ensembl VEP's key is compared as written (Utils/Sequence.pm:1207, 1232;
    /// the documented class `tabix_allele_case_insensitive_match`). A row in
    /// mixed case loses the bases it shares as letters.
    #[test]
    fn matched_variant_alleles_reads_a_row_in_either_letter_case() {
        let row = |start: u64, reference: &str, alt: &str| TabixRecord {
            chr: "21".to_string(),
            start,
            end: start + reference.len() as u64 - 1,
            columns: vec!["21".into(), start.to_string(), reference.into(), alt.into()],
            header: None,
        };
        let substitution = InputVariant::new("21".into(), 100, 100, b"G".to_vec(), b"A".to_vec());
        for (reference, alt) in [("g", "a"), ("G", "a"), ("g", "A"), ("G", "A")] {
            assert!(
                allele_match::matched_variant_alleles(
                    &row(100, reference, alt),
                    &substitution,
                    1,
                    2,
                    3
                ),
                "{reference}/{alt}"
            );
            assert_eq!(
                allele_match::matched_or_undefined(
                    &row(100, reference, alt),
                    &substitution,
                    1,
                    2,
                    3
                ),
                Some(true),
                "{reference}/{alt}"
            );
        }
        assert!(!allele_match::matched_variant_alleles(
            &row(100, "g", "t"),
            &substitution,
            1,
            2,
            3
        ));

        let insertion = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"AT".to_vec());
        for (reference, alt) in [("g", "gat"), ("G", "gat"), ("g", "GAT")] {
            assert!(
                allele_match::matched_variant_alleles(
                    &row(100, reference, alt),
                    &insertion,
                    1,
                    2,
                    3
                ),
                "{reference}/{alt}"
            );
        }
    }

    /// `get_matched_variant_alleles` trims the variant's alleles as well as the
    /// file's (Utils/Sequence.pm:1199-1208, 1223-1241): a substitution written
    /// with an unchanged leading base (`GG/GA` at 100) is the file's `G A` row at
    /// 101, which the one-sided match misses; an insertion written VCF-style in
    /// the file (`G GAT` at 100) is the variant `-/AT` at 101; a deletion inside
    /// a repeat (`ATA/A` at 100) matches a row trimmed in either order; a row at
    /// another position, or of other alleles, does not match; and a variant on
    /// the reverse strand compares through its forward-strand alleles.
    #[test]
    fn matched_variant_alleles_trims_both_sides_in_either_order() {
        use vep_core::coordinate::Strand;
        let row = |start: u64, reference: &str, alt: &str| TabixRecord {
            chr: "21".to_string(),
            start,
            end: start + reference.len() as u64 - 1,
            columns: vec!["21".into(), start.to_string(), reference.into(), alt.into()],
            header: None,
        };
        let variant = |start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            InputVariant::new("21".into(), start, end, reference.to_vec(), alt.to_vec())
        };

        let shared_base = variant(100, 101, b"GG", b"GA");
        assert!(allele_match::matched_variant_alleles(
            &row(101, "G", "A"),
            &shared_base,
            1,
            2,
            3
        ));
        assert!(!allele_match::matches_allele_with_position(
            &row(101, "G", "A"),
            &shared_base,
            1,
            2,
            3
        ));
        assert!(!allele_match::matched_variant_alleles(
            &row(100, "G", "A"),
            &shared_base,
            1,
            2,
            3
        ));
        assert!(!allele_match::matched_variant_alleles(
            &row(101, "G", "T"),
            &shared_base,
            1,
            2,
            3
        ));

        let snv = variant(101, 101, b"G", b"A");
        assert!(allele_match::matched_variant_alleles(
            &row(101, "G", "A"),
            &snv,
            1,
            2,
            3
        ));
        assert!(allele_match::matched_variant_alleles(
            &row(100, "GG", "GA"),
            &snv,
            1,
            2,
            3
        ));
        assert!(!allele_match::matched_variant_alleles(
            &row(102, "G", "A"),
            &snv,
            1,
            2,
            3
        ));

        let insertion = variant(101, 100, b"-", b"AT");
        assert!(allele_match::matched_variant_alleles(
            &row(100, "G", "GAT"),
            &insertion,
            1,
            2,
            3
        ));
        assert!(!allele_match::matched_variant_alleles(
            &row(100, "G", "GTA"),
            &insertion,
            1,
            2,
            3
        ));

        let deletion = variant(100, 102, b"ATA", b"A");
        assert!(allele_match::matched_variant_alleles(
            &row(101, "TA", "-"),
            &deletion,
            1,
            2,
            3
        ));
        assert!(allele_match::matched_variant_alleles(
            &row(100, "AT", "-"),
            &deletion,
            1,
            2,
            3
        ));
        assert!(allele_match::matched_variant_alleles(
            &row(99, "GATA", "GA"),
            &deletion,
            1,
            2,
            3
        ));
        assert!(!allele_match::matched_variant_alleles(
            &row(100, "TA", "-"),
            &deletion,
            1,
            2,
            3
        ));

        let mut reverse = variant(100, 101, b"CC", b"TC");
        reverse.strand = Strand::Reverse;
        assert!(allele_match::matched_variant_alleles(
            &row(101, "G", "A"),
            &reverse,
            1,
            2,
            3
        ));
        assert!(!allele_match::matched_variant_alleles(
            &row(100, "G", "A"),
            &reverse,
            1,
            2,
            3
        ));
    }

    #[test]
    fn test_lookup_handles_reversed_range_bounds() {
        let mut chr_positions: BTreeMap<u64, Vec<TabixRecord>> = BTreeMap::new();
        chr_positions.insert(
            150,
            vec![TabixRecord {
                chr: "1".to_string(),
                start: 150,
                end: 150,
                columns: vec![],
                header: None,
            }],
        );

        let mut records_by_pos: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>> = HashMap::new();
        records_by_pos.insert("1".to_string(), chr_positions);

        let batch = BatchQueryResult::from_records(records_by_pos);
        let hits = batch.lookup("1", 200, 100);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].start, 150);
    }

    /// A long record keyed by its start must be found from every position it
    /// covers, including after many shorter records that follow it. Expected
    /// hits are the records whose `start..=end` overlaps the range.
    #[test]
    fn lookup_finds_records_spanning_the_range_start() {
        let record = |start: u64, end: u64| TabixRecord {
            chr: "1".to_string(),
            start,
            end,
            columns: vec![],
            header: None,
        };
        let mut by_start: BTreeMap<u64, Vec<TabixRecord>> = BTreeMap::new();
        by_start.insert(100, vec![record(100, 5000)]);
        for start in (200..=4900).step_by(100) {
            by_start.insert(start, vec![record(start, start)]);
        }
        let mut records_by_pos = HashMap::new();
        records_by_pos.insert("1".to_string(), by_start);
        let batch = BatchQueryResult::from_records(records_by_pos);

        let starts = |hits: Vec<&TabixRecord>| hits.iter().map(|r| r.start).collect::<Vec<_>>();
        assert_eq!(starts(batch.lookup("1", 3000, 3000)), vec![100, 3000]);
        assert_eq!(starts(batch.lookup("1", 4950, 4999)), vec![100]);
        assert_eq!(starts(batch.lookup("1", 5001, 6000)), Vec::<u64>::new());
        assert_eq!(starts(batch.lookup("1", 50, 99)), Vec::<u64>::new());
        assert!(batch.lookup("2", 3000, 3000).is_empty());
    }

    /// The index's own recorded start column is the only thing distinguishing
    /// the two same-header copies of a dual-coordinate file (REVEL). If this
    /// reads the wrong field, the wrong-assembly guards built on it are inert.
    #[test]
    fn indexed_start_column_reports_the_column_the_tbi_was_built_on() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "#chr\thg19_pos\tgrch38_pos\tref\talt",
            "21\t9907250\t9068417\tC\tA",
        ];

        let on_col2 = write_tabix_fixture(
            dir.path(),
            "col2.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(2),
        );
        let on_col3 = write_tabix_fixture(
            dir.path(),
            "col3.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(3),
        );

        // 0-based: tabix's 1-based -b 2 / -b 3 become 1 / 2.
        assert_eq!(indexed_start_column(&on_col2).unwrap(), Some(1));
        assert_eq!(indexed_start_column(&on_col3).unwrap(), Some(2));

        // No index present is not an error: callers skip the check.
        let orphan = dir.path().join("no_index.tsv.gz");
        std::fs::write(&orphan, b"").unwrap();
        assert_eq!(indexed_start_column(&orphan).unwrap(), None);
    }

    /// A data file with no extension must still resolve its index: REVEL's
    /// upstream zip member unzips to a bare `revel_with_transcript_ids`, and
    /// replacing the extension yields `revel_with_transcript_ids..tbi`, a path
    /// that never exists.
    #[test]
    fn indexed_start_column_resolves_an_extensionless_data_file() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "#chr\thg19_pos\tgrch38_pos\tref\talt",
            "21\t9907250\t9068417\tC\tA",
        ];

        // Written with an extension, then renamed to the bare upstream shape so
        // both the data file and its `.tbi` carry no extension of their own.
        let staged = write_tabix_fixture(
            dir.path(),
            "staged.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(3),
        );
        let bare = dir.path().join("revel_with_transcript_ids");
        std::fs::rename(&staged, &bare).unwrap();
        std::fs::rename(
            staged.with_file_name("staged.tsv.gz.tbi"),
            dir.path().join("revel_with_transcript_ids.tbi"),
        )
        .unwrap();

        assert_eq!(indexed_start_column(&bare).unwrap(), Some(2));
    }

    /// A dual-coordinate file resolves the position column by its header name
    /// and is refused when its index was built on the other coordinate column;
    /// a header without the column is refused by name, and a file without an
    /// index is accepted on the header alone.
    #[test]
    fn resolve_indexed_position_column_checks_the_header_and_the_index() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "#chr\thg19_pos\tgrch38_pos\tref\talt",
            "21\t9907250\t9068417\tC\tA",
        ];
        let header: Vec<String> = ["chr", "hg19_pos", "grch38_pos", "ref", "alt"]
            .iter()
            .map(|name| name.to_string())
            .collect();
        let on_col2 = write_tabix_fixture(
            dir.path(),
            "col2.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(2),
        );
        let on_col3 = write_tabix_fixture(
            dir.path(),
            "col3.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(3),
        );

        assert_eq!(
            resolve_indexed_position_column(&on_col2, &header, "hg19_pos").unwrap(),
            1
        );
        assert_eq!(
            resolve_indexed_position_column(&on_col3, &header, "grch38_pos").unwrap(),
            2
        );

        let mismatch = resolve_indexed_position_column(&on_col2, &header, "grch38_pos")
            .expect_err("an index on hg19_pos must be refused for grch38_pos")
            .to_string();
        assert!(
            mismatch.contains("indexed on 'hg19_pos' (column 2)")
                && mismatch.contains("not on 'grch38_pos' (column 3)"),
            "{mismatch}"
        );

        let missing = resolve_indexed_position_column(&on_col2, &header, "pos")
            .expect_err("a header without the column must be refused")
            .to_string();
        assert!(
            missing.contains("has no 'pos' column") && missing.contains("hg19_pos,grch38_pos"),
            "{missing}"
        );

        let orphan = dir.path().join("no_index.tsv.gz");
        std::fs::write(&orphan, b"").unwrap();
        assert_eq!(
            resolve_indexed_position_column(&orphan, &header, "grch38_pos").unwrap(),
            2
        );

        // The check alone, for a plugin that resolves the column by its own
        // rule: the answer names the index's sequence-name column where the
        // index agrees and its start column where it does not.
        assert_eq!(
            indexed_position(&on_col2, 1).unwrap(),
            IndexedPosition::Agrees { sequence_column: 0 }
        );
        assert_eq!(
            indexed_position(&on_col2, 2).unwrap(),
            IndexedPosition::Disagrees { indexed_column: 1 }
        );
        assert_eq!(
            indexed_position(&orphan, 2).unwrap(),
            IndexedPosition::Unchecked
        );
    }

    /// A row whose start cell is a digit string at 2^63 - 1, 2^63 or 2^64 - 1
    /// starts past the last position a `.tbi` holds, so no read takes it for
    /// a row of the file's order: the per-window query, the sweep and the
    /// block read, on one thread, on four and in slices of one block, each
    /// pass over it and return the sound rows on both sides of it, whether
    /// the row after it shares its block, opens the next block or lies in
    /// another tile. `tabix` (htslib 1.9) indexes such a row only where it
    /// opens its chromosome, and prints the rows after it for every region
    /// that holds them: the last layout.
    #[test]
    fn sound_rows_read_alike_beside_a_start_at_the_integer_limit() {
        type Layout = (Vec<String>, usize, Vec<(u64, u64)>, &'static [&'static str]);
        let position = |block: u64| bgzf::VirtualPosition::try_from((block, 0)).unwrap();
        let names = |records: Vec<TabixRecord>| -> Vec<String> {
            records
                .into_iter()
                .map(|record| record.columns[2].clone())
                .collect()
        };
        let pool = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
        };
        let before = "21\t100\tb1\n21\t200\tb2\n21\t300\tb3\n";
        for limit in [i64::MAX as u64, 1 << 63, u64::MAX] {
            let limit_row = format!("21\t{limit}\tlimit\n");
            // The blocks of each file, the block the tile of 40000 opens with,
            // the windows asked for and the rows they hold.
            let layouts: [Layout; 3] = [
                (
                    vec![
                        format!("{before}{limit_row}21\t400\ta1\n"),
                        "21\t500\ta2\n21\t600\ta3\n".to_string(),
                    ],
                    0,
                    vec![(100, 100), (200, 300), (400, 400), (500, 600)],
                    &["b1", "b2", "b3", "a1", "a2", "a3"],
                ),
                (
                    vec![
                        format!("{before}{limit_row}"),
                        "21\t40000\ta1\n21\t40100\ta2\n21\t40200\ta3\n".to_string(),
                    ],
                    1,
                    vec![(100, 300), (40_000, 40_200)],
                    &["b1", "b2", "b3", "a1", "a2", "a3"],
                ),
                (
                    vec![format!(
                        "{limit_row}21\t400\ta1\n21\t500\ta2\n21\t600\ta3\n"
                    )],
                    0,
                    vec![(400, 400), (500, 600)],
                    &["a1", "a2", "a3"],
                ),
            ];
            for (blocks, far_tile_block, windows, sound) in layouts {
                let dir = tempfile::tempdir().unwrap();
                let data = dir.path().join("limit_row.tsv.gz");
                let blocks: Vec<&[u8]> = blocks.iter().map(String::as_bytes).collect();
                let offsets = write_blocks(&data, &blocks);
                write_tbi_bins(
                    &data,
                    vec![(
                        0,
                        vec![Chunk::new(position(0), position(offsets[blocks.len()]))],
                    )],
                    vec![position(0), position(0), position(offsets[far_tile_block])],
                );
                for ann in [
                    open(data.clone(), None),
                    open(data.clone(), None).reading_rows_as_tabix_indexes_them(),
                ] {
                    for asked in [windows.clone(), vec![(1, 50_000)]] {
                        let case = format!("{limit} in {} blocks, {asked:?}", blocks.len());
                        let mut by_query = Vec::new();
                        for &(start, end) in &asked {
                            by_query.extend(ann.query("21", start, end).unwrap());
                        }
                        assert_eq!(names(by_query), sound, "{case}");
                        assert_eq!(
                            names(ann.query_windows("21", &asked).unwrap()),
                            sound,
                            "{case}"
                        );
                        for threads in [1, 4] {
                            let by_blocks = pool(threads)
                                .install(|| ann.query_windows_parallel("21", &asked))
                                .unwrap();
                            assert_eq!(names(by_blocks), sound, "{case}, {threads} threads");
                        }
                        assert_eq!(
                            names(ann.query_windows_in_slices("21", &asked, 1).unwrap()),
                            sound,
                            "{case}"
                        );
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod key_position_tests {
    use super::allele_match::{self, KeyPosition};
    use super::TabixRecord;
    use vep_core::variant::InputVariant;

    fn row(start: u64, cell: &str, reference: &str, alt: &str) -> TabixRecord {
        TabixRecord {
            chr: "21".to_string(),
            start,
            end: start + reference.len() as u64 - 1,
            columns: vec!["21".into(), cell.into(), reference.into(), alt.into()],
            header: None,
        }
    }

    /// The key of a row whose alleles lose no leading base holds the cell as
    /// text (`"$ref\_$alt\_$pos"`, Utils/Sequence.pm:1232), equal to the
    /// variant's start only as its plain decimal: a fraction, an exponent, a
    /// trailing character, a sign, a leading zero or a blank is another text,
    /// and a Perl-false cell matches nothing (`:1173`).
    #[test]
    fn an_untrimmed_position_cell_is_compared_as_text() {
        let snv = InputVariant::new(
            "21".into(),
            30396092,
            30396092,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let matched = |cell: &str| {
            allele_match::matched_variant_alleles(&row(30396092, cell, "C", "A"), &snv, 1, 2, 3)
        };
        assert!(matched("30396092"));
        for cell in [
            "30396092.0",
            "30396092.5",
            "30396092x",
            "3.0396092e7",
            "abc",
            "+30396092",
            "030396092",
            " 30396092",
            "0",
            "",
        ] {
            assert!(!matched(cell), "{cell:?}");
        }
        assert!(KeyPosition::Cell {
            text: "100",
            moved: 0
        }
        .equals(KeyPosition::Number(100)));
        assert!(!KeyPosition::Cell {
            text: "100.0",
            moved: 0
        }
        .equals(KeyPosition::Number(100)));
    }

    /// A leading base shared by the row's alleles is a `$start++` on a scalar
    /// `$start + (length($ref) - 1)` numified (`:971, :1013`), so the key holds
    /// the cell's number plus the bases moved: a fraction of `.0`, an exponent
    /// and a trailing character reach the deletion's start, a fraction of `.5`
    /// and a cell nothing numeric leads do not.
    #[test]
    fn a_position_cell_moved_by_a_shared_base_is_compared_as_its_number() {
        let deletion = InputVariant::new(
            "21".into(),
            30881494,
            30881494,
            b"A".to_vec(),
            b"-".to_vec(),
        );
        let matched = |cell: &str| {
            allele_match::matched_variant_alleles(
                &row(30881493, cell, "GA", "G"),
                &deletion,
                1,
                2,
                3,
            )
        };
        for cell in [
            "30881493",
            "30881493.0",
            "3.0881493e7",
            "30881493x",
            " 30881493",
        ] {
            assert!(matched(cell), "{cell:?}");
        }
        for cell in ["30881493.5", "abc", "30881494", "0", ""] {
            assert!(!matched(cell), "{cell:?}");
        }
        // A Perl-false cell returns the match before an undefined allele would
        // fail it (`:1173` before `:1230`); a Perl-true one still fails it.
        let short = |cell: &str| TabixRecord {
            chr: "21".to_string(),
            start: 30881493,
            end: 30881493,
            columns: vec!["21".into(), cell.into(), "GA".into()],
            header: None,
        };
        assert_eq!(
            allele_match::matched_or_undefined(&short("0"), &deletion, 1, 2, 3),
            Some(false)
        );
        assert_eq!(
            allele_match::matched_or_undefined(&short("30881493"), &deletion, 1, 2, 3),
            None
        );
        let cell = KeyPosition::Cell {
            text: "30881493.0",
            moved: 0,
        };
        assert!(cell.moved().equals(KeyPosition::Number(30881494)));
        assert!(!cell.equals(KeyPosition::Number(30881493)));
        assert_eq!(KeyPosition::Number(5).moved(), KeyPosition::Number(6));
    }

    /// `$start++` on a cell nothing has numified (gnomADMt.pm:310): digits and
    /// letters step as a string, anything else as a number plus one.
    #[test]
    fn perl_increment_steps_a_string_as_perl_does() {
        for (cell, incremented) in [
            ("21", "22"),
            ("09", "10"),
            ("99", "100"),
            ("abc", "abd"),
            ("Az", "Ba"),
            ("zz", "aaa"),
            ("Zz", "AAa"),
            ("a9", "b0"),
            ("", "1"),
            ("21.0", "22"),
            ("42.5", "43.5"),
            ("6x", "7"),
            ("3.0881493e7", "30881494"),
            ("ab c", "1"),
            ("-1", "0"),
        ] {
            assert_eq!(allele_match::perl_increment(cell), incremented, "{cell:?}");
        }
    }
}
