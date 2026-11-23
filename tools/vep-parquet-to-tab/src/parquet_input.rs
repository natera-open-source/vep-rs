// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Parquet output: its files in path order, their footers, and the row groups as the
//! units the rows are printed by.

use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{Read, Seek, SeekFrom};
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::FileExt;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use arrow_schema::FieldRef;
use bytes::{Buf, Bytes};
use parquet::arrow::arrow_reader::{
    ArrowReaderMetadata, ArrowReaderOptions, ParquetRecordBatchReaderBuilder,
};
use parquet::arrow::ProjectionMask;
use parquet::errors::ParquetError;
use parquet::file::reader::{ChunkReader, Length};
use rayon::prelude::*;

use crate::cells::{stores_text, unsupported_type, Column};
use crate::Failure;

const KEY_COLUMNS: [&str; 6] = ["contig_index", "chrom", "pos", "end", "ref", "alt"];
/// What the writer stores a tab column named as a key column behind.
const TAB_COLUMN_PREFIX: &str = "tab_";
/// The footer key under which a file states its `--tab` column names, joined by tabs.
const TAB_COLUMNS_KEY: &str = "tab_columns";
const PARQUET_MAGIC: &[u8; 4] = b"PAR1";
/// The magic and the footer length that end a Parquet file.
const FOOTER_TAIL: u64 = 8;
/// Rows decoded at a time; a row group's text grows batch by batch.
const BATCH_ROWS: usize = 8192;

/// One file of the output: its path as messages print it, its footer and the columns it
/// prints.
struct PartFile {
    path: PathBuf,
    file: File,
    /// The file's length in bytes.
    length: u64,
    metadata: ArrowReaderMetadata,
    /// The value columns of the tab columns: a batch read under it holds the tab columns
    /// in their order.
    projection: ProjectionMask,
    /// Whether the field that prints under each tab column stores strings.
    prints_text: Vec<bool>,
}

/// A Parquet output opened for printing.
pub struct ParquetOutput {
    files: Vec<PartFile>,
    /// The tab columns: every field after the leading key columns, under the name the
    /// `--tab` format gives it (`tab_columns`).
    pub columns: Vec<String>,
    /// The longest text a row group has printed so far, which sizes the next one's.
    longest_text: AtomicUsize,
}

/// The bytes of one row group, read from its file in one piece; its pages are served from
/// memory.
struct RowGroupBytes {
    /// Where `bytes` start in the file.
    start: u64,
    bytes: Bytes,
    file_length: u64,
}

impl RowGroupBytes {
    fn range(&self, start: u64, length: usize) -> parquet::errors::Result<std::ops::Range<usize>> {
        start
            .checked_sub(self.start)
            .and_then(|from| usize::try_from(from).ok())
            .and_then(|from| Some(from..from.checked_add(length)?))
            .filter(|range| range.end <= self.bytes.len())
            .ok_or_else(|| ParquetError::General("a page lies outside its row group".into()))
    }
}

impl Length for RowGroupBytes {
    fn len(&self) -> u64 {
        self.file_length
    }
}

impl ChunkReader for RowGroupBytes {
    type T = bytes::buf::Reader<Bytes>;

    fn get_read(&self, start: u64) -> parquet::errors::Result<Self::T> {
        let from = self.range(start, 0)?.start;
        Ok(self.bytes.slice(from..).reader())
    }

    fn get_bytes(&self, start: u64, length: usize) -> parquet::errors::Result<Bytes> {
        Ok(self.bytes.slice(self.range(start, length)?))
    }
}

/// One row group of one file.
#[derive(Clone, Copy)]
pub struct Unit {
    file: usize,
    row_group: usize,
}

/// The entries under `directory` whose name ends in `.parquet` and that are files, as paths
/// relative to it split into their components. Directories are walked at any depth, hidden
/// ones included; a symbolic link to a directory is not followed, one to a file is a file;
/// a directory that cannot be listed holds nothing.
fn walk(directory: &Path, prefix: &[Vec<u8>], found: &mut Vec<Vec<Vec<u8>>>) {
    let Ok(entries) = fs::read_dir(directory) else {
        return;
    };
    for entry in entries.flatten() {
        let name = entry.file_name();
        let mut components = prefix.to_vec();
        components.push(name.as_bytes().to_vec());
        let path = entry.path();
        if entry.file_type().is_ok_and(|kind| kind.is_dir()) {
            walk(&path, &components, found);
        } else if name.as_bytes().ends_with(b".parquet") && path.is_file() {
            found.push(components);
        }
    }
}

/// The Parquet files of the output at `path`, which is a regular file or a directory: the
/// file itself, or every `.parquet` file under the directory ordered component by
/// component, so that `a/z.parquet` sorts before `a-c/x.parquet`.
fn parquet_files(path: &Path) -> Result<Vec<PathBuf>, Failure> {
    if path.is_file() {
        return Ok(vec![path.to_path_buf()]);
    }
    let mut found = Vec::new();
    walk(path, &[], &mut found);
    if found.is_empty() {
        return Err(Failure::at(path, "no .parquet file under the directory"));
    }
    found.sort();
    // Joined as pathlib joins: nothing is put before a part of the directory `.`.
    let base = (path.as_os_str() != ".").then_some(path);
    Ok(found
        .into_iter()
        .map(|components| {
            let relative = components.join(&b'/');
            let relative = Path::new(OsStr::from_bytes(&relative));
            base.map_or_else(|| relative.to_path_buf(), |base| base.join(relative))
        })
        .collect())
}

/// Opens one file's footer. The first three refusals carry the wording of the Arrow C++
/// reader, which a caller of this tool may already match on.
fn open_footer(path: &Path) -> Result<(File, u64, ArrowReaderMetadata), Failure> {
    let mut file = File::open(path).map_err(|error| {
        let errno = error.raw_os_error().unwrap_or(0);
        let text = error.to_string();
        let detail = text
            .strip_suffix(&format!(" (os error {errno})"))
            .unwrap_or(&text);
        let mut line = format!("[Errno {errno}] Failed to open local file '").into_bytes();
        line.extend_from_slice(path.as_os_str().as_bytes());
        line.extend_from_slice(format!("'. Detail: [errno {errno}] {detail}").as_bytes());
        Failure::Reported(line)
    })?;
    let size = file.metadata().map_err(|e| Failure::io(&e))?.len();
    if size == 0 {
        return Err(Failure::message("Parquet file size is 0 bytes"));
    }
    if size < FOOTER_TAIL {
        return Err(Failure::message(format!(
            "Parquet file size is {size} bytes, smaller than the minimum file footer (8 bytes)"
        )));
    }
    let mut tail = [0u8; FOOTER_TAIL as usize];
    file.seek(SeekFrom::Start(size - FOOTER_TAIL))
        .and_then(|_| file.read_exact(&mut tail))
        .map_err(|e| Failure::io(&e))?;
    if &tail[4..] != PARQUET_MAGIC {
        return Err(Failure::message(
            "Parquet magic bytes not found in footer. Either the file is corrupted or this is \
             not a parquet file.",
        ));
    }
    let footer = u64::from(u32::from_le_bytes([tail[0], tail[1], tail[2], tail[3]]));
    if size < footer + FOOTER_TAIL {
        return Err(Failure::message(format!(
            "Parquet file size is {size} bytes, smaller than the size reported by footer's \
             ({footer}bytes)"
        )));
    }
    // The types come from the Parquet schema alone, so that a string, a list or a number
    // reads as one Arrow type whatever Arrow type its writer declared for it.
    let options = ArrowReaderOptions::new().with_skip_arrow_metadata(true);
    let metadata =
        ArrowReaderMetadata::load(&file, options).map_err(|error| Failure::at(path, error))?;
    Ok((file, size, metadata))
}

/// Whether `name` is a key column's name behind any number of `TAB_COLUMN_PREFIX`, none
/// included.
fn names_a_key_column(name: &str) -> bool {
    let mut rest = name;
    while let Some(inner) = rest.strip_prefix(TAB_COLUMN_PREFIX) {
        rest = inner;
    }
    KEY_COLUMNS.contains(&rest)
}

/// The tab column a field after the key columns holds, by the name it is stored under, as
/// `tab_column_name` of the writer's module reads it. The writer stores a tab column named
/// as a key column, behind any number of the prefix, behind one more (`ref` as `tab_ref`,
/// `tab_ref` as `tab_tab_ref`), so that no two columns of a file share a name; every other
/// column is stored under its own where no earlier column carries that name in other
/// letters (`tab_columns`).
fn tab_column_name(stored: &str) -> &str {
    match stored.strip_prefix(TAB_COLUMN_PREFIX) {
        Some(tab_column) if names_a_key_column(tab_column) => tab_column,
        _ => stored,
    }
}

/// The tab columns of a file: every field after the leading key columns, under the name
/// the `--tab` format gives it.
///
/// The names are the ones the footer states under `tab_columns`, one for each such field:
/// the writer stores a name an earlier column carries in other letters behind the prefix
/// too, and such a name does not read back from the stored name. A file that states no
/// such list, or a list of another length, names each column by `tab_column_name` of its
/// field.
fn tab_columns(metadata: &ArrowReaderMetadata) -> Result<Vec<String>, Failure> {
    let fields = tab_fields(metadata);
    if fields.is_empty() {
        return Err(Failure::message("no tab columns in the Parquet schema"));
    }
    let stated: Vec<&str> = metadata
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .into_iter()
        .flatten()
        .find(|entry| entry.key == TAB_COLUMNS_KEY)
        .and_then(|entry| entry.value.as_deref())
        .map(|names| names.split('\t').collect())
        .unwrap_or_default();
    Ok(if stated.len() == fields.len() {
        stated.into_iter().map(str::to_string).collect()
    } else {
        fields
            .iter()
            .map(|field| tab_column_name(field.name()).to_string())
            .collect()
    })
}

/// The count of the key columns a file leads with.
fn leading_key_columns(metadata: &ArrowReaderMetadata) -> usize {
    metadata
        .schema()
        .fields()
        .iter()
        .take_while(|field| KEY_COLUMNS.contains(&field.name().as_str()))
        .count()
}

/// The fields a file's tab columns are read from, in tab-column order.
///
/// A tab column is read from its own field: the field at its place after the leading key
/// columns, whatever its name. So a part that holds a tab column under a key column's own
/// name, where the writer stores it behind `TAB_COLUMN_PREFIX`, prints the cells of that
/// column's own field.
fn tab_fields(metadata: &ArrowReaderMetadata) -> &[FieldRef] {
    &metadata.schema().fields()[leading_key_columns(metadata)..]
}

/// The value columns under a file's tab columns. A field below the root holds one value
/// column at least, since the reader takes a schema element without children for a value
/// column.
fn tab_value_columns(metadata: &ArrowReaderMetadata) -> Vec<usize> {
    let lead = leading_key_columns(metadata);
    let schema = metadata.parquet_schema();
    (0..schema.num_columns())
        .filter(|leaf| schema.get_column_root_idx(*leaf) >= lead)
        .collect()
}

impl ParquetOutput {
    /// Opens every file of the output at `path` and holds them to one list of tab columns.
    pub fn open(path: &Path) -> Result<Self, Failure> {
        let paths = parquet_files(path)?;
        let footers: Vec<Result<_, Failure>> =
            paths.par_iter().map(|path| open_footer(path)).collect();
        let footers = footers.into_iter().collect::<Result<Vec<_>, _>>()?;
        let mut files = Vec::with_capacity(paths.len());
        let mut columns: Option<Vec<String>> = None;
        for (path, (file, length, metadata)) in paths.iter().zip(footers) {
            let names = tab_columns(&metadata)?;
            match &columns {
                None => columns = Some(names),
                Some(first) if *first != names => {
                    let mut line = path.as_os_str().as_bytes().to_vec();
                    line.extend_from_slice(b": its tab columns differ from ");
                    line.extend_from_slice(paths[0].as_os_str().as_bytes());
                    line.extend_from_slice(b"'s");
                    return Err(Failure::Reported(line));
                }
                Some(_) => {}
            }
            files.push((path, file, length, metadata));
        }
        let columns = columns.unwrap_or_default();
        let mut parts = Vec::with_capacity(files.len());
        for (path, file, length, metadata) in files {
            let fields = tab_fields(&metadata);
            for field in fields {
                if let Some(data_type) = unsupported_type(field) {
                    return Err(Failure::at(
                        path,
                        format!(
                            "column {} is stored as {data_type}, which the cell contract has no \
                             tab text for",
                            field.name()
                        ),
                    ));
                }
            }
            let projection =
                ProjectionMask::leaves(metadata.parquet_schema(), tab_value_columns(&metadata));
            let prints_text = fields.iter().map(|field| stores_text(field)).collect();
            parts.push(PartFile {
                path: path.clone(),
                file,
                length,
                metadata,
                projection,
                prints_text,
            });
        }
        Ok(Self {
            files: parts,
            columns,
            longest_text: AtomicUsize::new(0),
        })
    }

    /// Whether the cell that `unit`'s file prints under each tab column is a stored string.
    pub fn prints_text(&self, unit: Unit) -> &[bool] {
        &self.files[unit.file].prints_text
    }

    /// The names of the tab columns the output stores as strings: every field of the name,
    /// in every file.
    pub fn text_columns(&self) -> Vec<&str> {
        // Every file names the tab columns of the output, each read from the field at its
        // place.
        let stored = || {
            self.files
                .iter()
                .flat_map(|part| self.columns.iter().zip(&part.prints_text))
        };
        let mut names: Vec<&str> = Vec::new();
        for (name, _) in stored().filter(|(_, text)| **text) {
            let all_text = stored().all(|(other, text)| other != name || *text);
            if all_text && !names.contains(&name.as_str()) {
                names.push(name);
            }
        }
        names
    }

    /// The row groups of the output in row order: file by file, group by group.
    pub fn units(&self) -> Vec<Unit> {
        self.files
            .iter()
            .enumerate()
            .flat_map(|(file, part)| {
                (0..part.metadata.metadata().num_row_groups())
                    .map(move |row_group| Unit { file, row_group })
            })
            .collect()
    }

    /// The bytes of a row group's columns, read from the file in one piece.
    fn row_group_bytes(part: &PartFile, row_group: usize) -> Result<RowGroupBytes, String> {
        let columns = part.metadata.metadata().row_group(row_group).columns();
        let mut span: Option<(u64, u64)> = None;
        for (leaf, column) in columns.iter().enumerate() {
            if !part.projection.leaf_included(leaf) {
                continue;
            }
            let start = column
                .dictionary_page_offset()
                .unwrap_or_else(|| column.data_page_offset());
            let (Ok(start), Ok(length)) = (
                u64::try_from(start),
                u64::try_from(column.compressed_size()),
            ) else {
                return Err("a column chunk has a negative offset or length".to_string());
            };
            let end = start
                .checked_add(length)
                .filter(|end| *end <= part.length)
                .ok_or("a column chunk ends past the end of the file")?;
            span = Some(span.map_or((start, end), |(first, last)| {
                (first.min(start), last.max(end))
            }));
        }
        let (start, end) = span.unwrap_or((0, 0));
        let mut bytes = vec![0u8; (end - start) as usize];
        part.file
            .read_exact_at(&mut bytes, start)
            .map_err(|error| error.to_string())?;
        Ok(RowGroupBytes {
            start,
            bytes: bytes.into(),
            file_length: part.length,
        })
    }

    /// Reads one row group and appends every row's `--tab` line to `out`, without a line
    /// end; after each row `row` is handed `out` and where that row's line starts in it.
    pub fn for_each_row(
        &self,
        unit: Unit,
        out: &mut Vec<u8>,
        mut row: impl FnMut(&mut Vec<u8>, usize) -> Result<(), Failure>,
    ) -> Result<(), Failure> {
        let part = &self.files[unit.file];
        let longest = self.longest_text.load(Ordering::Relaxed);
        out.reserve(longest + longest / 16);
        let before = out.len();
        // A reader that meets bytes it cannot follow may panic where it would report.
        let read = catch_unwind(AssertUnwindSafe(|| -> Result<(), Failure> {
            let failed = |error: &dyn std::fmt::Display| Failure::at(&part.path, error);
            let bytes = Self::row_group_bytes(part, unit.row_group).map_err(|e| failed(&e))?;
            let reader =
                ParquetRecordBatchReaderBuilder::new_with_metadata(bytes, part.metadata.clone())
                    .with_row_groups(vec![unit.row_group])
                    .with_projection(part.projection.clone())
                    .with_batch_size(BATCH_ROWS)
                    .build()
                    .map_err(|e| failed(&e))?;
            for batch in reader {
                let batch = batch.map_err(|e| failed(&e))?;
                // A batch holds the fields read in the file's own order: the tab columns.
                let columns: Vec<Column> = self
                    .columns
                    .iter()
                    .zip(batch.columns())
                    .map(|(name, column)| Column::of(name, column))
                    .collect();
                for index in 0..batch.num_rows() {
                    let start = out.len();
                    for (i, column) in columns.iter().enumerate() {
                        if i > 0 {
                            out.push(b'\t');
                        }
                        column.push(out, index);
                    }
                    row(out, start)?;
                }
            }
            Ok(())
        }));
        self.longest_text
            .fetch_max(out.len() - before, Ordering::Relaxed);
        read.unwrap_or_else(|_| {
            Err(Failure::at(
                &part.path,
                format!("row group {} cannot be read", unit.row_group),
            ))
        })
    }
}

/// Whether the regular file at `path` opens with the Parquet magic bytes.
pub fn starts_with_parquet_magic(file: &mut File) -> std::io::Result<bool> {
    let mut head = [0u8; 4];
    let mut filled = 0;
    while filled < head.len() {
        match file.read(&mut head[filled..])? {
            0 => break,
            n => filled += n,
        }
    }
    Ok(&head[..filled] == PARQUET_MAGIC)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A stored name reads back as the writer's `tab_column_name` reads it: one prefix
    /// goes where a key column's name stands behind it, and every other name stays.
    #[test]
    fn a_stored_name_reads_back_as_its_tab_column() {
        for key in KEY_COLUMNS {
            assert_eq!(tab_column_name(&format!("tab_{key}")), key);
            assert_eq!(
                tab_column_name(&format!("tab_tab_{key}")),
                format!("tab_{key}")
            );
            assert_eq!(tab_column_name(key), key);
        }
        for name in [
            "Allele",
            "REF_ALLELE",
            "Ref",
            "ALT",
            "reference",
            "ref_",
            "tab",
            "tab_",
            "tab_tab_",
            "tab_Allele",
            "tab_Ref",
            "TAB_ref",
            "xtab_ref",
            "",
        ] {
            assert_eq!(tab_column_name(name), name);
        }
    }

    /// The canonical form takes a column's kind from the field that prints it and, for a
    /// `--tab` output read beside the file, from the tab column's name: a string column
    /// stored as `tab_ref` is the tab column `ref` in the header and among the string
    /// columns, and an INT32 column is none of them.
    #[test]
    fn a_string_column_stored_behind_the_prefix_is_named_as_its_tab_column() {
        use arrow_array::{ArrayRef, Int32Array, RecordBatch, StringArray};
        use arrow_schema::{DataType, Field, Schema};
        use parquet::arrow::ArrowWriter;
        use std::sync::Arc;

        let dir = scratch("stored-text");
        let path = dir.join("part-00000.parquet");
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("Feature", DataType::Utf8, true),
            Field::new("tab_ref", DataType::Utf8, true),
            Field::new("TSL", DataType::Int32, true),
        ]));
        let columns: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from(vec!["21"])),
            Arc::new(StringArray::from(vec!["A"])),
            Arc::new(StringArray::from(vec!["ENST1"])),
            Arc::new(StringArray::from(vec!["T"])),
            Arc::new(Int32Array::from(vec![1])),
        ];
        let batch = RecordBatch::try_new(Arc::clone(&schema), columns).unwrap();
        let mut writer = ArrowWriter::try_new(File::create(&path).unwrap(), schema, None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();

        let output = ParquetOutput::open(&path).unwrap();
        assert_eq!(output.columns, ["Feature", "ref", "TSL"]);
        assert_eq!(output.text_columns(), ["Feature", "ref"]);
        let units = output.units();
        assert_eq!(output.prints_text(units[0]), [true, true, false]);
        let mut line = Vec::new();
        output
            .for_each_row(units[0], &mut line, |_, _| Ok(()))
            .unwrap();
        assert_eq!(line, b"ENST1\tT\t1");
        fs::remove_dir_all(&dir).unwrap();
    }

    fn scratch(name: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("vep-parquet-to-tab-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn the_files_of_a_directory_sort_component_by_component() {
        let dir = scratch("walk");
        for name in [
            "b.parquet",
            "a/z.parquet",
            "a-c/x.parquet",
            "a.d/x.parquet",
            ".hidden.parquet",
            ".parquet",
            ".hid/h.parquet",
            "UP.PARQUET",
            "deep/er/est/p.parquet",
            "dir.parquet/inner.parquet",
            "x.parquet.txt",
            "_vep_rs_manifest.json",
        ] {
            let path = dir.join(name);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, b"").unwrap();
        }
        let outside = scratch("walk-outside");
        fs::write(outside.join("o.parquet"), b"").unwrap();
        std::os::unix::fs::symlink(&outside, dir.join("linkdir")).unwrap();
        std::os::unix::fs::symlink(outside.join("o.parquet"), dir.join("linkfile.parquet"))
            .unwrap();
        std::os::unix::fs::symlink("nowhere", dir.join("dangling.parquet")).unwrap();
        let found: Vec<String> = parquet_files(&dir)
            .ok()
            .unwrap()
            .iter()
            .map(|path| path.strip_prefix(&dir).unwrap().display().to_string())
            .collect();
        assert_eq!(
            found,
            [
                ".hid/h.parquet",
                ".hidden.parquet",
                ".parquet",
                "a/z.parquet",
                "a-c/x.parquet",
                "a.d/x.parquet",
                "b.parquet",
                "deep/er/est/p.parquet",
                "dir.parquet/inner.parquet",
                "linkfile.parquet",
            ]
        );
        // A symbolic link to a directory is read when it is the output itself.
        let through_link = parquet_files(&dir.join("linkdir")).ok().unwrap();
        assert_eq!(through_link, [dir.join("linkdir").join("o.parquet")]);
        fs::remove_dir_all(&dir).unwrap();
        fs::remove_dir_all(&outside).unwrap();
    }
}
