// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The output directory's manifest, written last, and the check and removal of
//! a previous run's output under `--force_overwrite`, which take only what
//! vep-rs wrote.

use std::fs::{self, File};
use std::io;
use std::path::{Path, PathBuf};

use parquet::file::reader::FileReader;
use parquet::file::serialized_reader::SerializedFileReader;
use serde::{Deserialize, Serialize};

use super::parquet_error;
use crate::error::IoError;

/// The manifest's name under the output directory; the leading underscore
/// keeps it out of a reader's `*.parquet` glob and of Spark's listing.
pub const MANIFEST_NAME: &str = "_vep_rs_manifest.json";

/// The footer key every part vep-rs writes carries.
const VERSION_KEY: &str = "vep_rs_version";

/// The directory prefix of a `--parquet_partition chrom` directory.
const PARTITION_PREFIX: &str = "chrom=";

/// The output directory's manifest: the run that wrote it and its parts, in
/// write order, with their row counts.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Manifest {
    /// The footer's `vep_rs_version`.
    pub vep_rs_version: String,
    /// The footer's `run_id`.
    pub run_id: String,
    /// The rows of every part together.
    pub rows: u64,
    pub parts: Vec<PartEntry>,
}

/// One part of the output: its path under the directory and its rows.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PartEntry {
    pub path: String,
    pub rows: u64,
}

impl Manifest {
    /// A manifest with no parts, its identity read from the footer keys.
    pub fn new(kv_metadata: &[(String, String)]) -> Self {
        let value = |key: &str| {
            kv_metadata
                .iter()
                .find(|(k, _)| k == key)
                .map(|(_, v)| v.clone())
                .unwrap_or_default()
        };
        Self {
            vep_rs_version: value(VERSION_KEY),
            run_id: value("run_id"),
            rows: 0,
            parts: Vec::new(),
        }
    }

    /// Writes the manifest into `dir`, the row count summed from the parts.
    pub fn write(&mut self, dir: &Path) -> Result<(), IoError> {
        self.rows = self.parts.iter().map(|part| part.rows).sum();
        let text = serde_json::to_string_pretty(self).map_err(parquet_error)?;
        fs::write(dir.join(MANIFEST_NAME), text + "\n")?;
        Ok(())
    }

    /// Reads the manifest of `dir`.
    pub fn read(dir: &Path) -> Result<Self, IoError> {
        let text = fs::read_to_string(dir.join(MANIFEST_NAME))?;
        serde_json::from_str(&text).map_err(parquet_error)
    }
}

/// Whether `path` is a Parquet file whose footer carries `vep_rs_version`,
/// which every part vep-rs writes does.
pub fn is_vep_rs_part(path: &Path) -> bool {
    let Ok(file) = File::open(path) else {
        return false;
    };
    let Ok(reader) = SerializedFileReader::new(file) else {
        return false;
    };
    reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .is_some_and(|kv| kv.iter().any(|entry| entry.key == VERSION_KEY))
}

/// Why `--force_overwrite` left a previous output in place.
#[derive(Debug, thiserror::Error)]
pub enum PreviousOutput {
    /// The directory holds an entry vep-rs did not write: a file that is not a
    /// part with vep-rs's footer or the manifest, or a directory that is not a
    /// partition directory of such parts. Nothing is removed.
    #[error("output directory contains files vep-rs did not write")]
    Foreign(PathBuf),
    /// The path is a file, not a directory of parts, and not a part vep-rs
    /// wrote.
    #[error("output path {} is a file vep-rs did not write", .0.display())]
    ForeignFile(PathBuf),
    /// The path is a symbolic link, whose target a run neither reads as its
    /// previous output nor writes to.
    #[error("output path {} is a symbolic link", .0.display())]
    Link(PathBuf),
    /// An entry of the directory is a symbolic link. Nothing is removed.
    #[error("output directory entry {} is a symbolic link", .0.display())]
    LinkInside(PathBuf),
    #[error("could not remove the previous output: {0}")]
    Io(#[from] io::Error),
}

/// What a previous run left at an output path, every entry one vep-rs wrote.
#[derive(Default)]
struct PreviousEntries {
    files: Vec<PathBuf>,
    /// The partition directories, each holding only entries of `files`.
    directories: Vec<PathBuf>,
}

/// One entry of an output directory as its path and its type, the type read
/// without following a symbolic link. A link is the refusal.
fn non_link_entry(
    entry: io::Result<fs::DirEntry>,
) -> Result<(PathBuf, fs::FileType), PreviousOutput> {
    let entry = entry?;
    let file_type = entry.file_type()?;
    let path = entry.path();
    if file_type.is_symlink() {
        return Err(PreviousOutput::LinkInside(path));
    }
    Ok((path, file_type))
}

/// Reads what is at `path` without following a symbolic link there or below
/// it: nothing, a part vep-rs wrote, or a directory of what vep-rs writes (the
/// parts whose footer carries `vep_rs_version`, the manifest and the partition
/// directories holding such parts). Any other entry is the refusal, as is a
/// `path` whose metadata cannot be read. A `path` ending in a separator names
/// a link's target, so the caller passes none.
fn previous_entries(path: &Path) -> Result<PreviousEntries, PreviousOutput> {
    let mut previous = PreviousEntries::default();
    let metadata = match fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(previous),
        Err(error) => return Err(PreviousOutput::Io(error)),
    };
    let file_type = metadata.file_type();
    if file_type.is_symlink() {
        return Err(PreviousOutput::Link(path.to_path_buf()));
    }
    if !file_type.is_dir() {
        if file_type.is_file() && is_vep_rs_part(path) {
            previous.files.push(path.to_path_buf());
            return Ok(previous);
        }
        return Err(PreviousOutput::ForeignFile(path.to_path_buf()));
    }
    for entry in fs::read_dir(path)? {
        let (entry, file_type) = non_link_entry(entry)?;
        if file_type.is_dir() {
            let partition = entry
                .file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.starts_with(PARTITION_PREFIX));
            if !partition {
                return Err(PreviousOutput::Foreign(entry));
            }
            for part in fs::read_dir(&entry)? {
                let (part, file_type) = non_link_entry(part)?;
                if file_type.is_file() && is_vep_rs_part(&part) {
                    previous.files.push(part);
                } else {
                    return Err(PreviousOutput::Foreign(part));
                }
            }
            previous.directories.push(entry);
        } else if file_type.is_file()
            && (entry.file_name().is_some_and(|name| name == MANIFEST_NAME)
                || is_vep_rs_part(&entry))
        {
            previous.files.push(entry);
        } else {
            return Err(PreviousOutput::Foreign(entry));
        }
    }
    Ok(previous)
}

/// Checks that a run under `--force_overwrite` may replace what is at `path`,
/// removing nothing: nothing there, a part vep-rs wrote, or a directory whose
/// every entry is one vep-rs writes. A symbolic link at `path` or inside the
/// directory is refused by name, as is any other entry. A `path` whose
/// metadata cannot be read (one under a regular file or an unreadable
/// directory, or a name longer than the file system takes) passes: it holds
/// nothing to replace, and creating the output there fails with the output's
/// own refusal.
pub fn check_previous_output(path: &Path) -> Result<(), PreviousOutput> {
    if fs::symlink_metadata(path).is_err() {
        return Ok(());
    }
    previous_entries(path).map(|_| ())
}

/// Removes a previous run's output at `path`: the parts whose footer carries
/// `vep_rs_version`, the manifest and the partition directories holding such
/// parts, leaving the directory empty; a file at `path` goes when it is such a
/// part. Any entry [`check_previous_output`] refuses is refused here before
/// anything is removed, as is a `path` whose metadata cannot be read.
pub fn remove_previous_output(path: &Path) -> Result<(), PreviousOutput> {
    let previous = previous_entries(path)?;
    for file in previous.files {
        fs::remove_file(file)?;
    }
    for directory in previous.directories {
        fs::remove_dir(directory)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::fields::FieldOptions;
    use crate::output::parquet::{InputOrder, ParquetDirectoryWriter, ParquetOutputFormatter};
    use std::sync::Arc;

    /// A run's output with no rows into `dir`: one empty part at the root and
    /// the manifest.
    fn previous_output(dir: &Path) {
        let formatter = ParquetOutputFormatter::new(FieldOptions::default(), vec![], vec![]);
        let kv = vec![(VERSION_KEY.to_string(), "0.0.0".to_string())];
        ParquetDirectoryWriter::create(
            dir,
            Arc::clone(formatter.schema()),
            kv,
            InputOrder::Unchecked,
            None,
        )
        .unwrap()
        .close()
        .unwrap();
    }

    fn names(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = fs::read_dir(dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }

    #[test]
    fn the_manifest_round_trips_and_sums_the_rows() {
        let dir = tempfile::tempdir().unwrap();
        let mut manifest = Manifest::new(&[
            ("run_id".to_string(), "r".to_string()),
            (VERSION_KEY.to_string(), "1.2.3".to_string()),
        ]);
        manifest.parts = vec![
            PartEntry {
                path: "part-00000.parquet".into(),
                rows: 2,
            },
            PartEntry {
                path: "chrom=22/part-00001.parquet".into(),
                rows: 3,
            },
        ];
        manifest.write(dir.path()).unwrap();
        let read = Manifest::read(dir.path()).unwrap();
        assert_eq!(read, manifest);
        assert_eq!(read.rows, 5);
        assert_eq!(
            (read.vep_rs_version.as_str(), read.run_id.as_str()),
            ("1.2.3", "r")
        );
        assert!(Manifest::new(&[]).vep_rs_version.is_empty());
    }

    #[test]
    fn a_previous_output_is_removed_whole() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out");
        previous_output(&out);
        fs::copy(
            out.join("part-00000.parquet"),
            out.join("part-00001.parquet"),
        )
        .unwrap();
        let hive = out.join("chrom=22");
        fs::create_dir(&hive).unwrap();
        fs::copy(
            out.join("part-00000.parquet"),
            hive.join("part-00002.parquet"),
        )
        .unwrap();
        assert!(is_vep_rs_part(&out.join("part-00000.parquet")));
        remove_previous_output(&out).unwrap();
        assert_eq!(names(&out), Vec::<String>::new());
    }

    /// A foreign file, a foreign directory, or a file under a part's name that
    /// is not a Parquet file refuses the removal, and every entry stays.
    #[test]
    fn a_foreign_entry_refuses_the_removal_and_nothing_goes() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out");
        previous_output(&out);
        let before = names(&out);
        for foreign in ["notes.txt", "part-00001.parquet"] {
            fs::write(out.join(foreign), b"not a part").unwrap();
            let error = remove_previous_output(&out).unwrap_err();
            assert_eq!(
                error.to_string(),
                "output directory contains files vep-rs did not write"
            );
            assert!(matches!(&error, PreviousOutput::Foreign(p) if p == &out.join(foreign)));
            let mut expected = before.clone();
            expected.push(foreign.to_string());
            expected.sort();
            assert_eq!(names(&out), expected);
            fs::remove_file(out.join(foreign)).unwrap();
        }
        fs::create_dir(out.join("scratch")).unwrap();
        assert!(matches!(
            remove_previous_output(&out).unwrap_err(),
            PreviousOutput::Foreign(_)
        ));
        fs::remove_dir(out.join("scratch")).unwrap();
        fs::create_dir(out.join("chrom=1")).unwrap();
        fs::write(out.join("chrom=1").join("x.txt"), b"").unwrap();
        assert!(matches!(
            remove_previous_output(&out).unwrap_err(),
            PreviousOutput::Foreign(p) if p == out.join("chrom=1").join("x.txt")
        ));
        assert_eq!(names(&out).len(), before.len() + 1);
    }

    /// The check passes what the removal would take and what is not there,
    /// refuses what the removal refuses, and removes nothing either way.
    #[test]
    fn the_check_refuses_what_the_removal_refuses_and_removes_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out");
        check_previous_output(&out).unwrap();
        remove_previous_output(&out).unwrap();
        assert!(!out.exists());
        previous_output(&out);
        let hive = out.join("chrom=22");
        fs::create_dir(&hive).unwrap();
        fs::copy(
            out.join("part-00000.parquet"),
            hive.join("part-00001.parquet"),
        )
        .unwrap();
        let before = names(&out);
        check_previous_output(&out).unwrap();
        fs::write(hive.join("notes.txt"), b"mine").unwrap();
        assert!(matches!(
            check_previous_output(&out).unwrap_err(),
            PreviousOutput::Foreign(p) if p == hive.join("notes.txt")
        ));
        assert_eq!(names(&out), before);
        assert_eq!(names(&hive), ["notes.txt", "part-00001.parquet"]);
    }

    /// A path whose metadata cannot be read, one under a regular file or one
    /// whose name is longer than the file system takes, passes the check and
    /// is the removal's error.
    #[test]
    fn a_path_that_cannot_be_read_passes_the_check_and_fails_the_removal() {
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("file");
        fs::write(&file, b"text\n").unwrap();
        for unreadable in [file.join("out"), dir.path().join("n".repeat(300))] {
            check_previous_output(&unreadable).unwrap();
            assert!(matches!(
                remove_previous_output(&unreadable).unwrap_err(),
                PreviousOutput::Io(_)
            ));
        }
        assert_eq!(names(dir.path()), ["file"]);
    }

    /// A symbolic link at the path, dangling or not, or inside the directory
    /// is refused by name, and what it points at keeps its files.
    #[cfg(unix)]
    #[test]
    fn a_symbolic_link_is_refused_by_name_and_its_target_keeps_its_files() {
        use std::os::unix::fs::symlink;
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("target");
        previous_output(&target);
        let kept = names(&target);

        let link = dir.path().join("out");
        symlink(&target, &link).unwrap();
        for refusal in [check_previous_output(&link), remove_previous_output(&link)] {
            let error = refusal.unwrap_err();
            assert_eq!(
                error.to_string(),
                format!("output path {} is a symbolic link", link.display())
            );
            assert!(matches!(error, PreviousOutput::Link(p) if p == link));
        }
        let dangling = dir.path().join("dangling");
        symlink(dir.path().join("nowhere"), &dangling).unwrap();
        assert!(matches!(
            check_previous_output(&dangling).unwrap_err(),
            PreviousOutput::Link(p) if p == dangling
        ));

        let out = dir.path().join("real");
        previous_output(&out);
        let part = target.join("part-00000.parquet");
        let partition = dir.path().join("partition");
        fs::create_dir(&partition).unwrap();
        fs::copy(&part, partition.join("part-00001.parquet")).unwrap();
        for (inside, points_at) in [
            (out.join("part-00001.parquet"), &part),
            (out.join("chrom=22"), &partition),
            (out.join("chrom=21").join("part-00002.parquet"), &part),
        ] {
            fs::create_dir_all(inside.parent().unwrap()).unwrap();
            symlink(points_at, &inside).unwrap();
            let error = remove_previous_output(&out).unwrap_err();
            assert_eq!(
                error.to_string(),
                format!(
                    "output directory entry {} is a symbolic link",
                    inside.display()
                )
            );
            assert!(matches!(error, PreviousOutput::LinkInside(p) if p == inside));
            assert!(out.join("part-00000.parquet").exists());
            assert!(out.join(MANIFEST_NAME).exists());
            fs::remove_file(&inside).unwrap();
        }
        assert_eq!(names(&target), kept);
        assert_eq!(names(&partition), ["part-00001.parquet"]);
    }

    /// A file at the path goes when it is a part vep-rs wrote and refuses the
    /// removal otherwise.
    #[test]
    fn a_file_at_the_path_goes_only_when_vep_rs_wrote_it() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out");
        previous_output(&out);
        let own = dir.path().join("own.parquet");
        fs::rename(out.join("part-00000.parquet"), &own).unwrap();
        remove_previous_output(&own).unwrap();
        assert!(!own.exists());
        let text = dir.path().join("out.txt");
        fs::write(&text, b"## a text output\n").unwrap();
        let error = remove_previous_output(&text).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "output path {} is a file vep-rs did not write",
                text.display()
            )
        );
        assert!(text.exists());
    }
}
