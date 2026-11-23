// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The residual files of one engine.

use super::error::{Scored, ScorerError};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// The side of a tuple: the reference's or the engine's.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Side {
    /// Ordered first: a verdict line and a residual line carry `E` and `R`, and their
    /// order under one value is the bytes' order.
    Engine,
    Reference,
}

impl Side {
    /// The byte a verdict line carries for the side.
    pub fn letter(self) -> u8 {
        match self {
            Side::Engine => b'E',
            Side::Reference => b'R',
        }
    }

    /// The word a residual line opens with.
    pub fn name(self) -> &'static [u8] {
        match self {
            Side::Engine => b"engine",
            Side::Reference => b"reference",
        }
    }
}
/// A failed write of a residual file in the words CPython gives an `OSError` of a write,
/// `[Errno N] <reason>`, which names no file.
fn unnamed(error: &std::io::Error) -> ScorerError {
    ScorerError(super::paths::OsFailure::of(error, None).text())
}

/// The path of the residual file of one source, engine and projection.
pub fn residual_path(directory: &Path, source: &str, engine: &str, projection: &str) -> PathBuf {
    directory.join(format!("{source}.{engine}.{projection}.tsv"))
}

/// The residual lines of one projection over a run of keys: per source, in the order each
/// source's first line comes, one line of its side and its columns for every tuple left
/// unexplained.
#[derive(Debug, Default)]
pub struct ResidualLines(Vec<(String, Vec<u8>)>);

impl ResidualLines {
    /// Takes the `(value, side)` pairs left unexplained under `prefix`, the tuple's columns
    /// before its value, as lines of the source, by value then side.
    pub fn take(&mut self, source: &str, prefix: &[u8], unexplained: &mut [(&[u8], Side)]) {
        let at = match self.0.iter().position(|(held, _)| held == source) {
            Some(at) => at,
            None => {
                self.0.push((source.to_string(), Vec::new()));
                self.0.len() - 1
            }
        };
        unexplained.sort_unstable();
        let lines = &mut self.0[at].1;
        for (value, side) in unexplained.iter() {
            for part in [side.name(), b"\t", prefix, value, b"\n"] {
                lines.extend_from_slice(part);
            }
        }
    }

    /// Each source with its lines.
    pub fn by_source(&self) -> &[(String, Vec<u8>)] {
        &self.0
    }
}

/// The residual files of one engine: per source and projection, every tuple outside the
/// intersection that no applied class set aside, one line of its side and its columns.
pub struct Residuals {
    directory: PathBuf,
    engine: String,
    written: Vec<PathBuf>,
    files: Vec<ResidualFile>,
}

/// One residual file: `file` is `None` once it is ended.
struct ResidualFile {
    source: String,
    projection: String,
    file: Option<BufWriter<File>>,
}

impl Residuals {
    pub fn new(directory: &Path, engine: &str) -> Residuals {
        Residuals {
            directory: directory.to_path_buf(),
            engine: engine.to_string(),
            written: Vec::new(),
            files: Vec::new(),
        }
    }

    /// The paths of the files opened or left empty, in the order they were made.
    pub fn written(&self) -> &[PathBuf] {
        &self.written
    }

    fn path(&mut self, source: &str, projection: &str) -> PathBuf {
        let path = residual_path(&self.directory, source, &self.engine, projection);
        self.written.push(path.clone());
        path
    }

    /// Appends `lines`, whole residual lines, to the file of the source and projection,
    /// which is made at its first lines.
    pub fn append(&mut self, projection: &str, source: &str, lines: &[u8]) -> Scored<()> {
        let held = self
            .files
            .iter()
            .position(|held| held.source == source && held.projection == projection);
        let at = match held {
            Some(at) => at,
            None => {
                let path = self.path(source, projection);
                let file = File::create(&path).map_err(|error| ScorerError::os(&error, &path))?;
                self.files.push(ResidualFile {
                    source: source.to_string(),
                    projection: projection.to_string(),
                    file: Some(BufWriter::with_capacity(1 << 20, file)),
                });
                self.files.len() - 1
            }
        };
        let Some(file) = &mut self.files[at].file else {
            return Err(ScorerError(format!(
                "the residual file of {source} was closed before its last tuple"
            )));
        };
        file.write_all(lines).map_err(|error| unnamed(&error))
    }

    /// Ends every file opened, its last bytes written.
    pub fn close(&mut self) -> Scored<()> {
        let mut failed = None;
        for held in &mut self.files {
            if let Some(mut open) = held.file.take() {
                if let Err(error) = open.flush() {
                    failed.get_or_insert(unnamed(&error));
                }
            }
        }
        failed.map_or(Ok(()), Err)
    }

    /// Leaves an empty file for each scored `(source, projection)` without a residual.
    pub fn complete(&mut self, scored: &[(String, String)]) -> Scored<()> {
        for (source, projection) in scored {
            let held = self
                .files
                .iter()
                .any(|held| &held.source == source && &held.projection == projection);
            if !held {
                let path = self.path(source, projection);
                File::create(&path).map_err(|error| ScorerError::os(&error, &path))?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn scratch(name: &str) -> PathBuf {
        let directory = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-residuals-{name}-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&directory);
        fs::create_dir_all(&directory).unwrap();
        directory
    }

    #[test]
    fn unexplained_tuples_are_written_by_value_then_side_under_their_prefix() {
        let directory = scratch("lines");
        let mut residuals = Residuals::new(&directory, "engine");
        let mut lines = ResidualLines::default();
        let mut first: Vec<(&[u8], Side)> = vec![
            (b"b", Side::Reference),
            (b"a", Side::Reference),
            (b"b", Side::Engine),
        ];
        lines.take("Planted", b"1:100\tT\tENST1\tPLANT_score\t", &mut first);
        let mut other: Vec<(&[u8], Side)> = vec![(b"z", Side::Engine)];
        lines.take("Sown", b"1:100\tT\tENST1\tSOWN_score\t", &mut other);
        let mut second: Vec<(&[u8], Side)> = vec![(b"\xff", Side::Engine)];
        lines.take("Planted", b"1:200\tT\tENST1\tPLANT_score\t", &mut second);
        let sources: Vec<&str> = lines
            .by_source()
            .iter()
            .map(|(source, _)| source.as_str())
            .collect();
        assert_eq!(sources, ["Planted", "Sown"]);
        let planted = &lines.by_source()[0].1;
        // The lines of a source reach its file whatever the pieces they are appended in.
        let (opening, rest) = planted.split_at(planted.len() / 2);
        residuals.append("feature", "Planted", opening).unwrap();
        residuals.append("feature", "Planted", rest).unwrap();
        residuals.close().unwrap();
        residuals
            .complete(&[
                ("Planted".to_string(), "feature".to_string()),
                ("Planted".to_string(), "allele".to_string()),
            ])
            .unwrap();
        let feature = directory.join("Planted.engine.feature.tsv");
        let allele = directory.join("Planted.engine.allele.tsv");
        assert_eq!(
            fs::read(&feature).unwrap(),
            b"reference\t1:100\tT\tENST1\tPLANT_score\ta\n\
              engine\t1:100\tT\tENST1\tPLANT_score\tb\n\
              reference\t1:100\tT\tENST1\tPLANT_score\tb\n\
              engine\t1:200\tT\tENST1\tPLANT_score\t\xff\n"
        );
        assert_eq!(fs::read(&allele).unwrap(), b"");
        assert_eq!(residuals.written(), [feature, allele]);
        fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn file_that_cannot_be_opened_is_a_refusal_in_the_interpreters_words() {
        let directory = scratch("missing").join("not-made");
        let mut residuals = Residuals::new(&directory, "engine");
        let refused = residuals.append("feature", "Planted", b"");
        let path = directory.join("Planted.engine.feature.tsv");
        assert_eq!(
            refused,
            Err(ScorerError(format!(
                "[Errno 2] No such file or directory: '{}'",
                path.display()
            )))
        );
        fs::remove_dir_all(directory.parent().unwrap()).unwrap();
    }
}
