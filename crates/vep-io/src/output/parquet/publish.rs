// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Parquet output written beside its destination and moved into place once
//! the run has succeeded, so the destination holds a completed run or what it
//! held before the run.

use std::ffi::{OsStr, OsString};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use super::manifest::{check_previous_output, remove_previous_output, PreviousOutput};

/// The longest file name the common file systems take, in bytes.
const NAME_MAX: usize = 255;

/// Why a run's output did not replace what was at its destination.
#[derive(Debug, thiserror::Error)]
pub enum PublishError {
    /// Something is at the destination that the run may not replace. The
    /// destination is as it was.
    #[error("the output path exists")]
    Exists,
    /// What is at the destination is not a previous run's output. The
    /// destination is as it was.
    #[error(transparent)]
    Previous(#[from] PreviousOutput),
    /// A rename failed. The destination is as it was.
    #[error("could not move the output into place: {0}")]
    Io(#[source] io::Error),
    /// The output could not be moved into place and the previous output could
    /// not be moved back: nothing is at the destination and the previous output
    /// is whole at `previous`.
    #[error(
        "could not move the output into place ({error}); the previous output is whole at {}",
        previous.display()
    )]
    Unrestored {
        previous: PathBuf,
        #[source]
        error: io::Error,
    },
    /// The output is in place and complete; the previous output could not be
    /// removed, and what is left of it is at `previous`.
    #[error(
        "the output is complete, but the previous output is left at {}: {reason}",
        previous.display()
    )]
    Unremoved {
        previous: PathBuf,
        #[source]
        reason: PreviousOutput,
    },
}

/// The directories of one run's output: the destination the run names, the
/// temporary directory beside it that the run writes, and the name the previous
/// output takes while the temporary directory is moved into place. Both names
/// carry the destination's file name and the run's identifier, in the
/// destination's own directory, so each move is a rename within one directory
/// and two runs never share a name; a destination that cannot be renamed (a
/// mount point) ends a run under `replace` at publication with
/// [`PublishError::Io`], the destination as it was. A run that ends before
/// [`Self::publish`] has moved the temporary directory removes that directory
/// and leaves the destination alone: through [`Self::discard`], which returns
/// the removal's error, or when the value is dropped, which cannot.
#[derive(Debug)]
pub struct UnpublishedOutput {
    /// The output path as the run names it, which a trailing separator makes
    /// read differently from `destination`: a regular file is then not a
    /// directory, and a symbolic link is its target.
    output: PathBuf,
    destination: PathBuf,
    temporary: PathBuf,
    previous: PathBuf,
    /// Whether the temporary directory exists and is this run's to remove.
    created: bool,
}

impl UnpublishedOutput {
    /// The directories of a run writing to `output` under the identifier
    /// `run_id`, none of them created; `None` when `output` has no final name
    /// to place them beside (an empty path, `.`, `..` or a root).
    pub fn beside(output: &Path, run_id: &str) -> Option<Self> {
        let name = output.file_name()?;
        let parent = output.parent().unwrap_or(Path::new(""));
        Some(Self {
            output: output.to_path_buf(),
            destination: parent.join(name),
            temporary: parent.join(sibling_name(name, "tmp", run_id)),
            previous: parent.join(sibling_name(name, "previous", run_id)),
            created: false,
        })
    }

    /// The directory the run writes its parts and manifest into.
    pub fn path(&self) -> &Path {
        &self.temporary
    }

    /// The error of reading the output path as the run names it; `None` when
    /// the path reads or nothing is there.
    fn unreadable(&self) -> Option<io::Error> {
        fs::symlink_metadata(&self.output)
            .err()
            .filter(|error| error.kind() != io::ErrorKind::NotFound)
    }

    /// Checks, removing nothing, that the run may replace what is at the
    /// destination ([`check_previous_output`]). An output path that cannot be
    /// read as the run names it, for any reason but its absence (a regular
    /// file named with a trailing separator, a name longer than the file
    /// system takes), passes: it holds nothing to replace, and
    /// [`Self::create`] refuses it.
    pub fn check_destination(&self) -> Result<(), PreviousOutput> {
        if self.unreadable().is_some() {
            return Ok(());
        }
        check_previous_output(&self.destination)
    }

    /// Creates the temporary directory, which must not exist: a directory
    /// another run left is never taken over. An output path that cannot be
    /// read as the run names it, for any reason but its absence, is the
    /// error, before anything is created: the temporary name is made from the
    /// destination's file name, cut to fit, so its creation would succeed
    /// beside a regular file named with a trailing separator and beside a
    /// name longer than the file system takes.
    pub fn create(&mut self) -> io::Result<()> {
        if let Some(error) = self.unreadable() {
            return Err(error);
        }
        fs::create_dir(&self.temporary)?;
        self.created = true;
        Ok(())
    }

    /// Moves the temporary directory to the destination. With nothing there it
    /// is one rename. Otherwise, when `replace` allows it and the destination
    /// still passes [`check_previous_output`], the temporary directory takes
    /// the permission bits of a directory there, the previous output is
    /// renamed aside, the temporary directory renamed in, and the previous
    /// output removed through [`remove_previous_output`], which takes only
    /// what vep-rs wrote; a failed second rename moves the previous output
    /// back. Every error but [`PublishError::Unremoved`] and
    /// [`PublishError::Unrestored`] leaves the destination as it was, and
    /// every error but the first of those leaves the temporary directory for
    /// [`Self::discard`] or the drop to remove.
    pub fn publish(&mut self, replace: bool) -> Result<(), PublishError> {
        let previous = match fs::symlink_metadata(&self.destination) {
            Ok(metadata) => metadata,
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                fs::rename(&self.temporary, &self.destination).map_err(PublishError::Io)?;
                self.created = false;
                return Ok(());
            }
            Err(error) => return Err(PublishError::Io(error)),
        };
        if !replace {
            return Err(PublishError::Exists);
        }
        check_previous_output(&self.destination)?;
        let previous_is_directory = previous.is_dir();
        if previous_is_directory {
            set_directory_permissions(&self.temporary, previous.permissions())
                .map_err(PublishError::Io)?;
        }
        fs::rename(&self.destination, &self.previous).map_err(PublishError::Io)?;
        if let Err(error) = fs::rename(&self.temporary, &self.destination) {
            return Err(match fs::rename(&self.previous, &self.destination) {
                Ok(()) => PublishError::Io(error),
                Err(_) => PublishError::Unrestored {
                    previous: self.previous.clone(),
                    error,
                },
            });
        }
        self.created = false;
        remove_previous_output(&self.previous)
            .and_then(|()| {
                if previous_is_directory {
                    fs::remove_dir(&self.previous)?;
                }
                Ok(())
            })
            .map_err(|reason| PublishError::Unremoved {
                previous: self.previous.clone(),
                reason,
            })
    }

    /// Removes the temporary directory of a run that did not publish it,
    /// whatever it holds, a part without a footer included, and returns the
    /// removal's error; a directory already moved to the destination, or never
    /// created, is left alone.
    pub fn discard(mut self) -> io::Result<()> {
        if std::mem::take(&mut self.created) {
            fs::remove_dir_all(&self.temporary)?;
        }
        Ok(())
    }
}

impl Drop for UnpublishedOutput {
    fn drop(&mut self) {
        if self.created {
            let _ = fs::remove_dir_all(&self.temporary);
        }
    }
}

/// Gives the directory at `path` the bits of `permissions` without following
/// a symbolic link: a `path` that is not a directory itself is the error, and
/// what a link there points at keeps its bits. On Unix the bits go on a
/// handle, and only when the handle is on the directory `path` was read as: a
/// link put there after that read opens as what it points at.
fn set_directory_permissions(path: &Path, permissions: fs::Permissions) -> io::Result<()> {
    let named = fs::symlink_metadata(path)?;
    if !named.is_dir() {
        return Err(io::ErrorKind::NotADirectory.into());
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        let directory = fs::File::open(path)?;
        let opened = directory.metadata()?;
        if (opened.dev(), opened.ino()) != (named.dev(), named.ino()) {
            return Err(io::ErrorKind::NotADirectory.into());
        }
        directory.set_permissions(permissions)
    }
    #[cfg(not(unix))]
    fs::set_permissions(path, permissions)
}

/// The name of a directory beside an output named `name`:
/// `.<name>.vep-rs-<role>-<run_id>`. The leading dot keeps it out of a shell
/// glob over the output's directory, and out of the listing of a dataset
/// reader that skips hidden names, while the run writes it. A `name` too long
/// to leave room for the rest is cut at a character boundary, the run's
/// identifier keeping the result unique.
fn sibling_name(name: &OsStr, role: &str, run_id: &str) -> OsString {
    let suffix = format!(".vep-rs-{role}-{run_id}");
    let room = NAME_MAX.saturating_sub(1 + suffix.len());
    let mut sibling = OsString::from(".");
    if name.len() <= room {
        sibling.push(name);
    } else {
        let lossy = name.to_string_lossy();
        let mut end = room;
        while !lossy.is_char_boundary(end) {
            end -= 1;
        }
        sibling.push(&lossy[..end]);
    }
    sibling.push(suffix);
    sibling
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::fields::FieldOptions;
    use crate::output::parquet::{
        InputOrder, Manifest, ParquetDirectoryWriter, ParquetOutputFormatter, MANIFEST_NAME,
    };
    use std::sync::Arc;

    /// A run's output with no rows written into `dir`, its manifest naming
    /// the run `run_id`.
    fn write_output(dir: &Path, run_id: &str) {
        let formatter = ParquetOutputFormatter::new(FieldOptions::default(), vec![], vec![]);
        let kv = vec![
            ("vep_rs_version".to_string(), "0.0.0".to_string()),
            ("run_id".to_string(), run_id.to_string()),
        ];
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

    /// An output of the run `run_id` beside `destination`, written.
    fn unpublished(destination: &Path, run_id: &str) -> UnpublishedOutput {
        let mut unpublished = UnpublishedOutput::beside(destination, run_id).unwrap();
        unpublished.create().unwrap();
        write_output(unpublished.path(), run_id);
        unpublished
    }

    /// The entries of `dir` and, for the directories among them, of each, as
    /// paths under `dir` with their bytes, sorted.
    fn contents(dir: &Path) -> Vec<(String, Vec<u8>)> {
        let mut entries = Vec::new();
        for entry in fs::read_dir(dir).unwrap() {
            let entry = entry.unwrap();
            let name = entry.file_name().to_string_lossy().into_owned();
            if entry.file_type().unwrap().is_dir() {
                for (inner, bytes) in contents(&entry.path()) {
                    entries.push((format!("{name}/{inner}"), bytes));
                }
            } else {
                entries.push((name, fs::read(entry.path()).unwrap()));
            }
        }
        entries.sort();
        entries
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
    fn the_temporary_and_previous_names_carry_the_output_name_and_the_run() {
        let unpublished =
            UnpublishedOutput::beside(Path::new("results/out.parquet/"), "1a-2b-3c").unwrap();
        assert_eq!(unpublished.destination, Path::new("results/out.parquet"));
        assert_eq!(
            unpublished.path(),
            Path::new("results/.out.parquet.vep-rs-tmp-1a-2b-3c")
        );
        assert_eq!(
            unpublished.previous,
            Path::new("results/.out.parquet.vep-rs-previous-1a-2b-3c")
        );
        let bare = UnpublishedOutput::beside(Path::new("out"), "r").unwrap();
        assert_eq!(bare.destination, Path::new("out"));
        assert_eq!(bare.path(), Path::new(".out.vep-rs-tmp-r"));
        for nameless in ["", ".", "..", "/", "out/.."] {
            assert!(
                UnpublishedOutput::beside(Path::new(nameless), "r").is_none(),
                "{nameless:?}"
            );
        }
    }

    /// A name at the file system's limit is cut to leave room for the run's
    /// identifier, at a character boundary, and two runs still differ.
    #[test]
    fn a_long_output_name_is_cut_to_fit_the_temporary_name() {
        let long = "é".repeat(127) + "x";
        assert_eq!(long.len(), NAME_MAX);
        let name = sibling_name(OsStr::new(&long), "tmp", "1a-2b-3c");
        let name = name.to_str().unwrap();
        assert_eq!(name.len(), NAME_MAX);
        assert!(name.starts_with(".é") && name.ends_with("é.vep-rs-tmp-1a-2b-3c"));
        assert_ne!(
            sibling_name(OsStr::new(&long), "tmp", "1a-2b-3d"),
            OsString::from(name)
        );
        // An identifier one byte shorter leaves room that ends inside a
        // two-byte character, which is dropped whole.
        let odd = sibling_name(OsStr::new(&long), "tmp", "1a-2b-3");
        assert_eq!(odd.to_str().unwrap().len(), NAME_MAX - 1);
        let short = sibling_name(OsStr::new("out.parquet"), "previous", "r");
        assert_eq!(short, OsString::from(".out.parquet.vep-rs-previous-r"));
    }

    /// A destination name longer than the file system takes fails the
    /// creation, though the temporary name, cut to fit, could be created
    /// beside it, and no directory is made.
    #[test]
    fn a_destination_name_too_long_for_the_file_system_fails_the_creation() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("n".repeat(300));
        let mut unpublished = UnpublishedOutput::beside(&out, "r").unwrap();
        assert!(unpublished.path().file_name().unwrap().len() <= NAME_MAX);
        assert!(unpublished.create().is_err());
        assert_eq!(names(dir.path()), Vec::<String>::new());
    }

    /// A regular file named with a trailing separator is neither an output to
    /// replace nor a name to write beside, whether or not vep-rs wrote it: the
    /// check passes, the creation fails, no directory is made and the file
    /// keeps its bytes.
    #[test]
    fn a_file_named_with_a_trailing_separator_fails_the_creation() {
        let dir = tempfile::tempdir().unwrap();
        let written = dir.path().join("written");
        write_output(&written, "old");
        let part = dir.path().join("part.parquet");
        fs::copy(written.join("part-00000.parquet"), &part).unwrap();
        let text = dir.path().join("out.txt");
        fs::write(&text, b"text\n").unwrap();
        for file in [&part, &text] {
            let bytes = fs::read(file).unwrap();
            let mut unpublished = UnpublishedOutput::beside(&file.join(""), "r").unwrap();
            assert_eq!(unpublished.destination.as_os_str(), file.as_os_str());
            unpublished.check_destination().unwrap();
            assert!(unpublished.create().is_err());
            assert_eq!(names(dir.path()), ["out.txt", "part.parquet", "written"]);
            assert_eq!(fs::read(file).unwrap(), bytes);
        }
    }

    /// A symbolic link named with a trailing separator reads as its target,
    /// and the destination is still the link: a link to a directory, or to
    /// nothing, is refused by name, and a link to a regular file fails the
    /// creation as the file does. What the links point at stays as it was.
    #[cfg(unix)]
    #[test]
    fn a_link_named_with_a_trailing_separator_is_refused_by_name_or_fails_the_creation() {
        use std::os::unix::fs::symlink;
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("target");
        write_output(&target, "old");
        let before = contents(&target);
        let text = dir.path().join("out.txt");
        fs::write(&text, b"text\n").unwrap();

        let to_directory = dir.path().join("link");
        symlink(&target, &to_directory).unwrap();
        let to_nothing = dir.path().join("dangling");
        symlink(dir.path().join("nowhere"), &to_nothing).unwrap();
        for link in [&to_directory, &to_nothing] {
            let unpublished = UnpublishedOutput::beside(&link.join(""), "r").unwrap();
            let error = unpublished.check_destination().unwrap_err();
            assert_eq!(
                error.to_string(),
                format!("output path {} is a symbolic link", link.display())
            );
        }

        let to_file = dir.path().join("file-link");
        symlink(&text, &to_file).unwrap();
        let mut unpublished = UnpublishedOutput::beside(&to_file.join(""), "r").unwrap();
        unpublished.check_destination().unwrap();
        assert!(unpublished.create().is_err());

        assert_eq!(
            names(dir.path()),
            ["dangling", "file-link", "link", "out.txt", "target"]
        );
        assert_eq!(contents(&target), before);
        assert_eq!(fs::read(&text).unwrap(), b"text\n");
    }

    #[test]
    fn publishing_to_a_new_path_moves_the_temporary_directory_there() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        let mut unpublished = unpublished(&out, "new");
        let written = contents(unpublished.path());
        assert!(!out.exists(), "nothing is at the destination before");
        unpublished.publish(false).unwrap();
        assert_eq!(contents(&out), written);
        assert_eq!(names(dir.path()), ["out.parquet"]);
    }

    /// A run that ends before it publishes removes its temporary directory, a
    /// part without a footer included, and the destination is as it was.
    #[test]
    fn dropping_an_unpublished_output_removes_only_its_temporary_directory() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        write_output(&out, "old");
        let before = contents(&out);
        let failed = unpublished(&out, "failed");
        fs::write(failed.path().join("part-00001.parquet"), b"cut short").unwrap();
        drop(failed);
        assert_eq!(names(dir.path()), ["out.parquet"]);
        assert_eq!(contents(&out), before);

        let mut uncreated = UnpublishedOutput::beside(&out, "r").unwrap();
        let foreign = uncreated.path().to_path_buf();
        fs::create_dir(&foreign).unwrap();
        assert!(
            uncreated.create().is_err(),
            "an existing directory is not taken over"
        );
        drop(uncreated);
        assert!(foreign.exists(), "a directory the run did not make stays");
    }

    #[test]
    fn publishing_over_a_previous_output_replaces_it_whole() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        write_output(&out, "old");
        let partition = out.join("chrom=22");
        fs::create_dir(&partition).unwrap();
        fs::copy(
            out.join("part-00000.parquet"),
            partition.join("part-00001.parquet"),
        )
        .unwrap();
        let mut replacement = unpublished(&out, "new");
        let written = contents(replacement.path());
        replacement.publish(true).unwrap();
        assert_eq!(contents(&out), written);
        assert_eq!(Manifest::read(&out).unwrap().run_id, "new");
        assert_eq!(names(dir.path()), ["out.parquet"]);

        // A part vep-rs wrote at the path is replaced by the directory.
        let file = dir.path().join("part.parquet");
        fs::copy(out.join("part-00000.parquet"), &file).unwrap();
        unpublished(&file, "over-a-file").publish(true).unwrap();
        assert_eq!(Manifest::read(&file).unwrap().run_id, "over-a-file");
        assert_eq!(names(dir.path()), ["out.parquet", "part.parquet"]);
    }

    /// A run that did not publish removes its temporary directory through
    /// `discard`, which returns the error of a removal the directory's parent
    /// does not allow and leaves the destination alone either way.
    #[cfg(unix)]
    #[test]
    fn discarding_an_output_removes_its_temporary_directory_or_returns_the_error() {
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        let parent = dir.path().join("results");
        fs::create_dir(&parent).unwrap();
        let out = parent.join("out.parquet");
        write_output(&out, "old");
        let before = contents(&out);
        unpublished(&out, "failed").discard().unwrap();
        assert_eq!(names(&parent), ["out.parquet"]);

        let kept = unpublished(&out, "kept");
        let temporary = kept.path().to_path_buf();
        fs::set_permissions(&parent, fs::Permissions::from_mode(0o555)).unwrap();
        if fs::write(parent.join("probe"), b"").is_ok() {
            // This user's writes ignore the mode, so the removal cannot fail.
            fs::set_permissions(&parent, fs::Permissions::from_mode(0o755)).unwrap();
            return;
        }
        assert!(kept.discard().is_err());
        fs::set_permissions(&parent, fs::Permissions::from_mode(0o755)).unwrap();
        assert!(
            temporary.is_dir(),
            "the directory that could not be removed"
        );
        assert_eq!(contents(&out), before);
    }

    /// The directory published over a previous output directory takes that
    /// directory's permission bits; one published to a new path, or over a
    /// part file, keeps the bits it was created with.
    #[cfg(unix)]
    #[test]
    fn a_replacement_keeps_the_previous_directorys_permission_bits() {
        use std::os::unix::fs::PermissionsExt;
        let mode = |path: &Path| fs::metadata(path).unwrap().permissions().mode() & 0o7777;
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        write_output(&out, "old");
        let created = mode(&out);
        // A umask of 077 creates a directory with 0o700, so the bits to keep
        // differ from the created ones either way.
        let kept = if created == 0o700 { 0o750 } else { 0o700 };
        fs::set_permissions(&out, fs::Permissions::from_mode(kept)).unwrap();
        unpublished(&out, "new").publish(true).unwrap();
        assert_eq!(mode(&out), kept);
        assert_eq!(Manifest::read(&out).unwrap().run_id, "new");
        assert_eq!(names(dir.path()), ["out.parquet"]);

        let file = dir.path().join("part.parquet");
        fs::copy(out.join("part-00000.parquet"), &file).unwrap();
        fs::set_permissions(&file, fs::Permissions::from_mode(0o600)).unwrap();
        unpublished(&file, "over-a-file").publish(true).unwrap();
        assert_eq!(mode(&file), created);
        let new_path = dir.path().join("new.parquet");
        unpublished(&new_path, "new-path").publish(false).unwrap();
        assert_eq!(mode(&new_path), created);
    }

    /// The destination is checked again as the run publishes: an entry vep-rs
    /// did not write that arrived during the run, or a destination the run
    /// was not told to replace, refuses the move with everything in place.
    #[test]
    fn publishing_checks_the_destination_again_and_leaves_it_as_it_is() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        write_output(&out, "old");
        let before = contents(&out);

        let mut refused = unpublished(&out, "unasked");
        assert!(matches!(
            refused.publish(false).unwrap_err(),
            PublishError::Exists
        ));
        refused.discard().unwrap();
        assert_eq!(contents(&out), before);
        assert_eq!(names(dir.path()), ["out.parquet"]);

        let mut refused = unpublished(&out, "new");
        refused.check_destination().unwrap();
        fs::write(out.join("notes.txt"), b"mine\n").unwrap();
        let error = refused.publish(true).unwrap_err();
        refused.discard().unwrap();
        assert!(
            matches!(&error, PublishError::Previous(PreviousOutput::Foreign(p)) if p == &out.join("notes.txt")),
            "{error}"
        );
        assert_eq!(
            error.to_string(),
            "output directory contains files vep-rs did not write"
        );
        assert_eq!(names(dir.path()), ["out.parquet"]);
        let mut expected = before.clone();
        expected.push(("notes.txt".to_string(), b"mine\n".to_vec()));
        expected.sort();
        assert_eq!(contents(&out), expected);
    }

    /// The temporary directory gone, the publication fails and the previous
    /// output is whole under its own name: over a directory the permission
    /// step fails before any rename, and over a part file the second rename
    /// fails once the part is aside, which is moved back.
    #[test]
    fn a_failed_rename_into_place_restores_the_previous_output() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        write_output(&out, "old");
        let before = contents(&out);
        let file = dir.path().join("part.parquet");
        fs::copy(out.join("part-00000.parquet"), &file).unwrap();
        let part = fs::read(&file).unwrap();
        for destination in [&out, &file] {
            let mut unpublished = unpublished(destination, "new");
            fs::remove_dir_all(unpublished.path()).unwrap();
            let error = unpublished.publish(true).unwrap_err();
            assert!(matches!(error, PublishError::Io(_)), "{error}");
        }
        assert_eq!(contents(&out), before);
        assert_eq!(Manifest::read(&out).unwrap().run_id, "old");
        assert_eq!(fs::read(&file).unwrap(), part);
        assert_eq!(names(dir.path()), ["out.parquet", "part.parquet"]);
    }

    /// A symbolic link at the temporary name is not followed when the
    /// replacement takes the previous directory's bits: the publication
    /// fails, what the link points at keeps its bits, the destination is as
    /// it was, and the discard removes the link alone.
    #[cfg(unix)]
    #[test]
    fn a_link_at_the_temporary_name_keeps_its_targets_permission_bits() {
        use std::os::unix::fs::{symlink, PermissionsExt};
        let mode = |path: &Path| fs::metadata(path).unwrap().permissions().mode() & 0o7777;
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        write_output(&out, "old");
        fs::set_permissions(&out, fs::Permissions::from_mode(0o777)).unwrap();
        let before = contents(&out);
        let elsewhere = dir.path().join("elsewhere");
        fs::create_dir(&elsewhere).unwrap();
        fs::set_permissions(&elsewhere, fs::Permissions::from_mode(0o700)).unwrap();

        let mut replacement = unpublished(&out, "new");
        fs::rename(replacement.path(), dir.path().join("aside")).unwrap();
        symlink(&elsewhere, replacement.path()).unwrap();
        let error = replacement.publish(true).unwrap_err();
        assert!(matches!(error, PublishError::Io(_)), "{error}");
        assert_eq!(mode(&elsewhere), 0o700);
        assert_eq!(mode(&out), 0o777);
        assert_eq!(contents(&out), before);
        replacement.discard().unwrap();
        assert_eq!(names(dir.path()), ["aside", "elsewhere", "out.parquet"]);
    }

    /// A previous output whose entries cannot be removed once it is aside
    /// stays whole beside the new output, under the name the error states.
    #[cfg(unix)]
    #[test]
    fn a_previous_output_that_cannot_be_removed_is_left_beside_the_new_one() {
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        write_output(&out, "old");
        let before = contents(&out);
        let read_only = fs::Permissions::from_mode(0o555);
        fs::set_permissions(&out, read_only).unwrap();
        if fs::write(out.join("probe"), b"").is_ok() {
            // This user's writes ignore the mode, so the removal cannot fail.
            fs::remove_file(out.join("probe")).unwrap();
            return;
        }
        let mut unpublished = unpublished(&out, "new");
        let aside = unpublished.previous.clone();
        let error = unpublished.publish(true).unwrap_err();
        assert!(
            matches!(&error, PublishError::Unremoved { previous, .. } if previous == &aside),
            "{error}"
        );
        let stated = format!(
            "the output is complete, but the previous output is left at {}: could not remove the previous output: ",
            aside.display()
        );
        assert!(error.to_string().starts_with(&stated), "{error}");
        assert_eq!(Manifest::read(&out).unwrap().run_id, "new");
        assert_eq!(contents(&aside), before);
        assert_eq!(
            names(dir.path()),
            [aside.file_name().unwrap().to_str().unwrap(), "out.parquet"]
        );
        assert!(out.join(MANIFEST_NAME).exists());
        let mode = |path: &Path| fs::metadata(path).unwrap().permissions().mode() & 0o7777;
        assert_eq!(mode(&out), 0o555, "the new directory takes the bits");
        for directory in [&aside, &out] {
            fs::set_permissions(directory, fs::Permissions::from_mode(0o755)).unwrap();
        }
    }
}
