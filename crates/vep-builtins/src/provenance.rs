// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Content identity and completion metadata for binary annotation stores.
//!
//! A `.vpd.provenance.json` binds the complete compressed source, conversion
//! coordinates and both binary files. Publication stages files on the output
//! filesystem and replaces the manifest last. Inputs and published stores must
//! remain immutable while annotation runs use them.

use std::fs::{self, File, Metadata};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use serde::{Deserialize, Serialize};
use vep_core::cache_manifest::{file_digest, Digest, HashingReader};

use crate::tabix::TabixConfig;
use crate::PluginError;

const MANIFEST_VERSION: u32 = 1;
const MAX_MANIFEST_BYTES: u64 = 16 * 1024;

/// Coordinate interpretation used by the plugin-data converter.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConversionSettings {
    pub chr_col: usize,
    pub start_col: usize,
    pub end_col: Option<usize>,
    pub zero_based: bool,
}

impl From<&TabixConfig> for ConversionSettings {
    fn from(config: &TabixConfig) -> Self {
        Self {
            chr_col: config.chr_col,
            start_col: config.start_col,
            end_col: config.end_col,
            zero_based: config.zero_based,
        }
    }
}

/// Lowercase hexadecimal SHA256, suitable for logs and the sidecar.
pub fn digest_hex(digest: Digest) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

/// Path of the completion manifest for `data`.
pub fn manifest_path(data: &Path) -> PathBuf {
    let mut name = data.as_os_str().to_os_string();
    name.push(".provenance.json");
    name.into()
}

fn failure(message: impl Into<String>) -> PluginError {
    PluginError::Init(message.into())
}

fn at(path: &Path, error: impl std::fmt::Display) -> PluginError {
    failure(format!("{}: {error}", path.display()))
}

fn hash_file(path: &Path) -> Result<String, PluginError> {
    file_digest(path)
        .map(digest_hex)
        .map_err(|error| at(path, error))
}

fn hash_bytes(bytes: &[u8]) -> Result<String, PluginError> {
    HashingReader::new(bytes)
        .digest_to_end()
        .map(digest_hex)
        .map_err(PluginError::Io)
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Provenance {
    version: u32,
    converter: String,
    settings: ConversionSettings,
    pub(crate) source_sha256: String,
    pub(crate) data_sha256: String,
    pub(crate) index_sha256: String,
}

impl Provenance {
    pub(crate) fn load(data: &Path, config: &TabixConfig) -> Result<Self, PluginError> {
        let path = manifest_path(data);
        let mut bytes = Vec::new();
        File::open(&path)
            .and_then(|file| file.take(MAX_MANIFEST_BYTES + 1).read_to_end(&mut bytes))
            .map_err(|error| at(&path, error))?;
        if bytes.len() as u64 > MAX_MANIFEST_BYTES {
            return Err(at(&path, "provenance manifest exceeds the size limit"));
        }
        let manifest: Self = serde_json::from_slice(&bytes).map_err(|error| at(&path, error))?;
        if manifest.version != MANIFEST_VERSION {
            return Err(at(&path, "unsupported provenance version"));
        }
        for digest in [
            &manifest.source_sha256,
            &manifest.data_sha256,
            &manifest.index_sha256,
        ] {
            if digest.len() != 64
                || !digest
                    .bytes()
                    .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
            {
                return Err(at(&path, "invalid SHA256 identity"));
            }
        }
        if manifest.settings != ConversionSettings::from(config) {
            return Err(at(&path, "conversion settings mismatch"));
        }
        // The converter has no reference-allele fallback for a nonnumeric end.
        if config.ref_col.is_some() && config.end_col.is_some_and(|end| end != config.start_col) {
            return Err(at(
                &path,
                "reference-dependent end settings are not verifiable",
            ));
        }
        Ok(manifest)
    }

    /// Checks the exact bytes the binary loader will map and parse.
    pub(crate) fn verify(
        &self,
        source: &Path,
        data: &[u8],
        index: &[u8],
    ) -> Result<(), PluginError> {
        if hash_bytes(data)? != self.data_sha256 {
            return Err(failure("binary data SHA256 mismatch"));
        }
        if hash_bytes(index)? != self.index_sha256 {
            return Err(failure("binary index SHA256 mismatch"));
        }
        let before = FileState::read(source)?;
        if hash_file(source)? != self.source_sha256 {
            return Err(at(source, "source SHA256 mismatch"));
        }
        before.check(source)?;
        Ok(())
    }
}

#[derive(Debug, PartialEq, Eq)]
struct FileState {
    len: u64,
    modified: Option<SystemTime>,
    created: Option<SystemTime>,
    #[cfg(unix)]
    inode: (u64, u64, i64, i64),
}

impl FileState {
    fn of(metadata: &Metadata) -> Result<Self, PluginError> {
        if !metadata.is_file() {
            return Err(failure("provenance requires a regular source file"));
        }
        Ok(Self {
            len: metadata.len(),
            modified: metadata.modified().ok(),
            created: metadata.created().ok(),
            #[cfg(unix)]
            inode: {
                use std::os::unix::fs::MetadataExt;
                (
                    metadata.dev(),
                    metadata.ino(),
                    metadata.ctime(),
                    metadata.ctime_nsec(),
                )
            },
        })
    }

    fn read(path: &Path) -> Result<Self, PluginError> {
        Self::of(&fs::metadata(path).map_err(|error| at(path, error))?)
    }

    fn check(&self, path: &Path) -> Result<(), PluginError> {
        if *self != Self::read(path)? {
            return Err(at(path, "source changed during conversion or verification"));
        }
        Ok(())
    }
}

/// Exclusive, staged publication of a binary pair and its completion manifest.
///
/// Dropping this value before `publish` leaves any completed destination
/// untouched. A process killed while holding the lock leaves its
/// `.vpd.conversion-lock` directory for the operator to remove after confirming
/// no converter is using that destination.
pub struct Publication {
    source_path: PathBuf,
    source_file: File,
    source_state: FileState,
    output: PathBuf,
    settings: ConversionSettings,
    staging: tempfile::TempDir,
    lock: PublicationLock,
}

struct PublicationLock(PathBuf);

impl Drop for PublicationLock {
    fn drop(&mut self) {
        let _ = fs::remove_dir(&self.0);
    }
}

impl Publication {
    /// Reserve the destination and capture the source handle and settings.
    pub fn begin(
        source: &Path,
        output: &Path,
        settings: ConversionSettings,
    ) -> Result<Self, PluginError> {
        let source_path = fs::canonicalize(source).map_err(|error| at(source, error))?;
        let index = output.with_extension("vpdi");
        if index == output {
            return Err(at(output, "data and index must have distinct paths"));
        }
        for destination in [output.to_path_buf(), index, manifest_path(output)] {
            if fs::canonicalize(&destination).is_ok_and(|path| path == source_path) {
                return Err(at(&destination, "publication would replace its own source"));
            }
        }
        let parent = output
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let mut lock_name = output.as_os_str().to_os_string();
        lock_name.push(".conversion-lock");
        let lock_path = PathBuf::from(lock_name);
        fs::create_dir(&lock_path).map_err(|error| at(&lock_path, error))?;
        let lock = PublicationLock(lock_path);
        let staging = tempfile::Builder::new()
            .prefix(".vpd-conversion-")
            .tempdir_in(parent)
            .map_err(|error| at(parent, error))?;
        let source_file = File::open(source).map_err(|error| at(source, error))?;
        let source_state =
            FileState::of(&source_file.metadata().map_err(|error| at(source, error))?)?;
        source_state.check(source)?;
        Ok(Self {
            source_path: source.to_path_buf(),
            source_file,
            source_state,
            output: output.to_path_buf(),
            settings,
            staging,
            lock,
        })
    }

    /// Hash the complete compressed stream used by the parser, including its tail.
    ///
    /// Call once, consume through the parser, then call `digest_to_end` and pass
    /// its result to `publish`.
    pub fn source_reader(&self) -> Result<HashingReader<File>, PluginError> {
        self.source_file
            .try_clone()
            .map(HashingReader::new)
            .map_err(|error| at(&self.source_path, error))
    }

    /// Path for the writer and concordance verifier; never the live destination.
    pub fn data_path(&self) -> PathBuf {
        self.staging.path().join("data.vpd")
    }

    /// Publish only after the writer and requested concordance checks succeed.
    ///
    /// The digest must cover the bytes consumed by the parser. The source is
    /// rehashed in full here; metadata checks supplement, never replace, SHA256.
    pub fn publish(
        self,
        consumed_source: Digest,
        settings: ConversionSettings,
    ) -> Result<(), PluginError> {
        if settings != self.settings {
            return Err(failure("conversion settings changed before publication"));
        }
        self.source_state.check(&self.source_path)?;
        let source_sha256 = hash_file(&self.source_path)?;
        if source_sha256 != digest_hex(consumed_source) {
            return Err(at(
                &self.source_path,
                "source SHA256 changed during conversion",
            ));
        }
        self.source_state.check(&self.source_path)?;
        let data = self.data_path();
        let index = data.with_extension("vpdi");
        let manifest = Provenance {
            version: MANIFEST_VERSION,
            converter: format!("vep-cache-converter {}", env!("CARGO_PKG_VERSION")),
            settings,
            source_sha256,
            data_sha256: hash_file(&data)?,
            index_sha256: hash_file(&index)?,
        };
        let staged_manifest = manifest_path(&data);
        let mut file =
            File::create(&staged_manifest).map_err(|error| at(&staged_manifest, error))?;
        serde_json::to_writer_pretty(&mut file, &manifest)
            .map_err(|error| at(&staged_manifest, error))?;
        file.write_all(b"\n")
            .map_err(|error| at(&staged_manifest, error))?;
        for path in [&data, &index, &staged_manifest] {
            File::open(path)
                .and_then(|file| file.sync_all())
                .map_err(|error| at(path, error))?;
        }
        self.source_state.check(&self.source_path)?;
        let destination_manifest = manifest_path(&self.output);
        match fs::remove_file(&destination_manifest) {
            Ok(()) => {}
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(at(&destination_manifest, error)),
        }
        fs::rename(&data, &self.output).map_err(|error| at(&self.output, error))?;
        let destination_index = self.output.with_extension("vpdi");
        fs::rename(&index, &destination_index).map_err(|error| at(&destination_index, error))?;
        fs::rename(&staged_manifest, &destination_manifest)
            .map_err(|error| at(&destination_manifest, error))?;
        tracing::info!(
            backend = "binary",
            source = %self.source_path.display(),
            data = %self.output.display(),
            source_sha256 = %manifest.source_sha256,
            data_sha256 = %manifest.data_sha256,
            index_sha256 = %manifest.index_sha256,
            "published binary annotation provenance"
        );
        drop(self.lock);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_store::BinaryStoreWriter;
    use crate::tabix::TabixRecord;

    const SETTINGS: ConversionSettings = ConversionSettings {
        chr_col: 0,
        start_col: 1,
        end_col: None,
        zero_based: false,
    };

    fn stage(publication: &Publication) {
        BinaryStoreWriter::new(publication.data_path())
            .write(
                &[],
                &[(
                    "21".into(),
                    vec![TabixRecord {
                        chr: "21".into(),
                        start: 200,
                        end: 200,
                        columns: vec!["21".into(), "200".into()],
                        header: None,
                    }],
                )],
            )
            .unwrap();
    }

    fn setup() -> (tempfile::TempDir, PathBuf, PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join("data.tsv.gz");
        let output = dir.path().join("data.tsv.vpd");
        fs::write(&source, b"original").unwrap();
        (dir, source, output)
    }

    #[test]
    fn source_or_settings_changes_cannot_publish() {
        for change in ["source", "settings", "consumed-stream"] {
            let (_dir, source, output) = setup();
            let publication = Publication::begin(&source, &output, SETTINGS).unwrap();
            let mut digest = publication
                .source_reader()
                .unwrap()
                .digest_to_end()
                .unwrap();
            stage(&publication);
            let mut settings = SETTINGS;
            match change {
                "source" => fs::write(&source, b"replaced").unwrap(),
                "settings" => settings.zero_based = true,
                _ => digest[0] ^= 1,
            }
            assert!(publication.publish(digest, settings).is_err(), "{change}");
            assert!(!output.exists());
            assert!(!manifest_path(&output).exists());
        }
    }

    #[test]
    fn source_restored_after_reading_is_still_a_detected_change() {
        let (_dir, source, output) = setup();
        let publication = Publication::begin(&source, &output, SETTINGS).unwrap();
        let digest = publication
            .source_reader()
            .unwrap()
            .digest_to_end()
            .unwrap();
        stage(&publication);
        let modified = fs::metadata(&source).unwrap().modified().unwrap();
        fs::write(&source, b"replaced").unwrap();
        fs::write(&source, b"original").unwrap();
        fs::File::options()
            .write(true)
            .open(&source)
            .unwrap()
            .set_times(
                fs::FileTimes::new().set_modified(modified + std::time::Duration::from_secs(1)),
            )
            .unwrap();
        assert!(publication.publish(digest, SETTINGS).is_err());
        assert!(!manifest_path(&output).exists());
    }

    #[test]
    fn incomplete_conversion_preserves_the_completed_pair() {
        let (_dir, source, output) = setup();
        let first = Publication::begin(&source, &output, SETTINGS).unwrap();
        let digest = first.source_reader().unwrap().digest_to_end().unwrap();
        stage(&first);
        first.publish(digest, SETTINGS).unwrap();
        let original = [
            fs::read(&output).unwrap(),
            fs::read(output.with_extension("vpdi")).unwrap(),
            fs::read(manifest_path(&output)).unwrap(),
        ];

        let incomplete = Publication::begin(&source, &output, SETTINGS).unwrap();
        fs::write(incomplete.data_path(), b"partial").unwrap();
        assert!(incomplete.publish(digest, SETTINGS).is_err());
        assert_eq!(fs::read(&output).unwrap(), original[0]);
        assert_eq!(
            fs::read(output.with_extension("vpdi")).unwrap(),
            original[1]
        );
        assert_eq!(fs::read(manifest_path(&output)).unwrap(), original[2]);
    }

    #[test]
    fn failed_index_publication_removes_the_old_certificate() {
        let (_dir, source, output) = setup();
        let first = Publication::begin(&source, &output, SETTINGS).unwrap();
        let digest = first.source_reader().unwrap().digest_to_end().unwrap();
        stage(&first);
        first.publish(digest, SETTINGS).unwrap();
        let second = Publication::begin(&source, &output, SETTINGS).unwrap();
        stage(&second);
        fs::remove_file(output.with_extension("vpdi")).unwrap();
        fs::create_dir(output.with_extension("vpdi")).unwrap();
        assert!(second.publish(digest, SETTINGS).is_err());
        assert!(!manifest_path(&output).exists());
    }

    #[test]
    fn concurrent_publication_is_refused_and_abandoned_staging_is_cleaned() {
        let (dir, source, output) = setup();
        let first = Publication::begin(&source, &output, SETTINGS).unwrap();
        let staging = first.data_path().parent().unwrap().to_path_buf();
        assert!(Publication::begin(&source, &output, SETTINGS).is_err());
        drop(first);
        assert!(!staging.exists());
        let next = Publication::begin(&source, &output, SETTINGS).unwrap();
        drop(next);
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    #[test]
    fn publication_paths_do_not_alias_source_data_or_index() {
        let (_dir, source, output) = setup();
        assert!(Publication::begin(&source, &source, SETTINGS).is_err());
        assert!(Publication::begin(&source, &output.with_extension("vpdi"), SETTINGS).is_err());
        let index_source = output.with_extension("vpdi");
        fs::write(&index_source, b"source").unwrap();
        assert!(Publication::begin(&index_source, &output, SETTINGS).is_err());
        assert_eq!(fs::read(source).unwrap(), b"original");
        assert_eq!(fs::read(index_source).unwrap(), b"source");
    }
}
