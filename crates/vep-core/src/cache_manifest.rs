// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cache checksum manifest, `CHECKSUMS.sha256` in the cache directory
//! beside `info.json`: one line per file the cache holds, the manifest itself
//! excepted, in the text shape `sha256sum` writes and checks
//! (`<64 hex digits>  <path relative to the cache directory>`), the paths in
//! byte order with `/` between components. The cache tools write it last;
//! `sha256sum --check CHECKSUMS.sha256` run in the cache directory verifies
//! every file, and the annotation run verifies each file it reads whole as it
//! reads it. A file whose digest differs from its entry and a file the manifest
//! does not list are the same refusal, [`ManifestError::Mismatch`]. A cache
//! without the file is a manifest that is not present: every check passes.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Component, Path, PathBuf};

use sha2::{Digest as _, Sha256};

/// The manifest's file name, in the cache directory.
pub const MANIFEST_FILE: &str = "CHECKSUMS.sha256";

/// A 32-byte SHA256 digest, the value a manifest line records.
pub type Digest = [u8; 32];

/// Directory levels the writer descends below the cache directory before it
/// treats the tree as a symbolic-link loop; a cache is three levels deep.
const MAX_DEPTH: usize = 16;

/// Bytes read per step while a file is hashed.
const READ_CHUNK: usize = 1 << 20;

#[derive(Debug, thiserror::Error)]
pub enum ManifestError {
    /// The bytes of `relative` differ from its entry, or the manifest has no
    /// entry for it.
    #[error("checksum mismatch for {relative} in {dir}")]
    Mismatch { relative: String, dir: String },

    /// A line of the manifest is not `<64 hex digits>  <path>`, or repeats a
    /// path an earlier line listed.
    #[error("{path}: line {line} {reason}")]
    Malformed {
        path: String,
        line: usize,
        reason: String,
    },

    /// A path to verify does not lie under the cache directory.
    #[error("{path} is not under the cache directory {dir}")]
    Outside { path: String, dir: String },

    /// A file the manifest lists is not in the cache directory.
    #[error("{relative} is listed in {MANIFEST_FILE} but missing from {dir}")]
    Missing { relative: String, dir: String },

    /// The manifest or a file to hash could not be read or written.
    #[error("{path}: {source}")]
    Io { path: String, source: io::Error },

    /// The writer descended more than [`MAX_DEPTH`] directory levels.
    #[error("{path}: directories nest more than {MAX_DEPTH} levels deep")]
    TooDeep { path: String },
}

fn io_error(path: &Path) -> impl FnOnce(io::Error) -> ManifestError {
    let path = path.display().to_string();
    move |source| ManifestError::Io { path, source }
}

/// The manifest of one cache directory, loaded once per run.
#[derive(Debug, Clone)]
pub struct CacheManifest {
    dir: PathBuf,
    /// Digest by relative path; `None` when the directory has no manifest.
    entries: Option<BTreeMap<String, Digest>>,
}

impl CacheManifest {
    /// Reads `<dir>/CHECKSUMS.sha256`; a directory without the file gives a
    /// manifest that [`Self::is_present`] denies and whose checks all pass.
    pub fn load(dir: impl AsRef<Path>) -> Result<Self, ManifestError> {
        let dir = dir.as_ref().to_path_buf();
        let path = dir.join(MANIFEST_FILE);
        let text = match fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(Self::absent(dir)),
            Err(error) => return Err(io_error(&path)(error)),
        };
        let entries = parse(&text, &path)?;
        Ok(Self {
            dir,
            entries: Some(entries),
        })
    }

    /// A manifest that is not present: every check passes.
    pub fn absent(dir: impl AsRef<Path>) -> Self {
        Self {
            dir: dir.as_ref().to_path_buf(),
            entries: None,
        }
    }

    /// Whether the cache directory holds a manifest.
    pub fn is_present(&self) -> bool {
        self.entries.is_some()
    }

    /// The number of files the manifest lists; 0 when it is not present.
    pub fn len(&self) -> usize {
        self.entries.as_ref().map_or(0, BTreeMap::len)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Checks that every file the manifest lists is in the cache directory,
    /// naming the first one, in path order, that is not: a copy that lacks a
    /// listed shard or variation file would otherwise annotate as if the cache
    /// never held it. A directory in a listed file's place is missing; a stat
    /// that fails for a reason other than absence is [`ManifestError::Io`].
    pub fn require_complete(&self) -> Result<(), ManifestError> {
        let Some(entries) = &self.entries else {
            return Ok(());
        };
        for relative in entries.keys() {
            let path = self.dir.join(relative);
            let missing = || ManifestError::Missing {
                relative: relative.clone(),
                dir: self.dir.display().to_string(),
            };
            match fs::metadata(&path) {
                Ok(metadata) if metadata.is_file() => {}
                Ok(_) => return Err(missing()),
                Err(error) if error.kind() == io::ErrorKind::NotFound => return Err(missing()),
                Err(error) => return Err(io_error(&path)(error)),
            }
        }
        Ok(())
    }

    /// `path` relative to the cache directory, `/` between components, the
    /// form the manifest lists. `path` must lie under the directory lexically:
    /// one that does not start with it, or whose remainder holds a `..`, a
    /// root or a prefix component and so names a file elsewhere, is
    /// [`ManifestError::Outside`].
    pub fn relative(&self, path: &Path) -> Result<String, ManifestError> {
        let outside = || ManifestError::Outside {
            path: path.display().to_string(),
            dir: self.dir.display().to_string(),
        };
        let rest = path.strip_prefix(&self.dir).map_err(|_| outside())?;
        if !stays_under(rest) {
            return Err(outside());
        }
        Ok(relative_text(rest))
    }

    /// Checks `bytes`, the whole content of the file at `path`, against the
    /// manifest.
    pub fn verify_bytes(&self, path: &Path, bytes: &[u8]) -> Result<(), ManifestError> {
        if !self.is_present() {
            return Ok(());
        }
        self.verify_digest(path, &sha256(bytes))
    }

    /// Checks the digest of the whole content of the file at `path` against
    /// the manifest.
    pub fn verify_digest(&self, path: &Path, digest: &Digest) -> Result<(), ManifestError> {
        let Some(entries) = &self.entries else {
            return Ok(());
        };
        let relative = self.relative(path)?;
        match entries.get(&relative) {
            Some(expected) if expected == digest => Ok(()),
            _ => Err(self.mismatch(relative)),
        }
    }

    /// Checks that the manifest lists the file at `path`, for a file whose
    /// content is read by byte range and so never hashed whole by a run.
    pub fn require_listed(&self, path: &Path) -> Result<(), ManifestError> {
        let Some(entries) = &self.entries else {
            return Ok(());
        };
        let relative = self.relative(path)?;
        if entries.contains_key(&relative) {
            Ok(())
        } else {
            Err(self.mismatch(relative))
        }
    }

    /// The refusal for the file at `relative`, listed with another digest or
    /// not at all.
    pub fn mismatch(&self, relative: String) -> ManifestError {
        ManifestError::Mismatch {
            relative,
            dir: self.dir.display().to_string(),
        }
    }

    /// Hashes every regular file under `dir`, the manifest itself excepted,
    /// and writes `<dir>/CHECKSUMS.sha256`, replacing any earlier one; returns
    /// the number of files listed. Symbolic links are followed.
    pub fn write(dir: impl AsRef<Path>) -> Result<usize, ManifestError> {
        let dir = dir.as_ref();
        let mut files = Vec::new();
        collect_files(dir, 0, &mut files)?;
        let mut entries = BTreeMap::new();
        for path in files {
            let digest = file_digest(&path).map_err(io_error(&path))?;
            let relative = path
                .strip_prefix(dir)
                .map(relative_text)
                .expect("collected under dir");
            entries.insert(relative, digest);
        }
        let manifest = dir.join(MANIFEST_FILE);
        fs::write(&manifest, render(&entries)).map_err(io_error(&manifest))?;
        Ok(entries.len())
    }
}

/// Whether `rest`, a path relative to the cache directory, names a file under
/// it lexically: every component normal or `.`, so no `..`, root or prefix
/// component leads elsewhere.
fn stays_under(rest: &Path) -> bool {
    rest.components()
        .all(|component| matches!(component, Component::Normal(_) | Component::CurDir))
}

/// `rest`, a path below the cache directory of normal and `.` components, as
/// the manifest spells it.
fn relative_text(rest: &Path) -> String {
    let mut text = String::new();
    for component in rest.components() {
        if let Component::Normal(name) = component {
            if !text.is_empty() {
                text.push('/');
            }
            text.push_str(&name.to_string_lossy());
        }
    }
    text
}

/// The regular files under `dir`, through symbolic links, the manifest at
/// depth 0 left out.
fn collect_files(dir: &Path, depth: usize, files: &mut Vec<PathBuf>) -> Result<(), ManifestError> {
    if depth > MAX_DEPTH {
        return Err(ManifestError::TooDeep {
            path: dir.display().to_string(),
        });
    }
    let mut entries: Vec<PathBuf> = fs::read_dir(dir)
        .map_err(io_error(dir))?
        .map(|entry| entry.map(|e| e.path()))
        .collect::<io::Result<_>>()
        .map_err(io_error(dir))?;
    entries.sort();
    for path in entries {
        let metadata = fs::metadata(&path).map_err(io_error(&path))?;
        if metadata.is_dir() {
            collect_files(&path, depth + 1, files)?;
        } else if metadata.is_file() {
            if depth == 0 && path.file_name().is_some_and(|name| name == MANIFEST_FILE) {
                continue;
            }
            files.push(path);
        }
    }
    Ok(())
}

/// The SHA256 digest of the file at `path`, read in chunks.
pub fn file_digest(path: &Path) -> io::Result<Digest> {
    HashingReader::new(File::open(path)?).digest_to_end()
}

/// The SHA256 digest of `bytes`.
fn sha256(bytes: &[u8]) -> Digest {
    Sha256::digest(bytes).into()
}

/// The lower-case hexadecimal spelling of a digest, as `sha256sum` prints it.
fn to_hex(digest: &Digest) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

/// The digest a 64-character hexadecimal spelling denotes, in either case;
/// `None` for any other text.
fn from_hex(text: &str) -> Option<Digest> {
    let bytes = text.as_bytes();
    if bytes.len() != 64 {
        return None;
    }
    let mut digest = [0u8; 32];
    for (byte, pair) in digest.iter_mut().zip(bytes.chunks_exact(2)) {
        let high = (pair[0] as char).to_digit(16)?;
        let low = (pair[1] as char).to_digit(16)?;
        *byte = (high as u8) << 4 | low as u8;
    }
    Some(digest)
}

/// The manifest text for `entries`: one `sha256sum` line per path, in the
/// map's byte order.
pub fn render(entries: &BTreeMap<String, Digest>) -> String {
    let mut text = String::with_capacity(entries.len() * 100);
    for (relative, digest) in entries {
        text.push_str(&to_hex(digest));
        text.push_str("  ");
        text.push_str(relative);
        text.push('\n');
    }
    text
}

/// The entries of a manifest text. `sha256sum` separates the digest from the
/// path by two spaces in text mode and by a space and `*` in binary mode; both
/// are read. A path's leading `./` is dropped, so a line written as
/// `sha256sum ./info.json` names the file the run looks up; a path that leads
/// outside the directory through a `..`, root or prefix component is malformed,
/// as no file the run looks up can match it.
pub fn parse(text: &str, path: &Path) -> Result<BTreeMap<String, Digest>, ManifestError> {
    let malformed = |line: usize, reason: String| ManifestError::Malformed {
        path: path.display().to_string(),
        line,
        reason,
    };
    let mut entries = BTreeMap::new();
    for (index, line) in text.lines().enumerate() {
        let number = index + 1;
        if line.is_empty() {
            continue;
        }
        let (digest, mut relative) = parse_line(line)
            .ok_or_else(|| malformed(number, "is not a sha256sum line".to_string()))?;
        while let Some(stripped) = relative.strip_prefix("./") {
            relative = stripped;
        }
        if relative.is_empty() {
            return Err(malformed(number, "names no file".to_string()));
        }
        if !stays_under(Path::new(relative)) {
            return Err(malformed(
                number,
                format!("names {relative}, outside the cache directory"),
            ));
        }
        if entries.insert(relative.to_string(), digest).is_some() {
            return Err(malformed(number, format!("lists {relative} twice")));
        }
    }
    Ok(entries)
}

/// The digest and path of one `sha256sum` line; `None` for a line of another
/// shape.
fn parse_line(line: &str) -> Option<(Digest, &str)> {
    let (hex, rest) = line.split_once(' ')?;
    let digest = from_hex(hex)?;
    let relative = rest.strip_prefix(' ').or_else(|| rest.strip_prefix('*'))?;
    Some((digest, relative))
}

/// A reader that hashes every byte it hands out, for a cache file streamed
/// whole; [`Self::digest_to_end`] reads whatever its consumer left and gives
/// the digest of the whole stream.
pub struct HashingReader<R> {
    inner: R,
    hasher: Sha256,
}

impl<R: Read> HashingReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            hasher: Sha256::new(),
        }
    }

    /// Reads the rest of the stream and returns the digest of all of it.
    pub fn digest_to_end(mut self) -> io::Result<Digest> {
        let mut buffer = vec![0u8; READ_CHUNK];
        loop {
            let read = self.inner.read(&mut buffer)?;
            if read == 0 {
                return Ok(self.hasher.finalize().into());
            }
            self.hasher.update(&buffer[..read]);
        }
    }
}

impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(buf)?;
        self.hasher.update(&buf[..read]);
        Ok(read)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// A cache directory of three files in two levels.
    fn cache() -> tempfile::TempDir {
        let dir = tempfile::tempdir().unwrap();
        fs::create_dir_all(dir.path().join("transcripts/21")).unwrap();
        fs::write(dir.path().join("info.json"), b"{}\n").unwrap();
        fs::write(dir.path().join("transcripts/21/1-1000000.json"), b"[]").unwrap();
        fs::write(dir.path().join("chr_synonyms.txt"), b"21\tchr21\n").unwrap();
        dir
    }

    /// The writer lists every file but itself, by path, in `sha256sum`'s text
    /// shape, and writing again replaces the file rather than listing it.
    #[test]
    fn write_lists_every_file_but_the_manifest() {
        let dir = cache();
        assert_eq!(CacheManifest::write(dir.path()).unwrap(), 3);
        assert_eq!(CacheManifest::write(dir.path()).unwrap(), 3);
        let text = fs::read_to_string(dir.path().join(MANIFEST_FILE)).unwrap();
        let expected = format!(
            "{}  chr_synonyms.txt\n{}  info.json\n{}  transcripts/21/1-1000000.json\n",
            to_hex(&sha256(b"21\tchr21\n")),
            to_hex(&sha256(b"{}\n")),
            to_hex(&sha256(b"[]")),
        );
        assert_eq!(text, expected);
    }

    /// The digests FIPS 180-4's example vectors and `sha256sum` give, and
    /// their spelling read back in either case; a 64-character text that is
    /// not hexadecimal and a shorter one denote no digest.
    #[test]
    fn known_vectors_and_hex_round_trip() {
        let empty = sha256(b"");
        assert_eq!(
            to_hex(&empty),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        let abc = sha256(b"abc");
        assert_eq!(
            to_hex(&abc),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(from_hex(&to_hex(&abc)), Some(abc));
        assert_eq!(from_hex(&to_hex(&abc).to_uppercase()), Some(abc));
        assert_eq!(from_hex("abc"), None);
        assert_eq!(from_hex(&"g".repeat(64)), None);
    }

    /// Bytes equal to the file pass, a changed byte and an unlisted file fail
    /// with the same line, and a manifest-less directory checks nothing.
    #[test]
    fn verify_bytes_against_entries() {
        let dir = cache();
        CacheManifest::write(dir.path()).unwrap();
        let manifest = CacheManifest::load(dir.path()).unwrap();
        assert!(manifest.is_present());
        assert_eq!(manifest.len(), 3);
        let shard = dir
            .path()
            .join("transcripts")
            .join("21")
            .join("1-1000000.json");
        manifest.verify_bytes(&shard, b"[]").unwrap();
        let error = manifest.verify_bytes(&shard, b"[ ]").unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "checksum mismatch for transcripts/21/1-1000000.json in {}",
                dir.path().display()
            )
        );
        let extra = dir
            .path()
            .join("transcripts")
            .join("21")
            .join("2000001-3000000.json");
        let error = manifest.verify_bytes(&extra, b"[]").unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "checksum mismatch for transcripts/21/2000001-3000000.json in {}",
                dir.path().display()
            )
        );
        manifest
            .require_listed(&dir.path().join("info.json"))
            .unwrap();
        assert!(manifest.require_listed(&extra).is_err());
        let outside = manifest.verify_bytes(Path::new("/elsewhere/info.json"), b"{}\n");
        assert!(matches!(outside, Err(ManifestError::Outside { .. })));

        fs::remove_file(dir.path().join(MANIFEST_FILE)).unwrap();
        let none = CacheManifest::load(dir.path()).unwrap();
        assert!(!none.is_present());
        none.verify_bytes(&shard, b"anything").unwrap();
        none.require_listed(&extra).unwrap();
    }

    /// A listed file that is not in the directory is named, the first in path
    /// order, a directory in its place counting as missing; a listed path whose
    /// stat fails otherwise (a file where a directory of its path should be) is
    /// the I/O error; a complete directory and a manifest-less one pass.
    #[test]
    fn require_complete_names_the_first_missing_file() {
        let dir = cache();
        CacheManifest::write(dir.path()).unwrap();
        let manifest = CacheManifest::load(dir.path()).unwrap();
        manifest.require_complete().unwrap();
        fs::remove_file(dir.path().join("transcripts/21/1-1000000.json")).unwrap();
        fs::remove_file(dir.path().join("info.json")).unwrap();
        assert_eq!(
            manifest.require_complete().unwrap_err().to_string(),
            format!(
                "info.json is listed in CHECKSUMS.sha256 but missing from {}",
                dir.path().display()
            )
        );
        fs::remove_file(dir.path().join("chr_synonyms.txt")).unwrap();
        fs::create_dir(dir.path().join("chr_synonyms.txt")).unwrap();
        assert_eq!(
            manifest.require_complete().unwrap_err().to_string(),
            format!(
                "chr_synonyms.txt is listed in CHECKSUMS.sha256 but missing from {}",
                dir.path().display()
            )
        );
        CacheManifest::absent(dir.path())
            .require_complete()
            .unwrap();

        let dir = cache();
        CacheManifest::write(dir.path()).unwrap();
        let manifest = CacheManifest::load(dir.path()).unwrap();
        fs::remove_dir_all(dir.path().join("transcripts")).unwrap();
        fs::write(dir.path().join("transcripts"), b"").unwrap();
        let error = manifest.require_complete().unwrap_err();
        assert!(
            matches!(&error, ManifestError::Io { path, .. } if path.ends_with("1-1000000.json")),
            "{error}"
        );
    }

    /// The directory as given keeps its spelling in the refusal and in the
    /// relative path: a `./` prefix and a trailing `/` on the directory change
    /// neither.
    #[test]
    fn relative_paths_follow_the_directory_as_given() {
        let dir = cache();
        CacheManifest::write(dir.path()).unwrap();
        let spelled = format!("{}/", dir.path().display());
        let manifest = CacheManifest::load(&spelled).unwrap();
        let shard = Path::new(&spelled).join("transcripts/21/1-1000000.json");
        assert_eq!(
            manifest.relative(&shard).unwrap(),
            "transcripts/21/1-1000000.json"
        );
        manifest.verify_bytes(&shard, b"[]").unwrap();
        assert_eq!(
            manifest.verify_bytes(&shard, b"x").unwrap_err().to_string(),
            format!("checksum mismatch for transcripts/21/1-1000000.json in {spelled}")
        );
    }

    /// A path under the directory whose remainder climbs through `..` names a
    /// file elsewhere, so it is outside the manifest rather than mapped onto
    /// the listed path its normal components spell; every check refuses it.
    #[test]
    fn relative_refuses_a_path_that_climbs_out_of_the_directory() {
        let dir = cache();
        CacheManifest::write(dir.path()).unwrap();
        let manifest = CacheManifest::load(dir.path()).unwrap();
        let elsewhere = dir.path().join("transcripts/../21/1-1000000.json");
        let error = manifest.relative(&elsewhere).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "{} is not under the cache directory {}",
                elsewhere.display(),
                dir.path().display()
            )
        );
        for climbing in [
            elsewhere,
            dir.path().join("../info.json"),
            dir.path().join("transcripts/21/../../info.json"),
        ] {
            assert!(
                matches!(
                    manifest.verify_bytes(&climbing, b"[]"),
                    Err(ManifestError::Outside { .. })
                ),
                "{}",
                climbing.display()
            );
            assert!(
                matches!(
                    manifest.require_listed(&climbing),
                    Err(ManifestError::Outside { .. })
                ),
                "{}",
                climbing.display()
            );
        }
        assert_eq!(
            manifest
                .relative(&dir.path().join("transcripts/./21/1-1000000.json"))
                .unwrap(),
            "transcripts/21/1-1000000.json"
        );
    }

    /// Both `sha256sum` separators are read, `./` is dropped from a path, and
    /// a line of another shape, a repeated path, an empty path or a path that
    /// leads outside the directory names the line.
    #[test]
    fn parse_reads_sha256sum_lines() {
        let hex = to_hex(&sha256(b"[]"));
        let path = Path::new("CHECKSUMS.sha256");
        let entries = parse(
            &format!("{hex}  a.json\n{hex} *b.json\n\n{hex}  ./c/./d.json\n"),
            path,
        )
        .unwrap();
        assert_eq!(
            entries.keys().collect::<Vec<_>>(),
            ["a.json", "b.json", "c/./d.json"]
        );
        for (text, line, reason) in [
            ("not a line\n", 1, "is not a sha256sum line"),
            (&format!("{hex} x\n")[..], 1, "is not a sha256sum line"),
            (&format!("{hex}  a\n{hex}  a\n")[..], 2, "lists a twice"),
            (&format!("{hex}  ./\n")[..], 1, "names no file"),
            ("abc  a.json\n", 1, "is not a sha256sum line"),
            (
                &format!("{hex}  a.json\n{hex}  ../a.json\n")[..],
                2,
                "names ../a.json, outside the cache directory",
            ),
            (
                &format!("{hex}  b/../a.json\n")[..],
                1,
                "names b/../a.json, outside the cache directory",
            ),
            (
                &format!("{hex}  /a.json\n")[..],
                1,
                "names /a.json, outside the cache directory",
            ),
        ] {
            let error = parse(text, path).unwrap_err().to_string();
            assert_eq!(
                error,
                format!("CHECKSUMS.sha256: line {line} {reason}"),
                "{text:?}"
            );
        }
        let dir = cache();
        fs::write(dir.path().join(MANIFEST_FILE), "garbage\n").unwrap();
        let error = CacheManifest::load(dir.path()).unwrap_err().to_string();
        assert_eq!(
            error,
            format!(
                "{}: line 1 is not a sha256sum line",
                dir.path().join(MANIFEST_FILE).display()
            )
        );
    }

    /// The reader's digest covers every byte of the stream, those its consumer
    /// read and those `digest_to_end` reads after it stopped.
    #[test]
    fn hashing_reader_digests_the_whole_stream() {
        let data: Vec<u8> = (0..300_000u32).map(|i| (i % 253) as u8).collect();
        let mut reader = HashingReader::new(&data[..]);
        let mut head = vec![0u8; 1000];
        reader.read_exact(&mut head).unwrap();
        assert_eq!(reader.digest_to_end().unwrap(), sha256(&data));

        let mut whole = HashingReader::new(&data[..]);
        let mut sink = Vec::new();
        std::io::copy(&mut whole, &mut sink).unwrap();
        assert_eq!(sink, data);
        assert_eq!(whole.digest_to_end().unwrap(), sha256(&data));
    }

    /// A file hashed in chunks has the digest of its bytes.
    #[test]
    fn file_digest_matches_bytes() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("big");
        let data: Vec<u8> = (0..(READ_CHUNK * 2 + 17))
            .map(|i| (i % 251) as u8)
            .collect();
        File::create(&path).unwrap().write_all(&data).unwrap();
        assert_eq!(file_digest(&path).unwrap(), sha256(&data));
    }
}
