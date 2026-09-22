// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Indexed FASTA access wrapper for VEP.
//!
//! This crate provides fast, thread-safe random access into a reference
//! FASTA using a `.fai` index (samtools faidx format).
//!
//! Design notes:
//! - By default the FASTA file is mmapped once and offsets come from the `.fai`
//!   metadata.
//! - Set `VEP_FASTA_STORAGE=memory` or call [`IndexedFasta::from_path_in_memory`]
//!   to read the FASTA into owned RAM. This avoids runtime major page faults in
//!   service environments where random mmap reads hit remote or block-backed
//!   storage.
//! - All queries are 1-based genomic coordinates, matching VEP conventions.
//! - A small contig-alias map tolerates `chr`-prefixed references.

use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    sync::Arc,
};

use memmap2::Mmap;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum FastaError {
    #[error("FASTA file does not exist: {0}")]
    MissingFasta(String),

    #[error("FASTA index (.fai) does not exist: {0}")]
    MissingFai(String),

    #[error("failed to open FASTA: {0}")]
    OpenFasta(#[from] std::io::Error),

    #[error("failed to memory-map FASTA: {0}")]
    MmapFasta(std::io::Error),

    #[error("failed to read FASTA into memory: {0}")]
    ReadFasta(std::io::Error),

    #[error("failed to read FASTA index: {0}")]
    ReadFai(std::io::Error),

    #[error("invalid FASTA index line: {0}")]
    InvalidFaiLine(String),

    #[error("invalid VEP_FASTA_STORAGE={0:?}: expected 'mmap' or 'memory'")]
    InvalidStorageMode(String),
}

#[derive(Debug, Clone)]
struct FaiRecord {
    len: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// How [`IndexedFasta`] should store reference bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FastaStorageMode {
    /// File-backed mmap. Lower resident memory, but random access can fault at
    /// runtime if the file-backed pages are not resident.
    Mmap,
    /// Owned heap bytes. Higher resident memory, but no runtime file-backed
    /// page faults after startup.
    Memory,
}

impl FastaStorageMode {
    fn from_env() -> Result<Self, FastaError> {
        match std::env::var("VEP_FASTA_STORAGE") {
            Ok(value) => value.parse(),
            Err(std::env::VarError::NotPresent) => Ok(Self::Mmap),
            Err(std::env::VarError::NotUnicode(value)) => Err(FastaError::InvalidStorageMode(
                value.to_string_lossy().into_owned(),
            )),
        }
    }
}

impl std::str::FromStr for FastaStorageMode {
    type Err = FastaError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "mmap" | "map" | "mapped" => Ok(Self::Mmap),
            "memory" | "mem" | "ram" | "owned" => Ok(Self::Memory),
            _ => Err(FastaError::InvalidStorageMode(value.to_string())),
        }
    }
}

enum FastaStorage {
    Mmap(Mmap),
    Memory(Arc<[u8]>),
}

impl std::fmt::Debug for FastaStorage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Mmap(data) => f.debug_struct("Mmap").field("len", &data.len()).finish(),
            Self::Memory(data) => f.debug_struct("Memory").field("len", &data.len()).finish(),
        }
    }
}

impl FastaStorage {
    fn get(&self, idx: usize) -> Option<&u8> {
        match self {
            Self::Mmap(data) => data.get(idx),
            Self::Memory(data) => data.get(idx),
        }
    }

    fn mode(&self) -> FastaStorageMode {
        match self {
            Self::Mmap(_) => FastaStorageMode::Mmap,
            Self::Memory(_) => FastaStorageMode::Memory,
        }
    }
}

/// Thread-safe indexed reference FASTA.
///
/// Uses `.fai` entries to compute byte offsets for random access.
#[derive(Debug)]
pub struct IndexedFasta {
    fasta_path: PathBuf,
    storage: FastaStorage,
    records: HashMap<String, FaiRecord>,
    // Lowercased alias -> canonical contig name in `records`.
    aliases: HashMap<String, String>,
}

impl IndexedFasta {
    /// Open an indexed FASTA at `path` (expects a sibling `.fai`).
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, FastaError> {
        Self::from_path_with_storage(path, FastaStorageMode::from_env()?)
    }

    /// Open an indexed FASTA using a file-backed mmap.
    pub fn from_path_mmap(path: impl AsRef<Path>) -> Result<Self, FastaError> {
        Self::from_path_with_storage(path, FastaStorageMode::Mmap)
    }

    /// Open an indexed FASTA by reading the full file into owned memory.
    pub fn from_path_in_memory(path: impl AsRef<Path>) -> Result<Self, FastaError> {
        Self::from_path_with_storage(path, FastaStorageMode::Memory)
    }

    /// Open an indexed FASTA with an explicit storage mode.
    pub fn from_path_with_storage(
        path: impl AsRef<Path>,
        mode: FastaStorageMode,
    ) -> Result<Self, FastaError> {
        let fasta_path = path.as_ref().to_path_buf();
        if !fasta_path.exists() {
            return Err(FastaError::MissingFasta(fasta_path.display().to_string()));
        }

        let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
        if !fai_path.exists() {
            return Err(FastaError::MissingFai(fai_path.display().to_string()));
        }

        let storage = match mode {
            FastaStorageMode::Mmap => {
                let file = File::open(&fasta_path)?;
                // SAFETY: file is not modified while mapped.
                let mmap = unsafe { Mmap::map(&file).map_err(FastaError::MmapFasta)? };
                FastaStorage::Mmap(mmap)
            }
            FastaStorageMode::Memory => {
                let data = std::fs::read(&fasta_path).map_err(FastaError::ReadFasta)?;
                FastaStorage::Memory(Arc::from(data.into_boxed_slice()))
            }
        };

        let (records, aliases) = load_fai(&fai_path)?;
        Ok(Self {
            fasta_path,
            storage,
            records,
            aliases,
        })
    }

    /// Return the reference base at `chr:pos` (1-based).
    ///
    /// Returns `None` if the contig does not exist or `pos` is out of bounds.
    pub fn base(&self, chr: &str, pos: u64) -> Option<u8> {
        let rec = self.record_for(chr)?;
        if pos == 0 || pos > rec.len {
            return None;
        }

        let pos0 = pos - 1;
        let line = pos0 / rec.line_bases;
        let col = pos0 % rec.line_bases;
        let byte_offset = rec.offset + line * rec.line_width + col;
        let idx = byte_offset as usize;
        let b = *self.storage.get(idx)?;
        if b == b'\n' || b == b'\r' {
            return None;
        }
        Some(b.to_ascii_uppercase())
    }

    /// Extract a sequence of bases for a genomic range [start, end] (1-based, inclusive).
    ///
    /// Returns `None` if the contig does not exist or any position is out of bounds.
    ///
    /// Resolves the FAI record once for the contig, then indexes directly into
    /// the reference bytes for every requested position; looping through
    /// `base()` would re-hash the contig lookup for every base.
    pub fn sequence(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        if start == 0 || end < start {
            return None;
        }
        let rec = self.record_for(chr)?;
        if end > rec.len {
            return None;
        }
        let len = (end - start + 1) as usize;
        let mut seq = Vec::with_capacity(len);
        for pos in start..=end {
            let pos0 = pos - 1;
            let line = pos0 / rec.line_bases;
            let col = pos0 % rec.line_bases;
            let byte_offset = rec.offset + line * rec.line_width + col;
            let b = *self.storage.get(byte_offset as usize)?;
            if b == b'\n' || b == b'\r' {
                return None;
            }
            seq.push(b.to_ascii_uppercase());
        }
        Some(seq)
    }

    fn record_for(&self, chr: &str) -> Option<&FaiRecord> {
        if let Some(rec) = self.records.get(chr) {
            return Some(rec);
        }
        let key = chr.to_ascii_lowercase();
        let canon = self.aliases.get(&key)?;
        self.records.get(canon)
    }

    /// Returns the path to the underlying FASTA file.
    ///
    /// Available to consumers for diagnostics, such as logging the loaded
    /// reference.
    #[allow(dead_code)]
    pub fn fasta_path(&self) -> &Path {
        &self.fasta_path
    }

    /// Returns how this FASTA stores reference bytes.
    pub fn storage_mode(&self) -> FastaStorageMode {
        self.storage.mode()
    }
}

type FaiData = (HashMap<String, FaiRecord>, HashMap<String, String>);

fn load_fai(path: &Path) -> Result<FaiData, FastaError> {
    let file = File::open(path).map_err(FastaError::ReadFai)?;
    let reader = BufReader::new(file);

    let mut records: HashMap<String, FaiRecord> = HashMap::new();
    let mut aliases: HashMap<String, String> = HashMap::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.map_err(FastaError::ReadFai)?;
        if line.trim().is_empty() {
            continue;
        }

        // name \t len \t offset \t line_bases \t line_width
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            return Err(FastaError::InvalidFaiLine(format!(
                "line {}: expected >=5 tab-separated fields: {line}",
                line_no + 1
            )));
        }

        let name = parts[0].to_string();
        let len: u64 = parts[1]
            .parse()
            .map_err(|_| FastaError::InvalidFaiLine(format!("line {}: bad len", line_no + 1)))?;
        let offset: u64 = parts[2]
            .parse()
            .map_err(|_| FastaError::InvalidFaiLine(format!("line {}: bad offset", line_no + 1)))?;
        let line_bases: u64 = parts[3].parse().map_err(|_| {
            FastaError::InvalidFaiLine(format!("line {}: bad line_bases", line_no + 1))
        })?;
        let line_width: u64 = parts[4].parse().map_err(|_| {
            FastaError::InvalidFaiLine(format!("line {}: bad line_width", line_no + 1))
        })?;

        records.insert(
            name.clone(),
            FaiRecord {
                len,
                offset,
                line_bases,
                line_width,
            },
        );

        for alias in contig_aliases(&name) {
            aliases.insert(alias.to_ascii_lowercase(), name.clone());
        }
        aliases.insert(name.to_ascii_lowercase(), name.clone());
    }

    Ok((records, aliases))
}

fn contig_aliases(name: &str) -> Vec<String> {
    let mut out = Vec::new();

    if let Some(stripped) = name.strip_prefix("chr") {
        out.push(stripped.to_string());
    } else {
        out.push(format!("chr{name}"));
    }

    match name {
        "MT" => {
            out.push("M".to_string());
            out.push("chrM".to_string());
            out.push("chrMT".to_string());
        }
        "M" => {
            out.push("MT".to_string());
            out.push("chrM".to_string());
            out.push("chrMT".to_string());
        }
        "chrM" => {
            out.push("MT".to_string());
            out.push("M".to_string());
            out.push("chrMT".to_string());
        }
        "chrMT" => {
            out.push("MT".to_string());
            out.push("M".to_string());
            out.push("chrM".to_string());
        }
        _ => {}
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_indexed_fasta_base_lookup_single_line() {
        let dir = tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let fai = dir.path().join("ref.fa.fai");

        // Header ">1\n" is 3 bytes; sequence starts at offset 3.
        std::fs::write(&fa, b">1\nACACAC\n").unwrap();
        std::fs::write(&fai, b"1\t6\t3\t6\t7\n").unwrap();

        let idx = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(idx.base("1", 1), Some(b'A'));
        assert_eq!(idx.base("1", 2), Some(b'C'));
        assert_eq!(idx.base("1", 6), Some(b'C'));
        assert_eq!(idx.base("1", 7), None);

        assert_eq!(idx.base("chr1", 1), Some(b'A'));
    }

    #[test]
    fn test_indexed_fasta_in_memory_matches_mmap_lookup() {
        let dir = tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let fai = dir.path().join("ref.fa.fai");

        // Header ">1\n" is 3 bytes; sequence starts at offset 3.
        std::fs::write(&fa, b">1\nACACAC\n").unwrap();
        std::fs::write(&fai, b"1\t6\t3\t6\t7\n").unwrap();

        let mmap = IndexedFasta::from_path_mmap(&fa).unwrap();
        let in_memory = IndexedFasta::from_path_in_memory(&fa).unwrap();

        assert_eq!(mmap.storage_mode(), FastaStorageMode::Mmap);
        assert_eq!(in_memory.storage_mode(), FastaStorageMode::Memory);
        for pos in 1..=6 {
            assert_eq!(in_memory.base("1", pos), mmap.base("1", pos));
        }
        assert_eq!(in_memory.sequence("chr1", 2, 5), Some(b"CACA".to_vec()));
    }
}
