// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `.fai` and, for a BGZF file, the `.gzi` of a FASTA, written beside it:
//! [`index_if_missing`] writes the ones a FASTA given without them lacks, as
//! Ensembl VEP's reader does when it first opens such a file, and
//! [`write_fai`] and [`write_gzi`] write one index each, byte for byte the
//! files `samtools faidx` writes.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use noodles::{bgzf, fasta};

use super::{encoding_of, io_error, sibling, Encoding, FastaError, BGZF_EOF_LEN};

/// Writes the indexes the FASTA file at `fasta` lacks, the `.fai` and, for a
/// BGZF file, the `.gzi`, as Ensembl VEP's `_get_fasta_db`
/// (`Variation/Utils/FastaSequence.pm:308-358`) has them written: a
/// `<fasta>.vep.lock` left by an indexing run that did not finish removes the
/// indexes beside the file (`:315-321`), and when none of `.fai`, `.index` and
/// `.gzi` is there the lock is written for the time the indexes take
/// (`:323-338`, `:356`). A directory, a gzip file that is not BGZF and a file
/// with both its indexes are left as they are.
pub fn index_if_missing(fasta: &Path) -> Result<(), FastaError> {
    if !fasta.is_file() {
        return Ok(());
    }
    let lock = sibling(fasta, "vep.lock");
    if lock.exists() {
        for extension in ["fai", "index", "gzi", "vep.lock"] {
            let path = sibling(fasta, extension);
            if path.is_file() {
                std::fs::remove_file(&path)
                    .map_err(io_error("failed to remove FASTA index", &path))?;
            }
        }
    }
    let encoding = encoding_of(fasta)?;
    if encoding == Encoding::Gzip {
        return Ok(());
    }
    let fai = sibling(fasta, "fai");
    let gzi = sibling(fasta, "gzi");
    let write_gzi_too = encoding == Encoding::Bgzf && !gzi.is_file();
    if fai.is_file() && !write_gzi_too {
        return Ok(());
    }
    let index_exists = ["fai", "index", "gzi"]
        .iter()
        .any(|extension| sibling(fasta, extension).exists());
    if !index_exists {
        std::fs::write(&lock, "1\n")
            .map_err(|_| FastaError::LockFile(lock.display().to_string()))?;
    }
    tracing::debug!("Indexing {}", fasta.display());
    if write_gzi_too {
        write_gzi(fasta)?;
    }
    if !fai.is_file() {
        write_fai(fasta)?;
    }
    if !index_exists {
        std::fs::remove_file(&lock).map_err(io_error("failed to remove FASTA lock file", &lock))?;
    }
    Ok(())
}

/// Writes `<fasta>.fai`, the `samtools faidx` index of the plain or BGZF FASTA
/// at `fasta`, a BGZF file's offsets counted in its inflated stream.
pub fn write_fai(fasta: &Path) -> Result<(), FastaError> {
    let file = File::open(fasta).map_err(io_error("failed to open FASTA", fasta))?;
    let reader: Box<dyn BufRead> = match encoding_of(fasta)? {
        Encoding::Bgzf => Box::new(bgzf::io::Reader::new(BufReader::new(file))),
        _ => Box::new(BufReader::new(file)),
    };
    let mut indexer = fasta::io::Indexer::new(reader);
    let mut records = Vec::new();
    while let Some(record) = indexer
        .index_record()
        .map_err(|e| FastaError::Unindexable {
            path: fasta.display().to_string(),
            reason: e.to_string(),
        })?
    {
        records.push(record);
    }
    let index = fasta::fai::Index::from(records);
    let fai = sibling(fasta, "fai");
    let file = File::create(&fai).map_err(io_error("failed to create FASTA index", &fai))?;
    let mut writer = fasta::fai::io::Writer::new(BufWriter::new(file));
    writer
        .write_index(&index)
        .and_then(|()| writer.get_mut().flush())
        .map_err(io_error("failed to write FASTA index", &fai))
}

/// Writes `<fasta>.gzi`, the `bgzip` index of the BGZF file at `fasta`: the
/// compressed and inflated offsets of every non-empty block after the first,
/// read from the block headers and trailers without inflating the data.
pub fn write_gzi(fasta: &Path) -> Result<(), FastaError> {
    let unindexable = |reason: String| FastaError::Unindexable {
        path: fasta.display().to_string(),
        reason,
    };
    let mut file = File::open(fasta).map_err(io_error("failed to open FASTA", fasta))?;
    let length = file
        .metadata()
        .map_err(io_error("failed to read FASTA", fasta))?
        .len();
    let mut entries = Vec::new();
    let (mut compressed, mut inflated) = (0u64, 0u64);
    let mut header = [0u8; super::BGZF_HEADER_LEN];
    while compressed < length {
        file.seek(SeekFrom::Start(compressed))
            .and_then(|_| file.read_exact(&mut header))
            .map_err(|_| unindexable(format!("truncated block header at offset {compressed}")))?;
        let Some(size) = super::bgzf_header_len(&header) else {
            return Err(unindexable(format!(
                "gzip but not BGZF at offset {compressed}"
            )));
        };
        if size < BGZF_EOF_LEN {
            return Err(unindexable(format!(
                "a {size}-byte block at offset {compressed}, shorter than the {BGZF_EOF_LEN}-byte minimum of BGZF"
            )));
        }
        let size = size as u64;
        let mut trailer = [0u8; 4];
        file.seek(SeekFrom::Start(compressed + size - 4))
            .and_then(|_| file.read_exact(&mut trailer))
            .map_err(|_| unindexable(format!("truncated block at offset {compressed}")))?;
        let data_len = u64::from(u32::from_le_bytes(trailer));
        if data_len > 0 && compressed > 0 {
            entries.push((compressed, inflated));
        }
        compressed += size;
        inflated += data_len;
    }
    let gzi = sibling(fasta, "gzi");
    bgzf::gzi::fs::write(&gzi, &bgzf::gzi::Index::from(entries))
        .map_err(io_error("failed to write BGZF index", &gzi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IndexedFasta;
    use std::path::PathBuf;

    /// The GRCh37-hgvs corpus's reference, bgzipped, with the `.fai` and `.gzi`
    /// `samtools faidx` wrote for it.
    fn committed_reference() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh37-hgvs/reference.fa.gz")
    }

    const PLAIN: &[u8] = b">1 dna:chromosome\nACGTACGTAC\nGTACGTACGT\nACG\n>MT\nTTTTGGGG\nCC\n";
    const PLAIN_FAI: &str = "1\t23\t18\t10\t11\nMT\t10\t48\t8\t9\n";

    /// A plain FASTA given without its `.fai` gets the one `samtools faidx`
    /// writes, reads through it, and leaves no lock file behind.
    #[test]
    fn a_plain_fasta_without_its_fai_is_indexed() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, PLAIN).unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(
            std::fs::read_to_string(dir.path().join("ref.fa.fai")).unwrap(),
            PLAIN_FAI
        );
        assert!(!dir.path().join("ref.fa.vep.lock").exists());
        let fasta = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(fasta.sequence("1", 9, 13), Some(b"ACGTA".to_vec()));
        assert_eq!(fasta.sequence("MT", 8, 10), Some(b"GCC".to_vec()));
    }

    /// A BGZF FASTA given without its indexes gets the `.fai` and `.gzi`
    /// `samtools faidx` wrote for the committed reference, byte for byte; one
    /// with its `.fai` alone gets the `.gzi`.
    #[test]
    fn a_bgzf_fasta_without_its_indexes_is_indexed_as_samtools_indexes_it() {
        let reference = committed_reference();
        for kept in [None, Some("fai")] {
            let dir = tempfile::tempdir().unwrap();
            let fa = dir.path().join("reference.fa.gz");
            std::fs::copy(&reference, &fa).unwrap();
            if let Some(extension) = kept {
                std::fs::copy(sibling(&reference, extension), sibling(&fa, extension)).unwrap();
            }
            index_if_missing(&fa).unwrap();
            for extension in ["fai", "gzi"] {
                assert_eq!(
                    std::fs::read(sibling(&fa, extension)).unwrap(),
                    std::fs::read(sibling(&reference, extension)).unwrap(),
                    "{extension} with {kept:?} kept"
                );
            }
            assert!(!sibling(&fa, "vep.lock").exists());
            IndexedFasta::from_path(&fa).unwrap();
        }
    }

    /// A lock file left by an indexing run that did not finish removes the
    /// index beside the FASTA, which is written again (`FastaSequence.pm:315-321`).
    #[test]
    fn a_lock_file_left_behind_has_the_index_written_again() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, PLAIN).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "1\t5\t18\t10\t11\n").unwrap();
        std::fs::write(dir.path().join("ref.fa.vep.lock"), "1\n").unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(
            std::fs::read_to_string(dir.path().join("ref.fa.fai")).unwrap(),
            PLAIN_FAI
        );
        assert!(!dir.path().join("ref.fa.vep.lock").exists());
    }

    /// An index already beside the FASTA is read as it is, never rewritten.
    #[test]
    fn an_existing_fai_is_left_as_it_is() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, PLAIN).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "1\t5\t18\t10\t11\n").unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(
            std::fs::read_to_string(dir.path().join("ref.fa.fai")).unwrap(),
            "1\t5\t18\t10\t11\n"
        );
    }

    /// Where the lock file cannot be written, the run ends with Ensembl VEP's
    /// words (`FastaSequence.pm:334`).
    #[cfg(unix)]
    #[test]
    fn an_unwritable_directory_is_refused_with_the_lock_file_named() {
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, PLAIN).unwrap();
        std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o555)).unwrap();
        let result = index_if_missing(&fa);
        std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o755)).unwrap();
        if std::fs::write(dir.path().join("probe"), "").is_ok() {
            // Running as root, which writes into the directory anyway.
            return;
        }
        assert_eq!(
            result.unwrap_err().to_string(),
            format!(
                "Could not write to FASTA lock file {}",
                dir.path().join("ref.fa.vep.lock").display()
            )
        );
    }

    /// A FASTA whose sequence lines differ in length cannot be indexed, and the
    /// refusal names the file.
    #[test]
    fn a_fasta_with_uneven_lines_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, b">1\nACGT\nAC\nACGT\n").unwrap();
        let error = index_if_missing(&fa).unwrap_err().to_string();
        assert!(
            error.starts_with(&format!("failed to index {}: ", fa.display())),
            "{error}"
        );
    }

    /// A block stating a size shorter than the 28-byte minimum is refused with
    /// its offset.
    #[test]
    fn a_block_shorter_than_the_minimum_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        for bsize in [26u16, 2] {
            let mut header = [0u8; 18];
            header[..4].copy_from_slice(&[0x1f, 0x8b, 8, 0x04]);
            header[10..12].copy_from_slice(&6u16.to_le_bytes());
            header[12..16].copy_from_slice(&[b'B', b'C', 2, 0]);
            header[16..18].copy_from_slice(&bsize.to_le_bytes());
            let crafted = dir.path().join(format!("bsize{bsize}.fa.gz"));
            std::fs::write(&crafted, header).unwrap();
            let error = write_gzi(&crafted).unwrap_err().to_string();
            assert!(
                error.contains(&format!("{}-byte block at offset 0", bsize + 1))
                    && error.contains("28-byte minimum"),
                "{error}"
            );
        }
    }
}
