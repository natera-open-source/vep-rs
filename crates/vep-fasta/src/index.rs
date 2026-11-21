// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `.fai` and, for a BGZF file, the `.gzi` of a FASTA, written beside it:
//! [`index_if_missing`] writes the ones a FASTA given without them lacks, as
//! Ensembl VEP's reader does when it first opens such a file, [`open_indexed`]
//! opens the file once they are written, and [`write_fai`] and [`write_gzi`]
//! write one index each, byte for byte the files `samtools faidx` writes.

use std::collections::HashSet;
use std::fs::{File, OpenOptions};
use std::hash::{BuildHasher, RandomState};
use std::io::{BufRead, BufReader, BufWriter, ErrorKind, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use noodles::bgzf;

use super::{encoding_of, io_error, sibling, Encoding, FastaError, IndexedFasta, BGZF_EOF_LEN};

/// Writes the indexes the FASTA file at `fasta` lacks, the `.fai` and, for a
/// BGZF file, the `.gzi`, as Ensembl VEP's `_get_fasta_db`
/// (`Variation/Utils/FastaSequence.pm:308-358`) has them written: a
/// `<fasta>.vep.lock` left by an indexing run that did not finish removes the
/// indexes beside the file (`:315-321`), and when none of `.fai`, `.index` and
/// `.gzi` is there the lock is written for the time the indexes take
/// (`:323-338`, `:356`). A directory, a gzip file that is not BGZF and a file
/// with both its indexes are left as they are. A file that cannot be indexed
/// is refused with no lock and no part of an index left beside it.
pub fn index_if_missing(fasta: &Path) -> Result<(), FastaError> {
    write_missing(fasta, true)
}

/// How many times [`open_indexed`] writes the missing indexes and opens the
/// file before the open's refusal stands.
const OPEN_ATTEMPTS: usize = 32;

/// Opens the FASTA at `fasta` once [`index_if_missing`] has written the
/// indexes it lacks, as `vep` opens the file `--fasta` names. A run that
/// starts on the same file while this one indexes takes this run's lock for
/// one left behind and removes the indexes beside it, as this run does with a
/// lock it finds. An index that is gone at the open is therefore written
/// again, a lock found by then left as the other run's, and the file opened
/// again, up to `OPEN_ATTEMPTS` times: runs that meet on one FASTA wait for
/// no lock, and each ends with the file open or with the last open's refusal.
pub fn open_indexed(fasta: &Path) -> Result<IndexedFasta, FastaError> {
    let mut attempt = 1;
    loop {
        write_missing(fasta, attempt == 1)?;
        match IndexedFasta::from_path(fasta) {
            Err(error) if attempt < OPEN_ATTEMPTS && index_is_gone(&error) => attempt += 1,
            opened => return opened,
        }
    }
}

/// Whether `error` is the open's refusal of an index that is not there.
fn index_is_gone(error: &FastaError) -> bool {
    match error {
        FastaError::MissingFai(_) | FastaError::MissingGzi(_) => true,
        FastaError::Io { source, .. } => source.kind() == ErrorKind::NotFound,
        _ => false,
    }
}

/// [`index_if_missing`], with a lock found beside the FASTA taken for one left
/// behind, which removes the indexes, only when `lock_is_left_behind`.
fn write_missing(fasta: &Path, lock_is_left_behind: bool) -> Result<(), FastaError> {
    if !fasta.is_file() {
        return Ok(());
    }
    let lock = sibling(fasta, "vep.lock");
    if lock_is_left_behind && lock.exists() {
        // A removal that fails ends no run, as Ensembl VEP's `unlink` ends none
        // (`:319`): an index that stays is read, and one that is missing is
        // written or refused below.
        for extension in ["fai", "index", "gzi", "vep.lock"] {
            let _ = std::fs::remove_file(sibling(fasta, extension));
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
        write_lock(&lock)?;
    }
    tracing::debug!("Indexing {}", fasta.display());
    let written = (|| {
        if write_gzi_too {
            write_gzi(fasta)?;
        }
        if !fai.is_file() {
            write_fai(fasta)?;
        }
        Ok(())
    })();
    if !index_exists {
        // An index is put in place whole or not at all, so the lock marks
        // nothing once the writers have returned, whether or not they wrote;
        // one that cannot be removed ends no run (`:356`).
        let _ = std::fs::remove_file(&lock);
    }
    written
}

/// Writes the lock file `lock`, created new so that a link standing there is
/// not written through. A lock that is there is another run's, written since
/// this run looked for one, and is left as it is: Ensembl VEP writes over it
/// (`:334`), and neither run waits for the other.
fn write_lock(lock: &Path) -> Result<(), FastaError> {
    let refused = |source| FastaError::LockFile {
        lock: lock.display().to_string(),
        source,
    };
    match OpenOptions::new().write(true).create_new(true).open(lock) {
        Ok(mut file) => file.write_all(b"1\n").map_err(|error| {
            let _ = std::fs::remove_file(lock);
            refused(error)
        }),
        Err(error) if error.kind() == ErrorKind::AlreadyExists => Ok(()),
        Err(error) => Err(refused(error)),
    }
}

/// Writes `<fasta>.fai`, the `samtools faidx` index of the plain or BGZF FASTA
/// at `fasta`, a BGZF file's offsets counted in its inflated stream.
pub fn write_fai(fasta: &Path) -> Result<(), FastaError> {
    let file = File::open(fasta).map_err(io_error("failed to open FASTA", fasta))?;
    let file = BufReader::new(file);
    let records = match encoding_of(fasta)? {
        Encoding::Bgzf => fai_records(bgzf::io::Reader::new(file), fasta),
        _ => fai_records(file, fasta),
    }
    .map_err(|reason| FastaError::Unindexable {
        path: fasta.display().to_string(),
        reason,
    })?;
    let fai = sibling(fasta, "fai");
    write_new(&fai, |file| {
        for record in &records {
            file.write_all(&record.name)?;
            writeln!(
                file,
                "\t{}\t{}\t{}\t{}",
                record.bases, record.offset, record.line_bases, record.line_bytes
            )?;
        }
        Ok(())
    })
    .map_err(io_error("failed to write FASTA index", &fai))
}

/// One record of a `.fai`: a sequence's name, its bases, the offset of its
/// first base, and the bases and the bytes of its first line.
#[derive(Default)]
struct FaiLine {
    name: Vec<u8>,
    bases: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

/// The bytes of a FASTA, read one at a time or a line at a time, with the
/// offset of the next byte in the uncompressed text.
struct FastaBytes<R> {
    reader: R,
    offset: u64,
}

impl<R: BufRead> FastaBytes<R> {
    /// The bytes the reader holds, filled when it holds none; empty at the end
    /// of the text.
    fn buffered(&mut self) -> std::io::Result<&[u8]> {
        loop {
            match self.reader.fill_buf() {
                Ok(_) => break,
                Err(error) if error.kind() == ErrorKind::Interrupted => {}
                Err(error) => return Err(error),
            }
        }
        self.reader.fill_buf()
    }

    /// The next byte, `None` at the end of the text.
    fn next(&mut self) -> std::io::Result<Option<u8>> {
        let byte = self.buffered()?.first().copied();
        if byte.is_some() {
            self.reader.consume(1);
            self.offset += 1;
        }
        Ok(byte)
    }

    /// Reads through the next line feed, or to the end of the text: how many
    /// bytes stand before the line feed, and how many of them are printing
    /// characters.
    fn rest_of_line(&mut self) -> std::io::Result<(u64, u64)> {
        let (mut bytes, mut printing) = (0, 0);
        loop {
            let buffered = self.buffered()?;
            if buffered.is_empty() {
                return Ok((bytes, printing));
            }
            let line_feed = buffered.iter().position(|&byte| byte == b'\n');
            let line = &buffered[..line_feed.unwrap_or(buffered.len())];
            bytes += line.len() as u64;
            printing += line.iter().filter(|byte| byte.is_ascii_graphic()).count() as u64;
            let read = line.len() + usize::from(line_feed.is_some());
            self.reader.consume(read);
            self.offset += read as u64;
            if line_feed.is_some() {
                return Ok((bytes, printing));
            }
        }
    }
}

/// The `.fai` records of the FASTA text in `reader`, counted as htslib's
/// `fai_build_core` (`faidx.c`) counts them, so the index is the one
/// `samtools faidx` writes: a sequence's bases are its printing characters, a
/// line's bytes every byte through its line feed, and a last line without one
/// counts as if it had it; a name ends at its first white space, and white
/// space between the header mark and the name is passed over; blank lines may
/// follow a sequence; a header with no sequence line before another header
/// gets no record, and a name given again keeps its first record. The reason
/// is returned for a layout no index describes: a line longer than the first
/// of its sequence, text after a blank or a shorter line, text before the
/// first header, or a last header with no sequence line.
fn fai_records<R: BufRead>(reader: R, fasta: &Path) -> Result<Vec<FaiLine>, String> {
    enum At {
        /// Before the first header, or after a blank line or a line shorter
        /// than the first of its sequence: a header or a blank line is due.
        Between,
        /// Behind a header mark.
        Name,
        /// At the start of a line of the sequence, or of the next header.
        Sequence,
    }
    const ONE_WIDTH: &str = "wrap every line of a sequence but its last to one width";
    let unread = |error: std::io::Error| error.to_string();
    let named = |record: &FaiLine| String::from_utf8_lossy(&record.name).into_owned();
    let mut bytes = FastaBytes { reader, offset: 0 };
    let mut records: Vec<FaiLine> = Vec::new();
    let mut names: HashSet<Vec<u8>> = HashSet::new();
    let mut keep = |record: FaiLine| {
        if names.insert(record.name.clone()) {
            records.push(record);
        } else {
            tracing::warn!(
                "{} names the sequence {} again at byte {}: its first sequence of that name is the one indexed",
                fasta.display(),
                String::from_utf8_lossy(&record.name),
                record.offset
            );
        }
    };
    let mut at = At::Between;
    let mut line = 1u64;
    let mut record = FaiLine::default();
    let (mut opened, mut has_sequence) = (false, false);
    while let Some(byte) = bytes.next().map_err(unread)? {
        match at {
            At::Between => match byte {
                b'>' => {
                    opened = true;
                    at = At::Name;
                }
                b'\n' => line += 1,
                b'\r' => {
                    if bytes.next().map_err(unread)? != Some(b'\n') {
                        return Err(format!(
                            "line {line} holds a carriage return that no line feed follows"
                        ));
                    }
                    line += 1;
                }
                _ if opened => {
                    return Err(format!(
                        "line {line} follows a blank line or a line shorter than the first of sequence {}, where only a header line may stand: {ONE_WIDTH}",
                        named(&record)
                    ));
                }
                _ => {
                    return Err(format!(
                        "line {line} is not a header line: a FASTA opens each sequence with `>` and its name"
                    ));
                }
            },
            At::Name => {
                if has_sequence {
                    keep(std::mem::take(&mut record));
                    has_sequence = false;
                }
                record = FaiLine::default();
                let mut next = Some(byte);
                while let Some(byte) = next {
                    if !is_space(byte) {
                        record.name.push(byte);
                    } else if !record.name.is_empty() || byte == b'\n' {
                        break;
                    }
                    next = bytes.next().map_err(unread)?;
                }
                // htslib holds the name as a C string (`faidx.c:100`), so the
                // name ends before the first NUL byte in it.
                if let Some(nul) = record.name.iter().position(|&byte| byte == 0) {
                    record.name.truncate(nul);
                }
                match next {
                    None => {
                        return Err(format!(
                            "sequence {}, the last of the file, has no sequence line",
                            named(&record)
                        ));
                    }
                    Some(b'\n') => {}
                    Some(_) => {
                        bytes.rest_of_line().map_err(unread)?;
                    }
                }
                record.offset = bytes.offset;
                at = At::Sequence;
                line += 1;
            }
            At::Sequence => {
                if byte == b'\n' {
                    at = At::Between;
                    line += 1;
                    continue;
                }
                if byte == b'>' {
                    at = At::Name;
                    continue;
                }
                has_sequence = true;
                let (rest, printing) = bytes.rest_of_line().map_err(unread)?;
                let line_bytes = rest + 2;
                let line_bases = printing + u64::from(byte.is_ascii_graphic());
                record.bases += line_bases;
                if record.line_bytes == 0 {
                    record.line_bytes = line_bytes;
                    record.line_bases = line_bases;
                } else if line_bytes < record.line_bytes {
                    at = At::Between;
                } else if line_bytes > record.line_bytes {
                    return Err(format!(
                        "line {line} of sequence {} is longer than its first line: {ONE_WIDTH}",
                        named(&record)
                    ));
                }
                line += 1;
            }
        }
    }
    if !has_sequence {
        return Err(if opened {
            format!(
                "sequence {}, the last of the file, has no sequence line",
                named(&record)
            )
        } else {
            "the file holds no sequence".to_string()
        });
    }
    keep(record);
    Ok(records)
}

/// Whether `byte` is white space as C's `isspace` has it in its default
/// locale, the test htslib ends a sequence name on.
fn is_space(byte: u8) -> bool {
    byte == b' ' || (b'\t'..=b'\r').contains(&byte)
}

/// Writes `<fasta>.gzi`, the `bgzip` index of the BGZF file at `fasta`: the
/// compressed and inflated offsets of every block that holds data but the
/// first of them, read from the block headers and trailers without inflating
/// the data.
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
    let mut first_data_block_seen = false;
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
        if data_len > 0 {
            if first_data_block_seen {
                entries.push((compressed, inflated));
            }
            first_data_block_seen = true;
        }
        compressed += size;
        inflated += data_len;
    }
    let gzi = sibling(fasta, "gzi");
    let index = bgzf::gzi::Index::from(entries);
    write_new(&gzi, |file| {
        bgzf::gzi::io::Writer::new(file).write_index(&index)
    })
    .map_err(io_error("failed to write BGZF index", &gzi))
}

/// Writes `dst` through `write` into a file created new beside it, then
/// renames that file over `dst`, so a link or a reader at `dst` never sees a
/// partial index and no file outside the FASTA's directory is written through
/// a link. The file created is removed when the write or the rename fails.
fn write_new(
    dst: &Path,
    write: impl FnOnce(&mut BufWriter<File>) -> std::io::Result<()>,
) -> std::io::Result<()> {
    let (temporary, file) = create_beside(dst)?;
    let mut writer = BufWriter::new(file);
    let result = write(&mut writer)
        .and_then(|()| writer.flush())
        .and_then(|()| {
            drop(writer);
            std::fs::rename(&temporary, dst)
        });
    if result.is_err() {
        let _ = std::fs::remove_file(&temporary);
    }
    result
}

/// A file created new beside `dst`, as `<dst>.<token>.tmp`. The token is
/// drawn anew for every call from the process's randomly keyed hasher, so two
/// runs share no name, on one host or on two that share the directory; a name
/// that is taken all the same is passed over for the next, and the file it
/// names is left as it is.
fn create_beside(dst: &Path) -> std::io::Result<(PathBuf, File)> {
    static CALLS: AtomicU64 = AtomicU64::new(0);
    loop {
        let call = (std::process::id(), CALLS.fetch_add(1, Ordering::Relaxed));
        let token = RandomState::new().hash_one(call);
        let temporary = PathBuf::from(format!("{}.{token:016x}.tmp", dst.display()));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
        {
            Ok(file) => return Ok((temporary, file)),
            Err(error) if error.kind() == ErrorKind::AlreadyExists => {}
            Err(error) => return Err(error),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

    /// A link at the index's path is replaced by the index, and the file it
    /// points to is never written.
    #[cfg(unix)]
    #[test]
    fn a_link_at_the_index_path_is_replaced_not_followed() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, PLAIN).unwrap();
        let elsewhere = dir.path().join("elsewhere");
        std::os::unix::fs::symlink(&elsewhere, dir.path().join("ref.fa.fai")).unwrap();
        index_if_missing(&fa).unwrap();
        assert!(!elsewhere.exists());
        let fai = dir.path().join("ref.fa.fai");
        assert!(!std::fs::symlink_metadata(&fai)
            .unwrap()
            .file_type()
            .is_symlink());
        assert_eq!(std::fs::read_to_string(fai).unwrap(), PLAIN_FAI);
    }

    /// Where the lock file cannot be written, the run ends with Ensembl VEP's
    /// words (`FastaSequence.pm:334`), then the cause and what to do about it,
    /// and nothing is written beside the FASTA.
    #[cfg(unix)]
    #[test]
    fn an_unwritable_directory_is_refused_with_the_lock_file_named() {
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, PLAIN).unwrap();
        std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o555)).unwrap();
        // Root writes into the directory anyway, and the state cannot be made.
        let unwritable = std::fs::write(dir.path().join("probe"), "").is_err();
        let result = index_if_missing(&fa);
        std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o755)).unwrap();
        if !unwritable {
            return;
        }
        assert_eq!(
            result.unwrap_err().to_string(),
            format!(
                "Could not write to FASTA lock file {}: Permission denied (os error 13); index the FASTA with samtools faidx, or link it into a directory the run can write",
                dir.path().join("ref.fa.vep.lock").display()
            )
        );
        assert_eq!(names_in(dir.path()), ["ref.fa"]);
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

    /// `count` bases of a made-up four-base unit, repeated.
    fn made_up(count: usize) -> String {
        "ACGT".chars().cycle().take(count).collect()
    }

    /// The `.fai` written for a plain FASTA holding `text`.
    fn fai_of(text: &str) -> Result<String, FastaError> {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, text).unwrap();
        index_if_missing(&fa)?;
        Ok(std::fs::read_to_string(dir.path().join("ref.fa.fai")).unwrap())
    }

    /// White space before a line end is no base: a sequence's length and its
    /// bases per line count printing characters, its bytes per line every byte,
    /// as `samtools faidx` writes them, and each base reads as the file holds it.
    #[test]
    fn white_space_ending_a_sequence_line_is_not_counted_as_a_base() {
        let (full, short) = (made_up(8), made_up(5));
        for pad in [" ", "\t"] {
            let every_line = format!(">21\n{}", format!("{full}{pad}\n").repeat(3));
            assert_eq!(
                fai_of(&every_line).unwrap(),
                "21\t24\t4\t8\t10\n",
                "{pad:?} ends every line"
            );
            let first_record =
                format!(">decoy\n{full}\n{full}\n{short}{pad}\n>21\n{full}\n{full}\n");
            assert_eq!(
                fai_of(&first_record).unwrap(),
                "decoy\t21\t7\t8\t9\n21\t16\t36\t8\t9\n",
                "{pad:?} ends the last line of the first record"
            );
            let only_record = format!(">21\n{full}\n{full}\n{short}{pad}\n");
            assert_eq!(
                fai_of(&only_record).unwrap(),
                "21\t21\t4\t8\t9\n",
                "{pad:?} ends the last line of the only record"
            );
        }
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, format!(">21\n{}", format!("{full} \n").repeat(3))).unwrap();
        index_if_missing(&fa).unwrap();
        let fasta = IndexedFasta::from_path(&fa).unwrap();
        for (at, base) in made_up(24).bytes().enumerate() {
            assert_eq!(
                fasta.base("21", at as u64 + 1),
                Some(base),
                "base {}",
                at + 1
            );
        }
    }

    /// A name a FASTA gives two sequences is indexed once, under the first, as
    /// `samtools faidx` indexes it, and the file reads through that index.
    #[test]
    fn a_name_given_twice_is_indexed_under_its_first_sequence() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let text = format!(
            ">21\n{}\n>MT\n{}\n>21\n{}\n",
            made_up(8),
            made_up(4),
            "T".repeat(6)
        );
        std::fs::write(&fa, text).unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(
            std::fs::read_to_string(dir.path().join("ref.fa.fai")).unwrap(),
            "21\t8\t4\t8\t9\nMT\t4\t17\t4\t5\n"
        );
        let fasta = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(fasta.sequence("21", 1, 8), Some(made_up(8).into_bytes()));
        assert_eq!(fasta.sequence_length("21"), Some(8));
    }

    /// Blank lines after a sequence, a header with no sequence line before
    /// another header and white space after the header mark are indexed as
    /// `samtools faidx` indexes them.
    #[test]
    fn blank_lines_an_empty_sequence_and_a_spaced_header_are_indexed() {
        let (full, short, other) = (made_up(8), made_up(3), made_up(4));
        let layouts = [
            (
                "a blank line after a short last line",
                format!(">1\n{full}\n{full}\n{short}\n\n>2\n{other}\n"),
                "1\t19\t3\t8\t9\n2\t4\t29\t4\t5\n",
            ),
            (
                "two blank lines between two sequences",
                format!(">1\n{full}\n{full}\n\n\n>2\n{other}\n"),
                "1\t16\t3\t8\t9\n2\t4\t26\t4\t5\n",
            ),
            (
                "two blank lines at the end of the file",
                format!(">1\n{full}\n{full}\n\n\n"),
                "1\t16\t3\t8\t9\n",
            ),
            (
                "a blank first line",
                format!("\n>1\n{full}\n{full}\n"),
                "1\t16\t4\t8\t9\n",
            ),
            (
                "a header with no sequence line before another header",
                format!(">1\n{full}\n>empty\n>2\n{other}\n"),
                "1\t8\t3\t8\t9\n2\t4\t22\t4\t5\n",
            ),
            (
                "a space after the header mark",
                format!("> 1 made up\n{full}\n{short}\n"),
                "1\t11\t12\t8\t9\n",
            ),
        ];
        for (layout, text, fai) in layouts {
            let dir = tempfile::tempdir().unwrap();
            let fa = dir.path().join("ref.fa");
            std::fs::write(&fa, &text).unwrap();
            index_if_missing(&fa).unwrap_or_else(|error| panic!("{layout}: {error}"));
            assert_eq!(
                std::fs::read_to_string(dir.path().join("ref.fa.fai")).unwrap(),
                fai,
                "{layout}"
            );
            let fasta = IndexedFasta::from_path(&fa).unwrap();
            assert_eq!(fasta.base("1", 8), Some(b'T'), "{layout}");
        }
    }

    /// A sequence of one line that the file ends without a line end has the
    /// bytes per line of a whole line, its bases and one, as `samtools faidx`
    /// writes it.
    #[test]
    fn one_line_without_a_line_end_has_the_bytes_of_a_whole_line() {
        assert_eq!(
            fai_of(&format!(">1\n{}", made_up(8))).unwrap(),
            "1\t8\t3\t8\t9\n"
        );
        assert_eq!(
            fai_of(&format!(">1\n{} ", made_up(8))).unwrap(),
            "1\t8\t3\t8\t10\n"
        );
    }

    /// Windows line ends, a name followed by a description, a last line
    /// without its line end and a lone header mark ending the file are indexed
    /// as `samtools faidx` indexes them.
    #[test]
    fn windows_line_ends_and_a_last_line_without_its_end_are_indexed() {
        let (full, short, other) = (made_up(8), made_up(3), made_up(4));
        let windows = format!(">1 x\r\n{full}\r\n{full}\r\n{short}\r\n\r\n>2\r\n{other}\r\n");
        assert_eq!(
            fai_of(&windows).unwrap(),
            "1\t19\t6\t8\t10\n2\t4\t37\t4\t6\n"
        );
        assert_eq!(
            fai_of(&format!(">a b c\n{other}\n>x\ty\n{other}\n")).unwrap(),
            "a\t4\t7\t4\t5\nx\t4\t17\t4\t5\n"
        );
        assert_eq!(
            fai_of(&format!(">1\n{full}\n{short}")).unwrap(),
            "1\t11\t3\t8\t9\n"
        );
        assert_eq!(
            fai_of(&format!(">1\n{full}\n>")).unwrap(),
            "1\t8\t3\t8\t9\n"
        );
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, windows).unwrap();
        index_if_missing(&fa).unwrap();
        let fasta = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(
            fasta.sequence("1", 7, 10),
            Some(made_up(19).as_bytes()[6..10].to_vec())
        );
    }

    /// A layout no `.fai` describes is refused with the file, the line and the
    /// rule it breaks named, as `samtools faidx` refuses it.
    #[test]
    fn a_layout_no_index_describes_is_refused_with_its_line_named() {
        let (full, short) = (made_up(8), made_up(3));
        let one_width = "wrap every line of a sequence but its last to one width";
        let layouts = [
            (
                format!(">21\n{full}\n{full} \n{full}\n"),
                format!("line 3 of sequence 21 is longer than its first line: {one_width}"),
            ),
            (
                format!(">1\n{short}\n{full}\n"),
                format!("line 3 of sequence 1 is longer than its first line: {one_width}"),
            ),
            (
                format!(">1\n{full}\n{short}\n{full}\n"),
                format!(
                    "line 4 follows a blank line or a line shorter than the first of sequence 1, where only a header line may stand: {one_width}"
                ),
            ),
            (
                format!(">1\n\n{full}\n"),
                format!(
                    "line 3 follows a blank line or a line shorter than the first of sequence 1, where only a header line may stand: {one_width}"
                ),
            ),
            (
                format!("{full}\n>1\n{full}\n"),
                "line 1 is not a header line: a FASTA opens each sequence with `>` and its name"
                    .to_string(),
            ),
            (
                format!(">1\n{full}\n>2\n"),
                "sequence 2, the last of the file, has no sequence line".to_string(),
            ),
            (
                ">1".to_string(),
                "sequence 1, the last of the file, has no sequence line".to_string(),
            ),
            (String::new(), "the file holds no sequence".to_string()),
            (
                format!(">1\n{full}\n{short}\n\r{short}\n"),
                "line 4 holds a carriage return that no line feed follows".to_string(),
            ),
        ];
        for (text, reason) in layouts {
            let dir = tempfile::tempdir().unwrap();
            let fa = dir.path().join("ref.fa");
            std::fs::write(&fa, &text).unwrap();
            assert_eq!(
                index_if_missing(&fa).unwrap_err().to_string(),
                format!("failed to index {}: {reason}", fa.display()),
                "{text:?}"
            );
        }
    }

    /// Two BGZF files joined end to end get the `.fai` and the `.gzi`
    /// `samtools faidx` writes for the joined file, the end-of-file block of the
    /// first file without an entry, and both sequences read through them.
    #[test]
    fn two_bgzf_files_joined_are_indexed_and_read() {
        use crate::tests::{bases, bgzf_and_gzi, fasta_and_fai};
        let dir = tempfile::tempdir().unwrap();
        let (one, two) = (bases(21, 2_500), bases(22, 1_700));
        let (first, _) = fasta_and_fai(&[("1", &one)]);
        let (whole, fai) = fasta_and_fai(&[("1", &one), ("2", &two)]);
        let (first_bgzf, first_gzi) = bgzf_and_gzi(&first, 1_000);
        let (second_bgzf, second_gzi) = bgzf_and_gzi(&whole[first.len()..], 1_000);
        let (joined_at, text_at) = (first_bgzf.len() as u64, first.len() as u64);
        let mut pairs = first_gzi.as_ref().to_vec();
        pairs.push((joined_at, text_at));
        pairs.extend(
            second_gzi
                .as_ref()
                .iter()
                .map(|(compressed, text)| (joined_at + compressed, text_at + text)),
        );
        let mut gzi = (pairs.len() as u64).to_le_bytes().to_vec();
        for (compressed, text) in pairs {
            gzi.extend_from_slice(&compressed.to_le_bytes());
            gzi.extend_from_slice(&text.to_le_bytes());
        }
        let fa = dir.path().join("ref.fa.gz");
        std::fs::write(&fa, [first_bgzf, second_bgzf].concat()).unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(std::fs::read(sibling(&fa, "gzi")).unwrap(), gzi);
        assert_eq!(std::fs::read_to_string(sibling(&fa, "fai")).unwrap(), fai);
        let fasta = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(
            fasta.sequence("1", 1, 2_500),
            Some(one.to_ascii_uppercase())
        );
        assert_eq!(
            fasta.sequence("2", 1, 1_700),
            Some(two.to_ascii_uppercase())
        );
    }

    /// The names of the entries of `dir`, sorted.
    fn names_in(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = std::fs::read_dir(dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }

    /// A FASTA that cannot be indexed is left with no lock file and no part of
    /// an index beside it; a `.gzi` written whole before the `.fai` failed
    /// stays, and a `.fai` put beside the file afterwards is read as it is.
    #[test]
    fn a_fasta_that_cannot_be_indexed_keeps_no_lock_beside_it() {
        use crate::tests::bgzf_and_gzi;
        let uneven = format!(">1\n{}\n{}\n{}\n", made_up(8), made_up(3), made_up(8));
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, &uneven).unwrap();
        index_if_missing(&fa).unwrap_err();
        assert_eq!(names_in(dir.path()), ["ref.fa"]);
        let by_hand = "1\t11\t3\t8\t9\n";
        std::fs::write(dir.path().join("ref.fa.fai"), by_hand).unwrap();
        let fasta = open_indexed(&fa).unwrap();
        assert_eq!(fasta.base("1", 11), Some(b'G'));
        assert_eq!(
            std::fs::read_to_string(dir.path().join("ref.fa.fai")).unwrap(),
            by_hand
        );
        assert_eq!(names_in(dir.path()), ["ref.fa", "ref.fa.fai"]);

        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa.gz");
        std::fs::write(&fa, bgzf_and_gzi(uneven.as_bytes(), 16).0).unwrap();
        index_if_missing(&fa).unwrap_err();
        assert_eq!(names_in(dir.path()), ["ref.fa.gz", "ref.fa.gz.gzi"]);
    }

    /// An index that cannot be put in place takes its temporary file with it:
    /// with a directory standing at the `.fai` path the run is refused and
    /// nothing is left beside the FASTA but that directory.
    #[test]
    fn an_index_that_cannot_be_put_in_place_leaves_no_temporary_file() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, PLAIN).unwrap();
        std::fs::create_dir(dir.path().join("ref.fa.fai")).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai/kept"), "").unwrap();
        let error = index_if_missing(&fa).unwrap_err().to_string();
        assert!(
            error.starts_with(&format!(
                "failed to write FASTA index {}.fai: ",
                fa.display()
            )),
            "{error}"
        );
        assert_eq!(names_in(dir.path()), ["ref.fa", "ref.fa.fai"]);
        assert_eq!(names_in(&dir.path().join("ref.fa.fai")), ["kept"]);
    }

    /// A lock file left beside a FASTA and its indexes in a directory the run
    /// cannot write ends no run: nothing can be removed, and the indexes that
    /// are there are read, as Ensembl VEP reads them (`FastaSequence.pm:319`,
    /// whose `unlink` fails without a word).
    #[cfg(unix)]
    #[test]
    fn indexes_and_a_leftover_lock_in_an_unwritable_directory_are_read() {
        use crate::tests::bgzf_and_gzi;
        use std::os::unix::fs::PermissionsExt;
        for (name, bytes) in [
            ("ref.fa", PLAIN.to_vec()),
            ("ref.fa.gz", bgzf_and_gzi(PLAIN, 20).0),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let fa = dir.path().join(name);
            std::fs::write(&fa, bytes).unwrap();
            index_if_missing(&fa).unwrap();
            std::fs::write(sibling(&fa, "vep.lock"), "1\n").unwrap();
            let before = names_in(dir.path());
            std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o555)).unwrap();
            // Root writes into the directory anyway, and the state cannot be made.
            let unwritable = std::fs::write(dir.path().join("probe"), "").is_err();
            let result = open_indexed(&fa);
            let after = names_in(dir.path());
            std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o755)).unwrap();
            if !unwritable {
                continue;
            }
            let fasta = result.unwrap_or_else(|error| panic!("{name}: {error}"));
            assert_eq!(fasta.sequence("MT", 8, 10), Some(b"GCC".to_vec()), "{name}");
            assert_eq!(after, before, "{name}");
        }
    }

    /// A temporary file another run left beside the FASTA, under the name a
    /// run of this process id took, is neither in the way of the index nor
    /// removed.
    #[test]
    fn a_temporary_file_of_another_run_is_neither_in_the_way_nor_removed() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, PLAIN).unwrap();
        let left = format!("ref.fa.fai.{}.tmp", std::process::id());
        std::fs::write(dir.path().join(&left), "another run's").unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(
            std::fs::read_to_string(dir.path().join("ref.fa.fai")).unwrap(),
            PLAIN_FAI
        );
        assert_eq!(
            std::fs::read_to_string(dir.path().join(&left)).unwrap(),
            "another run's"
        );
        assert_eq!(
            names_in(dir.path()),
            ["ref.fa", "ref.fa.fai", left.as_str()]
        );
    }

    /// Runs that open one FASTA without its indexes at the same moment all
    /// read it, two or eight of them, plain or BGZF, with or without a lock
    /// file and a wrong index left by a run that did not finish; when the last
    /// has ended the indexes `samtools faidx` writes stand beside the file and
    /// no lock file or temporary file is left.
    #[test]
    fn runs_that_index_one_fasta_at_once_all_read_it() {
        use crate::tests::{bases, bgzf_and_gzi, fasta_and_fai};
        use std::sync::{Arc, Barrier};
        let (one, two) = (bases(31, 2_500), bases(32, 1_700));
        let (text, fai) = fasta_and_fai(&[("1", &one), ("2", &two)]);
        let (compressed, gzi) = bgzf_and_gzi(&text, 1_000);
        let mut gzi_bytes = (gzi.as_ref().len() as u64).to_le_bytes().to_vec();
        for (block, at) in gzi.as_ref() {
            gzi_bytes.extend_from_slice(&block.to_le_bytes());
            gzi_bytes.extend_from_slice(&at.to_le_bytes());
        }
        for (name, bytes) in [("ref.fa", &text), ("ref.fa.gz", &compressed)] {
            for leftover_lock in [false, true] {
                for width in [2, 8] {
                    for round in 0..10 {
                        let state =
                            format!("{name}, lock {leftover_lock}, {width} runs, round {round}");
                        let dir = tempfile::tempdir().unwrap();
                        let fa = dir.path().join(name);
                        std::fs::write(&fa, bytes).unwrap();
                        if leftover_lock {
                            std::fs::write(sibling(&fa, "fai"), "1\t5\t3\t60\t61\n").unwrap();
                            std::fs::write(sibling(&fa, "vep.lock"), "1\n").unwrap();
                        }
                        let barrier = Arc::new(Barrier::new(width));
                        let runs: Vec<_> = (0..width)
                            .map(|_| {
                                let (fa, barrier) = (fa.clone(), Arc::clone(&barrier));
                                std::thread::spawn(move || {
                                    barrier.wait();
                                    open_indexed(&fa).map(|fasta| fasta.sequence("2", 1, 1_700))
                                })
                            })
                            .collect();
                        for run in runs {
                            let read = run
                                .join()
                                .unwrap()
                                .unwrap_or_else(|error| panic!("{state}: {error}"));
                            assert_eq!(read, Some(two.to_ascii_uppercase()), "{state}");
                        }
                        assert_eq!(
                            std::fs::read_to_string(sibling(&fa, "fai")).unwrap(),
                            fai,
                            "{state}"
                        );
                        let mut left = vec![name.to_string(), format!("{name}.fai")];
                        if name.ends_with(".gz") {
                            assert_eq!(std::fs::read(sibling(&fa, "gzi")).unwrap(), gzi_bytes);
                            left.push(format!("{name}.gzi"));
                        }
                        assert_eq!(names_in(dir.path()), left, "{state}");
                    }
                }
            }
        }
    }

    /// A FASTA that is a link into a directory the run cannot write is
    /// indexed beside the link, where Ensembl VEP writes its lock and its
    /// index (`FastaSequence.pm:313`): nothing is written beside the file the
    /// link names.
    #[cfg(unix)]
    #[test]
    fn a_linked_fasta_is_indexed_beside_its_link() {
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        let (store, work) = (dir.path().join("store"), dir.path().join("work"));
        std::fs::create_dir(&store).unwrap();
        std::fs::create_dir(&work).unwrap();
        std::fs::write(store.join("ref.fa"), PLAIN).unwrap();
        std::os::unix::fs::symlink(store.join("ref.fa"), work.join("ref.fa")).unwrap();
        std::fs::set_permissions(&store, std::fs::Permissions::from_mode(0o555)).unwrap();
        let opened = open_indexed(&work.join("ref.fa"));
        let (in_store, in_work) = (names_in(&store), names_in(&work));
        std::fs::set_permissions(&store, std::fs::Permissions::from_mode(0o755)).unwrap();
        assert_eq!(opened.unwrap().sequence("MT", 8, 10), Some(b"GCC".to_vec()));
        assert_eq!(in_store, ["ref.fa"]);
        assert_eq!(in_work, ["ref.fa", "ref.fa.fai"]);
        assert_eq!(
            std::fs::read_to_string(work.join("ref.fa.fai")).unwrap(),
            PLAIN_FAI
        );
    }

    /// A BGZF FASTA with one of its two indexes gets the other, and the one
    /// it has is the file that was there, not written again: with its `.fai`
    /// alone the `.gzi` is written, with its `.gzi` alone the `.fai`.
    #[cfg(unix)]
    #[test]
    fn a_bgzf_fasta_with_one_index_gets_the_other_and_keeps_the_one_it_has() {
        use crate::tests::{bases, bgzf_and_gzi, fasta_and_fai};
        use std::os::unix::fs::MetadataExt;
        let (text, fai) = fasta_and_fai(&[("1", &bases(41, 2_500)), ("2", &bases(42, 700))]);
        let (compressed, gzi) = bgzf_and_gzi(&text, 1_000);
        for kept in ["fai", "gzi"] {
            let dir = tempfile::tempdir().unwrap();
            let fa = dir.path().join("ref.fa.gz");
            std::fs::write(&fa, &compressed).unwrap();
            match kept {
                "fai" => std::fs::write(sibling(&fa, "fai"), &fai).unwrap(),
                _ => bgzf::gzi::fs::write(sibling(&fa, "gzi"), &gzi).unwrap(),
            }
            let stamp = |path: &Path| {
                let file = std::fs::metadata(path).unwrap();
                (file.ino(), file.mtime(), file.mtime_nsec())
            };
            let before = stamp(&sibling(&fa, kept));
            index_if_missing(&fa).unwrap();
            assert_eq!(stamp(&sibling(&fa, kept)), before, "{kept} kept");
            assert_eq!(
                std::fs::read_to_string(sibling(&fa, "fai")).unwrap(),
                fai,
                "{kept} kept"
            );
            assert_eq!(
                bgzf::gzi::fs::read(sibling(&fa, "gzi")).unwrap().as_ref(),
                gzi.as_ref(),
                "{kept} kept"
            );
            assert_eq!(
                names_in(dir.path()),
                ["ref.fa.gz", "ref.fa.gz.fai", "ref.fa.gz.gzi"]
            );
            IndexedFasta::from_path(&fa).unwrap();
        }
    }

    /// A lock file left beside a BGZF FASTA removes both its indexes, and
    /// both are written again.
    #[test]
    fn a_lock_file_left_beside_a_bgzf_fasta_has_both_indexes_written_again() {
        use crate::tests::{bases, bgzf_and_gzi, fasta_and_fai};
        let (text, fai) = fasta_and_fai(&[("1", &bases(43, 2_500))]);
        let (compressed, gzi) = bgzf_and_gzi(&text, 1_000);
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa.gz");
        std::fs::write(&fa, compressed).unwrap();
        std::fs::write(sibling(&fa, "fai"), "1\t5\t3\t60\t61\n").unwrap();
        bgzf::gzi::fs::write(sibling(&fa, "gzi"), &bgzf::gzi::Index::default()).unwrap();
        std::fs::write(sibling(&fa, "vep.lock"), "1\n").unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(std::fs::read_to_string(sibling(&fa, "fai")).unwrap(), fai);
        assert_eq!(
            bgzf::gzi::fs::read(sibling(&fa, "gzi")).unwrap().as_ref(),
            gzi.as_ref()
        );
        assert_eq!(
            names_in(dir.path()),
            ["ref.fa.gz", "ref.fa.gz.fai", "ref.fa.gz.gzi"]
        );
    }

    /// A `<fasta>.index`, the index of Ensembl VEP's other FASTA reader, is
    /// kept beside the `.fai` written for the file, and a lock file found
    /// beside the FASTA removes it with the other indexes
    /// (`FastaSequence.pm:315-321`); it is not written again.
    #[test]
    fn an_index_of_another_reader_is_kept_and_removed_with_a_lock() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, PLAIN).unwrap();
        std::fs::write(dir.path().join("ref.fa.index"), "another reader's").unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(
            names_in(dir.path()),
            ["ref.fa", "ref.fa.fai", "ref.fa.index"]
        );
        assert_eq!(
            std::fs::read_to_string(dir.path().join("ref.fa.index")).unwrap(),
            "another reader's"
        );
        std::fs::write(dir.path().join("ref.fa.vep.lock"), "1\n").unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(names_in(dir.path()), ["ref.fa", "ref.fa.fai"]);
        assert_eq!(
            std::fs::read_to_string(dir.path().join("ref.fa.fai")).unwrap(),
            PLAIN_FAI
        );
    }

    /// A gzip file that is not BGZF is not indexed: nothing is written beside
    /// it, and the open refuses it as gzip, with or without an extra field in
    /// its header.
    #[test]
    fn a_gzip_file_that_is_not_bgzf_is_left_without_an_index() {
        use crate::tests::{gzip_with_extra, plain_gzip};
        for bytes in [
            plain_gzip(PLAIN),
            gzip_with_extra(PLAIN, &[b'X', b'Y', 2, 0, 0, 0]),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let fa = dir.path().join("ref.fa.gz");
            std::fs::write(&fa, bytes).unwrap();
            index_if_missing(&fa).unwrap();
            assert_eq!(names_in(dir.path()), ["ref.fa.gz"]);
            let error = open_indexed(&fa).unwrap_err();
            assert!(matches!(error, FastaError::NotBgzf(_)), "{error:?}");
            assert_eq!(names_in(dir.path()), ["ref.fa.gz"]);
        }
    }

    /// A BGZF file whose first block holds no data gets the `.gzi` `samtools
    /// faidx` writes for it, which leaves out the first block that holds data,
    /// and reads through it.
    #[test]
    fn a_first_bgzf_block_that_holds_no_data_is_stepped_over() {
        use crate::tests::{bases, bgzf_and_gzi, fasta_and_fai};
        let dir = tempfile::tempdir().unwrap();
        let one = bases(44, 2_500);
        let (text, fai) = fasta_and_fai(&[("1", &one)]);
        let (compressed, gzi) = bgzf_and_gzi(&text, 1_000);
        let empty = &compressed[compressed.len() - BGZF_EOF_LEN..];
        let mut expected = (gzi.as_ref().len() as u64).to_le_bytes().to_vec();
        for (block, at) in gzi.as_ref() {
            expected.extend_from_slice(&(block + BGZF_EOF_LEN as u64).to_le_bytes());
            expected.extend_from_slice(&at.to_le_bytes());
        }
        let fa = dir.path().join("ref.fa.gz");
        std::fs::write(&fa, [empty, compressed.as_slice()].concat()).unwrap();
        index_if_missing(&fa).unwrap();
        assert_eq!(std::fs::read(sibling(&fa, "gzi")).unwrap(), expected);
        assert_eq!(std::fs::read_to_string(sibling(&fa, "fai")).unwrap(), fai);
        let fasta = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(
            fasta.sequence("1", 1, 2_500),
            Some(one.to_ascii_uppercase())
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

    /// A name holding a NUL byte is indexed up to that byte, as `samtools
    /// faidx` indexes it: a name that opens with the byte is the empty name,
    /// and two names equal up to the byte are one name given twice.
    #[test]
    fn a_name_ends_before_a_nul_byte() {
        let (full, other) = (made_up(8), made_up(4));
        assert_eq!(
            fai_of(&format!(">a\0b desc\n{full}\n")).unwrap(),
            "a\t8\t10\t8\t9\n"
        );
        assert_eq!(
            fai_of(&format!(">\0ab\n{full}\n")).unwrap(),
            "\t8\t5\t8\t9\n"
        );
        assert_eq!(
            fai_of(&format!(">a\0b\n{full}\n>a\0c\n{other}\n>b\n{other}\n")).unwrap(),
            "a\t8\t5\t8\t9\nb\t4\t27\t4\t5\n"
        );
    }
}
