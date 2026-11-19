// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The genome FASTA a cache directory carries: a BGZF-compressed `.fa.gz` with
//! the `.fai` and `.gzi` indexes `vep` reads it through, written here from a
//! gzip or BGZF source without `bgzip` or `samtools`.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use tracing::info;

/// The path of an index beside `fasta`: `<fasta>.<suffix>`.
pub fn index_path(fasta: &Path, suffix: &str) -> PathBuf {
    let mut name = fasta.as_os_str().to_owned();
    name.push(".");
    name.push(suffix);
    PathBuf::from(name)
}

/// Places the gzip-compressed FASTA at `source` as the BGZF file `target` with
/// `target.fai` and `target.gzi` beside it, replacing anything at those paths. A
/// source that is already BGZF is copied as it is; a plain gzip source is
/// inflated and written again in BGZF blocks. The result is opened once, so an
/// index that does not describe the file is refused here, not at annotation
/// time.
pub fn install_fasta(source: &Path, target: &Path) -> Result<()> {
    if is_bgzf(source)? {
        info!("Copying the BGZF FASTA to {}", target.display());
        std::fs::copy(source, target).with_context(|| {
            format!(
                "Failed to copy {} to {}",
                source.display(),
                target.display()
            )
        })?;
    } else {
        info!("Recompressing {} to BGZF", source.display());
        recompress_to_bgzf(source, target)?;
    }
    info!("Indexing {}", target.display());
    write_fai(target)?;
    write_gzi(target)?;
    vep_fasta::IndexedFasta::from_path(target)
        .with_context(|| format!("The indexed FASTA {} does not open", target.display()))?;
    Ok(())
}

/// Whether the gzip file at `path` is BGZF: its first member carries the `BC`
/// extra subfield that states the block size.
pub fn is_bgzf(path: &Path) -> Result<bool> {
    let mut file =
        File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
    let mut header = [0u8; 18];
    match file.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(false),
        Err(e) => return Err(e).with_context(|| format!("Failed to read {}", path.display())),
    }
    Ok(block_size(&header).is_some())
}

/// The smallest BGZF block: an 18-byte header, the two bytes of an empty
/// deflate stream and the 8-byte trailer, the length of the EOF marker block.
const MIN_BLOCK_SIZE: usize = 28;

/// The gzip member header fields of BGZF: identification `1f 8b`, deflate, the
/// FEXTRA flag and a six-byte extra field whose `BC` subfield holds the block's
/// total size minus one; `None` for a header that is not BGZF.
fn block_size(header: &[u8; 18]) -> Option<usize> {
    let gzip_deflate = header[0] == 0x1f && header[1] == 0x8b && header[2] == 8;
    let fextra = header[3] & 0x04 != 0;
    let xlen = u16::from_le_bytes([header[10], header[11]]);
    let bc_subfield = header[12] == b'B'
        && header[13] == b'C'
        && u16::from_le_bytes([header[14], header[15]]) == 2;
    if !(gzip_deflate && fextra && xlen == 6 && bc_subfield) {
        return None;
    }
    Some(usize::from(u16::from_le_bytes([header[16], header[17]])) + 1)
}

/// Inflates the gzip file at `source` (every member) and writes it to `target`
/// as BGZF.
fn recompress_to_bgzf(source: &Path, target: &Path) -> Result<()> {
    let input =
        File::open(source).with_context(|| format!("Failed to open {}", source.display()))?;
    let mut decoder = flate2::read::MultiGzDecoder::new(BufReader::new(input));
    let output =
        File::create(target).with_context(|| format!("Failed to create {}", target.display()))?;
    let mut writer = noodles::bgzf::io::Writer::new(BufWriter::new(output));
    std::io::copy(&mut decoder, &mut writer)
        .with_context(|| format!("Failed to recompress {}", source.display()))?;
    writer
        .finish()
        .with_context(|| format!("Failed to finish {}", target.display()))?
        .flush()
        .with_context(|| format!("Failed to write {}", target.display()))?;
    Ok(())
}

/// Writes `<fasta>.fai`, the `samtools faidx` index of the BGZF FASTA at
/// `fasta`, each record's offset counted in the inflated stream.
fn write_fai(fasta: &Path) -> Result<()> {
    let file = File::open(fasta).with_context(|| format!("Failed to open {}", fasta.display()))?;
    let reader = noodles::bgzf::io::Reader::new(BufReader::new(file));
    let mut indexer = noodles::fasta::io::Indexer::new(reader);
    let mut records = Vec::new();
    while let Some(record) = indexer
        .index_record()
        .with_context(|| format!("Failed to index {}", fasta.display()))?
    {
        records.push(record);
    }
    let index = noodles::fasta::fai::Index::from(records);
    let fai = index_path(fasta, "fai");
    let mut writer = noodles::fasta::fai::io::Writer::new(BufWriter::new(
        File::create(&fai).with_context(|| format!("Failed to create {}", fai.display()))?,
    ));
    writer
        .write_index(&index)
        .with_context(|| format!("Failed to write {}", fai.display()))?;
    writer.into_inner().flush()?;
    Ok(())
}

/// Writes `<fasta>.gzi`, the `bgzip` index of the BGZF file at `fasta`: the
/// compressed and inflated offsets of every non-empty block after the first,
/// read from the block headers and trailers without inflating the data.
fn write_gzi(fasta: &Path) -> Result<()> {
    let mut file =
        File::open(fasta).with_context(|| format!("Failed to open {}", fasta.display()))?;
    let length = file.metadata()?.len();
    let mut entries = Vec::new();
    let (mut compressed, mut inflated) = (0u64, 0u64);
    let mut header = [0u8; 18];
    while compressed < length {
        file.seek(SeekFrom::Start(compressed))?;
        file.read_exact(&mut header)
            .with_context(|| format!("{}: truncated block header", fasta.display()))?;
        let Some(size) = block_size(&header) else {
            bail!(
                "{} is gzip but not BGZF at offset {compressed}",
                fasta.display()
            );
        };
        if size < MIN_BLOCK_SIZE {
            bail!(
                "{} states a {size}-byte block at offset {compressed}, shorter than the {MIN_BLOCK_SIZE}-byte minimum of BGZF",
                fasta.display()
            );
        }
        let size = size as u64;
        let mut trailer = [0u8; 4];
        file.seek(SeekFrom::Start(compressed + size - 4))?;
        file.read_exact(&mut trailer).with_context(|| {
            format!(
                "{}: truncated block at offset {compressed}",
                fasta.display()
            )
        })?;
        let data_len = u64::from(u32::from_le_bytes(trailer));
        if data_len > 0 && compressed > 0 {
            entries.push((compressed, inflated));
        }
        compressed += size;
        inflated += data_len;
    }
    let gzi = index_path(fasta, "gzi");
    noodles::bgzf::gzi::fs::write(&gzi, &noodles::bgzf::gzi::Index::from(entries))
        .with_context(|| format!("Failed to write {}", gzi.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn committed_reference() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh37-hgvs/reference.fa.gz")
    }

    /// A `bgzip`-written source is kept as it is, and the `.fai` and `.gzi`
    /// written for it are byte for byte the `samtools faidx` indexes committed
    /// beside the reference.
    #[test]
    fn bgzf_source_indexes_match_the_committed_reference() {
        let reference = committed_reference();
        assert!(is_bgzf(&reference).unwrap());
        let dir = tempfile::tempdir().unwrap();
        let target = dir
            .path()
            .join("Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz");
        install_fasta(&reference, &target).unwrap();
        assert_eq!(
            std::fs::read(&target).unwrap(),
            std::fs::read(&reference).unwrap()
        );
        for suffix in ["fai", "gzi"] {
            assert_eq!(
                std::fs::read(index_path(&target, suffix)).unwrap(),
                std::fs::read(index_path(&reference, suffix)).unwrap(),
                "{suffix}"
            );
        }
    }

    /// A plain gzip source is written again as BGZF, indexed, and reads back
    /// base for base through the indexes; the `.fai` is the plain file's.
    #[test]
    fn plain_gzip_source_is_recompressed_and_indexed() {
        let dir = tempfile::tempdir().unwrap();
        // Three sequences over 200 KiB, so the BGZF output spans several blocks
        // and the `.gzi` carries entries.
        let mut plain = Vec::new();
        for (name, base, lines) in [("1", b'A', 2000), ("2", b'C', 1500), ("MT", b'G', 3)] {
            plain.extend_from_slice(format!(">{name} dna:chromosome\n").as_bytes());
            for _ in 0..lines {
                plain.extend(std::iter::repeat_n(base, 60));
                plain.push(b'\n');
            }
            plain.extend_from_slice(b"TTTTT\n");
        }
        let source = dir.path().join("source.fa.gz");
        let mut encoder = flate2::write::GzEncoder::new(
            File::create(&source).unwrap(),
            flate2::Compression::default(),
        );
        encoder.write_all(&plain).unwrap();
        encoder.finish().unwrap();
        assert!(!is_bgzf(&source).unwrap());

        let target = dir.path().join("genome.fa.gz");
        install_fasta(&source, &target).unwrap();
        assert!(is_bgzf(&target).unwrap());
        let fai = std::fs::read_to_string(index_path(&target, "fai")).unwrap();
        assert_eq!(
            fai,
            "1\t120005\t18\t60\t61\n2\t90005\t122042\t60\t61\nMT\t185\t213567\t60\t61\n"
        );
        let gzi = noodles::bgzf::gzi::fs::read(index_path(&target, "gzi")).unwrap();
        assert!(gzi.as_ref().len() >= 3, "{gzi:?}");

        let fasta = vep_fasta::IndexedFasta::from_path(&target).unwrap();
        assert_eq!(fasta.base("1", 1), Some(b'A'));
        assert_eq!(fasta.base("1", 120_000), Some(b'A'));
        assert_eq!(fasta.base("1", 120_001), Some(b'T'));
        assert_eq!(fasta.base("2", 90_005), Some(b'T'));
        assert_eq!(fasta.base("MT", 180), Some(b'G'));
        assert_eq!(fasta.base("MT", 186), None);
        let mut inflated = Vec::new();
        noodles::bgzf::io::Reader::new(File::open(&target).unwrap())
            .read_to_end(&mut inflated)
            .unwrap();
        assert_eq!(inflated, plain);
    }

    /// A gzip header without the `BC` subfield, a short file and a non-gzip file
    /// are not BGZF.
    #[test]
    fn bgzf_detection_reads_the_extra_subfield() {
        let dir = tempfile::tempdir().unwrap();
        let short = dir.path().join("short");
        std::fs::write(&short, b"\x1f\x8b\x08").unwrap();
        assert!(!is_bgzf(&short).unwrap());
        let text = dir.path().join("text.fa");
        std::fs::write(&text, b">1\nACGT\n>2\nACGTACGTACGTACGT\n").unwrap();
        assert!(!is_bgzf(&text).unwrap());
        let mut header = [0u8; 18];
        header[..4].copy_from_slice(&[0x1f, 0x8b, 8, 0x04]);
        header[10..12].copy_from_slice(&6u16.to_le_bytes());
        header[12..16].copy_from_slice(&[b'B', b'C', 2, 0]);
        header[16..18].copy_from_slice(&27u16.to_le_bytes());
        assert_eq!(block_size(&header), Some(28));
        header[12] = b'X';
        assert_eq!(block_size(&header), None);
    }

    /// A BGZF header whose `BC` subfield states a block shorter than the
    /// 28-byte minimum (one byte short, and short enough that the trailer would
    /// sit before the block's own start) is refused with the file and offset
    /// named, by the index writer and by the installer.
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
            assert!(is_bgzf(&crafted).unwrap());

            let error = write_gzi(&crafted).unwrap_err().to_string();
            assert!(
                error.contains(&format!("{}-byte block at offset 0", bsize + 1))
                    && error.contains("28-byte minimum"),
                "{error}"
            );
            assert!(
                install_fasta(&crafted, &dir.path().join(format!("target{bsize}.fa.gz"))).is_err()
            );
        }
    }
}
