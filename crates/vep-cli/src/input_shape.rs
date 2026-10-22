// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The shape of a bgzipped, indexed input, read from four sampled windows: how
//! many bytes a record inflates to and how many records the file holds, for the
//! `auto` decompression thread count and `--resolved_config`.
//!
//! A BGZF file is recognised by its first block header, never by its name, with
//! the reader's own detection (`runner::bgzf_block_size`). The windows
//! open at the blocks nearest one, three, five and seven eighths of the file,
//! each resolved to a block start through the `.tbi` or `.csi` index's virtual
//! offsets, and each reads consecutive blocks until it holds
//! `WINDOW_MIN_LINES` complete data lines or `WINDOW_MAX_BLOCKS` blocks, so
//! a file whose records fill a block with a handful of lines is still measured
//! over tens of them.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use noodles::bgzf;
use noodles::csi::binning_index::index::reference_sequence::index::Index as ReferenceIndex;
use noodles::csi::binning_index::Index;

use crate::runner::{bgzf_block_size, GZIP_HEADER_BYTES};

/// Eighths of the file at which a window opens.
const SAMPLE_EIGHTHS: [u64; 4] = [1, 3, 5, 7];

/// Complete data lines a window reads at least, over consecutive blocks.
const WINDOW_MIN_LINES: u64 = 64;

/// Blocks a window reads at most, whatever it holds.
const WINDOW_MAX_BLOCKS: usize = 16;

/// What is known of an input before it is read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputShape {
    /// Standard input, a path that is not a regular file, a file that is not
    /// BGZF, a BGZF file with no `.tbi` or `.csi` beside it, or one whose
    /// windows cannot be read or hold no complete data line.
    Unknown,
    /// A bgzipped, indexed file.
    Indexed {
        /// Inflated bytes per complete data line over the windows, newline
        /// included.
        inflated_bytes_per_record: u64,
        /// The file size over the windows' compressed bytes per data line end.
        /// An estimate: the windows stand for the whole file, and a file of a
        /// few blocks carries its header's and end-of-file block's bytes as if
        /// they held records.
        records_estimate: u64,
        file_bytes: u64,
    },
}

impl InputShape {
    /// The shape of the input at `path`, `Unknown` for `STDIN` and for anything
    /// the module documentation lists.
    pub fn read(path: &str) -> InputShape {
        if path == "STDIN" {
            return InputShape::Unknown;
        }
        read_indexed(Path::new(path)).unwrap_or(InputShape::Unknown)
    }

    /// The shape as `--resolved_config` names it.
    pub fn name(&self) -> &'static str {
        match self {
            InputShape::Unknown => "unknown",
            InputShape::Indexed { .. } => "indexed",
        }
    }
}

/// What the windows hold between them.
#[derive(Default)]
struct Sample {
    /// Complete data lines and their bytes, newlines included.
    lines: u64,
    line_bytes: u64,
    /// Data lines ending inside the windows, the ones a window start cuts
    /// included.
    line_ends: u64,
    compressed_bytes: u64,
}

fn read_indexed(path: &Path) -> Option<InputShape> {
    // A named pipe is opened once, by the run: opening it here would take its
    // writer's bytes, and the open blocks until a writer arrives.
    let file_bytes = std::fs::metadata(path).ok().filter(|m| m.is_file())?.len();
    let mut file = File::open(path).ok()?;
    bgzf_block_size(&read_block_header(&mut file, 0)?)?;
    let offsets: Vec<u64> = block_offsets(path)?
        .into_iter()
        .filter(|&offset| offset < file_bytes)
        .collect();
    let mut sample = Sample::default();
    for eighth in SAMPLE_EIGHTHS {
        let target = file_bytes / 8 * eighth;
        // The first block at or past the target, else the last one before it:
        // a file of fewer blocks than windows still has a block to read, and a
        // chunk ending at a block boundary names the next block, which at the
        // end of the file is the empty end-of-file block.
        let at_or_past = offsets.iter().copied().find(|&offset| offset >= target);
        let before = offsets.iter().rev().copied().find(|&o| o < target);
        [at_or_past, before]
            .into_iter()
            .flatten()
            .find_map(|start| read_window(&mut file, start, &mut sample))?;
    }
    if sample.lines == 0 {
        return None;
    }
    let records_estimate = u64::try_from(
        u128::from(file_bytes) * u128::from(sample.line_ends) / u128::from(sample.compressed_bytes),
    )
    .ok()?;
    Some(InputShape::Indexed {
        inflated_bytes_per_record: sample.line_bytes / sample.lines,
        records_estimate,
        file_bytes,
    })
}

/// Reads the window opening at block `start` into `sample`: consecutive blocks
/// until the window holds [`WINDOW_MIN_LINES`] complete data lines,
/// [`WINDOW_MAX_BLOCKS`] blocks, a block that cannot be read or the empty
/// end-of-file block. `None` when the first block cannot be read.
fn read_window(file: &mut File, start: u64, sample: &mut Sample) -> Option<()> {
    let at_file_start = start == 0;
    let mut text = Vec::new();
    let mut compressed_bytes = 0u64;
    let (mut lines, mut line_bytes) = (0u64, 0u64);
    for _ in 0..WINDOW_MAX_BLOCKS {
        let Some((size, inflated)) = read_block(file, start + compressed_bytes) else {
            break;
        };
        if inflated.is_empty() {
            break;
        }
        text.extend_from_slice(&inflated);
        compressed_bytes += size;
        (lines, line_bytes) = complete_data_lines(&text, at_file_start);
        if lines >= WINDOW_MIN_LINES {
            break;
        }
    }
    if compressed_bytes == 0 {
        return None;
    }
    // A window opening past the file start cuts a line, whose end is a data
    // line's end as much as a complete line's; the trailing cut line has no end.
    let cut_line_end = u64::from(!at_file_start && text.contains(&b'\n'));
    sample.lines += lines;
    sample.line_bytes += line_bytes;
    sample.line_ends += lines + cut_line_end;
    sample.compressed_bytes += compressed_bytes;
    Some(())
}

/// The gzip header at `offset`, the fixed fields with the extra field.
fn read_block_header(file: &mut File, offset: u64) -> Option<Vec<u8>> {
    file.seek(SeekFrom::Start(offset)).ok()?;
    let mut header = vec![0u8; GZIP_HEADER_BYTES];
    file.read_exact(&mut header).ok()?;
    let extra_len = u16::from_le_bytes([header[10], header[11]]) as usize;
    header.resize(GZIP_HEADER_BYTES + extra_len, 0);
    file.read_exact(&mut header[GZIP_HEADER_BYTES..]).ok()?;
    Some(header)
}

/// The block at `offset`: its compressed size and its inflated bytes.
fn read_block(file: &mut File, offset: u64) -> Option<(u64, Vec<u8>)> {
    let header = read_block_header(file, offset)?;
    let size = bgzf_block_size(&header)?;
    let header_len = header.len();
    if size < header_len {
        return None;
    }
    let mut block = header;
    block.resize(size, 0);
    file.read_exact(&mut block[header_len..]).ok()?;
    let mut inflated = Vec::new();
    bgzf::io::Reader::new(&block[..])
        .read_to_end(&mut inflated)
        .ok()?;
    Some((size as u64, inflated))
}

/// Every block start the `.tbi` or `.csi` index beside `path` names, sorted
/// and deduplicated; `None` without a readable index.
fn block_offsets(path: &Path) -> Option<Vec<u64>> {
    let with_suffix = |suffix: &str| {
        let mut name = path.as_os_str().to_owned();
        name.push(suffix);
        PathBuf::from(name)
    };
    let tbi = with_suffix(".tbi");
    let csi = with_suffix(".csi");
    if tbi.is_file() {
        vep_builtins::tabix::read_tabix_index(&tbi)
            .ok()
            .map(|index| index_block_offsets(&index))
    } else if csi.is_file() {
        vep_builtins::tabix::read_csi_index(&csi)
            .ok()
            .map(|index| index_block_offsets(&index))
    } else {
        None
    }
}

/// The compressed offsets of every chunk start and end in `index`. A virtual
/// position's compressed part is a block start.
fn index_block_offsets<I: ReferenceIndex>(index: &Index<I>) -> Vec<u64> {
    let mut offsets: Vec<u64> = index
        .reference_sequences()
        .iter()
        .flat_map(|sequence| sequence.bins().values())
        .flat_map(|bin| bin.chunks())
        .flat_map(|chunk| [chunk.start().compressed(), chunk.end().compressed()])
        .collect();
    offsets.sort_unstable();
    offsets.dedup();
    offsets
}

/// The complete data lines of a window's inflated text: the newline-ended lines
/// that do not open with `#`, as a count and as bytes with their newlines. The
/// text after the last newline is the line the window's end cuts; the text
/// before the first newline is the line its start cuts, unless the window opens
/// the file (`at_file_start`), where that line is whole.
fn complete_data_lines(inflated: &[u8], at_file_start: bool) -> (u64, u64) {
    let Some(last_newline) = inflated.iter().rposition(|&b| b == b'\n') else {
        return (0, 0);
    };
    let mut lines = inflated[..=last_newline].split_inclusive(|&b| b == b'\n');
    if !at_file_start {
        lines.next();
    }
    lines
        .filter(|line| line.first() != Some(&b'#'))
        .fold((0, 0), |(count, bytes), line| {
            (count + 1, bytes + line.len() as u64)
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use noodles::core::Position;
    use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;

    fn golden_corpus(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116")
            .join(name)
            .join("variants.vcf")
    }

    /// A bgzipped copy of `vcf` at `out`, with its tabix index beside it and the
    /// same chunks as a CSI index at `<out>.csi.aside` (a file the reader does
    /// not look for, so a test places it as `.csi` itself); one block is closed
    /// every `records_per_block` records so a small file still spans several
    /// blocks, and `padding` bytes of a pseudo-random INFO field widen every
    /// record. Returns the data lines written, newlines included.
    fn bgzip_with_index(
        vcf: &Path,
        out: &Path,
        records_per_block: usize,
        padding: usize,
    ) -> Vec<String> {
        let text = std::fs::read_to_string(vcf).unwrap();
        let mut writer = bgzf::io::Writer::new(File::create(out).unwrap());
        let mut indexer = noodles::tabix::index::Indexer::default();
        indexer.set_header(HeaderBuilder::vcf().build());
        let mut csi_indexer = noodles::csi::binning_index::Indexer::<BinnedIndex>::default();
        let mut chromosomes: Vec<String> = Vec::new();
        let mut records = 0usize;
        let mut data_lines = Vec::new();
        // A linear congruential generator: letters that deflate like annotation
        // text rather than like a run of one byte.
        let mut state = 12_345u32;
        let mut pad = move || {
            (0..padding)
                .map(|_| {
                    state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                    b'A' + ((state >> 16) % 26) as u8
                })
                .map(char::from)
                .collect::<String>()
        };
        for line in text.lines() {
            if line.starts_with('#') {
                writeln!(writer, "{line}").unwrap();
                continue;
            }
            let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
            if padding > 0 {
                let pad = format!("PAD={}", pad());
                fields[7] = if fields[7] == "." {
                    pad
                } else {
                    format!("{};{pad}", fields[7])
                };
            }
            let line = fields.join("\t");
            let start = writer.virtual_position();
            writeln!(writer, "{line}").unwrap();
            data_lines.push(format!("{line}\n"));
            let pos: usize = fields[1].parse().unwrap();
            let end = pos + fields[3].len() - 1;
            let chunk = Chunk::new(start, writer.virtual_position());
            let (start_position, end_position) = (
                Position::try_from(pos).unwrap(),
                Position::try_from(end).unwrap(),
            );
            indexer
                .add_record(&fields[0], start_position, end_position, chunk)
                .unwrap();
            let chromosome_index = chromosomes
                .iter()
                .position(|c| *c == fields[0])
                .unwrap_or_else(|| {
                    chromosomes.push(fields[0].clone());
                    chromosomes.len() - 1
                });
            csi_indexer
                .add_record(
                    Some((chromosome_index, start_position, end_position, true)),
                    chunk,
                )
                .unwrap();
            records += 1;
            if records.is_multiple_of(records_per_block) {
                writer.flush().unwrap();
            }
        }
        writer.finish().unwrap();
        noodles::tabix::fs::write(format!("{}.tbi", out.display()), &indexer.build()).unwrap();
        noodles::csi::fs::write(
            format!("{}.csi.aside", out.display()),
            &csi_indexer.build(chromosomes.len()),
        )
        .unwrap();
        data_lines
    }

    /// Mean bytes of the data lines, newline included.
    fn mean_line_bytes<L: AsRef<[u8]>>(lines: &[L]) -> u64 {
        lines.iter().map(|l| l.as_ref().len() as u64).sum::<u64>() / lines.len() as u64
    }

    /// The golden corpora's plain VCFs, a gzip file without the `BC` subfield
    /// and a BGZF file without an index are all `Unknown`, as is standard input.
    #[test]
    fn plain_gzip_and_unindexed_inputs_are_unknown() {
        assert_eq!(InputShape::read("STDIN"), InputShape::Unknown);
        for corpus in ["GRCh37", "GRCh38", "GRCh38-run-control"] {
            let vcf = golden_corpus(corpus);
            assert_eq!(
                InputShape::read(vcf.to_str().unwrap()),
                InputShape::Unknown,
                "{corpus}"
            );
        }
        let tmp = tempfile::tempdir().unwrap();
        let text = std::fs::read(golden_corpus("GRCh38-run-control")).unwrap();

        let gzip = tmp.path().join("plain.vcf.gz");
        let mut encoder = flate2::write::GzEncoder::new(
            File::create(&gzip).unwrap(),
            flate2::Compression::default(),
        );
        encoder.write_all(&text).unwrap();
        encoder.finish().unwrap();
        assert_eq!(
            InputShape::read(gzip.to_str().unwrap()),
            InputShape::Unknown
        );

        let bgzf = tmp.path().join("unindexed.vcf.gz");
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf).unwrap());
        writer.write_all(&text).unwrap();
        writer.finish().unwrap();
        assert_eq!(
            InputShape::read(bgzf.to_str().unwrap()),
            InputShape::Unknown
        );

        assert_eq!(
            InputShape::read(tmp.path().join("absent.vcf.gz").to_str().unwrap()),
            InputShape::Unknown
        );
    }

    /// A bgzipped, tabix-indexed corpus reads as `Indexed`, with the inflated
    /// bytes per record within a fifth of the corpus's mean line length and the
    /// record estimate within a quarter of the record count, whether a block
    /// holds forty lines or four. The quarter is the corpora's own spread: the
    /// blocks of the GRCh37 corpus differ four-fold in compressed bytes per
    /// record (its structural-variant records against its point variants), and
    /// the small corpus's header is a sixth of its file, so four windows of
    /// either carry a share of that.
    #[test]
    fn indexed_corpus_reads_its_record_shape() {
        let tmp = tempfile::tempdir().unwrap();
        for (corpus, records_per_block) in [("GRCh37", 40), ("GRCh38-release116", 4)] {
            let vcf = golden_corpus(corpus);
            let out = tmp.path().join(format!("{corpus}.vcf.gz"));
            let lines = bgzip_with_index(&vcf, &out, records_per_block, 0);
            let shape = InputShape::read(out.to_str().unwrap());
            let InputShape::Indexed {
                inflated_bytes_per_record,
                records_estimate,
                file_bytes,
            } = shape
            else {
                panic!("{corpus}: {shape:?}");
            };
            let true_mean = mean_line_bytes(&lines);
            assert!(
                (true_mean * 4 / 5..=true_mean * 6 / 5).contains(&inflated_bytes_per_record),
                "{corpus}: {inflated_bytes_per_record} bytes per record against a mean of {true_mean}"
            );
            let records = lines.len() as u64;
            assert!(
                (records * 3 / 4..=records * 5 / 4).contains(&records_estimate),
                "{corpus}: {records_estimate} records estimated against {records}"
            );
            assert_eq!(file_bytes, std::fs::metadata(&out).unwrap().len());
            assert_eq!(shape.name(), "indexed");
        }
    }

    /// A `.csi` beside the input, with no `.tbi`, resolves the same block starts
    /// and so reads the same shape as the tabix index does.
    #[test]
    fn a_csi_index_reads_as_the_tabix_index_does() {
        let tmp = tempfile::tempdir().unwrap();
        let tbi_indexed = tmp.path().join("tbi.vcf.gz");
        bgzip_with_index(&golden_corpus("GRCh37"), &tbi_indexed, 40, 0);
        let csi_indexed = tmp.path().join("csi.vcf.gz");
        std::fs::copy(&tbi_indexed, &csi_indexed).unwrap();
        std::fs::copy(
            format!("{}.csi.aside", tbi_indexed.display()),
            format!("{}.csi", csi_indexed.display()),
        )
        .unwrap();
        let from_tbi = InputShape::read(tbi_indexed.to_str().unwrap());
        let from_csi = InputShape::read(csi_indexed.to_str().unwrap());
        assert!(
            matches!(from_tbi, InputShape::Indexed { .. }),
            "{from_tbi:?}"
        );
        assert_eq!(from_csi, from_tbi);
    }

    /// Records wide enough that a block holds four or five of them, the shape of
    /// an input whose records carry a long INFO field, are still counted within
    /// a tenth: a window reads consecutive blocks until it holds tens of lines,
    /// and every line end inside it is a record, the one its start cuts included.
    /// Counting complete lines alone over single blocks would lose one line in
    /// five here.
    #[test]
    fn wide_records_are_estimated_within_a_tenth() {
        let tmp = tempfile::tempdir().unwrap();
        let vcf = golden_corpus("GRCh38-input-compressed");
        let out = tmp.path().join("wide.vcf.gz");
        let lines = bgzip_with_index(&vcf, &out, usize::MAX, 14_000);
        let shape = InputShape::read(out.to_str().unwrap());
        let InputShape::Indexed {
            inflated_bytes_per_record,
            records_estimate,
            file_bytes,
        } = shape
        else {
            panic!("{shape:?}");
        };
        let records = lines.len() as u64;
        let true_mean = mean_line_bytes(&lines);
        assert!(
            file_bytes > 4 * 65_536,
            "{file_bytes} bytes: the input spans fewer blocks than windows"
        );
        assert!(
            (records * 9 / 10..=records * 11 / 10).contains(&records_estimate),
            "{records_estimate} records estimated against {records}"
        );
        assert!(
            (true_mean * 9 / 10..=true_mean * 11 / 10).contains(&inflated_bytes_per_record),
            "{inflated_bytes_per_record} bytes per record against a mean of {true_mean}"
        );
        assert!(inflated_bytes_per_record >= crate::config::INFLATE_BOUND_BYTES_PER_RECORD);
    }

    /// The committed compressed-input corpus, one BGZF block with its index
    /// beside it, reads its own shape: the bytes per record are the mean of its
    /// data lines exactly, since the one block holds every line whole, and the
    /// record estimate is the count plus the end-of-file block's share.
    #[test]
    fn committed_indexed_corpus_reads_its_shape() {
        let gz = golden_corpus("GRCh38-input-compressed").with_extension("vcf.gz");
        let shape = InputShape::read(gz.to_str().unwrap());
        let InputShape::Indexed {
            inflated_bytes_per_record,
            records_estimate,
            file_bytes,
        } = shape
        else {
            panic!("{shape:?}");
        };
        let lines: Vec<String> = std::fs::read_to_string(golden_corpus("GRCh38-input-compressed"))
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| format!("{l}\n"))
            .collect();
        assert_eq!(inflated_bytes_per_record, mean_line_bytes(&lines));
        let records = lines.len() as u64;
        assert!(
            (records..=records + 2).contains(&records_estimate),
            "{records_estimate} records estimated against {records}"
        );
        assert_eq!(file_bytes, std::fs::metadata(&gz).unwrap().len());
    }

    /// A named pipe is `Unknown` without being opened: opening a pipe with no
    /// writer blocks, and a pipe with one belongs to the run that reads it. The
    /// pipe has an index beside it so the name alone would read as indexed.
    #[cfg(unix)]
    #[test]
    fn a_named_pipe_is_unknown_and_never_opened() {
        let tmp = tempfile::tempdir().unwrap();
        let fifo = tmp.path().join("input.vcf.gz");
        let status = std::process::Command::new("mkfifo")
            .arg(&fifo)
            .status()
            .expect("mkfifo runs");
        assert!(status.success());
        std::fs::write(format!("{}.tbi", fifo.display()), b"").unwrap();
        assert_eq!(
            InputShape::read(fifo.to_str().unwrap()),
            InputShape::Unknown
        );
    }

    /// A damaged indexed input never panics the reader: the file emptied, its
    /// header's magic overwritten, its index replaced by noise and its index
    /// emptied read as `Unknown`; a file cut short reads the blocks that survive,
    /// so it is `Unknown` when a window's first block is cut and otherwise a
    /// shape whose `file_bytes` is the cut length.
    #[test]
    fn damaged_indexed_inputs_are_unknown() {
        let tmp = tempfile::tempdir().unwrap();
        let vcf = golden_corpus("GRCh37");
        let good = tmp.path().join("good.vcf.gz");
        bgzip_with_index(&vcf, &good, 40, 0);
        assert!(matches!(
            InputShape::read(good.to_str().unwrap()),
            InputShape::Indexed { .. }
        ));
        let bytes = std::fs::read(&good).unwrap();
        let index = std::fs::read(format!("{}.tbi", good.display())).unwrap();

        let damaged = |name: &str, data: &[u8], index_data: &[u8]| {
            let path = tmp.path().join(name);
            std::fs::write(&path, data).unwrap();
            std::fs::write(format!("{}.tbi", path.display()), index_data).unwrap();
            InputShape::read(path.to_str().unwrap())
        };
        let cut = bytes.len() / 2;
        match damaged("truncated.vcf.gz", &bytes[..cut], &index) {
            InputShape::Unknown => {}
            InputShape::Indexed { file_bytes, .. } => assert_eq!(file_bytes, cut as u64),
        }
        assert_eq!(damaged("empty.vcf.gz", &[], &index), InputShape::Unknown);
        let mut wrong_magic = bytes.clone();
        wrong_magic[0] = b'X';
        assert_eq!(
            damaged("wrong_magic.vcf.gz", &wrong_magic, &index),
            InputShape::Unknown
        );
        let noise: Vec<u8> = (0..200u32).map(|i| (i * 37 % 251) as u8).collect();
        assert_eq!(
            damaged("noise_index.vcf.gz", &bytes, &noise),
            InputShape::Unknown
        );
        assert_eq!(
            damaged("empty_index.vcf.gz", &bytes, &[]),
            InputShape::Unknown
        );
        // An index whose first reference claims more bins than the file holds
        // is refused by the count check before anything is allocated for them.
        let mut many_bins = b"TBI\x01".to_vec();
        for value in [1i32, 0, 1, 2, 3, i32::from(b'#'), 0, 3] {
            many_bins.extend_from_slice(&value.to_le_bytes());
        }
        many_bins.extend_from_slice(b"21\0");
        many_bins.extend_from_slice(&i32::MAX.to_le_bytes());
        many_bins.extend_from_slice(&[0; 4]);
        let mut encoder = bgzf::io::Writer::new(Vec::new());
        encoder.write_all(&many_bins).unwrap();
        let many_bins = encoder.finish().unwrap();
        assert_eq!(
            damaged("many_bins_index.vcf.gz", &bytes, &many_bins),
            InputShape::Unknown
        );
    }

    /// An indexed file whose one data block declares an ISIZE of 131072, over
    /// the most a BGZF block inflates to, reads as `Unknown`: the block reader
    /// refuses the window's block, and a window with no readable block is no
    /// sample. Eight records fit one block, so the index names that block alone
    /// and every window reads it; unpatched, the file reads as `Indexed`.
    #[test]
    fn a_block_declaring_an_isize_over_the_limit_is_an_unreadable_window() {
        let tmp = tempfile::tempdir().unwrap();
        let mut text =
            String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for record in 1..=8u32 {
            text.push_str(&format!(
                "21\t{}\t.\tA\tG\t.\t.\t.\n",
                10_000_000 + record * 1_000
            ));
        }
        let vcf = tmp.path().join("eight.vcf");
        std::fs::write(&vcf, text).unwrap();
        let good = tmp.path().join("good.vcf.gz");
        bgzip_with_index(&vcf, &good, 1_000, 0);
        let bytes = std::fs::read(&good).unwrap();
        let block = bgzf_block_size(&bytes).unwrap();
        assert_eq!(block + 28, bytes.len(), "one data block and the EOF block");
        assert!(matches!(
            InputShape::read(good.to_str().unwrap()),
            InputShape::Indexed { .. }
        ));
        let mut patched = bytes.clone();
        patched[block - 4..block].copy_from_slice(&131_072u32.to_le_bytes());
        let path = tmp.path().join("oversize_isize.vcf.gz");
        std::fs::write(&path, &patched).unwrap();
        std::fs::copy(
            format!("{}.tbi", good.display()),
            format!("{}.tbi", path.display()),
        )
        .unwrap();
        assert_eq!(
            InputShape::read(path.to_str().unwrap()),
            InputShape::Unknown
        );
    }

    /// A structural-variant fixture bgzipped and indexed by htslib reads as
    /// `Indexed`, with the bytes per record of its own data lines.
    #[test]
    fn htslib_indexed_fixture_reads_its_record_shape() {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/sv_validation/grch37/03_small_indels.vcf.gz");
        let shape = InputShape::read(fixture.to_str().unwrap());
        let InputShape::Indexed {
            inflated_bytes_per_record,
            records_estimate,
            file_bytes,
        } = shape
        else {
            panic!("{shape:?}");
        };
        let mut text = Vec::new();
        bgzf::io::Reader::new(File::open(&fixture).unwrap())
            .read_to_end(&mut text)
            .unwrap();
        let lines: Vec<&[u8]> = text
            .split_inclusive(|&b| b == b'\n')
            .filter(|l| l.first() != Some(&b'#'))
            .collect();
        let true_mean = mean_line_bytes(&lines);
        assert!(
            (true_mean * 4 / 5..=true_mean * 6 / 5).contains(&inflated_bytes_per_record),
            "{inflated_bytes_per_record} bytes per record against a mean of {true_mean}"
        );
        let records = lines.len() as u64;
        assert!(
            (records * 9 / 10..=records * 11 / 10).contains(&records_estimate),
            "{records_estimate} records estimated against {records}"
        );
        assert_eq!(file_bytes, std::fs::metadata(&fixture).unwrap().len());
    }

    #[test]
    fn complete_lines_exclude_the_cut_ends_and_header_lines() {
        assert_eq!(complete_data_lines(b"", false), (0, 0));
        assert_eq!(complete_data_lines(b"no newline", false), (0, 0));
        assert_eq!(complete_data_lines(b"one\n", false), (0, 0));
        assert_eq!(
            complete_data_lines(b"cut\nab\ncdef\n#x\ntail", false),
            (2, 8)
        );
        assert_eq!(complete_data_lines(b"\nab\n", false), (1, 3));
        // At the file start the first line is whole: a data line counts, a
        // header line does not.
        assert_eq!(complete_data_lines(b"one\n", true), (1, 4));
        assert_eq!(complete_data_lines(b"#h\nab\ntail", true), (1, 3));
        assert_eq!(complete_data_lines(b"no newline", true), (0, 0));
    }
}
