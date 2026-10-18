# vep-fasta

Indexed random access into a reference genome FASTA for vep-rs, the Rust reimplementation of the Ensembl Variant Effect Predictor (VEP). The `vep` binary (`vep-cli`) opens the file named by `--fasta` (or chosen through `--fasta_dir`, or found in the cache directory) through this crate and shares one handle with `vep-effects` and the built-in plugins of `vep-builtins`; `vep-cache-builder` opens the genome named by `--genome-fasta` the same way. Every query returns reference bases as upper-case bytes, or `None` when the contig or the position is not in the index.

## How it works

A query such as `fasta.sequence("chr1", 2, 5)` enters `IndexedFasta::sequence`, which checks the range (`start` at least 1, `end` not before `start`) and resolves the contig through `record_for`: an exact match in the `.fai` record map first, then a case-insensitive alias lookup that maps `chr1` to `1` (or `1` to `chr1`) and the four spellings of the mitochondrial contig (`MT`, `M`, `chrM`, `chrMT`) to whichever one the index carries. The matched `FaiRecord` holds the four numbers of a `samtools faidx` line (sequence length, byte offset of the first base, bases per line, bytes per line including the line terminator). The byte offsets of the first and last base, `offset + (pos - 1) / line_bases * line_width + (pos - 1) % line_bases`, bound one span of the uncompressed file; the crate reads that span from the storage, drops the line terminators inside it and upper-cases the rest, so no FASTA line is parsed at query time. A span whose base count disagrees with the index, or that runs past the end of the data, turns the whole result into `None`; otherwise the caller receives `Some(b"CACA")`.

The storage is one of three. A plain FASTA is memory-mapped (the default) or read into owned memory (`VEP_FASTA_STORAGE=memory`). A FASTA compressed with `bgzip` is mapped compressed: the `.gzi` index (compressed and uncompressed offset of every block after the first) names the block holding an uncompressed offset, the crate inflates that one block of at most 64 KiB with `noodles-bgzf`, and each thread keeps the sixteen blocks it inflated last, since a transcript's exons and its variants' bases sit in a handful of blocks that queries revisit; in the `memory` mode the whole file is inflated once instead. When a BGZF file is opened its `.gzi` is checked against the block layout the file's own headers give (the first block, a spread of entries and the last), so an index of another file, or of an earlier compression of the same sequence, is refused; a block that later inflates to a length the index does not give it is reported once and reads as `None`. The encoding is read from the file's first bytes, not its name: a gzip member carrying the `BC` extra subfield of the BGZF specification is BGZF, any other gzip stream is refused, since its `.fai` offsets cannot be reached without inflating the whole file.

A directory given as the path is read whole: every FASTA in it that has its indexes (a `.fa` with its `.fai`, or a `.fa.gz` with its `.fai` and `.gzi`) is opened and the sequences are served together, so a directory of one file per chromosome is one reference; two files naming the same sequence are refused.

## Modules

| Module | Role |
| --- | --- |
| `src/lib.rs` | The whole crate: `FastaError`, `FastaStorageMode` and `IndexedFasta`; the private `FaiRecord`, `FastaStorage` and `BgzfStorage` types and the per-thread block cache; the encoding sniff `encoding_of`; the `.fai` loader `load_fai`; the alias generator `contig_aliases`; the directory scans `indexed_fasta_files`, `first_indexed_fasta` and `fasta_dir_match`; the unit tests |

## Public API

| Item | Role |
| --- | --- |
| `IndexedFasta::from_path(path)` | Opens `path` with its sibling `<path>.fai` (and `<path>.gzi` for a BGZF file), or every indexed FASTA of a directory; the storage mode comes from `VEP_FASTA_STORAGE`, memory-mapped when the variable is unset |
| `IndexedFasta::from_path_mmap(path)` | Same, forcing the memory-mapped storage (blocks inflated on demand for a BGZF file) |
| `IndexedFasta::from_path_in_memory(path)` | Same, reading the whole file, inflated, into owned memory |
| `IndexedFasta::from_path_with_storage(path, mode)` | Same, with an explicit `FastaStorageMode` |
| `IndexedFasta::base(chr, pos)` | One reference base at a 1-based position, as `Option<u8>` |
| `IndexedFasta::sequence(chr, start, end)` | The bases of a 1-based inclusive range, as `Option<Vec<u8>>`; one contig lookup and one span read for the whole range |
| `IndexedFasta::fasta_path()` | The path the handle was opened from, a file or a directory, for log lines |
| `IndexedFasta::files()` | The FASTA files the handle reads, by name: one for a file, every indexed one for a directory |
| `IndexedFasta::storage_mode()` | Which `FastaStorageMode` the handle uses |
| `indexed_fasta_files(dir)` | The indexed `.fa` and `.fa.gz` files of a directory, by name |
| `first_indexed_fasta(dir)` | The first of those: what a cache directory without `--fasta` offers as the reference |
| `fasta_dir_match(dir, species, assembly)` | The file `--fasta_dir` selects: among the `.fa` and `.fa.gz` files whose names carry the species and the assembly, case-insensitively, the first by preference for `toplevel`, then for `.gz`, then by name; `Ok(None)` when none qualifies, `Err` when the directory cannot be read |
| `FastaStorageMode` | `Mmap` (file-backed pages, lower resident memory, page faults possible at query time) or `Memory` (owned heap bytes, no file reads after opening); `FromStr` accepts `mmap`, `map`, `mapped`, `memory`, `mem`, `ram` and `owned`, case-insensitively |
| `FastaError` | A missing FASTA, `.fai` or `.gzi`; a gzip file that is not BGZF; a `.gzi` that does not describe the file beside it, with the first disagreement named; a directory without an indexed FASTA, or whose files name one sequence twice; an open, map or read failure naming the file; a malformed `.fai` line, reported with its line number; an unrecognised `VEP_FASTA_STORAGE` value. Each message is the text the `vep` binary prints after `ERROR: [fasta] ` |

A caller opens the genome once and queries it by contig name and 1-based position:

```rust
use vep_fasta::{FastaStorageMode, IndexedFasta};

fn main() -> Result<(), vep_fasta::FastaError> {
    // Expects Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.fai and .gzi beside the FASTA,
    // as samtools faidx writes them for a bgzipped file.
    let fasta = IndexedFasta::from_path_with_storage(
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        FastaStorageMode::Mmap,
    )?;
    // Coordinates are 1-based; "1" and "chr1" name the same contig.
    let base = fasta.base("1", 1_000_000).map(char::from);
    let window = fasta
        .sequence("chr1", 999_996, 1_000_004)
        .map(|bases| String::from_utf8_lossy(&bases).into_owned());
    println!("{base:?} {window:?}");
    Ok(())
}
```

`IndexedFasta` is immutable after construction and holds a `memmap2::Mmap` (of the plain or the compressed bytes) or an `Arc<[u8]>`, so callers share one handle across threads behind an `Arc`: the `vep` binary opens the FASTA once and hands the same `Arc` to the effects engine and to the LoFTEE plugin. The memory-mapped mode requires that the file is not modified while the handle is alive.

The index is the format `samtools faidx` writes: one line per contig, five or more tab-separated fields, of which the crate reads the first five. For a bgzipped file the offsets are those of the uncompressed sequence, as `samtools faidx` writes them, and the `.gzi` beside it is the one `bgzip -r` or `samtools faidx` writes.

```text
name	length	offset	line_bases	line_width
1	6	3	6	7
```

Blank lines are skipped. A line with fewer than five fields, a non-numeric field or zero bases per line is reported as `FastaError::InvalidFaiLine` with its line number.

## Run it

```bash
cargo build -p vep-fasta
cargo test -p vep-fasta
```

The unit tests write small FASTA files into a temporary directory, plain and BGZF-compressed in several blocks, with the `.fai` and `.gzi` of each, and check base lookup, the `chr` alias, that the in-memory and memory-mapped modes return the same bases, that the compressed file reads as the plain one across block boundaries, the refusals (a plain gzip, a missing index, an index of another compression of the file, an empty or shifted index, a sequence named in two files of a directory), that a `.fai` longer than the file reads `None` past the data, the directory read as a whole and the `--fasta_dir` preference order. [CONTRIBUTING.md](../../CONTRIBUTING.md) has the workspace-wide build, lint and test commands.

## Related

Crates that depend on this one, and what they read:

| Crate | Use |
| --- | --- |
| [vep-cli](../vep-cli/README.md) | Opens the run's FASTA, refuses `--hgvs`, `--hgvsg` and `--shift_3prime 1` without one, and shares the handle with the crates below |
| [vep-effects](../vep-effects/README.md) | Reads single bases to walk an indel to its 3'-most position, and sequence ranges for Human Genome Variation Society (HGVS) notation and the untranslated regions of a transcript |
| [vep-builtins](../vep-builtins/README.md) | The LoFTEE plugin reads the splice-site dinucleotides of an intron, a nine-base window around a single-nucleotide variant, and the ancestral base from its own `human_ancestor_fa` file |
| [vep-cache-builder](../vep-cache-builder/README.md) | Reads each coding-sequence region of a transcript from the `--genome-fasta` genome to build its translateable sequence |

`vep-fasta` itself uses `memmap2`, `noodles` (its `bgzf` module) and `thiserror`; its `Cargo.toml` also declares `vep-core` and `tracing`, from which the source imports nothing.

Documentation:

- [docs/cli-reference.md](../../docs/cli-reference.md): the `--fasta`, `--fasta_dir`, `--no_fasta` and `--shift_3prime` flags and the `VEP_FASTA_STORAGE` variable.
- [docs/cache-setup.md](../../docs/cache-setup.md): the reference FASTA section (formats, directories, the HGVS refusal) and the `--genome-fasta` input of `vep-cache-builder`.
- [docs/plugins.md](../../docs/plugins.md): what the LoFTEE plugin reads from `--fasta` and from `human_ancestor_fa`.
- [Root README](../../README.md): the `--fasta` example under Usage.
