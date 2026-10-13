# vep-fasta

Indexed random access into a reference genome FASTA for vep-rs, the Rust reimplementation of the Ensembl Variant Effect Predictor (VEP). The `vep` binary (`vep-cli`) opens the file named by `--fasta` through this crate and shares one handle with `vep-effects` and the built-in plugins of `vep-builtins`; `vep-cache-builder` opens the genome named by `--genome-fasta` the same way. Every query returns reference bases as upper-case bytes, or `None` when the contig or the position is not in the index.

## How it works

A query such as `fasta.sequence("chr1", 2, 5)` enters `IndexedFasta::sequence`, which checks the range (`start` at least 1, `end` not before `start`) and resolves the contig through `record_for`: an exact match in the `.fai` record map first, then a case-insensitive alias lookup that maps `chr1` to `1` (or `1` to `chr1`) and the four spellings of the mitochondrial contig (`MT`, `M`, `chrM`, `chrMT`) to whichever one the index carries. The matched `FaiRecord` holds the four numbers of a `samtools faidx` line (sequence length, byte offset of the first base, bases per line, bytes per line including the line terminator). For each 1-based position the crate computes `offset + (pos - 1) / line_bases * line_width + (pos - 1) % line_bases` and reads that one byte from the storage, either the memory-mapped file or the owned in-memory copy, so no FASTA line is parsed at query time. A byte that lands on a line terminator, or an offset past the end of the data, turns the whole result into `None`; otherwise each byte is upper-cased and appended, and the caller receives `Some(b"CACA")`.

## Modules

| Module | Role |
| --- | --- |
| `src/lib.rs` | The whole crate: `FastaError`, `FastaStorageMode` and `IndexedFasta`; the private `FaiRecord` and `FastaStorage` types; the `.fai` loader `load_fai`; the alias generator `contig_aliases`; two unit tests |

## Public API

| Item | Role |
| --- | --- |
| `IndexedFasta::from_path(path)` | Opens `path` and its sibling `<path>.fai`; the storage mode comes from `VEP_FASTA_STORAGE`, memory-mapped when the variable is unset |
| `IndexedFasta::from_path_mmap(path)` | Same, forcing the memory-mapped storage |
| `IndexedFasta::from_path_in_memory(path)` | Same, reading the whole file into owned memory |
| `IndexedFasta::from_path_with_storage(path, mode)` | Same, with an explicit `FastaStorageMode` |
| `IndexedFasta::base(chr, pos)` | One reference base at a 1-based position, as `Option<u8>` |
| `IndexedFasta::sequence(chr, start, end)` | The bases of a 1-based inclusive range, as `Option<Vec<u8>>`; one contig lookup for the whole range |
| `IndexedFasta::fasta_path()` | The path the handle was opened from, for log lines |
| `IndexedFasta::storage_mode()` | Which `FastaStorageMode` the handle uses |
| `FastaStorageMode` | `Mmap` (file-backed pages, lower resident memory, page faults possible at query time) or `Memory` (owned heap bytes, no file reads after opening); `FromStr` accepts `mmap`, `map`, `mapped`, `memory`, `mem`, `ram` and `owned`, case-insensitively |
| `FastaError` | A missing FASTA or `.fai`; an open, map or read failure; a malformed `.fai` line, reported with its line number; an unrecognised `VEP_FASTA_STORAGE` value |

A caller opens the genome once and queries it by contig name and 1-based position:

```rust
use vep_fasta::{FastaStorageMode, IndexedFasta};

fn main() -> Result<(), vep_fasta::FastaError> {
    // Expects Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai beside the FASTA.
    let fasta = IndexedFasta::from_path_with_storage(
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
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

`IndexedFasta` is immutable after construction and holds either a `memmap2::Mmap` or an `Arc<[u8]>`, so callers share one handle across threads behind an `Arc`: the `vep` binary opens the `--fasta` file once and hands the same `Arc` to the effects engine and to the LoFTEE plugin. The memory-mapped mode requires that the file is not modified while the handle is alive.

The index is the format `samtools faidx` writes: one line per contig, five or more tab-separated fields, of which the crate reads the first five.

```text
name	length	offset	line_bases	line_width
1	6	3	6	7
```

Blank lines are skipped. A line with fewer than five fields, or a non-numeric field, is reported as `FastaError::InvalidFaiLine` with its line number.

## Run it

```bash
cargo build -p vep-fasta
cargo test -p vep-fasta
```

The two unit tests write a six-base FASTA and its index into a temporary directory and check base lookup, the `chr` alias, and that the in-memory and memory-mapped modes return the same bases. [CONTRIBUTING.md](../../CONTRIBUTING.md) has the workspace-wide build, lint and test commands.

## Related

Crates that depend on this one, and what they read:

| Crate | Use |
| --- | --- |
| [vep-cli](../vep-cli/README.md) | Opens the `--fasta` file, refuses `--shift_3prime` without it, and shares the handle with the crates below |
| [vep-effects](../vep-effects/README.md) | Reads single bases to walk an indel to its 3'-most position, and sequence ranges for Human Genome Variation Society (HGVS) notation and the untranslated regions of a transcript |
| [vep-builtins](../vep-builtins/README.md) | The LoFTEE plugin reads the splice-site dinucleotides of an intron, a nine-base window around a single-nucleotide variant, and the ancestral base from its own `human_ancestor_fa` file |
| [vep-cache-builder](../vep-cache-builder/README.md) | Reads each coding-sequence region of a transcript from the `--genome-fasta` genome to build its translateable sequence |

`vep-fasta` itself uses `memmap2` and `thiserror`; its `Cargo.toml` also declares `vep-core`, `noodles` and `tracing`, from which the source imports nothing.

Documentation:

- [docs/cli-reference.md](../../docs/cli-reference.md): the `--fasta` and `--shift_3prime` flags and the `VEP_FASTA_STORAGE` variable.
- [docs/plugins.md](../../docs/plugins.md): what the LoFTEE plugin reads from `--fasta` and from `human_ancestor_fa`.
- [docs/cache-setup.md](../../docs/cache-setup.md): the `--genome-fasta` input of `vep-cache-builder`.
- [Root README](../../README.md): the `--fasta` example under Usage.
