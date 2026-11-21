# vep-fasta

Indexed random access into a reference genome FASTA for vep-rs, the Rust reimplementation of the Ensembl Variant Effect Predictor (VEP). The `vep` binary (`vep-cli`) opens the file named by `--fasta` (or chosen through `--fasta_dir`, or found in the cache directory) through this crate and shares one handle with `vep-effects` and the built-in plugins of `vep-builtins`; `vep-cache-builder` opens the genome named by `--genome-fasta` the same way. A read returns bases in upper case, or in the file's case from `sequence_as_written`, and `None` when the contig or the position is not in the index.

## How it works

A query such as `fasta.sequence("chr1", 2, 5)` enters `IndexedFasta::sequence`, which checks the range (`start` at least 1, `end` not before `start`) and resolves the contig through `record`: an exact match in the `.fai` record map first, then a case-insensitive alias lookup that maps `chr1` to `1` (or `1` to `chr1`) and the four spellings of the mitochondrial contig (`MT`, `M`, `chrM`, `chrMT`) to whichever one the index carries, yielding the record with the name the FASTA gives it. When that name is `Y`, the handle was armed with `set_assembly("GRCh37")` or `set_assembly("GRCh38")` and the range overlaps one of the assembly's pseudoautosomal regions (`pars_for_assembly`), the bases of the region are read from the record `X` at the region's offset and the bases outside it from `Y`, as Ensembl VEP's reference reader serves a `Y` read. The matched `FaiRecord` holds the four numbers of a `samtools faidx` line (sequence length, byte offset of the first base, bases per line, bytes per line including the line terminator). The byte offsets of the first and last base, `offset + (pos - 1) / line_bases * line_width + (pos - 1) % line_bases`, bound one span of the uncompressed file; the crate reads that span from the storage, drops every byte inside it that is no printing character (the line ends, and white space before one) and upper-cases the rest, so no FASTA line is parsed at query time. A span whose base count disagrees with the index, or that runs past the end of the data, turns the whole result into `None`; otherwise the caller receives `Some(b"CACA")`.

The storage is one of three. A plain FASTA is memory-mapped (the default) or read into owned memory (`VEP_FASTA_STORAGE=memory`). A FASTA compressed with `bgzip` is mapped compressed: the `.gzi` index (compressed and uncompressed offset of every block that holds data after the first of them) names the block holding an uncompressed offset, the crate inflates that one block of at most 64 KiB with `noodles-bgzf`, and each thread keeps the sixteen blocks it inflated last, since a transcript's exons and its variants' bases sit in a handful of blocks that queries revisit; in the `memory` mode the whole file is inflated once instead. When a BGZF file is opened its `.gzi` is checked against the block layout the file itself gives: one pass reads the header and the stored uncompressed size of every block, and each index entry must sit where the block before it ends, past any block that holds no data (the end-of-file block of a BGZF file another was joined to has no entry), so an index of another file, or of an earlier compression of the same sequence, is refused. A block that declares more than 64 KiB of uncompressed bytes is refused at open as malformed; a block that later inflates to a length the index does not give it is reported once and reads as `None`. The encoding is read from the file's first bytes, not its name: a gzip member whose extra field is exactly the `BC` subfield of the BGZF specification is BGZF, any other gzip stream is refused, since its `.fai` offsets cannot be reached without inflating the whole file.

A file given without its indexes gets them from `index::index_if_missing`, which writes the `.fai` and, for a BGZF file, the `.gzi` beside it as `samtools faidx` writes them, under a `<fasta>.vep.lock` file as Ensembl VEP's reader does. The `.fai` is counted as htslib's `fai_build_core` counts it: a sequence's bases are its printing characters and a line's bytes every byte through its line feed, so white space ending a line is no base, blank lines may follow a sequence, and a name given twice keeps its first sequence. Each index is written into a file of its own name, `<index>.<token>.tmp`, and renamed into place, so no reader sees part of one. A lock found beside the file removes the indexes there before they are written again, no removal that fails ends a run, and a file that cannot be indexed is left with no lock beside it. `vep` opens a single file through `index::open_indexed`, which writes the missing indexes, opens the file, and does both again when a run started meanwhile removed an index; runs that meet on one file wait for no lock, and each reads it. `vep-cache-builder` writes the two indexes of the genome FASTA it installs with `index::write_fai` and `index::write_gzi`.

A directory given as the path is read whole: every FASTA in it that has its indexes (a `.fa` with its `.fai`, or a `.fa.gz` with its `.fai` and `.gzi`) is opened and the sequences are served together, so a directory of one file per chromosome is one reference; two files naming the same sequence are refused.

## Modules

| Module | Role |
| --- | --- |
| `src/lib.rs` | The reader: `FastaError`, `FastaStorageMode` and `IndexedFasta`; the pseudoautosomal tables `Par` and `pars_for_assembly`; the private `FaiRecord`, `FastaStorage` and `BgzfStorage` types and the per-thread block cache; the encoding sniff `encoding_of`; the `.fai` loader `load_fai`; the alias generator `contig_aliases`; the directory scans `indexed_fasta_files`, `first_indexed_fasta` and `fasta_dir_match`; the unit tests |
| `src/index.rs` | The index writers: `index_if_missing` and `open_indexed`, with the lock file and the leftover-lock clean-up; `write_fai` and its line counter `fai_records`; `write_gzi`, a walk of the block headers and trailers; `write_new`, which writes a file under a name of its own and renames it into place; the unit tests |

## Public API

| Item | Role |
| --- | --- |
| `IndexedFasta::from_path(path)` | Opens `path` with its sibling `<path>.fai` (and `<path>.gzi` for a BGZF file), or every indexed FASTA of a directory; the storage mode comes from `VEP_FASTA_STORAGE`, memory-mapped when the variable is unset |
| `IndexedFasta::from_path_mmap(path)` | Same, forcing the memory-mapped storage (blocks inflated on demand for a BGZF file) |
| `IndexedFasta::from_path_in_memory(path)` | Same, reading the whole file, inflated, into owned memory |
| `IndexedFasta::from_path_with_storage(path, mode)` | Same, with an explicit `FastaStorageMode` |
| `IndexedFasta::base(chr, pos)` | One reference base at a 1-based position, as `Option<u8>` |
| `IndexedFasta::sequence(chr, start, end)` | The bases of a 1-based inclusive range, as `Option<Vec<u8>>`; one contig lookup and one span read for the whole range, or one span per piece (the bases before the region, the region's from `X`, the bases after) for a `Y` read the pseudoautosomal redirect splits |
| `IndexedFasta::sequence_as_written(chr, start, end)` | `sequence` with every base in the case the file spells it, for a file whose case carries meaning: an ancestral-sequence FASTA marks a low-confidence call in lower case |
| `IndexedFasta::set_assembly(assembly)` | Arms the pseudoautosomal redirect of `GRCh37` or `GRCh38` (any other name arms nothing): a `base` or `sequence` read on the record `Y` inside a region is served from `X`, `N` when the FASTA has no `X` |
| `IndexedFasta::add_synonyms(synonyms)` | Registers, for each record, the names a `vep_core` `ChromosomeSynonyms` table gives it, so a FASTA named by accession is read under the cache's chromosome names; a record's own name and an alias already present win over a synonym |
| `IndexedFasta::record_name(chr)` | The name the FASTA gives the record `chr` reaches, by its own name, a `chr` spelling or a synonym; `None` when it reaches none |
| `IndexedFasta::sequence_length(chr)` | The number of bases of the record `chr` reaches, as its `.fai` line states it; `None` when it reaches none |
| `IndexedFasta::record_names()` | The records' names in `.fai` order, a directory's files in the order of `files()` |
| `IndexedFasta::sequence_names()` | Every record's name, as its `.fai` line spells it, sorted by byte |
| `IndexedFasta::records()` | Every record's name with the number of bases its `.fai` line states, in the order of `sequence_names()`: the chromosome lengths Ensembl VEP reads from a FASTA for its run statistics |
| `IndexedFasta::fasta_path()` | The path the handle was opened from, a file or a directory, for log lines |
| `IndexedFasta::files()` | The FASTA files the handle reads, by name: one for a file, every indexed one for a directory |
| `IndexedFasta::storage_mode()` | Which `FastaStorageMode` the handle uses |
| `IndexedFasta::corruption()` | The `FastaError` a compressed file's reads have earned, if any: a block that did not inflate as its `.gzi` says; the `vep` binary checks it when its reads are done and ends the run with that error |
| `index::index_if_missing(fasta)` | Writes the indexes the FASTA file lacks beside it, the `.fai` and, for a BGZF file, the `.gzi`, under a `<fasta>.vep.lock` when none of `.fai`, `.index` and `.gzi` is there; a lock found beside the file removes the indexes first. A directory, a gzip file that is not BGZF and a file with both indexes are left as they are |
| `index::open_indexed(fasta)` | `index_if_missing`, then `IndexedFasta::from_path`; an index another run removed in between is written again and the file opened again, up to 32 times, so runs that index one file at once all read it |
| `index::write_fai(fasta)` | Writes `<fasta>.fai` for a plain or BGZF FASTA, byte for byte the file `samtools faidx` writes; a layout no index describes is refused with its line |
| `index::write_gzi(fasta)` | Writes `<fasta>.gzi` for a BGZF file from its block headers and trailers, byte for byte the file `samtools faidx` writes: every block that holds data but the first of them |
| `indexed_fasta_files(dir)` | The indexed `.fa` and `.fa.gz` files of a directory, by name |
| `first_indexed_fasta(dir)` | The first of those: what a cache directory without `--fasta` offers as the reference |
| `fasta_dir_match(dir, species, assembly)` | The file `--fasta_dir` selects: among the `.fa` and `.fa.gz` files whose names carry the species and, when one is given, the assembly, case-insensitively, the first by preference for `toplevel`, then for `.gz`, then by name; `Ok(None)` when none qualifies, `Err` when the directory cannot be read |
| `Par`, `pars_for_assembly(assembly)` | A pseudoautosomal region of `Y` (`start`, `end`, and the `adj` added to reach `X`) and the table of an assembly, Ensembl VEP's `%PARS` |
| `FastaStorageMode` | `Mmap` (file-backed pages, lower resident memory, page faults possible at query time) or `Memory` (owned heap bytes, no file reads after opening); `FromStr` accepts `mmap`, `map`, `mapped`, `memory`, `mem`, `ram` and `owned`, case-insensitively |
| `FastaError` | A missing FASTA, `.fai` or `.gzi`; a gzip file that is not BGZF; a `.gzi` that does not describe the file beside it, with the first disagreement named; a BGZF block that declares more uncompressed bytes than a block holds, with its byte offset; a directory without an indexed FASTA, or whose files name one sequence twice; an open, map or read failure naming the file; a malformed `.fai` line, reported with its line number; a `.fai` record naming more bases than the file holds, with their count; an unrecognised `VEP_FASTA_STORAGE` value; a lock file that cannot be written beside a FASTA without an index, in Ensembl VEP's words, then the cause and the two ways on; a FASTA that cannot be indexed, with the line or the block at fault. Each message is the text the `vep` binary prints after `ERROR: [fasta] `. The type is `non_exhaustive`: a `match` on it outside the crate needs a wildcard arm |
| `MAX_LINE_BYTES` | The longest `.fai` line read, 64 MiB; a longer line is reported as an invalid index line |

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

`IndexedFasta` has two `&mut self` calls, `add_synonyms` and `set_assembly`, made before the handle is shared; after them it is immutable and holds a `memmap2::Mmap` (of the plain or the compressed bytes) or an `Arc<[u8]>`, so callers share one handle across threads behind an `Arc`: the `vep` binary opens the FASTA once, adds the cache's synonyms, arms the run's assembly and hands the same `Arc` to the effects engine and to the built-in plugins that read the reference (LoFTEE, MaxEntScan and UTRAnnotator). The memory-mapped mode requires that the file is not modified while the handle is alive.

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

The unit tests write small FASTA files into a temporary directory, plain and BGZF-compressed in several blocks, with the `.fai` and `.gzi` of each, and check base lookup, the `chr` alias, that the in-memory and memory-mapped modes return the same bases, that the compressed file reads as the plain one across block boundaries, the refusals (a plain gzip, a missing index, an index of another compression of the file, an empty or shifted index, a sequence named in two files of a directory), that a `.fai` naming more bases than the file holds is refused when the file opens, the directory read as a whole, the `--fasta_dir` preference order, that a synonym table reaches a record under the cache's name, and the pseudoautosomal redirect: the tables are Ensembl VEP's, a `Y` read under the armed assembly comes from `X` (from `Y` for the bases outside the region, as `N` when the FASTA has no `X`), and an unarmed handle or a `chrY` record reads as written. [CONTRIBUTING.md](../../CONTRIBUTING.md) has the workspace-wide build, lint and test commands.

## Related

Crates that depend on this one, and what they read:

| Crate | Use |
| --- | --- |
| [vep-cli](../vep-cli/README.md) | Opens the run's FASTA, refuses `--hgvs`, `--hgvsg` and `--shift_3prime 1` without one, and shares the handle with the crates below |
| [vep-effects](../vep-effects/README.md) | Reads single bases to walk an indel to its 3'-most position, and sequence ranges for Human Genome Variation Society (HGVS) notation and the untranslated regions of a transcript |
| [vep-builtins](../vep-builtins/README.md) | The LoFTEE plugin reads the splice-site dinucleotides of an intron, a nine-base window around a splice-acceptor variant covering one reference base, and the ancestral base from its own `human_ancestor_fa` file; MaxEntScan reads the windows at the transcript's splice sites and around the variant that it scores; UTRAnnotator reads the 5' untranslated region and the coding sequence it splices the allele into; AncestralAllele opens its own ancestral-sequence FASTA and reads the variant's span as the file spells it |
| [vep-cache-builder](../vep-cache-builder/README.md) | Reads each exon of a coding transcript from the `--genome-fasta` genome and splices them into the transcript sequence from which its translateable sequence and untranslated regions are cut |

`vep-fasta` itself uses `memmap2`, `noodles` (its `bgzf` module), `thiserror`, `tracing` (the one report of a compressed block that does not inflate as its index says) and `vep-core` (the `ChromosomeSynonyms` table `add_synonyms` reads).

Documentation:

- [docs/cli/cache-species-runtime.md](../../docs/cli/cache-species-runtime.md): the `--fasta`, `--fasta_dir`, `--no_fasta` and `--shift_3prime` flags; [docs/cli/compatibility.md](../../docs/cli/compatibility.md): the `VEP_FASTA_STORAGE` variable.
- [docs/cache-setup.md](../../docs/cache-setup.md): the reference FASTA section (formats, directories, the HGVS refusal) and the `--genome-fasta` input of `vep-cache-builder`.
- [docs/plugins/transcript-model.md](../../docs/plugins/transcript-model.md): what the MaxEntScan and UTRAnnotator plugins read from `--fasta`; [docs/plugins/loftee.md](../../docs/plugins/loftee.md): the LoFTEE splice motifs read from `--fasta` and the ancestral base from `human_ancestor_fa`; [docs/plugins/ancestral-allele.md](../../docs/plugins/ancestral-allele.md): the ancestral-sequence FASTA the AncestralAllele plugin reads in the case the file gives.
- [Root README](../../README.md): the Crates table, one row per crate.
