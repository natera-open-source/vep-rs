# vep-genome-features

A reader of a converted Ensembl VEP cache that states the features of the genome: the genetic
code of every coding transcript, the biotypes, the coding sequences cut off at an end, the
sequences, the synonyms, and the variant, SIFT and regulatory data the cache carries. It reads
one cache in the JSON layout `tools/vep-storable-to-json` writes, through the cache reader of
[vep-seeded-variants](../vep-seeded-variants/README.md), and writes the features as JSON.
A second subcommand gathers such files into one table, a row a genome.

```text
vep-genome-features extract <converted cache> --out <features.json>
    [--sequence-regions <file>] [--genome <key>] [--threads <n>]
vep-genome-features gather <directory> --out <features.tsv>
vep-genome-features --version
```

Build it with
`cargo build --release --locked --manifest-path tools/vep-genome-features/Cargo.toml`.

## extract

- The one argument is the converted cache directory, the one `--json_cache` takes.
- `--out` names the features file to write.
- `--sequence-regions` names the side file of sequence lengths. Without the option the file
  `sequence_regions.tsv` beside `info.json` is read where it exists.
- `--genome` is written into the file as the genome's key, for `gather`.
- `--threads` is the number of reading threads; the default is one a processor.

The sequences are read side by side, one shard in memory a thread, and the same cache gives
the same bytes whatever the number of threads.

| Key | What it states |
| --- | --- |
| `features_version` | the number of this file's shape, 2 |
| `genome` | the key `--genome` names, or null |
| `species`, `assembly` | the two values of `info.json` |
| `transcripts`, `coding_transcripts` | every transcript, a record counted once a sequence, and those with a coding sequence |
| `genetic_codes` | the coding transcripts by genetic code, with the sequences that carry each code and the first codon of the coding sequences |
| `biotypes` | the transcripts by biotype |
| `incomplete_coding` | the coding transcripts with the `cds_start_NF` or the `cds_end_NF` attribute: at the start, at the end, either, and their share of the coding transcripts to six decimal places |
| `sequences` | the directories under `transcripts`, how many of them hold a transcript, and from the side file its rows, the longest sequence and the circular ones |
| `synonyms` | whether `chr_synonyms.txt` names a synonym, and how many pairs |
| `variant_data` | a sequence directory under `variations` holds a variation file and `info.json` names the columns |
| `sift_data` | `info.json` carries the `sift` marker |
| `regulatory_data` | `info.json` sets `regulatory` |
| `exons_against_position` | the transcripts whose exons, in the cache's order, run against position on the transcript's strand |

A record that several region files of a sequence carry is counted once. Such a record spans
two regions, ends inside one region and stands again in the next region file, or is stored
across the origin of a circular sequence, which the cache files in the sequence's first and
last region file. A record is told apart by its stable id, start, end and `dbID`, as the cache
reader tells it: a stable id that stands at two places of a sequence, or at one place under
two `dbID`s, counts once for each record, and each record is tallied by its own biotype,
coding sequence, exons and peptide.

A transcript that names no genetic code, or names 0, counts under code 1. Under each code,
`opening_codons` tallies the first three bases of the cached coding sequence of every transcript
that starts complete (no `cds_start_NF` attribute and no start phase above zero), with how many
of those transcripts carry a cached peptide that opens with methionine and how many one that
opens with another residue.

A cache without a side file, without a synonyms file or without variant, SIFT or regulatory data
is no error: `sequences.side_file` is null, and the other keys read false or zero.

## gather

`gather` reads every file named `features.json` under the directory, at any depth, and writes a
tab-separated table: the header line, then a row a genome in bytewise order of the key. A
genome's key is the one its file names, or the name of the directory that holds the file where
it names none. The columns:

```text
genome  species  assembly  transcripts  coding_transcripts  genetic_codes
biotypes  incomplete_coding  incomplete_coding_share
sequences_with_features  sequences_with_transcripts  sequences_listed
longest_sequence  longest_sequence_length  circular_sequences
synonym_pairs  variant_data  sift_data  regulatory_data
exons_against_position
```

`genetic_codes` holds `code:count` pairs joined by commas, `biotypes` the number of biotypes,
and the three data columns 1 or 0. A cell without a value holds `-`: the four side-file columns
of a genome read without a side file, a species or an assembly `info.json` does not name, and
the code tally of a genome without a coding transcript.

## Exit codes

A run that cannot read its command line prints `ERROR: [vep-genome-features] <what>` and the
usage and exits 2. A run that fails prints one such line and exits 1: a cache or a shard that
cannot be read, and for `gather` a directory without a features file, a key that stands twice
or a file of another `features_version`.
