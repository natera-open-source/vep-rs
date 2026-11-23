# vep-seeded-variants

A seeded generator of input records for any converted Ensembl VEP cache. It reads one cache in
the JSON layout `tools/vep-storable-to-json` writes and writes a VCF of about thirty
records, one for each situation the genome check compares the two engines on
([docs/genome-checks.md](../../docs/genome-checks.md)). The same cache and seed give the same
bytes on any machine, so anyone can regenerate a genome's records.

```text
vep-seeded-variants --cache <dir> --seed <n> --output <file.vcf>
    [--records light|wide] [--per-category <n>] [--sequence <name>]...
    [--fasta <file>] [--sequence-regions <file>] [--summary <file.tsv>]
vep-seeded-variants --version
```

Build it with
`cargo build --release --locked --manifest-path tools/vep-seeded-variants/Cargo.toml`.

## Options

- `--cache` names the converted cache directory, the one `--json_cache` takes.
- `--seed` is a whole number from 0 to 2^64 - 1.
- `--records light`, the default, writes one record a category. `--records wide` writes each
  category twelve times over, spread over sequences and transcript biotypes: a few hundred
  records. `--per-category` changes the twelve.
- `--sequence` keeps the records to the named sequences of the cache; it may be given several
  times. The run reads every sequence of the cache all the same, since the rule for
  [organelle and plasmid sequences](#organelle-and-plasmid-sequences) is the cache's.
- `--fasta` names an indexed reference file: a plain FASTA with its `.fai`, or a block-gzip
  FASTA with its `.fai` and `.gzi`. Every reference allele is then read from that file.
- `--sequence-regions` names the side file of sequence lengths. Without the option the file
  `sequence_regions.tsv` beside `info.json` is read where it exists.
- `--summary` writes one line a category beside the records.
- `--version` prints the tool's version and the generator version.

A run that cannot read its command line prints `ERROR: [vep-seeded-variants] <what>` and the
usage and exits 2; `--summary` naming the file `--output` names is such a command line. A run
that fails prints one such line and exits 1. A reference file that cannot be read, or that does
not hold the bytes its index names, fails the run.

## The records

A record's ID is its category's name and an ordinal, so a difference found on a record names its
situation. The categories, in the order of the summary:

| Category | The record |
| --- | --- |
| `code<N>_changed_codon` | a substitution in a codon genetic code N reads differently from the standard code |
| `code<N>_start_codon` | a substitution in the first codon of a coding sequence on code N |
| `code<N>_stop_gained` | a substitution that turns a sense codon into a stop codon of code N |
| `code<N>_stop_lost` | a substitution that turns the closing stop codon into a sense codon |
| `code<N>_made_codon_synonymous` | a substitution that keeps a residue and makes a codon code N reads differently from the standard code |
| `code<N>_made_codon_missense` | a substitution that changes a residue and makes such a sense codon |
| `code<N>_made_codon_stop` | a substitution in the closing stop codon that makes such a codon |
| `first_sequence`, `last_sequence` | a substitution in an exon on the first and on the last sequence of the cache |
| `organelle` | a substitution in an exon on each organelle or plasmid sequence, three at most |
| `synonym` | a substitution in an exon, its sequence written by a synonym |
| `missense_forward`, `missense_reverse` | a substitution that changes a residue, on each strand |
| `synonymous_forward`, `synonymous_reverse` | a substitution that keeps a residue, on each strand |
| `splice_donor`, `splice_acceptor` | a substitution at one of the two intron bases next to an exon |
| `five_prime_utr`, `three_prime_utr` | a substitution in an untranslated region |
| `intron` | a substitution inside an intron, 20 bases or more from its ends |
| `upstream`, `downstream` | a substitution within 5,000 bases of a transcript's end, outside every transcript |
| `intergenic` | a substitution more than 5,000 bases from every transcript |
| `frameshift_insertion` | an insertion of one or two bases inside the coding part of one exon |
| `frameshift_deletion` | a deletion of one or two bases inside the coding part of one exon |
| `inframe_deletion` | a deletion of one whole codon inside the coding part of one exon |
| `first_base`, `last_base` | a substitution at the first and at the last base of a sequence |
| `origin` | a deletion written across the origin of a circular sequence |

The code categories are written once for each genetic code that a transcript with a coding
sequence names. Under a code that differs from the standard code in its start codons alone, as
code 11 does, `code<N>_changed_codon` is a substitution in the first codon of a coding sequence
whose start mark differs before or after the change. The standard code has no such record.

The three `made_codon` categories are written for a code that translates a codon differently
from the standard code, whether or not a coding sequence of the cache holds such a codon. The
codon before the change is one that both codes translate alike, and the codon after it is one
they translate differently, so each record tells the code from the standard code by itself.
In the closing stop codon the codon made is a sense codon or, under a code that reads a sense
codon of the standard code as a stop codon, a stop codon.

A code has each of the three that a substitution of one base of some codon can serve: code 12
has `code12_made_codon_missense` alone and code 23 `code23_made_codon_stop` alone. A code that
differs from the standard code in its start codons alone has none.

In a transcript with several exons the missense, synonymous, stop-gained, changed-codon,
untranslated-region and indel records, and the two made-codon records of a sense codon, keep
six bases from each end of their exon, which holds them clear of the splice region. The codon,
untranslated-region and indel categories read only a transcript whose cached coding sequence
is the genome's base for base: one without an incomplete-coding-sequence flag, a sequence edit
or a start phase above zero. A splice site record is placed in an intron of 20 bases or more.

The records of `first_sequence`, `last_sequence`, `organelle` and `synonym` stand in an exon of
a transcript of their sequence. A sequence none of whose transcripts gives such a place gets
its record within its first hundred bases, at its first base where no length is known.

The `origin` record is a deletion at the base before the last. Its reference allele is the last
two bases of the sequence followed by the first two, and its alternate allele is its first base.
Both engines read it as a deletion of three bases from the last base on. It needs the
sequence's length and a circular flag of `1`, and a sequence of four bases or more.

A sequence that holds a transcript stored across its origin, one whose start is below 1, takes
no record from its first base to 5,000 bases past that transcript's end, and none from 5,001
bases before that transcript's first base to the sequence's last base. The transcript's first
base is the lowest start among its mapper pairs that lies above its end, or the sequence's
length plus the stored start where the side file gives the length; the lower of the two where
both are given. Where neither is given the sequence takes no record.

### Organelle and plasmid sequences

A converted cache does not say which sequences are organelles or plasmids. The generator takes
for one every sequence one of whose transcripts names another genetic code than most
transcripts of the cache. Where no sequence does and the cache holds two to twelve sequences,
it takes every sequence but the longest, by the side file's lengths where it gives them all and
else by the highest transcript end.

The rule reads the whole cache whichever sequences `--sequence` names. A named sequence gets an
`organelle` record only where the rule takes it for an organelle or a plasmid among all the
sequences of the cache. The three at most are the first three such sequences among the named
ones.

### Synonyms

The `synonym` record is written on a sequence to which `chr_synonyms.txt` gives a name that
reaches it by the synonym alone: a name the cache holds no sequence by, and that the `chr`
prefix rule does not turn into the cache's own name.

## Reference alleles

A converted cache holds sequence only inside transcripts: the coding sequence and the two
untranslated regions. A record's reference base inside a coding sequence is the cached base,
turned to the forward strand, whichever category wrote the record. Elsewhere no base is known
and the reference base is drawn from the seed. Neither engine checks it under default flags
without a reference file.

With `--fasta` every reference allele is read from the file. A candidate whose bases there are
not `A`, `C`, `G` or `T`, or differ from the bases the cache gives, is passed over, and the
summary counts the candidates passed over.

## The determinism contract

- The random stream is SplitMix64 and a keyed FNV-1a hash, both written in this crate
  (`src/random.rs`), so a dependency's release cannot change a genome's records.
- Sequences are read in the order of `contig_order` in `info.json`, a sequence's shards by the
  start their names state, and a sequence's transcripts by start, then stable id. A sequence
  is a name that has a directory under `transcripts`; a directory `contig_order` does not list
  is read after the listed ones, in bytewise name order.
- A transcript is read once a sequence, from the first shard that holds it. An Ensembl VEP
  cache files a gene's transcripts in every region file the gene reaches, so the shards of a
  sequence repeat transcripts: a record under the stable id, start, end and `dbID` of a
  transcript already read is that transcript again. A repeat is no second candidate and does
  not count among the transcripts that name a genetic code. A record under the same stable id
  and another `dbID`, at the same place or at another, is another transcript.
- Every candidate has a priority, a hash of the seed, the category, the sequence and the
  transcript. A category keeps the candidates of the lowest priority, so the kept set does not
  depend on the machine, the thread count or the order a directory lists its files in.
- One sequence's transcripts are in memory at a time, and a category keeps a bounded number of
  candidates.
- `GENERATOR_VERSION` names the algorithm. The output's header and `--version` state it, and
  it changes whenever the same cache and seed would give other bytes. `cargo test` holds the
  light records of four committed cache slices for seed 1 to the files under `tests/expected/`,
  and with them the records and the summary of a made-up cache of fourteen sequences for seed 1
  under two `--sequence` options, of a made-up cache that files a gene's transcripts in two
  region files and of a made-up cache of one sequence under genetic code 4.

## The output

The VCF opens with `##fileformat`, a `##source` line that names the tool and the generator
version, the seed, the mode and where the reference alleles were read. The eight fixed columns
follow, and the records stand in the cache's sequence order, then by position.

## The side file

`sequence_regions.tsv` holds one row a sequence under the header `name`, `length`, `circular`;
`circular` is `1`, `0` or `-`. The cache converter writes it. Without it no length is known:
`last_base` and `origin` are left out, and `upstream`, `downstream` and `intergenic` records
stay below the highest transcript end of their sequence.

## The summary

One line a category under the header `category`, `wanted`, `written`, `passed_over`, `reason`.
`reason` is `-` where every wanted record was written. Else it is one of:

| Reason | Where |
| --- | --- |
| `the standard code` | `code1_changed_codon`, which wants no record |
| `no codon this code reads differently` | `code<N>_changed_codon` |
| `no substitution makes a codon this code reads differently` | the three `code<N>_made_codon` categories, where a coding sequence on the code is read |
| `no complete coding sequence` | `code<N>_start_codon`, `code<N>_stop_gained`, the indel categories, and the `code<N>_made_codon` categories where no coding sequence on the code is read |
| `no table for this code` | the categories of a code the annotator has no table for |
| `no complete stop codon` | `code<N>_stop_lost` |
| `one sequence` | `last_sequence` |
| `no organelle or plasmid sequence` | `organelle`: the cache holds none, or `--sequence` names none of them |
| `no synonym` | `synonym` |
| `no place on the sequence` | `first_sequence`, `last_sequence`, `organelle`, `synonym`, `first_base`, `last_base`, `origin`: the sequence takes no record there |
| `no transcript on this strand`, `no complete coding sequence on this strand` | the missense and synonymous categories |
| `no multi-exon transcript`, `no intron long enough` | the splice site categories and `intron` |
| `no untranslated region` | the two untranslated-region categories |
| `no room beside a transcript` | `upstream`, `downstream` |
| `no gap between genes` | `intergenic` |
| `no sequence length`, `no circular sequence` | `last_base`, `origin` |
| `reference masked` | any category whose every candidate the reference file passed over |
| `few candidates` | any category that wrote some records and fewer than it wanted |

## The library

The crate's library `vep_seeded_variants` carries the same run for a Rust caller. Its module
`cache` is the reader of a converted cache: `info.json`, the synonym pairs, the side file and
an iterator over one sequence's transcripts in shard order, each given once and read down to
its place, its exons and its coding sequence.
