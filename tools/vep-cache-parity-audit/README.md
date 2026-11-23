# vep-cache-parity-audit

Audits the conversion of an Ensembl VEP Storable cache into the JSON cache vep-rs annotates from
(`tools/vep-storable-to-json`). It takes a seeded sample of transcripts from the Storable cache,
finds each in the JSON cache, walks every field path of both and classes each Storable path against
the converter's rules, so that a field the conversion leaves out without a rule shows as a gap
before it shows as a difference in annotation.

```bash
vep-cache-parity-audit <storable_cache_dir> <json_cache_dir> [N] [SEED] > manifest.json
```

- `storable_cache_dir`: a directory above the cache's shards, the cache directory or a parent of it.
- `json_cache_dir`: the converter's output directory; every `*.json` under it is read.
- `N`: how many transcripts to sample (default 1000).
- `SEED`: the seed of the order the Storable files are taken in (default 1).

The tool reads the Storable images itself and needs no Perl. `VEP_CACHE_PARITY_AUDIT_THREADS` sets
the number of worker threads (default: the host's cores).

## The manifest

One JSON object on standard output, on one line without a final newline:

| Key | What it holds |
|---|---|
| `sampled_transcripts`, `sample_seed`, `storable_files_scanned` | the sample's size, its seed, the Storable files found |
| `json_transcripts_indexed`, `matched_in_json` | distinct transcript ids of the JSON cache; sampled transcripts found there |
| `unmatched_in_json_by_contig_class` | sampled transcripts not found, under `canonical_contig` (a directory of digits, `X`, `Y` or `MT`) and `non_canonical_contig` |
| `identical_field_counts` | paths present under the same path with the same kind |
| `renamed_field_counts` | paths present under the converter's renamed path with the same kind |
| `conversion_rules` | paths present with another kind, as `<path>:<kind>-><kind>` |
| `dropped_undef_counts` | paths absent from the JSON cache whose Storable value is undef |
| `dropped_by_converter` | paths absent by a named rule, as `<path>: <rule>` |
| `derived_on_json_side` | JSON paths no Storable path maps to and a rule produces |
| `parity_gaps` | `missing_in_json:<path>` and `only_in_json:<path>`: paths no rule accounts for |
| `parity_gap_examples` | up to three sampled ids for each gap |

Every count is the number of sampled transcripts that hold the path.

## Contract

`tests/fixtures/expected/` holds the reference bytes of 24 runs, made by perl 5.34.0 with Storable
3.23, JSON::PP 4.06 and File::Find 1.39, and `cargo test` holds the binary to them. The manifest of
a reference run is the one perl wrote, with the keys of every object in the order of their names
and, in a run whose sample spans several files, the number of the ids under `parity_gap_examples`
in the place of the ids; standard error and the exit code are fixed to the byte. Where a rule
belongs to a Perl module it is named here by its Perl name.

**Command line.** The two directories are required; an argument that is empty or `0` counts as
missing, as a false Perl string does, and the tool prints `Usage: <program> storable_cache_dir
json_cache_dir [N]` on standard error and exits 255. Arguments past the fourth are not read.

**The Storable files** are found as `File::Find::find` finds them: the files of a directory in the
order the system lists them, then each subdirectory in that order, whole, before the next; a
symbolic link to a directory is not entered, and a root that is one file is that file. A file is
taken when its name ends in `.gz` and not in `_reg.gz` or `_var.gz`. A root that does not exist is
`Can't stat <root>: <the system's reason>`, and a root without such a file ends the run with `No
.gz Storable chunks found under <root>`.

**The sample.** The files are put in the order `sort { rand(1) <=> rand(1) }` leaves them in after
`srand(SEED)`: Perl's merge sort (pairs first, longer runs where eight pairs in a row compare alike,
then merges two by two) with a comparison that draws two numbers of Perl's `drand48` each time. The
seed is the integer `SEED` opens with, cut to 32 bits; a text that is no number or has an exponent
seeds with all bits set. The files are then read in that order, each as `gunzip -c` writes it and
Storable's `fd_retrieve` reads that, until `N` transcripts are sampled (`N` as a Perl number: text
that is no number is 0, a fraction is rounded up). A file that is no image with an unblessed hash
at its root is passed over. In an image, each value of the root hash gives its transcripts: of a
list, the items that are references to a hash with a `stable_id` key; of an unblessed hash, itself
when it has that key. So the sample comes from the first files of the order that reach `N`, and is
no uniform draw over the cache.

**The JSON cache.** Every file under `json_cache_dir` whose name ends in `.json` is read, in the
byte order of their paths, as `JSON::PP`'s `decode_json` reads it; a file it dies on is passed over.
A file is a list of transcripts or an object of such lists; a transcript is an object whose
`stable_id` is a true value, and is indexed under that id (a number as Perl writes it, `true` as
`1`). A later transcript of an id replaces an earlier one. A sampled transcript whose `stable_id`
is a false value is counted in the sample and compared no further; one whose id the index lacks is
unmatched.

**Paths and kinds.** A path joins the keys from the transcript down with `.`; a key under a path
that is empty or `0` stands alone. The items of an array share the path `<array>[]`. A kind is
`scalar`, `undef`, `HASH` or `ARRAY`, a reference to anything else `SCALAR`, `REF`, `REGEXP` or
`VSTRING`, with the class before it for a blessed value (`Bio::EnsEMBL::Exon(HASH)`); a JSON `true`
or `false` is `JSON::PP::Boolean(SCALAR)`. A blessed hash that is being walked further up the path
is `<class>(HASH,cycle)` and is not walked again; a hash reached twice elsewhere is walked at each
place. Where two values share a path the last one's kind stands.

**The classes.** The Storable paths of a transcript are taken shortest first, then in code point
order. A path whose parent was written as a value that is no object or list (`conversion_rules`
with a kind that holds neither `HASH` nor `ARRAY`) is passed over with it. The JSON path of a
Storable path is the path with its first key renamed (`_gene_stable_id`, `_gene_symbol`,
`_gene_symbol_source`, `_gene_hgnc_id`, `_gene_phenotype`, `_gene_version`, `_ccds`, `_refseq`,
`_protein`, `_swissprot`, `_trembl`, `_uniparc`, `_uniprot_isoform`, `_source_cache`,
`_gene_symbol_id`, `_bam_edit_status` lose their underscore; `_trans_exon_array` is `exons`;
`_variation_effect_feature_cache` is `variation_effect_feature_cache`) and a prediction's `matrix`
as `predictions_data`. A path present there is identical, renamed or converted. One absent is
dropped as undef when its kind is `undef`, else dropped by the first rule that fits, else a gap:

| Paths | Rule |
|---|---|
| `slice` and below | the chromosome is the shard directory |
| `_gene` and below | carried by the `_gene_*` renames |
| `_vep_lazy_loaded` | a flag of the cache loader |
| `translation.start_exon`, `translation.end_exon` and below | derivable from the exons |
| `translation.transcript` | the back-reference to the transcript |
| the cache's `mapper.exon_coord_mapper` and below | replaced by `mapper.pairs` |
| `strand` of an exon, a sorted exon or an intron | the transcript's strand |
| an intron's `prev`, `next` and `slice` and below | the exons array; the shard directory |
| a sequence edit's `description` | read by nothing |
| a prediction's `matrix_compressed` | constant |
| below the cache's `five_prime_utr` and `three_prime_utr` | flattened to a sequence string |

A JSON path no Storable path accounted for is derived when it is the cache's `mapper.pairs` or below,
`mapper.pair_count`, `seq_edits` or a sequence edit's `length_diff`, else a gap. `src/audit.rs`
holds each rule with its text as the manifest writes it.

**Order and cores.** The Storable files are read ahead of the sample on the worker threads and
taken in the sample's order; the JSON files are read on the worker threads, one file a thread, and
taken in their order. The tool holds the field paths of the sampled transcripts and the ids of the
JSON cache, and a shard or a JSON file only while it is read: never a whole cache.

**Where the tool does not write what Perl writes.** A reference run that meets one of these cases
holds the tool's bytes for it, the ids of the second case excepted.

- The order of the keys of the manifest's objects. `encode_json` writes them in the order of a Perl
  hash, which changes from run to run; the tool writes them in the order of their names. The values
  are the same: counts as numbers, a seed that is given as a string and the seed of a run without
  one as the number `1`, an id that the cache stores as a number as that number.
- A sample over several files follows the order the system lists the cache's directories in, as
  `File::Find` reports them: the same cache copied to another file system can give another sample
  for one seed. The reference run of such a sample holds, under `parity_gap_examples`, the number
  of the ids of each gap in the place of the ids.
- A line that Perl ends with the place of a `die` or of a warning (` at <file> line <n>.`, also on
  a line of its own) is printed without it; Perl's own warnings (`Argument "..." isn't numeric`,
  `Integer overflow in srand`) and the two lines `gzip` writes for a file it cannot read are not
  printed.
- A run that dies exits 255. Perl exits with the error number its last failed system call left.
- An array that holds itself: the same walk in Perl has no end; the tool stops and says so.
- A root hash, and a JSON object of lists, that hold several lists: Perl takes them in the order of
  a hash; the tool takes them in the order the image or the file stores them.
- A `stable_id` that is a reference: Perl reads it as its address in memory. The tool counts such
  a transcript of the Storable cache as unmatched, and one of the JSON cache as one more id.
- A file name with a single quote in it is read as the file it names; a shell given `gunzip -c`
  and that name in single quotes does not read that file.
- A JSON file in UTF-16 or UTF-32, which `decode_json` transcodes, is passed over.

## Layout

- `src/lib.rs`: the modules, and those compiled from other tools by `#[path]`: the reader of Storable
  images, Perl's handling of a retrieved value, Storable's outcomes, `gzip -dc`'s output and the
  worker threads from `tools/vep-storable-to-json`; a retrieved shard and the entries of its hashes
  from `tools/vep-storable-digest`.
- `src/find.rs`: the order `File::Find` reports files in.
- `src/perlsort.rs`: Perl's merge sort, `srand` and `rand`.
- `src/jsonpp.rs`: JSON text as `JSON::PP` reads it, and a string as it writes one.
- `src/fields.rs`: the paths and kinds of a transcript on both sides.
- `src/audit.rs`: the sample, the index, the rules and the manifest; `src/main.rs`: the command line.
- `tests/golden.rs` with `tests/fixtures/storable/`, `tests/fixtures/json/` and
  `tests/fixtures/expected/`: the binary against the reference bytes of 24 runs, each on one thread
  and on four. `tests/fixtures/README.md` says what the caches hold.
- `tests/fixtures/perlsort/cases.txt` (126 orders Perl's randomising sort leaves),
  `tests/fixtures/jsonpp/cases.txt` (what `JSON::PP` takes, refuses and writes) and
  `tests/fixtures/rules/cases.tsv` (the mapped path and the rule of each path, as perl 5.34.0
  gives them): the unit tests of `src/perlsort.rs`, `src/jsonpp.rs` and `src/audit.rs` hold the
  code to them.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest.

```bash
cargo build --release --locked --manifest-path tools/vep-cache-parity-audit/Cargo.toml
cargo test --locked --manifest-path tools/vep-cache-parity-audit/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-cache-parity-audit/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-cache-parity-audit/Cargo.toml -- --check
```

The reference bytes are the contract, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the files it moves under `tests/fixtures/expected/<case>/` with
the binary's own output, so that the diff of the fixtures shows the change.
