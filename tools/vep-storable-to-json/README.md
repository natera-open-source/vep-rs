# vep-storable-to-json

Converts an Ensembl VEP Storable cache into the JSON cache vep-rs annotates from. The Storable
cache is the one Ensembl VEP's installer downloads: `info.txt`, one directory a chromosome, and in
each the gzipped Perl Storable shards `<start>-<end>.gz` (transcripts), `<start>-<end>_reg.gz`
(regulatory features) and the variation files. The tool reads the Storable images itself and needs
no Perl.

```bash
vep-storable-to-json <cache_dir> <output_dir> [chr1,chr2,...]
```

It writes, under `<output_dir>`:

| File | From |
|---|---|
| `transcripts/<chr>/<start>-<end>.json` | every `<start>-<end>.gz` of the chromosome directory |
| `regulatory/<chr>/<start>-<end>.json` | every `<start>-<end>_reg.gz` |
| `regulatory/binding_matrices.json` | the binding matrices the motifs of those shards name, each once |
| `variations/<chr>/` | `all_vars.gz` with its `.csi` or `.tbi` index, else every `<start>-<end>_var.gz`, copied byte for byte |
| `chr_synonyms.txt` | copied when the cache has one |
| `sequence_regions.tsv` | the slice of each chromosome's transcripts: its length and circular flag |
| `info.json` | `info.txt`, the cache directory's path, and what the run wrote |

The third argument limits the run to the named chromosome directories. `VEP_STORABLE_TO_JSON_THREADS`
sets the number of worker threads (default: the host's cores).

## Contract

The output is fixed to the byte: every file, every line on standard output and standard error, and
the exit code. `tests/fixtures/expected/` holds the reference bytes of 62 runs, made by perl 5.34.0
with Storable 3.23, JSON::XS 4.04 and zlib 1.2.11, and `cargo test` holds the binary to them. Where
a byte could be written in more than one way it follows the rule of the Perl module that defines the
format, named here by its Perl name.

**Command line.** `cache_dir` and `output_dir` are required; an argument that is empty or `0` counts
as missing, as a false Perl string does, and the tool prints `Usage: <program> cache_dir output_dir
[chr1,chr2,...]` on standard error and exits 255. The chromosome list is split at commas as Perl's
`split /,/` splits it (empty fields at the end are dropped, an empty field before a name is a name);
a list that is empty or `0` converts every chromosome. Arguments past the third are ignored.

**The run, in order.**

1. `info.txt` is read as Ensembl VEP's own reader takes it: a line that starts with `#` or has no
   tab is skipped; `source_<key>` goes under `source_versions`; `variation_cols` becomes a list with
   its empty items dropped; `cell_types` and `regulatory` are dropped; any other line is a key and a
   string value, unless the value is `-`; a later line replaces an earlier one. Without the file the
   run goes on after `WARN: no <cache_dir>/info.txt; info.json carries no source versions`.
2. `cache_version` is the number that opens the path's last component (`.../116_GRCh38`), and
   `source_type` is `refseq`, `merged` or `ensembl` by the suffix of the component before it
   (`homo_sapiens_refseq/116_GRCh38`). A path without those components declares neither.
3. Every directory under `cache_dir` whose name does not start with `.` is a chromosome, in byte
   order. A chromosome the third argument names and the cache lacks stops the run.
4. `info.json`, `chr_synonyms.txt`, `sequence_regions.tsv` and `regulatory/binding_matrices.json`
   of an earlier run are removed, so that a run that stops early leaves none of them.
5. Every regulatory shard is read once for the names of its features' activity maps: their sorted
   union is `cell_types`.
6. Chromosome by chromosome: the transcript shards, the variation files, the regulatory shards,
   then one line on standard output. A transcript shard that cannot be converted is removed,
   reported as `WARN: <shard>: <reason>` followed by an empty line, and fails the run after the last
   chromosome; a variation file that cannot be copied, an `all_vars.gz` without an index and a
   regulatory shard that cannot be converted stop the run at once.
7. `chr_synonyms.txt`, `regulatory/binding_matrices.json` when a regulatory shard was written,
   `sequence_regions.tsv` and `info.json`.

A fatal error is one line `ERROR: [storable_to_json] <message>`, followed by an empty line for a
regulatory shard that cannot be converted, and exit code 1. A directory that
cannot be created ends the run with the operating system's error number as the exit code.

**Order and cores.** Shards are converted on the worker threads, one shard a thread, and written by
one thread in the order above, so the files and the lines of a run are the same at any thread
count. A worker starts no shard more than twice the thread count ahead of the writer: a whole shard
and its JSON are held in memory for each shard in flight, and never a whole cache. Standard output
is held in a buffer and written at exit unless it is a terminal; standard error is written line by
line. The two streams are each fixed; how they interleave on one terminal or file is not.

**A transcript.** A shard is a hash of lists of transcripts (one list a chromosome), or a list of
transcripts. Each transcript object is written with the fields `stable_id`, `version`, `dbID`,
`start`, `end`, `strand`, `biotype`, `source`, `is_canonical`, `cdna_coding_start`,
`cdna_coding_end`; the identifiers Ensembl VEP keeps under underscore-prefixed keys, unprefixed
(`gene_stable_id`, `gene_symbol`, `gene_symbol_source`, `gene_hgnc_id`, `gene_phenotype`, `ccds`,
`refseq`, `gene_version`, `protein`, `swissprot`, `trembl`, `uniparc`, `uniprot_isoform`,
`source_cache`, `gene_symbol_id`, `bam_edit_status`); `exons` (from `exons`, else
`_trans_exon_array`); `translation` with its `seq` as a string; `attributes`; and
`variation_effect_feature_cache` with the codon table, the UTR sequences as strings, the cDNA and
peptide strings, `seq_edits`, `introns`, `sorted_exons`, `protein_features` (the analysis as its
display label), the `mapper` and the prediction matrices. A field that is undef is left out; a list
keeps only its items that are references. A `biotype` that is a reference to a hash, the shape of a
`Bio::EnsEMBL::Biotype` object, is written as the string of its `name`, the value Ensembl's
`Transcript::biotype` returns for it, and is left out when the name is undef or a reference. The
mapper's `pairs` come from the first source that has any: the mapper's own `pairs` list, the pairs
of its `exon_coord_mapper` object in cDNA order, the exons in transcript order. A prediction matrix
is written in base64
(`MIME::Base64::encode_base64($blob, '')`); a matrix that does not start with the gzip magic is
gzipped first as `Compress::Zlib::memGzip` gzips it (zlib's deflate at level 9, memory level 9,
operating system byte 255), which is why the tool links stock zlib.

**The sequence side file.** `sequence_regions.tsv` opens with the line `name`, `length`, `circular`,
tab-separated as every row is, and holds one row for each chromosome directory whose transcripts
state a length, in the order of `contig_order`. The cache dumper leaves the length and, where it was
primed before the dump, the circular flag on each transcript's `slice`, and every transcript of a
chromosome shares one slice. The row is read from the first transcript of the chromosome's shards,
in their order, whose `slice` is a reference to a hash and whose `seq_region_length` matches Perl's
`/^[1-9]\d*$/`. The length is printed as Perl prints the scalar; the flag is `1` or `0` by Perl's
truth of the slice's `circular`, and `-` where that is undef. A chromosome without such a transcript
has no row, a cache without one has the header alone, and no transcript record carries the slice.
The line `Wrote <output_dir>/sequence_regions.tsv (<n> sequences)` follows the file on standard
output. The pattern is perl 5.34's to the byte, but for the one case the list below names: one
newline may end the length and is printed with it; under Perl's UTF-8 flag `\d` is any decimal
digit of Unicode 13.0, and a row with a digit past ASCII is written as UTF-8 after the line
`Wide character in print` on standard error; a length under the flag whose bytes are not UTF-8
fails its shard with `Malformed UTF-8 character (fatal)` when the step Perl takes back from the
end of the string, or from the newline that ends it, lands past any continuation bytes on a byte
that starts no character; the newline counts only where the string's first byte does not announce
a character as long as the whole string. A shard is read for a length only while no shard before
it in its chromosome has stated one.

**A regulatory shard** is a hash chromosome => feature class => list. Regulatory features come
before motif features. A regulatory feature is `feature_type` (its feature class,
`RegulatoryFeature`), `stable_id`, `start`, `end`, `dbid`, `biotype` where the cache gives one and
`activity`: one letter per name of `cell_types`, `A` `I` `P` `R` `N` for `ACTIVE` `INACTIVE`
`POISED` `REPRESSED` `NA` and `-` where the feature's map lacks the name. A motif feature is
`feature_type` (`MotifFeature`), `stable_id`, `start`, `end`, `dbid`, `strand`, `score`, `seq` and
`matrix`, the stable id of its
binding matrix; `regulatory/binding_matrices.json` holds each matrix once (`frequencies` by position
as `[A, C, G, T]`, `length`, `name`, `source`, `threshold`, `transcription_factors`). A shard of
another shape, a feature without a required key, an activity outside the five values, a motif
whose sequence does not fill its span, and a matrix that differs from an earlier copy under the
same id each stop the run with a message that names the shard. After the shard, the message for an
activity, for a span and for a missing key a motif alone carries names the feature by its stable id.
For a missing key every feature carries, it names the feature class and the key; for a matrix, its
id and the shard that holds the earlier copy; for a shape, the part of another shape (the image, a
chromosome entry or a feature list).

**JSON text.** Shards and matrices are written as `JSON::XS` writes them under
`JSON->new->utf8->canonical`: one line, no spaces, keys in byte order, `"` and `\` escaped,
`\b \t \n \f \r` by letter, other control characters as `\u00xx`, `/` and DEL as themselves, every
other character as UTF-8 (a byte above 127 of a string without Perl's UTF-8 flag is the character
of that number). `info.json` is `JSON->new->canonical->pretty`: three spaces a level, ` : ` after a
key, a final newline, and its bytes are Perl's `print` of that text.

**Number or string.** A value copied from the cache is written as the kind the cache stores it as,
which is the kind Perl had it as when the cache was dumped: a scalar stored as a string is a JSON
string even when it is all digits (`"start":"1"` in a protein feature), one stored as an integer is
a JSON integer, one stored as a double is written with C's `%.15g`, and Perl's immortal true and
false are the strings `"1"` and `""`. A copied reference is written as `JSON::XS` writes it under
`allow_blessed`: an unblessed hash or array as itself, a blessed value as `null`, a reference to
the scalar `1` or `0` as `true` or `false`. A value the conversion computes is a number: the
coordinates of a sequence edit, of the mapper object's pairs and of exon-derived pairs, every number
of a regulatory feature and of a binding matrix, `pair_count`, `length_diff`, `cache_version`. Its
text is Perl's for `<value> + 0`: digits alone are an integer; a number written with an exponent is
an integer when it is whole; a number with a decimal point is a double (`"465.0"` is written `465`,
`"4.40"` is `4.4`); text after a number is ignored and text without a number is 0. A value the
conversion puts in quotes is a string whatever the cache holds: a regulatory `stable_id`,
`biotype`, `seq` and `matrix`, a matrix's `name`, `source` and `transcription_factors`, a sequence
edit's `alt_seq`, the `name` of a biotype object.

**Where the tool does not write what Perl writes.** A reference run that meets one of these cases
holds the tool's bytes for it.

- A line that Perl ends with the place of a `die` or `warn` (` at <file> line <n>.`, also chained
  as `, at <file> line <n>`) is printed without it: `WARN: <shard>: Not a HASH reference`. This
  covers the messages of Perl itself, of Storable and of JSON::XS; a message with the conversion's
  own wording never had the place.
- `gzip -dc`'s own two lines for a file it cannot read (an empty line, then `gzip: <file>:
  <reason>`) are not printed. The reason the run reports is the same: `gzip -dc failed (wait status
  256)` for a damaged or short stream, `(wait status 512)` for bytes after the last member that are
  neither a member nor zeros, `(wait status 13)` for more than 64 KiB after the image inside the
  stream.
- A number of 10^15 or more in a computed field. Perl keeps a double beside an integer once a
  statement has computed with it as a double, and then writes whole numbers of that size with an
  exponent or as digits depending on what the same statement computed earlier in the run. The tool
  writes a whole double below 2^53 as an integer and any other double with `%.15g`. Below 10^15
  both texts are the same digits.
- A sequence length under the UTF-8 flag whose bytes are not UTF-8, after many lengths of a run.
  Perl looks at the end of a length before it matches it, and stops looking for the rest of the run
  once it has tried 101 lengths with the flag, or 101 without it, that hold a digit from 1 to 9. A
  length that fails its shard before that point states nothing after it, and a run of the same
  steps in Perl can end with exit code 0. The tool fails the shard wherever in the run the length
  stands.
- Perl's warnings for an undef or a non-numeric string in arithmetic are not printed, nor is its
  line `Malformed UTF-8 character: ... in pattern match (m//)` for a sequence length under the
  UTF-8 flag whose bytes are not UTF-8 and whose match does not die.
- A reference where the conversion needs a number or a string: Perl uses its address in memory,
  which no other run repeats; the tool fails the shard and says so. The message of JSON::XS for a
  reference to a scalar other than 0 or 1 is printed without the address it carries.
- A shard whose hash of chromosomes has several lists of transcripts, and a mapper whose
  `_pair_<from>` hash has several lists that tie on the cDNA start: Perl walks a hash in an order
  of its own that changes from run to run; the tool walks it in the order the image stores it.
- A weak reference whose target nothing else holds is undef to Perl after the image is read; the
  tool still reads the target.
- An image of a native byte order or size other than little-endian with 8-byte integers and
  doubles, an image of a later binary major and the boolean opcodes of binary minor 12 are refused
  with Storable 3.23's message for them; an image of binary major 1, tied values, hooks, code
  references and an opcode the format does not define are refused with the reader's own message,
  which names the opcode and its offset.

## Layout

- `src/storable/`: the reader of Perl Storable images (binary major 2, network and native order).
  It uses nothing else of the package: another tool compiles it with
  `#[path = "../../vep-storable-to-json/src/storable/mod.rs"] mod storable;` and a dependency on
  `flate2`. `mod.rs` lists what it decodes and what it refuses.
- `src/shard.rs`: one shard file to an image, with the outcome of `gzip -dc` and `fd_retrieve` for a
  file that is damaged.
- `src/perl.rs`: what Perl does with a retrieved value (`defined`, `ref`, `reftype`, truth, a fetch,
  a number, a string).
- `src/number.rs`: Perl's numbers and their text; `src/json.rs`: the two JSON encoders, base64 and
  Perl's `print`; `src/gzip.rs`: `memGzip`.
- `src/transcripts.rs`, `src/regulatory.rs`, `src/info.rs`, `src/regions.rs`: the four kinds of
  output.
- `src/pool.rs`: the worker threads and the ordered hand-over; `src/run.rs`: the run;
  `src/main.rs`: the command line.
- `tools/vep-storable-digest` and `tools/vep-storable-dump` compile `src/storable/`, `src/json.rs`,
  `src/number.rs`, `src/perl.rs`, `src/pool.rs` and `src/shard.rs` by `#[path]`, so a change to one
  of them runs the test, clippy and fmt commands of those two tools as well.
- `tests/storable_opcodes.rs` with `tests/fixtures/storable/`: the reader against images written by
  perl 5.32.1 with Storable 3.37, one fixture an opcode family, each named with the line of Perl
  that wrote it.
- `tests/crafted_images.rs`: the conversion of a transcript on images built in the test, where a
  copied field nests past the encoder's limit or holds itself, where the biotype is an object, and
  where the slice's length is damaged or is a double.
- `tests/sequence_regions.rs`: the rows of `sequence_regions.tsv` for caches of the fixtures, and
  its absence after a run that fails.
- `tests/golden.rs` with `tests/fixtures/caches/` and `tests/fixtures/expected/`: the binary against
  the reference bytes of 62 runs (perl 5.34.0, Storable 3.23, JSON::XS 4.04, zlib 1.2.11).
- `tests/ensembl_vep_data.rs`: the binary on Ensembl VEP's own data, for a machine that holds it.
  With `VEP_ENSEMBL_VEP_DIR` naming an ensembl-vep checkout it converts the release 84 test caches
  under `t/testdata/cache` (RefSeq, merged, Ensembl) and holds the pairs of every transcript to the
  Storable mapper object's; with `VEP_REGULATORY_SHARD` naming the release 116 GRCh38 shard
  `21/44000001-45000000_reg.gz` it holds the values that shard converts to. Each test passes
  without its variable.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest. The build
compiles zlib from the source the `libz-sys` crate carries and needs a C compiler.

```bash
cargo build --release --locked --manifest-path tools/vep-storable-to-json/Cargo.toml
cargo test --locked --manifest-path tools/vep-storable-to-json/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-storable-to-json/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-storable-to-json/Cargo.toml -- --check
```

The reference bytes are the contract, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the files it moves under `tests/fixtures/expected/<case>/out/`
with the binary's own output, so that the diff of the fixtures shows the change.
