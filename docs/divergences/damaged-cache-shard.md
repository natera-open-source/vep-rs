# A transcript or regulatory shard the cache cannot read

Three differences in which records are written, on cache states no measured run has; nothing is set aside. They are a transcript or regulatory file of the cache cut short after its header, a cache holding no transcript file at all, and a variation directory whose files are not the form `var_type` names. On each, Ensembl VEP reads no feature or known variant where the cache is damaged, writes its rows without them and exits 0; vep-rs ends the run with an error, exit 1.

Ensembl VEP's transcript cache holds one Storable file per 1 Mb region, gzip-compressed. A region's file is opened as a region is first asked for, through `PerlIO::gzip`, the `gzip` binary or `Compress::Zlib` (`AnnotationSource/Cache/BaseSerialized.pm` lines 134, 140 and 146), and a file that cannot be opened ends the run with the open's error. The bytes are then handed to Storable's `fd_retrieve` (lines 135, 141 and 155), which dies on a stream whose header is not Storable's, so a file of foreign bytes ends the run too.

A file cut after its header is read differently. `fd_retrieve` wraps the C retrieval in an `eval` and dies only on an exception (`Storable.pm` lines 435 to 449 of the Storable 3.23 in the release 116.3 image). A stream that ends inside the data returns `undef` and raises none, through each of the three decompressors. The caller stores whatever came back as the region's transcripts (`BaseSerialized.pm` lines 261 to 265), and the transcript loader reads `undef` as an empty table (`AnnotationSource/Cache/Transcript.pm` line 236). The region is recorded as holding no transcript, and the run exits 0. When the `gzip` binary is the decompressor its `unexpected end of file` appears on standard error; the other two decompressors report nothing.

The regulatory tree is read through the same code. A regulatory region file cut after its header returns `undef` likewise, which the regulatory loader iterates as no feature (`AnnotationSource/Cache/RegFeat.pm` line 209). Under `--regulatory` a variant whose regulatory and motif features that file alone held is written without those rows, its `intergenic_variant` row kept, exit 0 and nothing reported.

A variant whose overlapping transcripts the cut file alone held is written as `intergenic_variant`. A variant within 5,000 bases of the region's edge can still meet a transcript that spans the boundary. The regions read for a variant extend by that distance on each side (`AnnotationSource.pm` line 190, `AnnotationType/Transcript.pm` lines 190 to 196). The cache stores a transcript in every region file it overlaps, which is why the loader drops a transcript already seen (`AnnotationType/Transcript.pm` lines 282 to 287).

vep-rs loads every shard of a chromosome when the first record on it arrives. A shard whose bytes cannot be read, inflated or parsed ends the run before any record of that chromosome is annotated, with the chromosome and the loader's own words for the shard, exit 1. The rows of batches on other chromosomes already written stay in the output:

```
ERROR: [cache] Cannot read the transcript cache for chromosome 21:
Failed to inflate <cache>/transcripts/21/42000001-43000000.json.gz:
incomplete deflate stream
```

The line is wrapped here and is one line on standard error. A `c.`, `n.` or `p.` notation whose lookup reaches the shard first gives the same refusal after `Cannot read the transcript cache for HGVS input: `. A cache that carries its [checksum manifest](../cache-checksums.md) refuses the cut shard earlier, as a checksum mismatch, before any byte of it is inflated.

A regulatory shard is refused the same way, as the chromosome's regulatory shards load, with the loader's own words for the shard, exit 1, or as a checksum mismatch under a manifest:

```
ERROR: [cache] Failed to inflate
<cache>/regulatory/21/46000001-47000000.json.gz: incomplete deflate stream
```

The line is wrapped here and is one line on standard error; the output file holds its header and no row. A `regulatory/<chromosome>/` directory that holds files but no `<start>-<end>.json` or `.json.gz` shard is refused as the cache opens, `ERROR: [cache] <dir> holds no <start>-<end>.json or .json.gz shard`, exit 1; a file that is not a shard beside the shards is ignored.

Ensembl VEP reads a region file by its computed name, finds none in a directory without `_reg.gz` files (`BaseSerialized.pm` line 263) and writes every record there without regulatory rows, exit 0. An empty regulatory chromosome directory is a chromosome with no regulatory feature on both engines. No measured run reads a damaged regulatory shard (tier 3 of [Intended divergences from Ensembl VEP](../intended-divergences.md); 54 flag sets run `--regulatory`).

To reproduce in Ensembl VEP 116.3: copy the GRCh38 cache and cut one region file of chromosome 21 short with `head -c` to half its bytes. Run [the command of the closing section](reproducing.md#reproducing-these-examples) over a record more than 5,000 bases inside that region; the record is written as `intergenic_variant`. In vep-rs: copy a JSON cache without its `CHECKSUMS.sha256`, cut one shard of chromosome 21 the same way, and run any record on chromosome 21 over the copy.

## A cache holding no transcript shard

Ensembl VEP's valid chromosomes are the cache's directories whose names do not start with a dot, whatever they hold (`CacheDir.pm` line 431), and a region whose file does not exist yields no feature (`AnnotationSource/Cache/BaseSerialized.pm` line 263, `if -e $file`). So a cache none of whose chromosome directories holds a transcript file runs to the end: every record on a directory's name is written as `intergenic_variant`, exit 0, and nothing reports the cache empty.

vep-rs lists the shards of every chromosome directory as the cache opens, before the input is read, and refuses a cache that holds none, a dot-prefixed directory's files not counted:

```
ERROR: [cache] no transcript shards under <cache>/transcripts (serialiser_type
storable)
```

The line is one line on standard error, exit 1, and names the serialiser `info.json` declares (`storable` when it declares none). One directory holding no shard beside others that do is a chromosome with no transcripts in both engines, and a record on it is `intergenic_variant` in both.

To reproduce in Ensembl VEP 116.3: a cache directory holding `info.txt` and an empty `21/`; [the command of the closing section](reproducing.md#reproducing-these-examples) writes any record on chromosome 21 as `intergenic_variant`. In vep-rs: a cache directory holding `info.json` and an empty `transcripts/21/`.

## A variation directory whose files are not the form `var_type` names

Ensembl VEP's cache holds its known variants in one of two forms: one indexed `all_vars.gz` per chromosome, written by `convert_cache.pl`, which then appends `var_type tabix` to `info.txt` (`convert_cache.pl` line 327), or one `<start>-<end>_var.gz` per region, as the dumper writes them. Under `--check_existing` the key selects the reader (`CacheDir.pm` lines 184 to 186), and each reader opens the one file name it knows only if the file exists (`AnnotationSource/Cache/VariationTabix.pm` lines 156 and 270, `AnnotationSource/Cache/Variation.pm` line 104). A chromosome directory holding only the other form therefore yields no variant. Every record on it is written with `Existing_variation` empty and every co-location column empty, exit 0, and nothing reports the files unread.

vep-rs lists each chromosome's variation files as the variation source opens, before any record is annotated, and refuses a cache whose `var_type` and files disagree, in either direction:

```
ERROR: [cache] <cache>/variations/21 holds region files (<start>-<end>_var.gz)
while <cache>/info.json says var_type tabix
```

```
ERROR: [cache] <cache>/variations/21 holds all_vars.gz while <cache>/info.json
says var_type (unset)
```

Each is one line on standard error, exit 1, naming the directory, the form it holds and the key's value. The rule reads `all_vars.gz` alone: a directory holding it under no key, or lacking it under `var_type tabix`, is refused, and a directory holding neither form is a chromosome with no known variant in both engines. A directory holding both forms is read through the key's form by both engines under `var_type tabix`; under no key Ensembl VEP reads its region files and vep-rs refuses it for its `all_vars.gz`.

The check runs only when a run reads co-located variants: `--check_existing` or a flag that implies it (`--everything`, the frequency flags, `--pubmed`, `--check_frequency`, `--filter_common` and `--var_synonyms`), `--format id` input, or `--recode` with `id` or `var_synonyms`. A run without them never opens the variation files. No tool writes the state: `convert_cache.pl` appends the key after writing `all_vars.gz`, `vep-storable-to-json` and `vep-cache-converter` carry the key and the files across together, and `vep-cache-builder` writes both. It comes from a hand-edited metadata file or a partial copy.

The column list is checked first. An `all_vars.gz` record starts with its chromosome, so a `var_type tabix` cache whose `variation_cols` does not begin with `chr` is refused as the source opens, exit 1, whatever its directories hold:

```
ERROR: [cache] variation_cols of <cache>/info.json does not begin with chr, the
first column of an all_vars.gz record
```

Ensembl VEP reads each record's fields against the list as written (`AnnotationSource/Cache/BaseCacheVariation.pm` lines 119 to 124), so every value lands under the column name before its own: the name under `failed`, the start under `end`, the alleles under `strand`. No known variant then matches a record. Each record is written with `Existing_variation` empty, exit 0, and `isn't numeric` and `Use of uninitialized value` lines from `AnnotationType/Variation.pm` reach standard error, or under `--fork` the warnings file, where the parent writes each child's captured standard error (`Runner.pm` lines 565 and 658 to 662). `convert_cache.pl` prepends `chr` to the list as it appends the key (lines 315 to 318 and 327), so this state too is a hand edit.

The list must also name `start`, the column each known variant is placed by: a list naming other columns but not `start` is refused before the `chr` and file checks, exit 1. An empty list is read as no column list, so the files are not opened:

```
ERROR: [cache] variation_cols of <cache>/info.json names no start column
```

Ensembl VEP reads such a list as written (`CacheDir.pm` line 204), so every known variant has no position and matches no record: each record is written with `Existing_variation` empty, exit 0, and nothing is reported.

To reproduce in Ensembl VEP 116.3: copy the GRCh38 cache's `info.txt` and chromosome 21 directory, remove `21/all_vars.gz` and its `.csi`, and add `21/5000001-6000000_var.gz` holding one line of single-space-delimited cells: `rs_craft`, two empty cells, `5025532`, an empty cell and `G/C`. Run [the command of the closing section](reproducing.md#reproducing-these-examples) with `--check_existing` over `21 5025532 craft1 G C`: the record is written with `Existing_variation` `-`. Remove the `var_type` line and the leading `chr,` of `variation_cols`, and the same run prints `rs_craft`. For the column list: keep `21/all_vars.gz`, its `.csi` and `var_type tabix`, remove the leading `chr,` alone, and the run writes `Existing_variation` `-` with the warning lines above. Rename `start` in the list instead, and it writes `Existing_variation` `-` and nothing else.

In vep-rs: a JSON cache holding the same `variations/21/` under an `info.json` with `var_type` `tabix` and `chr` first in `variation_cols` is refused with the first line above. Without the key and with `all_vars.gz` in place, it is refused with the second; with the key, `all_vars.gz` and no `chr`, with the third; with `start` renamed, with the fourth.
