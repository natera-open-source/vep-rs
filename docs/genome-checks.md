# Genome checks

The files that record, for each genome with a ready-made Ensembl VEP cache, whether vep-rs annotates a small generated input as Ensembl VEP does.

## The genome table

`tests/genomes/` holds four tab-separated files, each under one header line:

- `genomes.tsv`: one row a ready-made cache of the Ensembl and Ensembl Genomes release directories and of the human pangenome
- `bacteria.tsv.gz`: the same columns, gzip-compressed, one row a bacterial or archaeal cache
- `outcomes.tsv`: `outcome` and `meaning`, every `differs:` and `excluded:` value a row uses
- `sets.tsv`: `set`, `release`, `cache_version` and `caches`, how many ready-made caches each set holds

A row has nine columns, in this order.

| Column | Holds |
| --- | --- |
| `genome` | the cache's species name, without a `_refseq` or `_merged` suffix |
| `assembly` | the cache's assembly name |
| `transcripts` | `ensembl`, `refseq` or `merged` |
| `set` | `vertebrates`, `metazoa`, `plants`, `fungi`, `protists`, `pangenome`, `bacteria_<N>_collection`, `fungi_<name>_collection` or `protists_<name>_collection` |
| `cache_version` | the number Ensembl VEP's `--cache_version` takes |
| `seed` | the seed the genome's input records were generated from |
| `records` | how many input records were generated |
| `outcome` | `exact`, `differs:<name>` or `excluded:<name>` |
| `sha256` | the digest of vep-rs's output, or `-` where there is none |

### The genome key

A genome's key is the one string `<set>__<genome>__<assembly>__<transcripts>`. It takes all four fields: a species ships under up to three transcript sets, and a few genomes ship in two sets. The key is built from a row and never parsed back. No two rows share a key, and the rows of a file stand in bytewise key order. A row whose set starts with `bacteria_` belongs to `bacteria.tsv.gz`; every other row belongs to `genomes.tsv`.

### Outcomes and the digest

- `exact`: on the genome's records vep-rs matches Ensembl VEP exactly except for the documented divergences. The row carries a digest and a record count above zero.
- `differs:<name>`: the two outputs differ. `differs:unclassified` is a difference that has no name.
- `excluded:<name>`: the genome was not compared, for the reason `outcomes.tsv` gives.

A name is lower-case letters, digits and hyphens.

The digest is the SHA-256, as 64 lower-case hex digits, of vep-rs's default-format output with every line that starts with `##` left out, since those lines carry paths and times.

## The table tool

`tools/vep-genome-table` reads, checks and writes the table. Its [README](../tools/vep-genome-table/README.md) gives each subcommand's rules and exit codes.

```text
vep-genome-table check <dir> [--complete]
vep-genome-table row <dir> <genome key> [--field <column>]
vep-genome-table keys <dir> [--set <set>] [--outcome <prefix>]
vep-genome-table put <dir> <row file>...
vep-genome-table digest <output file>
```

`check` holds the four files to the rules above and, with `--complete`, each set to the `caches` count `sets.tsv` gives it. `put` replaces rows by key and rewrites both row files in key order, so equal rows give equal bytes. `digest` prints the digest of an output file, plain or gzip. Every digest in the table comes from it.

## The per-genome job

`scripts/genomes/check_genome.sh` checks one genome and writes its row. It takes:

- one cache: `--cache-tarball <file>`, `--cache-dir <extracted Ensembl cache directory>` or `--json-cache <converted cache directory>`
- the row's identity: `--genome`, `--assembly`, `--transcripts`, `--set` and `--cache-version`
- `--seed` and `--mode parity|stored|recheck`
- `--vep-image`, the Ensembl VEP image reference
- `--vep-rs`, `--generator`, `--comparator` and `--table-tool`: the paths of the `vep`, `vep-seeded-variants`, `vep-outputs-compare` and `vep-genome-table` binaries
- `--reference-dir` in `stored` mode and `--table-dir` in `recheck` mode
- `--stop-after extract|convert|generate`, which ends the run after that stage
- `--work` and `--out`, a work directory and an output directory

A converted cache may accompany `--cache-dir`, and the conversion is then not run. The conversion and Ensembl VEP run in the image under docker. The converter is the release build of `tools/vep-storable-to-json`, mounted read-only into the image: the build under `CARGO_TARGET_DIR` when that holds one, else the tool's own. The job builds nothing and stops, before it writes anything, when it finds no build. Both engines run offline under default flags in the default output format. The job reads and writes local paths alone, and `--help` prints every option.

A cache whose transcripts were edited from an alignment file names that file under `bam` in its `info.txt`, and the converted cache's `info.json` keeps the line. On such a cache Ensembl VEP and vep-rs both switch `--use_transcript_ref` on and then refuse a run that names no FASTA file. The job names none, so it gives both engines `--use_given_ref` on such a cache and on no other. Ensembl's human RefSeq and merged caches are edited so.

`parity` runs both engines and compares their outputs. `stored` runs vep-rs alone and compares its output with an Ensembl VEP output kept from a `parity` run. `recheck` runs vep-rs alone and compares its digest with the genome's row in the table.

A job that holds the Ensembl cache, not a converted cache alone, names it to the comparator, which then sets aside Ensembl VEP's rows of [a record under a name the synonyms file gives to several sequences](divergences/synonym-of-several-sequences-by-hash-order.md). On a cache whose synonyms file gives the name of the generator's `synonym` record to several sequences, that record is one vep-rs skips.

The job leaves these files in the output directory:

- `row.tsv`: the table's header line and the genome's row
- `variants.vcf`: the input records, and `variants.summary.tsv`: the generator's count of records a category
- `ensembl_vep.txt.gz` and `ensembl_vep.warnings.txt`: Ensembl VEP's output and its warnings, in `parity` mode
- `comparison/`: the report `vep-outputs-compare` wrote
- `job.json`: each step's seconds and peak memory, the tool versions, the image digest and the digest of the input records

`job.json` is diagnostics, and nothing in it is a table column. The job writes `exact`, `differs:unclassified` or an `excluded:` value, never another `differs:` name. The outputs are `exact` when the comparison leaves no row open. For a genome checked on its RefSeq or its merged cache the comparison names the converted cache. It so sets aside the rows of the documented difference [`merged_same_stable_id_collapse`](cache-conversion.md#a-transcript-model-the-cache-carries-at-two-loci-under---refseq-or---merged) and, on a merged cache, of [`merged_same_dbid_collapse`](divergences/merged-same-dbid-collapse.md).

The comparison names the converted cache for a genome checked on its Ensembl transcripts as well. A documented difference at a coding end, such as [`start_codon_overlap_without_start_codon`](divergences/start-codon-overlap-without-start-codon.md), is set aside only on the reading of the transcript's record there.

| Exit code | Meaning |
| --- | --- |
| 0 | `exact`, a `recheck` match, or the `--stop-after` stage reached |
| 1 | the outputs differ, or a `recheck` mismatch |
| 2 | a usage error, or a skip: no build of the converter, no docker, no image, or a `recheck` of a row that is not `exact` |
| 3 | a step failed; the row is written when the failure is an `excluded:` outcome |

## Checking one genome

`scripts/genomes/test_genome.sh` checks genomes of the table against their rows. For each genome it reads the row, takes the genome's cache archive, and runs the per-genome job on it with the row's seed.

```bash
bash scripts/genomes/test_genome.sh --mode recheck \
    --key vertebrates__saccharomyces_cerevisiae__R64-1-1__ensembl
bash scripts/genomes/test_genome.sh --mode parity \
    --key vertebrates__saccharomyces_cerevisiae__R64-1-1__ensembl
```

- `recheck`, the default, runs vep-rs alone. The genome passes when the digest and the record count equal the row's.
- `parity` runs Ensembl VEP and vep-rs. The genome passes when the job's outcome equals the row's. The job gives a difference no name, so a `differs:` row passes while the two outputs differ.

The genomes are named by `--key`, once a genome, by `--keys-file`, one key a line, or by `--set <set> --sample <n>`. `--seed` picks the sample, and the same seed picks the same keys. Each genome prints one line on standard output.

| Line | Meaning |
| --- | --- |
| `ok <key> <mode>` | the genome ran and gave what its row records |
| `FAIL <key> <mode>: <what>` | the run differs from the row, or the archive or a step failed; `<what>` gives both values |
| `skip <key>: <why>` | the genome did not run: its row is `excluded:<name>`, its archive is over `--max-bytes`, or the host lacks docker, the image, or a build of the cache converter |

The script exits 0 when every genome that ran is `ok`, 1 on a `FAIL`, and 2 when no genome ran. It takes the release builds of `vep` and of the three tools from the checkout, and `--help` lists the options that name others.

The conversion runs the release build of `tools/vep-storable-to-json`, which the per-genome job finds by the rule given above. `recheck` converts in the Ensembl VEP image when docker holds it and on the host otherwise. `parity` needs the image, and so does a run that names one with `--vep-image`. The script builds nothing: without a build it skips every genome on a line that gives the build command.

```bash
cargo build --release --locked \
    --manifest-path tools/vep-storable-to-json/Cargo.toml
```

### Where the cache comes from

The archive is downloaded from Ensembl's public site. Its directory follows from the row's set and from the set's release in `sets.tsv`.

| Set | Directory |
| --- | --- |
| `vertebrates` | `pub/release-<release>/variation/indexed_vep_cache/` on `ftp.ensembl.org` |
| `metazoa`, `plants`, `fungi`, `protists` | `pub/<set>/release-<release>/variation/indexed_vep_cache/` on `ftp.ensemblgenomes.ebi.ac.uk` |
| a collection | the folder of the set's name under the directory of `bacteria`, `fungi` or `protists` |
| `pangenome` | `pub/rapid-release/species/Homo_sapiens/<accession>/ensembl/variation/<release>/indexed_vep_cache/` on `ftp.ensembl.org` |

A download is held to the length the site states and to the directory's `CHECKSUMS` file. An archive that file does not list is held to gzip's own check. So is an archive whose directory has no such file, as every pangenome directory. A `CHECKSUMS` file the site fails to serve fails the genome. `--cache-root <dir>` takes the archives from a local directory, under the names Ensembl ships them, and downloads nothing.

### The workflow

The workflow `Genome check`, in `.github/workflows/genome-check.yml`, runs the script on a hosted runner. It starts by hand alone, from the repository's Actions page: choose `Genome check`, then `Run workflow`. Its inputs are `keys` (comma-separated genome keys), `set` and `sample` (used when `keys` is empty), `mode` and `max_bytes`. A genome whose archive is larger than `max_bytes` is skipped. Each genome's row, job record and comparison report are kept as the run's artifact.

## The sequence side file

A converted cache may carry `sequence_regions.tsv` beside `info.json`. The file is optional and tab-separated under a header line, one row a sequence. Its columns are `name` (the cache's own name of the sequence), `length` (in bases) and `circular` (`1`, `0`, or `-` where the Ensembl cache states neither).

The converted transcripts carry neither value. The cache converter writes the file from the slice each sequence's transcripts carry, so a sequence without a transcript has no row. The record generator reads the file to place records at the ends of a sequence and across the origin of a circular one. The annotator does not read it.

## The harness's genomes

`scripts/concordance/harness_genomes.tsv` holds the concordance harness's genomes as data, one row a genome under a header line. Its columns are `species`, `assembly`, `cache_version` (the number Ensembl VEP's `--cache_version` takes), `fasta` and `json_cache`. The last two are paths under the harness's data root.

`scripts/concordance/harness_genomes.sh` is a file to source. Its one function prints one field of a genome's row.

```text
harness_genome_field <table> <species> <assembly> <column>
```

The last row for a genome wins. A genome or a column the table lacks prints one `ERROR: [harness_genomes] <what>` line on standard error and returns 1. So does a row whose field count is not the header's. The function runs under bash 3.2.
