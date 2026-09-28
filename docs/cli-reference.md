# CLI Reference

The `vep` binary accepts the Perl VEP command-line interface. Most flags have the same name and behavior; a flag marked "Accepted for VEP compatibility" below is parsed and has no effect, so an existing VEP command line is accepted. The annotation cache is always `--json_cache`.

```bash
# Show all flags
cargo run --release -- --help

# Short form
vep -i input.vcf -o output.txt --offline --json_cache /path/to/json_cache
```

## Input

| Flag                 | Type   | Default | Description                                                                                                                                          |
| -------------------- | ------ | ------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input_file`, `-i` | String | `STDIN` | Input file path (use `STDIN` for standard input)                                                                                                     |
| `--input_data`       | String | --      | Accepted for VEP compatibility; not implemented (input is always read from `--input_file`)                                                           |
| `--format`           | String | `guess` | Input format: `vcf`, `ensembl`, `hgvs`, `region`. `guess` detects by file extension: `.hgvs` selects `hgvs`; anything else, including `STDIN`, is read as VCF |

## Output

| Flag                           | Type   | Default                     | Description                                                                                                                                                                                                                                                                   |
| ------------------------------ | ------ | --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--output_file`, `-o`          | String | `variant_effect_output.txt` | Output file path (use `STDOUT` for stdout)                                                                                                                                                                                                                                    |
| `--output_format`              | String | `vep`                       | Output format: `vep`, `vcf`, `json`, `tab`, `parquet`; any other value writes the default `vep` format. For `parquet`, the `duckdb` CLI must be on PATH (see Parquet section below).                                                                                          |
| `--vcf`                        | bool   | false                       | Shorthand for `--output_format vcf`                                                                                                                                                                                                                                           |
| `--json`                       | bool   | false                       | Shorthand for `--output_format json`                                                                                                                                                                                                                                          |
| `--tab`                        | bool   | false                       | Shorthand for `--output_format tab`                                                                                                                                                                                                                                           |
| `--parquet_shape`              | String | `nested`                    | Row shape for `--output_format=parquet`. `nested` = one row per (record, allele): `Uploaded_variation`, `Location`, `Allele` and `Existing_variation` stay scalar and every other `--tab` column is a LIST (VARCHAR elements; `DISTANCE` INTEGER, `STRAND` TINYINT). `flat` = one row per (allele × transcript-consequence) tuple, every column scalar. |
| `--compress_output`            | String | --                          | Accepted for VEP compatibility; not implemented (output is written uncompressed)                                                                                                                                                                                              |
| `--force_overwrite`, `--force` | bool   | false                       | Overwrite output file if it exists                                                                                                                                                                                                                                            |
| `--no_headers`                 | bool   | false                       | Suppress header lines in output (not for `parquet`, whose intermediate header is required)                                                                                                                                                                                    |
| `--fields`                     | String | --                          | Accepted for VEP compatibility; not implemented (the columns follow the boolean field flags below)                                                                                                                                                                            |
| `--vcf_info_field`             | String | `CSQ`                       | VCF INFO field name for annotations                                                                                                                                                                                                                                           |

### Output parity

For the golden-corpus flag set (`--offline --json_cache <dir> --species homo_sapiens
--assembly <assembly> --buffer_size 5000 --no_stats --quiet` plus one format flag),
the default, `--tab`, `--vcf` and `--json` data rows match Ensembl VEP 115.2 apart
from the consequence-term divergences the corpus manifests document; header lines
that carry a time, a path, the command line or Perl API component versions match
in shape. `--vcf` also switches on `--symbol`, `--biotype` and `--numbers`, as in
VEP, and `--json` output is unescaped (VEP's `--no_escape` behavior; vep-rs has no
`--no_escape` flag). Under `--everything`, vep-rs also switches on `--allele_number`
and `--hgvsg`, which VEP's `--everything` does not, so the `ALLELE_NUM` and `HGVSg`
header keys differ there. The formats, the header blocks and the golden-corpus
tests are described in [output-formats.md](output-formats.md).

### Parquet output

`--output_format=parquet` writes a TSV intermediate (one row per
(allele × transcript-consequence) tuple) and post-processes it to a
Hive-partitioned Parquet directory via a `duckdb` CLI subprocess. The
intermediate is deleted on success unless the hidden `--keep_intermediate` is
given. With `-o STDOUT` the TSV intermediate is streamed to stdout and no Parquet
is written.

**Requirements**: the `duckdb` CLI, 1.2 or newer, must be on `PATH`. Install via
`brew install duckdb` (macOS) or the CLI archive from
<https://duckdb.org/docs/installation/> (Linux: the `duckdb_cli-linux-*.zip`
release asset). vep-rs checks for it after the FASTA, plugins and variation
cache are loaded and before the first variant is read, and fails there if
`duckdb` is missing or too old, so no annotation is spent on a run it cannot
finalize.

**Schema**: five key columns (`chrom`, `pos`, `end`, `ref`, `alt`) then exactly the
`--tab` columns in `--tab` order; `pos`/`end` are BIGINT, `DISTANCE` INTEGER,
`STRAND` TINYINT, everything else VARCHAR with `-` stored as NULL, except the
allele columns (`ref`, `alt`, `Allele`), where `-` is a real allele (a deletion's
ALT, an insertion's REF) and is kept.

**Writer options** (all chosen for the reader): Parquet v2, ZSTD level 9, row
groups of 122,880 rows sorted by (chrom, pos, ref, alt), a dictionary and a
Bloom filter (0.1% false-positive rate) on every populated column of the `flat`
shape (asserted by the parity tests; a column NULL throughout a row group has
neither), the partition column stored inside the files, and the run's identity
(version, assembly, cache, command line, cache source versions) in the footer's
key-value metadata. The dictionary limit equals the row group size, which a
flat column's distinct values can never exceed; a LIST column of the `nested`
shape keeps its dictionary and Bloom filter only while its distinct elements in
a row group stay within that limit, so a column whose elements are mostly
distinct, such as `HGVSc`, carries neither. The hidden
`--parquet_row_group_size N` overrides the row group size (rounded up to a
multiple of 2,048, DuckDB's vector size, which the dictionary limit follows so
that a row group with a distinct value per row keeps its Bloom filter); it exists
so a small input can be written as several row groups under test.

**Shape** (`--parquet_shape`):

- `nested` (default): one row per (record, allele); the per-consequence columns
  are LIST columns aligned index-wise. Query with `UNNEST` or `list_contains`:

  ```sql
  SELECT chrom, pos, ref, alt
  FROM read_parquet('out.parquet/**/*.parquet', hive_partitioning=false)
  WHERE list_contains(Feature, 'ENST00000366667')
    AND list_contains(IMPACT, 'HIGH');
  ```

- `flat`: one row per (allele × transcript-consequence) tuple; every column
  scalar, so predicates push down directly:

  ```sql
  SELECT chrom, pos, ref, alt, Consequence, IMPACT
  FROM read_parquet('out.parquet/**/*.parquet', hive_partitioning=false)
  WHERE Feature = 'ENST00000366667' AND IMPACT = 'HIGH';
  ```

`hive_partitioning=false` keeps `chrom` as the stored VARCHAR; Hive type
inference would otherwise read `chrom=21` as a number.

## Cache / Offline

| Flag              | Type   | Default           | Description                                                                                                                                              |
| ----------------- | ------ | ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--cache`         | bool   | false             | Accepted for VEP compatibility; no effect (the annotation cache is `--json_cache`)                                                                       |
| `--offline`       | bool   | false             | Accepted for VEP compatibility; sets `--cache`, which has no effect. vep-rs always runs offline                                                          |
| `--dir`           | String | `~/.vep`          | Accepted for VEP compatibility; nothing is loaded from it. Heads the cache path printed in the output header when `--dir_cache` and `--json_cache` are absent |
| `--dir_cache`     | String | (same as `--dir`) | Accepted for VEP compatibility; nothing is loaded from it. Only names the cache path printed in the output header when `--json_cache` is absent          |
| `--cache_version` | u32    | 115               | Accepted for VEP compatibility; selects no cache. Only the version in the header's cache path when `--json_cache` is absent (defaults to the VEP release) |
| `--json_cache`    | String | --                | Path to the converted JSON cache directory: the only cache vep-rs reads. Without it no transcript is loaded and every variant is reported as `intergenic_variant` |
| `--fasta`, `--fa` | String | --                | Indexed FASTA file (`.fa` + `.fai`) for reference sequences                                                                                              |
| `--shift_3prime`  | bool   | false             | Enable 3' shifting of indels (requires `--fasta`)                                                                                                        |

The environment variable `VEP_FASTA_STORAGE` selects how the FASTA is held:
`mmap` (default) maps the file, `memory` reads it into RAM; any other value is an
error.

## Species / Assembly

| Flag               | Type   | Default        | Description                                |
| ------------------ | ------ | -------------- | ------------------------------------------ |
| `--species`, `-s`  | String | `homo_sapiens` | Accepted for VEP compatibility; only names the species in the header's cache path when `--json_cache` is absent. The data annotated is the `--json_cache` contents |
| `--assembly`, `-a` | String | --             | Genome assembly, `GRCh37` or `GRCh38` (`hg19`/`hg38` are accepted; an unrecognized value warns and leaves the plugin guards inert). Names the assembly in the JSON output's `assembly_name` and the Parquet footer metadata (without it the cache's own assembly is used there), and is **load-bearing for plugins**: it picks the coordinate column a dual-coordinate database is indexed on and refuses a data file whose own assembly contradicts it. Pass it whenever `--plugin` is used. |

## Runtime

| Flag                      | Type  | Default  | Description                                                                                                                                                                                                       |
| ------------------------- | ----- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--buffer_size`           | usize | 5000     | Number of variants per annotation batch                                                                                                                                                                           |
| `--max_sv_size`           | u64   | 10000000 | Structural variants spanning more than this many bases are carried through the VCF output without consequences and omitted from the JSON output (VEP's `--max_sv_size`); the default, `--tab` and `parquet` outputs list their rows       |
| `--fork`, `--threads`     | usize | 0        | Annotation worker threads. `0` = one per logical CPU; pass `1` for single-threaded annotation                                                                                                                      |
| `--decompression_threads` | usize | 0        | BGZF decompression threads for `.gz` VCF input. `0` = one per logical CPU, capped at 10; a positive value is used as given                                                                                        |
| `--distance`              | String | `5000`  | Upstream/downstream distance in bp (or `up,down`)                                                                                                                                                                 |
| `--transcript_index_impl` | String | `bin`   | Hidden. Transcript overlap index implementation: `bin`, `coitrees` or `sorted`. All three produce the same overlap sets; the flag exists for benchmarking                                                          |

## Output Control

These flags control which annotation fields appear in the output (`--coding_only` and `--no_intergenic` are accepted and have no effect). `--everything` switches on VEP's `--everything` set plus `--allele_number` and `--hgvsg`.

| Flag                 | Type   | Default | Description                                                                                                                                                                       |
| -------------------- | ------ | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--everything`, `-e` | bool   | false   | Switches on every flag in this table except `--xref_refseq`, `--coding_only` and `--no_intergenic`, every Allele Frequencies flag, `--mirna`, `--allele_number`, `--sift b` and `--polyphen b`. Implies `--regulatory`, so the not-implemented warning is logged on every `--everything` run |
| `--symbol`           | bool   | false   | Gene symbol (e.g., `BRCA1`)                                                                                                                                                       |
| `--biotype`          | bool   | false   | Transcript biotype                                                                                                                                                                |
| `--canonical`        | bool   | false   | Flag canonical transcripts                                                                                                                                                        |
| `--mane`             | bool   | false   | MANE Select and MANE Plus Clinical flags                                                                                                                                          |
| `--mane_select`      | bool   | false   | MANE Select flag only                                                                                                                                                             |
| `--tsl`              | bool   | false   | Transcript Support Level                                                                                                                                                          |
| `--appris`           | bool   | false   | APPRIS isoform annotation                                                                                                                                                         |
| `--ccds`             | bool   | false   | CCDS transcript ID                                                                                                                                                                |
| `--protein`          | bool   | false   | Ensembl protein ID                                                                                                                                                                |
| `--uniprot`          | bool   | false   | UniProt cross-references: `SWISSPROT` and `TREMBL`. The `UNIPARC` and `UNIPROT_ISOFORM` columns are added but never computed (always empty)                                        |
| `--xref_refseq`      | bool   | false   | RefSeq transcript cross-references                                                                                                                                                |
| `--hgvs`             | bool   | false   | HGVS coding and protein notations, and `HGVS_OFFSET` when an insertion or deletion was shifted 3' along the transcript for them                                                    |
| `--hgvsg`            | bool   | false   | Accepted for VEP compatibility; HGVSg is not computed (the `HGVSg` column is added but always empty)                                                                               |
| `--sift`             | String | --      | SIFT predictions (`p`=prediction, `s`=score, `b`=both)                                                                                                                            |
| `--polyphen`         | String | --      | PolyPhen predictions (`p`=prediction, `s`=score, `b`=both)                                                                                                                        |
| `--numbers`          | bool   | false   | Exon/intron numbers                                                                                                                                                               |
| `--domains`          | bool   | false   | Adds the `DOMAINS` column, which is never populated (see Compatibility with Perl VEP)                                                                                                                                                     |
| `--regulatory`       | bool   | false   | Accepted for VEP compatibility; regulatory-feature annotation is not implemented and the flag warns                                                                                |
| `--variant_class`    | bool   | false   | Variant class (SO term)                                                                                                                                                           |
| `--check_existing`   | bool   | false   | Check for co-located known variants                                                                                                                                               |
| `--coding_only`      | bool   | false   | Accepted for VEP compatibility; not implemented                                                                                                                                   |
| `--no_intergenic`    | bool   | false   | Accepted for VEP compatibility; not implemented                                                                                                                                   |

## Allele Frequencies

| Flag               | Type | Default | Description                                       |
| ------------------ | ---- | ------- | ------------------------------------------------- |
| `--af`             | bool | false   | Global allele frequency (1000 Genomes Phase 3)    |
| `--af_1kg`         | bool | false   | Continental frequencies (AFR, AMR, EAS, EUR, SAS) |
| `--af_gnomade`     | bool | false   | gnomAD exome frequencies                          |
| `--af_gnomadg`     | bool | false   | gnomAD genome frequencies                         |
| `--max_af`         | bool | false   | Highest frequency across all populations          |
| `--pubmed`         | bool | false   | PubMed IDs for co-located variants                |
| `--gene_phenotype` | bool | false   | Gene-phenotype associations                       |

## Filtering / Picking

| Flag                      | Type   | Default | Description                                                                                                                                                                                                     |
| ------------------------- | ------ | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--pick`                  | bool   | false   | Pick one consequence per variant                                                                                                                                                                                |
| `--pick_allele`           | bool   | false   | Pick one consequence per allele                                                                                                                                                                                 |
| `--per_gene`              | bool   | false   | Pick one consequence per gene                                                                                                                                                                                   |
| `--pick_allele_gene`      | bool   | false   | Pick one consequence per allele per gene                                                                                                                                                                        |
| `--most_severe`           | bool   | false   | Report only the single most severe consequence                                                                                                                                                                  |
| `--summary`               | bool   | false   | Report a summary of all consequences                                                                                                                                                                            |
| `--flag_pick`             | bool   | false   | Flag (but keep) non-picked consequences                                                                                                                                                                         |
| `--flag_pick_allele`      | bool   | false   | Flag non-picked per allele                                                                                                                                                                                      |
| `--flag_pick_allele_gene` | bool   | false   | Flag non-picked per allele per gene                                                                                                                                                                             |
| `--pick_order`            | String | --      | Accepted for VEP compatibility; not implemented. The pick order is fixed to VEP's default (`mane_select`, `mane_plus_clinical`, `canonical`, `appris`, `tsl`, `biotype`, `ccds`, `rank`); the final `length` tie-break is not applied |
| `--allele_number`         | bool   | false   | Add allele number from input                                                                                                                                                                                    |
| `--show_ref_allele`       | bool   | false   | Show reference allele in output                                                                                                                                                                                 |
| `--minimal`               | bool   | false   | Adds the `MINIMISED` column. Bi-allelic indels are always minimised, as VEP does regardless of the flag; per-allele minimisation of multi-allelic records is not implemented                                     |
| `--total_length`          | bool   | false   | Accepted for VEP compatibility; not implemented                                                                                                                                                                 |

## Plugins

| Flag            | Type                | Default | Description                                                                                                                                                                                                                                                                 |
| --------------- | ------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--plugin`      | String (repeatable) | --      | Load a plugin: `Name,param1,param2,...`. `Name` is a built-in (the names in [plugins.md](plugins.md), matched case-sensitively), the path of a C-ABI dylib, or a dylib name resolved as `lib<Name>.so` / `lib<Name>.dylib` under `--dir_plugins`. A name that is none of these is warned about and skipped (VEP's own default); a dylib that exists but fails to open, reports another ABI version or fails `init` aborts the run |
| `--custom`      | String (repeatable) | --      | Accepted for VEP compatibility; not implemented (custom annotation sources are not loaded)                                                                                                                                                                                  |
| `--dir_plugins` | String              | --      | Directory searched for `lib<Name>.so` / `lib<Name>.dylib` when `--plugin <Name>` is neither a built-in nor a path                                                                                                                                                           |

Dylib plugins implement ABI version 1 (`vep_plugin_abi_version`, `vep_plugin_create`,
`vep_plugin_init`, `vep_plugin_run`, `vep_plugin_free`, `vep_plugin_destroy`,
`vep_plugin_name`, `vep_plugin_version`). After the built-in plugins have run,
each is called once per transcript consequence with the consequence and the
variant as JSON and returns a JSON object whose entries become fields of the
default format's Extra column and of the `--json` output. The contract has no
header symbol, so a dylib field carries no `##` description line and no `--tab` or
`--vcf` column. A run-time failure in a dylib plugin is logged as an error and the
batch continues. `crates/vep-plugin/examples/echo_plugin.rs` is a working example
(built as a `cdylib`) and `crates/vep-cli/src/dylib_plugins.rs` is the host side.

## Transcript Sets

| Flag              | Type | Default | Description                                                                                       |
| ----------------- | ---- | ------- | ------------------------------------------------------------------------------------------------- |
| `--refseq`        | bool | false   | Accepted for VEP compatibility; not implemented (the transcript set is the `--json_cache` contents) |
| `--merged`        | bool | false   | Accepted for VEP compatibility; not implemented (the transcript set is the `--json_cache` contents) |
| `--gencode_basic` | bool | false   | Accepted for VEP compatibility; not implemented (the transcript set is the `--json_cache` contents) |

## Statistics

| Flag           | Type   | Default | Description                                                                                  |
| -------------- | ------ | ------- | -------------------------------------------------------------------------------------------- |
| `--no_stats`   | bool   | false   | Accepted for VEP compatibility; not implemented (no statistics file is written, with or without it) |
| `--stats_file` | String | --      | Accepted for VEP compatibility; not implemented (no statistics file is written)              |
| `--stats_text` | bool   | false   | Accepted for VEP compatibility; not implemented (no statistics file is written)              |

## Verbosity

| Flag              | Type | Default | Description                                                              |
| ----------------- | ---- | ------- | ------------------------------------------------------------------------ |
| `--verbose`, `-v` | bool | false   | Debug-level logging (to stderr)                                          |
| `--quiet`, `-q`   | bool | false   | Suppress status/progress messages and warnings; only errors are logged   |

## Miscellaneous

| Flag                  | Type | Default | Description                                                                                                                                                 |
| --------------------- | ---- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--dont_skip`         | bool | false   | Abort the run on a variant that fails parsing instead of skipping it (VEP keeps the variant); adds the `CHECK_REF` column, which is never computed (always empty) |
| `--allow_non_variant` | bool | false   | With `--vcf` output, write non-variant VCF records (ALT=`.`) through unchanged; ignored for the other output formats                                          |
| `--mirna`             | bool | false   | Accepted for VEP compatibility; miRNA secondary structure is not computed (the `miRNA` column is added but always empty)                                     |

## Compatibility with Perl VEP

What the flag tables above mark "Accepted for VEP compatibility" is the whole list of Perl VEP behaviour vep-rs does not implement; this section gathers it.

**Consequence terms never assigned.** Six regulatory-feature terms (`TFBS_ablation`, `TFBS_amplification`, `TF_binding_site_variant`, `regulatory_region_ablation`, `regulatory_region_amplification`, `regulatory_region_variant`) need a regulatory-feature cache that is not loaded: `--regulatory`, alone or through `--everything`, is accepted and warns that regulatory-feature annotation is not implemented. The root term `sequence_variant` has no predicate, and Perl VEP's consequence engine never assigns it either. The other 34 of Perl VEP's 41 Sequence Ontology terms are assigned.

**Input formats.** `--format` takes `vcf`, `ensembl` (the tab-delimited default format), `hgvs` and `region` (`chr:start-end/allele` and `chr:start-end:strand/allele`), or `guess`. HGVS input is genomic (`:g.`) only; `:c.`, `:n.`, `:p.` and `:r.` notation is rejected as not implemented. The `id`, `spdi` and `caid` formats are not read: Perl VEP resolves all three through its variation database and refuses them in `--offline` mode.

**No database mode.** vep-rs is offline-first and never contacts a database; `--database`, `--host`, `--port`, `--user` and `--password` do not exist and are rejected as unknown flags. The transcript set is the `--json_cache` contents ([cache-setup.md](cache-setup.md)), so `--refseq`, `--merged` and `--gencode_basic` select nothing.

**Accepted without effect.** `--custom` (custom annotation files are not processed); `--regulatory`; the statistics flags `--stats_file`, `--stats_text` and `--no_stats` (no statistics file is written with or without them); `--mirna` (the `miRNA` column is added and always empty); `--cache`, `--offline`, `--dir`, `--dir_cache` and `--cache_version` (they locate no cache; see the Cache table); `--coding_only` and `--pick_order`. Protein domains are never mapped: the `DOMAINS` field exists and protein features are read from the cache, but `--domains` leaves the column empty, and the native cache builder writes no domain features (the Ensembl GFF3 it reads carries none).

## Examples

### Basic offline annotation

```bash
vep -i input.vcf -o output.txt \
    --offline --cache \
    --json_cache /path/to/json_cache \
    --species homo_sapiens --assembly GRCh37
```

### VCF output with common fields

```bash
vep -i input.vcf -o annotated.vcf \
    --vcf --offline --cache \
    --json_cache /path/to/json_cache \
    --everything
```

### JSON output

```bash
vep -i input.vcf -o output.json \
    --json --offline --cache \
    --json_cache /path/to/json_cache \
    --everything
```

### Clinical plugins

```bash
vep -i input.vcf -o output.txt \
    --offline --cache --json_cache /path/to/json_cache \
    --assembly GRCh38 \
    --plugin CADD,snv=/path/to/cadd/whole_genome_SNVs.tsv.gz,indels=/path/to/cadd/gnomad.genomes.r4.0.indel.tsv.gz \
    --plugin REVEL,file=/path/to/revel/revel_scores.tsv.gz \
    --plugin SpliceAI,snv=/path/to/spliceai/spliceai_scores.raw.snv.vcf.gz \
    --everything --force
```

`--assembly` is not optional here. CADD and dbscSNV use it to reject a data file
built from the other genome; REVEL uses it to select which of its two coordinate
columns to read and rejects a copy whose tabix index was built on the other
column. Without it those guards are inert and REVEL falls back to the GRCh37
column with a warning, so a GRCh38-indexed file annotates nothing. SpliceAI also
accepts `indel=` and `cutoff=`; only the openly published GRCh38 SNV file is
validated, and the indel files, distributed through Illumina BaseSpace, are
unsupported data rather than a rejected parameter (see [plugins.md](plugins.md)).

### Parallel annotation

```bash
vep -i large_input.vcf -o output.txt \
    --offline --cache --json_cache /path/to/json_cache \
    --fork 8 --buffer_size 10000
```

### Filtered output (one consequence per allele per gene)

```bash
vep -i input.vcf -o output.txt \
    --offline --cache --json_cache /path/to/json_cache \
    --pick_allele_gene --canonical \
    --everything
```
