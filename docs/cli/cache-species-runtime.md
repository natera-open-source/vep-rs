# Cache, species and runtime

These flags name the annotation cache and the reference FASTA, the species and assembly, and the batch size, thread counts and distances of a run. Building or converting a cache is described under [Cache setup](../cache-setup.md).

## Cache and offline

| Flag | Default | What it does |
| --- | --- | --- |
| [`--cache`](#opt_cache) | off | No effect; the annotation cache is `--json_cache` |
| [`--offline`](#opt_offline) | off | No effect; vep-rs always runs offline |
| [`--dir`](#opt_dir) | `~/.vep` | Heads the cache path printed in the output header; nothing is loaded from it |
| [`--dir_cache`](#opt_dir_cache) | the value of `--dir` | Names the cache path printed in the output header; nothing is loaded from it |
| [`--cache_version`](#opt_cache_version) | `116` | Names the version in the header's cache path; selects no cache |
| [`--json_cache`](#opt_json_cache) | none | The converted JSON cache directory, the only cache vep-rs reads |
| [`--fasta`, `--fa`](#opt_fasta) | none | Indexed FASTA file for reference sequences |
| [`--shift_3prime`](#opt_shift_3prime) | off | Shifts indels 3'; requires `--fasta` |
| [`VEP_FASTA_STORAGE`](#env_vep_fasta_storage) | `mmap` | Environment variable: how the FASTA is held, `mmap` or `memory` |

## Species and assembly

| Flag | Default | What it does |
| --- | --- | --- |
| [`--species`, `-s`](#opt_species) | `homo_sapiens` | Names the species in the header's cache path; the data annotated is the `--json_cache` contents |
| [`--assembly`, `-a`](#opt_assembly) | none | Genome assembly, `GRCh37` or `GRCh38`; load-bearing for plugins |

## Runtime

| Flag | Default | What it does |
| --- | --- | --- |
| [`--buffer_size`](#opt_buffer_size) | `5000` | Variants per annotation batch |
| [`--max_sv_size`](#opt_max_sv_size) | 10000000 | Span above which a structural variant is oversize; `-1` lifts the limit |
| [`--fork`, `--threads`](#opt_fork) | `0` | Annotation worker threads; `0` is one per logical CPU, at most 32 |
| [`--decompression_threads`](#opt_decompression_threads) | `0` | BGZF decompression threads for `.gz` input; `0` is one per logical CPU, at most 10 |
| [`--distance`](#opt_distance) | `5000` | Upstream and downstream distance in bp, or `up,down` |
| [`--transcript_index_impl`](#opt_transcript_index_impl) | `bin` | Transcript overlap index implementation; hidden |

## Accepted for compatibility

The only cache vep-rs reads is `--json_cache`. The Ensembl VEP cache flags are accepted and only name the cache path the output header prints.

| Flag | What Ensembl VEP does with it |
| --- | --- |
| `--cache` | Reads annotation from a cache directory |
| `--offline` | Runs with no database connection, from the cache alone |
| `--dir` | Sets the base directory for the cache and plugins |
| `--dir_cache` | Sets the cache directory |
| `--cache_version` | Reads a cache of another version than the VEP release |
| `--species` | Selects the species annotated |

## Options

### `--cache` {#opt_cache}

**`--cache`**, boolean, default off. Ensembl VEP: accepted, no effect.

No effect: the annotation cache is `--json_cache`.

### `--offline` {#opt_offline}

**`--offline`**, boolean, default off. Ensembl VEP: accepted, no effect.

Sets `--cache`, which has no effect. vep-rs always runs offline.

### `--dir` {#opt_dir}

**`--dir <dir>`**, string, default `~/.vep`. Ensembl VEP: accepted, no effect.

Nothing is loaded from it. It heads the cache path printed in the output header when `--dir_cache` and `--json_cache` are absent.

### `--dir_cache` {#opt_dir_cache}

**`--dir_cache <dir>`**, string, default the value of `--dir`. Ensembl VEP: accepted, no effect.

Nothing is loaded from it. It only names the cache path printed in the output header when `--json_cache` is absent.

### `--cache_version` {#opt_cache_version}

**`--cache_version <n>`**, integer, default `116`. Ensembl VEP: accepted, no effect.

Selects no cache. It is only the version in the header's cache path when `--json_cache` is absent; the default is the Ensembl VEP release.

### `--json_cache` {#opt_json_cache}

**`--json_cache <dir>`**, string, no default. Ensembl VEP: vep-rs only.

Path to the converted JSON cache directory, the only cache vep-rs reads. Without it no transcript is loaded and every variant is reported as `intergenic_variant`. [Cache setup](../cache-setup.md) describes how to build or convert one.

### `--fasta`, `--fa` {#opt_fasta}

**`--fasta <file>`**, string, no default. Ensembl VEP: same flag and default.

Indexed FASTA file (`.fa` plus `.fai`) for reference sequences.

### `--shift_3prime` {#opt_shift_3prime}

**`--shift_3prime`**, boolean, default off. Ensembl VEP: same flag and default.

Enables 3' shifting of indels. Requires `--fasta`.

### `VEP_FASTA_STORAGE` {#env_vep_fasta_storage}

**`VEP_FASTA_STORAGE=<mode>`**, environment variable, default `mmap`. Ensembl VEP: vep-rs only.

Selects how the FASTA is held: `mmap` maps the file, `memory` reads it into RAM. Any other value is an error.

### `--species`, `-s` {#opt_species}

**`--species <name>`**, string, default `homo_sapiens`. Ensembl VEP: accepted, no effect.

Only names the species in the header's cache path when `--json_cache` is absent. The data annotated is the `--json_cache` contents.

### `--assembly`, `-a` {#opt_assembly}

**`--assembly <name>`**, string, no default. Ensembl VEP: same flag and default.

Genome assembly, `GRCh37` or `GRCh38`; `hg19` and `hg38` are accepted. An unrecognized value warns and leaves the plugin guards inert.

The flag names the assembly in the JSON output's `assembly_name` and in the Parquet footer metadata; without it the cache's own assembly is used there.

It is load-bearing for plugins: it picks the coordinate column a dual-coordinate database is indexed on and refuses a data file whose own assembly contradicts it. Pass it whenever `--plugin` is used.

When `--json_cache` is absent the flag also supplies the `_<assembly>` suffix of the header's cache path leaf, `<version>_GRCh38`. The leaf is the bare `<version>` when the flag is absent.

### `--buffer_size` {#opt_buffer_size}

**`--buffer_size <n>`**, integer, default `5000`. Ensembl VEP: same flag and default.

Number of variants per annotation batch.

### `--max_sv_size` {#opt_max_sv_size}

**`--max_sv_size <bases>`**, integer, default 10000000. Ensembl VEP: same flag and default.

A structural variant spanning more than this many bases is oversize; `-1` lifts the limit, as in VEP 116. A `<CNV:TR>` record whose repeat sequence is longer than this stays symbolic rather than reading a direction of change (under `-1` every one does), as in VEP 116.

Ensembl VEP 116 drops an oversize record before annotation and writes nothing for it in any format. vep-rs annotates it: rows in the default, `--tab` and `parquet` outputs, its VCF line without consequences, no JSON object. See [how structural variants are written](../output/structural-variants.md) and [class 5 on the Defect C page](../divergences/defect-c.md#5-a-record-ensembl-vep-drops-before-annotation).

### `--fork`, `--threads` {#opt_fork}

**`--fork <n>`**, alias `--threads`, integer, default `0`. Ensembl VEP: same flag, vep-rs default differs.

Annotation worker threads. `0` is one per logical CPU, at most 32. Pass `1` for single-threaded annotation, or any positive value to pin the count; a given value is not capped.

### `--decompression_threads` {#opt_decompression_threads}

**`--decompression_threads <n>`**, integer, default `0`. Ensembl VEP: vep-rs only.

BGZF decompression threads for `.gz` VCF input. `0` is one per logical CPU, capped at 10; a positive value is used as given.

### `--distance` {#opt_distance}

**`--distance <bp>`**, string, default `5000`. Ensembl VEP: same flag and default.

Upstream and downstream distance in base pairs, or `up,down` for two values.

### `--transcript_index_impl` {#opt_transcript_index_impl}

**`--transcript_index_impl <name>`**, string, default `bin`. Ensembl VEP: vep-rs only.

Hidden from `--help`. Transcript overlap index implementation: `bin`, `coitrees` or `sorted`. All three produce the same overlap sets; the flag exists for benchmarking.
