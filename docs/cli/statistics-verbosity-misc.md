# Statistics, verbosity and miscellaneous

These flags cover the statistics file Ensembl VEP writes and vep-rs does not, the logging level, and three flags that change how a record is handled.

## Statistics

| Flag | Default | What it does |
| --- | --- | --- |
| [`--no_stats`](#opt_no_stats) | off | No effect; no statistics file is written with or without it |
| [`--stats_file`](#opt_stats_file) | none | No effect; no statistics file is written |
| [`--stats_text`](#opt_stats_text) | off | No effect; no statistics file is written |

## Verbosity

| Flag | Default | What it does |
| --- | --- | --- |
| [`--verbose`, `-v`](#opt_verbose) | off | Debug-level logging to stderr |
| [`--quiet`, `-q`](#opt_quiet) | off | Logs errors only |

## Miscellaneous

| Flag | Default | What it does |
| --- | --- | --- |
| [`--dont_skip`](#opt_dont_skip) | off | Aborts the run on a variant that fails parsing |
| [`--allow_non_variant`](#opt_allow_non_variant) | off | Writes non-variant VCF records through unchanged in `--vcf` output |
| [`--mirna`](#opt_mirna) | off | Adds the `miRNA` column, which stays empty |

## Accepted for compatibility

| Flag | What Ensembl VEP does with it |
| --- | --- |
| `--no_stats` | Writes no statistics file |
| `--stats_file` | Names the statistics file |
| `--stats_text` | Writes the statistics as plain text instead of HTML |
| `--mirna` | Reports where in a miRNA's secondary structure the variant falls |

## Options

### `--no_stats` {#opt_no_stats}

**`--no_stats`**, boolean, default off. Ensembl VEP: accepted, no effect.

Not implemented: no statistics file is written, with or without it.

### `--stats_file` {#opt_stats_file}

**`--stats_file <path>`**, string, no default. Ensembl VEP: accepted, no effect.

Not implemented: no statistics file is written.

### `--stats_text` {#opt_stats_text}

**`--stats_text`**, boolean, default off. Ensembl VEP: accepted, no effect.

Not implemented: no statistics file is written.

### `--verbose`, `-v` {#opt_verbose}

**`--verbose`**, boolean, default off. Ensembl VEP: same flag and default.

Debug-level logging, to stderr.

### `--quiet`, `-q` {#opt_quiet}

**`--quiet`**, boolean, default off. Ensembl VEP: same flag and default.

Suppresses status and progress messages and warnings; only errors are logged.

### `--dont_skip` {#opt_dont_skip}

**`--dont_skip`**, boolean, default off. Ensembl VEP: same flag and default.

Aborts the run on a variant that fails parsing instead of skipping it; Ensembl VEP keeps the variant. Adds the `CHECK_REF` column, which is never computed and stays empty.

### `--allow_non_variant` {#opt_allow_non_variant}

**`--allow_non_variant`**, boolean, default off. Ensembl VEP: same flag and default.

With `--vcf` output, writes non-variant VCF records (ALT `.`) through unchanged. Ignored for the other output formats.

### `--mirna` {#opt_mirna}

**`--mirna`**, boolean, default off. Ensembl VEP: accepted, no effect.

The miRNA secondary structure is not computed: the `miRNA` column is added but always empty.
