# Input and output

These flags name the input file and its format, and the output file and its format. The formats themselves are described under [VEP default and tab](../output/vep-and-tab.md) and [Parquet](../output/parquet.md).

## Input

| Flag | Default | What it does |
| --- | --- | --- |
| [`--input_file`, `-i`](#opt_input_file) | `STDIN` | Input file path |
| [`--input_data`](#opt_input_data) | none | No effect; input is read from `--input_file` |
| [`--format`](#opt_format) | `guess` | Input format: `vcf`, `ensembl`, `hgvs`, `region` or `guess` |

## Output

| Flag | Default | What it does |
| --- | --- | --- |
| [`--output_file`, `-o`](#opt_output_file) | `variant_effect_output.txt` | Output file path |
| [`--output_format`](#opt_output_format) | `vep` | Output format: `vep`, `vcf`, `json`, `tab` or `parquet` |
| [`--vcf`](#opt_vcf) | off | Shorthand for `--output_format vcf` |
| [`--json`](#opt_json) | off | Shorthand for `--output_format json` |
| [`--tab`](#opt_tab) | off | Shorthand for `--output_format tab` |
| [`--parquet_shape`](#opt_parquet_shape) | `nested` | Row shape of the Parquet output: `nested` or `flat` |
| [`--parquet_row_group_size`](#opt_parquet_row_group_size) | `122880` | Rows per Parquet row group; hidden |
| [`--keep_intermediate`](#opt_keep_intermediate) | off | Keeps the TSV intermediate of a Parquet run; hidden |
| [`--compress_output`](#opt_compress_output) | none | No effect; the output is written uncompressed |
| [`--force_overwrite`, `--force`](#opt_force_overwrite) | off | Overwrites an existing output file |
| [`--no_headers`](#opt_no_headers) | off | Suppresses the header lines |
| [`--mask_header_cache_path`](#opt_mask_header_cache_path) | off | Prints the cache directory as `[PATH]/<leaf>` in the headers |
| [`--fields`](#opt_fields) | none | No effect; the columns follow the field flags |
| [`--vcf_info_field`](#opt_vcf_info_field) | `CSQ` | VCF INFO field name for the annotations |

## Accepted for compatibility

| Flag | What Ensembl VEP does with it |
| --- | --- |
| `--input_data` | Reads the input from a string given on the command line |
| `--compress_output` | Writes the output through gzip or bgzip |
| `--fields` | Selects the output columns from a comma-separated list |

## Options

### `--input_file`, `-i` {#opt_input_file}

**`--input_file <path>`**, string, default `STDIN`. Ensembl VEP: same flag and default.

Input file path. `STDIN` reads standard input.

### `--input_data` {#opt_input_data}

**`--input_data <string>`**, string, no default. Ensembl VEP: accepted, no effect.

Not implemented: input is always read from `--input_file`.

### `--format` {#opt_format}

**`--format <name>`**, string, default `guess`. Ensembl VEP: same flag and default.

Input format: `vcf`, `ensembl`, `hgvs` or `region`. `guess` detects the format by file extension: `.hgvs` selects `hgvs`; anything else, including `STDIN`, is read as VCF. The formats are described under [Input formats](compatibility.md#input-formats).

### `--output_file`, `-o` {#opt_output_file}

**`--output_file <path>`**, string, default `variant_effect_output.txt`. Ensembl VEP: same flag and default.

Output file path. `STDOUT` writes to standard output.

### `--output_format` {#opt_output_format}

**`--output_format <name>`**, string, default `vep`. Ensembl VEP: same flag and default.

Output format: `vep`, `vcf`, `json`, `tab` or `parquet`; any other value writes the default `vep` format. For `parquet`, the `duckdb` CLI must be on `PATH`; see [Parquet](../output/parquet.md).

### `--vcf` {#opt_vcf}

**`--vcf`**, boolean, default off. Ensembl VEP: same flag and default.

Shorthand for `--output_format vcf`. As in Ensembl VEP, it also switches on `--symbol`, `--biotype` and `--numbers`.

### `--json` {#opt_json}

**`--json`**, boolean, default off. Ensembl VEP: same flag and default.

Shorthand for `--output_format json`. The JSON output is unescaped, as under Ensembl VEP's `--no_escape`; vep-rs has no `--no_escape` flag.

### `--tab` {#opt_tab}

**`--tab`**, boolean, default off. Ensembl VEP: same flag and default.

Shorthand for `--output_format tab`.

### `--parquet_shape` {#opt_parquet_shape}

**`--parquet_shape <shape>`**, string, default `nested`. Ensembl VEP: vep-rs only.

Row shape for `--output_format parquet`. `nested` writes one row per (record, allele): `Uploaded_variation`, `Location`, `Allele` and `Existing_variation` stay scalar, and every other `--tab` column is a LIST (VARCHAR elements; `DISTANCE` INTEGER, `STRAND` TINYINT). `flat` writes one row per (allele × transcript-consequence) tuple, every column scalar. Query examples for both shapes are under [Parquet](../output/parquet.md).

### `--parquet_row_group_size` {#opt_parquet_row_group_size}

**`--parquet_row_group_size <rows>`**, integer, default `122880`. Ensembl VEP: vep-rs only.

Hidden from `--help`. Overrides the Parquet row group size. The value is rounded up to a multiple of 2,048, DuckDB's vector size, which the dictionary limit follows, so that a row group with a distinct value per row keeps its Bloom filter. The flag exists so a small input can be written as several row groups under test.

### `--keep_intermediate` {#opt_keep_intermediate}

**`--keep_intermediate`**, boolean, default off. Ensembl VEP: vep-rs only.

Hidden from `--help`. Parquet output is written through a TSV intermediate, which is deleted on success; this flag keeps it.

### `--compress_output` {#opt_compress_output}

**`--compress_output <method>`**, string, no default. Ensembl VEP: accepted, no effect.

Not implemented: the output is written uncompressed.

### `--force_overwrite`, `--force` {#opt_force_overwrite}

**`--force_overwrite`**, boolean, default off. Ensembl VEP: same flag and default.

Overwrites the output file if it exists.

### `--no_headers` {#opt_no_headers}

**`--no_headers`**, boolean, default off. Ensembl VEP: same flag and default.

Suppresses the header lines of the output. Not for `parquet`, whose intermediate header is required.

### `--mask_header_cache_path` {#opt_mask_header_cache_path}

**`--mask_header_cache_path`**, boolean, default off. Ensembl VEP: same flag and default.

Prints the cache directory as `[PATH]/<leaf>` in the `## Using cache in` header line of the default and `--tab` formats, in the `cache=` token of the VCF `##VEP` line, and in the Parquet footer. This is Ensembl VEP 116's flag; vep-rs writes no statistics page.

### `--fields` {#opt_fields}

**`--fields <list>`**, string, no default. Ensembl VEP: accepted, no effect.

Not implemented: the columns follow the boolean field flags under [Output control](output-control.md).

### `--vcf_info_field` {#opt_vcf_info_field}

**`--vcf_info_field <name>`**, string, default `CSQ`. Ensembl VEP: same flag and default.

VCF INFO field name for the annotations.
