# vep-cache-converter

`vep-cache-converter` is a command-line binary, shipped in every release beside `vep` and `vep-cache-builder`, that prepares data files for the annotator. Without a subcommand it copies a flat JSON (JavaScript Object Notation) transcript cache into the per-chromosome layout `vep-cli` loads, validating each region file; `convert-plugin-data` converts a tabix-indexed plugin data file into the binary annotation store (`.vpd` data plus `.vpdi` index) that the `vep-builtins` plugins memory-map at runtime. No crate depends on it; it depends on `vep-builtins` for its record and store types.

## How it works

`convert-plugin-data` follows one data line from the input file to a verified store. `read_tabix_records` opens the input through a BGZF (blocked gzip) reader, keeps the most recent `#` header line as the column names, and parses each data line into a `TabixRecord`: the chromosome from `--chr-col`, the start from `--start-col` (plus one when `--zero-based` is set), the end from `--end-col` or equal to the start when that flag is absent or the column does not parse, the raw columns, and a shared pointer to the header. A line with fewer columns than `--chr-col` or `--start-col` needs is skipped; a start that is not an integer is an error. Records are grouped by chromosome in a `BTreeMap`, so chromosomes come out in string order (`1`, `10`, `11`, `2`), and each group is sorted by `(start, end)`. `BinaryStoreWriter::write` from `vep-builtins` then writes the column table and the records to the `--output` path and the chromosome and position index to the same path with the extension `.vpdi`. Unless `--no-verify` is set, `verify_concordance` reopens the store through `BinaryAnnotator`, the reader the plugins use at runtime, queries every record's start position when the file holds at most 1,000 records or `--verify-all` is set, otherwise a random sample of 1,000 record positions, and compares record count, column count and every column value against the parsed source; the first ten mismatches are printed and any mismatch is a non-zero exit.

The cache mode follows one region file. `convert_cache` creates the output directory, then `copy_info_json` reads `info.json`, parses it, logs its `species` and `assembly`, and writes the text unchanged; a missing or unparseable `info.json` is an error. `convert_data_files` runs once for `transcripts/` and once for `variations/`: it lists the directory in file-name order, takes each `.json` file, and `parse_region_filename` splits the stem at its last underscore into chromosome and range, so `LRG_485_1-1000000` yields the LRG (Locus Reference Genomic) sequence name `LRG_485` as its chromosome; a stem that does not fit is skipped with a warning. The content must parse as a JSON array (anything else is skipped with a warning); its length is counted and the text is written unchanged to `<output>/transcripts/<chr>/<start>-<end>.json`, the shard layout `vep-cli` enumerates. A missing `variations/` directory is a warning, not an error. `--region-size` is parsed and ignored: the bounds come from the file names.

## Modules

| Module    | Role                                                                                 |
| --------- | ------------------------------------------------------------------------------------ |
| `main.rs` | The whole crate: the clap definitions, both conversion paths and the unit tests. |

The functions of `main.rs`:

| Function                | Role                                                                                                        |
| ----------------------- | ----------------------------------------------------------------------------------------------------------- |
| `main`                  | Installs the `tracing` subscriber on stderr at INFO level, parses the command line and dispatches.          |
| `resolve_json_dir`      | Requires `--json-dir` and checks that the path exists.                                                      |
| `convert_plugin_data`   | Drives read, write and verification for one plugin data file; logs record counts, sizes and timings.        |
| `read_tabix_records`    | Streams the BGZF input into `TabixRecord`s grouped by chromosome and sorted by position.                    |
| `verify_concordance`    | Reopens the store with `BinaryAnnotator` and compares sampled or all positions against the source records.  |
| `convert_cache`         | Creates the output directory, then runs `copy_info_json` and `convert_data_files` for both subdirectories.  |
| `copy_info_json`        | Parses `info.json`, logs `species` and `assembly`, copies it unchanged.                                     |
| `convert_data_files`    | Re-shelves each `{chr}_{start}-{end}.json` region file as `<chr>/<start>-<end>.json`, counting entries.     |
| `parse_region_filename` | Splits a region file stem at its last underscore into `(chromosome, start, end)`.                           |

## Command line

`vep-cache-converter --help` and `vep-cache-converter convert-plugin-data --help` print the two flag sets. Diagnostics go to stderr through `tracing` at INFO level; the mismatch listing of a failed verification is the one direct write.

Cache mode (no subcommand):

| Flag                   | Role                                                                                                        |
| ---------------------- | ----------------------------------------------------------------------------------------------------------- |
| `--json-dir <JSON_DIR>` | The flat cache: `{chr}_{start}-{end}.json` region files directly under `transcripts/` and `variations/`.    |
| `--output-dir <OUTPUT_DIR>` | Where the per-chromosome layout and the copied `info.json` are written. Required.                           |
| `--region-size <REGION_SIZE>` | Accepted and ignored; region bounds come from the file names. Default `1000000`.                             |

`convert-plugin-data`:

| Flag                  | Role                                                                                                         |
| --------------------- | ------------------------------------------------------------------------------------------------------------ |
| `--input <INPUT>`     | The tabix-indexed TSV (tab-separated values) or VCF (Variant Call Format) file, `.tsv.gz` or `.vcf.gz`.      |
| `--output <OUTPUT>`   | The `.vpd` data file; the `.vpdi` index is written beside it.                                                |
| `--chr-col <CHR_COL>` | 0-based column index of the chromosome. Default `0`.                                                         |
| `--start-col <START_COL>` | 0-based column index of the start position. Default `1`.                                                     |
| `--end-col <END_COL>` | 0-based column index of the end position; omitted, the end equals the start.                                 |
| `--zero-based`        | Positions are 0-based, as in BED (Browser Extensible Data); one is added. Default: 1-based.                  |
| `--verify-all`        | Verify every record after writing, not a sample of 1,000.                                                    |
| `--no-verify`         | Skip verification.                                                                                           |

The `vep` flags that consume the outputs, `--json_cache` for the cache and `--plugin` for a store, are in [docs/cli-reference.md](../../docs/cli-reference.md).

## Run it

```bash
cargo build -p vep-cache-converter
cargo test -p vep-cache-converter
target/debug/vep-cache-converter --help
```

```bash
# Re-shelve a flat JSON cache into the layout `--json_cache` reads
vep-cache-converter --json-dir /path/to/flat_cache --output-dir /path/to/json_cache

# Build a binary store beside a plugin data file. A plugin reads the store when it
# is named after the data file with `.vpd` in place of `.gz` or `.bgz`.
vep-cache-converter convert-plugin-data \
  --input whole_genome_SNVs.tsv.gz \
  --output whole_genome_SNVs.tsv.vpd
```

The unit tests cover `parse_region_filename`. The workspace's standard build, test, lint and format commands are in [the contributing guide](../../CONTRIBUTING.md).

## Related

- [`vep-builtins`](../vep-builtins/README.md): `tabix::TabixRecord`, `binary_store::BinaryStoreWriter` and `binary_store::BinaryAnnotator`, the types the plugin-data path writes and reads back through; its plugins read the `.vpd` store.
- [`vep-cli`](../vep-cli/README.md): enumerates the `transcripts/<chr>/<start>-<end>.json` shards the cache mode writes. Declared in `Cargo.toml`, not imported by `main.rs`.
- [`vep-core`](../vep-core/README.md): declared in `Cargo.toml`, not imported by `main.rs`.
- [`vep-cache-builder`](../vep-cache-builder/README.md): writes the per-chromosome layout directly, so a natively built cache needs no conversion.
- [docs/cache-setup.md](../../docs/cache-setup.md): the routes to the per-chromosome JSON transcript cache the flat copy targets.
- [docs/plugins.md § Data File Preparation](../../docs/plugins.md#data-file-preparation): how a plugin finds a `.vpd` store beside its data file.
- [Root README § Installation](../../README.md#installation): the release archives and the container image that ship this binary.
