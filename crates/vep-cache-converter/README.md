# vep-cache-converter

`vep-cache-converter` is a command-line binary, shipped in every release beside `vep` and `vep-cache-builder`, that prepares data files for the annotator. Without a subcommand it copies a flat JSON (JavaScript Object Notation) transcript cache into the per-chromosome layout `vep-cli` loads, validating each region file; `convert-plugin-data` converts a tabix-indexed plugin data file into the binary annotation store (`.vpd` data plus `.vpdi` index) that the `vep-builtins` plugins memory-map at runtime. No crate depends on it; it depends on `vep-builtins` for its record and store types and on `vep-cli` for the regulatory shard parser.

## How it works

`convert-plugin-data` follows one data line from the input file to a verified store. `read_tabix_records` opens the input through a BGZF (blocked gzip) reader, keeps the most recent `#` header line as the column names, and parses each data line into a `TabixRecord`: the chromosome from `--chr-col`, the start from `--start-col` (plus one when `--zero-based` is set), the end from `--end-col` or equal to the start when that flag is absent or the column does not parse, the raw columns, and a shared pointer to the header. A line with fewer columns than `--chr-col` or `--start-col` needs is skipped; a start that is not an integer is an error. Records are grouped by chromosome in a `BTreeMap`, so chromosomes come out in string order (`1`, `10`, `11`, `2`), and each group is sorted by `(start, end)`. `BinaryStoreWriter::write` from `vep-builtins` then writes the column table and the records to the `--output` path and the chromosome and position index to the same path with the extension `.vpdi`. Unless `--no-verify` is set, `verify_concordance` reopens the store through `BinaryAnnotator`, the reader the plugins use at runtime, queries every record's start position when the file holds at most 1,000 records or `--verify-all` is set, otherwise a random sample of 1,000 record positions, and compares record count, column count and every column value against the parsed source; the first ten mismatches are printed and any mismatch is a non-zero exit.

The cache mode follows one region file. `convert_cache` creates the output directory, then `copy_info_json` reads `info.json`, parses it, checks what the runtime's cache resolution reads (an object with a non-empty `species` and `assembly`, and a `contig_order` or `variation_cols` that is a list of strings when present), logs its `species` and `assembly`, writes it with the `regulatory` marker set when `copy_regulatory_files` wrote a shard and removed otherwise, and copies `chr_synonyms.txt` beside it when the source has one; a missing, unparseable or unusable `info.json` is an error. `convert_data_files` runs over `transcripts/`: it lists the directory in file-name order, takes each `.json` file, and `parse_region_filename` splits the stem at its last underscore into chromosome and range, so `LRG_485_1-1000000` yields the LRG (Locus Reference Genomic) sequence name `LRG_485` as its chromosome; a stem that does not fit is skipped with a warning. The content must parse as JSON (a parse failure is an error) and be an array (anything else is skipped with a warning); its length is counted and the text is written unchanged to `<output>/transcripts/<chr>/<start>-<end>.json`, the shard layout `vep-cli` enumerates. `copy_variation_files` runs over `variations/`, which holds the files the Ensembl VEP dumper writes: `all_vars.gz` with its `.csi` or `.tbi` index, or `{start}-{end}_var.gz` region files, under `variations/<chr>/`, or carrying a `{chr}_` prefix directly under `variations/`. `split_flat_variation_name` finds the chromosome in a flat name by trying each underscore until what follows is a variation file name, and each file is copied byte for byte to `<output>/variations/<chr>/`; any other name, flat or under `<chr>/`, is skipped with a warning. An `all_vars.gz` without its index is an error, as is a `.json` file under `variations/`: the annotator reads the dumper's files only. A missing `variations/` directory is a warning, not an error. `copy_regulatory_files` runs over `regulatory/<chr>/`, copying each `<start>-<end>.json` or `.json.gz` shard byte for byte once it parses as the annotator parses a shard (`vep_cli::regulatory::count_shard_features`); a shard that does not is an error, and a missing `regulatory/` directory writes no marker. `CacheManifest::write` then lists every file of the output directory, the manifest itself excepted, in `CHECKSUMS.sha256`, the file `vep-cli` checks each cache file it reads against. `--region-size` is parsed and ignored: the bounds come from the file names.

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
| `convert_cache`         | Creates the output directory, removes an earlier run's `info.json`, runs `convert_data_files` over `transcripts/`, `copy_variation_files` over `variations/` and `copy_regulatory_files` over `regulatory/`, then `copy_info_json`, and last writes `CHECKSUMS.sha256` over every file of the output directory (`vep_core::cache_manifest::CacheManifest::write`). |
| `copy_info_json`        | Parses and checks `info.json`, logs `species` and `assembly`, writes it with its `regulatory` marker set when regulatory shards were written and removed otherwise, and copies `chr_synonyms.txt`. |
| `copy_regulatory_files` | Copies each `regulatory/<chr>/<start>-<end>.json` or `.json.gz` shard byte for byte into `regulatory/<chr>/` once it parses as the annotator parses a shard; returns the count written. |
| `validate_info`         | Checks that `info.json` is an object with a non-empty `species` and `assembly`, and that `contig_order` and `variation_cols`, when present, are lists of strings. |
| `convert_data_files`    | Re-shelves each `{chr}_{start}-{end}.json` region file as `<chr>/<start>-<end>.json`, counting entries.     |
| `parse_region_filename` | Splits a region file stem at its last underscore into `(chromosome, start, end)`.                           |
| `copy_variation_files`  | Copies `all_vars.gz` with its index, or the `_var.gz` region files, into `variations/<chr>/`; refuses a `.json` and an unindexed `all_vars.gz`. |
| `variation_file_kind`   | The `VariationFileKind` a name denotes: `all_vars.gz`, its `.csi` or `.tbi` index, or a `<start>-<end>_var.gz` region file; none for any other name. |
| `json_shard_error`      | The error refusing a `.json` file under `variations/`.                                                      |
| `split_flat_variation_name` | Splits a flat `{chr}_<variation file>` name into the chromosome and the file's name under `variations/<chr>/`. |

## Command line

`vep-cache-converter --help` and `vep-cache-converter convert-plugin-data --help` print the two flag sets. Diagnostics go to stderr through `tracing` at INFO level; the mismatch listing of a failed verification is the one direct write.

Cache mode (no subcommand):

| Flag                   | Role                                                                                                        |
| ---------------------- | ----------------------------------------------------------------------------------------------------------- |
| `--json-dir <JSON_DIR>` | The flat cache: `{chr}_{start}-{end}.json` region files under `transcripts/`, the dumper's variation files under `variations/`. |
| `--output-dir <OUTPUT_DIR>` | Where the per-chromosome layout, the copied `info.json` and `chr_synonyms.txt`, then `CHECKSUMS.sha256` are written. Required. |
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

The unit tests cover `parse_region_filename`, `variation_file_kind`, `split_flat_variation_name`, `copy_variation_files` and `copy_info_json`. The workspace's standard build, test, lint and format commands are in [the contributing guide](../../CONTRIBUTING.md).

## Related

- [`vep-builtins`](../vep-builtins/README.md): `tabix::TabixRecord`, `binary_store::BinaryStoreWriter` and `binary_store::BinaryAnnotator`, the types the plugin-data path writes and reads back through; its plugins read the `.vpd` store.
- [`vep-cli`](../vep-cli/README.md): enumerates the `transcripts/<chr>/<start>-<end>.json` shards the cache mode writes. Declared in `Cargo.toml`, not imported by `main.rs`.
- [`vep-core`](../vep-core/README.md): `cache_manifest::CacheManifest`, whose `write` lists every file of the output directory in `CHECKSUMS.sha256`, and `coordinate::is_cache_path_component`, the check that a chromosome taken from a region file's stem or a `variations/` entry can name one directory of the output.
- [`vep-cache-builder`](../vep-cache-builder/README.md): writes the per-chromosome layout directly, so a natively built cache needs no conversion.
- [docs/cache-setup.md](../../docs/cache-setup.md): the routes to the per-chromosome JSON transcript cache the flat copy targets.
- [docs/plugins/data-files.md § Tabix indexes and `.vpd` stores](../../docs/plugins/data-files.md#tabix-indexes-and-vpd-stores): how a plugin finds a `.vpd` store beside its data file.
- [Root README § Installation](../../README.md#installation): the release archives and the container image that ship this binary.
