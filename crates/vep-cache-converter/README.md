# vep-cache-converter

`vep-cache-converter` is a command-line binary, shipped in every release beside `vep` and `vep-cache-builder`, that prepares data files for the annotator. Without a subcommand it copies a flat JSON (JavaScript Object Notation) transcript cache into the per-chromosome layout `vep-cli` loads, validating each region file; `convert-plugin-data` converts a tabix-indexed plugin data file into the binary annotation store (`.vpd` data plus `.vpdi` index) that the `vep-builtins` plugins memory-map at runtime. No crate depends on it; it depends on `vep-builtins` for its record and store types and on `vep-cli` for the regulatory shard parser.

## How it works

`convert-plugin-data` follows one data line from the input file to a verified store. An `--output` that is the input table or its `.tbi`, or whose `.vpdi` path is, is refused before anything is read or written, since the store's files replace whatever lies at their paths (`is_same_file`). `read_tabix_records` opens the input through a BGZF (blocked gzip) reader over a `HashingReader`, so the SHA256 and length it returns are of the bytes the conversion read, keeps the last `#` line before the first row as the column names, or, under no `#` line, the first row itself when its start cell leads with no integer, the bare header the plugins' tabix reader takes (`parse_header`; a `#` line after the first row is a comment), and parses each data line into a `TabixRecord`: the chromosome from `--chr-col`, the start from `--start-col` read to its leading integer as tabix reads it (plus one when `--zero-based` is set), the end by the tabix reader's span rule (`TabixConfig::record_end`: with `--ref-col`, the `--end-col` column's `END` key at or past the start, else the allele's span; without `--ref-col`, the `--end-col` cell's leading integer; the start alone without `--end-col`), the raw columns, and a shared pointer to the header. A line with fewer columns than `--chr-col` or `--start-col` needs is skipped; a start that no integer leads is an error. A line is split with its line end alone removed, as the plugins read it through tabix, so a last cell keeps its trailing whitespace and an empty last cell stays a cell ([Tabix indexes and `.vpd` stores](../../docs/plugins/data-files.md#tabix-indexes-and-vpd-stores)). Records are grouped by chromosome in a `BTreeMap`, so chromosomes come out in string order (`1`, `10`, `11`, `2`), and each group is sorted by start alone, a stable sort, so the rows at one start keep the file's order, the order `tabix` prints them in. `BinaryStoreWriter::write` from `vep-builtins` then writes the column table and the records to the `--output` path and, to the same path with the extension `.vpdi`, the layout flags the records were read by, the SHA256 and length of the input and of the data file, and the chromosome and position index; each file is written whole under another name beside its own and renamed into place, the data first and the index last, a path that is a symbolic link at the file it names, so a `vep` run that starts during a conversion finds no part-written store and one that holds the previous store open keeps it; a cell, column name or chromosome name over 65,535 bytes, or a header or row of more than 65,535 cells, is refused naming the row and the length, since the store's length fields are 16 bits wide, a position of more than 4,294,967,295 rows is refused naming the position, since an index entry counts rows in 32 bits, and a refused or failed write leaves neither file behind; a plugin whose table reads by other columns refuses the store, naming the flags that rebuild it, and reads the table through tabix, as does a plugin whose table, or whose store's data file, differs in length or SHA256 from what the index records, naming which. Unless `--no-verify` is set, `open_table` opens the input through `TabixAnnotator`, the reader the plugins use without a store, before anything is read or written, so an input without its index fails with no store left behind; after the write, `verify_concordance` reopens the store through `BinaryAnnotator`, the reader the plugins use at runtime, and checks every distinct record start when the file holds at most 1,000 of them or `--verify-all` is set, otherwise every n-th start from the first, n the smallest stride that keeps the sample at 1,000 starts or fewer, so a conversion checks the same starts on every run. `probe_positions` lists, for each checked start, the start and each of its records' end, midpoint, the position past its end and the one before its start; `verify_probes` queries each position through both readers and compares their records in the order each returns them, span and every cell (`record_line`), then asks the store for the position again through its batch query, alone and then with the positions beside it in batches of 64, in order and reversed, against the same answer of the table; the first ten mismatches are printed, the store and its index are removed (`remove_store`), and the exit is non-zero; a verification that does not run to its end (a query the table's index cannot answer) removes them the same way.

The cache mode follows one region file. `convert_cache` creates the output directory and, after the transcript and variation files, `read_info_json` reads `info.json`, parses it, checks what the runtime's cache resolution reads (an object with a non-empty `species` and `assembly`, a `contig_order`, `variation_cols` or `cell_types` that is a list of strings when present, a `regulatory_schema` that is a whole number when present, and a `bam`, the alignment file the cache's transcript models were edited against, or a `source_type`, the transcript set the cache holds, that is a string when present) and logs its `species` and `assembly`; after the regulatory files, `copy_info_json` writes it with the `regulatory` marker set when `copy_regulatory_files` wrote a shard, `cell_types` and `regulatory_schema` passing through, all three removed otherwise, and copies `chr_synonyms.txt` beside it when the source has one; a missing, unparseable or unusable `info.json` is an error. `convert_data_files` runs over `transcripts/`: it lists the directory in file-name order, takes each `.json` file, and `parse_region_filename` splits the stem at its last underscore into chromosome and range, so `LRG_485_1-1000000` yields the LRG (Locus Reference Genomic) sequence name `LRG_485` as its chromosome; a stem that does not fit is skipped with a warning. The content must parse as JSON (a parse failure is an error) and be an array (anything else is skipped with a warning); its length is counted and the text is written unchanged to `<output>/transcripts/<chr>/<start>-<end>.json`, the shard layout `vep-cli` enumerates. `copy_variation_files` runs over `variations/`, which holds the files the Ensembl VEP dumper writes: `all_vars.gz` with its `.csi` or `.tbi` index, or `{start}-{end}_var.gz` region files, under `variations/<chr>/`, or carrying a `{chr}_` prefix directly under `variations/`. `split_flat_variation_name` finds the chromosome in a flat name by trying each underscore until what follows is a variation file name, and each file is copied byte for byte to `<output>/variations/<chr>/`; any other name, flat or under `<chr>/`, is skipped with a warning. An `all_vars.gz` without its index is an error, as is a `.json` file under `variations/`: the annotator reads the dumper's files only. A missing `variations/` directory is a warning, not an error. `copy_regulatory_files` runs over `regulatory/`: it refuses a tree whose metadata lacks the `regulatory_schema` this build reads, loads and copies the binding matrices table `regulatory/binding_matrices.json`, then copies each `regulatory/<chr>/<start>-<end>.json` or `.json.gz` shard byte for byte once it parses as the annotator parses a shard (`vep_cli::regulatory::count_shard_features`, against the metadata's `cell_types` list and the table); a shard that does not is an error, any other entry is skipped with a warning, and a missing `regulatory/` directory writes no marker. `CacheManifest::write` then lists every file of the output directory, the manifest itself excepted, in `CHECKSUMS.sha256`, the file `vep-cli` checks each cache file it reads against. `--region-size` is parsed and ignored: the bounds come from the file names.

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
| `is_same_file`          | Whether two paths name one file, for the refusal of an `--output` that is the input table.                  |
| `read_tabix_records`    | Streams the BGZF input into `TabixRecord`s grouped by chromosome and sorted by position, with the SHA256 and length of the bytes read. |
| `open_table`            | Opens the input with `TabixAnnotator`, lines handed over as to a plugin, before the conversion reads or writes anything. |
| `verify_concordance`    | Reopens the store with `BinaryAnnotator` and compares its records with the table's at the probed positions.   |
| `probe_positions`       | The positions probed: every record start, or every n-th past 1,000 without `--verify-all`, and each of its records' end, midpoint and the positions on either side of its span. |
| `verify_probes`         | Compares the store's answers with the table's at each position: a single query, a batch of the position alone, and batches of 64 positions in order and reversed. |
| `record_line`           | A record as one comparable line: its span and its cells, the chromosome cell restored from the record's chromosome. |
| `remove_store`          | Removes the store and its index after a verification that mismatches or does not run to its end.            |
| `convert_cache`         | Creates the output directory, removes an earlier run's `info.json`, runs `convert_data_files` over `transcripts/`, `copy_variation_files` over `variations/` and `copy_regulatory_files` over `regulatory/`, then `copy_info_json`, and last writes `CHECKSUMS.sha256` over every file of the output directory (`vep_core::cache_manifest::CacheManifest::write`). |
| `read_info_json`        | Reads and parses `info.json`, checks it with `validate_info` and logs `species` and `assembly`; its `SourceInfo` carries the text and the parsed value. |
| `copy_info_json`        | Writes `info.json` with its `regulatory` marker set when regulatory shards were written (`cell_types` and `regulatory_schema` pass through) and the three keys removed otherwise, and copies `chr_synonyms.txt`. |
| `copy_regulatory_files` | Refuses a regulatory tree without the schema this build reads, copies `regulatory/binding_matrices.json`, then each `regulatory/<chr>/<start>-<end>.json` or `.json.gz` shard byte for byte once it parses as the annotator parses a shard; returns the count written. |
| `regulatory_shard_range` | The `(start, end)` range a regulatory shard's name carries, `<start>-<end>.json` or `.json.gz`; none for any other name. |
| `warn_of_stray_regulatory_file` | Warns of a file directly under `regulatory/` other than the binding matrices table. |
| `validate_info`         | Checks that `info.json` is an object with a non-empty `species` and `assembly`, that `contig_order`, `variation_cols` and `cell_types`, when present, are lists of strings, that `regulatory_schema`, when present, is a whole number, and that `bam` and `source_type`, when present, are strings. |
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
| `--output <OUTPUT>`   | The `.vpd` data file; the `.vpdi` index is written beside it. Neither path may be the input table or its `.tbi`. |
| `--chr-col <CHR_COL>` | 0-based column index of the chromosome. Default `0`.                                                         |
| `--start-col <START_COL>` | 0-based column index of the start position. Default `1`.                                                     |
| `--end-col <END_COL>` | 0-based column index of the end position; omitted, the end equals the start. With `--ref-col` the column is a VCF INFO column, read for its `END` key, the reference allele's span otherwise; without it the cell is read to its leading integer. |
| `--ref-col <REF_COL>` | 0-based column index of the reference allele (a VCF REF column), read beside `--end-col`: its length spans a record whose end column holds no `END` at or past the start. |
| `--zero-based`        | Positions are 0-based, as in BED (Browser Extensible Data); one is added. Default: 1-based.                  |
| `--verify-all`        | Verify every record after writing, not a sample of at most 1,000 record starts.                              |
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

The unit tests cover `convert_plugin_data`, `read_tabix_records`, `open_table`, `verify_concordance`, `probe_positions`, `verify_probes`, `remove_store`, `convert_cache`, `read_info_json`, `copy_info_json`, `copy_regulatory_files`, `regulatory_shard_range`, `validate_info`, `convert_data_files`, `parse_region_filename`, `copy_variation_files`, `variation_file_kind` and `split_flat_variation_name`. The workspace's standard build, test, lint and format commands are in [the contributing guide](../../CONTRIBUTING.md).

## Related

- [`vep-builtins`](../vep-builtins/README.md): `tabix::TabixRecord`, `binary_store::BinaryStoreWriter` and `binary_store::BinaryAnnotator`, the types the plugin-data path writes and reads back through; its plugins read the `.vpd` store.
- [`vep-cli`](../vep-cli/README.md): enumerates the `transcripts/<chr>/<start>-<end>.json` shards the cache mode writes; its `regulatory::count_shard_features` and `regulatory::BindingMatrices` are the shard parser and binding matrices table `copy_regulatory_files` checks each regulatory shard with.
- [`vep-core`](../vep-core/README.md): `cache_manifest::CacheManifest`, whose `write` lists every file of the output directory in `CHECKSUMS.sha256`, and `coordinate::is_cache_path_component`, the check that a chromosome taken from a region file's stem or a `variations/` entry can name one directory of the output.
- [`vep-cache-builder`](../vep-cache-builder/README.md): writes the per-chromosome layout directly, so a natively built cache needs no conversion.
- [docs/cache-setup.md](../../docs/cache-setup.md): the routes to the per-chromosome JSON transcript cache the flat copy targets.
- [docs/plugins/data-files.md § Tabix indexes and `.vpd` stores](../../docs/plugins/data-files.md#tabix-indexes-and-vpd-stores): how a plugin finds a `.vpd` store beside its data file.
- [Root README § Installation](../../README.md#installation): the release archives and the container image that ship this binary.
