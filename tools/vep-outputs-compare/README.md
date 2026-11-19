# vep-outputs-compare

A Rust build of `scripts/concordance/compare_vep_outputs.py`, the tuple comparator of the
concordance harness: the same command line in its three modes, the same report files byte for
byte (`summary.json`, `summary.md`, `discordant.tsv`, `discordant_open.tsv`, the per-class and
pooled CSVs, the per-term totals cache) and the same stdout lines. The Python stays the reference
implementation; this binary exists because the comparator is the longest step of a whole-genome
concordance run after the annotation itself.

```text
vep-outputs-compare --perl-dir DIR --rust-dir DIR --report-dir DIR [options]
vep-outputs-compare --by-consequence-class --perl-dir DIR --discordant-in TSV \
    --dataset NAME --class-csv CSV --sweep-id ID --binary-md5 DIGEST [options]
vep-outputs-compare --pool-consequence-class --class-csv CSV --pooled-csv CSV \
    [--per-term TSV] [options]
```

Every option of the Python is accepted with its default, `--flag=value` and an unambiguous prefix
included: `--glob`, `--assembly`, `--canonical-contigs` and `--no-canonical-contigs`,
`--reference-release {115.2,116.2}`, `--input-dir` and `--input-vcf`, `--summary-json`,
`--summary-md`, `--discordant-tsv`, `--fail-below-f1`, `--engine {vep-rs,fastvep}`,
`--assembly-label`, `--append`, `--min-tuples`, `--vep-totals-cache`, `--per-term` and
`--per-term-engine`.

## Contract

Every behaviour of the Python is reproduced on the bytes it produces, the Python-specific parts
included. Records end at `\n`, `\r\n` or a lone `\r`, as universal newlines read them, and are
decoded as UTF-8 with `errors="replace"`. Fields are stripped with `str.strip()`'s whitespace set.
Every path option goes through `str(pathlib.Path(arg))`, and the directories the summary names
through `Path.resolve()`. The JSON is `json.dumps(indent=2, ensure_ascii=True)` with the insertion
order of every object, and its precision, recall and F1 values are `float.__repr__`.

The markdown and stdout lines format with `f"{x:.6f}"`; the per-class CSVs round with
`Decimal.quantize` under `ROUND_HALF_UP` and are written as `csv.DictWriter` writes them, `\r\n`
rows and minimal quoting. The top consequence sets keep `Counter.most_common`'s tie order. The
intended-divergence rules are the Python's: the six pair rules in its classification order, the
row-aware pair rules of the 116.2 registry after them, the 115.2 transcript-selection class over
records Ensembl VEP marks `vep_skip`, and the fourteen 116.2 one-sided classes, with their
definitions and citations carried verbatim into `divergence_classes`. A row-aware pair rule
reads, beside the pair's two consequence sets and its allele, the first column, the Location,
the cDNA_position and the CDS_position of the rows that carry the pair's two tuples, fetched from
the two outputs for the pairs no other rule takes. A row among them whose first column spells no
alleles takes the reference bases its allele stands in for from the record of the input VCF behind
its output file (`--input-vcf`, or the file of `--input-dir` with the output's stem), which is
read once more where such a row asks; no cache or other program is read.

A divergence class that reaches the whole-genome tuples is added to the Python registry and to
the Rust class table in the same change. `cargo test` reads the Python's source and fails while
the two name different classes.

The keys of each side are compared in the bytewise order the Python gets from `sort -u` under
`LC_ALL=C`. On one thread a side's keys are sorted in memory, spilled as sorted runs under the
report directory past one GiB of keys and merged on finish. On more threads the keys of both
sides are routed by Location into partitions that cover that order end to end, written under the
report directory, and each partition is sorted and merged on its own; read in partition order the
keys, the counts and the discordant rows are those of the one sort. The merge holds every key
against the one before it on its side, and the first key of each partition against the last key
before it, and at the first that does not sort after it the run ends with exit code 1, writes no
summary, keeps the key files it stopped on in their directory under the report directory and
leaves the discordant rows written so far under the `--discordant-tsv` name with a `.partial`
suffix, naming both in the error: past such a key the merge would count a key both sides hold as
missing and as extra at once. The `--glob` pattern names
entries of the output directory itself (`*`, `?` and `[...]`); a pattern with a `/` is refused.
The one-sided classes of the 116.2 reference read a Location's coordinates and a breakend mate's
position, written in the digits `0` to `9`, whole at any length, as Python does; every other
coordinate and position is read as an integer up to 128 bits, where Python's are unbounded. Exit
codes: 2 on a usage error or an aggregate F1 below `--fail-below-f1`, 1 on a missing mode argument
(`error: mode <mode> requires ...` on stderr) or any other failure
(`ERROR: [vep-outputs-compare] <detail>` on stderr), 0 otherwise.

## Threads

The comparison runs across the host's available cores. `VEP_OUTPUTS_COMPARE_THREADS=<n>` bounds
the threads at work, `1` keeping the whole comparison on the calling thread; it is an environment
variable and not an option, so that the command line stays the Python's, and a value that is not a
whole number of 1 or more is refused. Every thread count writes the same bytes to every report file
and to stdout.

One thread reads an output file in runs of whole records while the others key them, and every
thread then sorts and merges partitions. The pass over `discordant.tsv` that applies the
intended-divergence rules, the read of an input VCF, the read of the record names behind the
skipped-record class, the read of the rows behind the row-aware pair rules and the
`--pool-consequence-class` mode stay on one thread.

## Layout

- `src/input.rs`: records, lossy decoding, the five-column tuple key with its Location
  normalised and the canonical-contig filter, the record names of the rows it keys, and the
  reader of an input VCF, plain or gzipped.
- `src/sort.rs`: the in-memory sort-unique with its spilled runs, and the sorted keys of one
  partition.
- `src/partition.rs`: the bounds between partitions sampled from an output file, the scan that
  routes a side's keys into partition files across threads, and the files themselves.
- `src/threads.rs`: the thread bound and the indexed fan-out over partitions.
- `src/merge.rs`: the streaming merge of two sorted key streams, the discordant rows and the
  per-set counters, feeding the selected release's one-sided collector.
- `src/rules/`: the pair rules (`pair.rs`), what a row-aware pair rule reads of an output row
  and of its record's alleles (`row.rs`), the Location and breakend readers (`breakend.rs`),
  the two collectors (`scope.rs`), the gVCF variant-site recognisers (`gvcf.rs`), the 116.2
  one-sided pass (`one_sided.rs`), the class table (`table.rs`) and the pass over
  `discordant.tsv` (`filter.rs`).
- `src/classes.rs`: the per-consequence-class counts, their pooling, the floor, the totals cache
  and the two CSV schemas.
- `src/report.rs`: the JSON emitter, the summary markdown and the per-file metrics.
- `src/compare.rs`, `src/per_class.rs`, `src/main.rs`: the three modes and the command line.
- `src/pyfmt.rs` and `src/csv.rs`: the CPython formatting and `csv` module rules listed above.
- `tests/golden.rs` and `tests/fixtures/`: synthetic output directories (`perl/`, `rust/`), an
  input VCF and the goldens the Python wrote over them in the three modes. The goldens carry
  `{{PERL_DIR}}`, `{{RUST_DIR}}`, `{{INPUT_DIR}}`, `{{REPORT_DIR}}`, `{{SCRATCH}}` and `{{GENERATED_AT}}` for
  the values that depend on the checkout and the run; `cargo test` runs the binary over the
  fixtures on one thread, on two and on the host's default, substitutes them and compares bytes,
  and compares the reports of a generated pair of directories across thread counts.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-outputs-compare/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

CI runs the test, clippy and fmt commands above in its Test, Clippy and Format jobs, builds the
package on the declared `rust-version` in MSRV and audits its lockfile. Its Comparator Tests job
regenerates the goldens with the Python and compares them with the committed files, so a change to
either side that moves the contract fails there first. That step of `.github/workflows/ci.yml` is
also the regeneration recipe: the same commands, their outputs written over
`tests/fixtures/expected_*` after the mask it applies.

## Harness

`scripts/concordance/run_clone_measurement.sh`, `scripts/concordance/run_concordance.sh` and
`scripts/validation/compare_cache_outputs.sh` score tuples with this binary by default. One
environment variable, `VEP_OUTPUTS_COMPARATOR`, selects the comparator in all three:

- Unset or empty: this tool. `run_clone_measurement.sh` and `compare_cache_outputs.sh` run this
  tree's release build, `tools/vep-outputs-compare/target/release/vep-outputs-compare`, and build
  nothing; `run_concordance.sh` builds it through cargo, as it builds the engine it runs.
- `python`: `scripts/concordance/compare_vep_outputs.py`, the reference implementation, which
  writes the same reports. No other value runs the Python.
- Any other value: a comparator binary, as a path or a command on `PATH`.

A comparator that is not an executable, or a build that fails, stops the script before its first
engine run with one line that names both remedies: build or name the binary, or select the Python.
`run_clone_measurement.sh` asks for a comparator only where it scores F1: not under
`--walltime-only`, and not for the Perl engine.
