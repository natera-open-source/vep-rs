# Contributing to vep-rs

## Building and Testing

```bash
# Run the full test suite (1,179 tests: 1,178 run, 1 #[ignore]d; the per-crate
# counts are in manuscript/data/code_metrics.csv)
cargo test --workspace

# Build a release binary (target/release/vep)
cargo build --release

# Lint (CI enforces -D warnings)
cargo clippy --workspace -- -D warnings

# Check formatting
cargo fmt --all -- --check
```

For faster iteration during development:

```bash
# Consequence changes
cargo test -p vep-effects consequences && cargo test -p vep-cli

# Plugin changes
cargo test -p vep-builtins
```

## Pre-Commit Hooks

The repo includes a `.pre-commit-config.yaml` that runs `cargo fmt` and `cargo clippy` before each commit. Install with:

```bash
pre-commit install
```

## Concordance Testing

Concordance testing compares Rust VEP output against Perl VEP (release 115) to verify semantic parity. Requires Docker and a Perl VEP cache.

```bash
scripts/concordance/run_concordance.sh --mode smoke --smoke-variants 5000 \
    --perl-cache-dir /path/to/vep_cache --json-cache-dir /path/to/json_cache \
    --fasta /path/to/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
```

`--perl-cache-dir` (the directory above `homo_sapiens/`), `--json-cache-dir` and `--fasta` (or an explicit `--no-fasta`) are required; the harness exits before annotating without them.

The data tree and the measurement entrypoint behind the published concordance figures are in [scripts/README.md](scripts/README.md#reproducing-the-published-concordance).

## Code of Conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md). By
participating, you are expected to uphold it.

## Security

Please report security vulnerabilities privately per [SECURITY.md](SECURITY.md),
not via a public issue.

## Golden corpora and output parity

`tests/golden/<release>/<assembly>/` holds, per Ensembl release and assembly, a
corpus of input records covering every consequence-set combination observed on
the benchmark datasets, Ensembl VEP's output for them in all four formats
(`expected/`, gzipped), a pruned JSON cache that reproduces vep-rs's full-cache
output on those records, and a `manifest.json` listing the exemplars per
combination and the consequence keys whose terms are documented to differ.
`cargo test -p vep-cli --test golden --test format_parity` annotates each corpus
and compares every column of every format (the Parquet round trip needs the
`duckdb` CLI; without it that one test skips). A failure lists every mismatch
grouped by column and consequence set.

To add a release or assembly, run the generator where the VEP reference outputs,
their inputs and the JSON cache are (`scripts/golden/build_golden_corpus.py
select`, `scripts/golden/prune_json_cache.py --gzip --drop sorted_exons --drop
protein_features --drop protein_function_predictions --perl-info <cache>/info.txt`,
VEP in the four formats, then `build_golden_corpus.py classify` against the
pinned vep-rs binary), record the VEP image digest and the expected files'
sha256 in `provenance.json`, and commit the new directory; older directories
stay until their release is dropped. Each corpus must stay under 15 MB on
disk (`each_corpus_fits_its_size_budget`); prefer narrower exemplars
(`--max-span`) over a larger cache.

## Branch and Commit Conventions

- Create feature branches off `main`.
- Write commit messages in imperative mood, focusing on what and why ("Add CADD plugin support", not "Added some plugin stuff").
- Keep the summary line under 72 characters.

## Sign your commits

Every commit must carry a `Signed-off-by` trailer. Pass `-s` to `git commit`, which
appends `Signed-off-by: Your Name <your.email@example.com>` from your git identity:

```bash
git commit -s -m "Add CADD plugin support"
```

The sign-off asserts, under the [Developer Certificate of Origin](https://developercertificate.org),
that you have the right to submit the code under this project's Apache-2.0 license.
A `DCO` workflow (`.github/workflows/dco.yml`) checks every commit of a pull request
for the trailer and fails the check if one is missing. To add a missing sign-off to
the commits on your branch, use `git rebase --signoff main` and force-push the branch.

## Pull Requests

Before submitting a PR:

1. All CI checks pass: `cargo test --workspace`, `cargo clippy --workspace -- -D warnings`, `cargo fmt --all -- --check`.
2. Add tests for new features or bug fixes.
3. Update the affected page under `docs/`: `cli-reference.md` for flags, `output-formats.md` for output layouts, `plugins.md` for plugins, `cache-setup.md` for the cache builder and converter.
4. For consequence or plugin changes, run concordance smoke to confirm no regressions.

## Comment and write-up discipline

A comment is one or two present-tense sentences stating a fact the code cannot: a contradicting engine or
domain behaviour, a constant's derivation, a rejected alternative, or the reference behaviour being
matched, cited (`Pick order follows Bio/EnsEMBL/VEP/Config.pm:306`). Plain `//` comments default to none.
Doc comments (`///`, `//!`) are the callable's contract as it renders on docs.rs and carry no dates or
history; a measured figure stays only where it explains a constant. Comments name nothing outside the
code: no organisation, infrastructure, account, issue key, person, date, or work-in-progress state; a
comment that cannot be reworded into a fact about the code is removed. Commented-out code is deleted.
Keep the two-line license header and the blank line after it on every source file.

## Code Conventions

- **Error handling**: `anyhow::Result` with `.context()` naming the file path in CLI and runner code; `thiserror` derives for typed errors in library crates.
- **Logging**: `tracing` for every diagnostic (`info!` for user-facing status, `debug!` behind `--verbose`, `warn!` for a recoverable condition such as a `--plugin` name that resolves to neither a built-in nor a dylib and is skipped). The only direct writes are the `vep` banner and end-of-run summary on stderr (`eprintln!` in `main.rs` and `runner.rs`, suppressed by `--quiet`), `vep-cache-converter`'s mismatch listing, and the `bench_plugins` tool's progress and report output.
- **Output ordering**: consequence terms within a `TranscriptConsequence` are ordered by SO rank; Extra field keys follow VEP's flag-field order (`@FLAG_FIELDS`, the fields each active flag adds, in VEP's group order; IMPACT, DISTANCE, STRAND and FLAGS are always active), then every other key, plugin keys included, alphabetically; within a record, rows are ordered by transcript id, then by ALT allele index, with one intergenic row per allele that overlapped nothing last.
- **Naming**: chromosome names are normalized for cache lookup when the `InputVariant` is built (`chr21` -> `21`, `chrM` -> `MT`), and `original_chr` keeps the input spelling for the Location column, as VEP does; struct fields follow Perl VEP's names where possible (`transcript_id`, `gene_symbol`); SO terms are Perl VEP's exact strings (`5_prime_UTR_variant`, `NMD_transcript_variant`).
- **Plugin data**: `IndexMap<String, String>` keyed by field name, variant-level in `InputVariant.plugin_data` and transcript-level in `TranscriptConsequence.plugin_data` (dylib plugin results merge there, one Extra field per key of the returned JSON object); output order does not depend on the map: the default-format Extra column and the JSON writer list plugin keys alphabetically after the flag fields, and the tab, VCF and Parquet writers take their column order from the plugins' declared field lists.
- **Shared strings**: `Arc<str>` for frequently cloned string fields (stable ids, gene symbols); a clone is a refcount increment.
