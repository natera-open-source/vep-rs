# Contributing to vep-rs

## Building and Testing

```bash
# Run the full test suite (the per-crate test counts are in
# manuscript/data/code_metrics.csv)
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

Concordance testing compares vep-rs output against Ensembl VEP release 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the image `run_concordance.sh` pins; `--perl-vep-docker-image` selects another) to verify semantic parity. Requires Docker and an Ensembl VEP cache of that release.

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
combination, the consequence keys whose terms are documented to differ, the
parity-matrix rows the corpus exercises (`covers`, written by `select --covers`
or `covers --add`; `covers --add`, `flag-set --covers` and `select --matrix` refuse
a name that is not a row of `tests/parity/matrix.tsv`) and, under
`flag_sets`, the further runs of the same records under other flags, each with
its own `flags`, input, formats, expected directory and `covers`, `unordered`
when Ensembl VEP emits the record's rows in hash order, and `perl_mode`
(`database`) when its reference came from Ensembl VEP connected to a database
rather than reading the cache, so its default and tab headers compare by their
column descriptions alone (`flag-set`, `flag-set --perl-mode`). `masks` lists the
mask classes applied to every run of the corpus before the comparison:
`regulatory_feature_rows` sets aside the reference's `RegulatoryFeature` and
`MotifFeature` rows for a corpus run under `--regulatory` or `--everything` against
a cache converted without its regulatory tree, where vep-rs writes none, and
`order_only` compares the `DOMAINS`,
`FREQS` and `ZYG` items of both sides as multisets, for a native cache that lists a
transcript's protein features in another order than the Ensembl cache and for the
sample order of `--individual_zyg`, which Ensembl VEP prints in hash order
(`scripts/concordance/compare_vep_formats.py --masks` applies the same two classes, its
`order_only` over the same three columns). A
flag set may carry its own `masks` for its run alone, kept when `flag-set` rewrites
the set's entry. The list is written into the manifest by hand (no `select` option
writes it), `select` keeps it when it regenerates the corpus, every subcommand
refuses a class other than those two, and `classify` applies the regulatory class
to the reference and compares the three columns as multisets under `order_only`.
`cargo test -p vep-cli --test golden --test format_parity` annotates each corpus
once per run and compares every column of every format (the Parquet round trip
needs the `duckdb` CLI; without it that one test skips). A failure lists every mismatch
grouped by column and consequence set. A record Ensembl VEP dropped before
annotation carries `reference_warning` (VEP's reason, from its warnings file,
recorded when the corpus is built); `classify` writes `reference_rows` 0 on it
and classes the keys only vep-rs emits for it `reference_skipped_record`, the
one class under which the harness accepts a record-level VCF line or JSON
object with no reference counterpart. A manifest may also carry
`reference_warnings`, the texts Ensembl VEP wrote to its warnings file about
the run rather than about one record (the notice a `<CNV:TR>` record draws
without a FASTA, for one); the harness expects them in the run's warnings file
beside the numbered reasons. A flag set may carry its own list for its run
alone; one without inherits the manifest's. Like `masks`, the list is written
by hand, kept by `select` and by a `flag-set` rewrite of the set, and refused
by `build_golden_corpus.py` and the golden tests alike when a text is empty.

To add a release or assembly, run the generator where the VEP reference outputs,
their inputs and the JSON cache are: select the exemplars, prune the cache, write
the four expected outputs with `generate_reference_output.sh --corpus`, then
classify against the pinned vep-rs binary.

```sh
scripts/golden/build_golden_corpus.py select
scripts/golden/prune_json_cache.py --gzip --drop sorted_exons \
    --drop protein_features --drop protein_function_predictions \
    --perl-info <cache>/info.txt
scripts/concordance/generate_reference_output.sh --corpus <dir> \
    --perl-cache-dir <perl cache>
scripts/golden/build_golden_corpus.py classify
```

`generate_reference_output.sh --corpus` runs one Ensembl VEP container per format
from the manifest's `assembly`, `flags` and `fasta`, records the image digest, the
vep command and the expected files' sha256 in `provenance.json` for you, and
`--force` regenerates an existing corpus. A manifest listing
`regulatory_feature_rows` under `masks` has its reference's `RegulatoryFeature` and
`MotifFeature` rows left out of the classification, as the golden test leaves them
out of the comparison. Commit the new directory; older directories stay until
their release is dropped. A corpus for Ensembl VEP's
per-sample flags (`--individual`, `--individual_zyg`, `--process_ref_homs`) is
selected with `--keep-samples N`, which keeps `FORMAT` and the first N sample
columns and adds records until every genotype class those samples show in the
source has exemplars (`check` refuses a corpus missing one); a class the corpus
must carry whether or not the source shows it is named with
`--require-genotype-class`, which `select` and `check` both enforce.
`--include sid=name` carries the record of suite `sid` that `name` names, by its
ID or by Ensembl VEP's auto-name, whatever the sampling chooses, and
`--keep-records` keeps every record of the manifest being regenerated. Each corpus
must stay under 15 MB on disk (`each_corpus_fits_its_size_budget`), every flag
set's expected directory included; prefer narrower exemplars (`--max-span`) and,
with `select --cache <json cache>`, the exemplars whose transcripts add the fewest
bytes to the pruned cache, over a larger cache.

## Parity matrix

`tests/parity/matrix.tsv` lists every Ensembl VEP parameter, feature, plugin and cache field with the
behaviour vep-rs targets for it and, once landed, the tests that prove it;
[tests/parity/README.md](tests/parity/README.md) explains the columns and how a row closes. A pull
request that adds a flag to `crates/vep-cli/src/args.rs` adds its row, and one that lands a row names
its tests, because CI runs `python3 scripts/parity/check_parity_matrix.py` (a flag without a row, a
landed row without its tests, or a named test or corpus that does not resolve fails) and
`crates/vep-cli/tests/parity_matrix.rs` (every long flag clap reports has a row).

## Branch and Commit Conventions

- Create feature branches off `main`.
- Write commit messages in imperative mood, focusing on what and why ("Add CADD plugin support", not "Added some plugin stuff").
- Keep the summary line under 72 characters.

## Changelog entries

`CHANGELOG.md` follows Keep a Changelog, and the release page either copies the version's section word
for word or restates its entries grouped by area, so each entry is written for the page: one bullet per
user-visible change, on one physical line; the subject first (the flag, format, input or behaviour,
identifiers in backticks), then what it now does; one sentence, a second only for a consequence the user
must act on; the measured effect inline where the change claims speed or memory; no type prefix, bold
lead-in or author handle; a trailing `(#N)` naming the pull request; at most 30 words. The pull request
number is written only after the pull request exists, never predicted. Refactors, CI and dependency
bumps appear only when behaviour changes. Under a kind heading, entries may be grouped by area: a bullet
`- <Area name>` at the margin (the area alone, at most 60 characters) with its entries nested under it
as `-` bullets indented two spaces, each in the grammar above with its own `(#N)`; two levels at most,
and every area has at least one entry.

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
3. Update the affected page under `docs/`: the format's page under `output/` for output layouts, the mechanism's page under `plugins/` and the plugin table in `plugins.md` for plugins, `cache-setup.md` for the cache builder and converter. A flag is documented by its help text in `crates/vep-cli/src/args.rs` and by the option tables, from which the CLI reference (`cli-reference.md` and the group pages under `cli/`) is generated: regenerate it with `VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference` and commit the result.
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
- **Logging**: `tracing` for every diagnostic (`info!` for user-facing status, `debug!` behind `--verbose`, `warn!` for a recoverable condition). In the `vep` binary every `warn!` event reaches the warning sink (`crates/vep-cli/src/warnings.rs`), which writes it to the warnings file as `WARNING: <text>` and echoes it on stderr unless `--quiet`; a record the run skips goes through the sink's `skipped_variant`, never through `warn!`, so it is counted for `--skipped_variants_file`. The only direct writes are the sink's, the `vep` banner, the `INFO:` lines of the setup checks and the end-of-run summary on stderr (`eprintln!` in `main.rs` and `runner.rs`, suppressed by `--quiet`), the repeated `-o` warning the `vep` script prints before any option is read (`main.rs`, whatever `--quiet` says), the one-line `ERROR: [<component>] <message>` a refusal or a failed run ends with, `vep-cache-converter`'s mismatch listing, and the `bench_plugins` tool's progress and report output.
- **Output ordering**: consequence terms within a `TranscriptConsequence` are ordered by SO rank; Extra field keys follow VEP's flag-field order (`@FLAG_FIELDS`, the fields each active flag adds, in VEP's group order; IMPACT, DISTANCE, STRAND and FLAGS are always active), then every other key, plugin keys included, alphabetically; within a record, rows are ordered by feature kind (transcript, regulatory, motif), then by stable id, then by ALT allele index, with one intergenic row per allele that overlapped no transcript last.
- **Naming**: the `vep` binary resolves each input chromosome name to the cache's before any lookup, as Ensembl VEP's `get_source_chr_name` does (the name itself when the cache holds it, else a synonym from `--synonyms` or the cache's `chr_synonyms.txt`, else the `chr` prefix removed or added: `chr21` -> `21`; `chrM` reaches `MT` through a synonym), and `original_chr` keeps the input spelling for the Location column, as VEP does; struct fields follow Ensembl VEP's names where possible (`transcript_id`, `gene_symbol`); SO terms are Ensembl VEP's exact strings (`5_prime_UTR_variant`, `NMD_transcript_variant`).
- **Plugin data**: `IndexMap<String, String>` keyed by field name, variant-level in `InputVariant.plugin_data` and transcript-level in `TranscriptConsequence.plugin_data` (dylib plugin results merge there, one Extra field per key of the returned JSON object); output order does not depend on the map: the default-format Extra column and the JSON writer list plugin keys alphabetically after the flag fields, and the tab, VCF and Parquet writers take their column order from the plugins' declared field lists.
- **Shared strings**: `Arc<str>` for frequently cloned string fields (stable ids, gene symbols); a clone is a refcount increment.
