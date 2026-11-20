# Contributing to vep-rs

How to contribute to vep-rs, a Rust reimplementation of Ensembl VEP: building and testing, concordance checks, signed commits and pull requests.

## Building and Testing

```bash
# Run every test (the per-crate test counts are in
# manuscript/data/code_metrics.csv)
cargo test --workspace --locked

# Build a release binary (target/release/vep)
cargo build --release

# Lint (CI enforces -D warnings; the feature brings in the bench_plugins binary)
cargo clippy --workspace --features vep-cli/bench -- -D warnings

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

When a release moves the Ensembl VEP reference, every entry of `docs/intended-divergences.md` is run again on the new release with `scripts/divergences/reevaluate.py`, an entry the release closes is retired with that release's number, and the page's `Verified against` column names the release last run.

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
its own `flags` (a word may carry `{corpus}`, the corpus directory, so a flag
names a file beside the manifest, such as the `--custom` sources
`scripts/golden/make_custom_fixtures.py` writes under `custom/`, and both engines
read the same file: vep-rs the directory itself, Ensembl VEP a read-only mount of
it; a word may also name a corpus file by a relative path, `file=data/scores.tsv.gz`,
since both engines run in the corpus directory, and a plugin whose header text
repeats the path it was given then writes the same header line on both sides;
or `{generated}`, the directory the manifest's `generated` block writes before
the run through `python3`, named by its `command` and the `sha256` of every file
it writes, for a fixture a corpus reproduces from a seed instead of committing,
such as the MaxEntScan matrices `scripts/fixtures/maxentscan_fixture_matrices.py`
writes), input, formats, expected directory and `covers`, `unordered`
when Ensembl VEP emits the record's rows in hash order, and `perl_mode`
(`database`) when its reference came from Ensembl VEP connected to a database
rather than reading the cache, so its default and tab headers compare by their
column descriptions alone (`flag-set`, `flag-set --perl-mode`), or (`gtf`, `gff`)
when Ensembl VEP read its transcripts from the file given to `--gtf` or `--gff`
while vep-rs reads the cache `vep-cache-builder` built from that file: the headers
compare the same way, and the source's `SOURCE` column and header lines, which no
cache carries, are set aside. A manifest may name its own run's `perl_mode` the
same way, and a flag set naming none inherits it. `masks` lists the
mask classes applied to every run of the corpus before the comparison:
`regulatory_feature_rows` sets aside the reference's `RegulatoryFeature` and
`MotifFeature` rows for a corpus run under `--regulatory` or `--everything` against
a cache converted without its regulatory tree, where vep-rs writes none, and
`order_only` compares the `DOMAINS`,
`FREQS`, `ZYG` and `NearestExonJB` items of both sides as multisets, and the
`5UTR_annotation` groups by their pairs, for a native cache that lists a
transcript's protein features in another order than the Ensembl cache and for the
orders Ensembl VEP prints from a hash: the samples of `--individual_zyg`, the exon
items of NearestExonJB and the pairs of UTRAnnotator
(`scripts/concordance/compare_vep_formats.py --masks` applies the same two classes, its
`order_only` over the same columns). Three further classes give `classify`
nothing to apply and belong to the golden harness: `unloaded_module_plugin` sets aside
the reference's warnings about a plugin module vep-rs cannot load,
`conservation_short_name_spelling` reads the reference's `phlyoP<N>way` key as
`phyloP<N>way` and `pc_description_text` reads the reference's `PC` header line with
vep-rs's description (the comparator applies the last two as well; both are on
`docs/intended-divergences.md`). A
flag set may carry its own `masks` for its run alone, kept when `flag-set` rewrites
the set's entry. The list is written into the manifest by hand (no `select` option
writes it), `select` keeps it when it regenerates the corpus, every subcommand
refuses a class outside `MASK_CLASSES` (those five and the classes of the corrected
behaviours, each on `docs/intended-divergences.md`; the golden harness and the
format comparator carry the same list), and `classify` applies the regulatory class
to the reference and compares `DOMAINS`, `FREQS` and `ZYG` as multisets under `order_only`
(the harness and the format comparator read `NearestExonJB` and the `5UTR_annotation`
groups that way as well).
`cargo test -p vep-cli --test golden --test format_parity` annotates each corpus
once per run and compares every column of every format (the Parquet round trip
reads the file back through the `parquet` crate, with no other tool). A failure lists every mismatch
grouped by column and consequence set. A record Ensembl VEP dropped before
annotation carries `reference_warning` (VEP's reason, from its warnings file,
recorded when the corpus is built); `classify` writes `reference_rows` 0 on it
and classes the keys only vep-rs emits for it `reference_skipped_record`, the
one class under which the harness accepts a record-level VCF line or JSON
object with no reference counterpart. A record Ensembl VEP's frequency filter
dropped (`--check_frequency`, `--filter_common`) carries no warning: `classify`
reads one neither engine wrote a row for as rowless by the run's nature, and
classes the keys only vep-rs emits for one it kept because an allele's
frequency is exactly 0 `frequency_filter_zero_value`. A manifest may also carry
`reference_warnings`, the texts Ensembl VEP wrote to its warnings file about
the run rather than about one record (the notice a `<CNV:TR>` record draws
without a FASTA, for one); the harness expects them in the run's warnings file
beside the numbered reasons. A flag set may carry its own list for its run
alone; one without inherits the manifest's. Like `masks`, the list is written
by hand, kept by `select` and by a `flag-set` rewrite of the set, and refused
by `build_golden_corpus.py` and the golden tests alike when a text is empty.

A `runs` entry that passes `--stats_text` may expect a `stats_file`, the
statistics report Ensembl VEP wrote for the same records
(`expected/<format>.summary.log`), which the harness compares with the run's as
`scripts/concordance/compare_vep_stats.py` does: the volatile rows by key alone,
the paths masked, every chart section as a set.

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
`--include sid=name` carries the record of dataset `sid` that `name` names, by its
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
request that adds a flag of vep-rs's own to `crates/vep-cli/src/args.rs` adds its row, one that adds
an Ensembl VEP flag lands the row the maintainers opened, and one that lands a row names its tests,
because CI runs `python3 scripts/parity/check_parity_matrix.py` (a flag without a row, a
landed row without its tests, or a named test or corpus that does not resolve fails; the plugin
support table of `docs/plugins.md` is held to the matrix's plugin rows, so a port's table status turns
supported in the pull request that lands its row) and `crates/vep-cli/tests/parity_matrix.rs`
(every long flag clap reports has a row).

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
lead-in or author handle; a trailing `(#N)` naming the pull request; at most 30 words in the bullet a
pull request adds, though a release's section may merge bullets into a longer one. A contributor's
entry takes its number once the pull request exists and is never predicted; a release's section is
numbered with the release's own pull request.

Refactors, CI and dependency bumps appear only when behaviour changes; the admission check asks every
crate code or manifest change for a changelog line, and a maintainer's `admission-override` label holds
the close for a refactor. Under a kind heading, entries may be grouped by area: a bullet `- <Area name>` at the
margin (the area alone, at most 60 characters) with its entries nested under it as `-` bullets indented
two spaces, each in the grammar above with its own `(#N)`; two levels at most, and every area has at
least one entry.

## Sign your commits

Every commit must carry a `Signed-off-by` trailer. Pass `-s` to `git commit`, which
appends `Signed-off-by: Your Name <your.email@example.com>` from your git identity:

```bash
git commit -s -m "Add CADD plugin support"
```

The sign-off asserts, under the [Developer Certificate of Origin](https://developercertificate.org),
that you have the right to submit the code under this project's Apache-2.0 license.
A `DCO` workflow (`.github/workflows/dco.yml`) checks every commit of a pull request
for the trailer and fails the check if one is missing. To add a missing sign-off before
the branch is pushed, use `git commit --amend -s` for the last commit or
`git rebase --signoff origin/main` for several.

A pushed branch with an unsigned commit is signed the same way on a new branch and opened
as a new pull request: never rebase or force-push a pushed branch, since maintainers may
push fixes onto it. Bring `main` in with `git fetch origin && git merge --signoff origin/main`;
from a fork clone, add this repository as the `upstream` remote, then
`git fetch upstream && git merge --signoff upstream/main`.

One commit of every pull request also declares the model that wrote the change, as a
second trailer beside the sign-off; the squash merge carries both trailers onto `main`.
The declaration takes one of two forms. The first names the model, any vendor's model
the same way: the first line below is the one a coding tool appends by default and the
second is its pattern. The second form, the last line, is for a change written without
a model:

```text
Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Co-Authored-By: <model name> <the vendor's noreply address>
Co-Authored-By: No LLM
```

The admission check reads a `Co-Authored-By` line as the declaration only when its name is
`No LLM`, or when it ends in an address whose mailbox is `noreply` or `no-reply` on a domain
other than `users.noreply.github.com`, where every address is a GitHub account's. Every
other such line credits a co-author and declares nothing.

## Issues

Every issue form carries a required choice about the intended divergences:
[docs/intended-divergences.md](docs/intended-divergences.md) lists every behaviour where vep-rs
deliberately differs from Ensembl VEP: Ensembl VEP defects it does not reproduce and differences in
what it writes or reads. A report or request touches an entry when it asks for a recorded behaviour
to change or go, argues against the entry's reasoning, or asks vep-rs to match Ensembl VEP where the
index says it deliberately does not; before reporting a difference, search the index for the flag,
field or shape. Say that the report touches none of the entries, or name the entries it touches
(the entry as the index's first column names it: the text of its row's first cell, inner backticks
dropped; with what should change and why the index's reasoning is wrong for this case), or propose
a new intended divergence: the input, what
Ensembl VEP and vep-rs each write, and why vep-rs's output is the one to keep. A proposal also states
the input's class with its basis, as `scripts/divergences/reproductions/README.md` defines the three:
`observed` in real data, `unobserved` there though well-formed, or `malformed`, with the rule of the
format's specification or of the documentation that the input breaks. Asking for an entry
to go needs the newer Ensembl VEP behaviour stated in plain words, with the record that shows it.
Additions are welcome; the maintainers write or review the entry's page. An issue opened or edited
without the choice, naming an entry the index lacks, naming existing entries with no entry line or
a line that says nothing after the entry, or proposing a divergence with an empty field, gets the
`needs-divergence-check` label and one comment; an edit of the description that supplies what was
missing removes the label and the comment reports the pass. An issue opened before the forms
carried the choice is not asked.

## Pull Requests

Before submitting a PR:

1. All CI checks pass on the toolchain `.github/workflows/ci.yml` pins: `cargo test --workspace --locked`, `cargo clippy --workspace --features vep-cli/bench -- -D warnings`, `cargo fmt --all -- --check`. The four tools under `tools/` are their own Cargo roots and run the same commands with `--manifest-path`; the Python tests under `scripts/` run with `python3 -m pytest -q <file>` after the `pip install` line the Comparator Tests job runs. The Advisory Reachability check holds the workspace's `Cargo.lock` with `scripts/check_advisory_reachability.sh`, which fails on an advisory reachable from the `vep` binary, and each tool's whole lockfile with `cargo audit --deny warnings --file tools/<name>/Cargo.lock`, which fails on any advisory.
2. Add tests for new features or defect fixes; the admission check below asks for a test change beside every code change. A test added to a file under `scripts/` raises that file's collected-count floor in `ci.yml`: name the new count under Verification for a maintainer to set.
3. Update the affected page under `docs/`: the format's page under `output/` for output layouts, the mechanism's page under `plugins/` and the plugin table in `plugins.md` for plugins, `cache-setup.md` for the cache builder and `cache-conversion.md` for the converter. A flag is documented by its help text in `crates/vep-cli/src/args.rs` and by the option tables, from which the CLI reference (`cli-reference.md` and the group pages under `cli/`) is generated: regenerate it with `VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference` and commit the result. The home page, `docs/index.md`, is written by `scripts/docs/render_home_page.py` from `scripts/docs/home_page.template.md` and the current release's records: edit the template, run the script and commit both. The plugin wall-time section of `docs/benchmarks/current-release.md` is written by `scripts/docs/render_plugin_wall_time.py` from its record under `docs/concordance-provenance/`: run the script and commit the page. A page added under `docs/` is listed in `docs/SUMMARY.md` under a title no other page has and opens with a paragraph that says what the page covers.
4. For consequence or plugin changes, run concordance smoke to confirm no regressions.

The description has the template's four sections: Summary, Verification, Upgrade impact and
Intended divergences. The last certifies the check of [Issues](#issues) for the change: one of its
three statements is checked, the first (the change touches none of the entries) excluding the other
two, and the second lists the entries touched one line each, the entry named in backticks as the
index's first column names it (the text of its row's first cell, inner backticks dropped) or by the
path of its page, and followed by what the change does to it and why the index's reasoning is wrong
for this case; the third describes a proposed entry. A change to the index's own text or order names
`docs/intended-divergences.md` itself under the second statement, a form only a change that edits
the index may use. A change that edits `docs/intended-divergences.md`, a page it links or the
reproduction registry under `scripts/divergences/reproductions/` (its README excepted) cannot claim
the first.
The section is the author's attestation, read by a maintainer: whether a code change that leaves
those files alone also leaves every recorded behaviour alone is not something the check can prove.

An admission check (the `Admission` workflow) runs on every pull request. It closes one that is
behind `main`, changes code without a test that gains a line (a benchmark or an example asks for
none), carries no model declaration, lacks a template section, or leaves the Intended divergences
section unchecked, inconsistent, naming an entry the index lacks, naming the index itself while it
leaves the index alone, or claiming none while it edits the index, a page it links or the
reproduction registry. It also closes one that changes a crate without a changelog line, changes a
file under `.github/workflows/`, fails the parity matrix or fails a required CI check. The comment
it leaves names each miss and its fix; push the fixes to the same branch and reopen it.
`python3 scripts/ci/admission.py --local --body-file <file>` applies the same rules to a saved
description before you open the pull request.

An unsigned commit is reported but never closes (the DCO check holds the merge; sign a pushed
branch by opening a new pull request from a signed branch). Maintainers' and bots' pull requests,
drafts and those opened before the check existed are commented but stay open, as is one opened
before the Intended divergences section was asked, for that miss alone; a maintainer's
`tests-not-needed` label waives the test requirement and `admission-override` holds the close for a
justified exception. The check judges a pull request by its author, so a maintainer sets
`admission-override` before pushing a workflow edit (a test floor in `ci.yml`) onto a contributor's
branch. `python3 scripts/ci/admission.py --local` applies the static rules to your branch before
you push (`--base upstream/main` from a fork clone).

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
- **Logging**: `tracing` for every diagnostic (`info!` for user-facing status, `debug!` behind `--verbose`, `warn!` for a recoverable condition). In the `vep` binary every `warn!` event reaches the warning sink (`crates/vep-cli/src/warnings.rs`), which writes it to the warnings file as `WARNING: <text>` and echoes it on stderr unless `--quiet`; a record the run skips goes through the sink's `skipped_variant`, never through `warn!`, so it is counted for `--skipped_variants_file`. The only direct writes are the sink's, the `vep` banner, the `INFO:` lines of the setup checks, the progress line and the end-of-run summary on stderr (in `main.rs`, `runner.rs` and `pipeline.rs`, suppressed by `--quiet`), the repeated `-o` warning the `vep` script prints before any option is read (`main.rs`, whatever `--quiet` says), the lines a plugin prints as its Ensembl VEP module's bare `warn` prints them, the transcript-set warning under `--refseq` or `--merged` and dbNSFP's column warning, file-name version warning and per-row `ERROR:` lines (`print_on_stderr` and `warn_on_stderr` in `crates/vep-builtins/src/traits.rs`, whatever `--quiet` says, never in the warnings file), the one-line `ERROR: [<component>] <message>` a refusal or a failed run ends with, `vep-cache-converter`'s mismatch listing, and the `bench_plugins` tool's progress and report output. The `vep` run's own writes go through `crates/vep-cli/src/stderr.rs`, where the `vep` binary switches on dropping a write that standard error refuses, so an unwritable stream never ends its run; `eprintln!` panics on one.
- **Output ordering**: consequence terms within a `TranscriptConsequence` are ordered by SO rank; Extra field keys follow VEP's flag-field order (`@FLAG_FIELDS`, the fields each active flag adds, in VEP's group order; IMPACT, DISTANCE, STRAND and FLAGS are always active), then every other key, plugin keys included, alphabetically; within a record, rows are ordered by feature kind (transcript, regulatory, motif; the transcript rows last on a record `--per_gene` or `--pick_allele_gene` reduced, `InputVariant::per_gene_pick`), then by Ensembl VEP's feature key (`TranscriptConsequence::ensembl_vep_key`: a transcript's stable id with `_` and its database id, a regulatory feature's stable id, a motif's database id), a structural variant's regulatory and motif rows instead in the order their features were met (`TranscriptConsequence::structural_feature_order`: the met rank, then the start and the stable id), then by ALT allele index, with one intergenic row per allele that overlapped no transcript last.
- **Naming**: the `vep` binary resolves each input chromosome name to the cache's before any lookup, as Ensembl VEP's `get_source_chr_name` does (the name itself when the cache holds it, else a synonym from `--synonyms` or the cache's `chr_synonyms.txt`, else the `chr` prefix removed or added: `chr21` -> `21`; `chrM` reaches `MT` through a synonym), and `original_chr` keeps the input spelling for the Location column, as VEP does; struct fields follow Ensembl VEP's names where possible (`transcript_id`, `gene_symbol`); SO terms are Ensembl VEP's exact strings (`5_prime_UTR_variant`, `NMD_transcript_variant`).
- **Plugin data**: `IndexMap<String, String>` keyed by field name, variant-level in `InputVariant.plugin_data` and transcript-level in `TranscriptConsequence.plugin_data` (dylib plugin results merge there, one Extra field per key of the returned JSON object); output order does not depend on the map: the default-format Extra column and the JSON writer list plugin keys alphabetically after the flag fields, and the tab, VCF and Parquet writers take their column order from the plugins' declared field lists.
- **Shared strings**: `Arc<str>` for frequently cloned string fields (stable ids, gene symbols); a clone is a refcount increment.
