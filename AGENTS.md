# Agent guide

vep-rs is a Rust port of Ensembl VEP release 116.3. Ensembl VEP's behaviour is the specification;
[CONTRIBUTING.md](CONTRIBUTING.md) is the source of truth for conventions, and this file caches the
parts that fail CI or review when missed.

## Where the truth lives

- **Spec**: Ensembl VEP's Perl source, read from the pinned image
  `ensemblorg/ensembl-vep:release_116.3` (step 2 of the close-parity-row skill gives the paths).
  Cite what you match as `Module.pm:line`, in the test or comment that depends on it.
- **Matrix**: `tests/parity/matrix.tsv`, one row per Ensembl VEP flag, feature, plugin and cache
  field; [tests/parity/README.md](tests/parity/README.md) defines the columns. Every row whose
  target vep-rs ports has landed; the rows still `open` are the `excluded` ones, which never land.
  New rows come from the maintainers, who open them when an Ensembl VEP release adds a parameter,
  plugin or cache field; a contributor who finds an item without a row proposes it on an issue with
  the parity-row form. [.claude/skills/close-parity-row/SKILL.md](.claude/skills/close-parity-row/SKILL.md)
  is the recipe for such a row once the maintainers have opened it.
- **Recipes**: the three kinds of change a pull request brings, each a recipe under
  `.claude/skills/`: [fix-a-divergence](.claude/skills/fix-a-divergence/SKILL.md) for a row where the
  engines disagree, [performance-change](.claude/skills/performance-change/SKILL.md) for a change
  made for speed or memory, and [port-a-plugin](.claude/skills/port-a-plugin/SKILL.md) for an
  Ensembl VEP plugin the support table lists as not ported.
- **Divergences**: [docs/intended-divergences.md](docs/intended-divergences.md) indexes every
  intended disagreement with Ensembl VEP, each with a page (or a section of its subject's page) and
  a reproduction under `scripts/divergences/reproductions/`, whose `implementation` table names the
  code, the tests and the comparison rules that implement the entry; `cargo run --release --locked
  --manifest-path tools/vep-divergence-reevaluate/Cargo.toml -- --check` holds the index and the
  reproductions to each other and to the reference release, each `implementation` table to the
  files and items the tree holds, and the index's `cache transcripts` cells to one total. A
  disagreement outside the index is charged to vep-rs, whichever engine is wrong.
- **Plugin support**: the support table of [docs/plugins.md](docs/plugins.md), one row per Ensembl
  VEP plugin with one of three statuses; `check_parity_matrix.py` holds it to the matrix's plugin
  rows, so a port's status turns supported in the pull request that lands its row.
- **Oracle**: the golden corpora under `tests/golden/` (Ensembl VEP's output, every column of every
  format, for human on both assemblies and for the other species
  [docs/cache-species.md](docs/cache-species.md) tables) and the concordance harness in
  `scripts/concordance/`.

## Gotchas

- Sign every commit: `git commit -s`. The DCO check fails a pull request with an unsigned commit.
- Declare the model that wrote the change: one commit carries a trailer beside its sign-off, in
  one of the two forms below, and the squash body carries both trailers onto `main`. The first
  form names the model: the first line is the one a coding tool appends by default and the
  second is its pattern for any vendor's model. The second form, the last line, is for a change
  written without a model.

  ```text
  Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
  Co-Authored-By: <model name> <the vendor's noreply address>
  Co-Authored-By: No LLM
  ```

  Any other `Co-Authored-By` line is a co-author's credit and declares nothing: the admission
  check reads the name `No LLM`, or a name before an address whose mailbox is `noreply` or
  `no-reply`, and a person's `<login>@users.noreply.github.com` is neither.

- A new `.rs` file opens with the two-line license header, a blank line, then its `//!` module doc
  (`scripts/check_module_docs.sh` checks the blank line).
- A doc comment links only to an item its built page can reach: a private item is named in a code
  span, a URL is written `<https://...>`, and a placeholder such as `<name>` goes in a code span.
  `scripts/check_api_docs.sh` builds the API documentation with rustdoc's warnings denied.
- A documented public item of a library crate that is added, removed or changed, a cargo feature
  or a default feature set that changes, and a library target a package gains fail the Module Docs
  job's `python3 scripts/check_public_interface.py --offline`, which names each one. Where the
  change is meant, write the record again with `python3 scripts/check_public_interface.py --update`,
  commit the diff of `tests/public-interface/` and name it in the pull request; an item no page
  offers takes `#[doc(hidden)]` or `pub(crate)` instead. A signature that names a type of a crate
  outside the workspace needs that crate's line, with its reason, in
  `tests/public-interface/outside-crates.tsv`.
- A flag added or changed in `crates/vep-cli/src/args.rs` has its matrix row in the same change (a
  flag of vep-rs's own gets a new row there; an Ensembl VEP flag lands the row the maintainers
  opened), and the generated CLI reference is rewritten and committed:
  `VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference`.
- A comment is one or two present-tense sentences stating a fact the code cannot, usually the Ensembl
  VEP behaviour matched, cited. It names nothing outside the code.
- Diagnostics go through `tracing`; a refusal ends the run with one `ERROR: [<component>] <message>`
  line. CONTRIBUTING.md's Code Conventions lists the few direct writes. A plugin or a `--custom`
  source names the input record a warning is about through the event field `input_record`, as in
  `warn!(input_record = variant.input_record, "...")`, so the warnings file lists the line with
  that record's warnings.
- `CHANGELOG.md` entries follow the grammar in CONTRIBUTING.md's Changelog entries. The trailing
  `(#N)` is added in a commit after the pull request exists.
- Before reporting or changing a difference from Ensembl VEP, read `docs/intended-divergences.md`
  and search its pages for the flag, field or shape: the entries are Ensembl VEP defects vep-rs
  deliberately does not reproduce and differences in what it writes or reads, and a change that
  makes vep-rs match one is refused. Cite the entry as the index's first column names it (the text
  of its row's first cell, inner backticks dropped) or by the path of its page; a change to the
  index's own text or order cites `docs/intended-divergences.md` itself. Removing an entry needs the
  newer Ensembl VEP behaviour stated in plain words, with the record that shows it. Fill the
  template's fourth section, Intended divergences, truthfully (the first statement excludes the
  other two), and run `python3 scripts/ci/admission.py --local --body-file <description>` on the
  saved description: the admission check closes a pull request whose section is missing, unchecked,
  inconsistent, names an entry the index lacks, cites the index itself while leaving it alone, or
  claims none while editing the index, a page it links or the reproduction registry under
  `scripts/divergences/reproductions/` (its README excepted), or while changing a line of an item
  a record there lists in its `implementation` table: such a change names the entry under the
  second statement and says what it does to it.
- A change that can alter annotation output runs `cargo test -p vep-cli --test golden --test
  format_parity`, which compares every column of every format over the committed corpora; the
  Parquet round trip reads the files back through the `parquet` crate, so no other tool is needed.
  A consequence or plugin change also runs the concordance smoke (`just concordance-smoke
  <perl_cache> <json_cache> <fasta>`), which needs Docker and an Ensembl VEP 116 cache; on a host
  without them, the pull request's Verification section says the smoke was not run and why.
- A test file under `scripts/` runs in CI only when `.github/workflows/ci.yml` lists it, with a floor
  on the tests it collects (pytest) or the `ok` lines it prints (bash). Adding or deleting a test
  moves that floor, and a new test file needs its line there: name the new count under Verification
  and a maintainer edits `ci.yml`, since a change to a workflow file from outside the maintainers
  closes the pull request. The counts are measured, never estimated:

  ```bash
  python3 -m pytest --collect-only -q <file> | grep -cE '^scripts/.*::'   # a pytest file's floor
  bash <test file> 2>&1 | grep -c '^  ok '                                  # a bash test file's floor
  ```

  A floor line reads `"<file> <count>"`; the one for this guide's own test is
  `"scripts/docs/test_contributor_guide.py 20"`, and that test holds the commands, versions and
  floors quoted on this page to what `ci.yml` states.
- Before deleting or renaming a test, grep `tests/parity/matrix.tsv` for its name: a cell naming a
  test the tree lacks fails `check_parity_matrix.py`.
- A test file under `crates/vep-cli/tests` starts `vep` through `common::isolated_vep` or
  `common::isolated_vep_anywhere`, which keep the machine's `VEP_*` variables and its `vep.ini`
  from the run, never by the binary's path (`CARGO_BIN_EXE_vep`). The test
  `no_other_test_file_starts_the_binary_by_its_path` of `cli_argv.rs` holds every other file of
  that directory to it.
- A docs edit runs `python3 scripts/docs/site_metrics.py --gate` and, after
  `python3 scripts/docs/site_build.py`, `python3 scripts/docs/check_site_links.py` (CI adds
  `--tree HEAD`, which reads the committed tree, so run it that way after committing).
  `site_build.py` runs `mdbook build` and exits 1 on a warning of the builder: mdBook warns of an
  HTML tag a page leaves unclosed and writes it into the page, so a placeholder such as `<dir>`
  stands in a code span that opens on the placeholder's own line. A page over 4,000 words (4,500 for
  `docs/intended-divergences.md`), a paragraph over 120, a sentence
  over 45, a table over 4 columns, a code line over 80 characters, a broken relative link, a
  redirect in `book.toml` to a page or heading the site lacks, a page under `docs/` the built site lacks (one with no line in `docs/SUMMARY.md`), or a
  repository link naming a missing path fails the Pages job. So does a chapter title another page has, or an opening paragraph that begins with the same 160 characters as another page's. A page added to `docs/SUMMARY.md` has a title no other page has and opens with a paragraph of 80 to 160 characters that says what the page covers; the site shows that paragraph as the page's description in a search result. `CHANGELOG.md` and `CONTRIBUTING.md`
  are site pages exempt from the paragraph and sentence caps (`CHANGELOG.md` from the page cap too);
  of `README.md` only the sections the site includes are measured, inside the pages that include
  them.
- A `Cargo.toml` change is checked with `cargo package --list -p <crate>` for every workspace crate,
  reading the list for what the package ships. `cargo package` takes the files under the crate's own
  directory plus the file its manifest's `license-file` names, so the repository's `LICENSE` and
  `NOTICE` ship only as a copy in the crate or, for `LICENSE`, through `license-file`. A scratch
  `cargo install --path crates/vep-cli --root <dir>` shows which binaries a user gets.
- A crate added to the workspace carries a `description` that names vep-rs and Ensembl VEP, with `keywords`, `categories`, `repository` and `homepage`; a crate with no library target also sets `documentation`, since docs.rs builds no page for it.
  A crate with a `tests/` directory carries the line `exclude = ["tests/"]`, as every such crate of
  the workspace does, so its package ships without that directory. A crate the workspace publishes
  joins the list `CRATES` of `.github/workflows/crates-io.yml` in dependency order, after every
  workspace crate it depends on. Its first publication is a maintainer's by hand, since crates.io
  trusted publishing cannot create a crate.

## Before pushing

Run these on CI's toolchain, Rust 1.97.1 (the `toolchain:` input in `.github/workflows/ci.yml`), and
quote the result in the pull request's Verification section. The MSRV job builds the workspace on
Rust 1.89 (`rust-version` in `Cargo.toml`); `rustup toolchain install 1.97.1 1.89` installs both.

```bash
cargo fmt --all -- --check && cargo clippy --workspace --features vep-cli/bench -- -D warnings
cargo test --workspace --locked --no-fail-fast
cargo test --locked -p vep-cli --features bench --bin bench_plugins
python3 scripts/parity/check_parity_matrix.py && scripts/check_module_docs.sh
scripts/check_api_docs.sh workspace
cargo +1.89 build --workspace --locked --features vep-cli/bench
```

CI builds and tests with `--locked`: a dependency change commits `Cargo.lock` in the same change. A
dependency whose newer release declares a higher `rust-version` fails the MSRV job, so a bump is
built on 1.89 first. The clippy feature `vep-cli/bench` puts the `bench_plugins` binary, which
builds only under that feature, under the lint. Every directory under `tools/` is its own Cargo
root: `tools/vep-cache-diff`, `tools/vep-cache-output-gates`, `tools/vep-cache-parity-audit`,
`tools/vep-discordant-classes`, `tools/vep-divergence-reevaluate`, `tools/vep-fastvep-normalize`,
`tools/vep-fastvep-plugin-tuples`, `tools/vep-fields-compare`, `tools/vep-file-census`,
`tools/vep-formats-compare`, `tools/vep-genome-features`, `tools/vep-genome-table`,
`tools/vep-golden-corpus`, `tools/vep-gtf-tags-compare`, `tools/vep-intermediates-compare`,
`tools/vep-loftee-synthetic-conservation`,
`tools/vep-outputs-compare`, `tools/vep-parquet-to-tab`, `tools/vep-plugin-pairs`,
`tools/vep-prepare-benchmark-vcfs`, `tools/vep-seeded-variants`, `tools/vep-stats-compare`,
`tools/vep-storable-digest`, `tools/vep-storable-digest-join`, `tools/vep-storable-dump`,
`tools/vep-storable-to-json`,
`tools/vep-storable-transcript-classes` and `tools/vep-sv-compare`. A change under any of them runs
the test, clippy (with `--all-targets`) and fmt commands and the 1.89 build again with
`--manifest-path tools/<name>/Cargo.toml`, and `scripts/check_module_docs.sh tools` and
`scripts/check_api_docs.sh tools`, the two commands the Module Docs job runs over the tool roots.
Some tools compile a source file of another by `#[path]`
(`grep -rnE '#\[path = "(\.\./)+vep-' tools/*/src` lists each), so a change to such a file runs the
commands of the tools that compile it as well. The Advisory Reachability job checks the workspace's
`Cargo.lock` with `scripts/check_advisory_reachability.sh`, which fails on an advisory reachable
from the `vep` binary unless its allowlist of dated entries
(`scripts/advisory-reachability-allow.txt`) carries one for it that has not expired, and each tool
root's whole lockfile with `cargo audit --deny warnings --file tools/<name>/Cargo.lock`, which
fails on any advisory but one its row names with `--ignore`, an advisory that allowlist carries for
the workspace.

The Python tests need `pip install pytest pyyaml`, and
`scripts/release/test_crate_publish_order.py` needs cargo on `PATH` as well, since one of its tests
reads the workspace through `cargo metadata`. Run `python3 -m pytest -q <file>` for
every `scripts/**/test_*.py` the change touches, and after `cargo test` run
`cargo test --locked --manifest-path tools/vep-parquet-to-tab/Cargo.toml` when the Parquet writer or
its reader changes, since its corpus round trips read the `vep` binary under `target/` and skip
without one, where CI holds them to run. `scripts/concordance/test_feature_oracle.sh` compares through the
release build of `tools/vep-formats-compare` and stops without it; build it with
`cargo build --release --locked --manifest-path tools/vep-formats-compare/Cargo.toml`. The tests
of `scripts/concordance/test_reproduce_plugin_concordance.py` that score run a build of
`tools/vep-plugin-pairs`, release or debug, and skip without one. The tests of
`tools/vep-golden-corpus` that classify with SpliceAI score files ask a build of
`tools/vep-formats-compare`, release or debug, and fail without one. The checks
of the bash test files under `scripts/concordance/` and `scripts/validation/` that run a script on
the release build of a tool it calls skip without that build: `tools/vep-outputs-compare`,
`tools/vep-fields-compare`, `tools/vep-fastvep-normalize`, `tools/vep-file-census`,
`tools/vep-cache-diff` and `tools/vep-cache-output-gates`; build each with the same command and its
own manifest. A bash test file under `scripts/concordance/`,
`scripts/docker/`, `scripts/validation/` or `scripts/data/` runs as `bash <test file>` and exits 1
on any failing check. The floors are read
from `ci.yml` itself: the Comparator Tests job holds each file of its `tests=(...)` list to the
floor beside it before running them all, and the Harness Tests job holds each of its bash test files
to the floor of its `check <test file> <floor>` line (the plugin data layout test file is held the
same way inside the Comparator Tests job).

A failure that reproduces on a clean checkout of `origin/main` is pre-existing. Name it under
Verification with that reproduction, and leave it out of this change. A newer local toolchain adds
clippy lints and can move type sizes, so compare like with like.

Cut the branch from the repository's freshly fetched `main`, and before pushing run
`python3 scripts/ci/admission.py --local --body-file <description>` on the saved description (from
a fork clone, with the repository added as the `upstream` remote, `--base upstream/main`). It
applies the admission check's static rules to the local branch against that base: current with
`main`, a test gaining a line beside every code change, every commit signed off, the model declared
with a trailer form under Gotchas, the changelog line, the parity matrix, and the description's
Intended divergences section (without `--body-file` the run certifies none on its own, which a
change to an item a reproduction's `implementation` table lists then fails).
CONTRIBUTING.md's Pull Requests says what the `Admission` check closes and how a closed pull request
reopens.

On a pull request from a fork, the CI and DCO workflows, and Pages or Docker when their paths
change, do not start until a maintainer presses "Approve and run", after every push, and those
checks read "no checks reported" until then. Verification therefore carries the local results; the
admission check runs at once. `main` requires the branch to be current with it and every review
thread resolved before merge. Bring `main` in with
`git fetch origin && git merge --signoff origin/main` (from a fork clone,
`git fetch upstream && git merge --signoff upstream/main`); never rebase or force-push a pushed
branch, since maintainers may push fixes onto it.

Fill in [.github/pull_request_template.md](.github/pull_request_template.md) for the description:
its four sections, the fourth being the Intended divergences statement under Gotchas, and nothing
else.

## Scope

- A row for an Ensembl VEP item is opened by the maintainers, never by the pull request that lands
  it: when an Ensembl VEP release adds a parameter, plugin or cache field, they open its row with
  the target it gets; a contributor who finds such an item without a row proposes it on an issue
  with the parity-row form, and closes the row only when it is open. A flag vep-rs adds beyond
  Ensembl VEP's gets its row in the change that declares it.
- Before taking a divergence, a plugin or a performance item, ask on its issue whether a release
  branch already carries it; in-flight maintainer work closes a duplicate with credit rather than
  merging it.
- One issue per branch, cut from `origin/main`: one divergence entry, one plugin, one performance
  change. Items that share one mechanism, one message shape or one docs sentence (`gff` and `gtf`; two
  plugins of one data layout) are a **family** and land together in one branch.
- A pull request is squash-merged: its body is written at merge from the commit messages, so wording
  fixes to commits need no force-push.
- A row's `target` is the maintainers' call. When the Perl source or the code shows the target is
  wrong, or another row contradicts it, stop and report it as a **target decision**, with the
  evidence, and keep the row as it is.
- Parallel branches collide on `tests/parity/matrix.tsv`, `CHANGELOG.md`,
  `crates/vep-cli/src/args.rs`, `crates/vep-cli/src/config.rs`, `docs/cli-reference.md`,
  `docs/cli/`, `docs/intended-divergences.md` and the support table of `docs/plugins.md`. Rebase onto
  `origin/main` before the first push; once the pull request exists, merge `origin/main` in instead.
  Regenerate the CLI reference after either.
- Releases and dependency bumps stay with the maintainers: the version in `Cargo.toml`, tags,
  crates.io publishing, every workflow under `.github/workflows/`, and a bump in `Cargo.lock` that
  no change of yours needs.
