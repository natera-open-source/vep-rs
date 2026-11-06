# Agent guide

vep-rs is a Rust port of Ensembl VEP release 116.2. Ensembl VEP's behaviour is the specification;
[CONTRIBUTING.md](CONTRIBUTING.md) is the source of truth for conventions, and this file caches the
parts that fail CI or review when missed.

## Where the truth lives

- **Spec**: Ensembl VEP's Perl source, read from the pinned image
  `ensemblorg/ensembl-vep:release_116.2` (Step 2 of the close-parity-row skill gives the paths).
  Cite what you match as `Module.pm:line`, in the test or comment that depends on it.
- **Backlog**: `tests/parity/matrix.tsv`, one row per Ensembl VEP flag, feature, plugin and cache
  field; [tests/parity/README.md](tests/parity/README.md) defines the columns. To close a row, follow
  [.claude/skills/close-parity-row/SKILL.md](.claude/skills/close-parity-row/SKILL.md).
- **Oracle**: the golden corpora under `tests/golden/` (Ensembl VEP's output, every column of every
  format) and the concordance harness in `scripts/concordance/`.

## Gotchas

- Sign every commit: `git commit -s`. The DCO check fails a pull request with an unsigned commit.
- Declare the model that wrote the change: one commit carries a trailer beside its sign-off, in
  one of the three forms below, and the squash body carries both trailers onto `main`. The first
  form is the line a coding tool appends by default; the last is for a change written without a
  model.

```text
Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Co-Authored-By: <model name> <the vendor's noreply address>
Co-Authored-By: No LLM
```

- A new `.rs` file opens with the two-line license header, a blank line, then its `//!` module doc
  (`scripts/check_module_docs.sh` checks the blank line).
- A flag added or changed in `crates/vep-cli/src/args.rs` gets its matrix row in the same change, and
  the generated CLI reference is rewritten and committed:
  `VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference`.
- A comment is one or two present-tense sentences stating a fact the code cannot, usually the Ensembl
  VEP behaviour matched, cited. It names nothing outside the code.
- Diagnostics go through `tracing`; a refusal ends the run with one `ERROR: [<component>] <message>`
  line. CONTRIBUTING.md's Code Conventions lists the few direct writes.
- `CHANGELOG.md` entries follow the grammar in CONTRIBUTING.md's Changelog entries. The trailing
  `(#N)` is added in a commit after the pull request exists.
- A change that can alter annotation output runs `cargo test -p vep-cli --test golden --test
  format_parity`. A consequence or plugin change also runs the concordance smoke (`just
  concordance-smoke <perl_cache> <json_cache> <fasta>`), which needs Docker and an Ensembl VEP 116
  cache; on a host without them, the pull request's Verification section says the smoke was not run
  and why.
- A test file under `scripts/` runs in CI only when `.github/workflows/ci.yml` lists it, with a floor
  on the tests it collects (pytest) or the `ok` lines it prints (bash). Adding or deleting a test
  moves that floor in the same change; a test file `ci.yml` does not list never runs.
- Before deleting or renaming a test, grep `tests/parity/matrix.tsv` for its name: a cell naming a
  test the tree lacks fails `check_parity_matrix.py`.
- A docs edit runs `python3 scripts/docs/site_metrics.py --gate` and, after `mdbook build`,
  `python3 scripts/docs/check_site_links.py` (CI adds `--tree HEAD`, which reads the committed tree,
  so run it that way after committing). A page over 4,000 words, a paragraph over 120, a sentence
  over 45, a table over 4 columns, a code line over 80 characters or a broken relative link fails the
  Pages job. It also runs for `README.md`, `CHANGELOG.md` and `CONTRIBUTING.md`, where the table,
  code-line and link checks apply; the page cap applies to all but `CHANGELOG.md`.
- A `Cargo.toml` change is checked with `cargo package --list -p <crate>` for every workspace crate,
  reading the list for what the package ships. `cargo package` takes the files under the crate's own
  directory plus the file its manifest's `license-file` names, so the repository's `LICENSE` and
  `NOTICE` ship only as a copy in the crate or, for `LICENSE`, through `license-file`. A scratch
  `cargo install --path crates/vep-cli --root <dir>` shows which binaries a user gets.

## Before pushing

Run these on CI's toolchain, Rust 1.97.1 (the `toolchain:` input in `.github/workflows/ci.yml`), and
quote the result in the pull request's Verification section. The MSRV job builds the workspace on
Rust 1.89 (`rust-version` in `Cargo.toml`); `rustup toolchain install 1.97.1 1.89` installs both.

```bash
cargo fmt --all -- --check && cargo clippy --workspace -- -D warnings
cargo test --workspace --locked --no-fail-fast
python3 scripts/parity/check_parity_matrix.py && scripts/check_module_docs.sh
cargo +1.89 build --workspace --locked
```

CI builds and tests with `--locked`: a dependency change commits `Cargo.lock` in the same change. A
dependency whose newer release declares a higher `rust-version` fails the MSRV job, so a bump is
built on 1.89 first. `tools/vep-fields-compare` is its own Cargo root: a change there runs the test,
clippy (with `--all-targets`) and fmt commands and the 1.89 build again with
`--manifest-path tools/vep-fields-compare/Cargo.toml`, and `scripts/check_module_docs.sh tools`.

A failure that reproduces on a clean checkout of `origin/main` is pre-existing. Name it under
Verification with that reproduction, and leave it out of this change. A newer local toolchain adds
clippy lints and can move type sizes, so compare like with like.

Cut the branch from the repository's freshly fetched `main`, and before pushing run
`python3 scripts/ci/admission.py --local` (from a fork clone, with the repository added as the
`upstream` remote, `--base upstream/main`). It applies the admission check's static rules to the
local branch against that base: current with `main`, a new or changed test with every code change,
every commit signed off, the model declared with a trailer form under Gotchas, the changelog line,
the parity matrix. CONTRIBUTING.md's Pull Requests says what the `Admission` check closes and how a
closed pull request reopens.

On a pull request from a fork, the CI, DCO, Pages and Docker workflows do not start until a
maintainer presses "Approve and run", after every push, and those checks read "no checks reported"
until then. Verification therefore carries the local results; the admission check runs at once.
`main` requires the branch to be current with it and every review thread resolved before merge.
Bring `main` in with `git fetch origin && git merge --signoff origin/main` (from a fork clone,
`git fetch upstream && git merge --signoff upstream/main`); never rebase or force-push a pushed
branch, since maintainers may push fixes onto it.

Fill in [.github/pull_request_template.md](.github/pull_request_template.md) for the description.

## Scope

- Before taking a row, a plugin or a performance item, ask on its issue whether a release branch
  already carries it; in-flight maintainer work closes a duplicate with credit rather than merging
  it.
- One matrix row or one issue per branch, cut from `origin/main`. Rows that share one mechanism, one
  message shape or one docs sentence (`gff` and `gtf`; the unported plugins) are a **family** and land
  together in one branch.
- A pull request is squash-merged: its body is written at merge from the commit messages, so wording
  fixes to commits need no force-push.
- A row's `target` is the maintainers' call. When the Perl source or the code shows the target is
  wrong, or another row contradicts it, stop and report it as a **target decision**, with the
  evidence, and keep the row as it is.
- Parallel branches collide on `tests/parity/matrix.tsv`, `CHANGELOG.md`,
  `crates/vep-cli/src/args.rs`, `crates/vep-cli/src/config.rs`, `docs/cli-reference.md` and
  `docs/cli/`. Rebase onto `origin/main` before the first push; once the pull request exists, merge
  `origin/main` in instead. Regenerate the CLI reference after either.
- Releases and dependency bumps stay with the maintainers: the version in `Cargo.toml`, tags,
  crates.io publishing, the `release.yml`, `docker.yml` and `pages.yml` workflows, and a bump in
  `Cargo.lock` that no change of yours needs.
