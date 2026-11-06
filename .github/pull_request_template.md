<!-- Title: imperative, plain English, no type prefix. Example: "Inflate bgzip input through libdeflate". -->
<!-- The body is these three sections and nothing else. A question for the maintainers goes on the issue; an idea beyond this change opens its own issue. -->

## Summary
<!-- One or two sentences: what changes for a user and why. Link an issue with "Closes #N" when one exists. Docs or data: which files change and what they now say or show. -->

## Verification
<!-- What was run and what it showed. Code: tests, clippy, fmt, and the concordance smoke whenever annotation output can change (state "byte-identical" or the F1). Docs or data: the source the text or figures were checked against (the script, the record). -->
<!-- A figure comes with its command, the image digest, the cache and the input. A speed claim names the machine, the input, the run count, the median and the spread, and shows output byte-identical to `main` in every format. -->
<!-- Without Docker and an Ensembl VEP 116 cache, say the smoke was not run and why; a maintainer runs it. -->

## Upgrade impact
<!-- "None", or the breaking change or changed default and what a user must do; on the release page this becomes the leading Breaking: entry of its area. -->

- [ ] Every commit is signed off (`git commit -s`); the DCO check fails without it.
- [ ] One commit declares the model that wrote the change, in a trailer in the form CONTRIBUTING.md's Sign your commits gives, or declares that none did.
- [ ] The branch contains `main`'s tip, brought in by a merge (never a rebase once pushed), and `python3 scripts/ci/admission.py --local` passes.
- [ ] The CI gate passed locally on the toolchain `ci.yml` pins (`cargo test --workspace --locked --no-fail-fast`, clippy with `-D warnings`, fmt, `check_parity_matrix.py`, `check_module_docs.sh`) and the MSRV build passed, with the output quoted under Verification, or no code changed.
- [ ] `CHANGELOG.md` has an entry under Unreleased in the release-notes grammar (one line, the subject first, at most 30 words, a trailing `(#N)` added once the PR number exists), or no crate code or manifest under `crates/` changed. The admission check asks every crate change for a line.
- [ ] Tests cover the change, added or updated in this PR (a `tests/` file, a `#[test]`, a `test_*.py` or `test_*.sh`, or the golden expectation or corpus `manifest.json` the change moves), or it changes no code.
- [ ] A test file under `scripts/` that gained or lost tests has its floor in `.github/workflows/ci.yml` moved to the new count, or no such file changed.
- [ ] The affected page under `docs/` is updated and `python3 scripts/docs/site_metrics.py --gate` passes, or no page is affected.
- [ ] Consequence or plugin change: concordance smoke passes (summary under Verification), or Verification says it was not run and why, or this PR does not touch consequences or plugins.
