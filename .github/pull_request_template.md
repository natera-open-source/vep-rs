<!-- Title: imperative, plain English, no type prefix. Example: "Inflate bgzip input through libdeflate". -->

## Summary
<!-- One or two sentences: what changes for a user and why. Link an issue with "Closes #N" when one exists. Docs or data: which files change and what they now say or show. -->

## Verification
<!-- What was run and what it showed. Code: tests, clippy, fmt, and the concordance smoke whenever annotation output can change (state "byte-identical" or the F1). Docs or data: the source the text or figures were checked against (the script, the record). -->

## Upgrade impact
<!-- "None", or the breaking change or changed default and what a user must do; this line becomes the release page's Upgrade notes bullet. -->

- [ ] Every commit is signed off (`git commit -s`); the DCO check fails without it.
- [ ] `CHANGELOG.md` has an entry under Unreleased in the release-notes grammar (one line, the subject first, at most 30 words, a trailing `(#N)` added once the PR number exists), or the change is not user-visible.
- [ ] Tests cover the change, or it changes no code.
- [ ] The affected page under `docs/` is updated, or no page is affected.
- [ ] Consequence or plugin change: concordance smoke passes (summary under Verification), or this PR does not touch consequences or plugins.
