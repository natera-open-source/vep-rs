# Test 2 and Test 3: golden corpus and concordance

Reached from step 5 of [SKILL.md](SKILL.md), for `implemented` and `documented_subset` rows. The
corpus format, its masks and its size budget are defined in CONTRIBUTING.md's Golden corpora and
output parity section; this file is the order of operations.

## Choose the corpus

Use the narrowest change that gives the row records exercising it, in this order:

1. **An existing run already exercises it.** The records and flags of a corpus under
   `tests/golden/116/` already produce the item's output. Add the row to its `covers`:
   `scripts/golden/build_golden_corpus.py covers <corpus> --add <item>`, with `--flag-set <name>` for
   a flag set's list.
2. **The records exist but the flags differ.** Add a flag set to the corpus whose records exercise the
   item: `build_golden_corpus.py flag-set <corpus> --name <name> --flag=--<flag> --covers <item>`.
   Its expected outputs need the reference run below.
3. **No corpus has the records.** Select a new corpus (`build_golden_corpus.py select`, the recipe in
   CONTRIBUTING.md). This needs the reference run and the Ensembl VEP outputs the selection reads.

Done when `build_golden_corpus.py covers <corpus> --list` (with `--flag-set <name>` for a set) names
the row and `python3 scripts/parity/check_parity_matrix.py` passes.

## Reference run (options 2 and 3)

Ensembl VEP's expected outputs need a host with Docker and an Ensembl VEP 116 cache for the corpus's
assembly. `scripts/concordance/generate_reference_output.sh --corpus <dir>` writes a corpus's
top-level `expected/` and its `provenance.json`. A flag set's `expected_<name>/` has no generator in
the repository: run Ensembl VEP under the set's `flags` the way `provenance.json` records the
top-level run, one run per format, and say in the pull request how the files were made. Then
`build_golden_corpus.py classify` (with `--flag-set <name>` for a flag set) records the divergence
keys against the vep-rs binary.

On a host without Docker or the cache, the reference run is **blocked**. Write Test 1 and the
implementation, leave the row `open`, and open a draft pull request whose Verification section names
the corpus, the flags and the reference run it is waiting for.

## Compare

Run `cargo test -p vep-cli --test golden --test format_parity`. Done when every column of every
format of every run of the corpus passes. A mismatch is fixed in the code or recorded as a documented
divergence under `docs/divergences/`; the expected outputs stay as Ensembl VEP wrote them.

Keep the corpus under its 15 MB budget (`each_corpus_fits_its_size_budget`): narrow the exemplars
with `--max-span` before growing the cache.

## Test 3

For a consequence or plugin change, run the concordance smoke (`just concordance-smoke <perl_cache>
<json_cache> <fasta>`) and quote its F1 summary in the pull request, with the command, the image
digest, the cache and the input it ran on. Done when the smoke's F1 is not below `main`'s on the
same input, or the pull request states that it was not run and why.
