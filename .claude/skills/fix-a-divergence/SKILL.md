---
name: fix-a-divergence
description: 'Resolve a disagreement between vep-rs and Ensembl VEP on one record: trace it to the engine that is wrong, fix vep-rs when it is, and document an Ensembl VEP defect as an intended divergence with its reproduction. Use when a concordance report, a golden-corpus mismatch or an "Output differs from Ensembl VEP" issue shows the two engines writing different rows.'
---

# Fix a divergence

The two engines may differ only where
[docs/intended-divergences.md](../../../docs/intended-divergences.md) says so; every other
difference is charged to vep-rs, whichever engine is wrong. A divergence is resolved by one of two
changes: a vep-rs defect is fixed, or an Ensembl VEP defect is corrected in vep-rs and documented as
an intended divergence. An Ensembl VEP defect is never reproduced to make the rows agree, and a
correctness fix is never held back for its speed cost. Branching, signing and the pull request
follow [AGENTS.md](../../../AGENTS.md); the Ensembl VEP source paths are in step 2 of
[close-parity-row](../close-parity-row/SKILL.md).

## 1. Reproduce it on one record

Take one record from a public dataset, or write one and print it in full. Run both engines on it
with the same cache, assembly, flags and FASTA, as
[docs/divergences/reproducing.md](../../../docs/divergences/reproducing.md) runs every example.
Ensembl VEP runs from the pinned image `ensemblorg/ensembl-vep:release_116.3`. Without Docker and an
Ensembl VEP 116 cache, start from a committed corpus instead: `cargo test -p vep-cli --test golden
--test format_parity` lists every mismatching column of every format over `tests/golden/`, and the
corpus's `expected/` directory is Ensembl VEP's output.

Check the index before classing the difference. Read
[docs/intended-divergences.md](../../../docs/intended-divergences.md) and search its sections and the
pages they link for the flag, the field or the shape of the row. A difference that is an entry there
is intended, and a change that makes vep-rs match Ensembl VEP on it is refused. A row that differs
from its entry's recorded shape is a `changed` verdict of `tools/vep-divergence-reevaluate`, which
re-traces the entry rather than opening a new one. Cite an entry as the index's first column names
it (the text of its row's first cell, inner backticks dropped) or by the path of its page. Removing
an entry needs the newer Ensembl VEP behaviour stated in plain words, with the record that shows it.

Done when you have both engines' rows for the record, the exact commands, the image digest and the
cache version, and either the index has no entry for the difference or you can cite the one it has.

## 2. Trace the mechanism

Read the Ensembl VEP code that writes its row and cite every claim as `Module.pm:line` at
`release/116.3`. Then derive the correct row from the record's own coordinates and the reference
sequence, not from either engine's output: which bases change, which transcript positions they
touch, which Sequence Ontology terms follow.

Done when you can name the engine that is wrong and the line that makes it so, and the derived row
equals one engine's row.

## 3. A vep-rs defect: fix it

Write the test first and see it fail on the wrong row. The test pins the literal expected row,
whole, taken from Ensembl VEP's output or derived in step 2, never from the function under test. A
change you cannot make fail before applying it is untested. Then fix the code until the test and the
golden harness pass. The committed `expected/` outputs are Ensembl VEP's and do not move for a
vep-rs fix.

Add a `CHANGELOG.md` line under Unreleased in CONTRIBUTING.md's grammar, under Fixed. Done when
`cargo test -p vep-cli --test golden --test format_parity` passes with no new mask.

## 4. An Ensembl VEP defect: document it

vep-rs writes the correct row and the difference becomes an entry of the index. One entry stands for
one mechanism, and the entry's name is the class the comparators count its rows under. The pull
request carries every part:

- **The class name**, in `CORRECTED_BEHAVIOUR_CLASSES` of `tools/vep-golden-corpus/src/constants.rs`
  and of `crates/vep-cli/tests/common/mod.rs`, and as a built-in class of
  `tools/vep-formats-compare`, each recognising the shape of the disagreement from
  the rows alone. The corpus that carries the record names the class in its `masks` list, and
  `vep-golden-corpus classify` records it on the explained rows as `expected_divergence`.
- **The page**: `docs/divergences/<entry-with-hyphens>.md` (the entry's name is written with
  underscores, its page's with hyphens, and the reproduction's `page` key names the page), or a
  section of its subject's page where one exists (a plugin's divergences live on the plugin's page).
  It states the difference in one sentence, what each engine does and where in Ensembl VEP's source
  its behaviour is, both engines' rows for the record, why vep-rs's row is the correct one from the
  record's own coordinates, and how a reader reproduces it. A paragraph stays under 120 words and a
  sentence under 45.
- **The index row**, in the section of `docs/intended-divergences.md` for the part of a run the
  entry touches, or under the same category of its last section of entries, `Differences only
  malformed input reaches`, when the entry's input class is `malformed`: the entry linked to its
  page, the difference in one sentence, its `Tier and count`, `Verified against` 116.3, and the
  page. The section's rules place the row: the tier from the input the divergence needs, then the
  counts. Count what the tree lets you count (the committed flag sets whose command line reaches
  the divergence, the corpora whose records carry the shape) and name the counts under
  Verification. The maintainers measure the rate over the datasets the page names and hold the
  order with `scripts/divergences/check_divergence_order.py`.
- **The reproduction**: `scripts/divergences/reproductions/<entry>.toml` in the shape the
  directory's README defines, `verified_against = "116.3"`, the image by digest, the input class
  with its basis (`observed`, `unobserved` or `malformed`, as that README defines them), the record,
  the committed cache, the compared field or rows, both engines' recorded values and the
  `implementation` table: the functions that hold the entry's behaviour, the tests that pin
  vep-rs's row and the comparison rules that recognise the difference, each as `<path>::<item>`.
  `cargo run
  --release --locked --manifest-path tools/vep-divergence-reevaluate/Cargo.toml -- --check` passes,
  and `cargo test --locked --manifest-path tools/vep-divergence-reevaluate/Cargo.toml` and `python3
  -m pytest -q scripts/divergences/test_check_divergence_order.py` stay green.
- **The changelog line** under Unreleased, naming the entry and the page.

Done when the harness passes with the class named, the index row and the reproduction agree, and
`--check` reports every entry verified against 116.3.

## 5. Gate

Commit with `-s` and the model declaration, then run the full set in
[AGENTS.md](../../../AGENTS.md#before-pushing) (the admission check included, on the saved
description). A consequence change also runs the concordance smoke (`just concordance-smoke
<perl_cache> <json_cache> <fasta>`); its F1 is not below `main`'s on the same input, or Verification
says the smoke was not run and why. Open the pull request with the template's four sections and
nothing else; Verification carries the record, both engines' rows, the citations and the commands.
Done when every command passes and its output is quoted.

The fourth section, Intended divergences, states truthfully which of its three statements holds. A
vep-rs fix of step 3 that leaves every entry alone checks the first. A change that alters or removes
an entry's recorded behaviour checks the second, one line per entry: the entry as the index names
it, what the change does to it, and why the index's reasoning is wrong for this case. The entry of
step 4 checks the third, with its input, both engines' rows and why vep-rs's row is the one to keep.
The first statement excludes the other two, and a change that edits the index, a page it links or
the reproduction registry under `scripts/divergences/reproductions/` (its README excepted) cannot
claim it; nor can a change to a line of an item a record's `implementation` table lists, which
names that record's entry under the second statement.
