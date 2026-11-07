---
name: close-parity-row
description: Close a row of tests/parity/matrix.tsv, moving it from `open` to `landed:<target>` with the tests its target requires. Use when porting an Ensembl VEP flag, plugin, input format or cache field, or when picking up a parity-row issue.
---

# Close a parity row

A row lands in the pull request that brings its behaviour to its `target`, with the tests the target
requires. [tests/parity/README.md](../../../tests/parity/README.md) defines the columns and the tests
each target needs; read it before step 1. Branching, signing and the pull request follow
[AGENTS.md](../../../AGENTS.md).

## 1. Read the row

```bash
awk -F'\t' 'NR==1 || $1=="<item>"' tests/parity/matrix.tsv
```

Done when you can state the row's `target`, its `perl_reference` and its `interim` test, if it has one.

## 2. Read the reference

Open the `perl_reference` module and find the code that implements the item. Read it from GitHub at
the pinned tag, `https://raw.githubusercontent.com/Ensembl/ensembl-vep/release/116.2/<path>`, or from
the image with `docker run --rm ensemblorg/ensembl-vep:release_116.2 cat /opt/vep/src/ensembl-vep/<path>`:

| `perl_reference` | `<path>` |
| --- | --- |
| a VEP module (`Config.pm`, `Parser/VCF.pm`) | `modules/Bio/EnsEMBL/VEP/<module>` |
| `vep`, `INSTALL.pl` | the name itself |
| `Variation/<module>` | `Bio/EnsEMBL/Variation/<module>` in the image; in `Ensembl/ensembl-variation` at `modules/Bio/EnsEMBL/Variation/<module>`, branch `release/116` |
| a plugin (`<Name>.pm`) | `<Name>.pm` in `Ensembl/VEP_plugins`, branch `release/116` |

Write down the behaviour vep-rs must match, each claim cited as `Module.pm:line`: for an `error` or
`accepted_noop` row, where Ensembl VEP declares the flag and what it does with it; for the other
targets, every behaviour the target covers.

When the Perl source shows the `target` is wrong (an `error` row Ensembl VEP's own behaviour makes an
`accepted_noop`, an `implemented` row that needs the Ensembl database), stop: report a **target
decision** to the maintainers with the citations, and leave the row unchanged.

Before classing a difference between the two engines as a defect of either, search
`docs/intended-divergences.md` and its pages for the flag, field or shape: a difference the index
records is deliberate, and a row that would make vep-rs match Ensembl VEP there is a target decision,
not a port. Cite the entry as the index's first column names it (the text of its row's first cell,
inner backticks dropped) in the pull request's fourth section; asking
for an entry to go needs the newer Ensembl VEP behaviour stated in plain words, with the record that
shows it; a new difference the row uncovers is proposed there, with the input, both engines' rows
and why vep-rs's row is the one to keep.

Done when every claim has a citation you read yourself at `release/116.2`.

## 3. Write Test 1 red

The test lives in the crate the behaviour lives in. Follow the pattern already used for the target:

| `target` | Test 1 | Pattern |
| --- | --- | --- |
| `error` | The refusal fires with its message, exit code 1 | `assert_refused` in `crates/vep-cli/tests/cli_errors.rs` |
| `accepted_noop` | Output is byte-identical with and without the flag | `assert_no_effect` in `crates/vep-cli/tests/cli_noops.rs` |
| `removed` | The flag is rejected and `args.rs` does not declare it | `crates/vep-cli/tests/cli_errors.rs` |
| `implemented` | The effect on a hand-built variant and transcript, the expected value taken from Ensembl VEP's output on the same input and cited | the owning crate's unit tests |
| `documented_subset` | The subset's effect, plus a test that the page states the exclusion | the `implemented` pattern, plus `crates/vep-cli/tests/documented_subsets.rs` |

A test pins the literal expected value, the whole row, taken from Ensembl VEP's output or the built
artefact, never from the function under test. A test that asserts a substring, or a value the code
computed, proves nothing.

Done when the test fails for the reason the row describes, and you have seen it fail.

## 4. Make it green

Implement until Test 1 passes.

- **Command-line refusals and no-ops**: a flag refused, or accepted without effect, is an entry of
  `ERROR_TABLE` in `crates/vep-cli/src/config.rs`. Its doc comments and the existing entries are the
  contract, and `every_error_table_entry_names_flag_need_and_route` in `cli_errors.rs` checks the
  message's shape. The refusal message, help text and CLI reference all read the entry. A flag in
  Ensembl VEP's `@VEP_PARAMS` takes `refused(...)` or `noop(...)`; a flag outside it takes the
  inline `declared: false` form, caught before parsing.
- **Flags**: a flag in `@VEP_PARAMS` is declared in `crates/vep-cli/src/args.rs`. A refused or no-op
  flag carries no `///` doc comment, because its help text comes from the error table. Give it the
  `help_heading` of the vep-rs flag it sits beside in `Config.pm` (the `--custom` shortcuts beside
  `--custom`), because the heading picks the `docs/cli/` page it renders on. After any `args.rs` or
  error-table change, regenerate the CLI reference: `VEP_UPDATE_CLI_REFERENCE=1 cargo test -p
  vep-cli --test cli_reference`.
- **Examples that name the flag**: tests and docs use undeclared Ensembl VEP flags as examples of
  "an option this version lacks". Declaring the flag breaks them. Find them with `grep -rnw --
  '--<flag>' crates docs` and `grep -rn '"<flag>"' crates`, move each to another flag that is still
  undeclared, and name the new example in the pull request.
- **The interim test**: once the real behaviour replaces the stop-gap, delete the test the row's
  `interim` names, and fix every doc comment, module doc or intra-doc link that names it, including
  in tests that belong to other rows. Before deleting or renaming any test, grep
  `tests/parity/matrix.tsv` for its name: another row's cell may name it.
- **New row kinds**: a new kind of output row (a feature type, an intergenic or companion row) is
  driven through every writer, the `--pick` family, the row filters and the plugin hosts, each
  compared whole-row with Ensembl VEP's output for the same records. A row selector or mask is a
  reporting aid, never the gate: the gate is every row of every format.

Done when Test 1 passes and `cargo test -p <crate>` is green.

## 5. Test 2 and Test 3

`implemented` and `documented_subset` rows need a golden corpus that covers the row and a comparison
against Ensembl VEP: follow [corpus.md](corpus.md). Other targets skip this step.

## 6. Land the row

Edit the row:

- `state`: `landed:<target>`.
- `test1`: the tests written for this row, as `<path>.rs::<fn>` separated by `;`.
- `test2`: the corpus, for `implemented` and `documented_subset`; empty for the other targets.
- `interim`: empty.

Find every hand-written page that states the item's old behaviour with `grep -rnw -- '<flag or
item>' docs/` (`-w` keeps `--gff` from matching `--gff3`) and bring it in line; `docs/cli-reference.md` and `docs/cli/` are generated in step 4. A
`documented_subset` exclusion is stated on its page in the sentences Test 1 asserts. Add a
`CHANGELOG.md` entry under Unreleased when the change is visible to a user; the admission check asks
every crate code or manifest change for a changelog line, and a maintainer's `admission-override`
label holds the close for a refactor.

Done when `python3 scripts/parity/check_parity_matrix.py` and `cargo test -p vep-cli --test
parity_matrix` pass.

## 7. Gate

Commit with `-s` and the model declaration, run the full set in
[AGENTS.md](../../../AGENTS.md#before-pushing) (the admission check included, on the saved
description with `--body-file`), then open the pull request with the template's four sections and
nothing else, the Intended divergences statement checked truthfully. A test added to a file under `scripts/` moves that file's floor in
`.github/workflows/ci.yml`: name the new count under Verification for a maintainer to set. Done when every command passes, its output is quoted in the pull
request's Verification section, and the smoke's F1 is quoted there or stated as not run and why.
