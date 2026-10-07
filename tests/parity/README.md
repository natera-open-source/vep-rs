# Parity matrix

`matrix.tsv` is the register of every Ensembl VEP parameter, feature, plugin and cache field that
vep-rs accounts for: one row per item, with the behaviour vep-rs targets for it and the tests that
prove it once it has landed. `scripts/parity/check_parity_matrix.py` holds the file to the tree on
every pull request.

## Columns

| Column           | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `item`           | The row's stable name: the `Config.pm` parameter name for a flag (`max_sv_size`, `phyloP`), the plugin or cache field name (`AlphaMissense`, `_gene_symbol`), or a fixed kebab-case slug for a behaviour (`symbolic-sv-classification`). A flag `--foo_bar` or `--foo-bar` in `crates/vep-cli/src/args.rs` is the row `foo_bar`; flags vep-rs adds beyond Ensembl VEP's have rows too.                                                       |
| `perl_reference` | The Ensembl VEP module that implements the item, as a path below `modules/Bio/EnsEMBL/VEP/` (ensembl-variation modules under `Variation/`), with the consuming sub after `::` where one is the reference (`BaseVEP.pm::status_msg`); `vep` and `INSTALL.pl` name the scripts; `-` when there is none (an item vep-rs alone has, or a distribution item).                                                                                     |
| `target`         | What vep-rs does for the item once the row has landed. `implemented`: the behaviour of Ensembl VEP. `documented_subset`: part of it, with the exclusion documented. `accepted_noop`: the flag parses and changes nothing, as in Ensembl VEP. `error`: the flag or input is refused with a message naming it unsupported. `removed`: a flag or cache field vep-rs alone had, and drops. `excluded`: outside the port; such a row never lands. |
| `state`          | `open`, or `landed:<target>` naming the target the row landed at, which must equal `target`: a row that proved more or less than its target is changed by a target decision, not landed under another target; a bare `landed` is an error, and an `excluded` row never lands.                                                                                                                                                                |
| `test1`          | Test 1, the unit or behaviour test, as `<path>.rs::<fn>` under `crates/` (`crates/vep-effects/src/coding.rs::test_snv_stop_gained`); the function carries a test attribute (`#[test]`, `#[tokio::test]`, `#[rstest]`). Several entries are separated by `;`.                                                                                                                                                                                 |
| `test2`          | Test 2, the golden corpus whose records exercise the item, as `tests/golden/<set>/<corpus>`. The corpus's `manifest.json` carries a `covers` list naming every row its records exercise; a row's `test2` resolves only when that list names the row, so a corpus regenerated without the row's records cannot certify it.                                                                                                                    |
| `interim`        | A refusal or stop-gap shipped before the row lands, as the `<path>.rs::<fn>` test that proves it; empty otherwise, and always empty once the row has landed.                                                                                                                                                                                                                                                                                 |

## How a row closes

A row moves from `open` to `landed:<target>` in the pull request that brings its behaviour to its
target, with the tests its target requires:

1. **Test 1**, in the crate the feature lives in. `implemented`: the flag's effect on a hand-built
   variant and transcript, with the expected value taken from Ensembl VEP's output on the same input
   and cited in the test. `error`: the error fires with its message. `removed`: the flag is rejected.
   `accepted_noop`: the output is byte-identical with and without the flag. `documented_subset`: the
   subset's effect, plus a test that the documented exclusion is stated.
2. **Test 2**, for `implemented` and `documented_subset` rows: a golden corpus under `tests/golden/`
   holding records that exercise the item, Ensembl VEP's output for them in all four formats, and a
   manifest whose `covers` list names the row. `cargo test -p vep-cli --test golden --test format_parity`
   compares every column of every format.
3. **Test 3**, for the same two targets: a comparison against Ensembl VEP's output on the same input,
   cache and flags, every compared column agreeing exactly after the documented divergences are masked.

`python3 scripts/parity/check_parity_matrix.py` fails a pull request whose `args.rs` declares a flag
without a row, whose landed row lacks a required test or keeps an `interim`, whose `removed` row lands
while `args.rs` still declares the flag, or whose named test or corpus does not resolve; it reads the
flags from the attributes of `args.rs`, and `cargo test -p vep-cli --test parity_matrix` asks clap for
the same list, hidden and flattened flags included.
