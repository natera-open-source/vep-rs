# Parity matrix

`matrix.tsv` is the register of every Ensembl VEP parameter, feature, plugin and cache field that
vep-rs accounts for: one row per item, with the behaviour vep-rs targets for it and the tests that
prove it once it has landed. `scripts/parity/check_parity_matrix.py` holds the file to the tree on
every pull request.

## Columns

| Column           | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `item`           | The row's stable name: the `Config.pm` parameter name for a flag (`max_sv_size`, `phyloP`), the plugin or cache field name (`AlphaMissense`, `_gene_symbol`), or a fixed kebab-case slug for a behaviour (`symbolic-sv-classification`). A flag `--foo_bar` or `--foo-bar` in `crates/vep-cli/src/args.rs` is the row `foo_bar`; flags vep-rs adds beyond Ensembl VEP's have rows too.                                                       |
| `perl_reference` | The file of Ensembl VEP's source that implements the item. A module of `ensembl-vep` is its path below `modules/Bio/EnsEMBL/VEP/` (`Parser/VCF.pm`); a file under `modules/Bio/EnsEMBL/Variation/` of `ensembl-variation` or under `modules/Bio/EnsEMBL/IO/` of `ensembl-io` is its path below `modules/Bio/EnsEMBL/` (`Variation/Utils/Sequence.pm`, `IO/Parser/BaseVCF4.pm`); a script of `ensembl-vep` is its name (`vep`, `INSTALL.pl`, `filter_vep`). Each is a line of [`reference_files.txt`](#the-reference-files), which also lists the BioPerl package `Bio::Tools::CodonTable`. A plugin's module is `<Name>.pm`, `<Name>` a plugin of the support table, except LoFTEE's, which is `LoF.pm`. The consuming sub, or the parameter the module sets, follows `::` where one is the reference (`BaseVEP.pm::status_msg`, `VariantRecoder.pm::is_vr`). `-` when there is none (an item vep-rs alone has, or a distribution item). |
| `target`         | What vep-rs does for the item once the row has landed. `implemented`: the behaviour of Ensembl VEP. `documented_subset`: part of it, with the exclusion documented. `accepted_noop`: the flag parses and changes nothing but the command line the run records, as in Ensembl VEP where it reads nothing from the flag; `--debug` also prints one `INFO:` line. `error`: the flag or input is refused with a message naming it unsupported. `removed`: a flag or cache field vep-rs alone had, and drops. `excluded`: outside the port; such a row never lands. |
| `state`          | `open`, or `landed:<target>` naming the target the row landed at, which must equal `target`: a row that proved more or less than its target is changed by a target decision, not landed under another target; a bare `landed` is an error, and an `excluded` row never lands.                                                                                                                                                                |
| `test1`          | Test 1, the unit or behaviour test, as `<path>.rs::<fn>` under `crates/` (`crates/vep-effects/src/coding.rs::test_snv_stop_gained`); the function carries a test attribute (`#[test]`, `#[tokio::test]`, `#[rstest]`). A process row, which registers a release procedure rather than a behaviour of Ensembl VEP's, names a pytest function instead, as `scripts/<dir>/test_<name>.py::<fn>`, a named check of a shell suite, as `scripts/[<dir>/]test_<name>.sh::<check>`, or a test of a tool's own Cargo root, as `tools/<tool>/<path>.rs::<fn>`; its `perl_reference` is `-`. A shell check is named by the text the suite gives, whole and quoted, as the first argument of `ok` or of a function of its own that hands its first argument on to `ok`, in a command outside comments, quotes and here-documents; a name holding `$`, a backtick or a backslash is taken for one the suite computes and names no check, and neither does one holding `;`, which separates entries. Several entries are separated by `;`. |
| `test2`          | Test 2, the golden corpus whose records exercise the item, as `tests/golden/<set>/<corpus>`. The corpus's `manifest.json` carries a `covers` list naming every row its records exercise, and each of its `flag_sets` entries (a run of the corpus under other flags) carries one of its own; a row's `test2` resolves only when one of those lists names the row, so a corpus regenerated without the row's records cannot certify it, and every name in them is a row. |
| `interim`        | A refusal or stop-gap shipped before the row lands, as the `<path>.rs::<fn>` test that proves it; empty otherwise, and always empty once the row has landed.                                                                                                                                                                                                                                                                                 |

## How a row closes

A row moves from `open` to `landed:<target>` in the pull request that brings its behaviour to its
target, with the tests its target requires:

1. **Test 1**, in the crate the feature lives in. `implemented`: the flag's effect on a hand-built
   variant and transcript, with the expected value taken from Ensembl VEP's output on the same input
   and cited in the test. `error`: the error fires with its message. `removed`: the flag is rejected.
   `accepted_noop`: the output is the same with and without the flag, apart from the header lines
   that record the command line and the time of the run. `documented_subset`: the
   subset's effect, plus a test that the documented exclusion is stated. The check asks a landed
   `documented_subset` row for that test by what a test function reads: its body holds a string
   literal that is a Markdown file's path, with no white space in it and the file's name at its end
   (the page or the README whose sentences it asserts), or asks clap for a help text (`get_help`,
   `get_long_help`, `render_help`, `render_long_help`), or calls, by its bare name, a function of
   its own file whose body does. A message that names a page is no such path. The file a test
   stands in and the name it carries show nothing of what it reads.
2. **Test 2**, for `implemented` and `documented_subset` rows: a golden corpus under `tests/golden/`
   holding records that exercise the item, Ensembl VEP's output for them in all four formats, and a
   manifest whose `covers` list names the row, in the manifest's own list or in that of a `flag_sets`
   entry, a run of the same records under other flags with expected outputs of its own.
   `cargo test -p vep-cli --test golden --test format_parity` compares every column of every format
   of every run. `tools/vep-golden-corpus` writes the lists (`select --covers`, `covers
   --add`, `flag-set --covers`); `covers --add`, `flag-set --covers`, `select --matrix` and `check
   --matrix` refuse a name that is not a row.
3. **Test 3**, for the same two targets: a comparison against Ensembl VEP's output on the same input,
   cache and flags, every compared column agreeing exactly after the documented divergences are masked.

A process row (`flag_matrix_concordance`) registers a release procedure rather than a behaviour of
Ensembl VEP's: its `perl_reference` is `-`, its Test 1 is the test that proves the procedure's
mechanism, a pytest function named `scripts/<dir>/test_<name>.py::<fn>`, a check of a shell suite
named `scripts/[<dir>/]test_<name>.sh::<check>` or a test of a tool's own Cargo root named
`tools/<tool>/<path>.rs::<fn>`, and its Test 2 is a corpus whose `covers` names the row, as for any
other row.

## The reference files

`reference_files.txt` lists, one a line and in byte order, every file a `perl_reference` cell names
outside the plugins: the six scripts at the root of `ensembl-vep` and every file under its
`modules/Bio/EnsEMBL/VEP/` at the tag `release/116.3`, every file under
`modules/Bio/EnsEMBL/Variation/` of `ensembl-variation` and under `modules/Bio/EnsEMBL/IO/` of
`ensembl-io` at the commits the image `ensemblorg/ensembl-vep:release_116.3` holds, and the BioPerl
package that holds the genetic codes. A name that is no line of the list and no plugin's module is
refused, as `VCF.pm` is, the last component of `Parser/VCF.pm`. The list is written from a checkout
of each repository and moves with the reference release:

```bash
{
  git -C ensembl-vep ls-tree -r --name-only release/116.3 -- \
    INSTALL.pl convert_cache.pl filter_vep haplo variant_recoder vep \
    modules/Bio/EnsEMBL/VEP | sed 's#^modules/Bio/EnsEMBL/VEP/##'
  git -C ensembl-variation ls-tree -r --name-only \
    db7394a02523a6e519cb0d89f84b2889bac12737 -- \
    modules/Bio/EnsEMBL/Variation | sed 's#^modules/Bio/EnsEMBL/##'
  git -C ensembl-io ls-tree -r --name-only \
    6afb5dc27a5ae6881d75959153fbe6e9a4a7e788 -- \
    modules/Bio/EnsEMBL/IO | sed 's#^modules/Bio/EnsEMBL/##'
  echo 'Bio::Tools::CodonTable'
} | LC_ALL=C sort > tests/parity/reference_files.txt
```

## The plugin support table

A plugin row is one whose `perl_reference` is `<item>.pm`, the plugin's own module, a file of the
Ensembl VEP plugin collection unless the plugin's row of the support table says the collection lacks
it. Each has a row in the support table of `docs/plugins.md`, the six-column
table under its `## Support table` heading, whose first cell is the plugin name (linked or bare) and
whose fifth is one of three status phrases. The status fixes the matrix row's `target`: `supported,
verified on real data` and `supported on the user's licensed file, verified on synthetic fixtures` are
`implemented`, and a supported status also needs the row landed; `not ported` is `error`, and fits any
open row, a port decided but not landed included. Every table row names a matrix row by its `item`,
so a plugin whose module is named otherwise (LoFTEE, `LoF.pm`) is checked through its table row
alone, and the rows keep alphabetical order of the plugin name, case ignored. A port's table status
turns supported in the pull request that lands its row.

`python3 scripts/parity/check_parity_matrix.py` fails a pull request whose `args.rs` declares a flag
without a row, whose landed row lacks a required test or keeps an `interim`, whose `removed` row lands
while `args.rs` still declares the flag, whose process row names an Ensembl VEP module, whose
`perl_reference` is in none of the forms its column gives (a file the reference list lacks, a
`.pm` name that is no module of a plugin of the support table, or text after `::` that is not one
identifier), whose named test or corpus does not resolve, or whose corpus `covers` names something
that is not a row; it reads
the flags from the attributes of `args.rs`, and `cargo test -p vep-cli --test parity_matrix` asks clap
for the same list, hidden and flattened flags included. It also fails a pull request whose support
table names a plugin without a matrix row, gives a status that contradicts the row's `target` or a
supported status to a row still open, leaves alphabetical order, lacks a plugin row of the matrix, or
is malformed: a row without six cells, a status outside the three phrases, or a name given twice.
