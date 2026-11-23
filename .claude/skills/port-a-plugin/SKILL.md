---
name: port-a-plugin
description: Port an Ensembl VEP plugin the support table of docs/plugins.md lists as not ported, as a built-in of the vep-builtins crate matched field for field to the plugin's own Perl module, with its fixtures, golden corpus and concordance case. Use when picking up a plugin request or a plugin's matrix row.
---

# Port a plugin

A plugin's specification is its own Ensembl VEP module, `<Name>.pm` in the Ensembl VEP plugin
collection at `release/116`: every field name, header description, parameter, allele-matching rule
and warning vep-rs prints is that module's, cited as `Name.pm:line`. The port is a built-in of the
`vep-builtins` crate
([docs/developers/writing-a-plugin.md](../../../docs/developers/writing-a-plugin.md), Adding a new
plugin), and it lands in one pull request with its matrix row, its support-table row, its fixtures,
its golden corpus and its page. Branching, signing and the pull request follow
[AGENTS.md](../../../AGENTS.md); the Ensembl VEP source paths are in step 2 of
[close-parity-row](../close-parity-row/SKILL.md).

## 1. Ask first

Open a Plugin request, or comment on the one that exists, before writing code: a release branch may
carry the plugin, and in-flight maintainer work closes a duplicate with credit rather than merging
it. The plugin's matrix row has `target` `error` while the support table says `not ported`, and a
row's `target` is the maintainers' call. The port therefore needs their word on the issue that the
target moves to `implemented`, answering the reason the table's Reason cell gives for the plugin
being unported (a GRCh37-only file, a web service, a database, a licence).

Done when a maintainer has agreed the target on the issue.

## 2. Read the module

Read `<Name>.pm` whole and write down, each with its line:

- the parameters `new` reads and their grammar;
- the data file it opens and how it is indexed, and the released version of that file the module
  documents;
- `get_header_info`, the field names and descriptions as printed;
- `variant_feature_types` and `feature_types`, which decide the row kinds the plugin writes on;
- how `run` matches a record to a data row (position, alleles, strand, transcript or gene, amino
  acids), and what it writes for a record with no match;
- every warning the module prints.

Done when every behaviour has a citation you read yourself, and you can state the rows the plugin
writes on and the fields it writes.

## 3. The fixture

The corpus commits the data the plugin reads, so both engines read one copy of the same bytes;
[scripts/fixtures/README.md](../../../scripts/fixtures/README.md) decides the kind by the data's
terms. A file that may be redistributed is a cut of the real file over the corpus's regions
(`scripts/fixtures/cut_plugin_regions.py`). A file that is not redistributed here is a synthetic
file with fabricated scores at the corpus's own sites, written by a seeded
script beside the existing ones (`plugin_score_fixtures.py`, `make_plugin_fixture.py`). The support
table's status for such a plugin is `supported on the user's licensed file, verified on synthetic
fixtures`. The corpus manifest names the files under `plugin_files`, and its flags reach them through
`{corpus}` and `{plugins}`.

Done when the fixture follows the layout the Ensembl VEP plugin's source and documentation state, in its preamble lines, header and column order,
and a test of the script
that writes it (a `scripts/fixtures/test_*.py`, floored in `ci.yml`) covers the layout.

## 4. Write the tests red

Unit tests in the plugin's own module, each pinning the literal field values Ensembl VEP's plugin
writes for a record of the fixture, whole, and cited. They cover the matched row, the unmatched
record, each parameter's effect and the refusal of a malformed parameter: the built-in grammar is
all positional or all `key=value`, and a mixed, repeated or undeclared parameter is refused at
initialisation naming it. `crates/vep-cli/tests/plugin_headers.rs` loads a stack of the built-ins on
committed data and holds each one's header descriptions to the three headed formats: add the port to
that stack with its fixture, so its `header_info` is held to what `get_header_info` prints.

Done when the tests fail for want of the port, and you have seen them fail.

## 5. Make it green

Implement the module under `crates/vep-builtins/src/plugins/`, register it in `create_builtin()` of
`crates/vep-builtins/src/registry.rs` and remove its entry from `UNPORTED_ENSEMBL_PLUGINS` there:
the refusal test (`crates/vep-cli/tests/cli_plugins_unported.rs`) and the partition test
(`unported_and_builtin_names_partition_the_collection`) hold the two lists to the collection. A
tabix-backed plugin implements the two-phase `prefetch` and `annotate` protocol the page describes
and reads through `open_best_store_with`, so a `.vpd` store serves it as the tabix file does.

Done when the unit tests pass and `cargo test -p vep-builtins -p vep-cli` is green.

## 6. The corpus and the concordance case

The golden corpus is `tests/golden/116/<assembly>-plugin-<name>` (a corpus over synthetic data
takes a `-fixture` suffix), with the plugin's flags as a flag set and the fixture under
`plugin_data/`. Ensembl VEP's expected outputs come from
`scripts/concordance/generate_reference_output.sh --corpus <dir> --flag-set <name> --dir-plugins
<modules>`, which needs Docker and an Ensembl VEP 116 cache; follow
[corpus.md](../close-parity-row/corpus.md) for the order of operations. On a host without them, open
the pull request as a draft whose Verification names the reference run it waits for. `cargo test -p
vep-cli --test golden --test format_parity` then compares every field of every format.

A plugin whose open data file the harness prepares has a row in `scripts/data/plugin_data_layout.sh`
(`scripts/data/test_plugin_data_layout.sh` proves it, with its floor in `ci.yml`) and runs through
`scripts/concordance/run_concordance.sh --plugins <Name>`, whose `vep-plugin-pairs compare` report
is the concordance case; the pull request quotes it or states that it was not run and why.

Done when the golden tests pass over the new corpus and the matrix row names its tests.

## 7. Land the row and the pages

Edit the matrix row: `target` `implemented`, `state` `landed:implemented`, `test1` the unit tests as
`<path>.rs::<fn>` separated by `;`, `test2` the corpus. Edit the plugin's row of the support table
of `docs/plugins.md`: the status phrase (`supported, verified on real data`, or the synthetic-fixtures
phrase of step 3), the data file and provider, the terms and the assemblies, the rows kept in
alphabetical order. `python3 scripts/parity/check_parity_matrix.py` holds the table to the row.

Write the plugin's section on the mechanism's page under `docs/plugins/` (a tabix table on
`annotation-store.md`, a flat file on `flat-file.md`, a transcript-model plugin on
`transcript-model.md`): the parameters, the fields and the data file's preparation. When the harness
prepares the file, add its preparation to
[docs/plugins/data-files.md](../../../docs/plugins/data-files.md). A difference from the module that
vep-rs keeps on purpose is an intended divergence, documented through
[fix-a-divergence](../fix-a-divergence/SKILL.md), never left unexplained.

Before classing a difference between the port's rows and the module's, read
[docs/intended-divergences.md](../../../docs/intended-divergences.md) and search its sections and
the pages they link for the plugin, the field or the shape. The Plugins section has a subsection for
shapes several plugins share and one per plugin with entries of its own. An entry there is intended
and the port keeps it, cited as the index's first column names it (the text of its row's first cell,
inner backticks dropped) or by the path of its page. Removing an entry needs the newer Ensembl VEP
behaviour stated in plain words, with the record that shows it.

Done when `check_parity_matrix.py` and `python3 scripts/docs/site_metrics.py --gate` pass.

## 8. Gate

Commit with `-s` and the model declaration, then run the full set in
[AGENTS.md](../../../AGENTS.md#before-pushing) (the admission check included, on the saved
description). A new `scripts/` test file or a test added to one moves its floor in `ci.yml`: name the
collected count under Verification for a maintainer to set. The `CHANGELOG.md` line under Added names
the plugin and its parameters in CONTRIBUTING.md's grammar. Open the pull request with the template's
four sections and nothing else, Verification carrying the citations, the golden result, the
concordance case or the reason it was not run, and the floors. Done when every command passes and its
output is quoted.

The fourth section, Intended divergences, cannot claim its first statement for a port. The support
table of `docs/plugins.md` and the mechanism pages under `docs/plugins/` carry entries of the index,
and the admission check refuses that statement to a change that edits a page the index links. Check
the second and name, one line each, the entries whose index rows link a page the port edits: the
entry as the index names it, that the port keeps it, and the row or section the edit adds beside it.
`grep -n '](plugins.md' docs/intended-divergences.md` lists the rows for the support table's page,
and the same grep with `plugins/<page>.md` those for a mechanism page. A difference of step 7 that
the port keeps on purpose checks the third as well, with the input, both engines' rows and why
vep-rs's row is the one to keep.
