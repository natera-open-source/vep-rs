# Versioning and compatibility

What a vep-rs version number covers, what may change in any release, and what a release asks of caches, plugin stores and plugins an earlier one built.

**In short.** One version number covers the nine crates and the programs built from them. In a release whose major number is 1 or more it makes a promise for five things: the documented flags and what they do, the output formats, the exit codes, the plugin interface and the documented part of the libraries. Four things may change in any release: hidden tuning flags, abbreviations, message wording and library internals. A release that corrects an annotation row is a fix, never a breaking change. A release inside a major version may need caches and plugin stores built again, and its release notes say when.

## What the version number covers

vep-rs has one version number. The nine crates of the workspace and the programs built from them are released together under it: `vep` with its `filter` and `recode` commands, `vep-cache-builder` and `vep-cache-converter`. `vep --version` prints the number. The number follows [Semantic Versioning](https://semver.org/), as the [Changelog](../changelog.md) declares.

The version number covers five things: the documented flags and what they do, the output formats, the exit codes, the plugin interface and the documented part of the libraries. In a release whose major number is 1 or more, a change that breaks one of the five raises the major number. Every rule on this page binds such a release. A release numbered 0.y.z stands outside them: it raises its minor number for a change that breaks one of the five.

| Covered | Where it is documented |
| --- | --- |
| The documented flags and what they do | The [CLI reference](../cli-reference.md): a row for each flag of `vep`, `vep filter` and `vep recode`, with its names, its value type, its default and its description. [Cache setup](../cache-setup.md) and the pages beside it document the flags of the two cache tools |
| The output formats | [Output formats](../output-formats.md) and the page of each format: [default and tab](../output/vep-and-tab.md), [VCF](../output/vcf.md), [JSON](../output/json.md) and [Parquet](../output/parquet.md) |
| The exit codes | [Exit status and errors](../cli/compatibility.md#exit-status-and-errors): `0`, `1` and `2` |
| The plugin interface | [Writing a plugin](../developers/writing-a-plugin.md); [The plugin interface](#the-plugin-interface) below says which parts |
| The documented part of the libraries | Each crate's API documentation and [Embedding vep-rs as a library](../developers/embedding.md); [The libraries](#the-libraries) below says which part |

## What may change in any release

Four things may change in any release: hidden tuning flags, abbreviations, message wording and library internals. A script or a program that must keep working across releases builds on none of them.

- **Hidden tuning flags.** The CLI reference marks such a flag `(hidden)`, and `--help` leaves it out. Its name, its values, its default and its effect may change, and the flag may be removed, however it is given ([Cache, species and runtime](../cli/cache-species-runtime.md#runtime)).
- **Abbreviations.** `vep`, `vep filter` and `vep recode` read a flag from a prefix of its name. A prefix is resolved against the flags of the release that reads it, so a release that adds a flag can give the prefix another meaning or have it refused. The pages that offer abbreviations say so: [Option spellings](../cli/compatibility.md#option-spellings-configpm62), [vep filter](../cli/filter.md) and [vep recode](../cli/recode.md). A script writes each flag's full name.
- **Message wording.** The words of an error, a warning or a notice that a program prints. A script reads the exit code, not the words.
- **Library internals.** An item that a crate's API documentation does not show ([The libraries](#the-libraries)).

## Annotation output

The output formats are covered; the annotation a row carries is measured against the reference Ensembl VEP release and the [Intended divergences](../intended-divergences.md) index, never against an earlier vep-rs release.

A release that corrects an annotation row is a fix and never a breaking change. A user reads of such a row in two places. The [Changelog](../changelog.md) lists the correction in the release's section. The [Intended divergences](../intended-divergences.md) index gains, changes or loses an entry where the row differs from Ensembl VEP's on purpose.

A move of the reference Ensembl VEP release is a minor release, and the rows that change with it are not a breaking change. On the `##VEP=` line of a VCF output, the reference release, its API version and its five module versions move with the reference, and the `##vep-rs=` version moves with each release; the keys of `##VEP=` and their order belong to the output format ([VCF output](../output/vcf.md#header)).

## Converted caches and plugin stores

A release inside a major version may need caches and plugin stores built again, and its release notes say when. This holds for a cache converted from an Ensembl VEP cache, for a cache `vep-cache-builder` built and for a plugin store `vep-cache-converter convert-plugin-data` wrote.

The list below says what the programs of this release do when they meet a file an earlier release wrote. None of it is a promise for a later release.

- **The transcript cache** carries no mark of its layout. Its `info.json` names the Ensembl release under `cache_version`, and a key the program does not know is passed over.
- **A regulatory tree** carries the mark `regulatory_schema` in `info.json`. A run that reads the tree refuses one whose mark is missing or is not the one the program reads, and the message names the conversion ([Regulatory and motif features](../regulatory.md#the-cache)).
- **A cache with no regulatory tree**, which is what the converters of releases 0.1.0 to 0.5.0 wrote, runs. Under `--regulatory` it draws one warning and no regulatory or motif row.
- **A plugin store** (`.vpd`) carries an index version. A plugin refuses a store whose index version is not the one it reads and names the conversion that rebuilds it ([When a plugin refuses a store](../plugins/data-files.md#when-a-plugin-refuses-a-store)).
- **A cache whose `info.json` has no `left_out` key** loads and runs. A build or a conversion lists under that key what it left out, and a run draws from the key the lines that say so; without the key it draws none of them ([What a native cache lacks](../cache-variation-data.md#what-a-native-cache-lacks)).

## The plugin interface

A plugin is a shared library that exports the C interface [Writing a plugin](../developers/writing-a-plugin.md#the-c-abi-version-2) documents. The version number covers three parts of it:

- the interface number, which the library's `vep_plugin_abi_version` returns and the host compares with its own
- the exported symbols with their signatures
- the two JSON documents `vep_plugin_run` receives, the consequence document and the variant document ([What `run` receives](../developers/writing-a-plugin.md#what-run-receives))

A plugin library built against the interface number a release documents loads through that whole major version. A release may add a key to a document; inside a major version it renames or removes none.

## The libraries

The nine crates carry one version: `vep-core`, `vep-fasta`, `vep-io`, `vep-effects`, `vep-plugin`, `vep-builtins`, `vep-cli`, `vep-cache-builder` and `vep-cache-converter`. The last is a program with no library.

The documented part of a library is what the crate's API documentation shows: the pages `cargo doc` builds from the crate. [Embedding vep-rs as a library](../developers/embedding.md) walks through the documented part of `vep-cli`.

Everything else belongs to the library internals. An item hidden from the API documentation belongs to them, though it is `pub` for the workspace's own programs and tests. A binary and a cargo feature that no page documents belong to them as well: the development binary `bench_plugins` of `vep-cli`, behind its default-off `bench` feature, is such a binary.

Some documented signatures name a type of another project: `anyhow::Result`, `indexmap::IndexMap`, `rayon::ThreadPool`, `serde_json::Value`, `smallvec::SmallVec` and the `serde` traits `Serialize` and `Deserialize`.

## The minimum supported Rust version

The minimum supported Rust version is 1.89. It is declared once, as `rust-version` under `[workspace.package]` of the workspace's [`Cargo.toml`](../../Cargo.toml), and the nine crates inherit it. The `MSRV` job of the CI workflow (`.github/workflows/ci.yml`) builds the workspace on that toolchain and fails where its toolchain and the declared figure differ.

The minimum Rust version may rise in a minor release, and that release's notes name the new figure. The [Changelog](../changelog.md) records the rise.

## Which release receives fixes

Fixes go to the newest release only. A fix is released as a new version above the newest one, and no earlier release receives it.

## How a change is announced

Every release has a section of the [Changelog](../changelog.md), written as [Changelog entries](../contributing.md#changelog-entries) says: one entry for each user-visible change, under its kind. The entry of a breaking change says so in its own sentence, as `a breaking change of the libraries` or `a breaking change of the command line`, and then says what a user does. The release notes, the text of the release's page, copy the section or restate its entries.

A documented flag or library item that is retired is marked for at least one minor release, and it is removed only in a release that raises the major number.
