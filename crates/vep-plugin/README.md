# vep-plugin

`vep-plugin` defines the C application binary interface (ABI) that a dynamically loaded
[Variant Effect Predictor (VEP)](../../README.md) plugin exports, opens such a library at run
time, and keeps the loaded plugins in a registry. `vep-cli`, the only workspace crate
that depends on it, resolves each `--plugin` name against [`vep-builtins`](../vep-builtins/README.md)
first and hands the names no built-in claims to this crate. Each loaded plugin returns one
JavaScript Object Notation (JSON) object per consequence row of the feature types it declares
(transcript rows alone without `vep_plugin_feature_types`), merged by `vep-cli` into that row's
`plugin_data` fields; a plugin that declares `Intergenic` returns one as well for the
intergenic row of a variant that overlaps no transcript, merged into the variant's
`plugin_data`.

The plugin support table is in [docs/plugins.md](../../docs/plugins.md), the data-file
preparation steps are in [docs/plugins/data-files.md](../../docs/plugins/data-files.md) and
the guide to writing a plugin is
[docs/developers/writing-a-plugin.md](../../docs/developers/writing-a-plugin.md); this page
describes the crate.

## How it works

Follow `--plugin target/debug/examples/libecho_plugin.so,hello`, the crate's own example
plugin passed by path with one parameter. `PluginRegistry::load_from_args` splits the argument
on every comma, as Perl's `split` in Ensembl VEP reads it (a trailing empty field is no parameter), into the
name or path and the parameter list `["hello"]`. The private
`resolve_plugin_path` keeps the first token because it is written as a file (`names_a_file`: a
path separator or a suffix) and a file exists at it; a file name without a separator is handed
over as `./<name>`, the form the loader opens as a file. A bare name is never read against the
working directory and is instead searched under `--dir_plugins` as `lib<name>.so` and then
`lib<name>.dylib`, each a regular file, and a name that resolves to neither is a
`PluginError::Load` reading `plugin '<name>' not found`. The resolved
path and the parameter list go to `LoadedPlugin::load`, which opens the library through
`libloading`, calls `vep_plugin_abi_version` and refuses any value other than
`VEP_PLUGIN_ABI_VERSION` with `PluginError::AbiMismatch`, resolves the other eight required
symbols (a missing one is a `PluginError::Load` naming the symbol) and the optional
`vep_plugin_feature_types`, calls `vep_plugin_create`, and reads the instance's name and
version strings. It then converts `hello` to a C string and passes it to `vep_plugin_init`,
which it calls on every load, with an empty list when the argument carried no parameter; a
non-zero return is a `PluginError::Init`. After `init` it reads `vep_plugin_header_info`, a
JSON array of `[key, description]` pairs it sorts by key, or `null` for none, and, when the
library exports it, `vep_plugin_feature_types`, a JSON array of row kinds it parses into
`FeatureType`s; text of another shape or an unknown kind is a `PluginError::Load`. The registry
logs the plugin's name
and version and keeps the `LoadedPlugin`. `vep-cli` reads `PluginRegistry::all_header_info`
when it opens the output, so the plugin's fields get their header lines and their `--tab`, VCF
CSQ and Parquet columns beside the built-in plugins'. During annotation, after the built-in
plugins have run on a batch, `vep-cli` asks `LoadedPlugin::runs_on` which rows each plugin
takes, serialises the variant once and each consequence some plugin takes once, converts each
document to a C string once, and calls `LoadedPlugin::run` with the two C strings for every
plugin that takes the row, in `--plugin` order; on the intergenic row of a variant that overlaps
no transcript the consequence document is `null`. `run` hands the strings to `vep_plugin_run`,
copies the returned text into a Rust `String`, releases the plugin's buffer through
`vep_plugin_free`, and returns the copy. `vep-cli` parses the copy as a JSON object and inserts
each entry into the consequence's `plugin_data`, or the variant's on an intergenic row, which
is where the example's `ECHO_FEATURE` and `ECHO_PARAM` fields reach the output's Extra column
and plugin columns. When the registry is dropped, each `LoadedPlugin` calls
`vep_plugin_destroy` on its instance before the library handle is released.

## Modules

| Module            | Role                                                                                                                                                                                       |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `src/lib.rs`      | Crate root: `PluginError` and the `FeatureType` row kind (`vep-core`'s `RowKind`, re-exported); re-exports `LoadedPlugin` and `PluginRegistry`.                                            |
| `src/abi.rs`      | The contract a plugin library exports: `VEP_PLUGIN_ABI_VERSION`, the ten `vep_plugin_*` symbol names in `symbol_names`, and their function-pointer types in `ffi`.                        |
| `src/loader.rs`   | `LoadedPlugin`: opens one library, checks the ABI version, resolves the symbols, creates and initialises the instance, reads its header info and feature types, and wraps `run` and destruction in safe methods. |
| `src/registry.rs` | `PluginRegistry`: parses `--plugin` arguments, resolves a name or path to a file, loads each plugin, lists their header info, answers which row kinds they take, and runs them all on one consequence and variant pair. |

Beside `src/`, [`examples/echo_plugin.rs`](examples/echo_plugin.rs) is a complete plugin that
the manifest builds as a `cdylib`, declaring two fields and the `Transcript` feature type, and
[`examples/echo_all_rows_plugin.rs`](examples/echo_all_rows_plugin.rs) is one that lists all
four feature types and so runs on every row, and [`examples/faulty_plugin.rs`](examples/faulty_plugin.rs)
one that fails on purpose (a non-object result, then a NULL) before it annotates, so the host's
per-plugin warnings can be tested. [`tests/echo_roundtrip.rs`](tests/echo_roundtrip.rs)
loads the two echo plugins through the loader and the registry, by path and by bare name under
a plugin directory, where only a regular `lib<name>` file answers to a bare name (a directory
or a plain entry of that name is passed over), compiles the first with `--cfg fixture_abi_one`
to check that a library of interface number 1 is refused and with `--cfg fixture_missing_free` to check that
a library lacking a required symbol is refused with the dynamic loader's reason, and the second
with `--cfg fixture_no_feature_types` to check that a library without the export runs on
transcript rows alone. The manifest lists the three `cfg` names under `[lints.rust]`; none is a
feature of the crate. [`tests/interface_contract.rs`](tests/interface_contract.rs) holds the
ABI version, the ten symbol names and their C signatures as literals of its own.

## Public API

### The C ABI

A plugin library exports ten symbols, nine required and `vep_plugin_feature_types` optional.
`VEP_PLUGIN_ABI_VERSION` is `2`; the loader refuses any other value before it touches the
remaining symbols. The signatures and the field lists of the two JSON documents are in
docs/developers/writing-a-plugin.md under
[The C ABI](../../docs/developers/writing-a-plugin.md#the-c-abi-version-2) and
[What `run` receives](../../docs/developers/writing-a-plugin.md#what-run-receives); this table states when the
host calls each symbol.

| Symbol                                   | When the host calls it                                                                                                  |
| ---------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `vep_plugin_abi_version`                 | Once at load, before any other symbol; a value other than `2` fails the load with `PluginError::AbiMismatch`.           |
| `vep_plugin_create`                      | Once at load; a null pointer fails the load with `PluginError::Init`.                                                   |
| `vep_plugin_name`, `vep_plugin_version`  | Once at load; a null pointer or text that is not valid UTF-8 fails the load; the pointers stay valid for the instance's lifetime. |
| `vep_plugin_init`                        | Once after `create`, on every load, with the `--plugin` parameters after the name, zero or more; a non-zero return is `PluginError::Init`. |
| `vep_plugin_header_info`                 | Once after `init`; the JSON array of `[key, description]` pairs is sorted by key, and the JSON text `null` declares no field; a null pointer or text of another shape is `PluginError::Load`. |
| `vep_plugin_feature_types`               | Once after `init`, only when the library exports it; an unknown kind in the JSON array is `PluginError::Load`; absent, the plugin runs on transcript rows alone. |
| `vep_plugin_run`                         | Once per row the plugin takes; a null return or text that is not valid UTF-8 is `PluginError::Run`.                     |
| `vep_plugin_free`                        | After every `vep_plugin_run`, `vep_plugin_header_info` or `vep_plugin_feature_types` result has been copied into a Rust `String`. |
| `vep_plugin_destroy`                     | When the `LoadedPlugin` is dropped, before the library is unloaded.                                                     |

### Rust types

- `PluginRegistry` (`src/registry.rs`) holds the loaded plugins. `new` makes an empty
  registry; `load_from_args(plugin_args, plugin_dirs)` loads one plugin per `--plugin` string,
  resolving a bare name under the `plugin_dirs` directories in their order; `plugins`,
  `plugin_count` and `is_empty` expose the set; `all_header_info` lists every plugin's
  `(key, description)` pairs, plugin by plugin in load order; `runs_on(row)` says whether any
  plugin takes a row of that kind; `plugins_wanting_regulatory_features` names the loaded
  plugins that list `RegulatoryFeature` or `MotifFeature`, each of which switches regulatory
  annotation on for the run, as it does in Ensembl VEP; `field_row_kinds` pairs each plugin's
  fields with the row kinds the plugin takes, for the writers that print a field on those rows
  alone; `run_all(consequence_json, variant_json)` runs every plugin and returns
  `(plugin_name, result_json)` pairs. `registry::names_a_file(token)` says whether a `--plugin`
  token is written as a file (a path separator or a suffix) rather than as a bare plugin name;
  `vep-cli` applies it before its Perl-module check.
- `LoadedPlugin` (`src/loader.rs`) is one opened library and its instance. `load(path, params)`
  performs the ABI check, symbol resolution, `create`, `init` and the header and feature-type
  reads; `run(consequence_json, variant_json)` takes the two documents as `&CStr` and wraps the
  C call; `name`, `version`, `header_info` and `feature_types` return what was read at load,
  and `runs_on(row)` applies the feature-type list. The type is marked `Send + Sync`, so a
  plugin instance must tolerate calls from any thread.
- `PluginError` (`src/lib.rs`) has four variants: `Init`, `Run`, `Load` and
  `AbiMismatch { expected, got }`. It is `#[non_exhaustive]`, so a `match` on it outside the
  crate ends in a wildcard arm. `vep-cli` matches the `Load` text
  `plugin '<name>' not found` to tell a name that resolves to no file from a file that failed
  to load. The name is settled by itself: an Ensembl VEP plugin vep-rs does not provide, or a
  `<name>.pm` under the plugins directories, is refused with or without `--safe`, and any other
  name is warned about and skipped, Ensembl VEP's default, or refused under `--safe`. A file
  that failed to load is dropped with a warning and the run goes on without it,
  `Failed to compile plugin <name>` for `Load` and `AbiMismatch` and
  `Failed to instantiate plugin <name>` for `Init`; under `--safe` the warning is the run's
  error.
- `abi::symbol_names` and `abi::ffi` (`src/abi.rs`) are the null-terminated symbol names the
  loader looks up and the `unsafe extern "C"` function-pointer types it casts them to.
- `FeatureType` (`src/lib.rs`) is the kind of row a plugin takes, `Transcript`,
  `RegulatoryFeature`, `MotifFeature` or `Intergenic`: `vep-core`'s `RowKind` under the name
  Ensembl VEP gives the list, with `parse` and `as_str` for the names a plugin lists and a
  conversion from the consequence feature type, the one type the built-in plugins' `feature_types`
  return as well.

## Run it

```bash
cargo build -p vep-plugin
cargo test -p vep-plugin
cargo build --example echo_plugin -p vep-plugin
cargo doc --no-deps -p vep-plugin
```

`cargo test -p vep-plugin` builds the three example libraries into `target/debug/examples/`
before `tests/echo_roundtrip.rs` loads the two echo plugins; a
`cargo test -p vep-plugin --test echo_roundtrip` filter alone does not build the examples, and
the test fails with a message saying so. The example build writes `libecho_plugin.so` on Linux
and `libecho_plugin.dylib` on macOS, the file the worked example in
docs/developers/writing-a-plugin.md passes to `vep --plugin`. The refusal test also runs
`rustc` on `examples/echo_plugin.rs` with `--cfg fixture_abi_one`, so the toolchain's `rustc`
must be on the path. The repository's standard commands (`cargo test --workspace`,
`cargo clippy`, `cargo fmt`) are in the [contributing guide](../../CONTRIBUTING.md).

## Related

- [`vep-core`](../vep-core/README.md): the `TranscriptConsequence` and `InputVariant` types
  serialised into the two JSON documents a plugin receives.
- [`vep-cli`](../vep-cli/README.md): the only crate that depends on `vep-plugin`. Its
  [`dylib_plugins` module](../vep-cli/src/dylib_plugins.rs) loads the `--plugin` names no
  built-in claimed, runs each plugin on the rows it takes and merges each `run` result into
  `plugin_data`.
- [`vep-builtins`](../vep-builtins/README.md): the compiled-in plugins the runner resolves
  before falling back to this crate.
- [docs/plugins.md](../../docs/plugins.md): the plugin support table.
- [docs/developers/writing-a-plugin.md](../../docs/developers/writing-a-plugin.md): the ABI
  signatures, the JSON field lists and the worked example.
- [docs/cli/plugins.md](../../docs/cli/plugins.md): the `--plugin` and `--dir_plugins` flags.
- [Root README](../../README.md).
