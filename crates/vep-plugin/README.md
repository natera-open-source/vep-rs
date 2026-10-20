# vep-plugin

`vep-plugin` defines the C application binary interface (ABI) that a dynamically loaded
[Variant Effect Predictor (VEP)](../../README.md) plugin exports, opens such a library at run
time, and keeps the loaded plugins in a registry. `vep-cli`, the only workspace crate
that depends on it, resolves each `--plugin` name against [`vep-builtins`](../vep-builtins/README.md)
first and hands the names no built-in claims to this crate. Each loaded plugin returns one
JavaScript Object Notation (JSON) object per transcript consequence, merged by `vep-cli` into
that consequence's `plugin_data` fields.

The plugin support matrix, the data-file preparation steps and the guide to writing a plugin
live in [docs/plugins.md](../../docs/plugins.md); this page describes the crate.

## How it works

Follow `--plugin target/debug/examples/libecho_plugin.so,hello`, the crate's own example
plugin passed by path with one parameter. `PluginRegistry::load_from_args` splits the argument
at its first comma into the name or path and the parameter list `["hello"]`. The private
`resolve_plugin_path` keeps the first token because a file exists at it; a bare name is instead
searched under `--dir_plugins` as `lib<name>.so` and then `lib<name>.dylib`, and a name that
resolves to neither is a `PluginError::Load` reading `plugin '<name>' not found`. The resolved
path goes to `LoadedPlugin::load`, which opens the library through `libloading`, calls
`vep_plugin_abi_version` and refuses any value other than `VEP_PLUGIN_ABI_VERSION` with
`PluginError::AbiMismatch`, resolves the other seven symbols (a missing one is a
`PluginError::Load` naming the symbol), calls `vep_plugin_create`, and reads the instance's
name and version strings. Because the parameter list is not empty, the registry calls
`LoadedPlugin::init`, which converts `hello` to a C string and passes it to `vep_plugin_init`;
a non-zero return is a `PluginError::Init`. The registry logs the plugin's name and version and
keeps the `LoadedPlugin`. During annotation, after the built-in plugins have run on a batch,
`vep-cli` serialises the variant once and each of its transcript consequences, and calls
`LoadedPlugin::run` with the two JSON strings for every consequence and every loaded plugin
in `--plugin` order. `run` hands the strings to `vep_plugin_run` as C strings, copies the
returned text into a Rust `String`, releases the plugin's buffer through `vep_plugin_free`, and
returns the copy. `vep-cli` parses the copy as a JSON object and inserts each entry into the
consequence's `plugin_data`, which is where the example's `ECHO_FEATURE` and `ECHO_PARAM`
fields reach the output's Extra column. When the registry is dropped, each `LoadedPlugin`
calls `vep_plugin_destroy` on its instance before the library handle is released.

## Modules

| Module            | Role                                                                                                                                                                             |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/lib.rs`      | Crate root: `PluginError`, the `FeatureType` enum and the `VepPlugin` trait; re-exports `LoadedPlugin` and `PluginRegistry`.                                                     |
| `src/abi.rs`      | The contract a plugin library exports: `VEP_PLUGIN_ABI_VERSION`, the eight `vep_plugin_*` symbol names in `symbol_names`, and their function-pointer types in `ffi`.             |
| `src/loader.rs`   | `LoadedPlugin`: opens one library, checks the ABI version, resolves the symbols, creates the instance, and wraps `init`, `run` and destruction in safe methods.                   |
| `src/registry.rs` | `PluginRegistry`: parses `--plugin` arguments, resolves a name or path to a file, loads and initialises each plugin, and runs them all on one consequence and variant pair.       |

Beside `src/`, [`examples/echo_plugin.rs`](examples/echo_plugin.rs) is a complete plugin that
the manifest builds as a `cdylib`, and [`tests/echo_roundtrip.rs`](tests/echo_roundtrip.rs)
loads it through the loader and the registry, by path and by bare name under a plugin
directory.

## Public API

### The C ABI

A plugin library exports eight symbols. `VEP_PLUGIN_ABI_VERSION` is `1`; the loader refuses
any other value before it touches the remaining seven. The signatures and the field lists of
the two JSON documents are in docs/developers/writing-a-plugin.md under
[The C ABI](../../docs/developers/writing-a-plugin.md#the-c-abi-version-1) and
[What `run` receives](../../docs/developers/writing-a-plugin.md#what-run-receives); this table states when the
host calls each symbol.

| Symbol                                   | When the host calls it                                                                                                  |
| ---------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `vep_plugin_abi_version`                 | Once at load, before any other symbol; a value other than `1` fails the load with `PluginError::AbiMismatch`.           |
| `vep_plugin_create`                      | Once at load; a null pointer fails the load with `PluginError::Init`.                                                   |
| `vep_plugin_name`, `vep_plugin_version`  | Once at load; a null pointer or text that is not valid UTF-8 fails the load; the pointers stay valid for the instance's lifetime. |
| `vep_plugin_init`                        | Once after load, only when the `--plugin` argument carried at least one parameter; a non-zero return is `PluginError::Init`. |
| `vep_plugin_run`                         | Once per transcript consequence per loaded plugin; a null return or text that is not valid UTF-8 is `PluginError::Run`. |
| `vep_plugin_free`                        | After every `vep_plugin_run` result has been copied into a Rust `String`.                                               |
| `vep_plugin_destroy`                     | When the `LoadedPlugin` is dropped, before the library is unloaded.                                                     |

### Rust types

- `PluginRegistry` (`src/registry.rs`) holds the loaded plugins. `new` makes an empty
  registry; `load_from_args(plugin_args, plugin_dir)` loads one plugin per `--plugin` string,
  resolving a bare name under `plugin_dir`; `plugins`, `plugin_count` and `is_empty` expose
  the set; `run_all(consequence_json, variant_json)` runs every plugin and returns
  `(plugin_name, result_json)` pairs.
- `LoadedPlugin` (`src/loader.rs`) is one opened library and its instance. `load(path)`
  performs the ABI check, symbol resolution and `create`; `init(params)` and
  `run(consequence_json, variant_json)` wrap the C calls; `name` and `version` return the
  strings read at load. The type is marked `Send + Sync`, so a plugin instance must tolerate
  calls from any thread.
- `PluginError` (`src/lib.rs`) has four variants: `Init`, `Run`, `Load` and
  `AbiMismatch { expected, got }`. `vep-cli` matches the `Load` text
  `plugin '<name>' not found` to warn and skip a name that resolves to no file, Ensembl VEP's
  default without `--safe`; any other load error aborts the run before annotation starts.
- `abi::symbol_names` and `abi::ffi` (`src/abi.rs`) are the null-terminated symbol names the
  loader looks up and the `unsafe extern "C"` function-pointer types it casts them to.
- `VepPlugin` and `FeatureType` (`src/lib.rs`) describe an in-process plugin interface with
  the same `run(consequence, variant)` shape over native `vep-core` types. The workspace's
  built-in plugins implement `BuiltinPlugin` from `vep-builtins` instead, and a loaded library
  exports the C ABI, so no type in this repository implements `VepPlugin`.

## Run it

```bash
cargo build -p vep-plugin
cargo test -p vep-plugin
cargo build --example echo_plugin -p vep-plugin
cargo doc --no-deps -p vep-plugin
```

`cargo test -p vep-plugin` builds the example library into `target/debug/examples/` before
`tests/echo_roundtrip.rs` loads it; a `cargo test -p vep-plugin --test echo_roundtrip` filter
alone does not build the example, and the test fails with a message saying so. The example
build writes `libecho_plugin.so` on Linux and `libecho_plugin.dylib` on macOS, the file the
worked example in docs/plugins.md passes to `vep --plugin`. The repository's standard
commands (`cargo test --workspace`, `cargo clippy`, `cargo fmt`) are in the
[contributing guide](../../CONTRIBUTING.md).

## Related

- [`vep-core`](../vep-core/README.md): the `TranscriptConsequence` and `InputVariant` types
  named by the `VepPlugin` trait and serialised into the two JSON documents a plugin
  receives.
- [`vep-cli`](../vep-cli/README.md): the only crate that depends on `vep-plugin`. Its
  [`dylib_plugins` module](../vep-cli/src/dylib_plugins.rs) loads the `--plugin` names no
  built-in claimed and merges each `run` result into `plugin_data`.
- [`vep-builtins`](../vep-builtins/README.md): the compiled-in plugins the runner resolves
  before falling back to this crate.
- [docs/plugins.md](../../docs/plugins.md): the plugin matrix, the ABI table, the JSON field
  lists and the worked example.
- [docs/cli-reference.md](../../docs/cli-reference.md): the `--plugin` and `--dir_plugins`
  flags.
- [Root README](../../README.md).
