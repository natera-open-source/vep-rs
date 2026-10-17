# Plugins

These flags load the built-in plugins and dynamic-library plugins, and name the directory the latter are searched in. The built-in plugins, their parameters and their data files are described under [Plugins](../plugins.md); the C ABI a dynamic-library plugin implements is described under [Writing a plugin](../developers/writing-a-plugin.md).

## Summary

| Flag | Default | What it does |
| --- | --- | --- |
| [`--plugin`](#opt_plugin) | none | Loads a plugin by built-in name, dylib path or dylib name; repeatable |
| [`--custom`](#opt_custom) | none | No effect; custom annotation sources are not loaded |
| [`--dir_plugins`](#opt_dir_plugins) | none | Directory searched for a dylib named by `--plugin` |

## Accepted for compatibility

`--custom_suppress_filter` is refused with an error and is described under [Refused flags](compatibility.md#refused-flags).

| Flag | What Ensembl VEP does with it |
| --- | --- |
| `--custom` | Adds annotation from a tabix-indexed or bigWig file |

## Options

### `--plugin` {#opt_plugin}

**`--plugin <Name,param1,param2,...>`**, string, repeatable, no default. Ensembl VEP: same flag and default.

Loads a plugin. `Name` is a built-in (the names under [Plugins](../plugins.md), matched case-sensitively), the path of a C-ABI dylib, or a dylib name resolved as `lib<Name>.so` or `lib<Name>.dylib` under `--dir_plugins`. A name that is none of these is warned about and skipped, which is Ensembl VEP's own default. A dylib that exists but fails to open, reports another ABI version or fails `init` aborts the run.

### `--custom` {#opt_custom}

**`--custom <spec>`**, string, repeatable, no default. Ensembl VEP: accepted, no effect.

Not implemented: custom annotation sources are not loaded.

### `--dir_plugins` {#opt_dir_plugins}

**`--dir_plugins <dir>`**, string, no default. Ensembl VEP: same flag, vep-rs default differs: no directory is searched without the flag.

Directory searched for `lib<Name>.so` or `lib<Name>.dylib` when `--plugin <Name>` is neither a built-in nor a path.
