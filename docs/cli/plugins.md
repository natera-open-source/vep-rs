# Plugins

These flags load the built-in plugins and dynamic-library plugins, and name the directory the latter are searched in. The built-in plugins, their parameters and their data files are described under [Plugins](../plugins.md); the C ABI a dynamic-library plugin implements is described under [Writing a plugin](../developers/writing-a-plugin.md).

## Flags

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--plugin` | String (repeatable) | none | Load plugin(s) - can be specified multiple times |
| `--custom` | String (repeatable) | none | Custom annotation source(s) (accepted for VEP compatibility; not implemented, sources are not loaded) |
| `--custom_suppress_filter` | bool | false | Refused with an error: `--custom_suppress_filter` is not supported: it would keep a `--custom` VCF source's FILTER column out of the output unless the `--custom` source's fields argument names it, and custom annotation is not supported |
| `--dir_plugins` | String | none | Directory the plugin libraries are looked up in; default \<`--dir`\>/Plugins. A directory that does not exist is warned about when `--plugin` names a plugin (an error under `--safe`). A value holding a glob pattern (`*`, `?`, `[...]`, `{a,b}`) or a `~user` prefix is read as one literal path; Ensembl VEP expands it with `glob` |
| `--safe` | bool | false | Make a plugin that cannot be loaded, or a missing plugins directory, an error instead of a warning, and skip the check of incompatible flag pairs |

