# Plugins

These flags load the built-in plugins and dynamic-library plugins, and name the directory the latter are searched in. The built-in plugins, their parameters and their data files are described under [Plugins](../plugins.md); the C ABI a dynamic-library plugin implements is described under [Writing a plugin](../developers/writing-a-plugin.md). A `--plugin` naming an Ensembl VEP plugin vep-rs does not provide ends the run with an error naming the reason, with or without `--safe`, as does a Perl module (`Name.pm` by path, or `Name` with `Name.pm` under the plugins directory). The [support table](../plugins.md) lists every Ensembl VEP plugin and its status.

The flags under [Custom annotation](#custom-annotation) name the `--custom` sources, their bigWig and conservation-score shortcuts, and how a VCF source's values are written. The sources themselves are described under [Custom annotation](../custom-annotation.md).

## Loading plugins

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--plugin` | String (repeatable) | none | Load plugin(s) - can be specified multiple times |
| `--dir_plugins` | String | none | Directory the plugin libraries are looked up in; default \<`--dir`\>/Plugins. The value is expanded as Ensembl VEP's `glob` expands it: `*`, `?`, `[...]` and `{a,b}` match as in a shell, `~` and `~user` name a home directory, and a backslash quotes the character after it. Every match that is a directory is searched, the last match first; a match that is not a directory is warned about when `--plugin` names a plugin (an error under `--safe`). Spaces are part of the value, not separators |
| `--safe` | bool | false | Make a plugin that cannot be loaded, or a missing plugins directory, an error instead of a warning, and skip the check of incompatible flag pairs |
| `--debug` | bool | false | Accepted without effect: Ensembl VEP reads it only to print the plugin directories and modules it loads; the annotation output is unchanged |

## Custom annotation

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--custom` | String (repeatable) | none | Custom annotation source, as Ensembl VEP reads it: file=\<f\>,format=\<vcf\|bed\|bigwig\|gff\|gtf\>, short\_name=\<n\>,type=\<overlap\|within\|surrounding\|exact\>,coords=\<0\|1\>,fields=\<A%B\>, overlap\_cutoff=\<pct\>,distance=\<bp\>,same\_type=\<0\|1\>,reciprocal=\<0\|1\>,num\_records=\<N\|all\>, summary\_stats=\<min%max%mean%count%sum\>,gff\_type=\<transcript\|gencode\_promoter\>, or the positional \<file\>,\<short\_name\>,\<format\>,\<type\>,\<coords\>,\<field\>,...; a file holding ###CHR### stands for one source per chromosome whose file exists. The file is a local path (a remote URL is refused). A VCF, a BED or a GFF in gencode\_promoter mode is bgzipped and tabix-indexed; a bigWig needs neither. A GTF, or a GFF in its default transcript mode, is a transcript model read whole, plain, or gzipped or bgzipped under a name ending .gz, with no index, and annotated beside the cache (see `--gff`). May be given more than once |
| `--bigwig` | String | none | Shortcut for `--custom` \<FILE\>,,bigwig,exact: a bigWig whose value over exactly the variant's span prints under the file's name |
| `--phyloP` | String (repeatable) | none | Shortcut for the UCSC phyloP conservation scores of an N-way alignment, read from the local mirror `--ucsc_data_root` names as \<root\>/\<assembly\>/phyloP\<N\>way/\<assembly\>.phyloP\<N\>way.bw with `--custom` ...,phyloP\<N\>way,bigwig,exact, so the column is phyloP\<N\>way. Needs `--ucsc_assembly`. May be given more than once |
| `--phastCons` | String (repeatable) | none | Shortcut for the UCSC phastCons conservation scores of an N-way alignment, read from the local mirror `--ucsc_data_root` names as \<root\>/\<assembly\>/phastCons\<N\>way/\<assembly\>.phastCons\<N\>way.bw with `--custom` ...,phastCons\<N\>way,bigwig,exact. Needs `--ucsc_assembly`. May be given more than once |
| `--ucsc_assembly` | String | none | The UCSC assembly name `--phyloP` and `--phastCons` build their paths from: hg38 for GRCh38, hg19 for GRCh37 |
| `--ucsc_data_root` | String | `http://hgdownload.cse.ucsc.edu/goldenpath/` | The directory `--phyloP` and `--phastCons` read their bigWigs from, laid out as UCSC's goldenPath tree (\<root\>/\<assembly\>/phyloP\<N\>way/\<assembly\>.phyloP\<N\>way.bw for `--phyloP`, \<root\>/\<assembly\>/phastCons\<N\>way/\<assembly\>.phastCons\<N\>way.bw for `--phastCons`). The default is Ensembl VEP's remote root, which is refused: a local mirror must be named |
| `--custom_multi_allelic` | bool | false | Keep a comma-separated INFO value of a `--custom` VCF source whole on every allele of a multi-allelic record; without it a value with one entry per alternate allele is split, each allele taking its own entry |
| `--custom_suppress_filter` | bool | false | Leave the FILTER column of a `--custom` VCF source out of its records unless the source's fields name FILTER; without it every source with fields adds \<short\_name\>\_FILTER |
| `--gff` | String | none | Shortcut for `--custom` file=\<FILE\>,format=gff: a GFF3 transcript model whose transcripts are annotated beside the cache's, each of their rows carrying SOURCE=\<the file's name\>. Needs `--fasta`, which supplies the transcripts' sequence |
| `--gtf` | String | none | Shortcut for `--custom` file=\<FILE\>,format=gtf: a GTF transcript model whose transcripts are annotated beside the cache's, each of their rows carrying SOURCE=\<the file's name\>. Needs `--fasta`, which supplies the transcripts' sequence |

