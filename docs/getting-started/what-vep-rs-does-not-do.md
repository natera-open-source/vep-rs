# What vep-rs does not do

What to expect when a command line brought from Ensembl VEP asks vep-rs for a feature it does not have: which features those are and what the run does.

Every flag Ensembl VEP accepts is honoured, refused with an error that names it, or accepted without effect with its help text saying so. The [CLI reference](../cli-reference.md) and its [Compatibility with Ensembl VEP](../cli/compatibility.md) page carry the per-flag detail. Each family below gets its flags and what a run does with them.

## RefSeq and merged caches

The RefSeq and merged caches that `--refseq` and `--merged` annotate are Ensembl's, converted ([Cache setup](../cache-setup.md#refseq-and-merged-caches)): vep-rs builds no RefSeq model natively and produces no merged cache. `--bam`, which corrects RefSeq transcript models against an alignment file as a cache, GFF, GTF or database source is read, is refused with an error naming the route (exit 1); a converted RefSeq or merged cache carries the edits Ensembl applied when it was dumped.

## Custom annotation sources

`--custom` reads local VCF, BED, GFF promoter and bigWig sources ([Custom annotation](../custom-annotation.md)). Two things Ensembl VEP reads through it are not read. A remote file (`http://`, `ftp://` or any URL) is refused with an error naming a local copy, where Ensembl VEP reads it over the network. A GFF or GTF transcript model is read beside a cache in the same run (`--gff`, `--gtf`); a run with the file alone is refused with the two-step cache recipe, where Ensembl VEP reads it in place of its cache.

A BAM is not read. Ensembl VEP corrects such a model against one named by `--bam`. `--bam` is refused, and `format=bam` is refused as an unsupported format (exit 1), as Ensembl VEP refuses it.

A `--custom` run still reads a cache, where Ensembl VEP accepts the source as the run's only annotation source. `--gff` and `--gtf`, Ensembl VEP's shortcuts to a `--custom` transcript source, read the file beside the cache. Without a cache they are refused with an error that prints the cache recipe (exit 1): build a cache from the file with `vep-cache-builder` and annotate against it ([A GTF or GFF3 as the transcript model](../cache-from-gtf-or-gff3.md)). `--phyloP` and `--phastCons` read their bigWigs from a local mirror alone: the remote root Ensembl VEP reads by default is refused ([Conservation shortcuts](../output/custom-bigwig.md#conservation-shortcuts)).

`--regulatory_gff` reads local files alone as well. A file whose name opens `http://`, `https://` or `ftp://` is refused with an error that says so, where Ensembl VEP reads an `http://` or `ftp://` file over the network ([Regulatory and motif features](../regulatory.md#files-and-records-that-are-refused)).

## Plugins

Ensembl VEP's plugins are Perl modules, which vep-rs cannot load. A `--plugin` name resolves to one of the built-in plugins listed on [Plugins](../plugins.md), else to a shared library under `--dir_plugins`. The rest of Ensembl VEP's plugin catalogue has no vep-rs counterpart.

What a run does with any other name:

- With the name of a plugin from that catalogue and no shared library of that name, `--plugin` ends the run (exit 1), with or without `--safe`. The error begins `ERROR: [plugin] <Name> is an Ensembl VEP plugin vep-rs does not provide:` and states the reason.
- A `.pm` file, by path or under `--dir_plugins`, ends the run with an error that begins `ERROR: [plugin] <path> is an Ensembl VEP Perl plugin`.
- A name that is none of these is skipped with a warning, or refused under `--safe` (an error that begins `ERROR: [plugin] Failed to load plugin <Name>`).
- Under `--safe` a plugins directory that does not exist ends the run before any of these (`ERROR: [plugin] Plugins directory '<dir>' not found.`).

A built-in plugin's data file is a local path. A URL (`http://`, `ftp://` or any
other scheme) fails the plugin's instantiation with `remote data files are not
supported; download <url> and index it locally`, where Ensembl VEP's tabix
plugins stream an `http://` or `ftp://` file over the network
(`BaseVepTabixPlugin.pm` lines 271-274). LoFTEE's `key:value` word holding a
URL is dropped and logged instead, as Ensembl VEP's plugin drops it
([LoFTEE](../plugins/loftee.md)).

## Parquet columns

`--output_format parquet` writes a fixed column set, six key columns and then the default `--tab` columns, each under a declared type ([Parquet](../output/parquet.md)). `--fields` beside it is refused with an error: `ERROR: [output] --fields applies to --tab and --vcf output; the Parquet column set is fixed`. Parquet is vep-rs's own format, with no Ensembl VEP counterpart; a column set chosen with `--fields` is what it lacks.

## The companion tools

Ensembl VEP ships three companion tools:

- `filter_vep`, which filters an annotated output by field expressions;
- `haplo`, which computes transcript haplotypes from phased genotypes;
- `variant_recoder`, which translates a variant between its identifiers and notations.

The release archives and the container image carry no `haplo`.

`filter_vep` is `vep filter` ([its page](../cli/filter.md)). Its `--ontology` is answered from the built-in Sequence Ontology table instead of an Ensembl database. The database connection flags (`--host`, `--port`, `--user`, `--pass`, `--version`, `--registry`) are refused with an error naming the flag.

`variant_recoder` is `vep recode` ([its page](../cli/recode.md)): the same options over the cache alone. The synonyms and the co-located variants come from the cache's variation data, where `variant_recoder` reads the Ensembl database. So the merged dbSNP identifiers the database lists as synonyms (`Archive dbSNP`) are not written. `--lrg` and `--ambiguous_hgvs` are refused with the database connection flags.

A RefSeq notation is written from a merged cache, and a run on a cache read as the Ensembl set says so. A ClinGen allele identifier is not an input, and `ga4gh_vrs` is written only for a chromosome whose synonyms name an `NC_` accession. For any other chromosome the `ga4gh_vrs` key is still present under `--ga4gh_vrs`, as it is in `variant_recoder`, and holds an empty list.

## Debug output

`--debug` is accepted without effect. Ensembl VEP reads it only to print the plugin directories and modules it loads; that detail is not printed, and one `INFO:` line on standard error says so unless `--quiet` is given. The flag changes no annotation.

## Ensembl database mode

`--database`, `--host`, `--port`, `--user`, `--password`, `--registry`, `--db_version`, `--genomes` and `--is_multispecies` are each refused with an error naming the flag. `--lrg` and `--ambiguous_hgvs`, which need the database, are refused with them. The error for `--database` begins `ERROR: [cli] --database is not supported: vep-rs reads a local cache only and never connects to an Ensembl database`.

The refusal holds under `--offline` too. Ensembl VEP accepts `--offline --database`, since its `--offline` option set switches `--database` off before the incompatibility check runs (`Config.pm` line 766, then line 769), and runs offline. vep-rs refuses the pair.

## 1000 Genomes pilot coordinates

`--gp` is refused with an error naming it, which begins `ERROR: [cli] --gp is not supported: it reads a GRCh37 chr:pos from the INFO/GP field that some 1000 Genomes VCFs carry`. In Ensembl VEP the flag reads a record's chromosome and position from the `GP` field of its INFO column and skips a record without that field (`Config.pm` line 120, `Parser/VCF.pm` lines 271-285). That field is the convention of the 1000 Genomes pilot VCFs.

## Very long lines

The limits, by what reads the line:

- An input, cache variation or FASTA index line longer than 64 MiB is refused under the `input`, `cache` or `fasta` component (exit 1), with `longer than the 67108864-byte limit for one line` in the message. Ensembl VEP reads such a line whole.
- Such a line drops a plugin at load (the `Failed to instantiate plugin` warning, or the `plugin` error under `--safe`).
- A plugin query that reads it costs the queried variants that plugin's values, with a `Plugin '<Name>' went wrong` warning, `--safe` or not.
- `vep-cache-builder` refuses a GFF3, GTF, protein FASTA or core-dump line longer than 16 MiB under the `cache-builder` component, with `16777216-byte` in the message.
- `vep-cache-converter` refuses a plugin data line over 64 MiB.

## Bioconda

The Bioconda channel carries vep-rs as the package `vep-rs`, and the version it serves can trail the newest release: the channel builds a release once its own copy of the recipe names it ([Installation](installation.md#bioconda)). `vep --version` prints the release the installed package holds. Every release is published as pre-built binaries on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page, as source and as the container image ([Installation](installation.md)). A pipeline that installs Ensembl VEP from Bioconda does not get vep-rs by that package: its environment names `vep-rs` in place of `ensembl-vep`, or it runs vep-rs from the container image or the release archive ([Workflow managers](workflow-managers.md)). The repository holds a Bioconda recipe, which builds the binaries from a published release's source tarball ([Installation](installation.md#bioconda)).

To ask for a feature, or to report one this page leaves out, open an issue on the [issue tracker](https://github.com/natera-open-source/vep-rs/issues).
