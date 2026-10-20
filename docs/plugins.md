# Plugins {#plugin-guide}

vep-rs registers 11 built-in annotation plugins and loads C-ABI shared-library
plugins built against the `vep-plugin` crate without a rebuild. A plugin adds
fields to the output: the Extra column in VEP format, CSQ in VCF, additional
fields in JSON. This page lists the plugins, shows how one is enabled and points
to the page that documents each.

Support is per reference genome, because each annotation database decides which
builds it publishes. Eight plugins are supported on both GRCh37 and GRCh38. GWAS
is GRCh38-only, SpliceAI is GRCh38-SNV-only, and dbNSFP is not supported. Each
restriction is a property of the upstream database rather than of vep-rs; the
plugin's page gives the reason.

Pass `--assembly GRCh37` or `--assembly GRCh38` whenever plugins are enabled. The
flag is load-bearing, not just cache selection. REVEL uses it to choose which
coordinate column its index was built on. CADD, REVEL, dbscSNV and AlphaMissense
use it to refuse a data file whose own assembly contradicts the run, rather than
silently annotating nothing. Those guards are inert without it, and REVEL falls
back to the GRCh37 column with only a warning.

## Available plugins

The pages group the built-in plugins by mechanism, as
`crates/vep-builtins/README.md` does: how a plugin reads its data and matches a
variant decides which data file it needs and where its fields land. The Builds
column states which assemblies each upstream database publishes, and each
plugin name links to its entry.

| Plugin                                                     | Mechanism                                                            | Builds         |
| ---------------------------------------------------------- | -------------------------------------------------------------------- | -------------- |
| [CADD](plugins/annotation-store.md#cadd)                   | annotation store; position and allele                                | GRCh37, GRCh38 |
| [REVEL](plugins/annotation-store.md#revel)                 | annotation store; position, allele and transcript                    | GRCh37, GRCh38 |
| [gnomADc](plugins/annotation-store.md#gnomadc)             | annotation store; position only                                      | GRCh37, GRCh38 |
| [AlphaMissense](plugins/annotation-store.md#alphamissense) | annotation store; position and allele                                | GRCh37, GRCh38 |
| [dbscSNV](plugins/annotation-store.md#dbscsnv)             | annotation store; position and allele, per transcript                | GRCh37, GRCh38 |
| [LoFTEE](plugins/loftee.md)                                | consequence filter over the transcript model, splice motifs and GERP | GRCh37, GRCh38 |
| [LoFtool](plugins/flat-file.md#loftool)                    | flat file; gene symbol                                               | GRCh37, GRCh38 |
| [pLI](plugins/flat-file.md#pli)                            | flat file; gene symbol                                               | GRCh37, GRCh38 |
| [GWAS](plugins/flat-file.md#gwas)                          | flat file; position and risk allele                                  | GRCh38         |
| [SpliceAI](plugins/spliceai.md)                            | annotation store; position, allele and gene symbol                   | GRCh38 SNV     |
| [dbNSFP](#not-supported-dbnsfp)                            | annotation store; position, allele and transcript                    | not supported  |

Preparing the data files, and the `.vpd` store that replaces a tabix index, is
on [Preparing data files](plugins/data-files.md).

## Usage

Plugins are loaded via `--plugin Name,param1,param2,...` and can be specified
multiple times:

```bash
SNV=/path/to/cadd/whole_genome_SNVs.tsv.gz
INDELS=/path/to/cadd/gnomad.genomes.r4.0.indel.tsv.gz
vep -i input.vcf -o output.txt --offline --cache \
    --json_cache /path/to/json_cache --assembly GRCh38 \
    --plugin CADD,snv=$SNV,indels=$INDELS \
    --plugin REVEL,file=/path/to/revel/new_tabbed_revel_grch38.tsv.gz \
    --plugin gnomADc,/path/to/gnomad/gnomad.exomes.coverage.tsv.gz
```

Each plugin's parameters are on its page.

## Not supported {#not-supported-dbnsfp}

`--plugin dbNSFP` resolves to a built-in plugin, but dbNSFP is **not a supported
plugin**: its license grants free use only for academic, non-commercial purposes
and requires a paid commercial agreement for any other use, so vep-rs neither
distributes the data nor validates against it. A user who holds a dbNSFP license
can point the plugin at their own tabix-indexed file and receive user-selected
columns prefixed `dbNSFP_`. Native SIFT and PolyPhen decoding needs no dbNSFP
file (see [Cache setup](cache-setup.md#sift-and-polyphen-predictions)).

The parameters are the file path followed by the column names to extract:

```bash
--plugin dbNSFP,/path/to/dbNSFP.gz,\
SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred
```

## Writing a plugin

How the two plugin systems work, the C ABI a shared-library plugin exports, and
how a new built-in plugin is added and tested are on
[Writing a plugin](developers/writing-a-plugin.md).
