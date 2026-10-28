# CLI reference

The `vep` binary accepts the Ensembl VEP command line. Every flag is Ensembl VEP's and does what Ensembl VEP's does, or is refused with an error naming what it would do; the few flags vep-rs adds are marked as its own. A flag the error table governs (refused, or accepted without effect) carries that table's text as its help.

These pages are generated from the command-line definition and the option tables by the `cli_reference` test of `crates/vep-cli`. Edit the sources, then run `VEP_UPDATE_CLI_REFERENCE=1 cargo test -p vep-cli --test cli_reference` to rewrite them.

```bash
# Show all flags
vep --help

# Short form
vep -i input.vcf -o output.txt --offline --json_cache /path/to/json_cache
```

Worked examples are on [Examples](examples.md); the output formats and their parity with Ensembl VEP on [Output formats](output-formats.md); the divergences from Ensembl VEP that a flag's description names on [Intended divergences](intended-divergences.md).

## Option groups

Every flag has a row on one of these pages, grouped as `--help` groups them: its names, its value type, its default and its description. A default of `false` is a switch that is off; `none` marks a flag with no declared default, and where a value applies in its absence the description names it.

| Group | Covers |
| --- | --- |
| [Input and output](cli/input-output.md) | The input file, its format and how it is read, the output file, its format and compression, the header lines, the warnings and skipped-variants files. |
| [Cache, species and runtime](cli/cache-species-runtime.md) | The cache directory and its layout, the reference FASTA, species and assembly, the batch size, thread counts and distances. |
| [Output control](cli/output-control.md) | The field flags, `--everything` among them, that choose which annotation columns appear. |
| [Allele frequencies and co-located variants](cli/allele-frequencies.md) | Co-located known variants, population frequencies, the frequency filter, PubMed identifiers and gene-phenotype associations. |
| [Filtering, picking and transcript sets](cli/filtering-picking.md) | The pick and flag-pick families, the allele columns, the transcript filter and the transcript-set flags. |
| [Plugins](cli/plugins.md) | Loading built-in and dynamic-library plugins. |
| [Statistics, verbosity and miscellaneous](cli/statistics-verbosity-misc.md) | The statistics flags, logging, and the flags that change how a record is handled. |
| [Compatibility with Ensembl VEP](cli/compatibility.md) | The database flags, the flags refused or accepted without effect, where `Config.pm` is not followed, the configuration sources, the option spellings and the exit status. |

A repeated single-value flag keeps its last value, as Getopt::Long does. `-o` accumulates and the last wins with a warning, as the `vep` script does (vep:32-33). A long option is also read in any case, as a unique prefix and after one dash, and from the `--config` file, `vep.ini` and the `VEP_*` variables ([configuration sources](cli/compatibility.md#configuration-sources)).

`--help` prints the flag list and exits 0, as does a command line with no argument at all (vep:29); `--version` prints `vep <version>` and exits 0.

## Option rules

Ensembl VEP's option rules are applied as data, from `Config.pm:327-645` at release/116.2, in this order after parsing. A flag the error table refuses is refused first. Then the option sets switch on the flags they imply, in the order below, a set's value overwriting an explicit one (`apply_option_sets`, `Config.pm:816-884`). Then the checks of `check_config` (`Config.pm:944-985`) run: valid values, incompatible pairs, requirements, deprecations, each refused with Ensembl VEP's own message. A rule naming a flag vep-rs does not declare cannot fire and is listed so the table reads as Ensembl VEP's.

### Option sets (`Config.pm:351-593`)

| Flag | Switches on |
| --- | --- |
| `--everything` | `--sift b`, `--polyphen b`, `--ccds`, `--hgvs`, `--symbol`, `--numbers`, `--domains`, `--regulatory`, `--canonical`, `--protein`, `--biotype`, `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`, `--uniprot`, `--mane`, `--tsl`, `--appris`, `--variant_class`, `--gene_phenotype`, `--mirna` |
| `--hgvs` | `--hgvsc`, `--hgvsp` |
| `--genomes` | `--host mysql-eg-publicsql.ebi.ac.uk`, `--port 4157` |
| `--grch37` | `--port 3337`, `--assembly GRCh37` |
| `--refseq` | `--core_type otherfeatures` |
| `--individual` | `--allow_non_variant` |
| `--individual_zyg` | `--allow_non_variant` |
| `--cell_type` | `--regulatory` |
| `--regulatory_gff` | `--regulatory` |
| `--filter_common` | `--check_frequency` |
| `--check_frequency`, `--af`, `--af_1kg`, `--af_gnomad`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed` | `--check_existing` |
| `--show_cache_info` | `--offline`, `--format ensembl` |
| `--offline` | `--cache`, `--database 0` |
| `--cache` | `--no_slice_cache` |
| `--convert` | `--no_stats` |
| `--json` | `--rest`, `--output_format json` |
| `--rest` | `--no_escape` |
| `--vcf` | `--output_format vcf`, `--symbol`, `--biotype`, `--numbers` |
| `--tab` | `--output_format tab` |
| `--humdiv` | `--polyphen_analysis humdiv` |
| `--minimal` | `--allele_number` |
| `--gff` | `--custom file=%gff%,format=gff` |
| `--gtf` | `--custom file=%gtf%,format=gtf` |
| `--bigwig` | `--custom %bigwig%,,bigwig,exact` |
| `--phyloP` | `--custom %ucsc_data_root%%ucsc_assembly%/phyloP%phyloP%way/%ucsc_assembly%.phyloP%phyloP%way.bw,phlyoP%phyloP%way,bigwig,exact` |
| `--phastCons` | `--custom %ucsc_data_root%%ucsc_assembly%/phastCons%phastCons%way/%ucsc_assembly%.phastCons%phastCons%way.bw,phastCons%phastCons%way,bigwig,exact` |
| `--exclude_predicted` | `--transcript_filter not stable_id match ^X._` |
| `--bam` | `--use_transcript_ref`, `--bam_edited` |
| `--mane` | `--mane_select`, `--mane_plus_clinical` |
| `--use_given_ref` | `--use_transcript_ref 0` |

`--everything` therefore reaches `--hgvsc` and `--hgvsp` through `--hgvs`, `--mane_select` and `--mane_plus_clinical` through `--mane`, and `--check_existing` through the frequency flags; it does not switch on `--allele_number` or `--hgvsg`. `--regulatory` is switched on and emits RegulatoryFeature and MotifFeature rows when the converted cache carries regulatory data. Motif scores, high-information positions and cell-type activity are omitted.

### Incompatible flags (`Config.pm:615-632`)

`Can't use --<flag> and --<other> together`, naming the first pair in this table's order (Ensembl VEP iterates a hash there, so which pair it names varies from run to run); Ensembl VEP skips the check under `--safe`.

| Flag | Cannot be combined with |
| --- | --- |
| `--most_severe` | `--biotype`, `--no_intergenic`, `--protein`, `--symbol`, `--sift`, `--polyphen`, `--coding_only`, `--ccds`, `--mane`, `--mane_select`, `--canonical`, `--xref_refseq`, `--numbers`, `--domains`, `--tsl`, `--appris`, `--uniprot`, `--pick`, `--flag_pick`, `--pick_allele`, `--flag_pick_allele`, `--vcf`, `--summary` |
| `--summary` | `--biotype`, `--no_intergenic`, `--protein`, `--symbol`, `--sift`, `--polyphen`, `--coding_only`, `--ccds`, `--mane`, `--mane_select`, `--canonical`, `--xref_refseq`, `--numbers`, `--domains`, `--tsl`, `--appris`, `--uniprot`, `--pick`, `--flag_pick`, `--pick_allele`, `--flag_pick_allele`, `--vcf`, `--most_severe` |
| `--database` | `--af_1kg`, `--af_gnomad`, `--af_gnomade`, `--af_gnomadg`, `--af_exac`, `--max_af`, `--pubmed`, `--var_synonyms`, `--offline`, `--cache` |
| `--af_gnomade` | `--af_gnomad` |
| `--quiet` | `--verbose` |
| `--refseq` | `--gencode_basic`, `--merged` |
| `--gencode_primary` | `--gencode_basic`, `--refseq` |
| `--json` | `--vcf`, `--tab` |
| `--vcf` | `--json`, `--tab` |
| `--tab` | `--vcf`, `--json` |
| `--individual` | `--minimal`, `--individual_zyg` |
| `--check_ref` | `--lookup_ref` |
| `--check_svs` | `--offline` |
| `--ga4gh_vrs` | `--vcf` |
| `--shift_hgvs` | `--shift_3prime`, `--shift_genomic` |
| `--lrg` | `--offline` |

Pairs Ensembl VEP refuses that are accepted here:

| Flag | Beside | Why the pair is accepted |
| --- | --- | --- |
| `--check_svs` | `--offline` | Ensembl VEP reads the overlapping structural variants from its database, which `--offline` forgoes (`AnnotationSourceAdaptor.pm:203-206`); `--check_svs <file>` reads them from Ensembl's structural-variant VCF and needs no database, so the pair is accepted and the file is read offline |

### Required flags (`Config.pm:607-612`)

`You must set --<required> to use --<flag>`.

| Flag | Needs |
| --- | --- |
| `--original` | `--filters` |
| `--phyloP` | `--ucsc_assembly` |
| `--phastCons` | `--ucsc_assembly` |
| `--custom_multi_allelic` | `--custom` |

### Valid values (`Config.pm:596-604`)

A declared flag carries its list as a value enumeration, so a value outside it is an argument syntax error (exit 2) listing the possible values.

| Flag | Values |
| --- | --- |
| `--format` | `ensembl`, `vcf`, `hgvs`, `id`, `spdi`, `region`, `guess` |
| `--terms` | `SO`, `display`, `NCBI` |
| `--sift` | `s`, `p`, `b` |
| `--polyphen` | `s`, `p`, `b` |
| `--pick_order` | `mane_select`, `mane_plus_clinical`, `canonical`, `appris`, `tsl`, `biotype`, `ccds`, `rank`, `length`, `ensembl`, `refseq` |
| `--nearest` | `transcript`, `gene`, `symbol` |
| `--compress_output` | `gzip`, `bgzip` |

### Deprecated flags (`Config.pm:636-638`)

- `--convert`: `--convert has been deprecated`

### Lists (`Config.pm:327-346`)

- Comma-separated values: `--individual`, `--individual_zyg`, `--pick_order`, `--fields`, `--chr`, `--distance`. `--cell_type`, comma-separated in Ensembl VEP, is an option this version lacks, refused on the command line and in a configuration file (exit 2). `--dont_export`, comma-separated in Ensembl VEP's list, is an option neither its `vep` script (`Config.pm:73-276`) nor this version declares: unknown on the command line (exit 2) and dropped from a configuration file, as in Ensembl VEP.
- Repeatable, each occurrence appending: `--plugin`, `--custom`, `--transcript_filter`. `--phyloP` and `--phastCons`, repeatable in Ensembl VEP, are options this version lacks, refused on the command line and in a configuration file (exit 2).
- Input formats Ensembl VEP exempts from its sorted-input check (`Config.pm:640-645`): `id`, `hgvs`, `spdi`, `region`. vep-rs performs no such check in any format: every record is annotated in input order, and an input whose chromosomes interleave keeps more transcript indexes in memory.

### Rules of vep-rs's own

Checked after Ensembl VEP's rules on the same state, for behaviour Ensembl VEP has no flag for; each is refused with the line below.

| Component | Error |
| --- | --- |
| `output` | `ERROR: [output] --output_format parquet needs a file path: Parquet cannot be written to STDOUT` |
| `output` | `ERROR: [output] --fields applies to --tab and --vcf output; the Parquet column set is fixed` |
| `output` | `ERROR: [output] --compress_output applies to the default, tab, VCF and JSON outputs; Parquet output is compressed by its own codec` |

