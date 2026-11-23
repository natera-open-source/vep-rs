# `vep recode --config`: a `vep` option the recoder never declares

**The defect in one sentence.** `variant_recoder` reads a configuration source through `vep`'s option list, so a line naming a `vep` option its own command line never declares, `symbol 1` among them, is read and applied: the recoder computes what the option asks for, and `--fields` can print it.

**Why this is a defect.** A configuration file stands for the command line it replaces, and the recoder's command line is a fixed list of options (`variant_recoder`, lines 36-95), none of them `--symbol`. The file is instead tested against `vep`'s list of about two hundred options, which admits every `vep` setting the recoder's interface does not offer. A file written for `vep` therefore changes a recoder run in ways no documented option can, and `--fields gene_symbol` prints a field the tool's command line cannot ask for.

**What vep-rs does instead.** `vep recode` reads every configuration source as its own command line: an option of [its page](../cli/recode.md) may be set there, the first source naming it wins, and a name the page lacks is dropped from every source. The recoder's own presets are read alike on both tools. `fields`, `buffer_size`, `failed`, `merged`, `lrg` and `ambiguous_hgvs` keep the values the recoder settles before a source is read; `vcf_string`, `var_synonyms`, `mane_select` and `ga4gh_vrs` add no field from a source; `pretty` and `vcf_string` are names the file reader does not know.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):

- `variant_recoder`, lines 36-95: the recoder's own option list, with `config` (41); line 108 hands the parsed options to the recoder.
- `Bio/EnsEMBL/VEP/VariantRecoder.pm`, lines 96-109: the presets; line 111: the default field list; lines 131-150: the field flags append their fields from the command line alone; line 152: the configuration object is built after all of this.
- `Bio/EnsEMBL/VEP/Config.pm`, line 701: the `--config` file is read inside that constructor. Lines 1068-1110, `read_config_from_file`: each line's name is tested with `is_valid_param` (1100), which parses it against `@VEP_PARAMS` (lines 73-276, `vep`'s own options). An admitted scalar is stored with `||=` (1101), so a value already set stands and any other admitted name is set.

**What correcting it would change.** The file would set the recoder's own options and no others, as `vep recode` reads it.

**The class it produces.** `recode_config_file_vep_option`: the output of a `vep recode` run whose configuration source names a `vep` option the recoder's command line lacks, where `variant_recoder` applies the option and `vep recode` drops the line. The class is documented only: no committed recoder run names `--config`, and a comparison case that does names the class from its documentation. The concordance figures never see a recoder run.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A file naming `symbol`

**In short.** A one-line file, `symbol 1`, given to both tools with `--fields hgvsg,gene_symbol` and the identifier `rs781261918`. `variant_recoder` prints a `gene_symbol` list beside `hgvsg`: the file switched on `vep`'s `--symbol`, which its command line never offers. `vep recode` prints `hgvsg` and `input` alone: the line names no option of its page and is dropped, and a field nothing computes prints nothing.

### 1.1 What a reader needs to know first

Both tools translate a variant between its identifiers and notations and print one JSON array, one object per input line, keyed by allele, each allele's object holding `input` and the fields asked for. `--fields` may name any key of the annotation's JSON output, `gene_symbol` among them, and a key nothing computes prints no field. `symbol` is a `vep` option that writes gene symbols into that JSON; neither tool's command line declares it. `--config` names a file of `<name> <value>` lines, one option per line without the dashes.

A file naming the recoder's own options is read alike by both tools. `fields hgvsg`, `pretty 1` and `vcf_string 1` in a file leave both outputs as the command line alone gives them, the five default fields on one line. `mane_select 1` and `buffer_size 7` add nothing on either, since the recoder settles those options before a source is read.

`vep recode` reads the cache and a reference FASTA. `variant_recoder` reads the Ensembl database over the network and has no offline mode: its `lrg` preset is incompatible with `--offline`, so a file setting `offline 1` ends the run with `Can't use --lrg and --offline together`. The Ensembl VEP outputs below are database runs, and their notation lists carry the database's RefSeq models and co-located identifiers beside the cache's Ensembl ones.

### 1.2 The record

The identifier `rs781261918`, which the `GRCh38-recoder` golden corpus's cache (`tests/golden/116/GRCh38-recoder/`) resolves, `vep recode` reading that cache and the corpus's FASTA, and the file:

```
symbol 1
```

Both tools run it as `--id rs781261918 --config recode.ini --fields hgvsg,gene_symbol`, then, as the control, without the file.

### 1.3 What each engine writes

`variant_recoder` with the file:

```
[{"T":{"input":"rs781261918","gene_symbol":["JAM2","ATP5PF","LOC124905002"],"hgvsg":["NC_000021.9:g.25702257C>T"]}}]
```

`vep recode` with the file:

```
[{"T":{"hgvsg":["NC_000021.9:g.25702257C>T"],"input":"rs781261918"}}]
```

Without the file both tools print `hgvsg` and `input` alone; the three symbols are the genes of the transcripts the database annotates the variant on, which `--symbol` writes into the JSON and `--fields gene_symbol` collects.

### 1.4 Why Ensembl VEP writes that

`variant_recoder` parses its command line into a hash (lines 36-95) and hands it to its constructor (line 108), which builds the configuration object last (`VariantRecoder.pm` line 152). That object's constructor reads the `--config` file (`Config.pm` line 701). Each line's name goes through `is_valid_param` (line 1100), which parses `--<name> <value>` against `@VEP_PARAMS`, the option list of `vep` (lines 73-276), and stores an admitted value with `||=` (line 1101). `symbol` is on that list and nothing preset it, so the file sets it, the annotation run writes gene symbols into its JSON, and the recoder's field search finds `gene_symbol` there.

The same reader explains why the recoder's own options read alike on both tools. `pretty` and `vcf_string` are not on `vep`'s list and are skipped; `fields`, `buffer_size`, `failed`, `merged` and `lrg` were preset before the file was read (`VariantRecoder.pm` lines 96-111), so `||=` keeps them. `ambiguous_hgvs` is preset there too; `variant_recoder`'s own command line does not declare it, and `vep recode`'s does.

Source (`ensembl-vep`):

- `variant_recoder`, lines 36-95 and 108: the option list and the constructor call.
- `Bio/EnsEMBL/VEP/VariantRecoder.pm`, lines 96-111 and 152: the presets, the field list and the constructor call.
- `Bio/EnsEMBL/VEP/Config.pm`, lines 73-276, 701, 1035-1052 and 1068-1110: `@VEP_PARAMS`, the file read, `is_valid_param` and the `||=` store.

### 1.5 Why vep-rs's row is the correct one

A configuration file stands for the command line it replaces, and the recoder's command line has no `--symbol`. `vep recode` reads a source into the same parser as the command line, so a file can set what the command line can set and nothing else. A `vep` setting the tool never offers is dropped, as `vep` itself drops a name its own list lacks. The rule the two tools share, the command line winning over the file, and the recoder's presets, read alike on both, hold in both.

### 1.6 Try it yourself

Write the one-line file of 1.2 as `recode.ini` and run both tools on the identifier. `vep recode` reads the corpus's cache; `variant_recoder` reads the Ensembl database (`--db_version 116`):

```
vep recode --id rs781261918 --config recode.ini --fields hgvsg,gene_symbol \
  --json_cache <corpus>/json_cache --fasta <corpus>/reference.fa.gz \
  --synonyms <corpus>/chr_synonyms.txt
```

| Run | `variant_recoder` | `vep recode` |
| --- | --- | --- |
| `--config recode.ini` | `hgvsg`, `gene_symbol`, `input` | `hgvsg`, `input` |
| no file | `hgvsg`, `input` | the same |

### 1.7 How the comparator treats this class

The class is documented only. No committed recoder run names `--config`, so the tool comparison (`tools/vep-formats-compare` over the recoder cases) meets no output shaped by a file, and a comparison case that names one records the class from its documentation.
