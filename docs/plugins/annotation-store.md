# Annotation-store plugins

The parameters, commands and data preparation of CADD, REVEL, gnomADc,
AlphaMissense, dbscSNV, AVI, OpenTargets, SpliceVault, TranscriptAnnotator,
gnomADMt and LOEUF.

An annotation-store plugin reads a bgzipped, tabix-indexed data file, or the
`.vpd` binary store built from it. It matches records by position and then by
whatever its Matches row adds.

Pass `--assembly` whenever one of these plugins is enabled. CADD, REVEL, dbscSNV
and AlphaMissense use it to refuse a data file built on the other assembly.
REVEL and dbscSNV use it to pick the coordinate column their index was built on.

## CADD {#cadd}

CADD adds Combined Annotation Dependent Depletion scores to a variant. The SNV
and indel files score sequence variants, and the CADD-SV file scores structural
variants (`<DEL>`, `<INS>` and `<DUP>` records) by their class and start.

| Property      | Value                                                                                        |
| ------------- | -------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38 (CADD-SV on GRCh38)                                                           |
| Matches       | position and allele, after the bases the alleles share at either end are trimmed; a structural variant by its class word (`DEL`, `INS`, `DUP`) and the base before its start |
| Output fields | `CADD_PHRED`, `CADD_RAW`, printed on every row of the record                                 |
| Data source   | [CADD download](https://cadd.gs.washington.edu/download); [CADD-SV download](https://kircherlab.bihealth.org/download/CADD-SV/v1.1/) |

| Parameter          | Meaning                                                                                                                                                |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `snv=<path>`       | The SNV score file.                                                                                                                                    |
| `indels=<path>`    | The indel score file.                                                                                                                                  |
| `sv=<path>`        | The CADD-SV file.                                                                                                                                      |
| `force_annotate=1` | Annotate structural variants from the SNV and indel files too, and annotate a variant whose region the files answer with more than 100,000 rows. Off by default. |

The parameters are either every one `key=value` or every one a file path by
position, where a trailing `1` is `force_annotate`, as Ensembl VEP's plugin
reads them. A structural variant is not annotated while a file is loaded under
`snv=` or `indels=` unless `force_annotate=1` is set. A run that sets it prints
Ensembl VEP's warning that sequence files slow a structural-variant run. A
region the files answer with more than 100,000 rows is left unannotated with a
warning unless `force_annotate=1` is set. A file carries the fields its header
names (`CADD.pm` lines 81-98 and 171-175): `PHRED` and `RawScore`, or
`CADD-SV_PHRED-score` and `CADD-SV_Raw-score`; the first file whose table row
matches writes them.

```bash
--plugin CADD,snv=/path/to/snv.tsv.gz,indels=/path/to/indels.tsv.gz
--plugin CADD,sv=/path/to/1000G_phase3_SVs.tsv.gz
--plugin CADD,snv=/path/to/snv.tsv.gz,indels=/path/to/indels.tsv.gz,\
force_annotate=1
--plugin CADD,/path/to/snv.tsv.gz,/path/to/indels.tsv.gz,1
```

**Preparing the data.** Ensembl VEP's `CADD.pm` names two downloads (lines
45-49): the SNV and indel data files with their tabix indexes at
[CADD download](https://cadd.gs.washington.edu/download), and the
structural-variant data files with theirs at [CADD-SV
download](https://kircherlab.bihealth.org/download/CADD-SV/v1.1/). Its examples
name `whole_genome_SNVs.tsv.gz`, `InDels.tsv.gz` and `1000G_phase3_SVs.tsv.gz`
(lines 31-33). vep-rs reads a structural-variant file's index as the index
declares it.

Under `--assembly`, a file whose header lines do not name the assembly is
refused, as Ensembl VEP's plugin refuses it (`CADD.pm` lines 151-157).

## REVEL {#revel}

REVEL adds the Rare Exome Variant Ensemble Learner score to a missense row
whose amino-acid change the table names.

| Property      | Value                                                                                                        |
| ------------- | ------------------------------------------------------------------------------------------------------------ |
| Builds        | GRCh37, GRCh38                                                                                               |
| Matches       | missense rows only: position and allele, the record's `aaref` and `aaalt` equal, as written, to the row's reference and alternate amino acids, and one of the record's transcript ids equal to the row's transcript unless `no_match` is set or the file carries no transcript ids. Ensembl VEP's `REVEL.pm` tests the alternate amino acid alone: [the score of a table row of another reference amino acid](../divergences/revel-other-reference-amino-acid-rows.md) |
| Output fields | `REVEL`                                                                                                      |
| Data source   | [REVEL download](https://sites.google.com/site/revelgenomics/downloads)                                     |

The parameters are read by key or by position, never mixed: the file first, then
`no_match`.

| Parameter          | Meaning                                                                                                            |
| ------------------ | ------------------------------------------------------------------------------------------------------------------ |
| `file=<path>`      | The REVEL score file (or the first parameter).                                                                     |
| `no_match=<value>` | Write the score without matching the record's transcripts to the row's (or the second parameter); every value but an empty one and `0` sets it. |

```bash
--plugin REVEL,file=/path/to/revel.tsv.gz
--plugin REVEL,file=/path/to/revel.tsv.gz,no_match=1
--plugin REVEL,/path/to/revel.tsv.gz,1
```

The header's column count fixes the layout, whatever the header names. Seven
columns are the GRCh37-only file, eight and nine the file with both assemblies'
coordinates, a row's ninth cell its transcript ids. Any other count is refused,
as is a header without the run's position column (`hg19_pos` under GRCh37,
`grch38_pos` under GRCh38) and a copy indexed on the other assembly's column. A
header line lacking its `#` is refused, as Ensembl VEP refuses it (`REVEL.pm`
lines 116-125); the recipe below adds it.

The score is written as the file spells it. A run under `--refseq` or `--merged`
without `no_match` is warned on standard error alone, under `--quiet` too, that
transcript matching may fail, as Ensembl VEP does.

**Preparing the data.** Download from
[REVEL download](https://sites.google.com/site/revelgenomics/downloads) and
process:

```bash
unzip revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv

# GRCh37:
tabix -f -s 1 -b 2 -e 2 new_tabbed_revel.tsv.gz

# GRCh38:
zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v '^#chr' new_tabbed_revel.tsv.gz | awk '$3 != "." ' \
    | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz
```

One download serves both assemblies, indexed on column 2 for GRCh37 and column 3
for GRCh38. The run's assembly selects the column, from `--assembly` or else the
cache, and the plugin is refused when neither names one. Ensembl VEP 116.2
builds its plugins before it reads the cache (`Runner.pm` lines 125-128). A run
without `--assembly` therefore drops REVEL there with a warning (`REVEL.pm`
lines 132-133), or ends under `--safe`.

## gnomADc {#gnomadc}

gnomADc adds gnomAD coverage statistics to a variant.

| Property      | Value                                                                                                   |
| ------------- | ------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                          |
| Matches       | position only                                                                                           |
| Output fields | `gnomAD_mean_cov`, `gnomAD_median_cov` or `gnomAD_median_approx_cov` and `gnomAD_total_DP_cov`, and `gnomAD_<N>x_cov` for N in 1, 5, 10, 15, 20, 25, 30, 50, 100 |
| Data source   | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads)                                        |

| Parameter | Meaning                                  |
| --------- | ---------------------------------------- |
| `<path>`  | The coverage file, as the only parameter. |

```bash
--plugin gnomADc,/path/to/gnomad.coverage.tsv.gz
```

The fields are Ensembl VEP's, keyed by the file's column count, not its header
names. Fourteen columns, the split GRCh38 v3 and v4 layout, yield
`gnomAD_median_approx_cov` and `gnomAD_total_DP_cov`. Any other count, GRCh37
v2.1's thirteen among them, yields `gnomAD_median_cov` and no total depth.
Values are read by position from the third column on.

Every value is the mean over the positions the variant spans, printed with four
decimals. A single-base variant prints its own table row's values, `1.1177e-01`
as `0.1118` and a depth of `81695` as `81695.0000`. A deletion averages the
table rows of its deleted bases; an insertion reads the table row at the base
after the insertion point. A variant with no table row in the file prints no
field. A record written `Chr21` reads the file's `chr21` contig, as every
tabix-backed plugin's lookup does here; Ensembl VEP's finds nothing for it ([the
capitalised
prefix](../plugins.md#4-a-chromosome-written-with-a-capitalised-chr-prefix-in-a-plugins-file-lookup)).

Under `--minimal`, a split allele reads its own trimmed span, the one its row's
alleles describe; Ensembl VEP's rejoin hands the plugin the record's span
([a documented class](../divergences/minimal-plugin-sees-the-rejoined-record.md)).

**Preparing the data.** Download coverage files from [gnomAD
downloads](https://gnomad.broadinstitute.org/downloads). Split the GRCh38 v3 and
v4 files' single `locus` column (`chr1:11819`) into chromosome and position
first (`sed "1s/locus/chr\tpos/; s/:/\t/g"`). Build the unpublished tabix index
on the chromosome and position columns (`bgzip`, then `tabix -s 1 -b 2 -e 2`).
The header line works with or without a leading `#`.

## AlphaMissense {#alphamissense}

AlphaMissense adds DeepMind's AlphaMissense pathogenicity prediction to the
missense transcript consequences of a variant, as Ensembl VEP's
`AlphaMissense.pm` does.

| Property      | Value                                                                                                                   |
| ------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                          |
| Matches       | missense rows only: a record from two bases before the variant's start to its end whose alleles, trimmed of the bases they share, are the variant's at its position, and the first and last letters of the record's `protein_variant` equal to the row's reference and alternate amino acids; the record's transcript equal to the row's under `transcript_match=1`. On a multi-allelic record each allele's rows are matched against that allele alone, where Ensembl VEP 116.2 matches every allele of the record ([The AlphaMissense prediction of a sibling allele on a multi-allelic record](../divergences/alphamissense-sibling-allele-rows.md)) |
| Output fields | `am_pathogenicity`, `am_class` by default; any column of the file under `cols`, each named with the `am_` prefix          |
| Data source   | [AlphaMissense data](https://console.cloud.google.com/storage/browser/dm_alphamissense)                                |

Every parameter is `key=value`; a file given by position is refused.

| Parameter                  | Meaning                                                                                                          |
| -------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `file=<path>`              | The prediction file (mandatory). The published tables name the mitochondrial contig `chrM`; vep-rs reads them for a record spelt `MT`, `chrM` or `chrMT`, where Ensembl VEP 116.2 reaches them from a record spelt `chrM` or `M` alone ([A mitochondrial record against a plugin's file naming the contig under another spelling](../divergences/plugin-contig-mitochondrial-spelling.md)). |
| `cols=<name>[:<name>...]`  | The columns written, by their output names or the file's own; `cols=all` writes every column. A name the header lacks is refused with the valid names listed; Ensembl VEP 116.2 tests each name as a pattern ([a documented class](../divergences/plugin-cols-names-matched-exactly.md)). |
| `transcript_match=<value>` | Write a prediction only when its transcript is the row's; every value but an empty one and `0` sets it.          |

```bash
--plugin AlphaMissense,file=/path/to/AlphaMissense_hg38.tsv.gz
--plugin AlphaMissense,file=/path/to/AlphaMissense_hg38.tsv.gz,cols=all
--plugin AlphaMissense,file=/path/to/AlphaMissense_hg38.tsv.gz,\
    transcript_match=1
```

The header describes each column as `column from <file>`, the published columns
with their own text first. In the JSON format the fields nest under
`alphamissense`. The transcript comparison of `transcript_match` sets the
version suffix aside on both sides, so a prediction on `ENST00000000233.10`
matches the transcript `ENST00000000233`. Ensembl VEP's plugin removes one digit
of the version and never matches a two-digit one.

A run under `--refseq` or `--merged` with `transcript_match=1` is warned that
the match may fail, since the file names Ensembl transcripts. The warning is
printed on standard error, under `--quiet` too, and never written to the
warnings file, as Ensembl VEP prints it.

**Preparing the data.** AlphaMissense publishes one file per assembly, and every
data row names its assembly in the `genome` column. A file whose rows name the
other assembly is refused under `--assembly`.

## dbscSNV {#dbscsnv}

dbscSNV adds AdaBoost and random-forest splice-site predictions to a transcript
consequence.

| Property      | Value                                                                                                                     |
| ------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                            |
| Matches       | position and alternate allele of a single-base substitution; the scores print on every transcript row of the variant   |
| Output fields | `ada_score`, `rf_score`                                                                                                   |
| Data source   | [dbscSNV 1.1](https://sites.google.com/site/jpopgen/dbNSFP), the `dbscSNV1.1.zip` download prepared per assembly as the Ensembl VEP plugin documents (`dbscSNV.pm` lines 53-59)          |

| Parameter   | Meaning                                                                                  |
| ----------- | ---------------------------------------------------------------------------------------- |
| `<path>`    | The prepared file, as the first parameter.                                               |
| `<assembly>` | `GRCh37` or `GRCh38` as the second parameter: the assembly the file was prepared for, in place of the run's. An empty second parameter leaves the run's, as Ensembl VEP's plugin reads it; a later positional parameter is not read, and is noted in the warnings file. |

```bash
--plugin dbscSNV,/path/to/dbscSNV1.1_GRCh38.txt.gz
--plugin dbscSNV,/path/to/dbscSNV1.1_GRCh37.txt.gz,GRCh37
```

**Preparing the data.** Prepare the download as Ensembl VEP's `dbscSNV.pm`
documents (lines 45-62). The file holds both assemblies' coordinates (lines
46-47), and tabix is told that its header line starts with `c` (lines 61-62).
Ensembl VEP's plugin strips a carriage return from a row's end (line 191).

For GRCh37, concatenate the `dbscSNV1.1.chr*` files under the header and index
on columns 1 and 2 (`tabix -s 1 -b 2 -e 2 -c c`). For GRCh38, sort them on
columns 5 and 6 and index on those (`tabix -s 5 -b 6 -e 6 -c c`). The plugin
reads `pos` under GRCh37 and `hg38_pos` under GRCh38 (`dbscSNV.pm` lines
215-219), and `ada_score` and `rf_score` by header name.

The two preparations share the header and differ in their index (`dbscSNV.pm`
lines 51-59). vep-rs therefore checks that the index was built on the
assembly's position column and refuses the other preparation, naming both
columns. A file that keeps only the columns the plugin reads
(`chr pos ref alt ada_score rf_score`, or `hg38_chr hg38_pos ...` for GRCh38)
is read the same way.

The assembly comes from the second parameter, else `--assembly`, else the cache;
a run naming none is refused. Ensembl VEP 116.2 builds the plugin before it
reads the cache (`Runner.pm` lines 125-128). A run without the flag or the
parameter is therefore refused there (`dbscSNV.pm` line 226) and the plugin
dropped with a warning.

## AVI {#avi}

AVI adds the AlphaGenome Variant Impact score of a single-base substitution to
a variant. The values print on every row of the variant, transcript, regulatory
and intergenic rows alike, as the Ensembl VEP plugin's `Feature` and
`Intergenic` feature types do.

| Property      | Value                                                                                   |
| ------------- | --------------------------------------------------------------------------------------- |
| Builds        | GRCh38                                                                                  |
| Matches       | position and allele, the table's letter case disregarded; one base replaced by one base only; under `--minimal` a split allele over its own trimmed span, where Ensembl VEP reads the rejoined record and scores no allele of it ([a documented class](../divergences/minimal-plugin-sees-the-rejoined-record.md)) |
| Output fields | `AVI_PHRED`; `AVI_RAW` under `raw=1`; `AVI_LINK` under `links=1`                        |
| Data source   | [AlphaGenome downloads](https://deepmind.google.com/science/alphagenome/downloads)      |

| Parameter     | Meaning                                                                        |
| ------------- | ------------------------------------------------------------------------------ |
| `file=<path>` | The score table; the path may also be the first positional parameter. A value Perl reads as false (`0` or the empty string) is no file, and the plugin is refused without one. |
| `raw=1`       | Also write the raw score. `1`, `true` and `yes` in any case turn a flag on.   |
| `links=1`     | Also write a link to the variant in the AlphaGenome Atlas.                     |

```bash
--plugin AVI,file=/path/to/avi_scores.tsv.gz
--plugin AVI,file=/path/to/avi_scores.tsv.gz,raw=1,links=1
```

`AVI_PHRED` is the PHRED-scaled quantile of the raw score, 0 to about 80. A
value of 20 places the substitution among the top 1 percent most deleterious
substitutions genome-wide. `AVI_RAW` runs from about -1.3 to +4.0; higher means
more deleterious in both (`AVI.pm` lines 53-58). `AVI_LINK` names the matching
table row's own chromosome, position and alleles in the Atlas URL, with `chr`
prefixed when the table's chromosome names lack it. An insertion, a deletion, a
substitution of more than one base, or one whose reference or alternate base is
`N` takes no score.

**Preparing the data.** The AVI SNV score table is a permissive-use artifact
under the AlphaGenome terms: commercial and non-commercial use. Ensembl VEP's
plugin reads a tabix-indexed, bgzip-compressed TSV (`AVI.pm` lines 37-38) by
column position, six cells to a row (line 161), so no header is needed:

| Column | Content |
| --- | --- |
| 1 | chromosome |
| 2 | position |
| 3 | reference allele |
| 4 | alternate allele |
| 5 | raw score |
| 6 | PHRED score |

The plugin's build is GRCh38, and none of the six cells the plugin reads names
an assembly. A GRCh37 run is not refused, only warned, and matches table rows
by position alone.

## OpenTargets {#opentargets}

OpenTargets adds the Open Targets Platform's GWAS and QTL associations of a
variant. These are the locus-to-gene score, gene and disease of each GWAS
credible set the variant leads, and the gene and biosample of each QTL.

| Property      | Value                                                                                                                   |
| ------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh38                                                                                                                  |
| Matches       | position and allele, at the variant level; under `--minimal` a split allele by its own trimmed form, where Ensembl VEP reads the rejoined record and writes a matched record's fields on every allele's rows ([a documented class](../divergences/minimal-plugin-sees-the-rejoined-record.md)) |
| Output fields | `OpenTargets_<column>` for every selected column; the defaults are listed under `cols`                                  |
| Data source   | [Open Targets Platform downloads](https://platform.opentargets.org/downloads), the VEP export under `etc/export/vep_data/` (CC0 1.0) |

| Parameter                | Meaning                                                                                                                                                                                                      |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `file=<path>`            | The tabix-indexed export, with its header line commented (`#`).                                                                                                                                              |
| `cols=<a:b:c>`           | The header columns to write, `:`-separated, or `all`. Default: `gwasLocusToGeneScore`, `gwasGeneId` and `gwasDiseases` when GWAS study types are selected, `qtlGeneId` and `qtlBiosampleName` when QTL types are. A name the header lacks is refused with the valid names listed; Ensembl VEP 116.2 tests each name as a pattern, so `cols=gwasGeneI.` loads there with every value `NA` ([a documented class](../divergences/plugin-cols-names-matched-exactly.md)). |
| `study_types=<list>`     | `all` (default), `gwas`, `QTL` (every study type ending in `qtl`) or a `:`-separated list of such names, `eqtl:pqtl` among them. A list naming a QTL type beside `QTL` is refused. A list naming `QTL` twice, which Ensembl VEP 116.2 dies on, is [a documented class](../divergences/opentargets-duplicate-study-type-and-missing-islead.md). |
| `lead_variants_only=0|1` | `1` (default) keeps the records whose `isLead` is `true`. An export whose header lacks `isLead` (or lacks `studyType` under a `study_types` filter) is refused at load, where Ensembl VEP 116.2 prints the header lines and annotates nothing: [the same class](../divergences/opentargets-duplicate-study-type-and-missing-islead.md). |

```bash
--plugin OpenTargets,file=/path/to/open_targets_vep_full.tsv.gz
--plugin OpenTargets,file=/path/to/open_targets_vep_gwas.tsv.gz,\
study_types=gwas,cols=gwasGeneId:gwasLocusToGeneScore:gwasDiseases:studyId
```

A variant that matches several table rows prints each field's values `,`-joined
in file order. `NA` stands where a table row has no value in that column (a QTL
table row under a GWAS column, and the reverse). The JSON output instead carries
one object per table row under `opentargets`, with the selected columns as keys.
The header line describes every field as a column of the file, the locus-to-gene
score with its meaning first. A table row's VCF-style insertion or deletion (an
anchor base shared by both alleles) matches the variant as the input states it.

An export holding one family alone, the Platform's `open_targets_vep_gwas.tsv.gz`
or `open_targets_vep_qtl.tsv.gz` or a file split as Ensembl VEP's plugin
documents, loads under the default `study_types` with that family's default
columns. Ensembl VEP 116.2's OpenTargets plugin selects both families
whatever the header carries and asks for both families' default columns
(`OpenTargets.pm`, lines 202-205 and 273-278 of the `VEP_plugins` release/116
branch). It drops the plugin, writing
`Failed to instantiate plugin OpenTargets: ERROR: The following columns were not found in file header`
to its warnings file (lines 227-228) and every row without an `OpenTargets_`
field; there, `study_types=all` or the family's own name loads it.

`lead_variants_only` takes `1` or `0`, and `0` switches the filter off in both
engines. vep-rs drops the plugin at load on any other value, with the warning
`lead_variants_only must be 0 or 1, not '<value>'`. Ensembl VEP's plugin reads
the value by Perl truth (`OpenTargets.pm`, lines 262-267 and 352).
`lead_variants_only=false` is therefore a true string there and keeps the filter
on without a warning. Neither shape is on a committed dataset (tier 3 of
[Intended divergences from Ensembl VEP](../intended-divergences.md); 4 flag sets
run the plugin).

**Preparing the data.** The Platform publishes `open_targets_vep_full.tsv.gz`
(every association), `open_targets_vep_gwas.tsv.gz` and
`open_targets_vep_qtl.tsv.gz`, each with its `.tbi`, under
`etc/export/vep_data/` of a release directory. The published header line is not
commented, which Ensembl VEP's OpenTargets plugin cannot read; comment it and
index the file again, as that plugin's own instructions say:

```bash
zcat open_targets_vep_full.tsv.gz | sed '1 s/.*/#&/' \
  | bgzip -c > open_targets_vep.tsv.gz
tabix -s 1 -b 2 -e 2 open_targets_vep.tsv.gz
```

vep-rs reads the header with or without the `#`.

## SpliceVault {#splicevault}

SpliceVault adds the most common mis-splicing events around a splice site a
variant is predicted to disrupt, from the SpliceVault table the user passes
it, as Ensembl VEP's `SpliceVault.pm` does.

| Property      | Value                                                                                                                      |
| ------------- | -------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                             |
| Matches       | transcript rows only: position and allele, the record's alleles and the row's each trimmed of the bases they share at either end, as Ensembl VEP trims them; under `transcript_match=1` the table row must also name the row's transcript |
| Output fields | `SpliceVault_top_events`, `SpliceVault_site_sample_count`, `SpliceVault_site_max_depth`, `SpliceVault_out_of_frame_events`, `SpliceVault_site_pos`, `SpliceVault_site_type`, `SpliceVault_SpliceAI_delta` |
| Data source   | `SpliceVault_data_GRCh38.tsv.gz` and `SpliceVault_data_hg19.tsv.gz` |

Every parameter is `key=value`; every key but `transcript_match` names a table,
read in the order written (`SpliceVault.pm` lines 119-122).

| Parameter                  | Meaning                                                                                                                              |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `file=<path>`              | A table; an empty or `0` path is refused.                                                                                            |
| `transcript_match=<value>`  | `1`, the default, annotates a transcript row from the table rows naming its stable id; `0` drops that test, so every transcript row of the record takes the first table row at the site whose allele matches, whatever transcript it names. Every value but an empty one and `0` keeps the test, as Ensembl VEP's truth test reads it (`SpliceVault.pm` line 161). |

```bash
--plugin SpliceVault,file=/path/to/SpliceVault_data_GRCh38.tsv.gz
--plugin SpliceVault,file=/path/to/SpliceVault_data_GRCh38.tsv.gz,\
transcript_match=0
```

The plugin reads a table row's twelve cells in this order (`SpliceVault.pm` line
201):

| Column | Content |
| --- | --- |
| 1 | chromosome |
| 2 | position |
| 3 | reference allele |
| 4 | alternate allele |
| 5 | transcript |
| 6 | site type |
| 7 | site position |
| 8 | SpliceAI delta |
| 9 | out-of-frame fraction |
| 10 | top events |
| 11 | sample count |
| 12 | maximum depth |

The table rows from two bases before the record's start to its end are read; the
first whose alleles match the record's is the row's. `SpliceVault_top_events`
lists the events `rank:type:impact:percent:frame`, joined by `,` in the default
and tab formats and `&` in the VCF. In JSON the seven fields nest under
`splicevault`, the events keyed by rank.

A run under `--refseq` or `--merged` with `transcript_match=1` is warned on
standard error, under `--quiet` too and never in the warnings file, that the
match may fail, as Ensembl VEP warns (`SpliceVault.pm` lines 129-131).

## TranscriptAnnotator {#transcriptannotator}

TranscriptAnnotator appends the columns of a table of your own to the
transcript rows the table names, as Ensembl VEP's `TranscriptAnnotator.pm`
does.

| Property      | Value                                                                                                                             |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                                    |
| Matches       | transcript rows only: position and allele, and the row's transcript id equal to the transcript's stable id or one of its RefSeq cross-references |
| Output fields | one `<prefix>_<column>` field per selected column                                                                                  |
| Data source   | your own bgzipped, tabix-indexed table                                                                                            |

| Parameter       | Meaning                                                                                                                             |
| --------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `file=<path>`   | The table (mandatory): chromosome, position, reference allele, alternate allele and transcript id, then the columns to append, named by the table's last header line. A header line without a leading `#` (a table indexed with `tabix -S 1`) names the columns too; Ensembl VEP 116.2 reads the header through `tabix -h`, which prints `#` lines alone, so such a table gives it no column names and no output fields (`TranscriptAnnotator.pm` lines 99-105). |
| `cols=<a:b>`    | The columns written, `:`-separated; every column after the fifth by default. A name is matched to the header names exactly; Ensembl VEP 116.2 tests each name as a pattern and keeps it as typed, so `SIFT_.core` heads an empty column there ([a documented class](../divergences/plugin-cols-names-matched-exactly.md)). |
| `prefix=<text>` | The field prefix; the file's name without its extensions by default, and `prefix=0` writes the column names bare.                   |
| `trim=<0 or 1>` | `1`, the default, removes the whitespace around every header name and cell before it is read; `0` keeps the cells as written. |

```bash
--plugin TranscriptAnnotator,file=/path/to/annotations.tsv.gz
--plugin TranscriptAnnotator,file=/path/to/annotations.tsv.gz,\
cols=SIFT_score:Comment,prefix=my
```

`trim` and `cols` are read by their truth, as Ensembl VEP 116.2 reads them
(`TranscriptAnnotator.pm` lines 132 and 164). Only `0` and an empty value are
false. `trim=yes` trims, `trim=` does not, and `cols=0` writes every column.

An empty cell is an empty value and a missing column an absent one
([data files](data-files.md#transcriptannotator-cells)).

The rows from two bases before the record's span to its end are read. The span
runs from its lower to its higher coordinate, so that an insertion's start lies
inside it. The first row whose transcript and alleles match the transcript row's
is written on it. Ensembl VEP 116.2 reads from two bases before the record's
start to its end (`TranscriptAnnotator.pm` line 218), a two-base window for an
insertion. 996,643 of the chromosome 21 population cuts' records have that shape
([the insertion-window class, section
10](#10-transcriptannotator-a-table-row-written-at-an-insertions-start)).

## gnomADMt {#gnomadmt}

gnomADMt adds gnomAD's mitochondrial allele frequencies and haplogroup
statistics to every row of a record on `MT`, from the gnomAD v3.1 chrM sites
VCF, as `gnomADMt.pm` does.

| Property      | Value                                                                                                                                 |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh38                                                                                                                                |
| Matches       | position and allele on `MT`; the values print on every row of the record, intergenic rows included                                   |
| Output fields | `gnomAD_<field>` for each selected INFO field of the file: by default the five haplogroup fields `gnomAD_hap_AC_hom`, `gnomAD_hap_AC_het`, `gnomAD_hap_AF_hom`, `gnomAD_hap_AF_het` and `gnomAD_hap_AN`; every field under `fields=all` |
| Data source   | `gnomad.genomes.v3.1.sites.chrM.vcf.bgz` with its index, from the gnomAD downloads                                                     |

| Parameter             | Meaning                                                                                                                        |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `file=<path>`         | The sites VCF (mandatory); a `##fileformat` not naming major version 4 is refused. Its contig may be named `MT`, `M`, `chrM` or `chrMT`; Ensembl VEP 116.2 reads a `chrM`- or `M`-named one only ([A gnomADMt data file whose contig is named `MT`](../divergences/gnomadmt-mt-named-contig.md)). |
| `fields=<a:b>`        | The INFO fields written, `:`-separated by their names in the file, or `all`.                                                   |
| `hap_filter=<a:b>`    | The haplogroups kept, in the order named, in every per-haplogroup field; `all` keeps every one in the file's order. On a file whose `hap_AN` description states no order list the filter is refused at plugin initialization, `all` included; Ensembl VEP 116.2 refuses a named filter but loads `all` with a warning and writes every per-haplogroup field empty (`gnomADMt.pm` lines 82, 114 and 347). The gnomAD chrM sites file states the list. |
| `pop_filter=<a:b>`    | The populations kept, in the order named, in every per-population field; `all` keeps every one in the file's order. A file whose `pop_AN` description states no order list is read as `hap_filter` reads `hap_AN`. |

```bash
--plugin gnomADMt,file=/path/to/gnomad.genomes.v3.1.sites.chrM.vcf.bgz
--plugin gnomADMt,file=/path/to/gnomad.genomes.v3.1.sites.chrM.vcf.bgz,\
fields=all
```

Each field's header line carries its `##INFO` description. A per-haplogroup or
per-population field is a `|`-joined list in the order the file's `hap_AN` or
`pop_AN` description states. The `hap_hl_hist` and `pop_hl_hist` histograms are
a `,`-joined list in that order. A filter keeps the entries named, or every
entry in file order under `all`, and rewrites the description's order list to
the entries kept.

## LOEUF {#loeuf}

LOEUF adds gnomAD's loss-of-function observed/expected upper bound fraction of
a transcript, or the minimum over its gene's transcripts, as Ensembl VEP's
`LOEUF.pm` does.

| Property      | Value                                                                                                                                  |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37 (the gnomAD v2.1.1 table by transcript), GRCh38 (the gnomAD v4 table with gene coordinates appended)                              |
| Matches       | transcript rows of sequence variants, never an `upstream_gene_variant` or `downstream_gene_variant` row: the table rows whose gene span covers the variant, then the row naming the transcript, or every row naming the gene |
| Output fields | `LOEUF`                                                                                                                                |
| Data source   | the gnomAD constraint tables, tabix-indexed on their gene coordinate columns as Ensembl VEP's plugin documents; the prepared GRCh38 table is on the Ensembl FTP site under `pub/current/variation/LOEUF/` |

Both parameters are required.

| Parameter          | Meaning                                                                                                              |
| ------------------ | -------------------------------------------------------------------------------------------------------------------- |
| `file=<path>`      | The table.                                                                                                           |
| `match_by=<mode>`  | `transcript`: the value of the table row naming the row's transcript, as written. `gene`: the smallest numeric value over the rows naming the row's gene, as the table wrote it. Another word is refused, where Ensembl VEP loads the plugin and writes no value ([a documented divergence](../divergences/plugin-unknown-mode-word-refused.md)). |

```bash
--plugin LOEUF,file=/path/to/loeuf_dataset_grch38.tsv.gz,match_by=gene
--plugin LOEUF,file=/path/to/loeuf_dataset_grch38.tsv.gz,match_by=transcript
```

The two gnomAD layouts are told apart by the column names at fixed positions
(`transcript`, `oe_lof_upper` and `gene_id` in the v2 table; `transcript`,
`lof.oe_ci.upper` and `gene_id` in the v4 table). A table with other names there
is refused with the missing names listed. The coordinate columns are the ones
the tabix index declares.

## Where these plugins differ from Ensembl VEP

Eight behaviours of Ensembl VEP 116.2's plugins are not reproduced, each a
walkthrough of one record on
[Where the annotation-store plugins differ from Ensembl VEP](annotation-store-divergences.md):

1. [AlphaMissense: a two-digit transcript version under `transcript_match=1`](annotation-store-divergences.md#1-alphamissense-a-two-digit-transcript-version-under-transcript_match1)
2. [CADD: structural alleles outside insertion, deletion and duplication](annotation-store-divergences.md#2-cadd-structural-alleles-outside-insertion-deletion-and-duplication)
3. [SpliceVault: the table row of a sibling allele](annotation-store-divergences.md#3-splicevault-the-table-row-of-a-sibling-allele)
4. [SpliceVault: the JSON object of a record's later transcript rows](annotation-store-divergences.md#4-splicevault-the-json-object-of-a-records-later-transcript-rows)
5. [TranscriptAnnotator: the table row of a sibling allele](annotation-store-divergences.md#5-transcriptannotator-the-table-row-of-a-sibling-allele)
6. [gnomADMt: the `##INFO` descriptions of the file's fields](annotation-store-divergences.md#6-gnomadmt-the-info-descriptions-of-the-files-fields)
7. [LOEUF: a gene's minimum of `0`](annotation-store-divergences.md#7-loeuf-a-genes-minimum-of-0)
8. [OpenTargets: a cell of `0`](annotation-store-divergences.md#8-opentargets-a-cell-of-0)

The comparison tools set six aside; the gnomADMt and OpenTargets ones are
documented only. Eight more follow, three below, five on their own pages:

- [the `MT`-named gnomADMt file](../divergences/gnomadmt-mt-named-contig.md)
- [the TranscriptAnnotator record at position 1 of contig
  `1`](../divergences/transcript-annotator-header-from-record.md)
- [AlphaMissense's prediction of a sibling
  allele](../divergences/alphamissense-sibling-allele-rows.md)
- [REVEL's reference amino
  acid](../divergences/revel-other-reference-amino-acid-rows.md)
- [OpenTargets' repeated `QTL` and export without
  `isLead`](../divergences/opentargets-duplicate-study-type-and-missing-islead.md)

Each names its rows; the format comparator sets aside the three below,
AlphaMissense's and REVEL's. [Intended divergences from Ensembl
VEP](../intended-divergences.md) indexes all sixteen.

Each of the three carries its tier under that page's rule. Tier 1 is a shape the
chromosome 21 population cuts carry with no flag. Tier 2 is one that needs a
flag or plugin on ordinary data. Tier 3 is one ordinary data does not carry.

The counts are over the cuts' 15,646,481 records, the 111 golden corpora and
their 392 flag sets. Every row below is in the default output format.

### 9. REVEL: the score of a sibling allele on a multi-allelic record

**The defect in one sentence.** Ensembl VEP's REVEL plugin matches the score
table against every alternate allele of a record and returns the first row
whose amino acid equals the row's own, so a multi-allelic record's allele takes
its sibling's score when both give the same amino acid.

**Why this is a defect.** Each table row scores one substitution: a position, a
reference base, an alternate base, the amino acid change and a transcript. The
plugin tests the row's own alternate amino acid and transcript, so its contract
is this allele's score on this transcript. On a one-allele record it reads
exactly that row. Only the record's allele list widens the match to a sibling's
row.

**What vep-rs does instead.** vep-rs matches each allele's rows against that
allele alone, so every row carries its own allele's score or none.

**Where it is in Ensembl VEP's source** (`REVEL.pm`, `VEP_plugins` at
release/116):

- `REVEL.pm`, line 165, `run`: the variant side of the match is `alt_alleles`,
  every alternate allele of the record, not the allele being annotated.
- `REVEL.pm`, lines 176-199: each table row at the position is matched against
  that list. The first row whose alternate amino acid equals the row's peptide
  is returned: line 193 with a transcript match, line 197 under `no_match`.

**What correcting it would change.** Each allele's rows would carry the table's
score for that allele, as they do on one-allele records.

**The class it produces.** `revel_sibling_allele_rows`: a `REVEL` field (the
JSON `revel` member) where Ensembl VEP's value is the one its row of a sibling
allele carries on that feature beside the same amino acids, vep-rs's being empty
or different. The format comparator sets it aside where a comparison case names
the class; any other difference stays one. Each REVEL golden corpus holds one
two-allele record, `G A,T`, whose alleles give different amino acids. The class
therefore applies to 0 rows of their 24 compared outputs (6 flag sets, four
formats) and of the four comparison cases that run them.

**How often the shape occurs.** Tier 2: 6,347 of the population cuts' records
(406 per million) carry two or more sequence alleles. The codon decides whether
two give one amino acid. 64 corpora carry the shape and 6 flag sets run the
plugin.

#### 9.1 Two alleles of one codon that both give leucine

**In short.** The record `21 13610388 A C,T` on the forward-strand transcript
`ENST00000299443`, against the GRCh38 REVEL golden corpus's table, whose rows at
the position read `A C M L 0.025`, `A G M V 0.043` and `A T M L 0.027`.

Both engines write `REVEL=0.025` on the `C` row. On the `T` row Ensembl VEP
writes `REVEL=0.025`, the `C` row's score. The `A C` table row matches the
record's `C` allele, its amino acid `L` equals the `T` allele's, and it comes
first. vep-rs writes `REVEL=0.027`, the table's row for `A T`. Run as two
one-allele records, Ensembl VEP writes vep-rs's rows.

### 10. TranscriptAnnotator: a table row written at an insertion's start

**The defect in one sentence.** Ensembl VEP's TranscriptAnnotator plugin reads
the table rows from two bases before the variant's start to its end, and an
insertion's end precedes its start, so the window is two bases wide and never
reaches a row written at the insertion's start.

**Why this is a defect.** Ensembl VEP numbers an insertion between two bases as
start = end + 1. A window from `start - 2` to `end` therefore covers the two
bases before the insertion point and no more. A row written at the insertion's
own start, the coordinate Ensembl VEP prints for the allele, lies outside it and
is never offered to the allele match.

Three other plugins of the same collection read the corrected span. gnomADMt and
PromoterAI set the end to the start when the start is past it. MaveDB adds one
to the end "to account for insertions".

**What vep-rs does instead.** vep-rs reads from two bases before the lower of
the two coordinates to the higher. A row at the insertion's start is therefore a
candidate and is written when its alleles match.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `TranscriptAnnotator.pm`, line 218, `run`:
  `get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})`, the rows with a position
  in that closed span.
- `gnomADMt.pm`, line 265, and `PromoterAI.pm`, line 232:
  `$vf_end = $vf_start if $vf_start > $vf_end;` before the same read;
  `MaveDB.pm`, line 283: the read ends at `$vf->{end} + 1`.

**What correcting it would change.** An insertion would take a table row
written at its start, as a substitution and a deletion take theirs.

**The class it produces.** `transcript_annotator_insertion_window`, applied from
a corpus manifest: the fields of the table's prefix on an insertion's row that
vep-rs alone fills, from a row at the insertion's start or two bases before its
end. The `GRCh38-plugin-transcript-annotator` corpus holds no insertion, and the
one insertion of the `GRCh37-plugin-transcript-annotator` corpus has its row at
the record's own position, inside both windows. No row of their five flag sets
is therefore set aside. No comparison case loads the plugin.

**How often the shape occurs.** Tier 3: the shape is a table row written at an
insertion's start with reference `-`. None of the 18 rows of the two committed
tables has one, and 5 flag sets run the plugin. Pure insertions themselves are
common, 996,643 of the population cuts' records (63,698 per million). A table of
the user's own therefore reaches the shape as soon as it names one that way.

#### 10.1 A one-base insertion against a row written at its start

**In short.** For the insertion `21 10541124 insA T TA` on GRCh38, both engines
print the allele at `21:10541124-10541125`. A table row
`21 10541125 - A ENST00000618007 0.42 Tolerated VEPstyleIns` names that start.
vep-rs writes `ins_vepstyle_SIFT_score=0.42`, `ins_vepstyle_SIFT_pred=Tolerated`
and `ins_vepstyle_Comment=VEPstyleIns` on the `ENST00000618007` row; Ensembl
VEP 116.2 writes nothing on it, its window being 10541123-10541124. The same row
written as `21 10541124 T TA ENST00000618007 0.43 Tolerated VCFstyleIns` is
inside both windows, and both engines write its three fields.

### 11. gnomADMt: an INFO value that holds `=` is cut at that sign

**The defect in one sentence.** Ensembl VEP's gnomADMt plugin reads each
`key=value` item of a file record's INFO column with a two-variable `split`,
which Perl runs with a limit one more than the number of variables, so a value
holding an `=` of its own is cut there.

**Why this is a defect.** Perl's `split` documentation states the implicit
limit: a list assignment of two variables splits into at most three parts and
discards the third. The VCF specification lets an INFO value carry `=`, and the
gnomAD mitochondrial file does: in its `vep` field the HGVSp of every synonymous
change ends in `=`. The plugin's documentation offers every INFO field by name
and all of them under `fields=all`. Yet the `vep` value it prints for such a
record stops at that sign, its remaining 33 columns gone.

**What vep-rs does instead.** vep-rs splits each INFO item at its first `=`
only and prints the file's whole value.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `gnomADMt.pm`, lines 327-330, `_parse_vcf`:
  `my ($key, $value) = split /=/, $field` for each `;`-separated item, the value
  stored under the key; line 373, `parse_data`, and line 293, `run`, hand that
  value to the row.

**What correcting it would change.** Every selected field would print the
file's value, as vep-rs prints it; the five default haplogroup fields carry no
`=` and are unchanged.

**The class it produces.** `gnomadmt_info_value_second_equals`: a
`gnomAD_<field>` (JSON `gnomad_<field>`) both engines fill, where the
reference's value is the candidate's cut before the candidate's first `=`.

No row on the committed corpora or their comparison cases. In the committed chrM
slice 47 of 373 records carry such a value, all in `vep`. None of the 11 records
the corpus's variants match does, so the `fields=all` flag set has 0 rows in its
JSON format. The four other flag sets select no field that carries one, and the
population cases select the default fields.

**How often the shape occurs.** Tier 2: the shape is a file record whose INFO
item carries a second `=` inside its value. 47 of the 373 records of the
committed gnomAD mitochondrial slice do (126,005 per million), every one in the
`vep` field. 5 flag sets run the plugin.

#### 11.1 chrM:3441 A>G under `fields=vep`

**In short.** The file's record `chrM 3441 A G` carries a `vep` value of 362
characters in 45 columns,
`G|synonymous_variant|LOW|MT-ND1|ENSG00000198888|Transcript|ENST00000361390|protein_coding|1/1||ENST00000361390.2:c.135A>G|ENSP00000354687.2:p.Leu45=|135|135|45|L|ctA/ctG|...`.
On the input `MT 3441 A G` under `fields=vep`, vep-rs writes `gnomAD_vep=`
followed by all 362 characters. Ensembl VEP's plugin, handed the file's line,
stores the 147 characters before the `=` of `p.Leu45=`. Its rows therefore end
the value at `p.Leu45` and carry no CDS position, codons, HGNC id or protein
domains.

## Where next

- [Preparing data files](data-files.md): how to bgzip, index and convert the
  data files.
- [SpliceAI](spliceai.md) reads its data likewise and has its own page.
- [Where the annotation-store plugins differ from Ensembl
  VEP](annotation-store-divergences.md): a walkthrough of one record for each of
  the eight classes.
- [Intended divergences from Ensembl VEP](../intended-divergences.md): the index
  of every class.
