# Annotation-store plugins

An annotation-store plugin reads a bgzipped, tabix-indexed data file, or the
`.vpd` binary store built from it, and matches records by position and then by
whatever its Matches row adds. This page documents CADD, REVEL, gnomADc,
AlphaMissense and dbscSNV. [SpliceAI](spliceai.md) reads its data the same way
and has its own page because its open data covers GRCh38 SNVs only. How to
bgzip, index and convert the files is on [Preparing data files](data-files.md).

Pass `--assembly` whenever one of these plugins is enabled: CADD, REVEL, dbscSNV
and AlphaMissense use it to refuse a data file built on the other assembly, and
REVEL uses it to pick the coordinate column its index was built on.

## CADD {#cadd}

CADD adds Combined Annotation Dependent Depletion scores to a variant.

| Property      | Value                                                     |
| ------------- | --------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                            |
| Matches       | position and allele                                       |
| Output fields | `CADD_PHRED`, `CADD_RAW`                                  |
| Data source   | [CADD download](https://cadd.gs.washington.edu/download) |

| Parameter          | Meaning                                                                                                          |
| ------------------ | ---------------------------------------------------------------------------------------------------------------- |
| `snv=<path>`       | The SNV score file.                                                                                              |
| `indels=<path>`    | The indel score file.                                                                                            |
| `force_annotate=1` | Let a later file overwrite the scores an earlier file set. Off by default, so the first file that matches wins. |

```bash
--plugin CADD,snv=/path/to/snv.tsv.gz,indels=/path/to/indels.tsv.gz
--plugin CADD,snv=/path/to/snv.tsv.gz,indels=/path/to/indels.tsv.gz,\
force_annotate=1
```

**Preparing the data.** Download from
[CADD download](https://cadd.gs.washington.edu/download):

- `whole_genome_SNVs.tsv.gz` + `.tbi` (SNVs)
- the release's indel file + `.tbi`. In v1.7 it is the first name below on
  GRCh38 and the second on GRCh37:

  ```text
  gnomad.genomes.r4.0.indel.tsv.gz
  gnomad.genomes-exomes.r4.0.indel.tsv.gz
  ```

  v1.6 GRCh37 and earlier releases name it `InDels.tsv.gz`.

A file whose first line names the other assembly is refused under `--assembly`.

## REVEL {#revel}

REVEL adds the Rare Exome Variant Ensemble Learner missense pathogenicity score
to a transcript consequence.

| Property      | Value                                                                        |
| ------------- | ---------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                               |
| Matches       | position and allele, then transcript unless `no_match=1`; missense only      |
| Output fields | `REVEL_score`                                                                |
| Data source   | [REVEL download](https://sites.google.com/site/revelgenomics/downloads)     |

| Parameter     | Meaning                                                                |
| ------------- | ---------------------------------------------------------------------- |
| `file=<path>` | The REVEL score file.                                                  |
| `no_match=1`  | Report the score even when the record's transcript is not this one.   |

```bash
--plugin REVEL,file=/path/to/revel.tsv.gz
--plugin REVEL,file=/path/to/revel.tsv.gz,no_match=1
```

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

One download serves both assemblies: the GRCh37 index is built on column 2 and
the GRCh38 index on column 3. `--assembly` selects the column; without it REVEL
falls back to the GRCh37 column with only a warning.

## gnomADc {#gnomadc}

gnomADc adds gnomAD coverage statistics to a variant.

| Property      | Value                                                                                                   |
| ------------- | ------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                          |
| Matches       | position only                                                                                           |
| Output fields | `gnomADc_<column>` for every header column after chromosome and position (`gnomADc_mean`, `gnomADc_over_1`, ...) |
| Data source   | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads)                                        |

| Parameter | Meaning                                  |
| --------- | ---------------------------------------- |
| `<path>`  | The coverage file, as the only parameter. |

```bash
--plugin gnomADc,/path/to/gnomad.coverage.tsv.gz
```

The field names follow the file's own header, so the GRCh37 v2.1 coverage file
yields `gnomADc_median` and the GRCh38 v4 file `gnomADc_median_approx` and
`gnomADc_total_DP`. These names are vep-rs's own: the Perl gnomADc plugin
writes `gnomAD_mean_cov`, `gnomAD_median_cov` and `gnomAD_<N>x_cov` instead.

**Preparing the data.** Download coverage files from
[gnomAD downloads](https://gnomad.broadinstitute.org/downloads). No tabix index
is published, so build one on the chromosome and position columns (`bgzip`,
then `tabix -s 1 -b 2 -e 2`). The GRCh38 v3 and v4 files carry a single `locus`
column (`chr1:11819`); split it into chromosome and position first
(`sed "1s/locus/chr\tpos/; s/:/\t/g"`). The header row works with or without a
leading `#`.

## AlphaMissense {#alphamissense}

AlphaMissense adds DeepMind's AlphaMissense pathogenicity prediction to a
variant.

| Property      | Value                                                       |
| ------------- | ----------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                              |
| Matches       | position and allele                                         |
| Output fields | `am_class`, `am_pathogenicity`                              |
| Data source   | [AlphaMissense data](https://zenodo.org/records/8208688)   |

| Parameter | Meaning                                     |
| --------- | ------------------------------------------- |
| `<path>`  | The prediction file, as the only parameter. |

```bash
--plugin AlphaMissense,/path/to/AlphaMissense_hg38.tsv.gz
```

**Preparing the data.** AlphaMissense publishes one file per build, and every
data row names its build in the `genome` column. A file whose rows name the
other build is refused under `--assembly`.

## dbscSNV {#dbscsnv}

dbscSNV adds AdaBoost and Random Forest splice-site predictions to a transcript
consequence.

| Property      | Value                                                                                                                                                                                |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Builds        | GRCh37, GRCh38                                                                                                                                                                       |
| Matches       | position and allele, per transcript                                                                                                                                                  |
| Output fields | `dbscSNV_ADA_SCORE`, `dbscSNV_RF_SCORE`                                                                                                                                              |
| Data source   | [dbscSNV](https://sites.google.com/site/jpopgen/dbNSFP); the plugin reads the 6-column `dbscSNV1.1_<assembly>.txt.gz` redistributed by [CADD](https://kircherlab.bihealth.org/download/CADD/v1.7/) |

| Parameter | Meaning                                    |
| --------- | ------------------------------------------ |
| `<path>`  | The 6-column file, as the only parameter.  |

```bash
--plugin dbscSNV,/path/to/dbscSNV1.1_GRCh38.txt.gz
```

**Preparing the data.** Use the 6-column `dbscSNV1.1_GRCh37.txt.gz` or
`dbscSNV1.1_GRCh38.txt.gz` (+ `.tbi`) redistributed under
[CADD v1.7](https://kircherlab.bihealth.org/download/CADD/v1.7/) (`GRCh37/`,
`GRCh38/`). The plugin reads chromosome, position, ref and alt from the first
four columns and `ada_score` and `rf_score` by header name. The header may lack
its leading `#` and the lines may end in CRLF; vep-rs accepts both. The GRCh38
file is recognised by its `hg38_chr` first column.

The preparation the Perl dbscSNV plugin documents for the full dbscSNV1.1 zip
keeps every column with `chr` first and, on GRCh38, indexes columns 5 and 6.
The plugin does not read those columns, and under `--assembly GRCh38` it refuses
such a file as a GRCh37 one.
