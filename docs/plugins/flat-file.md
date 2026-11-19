# Flat-file plugins

LoFtool, pLI, DosageSensitivity and GWAS add gene intolerance scores, dosage
sensitivity scores and GWAS Catalog associations to a variant from flat files.

A flat-file plugin loads its whole data file at init and looks records up in
memory, by gene symbol or by position. LoFtool, pLI and DosageSensitivity are
keyed by gene symbol, so one data file serves both assemblies. GWAS is GRCh38-only because
the NHGRI-EBI catalog publishes `CHR_POS` in GRCh38 coordinates only.

## LoFtool {#loftool}

LoFtool adds a gene's intolerance to functional variation to a transcript
consequence, as Ensembl VEP's plugin does.

| Property      | Value                                                                                                                                                       |
| ------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                                                              |
| Matches       | gene symbol lookup, case-insensitive                                                                                                                        |
| Output fields | `LoFtool`                                                                                                                                                   |
| Data source   | `LoFtool_scores.txt` from the [Ensembl VEP_plugins](https://github.com/Ensembl/VEP_plugins) repository; a copy is committed under `tests/fixtures/plugins/` |

| Parameter | Meaning                                                                                                        |
| --------- | -------------------------------------------------------------------------------------------------------------- |
| `<path>`  | The scores file. Without it, or with a value Perl reads as false (`0` or the empty string), `LoFtool_scores.txt` is read from the first `--dir_plugins` directory holding it. |

```bash
--plugin LoFtool,/path/to/LoFtool_scores.txt
--dir_plugins /path/to/VEP_plugins --plugin LoFtool
```

The file has two whitespace-separated columns, a gene symbol and its
percentile, under a header line. A score prints as Perl's `%g` writes it, so
the file's `6.89E-05` reaches the output as `6.89e-05` and its `1.38E-04` as
`0.000138`.

A cell is read as `sprintf` numifies it (`LoFtool.pm` line 89): its
leading number, so `1.5abc` scores `1.5`. A gene listed twice takes its
later line. A cell with no leading number (`NA`) scores nothing, warned once at
load. Ensembl VEP 116.2 formats it as `0`, which its truth test drops, so
neither engine prints it.

**Preparing the data.** Fetch `LoFtool_scores.txt` from the VEP_plugins
repository and pass its path, or point `--dir_plugins` at a checkout, where the
file sits beside the plugin module.

## pLI {#pli}

pLI adds a gene's probability of loss-of-function intolerance to a transcript
consequence, as Ensembl VEP's pLI plugin does.

| Property      | Value                                                                                                                                                                        |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                                                                               |
| Matches       | gene symbol lookup, or transcript stable id lookup in `transcript` mode; both case-insensitive                                                                               |
| Output fields | `pLI_gene_value`, or `pLI_transcript_value` in `transcript` mode                                                                                                             |
| Data source   | `pLI_values.txt` from the [Ensembl VEP_plugins](https://github.com/Ensembl/VEP_plugins) repository (ExAC r0.3 gene pLI); a copy is committed under `tests/fixtures/plugins/` |

| Parameter | Meaning                                                                                                                    |
| --------- | -------------------------------------------------------------------------------------------------------------------------- |
| `<path>`  | The values file, first. Without it, or with a value Perl reads as false (`0` or the empty string), `pLI_values.txt` is read from the first `--dir_plugins` directory holding it. |
| `<mode>`  | `gene` (the default) or `transcript`, second: whether the file's first column holds gene symbols or transcript stable ids. A later positional parameter is not read, as Ensembl VEP's plugin reads it, and is noted in the warnings file. |
| `sigfig`  | Anywhere in the list: print a score with three significant figures (`%.3g`) instead of two decimals (`%.2f`).              |

```bash
--plugin pLI,/path/to/pLI_values.txt
--plugin pLI,/path/to/pLI_transcript_values.txt,transcript,sigfig
--dir_plugins /path/to/VEP_plugins --plugin pLI
```

The file has two whitespace-separated columns under a header, the last line
containing `gene` or `transcript` in lower case (`pLI.pm` lines 125 to 131). A
header without the mode's column is refused, as is a mode other than the two
([a documented divergence](../divergences/plugin-unknown-mode-word-refused.md)).
A transcript file keys its rows by
the stable id without a version, the form the lookup uses. The extract the
plugin's documentation describes matches once the version suffix is dropped.

In `transcript` mode under `--refseq` or `--merged` the plugin warns that
RefSeq ids never match an Ensembl stable id. The warning is printed on standard
error, under `--quiet` too, and never written to the warnings file, as Ensembl
VEP prints it. Ensembl VEP 116.2 tests the second parameter as given, so
`pLI,<path>,sigfig,transcript` warns here and not there (`pLI.pm` lines 109
and 187).

A cell is read as `sprintf` numifies it (`pLI.pm` line 145): its leading
number, so `1.5abc` scores `1.50`. A key listed twice takes its later line.
A cell with no leading number (`NA`) scores nothing, warned once at load.
Ensembl VEP 116.2 prints the `0.00` its `sprintf` makes of it on every gene row.
No distributed values file carries one, and the comparison tools set nothing
aside for it.

## GWAS {#gwas}

GWAS adds an NHGRI-EBI GWAS Catalog association to every row of a variant, from
the catalog's associations download or one study's harmonised summary-statistics
file.

| Property      | Value                                                                                                                                                                                                   |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh38                                                                                                                                                                                                  |
| Matches       | the site, within two bases ahead of the record's start, and the risk allele against the record's alternate alleles, each pair of reference and allele trimmed to its minimal form as Ensembl VEP trims them. Under `--minimal` a split allele is matched by its own trimmed alleles over its own span; Ensembl VEP matches the rejoined record's whole allele list over its span ([a documented class](../divergences/minimal-plugin-sees-the-rejoined-record.md)) |
| Output fields | `GWAS_accessions`, `GWAS_associated_gene`, `GWAS_beta_coef`, `GWAS_odds_ratio`, `GWAS_p_value`, `GWAS_pmid`, `GWAS_risk_allele`, `GWAS_study`, the field set of Ensembl VEP's GWAS plugin                  |
| Data source   | [GWAS Catalog downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads) and the catalog's [summary statistics](https://www.ebi.ac.uk/gwas/docs/summary-statistics)                                      |

Every parameter is `key=value`; a bare path is refused.

| Parameter       | Meaning                                                                                                                                      |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `file=<path>`   | The data file. Required. A catalog that is empty or lacks `CHR_ID`, `CHR_POS`, `SNPS`, `DISEASE/TRAIT` or `STRONGEST SNP-RISK ALLELE` is refused at plugin initialization, where Ensembl VEP 116.2 loads it and annotates nothing. A row whose SNP count differs from its location count puts every SNP at its one location, and is skipped when it names several; Ensembl VEP writes every SNP at every location of it, the shape `gwas_haplotype_cross_product` sets aside (the associations download holds one such row, `GCST002273`, three SNPs over two locations on chromosome 8). A chromosome written with a `chr` prefix or as `M` is read as the engine names it (`21`, `MT`); Ensembl VEP's plugin matches a `chr` prefix only, and the catalog writes neither. |
| `type=<kind>`   | `curated` (the default), the associations download, or `sstate`, a harmonised summary-statistics file named `<PubMed id>-<study accession>-<trait accession>.<suffix>`, the catalog's `.h.tsv.gz` or the uncompressed `.h.tsv`. A file whose name lacks one of the three parts or the suffix is refused at plugin initialization, where Ensembl VEP 116.2 loads it and writes `GWAS_pmid`, `GWAS_study` and `GWAS_accessions` empty (`GWAS.pm` lines 139-142); a harmonised file without `hm_chrom` and `hm_pos` is refused by both engines. The file is read by its tab-separated columns throughout; Ensembl VEP 116.2's filter and sort count blank-separated fields, so a row with an empty cell before `hm_pos` can fail its index build and drop the plugin ([a documented class](../divergences/gwas-sstate-index-failure.md)). |
| `verbose=1`     | Name every catalog row the plugin skips: one without a trait or without an `rs` id. A risk allele of `?` is read as none, with no message. |
| `db=1`          | Refused: Ensembl VEP's plugin reads each SNP's locations from the Ensembl database under it, and vep-rs has no database connection.          |

```bash
--plugin GWAS,file=/path/to/gwas_catalog_associations.tsv
--plugin GWAS,type=sstate,file=/path/to/17463246-GCST000028-EFO_0001360.h.tsv
```

**Preparing the data.** Download the associations TSV, or one study's
harmonised summary-statistics file, from the catalog. The plugin reads either
plain or gzipped, told apart by content, and locates its columns by header
name, so no reformatting is needed. Ensembl VEP 116.2 tells a compressed catalog
by its `gz` suffix (`GWAS.pm` line 366) and reads a summary-statistics file
uncompressed only (`GWAS.pm` lines 437 and 460). A summary-statistics file
keeps the name the catalog gives it,
`<PubMed id>-<study accession>-<trait accession>.h.tsv.gz`: the three fields
reach the output as `GWAS_pmid`, `GWAS_study` and `GWAS_accessions`.

**What the fields carry.** The catalog's text reaches the output as Ensembl
VEP's plugin writes it. The reported genes lose their whitespace and the trait
URIs stay as the catalog lists them. `OR or BETA` becomes `GWAS_beta_coef`,
with the unit of the confidence-interval column, when that column names an
increase or a decrease, and `GWAS_odds_ratio` otherwise. A field the row lacks
is the empty string: an empty cell in the tab format and `""` in JSON.

**Columns the header may lack.** A header lacking any column but the five the
`file=` row names loads with that field empty (`GWAS.pm` lines 237, 272, 330
and 395). The v1.0 download has no `STUDY ACCESSION` or `MAPPED_TRAIT_URI`, so
`GWAS_study` and `GWAS_accessions` are empty, as Ensembl VEP writes them.

**Several records at one site.** The record reported is the first in the order
of the plugin's processed file, which GNU `sort` writes in the English UTF-8
locale. Letters and digits come first, without case and with punctuation
ignored, then case, then punctuation.

## DosageSensitivity {#dosagesensitivity}

DosageSensitivity adds a gene's probabilities of haploinsufficiency and
triplosensitivity from the dosage sensitivity catalogue of Collins et al., 2022,
to a transcript consequence, as Ensembl VEP's `DosageSensitivity.pm` does.

| Property      | Value                                                                                                                       |
| ------------- | --------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                              |
| Matches       | gene symbol lookup on transcript rows, of sequence and structural variants alike                                            |
| Output fields | `pHaplo`, the probability that the gene is haploinsufficient; `pTriplo`, the probability that it is triplosensitive           |
| Data source   | `Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz`, Zenodo record 6347673                                                  |

Every parameter is `key=value`.

| Parameter     | Meaning                                                                                                               |
| ------------- | --------------------------------------------------------------------------------------------------------------------- |
| `file=<path>` | The scores table (mandatory): three columns, the gene symbol, `pHaplo` and `pTriplo`, read whole at start-up, gzip-compressed or plain; Ensembl VEP 116.2 reads the compressed form only (`Utils.pm` line 434). |
| `cover=1`     | Write the scores only on a row whose variant spans the whole transcript, such as a deletion or duplication of the gene. `cover=0` is the default; every value but the string `1` is `0`, as Ensembl VEP 116.2 reads it (`DosageSensitivity.pm` line 101). |

```bash
--plugin DosageSensitivity,\
file=/path/to/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
--plugin DosageSensitivity,\
file=/path/to/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz,cover=1
```

## Where these plugins differ from Ensembl VEP

Nine behaviours of Ensembl VEP 116.2's plugins are not reproduced. The first
two are classes the comparison tools set aside. The other seven are documented
classes no tool codes, since no committed row falls under them. Each section
names the rows a comparison would set aside.

[Intended divergences from Ensembl VEP](../intended-divergences.md) lists all
nine; its
[How to read the counts](../intended-divergences.md#how-to-read-the-counts)
defines the format comparator and the warnings comparison the sections name.
The concordance figures score consequence sets and never read a plugin field,
so nothing is set aside in them. Every row below is in the default output
format, against the Ensembl 116 GRCh38 cache.

Four of the documented classes are sections of this page. The GWAS plugin's
three others are:

- [A GWAS record's match that depends on the record before it at the same site](../divergences/gwas-shared-site-previous-record-reference.md)
- [No GWAS annotation at position 1 or 2 of a chromosome](../divergences/gwas-window-below-position-one.md)
- [A GWAS summary-statistics file with an empty cell before `hm_pos`](../divergences/gwas-sstate-index-failure.md)

Sections 3 to 6 carry their tier under the rule of [Intended divergences from
Ensembl VEP](../intended-divergences.md), which lists every class's tier. Tier 1
is a shape the chromosome 21 population cuts carry with no flag. Tier 2 is one
that needs a flag or plugin on ordinary data. Tier 3 is one ordinary data does
not carry.
The counts are over the cuts' 15,646,481 records, the 111 golden corpora and
their 392 flag sets.

### 1. GWAS: every SNP of a haplotype row at every location of the row

**In short.** A GWAS Catalog row for a haplotype or a SNP x SNP interaction
lists several SNPs and as many locations. Ensembl VEP's GWAS plugin hands every
location of such a row to every SNP, so each site carries every SNP's risk
allele and study. vep-rs pairs the first SNP with the first location, the second
with the second, and so on. On the gnomAD v4.1 chromosome 21 sites against the
full catalog, 590 rows carry an annotation from Ensembl VEP alone and 139 a
different study or p-value on the two engines, every one at a site of such a
row.

The first SNP of an interaction row is the one exception in Ensembl VEP 116.2.
Its risk allele keeps a trailing space (`GWAS.pm` lines 346-348) and never
matches, so no site carries it ([section 3](#3-gwas-the-first-snp-of-an-interaction-row-never-matches-its-own-risk-allele)).

#### 1.1 What a reader needs to know first

The catalog's `SNPS`, `CHR_ID`, `CHR_POS` and `STRONGEST SNP-RISK ALLELE`
columns hold `;`-separated lists of equal length on a haplotype row, one entry
per SNP. Each risk allele is written `rs<id>-<base>`. A record is annotated
with a row when its position is one of the row's locations and its allele is
the risk allele there. The record then carries the row's study, p-value and
trait.

#### 1.2 The record

Two gnomAD v4.1 chromosome 21 records, `21 45427844 C T` and
`21 38297596 C T`, against the full associations download. The first lies at
the first location of a nine-SNP haplotype row of study `GCST008413`:

```
CHR_POS=45427844;45428241;45444973;45451558;45463239;45466796;45472737;45474107;45477179
SNPS=rs2838925; rs11911327; rs10854470; rs2838932; rs4819119; rs2838933; rs8126650; rs2236458; rs2236461
STRONGEST SNP-RISK ALLELE=rs2838925-C; rs11911327-C; rs10854470-C; rs2838932-C; rs4819119-A; rs2838933-T; rs8126650-T; rs2236458-C; rs2236461-T
```

The second lies at the third location of two six-SNP haplotype rows of the
same study, whose SNP there is `rs2236607` with the risk allele `C`, and at
the single location of a row of study `GCST90245845`, `rs2236607-T`.

#### 1.3 What each engine writes

The first record on `ENST00000397787`. Ensembl VEP:

```
21_45427844_C/T	21:45427844	T	ENSG00000183535	ENST00000397787	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2804;STRAND=-1;GWAS_accessions=http://purl.obolibrary.org/obo/MONDO_0018874;GWAS_associated_gene=;GWAS_beta_coef=;GWAS_odds_ratio=;GWAS_p_value=3E-8;GWAS_pmid=27903959;GWAS_risk_allele=T;GWAS_study=GCST008413
```

vep-rs:

```
21_45427844_C/T	21:45427844	T	ENSG00000183535	ENST00000397787	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2804;STRAND=-1
```

The second record on `ENST00000328656`. Ensembl VEP:

```
21_38297596_C/T	21:38297596	T	ENSG00000157551	ENST00000328656	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;GWAS_accessions=http://purl.obolibrary.org/obo/MONDO_0018874;GWAS_associated_gene=;GWAS_beta_coef=;GWAS_odds_ratio=;GWAS_p_value=3E-15;GWAS_pmid=27903959;GWAS_risk_allele=T;GWAS_study=GCST008413
```

vep-rs:

```
21_38297596_C/T	21:38297596	T	ENSG00000157551	ENST00000328656	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;GWAS_accessions=http://purl.obolibrary.org/obo/OBA_VT0001253;GWAS_associated_gene=;GWAS_beta_coef=0.027 unit decrease;GWAS_odds_ratio=;GWAS_p_value=3E-25;GWAS_pmid=36224396;GWAS_risk_allele=T;GWAS_study=GCST90245845
```

#### 1.4 Why Ensembl VEP writes that

`GWAS.pm` splits a row's `CHR_ID`, `CHR_POS` and risk-allele lists (lines
233-235). For each SNP id of the row in turn (lines 334-355) it builds the
SNP's entry. The entry holds every location of the row, since line 338 calls
`get_vfs_from_file` with the whole lists, and that SNP's risk allele (lines
343-350).

The first record's site, `45427844`, therefore carries nine entries, one per SNP,
with the risk alleles `C`, `C`, `C`, `C`, `A`, `T`, `T`, `C` and
`T`. The record's `T` matches three of them and the row's study is written. At
the second record's site the same mechanism gives the `T` allele the haplotype
rows' study through another SNP's risk allele. Of the matching rows the
plugin reports the first in its processed file, a haplotype row, not the site's
own single-SNP row.

#### 1.5 Why vep-rs's row is the correct one

The catalog states one risk allele and one location per SNP. The first
record's site is `rs2838925`, whose risk allele is `C`, not `T`, so no
association is reported there. At the second site the SNP is `rs2236607`, whose
risk allele in the haplotype rows is `C`. The record's `T` matches only the
single-SNP row `rs2236607-T` of `GCST90245845`, which vep-rs reports.
The `T` alleles Ensembl VEP matched belong to SNPs at other positions.

#### 1.6 Try it yourself

```
vep -i records.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin GWAS,file=gwas_catalog_associations.tsv
```

`records.vcf` holds the two records of 1.2, `<cache>` is the Ensembl 116
GRCh38 cache and the TSV is the full associations download.

| Record | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `21 45427844 C T` | `GWAS_study=GCST008413`, `GWAS_risk_allele=T` | no `GWAS_` field |
| `21 38297596 C T` | `GWAS_study=GCST008413`, `GWAS_p_value=3E-15` | `GWAS_study=GCST90245845`, `GWAS_p_value=3E-25` |

Behind both are the catalog's `GCST008413` lines holding `45427844` and
`38297596` in `CHR_POS` and its `GCST90245845` line holding `rs2236607`.

#### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies
`gwas_haplotype_cross_product` when a comparison case names it. A `GWAS_` field
the reference fills, or writes empty beside its filled fields, and the candidate
leaves empty, or both fill differently, is set aside on covered rows. Covered
is a row carrying, at one location of a catalog row (`--gwas-catalog`), the
risk allele of another of its SNPs, or of any where SNP and location counts
differ, and whose reference `GWAS_study` names that catalog row's
study. A `GWAS_` field the candidate alone fills is set aside there only where
both engines annotate the row under different studies.

The class is silent where the engines agree. The `GRCh38-plugin-gwas` corpus
has records at three such locations (21:34300618, 21:39375657 and 21:45444973),
each with its own SNP's risk allele, a base other SNPs of its one catalog row
may share.

### 2. DosageSensitivity: a warning per row kind on regulatory and motif rows

**In short.** Under `--regulatory`, Ensembl VEP's DosageSensitivity plugin is
run on regulatory-feature and motif-feature rows. There it reads a transcript
the row does not have and dies on every such row. Ensembl VEP catches each
exception and writes a `Plugin 'DosageSensitivity' went wrong` warning, one
line per distinct message, so one per row kind in a run. The rows themselves
are written without a field by both engines, and vep-rs writes no warning.

#### 2.1 What a reader needs to know first

A plugin names the row kinds it runs on. DosageSensitivity's module names
`Feature`, which Ensembl VEP reads as every kind, so the plugin is run on
regulatory and motif rows as well as transcript rows. Ensembl VEP writes each
distinct warning text once per run (`BaseVEP.pm`, lines 809-813, keep the
texts already written). The message names the row's class, so ten failing
rows of two kinds give two lines.

#### 2.2 The record

The first five records of the `GRCh38-regulatory` golden corpus
(`tests/golden/116/GRCh38-regulatory/variants.vcf`), which overlap regulatory
features and motifs, with `--regulatory --plugin DosageSensitivity,file=<file>`
and the published `Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz`.

#### 2.3 What each engine writes

Both engines write the same ten regulatory and motif rows, the first:

```
del_covers_rf	21:5014322-5014344	-	-	ENSR21_536SRH	RegulatoryFeature	regulatory_region_ablation,regulatory_region_variant	-	-	-	-	-	-	IMPACT=MODIFIER;BIOTYPE=CTCF_binding_site
```

Ensembl VEP's warnings file carries two lines, one per row kind:

```
WARNING: 11 : WARNING: Plugin 'DosageSensitivity' went wrong: Can't locate object method "transcript" via package "Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele" at /plugins/DosageSensitivity.pm line 169, <__ANONIO__> line 7.
WARNING: Plugin 'DosageSensitivity' went wrong: Can't locate object method "transcript" via package "Bio::EnsEMBL::Variation::MotifFeatureVariationAllele" at /plugins/DosageSensitivity.pm line 169, <__ANONIO__> line 7.
```

vep-rs's warnings file has no line for the records.

#### 2.4 Why Ensembl VEP writes that

`DosageSensitivity.pm`, line 169, reads `$tva->transcript->{_gene_symbol}` on
every row it is run on. A regulatory or motif row's allele object has no
`transcript` method, so the call dies. `OutputFactory.pm`, line 2244,
catches the exception from the plugin's `run` and writes it as the warning.
The row is then written without the plugin's fields.

#### 2.5 Why vep-rs's row is the correct one

A regulatory feature has no gene, so the plugin has nothing to look up on such a
row and nothing to report. Writing the row without a field and without a
warning is what the plugin's own `Feature` declaration describes for a feature
without a gene symbol. The values on every transcript row agree.

#### 2.6 Try it yourself

```
vep -i variants.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 --regulatory \
  --plugin DosageSensitivity,file=Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
```

`<corpus>` is `tests/golden/116/GRCh38-regulatory`, and `variants.vcf` holds
the five records of 2.2.

Compare `out.txt_warnings.txt`: Ensembl VEP's carries the two lines above,
vep-rs's none.

#### 2.7 How the comparator treats this class

The format comparator applies `dosage_sensitivity_regulatory_row_warnings` in
`--warnings` mode when a comparison case names it. The reference's `Plugin
'DosageSensitivity' went wrong` lines are left out of the comparison, each
counted under the class. The output rows compare as usual. No golden corpus
runs the plugin under `--regulatory`, so the class applies to none of their
warnings files.

### 3. GWAS: the first SNP of an interaction row never matches its own risk allele

**The defect in one sentence.** On a `SNP x SNP` interaction row, Ensembl VEP's
GWAS plugin keeps the space before the ` x ` in every risk allele but the last,
so `G ` matches no variant and its SNP's site gets no `GWAS_` field.

**Why this is a defect.** The plugin splits the row's `CHR_ID`, `CHR_POS` and
`STRONGEST SNP-RISK ALLELE` lists on the same `;`, `,` and `x` separators. It
strips the whitespace from the first two lists but not from the third. The
catalog writes a risk allele as `rs<id>-<base>`; the space after the base is
list punctuation, not allele text. The same row's last SNP, whose token the
` x ` precedes, is matched by its allele as written.

**What vep-rs does instead.** vep-rs trims the whitespace around each risk
allele, so every SNP of an interaction row is matched by its own base at its own
location.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116,
`ensembl-variation` at release/116):

- `GWAS.pm`, lines 343-348, `parse_curated_file`: the risk-allele list is split
  on `/[;,x]/` and the second `-`-separated field of the token naming the SNP is
  taken as written, so `rs10830963-G x rs4818032-C` gives `G ` and `C`. Lines
  233-235 strip the whitespace of the location lists the same split produces.
- `GWAS.pm`, lines 395-418 and 487: the allele is written into the processed
  file and read back with its space.
  `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 1207 and 1232,
  `get_matched_variant_alleles`: the key `<ref>_<alt>_<pos>` is an exact string.
  `trim_sequences` (lines 965-1038) removes only bases the reference allele
  shares, never the space.

**What correcting it would change.** The first SNP of every interaction row
would annotate variants carrying its risk allele at its own location, as the
last SNP of the row does in both engines.

**The class it produces.** `gwas_interaction_first_snp_allele_whitespace`: a
`GWAS_` field (the JSON `gwas_` keys) that vep-rs alone fills, on a row at the
site of an interaction row's first SNP. No committed comparison carries such a
row. The `GRCh38-plugin-gwas` corpus's one interaction row names chromosome 11
first, and the full associations download's 127 interaction rows with stated
alleles have no first SNP on chromosome 21. The gnomAD v4.1 chromosome 21
comparison therefore sets aside 0 rows under it.

**How often the shape occurs.** Tier 2: the shape is a catalog row, not a
record. 1 of the 790 rows of the committed catalog slice is an interaction row
joining risk alleles with ` x `. 2 flag sets run the plugin. A site carries
the defect when such a row's first SNP names it.

#### 3.1 The chromosome 21 SNP of an interaction row, written first

**In short.** The corpus record `21 39375657 rs4818032 T C`, run against the
catalog row `rs10830963-G x rs4818032-C` with its two SNPs written in the other
order (`rs4818032-C x rs10830963-G`, `21 x 11`). Ensembl VEP writes no `GWAS_`
field on any of the record's 24 rows: at `21:39375657` its processed file holds
the alleles `C `, with the space, and `G`, neither equal to the record's `C`.
vep-rs writes
`GWAS_p_value=3E-7;GWAS_pmid=29691896;GWAS_risk_allele=C;GWAS_study=GCST006494`
among the row's `GWAS_` fields on each. The catalog states `C` as rs4818032's risk allele at that position, the
record's own allele, and Ensembl VEP itself matches the row's last SNP by the
same reading.

### 4. GWAS: the last column of a summary-statistics header is never read

**The defect in one sentence.** Under `type=sstate`, Ensembl VEP's GWAS plugin
reads the header of a harmonised summary-statistics file with the line end
still on its last column name, so the last column maps to no field and the
field it would fill prints empty.

**Why this is a defect.** The plugin looks up six column names (`hm_chrom`,
`hm_pos`, `p_value`, `hm_beta`, `hm_odds_ratio`, `hm_effect_allele`) by exact
match against the header. A file may then order its columns as the catalog's
harmonisation leaves them. The last name read ends in a newline the code never
removes, and `p_value` followed by a newline equals no name. A file whose last
column is `p_value` therefore annotates every record with `GWAS_p_value=`.

**What vep-rs does instead.** vep-rs reads the header line without its line end
and maps every column by name, the last one included.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `GWAS.pm`, lines 437-448, `parse_sstate_header`: the header comes from
  `head -n 1` with its trailing newline, and `split("\t", ...)` leaves it on
  the last name. `grep {/^$cols[$index]$/}` over the six required names matches
  a pattern carrying that newline against none of them.
- `GWAS.pm`, lines 496-509, `parse_data`: a column with no mapping is never
  read, and `|| ""` prints the field empty. With `hm_chrom` or `hm_pos` last,
  the plugin has no position column and fails to instantiate (lines 127-136 and
  460), annotating nothing.

**What correcting it would change.** The field named by the header's last column
would print the file's value, as the other five do.

**The class it produces.** `gwas_sstats_last_header_column_unmapped`: under
`type=sstate`, the one `GWAS_` field the header's last column names
(`GWAS_p_value` for `p_value`, `GWAS_beta_coef` for `hm_beta`,
`GWAS_odds_ratio` for `hm_odds_ratio`) that Ensembl VEP leaves empty and vep-rs
fills, the row's other `GWAS_` fields agreeing. The committed file
`17463246-GCST000028-EFO_0001360.h.tsv` ends in `ci_upper`, so the
`GRCh38-plugin-gwas` corpus's `sstate` run has 0 such rows in each format. No
other committed comparison runs a summary-statistics file.

**How often the shape occurs.** Tier 3: the shape is a header whose last column
is one of the six the plugin maps. Neither of the 2 committed summary-statistics
files has one, and 2 flag sets run the plugin.

#### 4.1 A nine-column file ending in `p_value`

**In short.** A harmonised file whose header is `hm_variant_id hm_rsid hm_chrom
hm_pos hm_other_allele hm_effect_allele hm_beta hm_odds_ratio p_value`, with the
one row `21_34036428_A_G rs4817614 21 34036428 A G NA NA 0.1965`, and the
record `21 34036428 rs4817614 A G`. Ensembl VEP writes `GWAS_p_value=` and
vep-rs `GWAS_p_value=0.1965`; the other seven `GWAS_` fields agree,
`GWAS_risk_allele=G` and `GWAS_study=GCST000028` among them. The file states
0.1965 in the column the plugin itself names as the one it reads for the
p-value.

### 5. LoFtool and pLI: a score of `0`

**The defect in one sentence.** Ensembl VEP's LoFtool plugin, and its pLI plugin
under `sigfig`, format a gene's score and then test the text for truth before
writing it, so a score of exactly `0`, formatted as `0`, is written on no row.

**Why this is a defect.** Both plugins format the score with `sprintf` and
return the field only when the formatted string is true in Perl. `0` is one
of the two strings Perl treats as false. The same cell prints under pLI's
default format, where `%.2f` writes it as `0.00`, a true string. A LoFtool
percentile of `0` is the most intolerant rank the file can state; a score the
file carries is a value, not an absence.

**What vep-rs does instead.** vep-rs writes every score the file carries for the
row's gene or transcript, `0` included, in the format the parameters select.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `LoFtool.pm`, line 89: `sprintf("%g", $score)` writes a zero as `0`. Line
  118: the score is returned through
  `$self->{scores}->{lc($symbol)} ? { LoFtool => ... } : {}`, a truth test on
  that string.
- `pLI.pm`, line 110: `%.3g` under `sigfig`, else `%.2f`. Lines 145, 156, 169
  and 182: every score formatted at load. Lines 221 and 227: the same truth
  test on the formatted string, in gene and in transcript mode.

**What correcting it would change.** Every gene whose file score is `0` would
carry `LoFtool=0`, or under `sigfig` `pLI_gene_value=0` or
`pLI_transcript_value=0`, as it carries `0.00` without `sigfig`.

**The class it produces.** `gene_score_zero_value`: a `LoFtool`,
`pLI_gene_value` or `pLI_transcript_value` field (the JSON `loftool`,
`pli_gene_value` or `pli_transcript_value`) the candidate alone fills, with `0`.
No committed row: the `LoFtool_scores.txt` (14,515 entries) and
`pLI_values.txt` (18,225 entries) the plugins are distributed with, committed
under `tests/fixtures/plugins/`, and the score files of the four committed
LoFtool and pLI corpora carry no zero score. The smallest pLI value of the four
corpora's tables, `5.11e-31`, prints as itself under `sigfig`.

**How often the shape occurs.** Tier 3: the shape is a score cell of exactly
zero. None of the 690 rows of the golden corpora's score tables carries one,
and 14 flag sets run one of the two plugins.

#### 5.1 A percentile of `0` for SOD1

**In short.** Against a scores file that gives SOD1 the percentile `0.00E+00`,
the GRCh38 record `21 31667290 14766 A C`
(`tests/golden/116/GRCh38-plugin-loftool/variants.vcf`) has 26 SOD1 transcript
rows. Ensembl VEP 116.2 writes no `LoFtool` on any of them. vep-rs writes
`LoFtool=0` on all 26, as both engines write `LoFtool=0.126` against the bundled
file's `1.26E-01`. With pLI, a cell `0` for SOD1 prints `pLI_gene_value=0.00` on
both engines by default and under `sigfig` on vep-rs alone, as
`pLI_gene_value=0`.

### 6. pLI: a table cell that is not a number is printed as `0.00` {#6-pli-a-table-cell-that-is-not-a-number}

**The defect in one sentence.** Ensembl VEP's pLI plugin formats every score
cell with `sprintf` as it loads the table, and Perl formats text that is not a
number as 0, so a gene whose cell reads `NA` or is missing gets
`pLI_gene_value=0.00` on each of its transcript rows.

**Why this is a defect.** The plugin's documentation describes the table as each
gene and its pLI, the probability that the gene is intolerant of loss of
function. `NA` in a constraint table states that no estimate exists, and the
gnomAD v4 constraint table the documentation points at writes `NA` in its pLI
column for such genes. A printed `0.00` claims the gene is as tolerant as a gene
can be. The same cell under `sigfig` formats as `0`, which the plugin's own
truth test drops, so Ensembl VEP itself prints nothing for it there.

**What vep-rs does instead.** vep-rs reads each cell as its leading number. A
cell without one, empty or `NA`, gives its gene or transcript no field in both
modes, with or without `sigfig`. The warnings file names the gene and the cell
once.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `pLI.pm`, lines 145, 156, 169 and 182, `new`: `sprintf($score_format,
  $score)` on every cell, `%.2f` by default and `%.3g` under `sigfig`. Lines
  221 and 227, `run`: the formatted string is written when it is true, so
  `0.00` prints and `0` does not.
- `LoFtool.pm`, line 89, runs the same `sprintf("%g", $score)`, which makes `0`
  of such a cell. Its truth test at line 118 drops it, so LoFtool's rows agree
  with vep-rs's.

**What correcting it would change.** A cell that is not a number would be
skipped, or refused with its line, and no row would carry a score the table
never stated.

**The class it produces.** `pli_non_numeric_score_prints_zero`: a
`pLI_gene_value` or `pLI_transcript_value` (JSON `pli_gene_value`,
`pli_transcript_value`) the reference alone fills, with `0.00`, on a row whose
gene or transcript has a cell that is empty or has no leading number. No row on the
committed corpora or their comparison cases: the four tables of the GRCh37 and
GRCh38 pLI corpora and the bundled `pLI_values.txt` hold numeric cells only.
Every pLI comparison case therefore has 0 rows in each format. The LoFtool corpora, which
the class does not reach, have 0.

**How often the shape occurs.** Tier 3: the shape is a pLI table cell that is
not a number. None of the 643 rows of the committed pLI tables carries one, and
10 flag sets run the plugin.

#### 6.1 CXADR with its cell written as `NA`

**In short.** The record `4453310`, `21 17561425 G T` of the GRCh38 pLI corpus,
lies in CXADR, whose cell in the corpus table is `0.342891274842555`. Both
engines write `pLI_gene_value=0.34` on its 33 transcript rows. With the cell
changed to `NA` in a copy of the table, vep-rs writes the 33 rows without the
field. Ensembl VEP writes `pLI_gene_value=0.00` on each of them, a score the
table does not state, and warns `Argument "NA" isn't numeric in sprintf` on
standard error.

## Where next

- [Plugins](../plugins.md): every plugin's support status and how one is
  enabled.
- [Intended divergences from Ensembl VEP](../intended-divergences.md): the
  index that lists all nine classes.
- [Reproducing these examples](../divergences/reproducing.md): the engines,
  inputs and commands behind the examples.
