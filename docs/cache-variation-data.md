# Variation data in a native cache

A native cache is one `vep-cache-builder` builds. This page covers the variation files it holds, the form of their records and what they lack.

Variation data, which `--check_existing` and the `--af*` flags read, comes with an `--auto` installation of a genome the release publishes variation VCFs for, and with `--include-variations` in a build with `--output-dir`. The builder downloads Ensembl's variation VCFs from this directory (`pub/grch37/release-116/` for GRCh37), or reads local copies under the same names from `--variation-vcf-dir`:

```
https://ftp.ensembl.org/pub/release-116/variation/vcf/homo_sapiens/
```

For `homo_sapiens` the files are the per-chromosome germline files `homo_sapiens-chr<N>.vcf.gz` and four genome-wide files, each with its `.csi`. The download covers chromosomes 1-22, X, Y and MT unless `--chromosomes` narrows it. Any other species has one file, `<species>.vcf.gz`, read for every sequence with records; `--chromosomes` narrows the read. From a local directory every `*-chr*.vcf.gz` is read, else the one `<species>.vcf.gz`; a `##fileformat` outside major version 4 stops the build.

The four genome-wide files:

| Genome-wide file | Contents |
| --- | --- |
| `homo_sapiens_somatic.vcf.gz` | every record flagged `somatic` |
| `homo_sapiens_phenotype_associated.vcf.gz` | records the per-chromosome files already carry at release 116; a record found only there is added |
| `homo_sapiens_clinically_associated.vcf.gz` | records the per-chromosome files already carry at release 116; a record found only there is added |
| `1000GENOMES-phase_3.vcf.gz` | the 1000 Genomes frequencies |

A genome-wide file the directory or the site lacks is left out with a warning naming what the cache then lacks. An installation of a genome whose tree publishes none of the files holds no variation data, and its last lines say so, naming the listing. A run that asks such a cache for known variants warns once, naming the flags that ask, and says that a frequency filter then keeps every record. `vep recode` over such a cache names in one `INFO:` line the fields among `id` and `var_synonyms` that stay unwritten. A human listing that names some per-chromosome files and lacks another, or lacks one `--chromosomes` names, stops the installation before any download.

The records are written per chromosome as `variations/<chr>/all_vars.gz` with a `.csi` on the start column, the file Ensembl VEP's own cache carries. So `--check_existing` and the frequency flags read a native cache exactly as a converted one.

## What an installation reads beside the VCFs

An `--auto` installation fills more of each record from the release's other files, each where the genome's tree lists it:

| What the installation adds | Read from |
| --- | --- |
| `pubmed`, the PubMed ids of the publications that cite a variant | The variation database's `variation`, `variation_citation` and `publication` table dumps |
| `var_synonyms` | Its `variation`, `variation_synonym` and `source` table dumps; a source whose name holds `dbSNP` is left out, as in Ensembl VEP's dump |
| `failed`, and the failed variants the VCFs leave out | `<species>_failed.gvf.gz` of the release's GVF directory, with the `failed_variation` table dump where the tree lists it |
| The known variants of a sequence on which the one VCF of a species other than human holds no record | The displayed rows of the `variation_feature` table dump, fetched only when the cache has such a sequence |

A run leaves a failed variant out unless `--failed 1`, as on Ensembl VEP's cache. The GRCh37 tree lists the failed-variants file and no `variation` table, so a GRCh37 installation holds the failed variants and no `pubmed` or `var_synonyms`. The installation's last lines name each cell whose file the tree lacks, with the listing. A build with `--include-variations` reads none of these files.

The table dumps add about 10 GB to the downloads of a GRCh38 installation, `variation` 7.5 GB and `variation_synonym` 2.5 GB of it; a GRCh37 installation fetches neither. `--chromosomes` narrows the variation VCFs and not the table dumps. The `variation_feature` dump is fetched whole: about 4 MB for yeast, 213 MB for Arabidopsis, 357 MB for zebrafish and 2.6 GB for mouse.

### Where an installed cache differs from Ensembl VEP's

The published files cause four differences from Ensembl VEP's own cache:

- The published `variation_synonym` table does not record the order of one source's names in Ensembl's dump: of chromosome 21's 92,226 `var_synonyms` cells on GRCh38, 7,101 hold a source's names in another order and the other 85,125 are equal as stored.
- A `pubmed` or `var_synonyms` list longer than 1,024 characters, which Ensembl VEP's cache [cuts there](divergences/colocated-list-cut-at-cache-cap.md), is stored at its whole items within that length, as a run prints it: 2 of chromosome 21's 10,539 `pubmed` cells are equal as printed and not as stored.
- Chromosome 21 of GRCh38 holds all 178 failed records of Ensembl VEP's cache, 169 of them added from the failed-variants file, which gives no 1000 Genomes frequency where that cache has one on 9 and no clinical significance where it has one on 1.
- The `variation_feature` table is read for a sequence without a VCF record alone, so yeast's `IX`, which has VCF records, lacks the one displayed row no VCF or GVF holds: 261,055 of the 261,056 nuclear records of Ensembl VEP's cache.

## The form of a record

Each record takes the form the variation database holds, as far as the export (Ensembl's variation VCFs) preserves it. A deletion's alleles have the anchor base the export prepended removed, with `-` for the emptied allele. An insertion the export marks `TSA=insertion` sits two bases after the base its VCF line names (`start` one past `end`), where the Ensembl cache places it.

An indel (`TSA=indel`, every allele non-empty) keeps the alleles the export wrote, at the line's `POS`. The export prepends an anchor only to alleles that do not share a first base. So an `A/AA` record stays `A/AA`, while a `CC/CTTCT` line, whose alleles differ once the shared base goes, is the database's `C/TTCT` one base on.

A COSMIC variant without alleles is written as `COSMIC_MUTATION` on its span. Every HGMD-PUBLIC variant, alleles or not, is written as `HGMD_MUTATION`, the string the variation database holds for that source. `clin_sig` comes from the record's `CLIN_<term>` flags (a term the export declares no flag for, `likely_risk_allele`, is absent) and `phenotype_or_disease` from `E_Phenotype_or_Disease`.

`clin_sig_allele`, which a build fills from the ClinVar records of the variation database's tables, lists their `allele:term` pairs in the order of the records' identifiers. Ensembl VEP's cache holds its pairs in [Perl's hash order](divergences/clin-sig-allele-term-order.md), so `CLIN_SIG` can list one allele's terms in another order on the two caches. A record whose significance holds several classifications, such as `uncertain significance; risk factor`, gives a pair for each. Ensembl VEP's cache holds [a pair for the first alone](divergences/clin-sig-allele-value-cut-at-semicolon.md), so `CLIN_SIG` can list a term more on a native cache.

## An indel that can shift along a run

The Ensembl cache keeps dbSNP's alleles over the whole run an indel can shift along (`CCCCC/CCCC` at 5030241-5030245), which the export trims to `CC/C` at 5030241. Such records are 5 of every 100 of chromosome 21's `indel` records. A native cache holds the trimmed form, which minimises from the right end to the same alleles and position as the Ensembl cache's record. So an input spelled at the left end of the run (the VCF convention, `AC/A` at 5030240) co-locates with it exactly as Ensembl VEP co-locates it.

An input anchored on the run's first base (`CC/C` at 5030241) minimises to the record's second base. The trimmed form also reaches that base and the Ensembl cache's record does not, so the input co-locates with the native record alone. An input spelled at the run's right end co-locates with neither, its start lying outside the one-base window the start-indexed files answer.

## What a native cache lacks

No native cache carries a gnomAD frequency. Ensembl publishes no gnomAD VCF for release 116, so gnomAD frequencies come from the Ensembl variation cache, by conversion.

The export leaves the failed variants out, the `dbSNP_novariation` records among them, and carries no `pubmed` or `var_synonyms`. An installation brings the three from the release's other files ([above](#what-an-installation-reads-beside-the-vcfs)). A build with `--include-variations` reads the export alone: its cache carries no `pubmed` or `var_synonyms` and no `failed` flag, and the few failed variants the export keeps are written unflagged. A native cache built by a vep-rs release before 0.5.0 holds JSON variation files that nothing reads; build it again.

The dbSNP, ClinVar, COSMIC and gnomAD versions an Ensembl cache's `info.txt` lists come from the variation database. The builder takes no version from it, so a native cache's headers print none.

A native cache holds a known variant only on a sequence the builder read a variation record on. The one file of a species other than human can hold no record on a sequence whose known variants Ensembl VEP's own cache carries. The yeast file of Ensembl Genomes release 63, `saccharomyces_cerevisiae.vcf.gz`, holds none on `Mito`, where the release's GVF file holds 2,481. An installation reads such a sequence from the `variation_feature` table, 2,481 known variants on `Mito`; a build with `--include-variations` reads the file alone.

For `homo_sapiens` the release publishes a file for each of chromosomes 1-22, X, Y and MT and for no scaffold, and Ensembl VEP's own caches of release 116 hold no known variant on a scaffold either.

An installation's last lines name the sequences of the cache it read no variation record on, the first ten by name and the rest by their count. A conversion of Ensembl VEP's cache carries any known variant that cache holds there.

The cache's `info.json` lists under `left_out` what its build left out. A run names once, in one line, the columns its command asks for and the cache does not carry: those names, and any frequency column `variation_cols` does not declare. A cache whose `info.json` has no `left_out` key draws no line for a name.

## The 1000 Genomes frequencies

The 1000 Genomes file exists for GRCh38 only, so a GRCh37 native cache has no frequency column. Where it is read, the `AFR`, `AMR`, `EAS`, `EUR` and `SAS` columns hold each allele's continental frequency as the file gives it, and `AF` the whole panel's. The file states no panel-wide frequency, so `AF` is the allele's count over the panel's 5,008 chromosomes. Each population's count is recovered from its frequency and sample count with two chromosomes per sample (an approximation on X and Y). A population the file gives no value for drops out of both the count and the total. The build's last lines say that `AF` is recovered, as does a run that reads it.

`MAX_AF` and `MAX_AF_POPS` range over the populations a cache holds, on a native cache the 1000 Genomes populations alone, which a run under `--max_af` says in one line. On a cache with known variants and no frequency column, a native GRCh37 cache among them, the line of the empty columns also says what a frequency filter then does.

An allele whose frequency is 0 in every population gets no entry. The release file lists no frequency for an allele absent from its call set, and the export does not tell such an allele from one the release lists with a count of 0. For the latter the Ensembl cache prints `0` (among the 1,106,899 chromosome 21 records the Ensembl cache holds frequencies for). The release's own panel-wide frequency, which the Ensembl cache prints, differs from the recovered count's for 443 of those records, 315 of them by one unit in the last significant digit.

## Where next

- [Cache setup](cache-setup.md#build-a-cache): the builder's command and its other inputs.
- [Cache conversion](cache-conversion.md#variation-files): the variation files a converted Ensembl VEP cache carries.
- [Co-located variants, frequencies and known structural variants](output/known-variants.md): what a run prints from the variation files.
- [Intended divergences](intended-divergences.md): where the two engines differ, and why.
