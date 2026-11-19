# Plugins {#plugin-guide}

vep-rs provides built-in ports of Ensembl VEP plugins, loaded by their Ensembl VEP names and parameters, and loads C-ABI shared-library plugins.

Built against the `vep-plugin` crate, they load without a rebuild. A plugin adds fields to the output:
the Extra column in VEP format, CSQ in VCF, fields in JSON. This page lists every
Ensembl VEP plugin with its support status, shows how one is enabled and points
to the page documenting each.

## Usage

Plugins are loaded with `--plugin Name,param1,param2,...`, given as many times
as needed:

```bash
SNV=/path/to/your/snv.tsv.gz
INDELS=/path/to/your/indels.tsv.gz
vep -i input.vcf -o output.txt --offline --cache \
    --json_cache /path/to/json_cache --assembly GRCh38 \
    --plugin CADD,snv=$SNV,indels=$INDELS \
    --plugin REVEL,file=/path/to/revel/new_tabbed_revel_grch38.tsv.gz \
    --plugin gnomADc,/path/to/gnomad/gnomad.exomes.coverage.tsv.gz
```

`--json_cache` names the cache directory, the one holding `info.json` and
`transcripts/` ([Cache setup](cache-setup.md)).
Preparing the data files, and the `.vpd` store that replaces a tabix index, is
on [Preparing data files](plugins/data-files.md).

Pass `--assembly GRCh37` or `--assembly GRCh38` whenever plugins are enabled. The
flag is load-bearing, not just cache selection.

REVEL, dbNSFP and dbscSNV use it
to choose which coordinate column their index was built on. REVEL refuses to
load when neither the flag nor the cache names an assembly. CADD, REVEL, dbNSFP,
dbscSNV and AlphaMissense use it to refuse a data file whose own assembly
contradicts the run, rather than silently annotating nothing. Those guards are
inert when neither names an assembly.

Without the flag, the cache's own assembly serves every plugin that reads one.
Ensembl VEP 116.2 builds its plugins before it reads the cache (`Runner.pm`
lines 125-128). There a plugin that reads the assembly as it loads, REVEL,
dbscSNV and PrimateAI among them, refuses a run without `--assembly` and is
dropped with a warning. CADD's header check (`CADD.pm` line 141) then has no assembly to test.
A GRCh37 file loads on a GRCh38 cache and its scores are written, where vep-rs
refuses the file from the cache's assembly.

Each plugin's parameters are on its page. A parameter list is all positional or
all `key=value`. A value is the second `=`-field alone, as Ensembl VEP 116.2
splits it (`BaseVepPlugin.pm` line 126): `file=a=b` is `a`, `file==b` and
`file=` the empty value. LoFTEE's own rules are on [its page](plugins/loftee.md).

A key given twice is refused, an undeclared one too, except by CADD and
SpliceVault. Ensembl VEP keeps the later value (line 128) and never reads an
undeclared key, with three exceptions. dbNSFP refuses one (`dbNSFP.pm` line
601); SpliceVault and CADD open every other key as a data file
(`SpliceVault.pm` lines 119-122, `CADD.pm` 124-129).

## Support table

Support is per reference genome, because each annotation database decides which
builds it publishes, and per data licence, because a database whose file is not
redistributed here is compared with Ensembl VEP on a synthetic file. The table
states both; each restriction is the upstream database's, and the plugin's page
gives the reason.

One row per Ensembl VEP plugin: the plugin collection of releases 115 and 116,
plus LoFTEE and Gwava, which Ensembl VEP loads from outside the collection. Rows
are in alphabetical order of the plugin name, case ignored.

The status is one of three phrases.

- "Supported, verified on real data" means
  vep-rs's output is compared with Ensembl VEP's plugin on the published data file,
  on public variant records when the plugin reads no data file, or on a table of the
  documented layout when the file is the user's own.
- "Supported on the user's licensed file, verified on synthetic fixtures" means the
  plugin reads the file the user holds under their own licence. The comparison
  with Ensembl VEP then runs on a synthetic file in the layout the plugin's
  documentation gives, because the file is not redistributed here.
- "Not ported" means vep-rs does not provide the
  plugin, for the reason the row gives; a `--plugin` naming one stops the run with
  that reason.

The Terms column states the data provider's licence and whether it permits
commercial use, as far as the plugin or the provider states it; "terms not
stated by the plugin" means neither names any. The staged AlphaMissense cuts
keep the file's own header, with its CC BY-NC-SA 4.0 notice, verbatim, and the
corpus provenance record states both licences. A linked plugin name leads to the
page that documents it; the pages group the built-in plugins by mechanism, as
`crates/vep-builtins/README.md` does. A supported plugin without a page of its
own takes the parameters its Reason cell names, each with the meaning Ensembl
VEP's plugin gives it.

| Plugin | Data file and provider | Terms | Assemblies | Status | Reason |
| --- | --- | --- | --- | --- | --- |
| [AlphaMissense](plugins/annotation-store.md#alphamissense) | `AlphaMissense_hg19.tsv.gz` and `AlphaMissense_hg38.tsv.gz`; Google DeepMind, the 2023 release | CC BY 4.0 (the publisher's README; the file header carries the older CC BY-NC-SA 4.0 notice) | GRCh37, GRCh38 | supported, verified on real data | open data; one file per assembly, each row naming its build; missense rows whose amino-acid change the prediction names, with `cols=` and `transcript_match=` as Ensembl VEP's plugin reads them; under `transcript_match=1` Ensembl VEP matches no prediction on a two-digit transcript version and vep-rs does ([the class](plugins/annotation-store-divergences.md#1-alphamissense-a-two-digit-transcript-version-under-transcript_match1)); the published tables name the mitochondrial contig `chrM`, which a record written `MT` reaches in vep-rs alone ([the class](divergences/plugin-contig-mitochondrial-spelling.md)) |
| [AncestralAllele](plugins/ancestral-allele.md) | `homo_sapiens_ancestor_GRCh38.tar.gz` (Ensembl, release 116 FTP) and `homo_sapiens_ancestor_GRCh37_e71.tar.bz2` (release 75) | Ensembl data, no restriction on access or use | GRCh37, GRCh38 | supported, verified on real data | the record's span read from the ancestral FASTA as written; `-` for an insertion, `?` for a chromosome the FASTA lacks |
| AVADA | `avada_v1.00_2016.vcf.gz`; Bejerano laboratory, Stanford | an academic resource; terms not stated by the plugin | GRCh37 | not ported | literature-evidence PubMed identifiers from one GRCh37-only file; absent from the release 116 web manifest |
| [AVI](plugins/annotation-store.md#avi) | `avi_scores.tsv.gz` in the plugin's synopsis, a tabix-indexed, bgzip-compressed TSV file of AlphaGenome Variant Impact scores for single nucleotide variants (`AVI.pm` lines 31 and 37-38), available from the address its documentation gives (lines 40-41); Google DeepMind | a permissive-use artifact under the AlphaGenome terms: commercial and non-commercial use | GRCh38 | supported on the user's licensed file, verified on synthetic fixtures | the table is not redistributed here, so the comparison with Ensembl VEP runs on a synthetic table of the six columns the plugin reads (`AVI.pm` line 161), in all four output formats; `AVI_PHRED`, with `AVI_RAW` under `raw=1` and `AVI_LINK` under `links=1` |
| BayesDel | the `BayesDel_170824_addAF` per-chromosome files concatenated, sorted, bgzip-compressed and tabix-indexed as the Ensembl VEP plugin documents (`BayesDel.pm` lines 53-62); the authors' download | no licence text published; commercial use not stated | GRCh37; GRCh38 on a file remapped as the plugin's documentation describes (`BayesDel.pm` lines 64-65) | supported on the user's licensed file, verified on synthetic fixtures | with no licence text published the file is not redistributed here, so the comparison with Ensembl VEP runs on a synthetic file of the five columns the plugin reads (`BayesDel.pm` line 151); `file=`; position, reference and alternate allele on every row of the record, intergenic rows included; a reverse-strand Ensembl-format record: [section 2 below](#2-bayesdel-a-record-given-on-the-reverse-strand) |
| [Blosum62](plugins/transcript-model.md#blosum62) | none; the BLOSUM62 matrix is embedded in the plugin | Apache 2.0 plugin code; no data | GRCh37, GRCh38 | supported, verified on real data | the score of the row's amino-acid change, compared with Ensembl VEP on public ClinVar records |
| [CADD](plugins/annotation-store.md#cadd) | `whole_genome_SNVs.tsv.gz` and the release's indel file, CADD, University of Washington; the CADD-SV v1.1 structural-variant files, Berlin Institute of Health | free for academic and non-commercial use; commercial use needs a licence from UW CoMotion | GRCh37, GRCh38 (CADD-SV on GRCh38) | supported on the user's licensed file, verified on synthetic fixtures | the files are not redistributed here, so the comparison with Ensembl VEP runs on synthetic files whose headers name the columns the plugin reads (`CADD.pm` lines 81-98 and 265-268), the structural-variant file included; Ensembl VEP gives a structural allele outside insertion, deletion and duplication a `Plugin 'CADD' went wrong` warning where both engines write its row without a score, and matches no reverse-strand record of Ensembl's input format, where vep-rs matches the class word as written and scores it ([the class](plugins/annotation-store-divergences.md#2-cadd-structural-alleles-outside-insertion-deletion-and-duplication)) |
| CAPICE | `capice_v1.0_build37_snvs.tsv.gz` and `capice_v1.0_build37_indels.tsv.gz`, Zenodo record 3928295 | terms not stated by the plugin | GRCh37 | not ported | a GRCh37-only pathogenicity score from a single data release |
| Carol | none; the cache's SIFT and PolyPhen scores | no data | GRCh37, GRCh38 | not ported | a combined SIFT and PolyPhen score; REVEL and AlphaMissense cover the use |
| CCDSFilter | none | no data | GRCh37, GRCh38 | not ported | an example filter plugin absent from the release 116 web manifest; offline, the cached transcripts carry no database cross-references for it to filter on, and the `vep filter` tool covers the use |
| ClinPred | `ClinPred.txt.gz` (GRCh37) and `ClinPred_hg38.txt.gz` (GRCh38), tab-separated, bgzip-compressed and tabix-indexed as the Ensembl VEP plugin documents (`ClinPred.pm` lines 36-53); the ClinPred download site | non-commercial use under the authors' terms | GRCh37, GRCh38 | supported on the user's licensed file, verified on synthetic fixtures | the file is not redistributed here, so the comparison with Ensembl VEP runs on a synthetic file of the five columns the plugin reads (`ClinPred.pm` line 127); `file=` or the first positional parameter; single-base substitutions on missense, stop-gained, stop-lost and start-lost rows; a reverse-strand Ensembl-format record: [section 1 below](#1-clinpred-a-record-given-on-the-reverse-strand); a mitochondrial record against a file naming the contig otherwise: [its page](divergences/plugin-contig-mitochondrial-spelling.md) |
| Condel | the `config/Condel/` distributions bundled with the Ensembl VEP plugin collection | Apache 2.0, bundled with the collection | GRCh37, GRCh38 | not ported | a consensus of SIFT and PolyPhen scores that REVEL and AlphaMissense cover; the `--condel` flag is accepted without effect |
| [Conservation](plugins/conservation.md) | `gerp_conservation_scores.homo_sapiens.GRCh38.bw`; Ensembl Compara, the GERP scores of the 92-mammal alignment, release 116 | Ensembl data, no restriction on access or use | GRCh38 | supported, verified on real data | open data; Ensembl publishes the bigWig for GRCh38 alone, so on GRCh37 the plugin reads a bigWig of your own; the database and FTP-lookup forms of the Ensembl VEP plugin are refused by name; under `MAX` Ensembl VEP reports the score of an all-negative span that sorts last as text, which can lie below the numeric maximum vep-rs reports ([the class](plugins/conservation.md#where-max-differs-from-ensembl-vep)) |
| CSN | none; computes HGVS itself and needs `--fasta` offline | no data | GRCh37, GRCh38 | not ported | Clinical Sequencing Nomenclature strings computed from the transcript model; rarely used, and absent from the release 116 web manifest |
| DAS | a DAS server queried per variant | a web service; terms not stated by the plugin | any | not ported | an obsolete remote protocol queried per variant; absent from the release 116 web manifest |
| [dbNSFP](#dbnsfp) | the `_variant` download prepared per assembly as the Ensembl VEP plugin documents (`dbNSFP<version>_grch38.gz`, `_grch37.gz`), or a dbNSFP 5.x file; Liu et al. | free for academic use; commercial use needs a paid agreement | GRCh37, GRCh38 | supported on the user's licensed file, verified on synthetic fixtures | the file is not redistributed here, so the comparison with Ensembl VEP runs on a synthetic file named for dbNSFP 4.9a, in both preparations the plugin's documentation gives (`dbNSFP.pm` lines 68-74) |
| [dbscSNV](plugins/annotation-store.md#dbscsnv) | `dbscSNV1.1_GRCh37.txt.gz` and `dbscSNV1.1_GRCh38.txt.gz`, the download prepared per assembly as the Ensembl VEP plugin documents; dbscSNV 1.1, Jian, Boerwinkle and Liu | no licence text published; commercial use not stated | GRCh37, GRCh38 | supported on the user's licensed file, verified on synthetic fixtures | with no licence text published the file is not redistributed here, so the comparison with Ensembl VEP runs on a synthetic file in both preparations the plugin's documentation gives (`dbscSNV.pm` lines 53-59) |
| DeNovo | a PED pedigree file and the record's genotypes | the user's own data | GRCh37, GRCh38 | not ported | a filter plugin that writes per-run report files from a pedigree and the record's genotypes; a vep-rs plugin drops no output line and writes no file |
| [DosageSensitivity](plugins/flat-file.md#dosagesensitivity) | `Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz`, Zenodo record 6347673; Collins et al., 2022 | CC BY 4.0; commercial use permitted | GRCh37, GRCh38 | supported, verified on real data | open data keyed on gene symbol, so one file serves both assemblies; `file=` and `cover=0\|1`; `pHaplo` and `pTriplo` on every transcript row, or under `cover=1` only where the variant spans the transcript, structural variants included; Ensembl VEP's module dies on a regulatory or motif row under `--regulatory` and writes a warning per row kind; vep-rs writes none ([the class](plugins/flat-file.md#2-dosagesensitivity-a-warning-per-row-kind-on-regulatory-and-motif-rows)) |
| [Downstream](plugins/transcript-model.md#downstream) | none; the cache's coding sequence, 3' UTR and protein | Apache 2.0 plugin code; no data | GRCh37, GRCh38 | supported, verified on real data | the changed coding sequence translated from the frameshift on the record's own span, with and without `--hgvs`, where Ensembl VEP's plugin reads the HGVS-shifted span ([the class](divergences/plugin-hgvs-shifted-cds-coordinates.md)) |
| Draw | none; the GD graphics library | no data | any | not ported | writes a PNG image per output row through the GD library |
| Enformer | `Enformer_grch37.vcf.gz` and `Enformer_grch38.vcf.gz`, the Enformer variant scores lifted over; Ensembl FTP | redistributed by Ensembl; the score provider's terms not stated by the plugin | GRCh37, GRCh38 | not ported | over five thousand regulatory expression-impact features per variant from a lifted-over score file; rarely used |
| EVE | `eve_merged.vcf.gz`, built by the user from the evemodel.org bulk VCFs, and the popEVE file `grch38_popEVE_ukbb_20250715.vcf.gz` | MIT code; no data terms stated | GRCh38 | not ported | missense-only evolutionary scores, popEVE included since release 116, distributed without stated data terms; AlphaMissense covers missense pathogenicity |
| FATHMM | the `fathmm.py` script against a MySQL load of `fathmm.v2.3.SQL` | terms not stated by the plugin | any | not ported | runs an external Python script against a MySQL database |
| FATHMM_MKL | `fathmm-MKL_Current.tab.gz` | terms not stated by the plugin | GRCh37 | not ported | GRCh37-only coding and non-coding scores |
| FlagLRG | `list_LRGs_transcripts_xrefs.txt`; EMBL-EBI, the Locus Reference Genomic resource | EMBL-EBI terms | any; keyed on transcript id | not ported | flags Locus Reference Genomic transcripts; the LRG resource is retired and the `--lrg` flag is refused |
| FunMotifs | tissue BED files served from the FunMotifs data server | terms not stated by the plugin | GRCh37 | not ported | tissue BED files served from a preprint's data server, GRCh37 only; rarely used |
| G2P | a G2P, PanelApp or GenCC CSV, optional VCF frequency collections and an include-list VCF; EMBL-EBI | EMBL-EBI terms | GRCh37, GRCh38 | not ported | diagnostic allelic-requirement filtering that aggregates per individual, reads VCF frequency collections and writes HTML and text reports |
| GeneBe | the `genebe.net` web service | academic use only; rate-limited | GRCh38 | not ported | queries the GeneBe web service per variant under academic-only, rate-limited terms |
| GeneSplicer | the GeneSplicer binary and its training directory, Johns Hopkins University, plus `--fasta` | the GeneSplicer package's own terms | GRCh37, GRCh38 | not ported | runs the external GeneSplicer binary per variant; MaxEntScan covers splice-site scoring |
| Geno2MP | the Geno2MP VCF; University of Washington Center for Mendelian Genomics | the Geno2MP data-use terms, not stated by the plugin | GRCh37 | not ported | rare-variant genotype counts linked to HPO terms in GRCh37 coordinates; rarely used |
| [gnomADc](plugins/annotation-store.md#gnomadc) | the gnomAD coverage summary files, v2.1 on GRCh37 and v4 on GRCh38; Broad Institute | CC0; commercial use permitted | GRCh37, GRCh38 | supported, verified on real data | open data; the `gnomAD_<key>_cov` fields of Ensembl VEP's plugin, their keys chosen by the file's column count |
| [gnomADMt](plugins/annotation-store.md#gnomadmt) | `gnomad.genomes.v3.1.sites.chrM.vcf.bgz`, the gnomAD v3.1 mitochondrial sites VCF; Broad Institute | CC0, the gnomAD terms; commercial use permitted | GRCh38 | supported, verified on real data | open data; position and allele on MT, intergenic rows included; the fields are the file's INFO fields, the five haplogroup fields by default, any under `fields=`, reordered by `hap_filter=` and `pop_filter=`; Ensembl VEP's header cuts five of the file's `##INFO` descriptions at their first `=`; vep-rs writes the file's text ([the class](plugins/annotation-store-divergences.md#6-gnomadmt-the-info-descriptions-of-the-files-fields)) |
| GO | a Gene Ontology GFF generated from the Ensembl core database on first use, or a pre-built copy | Ensembl data | GRCh37, GRCh38 | not ported | the gene ontology file is built from the Ensembl core database on first use; vep-rs connects to no database |
| [GWAS](plugins/flat-file.md#gwas) | the NHGRI-EBI GWAS Catalog associations download, released continuously, and a study's harmonised summary-statistics file; EMBL-EBI | EMBL-EBI terms of use; commercial use permitted | GRCh38 | supported, verified on real data | open data; the catalog publishes GRCh38 coordinates only; both file types read under the `key=value` grammar of Ensembl VEP's plugin; Ensembl VEP hands every location of a haplotype or interaction row to every SNP of the row and vep-rs pairs each SNP with its own location ([the class](plugins/flat-file.md#1-gwas-every-snp-of-a-haplotype-row-at-every-location-of-the-row)) |
| Gwava | GRCh37 BED score files from the Sanger Institute FTP; the module sits outside the Ensembl VEP plugin collection | retired data source | GRCh37 | not ported | an external plugin outside the Ensembl VEP collection whose GRCh37 data source is retired |
| GXA | the Gene Expression Atlas web service, queried per gene | a web service; terms not stated by the plugin | any | not ported | queries the Gene Expression Atlas web service per gene; absent from the release 116 web manifest |
| [HGVSIntronOffset](plugins/transcript-model.md#hgvsintronoffset) | none; the HGVSc coordinates of `--hgvs` | Apache 2.0 plugin code; no data | GRCh37, GRCh38 | supported, verified on real data | the two intron offsets of the HGVSc coordinates; a run without `--hgvs` is refused |
| HGVSReferenceBase | none; `--hgvs` required | no data | GRCh37, GRCh38 | not ported | the reference base of the long HGVS deletion or insertion form, computed from the transcript model; absent from the release 116 web manifest |
| IntAct | `mutation_gc_map.txt.gz` and `mutations.tsv`; EMBL-EBI IntAct | EMBL-EBI terms | GRCh37, GRCh38 | not ported | molecular-interaction effects from the IntAct mutation tables across seven species; rarely used |
| LD | the Ensembl variation database plus 1000 Genomes genotype VCFs | Ensembl data | GRCh37, GRCh38 | not ported | linkage disequilibrium computed through the Ensembl variation database and genotype VCFs; vep-rs connects to no database |
| LocalID | a SQLite database built from the Ensembl variation database or the cache | Ensembl data | any | not ported | serves identifier input from a SQLite database built from the variation database; vep-rs reads identifier input from the variation cache with `--format id` |
| [LOEUF](plugins/annotation-store.md#loeuf) | the gnomAD constraint tables tabix-indexed on their gene coordinate columns as the Ensembl VEP plugin documents: v2.1.1 by transcript on GRCh37, v4.1.1 on GRCh38 as Ensembl publishes it with gene coordinates appended; Broad Institute | CC0, the gnomAD terms; commercial use permitted | GRCh37, GRCh38 | supported, verified on real data | open data; `file=` and `match_by=gene\|transcript`, both required; `LOEUF` on transcript rows, the minimum over the gene's transcripts under `match_by=gene`, none on upstream or downstream rows; Ensembl VEP's module drops a gene's minimum written `0`; vep-rs prints it ([the class](plugins/annotation-store-divergences.md#7-loeuf-a-genes-minimum-of-0)) |
| [LoFTEE](plugins/loftee.md) | `human_ancestor.fa.gz` and `gerp_conservation_scores.homo_sapiens.GRCh38.bw` of the LoFTEE `grch38` data bundle; Ensembl ancestral alleles and Compara GERP scores | MIT code on Ensembl data; commercial use permitted | GRCh38 | supported, verified on real data | loads as `--plugin LoF` or `--plugin LoFTEE`; the GRCh37 bundle is read, but its GERP table carries no licence text, so GRCh37 verdicts are unverified; a `conservation_file` naming a database is refused, since vep-rs does not compute the PhyloCSF items |
| [LoFtool](plugins/flat-file.md#loftool) | `LoFtool_scores.txt` beside the Ensembl VEP plugin; Fadista et al. | Apache 2.0; commercial use permitted | GRCh37, GRCh38 | supported, verified on real data | open data keyed on gene symbol, so one file serves both assemblies; given no file, read from the `--dir_plugins` directory as Ensembl VEP's plugin reads its bundled copy |
| LOVD | the `lovd.nl` web service, queried per variant | a web service; terms not stated by the plugin | GRCh37, GRCh38 | not ported | queries the LOVD web service per variant, which the plugin itself asks users not to do intensively |
| Mastermind | `mastermind_cited_variants_reference-*.vcf.gz`; Genomenon | Genomenon's terms, behind a free registration | GRCh37, GRCh38 | not ported | literature citation counts and identifiers from a vendor file behind registration, under commercial terms |
| MaveDB | `MaveDB_variants.tsv.gz`; Ensembl FTP | redistributed by Ensembl; each MaveDB dataset carries its own licence | GRCh38 | not ported | deep mutational scanning scores whose licence differs per dataset; rarely used |
| [MaxEntScan](plugins/transcript-model.md#maxentscan) | the directory of matrix files the Ensembl VEP plugin opens (`MaxEntScan.pm` lines 104-105 and 824-831), read as a positional parameter | no licence text published; commercial use not stated | GRCh37, GRCh38 | supported on the user's licensed file, verified on synthetic fixtures | with no licence text published the matrices are not redistributed here, so the comparison with Ensembl VEP runs on a synthetic matrix set under the file names the plugin opens (`MaxEntScan.pm` lines 104-105 and 824-831); splice-site windows read from `--fasta`; `SWA`, `NCSS` and `verbose` as Ensembl VEP's plugin reads them; a substitution in an exon of at most five bases lies in two canonical windows, which Ensembl VEP resolves in the order its interval tree returns the introns and vep-rs in transcript order, the upstream intron's acceptor window first ([the class](plugins/transcript-model.md#2-maxentscan-two-splice-site-windows-over-one-substitution)) |
| MechPredict | `MechPredict_input.tsv`, assembled from three OSF datasets; Badonyi, 2024 | terms not stated by the plugin | any; keyed on gene | not ported | gene-level probabilities of dominant-negative, gain-of-function and loss-of-function mechanisms, assembled from three OSF datasets; rarely used |
| MPC | `fordist_constraint_official_mpc_values.txt.gz`; Broad Institute, the ExAC release | terms not stated by the plugin | GRCh37 | not ported | GRCh37-only regional missense constraint scores |
| MTR | `mtrflatfile_2.0`; University of Melbourne | terms not stated by the plugin | GRCh37 | not ported | a GRCh37-only missense tolerance ratio |
| mutfunc | `mutfunc_data.db`, a SQLite database, plus `--fasta`; Ensembl FTP | Ensembl data | any; keyed on protein position | not ported | protein destabilisation and motif effects from a SQLite database; vep-rs embeds no SQLite |
| [NearestExonJB](plugins/transcript-model.md#nearestexonjb) | none; the cache's exon model | Apache 2.0 plugin code; no data | GRCh37, GRCh38 | supported, verified on real data | the nearest junction of each overlapped exon, `max_range=` and `intronic=1`; items in transcript order where Ensembl VEP's follow hash order; `gff3=<file>` is refused |
| NearestGene | the Ensembl core database by default; `gff3=` names a tabix-indexed GFF3 of gene locations | Ensembl data | GRCh37, GRCh38 | not ported | reads the Ensembl core database by default; the `--nearest` flag covers the nearest-gene lookup from the cache |
| neXtProt | the neXtProt SPARQL endpoint, queried per variant, plus a database connection | a web service; terms not stated by the plugin | any | not ported | queries the neXtProt SPARQL endpoint per variant and needs a database connection; absent from the release 116 web manifest |
| NMD | none; the transcript's exon model and CDS coordinates | Apache 2.0 plugin code; no data | GRCh37, GRCh38 | supported, verified on real data | no data file; `NMD_escaping_variant` on stop-gained, frameshift and splice-site rows that escape nonsense-mediated decay under Ensembl VEP's four exon-position rules; takes no parameter; under `--minimal` a split allele is judged over its own trimmed span, where Ensembl VEP's plugin reads the rejoined record's ([the class](divergences/nmd-minimal-rejoined-end.md)) |
| NonSynonymousFilter | none | no data | any | not ported | an example filter plugin absent from the release 116 web manifest; a vep-rs plugin drops no output line, and the `vep filter` tool covers the use |
| [OpenTargets](plugins/annotation-store.md#opentargets) | `open_targets_vep_full.tsv.gz` with its index, the Open Targets Platform 26.03 VEP export from the EMBL-EBI FTP, its header line commented and re-indexed | CC0 1.0; commercial use permitted | GRCh38 | supported, verified on real data | open data; locus-to-gene scores for GWAS and QTL loci with `cols=`, `study_types=` and `lead_variants_only=` as Ensembl VEP's plugin reads them |
| Paralogues | a `matches=` table from the Ensembl FTP, the Compara database or a web service | Ensembl data | GRCh38 (the matches table) | not ported | variants at paralogous aligned positions, read from the Compara database, a web service or a GRCh38 matches file; rarely used |
| PhenotypeOrthologous | `PhenotypesOrthologous_homo_sapiens_<release>_GRCh38.gff3.gz`; Ensembl FTP | Ensembl data | GRCh38 | not ported | mouse and rat orthologue phenotypes from a GRCh38 GFF3; rarely used |
| Phenotypes | a phenotype GVF generated from the Ensembl variation database on first use | Ensembl data | GRCh37, GRCh38 | not ported | the phenotype file is generated from the Ensembl variation database on first use; vep-rs connects to no database |
| [pLI](plugins/flat-file.md#pli) | `pLI_values.txt` beside the Ensembl VEP plugin, cut from the ExAC r0.3 constraint table; Samocha and MacArthur | Apache 2.0; commercial use permitted | GRCh37, GRCh38 | supported, verified on real data | open data keyed on gene symbol, or on transcript stable id in `transcript` mode, so one file serves both assemblies; given no file, read from the `--dir_plugins` directory as Ensembl VEP's plugin reads its bundled copy |
| PolyPhen_SIFT | a SQLite database of predictions built from the Ensembl database, or the pangenome databases on the Ensembl FTP | Ensembl data | GRCh37, GRCh38 | not ported | SIFT and PolyPhen predictions from a SQLite database built from the Ensembl database; the vep-rs cache carries both natively |
| PON_P2 | its web service through `ponp2.py`, or a pre-computed COSMIC-derived GRCh37 file | COSMIC terms, restricted, for the file mode | GRCh37 (the file mode) | not ported | queries its web service per variant, or reads a pre-computed GRCh37 file derived from COSMIC under its restricted terms |
| PostGAP | none; the module was removed from the Ensembl VEP plugin collection in release 116 | none | none | not ported | removed from the Ensembl VEP plugin collection in release 116 |
| PrimateAI | `PrimateAI_scores_v0.2.tsv.gz` (GRCh37) and `PrimateAI_scores_v0.2_hg38.tsv.gz` (GRCh38), sorted, bgzip-compressed and tabix-indexed as the Ensembl VEP plugin documents (`PrimateAI.pm` lines 65-69); Illumina | non-commercial use under Illumina's terms | GRCh37, GRCh38 | supported on the user's licensed file, verified on synthetic fixtures | the file is not redistributed here, so the comparison with Ensembl VEP runs on a synthetic file with the header columns the plugin requires (`PrimateAI.pm` lines 136-138); the file is the first positional parameter, and a file name carrying an assembly must match the run; single-base substitutions on every row of the record, intergenic rows included |
| PromoterAI | `promoterAI.tsv.bgz`; Illumina | Illumina licence agreement, restricted | GRCh38 | not ported | promoter expression-impact scores under a commercial licence agreement |
| ProteinSeqs | none; writes two FASTA files | no data | GRCh37, GRCh38 | not ported | writes two FASTA files of reference and mutated protein sequences as a side effect |
| ProtVar | `ProtVar_data.db`, a SQLite database, plus `--fasta`; EMBL-EBI | EMBL-EBI terms | any; keyed on protein position | not ported | structure stability, pocket and interface context from a SQLite database; vep-rs embeds no SQLite |
| RankFilter | none | no data | any | not ported | an example filter plugin absent from the release 116 web manifest; `--most_severe`, `--pick` and the `vep filter` tool cover the use |
| ReferenceQuality | GRC issue and clone-problem GFF3 files merged by the user; NCBI | NCBI, public domain | GRCh37, GRCh38 | not ported | reference-assembly issue annotations from GRC GFF3 files merged by the user; rarely used |
| RefSeqHGVS | the Ensembl otherfeatures database | Ensembl data | any | not ported | reads the Ensembl otherfeatures database; `--refseq` with `--hgvs` yields RefSeq HGVS notation directly, and the plugin is absent from the release 116 web manifest |
| [REVEL](plugins/annotation-store.md#revel) | `revel-v1.3_all_chromosomes.zip`, tabbed, its header line commented, bgzipped and indexed per assembly; REVEL, v1.3 with transcript ids | Open Database License 1.0; commercial use permitted | GRCh37, GRCh38 | supported, verified on real data | open data; one file carries both assemblies' coordinates, indexed per assembly; the field is `REVEL`, written on missense rows whose reference and alternate amino acids and, unless `no_match`, transcript the record names |
| RiboseqORFs | `Ribo-seq_ORFs.bed.gz` plus `--fasta`; GENCODE | GENCODE, EMBL-EBI terms | GRCh38 | not ported | recomputes consequences over Ribo-seq open reading frame models; absent from the release 116 web manifest |
| SameCodon | the Ensembl variation database | Ensembl data | any | not ported | reads the Ensembl variation database for other variants in the same codon |
| satMutMPRA | `satMutMPRA_GRCh37_ALL` and `satMutMPRA_GRCh38_ALL`; University of Washington | terms not stated by the plugin | GRCh37, GRCh38 | not ported | saturation mutagenesis MPRA effects for twenty-one regulatory elements; rarely used |
| SingleLetterAA | none | no data | GRCh37, GRCh38 | not ported | rewrites `HGVSp` with one-letter amino-acid codes, a form of the protein notation vep-rs does not compute |
| [SpliceAI](plugins/spliceai.md) | `spliceai_scores.raw.snv.hg38.vcf.gz`, `spliceai_scores.raw.indel.hg38.vcf.gz` and their hg19 twins; Illumina, `genome_scores_v1.3` (`SpliceAI.pm` line 109); the Ensembl MANE SNV file `spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz` | academic use only under Illumina's terms; no licence text published for the Ensembl-derived SNV file | GRCh37, GRCh38 | supported on the user's licensed file, verified on synthetic fixtures | the files are not redistributed here, so the comparison with Ensembl VEP runs on synthetic SNV and indel files in the VCF layout the plugin reads (`SpliceAI.pm` lines 257 and 326-332), in every output format and under `cutoff` and `split_output` |
| SpliceRegion | none | no data | GRCh37, GRCh38 | not ported | its `splice_donor_5th_base_variant`, `splice_donor_region_variant` and `splice_polypyrimidine_tract_variant` terms are core consequence terms of both engines; only its `Extended` mode adds a key of its own |
| [SpliceVault](plugins/annotation-store.md#splicevault) | `SpliceVault_data_GRCh38.tsv.gz` and `SpliceVault_data_hg19.tsv.gz`; Dawes et al., 2023 | no licence text published; commercial use not stated | GRCh37, GRCh38 | supported on the user's licensed file, verified on synthetic fixtures | with no licence text published the tables are not redistributed here, so the comparison with Ensembl VEP runs on a synthetic file of the twelve columns the plugin reads (`SpliceVault.pm` line 201); `file=`; seven `SpliceVault_` fields on the rows naming the transcript under `transcript_match=1`, and on every transcript row from the first allele match under `transcript_match=0`; on a multi-allelic record Ensembl VEP's module gives every allele the first allele's table row; vep-rs matches each allele's own ([the class](plugins/annotation-store-divergences.md#3-splicevault-the-table-row-of-a-sibling-allele)); in the JSON format under `transcript_match=0` Ensembl VEP's module writes the `splicevault` object on a record's first transcript row alone; vep-rs writes it on every row ([the class](plugins/annotation-store-divergences.md#4-splicevault-the-json-object-of-a-records-later-transcript-rows)) |
| StructuralVariantOverlap | none; the module was removed from the Ensembl VEP plugin collection in release 116 | none | none | not ported | removed from the Ensembl VEP plugin collection in release 116; a `--custom` structural-variant VCF source covers the use |
| SubsetVCF | any tabix-indexed VCF named by the user | the user's own file | GRCh37, GRCh38 | not ported | all alternate alleles and INFO fields of the overlapping records of a VCF; a `--custom` VCF source covers the common case, and the plugin is absent from the release 116 web manifest |
| [TranscriptAnnotator](plugins/annotation-store.md#transcriptannotator) | a tabix-indexed table of your own: chromosome, position, reference and alternate allele, transcript id, then the columns to append, named by its last header line | the user's own file; no third-party terms | GRCh37, GRCh38 | supported, verified on real data | the table is the user's, so the comparison with Ensembl VEP runs on a table of the documented layout; `file=`, `cols=`, `prefix=` and `trim=`; one `<prefix>_<column>` field per selected column on the rows whose transcript id or RefSeq cross-reference the table names; on a multi-allelic record Ensembl VEP's module gives an allele without a row a sibling allele's; vep-rs matches each allele's own ([the class](plugins/annotation-store-divergences.md#5-transcriptannotator-the-table-row-of-a-sibling-allele)) |
| [TSSDistance](plugins/transcript-model.md#tssdistance) | none; the transcript's span | Apache 2.0 plugin code; no data | GRCh37, GRCh38 | supported, verified on real data | the distance to the transcription start site, upstream rows by default and every row under `both_direction=1`, structural variants included |
| [UTRAnnotator](plugins/transcript-model.md#utrannotator) | `uORF_5UTR_GRCh37_PUBLIC.txt` and `uORF_5UTR_GRCh38_PUBLIC.txt` of the Ensembl UTRannotator repository; Whiffin et al., 2020 | MIT, the repository's licence; commercial use permitted | GRCh37, GRCh38 | supported, verified on real data | open data; a 5' UTR analysis over the exon model and `--fasta`, which the plugin needs; `file=` names the uORF evidence table, without which every `Evidence` is `NA`; `max_overlap=` defaults to 100; Ensembl VEP's gained-start window is sized by the reference allele, so an insertion's gained `ATG` can be missed; vep-rs sizes it by the alternate ([the class](plugins/transcript-model.md#1-utrannotator-the-window-of-the-gained-start-test)) |
| VARITY | `varity_all_predictions`; varianteffect.org | terms not stated by the plugin | GRCh37 | not ported | a GRCh37-only rare-missense predictor |

## Rows a plugin writes on

A plugin runs on the rows whose kind its Ensembl VEP module's `feature_types`
names, as `OutputFactory.pm` `run_plugins` checks it (line 2223): a transcript
row, a regulatory-feature row, a motif-feature row or the intergenic row of a
record that overlaps no transcript. CADD, gnomADc, GWAS, BayesDel, PrimateAI,
gnomADMt, AncestralAllele, AVI, OpenTargets and Conservation declare every kind.
Their fields therefore print on every row of a record, the regulatory, motif and
intergenic rows of a `--regulatory` run included. ClinPred, LOEUF, NMD,
TranscriptAnnotator and UTRAnnotator declare `Transcript`, and
DosageSensitivity's `Feature` reaches no row without a gene symbol.

REVEL, AlphaMissense, dbscSNV, SpliceAI, dbNSFP, LoFTEE, LoFtool and pLI declare
`Transcript` too, so none of their fields prints on a regulatory, motif or
intergenic row. A shared-library plugin names its kinds in
`vep_plugin_feature_types`, `Transcript` alone when it exports none. One naming
`RegulatoryFeature` or `MotifFeature` switches regulatory annotation on
([Writing a plugin](developers/writing-a-plugin.md)).

A plugin also runs only on the variant kinds its module's `variant_feature_types`
names, checked before the row kind (`OutputFactory.pm` line 2220). The default is
sequence variants alone (`BaseVepPlugin.pm` line 100). Of the built-in plugins,
CADD, DosageSensitivity and UTRAnnotator declare `StructuralVariationFeature` as
well and TSSDistance declares `BaseVariationFeature`. These four annotate the
rows of a structural variant and every other built-in writes nothing on them, as
in Ensembl VEP. A shared-library plugin sees both kinds and reads
`is_structural` in the variant document when it wants one.

## dbNSFP {#dbnsfp}

`--plugin dbNSFP` reads the user's own tabix-indexed dbNSFP file; vep-rs ships
no data and compares with Ensembl VEP on a synthetic file prepared as
`dbNSFP.pm` lines 62-74 document. Native SIFT and PolyPhen decoding needs
no dbNSFP file ([SIFT and PolyPhen predictions](cache-predictions.md)). The
parameters are those of `dbNSFP.pm`, `key=value` options anywhere in the list.

| Parameter | Meaning |
| --- | --- |
| first positional | the bgzip-compressed, tabix-indexed dbNSFP file |
| second positional, when it names an existing file | a replacement-logic file (`colname from to` lines with `default` rows) in place of the built-in `;` to `,` and `\|` to `&` |
| remaining positionals | the columns to write under their own names; `ALL` writes every column; a name the header lacks prints `invalid_field` |
| `consequence=` | `ALL`, or SO terms joined by `&`; the default annotates rows whose consequences include `missense_variant`, `stop_lost`, `stop_gained` or `start_lost` |
| `output_prefix=1` | prefixes every output name with `dbNSFP_` |
| `pep_match=0` | matches a record on position and allele alone; by default its `aaref` and `aaalt` must equal the row's amino acid change |
| `transcript_match=1` | keeps, in a `;`-joined transcript-specific column, the entry of the row's transcript alone |

Two `--plugin` options as they stand on a `vep` command line:

```bash
--plugin dbNSFP,/path/to/dbNSFP4.9a_grch38.gz,\
SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred
--plugin dbNSFP,/path/to/dbNSFP4.9a_grch38.gz,\
transcript_match=1,output_prefix=1,SIFT_score,GERP++_RS
```

### Matching and values

A record matches a single-base substitution (`A`, `C`, `G` or `T`) on its
position and alternate allele. The first matching table row of a position is the
row's, and a mitochondrial record is read under the name `M` (`dbNSFP.pm` line
313). A `.` value is
not written, nor the transcript's own `.` list entry under `transcript_match=1`.
An empty cell the row has is written as the empty value.

`;` and `|` in a value become `,` and `&` unless the
replacement-logic file says otherwise for that column. A logic line without
three blank-separated values, one opening with a blank included, fails the
plugin's load as in Ensembl VEP.
A `dbNSFP_replacement_logic` file in the working directory is not read;
Ensembl VEP 116.2 reads one when no second positional names a file
(`dbNSFP.pm` lines 548-549).

### Core fields and warnings

A selected column named as a core field, `TSL` or `APPRIS`, prints
the file's value on the plugin's rows over the cache's, with or without `--tsl`
and `--appris`, as Ensembl VEP does (`OutputFactory.pm` line 2231). Other rows
print the cache's value under the flag alone. Under `--output_format
parquet` the INT32 `TSL` column stores a list of support levels as NULL
([Parquet output](output/parquet.md)).

The `--refseq` or `--merged` warning under `transcript_match=1` and the warning
naming a selected column the header lacks go to standard error alone, under
`--quiet` too, as Ensembl VEP prints them. A `;`-joined column whose entries do
not pair with the row's transcript ids keeps its whole value, reported on
standard error for every such row. Under `--fork` Ensembl VEP also copies the
report into the warnings file (`Runner.pm` lines 565 and 659-662), vep-rs does
not.

### Preparing the dbNSFP file

**Preparing the data.** Prepare the `_variant` download as the Ensembl VEP
plugin documents (`dbNSFP.pm` lines 68-74): for GRCh38 sorted and indexed on `#chr` and `pos(1-based)`
(`tabix -s 1 -b 2 -e 2`), for GRCh37 on column 8 and `hg19_pos(1-based)`
(`tabix -s 8 -b 9 -e 9`). The two share the header and differ in their index.
vep-rs therefore resolves the position column by name for the run's assembly and
refuses a file indexed on the other column, naming both. A dbNSFP 2.9 file
(`pos(1-coor)`) serves a GRCh37 run (`dbNSFP.pm` lines 317-319), a file named
as 4.0b1 a GRCh38 one (lines 329-330).

Transcript-specific columns are detected in `1:1008170-1082927` under the bare
name `1`, as Ensembl VEP queries it. `transcript_match=1` is refused on a file
naming its contigs `chr1`, whose rows are annotated without it. The README is
any entry beside the data file whose name carries `dbNSFP` and then
`readme.txt` in any letter case, a `.bak` copy included. It supplies the header
descriptions and marks the transcript-specific columns its text calls
"separated by" or ties to an "Ensembl gene id". Of several, Ensembl VEP reads
the first in directory order, vep-rs the first in name order.

### File name and README descriptions

The data file's name must carry a dbNSFP version, read as Ensembl VEP's plugin
reads it (`dbNSFP.pm` lines 197-211). A name carrying none fails the plugin's
load, `Failed to instantiate plugin dbNSFP` in the warnings file or a stop under
`--safe`. A bare `4.9` loads with the suffix warning on standard error. The
version decides acceptance alone; the layout comes from the header.

A description is `(from <data file>)`, the name without `.gz`, a blank and the
README's text. ASCII whitespace is collapsed to single spaces, none at either end.
Non-UTF-8 README or logic-file bytes print as the replacement character where
Ensembl VEP prints them. Ensembl VEP 116.2 keeps one trailing space when the
text's last line ends in whitespace. It prints `<column> from dbNSFP file`, its
text for an undescribed column, when no space follows the colon (`dbNSFP.pm`
lines 520, 537 and 278).

A `0` name or text is kept. Ensembl VEP drops it (line 523), printing
`SIFT_score from dbNSFP file` for `SIFT_score: 0` where vep-rs prints
`(from <data file>) 0`. No committed corpus carries a dbNSFP README (tier 3 of
[Intended divergences from Ensembl VEP](intended-divergences.md); 14 flag
sets run the plugin).

A numbered line opens a column's description when a name, a colon, a blank and
one more character follow its number, the name ending at the first such colon,
or when non-blank text follows the colon. Any other numbered line describes no
column; its continuation lines are dropped.

## Where the tabix-backed plugins differ from Ensembl VEP

Beside those documented with each plugin, seven behaviours of Ensembl VEP
116.2's tabix-backed plugins are not reproduced. Four are set out here: two of
ClinPred and BayesDel, which have no page of their own, and two several
plugins share.

Each is a documented
class: the section names the rows a comparison would set aside. No
comparison tool codes the class, since no committed row falls under it. The
index's
[How to read the counts](intended-divergences.md#how-to-read-the-counts)
defines the format comparator and the comparison case.

The concordance figures
score consequence sets, not plugin fields, so nothing is set aside in them. Rows
below are in the default output format.

The fifth is
[a mitochondrial record against a plugin's file naming the contig under another spelling](divergences/plugin-contig-mitochondrial-spelling.md).
The sixth is
[no annotation from seven tabix-backed plugins at the start of a chromosome](divergences/tabix-window-below-position-one.md).
The seventh is
[a VCF table row with `END` below its position](divergences/tabix-row-with-an-end-below-its-start.md).

Each class carries its tier under the index's rule. Tier 1 is a shape the
chromosome 21 population cuts carry with no flag, tier 2 one that needs a flag
or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts
are over the cuts' 15,646,481 records, the 111 golden corpora and their 392 flag
sets.

### 1. ClinPred: a record given on the reverse strand

**The defect in one sentence.** For an Ensembl-format or region record given on
the reverse strand, Ensembl VEP's ClinPred plugin compares the allele as the
input wrote it with the alternate-allele column of its score table, so the
row carries another substitution's score, or none.

**Why this is a defect.** The synthetic tables of the committed corpora list
each position's substitutions on the chromosome's forward strand. An
allele written on the reverse strand names the same substitution as its
complement on the forward strand. Ensembl VEP's own allele matcher
reverse-complements a variant's alleles before comparing them with a file whose
strand differs.

The plugin skips that step. A reverse-strand `A/G` is thus read as
a forward-strand `G` and scored from the table's `T G` row, where the record
describes `T>C`.

**What vep-rs does instead.** vep-rs compares the record's forward-strand
alleles with the table, so a reverse-strand record scores the row its
forward-strand twin scores.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116,
`ensembl-variation` at release/116):

- `ClinPred.pm`, line 109, `run`: the allele is read through
  `variation_feature_seq`, the sequence relative to the variation feature, on
  the record's own strand. Line 117: it is compared with a table row's `alt` by
  string equality, with no complement.
- `Bio/EnsEMBL/Variation/VariationFeatureOverlapAllele.pm`, lines 276-281:
  `variation_feature_seq` returns the allele as written; lines 246-264:
  `feature_seq` alone reverse-complements, and the plugin never calls it.
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 1181-1182 and 1195-1196,
  `get_matched_variant_alleles`: the engine's matcher reverse-complements a
  variant's alleles when its strand differs from the file's.

**What correcting it would change.** A reverse-strand record would carry the
score of the substitution it describes, the value its forward-strand twin gets
from both engines.

**The class it produces.** `reverse_strand_record_as_written_allele`: a
`ClinPred` field (the JSON `clinpred`) whose two values differ, or that one
engine alone fills, on a row of a record given on the reverse strand in Ensembl
or region format. A VCF record is always on the forward strand, so the shape has
no row on the committed corpora or the comparison cases.

**How often the shape occurs.** Tier 3: the shape is an Ensembl-format line on
the reverse strand, which no population cut uses (VCF is forward-strand). 107 of
the 217 Ensembl-format lines of the committed corpora are reverse-strand, in two
corpora, and 4 flag sets run the plugin; the table's rows at the position decide
which row is wrongly taken.

#### 1.1 A reverse-strand missense record in KRTAP24-1

**In short.** `21 30282172 30282172 A/G - rev_missense` is the reverse-strand
spelling of the substitution `T>C` at `21:30282172`, which the synthetic score
table of the GRCh38 ClinPred corpus
(`tests/golden/116/GRCh38-plugin-clinpred-fixture/`) scores `0.163359` (`T C`),
beside `0.507213` (`T A`) and `0.375971` (`T G`). Ensembl VEP writes
`ClinPred=0.375971`, the `T G` row, because it matched the written `G`; vep-rs
writes `ClinPred=0.163359`. The forward-strand twin
`21 30282172 30282172 T/C + fwd_missense` and the VCF record
`21 30282172 . T C` get `0.163359` from both engines, and for each record the
two engines' rows agree in every other column, the codon change `tAt/tGt` on
`ENST00000340345` included.

### 2. BayesDel: a record given on the reverse strand

**The defect in one sentence.** For an Ensembl-format or region record given on
the reverse strand, Ensembl VEP's BayesDel plugin calls a `reverse_comp` helper
its module never imports, so the call dies with a `Plugin 'BayesDel' went
wrong` warning and the record's rows carry no `BayesDel` field.

**Why this is a defect.** The plugin's own lines state the rule. When the
record's strand is -1, reverse-complement its reference and alternate alleles,
then match them against the table's forward-strand columns. The helper they
name lives in `Bio::EnsEMBL::Utils::Sequence`, which sibling plugins such as
`SpliceAI.pm` import and `BayesDel.pm` does not. Perl therefore finds no
`BayesDel::reverse_comp` and the match never runs.

The plugin runner catches
the death, prints one warning per row and leaves the row without the field,
exit 0.

**What vep-rs does instead.** vep-rs complements the alleles of a reverse-strand
record and matches the table on the forward strand, the rule the lines state.
The record therefore scores the row its forward-strand twin scores.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116,
`ensembl-vep` at release/116.2):

- `BayesDel.pm`, lines 135-138, `run`: `reverse_comp(\$ref_allele)` and
  `reverse_comp(\$alt_allele)` under `if ($vf->{strand} == -1)`. Lines 74-79,
  the module's `use` lines, load `strict`, `warnings` and the tabix base class
  alone, and the file defines no `sub reverse_comp`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2225-2245, `run_plugins`: the
  plugin runs inside `eval`, and its death becomes the warning
  `Plugin 'BayesDel' went wrong: Undefined subroutine &BayesDel::reverse_comp
  called at ... line 136`.

**What correcting it would change.** One import line would make the plugin
score a reverse-strand record as vep-rs does. A VCF record is never on the
reverse strand, so no VCF run changes.

**The class it produces.** `bayesdel_reverse_strand_ensembl_input`: a
`BayesDel` field (the JSON `bayesdel`) the candidate alone fills on a row of a
record given on the reverse strand in Ensembl or region format, and in warnings
mode the reference's `Plugin 'BayesDel' went wrong` lines. No committed
corpus or comparison case runs the plugin on Ensembl-format or region input.

**How often the shape occurs.** Tier 3: the shape is an Ensembl-format line on
the reverse strand, which no population cut uses. 107 of the 217 Ensembl-format
lines of the committed corpora are reverse-strand, in two corpora, and 2 flag
sets run the plugin; the plugin dies before any lookup.

#### 2.1 A reverse-strand missense record in a forward-strand transcript

**In short.** `21 30396091 30396091 G/A - missense_fwd_minus` is the
reverse-strand spelling of `C>T` at `21:30396091`, the GRCh38 BayesDel corpus's
own record `missense_fwd`
(`tests/golden/116/GRCh38-plugin-bayesdel-fixture/`), whose synthetic table
holds `21 30396091 C T -1.18855`. Ensembl VEP writes the row
`missense_fwd_minus 21:30396091 A ... missense_variant ... IMPACT=MODERATE;STRAND=1`
with no `BayesDel`, and its warnings file names line 136. vep-rs ends the same
row with `BayesDel=-1.18855`. Both engines write that value for the corpus's
VCF record `missense_fwd`, `21 30396091 C T`, and Ensembl VEP writes it for the
forward-strand twin `21 30396091 30396091 C/T +`.

### 3. A score-table row whose alleles are written in lower case

**The defect in one sentence.** Ensembl VEP's AVI plugin, like every plugin that
matches a table row by allele, compares the table's alleles as written while the
engine has upper-cased the record's, so a table row whose alleles are in lower
case matches no record.

**Why this is a defect.** Letter case carries no meaning in a nucleotide allele.
Ensembl VEP upper-cases every input allele before annotation, and the VCF
specification, which the SpliceAI and gnomADMt fixture tables follow, declares REF and
ALT bases case-insensitive.

The engine's matcher builds `ref_alt_pos` keys from
the strings it is given and compares them as text, so `c_t_5030633` never
equals `C_T_5030633`. The plugins hand it the table's cells unchanged. The
same substitution, written in two cases, is read as two.

**What vep-rs does instead.** vep-rs compares a table row's alleles with the
record's without regard to letter case, in the matcher its AVI, gnomADMt,
SpliceAI, SpliceVault and TranscriptAnnotator plugins share. A lower-case row therefore
scores the record it describes. The row's own text is printed as written.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2,
`ensembl-variation` at release/116, `VEP_plugins` at release/116):

- `Bio/EnsEMBL/VEP/Parser.pm`, line 554, `validate_vf`: the record's allele
  string is upper-cased (`tr/[a-z]/[A-Z]/`) before any plugin runs.
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 1207 and 1232,
  `get_matched_variant_alleles`: the `"$ref\_$alt\_$pos"` keys of the two
  variants are compared as strings, and nothing in lines 1144-1247 folds case.
- `AVI.pm`, lines 135-136 and 139-151, `run`: the record's alleles pass an
  upper-case `/^[ACGT]$/` gate and are matched through that function. Lines 161
  and 176-181, `parse_data`: the row's `ref` and `alt` are split from the line
  and stored as written. `gnomADMt.pm`, line 278, `SpliceVault.pm`, line 165,
  and `TranscriptAnnotator.pm`, line 224, call the same function on their rows.
  `SpliceAI.pm`, line 272, compares its row's alleles with `eq`.

**What correcting it would change.** A table row in lower case would score the
record it describes, as its upper-case twin does, in every plugin that matches a
row by allele.

**The class it produces.** `tabix_allele_case_insensitive_match`: a field of
the AVI, gnomADMt, SpliceAI or SpliceVault plugin (the JSON `avi_`, `gnomad_`,
`spliceai` and `splicevault` keys), or a TranscriptAnnotator field a manifest
names, that the candidate alone fills on a row whose table row carries a
lower-case allele. Every committed table carries upper-case alleles alone, the
AVI, SpliceAI and SpliceVault tables being synthetic fixtures: 0 of
225 AVI rows, 0 of 373 gnomADMt rows, 0 of 141 SpliceAI rows, 0 of 92 SpliceVault
rows and 0 of 18 TranscriptAnnotator rows. No committed corpus and no comparison
case reaches the shape.

**How often the shape occurs.** Tier 3: the shape is a table row whose allele
cell holds a lower-case base. None of the 849 rows of the ten corpus tables
has one, and 23 flag sets run one of the five plugins.

#### 3.1 The AVI row for 21:5030633 with its alleles lower-cased

**In short.** The GRCh38 AVI corpus's synthetic table
(`tests/golden/116/GRCh38-plugin-avi/`)
holds the row `chr21 5030633 C T 1.6133 3.465`, and for the record
`21 5030633 . C T` under `raw=1,links=1` both engines write
`AVI_PHRED=3.465;AVI_RAW=1.6133` and the Atlas link ending
`chr21:5030633:C>T&m=variant`. Against a copy of the table in which that row
reads `chr21 5030633 c t 1.6133 3.465`, re-indexed with `tabix -s 1 -b 2 -e 2`,
Ensembl VEP writes the intergenic row
`21_5030633_C/T 21:5030633 T ... IMPACT=MODIFIER` with no `AVI_` field. vep-rs
writes the three fields, its link ending `chr21:5030633:c>t&m=variant`, the
row's text as written.

### 4. A chromosome written with a capitalised `chr` prefix in a plugin's file lookup

**The defect in one sentence.** Ensembl VEP annotates a record written `Chr21`
from the cache, which holds `21`, but the base class every tabix-backed plugin
reads its file through maps that spelling one step only, so against a file
whose contigs are named `chr21` the plugin silently finds nothing.

**Why this is a defect.** The engine keeps the record because its name, with
the prefix removed in any letter case, is a cache chromosome. By the engine's own
rule, `Chr21`, `chr21` and `21` name one chromosome.

The plugin base class
applies that rule once, trying the stripped `21` against the file's names. It
never tries the `chr21` the file has, although adding `chr` to a bare name is
its own next rule. The record is thus annotated from the cache and skipped by
the plugin, with no warning. gnomADc's own strip of the prefix is case-sensitive
and leaves `Chr21` untouched, so it does not rescue the lookup.

**What vep-rs does instead.** vep-rs names the record's chromosome as the cache
names it before any lookup, `21` here. It queries a plugin's file under that
name and its `chr`-prefixed form, so the lookup reaches `chr21`. The output
keeps the input's spelling.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2,
`ensembl-variation` at release/116, `VEP_plugins` at release/116):

- `Bio/EnsEMBL/VEP/Parser.pm`, lines 788-793, `_have_chr`: a record is kept
  when its name with a leading `chr` removed in any case, or with `chr` added,
  is a cache chromosome; the record's own name is not rewritten.
- `Bio/EnsEMBL/Variation/Utils/BaseVepTabixPlugin.pm`, lines 693-708,
  `_get_source_chr_name`: a name the file lacks is tried with its `chr` prefix
  removed in any case (lines 699-704) or, only when it has no such prefix, with
  `chr` added (lines 705-707). `Chr21` yields `21`, which a file with bare names
  holds and a `chr21` file lacks, so the query keeps `Chr21`.
- `gnomADc.pm`, line 158, `run`: `s/^chr//` strips a lower-case prefix alone
  before the base class is called.

**What correcting it would change.** A `Chr21` record would carry the plugin
fields its `21` and `chr21` spellings carry, in every tabix-backed plugin run
against a `chr`-prefixed file.

**The class it produces.** `plugin_contig_case_prefix`: a field of a
tabix-backed plugin that the candidate alone fills on a row whose chromosome is
written with a `chr` prefix not in lower case. No committed input carries such a
spelling (0 records across every corpus; the prefixed names present are 50
`chr21` and 4 `chrM`), and no comparison case does, so the shape has no row on
them.

**How often the shape occurs.** Tier 3: the shape is a chromosome name whose
`chr` prefix is not all lower-case. None of the 15,646,481 records of the
population cuts and no corpus record carries one; 85 flag sets run a
tabix-backed plugin.

#### 4.1 A gnomADc record spelt `Chr21`

**In short.** The GRCh38 gnomADc corpus (`tests/golden/116/GRCh38-plugin-gnomadc/`)
names its coverage file's contig `chr21`, and for
`21 32699959 rs2042341561 G A` both engines write twelve coverage fields,
`gnomAD_mean_cov=31.1000` and `gnomAD_total_DP_cov=22732808.0000` among them, on
each of the record's 14 transcript rows. The same record written `Chr21 32699959 cap_chr G A`
gets those fields from vep-rs on all 14 rows and none from Ensembl VEP, whose 14
rows are otherwise identical, column for column.

## Writing a plugin

How the two plugin systems work, the C ABI a shared-library plugin exports, and
how a new built-in plugin is added and tested are on
[Writing a plugin](developers/writing-a-plugin.md).
