# Plugin Guide

vep-rs registers 11 built-in annotation plugins covering common clinical and research use cases, and loads C-ABI shared-library plugins built against the `vep-plugin` crate without a rebuild (see [Built-in vs Dynamic Plugins](#built-in-vs-dynamic-plugins)). Plugins add fields to the output (Extra column in VEP format, CSQ in VCF, additional fields in JSON).

**Support is per reference genome**, because each annotation database decides which
builds it publishes: 8 plugins are supported on both GRCh37 and GRCh38 (CADD, REVEL,
gnomADc, AlphaMissense, dbscSNV, LoFTEE, LoFtool, pLI); GWAS is GRCh38-only and
SpliceAI is GRCh38-SNV-only; dbNSFP is not supported. The
[support matrix](#available-plugins) below gives the reasons, each of which is a property
of the upstream database rather than of vep-rs.

**Pass `--assembly GRCh37` or `--assembly GRCh38` whenever plugins are enabled.** It
is load-bearing, not just cache selection: REVEL uses it to choose which coordinate
column its index was built on, and CADD, REVEL, dbscSNV and AlphaMissense use it to
**refuse a data file whose own assembly contradicts the run** rather than silently
annotating nothing. Those guards are inert without it, and REVEL falls back to the
GRCh37 column with only a warning.

## Available Plugins

The tables below group the built-in plugins by mechanism, as `crates/vep-builtins/README.md`
does: how a plugin reads its data and matches a variant decides which data file it needs and
where its fields land. The GRCh37 and GRCh38 columns state which builds each upstream
database publishes.

### Annotation-store plugins

Each reads a bgzipped, tabix-indexed data file, or the `.vpd` binary store built from it
([Data File Preparation](#data-file-preparation)), and matches records by position and then
by whatever the Mechanism column adds.

| Plugin        | Description                                                  | Mechanism                                                                    | GRCh37 | GRCh38 | Output Fields                                                                                                    | Data Source                                                                                                                                                                                                                      |
| ------------- | ------------------------------------------------------------ | ---------------------------------------------------------------------------- | ------ | ------ | ---------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| CADD          | Combined Annotation Dependent Depletion scores               | position + allele match                                                      | Yes    | Yes    | `CADD_PHRED`, `CADD_RAW`                                                                                         | [CADD download](https://cadd.gs.washington.edu/download)                                                                                                                                                                         |
| REVEL         | Rare Exome Variant Ensemble Learner (missense pathogenicity) | position + allele match, transcript match unless `no_match=1`; missense only | Yes    | Yes    | `REVEL_score`                                                                                                    | [REVEL download](https://sites.google.com/site/revelgenomics/downloads)                                                                                                                                                          |
| gnomADc       | gnomAD coverage statistics                                   | position only                                                                | Yes    | Yes    | `gnomADc_<column>` for every header column after chromosome and position (`gnomADc_mean`, `gnomADc_over_1`, ...) | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads)                                                                                                                                                                  |
| AlphaMissense | DeepMind AlphaMissense pathogenicity predictions             | position + allele match                                                      | Yes    | Yes    | `am_class`, `am_pathogenicity`                                                                                   | [AlphaMissense data](https://zenodo.org/records/8208688)                                                                                                                                                                         |
| dbscSNV       | Splice site predictions (AdaBoost + Random Forest)           | position + allele match, per transcript                                      | Yes    | Yes    | `dbscSNV_ADA_SCORE`, `dbscSNV_RF_SCORE`                                                                          | [dbscSNV](https://sites.google.com/site/jpopgen/dbNSFP); the plugin reads the 6-column `dbscSNV1.1_<assembly>.txt.gz` redistributed by [CADD](https://kircherlab.bihealth.org/download/CADD/v1.7/) (see [preparation](#dbscsnv)) |

SpliceAI (position + allele + gene symbol match) and dbNSFP (position + allele + transcript
match, per transcript) are annotation-store plugins too; their support is restricted, see
[Restricted: SpliceAI](#restricted-spliceai-grch38-snv-only) and
[Not supported: dbNSFP](#not-supported-dbnsfp).

gnomADc names its fields after the file's own header, so the GRCh37 v2.1 coverage
file yields `gnomADc_median` and the GRCh38 v4 file `gnomADc_median_approx` and
`gnomADc_total_DP`. These names are vep-rs's own: the Perl gnomADc plugin writes
`gnomAD_mean_cov`, `gnomAD_median_cov` and `gnomAD_<N>x_cov` instead.

### Consequence filter: LoFTEE

LoFTEE is not a score lookup. It filters the consequences vep-rs has already assigned
(`stop_gained`, `frameshift_variant`, `splice_acceptor_variant` and `splice_donor_variant`
on protein-coding transcripts) against the transcript model, the `--fasta` splice motifs
and GERP, and emits a loss-of-function confidence call; the parameters are under
[Plugin Parameter Formats](#plugin-parameter-formats).

| Plugin | Description                                        | Mechanism                                                                   | GRCh37 | GRCh38 | Output Fields                                | Data Source                                                                                                                                                            |
| ------ | -------------------------------------------------- | --------------------------------------------------------------------------- | ------ | ------ | -------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| LoFTEE | Loss-of-function confidence (HC/LC) + filter flags | consequence filter over the transcript model, `--fasta` splice motifs, GERP | Yes    | Yes    | `LoF`, `LoF_filter`, `LoF_flags`, `LoF_info` | [LoFTEE data](https://personal.broadinstitute.org/konradk/loftee_data/) -- per assembly: `GRCh38/` (GERP bigWig) or `GRCh37/` (GERP tabix TSV), plus human_ancestor.fa |

### Flat-file plugins

Each loads its whole data file at init and looks records up in memory, by gene symbol or
by position.

| Plugin  | Description                                 | Mechanism                    | GRCh37 | GRCh38 | Output Fields                                                                                                                                                                    | Data Source                                                                                                                                                                                    |
| ------- | ------------------------------------------- | ---------------------------- | ------ | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| LoFtool | Gene intolerance to functional variation    | gene symbol lookup           | Yes    | Yes    | `LoFtool_percentile`                                                                                                                                                             | `LoFtool_scores.txt` from the [Ensembl VEP_plugins](https://github.com/Ensembl/VEP_plugins) repository; vep-rs does not bundle it                                                              |
| pLI     | Probability of Loss-of-function Intolerance | gene symbol lookup           | Yes    | Yes    | `pLI_gene_value`                                                                                                                                                                 | `pLI_values.txt` from the [Ensembl VEP_plugins](https://github.com/Ensembl/VEP_plugins) repository (ExAC r0.3 gene pLI); any uncompressed TSV with `gene` and `pLI` header columns is accepted |
| GWAS    | NHGRI-EBI GWAS Catalog associations         | position + risk allele match | No     | Yes    | `GWAS_accessions`, `GWAS_associated_gene`, `GWAS_beta_coef`, `GWAS_odds_ratio`, `GWAS_p_value`, `GWAS_pmid`, `GWAS_risk_allele`, `GWAS_study` (the Perl GWAS plugin's field set) | [GWAS Catalog downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads)                                                                                                                       |

LoFtool and pLI are keyed by gene symbol, so one data file serves both assemblies. GWAS is
GRCh38-only because the NHGRI-EBI catalog publishes `CHR_POS` in GRCh38 coordinates only.

### Restricted: SpliceAI (GRCh38 SNV only)

| Plugin   | Description                           | Output Fields                                                                                                                                                                                                  | Data Source                                                                                                          |
| -------- | ------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| SpliceAI | Deep learning splice site predictions | `SpliceAI_pred_DS_AG`, `SpliceAI_pred_DS_AL`, `SpliceAI_pred_DS_DG`, `SpliceAI_pred_DS_DL`, `SpliceAI_pred_DP_AG`, `SpliceAI_pred_DP_AL`, `SpliceAI_pred_DP_DG`, `SpliceAI_pred_DP_DL`, `SpliceAI_pred_SYMBOL` | [Ensembl MANE release](https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/) (GRCh38 SNV) |

The supported SpliceAI data is the openly published GRCh38 SNV score file from
the Ensembl MANE release. GRCh37 scores and the indel files for both assemblies
are distributed solely through Illumina BaseSpace, which requires an account
signup, so they are not supported out of the gate; they drop in unchanged if
you hold credentials (same parameters, same tabix contract).

### Not supported: dbNSFP

`--plugin dbNSFP` resolves to a built-in plugin, but dbNSFP is
**not a supported plugin**: its license grants free use only for academic,
non-commercial purposes and requires a paid commercial agreement for any other
use, so vep-rs neither distributes the data nor validates against it. A user who
holds a dbNSFP license can point the plugin at their own file and receive
user-selected columns prefixed `dbNSFP_`. Native SIFT and PolyPhen decoding needs
no dbNSFP file (see
[cache-setup.md](cache-setup.md#sift-and-polyphen-predictions)).

## Usage

Plugins are loaded via `--plugin Name,param1,param2,...` and can be specified multiple times:

```bash
vep -i input.vcf -o output.txt --offline --cache \
    --json_cache /path/to/json_cache --assembly GRCh38 \
    --plugin CADD,snv=/path/to/cadd/whole_genome_SNVs.tsv.gz,indels=/path/to/cadd/gnomad.genomes.r4.0.indel.tsv.gz \
    --plugin REVEL,file=/path/to/revel/new_tabbed_revel_grch38.tsv.gz \
    --plugin gnomADc,/path/to/gnomad/gnomad.exomes.coverage.tsv.gz
```

### Plugin Parameter Formats

**CADD** -- named parameters for file types:

```bash
--plugin CADD,snv=/path/to/snv.tsv.gz,indels=/path/to/indels.tsv.gz
--plugin CADD,snv=/path/to/snv.tsv.gz,indels=/path/to/indels.tsv.gz,force_annotate=1
```

**REVEL** -- named file parameter with optional transcript matching control:

```bash
--plugin REVEL,file=/path/to/revel.tsv.gz
--plugin REVEL,file=/path/to/revel.tsv.gz,no_match=1
```

**SpliceAI** -- separate files for SNVs and indels, optional cutoff. Only the
GRCh38 SNV file is supported from open data; `indel=` takes a BaseSpace-gated
file if you have one:

```bash
--plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz
--plugin SpliceAI,snv=/path/to/spliceai_snv.vcf.gz,indel=/path/to/spliceai_indel.vcf.gz,cutoff=0.5
```

**dbNSFP** (not supported; for license holders only) -- file path followed by
column names to extract:

```bash
--plugin dbNSFP,/path/to/dbNSFP.gz,SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred
```

**gnomADc** -- single file path:

```bash
--plugin gnomADc,/path/to/gnomad.coverage.tsv.gz
```

**AlphaMissense** -- single file path:

```bash
--plugin AlphaMissense,/path/to/AlphaMissense_hg38.tsv.gz
```

**dbscSNV** -- single file path:

```bash
--plugin dbscSNV,/path/to/dbscSNV1.1_GRCh38.txt.gz
```

**LoFtool** -- scores file path (required). vep-rs does not bundle
`LoFtool_scores.txt`; fetch it from the VEP_plugins repository and pass its path:

```bash
--plugin LoFtool,/path/to/LoFtool_scores.txt
```

**pLI** -- gene constraint file:

```bash
--plugin pLI,/path/to/pLI_values.txt
```

**LoFTEE** -- loss-of-function confidence flagging (consequence-filtering, not a
lookup). Targets parity with the upstream LoFTEE `grch38` branch defaults
(splice-prediction layer off). Built-in names are matched exactly as written in the
matrix: Perl VEP loads LoFTEE as `--plugin LoF,...` and vep-rs as `--plugin LoFTEE,...`.
Named `key=value` parameters, all optional (LoFTEE's native `key:value` form is
accepted alongside):

```bash
# GRCh38 -- upstream ships GERP as a bigWig
--plugin LoFTEE,human_ancestor_fa=/path/to/loftee/grch38/human_ancestor.fa,gerp_bigwig=/path/to/loftee/grch38/gerp_conservation_scores.homo_sapiens.GRCh38.bw

# GRCh37 -- upstream ships GERP as a tabix per-base TSV instead
--plugin LoFTEE,human_ancestor_fa=/path/to/loftee/grch37/human_ancestor.fa,gerp_tabix=/path/to/loftee/grch37/GERP_scores.final.sorted.txt.gz
```

The two upstream LoFTEE bundles are **structurally different**, which is why the
GERP parameter differs per assembly. vep-rs reads both formats natively.

- `human_ancestor_fa` -- PLAIN (uncompressed) indexed FASTA (`.fa` + `.fai`).
  Enables the `ANC_ALLELE` filter (SNP-only). Decompress the upstream
  `human_ancestor.fa.gz` first; vep-rs does not read bgzip here.
- `gerp_bigwig` -- GERP conservation bigWig, as shipped in the upstream **grch38**
  bundle. Enables the default GERP-weighted `END_TRUNC` filter.
- `gerp_tabix` (alias `gerp_scores`) -- GERP as a tabix-indexed per-base TSV, as
  shipped in the upstream **master/GRCh37** bundle
  (`GERP_scores.final.sorted.txt.gz` + `.tbi`). Same filter, different on-disk
  format.
- `gerp` -- assembly-agnostic spelling; classifies the source by file extension.
- Supply exactly one of the two. With neither, set `use_gerp_end_trunc=false` to
  use the 5%-percentile `END_TRUNC` fallback (no data file); an unreadable GERP
  path is a hard error at plugin init rather than a silent degradation.
- `min_intron_size` (default 15), `filter_position` (default 0.05),
  `use_gerp_end_trunc` (default true), `check_complete_cds` (default false) --
  tunables matching the Perl plugin.
- `gerp_end_trunc_cutoff` -- **derived from `--assembly`** unless set explicitly,
  because the two upstream branches use cutoffs of opposite sign: `-58` on grch38,
  `+180` on master/GRCh37. Setting it pins that value for both.
- `--fasta` (the pipeline reference genome) is required for the splice-motif
  checks (the `GC_TO_GT_DONOR` filter and the `NON_CAN_SPLICE` and `NAGNAG_SITE`
  flags); without it those checks are skipped.

LoFTEE only annotates `stop_gained`, `frameshift_variant`,
`splice_acceptor_variant`, and `splice_donor_variant` on protein-coding
transcripts; everything else gets no `LoF` field. The MaxEntScan +
logistic-regression + de-novo-donor-SVM splice-prediction layer and the `OS`
confidence class (off by default on the `grch38` branch) are not implemented.
The PhyloCSF `LoF_flags` (`PHYLOCSF_WEAK`/`PHYLOCSF_UNLIKELY_ORF`) are a
documented gap: they are flags, not filters, so they never change the HC/LC
call, and they would require a SQLite dependency; the `conservation_file`
parameter is accepted for arg parity but is a no-op.

## Data File Preparation

All tabix plugins require bgzipped, tabix-indexed data files (`.tsv.gz` + `.tsv.gz.tbi` or `.vcf.gz` + `.vcf.gz.tbi`). The one exception is a `.vpd` binary store built from the data file by `vep-cache-converter convert-plugin-data` (its `.vpdi` index beside it): when a store named after the data file with `.vpd` in place of `.gz` sits beside it (`whole_genome_SNVs.tsv.vpd` beside `whole_genome_SNVs.tsv.gz`), the plugin reads its records from the store instead of through tabix.

### CADD

Download from https://cadd.gs.washington.edu/download:

- `whole_genome_SNVs.tsv.gz` + `.tbi` (SNVs)
- the release's indel file + `.tbi`: `gnomad.genomes.r4.0.indel.tsv.gz` (v1.7 GRCh38) or `gnomad.genomes-exomes.r4.0.indel.tsv.gz` (v1.7 GRCh37); v1.6 GRCh37 and earlier releases name it `InDels.tsv.gz`

### REVEL

Download from https://sites.google.com/site/revelgenomics/downloads and process:

```bash
unzip revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv

# GRCh37:
tabix -f -s 1 -b 2 -e 2 new_tabbed_revel.tsv.gz

# GRCh38:
zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz
```

### SpliceAI

Supported data: the GRCh38 SNV score VCF (bgzipped, with `.tbi`) from the Ensembl MANE release at https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/. Remote tabix works against it, so a chr21 subset can be sliced without the full download. GRCh37 scores and both assemblies' indel files are available only through Illumina BaseSpace and are not supported.

### dbNSFP

Not supported (free use is limited to academic, non-commercial purposes; any other use requires a paid commercial license); no data file is distributed or validated. License holders can use their own tabix-indexed dbNSFP file directly.

### dbscSNV

Use the 6-column `dbscSNV1.1_GRCh37.txt.gz` or `dbscSNV1.1_GRCh38.txt.gz` (+ `.tbi`) redistributed under https://kircherlab.bihealth.org/download/CADD/v1.7/ (`GRCh37/`, `GRCh38/`). The plugin reads chromosome, position, ref and alt from the first four columns and `ada_score` and `rf_score` by header name; the header may lack its leading `#` and the lines may end in CRLF, both of which vep-rs accepts. The GRCh38 file is recognised by its `hg38_chr` first column. The preparation the Perl dbscSNV plugin documents for the full dbscSNV1.1 zip keeps every column with `chr` first and, on GRCh38, indexes columns 5 and 6; the plugin does not read those columns, and under `--assembly GRCh38` it refuses such a file as a GRCh37 one.

### gnomADc

Download coverage files from https://gnomad.broadinstitute.org/downloads. No tabix index is published, so build one on the chromosome and position columns (`bgzip`, then `tabix -s 1 -b 2 -e 2`). The GRCh38 v3 and v4 files carry a single `locus` column (`chr1:11819`); split it into chromosome and position first (`sed "1s/locus/chr\tpos/; s/:/\t/g"`). The header row works with or without a leading `#`.

### Preparing data files with `scripts/data/setup_plugin_data.sh`

Stage plugin data files with `scripts/data/setup_plugin_data.sh`. It does not
download from the upstream databases itself. `--source upstream`, the default,
prints, per plugin, the upstream URL and the preparation each file needs (header
fixes, per-assembly tabix indexes, decompression) for you to run; `--source s3`
syncs an already-prepared copy from an S3 bucket of your own, named by `--bucket`
or `$VEP_BENCHMARK_BUCKET`, that you seeded from those steps:

```bash
scripts/data/setup_plugin_data.sh --assembly GRCh38 --tier quick --source upstream            # print sources for chr21 subsets
scripts/data/setup_plugin_data.sh --assembly GRCh38 --tier full --source s3 --bucket my-bucket   # sync genome-wide files (roughly 135 GB for GRCh38, 87 GB of it CADD) from your bucket
```

`--assembly` is required: it selects which files are staged and which plugins
apply, and an unsupported combination is skipped with its reason rather than
substituted from the wrong build. Files land under
`~/.vep/plugin_data/<assembly>/<tier>/` by default, in per-plugin subdirectories
(`cadd/whole_genome_SNVs.tsv.gz`, `revel/revel.tsv.gz`, ...), the layout
`scripts/data/plugin_data_layout.sh` defines for the data-preparation and concordance scripts
alike.

## Architecture

### Built-in vs Dynamic Plugins

vep-rs has two plugin systems:

**Built-in plugins** (`vep-builtins` crate) are compiled into the binary. They operate on native Rust types with no serialization overhead. This is the recommended approach for canonical VEP plugins.

**Dynamic plugins** (`vep-plugin` crate) are compiled as separate `cdylib` shared libraries loaded at runtime via C FFI. Communication uses JSON serialization. This is the route for extensions that ship separately from vep-rs. The host side is `crates/vep-cli/src/dylib_plugins.rs`.

When `--plugin Name,...` is specified, the runner checks the built-in registry first. A name no built-in claims is resolved as a dylib: as a file path if one exists at `Name`, otherwise as `libName.so` and then `libName.dylib` under `--dir_plugins`. The name is the library file name, so `--dir_plugins /opt/vep/plugins --plugin echo_plugin` loads `/opt/vep/plugins/libecho_plugin.so`; the string a plugin returns from `vep_plugin_name` appears only in log and error messages. A name that resolves to neither is reported (`plugin 'Name' not found as a built-in or as a dylib under --dir_plugins; skipped`) and the run continues without it, which is Ensembl VEP's default without `--safe`. A file that exists but does not open as a shared library, lacks one of the eight symbols, reports an ABI version other than 1, or returns non-zero from `vep_plugin_init` aborts the run before annotation starts.

#### The C ABI (version 1)

A plugin exports these eight symbols with these signatures (`crates/vep-plugin/src/abi.rs`):

| Symbol                   | Signature                                                              | Contract                                                                                                                                                         |
| ------------------------ | ---------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `vep_plugin_abi_version` | `extern "C" fn() -> u32`                                               | Returns `1`; any other value is a load error.                                                                                                                    |
| `vep_plugin_create`      | `extern "C" fn() -> *mut c_void`                                       | Returns an opaque instance pointer; NULL is a load error.                                                                                                        |
| `vep_plugin_destroy`     | `extern "C" fn(*mut c_void)`                                           | Frees the instance; called when the host drops the plugin.                                                                                                       |
| `vep_plugin_name`        | `extern "C" fn(*const c_void) -> *const c_char`                        | NUL-terminated name; the pointer must stay valid for the instance's lifetime.                                                                                    |
| `vep_plugin_version`     | `extern "C" fn(*const c_void) -> *const c_char`                        | NUL-terminated version string, same lifetime rule.                                                                                                               |
| `vep_plugin_init`        | `extern "C" fn(*mut c_void, *const *const c_char, usize) -> i32`       | Receives the `--plugin` parameters after the name as NUL-terminated strings; returns 0 on success. Called only when at least one parameter was given.            |
| `vep_plugin_run`         | `extern "C" fn(*const c_void, *const c_char, *const c_char) -> *mut c_char` | `(instance, consequence JSON, variant JSON)`; returns a NUL-terminated JSON object that the host frees with `vep_plugin_free`. NULL signals an error.           |
| `vep_plugin_free`        | `extern "C" fn(*mut c_char)`                                           | Frees a string returned by `vep_plugin_run`.                                                                                                                     |

#### What `run` receives

After the built-in plugins have run on a batch, the host calls `vep_plugin_run` once per entry in each variant's `transcript_consequences` (transcript, regulatory-feature and motif-feature consequences alike), for every loaded dylib in `--plugin` order. A variant with no such entries, such as an intergenic one, is never passed to a plugin. The two JSON documents are the serde serialization of `TranscriptConsequence` (`crates/vep-core/src/consequence.rs`) and `InputVariant` (`crates/vep-core/src/variant.rs`); the Rust field names are the JSON keys, and those two structs are the authoritative schema. Absent optional values are `null`. The fields a plugin will find:

- Consequence: `transcript_id`, `gene_id`, `gene_symbol`, `gene_symbol_source`, `hgnc_id`, `consequences` (an array of the `Consequence` enum's variant names, `"MissenseVariant"`, `"StopGained"`, not the SO terms), `impact` (`"HIGH"`, `"MODERATE"`, `"LOW"`, `"MODIFIER"`), `biotype`, `canonical`, `cdna_position`, `cds_position`, `protein_position` (strings such as `"428"` or `"301-303"`), `amino_acids`, `codons`, `protein_id`, `distance`, `strand` (`1` or `-1`), `exon`, `intron`, `hgvsc`, `hgvsp`, `sift`, `polyphen`, `domains` (an array of `[source, id]` pairs), `feature_type` (`"Transcript"`, `"RegulatoryFeature"`, `"MotifFeature"`), `flags`, `tsl`, `mane_select`, `mane_plus_clinical`, `appris`, `ccds`, `swissprot`, `trembl`, `refseq`, `feature_start`, `feature_end`, and `plugin_data`: the fields the built-in plugins have already written to this consequence (`REVEL_score`, `LoF`, ...), present only when non-empty. `loftee_ctx` (the exon/intron model) is present only when LoFTEE is enabled.
- Variant: `chr` (normalized, `"21"`), `original_chr` (as input, `"chr21"`), `start`, `end` (1-based inclusive; `end` is `start - 1` for an insertion), `strand` (`"Forward"` or `"Reverse"`), `allele_string` (`"A/G"`, `"AC/-"`, `"-/TT"`), `ref_allele` and `alt_alleles` (serialized as arrays of byte values, so read `allele_string` for the sequences), `allele_index`, `id`, `variant_class` (`"Snv"`, `"Insertion"`, `"Deletion"`, `"StructuralDeletion"`, ...), `transcript_consequences` (the variant's whole list, so the variant document repeats every consequence), `most_severe_consequence`, `colocated_variants`, `existing_variation`, `plugin_data` (variant-level built-in fields such as `CADD_PHRED`, present only when non-empty), the structural-variant fields (`is_structural`, `sv_end`, `sv_type`, `sv_len`, `ci_pos`, `ci_end`, `mate_id`, `mate_chr`, `mate_pos`, `is_single_breakend`), `minimised`, `original_allele_string`, `original_start`, `original_end`, `raw_input`, `input_record` and the `Uploaded_variation` helpers `uploaded_allele_string` and `record_allele_string_multi`.

#### What `run` returns, and where it lands

The returned string must be a JSON object. Each entry becomes one field on that consequence: a string value is stored verbatim, any other value as its JSON text (`0.42`, `null`), and a key already present, from a built-in or an earlier dylib, is overwritten. The fields appear in the Extra column of the default output format and in the per-consequence objects of `--json`, where, as for every other field, the key is lowercased and a numeric-looking value is written as a number. The ABI has no header-declaration symbol, so a dylib field gets no `##` description line in the header and no column in `--tab`, `--vcf` (CSQ) or `--output_format parquet` output; those writers take their plugin columns from the built-in plugins' header info. A `run` call that returns NULL, invalid JSON, or a JSON value that is not an object aborts the run.

#### Worked example

`crates/vep-plugin/examples/echo_plugin.rs` is a complete plugin: it echoes the consequence's `transcript_id` as `ECHO_FEATURE` and its first parameter, when one was given, as `ECHO_PARAM`. It reads the one field it needs with a string scan, so it adds no dependency. Build it as the shared library the loader opens, then pass it by path or by name:

```bash
cargo build --example echo_plugin -p vep-plugin
# -> target/debug/examples/libecho_plugin.so (Linux) or libecho_plugin.dylib (macOS)

vep -i input.vcf -o output.txt --offline --cache --json_cache /path/to/json_cache --assembly GRCh38 \
    --plugin target/debug/examples/libecho_plugin.so,hello

vep -i input.vcf -o output.txt --offline --cache --json_cache /path/to/json_cache --assembly GRCh38 \
    --dir_plugins target/debug/examples --plugin echo_plugin,hello
```

Every consequence line then carries `ECHO_FEATURE=<transcript_id>;ECHO_PARAM=hello` in its Extra column. `cargo test -p vep-plugin` builds the example and runs `tests/echo_roundtrip.rs`, which loads it through the same loader, by path and by name.

### Batch Query Strategy

Tabix-based plugins use the `TabixAnnotator` engine (or the `BinaryAnnotator` when a `.vpd` store sits beside the data file), which batches I/O across the entire variant buffer (`--buffer_size`, default 5000):

1. Collect query regions from all variants in the batch
2. Group the regions by chromosome and sort each chromosome's regions by position
3. Merge overlapping regions, bridging gaps of up to 1000 bp
4. Issue one tabix query per merged region rather than one per variant
5. Index returned records by position
6. Match each variant against the indexed records

Perl VEP's `BaseVepTabixPlugin` runs one tabix query per variant whenever a plugin sets no region expansion, as its CADD and gnomADc plugins do. The batched strategy replaces that with one query per merged region. Its contribution to run time is not separately measured: the published wall-time measurements were taken with no plugin enabled.

### Plugin Execution in the Pipeline

Plugins execute after consequence assignment:

1. Input batch parsed
2. Consequences calculated (parallel if `--fork > 1`)
3. Co-located variants matched
4. **`registry.run_batch(&mut variants)`**: every built-in plugin's `prefetch()` runs in parallel across plugins, then each plugin's `annotate()` (or `run_batch()`, for a plugin without `prefetch`) runs sequentially on the full batch
5. Dynamic (dylib) plugins run on every transcript consequence
6. Output written

Variant-level plugins (CADD, gnomADc, AlphaMissense, GWAS) write to `variant.plugin_data`. Transcript-level plugins write to `consequence.plugin_data`: REVEL, SpliceAI, dbNSFP and dbscSNV do so in their batch `annotate()` after a tabix `prefetch()`, and LoFtool, pLI and LoFTEE through the per-consequence `run()` method. LoFTEE is the only consequence-_filtering_ plugin: it inspects the transcript exon/intron model (carried on the consequence via `LofteeContext`, populated only when LoFTEE is active) rather than looking up a precomputed score by position.

## Adding a New Plugin

1. Create a new module in `crates/vep-builtins/src/plugins/my_plugin.rs`
2. Implement the `BuiltinPlugin` trait:

```rust
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::variant::InputVariant;
use vep_core::consequence::TranscriptConsequence;
use indexmap::IndexMap;

pub struct MyPlugin { /* fields */ }

impl BuiltinPlugin for MyPlugin {
    fn name(&self) -> &str { "MyPlugin" }
    fn header_info(&self) -> Vec<(String, String)> {
        vec![("MY_FIELD".into(), "Description".into())]
    }
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        // Parse params, open data files
        Ok(())
    }
    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Return key-value pairs for the Extra column
        Ok(IndexMap::new())
    }
    // Implement prefetch() + annotate() for batch I/O (tabix plugins); see below
}
```

3. Add the module to `crates/vep-builtins/src/plugins/mod.rs`
4. Register the name in `crates/vep-builtins/src/registry.rs` (`create_builtin()`)
5. Add tests

### For Tabix-Based Plugins

Open the data file with `open_best_store` (tabix, or the `.vpd` store when present) and implement the two-phase protocol the built-in tabix plugins use: `prefetch()` takes a shared reference to the batch, so the registry runs it for every plugin in parallel, and `annotate()` applies the fetched records sequentially. `run_batch()` is the single-phase fallback for a plugin that returns `None` from `prefetch()`.

```rust
use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::{BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};

fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
    let Some(store) = &self.store else { return Ok(None) };
    let regions: Vec<(&str, u64, u64)> = variants.iter()
        .map(|v| (v.chr.as_str(), v.start, v.end))
        .collect();
    Ok(Some(PrefetchData::new(store.query_batch(&regions)?)))
}

fn annotate(&self, variants: &mut [InputVariant], data: PrefetchData) -> Result<(), PluginError> {
    let batch: BatchQueryResult = data.downcast()?;
    for variant in variants.iter_mut() {
        let records = batch.lookup(&variant.chr, variant.start, variant.end);
        // Match alleles, then insert into variant.plugin_data or a consequence's plugin_data
    }
    Ok(())
}
```

## Plugin Concordance Testing

`scripts/concordance/run_concordance.sh` runs Perl VEP (in Docker) and vep-rs on the same inputs with the same plugins and compares each plugin's fields. Every run needs `--perl-cache-dir` (the directory that contains `homo_sapiens/<version>_<assembly>/`), `--json-cache-dir` (a vep-rs JSON cache directory containing `info.json`), `--fasta` (an indexed reference FASTA given to both engines; `--no-fasta` is the explicit opt-out and is recorded in `provenance.json`) and a benchmark directory (`--benchmark-dir`, or `$VEP_BENCHMARK_DIR`); it exits 1 before running anything when one is missing. The plugin flags:

| Flag                          | Default                                | Description                                                                                                     |
| ----------------------------- | -------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| `--plugins <list>`            | (none)                                 | Comma-separated plugin names to test                                                                            |
| `--plugin-data-dir <path>`    | `~/.vep/plugin_data/<assembly>/<tier>` | Plugin data directory (lower-case assembly, from `--assembly`); a path that names the other assembly is refused |
| `--plugin-tier <quick\|full>` | `quick`                                | Subdirectory used by the default data directory                                                                 |

Two data tiers share the layout of `scripts/data/plugin_data_layout.sh`. `quick` holds chromosome-21 slices of the full files, enough for the chr21 benchmark VCFs that `scripts/data/download_real_world_vcfs.sh` seeds (a genome-wide input receives plugin values only on its chr21 records; where the upstream host serves a working remote tabix index, `tabix <url> 21 | bgzip` produces the slice without a full download). `full` holds the genome-wide files. `setup_plugin_data.sh --assembly <asm> --tier <tier>` stages either, as described above.

```bash
scripts/concordance/run_concordance.sh \
  --assembly GRCh38 \
  --perl-cache-dir /path/to/vep_cache \
  --json-cache-dir ~/.vep/json_cache/homo_sapiens/115_GRCh38 \
  --fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --benchmark-dir .vep/vcf/grch38/all_variants \
  --plugins CADD,REVEL,gnomADc \
  --plugin-data-dir ~/.vep/plugin_data/grch38/quick \
  --mode smoke
```

With `--plugins` set, the harness maps each plugin to a `--plugin Name,params` flag for both engines, with file paths from the layout table (Perl's remapped under the container mount `/work/plugin_data`); mounts the data directory into the Perl container, plus a `VEP_plugins/` checkout when one exists at `$VEP_PLUGINS_DIR` (default: a `VEP_plugins/` sibling of the repository; otherwise the `ensemblorg/ensembl-vep` image's own `/plugins` copies serve); keys its Perl output cache on the plugin list, the data directory and a size-and-mtime fingerprint of its files, so cached outputs from another plugin configuration are not reused; and, after the consequence comparison, runs `compare_plugin_outputs.py` per plugin and writes `reports/plugin_concordance.json` (per plugin: `match`, `mismatch`, unrounded `concordance_pct`, plus `tolerance`, `min_annotated`, `file_count` and per-file `results`). A failing plugin makes the harness exit non-zero. Unless `--no-perf` is passed, a baseline vep-rs run without plugins is also timed on the first prepared input and `reports/perf/plugin_timing.json` records the plugin overhead and the Perl-with-plugins over vep-rs-with-plugins wall time.

`compare_plugin_outputs.py` also runs standalone on two default-format outputs:

```bash
python3 scripts/concordance/compare_plugin_outputs.py \
    --perl-output perl_vep_with_cadd.txt \
    --rust-output rust_vep_with_cadd.txt \
    --plugins CADD \
    --fields CADD_PHRED,CADD_RAW
```

**LoFTEE.** A docker-mode run with LoFTEE also needs `LOFTEE_PERL_PATH` set to a checkout of https://github.com/konradjk/loftee (branch `master` for GRCh37, `grch38` for GRCh38), which the harness mounts into the Perl container; `LOFTEE_DATA_DIR` overrides the default data location, `<plugin-data-dir>/loftee/<assembly>/`. vep-rs reads the plain `human_ancestor.fa`, so decompress the upstream bgzipped file and re-index it (`bgzip -d`, then `samtools faidx`); the harness hands Perl's `LoF` plugin the bgzipped `human_ancestor.fa.gz`, and `setup_plugin_data.sh` reports the directory READY only when the plain file, the compressed copy and all three indexes are present. The `END_TRUNC` cutoff differs in sign between the two upstream branches (`-58` on `grch38`, `+180` on `master`); vep-rs derives it from `--assembly`, so pass that flag rather than `gerp_end_trunc_cutoff=`. Compare `LoF` and `LoF_filter` only (`--fields LoF,LoF_filter --tolerance 0`): `LoF_info` and `LoF_flags` are informational and differ in formatting between the two implementations. Three things to know when reading a LoFTEE comparison: on GRCh38, Perl's plugin reads the GERP bigWig through lossy zoom-level summaries while vep-rs sums exact per-base values, so a borderline variant can land on either side of the `END_TRUNC` cutoff; on GRCh37, Perl's `master` plugin takes whole intermediate exons from the `gerp_exons` table of its conservation database, which the harness does not pass, so a record whose exon walk needs a whole-exon total dies on the Perl side and is counted as vep-rs-only; and Perl's `ANC_ALLELE` check shells out to `samtools faidx`, which the `ensemblorg/ensembl-vep` image does not ship, so that filter is never applied on the Perl side in Docker.

No per-plugin agreement figure is published: the released concordance results are core-annotation F1 with no plugin enabled, and a plugin comparison here is a check on your own data and binary, not a reproduction of a published number.
