# vep-rs

[![CI](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![MSRV](https://img.shields.io/badge/MSRV-1.89-orange.svg)](#installation)

Documentation: https://natera-open-source.github.io/vep-rs/

vep-rs is a fast, memory-efficient Rust reimplementation of Ensembl VEP, the Variant Effect Predictor, for VCF variant annotation. It is a drop-in replacement for Ensembl VEP: it takes VEP's inputs and flags, writes VEP's consequence vocabulary in VEP's output formats, and annotates population-scale callsets up to 432 times faster than Ensembl VEP.

In this README, "Ensembl VEP" is the reference implementation, `ensemblorg/ensembl-vep`, written in Perl; unqualified, it means release 116.2. The accompanying paper calls the same engine "VEP".

## Why vep-rs instead of Ensembl VEP or fastVEP

|                                  | vep-rs v0.6.0                                                                       |
| -------------------------------- | ----------------------------------------------------------------------------------- |
| **Speed vs Ensembl VEP 116.2**   | 429x faster on ARM (Graviton5), 315x on x86; 275x to 432x per population dataset    |
| **Speed vs fastVEP v0.5.0**      | 10.7x faster on ARM, 8.7x on x86; 5.0x to 18.9x per population dataset              |
| **All of gnomAD v4.1**           | 759,336,320 variants, 24 chromosomes: 591 s on one ARM machine, 799 s on x86        |
| **Peak memory**                  | At most 2.8 GB RSS on any chromosome                                                |
| **Concordance with Ensembl VEP** | Adjusted F1 1.000000000 on all six whole-genome datasets, raw F1 0.999996998 pooled |
| **Plugins**                      | 31 of the 87 Ensembl VEP plugins built in, verified on real data or, where the data cannot be redistributed, on synthetic fixtures in the layout each plugin's documentation gives ([the support table](docs/plugins.md#support-table)); C-ABI dynamic-library plugins |
| **Output formats**               | VEP default, VCF, JSON, Tab, typed Parquet                                          |
| **Threading**                    | `--fork N`; the default uses every logical CPU up to 32, less an inflate-bound input's decompression threads |

The speed figures cover all of four population datasets (1,179,639,224 variants) under one protocol: the sum of per-chromosome medians at 16 threads, over 5 machines for vep-rs and fastVEP and over Ensembl VEP's one reference run per chromosome. The concordance figures cover those four plus the two genome-wide ClinVar releases, six whole-genome datasets in all (1,188,425,523 variants, 11,838,755,797 annotation tuples, 11,838,727,681 matched), against Ensembl VEP 116.2. The [Current release](#current-release-v060) section has the tables and the method.

The adjusted F1 sets aside only the documented cases where Ensembl VEP's own output is wrong, or where the two engines represent a record differently. Each case is traced to the mechanism in Ensembl VEP's source and reproducible from the record ([docs/intended-divergences.md](docs/intended-divergences.md)). Outside those cases every tuple matches, so vep-rs's remaining disagreements with Ensembl VEP are deliberate, documented corrections rather than defects.

vep-rs is verified by concordance against Ensembl VEP 116.2 on human GRCh37 and GRCh38, the whole-genome datasets of the current release, and on mouse GRCm39, fruit fly BDGP6.54, Arabidopsis thaliana TAIR10 and yeast R64-1-1. Each of the four has a golden corpus of its own ([Other species](docs/cache-species.md)).

Any other species in Ensembl release 116 or Ensembl Genomes release 63 with a GFF3 and a peptide FASTA, or with an Ensembl VEP cache to convert, builds a cache ([Cache setup](docs/cache-setup.md)) and annotates under the codon tables the cache names, without a concordance measurement. `--sift` and `--polyphen` refuse on a cache without their data, and `--af_gnomad` on a cache that declares variation columns without the gnomAD exome column; every other flag that reads data Ensembl ships for human alone runs there as [the flag table](docs/cli/cache-species-runtime.md#flags-that-read-data-ensembl-ships-for-human-alone) states. The Peak memory row is the largest per-chromosome median of the release measurement; [Installation](docs/getting-started/installation.md#memory) states what a run holds in memory and how the figure grows with the chromosomes one input file touches.

The built-in plugins are ports of Ensembl VEP's, loaded by the same names and parameters. [The plugin support table](docs/plugins.md#support-table) gives every Ensembl VEP plugin's data file, terms, assemblies and status; a `--plugin` naming one vep-rs does not provide stops the run with the reason. Dynamic-library plugins load from `--dir_plugins` ([docs/plugins.md](docs/plugins.md)). Parquet output is written natively, with no external tool.

### What "drop-in" means

vep-rs takes the same input files and flag names as Ensembl VEP and writes the same output layouts. Two things differ:

- **The transcript cache.** The cache is a directory built by `vep-cache-builder` or converted from an Ensembl VEP cache ([docs/cache-setup.md](docs/cache-setup.md)), found under `--dir_cache`, `--species` and `--cache_version` in Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout or named directly with `--json_cache`; a run with no cache stops with an error, and vep-rs never contacts a database.
- **Some flags are refused or accepted without effect.** The database-connection flags, `--database` first among them, are refused with an error naming the flag, as are `--gff` and `--gtf`, which read a transcript model from a file in place of the cache: the file is built into a cache first ([docs/cache-from-gtf-or-gff3.md](docs/cache-from-gtf-or-gff3.md)). `--debug` is accepted without effect. `--regulatory`, alone or through `--everything`, annotates regulatory and motif features from a converted Ensembl VEP cache, and `--regulatory_gff` from Ensembl's regulatory GFF files ([docs/regulatory.md](docs/regulatory.md)). `--refseq` and `--merged` select the `_refseq` or `_merged` species directory, and `--refseq` also resolves an HGVS gene name among RefSeq transcripts first. [The CLI reference](docs/cli-reference.md) and [limitations page](docs/getting-started/what-vep-rs-does-not-do.md) detail supported flags and omissions.

## Installation

Pre-built binaries for Linux x86_64, Linux aarch64, macOS aarch64 and macOS x86_64 are on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page; each Linux architecture has a measured archive and a generic one, and both run on glibc 2.17 or later. Building from source needs Rust 1.89+ ([rustup](https://rustup.rs/)):

```bash
git clone https://github.com/natera-open-source/vep-rs.git
cd vep-rs
cargo build --release
# Binary at target/release/vep
```

Each release from 0.2.0 on is also a container image, `ghcr.io/natera-open-source/vep-rs`, multi-platform from 0.5.0: the `X.Y.Z`, `X.Y` and `latest` tags carry the generic Linux binaries for linux/amd64 and linux/arm64 and run on Apple silicon and Graviton2 without a platform flag; `X.Y.Z-x86-64-v3` and `X.Y.Z-neoverse-v1` carry the measured binaries. The [Installation](docs/getting-started/installation.md) page covers each archive's CPU floor, how much slower the generic archives are (each release page reports it, or that the release did not time them), how to build for an older CPU, and the image's tags, platforms and attestation.

## Quick start

### 1. Try the bundled example

This needs no download. It annotates `examples/homo_sapiens_GRCh37.vcf`, 175 records cut from the GRCh37 golden corpus, against the pruned cache that ships beside that corpus in `tests/golden/` ([examples/README.md](examples/README.md) describes both example files):

```bash
target/release/vep \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --assembly GRCh37 \
  -i examples/homo_sapiens_GRCh37.vcf -o output.txt
```

### 2. Annotate your own data

Real data needs a JSON transcript cache, built natively with `vep-cache-builder` or converted from an Ensembl VEP cache. Then run the same VEP command you run, with `--json_cache` naming the cache directory. The [Quick start](docs/getting-started/quick-start.md) page gives the cache-builder command, the annotation command and the plugin and HGVS options. [Cache setup](docs/cache-setup.md) explains each input of the cache build.

## Current release: v0.6.0

These are the released version's own figures, measured against Ensembl VEP 116.2 under the release protocol stated under Method. The figures the paper reports, measured at its pinned build, are in [docs/published-figures.md](docs/published-figures.md) and stay as published. How each number was produced is under [Method](#method).

![Four panels comparing vep-rs v0.6.0 and fastVEP v0.5.0 against Ensembl VEP 116.2. (a) Raw F1 per dataset as a bar for each of vep-rs and fastVEP, with each tool's adjusted F1 as a marker joined to the end of the bar by a line where the two differ. (b) Wall time per dataset and architecture of vep-rs, fastVEP and Ensembl VEP, a circle each, filled for ARM and open for x86; from left to right, the circles of every dataset and architecture are those of vep-rs, fastVEP and Ensembl VEP; each wall-time cell of Ensembl VEP holds 1 run, against 5 for vep-rs and fastVEP. (c) Raw F1 per plugin (plugins verified on real data; TranscriptAnnotator on its committed table; plugins whose data cannot be redistributed on synthetic fixtures, marked synthetic) as a bar for each of vep-rs and fastVEP, with each tool's adjusted F1 as a marker joined to the end of the bar by a line where the two differ: fastVEP has a figure on 14 of the 57 rows, its marker on 12 of them, and a cross (no such source) on 43; vep-rs's adjusted F1 is exactly 1 on 57 of the 57 rows. (d) The plugins read by vep-rs and fastVEP, a bar each: 31 of 31 and 7 of 31. A marker standing at the end of its bar means that the tool's raw and adjusted F1 are equal as drawn. On a row of a plugin verified on real data, vep-rs's counts are summed over the genome-wide ClinVar dataset of the row's genome build and whole-genome population datasets of that build (GRCh37: gnomAD v2.1.1 and 1000 Genomes Phase 3; GRCh38: gnomAD v4.1 and 1000 Genomes high-coverage), except on gnomADMt GRCh38 (ClinVar alone); fastVEP's counts are summed over the genome-wide ClinVar dataset of the row's genome build and whole-genome population datasets of that build (GRCh37: gnomAD v2.1.1 and 1000 Genomes Phase 3; GRCh38: gnomAD v4.1 and 1000 Genomes high-coverage). The plugin panel's 57 rows carry 32 row names for 31 plugins: CADD is scored under CADD and CADD (structural variants). A row's figure is summed over the sets of the row that the record holds for the tool, and the tools differ in them on 4 rows: dbNSFP GRCh37 synthetic (vep-rs 5, fastVEP 1), dbNSFP GRCh38 synthetic (vep-rs 9, fastVEP 1), SpliceAI GRCh37 synthetic (vep-rs 4, fastVEP 1) and SpliceAI GRCh38 synthetic (vep-rs 4, fastVEP 1). An F1 of exactly 1 has no place on a scale of 1 - F1: it is drawn where the axis reads 1, past a break in the axis line, and 66 of the 66 adjusted markers of vep-rs stand there. The plugin panel runs over both columns: its groups from Committed table on continue in the second column, under the same title.](docs/assets/comparison.png)

**vep-rs v0.6.0 against Ensembl VEP 116.2, with fastVEP v0.5.0 as comparator.** (a) Raw F1 per dataset on a log(1 - F1) axis; a marker is its tool's adjusted F1, joined to the end of its bar by a line where the two differ. All six is the pool of the six whole-genome datasets, the ClinVar and the population datasets. (b) Wall time per dataset and architecture: a ClinVar line is the median of its runs, any other the sum of its dataset's cell medians (the record's methodology note: "a dataset's figure is the sum of its per-chromosome medians"); each wall-time cell of Ensembl VEP holds 1 run, against 5 for vep-rs and fastVEP. (c) Raw F1 per plugin: plugins verified on real data; TranscriptAnnotator on its committed table; plugins whose data cannot be redistributed on synthetic fixtures, marked synthetic. On a row of a plugin verified on real data, vep-rs's counts are summed over the genome-wide ClinVar dataset of the row's genome build and whole-genome population datasets of that build (GRCh37: gnomAD v2.1.1 and 1000 Genomes Phase 3; GRCh38: gnomAD v4.1 and 1000 Genomes high-coverage), except on gnomADMt GRCh38 (ClinVar alone); fastVEP's counts are summed over the genome-wide ClinVar dataset of the row's genome build and whole-genome population datasets of that build (GRCh37: gnomAD v2.1.1 and 1000 Genomes Phase 3; GRCh38: gnomAD v4.1 and 1000 Genomes high-coverage). The bar is each tool's per-feature raw F1, the marker its tool's adjusted F1, joined to the end of its bar by a line where the two differ, as in (a); a cross marks a plugin for which the comparator has no counterpart ("no such source"). LoFTEE's reference ran with stand-ins for the MaxEntScan files its module bundles, on the Ensembl VEP image with samtools added. (d) The plugins each tool reads, of the 31 that vep-rs supports: vep-rs 31 and fastVEP 7; plugins vep-rs does not port are not counted and are listed by reason in the README. The plugin panel's 57 rows carry 32 row names for 31 plugins: CADD is scored under CADD and CADD (structural variants). A row's figure is summed over the sets of the row that the record holds for the tool, and the tools differ in them on 4 rows: dbNSFP GRCh37 synthetic (vep-rs 5, fastVEP 1), dbNSFP GRCh38 synthetic (vep-rs 9, fastVEP 1), SpliceAI GRCh37 synthetic (vep-rs 4, fastVEP 1) and SpliceAI GRCh38 synthetic (vep-rs 4, fastVEP 1). On an F1 axis a value is a place, not a length: every bar starts at the left edge of its plot, left of where an F1 of 0 stands, and ends at the place of its F1. A marker standing at the end of its bar means that the tool's raw and adjusted F1 are equal as drawn. On a comparator's adjusted F1 the record's note reads: "The adjusted figure sets aside the same documented classes as for vep-rs; where the engine produces none of their shapes it equals the raw one." fastVEP's adjusted F1 is left out where no record holds its adjusted counts. An F1 of exactly 1 has no place on a scale of 1 - F1: it is drawn where the axis reads 1, past a break in the axis line, and 66 of the 66 adjusted markers of vep-rs stand there. vep-rs matches Ensembl VEP exactly except for the documented divergences. Every value is in the records listed under Provenance records.

### Whole-genome concordance

Against Ensembl VEP 116.2. A population dataset is scored per chromosome and its shards pooled; each ClinVar release is one genome-wide file. The per-chromosome tables are in the [provenance records](#provenance-records).

| Engine | Dataset | Assembly | Shards | Raw F1 | Adjusted F1 | Variants | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.6.0 | gnomAD v4.1 | GRCh38 | 24 | 0.999998755 | 1.000000000 | 759,336,320 | 9,601,562,823 | 9,601,562,823 | 9,601,550,872 |
| vep-rs v0.6.0 | gnomAD v2.1.1 | GRCh37 | 23 | 0.999999005 | 1.000000000 | 261,942,336 | 885,664,902 | 885,664,902 | 885,664,021 |
| vep-rs v0.6.0 | 1000 Genomes Phase 3 | GRCh37 | 25 | 0.999999400 | 1.000000000 | 84,805,772 | 286,649,478 | 286,649,478 | 286,649,306 |
| vep-rs v0.6.0 | 1000 Genomes high-coverage | GRCh38 | 23 | 0.999990864 | 1.000000000 | 73,554,796 | 905,822,094 | 905,836,944 | 905,821,243 |
| vep-rs v0.6.0 | ClinVar | GRCh37 | 1 | 0.999929250 | 1.000000000 | 4,388,172 | 33,483,881 | 33,483,881 | 33,481,512 |
| vep-rs v0.6.0 | ClinVar | GRCh38 | 1 | 0.999905298 | 1.000000000 | 4,398,127 | 125,572,619 | 125,572,619 | 125,560,727 |
| vep-rs v0.6.0 | All six | GRCh37 and GRCh38 | 97 | 0.999996998 | 1.000000000 | 1,188,425,523 | 11,838,755,797 | 11,838,770,647 | 11,838,727,681 |
| fastVEP v0.5.0 | gnomAD v4.1 | GRCh38 | 24 | 0.999997730 | 0.999998149 | 759,336,320 | 9,601,562,823 | 9,601,562,823 | 9,601,541,024 |
| fastVEP v0.5.0 | gnomAD v2.1.1 | GRCh37 | 23 | 0.992171351 | 0.992171764 | 261,942,336 | 885,664,902 | 889,340,923 | 880,554,964 |
| fastVEP v0.5.0 | 1000 Genomes Phase 3 | GRCh37 | 25 | 0.991905131 | 0.991905307 | 84,805,772 | 286,649,478 | 287,933,151 | 284,965,729 |
| fastVEP v0.5.0 | 1000 Genomes high-coverage | GRCh38 | 23 | 0.998446204 | 0.998454606 | 73,554,796 | 905,822,094 | 905,859,484 | 904,433,297 |
| fastVEP v0.5.0 | ClinVar | GRCh37 | 1 | 0.994836380 | 0.994842092 | 4,388,172 | 33,483,881 | 33,742,985 | 33,439,866 |
| fastVEP v0.5.0 | ClinVar | GRCh38 | 1 | 0.999982066 | 0.999990699 | 4,398,127 | 125,572,619 | 125,572,619 | 125,570,367 |
| fastVEP v0.5.0 | All six | GRCh37 and GRCh38 | 97 | 0.999081303 | 0.999082429 | 1,188,425,523 | 11,838,755,797 | 11,844,011,985 | 11,830,505,247 |

### Whole-genome wall time

Seconds: a population dataset as the sum of its per-chromosome medians, each ClinVar release as its whole-file median. The figure in parentheses beside a comparator's time is that time divided by vep-rs's, taken only where both were measured under this protocol. The two right-hand columns are vep-rs on its largest chromosome, or on the whole file for ClinVar. With 5 measurements per cell the P5 and P95 are the cell's minimum and maximum.

- vep-rs and fastVEP ran the population datasets on `c9gd.8xlarge` (ARM Graviton5) and `c8id.8xlarge` (x86 Intel), the newest available class per architecture (Graviton5 against 8th-generation Intel), 32 vCPU with local NVMe.
- Ensembl VEP 116.2's figures are its single reference run per chromosome of the same inputs under the same protocol, on `c9gd.8xlarge` and `c8id.8xlarge` (the same processors with 64 GiB), the newest available class per architecture (Graviton5 against 8th-generation Intel), except gnomAD v4.1 chromosomes 1, 2 and 3 on `c9gd.12xlarge` (ARM, 48 vCPU / 96 GiB; peak memory 77.9, 76.4 and 80.2 GiB) and `c8id.12xlarge` (x86, 48 vCPU / 96 GiB; peak memory 73.2, 75.7 and 77.3 GiB), over the standard class's 64 GiB. The concordance tables compare against this run's output.
- The ClinVar rows are the released version's whole-file medians on `c9gd.8xlarge` and `c8id.8xlarge`, the newest available class per architecture (Graviton5 against 8th-generation Intel); the comparators' times on them were measured on the same classes under the same protocol, fastVEP v0.5.0 as the median over 5 machines and Ensembl VEP 116.2 as one machine's timed run.

| Dataset | Assembly | Variants | Arch | vep-rs v0.6.0 (s) | fastVEP v0.5.0 (s) | Ensembl VEP 116.2 (s) | vep-rs, largest cell median (s) | P5 to P95 (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gnomAD v4.1 | GRCh38 | 759,336,320 | ARM | 591 | 4,902 (8.3x) | 255,394 (432x) | 55.32 (chromosome 2) | 50.71 to 63.31 |
| gnomAD v4.1 | GRCh38 | 759,336,320 | x86 | 799 | 6,079 (7.6x) | 265,212 (332x) | 76.01 (chromosome 2) | 71.54 to 79.49 |
| gnomAD v2.1.1 | GRCh37 | 261,942,336 | ARM | 207 | 3,918 (18.9x) | 88,002 (425x) | 17.10 (chromosome 2) | 16.77 to 18.00 |
| gnomAD v2.1.1 | GRCh37 | 261,942,336 | x86 | 377 | 4,397 (11.7x) | 103,809 (275x) | 31.88 (chromosome 2) | 30.92 to 31.99 |
| 1000 Genomes Phase 3 | GRCh37 | 84,805,772 | ARM | 18 | 113 (6.4x) | 7,209 (407x) | 1.46 (chromosome 2) | 1.45 to 1.46 |
| 1000 Genomes Phase 3 | GRCh37 | 84,805,772 | x86 | 22 | 142 (6.4x) | 7,039 (317x) | 1.84 (chromosome 2) | 1.69 to 1.86 |
| 1000 Genomes high-coverage | GRCh38 | 73,554,796 | ARM | 39 | 193 (5.0x) | 15,817 (410x) | 3.40 (chromosome 1) | 3.21 to 3.43 |
| 1000 Genomes high-coverage | GRCh38 | 73,554,796 | x86 | 52 | 262 (5.0x) | 17,660 (340x) | 4.39 (chromosome 1) | 4.29 to 4.41 |
| ClinVar | GRCh37 | 4,388,172 | ARM | 1.90 | 6.85 (3.6x) | 719.22 (379x) | 1.90 (whole file) | 1.87 to 1.91 |
| ClinVar | GRCh37 | 4,388,172 | x86 | 2.78 | 9.84 (3.5x) | 837.80 (301x) | 2.78 (whole file) | 2.75 to 2.80 |
| ClinVar | GRCh38 | 4,398,127 | ARM | 8.02 | 12.84 (1.6x) | 2,543.75 (317x) | 8.02 (whole file) | 7.48 to 8.93 |
| ClinVar | GRCh38 | 4,398,127 | x86 | 9.26 | 17.97 (1.9x) | 3,093.35 (334x) | 9.26 (whole file) | 9.23 to 9.59 |
| All four population datasets | GRCh37 and GRCh38 | 1,179,639,224 | ARM | 854 | 9,125 (10.7x) | 366,422 (429x) | | |
| All four population datasets | GRCh37 and GRCh38 | 1,179,639,224 | x86 | 1,251 | 10,880 (8.7x) | 393,720 (315x) | | |

### Structural variants

The structural-variant per-VCF sets (16 files per assembly, the files' counts pooled), against Ensembl VEP 116.2. fastVEP v0.5.0's rows are its scoring against the same reference, its adjusted figure differing from its raw one where the documented classes reach its rows (GRCh37 and GRCh38).

| Engine | Assembly | Raw F1 | Adjusted F1 | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.6.0 | GRCh37 | 0.963822119 | 1.000000000 | 317,047 | 337,773 | 315,565 |
| vep-rs v0.6.0 | GRCh38 | 0.905536410 | 1.000000000 | 1,455,318 | 1,734,800 | 1,444,384 |
| fastVEP v0.5.0 | GRCh37 | 0.245961547 | 0.246102670 | 317,047 | 286,340 | 74,205 |
| fastVEP v0.5.0 | GRCh38 | 0.228874968 | 0.228904209 | 1,455,318 | 1,315,802 | 317,120 |

### Structural-variant wall time

The released version on `c9gd.8xlarge` and `c8id.8xlarge`, the newest available class per architecture (Graviton5 against 8th-generation Intel), 16 files per assembly run as 16 processes: a start-up-bound regime, excluded from every speed-up above. With 5 measurements per cell the P5 and P95 are the cell's minimum and maximum.

| Arch | Dataset | Assembly | v0.6.0 median (s) | P5 to P95 (s) |
| --- | --- | --- | --- | --- |
| ARM | SV per-VCF (16 files) | GRCh37 | 0.56 | 0.54 to 0.58 |
| ARM | SV per-VCF (16 files) | GRCh38 | 1.41 | 1.39 to 1.48 |
| x86 | SV per-VCF (16 files) | GRCh37 | 0.91 | 0.90 to 0.99 |
| x86 | SV per-VCF (16 files) | GRCh38 | 2.14 | 2.11 to 2.22 |

### Plugins

Against Ensembl VEP 116.2, on the plugin tuple stated under [Method](#method): the plugins verified on real data, each on the genome-wide ClinVar dataset of the assembly and on the whole-genome population datasets its Data cell names, the row's counts summed over those datasets, and TranscriptAnnotator on its committed tables (`GRCh37-plugin-transcript-annotator` and `GRCh38-plugin-transcript-annotator` under `tests/golden/116/`).

| Engine | Plugin or source | Data | Assembly | Raw F1 | Adjusted F1 | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.6.0 | AlphaMissense | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.999898999 | 1.000000000 | 27,293,920 | 27,299,434 | 27,293,920 |
| fastVEP v0.5.0 (per feature) | AlphaMissense | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.477551794 | not applied: `alphamissense_two_digit_transcript_version` | 27,293,920 | 34,675,396 | 14,796,779 |
| fastVEP v0.5.0 (per allele) | AlphaMissense | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.698912693 | 0.699069738 | 8,724,832 | 4,785,153 | 4,721,150 |
| vep-rs v0.6.0 | AlphaMissense | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.999989474 | 1.000000000 | 182,213,448 | 182,217,284 | 182,213,448 |
| fastVEP v0.5.0 (per feature) | AlphaMissense | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.678869469 | not applied: `alphamissense_two_digit_transcript_version` | 182,213,448 | 350,983,024 | 180,985,403 |
| fastVEP v0.5.0 (per allele) | AlphaMissense | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.991301689 | not applied: `alphamissense_two_digit_transcript_version` | 13,936,401 | 13,986,606 | 13,840,062 |
| vep-rs v0.6.0 | AncestralAllele | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 1.000000000 | 1.000000000 | 1,138,649,670 | 1,138,649,670 | 1,138,649,670 |
| fastVEP v0.5.0 | AncestralAllele | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | AncestralAllele | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 1.000000000 | 1.000000000 | 9,601,363,499 | 9,601,363,499 | 9,601,363,499 |
| fastVEP v0.5.0 | AncestralAllele | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | Blosum62 | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.999999933 | 1.000000000 | 14,821,388 | 14,821,386 | 14,821,386 |
| fastVEP v0.5.0 | Blosum62 | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | Blosum62 | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 1.000000000 | 1.000000000 | 96,829,002 | 96,829,002 | 96,829,002 |
| fastVEP v0.5.0 | Blosum62 | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | Conservation | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 1.000000000 | 1.000000000 | 9,867,420,929 | 9,867,420,929 | 9,867,420,929 |
| fastVEP v0.5.0 | Conservation | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | DosageSensitivity | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 1.000000000 | 1.000000000 | 1,798,317,572 | 1,798,317,572 | 1,798,317,572 |
| fastVEP v0.5.0 | DosageSensitivity | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | DosageSensitivity | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.999999663 | 1.000000000 | 14,187,431,686 | 14,187,441,252 | 14,187,431,686 |
| fastVEP v0.5.0 | DosageSensitivity | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | Downstream | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.788775391 | 1.000000000 | 1,461,951 | 1,471,817 | 1,157,042 |
| fastVEP v0.5.0 | Downstream | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | Downstream | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.788627126 | 1.000000000 | 10,273,147 | 10,339,233 | 8,127,741 |
| fastVEP v0.5.0 | Downstream | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | gnomADc | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 1.000000000 | 1.000000000 | 13,255,991,683 | 13,255,991,683 | 13,255,991,683 |
| fastVEP v0.5.0 | gnomADc | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | gnomADc | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.999974559 | 1.000000000 | 15,951,642,432 | 15,952,454,100 | 15,951,642,432 |
| fastVEP v0.5.0 | gnomADc | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | gnomADMt | ClinVar | GRCh38 | 1.000000000 | 1.000000000 | 297,545 | 297,545 | 297,545 |
| fastVEP v0.5.0 | gnomADMt | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | GWAS | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.992879519 | 1.000000000 | 47,507,412 | 46,907,131 | 46,871,133 |
| fastVEP v0.5.0 | GWAS | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | HGVSIntronOffset | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.999994243 | 1.000000000 | 1,645,147,996 | 1,645,157,824 | 1,645,143,439 |
| fastVEP v0.5.0 | HGVSIntronOffset | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | HGVSIntronOffset | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.999989214 | 1.000000000 | 16,960,725,510 | 16,960,830,466 | 16,960,595,049 |
| fastVEP v0.5.0 | HGVSIntronOffset | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | LOEUF | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 1.000000000 | 1.000000000 | 736,165,178 | 736,165,178 | 736,165,178 |
| fastVEP v0.5.0 (per feature) | LOEUF | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.372314570 | 0.372314570 | 736,165,178 | 935,154,519 | 311,128,337 |
| fastVEP v0.5.0 (per allele) | LOEUF | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.503374478 | 0.503374478 | 159,738,518 | 185,599,022 | 86,917,052 |
| vep-rs v0.6.0 | LOEUF | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 1.000000000 | 1.000000000 | 6,310,931,246 | 6,310,931,246 | 6,310,931,246 |
| fastVEP v0.5.0 (per feature) | LOEUF | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.126740983 | 0.126740983 | 6,310,931,246 | 7,817,653,301 | 895,335,349 |
| fastVEP v0.5.0 (per allele) | LOEUF | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.152210193 | 0.152210193 | 359,558,712 | 437,457,703 | 60,657,011 |
| vep-rs v0.6.0 | LoFTEE | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.701585606 | 1.000000000 | 22,020,844 | 21,990,039 | 15,438,701 |
| fastVEP v0.5.0 | LoFTEE | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | LoFtool | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 1.000000000 | 1.000000000 | 772,751,242 | 772,751,242 | 772,751,242 |
| fastVEP v0.5.0 | LoFtool | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | LoFtool | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 1.000000000 | 1.000000000 | 6,140,527,244 | 6,140,527,244 | 6,140,527,244 |
| fastVEP v0.5.0 | LoFtool | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | NearestExonJB | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.999984102 | 1.000000000 | 48,368,873 | 48,370,079 | 48,368,707 |
| fastVEP v0.5.0 | NearestExonJB | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | NearestExonJB | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.999998706 | 1.000000000 | 445,543,072 | 445,543,287 | 445,542,603 |
| fastVEP v0.5.0 | NearestExonJB | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | NMD | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.998091399 | 1.000000000 | 318,050 | 318,018 | 317,427 |
| fastVEP v0.5.0 | NMD | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | NMD | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.997522027 | 1.000000000 | 2,074,273 | 2,074,279 | 2,069,136 |
| fastVEP v0.5.0 | NMD | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | OpenTargets | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 1.000000000 | 1.000000000 | 237,997,475 | 237,997,475 | 237,997,475 |
| fastVEP v0.5.0 | OpenTargets | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | pLI | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 1.000000000 | 1.000000000 | 896,438,668 | 896,438,668 | 896,438,668 |
| fastVEP v0.5.0 (per feature) | pLI | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.335976085 | 0.335976085 | 896,438,668 | 905,826,126 | 302,758,935 |
| fastVEP v0.5.0 (per allele) | pLI | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.343099987 | 0.343099987 | 173,189,135 | 175,944,985 | 59,893,956 |
| vep-rs v0.6.0 | pLI | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 1.000000000 | 1.000000000 | 7,149,538,152 | 7,149,538,152 | 7,149,538,152 |
| fastVEP v0.5.0 (per feature) | pLI | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.321559950 | 0.321559950 | 7,149,538,152 | 7,345,403,197 | 2,330,496,309 |
| fastVEP v0.5.0 (per allele) | pLI | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.336761364 | 0.336761364 | 398,673,174 | 402,875,278 | 134,965,275 |
| vep-rs v0.6.0 | REVEL | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.999667478 | 1.000000000 | 12,333,981 | 12,326,039 | 12,325,910 |
| fastVEP v0.5.0 (per feature) | REVEL | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.491066649 | 0.491066649 | 12,333,981 | 36,325,361 | 11,947,490 |
| fastVEP v0.5.0 (per allele) | REVEL | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.955212957 | 0.955212957 | 4,673,814 | 4,855,192 | 4,551,115 |
| vep-rs v0.6.0 | REVEL | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.999615257 | 1.000000000 | 15,521,504 | 15,509,565 | 15,509,565 |
| fastVEP v0.5.0 (per feature) | REVEL | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.141474058 | 0.141474058 | 15,521,504 | 197,337,209 | 15,056,993 |
| fastVEP v0.5.0 (per allele) | REVEL | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.925996734 | 0.925996734 | 6,893,122 | 7,616,035 | 6,717,716 |
| vep-rs v0.6.0 | TSSDistance | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 1.000000000 | 1.000000000 | 125,704,335 | 125,704,335 | 125,704,335 |
| fastVEP v0.5.0 | TSSDistance | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | TSSDistance | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 1.000000000 | 1.000000000 | 988,520,633 | 988,520,633 | 988,520,633 |
| fastVEP v0.5.0 | TSSDistance | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | UTRAnnotator | ClinVar, gnomAD v2.1.1 and 1000 Genomes Phase 3 | GRCh37 | 0.968210791 | 1.000000000 | 9,481,803 | 9,478,279 | 9,178,678 |
| fastVEP v0.5.0 | UTRAnnotator | ClinVar | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | UTRAnnotator | ClinVar, gnomAD v4.1 and 1000 Genomes high-coverage | GRCh38 | 0.972618929 | 1.000000000 | 128,606,861 | 128,571,761 | 125,068,398 |
| fastVEP v0.5.0 | UTRAnnotator | ClinVar | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | TranscriptAnnotator | committed table (1 set) | GRCh37 | 1.000000000 | 1.000000000 | 18 | 18 | 18 |
| fastVEP v0.5.0 | TranscriptAnnotator | committed table | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | TranscriptAnnotator | committed table (4 sets) | GRCh38 | 1.000000000 | 1.000000000 | 117 | 117 | 117 |
| fastVEP v0.5.0 | TranscriptAnnotator | committed table | GRCh38 | no such source | | | | |

#### Synthetic fixtures (data that cannot be redistributed)

Every row is marked synthetic: the plugin's data cannot be redistributed, so it is scored on the synthetic fixture of its assembly under `tests/golden/116/`, its counts summed over the fixture's reference sets.

| Engine | Plugin or source | Data | Assembly | Raw F1 | Adjusted F1 | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.6.0 | AVI | synthetic fixture (1 set) | GRCh38 | 1.000000000 | 1.000000000 | 414 | 414 | 414 |
| fastVEP v0.5.0 | AVI | synthetic fixture | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | BayesDel | synthetic fixture (1 set) | GRCh37 | 1.000000000 | 1.000000000 | 118 | 118 | 118 |
| fastVEP v0.5.0 | BayesDel | synthetic fixture | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | BayesDel | synthetic fixture (1 set) | GRCh38 | 1.000000000 | 1.000000000 | 14 | 14 | 14 |
| fastVEP v0.5.0 | BayesDel | synthetic fixture | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | CADD | synthetic fixture (3 sets) | GRCh37 | 1.000000000 | 1.000000000 | 708 | 708 | 708 |
| fastVEP v0.5.0 | CADD | synthetic fixture | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | CADD | synthetic fixture (4 sets) | GRCh38 | 1.000000000 | 1.000000000 | 112 | 112 | 112 |
| fastVEP v0.5.0 | CADD | synthetic fixture | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | CADD (structural variants) | synthetic fixture (3 sets) | GRCh38 | 1.000000000 | 1.000000000 | 120 | 120 | 120 |
| fastVEP v0.5.0 | CADD (structural variants) | synthetic fixture | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | ClinPred | synthetic fixture (2 sets) | GRCh37 | 1.000000000 | 1.000000000 | 78 | 78 | 78 |
| fastVEP v0.5.0 | ClinPred | synthetic fixture | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | ClinPred | synthetic fixture (2 sets) | GRCh38 | 1.000000000 | 1.000000000 | 14 | 14 | 14 |
| fastVEP v0.5.0 | ClinPred | synthetic fixture | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | dbNSFP | synthetic fixture (5 sets) | GRCh37 | 1.000000000 | 1.000000000 | 697 | 697 | 697 |
| fastVEP v0.5.0 (per feature) | dbNSFP | synthetic fixture (1 set) | GRCh37 | 0.012578616 | 0.012578616 | 88 | 230 | 2 |
| fastVEP v0.5.0 (per allele) | dbNSFP | synthetic fixture (1 set) | GRCh37 | 0.105263158 | 0.105263158 | 16 | 22 | 2 |
| vep-rs v0.6.0 | dbNSFP | synthetic fixture (9 sets) | GRCh38 | 1.000000000 | 1.000000000 | 4,670 | 4,670 | 4,670 |
| fastVEP v0.5.0 (per feature) | dbNSFP | synthetic fixture (1 set) | GRCh38 | 0.062717770 | 0.062717770 | 126 | 448 | 18 |
| fastVEP v0.5.0 (per allele) | dbNSFP | synthetic fixture (1 set) | GRCh38 | 0.500000000 | 0.500000000 | 38 | 34 | 18 |
| vep-rs v0.6.0 | dbscSNV | synthetic fixture (2 sets) | GRCh37 | 1.000000000 | 1.000000000 | 188 | 188 | 188 |
| fastVEP v0.5.0 | dbscSNV | synthetic fixture | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | dbscSNV | synthetic fixture (2 sets) | GRCh38 | 1.000000000 | 1.000000000 | 1,376 | 1,376 | 1,376 |
| fastVEP v0.5.0 | dbscSNV | synthetic fixture | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | MaxEntScan | synthetic fixture (3 sets) | GRCh37 | 0.995277957 | 1.000000000 | 4,659 | 4,659 | 4,637 |
| fastVEP v0.5.0 | MaxEntScan | synthetic fixture | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | MaxEntScan | synthetic fixture (3 sets) | GRCh38 | 0.999236164 | 1.000000000 | 14,401 | 14,401 | 14,390 |
| fastVEP v0.5.0 | MaxEntScan | synthetic fixture | GRCh38 | no such source | | | | |
| vep-rs v0.6.0 | PrimateAI | synthetic fixture (1 set) | GRCh37 | 1.000000000 | 1.000000000 | 88 | 88 | 88 |
| fastVEP v0.5.0 (per feature) | PrimateAI | synthetic fixture (1 set) | GRCh37 | 1.000000000 | 1.000000000 | 88 | 88 | 88 |
| fastVEP v0.5.0 (per allele) | PrimateAI | synthetic fixture (1 set) | GRCh37 | 1.000000000 | 1.000000000 | 11 | 11 | 11 |
| vep-rs v0.6.0 | PrimateAI | synthetic fixture (1 set) | GRCh38 | 1.000000000 | 1.000000000 | 11 | 11 | 11 |
| fastVEP v0.5.0 (per feature) | PrimateAI | synthetic fixture (1 set) | GRCh38 | 1.000000000 | 1.000000000 | 11 | 11 | 11 |
| fastVEP v0.5.0 (per allele) | PrimateAI | synthetic fixture (1 set) | GRCh38 | 1.000000000 | 1.000000000 | 11 | 11 | 11 |
| vep-rs v0.6.0 | SpliceAI | synthetic fixture (4 sets) | GRCh37 | 0.928977273 | 1.000000000 | 2,002 | 2,222 | 1,962 |
| fastVEP v0.5.0 (per feature) | SpliceAI | synthetic fixture (1 set) | GRCh37 | 0.810256410 | 0.882681564 | 91 | 104 | 79 |
| fastVEP v0.5.0 (per allele) | SpliceAI | synthetic fixture (1 set) | GRCh37 | 0.769230769 | 0.833333333 | 26 | 26 | 20 |
| vep-rs v0.6.0 | SpliceAI | synthetic fixture (4 sets) | GRCh38 | 0.935861684 | 1.000000000 | 1,738 | 1,848 | 1,678 |
| fastVEP v0.5.0 (per feature) | SpliceAI | synthetic fixture (1 set) | GRCh38 | 0.759036145 | 0.887323944 | 79 | 87 | 63 |
| fastVEP v0.5.0 (per allele) | SpliceAI | synthetic fixture (1 set) | GRCh38 | 0.760000000 | 0.826086957 | 25 | 25 | 19 |
| vep-rs v0.6.0 | SpliceVault | synthetic fixture (2 sets) | GRCh37 | 1.000000000 | 1.000000000 | 189 | 189 | 189 |
| fastVEP v0.5.0 | SpliceVault | synthetic fixture | GRCh37 | no such source | | | | |
| vep-rs v0.6.0 | SpliceVault | synthetic fixture (2 sets) | GRCh38 | 1.000000000 | 1.000000000 | 245 | 245 | 245 |
| fastVEP v0.5.0 | SpliceVault | synthetic fixture | GRCh38 | no such source | | | | |

#### Plugins each tool reads

A tool reads a plugin when it writes at least one per-feature tuple of it. The counts are of the 31 plugins vep-rs supports; the plugins it does not port are not counted.

| Engine | Plugins read, of 31 | On ClinVar | On synthetic fixtures | On a committed table |
| --- | --- | --- | --- | --- |
| vep-rs v0.6.0 | 31 | 20 | 10 | 1 |
| fastVEP v0.5.0 | 7 | 4 | 3 | 0 |

#### Plugins not ported, by reason

The support table of [docs/plugins.md](docs/plugins.md) lists 87 Ensembl VEP plugins; vep-rs does not port these 56.

| Reason | Plugins | Count |
| --- | --- | --- |
| removed from or absent from Ensembl VEP's release 116 collection, or resting on a retired resource | AVADA, CCDSFilter, CSN, DAS, FlagLRG, Gwava, GXA, HGVSReferenceBase, neXtProt, NonSynonymousFilter, PostGAP, RankFilter, RefSeqHGVS, RiboseqORFs, StructuralVariantOverlap, SubsetVCF | 16 |
| needing a live database, a web service or an external program | FATHMM, GeneSplicer, GO, LD, LocalID, LOVD, mutfunc, NearestGene, Paralogues, Phenotypes, PolyPhen_SIFT, ProtVar, SameCodon | 13 |
| data for GRCh37 only | CAPICE, FATHMM_MKL, FunMotifs, Geno2MP, MPC, MTR, VARITY | 7 |
| specialised data, rarely used | Enformer, IntAct, MaveDB, MechPredict, PhenotypeOrthologous, ReferenceQuality, satMutMPRA | 7 |
| writing files or reports, or rewriting or filtering output | DeNovo, Draw, G2P, ProteinSeqs, SingleLetterAA | 5 |
| commercial or restricted licence | GeneBe, Mastermind, PON_P2, PromoterAI | 4 |
| covered by something vep-rs has | Carol, Condel, EVE, SpliceRegion | 4 |
| Total | | 56 |

### Method

- **Terms.** A variant is one VCF record. A tuple is one output row reduced to its Location, Allele, Feature, Feature_type and Consequence set. F1 pools a dataset's chromosomes by summed counts. HGVS notation is outside every F1 above.
- **Concordance.** One run per dataset of a build of the released version, against Ensembl VEP 116.2's output on the same inputs, with the same comparator, cache and reference FASTA. F1 is deterministic per binary and independent of architecture.
- **Adjusted F1.** The tuples it sets aside are the documented classes of [docs/intended-divergences.md](docs/intended-divergences.md). Every tuple outside those classes is matched.
- **fastVEP v0.5.0** was scored on the same Ensembl VEP output, with the same comparator, cache and reference and the same classes set aside. It was run with `--exclude-readthrough --exclude-artifact`, the two flags its documentation gives for a run comparable with a run on Ensembl VEP's cache. Its rows state its own tuple counts. Its adjusted figure differs from its raw one where the documented classes reach its whole-genome rows (gnomAD v4.1 on GRCh38; gnomAD v2.1.1 on GRCh37; 1000 Genomes Phase 3 on GRCh37; 1000 Genomes high-coverage on GRCh38; ClinVar on GRCh37; ClinVar on GRCh38; All six on both assemblies).
- **vep-rs v0.5.0, re-scored.** vep-rs v0.5.0 was re-scored against the same Ensembl VEP 116.2 reference (`previous_release_concordance`); on ClinVar GRCh37, ClinVar GRCh38, SV per-VCF (16 files) GRCh37 and SV per-VCF (16 files) GRCh38 vep-rs v0.6.0 matches 1,668, 8,432, 72 and 358 fewer tuples than v0.5.0's re-scoring (33,481,512, 125,560,727, 315,565 and 1,444,384 against 33,483,180, 125,569,159, 315,637 and 1,444,742).
- **Wall time.** For vep-rs and fastVEP, the median of 5 independent machine measurements per cell on both architectures, with `--fork 16` on sites-only inputs; for Ensembl VEP, its reference run per chromosome on the population datasets and, on ClinVar, a timed run on the same classes. Each measurement runs on a fresh instance after a discarded warmup, whose output is flushed and deleted before the timed run.
- **Speed against vep-rs v0.5.0.** vep-rs v0.5.0's archived binaries were timed on the four population datasets on the same day as vep-rs v0.6.0, on the same instance classes and with the same number of machines per cell (`previous_release_wall_time`). vep-rs v0.6.0's summed per-chromosome medians divided by vep-rs v0.5.0's: gnomAD v4.1, 0.85x on ARM and 0.76x on x86; gnomAD v2.1.1, 0.97x on ARM and 0.93x on x86; 1000 Genomes Phase 3, 0.66x on ARM and 0.51x on x86; 1000 Genomes high-coverage, 0.64x on ARM and 0.57x on x86. A figure below 1 is a shorter wall time for vep-rs v0.6.0.
- **Plugins.** A plugin tuple is one populated plugin field of one output row: its Location, Allele, Feature, field name and value. Raw F1 compares a plugin's tuples as the two engines wrote them; adjusted F1 sets aside the tuples of the documented classes of [docs/intended-divergences.md](docs/intended-divergences.md) and nothing else. The 31 plugins vep-rs supports are scored: a plugin verified on real data on the genome-wide ClinVar dataset of each assembly it supports and on the whole-genome population datasets of that assembly its row names, the row's counts summed over those datasets, with the publishers' files; a plugin whose data cannot be redistributed on the synthetic fixtures under `tests/golden/116/`, each such row marked synthetic and scored on every reference set of its fixture; TranscriptAnnotator, which reads a table of the user's own, on the committed table of the documented layout. The tuple counts stand beside each F1. The Ensembl VEP reference for LoFTEE ran LoFTEE's module with stand-ins for the MaxEntScan files it bundles, splice predictions off, on the Ensembl VEP image with samtools added. Plugins that vep-rs does not port are not counted; [Plugins not ported, by reason](#plugins-not-ported-by-reason) lists them. fastVEP v0.5.0 ran as `run plugins: --output-format tab --exclude-readthrough --exclude-artifact --sa-dir <databases>/plugins; run pli: --output-format tab --exclude-readthrough --exclude-artifact --sa-dir <databases>/pli` on databases it built from the data Ensembl VEP read: for its rows not marked synthetic the publishers' files (REVEL from the publisher's comma-separated archive, LOEUF and pLI from the constraint tables the plugins read), and for its rows marked synthetic the same synthetic fixtures. A plugin for which fastVEP has no counterpart is marked "no such source". fastVEP attaches a value to an allele and writes it on every feature of that allele, so its rows state a per-feature figure, on the tuple above, and a per-allele figure, on the tuple without the Feature. Where a documented class that is not applied would set a tuple of one of fastVEP's rows aside, the row states no adjusted F1 and names the class in its place. Two numbers are equal when they differ by at most 0.000001 of the larger of the two and by at most 0.000001 outright, both conditions; the record also counts the values that are equal under that rule and differ as text.

### Provenance records

Every value above is recorded in these files:

- [`docs/concordance-provenance/2026-11-19-release-v0.6.0.json`](docs/concordance-provenance/2026-11-19-release-v0.6.0.json): the ClinVar and structural-variant cells, the per-column field agreement, the compiler and build inputs of the measured binaries, the reference sets, the same-day comparison against v0.5.0, and the whole-genome concordance per chromosome.
- [`docs/concordance-provenance/2026-11-19-population-v0.6.0.json`](docs/concordance-provenance/2026-11-19-population-v0.6.0.json): every chromosome's wall-time median, P5 and P95 for all three engines, with the 5 wall times behind vep-rs's and fastVEP's, and fastVEP's concordance per chromosome; fastVEP ran with `--exclude-readthrough --exclude-artifact`.
- [`docs/concordance-provenance/2026-11-18-comparator-fastvep-v0.5.0.json`](docs/concordance-provenance/2026-11-18-comparator-fastvep-v0.5.0.json): fastVEP v0.5.0's ClinVar and structural-variant figures with `--exclude-readthrough --exclude-artifact`, scored on the same Ensembl VEP output with the same comparator, cache and reference, the same classes set aside.
- [`docs/concordance-provenance/2026-11-12-comparator-walltime-v0.6.0.json`](docs/concordance-provenance/2026-11-12-comparator-walltime-v0.6.0.json): the comparators' ClinVar wall times, measured on the same classes as vep-rs under the same protocol; fastVEP ran with `--exclude-readthrough --exclude-artifact`.
- [`docs/concordance-provenance/2026-11-19-plugins-v0.6.0.json`](docs/concordance-provenance/2026-11-19-plugins-v0.6.0.json): every plugin's tuple counts per engine, assembly and dataset or fixture reference set, the Ensembl VEP references they were scored against, each tool's version, binaries and arguments, the plugins each tool reads and the plugins not ported, by reason.
- [`docs/concordance-provenance/2026-10-25-release-v0.5.0.json`](docs/concordance-provenance/2026-10-25-release-v0.5.0.json): vep-rs v0.5.0's own record, with the figures it published.
- [`docs/concordance-provenance/2026-10-25-population-v0.5.0.json`](docs/concordance-provenance/2026-10-25-population-v0.5.0.json): vep-rs v0.5.0's whole-genome measurement.

### Reproducing the figures

- `scripts/concordance/run_clone_measurement.sh` re-measures the ClinVar and structural-variant cells on your own machines.
- `scripts/concordance/run_concordance.sh` times both engines on a directory of your own VCFs.
- [scripts/README.md](scripts/README.md#reproducing-the-published-concordance) is the runbook for reproducing the concordance numbers.

## Documentation

The documentation site, https://natera-open-source.github.io/vep-rs/, describes the latest release; its pages are the `docs/` tree of this repository, in six sections.

- **Getting started**: [Installation](docs/getting-started/installation.md), the transcript cache and [the move from Ensembl VEP](docs/getting-started/coming-from-ensembl-vep.md); [Examples](docs/examples.md), the common shapes of a run.
- **Reference**: [every flag](docs/cli-reference.md), with the option rules, the flags refused or accepted without effect and the exit-status contract; [every output format](docs/output-formats.md); [every plugin](docs/plugins.md).
- **Concordance with Ensembl VEP**: [the divergence classes](docs/intended-divergences.md) the adjusted F1 sets aside.
- **Benchmarks**: the [Current release](#current-release-v060) figures above, also on the site's [Benchmarks](docs/benchmarks/current-release.md) page, beside [the paper's figures](docs/published-figures.md).
- **Developers**: [Writing a plugin](docs/developers/writing-a-plugin.md).
- **Project**: [Citing vep-rs](docs/project/citing.md) and the [changelog](CHANGELOG.md).

## Crates

One row per crate under `crates/`, in directory order; each links the crate's own README.

| Crate | Purpose | Kind |
| --- | --- | --- |
| [vep-builtins](crates/vep-builtins/README.md) | The built-in plugins: the tabix annotation engine, allele matching and the clinical plugins such as LoFTEE | library |
| [vep-cache-builder](crates/vep-cache-builder/README.md) | Builds the JSON transcript cache from an Ensembl GFF3, protein FASTA and variation VCF release | binary |
| [vep-cache-converter](crates/vep-cache-converter/README.md) | Copies a flat JSON transcript cache into the per-chromosome runtime layout, validating each region, and converts tabix plugin data to the binary annotation format | binary |
| [vep-cli](crates/vep-cli/README.md) | The `vep` command line: argument parsing, the runner and the output pipeline | library and binary |
| [vep-core](crates/vep-core/README.md) | The shared types: variants, transcripts, consequences and coordinates | library |
| [vep-effects](crates/vep-effects/README.md) | The consequence-calculation engine, with HGVS and protein prediction | library |
| [vep-fasta](crates/vep-fasta/README.md) | Indexed FASTA access for the reference genome | library |
| [vep-io](crates/vep-io/README.md) | Input parsers (Ensembl default, HGVS, SPDI, region, variant identifier) with format detection, and the output writers (VEP, VCF, JSON, Tab, Parquet) | library |
| [vep-plugin](crates/vep-plugin/README.md) | The plugin trait, the dynamic-library loader and the registry | library |

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for build, test and pull-request guidelines, and [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) for community expectations. To report a security vulnerability, follow [SECURITY.md](SECURITY.md) rather than filing a public issue.

## Authors

vep-rs was created and is maintained by:

- **Matthew Porter**, Natera, Inc.
- **Robert Borkowski**, Natera, Inc.

## Citing vep-rs

If you use vep-rs in published work, cite the paper:

> Porter M, Borkowski R. _vep-rs: high-throughput Rust variant annotation with population-scale concordance to Ensembl VEP._ bioRxiv (2026). [doi:10.64898/2026.09.22.753614](https://doi.org/10.64898/2026.09.22.753614).

The paper's figures are in [docs/published-figures.md](docs/published-figures.md) and re-derive from `manuscript/data/*.csv` in this tree ([manuscript/data/README.md](manuscript/data/README.md)). The current release's own measurements are in [Current release](#current-release-v060).

The software itself is archived at Zenodo, and that archive is the citation for the code and the released data files:

> Porter M, Borkowski R. _vep-rs_, version 0.1.0. Zenodo (2026). [doi:10.5281/zenodo.22837897](https://doi.org/10.5281/zenodo.22837897)

`CITATION.cff` at the repository root carries the software citation in machine-readable form, with the all-versions DOI ([10.5281/zenodo.22837896](https://doi.org/10.5281/zenodo.22837896)) and the archived version's own. GitHub's "Cite this repository" button and citation managers pick it up directly.

## Acknowledgments

vep-rs is an independent Rust implementation of the variant effect prediction algorithm. It is not affiliated with or endorsed by EMBL-EBI or the Ensembl project.

The algorithm is described in:

> McLaren W, Gil L, Hunt SE, et al. _The Ensembl Variant Effect Predictor._ Genome Biology 17, 122 (2016). [doi:10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)

## License

Apache-2.0. See [LICENSE](LICENSE).

## Disclaimer

vep-rs is distributed under the Apache License 2.0 on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied, including, without limitation, any warranties
or conditions of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE.
You are solely responsible for determining the appropriateness of using vep-rs and assume any
risks associated with that use. See sections 7 and 8 of the LICENSE for the complete disclaimer
of warranty and limitation of liability.

vep-rs is released for research and informational use. It has not been validated, cleared, or
approved by any regulatory authority for clinical use, and its output must not be relied upon
for patient care. Any clinical or diagnostic application of this software is the sole
responsibility of the party undertaking it.
