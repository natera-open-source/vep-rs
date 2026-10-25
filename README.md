# vep-rs

[![CI](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![MSRV](https://img.shields.io/badge/MSRV-1.89-orange.svg)](#installation)

Documentation: https://natera-open-source.github.io/vep-rs/

A fast, memory-efficient reimplementation of the Ensembl Variant Effect Predictor (VEP) in Rust. It takes VEP's inputs and flags, writes VEP's consequence vocabulary in VEP's output formats, and annotates population-scale callsets two orders of magnitude faster than the Perl implementation.

In this README, "Ensembl VEP" is the reference implementation, `ensemblorg/ensembl-vep`, written in Perl; unqualified, it means release 116.2. The accompanying paper calls the same engine "VEP".

## Why vep-rs instead of Ensembl VEP or fastVEP

|                                  | vep-rs v0.5.0                                                                       |
| -------------------------------- | ----------------------------------------------------------------------------------- |
| **Speed vs Ensembl VEP 116.2**   | 396x faster on ARM (Graviton5), 254x on x86; 162x to 411x per dataset               |
| **Speed vs fastVEP v0.4.0**      | 19.6x faster on ARM, 11.8x on x86; 9.2x to 22.4x per dataset                        |
| **All of gnomAD v4.1**           | 759,336,320 variants, 24 chromosomes: 624 s on one ARM machine, 1,012 s on x86      |
| **Peak memory**                  | At most 4.0 GB RSS on any chromosome                                                |
| **Concordance with Ensembl VEP** | Adjusted F1 1.000000000 on all six whole-genome datasets, raw F1 0.999998408 pooled |
| **Plugins**                      | 10 built-in databases, plus C-ABI dynamic-library plugins                           |
| **Output formats**               | VEP default, VCF, JSON, Tab, Parquet                                                |
| **Threading**                    | `--fork N`; the default uses every logical CPU up to 32                             |

The speed figures cover all of four population datasets (1,179,639,224 variants) under one protocol: the sum of per-chromosome medians at 16 threads, over 5 machines for vep-rs and fastVEP and over Ensembl VEP's one reference run per chromosome. The concordance figures cover those four plus the two genome-wide ClinVar releases, six whole-genome datasets in all (1,188,425,523 variants, 11,838,755,797 annotation tuples, 11,838,744,378 matched), against Ensembl VEP 116.2. The [Current release](#current-release-v050) section has the tables and the method.

The adjusted F1 sets aside only the documented cases where Ensembl VEP's own output is wrong, or where the two engines represent a record differently. Each case is traced to the mechanism in Ensembl VEP's source and reproducible from the record ([docs/intended-divergences.md](docs/intended-divergences.md)). Outside those cases every tuple matches, so vep-rs's remaining disagreements with Ensembl VEP are deliberate, documented corrections rather than defects.

The ten built-in plugin databases are CADD, REVEL, gnomADc, AlphaMissense, dbscSNV, LoFTEE, LoFtool and pLI on both GRCh37 and GRCh38, plus GWAS and SpliceAI on GRCh38. Dynamic-library plugins load from `--dir_plugins` ([docs/plugins.md](docs/plugins.md)). Parquet output needs the `duckdb` CLI on PATH.

### What "drop-in" means

vep-rs takes the same input files and flag names as Ensembl VEP and writes the same output layouts. Two things differ:

- **The transcript cache.** The cache is a directory built by `vep-cache-builder` or converted from an Ensembl VEP cache ([docs/cache-setup.md](docs/cache-setup.md)), found under `--dir_cache`, `--species` and `--cache_version` in Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout or named directly with `--json_cache`; a run with no cache stops with an error, and vep-rs never contacts a database.
- **Some flags are accepted but not implemented.** Among them are `--regulatory`, `--custom` and the statistics flags (vep-rs writes no summary file); `--refseq` and `--merged` select the `_refseq` or `_merged` species directory, and `--refseq` also resolves an HGVS gene name among RefSeq transcripts first. [docs/cli-reference.md](docs/cli-reference.md) marks every such flag. [docs/getting-started/what-vep-rs-does-not-do.md](docs/getting-started/what-vep-rs-does-not-do.md) says in plain words what a run does when asked for each missing feature, and which of them the next release, 0.6.0, adds.

## Current release: v0.5.0

These are the released version's own figures, measured against Ensembl VEP 116.2 under the release protocol stated under Method. The figures the paper reports, measured at its pinned build, are in [docs/published-figures.md](docs/published-figures.md) and stay as published. How each number was produced is under [Method](#method).

### Whole-genome concordance

Against Ensembl VEP 116.2. A population dataset is scored per chromosome and its shards pooled; each ClinVar release is one genome-wide file. The per-chromosome tables are in the [provenance records](#provenance-records).

| Engine | Dataset | Assembly | Shards | Raw F1 | Adjusted F1 | Variants | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.5.0 | gnomAD v4.1 | GRCh38 | 24 | 0.999999353 | 1.000000000 | 759,336,320 | 9,601,562,823 | 9,601,562,823 | 9,601,556,614 |
| vep-rs v0.5.0 | gnomAD v2.1.1 | GRCh37 | 23 | 0.999999521 | 1.000000000 | 261,942,336 | 885,664,902 | 885,664,902 | 885,664,478 |
| vep-rs v0.5.0 | 1000 Genomes Phase 3 | GRCh37 | 25 | 0.999999714 | 1.000000000 | 84,805,772 | 286,649,478 | 286,649,478 | 286,649,396 |
| vep-rs v0.5.0 | 1000 Genomes high-coverage | GRCh38 | 23 | 0.999991204 | 1.000000000 | 73,554,796 | 905,822,094 | 905,836,944 | 905,821,551 |
| vep-rs v0.5.0 | ClinVar | GRCh37 | 1 | 0.999979065 | 1.000000000 | 4,388,172 | 33,483,881 | 33,483,881 | 33,483,180 |
| vep-rs v0.5.0 | ClinVar | GRCh38 | 1 | 0.999972446 | 1.000000000 | 4,398,127 | 125,572,619 | 125,572,619 | 125,569,159 |
| vep-rs v0.5.0 | All six | GRCh37 and GRCh38 | 97 | 0.999998408 | 1.000000000 | 1,188,425,523 | 11,838,755,797 | 11,838,770,647 | 11,838,744,378 |
| vep-rs v0.4.0 (re-scored) | ClinVar | GRCh37 | 1 | 0.999979065 | 1.000000000 | 4,388,172 | 33,483,881 | 33,483,881 | 33,483,180 |
| vep-rs v0.4.0 (re-scored) | ClinVar | GRCh38 | 1 | 0.999972446 | 1.000000000 | 4,398,127 | 125,572,619 | 125,572,619 | 125,569,159 |
| fastVEP v0.4.0 | gnomAD v4.1 | GRCh38 | 24 | 0.923598441 | 0.923598441 | 759,336,320 | 9,601,562,823 | 9,646,230,422 | 8,888,615,919 |
| fastVEP v0.4.0 | gnomAD v2.1.1 | GRCh37 | 23 | 0.939631987 | 0.939631987 | 261,942,336 | 885,664,902 | 889,340,923 | 833,926,125 |
| fastVEP v0.4.0 | 1000 Genomes Phase 3 | GRCh37 | 25 | 0.977591321 | 0.978015902 | 84,805,772 | 286,649,478 | 287,933,154 | 280,853,497 |
| fastVEP v0.4.0 | 1000 Genomes high-coverage | GRCh38 | 23 | 0.934882830 | 0.935614047 | 73,554,796 | 905,822,094 | 910,030,911 | 848,804,898 |
| fastVEP v0.4.0 | ClinVar | GRCh37 | 1 | 0.971890563 | 0.971919461 | 4,388,172 | 33,483,881 | 33,748,452 | 32,671,235 |
| fastVEP v0.4.0 | ClinVar | GRCh38 | 1 | 0.974289466 | 0.974360417 | 4,398,127 | 125,572,619 | 126,203,713 | 122,651,514 |
| fastVEP v0.4.0 | All six | GRCh37 and GRCh38 | 97 | 0.927642871 | 0.927708892 | 1,188,425,523 | 11,838,755,797 | 11,893,487,575 | 11,007,523,188 |

### Whole-genome wall time

Seconds: a population dataset as the sum of its per-chromosome medians, each ClinVar release as its whole-file median. The figure in parentheses beside a comparator's time is that time divided by vep-rs's, taken only where both were measured under this protocol. The two right-hand columns are vep-rs on its largest chromosome, or on the whole file for ClinVar. With 5 measurements per cell the P5 and P95 are the cell's minimum and maximum.

- vep-rs and fastVEP ran the population datasets on `c9gd.8xlarge` (ARM Graviton5) and `c8id.8xlarge` (x86 Intel), the newest available class per architecture (Graviton5 against 8th-generation Intel), 32 vCPU with local NVMe.
- Ensembl VEP 116.2's figures are its single reference run per chromosome of the same inputs under the same protocol, on `c9gd.8xlarge` and `c8id.8xlarge` (the same processors with 64 GiB), the newest available class per architecture (Graviton5 against 8th-generation Intel), except gnomAD v4.1 chromosomes 1, 2 and 3 on `c9gd.12xlarge` (ARM, 48 vCPU / 96 GiB; peak memory 77.9, 76.4 and 80.2 GiB) and `c8id.12xlarge` (x86, 48 vCPU / 96 GiB; peak memory 73.2, 75.7 and 77.3 GiB), over the standard class's 64 GiB. The concordance tables compare against this run's output.
- The ClinVar rows are the released version's whole-file medians on `c9gd.8xlarge` and `c8id.8xlarge`, the newest available class per architecture (Graviton5 against 8th-generation Intel); no comparator time is stated on them.

| Dataset | Assembly | Variants | Arch | vep-rs v0.5.0 (s) | fastVEP v0.4.0 (s) | Ensembl VEP 116.2 (s) | vep-rs, largest cell median (s) | P5 to P95 (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gnomAD v4.1 | GRCh38 | 759,336,320 | ARM | 624 | 12,018 (19.3x) | 255,394 (410x) | 53.44 (chromosome 2) | 53.35 to 53.95 |
| gnomAD v4.1 | GRCh38 | 759,336,320 | x86 | 1,012 | 11,924 (11.8x) | 265,212 (262x) | 85.12 (chromosome 2) | 84.27 to 85.54 |
| gnomAD v2.1.1 | GRCh37 | 261,942,336 | ARM | 214 | 4,783 (22.4x) | 88,002 (411x) | 18.06 (chromosome 2) | 17.47 to 19.70 |
| gnomAD v2.1.1 | GRCh37 | 261,942,336 | x86 | 405 | 5,046 (12.5x) | 103,809 (256x) | 33.70 (chromosome 2) | 33.21 to 33.98 |
| 1000 Genomes Phase 3 | GRCh37 | 84,805,772 | ARM | 27 | 398 (14.9x) | 7,209 (270x) | 2.15 (chromosome 2) | 2.14 to 2.17 |
| 1000 Genomes Phase 3 | GRCh37 | 84,805,772 | x86 | 43 | 397 (9.2x) | 7,039 (162x) | 3.54 (chromosome 2) | 3.53 to 3.55 |
| 1000 Genomes high-coverage | GRCh38 | 73,554,796 | ARM | 60 | 951 (15.8x) | 15,817 (264x) | 5.02 (chromosome 1) | 4.99 to 5.02 |
| 1000 Genomes high-coverage | GRCh38 | 73,554,796 | x86 | 92 | 917 (10.0x) | 17,660 (193x) | 7.62 (chromosome 2) | 7.42 to 7.70 |
| ClinVar | GRCh37 | 4,388,172 | ARM | 3.43 |  |  | 3.43 (whole file) | 3.42 to 4.04 |
| ClinVar | GRCh37 | 4,388,172 | x86 | 5.15 |  |  | 5.15 (whole file) | 5.08 to 5.22 |
| ClinVar | GRCh38 | 4,398,127 | ARM | 11.34 |  |  | 11.34 (whole file) | 11.30 to 11.44 |
| ClinVar | GRCh38 | 4,398,127 | x86 | 15.43 |  |  | 15.43 (whole file) | 15.30 to 15.57 |
| All four population datasets | GRCh37 and GRCh38 | 1,179,639,224 | ARM | 924 | 18,149 (19.6x) | 366,422 (396x) | | |
| All four population datasets | GRCh37 and GRCh38 | 1,179,639,224 | x86 | 1,551 | 18,284 (11.8x) | 393,720 (254x) | | |

### Structural variants

The structural-variant per-VCF sets (16 files per assembly, the files' counts pooled), against Ensembl VEP 116.2. vep-rs v0.4.0's rows are its re-scoring against the same reference.

| Engine | Assembly | Raw F1 | Adjusted F1 | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.5.0 | GRCh37 | 0.964042027 | 1.000000000 | 317,047 | 337,773 | 315,637 |
| vep-rs v0.5.0 | GRCh38 | 0.905760853 | 1.000000000 | 1,455,318 | 1,734,800 | 1,444,742 |
| vep-rs v0.4.0 (re-scored) | GRCh37 | 0.959863779 | 1.000000000 | 317,047 | 337,773 | 314,269 |
| vep-rs v0.4.0 (re-scored) | GRCh38 | 0.903409216 | 1.000000000 | 1,455,318 | 1,734,800 | 1,440,991 |

### Structural-variant wall time

The released version on `c9gd.8xlarge` and `c8id.8xlarge`, the newest available class per architecture (Graviton5 against 8th-generation Intel), 16 files per assembly run as 16 processes: a start-up-bound regime, excluded from every speed-up above. With 5 measurements per cell the P5 and P95 are the cell's minimum and maximum.

| Arch | Dataset | Assembly | v0.5.0 median (s) | P5 to P95 (s) |
| --- | --- | --- | --- | --- |
| ARM | SV per-VCF (16 files) | GRCh37 | 1.24 | 1.23 to 1.25 |
| ARM | SV per-VCF (16 files) | GRCh38 | 3.39 | 3.32 to 3.46 |
| x86 | SV per-VCF (16 files) | GRCh37 | 1.85 | 1.84 to 1.87 |
| x86 | SV per-VCF (16 files) | GRCh38 | 4.82 | 4.77 to 4.87 |

### Method

- **Terms.** A variant is one VCF record. A tuple is one output row reduced to its Location, Allele, Feature, Feature_type and Consequence set. F1 pools a dataset's chromosomes by summed counts. HGVS notation is outside every F1 above.
- **Concordance.** One run per dataset of a build of the released version, against Ensembl VEP 116.2's output on the same inputs, with the same comparator, cache and reference FASTA. F1 is deterministic per binary and independent of architecture.
- **Adjusted F1.** The tuples it sets aside are the documented classes of [docs/intended-divergences.md](docs/intended-divergences.md). Every tuple outside those classes is matched.
- **fastVEP v0.4.0** was scored on the same Ensembl VEP output, with the same comparator, cache and reference and the same classes set aside. Its rows state its own tuple counts. Its adjusted figure differs from its raw one where the documented classes reach its whole-genome rows (1000 Genomes Phase 3 on GRCh37; 1000 Genomes high-coverage on GRCh38; ClinVar on GRCh37; ClinVar on GRCh38; All six on both assemblies). It was scored on the ClinVar cells and not on the structural-variant cells, whose tables carry no fastVEP v0.4.0 row.
- **Previous release.** vep-rs v0.4.0's rows are its re-scoring against the same Ensembl VEP 116.2 reference on the cells this release's record pairs (`previous_release_concordance`), measured under the same protocol; tuple counts move against v0.4.0 on SV per-VCF (16 files) GRCh37, SV per-VCF (16 files) GRCh38; its published figures were scored against the same reference and equal the re-scoring.
- **Wall time.** For vep-rs and fastVEP, the median of 5 independent machine measurements per cell on both architectures, with `--fork 16` on sites-only inputs; for Ensembl VEP, its one reference run per cell. Each measurement runs on a fresh instance after a discarded warmup, whose output is flushed and deleted before the timed run.

### Provenance records

Every value above is recorded in these files:

- [`docs/concordance-provenance/2026-10-25-release-v0.5.0.json`](docs/concordance-provenance/2026-10-25-release-v0.5.0.json): the ClinVar and structural-variant cells, the per-column field agreement, the compiler and build inputs of the measured binaries, the reference sets, the same-day wall-time control and the re-scoring of v0.4.0, and the whole-genome concordance per chromosome.
- [`docs/concordance-provenance/2026-10-25-population-v0.5.0.json`](docs/concordance-provenance/2026-10-25-population-v0.5.0.json): every chromosome's wall-time median, P5 and P95 for all three engines, with the 5 wall times behind vep-rs's and fastVEP's, and fastVEP's concordance per chromosome.
- [`docs/concordance-provenance/2026-10-24-comparator-fastvep-v0.4.0.json`](docs/concordance-provenance/2026-10-24-comparator-fastvep-v0.4.0.json): fastVEP v0.4.0's ClinVar figures, scored on the same Ensembl VEP output with the same comparator, cache and reference, the same classes set aside.
- [`docs/concordance-provenance/2026-10-16-release-v0.4.0.json`](docs/concordance-provenance/2026-10-16-release-v0.4.0.json): vep-rs v0.4.0's own record, with the figures it published.
- [`docs/concordance-provenance/2026-10-16-population-v0.4.0.json`](docs/concordance-provenance/2026-10-16-population-v0.4.0.json): vep-rs v0.4.0's whole-genome measurement.

### Reproducing the figures

- `scripts/concordance/run_clone_measurement.sh` re-measures the ClinVar and structural-variant cells on your own machines.
- `scripts/concordance/run_concordance.sh` times both engines on a directory of your own VCFs.
- [scripts/README.md](scripts/README.md#reproducing-the-published-concordance) is the runbook for reproducing the concordance numbers.

## Installation

Pre-built binaries for Linux x86_64, Linux aarch64 and macOS aarch64 are on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page; each Linux architecture has a measured archive and a generic one, and both run on glibc 2.17 or later. Building from source needs Rust 1.89+ ([rustup](https://rustup.rs/)):

```bash
git clone https://github.com/natera-open-source/vep-rs.git
cd vep-rs
cargo build --release
# Binary at target/release/vep
```

Each release from 0.2.0 on is also a container image, `ghcr.io/natera-open-source/vep-rs`, multi-platform from 0.5.0: the `X.Y.Z`, `X.Y` and `latest` tags carry the generic Linux binaries for linux/amd64 and linux/arm64 and run on Apple silicon and Graviton2 without a platform flag; `X.Y.Z-x86-64-v3` and `X.Y.Z-neoverse-v1` carry the measured binaries. The [Installation](docs/getting-started/installation.md) page covers each archive's CPU floor, how much slower the generic archives are (each release page reports it, or that the release did not time them), how to build for an older CPU, and the image's tags, platforms and attestation.

## Quick start

### 1. Try the bundled corpus

This needs no download. It runs the bundled 1,008-record excerpt of the benchmark datasets against the pruned GRCh37 cache that ships in `tests/golden/`:

```bash
target/release/vep \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --assembly GRCh37 \
  -i tests/golden/116/GRCh37/variants.vcf -o output.txt
```

### 2. Annotate your own data

Real data needs a JSON transcript cache, built natively with `vep-cache-builder` or converted from an Ensembl VEP cache. Then run the same VEP command you run, with `--json_cache` naming the cache directory. The [Quick start](docs/getting-started/quick-start.md) page gives the cache-builder command, the annotation command and the plugin and HGVS options. [Cache setup](docs/cache-setup.md) explains each input of the cache build.

## Documentation

The documentation site, https://natera-open-source.github.io/vep-rs/, describes the latest release; its pages are the `docs/` tree of this repository, in six sections. [Getting started](docs/getting-started/installation.md) covers installation, the transcript cache, the move from Ensembl VEP and, on the [Examples](docs/examples.md) page, the common shapes of a run. The reference pages cover [every flag](docs/cli-reference.md), with the option rules, the flags refused or accepted without effect and the exit-status contract, [every output format](docs/output-formats.md) and [every plugin](docs/plugins.md). [Concordance with Ensembl VEP](docs/intended-divergences.md) documents each divergence class the adjusted F1 sets aside, the [Current release](#current-release-v050) section above states the release's figures beside [the paper's](docs/published-figures.md), [Writing a plugin](docs/developers/writing-a-plugin.md) is the developer guide, and [Project](docs/project/citing.md) holds the citation and the changelog.

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

The journal DOI will be added here on acceptance. The paper's figures are in [docs/published-figures.md](docs/published-figures.md) and re-derive from `manuscript/data/*.csv` in this tree ([manuscript/data/README.md](manuscript/data/README.md)). The current release's own measurements are in [Current release](#current-release-v050).

The software itself is archived at Zenodo, and that archive is the citation for the code and the released data files:

> Porter M, Borkowski R. _vep-rs_, version 0.1.0. Zenodo (2026). [doi:10.5281/zenodo.22837897](https://doi.org/10.5281/zenodo.22837897)

`CITATION.cff` at the repository root carries the software citation in machine-readable form, with the all-versions DOI ([10.5281/zenodo.22837896](https://doi.org/10.5281/zenodo.22837896)) and the archived version's own. GitHub's "Cite this repository" button and citation managers pick it up directly. The article is added there on acceptance.

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
