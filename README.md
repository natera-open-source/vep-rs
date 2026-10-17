# vep-rs

[![CI](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![MSRV](https://img.shields.io/badge/MSRV-1.88-orange.svg)](#installation)

Documentation: https://natera-open-source.github.io/vep-rs/

A fast, memory-efficient reimplementation of the Ensembl Variant Effect Predictor (VEP) in Rust. It takes VEP's inputs and flags, writes VEP's consequence vocabulary in VEP's output formats, and annotates population-scale callsets two orders of magnitude faster than the Perl implementation.

In this README, "Perl VEP" is the reference implementation, `ensemblorg/ensembl-vep` release 116.2. The accompanying paper calls the same engine "VEP".

## Why vep-rs instead of Ensembl VEP or fastVEP

|                                  | vep-rs v0.4.0                                                                       |
| -------------------------------- | ----------------------------------------------------------------------------------- |
| **Speed vs Ensembl VEP 115.2**   | 343x faster on ARM (Graviton4), 187x on x86; 96.2x to 421x per dataset              |
| **Speed vs fastVEP v0.4.0**      | 16.2x faster on ARM, 8.2x on x86; 5.2x to 21.6x per dataset                         |
| **All of gnomAD v4.1**           | 759,336,320 variants, 24 chromosomes: 881 s on one ARM machine, 1,365 s on x86      |
| **Peak memory**                  | At most 2.6 GB RSS on any chromosome                                                |
| **Concordance with Ensembl VEP** | Adjusted F1 1.000000000 on all six whole-genome datasets, raw F1 0.999998408 pooled |
| **Plugins**                      | 10 built-in databases, plus C-ABI dynamic-library plugins                           |
| **Output formats**               | VEP default, VCF, JSON, Tab, Parquet                                                |
| **Threading**                    | `--fork N`; the default uses every logical CPU up to 32                             |

The speed figures cover all of four population datasets (1,179,639,224 variants) under one protocol: the sum of per-chromosome medians over 20 machines at 16 threads. The concordance figures cover those four plus the two genome-wide ClinVar releases, six whole-genome datasets in all (1,188,425,523 variants, 11,838,755,797 annotation tuples, 11,838,744,378 matched), against Ensembl VEP 116.2. The [Current release](#current-release-v040) section has the tables and the method.

The adjusted F1 sets aside only the documented cases where Ensembl VEP's own output is wrong, or where the two engines represent a record differently. Each case is traced to the mechanism in Ensembl VEP's source and reproducible from the record ([docs/intended-divergences.md](docs/intended-divergences.md)). Outside those cases every tuple matches, so vep-rs's remaining disagreements with Ensembl VEP are deliberate, documented corrections rather than defects.

The ten built-in plugin databases are CADD, REVEL, gnomADc, AlphaMissense, dbscSNV, LoFTEE, LoFtool and pLI on both GRCh37 and GRCh38, plus GWAS and SpliceAI on GRCh38. Dynamic-library plugins load from `--dir_plugins` ([docs/plugins.md](docs/plugins.md)). Parquet output needs the `duckdb` CLI on PATH.

### What "drop-in" means

vep-rs takes the same input files and flag names as Perl VEP and writes the same output layouts. Two things differ:

- **The transcript cache.** `--json_cache` names a directory built by `vep-cache-builder` or converted from a Perl VEP cache ([docs/cache-setup.md](docs/cache-setup.md)). `--cache`, `--offline`, `--dir_cache`, `--cache_version` and `--species` are accepted but don't locate a cache, and vep-rs never contacts a database.
- **Some flags are accepted but not implemented.** Among them are `--regulatory`, `--custom`, `--refseq`, `--merged`, `--coding_only`, `--pick_order` and the statistics flags (vep-rs writes no summary file). [docs/cli-reference.md](docs/cli-reference.md) marks every such flag.

## Current release: v0.4.0

These are the released version's own figures, measured under the paper's protocol against Ensembl VEP 116.2. The figures the paper reports, measured at its pinned build, are in [docs/published-figures.md](docs/published-figures.md) and stay as published. How each number was produced is under [Method](#method).

### Whole-genome concordance

Against Ensembl VEP 116.2. A population dataset is scored per chromosome and its shards pooled; each ClinVar release is one genome-wide file. The per-chromosome tables are in the [provenance records](#provenance-records).

| Engine | Dataset | Assembly | Shards | Raw F1 | Adjusted F1 | Variants | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.4.0 | gnomAD v4.1 | GRCh38 | 24 | 0.999999353 | 1.000000000 | 759,336,320 | 9,601,562,823 | 9,601,562,823 | 9,601,556,614 |
| vep-rs v0.4.0 | gnomAD v2.1.1 | GRCh37 | 23 | 0.999999521 | 1.000000000 | 261,942,336 | 885,664,902 | 885,664,902 | 885,664,478 |
| vep-rs v0.4.0 | 1000 Genomes Phase 3 | GRCh37 | 25 | 0.999999714 | 1.000000000 | 84,805,772 | 286,649,478 | 286,649,478 | 286,649,396 |
| vep-rs v0.4.0 | 1000 Genomes high-coverage | GRCh38 | 23 | 0.999991204 | 1.000000000 | 73,554,796 | 905,822,094 | 905,836,944 | 905,821,551 |
| vep-rs v0.4.0 | ClinVar | GRCh37 | 1 | 0.999979065 | 1.000000000 | 4,388,172 | 33,483,881 | 33,483,881 | 33,483,180 |
| vep-rs v0.4.0 | ClinVar | GRCh38 | 1 | 0.999972446 | 1.000000000 | 4,398,127 | 125,572,619 | 125,572,619 | 125,569,159 |
| vep-rs v0.4.0 | All six | GRCh37 and GRCh38 | 97 | 0.999998408 | 1.000000000 | 1,188,425,523 | 11,838,755,797 | 11,838,770,647 | 11,838,744,378 |
| vep-rs v0.3.2 (re-scored) | ClinVar | GRCh37 | 1 | 0.999883317 | 0.999904250 | 4,388,172 | 33,483,881 | 33,483,881 | 33,479,974 |
| vep-rs v0.3.2 (re-scored) | ClinVar | GRCh38 | 1 | 0.999844146 | 0.999871696 | 4,398,127 | 125,572,619 | 125,572,619 | 125,553,048 |
| fastVEP v0.4.0 | gnomAD v4.1 | GRCh38 | 24 | 0.923598441 | 0.923598441 | 759,336,320 | 9,601,562,823 | 9,646,230,422 | 8,888,615,919 |
| fastVEP v0.4.0 | gnomAD v2.1.1 | GRCh37 | 23 | 0.939631987 | 0.939631987 | 261,942,336 | 885,664,902 | 889,340,923 | 833,926,125 |
| fastVEP v0.4.0 | 1000 Genomes Phase 3 | GRCh37 | 25 | 0.977591321 | 0.978015902 | 84,805,772 | 286,649,478 | 287,933,154 | 280,853,497 |
| fastVEP v0.4.0 | 1000 Genomes high-coverage | GRCh38 | 23 | 0.934882830 | 0.935614047 | 73,554,796 | 905,822,094 | 910,030,911 | 848,804,898 |
| fastVEP v0.4.0 | ClinVar | GRCh37 | 1 | 0.971890563 | 0.971919461 | 4,388,172 | 33,483,881 | 33,748,452 | 32,671,235 |
| fastVEP v0.4.0 | ClinVar | GRCh38 | 1 | 0.974289466 | 0.974360417 | 4,398,127 | 125,572,619 | 126,203,713 | 122,651,514 |
| fastVEP v0.4.0 | All six | GRCh37 and GRCh38 | 97 | 0.927642871 | 0.927708892 | 1,188,425,523 | 11,838,755,797 | 11,893,487,575 | 11,007,523,188 |

### Whole-genome wall time

Seconds: a population dataset as the sum of its per-chromosome medians, each ClinVar release as its whole-file median. The figure in parentheses beside a comparator's time is that time divided by vep-rs's, taken only where both were measured under this protocol. The two right-hand columns are vep-rs on its largest chromosome, or on the whole file for ClinVar.

- vep-rs and fastVEP ran the population datasets on `c8gd.8xlarge` (ARM Graviton4) and `c8id.8xlarge` (x86 Intel), 32 vCPU with local NVMe.
- Ensembl VEP 115.2's figures are its reference run of the same inputs under the same protocol, on `r8gd.8xlarge` and `r8id.8xlarge` (the same processors with 256 GiB), carried forward from the earlier reference generation. The concordance tables compare against Ensembl VEP 116.2's output.
- The ClinVar rows are the released version's whole-file medians on the paper's instance types (`r8gd.8xlarge`, `r8id.8xlarge`); the comparators' times on them are carried forward from an earlier measurement and labelled with the engine version they were measured at, and no ratio is taken against them.

| Dataset | Assembly | Variants | Arch | vep-rs v0.4.0 (s) | fastVEP v0.4.0 (s) | Ensembl VEP 115.2 (s) | vep-rs, largest cell median (s) | P5 to P95 (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gnomAD v4.1 | GRCh38 | 759,336,320 | ARM | 881 | 13,437 (15.2x) | 296,805 (337x) | 75.57 (chromosome 2) | 74.86 to 78.00 |
| gnomAD v4.1 | GRCh38 | 759,336,320 | x86 | 1,365 | 10,465 (7.7x) | 251,614 (184x) | 115.44 (chromosome 2) | 113.05 to 115.95 |
| gnomAD v2.1.1 | GRCh37 | 261,942,336 | ARM | 282 | 6,111 (21.6x) | 119,011 (421x) | 23.63 (chromosome 2) | 23.21 to 24.22 |
| gnomAD v2.1.1 | GRCh37 | 261,942,336 | x86 | 472 | 5,037 (10.7x) | 104,430 (221x) | 39.03 (chromosome 2) | 38.46 to 39.74 |
| 1000 Genomes Phase 3 | GRCh37 | 84,805,772 | ARM | 49 | 504 (10.4x) | 9,521 (196x) | 4.06 (chromosome 2) | 3.91 to 4.33 |
| 1000 Genomes Phase 3 | GRCh37 | 84,805,772 | x86 | 76 | 396 (5.2x) | 7,287 (96.2x) | 6.33 (chromosome 2) | 6.12 to 6.57 |
| 1000 Genomes high-coverage | GRCh38 | 73,554,796 | ARM | 85 | 991 (11.6x) | 19,591 (229x) | 7.08 (chromosome 2) | 6.76 to 7.76 |
| 1000 Genomes high-coverage | GRCh38 | 73,554,796 | x86 | 116 | 752 (6.5x) | 15,840 (137x) | 9.75 (chromosome 2) | 9.56 to 9.98 |
| ClinVar | GRCh37 | 4,388,172 | ARM | 4.54 | 43.65 (fastVEP 0.3.0) | 961.08 | 4.54 (whole file) | 4.36 to 4.86 |
| ClinVar | GRCh37 | 4,388,172 | x86 | 5.73 | 33.16 (fastVEP 0.3.0) | 849.97 | 5.73 (whole file) | 5.62 to 5.99 |
| ClinVar | GRCh38 | 4,398,127 | ARM | 13.96 | 121.13 (fastVEP 0.3.0) | 2,722.01 | 13.96 (whole file) | 13.53 to 14.51 |
| ClinVar | GRCh38 | 4,398,127 | x86 | 14.39 | 90.26 (fastVEP 0.3.0) | 2,378.35 | 14.39 (whole file) | 14.29 to 15.39 |
| All four population datasets | GRCh37 and GRCh38 | 1,179,639,224 | ARM | 1,298 | 21,043 (16.2x) | 444,929 (343x) | | |
| All four population datasets | GRCh37 and GRCh38 | 1,179,639,224 | x86 | 2,028 | 16,649 (8.2x) | 379,170 (187x) | | |

### Structural variants

The structural-variant per-VCF sets (16 files per assembly, the files' counts pooled), against Ensembl VEP 116.2. vep-rs v0.3.2's rows are its re-scoring against the same reference.

| Engine | Assembly | Raw F1 | Adjusted F1 | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.4.0 | GRCh37 | 0.959863779 | 1.000000000 | 317,047 | 337,773 | 314,269 |
| vep-rs v0.4.0 | GRCh38 | 0.903409216 | 1.000000000 | 1,455,318 | 1,734,800 | 1,440,991 |
| vep-rs v0.3.2 (re-scored) | GRCh37 | 0.959120631 | 0.999242374 | 317,047 | 337,513 | 313,901 |
| vep-rs v0.3.2 (re-scored) | GRCh38 | 0.901958840 | 0.998501112 | 1,455,318 | 1,731,619 | 1,437,243 |

### Structural-variant wall time

The released version on the paper's instance types (`r8gd.8xlarge`, `r8id.8xlarge`), 16 files per assembly run as 16 processes: a start-up-bound regime, excluded from every speed-up above.

| Arch | Dataset | Assembly | v0.4.0 median (s) | P5 to P95 (s) |
| --- | --- | --- | --- | --- |
| ARM | SV per-VCF (16 files) | GRCh37 | 0.55 | 0.53 to 0.60 |
| ARM | SV per-VCF (16 files) | GRCh38 | 1.71 | 1.65 to 1.87 |
| x86 | SV per-VCF (16 files) | GRCh37 | 0.74 | 0.71 to 0.83 |
| x86 | SV per-VCF (16 files) | GRCh38 | 1.98 | 1.95 to 2.18 |

### Method

- **Terms.** A variant is one VCF record. A tuple is one output row reduced to its Location, Allele, Feature, Feature_type and Consequence set. F1 pools a dataset's chromosomes by summed counts. HGVS notation is outside every F1 above.
- **Concordance.** One run per dataset of a build of the released version, against Ensembl VEP 116.2's output on the same inputs, with the same comparator, cache and reference FASTA. F1 is deterministic per binary and independent of architecture.
- **Adjusted F1.** The tuples it sets aside are the documented classes of [docs/intended-divergences.md](docs/intended-divergences.md). Every tuple outside those classes is matched.
- **fastVEP v0.4.0** was scored on the same Ensembl VEP output, with the same comparator, cache and reference and the same classes set aside. Its rows state its own tuple counts. Its adjusted figure differs from its raw one where the documented classes reach its whole-genome rows (1000 Genomes Phase 3 on GRCh37; 1000 Genomes high-coverage on GRCh38; ClinVar on GRCh37; ClinVar on GRCh38; All six on both assemblies). It was scored on the ClinVar cells and not on the structural-variant cells, whose tables carry no fastVEP v0.4.0 row.
- **Previous release.** vep-rs v0.3.2's rows are its re-scoring against the same Ensembl VEP 116.2 reference on the cells this release's record pairs (`previous_release_concordance`), measured the same day under the same protocol; tuple counts move against v0.3.2 on ClinVar GRCh37, ClinVar GRCh38, SV per-VCF (16 files) GRCh37, SV per-VCF (16 files) GRCh38; its published figures were scored against an earlier reference and are not compared across references.
- **Wall time.** The median of 20 independent machine measurements per cell on both architectures, with `--fork 16` on sites-only inputs. Each measurement runs on a fresh instance after a discarded warmup, whose output is flushed and deleted before the timed run.

### Provenance records

Every value above is recorded in these files:

- [`docs/concordance-provenance/2026-10-16-release-v0.4.0.json`](docs/concordance-provenance/2026-10-16-release-v0.4.0.json): the ClinVar and structural-variant cells, the per-column field agreement, the compiler and build inputs of the measured binaries, the reference sets, the same-day comparison against v0.3.2, and the whole-genome concordance per chromosome.
- [`docs/concordance-provenance/2026-10-16-population-v0.4.0.json`](docs/concordance-provenance/2026-10-16-population-v0.4.0.json): every chromosome's wall-time median, P5 and P95 with the 20 wall times behind them, for both engines, and fastVEP's concordance per chromosome.
- [`docs/concordance-provenance/2026-10-15-comparator-fastvep-v0.4.0.json`](docs/concordance-provenance/2026-10-15-comparator-fastvep-v0.4.0.json): fastVEP v0.4.0's ClinVar figures, scored on the same Ensembl VEP output with the same comparator, cache and reference, the same classes set aside.
- [`docs/concordance-provenance/2026-09-20-n20-parallel-clones.json`](docs/concordance-provenance/2026-09-20-n20-parallel-clones.json): the comparators' ClinVar wall times, carried forward and labelled with the engine version they were measured at.
- [`docs/concordance-provenance/2026-10-06-release-v0.3.2.json`](docs/concordance-provenance/2026-10-06-release-v0.3.2.json): vep-rs v0.3.2's own record, with the figures it published.
- [`docs/concordance-provenance/2026-10-06-population-v0.3.2.json`](docs/concordance-provenance/2026-10-06-population-v0.3.2.json): vep-rs v0.3.2's whole-genome measurement.

### Reproducing the figures

- `scripts/concordance/run_clone_measurement.sh` re-measures the ClinVar and structural-variant cells on your own machines.
- `scripts/concordance/run_concordance.sh` times both engines on a directory of your own VCFs.
- [scripts/README.md](scripts/README.md#reproducing-the-published-concordance) is the runbook for reproducing the concordance numbers.

## Installation

Pre-built binaries for Linux x86_64, Linux aarch64 and macOS aarch64 are on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page. The x86_64 binary needs AVX2. Building from source needs Rust 1.88+ ([rustup](https://rustup.rs/)):

```bash
git clone https://github.com/natera-open-source/vep-rs.git
cd vep-rs
cargo build --release
# Binary at target/release/vep
```

Each release is also a container image, `ghcr.io/natera-open-source/vep-rs`, holding the release's x86_64 Linux binaries and the `duckdb` CLI. The [Installation](docs/getting-started/installation.md) page covers the CPU floor every build targets, how to build for an older CPU, and the image's tags, platform and attestation.

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

Real data needs a JSON transcript cache, built natively with `vep-cache-builder` or converted from a Perl VEP cache. Then run the same VEP command you run today, with `--json_cache` naming the cache directory. The [Quick start](docs/getting-started/quick-start.md) page gives the cache-builder command, the annotation command and the plugin and HGVS options. [Cache setup](docs/cache-setup.md) explains each input of the cache build.

## Documentation

The documentation site, https://natera-open-source.github.io/vep-rs/, describes the latest release; its pages are the `docs/` tree of this repository, in six sections. [Getting started](docs/getting-started/installation.md) covers installation, the transcript cache and the move from Ensembl VEP, and the reference pages cover [every flag](docs/cli-reference.md), [every output format](docs/output-formats.md) and [every plugin](docs/plugins.md). [Concordance with Ensembl VEP](docs/intended-divergences.md) documents each divergence class the adjusted F1 sets aside, the [Current release](#current-release-v040) section above states the release's figures beside [the paper's](docs/published-figures.md), [Writing a plugin](docs/developers/writing-a-plugin.md) is the developer guide, and [Project](docs/project/citing.md) holds the citation and the changelog.

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
| [vep-io](crates/vep-io/README.md) | Input parsers (VCF, Ensembl default, HGVS, region) and the output writers (VEP, VCF, JSON, Tab, Parquet) | library |
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

The journal DOI will be added here on acceptance. The paper's figures are in [docs/published-figures.md](docs/published-figures.md) and re-derive from `manuscript/data/*.csv` in this tree ([manuscript/data/README.md](manuscript/data/README.md)). The current release's own measurements are in [Current release](#current-release-v040).

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
