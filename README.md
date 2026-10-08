# vep-rs

A fast, memory-efficient variant effect predictor written in Rust: a from-scratch reimplementation of the Ensembl Variant Effect Predictor (VEP) that takes the same inputs and flags, emits VEP's consequence vocabulary in VEP's output formats, and annotates population-scale callsets two orders of magnitude faster than the Perl implementation.

[![CI](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![MSRV](https://img.shields.io/badge/MSRV-1.88-orange.svg)](#installation)

vep-rs reads a JSON transcript cache, built from Ensembl's release files by `vep-cache-builder` or converted from an existing Perl VEP cache, and writes VEP's default, VCF, JSON and Tab formats plus Parquet. Throughout this README "Perl VEP" is the reference implementation, `ensemblorg/ensembl-vep` release 115.2; the accompanying paper calls the same engine "VEP".

## Why vep-rs instead of Ensembl VEP or fastVEP

|                    |                                                                                                                                                                                                                                 |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Performance**    | 397x faster than Ensembl VEP 115.2 on ARM Graviton4 and 209x on x86, and 18.8x and 9.2x faster than fastVEP v0.4.0, on the whole of four population datasets (1,179,639,224 variants) under one protocol, the sum of per-chromosome medians of 20 machines at 16 threads; per dataset the fold runs 98.0x to 422x over Ensembl VEP and 5.3x to 21.6x over fastVEP. The whole of gnomAD v4.1 (759,336,320 variants, 24 chromosomes) takes 724 s on one ARM machine and 1,169 s on x86; peak RSS at most 2.4 GB on any chromosome ([Current release](#current-release-v032)) |
| **Concordance**    | Adjusted F1 1.000000000 (raw F1 0.999998884) with Ensembl VEP 115.2 on every dataset: 1,179,639,224 variants, 9,281,654,995 annotation tuples, 9,281,650,663 matched ([Current release](#current-release-v032)). The tuples the adjusted figure sets aside are the documented cases in which Ensembl VEP's own output is wrong, or the two engines represent a record differently, each traced to the mechanism in Ensembl VEP's source and reproducible from the record ([docs/intended-divergences.md](docs/intended-divergences.md)); outside them every tuple is matched, so vep-rs's remaining disagreements with Ensembl VEP are deliberate, documented corrections rather than defects |
| **Plugins**        | 10 supported built-in databases: 8 on both GRCh37 and GRCh38 (CADD, REVEL, gnomADc, AlphaMissense, dbscSNV, LoFTEE, LoFtool, pLI) plus GWAS and SpliceAI on GRCh38, and C-ABI dynamic-library plugins loaded from `--dir_plugins`; see [docs/plugins.md](docs/plugins.md) |
| **Output formats** | VEP default, VCF, JSON, Tab, Parquet (Parquet needs the `duckdb` CLI on PATH)                                                                                                                                                    |
| **Threading**      | Parallel annotation via `--fork N`; the default uses every logical CPU up to 32                                                                                                                                                  |

**What drop-in means here.** The same input files, the same flag names, the same output layouts. Two things differ from Perl VEP. The transcript cache is `--json_cache`, a directory built by `vep-cache-builder` or converted from a Perl VEP cache ([docs/cache-setup.md](docs/cache-setup.md)); `--cache`, `--offline`, `--dir_cache`, `--cache_version` and `--species` are accepted, but they do not locate a cache, and vep-rs never contacts a database. And some Perl VEP flags are accepted for compatibility without being implemented, among them `--regulatory`, `--custom`, `--refseq`, `--merged`, `--coding_only`, `--pick_order` and the statistics flags (vep-rs writes no summary file); [docs/cli-reference.md](docs/cli-reference.md) marks every such flag.

## Current release: v0.3.2

The figures here are the released version's own, measured under the paper's protocol; the figures the paper reports, measured at the paper's pinned build, are in [docs/published-figures.md](docs/published-figures.md) and stay as published. Concordance is one run of a build of the released version per dataset (F1 is deterministic per binary and independent of architecture) against Ensembl VEP 115.2's output on the same inputs, with the same comparator, cache and reference FASTA; a variant is one VCF record, a tuple is one output row reduced to its Location, Allele, Feature, Feature_type and Consequence set, and F1 pools a dataset's chromosomes by summed counts. The tuples the adjusted figure sets aside are the documented classes of [docs/intended-divergences.md](docs/intended-divergences.md); every tuple outside those classes is matched. fastVEP v0.4.0 was scored on the same Ensembl VEP output with the same comparator, cache and reference, the same classes set aside; its rows state its own tuple counts, and its adjusted figure equals its raw one because the classes set aside are shapes it does not produce. Wall time is the median of 20 independent machine measurements per cell on both architectures, each on a fresh instance after a discarded warmup whose output is flushed and deleted before the timed run, `--fork 16`, on sites-only inputs. Every value below is recorded in [`docs/concordance-provenance/2026-10-06-release-v0.3.2.json`](docs/concordance-provenance/2026-10-06-release-v0.3.2.json) (the chromosome 21 cells, the per-column field agreement, the compiler and build inputs of the measured binaries, the same-day comparison against v0.3.1, the whole-genome concordance per chromosome) and [`docs/concordance-provenance/2026-10-06-population-v0.3.2.json`](docs/concordance-provenance/2026-10-06-population-v0.3.2.json) (every chromosome's wall-time median, P5 and P95 with the 20 wall times behind them, for both engines, and fastVEP's concordance per chromosome). vep-rs v0.3.1's whole-genome measurement is [`docs/concordance-provenance/2026-10-04-population-v0.3.1.json`](docs/concordance-provenance/2026-10-04-population-v0.3.1.json). fastVEP v0.4.0's ClinVar and structural-variant figures are [`docs/concordance-provenance/2026-10-07-comparator-fastvep-v0.4.0.json`](docs/concordance-provenance/2026-10-07-comparator-fastvep-v0.4.0.json), scored on the same Ensembl VEP output with the same comparator, cache and reference, the same classes set aside.

Whole-genome concordance against Ensembl VEP 115.2 (the per-chromosome tables are in the records):

| Engine | Dataset | Assembly | Chromosomes | Raw F1 | Adjusted F1 | Variants | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.3.2 | gnomAD v4.1 | GRCh38 | 24 | 0.999999522 | 1.000000000 | 759,336,320 | 7,409,818,550 | 7,409,818,550 | 7,409,815,011 |
| vep-rs v0.3.2 | gnomAD v2.1.1 | GRCh37 | 23 | 0.999999521 | 1.000000000 | 261,942,336 | 885,664,902 | 885,664,902 | 885,664,478 |
| vep-rs v0.3.2 | 1000 Genomes Phase 3 | GRCh37 | 25 | 0.999999714 | 1.000000000 | 84,805,772 | 286,649,478 | 286,649,478 | 286,649,396 |
| vep-rs v0.3.2 | 1000 Genomes high-coverage | GRCh38 | 23 | 0.999990975 | 1.000000000 | 73,554,796 | 699,522,065 | 699,534,118 | 699,521,778 |
| vep-rs v0.3.2 | All four | GRCh37 and GRCh38 | 95 | 0.999998884 | 1.000000000 | 1,179,639,224 | 9,281,654,995 | 9,281,667,048 | 9,281,650,663 |
| fastVEP v0.4.0 | gnomAD v4.1 | GRCh38 | 24 | 0.923864433 | 0.923864433 | 759,336,320 | 7,409,818,550 | 7,454,247,630 | 6,866,191,040 |
| fastVEP v0.4.0 | gnomAD v2.1.1 | GRCh37 | 23 | 0.939632510 | 0.939632510 | 261,942,336 | 885,664,902 | 889,340,923 | 833,926,589 |
| fastVEP v0.4.0 | 1000 Genomes Phase 3 | GRCh37 | 25 | 0.977591449 | 0.977591449 | 84,805,772 | 286,649,478 | 287,933,154 | 280,853,534 |
| fastVEP v0.4.0 | 1000 Genomes high-coverage | GRCh38 | 23 | 0.935064819 | 0.935064819 | 73,554,796 | 699,522,065 | 703,703,820 | 656,053,579 |
| fastVEP v0.4.0 | All four | GRCh37 and GRCh38 | 95 | 0.927870245 | 0.927870245 | 1,179,639,224 | 9,281,654,995 | 9,335,225,527 | 8,637,024,742 |

Whole-genome wall time, in seconds, as the sum of the per-chromosome medians; the parenthesised figure beside a comparator's time is that time over vep-rs's, and the largest chromosome's own median and P5 to P95 stand at the right. vep-rs and fastVEP ran on `c8gd.8xlarge` (ARM Graviton4) and `c8id.8xlarge` (x86 Intel), 32 vCPU with local NVMe; Ensembl VEP 115.2's figures are its reference run of the same inputs under the same protocol on `r8gd.8xlarge` and `r8id.8xlarge` (the same processors with 256 GiB), the run whose output the concordance tables compare against:

| Dataset | Assembly | Variants | Arch | vep-rs v0.3.2 (s) | fastVEP v0.4.0 (s) | Ensembl VEP 115.2 (s) | vep-rs, chromosome 2 median (s) | P5 to P95 (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gnomAD v4.1 | GRCh38 | 759,336,320 | ARM | 724 | 13,437 (18.6x) | 296,805 (410x) | 62.48 | 62.00 to 63.73 |
| gnomAD v4.1 | GRCh38 | 759,336,320 | x86 | 1,169 | 10,465 (9.0x) | 251,614 (215x) | 99.34 | 97.59 to 100.30 |
| gnomAD v2.1.1 | GRCh37 | 261,942,336 | ARM | 282 | 6,111 (21.6x) | 119,011 (422x) | 23.61 | 23.12 to 24.30 |
| gnomAD v2.1.1 | GRCh37 | 261,942,336 | x86 | 476 | 5,037 (10.6x) | 104,430 (219x) | 39.65 | 38.83 to 40.19 |
| 1000 Genomes Phase 3 | GRCh37 | 84,805,772 | ARM | 46 | 504 (11.0x) | 9,521 (208x) | 3.82 | 3.63 to 4.02 |
| 1000 Genomes Phase 3 | GRCh37 | 84,805,772 | x86 | 74 | 396 (5.3x) | 7,287 (98.0x) | 6.17 | 6.03 to 6.49 |
| 1000 Genomes high-coverage | GRCh38 | 73,554,796 | ARM | 69 | 991 (14.3x) | 19,591 (283x) | 5.84 | 5.70 to 6.13 |
| 1000 Genomes high-coverage | GRCh38 | 73,554,796 | x86 | 99 | 752 (7.6x) | 15,840 (160x) | 8.38 | 8.09 to 8.63 |
| All four | GRCh37 and GRCh38 | 1,179,639,224 | ARM | 1,121 | 21,043 (18.8x) | 444,929 (397x) | | |
| All four | GRCh37 and GRCh38 | 1,179,639,224 | x86 | 1,818 | 16,649 (9.2x) | 379,170 (209x) | | |

ClinVar, full set, and the structural-variant per-VCF sets (16 files per assembly, the files' counts pooled), both engines against Ensembl VEP 115.2; fastVEP v0.4.0's rows are from [`docs/concordance-provenance/2026-10-07-comparator-fastvep-v0.4.0.json`](docs/concordance-provenance/2026-10-07-comparator-fastvep-v0.4.0.json), and its adjusted figure differs from its raw one where the documented classes reach its rows (the structural-variant sets on GRCh37 and GRCh38):

| Engine | Dataset | Assembly | Raw F1 | Adjusted F1 | Ensembl VEP tuples | Engine tuples | Matched |
| --- | --- | --- | --- | --- | --- | --- | --- |
| vep-rs v0.3.2 | ClinVar full | GRCh37 | 0.999979065 | 1.000000000 | 33,483,881 | 33,483,881 | 33,483,180 |
| vep-rs v0.3.2 | ClinVar full | GRCh38 | 0.999974091 | 1.000000000 | 91,552,813 | 91,552,813 | 91,550,441 |
| vep-rs v0.3.2 | SV per-VCF (16 files) | GRCh37 | 0.969197320 | 1.000000000 | 329,864 | 337,513 | 323,410 |
| vep-rs v0.3.2 | SV per-VCF (16 files) | GRCh38 | 0.902664038 | 1.000000000 | 1,117,037 | 1,274,145 | 1,079,217 |
| fastVEP v0.4.0 | ClinVar full | GRCh37 | 0.971925190 | 0.971925190 | 33,483,881 | 33,748,452 | 32,672,399 |
| fastVEP v0.4.0 | ClinVar full | GRCh38 | 0.973008467 | 0.973008467 | 91,552,813 | 92,175,395 | 89,384,551 |
| fastVEP v0.4.0 | SV per-VCF (16 files) | GRCh37 | 0.163918745 | 0.164480344 | 329,864 | 301,583 | 51,753 |
| fastVEP v0.4.0 | SV per-VCF (16 files) | GRCh38 | 0.130548022 | 0.131981905 | 1,117,037 | 982,788 | 137,064 |

Chromosome 21 suites, the paper's cells measured on the released version (the v0.3.1 columns are the concordance of [`docs/concordance-provenance/2026-10-02-release-v0.3.1.json`](docs/concordance-provenance/2026-10-02-release-v0.3.1.json); raw F1 equals v0.3.1's on every suite because this version changes no annotation: its output is byte-identical to the previous release's on every gated input):

| Dataset | Assembly | v0.3.2 Raw F1 | v0.3.2 Adj F1 | v0.3.1 Raw F1 | v0.3.1 Adj F1 |
| --- | --- | --- | --- | --- | --- |
| ClinVar full | GRCh37 | 0.999979065 | 1.000000000 | 0.999979065 | 1.000000000 |
| ClinVar full | GRCh38 | 0.999974091 | 1.000000000 | 0.999974091 | 1.000000000 |
| gnomAD v2.1.1 chr21 | GRCh37 | 0.999999729 | 1.000000000 | 0.999999729 | 1.000000000 |
| gnomAD v4.1 chr21 | GRCh38 | 0.999999303 | 1.000000000 | 0.999999303 | 1.000000000 |
| 1KG Phase 3 chr21 | GRCh37 | 1.000000000 | 1.000000000 | 1.000000000 | 1.000000000 |
| 1KG high-cov chr21 | GRCh38 | 0.999999904 | 1.000000000 | 0.999999904 | 1.000000000 |
| SV per-VCF (16 files) | GRCh37 | 0.969197320 | 1.000000000 | 0.969197320 | 1.000000000 |
| SV per-VCF (16 files) | GRCh38 | 0.902664038 | 1.000000000 | 0.902664038 | 1.000000000 |

Chromosome 21 wall time of the released version on the paper's instance types (`r8gd.8xlarge`, `r8id.8xlarge`):

| Arch | Dataset | Assembly | v0.3.2 median (s) | P5 to P95 (s) |
| --- | --- | --- | --- | --- |
| ARM | ClinVar full | GRCh37 | 4.31 | 4.18 to 4.40 |
| ARM | ClinVar full | GRCh38 | 10.43 | 10.24 to 10.61 |
| ARM | gnomAD v2.1.1 chr21 | GRCh37 | 3.32 | 3.27 to 3.40 |
| ARM | gnomAD v4.1 chr21 | GRCh38 | 9.66 | 9.37 to 10.07 |
| ARM | 1KG Phase 3 chr21 | GRCh37 | 0.56 | 0.54 to 0.60 |
| ARM | 1KG high-cov chr21 | GRCh38 | 0.95 | 0.92 to 0.98 |
| ARM | SV per-VCF (16 files) | GRCh37 | 0.53 | 0.52 to 0.53 |
| ARM | SV per-VCF (16 files) | GRCh38 | 1.27 | 1.25 to 1.34 |
| x86 | ClinVar full | GRCh37 | 5.64 | 5.50 to 5.74 |
| x86 | ClinVar full | GRCh38 | 11.46 | 11.26 to 11.88 |
| x86 | gnomAD v2.1.1 chr21 | GRCh37 | 5.86 | 5.55 to 6.15 |
| x86 | gnomAD v4.1 chr21 | GRCh38 | 15.63 | 15.51 to 15.90 |
| x86 | 1KG Phase 3 chr21 | GRCh37 | 0.97 | 0.95 to 1.04 |
| x86 | 1KG high-cov chr21 | GRCh38 | 1.42 | 1.37 to 1.49 |
| x86 | SV per-VCF (16 files) | GRCh37 | 0.72 | 0.70 to 0.78 |
| x86 | SV per-VCF (16 files) | GRCh38 | 1.53 | 1.50 to 1.63 |

HGVS notation is outside every F1 above. `scripts/concordance/run_clone_measurement.sh` re-measures the chromosome 21 cells on your own machines and `scripts/concordance/run_concordance.sh` times both engines on a directory of your own VCFs; [scripts/README.md](scripts/README.md#reproducing-the-published-concordance) is the runbook for reproducing the concordance numbers.

## Installation

Pre-built binaries for Linux x86_64, Linux aarch64 and macOS aarch64 are on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page. Building from source needs Rust 1.88+ ([rustup](https://rustup.rs/)):

```bash
git clone https://github.com/natera-open-source/vep-rs.git
cd vep-rs
cargo build --release
# Binary at target/release/vep
```

**CPU floor.** `.cargo/config.toml` compiles every build for a fixed CPU baseline: `-C target-cpu=x86-64-v3` on x86_64 Linux (AVX2 required), `-C target-cpu=neoverse-v1` on aarch64 Linux and `-C target-cpu=apple-m1` on Apple-silicon macOS. Release binaries are built the same way, so they need AVX2 on x86_64. A binary built or downloaded this way aborts with an illegal-instruction fault on an older CPU rather than running slowly; to build for the machine you are on, override the flag for that build with `RUSTFLAGS="-C target-cpu=native" cargo build --release` (`RUSTFLAGS` replaces the file's setting, it does not add to it), or edit the `rustflags` line for your target in `.cargo/config.toml`.

### Container image

Each published release from 0.2.0 on is also an image on the GitHub Container Registry, holding the release's own x86_64 Linux binaries, `vep-cache-builder`, `vep-cache-converter`, and the `duckdb` CLI that Parquet output needs:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/tests/golden/116/GRCh37:/corpus:ro" -v "$PWD:/data" \
  ghcr.io/natera-open-source/vep-rs:0.2.0 \
  vep --json_cache /corpus/json_cache --assembly GRCh37 -i /corpus/variants.vcf -o /data/output.txt
```

Tags are `X.Y.Z`, `X.Y` (the newest patch of that line) and `latest`. The image is linux/amd64 only and has the same x86-64-v3 CPU floor as the release binary. On Apple silicon, run it with `--platform linux/amd64` under Rosetta. That works, but it is a translated binary, so don't use it for timing. Each image carries a build-provenance attestation: `gh attestation verify oci://ghcr.io/natera-open-source/vep-rs:0.2.0 --repo natera-open-source/vep-rs`.

## Quick start

Run it first on the bundled corpus, which needs no download:

```bash
# The bundled 1,008-record excerpt of the benchmark datasets against the pruned
# GRCh37 cache that ships in tests/golden/
target/release/vep \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --assembly GRCh37 \
  -i tests/golden/116/GRCh37/variants.vcf -o output.txt
```

For real data you need a JSON transcript cache. Build one natively with `vep-cache-builder` (GRCh38 shown; [docs/cache-setup.md](docs/cache-setup.md) explains each input and gives the GRCh37 command):

```bash
# Ensembl release 116, GRCh38
wget https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

cargo build --release -p vep-cache-builder
target/release/vep-cache-builder \
  --species homo_sapiens --assembly GRCh38 --release 116 \
  --gtf Homo_sapiens.GRCh38.116.gtf.gz \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/json_cache/homo_sapiens/116_GRCh38
```

or convert an existing Perl VEP cache with `scripts/data/storable_to_json.pl` ([scripts/README.md](scripts/README.md)). Then run the same VEP command you run today, with `--json_cache` naming the cache directory (the one that holds `info.json` and `transcripts/`, so `116_GRCh38` itself, not its parent):

```bash
target/release/vep \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/116_GRCh38 \
  --assembly GRCh38 \
  -i input.vcf -o output.txt
```

## Usage

```bash
# With annotation plugins, on data staged by scripts/data/setup_plugin_data.sh
target/release/vep \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/116_GRCh38 \
  --assembly GRCh38 \
  --plugin CADD,snv=$HOME/.vep/plugin_data/grch38/full/cadd/whole_genome_SNVs.tsv.gz \
  --plugin REVEL,file=$HOME/.vep/plugin_data/grch38/full/revel/revel.tsv.gz \
  -i input.vcf -o output.txt --force_overwrite

# JSON output with HGVS notation
target/release/vep \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/116_GRCh38 \
  --assembly GRCh38 --json \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --hgvs --everything \
  -i input.vcf -o output.json
```

- **Output formats.** VEP default (the tab-delimited text format with the `Extra` column), `--vcf` (annotations in the `CSQ` INFO field), `--json` (one object per variant), `--tab`, and `--parquet`, a vep-rs extension written through the `duckdb` CLI (1.2 or newer on PATH). The four VEP formats match Perl VEP's layouts; field-level specs in [docs/output-formats.md](docs/output-formats.md).
- **Plugins.** `--plugin <Name>,<params>` resolves one of the eleven built-in plugins first (ten with supported data; support is per reference genome, and `--assembly` selects the right column and refuses a data file whose own assembly contradicts it) and otherwise loads a C-ABI dynamic library from `--dir_plugins`. The support matrix, the data files and how to write a plugin are in [docs/plugins.md](docs/plugins.md).
- **Flags.** `--fasta` takes an uncompressed, `faidx`-indexed genome and gives `--hgvs` the reference context for the 3' shifting of indels. [docs/cli-reference.md](docs/cli-reference.md) is the complete flag reference, including the flags accepted for Perl VEP compatibility but not implemented.

## Documentation

| Document                                     | Description                                                                 |
| -------------------------------------------- | --------------------------------------------------------------------------- |
| [CLI Reference](docs/cli-reference.md)       | Complete CLI flag reference with examples, and what is not ported from Perl VEP |
| [Output Formats](docs/output-formats.md)     | VEP, VCF, JSON, Tab, and Parquet output specs                               |
| [Plugins](docs/plugins.md)                   | Plugin support matrix, data files, preparation, development, concordance testing |
| [Cache setup](docs/cache-setup.md)           | Building a JSON transcript cache natively or converting a Perl VEP cache; SIFT and PolyPhen matrices |
| [Intended divergences](docs/intended-divergences.md) | The five Perl VEP defects behind the classes the adjusted concordance sets aside, and the one representation difference: for each class the record, both engines' rows, the mechanism in Perl VEP's source, how to reproduce it |
| [Published figures](docs/published-figures.md) | The paper's concordance and wall-time tables, as published, and the released data they re-derive from |
| [Scripts](scripts/README.md)                 | Utility script index and the runbook for reproducing the published concordance |
| [Released data](manuscript/data/README.md)   | The measurement CSVs every published figure re-derives from                 |
| [Contributing](CONTRIBUTING.md)              | Build, test, code conventions, sign-off and pull-request guidelines         |

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for build, test, and pull-request guidelines, and [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) for community expectations. To report a security vulnerability, follow [SECURITY.md](SECURITY.md) rather than filing a public issue.

## Authors

vep-rs was created and is maintained by:

- **Matthew Porter**, Natera, Inc.
- **Robert Borkowski**, Natera, Inc.

## Citing vep-rs

If you use vep-rs in published work, cite the paper:

> Porter M, Borkowski R. _vep-rs: high-throughput Rust variant annotation with population-scale concordance to Ensembl VEP._ bioRxiv (2026). [doi:10.64898/2026.09.22.753614](https://doi.org/10.64898/2026.09.22.753614). The journal DOI will be added here on acceptance.

The figures the paper reports are in [docs/published-figures.md](docs/published-figures.md); the current release's own measurements are in the [Current release](#current-release-v032) section.

The software itself is archived at Zenodo, and that archive is the citation for the code and the released data files:

> Porter M, Borkowski R. _vep-rs_, version 0.1.0. Zenodo (2026). [doi:10.5281/zenodo.22837897](https://doi.org/10.5281/zenodo.22837897)

`CITATION.cff` at the repository root carries the software citation in machine-readable form with the all-versions DOI ([10.5281/zenodo.22837896](https://doi.org/10.5281/zenodo.22837896)) and the archived version's own, so GitHub's "Cite this repository" button and citation managers pick it up directly; the article is added there on acceptance. The concordance and wall-time figures the paper reports re-derive from `manuscript/data/*.csv` in this tree; see [manuscript/data/README.md](manuscript/data/README.md).

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
