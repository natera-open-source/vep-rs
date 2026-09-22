# vep-rs

A fast, memory-efficient variant effect predictor written in Rust: a from-scratch reimplementation of the Ensembl Variant Effect Predictor (VEP) that takes the same inputs and flags, emits VEP's consequence vocabulary in VEP's output formats, and annotates population-scale callsets two orders of magnitude faster than the Perl implementation.

[![CI](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/natera-open-source/vep-rs/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![MSRV](https://img.shields.io/badge/MSRV-1.88-orange.svg)](#installation)

vep-rs reads a JSON transcript cache, built from Ensembl's release files by `vep-cache-builder` or converted from an existing Perl VEP cache, and writes VEP's default, VCF, JSON and Tab formats plus Parquet. Throughout this README "Perl VEP" is the reference implementation, `ensemblorg/ensembl-vep` release 115.2; the accompanying paper calls the same engine "VEP".

## Why vep-rs instead of Ensembl VEP or fastVEP

|                    |                                                                                                                                                                                                                                 |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Performance**    | Faster than Perl VEP on every cell (per-cell 101x to 284x ARM and 78.3x to 192x x86, N=20 on ARM Graviton4 and x86 Intel) and faster than fastVEP on every cell (geomean 9.07x ARM / 6.55x x86); peak RSS 0.62 to 4.56 GiB, lighter than fastVEP on four of the ten measured cells on ARM and five on x86 |
| **Concordance**    | Adjusted F1 = 1.000000 with Perl VEP on all six SNP/indel datasets (GRCh37 + GRCh38); raw F1 ≥ 0.999974, with exactly Perl VEP's consequence-tuple count on every dataset                                                             |
| **Plugins**        | 10 supported built-in databases: 8 on both GRCh37 and GRCh38 (CADD, REVEL, gnomADc, AlphaMissense, dbscSNV, LoFTEE, LoFtool, pLI) plus GWAS and SpliceAI on GRCh38, and C-ABI dynamic-library plugins loaded from `--dir_plugins`; see [docs/plugins.md](docs/plugins.md) |
| **Output formats** | VEP default, VCF, JSON, Tab, Parquet (Parquet needs the `duckdb` CLI on PATH)                                                                                                                                                    |
| **Threading**      | Parallel annotation via `--fork N`; the default uses every logical CPU                                                                                                                                                          |

**What drop-in means here.** The same input files, the same flag names, the same output layouts. Two things differ from Perl VEP. The transcript cache is `--json_cache`, a directory built by `vep-cache-builder` or converted from a Perl VEP cache ([docs/cache-setup.md](docs/cache-setup.md)); `--cache`, `--offline`, `--dir_cache`, `--cache_version` and `--species` are accepted, but they do not locate a cache, and vep-rs never contacts a database. And some Perl VEP flags are accepted for compatibility without being implemented, among them `--regulatory`, `--custom`, `--refseq`, `--merged`, `--coding_only`, `--pick_order` and the statistics flags (vep-rs writes no summary file); [docs/cli-reference.md](docs/cli-reference.md) marks every such flag.

## Installation

Pre-built binaries for Linux x86_64, Linux aarch64 and macOS aarch64 are on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page. Building from source needs Rust 1.88+ ([rustup](https://rustup.rs/)):

```bash
git clone https://github.com/natera-open-source/vep-rs.git
cd vep-rs
cargo build --release
# Binary at target/release/vep
```

**CPU floor.** `.cargo/config.toml` compiles every build for a fixed CPU baseline: `-C target-cpu=x86-64-v3` on x86_64 Linux (AVX2 required), `-C target-cpu=neoverse-v1` on aarch64 Linux and `-C target-cpu=apple-m1` on Apple-silicon macOS. Release binaries are built the same way, so they need AVX2 on x86_64. A binary built or downloaded this way aborts with an illegal-instruction fault on an older CPU rather than running slowly; to build for the machine you are on, override the flag for that build with `RUSTFLAGS="-C target-cpu=native" cargo build --release` (`RUSTFLAGS` replaces the file's setting, it does not add to it), or edit the `rustflags` line for your target in `.cargo/config.toml`.

## Quick start

Run it first on the bundled corpus, which needs no download:

```bash
# The bundled 1,008-record excerpt of the benchmark datasets against the pruned
# GRCh37 cache that ships in tests/golden/
target/release/vep \
  --json_cache tests/golden/115/GRCh37/json_cache \
  --assembly GRCh37 \
  -i tests/golden/115/GRCh37/variants.vcf -o output.txt
```

For real data you need a JSON transcript cache. Build one natively with `vep-cache-builder` (GRCh38 shown; [docs/cache-setup.md](docs/cache-setup.md) explains each input and gives the GRCh37 command):

```bash
# Ensembl release 115, GRCh38
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

cargo build --release -p vep-cache-builder
target/release/vep-cache-builder \
  --species homo_sapiens --assembly GRCh38 --release 115 \
  --gtf Homo_sapiens.GRCh38.115.gtf.gz \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/json_cache/homo_sapiens/115_GRCh38
```

or convert an existing Perl VEP cache with `scripts/data/storable_to_json.pl` ([scripts/README.md](scripts/README.md)). Then run the same VEP command you run today, with `--json_cache` naming the cache directory (the one that holds `info.json` and `transcripts/`, so `115_GRCh38` itself, not its parent):

```bash
target/release/vep \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/115_GRCh38 \
  --assembly GRCh38 \
  -i input.vcf -o output.txt
```

## Usage

```bash
# With annotation plugins, on data staged by scripts/data/setup_plugin_data.sh
target/release/vep \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/115_GRCh38 \
  --assembly GRCh38 \
  --plugin CADD,snv=$HOME/.vep/plugin_data/grch38/full/cadd/whole_genome_SNVs.tsv.gz \
  --plugin REVEL,file=$HOME/.vep/plugin_data/grch38/full/revel/revel.tsv.gz \
  -i input.vcf -o output.txt --force_overwrite

# JSON output with HGVS notation
target/release/vep \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/115_GRCh38 \
  --assembly GRCh38 --json \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --hgvs --everything \
  -i input.vcf -o output.json
```

- **Output formats.** VEP default (the tab-delimited text format with the `Extra` column), `--vcf` (annotations in the `CSQ` INFO field), `--json` (one object per variant), `--tab`, and `--parquet`, a vep-rs extension written through the `duckdb` CLI (1.2 or newer on PATH). The four VEP formats match Perl VEP's layouts; field-level specs in [docs/output-formats.md](docs/output-formats.md).
- **Plugins.** `--plugin <Name>,<params>` resolves one of the eleven built-in plugins first (ten with supported data; support is per reference genome, and `--assembly` selects the right column and refuses a data file whose own assembly contradicts it) and otherwise loads a C-ABI dynamic library from `--dir_plugins`. The support matrix, the data files and how to write a plugin are in [docs/plugins.md](docs/plugins.md).
- **Flags.** `--fasta` takes an uncompressed, `faidx`-indexed genome and gives `--hgvs` the reference context for the 3' shifting of indels. [docs/cli-reference.md](docs/cli-reference.md) is the complete flag reference, including the flags accepted for Perl VEP compatibility but not implemented.

## Results

Perl VEP release 115.2 is the reference engine in every table; fastVEP is the other Rust VEP port. Every figure here, and every figure in the table at the top, re-derives from the released measurement data in [`manuscript/data/`](manuscript/data/README.md). Concordance is measured on the complete consequence output of each dataset, tuple by tuple on matched caches, and F1 is architecture-independent, so one value per dataset covers both ARM and x86. **Adjusted F1** sets aside the tuples of five documented divergence classes, four of them Perl VEP defects vep-rs does not reproduce; on the six SNP/indel datasets the mask removes matched pairs from both engines' sides in equal numbers, the covered `splice_region_variant` swap (all 1,884 pairs on the two ClinVar datasets) and the `start_lost`/`start_retained_variant` co-emission. The classes, the mask and the structural-variant results per input file are in the accompanying paper's supplement (S2, S4.3 and Table S5), and every value in them re-derives from [`manuscript/data/`](manuscript/data/README.md).

| Dataset             | Assembly | vep-rs Raw F1 | vep-rs Adj F1 | fastVEP Raw F1 |
| ------------------- | -------- | ------------- | ------------- | -------------- |
| ClinVar full        | GRCh37   | 0.999979      | 1.000000      | 0.825597       |
| ClinVar full        | GRCh38   | 0.999974      | 1.000000      | 0.969701       |
| gnomAD v2.1.1 chr21 | GRCh37   | 0.999999      | 1.000000      | 0.762022       |
| gnomAD v4.1 chr21   | GRCh38   | 0.999999      | 1.000000      | 0.919675       |
| 1KG Phase 3 chr21   | GRCh37   | 1.000000      | 1.000000      | 0.790928       |
| 1KG high-cov chr21  | GRCh38   | 1.000000      | 1.000000      | 0.938604       |
| SV per-VCF (16 files) | GRCh37 | 0.975395      | 0.998524      | 0.153626       |
| SV per-VCF (16 files) | GRCh38 | 0.909998      | 0.998031      | 0.133436       |

Wall time is the median of 20 independent machine measurements per cell, each on a fresh instance after a discarded warmup, on sites-only inputs; vep-rs and Perl VEP run `--fork 16`, fastVEP parallelizes internally at 16 threads. vep-rs annotates full ClinVar GRCh37 (4,388,172 input records) in 6.00 s on ARM and 6.20 s on x86, and its peak resident memory spans 0.62 to 4.56 GiB across every measured cell.

| Arch | Dataset             | Perl (s) | fastVEP (s) | vep-rs (s) | vep-rs vs Perl | vep-rs vs fastVEP |
| ---- | ------------------- | -------- | ----------- | ---------- | -------------- | ----------------- |
| ARM  | ClinVar GRCh37      | 961.08   | 43.65       | 6.00       | 160×           | 7.28×             |
| ARM  | ClinVar GRCh38      | 2,722.01 | 121.13      | 12.35      | 220×           | 9.81×             |
| ARM  | gnomAD v2.1.1 chr21 | 1,500.91 | 68.24       | 6.84       | 219×           | 9.98×             |
| ARM  | gnomAD v4.1 chr21   | 3,953.00 | 193.44      | 13.92      | 284×           | 13.9×             |
| ARM  | 1KG Phase 3 chr21   | 106.29   | 6.26        | 1.05       | 101×           | 5.96×             |
| ARM  | 1KG high-cov chr21  | 249.15   | 17.45       | 1.85       | 135×           | 9.43×             |
| x86  | ClinVar GRCh37      | 849.97   | 33.16       | 6.20       | 137×           | 5.35×             |
| x86  | ClinVar GRCh38      | 2,378.35 | 90.26       | 12.56      | 189×           | 7.19×             |
| x86  | gnomAD v2.1.1 chr21 | 1,341.62 | 57.69       | 10.16      | 132×           | 5.68×             |
| x86  | gnomAD v4.1 chr21   | 3,276.52 | 158.24      | 17.06      | 192×           | 9.28×             |
| x86  | 1KG Phase 3 chr21   | 85.31    | 5.10        | 1.09       | 78.3×          | 4.68×             |
| x86  | 1KG high-cov chr21  | 205.16   | 14.60       | 1.75       | 118×           | 8.36×             |

Across those cells vep-rs is faster than Perl VEP by a geometric mean of 176× (ARM) / 135× (x86) and faster than fastVEP by 9.07× (ARM) / 6.55× (x86). vep-rs was timed in one campaign and the comparators' SNP and indel cells in another, on the same instance types against the same dataset inventory, so every ratio in the table divides medians from separate machines and campaigns; the structural-variant cells of all three engines come from the vep-rs campaign. The structural-variant sets are annotated as 16 separate per-file invocations, so per-invocation start-up weighs far more there; vep-rs is faster than Perl VEP on those cells by 52.9× (GRCh37) and 23.1× (GRCh38) on ARM and 37.5× and 22.5× on x86 (0.61 s, 2.81 s, 0.79 s and 2.63 s against Perl's 32.28 s, 64.92 s, 29.62 s and 59.29 s). No fastVEP structural-variant speedup is given, because fastVEP emits 62% and 88% of Perl VEP's tuple volume on those sets and recovers only 12.5% and 12.7% of Perl VEP's tuples, so its wall time does not buy comparable annotations. Plugin output and HGVS notation are outside every F1 above.

`scripts/concordance/run_clone_measurement.sh` re-measures these cells on your own machines and `scripts/concordance/run_concordance.sh` times both engines on a directory of your own VCFs; [scripts/README.md](scripts/README.md#reproducing-the-published-concordance) is the runbook for reproducing the concordance numbers.

## Documentation

| Document                                     | Description                                                                 |
| -------------------------------------------- | --------------------------------------------------------------------------- |
| [CLI Reference](docs/cli-reference.md)       | Complete CLI flag reference with examples, and what is not ported from Perl VEP |
| [Output Formats](docs/output-formats.md)     | VEP, VCF, JSON, Tab, and Parquet output specs                               |
| [Plugins](docs/plugins.md)                   | Plugin support matrix, data files, preparation, development, concordance testing |
| [Cache setup](docs/cache-setup.md)           | Building a JSON transcript cache natively or converting a Perl VEP cache; SIFT and PolyPhen matrices |
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

> Porter M, Borkowski R. _vep-rs: high-throughput Rust variant annotation with population-scale concordance to Ensembl VEP._ The journal DOI will be added here on acceptance.

The software itself is archived at Zenodo, and that archive is the citation for the code and the released data files:

> Porter M, Borkowski R. _vep-rs_, version 0.1.0. Zenodo (2026). [doi:10.5281/zenodo.22837897](https://doi.org/10.5281/zenodo.22837897)

`CITATION.cff` at the repository root carries the software citation and its DOI in machine-readable form, so GitHub's "Cite this repository" button and citation managers pick it up directly; the article is added there on acceptance. The concordance and wall-time figures the paper reports re-derive from `manuscript/data/*.csv` in this tree; see [manuscript/data/README.md](manuscript/data/README.md).

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
