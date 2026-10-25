# Installation

This page covers the three ways to get a `vep` binary: a pre-built release, a build from source, and the container image.

## Pre-built binaries

Binaries for Linux x86_64, Linux aarch64 and macOS aarch64 are on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page. Each Linux architecture has two archives: the measured one, built for the CPU floor below, and a generic one that runs on any 64-bit Arm CPU or almost every 64-bit x86 CPU; both run on glibc 2.17 or later. See [CPU floor](#cpu-floor).

## Build from source

Building needs Rust 1.89+ ([rustup](https://rustup.rs/)):

```bash
git clone https://github.com/natera-open-source/vep-rs.git
cd vep-rs
cargo build --release
# Binary at target/release/vep
```

## CPU floor

`.cargo/config.toml` compiles every build for a fixed CPU baseline, and the measured release binaries are built the same way; the generic Linux archives are built for a CPU floor almost every machine of their architecture meets:

| Archive | `-C target-cpu` | Runs on |
| --- | --- | --- |
| `x86_64-unknown-linux-gnu` | `x86-64-v3` | AVX2, BMI2 and FMA: Intel Haswell (2013) or later, AMD Zen or later |
| `x86_64-unknown-linux-gnu-v2` | `x86-64-v2` | SSE4.2 and POPCNT: Intel Nehalem (2008) or later, AMD Bulldozer (2011) or later, VMs that hide AVX2 included |
| `aarch64-unknown-linux-gnu` | `neoverse-v1` | SVE: AWS Graviton3 and later, Arm Neoverse V1, V2 and N2 |
| `aarch64-unknown-linux-gnu-generic` | `generic` | any 64-bit Arm: Graviton2, Ampere Altra, Raspberry Pi, Apple silicon under Linux |
| `aarch64-apple-darwin` | `apple-m1` | any Apple M-series Mac |

- On an older CPU, a binary built or downloaded this way aborts with an illegal-instruction fault rather than running slowly. The generic archives are slower than the measured ones; each release page reports the measured difference, or that the release did not time them.
- To build for the machine you are on, override the flag for that build. `RUSTFLAGS` replaces the file's setting; it does not add to it. You can also edit the `rustflags` line for your target in `.cargo/config.toml`.

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## Container image

Each published release from 0.2.0 on is also an image on the GitHub Container Registry; from 0.5.0 it is a multi-platform image. The `X.Y.Z`, `X.Y` and `latest` tags carry the generic Linux binaries (`x86-64-v2` on linux/amd64, generic aarch64 on linux/arm64), so `docker run` works on any 64-bit Arm host and almost every 64-bit x86 host, Apple silicon and Graviton2 included, with no `--platform` flag. The `X.Y.Z-x86-64-v3` and `X.Y.Z-neoverse-v1` tags carry the measured binaries for hosts that meet their CPU floors. Every image holds the release's own `vep`, `vep-cache-builder` and `vep-cache-converter`, and the `duckdb` CLI that Parquet output needs. Images from 0.5.0 add `bgzip` and `tabix` from Debian's `tabix` package, so `--compress_output bgzip` VCF output can be indexed in the container.

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/tests/golden/116/GRCh37:/corpus:ro" -v "$PWD:/data" \
  ghcr.io/natera-open-source/vep-rs:0.5.0 \
  vep --json_cache /corpus/json_cache --assembly GRCh37 \
  -i /corpus/variants.vcf -o /data/output.txt
```

- Tags are `X.Y.Z`, `X.Y` (the newest patch of that line) and `latest`, each a manifest list over linux/amd64 and linux/arm64, plus the single-platform `X.Y.Z-x86-64-v3` (linux/amd64) and `X.Y.Z-neoverse-v1` (linux/arm64).
- The multi-platform tags run on any 64-bit Arm CPU and almost every 64-bit x86 CPU. Their generic binaries are slower than the measured ones; each release page reports the measured difference, or that the release did not time them. The `-x86-64-v3` tag needs AVX2 and the `-neoverse-v1` tag SVE (AWS Graviton3 and later); below those floors they abort with an illegal-instruction fault.
- On Apple silicon the multi-platform tags run natively as linux/arm64. Their speed there is not measured, so don't cite it for timing.

Each image carries a build-provenance attestation:

```bash
gh attestation verify oci://ghcr.io/natera-open-source/vep-rs:0.5.0 \
  --repo natera-open-source/vep-rs
```

The next step is the [Quick start](quick-start.md), which runs the bundled corpus and then builds a transcript cache.
