# Installation

This page covers the three ways to get a `vep` binary: a pre-built release, a build from source, and the container image.

## Pre-built binaries

Binaries for Linux x86_64, Linux aarch64 and macOS aarch64 are on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page. The x86_64 binary needs AVX2; see [CPU floor](#cpu-floor).

## Build from source

Building needs Rust 1.88+ ([rustup](https://rustup.rs/)):

```bash
git clone https://github.com/natera-open-source/vep-rs.git
cd vep-rs
cargo build --release
# Binary at target/release/vep
```

## CPU floor

`.cargo/config.toml` compiles every build for a fixed CPU baseline, and release binaries are built the same way:

| Target | `-C target-cpu` |
| --- | --- |
| x86_64 Linux | `x86-64-v3` (AVX2 required) |
| aarch64 Linux | `neoverse-v1` |
| Apple-silicon macOS | `apple-m1` |

- On an older CPU, a binary built or downloaded this way aborts with an illegal-instruction fault rather than running slowly.
- To build for the machine you are on, override the flag for that build. `RUSTFLAGS` replaces the file's setting; it does not add to it. You can also edit the `rustflags` line for your target in `.cargo/config.toml`.

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## Container image

Each published release from 0.2.0 on is also an image on the GitHub Container Registry. It holds the release's own x86_64 Linux binaries, `vep-cache-builder`, `vep-cache-converter`, and the `duckdb` CLI that Parquet output needs.

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/tests/golden/116/GRCh37:/corpus:ro" -v "$PWD:/data" \
  ghcr.io/natera-open-source/vep-rs:0.4.0 \
  vep --json_cache /corpus/json_cache --assembly GRCh37 \
  -i /corpus/variants.vcf -o /data/output.txt
```

- Tags are `X.Y.Z`, `X.Y` (the newest patch of that line) and `latest`.
- The image is linux/amd64 only and has the same x86-64-v3 CPU floor as the release binary.
- On Apple silicon, run it with `--platform linux/amd64` under Rosetta. That works, but it is a translated binary, so don't use it for timing.

Each image carries a build-provenance attestation:

```bash
gh attestation verify oci://ghcr.io/natera-open-source/vep-rs:0.4.0 \
  --repo natera-open-source/vep-rs
```

The next step is the [Quick start](quick-start.md), which runs the bundled corpus and then builds a transcript cache.
