# Installation

Install vep-rs in one of four ways: a pre-built `vep` binary for Linux or macOS, a build from source, the container image, or the Bioconda package.

## Pre-built binaries

Binaries for Linux x86_64, Linux aarch64, macOS aarch64 and macOS x86_64 are on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page. Each Linux architecture has two archives: the measured one, built for the CPU floor below, and a generic one that runs on any 64-bit Arm CPU or almost every 64-bit x86 CPU; both run on glibc 2.17 or later. Each macOS architecture has one archive. See [CPU floor](#cpu-floor).

## Build from source

Building needs Rust 1.89+ ([rustup](https://rustup.rs/)):

```bash
git clone https://github.com/natera-open-source/vep-rs.git
cd vep-rs
cargo build --release
# Binary at target/release/vep
```

The same build writes `target/release/vep-cache-builder`, which installs a release's transcript cache and its indexed genome FASTA in one command. That command is `vep-cache-builder --auto homo_sapiens GRCh38 --release 116` ([Cache setup](../cache-setup.md)). On Linux the builder links OpenSSL, so the build needs the OpenSSL development package.

### From crates.io

The workspace's crates are on crates.io from 0.6.0 on, each under the release's version: the release workflow publishes them when a maintainer publishes the release page, so a version appears there minutes after its page. `vep-cli` is the crate of the `vep` command, and `cargo install` then builds that release's `vep` from source with Rust 1.89+:

```bash
cargo install vep-cli --locked
```

The install puts the `vep` binary alone into Cargo's bin directory (`~/.cargo/bin`). The crate named `vep` on crates.io is an unrelated crate: `cargo install vep` does not install vep-rs. `vep-cache-builder` and `vep-cache-converter` are crates of their own, installed the same way; the builder links OpenSSL on Linux, as in a build from a checkout.

A published crate carries neither the repository's `.cargo/config.toml` nor its release profile, so the install compiles for the generic CPU of its target under Cargo's default release profile, not the configuration the release archives are built with. `RUSTFLAGS="-C target-cpu=native"` before the command builds for the machine it runs on ([CPU floor](#cpu-floor)).

## CPU floor

In a checkout of the repository, `.cargo/config.toml` compiles every build for a fixed CPU floor, and the measured release binaries are built the same way. The generic Linux archives are built for a CPU floor almost every machine of their architecture meets:

| Archive | `-C target-cpu` | Runs on |
| --- | --- | --- |
| `x86_64-unknown-linux-gnu` | `x86-64-v3` | AVX2, BMI2 and FMA: Intel Haswell (2013) or later, AMD Zen or later |
| `x86_64-unknown-linux-gnu-v2` | `x86-64-v2` | SSE4.2 and POPCNT: Intel Nehalem (2008) or later, AMD Bulldozer (2011) or later, VMs that hide AVX2 included |
| `aarch64-unknown-linux-gnu` | `neoverse-v1` | SVE: AWS Graviton3 and later, Arm Neoverse V1, V2 and N2 |
| `aarch64-unknown-linux-gnu-generic` | `generic` | any 64-bit Arm: Graviton2, Ampere Altra, Raspberry Pi, Apple silicon under Linux |
| `aarch64-apple-darwin` | `apple-m1` | any Apple M-series Mac |
| `x86_64-apple-darwin` | `x86-64-v3` | AVX2, BMI2 and FMA: every Intel Mac macOS 26 supports |

- On an older CPU, a binary built or downloaded this way aborts with an illegal-instruction fault rather than running slowly. The generic archives are slower than the measured ones; each release page reports the measured difference, or that the release did not time them.
- To build for the machine you are on, override the flag for that build. `RUSTFLAGS` replaces the file's setting; it does not add to it. You can also edit the `rustflags` line for your target in `.cargo/config.toml`.

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## Container image

Each published release from 0.2.0 on is also an image on the GitHub Container Registry; from 0.5.0 it is a multi-platform image. The `X.Y.Z`, `X.Y` and `latest` tags carry the generic Linux binaries (`x86-64-v2` on linux/amd64, generic aarch64 on linux/arm64), so `docker run` works on any 64-bit Arm host and almost every 64-bit x86 host, Apple silicon and Graviton2 included, with no `--platform` flag. The `X.Y.Z-x86-64-v3` and `X.Y.Z-neoverse-v1` tags carry the measured binaries for hosts that meet their CPU floors.

Every image holds the release's own `vep`, `vep-cache-builder` and `vep-cache-converter`; Parquet output needs no other tool. Images from 0.5.0 add `bgzip` and `tabix` from Debian's `tabix` package, so `--compress_output bgzip` VCF output can be indexed in the container. Images from 0.6.0 add `ps` from Debian's `procps` package, which Nextflow's task wrapper calls inside the container.

The command below runs from a checkout of the repository, whose `tests/golden/116/GRCh37` directory holds the pruned cache and the VCF it reads:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/tests/golden/116/GRCh37:/corpus:ro" -v "$PWD:/data" \
  ghcr.io/natera-open-source/vep-rs:0.6.0 \
  vep --json_cache /corpus/json_cache --assembly GRCh37 \
  -i /corpus/variants.vcf -o /data/output.txt
```

- Tags are `X.Y.Z`, `X.Y` (the newest patch of that line) and `latest`, each a manifest list over linux/amd64 and linux/arm64, plus the single-platform `X.Y.Z-x86-64-v3` (linux/amd64) and `X.Y.Z-neoverse-v1` (linux/arm64).
- The `-x86-64-v3` tag needs AVX2 and the `-neoverse-v1` tag SVE (AWS Graviton3 and later); below those floors they abort with an illegal-instruction fault.
- On Apple silicon the multi-platform tags run natively as linux/arm64. Their speed there is not measured, so don't cite it for timing.

Each image carries a build-provenance attestation:

```bash
gh attestation verify oci://ghcr.io/natera-open-source/vep-rs:0.6.0 \
  --repo natera-open-source/vep-rs
```

### Singularity and Apptainer

Apptainer and Singularity run the same image from a single file. `apptainer pull docker://ghcr.io/natera-open-source/vep-rs:X.Y.Z` writes `vep-rs_X.Y.Z.sif` from the image for the host's architecture: the generic binary under the multi-platform tags, the measured one under `X.Y.Z-x86-64-v3` or `X.Y.Z-neoverse-v1` on a host that meets its CPU floor. The container runs as the invoking user with the host's current directory, home directory and `/tmp` mounted, so the image's own user and `/data` working directory do not apply. A cache or output directory anywhere else is mounted with `--bind`:

```bash
apptainer pull docker://ghcr.io/natera-open-source/vep-rs:0.6.0
apptainer exec --bind /shared/vep-cache:/cache:ro --bind "$PWD/results:/out" \
  vep-rs_0.6.0.sif \
  vep --json_cache /cache/json_cache --assembly GRCh37 \
  -i variants.vcf -o /out/output.txt
```

`singularity pull` and `singularity exec` take the same arguments.

## Bioconda

The Bioconda channel carries vep-rs as the package `vep-rs`, built for Linux and macOS on x86_64 and arm64. It installs `vep`, `vep-cache-builder` and `vep-cache-converter`:

```bash
conda install -c conda-forge -c bioconda vep-rs
```

The channel builds a release from its own copy of the recipe, `recipes/vep-rs` in [bioconda-recipes](https://github.com/bioconda/bioconda-recipes), once that copy names the release, so the version it serves can trail the newest release. `conda search -c bioconda vep-rs` lists the versions the channel holds.

The repository holds a Bioconda recipe, [`bioconda/vep-rs`](../../bioconda/vep-rs/meta.yaml). It builds `vep`, `vep-cache-builder` and `vep-cache-converter` from the source tarball of the published release its `version` names (`vep-rs-X.Y.Z.tar.gz`, the `git archive` the release workflow cuts) with `cargo build --release --locked`.

Linux x86_64 compiles for `x86-64-v2` and Linux aarch64 for a generic CPU, the floors of the generic release archives. `.cargo/config.toml` would otherwise compile for `x86-64-v3` and Neoverse V1, which abort with an illegal-instruction fault below those floors. On Linux the build links conda-forge OpenSSL, the same library the release binaries link.

[What vep-rs does not do](what-vep-rs-does-not-do.md#bioconda) says what a pipeline that installs Ensembl VEP from Bioconda changes.

## Operating systems and architectures

vep-rs runs where a release archive or the container image runs. Every archive holds the same three binaries; the measured archives are the ones the release measurement times, and the generic ones trade speed for a lower CPU floor ([CPU floor](#cpu-floor)).

| Platform | Archive or image | Runs on |
| --- | --- | --- |
| Linux x86_64 | `x86_64-unknown-linux-gnu` (measured, `x86-64-v3`) and `x86_64-unknown-linux-gnu-v2` (generic, `x86-64-v2`) | glibc 2.17 or later |
| Linux aarch64 | `aarch64-unknown-linux-gnu` (measured, `neoverse-v1`) and `aarch64-unknown-linux-gnu-generic` (generic, Armv8.0-A) | glibc 2.17 or later |
| macOS on Apple silicon | `aarch64-apple-darwin` | macOS 11 or later |
| macOS on Intel | `x86_64-apple-darwin` (`x86-64-v3`) | an Intel Mac with AVX2 (every Intel Mac macOS 26 supports) |
| Container image | `ghcr.io/natera-open-source/vep-rs` on linux/amd64 and linux/arm64 | Docker, Podman, Singularity or Apptainer |
| Windows | none | the container image under Docker, or the Linux x86_64 archive under WSL2; neither is tested on Windows |
| FreeBSD, Linux with musl | none | not built |

There is no Homebrew formula: macOS users install from the release archive or build from source. The two macOS archives are not timed by the release measurement, so no page states their speed.

## Memory

The transcript index is loaded one chromosome at a time, and every chromosome a run touches stays resident until the run ends. Variation files are read per chromosome as needed, and an identifier input holds every name of the list while it reads the variation files once. So an input on one chromosome peaks at that chromosome's index plus the batches in flight (`--buffer_size` records each, 5,000 by default). The index is the part of the peak that grows with the chromosomes an input touches.

The README's Peak memory row is the largest per-chromosome median of the release measurement, one chromosome per input file at 16 threads. A whole genome in one input file holds every chromosome's index at once. The release page states the measured peak of gnomAD v4.1 GRCh38 as a single file per architecture, beside the per-dataset peaks of its wall-time table; the README's provenance records carry each cell's peak.

The next step is the [Quick start](quick-start.md), which runs the bundled example and then builds a transcript cache.
