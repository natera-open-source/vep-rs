# vep-rs

A Rust reimplementation of the Ensembl Variant Effect Predictor (VEP): the same inputs, flags and output formats, two orders of magnitude faster.

This site describes vep-rs {{vep-rs-ref}}; the latest release is {{vep-rs-release}}.

{{#include assets/benchmark.svg}}

Wall time on ARM under the protocol the chart's title names; [Benchmarks](benchmarks/current-release.md) has the x86 figures and the method.

```bash
https://github.com/natera-open-source/vep-rs/releases   # Binary
docker pull ghcr.io/natera-open-source/vep-rs:latest    # Container
cargo build --release                                   # Source
```

The bundled corpus in `tests/golden/` runs with no download:

```bash
target/release/vep \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --assembly GRCh37 \
  -i tests/golden/116/GRCh37/variants.vcf -o output.txt
```

[Installation](getting-started/installation.md) covers the CPU floor, the container tags and building for your own machine. The [Quick start](getting-started/quick-start.md) continues with a transcript cache and your own VCF. [What vep-rs does not do](getting-started/what-vep-rs-does-not-do.md) lists the Ensembl VEP features vep-rs does not have, what a run does when asked for one, and which of them the next release, 0.6.0, adds.

| [Getting started](getting-started/installation.md) | [Reference](cli-reference.md) |
| --- | --- |
| Install, run the bundled corpus, build a cache, and move over from Ensembl VEP. | Every flag, output format and plugin, one page per group. |

| [Concordance with Ensembl VEP](intended-divergences.md) | [Benchmarks](benchmarks/current-release.md) |
| --- | --- |
| Where the two engines' outputs are known to differ, and why each difference is intended. | The current release's speed and concordance figures, measured under the paper's protocol. |

<!-- readme-section: Why vep-rs | fragments: benchmarks/current-release.md -->

If you use vep-rs in published work, cite the paper and the software archive as described on the [Citing vep-rs](project/citing.md) page.
