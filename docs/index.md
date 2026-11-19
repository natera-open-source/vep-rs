# vep-rs

vep-rs is a drop-in replacement for Ensembl VEP in Rust, with the same flags and output formats: it runs up to 411 times faster than Ensembl VEP and matches Ensembl VEP exactly except for [the documented divergences](intended-divergences.md).

It does the work of Ensembl VEP, the Variant Effect Predictor: for each variant of a VCF file it writes the genes and transcripts affected and the consequence on each. [Install](#install) a binary, a container or the source.

## Why vep-rs

<div class="home-reasons">

- **Fast on whole genomes.** All of gnomAD v4.1, 759,336,320 variants, annotates in 624 s on one ARM machine.
- **Scored against Ensembl VEP.** Release 116.2, on six whole-genome datasets.
- **The same command line.** Ensembl VEP's input files, flag names and output layouts, plus typed Parquet.

</div>

## Install

```bash
https://github.com/natera-open-source/vep-rs/releases   # Binary
docker pull ghcr.io/natera-open-source/vep-rs:latest    # Container
cargo build --release                                   # Source
```

[Installation](getting-started/installation.md) covers the CPU floor, the container tags, the Bioconda recipe and building for your own machine.

## Try it

The bundled example in `examples/` runs against the pruned cache in `tests/golden/` with no download:

```bash
target/release/vep \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --assembly GRCh37 \
  -i examples/homo_sapiens_GRCh37.vcf -o output.txt
```

## Where next

- [Quick start](getting-started/quick-start.md): a transcript cache, your own data and plugins.
- [Coming from Ensembl VEP](getting-started/coming-from-ensembl-vep.md): what stays the same and what changes.
- [CLI reference](cli-reference.md): every flag, output format and plugin.
- [Benchmarks](benchmarks/current-release.md): the release's speed and concordance tables.
- [Intended divergences](intended-divergences.md): where the two engines differ, and why.
- [What vep-rs does not do](getting-started/what-vep-rs-does-not-do.md): the Ensembl VEP features it lacks.

This site describes vep-rs {{vep-rs-ref}}; the latest release is {{vep-rs-release}}.

If you use vep-rs in published work, cite the paper and the software archive as described on the [Citing vep-rs](project/citing.md) page.
