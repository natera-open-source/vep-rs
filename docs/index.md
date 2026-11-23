# vep-rs

vep-rs is a drop-in replacement for Ensembl VEP in Rust, with the same flags and output formats: it runs up to 432 times faster than Ensembl VEP and matches Ensembl VEP exactly except for [the documented divergences](intended-divergences.md).

It does the work of Ensembl VEP, the Variant Effect Predictor: for each variant of a VCF file it writes the genes and transcripts affected and the consequence on each. [Install](#install) a binary, a container, the Bioconda package or the source.

<figure class="home-figure">
<a href="assets/comparison.png">

![Counts are of the 31 plugins vep-rs supports, of Ensembl VEP's 87. Plugins vep-rs does not port are not counted; the README lists them by reason. Four panels comparing vep-rs v0.6.0 and fastVEP v0.5.0 against Ensembl VEP 116.2: raw F1 per dataset, wall time per dataset and architecture, raw F1 per plugin and the plugins each tool reads.](assets/comparison.png)

</a>
<figcaption>

Counts are of the 31 plugins vep-rs supports, of Ensembl VEP's 87. Plugins vep-rs does not port are not counted; the README lists them by reason. vep-rs v0.6.0 against Ensembl VEP 116.2, with fastVEP v0.5.0 as comparator.

</figcaption>
</figure>

## Why vep-rs

<div class="home-reasons">

- **Fast on whole genomes.** All of gnomAD v4.1, 759,336,320 variants, annotates in 591 s on one ARM machine.
- **Scored against Ensembl VEP.** Release 116.2, on six whole-genome datasets and [31 plugins](plugins.md#support-table).
- **The same command line.** Ensembl VEP's input files, flag names and output layouts, plus typed Parquet.

</div>

## Install

```bash
https://github.com/natera-open-source/vep-rs/releases   # Binary
docker pull ghcr.io/natera-open-source/vep-rs:latest    # Container
conda install -c conda-forge -c bioconda vep-rs         # Bioconda
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
