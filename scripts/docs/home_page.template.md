<<panel a>>raw F1 per dataset
<<panel b>>wall time per dataset and architecture
<<panel c>>raw F1 per plugin
<<panel d>>the plugins each tool reads
<<panel e>>the ready-made Ensembl VEP caches that are exact, of those checked
<<page>>
# vep-rs

vep-rs is a drop-in replacement for Ensembl VEP in Rust, with the same flags and output formats: it runs up to <<speed-up>> times faster than Ensembl VEP and matches Ensembl VEP exactly except for [the documented divergences](intended-divergences.md).

It does the work of Ensembl VEP, the Variant Effect Predictor: for each variant of a VCF file it writes the genes and transcripts affected and the consequence on each. [Install](#install) a binary, a container, the Bioconda package or the source.

<<if figure>>
<figure class="home-figure">
<a href="assets/comparison.png">

![<<figure-notes>><<panel-count>> comparing <<compared>> against <<reference>>: <<panels>>.](assets/comparison.png)

</a>
<figcaption>

<<figure-notes>><<released>> against <<reference>><<if one-comparator>>, with <<comparators>> as comparator<<end one-comparator>><<if several-comparators>>, with <<comparators>> as comparators<<end several-comparators>>.

</figcaption>
</figure>

<<end figure>>
## Why vep-rs

<div class="home-reasons">

- **Fast on whole genomes.** All of <<gnomad-label>>, <<gnomad-variants>> variants, annotates in <<gnomad-arm-seconds>> s on one ARM machine.
- **Scored against Ensembl VEP.** Release <<reference-release>>, on <<whole-genome-datasets>> whole-genome datasets<<if plugins>> and [<<plugins-read>> plugins](plugins.md#support-table)<<end plugins>>.
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
