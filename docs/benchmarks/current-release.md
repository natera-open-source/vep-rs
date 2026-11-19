# Benchmarks against Ensembl VEP and fastVEP

vep-rs benchmarks for the current release: wall time against Ensembl VEP and fastVEP and concordance with Ensembl VEP, measured on whole-genome datasets.

<!-- readme-section: Why vep-rs | fragments: benchmarks/current-release.md -->

<!-- readme-section: Current release -->

The figures the paper reports, measured at its pinned build, are on [the paper's figures page](../published-figures.md) and stay as published.

Each release page on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page states each dataset's peak RSS per architecture beside its wall time. That figure is the largest per-chromosome median of the release measurement. Where the release measured it, the page also states the peak of a whole genome annotated as one input file. [Installation](../getting-started/installation.md#memory) states what those figures cover and how a run's memory grows with the chromosomes an input touches.
