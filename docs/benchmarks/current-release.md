# Benchmarks against Ensembl VEP and fastVEP

vep-rs benchmarks for the current release: wall time against Ensembl VEP and fastVEP and concordance with Ensembl VEP, measured on whole-genome datasets.

<!-- readme-section: Why vep-rs | fragments: benchmarks/current-release.md -->

<!-- readme-section: Current release -->

The figures the paper reports, measured at its pinned build, are on [the paper's figures page](../published-figures.md) and stay as published.

Each release page on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page states each dataset's peak RSS per architecture beside its wall time. That figure is the largest per-chromosome median of the release measurement. Where the release measured it, the page also states the peak of a whole genome annotated as one input file. [Installation](../getting-started/installation.md#memory) states what those figures cover and how a run's memory grows with the chromosomes an input touches.

Two further records under `docs/concordance-provenance/`, their names ending `-no-extra-args.json`, hold fastVEP v0.6.2 run without `--exclude-readthrough --exclude-artifact`: one its ClinVar and structural-variant concordance, the other its ClinVar wall times. No table above reads them.

<!-- plugin-wall-time: begin -->

## Wall time with each plugin

Seconds for one annotation run of the genome-wide ClinVar dataset of each assembly: 4,388,172 variants on GRCh37 and 4,398,127 variants on GRCh38. A plugin's row is the run of the first row, `No plugin`, with that one plugin added. HGVSIntronOffset, which reads `--hgvs`, is added to the run under `--hgvs`. The 20 plugins timed are verified on real data, and the [Plugins](#plugins) tables score them on ClinVar: a plugin whose data cannot be redistributed is not timed, and neither is TranscriptAnnotator, which reads a table of the user's own.

| Plugin | Assembly | Ensembl VEP 116.3, ARM (s) | vep-rs v1.0.0, ARM (s) | Ensembl VEP 116.3, x86 (s) | vep-rs v1.0.0, x86 (s) |
| --- | --- | --- | --- | --- | --- |
| No plugin | GRCh37 | 1,271.36 | 2.06 | 1,026.69 | 3.07 |
| No plugin | GRCh38 | 3,708.45 | 7.69 | 3,516.17 | 9.99 |
| No plugin, `--hgvs` | GRCh37 | 1,586.09 | 6.30 | 1,294.86 | 9.66 |
| No plugin, `--hgvs` | GRCh38 | 4,222.00 | 23.88 | 4,688.72 | 36.63 |
| AlphaMissense | GRCh37 | 1,380.81 | 41.89 | 1,212.23 | 55.70 |
| AlphaMissense | GRCh38 | 3,746.75 | 47.26 | 3,760.73 | 62.10 |
| AncestralAllele | GRCh37 | 1,444.06 | 6.24 | 1,079.60 | 8.59 |
| AncestralAllele | GRCh38 | 3,460.03 | 15.92 | 3,651.43 | 20.49 |
| Blosum62 | GRCh37 | 1,470.12 | 4.83 | 1,081.77 | 8.33 |
| Blosum62 | GRCh38 | 3,458.85 | 18.06 | 3,625.24 | 24.97 |
| Conservation | GRCh38 | 5,626.90 | 16.65 | 5,799.53 | 20.66 |
| DosageSensitivity | GRCh37 | 1,544.92 | 8.85 | 1,053.31 | 11.44 |
| DosageSensitivity | GRCh38 | 3,849.89 | 31.35 | 3,727.71 | 39.41 |
| Downstream | GRCh37 | 1,473.06 | 21.27 | 1,100.73 | 50.39 |
| Downstream | GRCh38 | 3,858.43 | 70.50 | 3,856.49 | 191.83 |
| gnomADc | GRCh37 | 2,181.00 | 210.10 | 2,206.56 | 238.68 |
| gnomADc | GRCh38 | 4,869.35 | 80.84 | 4,416.30 | 99.10 |
| gnomADMt | GRCh38 | 4,050.22 | 10.30 | 3,541.59 | 14.24 |
| GWAS | GRCh38 | 4,079.26 | 12.13 | 3,592.02 | 17.73 |
| HGVSIntronOffset | GRCh37 | 1,547.11 | 13.75 | 1,347.29 | 19.98 |
| HGVSIntronOffset | GRCh38 | 4,957.34 | 53.53 | 5,103.18 | 77.89 |
| LOEUF | GRCh37 | 1,650.07 | 186.11 | 1,487.94 | 202.81 |
| LOEUF | GRCh38 | 4,365.30 | 91.77 | 4,068.43 | 99.78 |
| LoFTEE | GRCh38 | 10,141.20 | 1,085.49 | 7,083.69 | 1,210.78 |
| LoFtool | GRCh37 | 1,433.72 | 8.49 | 1,016.26 | 10.66 |
| LoFtool | GRCh38 | 3,638.60 | 29.83 | 3,595.54 | 36.41 |
| NearestExonJB | GRCh37 | 1,497.79 | 11.39 | 1,033.45 | 15.41 |
| NearestExonJB | GRCh38 | 3,705.35 | 41.71 | 3,613.16 | 59.91 |
| NMD | GRCh37 | 1,457.15 | 4.40 | 1,057.35 | 7.42 |
| NMD | GRCh38 | 4,000.75 | 15.12 | 3,596.98 | 23.18 |
| OpenTargets | GRCh38 | 4,225.77 | 42.86 | 4,057.12 | 51.80 |
| pLI | GRCh37 | 1,493.46 | 9.13 | 1,069.89 | 13.17 |
| pLI | GRCh38 | 4,037.04 | 31.65 | 3,570.55 | 39.05 |
| REVEL | GRCh37 | 1,535.44 | 44.23 | 1,191.70 | 55.63 |
| REVEL | GRCh38 | 4,108.96 | 50.19 | 3,805.47 | 63.33 |
| TSSDistance | GRCh37 | 1,472.09 | 4.10 | 1,047.73 | 6.86 |
| TSSDistance | GRCh38 | 3,960.81 | 12.93 | 3,524.43 | 19.05 |
| UTRAnnotator | GRCh37 | 1,714.13 | 4.26 | 1,248.34 | 6.47 |
| UTRAnnotator | GRCh38 | 3,722.67 | 13.38 | 3,908.05 | 18.97 |

**Method.** For a row and architecture both engines ran on one machine, `c9gd.8xlarge` (ARM Graviton5) or `c8id.8xlarge` (x86 Intel), 32 vCPU with local NVMe, each with `--fork 16` and `--buffer_size 5000`; the timed span is the engine's process alone. Ensembl VEP's figure is one timed run. vep-rs's is the median of 5 timed runs, after one discarded full run.

Ensembl VEP ran with `--stats_text`, as its concordance reference runs do, where the ClinVar rows of [Whole-genome wall time](#whole-genome-wall-time) time it with `--no_stats`. Ensembl VEP ran LoFTEE's module with stand-ins for the MaxEntScan files it bundles, on its image with samtools added. Every timed run is in [`docs/concordance-provenance/plugin-walltime-v1.0.0.json`](../concordance-provenance/plugin-walltime-v1.0.0.json).

<!-- plugin-wall-time: end -->
