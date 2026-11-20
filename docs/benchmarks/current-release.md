# Benchmarks against Ensembl VEP and fastVEP

vep-rs benchmarks for the current release: wall time against Ensembl VEP and fastVEP and concordance with Ensembl VEP, measured on whole-genome datasets.

<!-- readme-section: Why vep-rs | fragments: benchmarks/current-release.md -->

<!-- readme-section: Current release -->

The figures the paper reports, measured at its pinned build, are on [the paper's figures page](../published-figures.md) and stay as published.

Each release page on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page states each dataset's peak RSS per architecture beside its wall time. That figure is the largest per-chromosome median of the release measurement. Where the release measured it, the page also states the peak of a whole genome annotated as one input file. [Installation](../getting-started/installation.md#memory) states what those figures cover and how a run's memory grows with the chromosomes an input touches.

Two further records under `docs/concordance-provenance/`, their names ending `-no-extra-args.json`, hold fastVEP v0.5.0 run without `--exclude-readthrough --exclude-artifact`: one its ClinVar and structural-variant concordance, the other its ClinVar wall times. No table above reads them.

<!-- plugin-wall-time: begin -->

## Wall time with each plugin

Seconds for one annotation run of the genome-wide ClinVar dataset of each assembly: 4,388,172 variants on GRCh37 and 4,398,127 variants on GRCh38. A plugin's row is the run of the first row, `No plugin`, with that one plugin added. HGVSIntronOffset, which reads `--hgvs`, is added to the run under `--hgvs`. The 20 plugins are those the [Plugins](#plugins) tables score on ClinVar: a plugin whose data cannot be redistributed is not timed, and neither is TranscriptAnnotator, which reads a table of the user's own.

| Plugin | Assembly | Ensembl VEP 116.2, ARM (s) | vep-rs v0.6.0, ARM (s) | Ensembl VEP 116.2, x86 (s) | vep-rs v0.6.0, x86 (s) |
| --- | --- | --- | --- | --- | --- |
| No plugin | GRCh37 | 1,047.51 | 1.96 | 1,034.36 | 2.96 |
| No plugin | GRCh38 | 4,027.45 | 7.51 | 3,622.69 | 9.92 |
| No plugin, `--hgvs` | GRCh37 | 1,161.92 | 6.16 | 1,271.86 | 9.37 |
| No plugin, `--hgvs` | GRCh38 | 4,127.65 | 23.48 | 4,741.23 | 35.71 |
| AlphaMissense | GRCh37 | 1,086.38 | 41.46 | 1,197.08 | 53.89 |
| AlphaMissense | GRCh38 | 4,148.27 | 47.59 | 3,884.30 | 61.66 |
| AncestralAllele | GRCh37 | 1,198.21 | 6.16 | 1,122.80 | 8.54 |
| AncestralAllele | GRCh38 | 3,708.30 | 14.24 | 3,772.67 | 20.36 |
| Blosum62 | GRCh37 | 1,108.46 | 5.21 | 1,111.93 | 7.47 |
| Blosum62 | GRCh38 | 4,050.02 | 18.12 | 3,636.90 | 25.33 |
| Conservation | GRCh38 | 5,698.72 | 15.76 | 6,305.36 | 22.09 |
| DosageSensitivity | GRCh37 | 1,094.55 | 8.78 | 1,125.38 | 11.40 |
| DosageSensitivity | GRCh38 | 4,020.92 | 30.80 | 3,646.24 | 43.32 |
| Downstream | GRCh37 | 1,012.72 | 21.08 | 1,110.54 | 50.95 |
| Downstream | GRCh38 | 3,906.16 | 68.97 | 3,823.95 | 198.47 |
| gnomADc | GRCh37 | 2,276.60 | 222.55 | 2,201.99 | 241.08 |
| gnomADc | GRCh38 | 4,688.07 | 82.75 | 4,427.76 | 99.36 |
| gnomADMt | GRCh38 | 3,665.99 | 10.08 | 3,766.90 | 15.07 |
| GWAS | GRCh38 | 4,063.44 | 12.08 | 3,764.84 | 16.70 |
| HGVSIntronOffset | GRCh37 | 1,700.94 | 13.26 | 1,381.82 | 20.38 |
| HGVSIntronOffset | GRCh38 | 4,800.54 | 51.62 | 4,976.86 | 75.96 |
| LOEUF | GRCh37 | 1,782.91 | 192.24 | 1,548.06 | 204.83 |
| LOEUF | GRCh38 | 4,373.12 | 93.77 | 4,188.71 | 100.42 |
| LoFTEE | GRCh38 | 10,969.10 | 1,087.50 | 7,252.10 | 1,265.82 |
| LoFtool | GRCh37 | 1,610.80 | 8.91 | 1,077.17 | 10.38 |
| LoFtool | GRCh38 | 4,028.74 | 31.40 | 3,620.39 | 36.58 |
| NearestExonJB | GRCh37 | 1,516.26 | 11.37 | 1,104.09 | 15.25 |
| NearestExonJB | GRCh38 | 4,093.40 | 42.64 | 3,644.92 | 59.50 |
| NMD | GRCh37 | 1,497.12 | 4.54 | 1,062.37 | 7.69 |
| NMD | GRCh38 | 3,987.93 | 16.15 | 3,682.05 | 24.63 |
| OpenTargets | GRCh38 | 4,303.80 | 44.43 | 3,940.35 | 52.15 |
| pLI | GRCh37 | 1,633.74 | 9.20 | 1,074.19 | 12.62 |
| pLI | GRCh38 | 3,067.19 | 31.39 | 3,582.17 | 40.28 |
| REVEL | GRCh37 | 1,575.83 | 45.64 | 1,167.42 | 55.53 |
| REVEL | GRCh38 | 3,389.23 | 51.20 | 3,947.92 | 65.35 |
| TSSDistance | GRCh37 | 1,089.83 | 4.11 | 1,025.89 | 7.01 |
| TSSDistance | GRCh38 | 3,179.11 | 13.27 | 3,609.68 | 20.44 |
| UTRAnnotator | GRCh37 | 1,410.60 | 4.25 | 1,267.13 | 6.49 |
| UTRAnnotator | GRCh38 | 3,725.09 | 13.62 | 3,924.28 | 18.76 |

**Method.** For a row and architecture both engines ran on one machine, `c9gd.8xlarge` (ARM Graviton5) or `c8id.8xlarge` (x86 Intel), 32 vCPU with local NVMe, each with `--fork 16` and `--buffer_size 5000`; the timed span is the engine's process alone. Ensembl VEP's figure is one timed run. vep-rs's is the median of 5 timed runs on the rows without a plugin and on LoFTEE, and one timed run on every other row, after one discarded full run.

Ensembl VEP ran with `--stats_text`, as its concordance reference runs do, where the ClinVar rows of [Whole-genome wall time](#whole-genome-wall-time) time it with `--no_stats`. Ensembl VEP ran LoFTEE's module with stand-ins for the MaxEntScan files it bundles, on its image with samtools added. Every timed run is in [`docs/concordance-provenance/plugin-walltime-v0.6.0.json`](../concordance-provenance/plugin-walltime-v0.6.0.json).

<!-- plugin-wall-time: end -->
