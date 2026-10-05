# The figures the paper reports

These are the concordance and wall-time tables of Porter and Borkowski, _vep-rs: high-throughput Rust variant annotation with population-scale concordance to Ensembl VEP_ (bioRxiv 2026, [doi:10.64898/2026.09.22.753614](https://doi.org/10.64898/2026.09.22.753614)), measured at the paper's pinned build of vep-rs, fastVEP v0.3.0 and Ensembl VEP release 115.2 on the chromosome 21 suites and the two ClinVar releases. They stay as published; the current release's own figures are in the [README](../README.md#current-release-v031).

Perl VEP release 115.2 is the reference engine in every table; fastVEP is the other Rust VEP port. Every figure here re-derives from the released measurement data in [`manuscript/data/`](../manuscript/data/README.md). Concordance is measured on the complete consequence output of each dataset, tuple by tuple on matched caches, and F1 is architecture-independent, so one value per dataset covers both ARM and x86. **Adjusted F1** sets aside the tuples of documented divergence classes: each is a shape of disagreement traced to a mechanism in one engine's source and reproduced on a public record with a command a reader can run. [docs/intended-divergences.md](intended-divergences.md) documents every class the comparators in this repository set aside: eleven classes, ten of them Perl VEP defects vep-rs does not reproduce and one (`<CNV:TR>`) a difference of representation in which neither engine is wrong. The paper's figures below set aside the five classes its supplement documents (S2, S4.3 and Table S5), four of them Perl VEP defects; on the six SNP/indel datasets that mask removes matched pairs from both engines' sides in equal numbers, the covered `splice_region_variant` swap (all 1,884 pairs on the two ClinVar datasets) and the `start_lost`/`start_retained_variant` co-emission, and the structural-variant results per input file are in the supplement, every value re-deriving from [`manuscript/data/`](../manuscript/data/README.md).

| Dataset             | Assembly | vep-rs Raw F1 | vep-rs Adj F1 | fastVEP Raw F1 |
| ------------------- | -------- | ------------- | ------------- | -------------- |
| ClinVar full        | GRCh37   | 0.999979      | 1.000000      | 0.825597       |
| ClinVar full        | GRCh38   | 0.999974      | 1.000000      | 0.969701       |
| gnomAD v2.1.1 chr21 | GRCh37   | 0.999999      | 1.000000      | 0.762022       |
| gnomAD v4.1 chr21   | GRCh38   | 0.999999      | 1.000000      | 0.919675       |
| 1KG Phase 3 chr21   | GRCh37   | 1.000000      | 1.000000      | 0.790928       |
| 1KG high-cov chr21  | GRCh38   | 1.000000      | 1.000000      | 0.938604       |
| SV per-VCF (16 files) | GRCh37 | 0.975395      | 0.998524      | 0.153626       |
| SV per-VCF (16 files) | GRCh38 | 0.909998      | 0.998031      | 0.133436       |

Wall time is the median of 20 independent machine measurements per cell, each on a fresh instance after a discarded warmup, on sites-only inputs; vep-rs and Perl VEP run `--fork 16`, fastVEP parallelizes internally at 16 threads. vep-rs annotates full ClinVar GRCh37 (4,388,172 input records) in 6.00 s on ARM and 6.20 s on x86, and its peak resident memory spans 0.62 to 4.56 GiB across every measured cell.

| Arch | Dataset             | Perl (s) | fastVEP (s) | vep-rs (s) | vep-rs vs Perl | vep-rs vs fastVEP |
| ---- | ------------------- | -------- | ----------- | ---------- | -------------- | ----------------- |
| ARM  | ClinVar GRCh37      | 961.08   | 43.65       | 6.00       | 160×           | 7.28×             |
| ARM  | ClinVar GRCh38      | 2,722.01 | 121.13      | 12.35      | 220×           | 9.81×             |
| ARM  | gnomAD v2.1.1 chr21 | 1,500.91 | 68.24       | 6.84       | 219×           | 9.98×             |
| ARM  | gnomAD v4.1 chr21   | 3,953.00 | 193.44      | 13.92      | 284×           | 13.9×             |
| ARM  | 1KG Phase 3 chr21   | 106.29   | 6.26        | 1.05       | 101×           | 5.96×             |
| ARM  | 1KG high-cov chr21  | 249.15   | 17.45       | 1.85       | 135×           | 9.43×             |
| x86  | ClinVar GRCh37      | 849.97   | 33.16       | 6.20       | 137×           | 5.35×             |
| x86  | ClinVar GRCh38      | 2,378.35 | 90.26       | 12.56      | 189×           | 7.19×             |
| x86  | gnomAD v2.1.1 chr21 | 1,341.62 | 57.69       | 10.16      | 132×           | 5.68×             |
| x86  | gnomAD v4.1 chr21   | 3,276.52 | 158.24      | 17.06      | 192×           | 9.28×             |
| x86  | 1KG Phase 3 chr21   | 85.31    | 5.10        | 1.09       | 78.3×          | 4.68×             |
| x86  | 1KG high-cov chr21  | 205.16   | 14.60       | 1.75       | 118×           | 8.36×             |

Across those cells vep-rs is faster than Perl VEP by a geometric mean of 176× (ARM) / 135× (x86) and faster than fastVEP by 9.07× (ARM) / 6.55× (x86). vep-rs was timed in one campaign and the comparators' SNP and indel cells in another, on the same instance types against the same dataset inventory, so every ratio in the table divides medians from separate machines and campaigns; the structural-variant cells of all three engines come from the vep-rs campaign. The structural-variant sets are annotated as 16 separate per-file invocations, so per-invocation start-up weighs far more there; vep-rs is faster than Perl VEP on those cells by 52.9× (GRCh37) and 23.1× (GRCh38) on ARM and 37.5× and 22.5× on x86 (0.61 s, 2.81 s, 0.79 s and 2.63 s against Perl's 32.28 s, 64.92 s, 29.62 s and 59.29 s). No fastVEP structural-variant speedup is given, because fastVEP emits 62% and 88% of Perl VEP's tuple volume on those sets and recovers only 12.5% and 12.7% of Perl VEP's tuples, so its wall time does not buy comparable annotations. Plugin output and HGVS notation are outside every F1 above.
