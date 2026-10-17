# Consequence filter: LoFTEE

LoFTEE is not a score lookup. It filters the consequences vep-rs has already
assigned against the transcript model, the `--fasta` splice motifs and GERP, and
emits a loss-of-function confidence call. This page lists its parameters, the
two upstream data bundles and what the port leaves out.

LoFTEE targets parity with the upstream LoFTEE `grch38` branch defaults, with
the splice-prediction layer off. On GRCh37 the matching Perl run is the
`master` branch with `run_splice_predictions:0`, because `master` turns that
layer on by default. Built-in names are matched exactly as written: Perl VEP
loads LoFTEE as `--plugin LoF,...` and vep-rs as `--plugin LoFTEE,...`.

## Parameters and data {#loftee}

| Property      | Value                                                                                                                                                                  |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                                                                         |
| Matches       | consequence filter over the transcript model, `--fasta` splice motifs, GERP                                                                                            |
| Output fields | `LoF`, `LoF_filter`, `LoF_flags`, `LoF_info`                                                                                                                           |
| Data source   | [LoFTEE data](https://personal.broadinstitute.org/konradk/loftee_data/), per assembly: `GRCh38/` (GERP bigWig) or `GRCh37/` (GERP tabix TSV), plus human_ancestor.fa |

Every parameter is a named `key=value` and optional; LoFTEE's native `key:value`
form is accepted alongside.

| Parameter                | Meaning                                                                                                                                                                                                      |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `human_ancestor_fa`      | Plain (uncompressed) indexed FASTA (`.fa` + `.fai`). Enables the `ANC_ALLELE` filter (SNP-only). Decompress the upstream `human_ancestor.fa.gz` first; vep-rs does not read bgzip here.                       |
| `gerp_bigwig`            | GERP conservation bigWig, as shipped in the upstream `grch38` bundle. Enables the default GERP-weighted `END_TRUNC` filter.                                                                                   |
| `gerp_tabix`             | GERP as a tabix-indexed per-base TSV, as shipped in the upstream `master`/GRCh37 bundle (`GERP_scores.final.sorted.txt.gz` + `.tbi`). Same filter, different on-disk format. Alias `gerp_scores`.           |
| `gerp`                   | Assembly-agnostic spelling of the GERP source; classifies the file by its extension.                                                                                                                         |
| `use_gerp_end_trunc`     | Default `true`. With no GERP source, `use_gerp_end_trunc=false` selects the 5%-percentile `END_TRUNC` fallback, which needs no data file.                                                                      |
| `gerp_end_trunc_cutoff`  | Derived from `--assembly` unless set, because the two upstream branches use cutoffs of opposite sign: `-58` on `grch38`, `+180` on `master`/GRCh37. Setting it pins that value for both.                      |
| `min_intron_size`        | Default `15`. Tunable matching the Perl plugin.                                                                                                                                                              |
| `filter_position`        | Default `0.05`. Tunable matching the Perl plugin.                                                                                                                                                            |
| `check_complete_cds`     | Default `false`. Tunable matching the Perl plugin.                                                                                                                                                           |
| `run_splice_predictions` | Accepted so the key does not trip the unknown-parameter warning when a Perl `LoF` argument string is reused. `0` or an empty value is accepted silently; any other value (true by Perl's string truth) logs a warning, and the run continues without the layer. |
| `conservation_file`      | Accepted for argument parity; a no-op.                                                                                                                                                                       |

Supply exactly one of `gerp_bigwig` and `gerp_tabix`. An unreadable GERP path is
a hard error at plugin init rather than a silent degradation. The two upstream
bundles are structurally different, which is why the GERP parameter differs per
assembly; vep-rs reads both formats natively.

`--fasta`, the pipeline reference genome, is required for the splice-motif
checks: the `GC_TO_GT_DONOR` filter and the `NON_CAN_SPLICE` and `NAGNAG_SITE`
flags. Without it those checks are skipped.

`$LOFTEE` below is the unpacked bundle for the run's assembly,
`/path/to/loftee/grch38` or `/path/to/loftee/grch37`:

```bash
# GRCh38: upstream ships GERP as a bigWig
--plugin LoFTEE,human_ancestor_fa=$LOFTEE/human_ancestor.fa,\
gerp_bigwig=$LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw

# GRCh37: upstream ships GERP as a tabix per-base TSV instead
--plugin LoFTEE,human_ancestor_fa=$LOFTEE/human_ancestor.fa,\
gerp_tabix=$LOFTEE/GERP_scores.final.sorted.txt.gz
```

**Preparing the data.** Download the bundle for the run's assembly from
[LoFTEE data](https://personal.broadinstitute.org/konradk/loftee_data/). vep-rs
reads the plain `human_ancestor.fa`, so decompress the upstream bgzipped file
and re-index it (`bgzip -d`, then `samtools faidx`). The GERP file is used as
shipped: the bigWig on GRCh38, the tabix-indexed TSV with its `.tbi` on GRCh37.

## What LoFTEE annotates

LoFTEE only annotates `stop_gained`, `frameshift_variant`,
`splice_acceptor_variant` and `splice_donor_variant` on protein-coding
transcripts; everything else gets no `LoF` field.

The MaxEntScan + logistic-regression + de-novo-donor-SVM splice-prediction layer
and the `OS` confidence class, off by default on the `grch38` branch, are not
implemented. The PhyloCSF `LoF_flags` (`PHYLOCSF_WEAK`/`PHYLOCSF_UNLIKELY_ORF`)
are a documented gap: they are flags, not filters, so they never change the
HC/LC call, and they would require a SQLite dependency.
