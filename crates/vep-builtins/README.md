# vep-builtins

Built-in VEP plugin implementations compiled directly into the `vep` binary.

## Architecture

Instead of compiling each plugin as a separate `cdylib` shared library (the C FFI
approach in `vep-plugin`), built-in plugins operate on native Rust types with no
JSON serialization or FFI overhead.

```
vep-cli (runner.rs)
  |
  +-- vep-builtins (this crate)
  |     |
  |     +-- BuiltinPlugin trait
  |     |     +-- prefetch() + annotate()  <-- two-phase: parallel I/O, then sequential writes
  |     |     +-- run_batch()              <-- single-phase fallback, once per buffer
  |     |     +-- run()                    <-- per-transcript-consequence fallback
  |     |
  |     +-- AnnotationStore trait (one query interface over two backends)
  |     |     +-- TabixAnnotator: noodles::tabix + CSI binning index over `.tsv.gz`
  |     |     +-- BinaryAnnotator: mmap-backed `.vpd` data + `.vpdi` index, chosen
  |     |     |   when a `.vpd` sits beside the `.tsv.gz` (`vep-cache-converter
  |     |     |   convert-plugin-data` writes it)
  |     |     +-- Batch region merging (1000 bp gap)
  |     |     +-- Allele matching & normalization
  |     |
  |     +-- BuiltinRegistry
  |           +-- Name-based plugin resolution
  |           +-- Returns the names no built-in claims to the caller
  |
  +-- vep-plugin (C FFI dylib loader, for the names no built-in claims)
```

## Plugins

Support is per reference genome, as in the root README's plugin matrix. The
Mechanism column names how the plugin reads its data and matches a variant.

| Plugin        | Mechanism                                                                                      | GRCh37 | GRCh38   | Output fields                                                                                                                                 |
| ------------- | ---------------------------------------------------------------------------------------------- | ------ | -------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| CADD          | Annotation store, position + allele match                                                      | Yes    | Yes      | `CADD_PHRED`, `CADD_RAW`                                                                                                                      |
| REVEL         | Annotation store, position + allele match, transcript match unless `no_match=1`; missense only | Yes    | Yes      | `REVEL_score`                                                                                                                                 |
| gnomADc       | Annotation store, position only                                                                | Yes    | Yes      | `gnomADc_*` (one field per coverage column in the file header)                                                                                |
| AlphaMissense | Annotation store, position + allele match                                                      | Yes    | Yes      | `am_class`, `am_pathogenicity`                                                                                                                |
| dbscSNV       | Annotation store, position + allele match, per transcript                                      | Yes    | Yes      | `dbscSNV_ADA_SCORE`, `dbscSNV_RF_SCORE`                                                                                                       |
| LoFTEE        | Consequence filter over the transcript model, `--fasta` splice motifs, GERP                    | Yes    | Yes      | `LoF`, `LoF_filter`, `LoF_flags`, `LoF_info`                                                                                                  |
| LoFtool       | Flat file loaded at init, gene symbol lookup                                                   | Yes    | Yes      | `LoFtool_percentile`                                                                                                                          |
| pLI           | Flat file loaded at init, gene symbol lookup                                                   | Yes    | Yes      | `pLI_gene_value`                                                                                                                              |
| GWAS          | Flat file loaded at init, position + risk allele match                                         | No     | Yes      | `GWAS_accessions`, `GWAS_associated_gene`, `GWAS_beta_coef`, `GWAS_odds_ratio`, `GWAS_p_value`, `GWAS_pmid`, `GWAS_risk_allele`, `GWAS_study` |
| SpliceAI      | Annotation store, position + allele + gene symbol match                                        | No     | SNV only | `SpliceAI_pred_*` (9 fields)                                                                                                                  |
| dbNSFP        | Annotation store, position + allele + transcript match, per transcript                         | No     | No       | User-selected columns; not supported (its academic license excludes commercial use)                                                          |

LoFtool and pLI are keyed by gene symbol, so one data file serves both assemblies.
LoFTEE reads GERP from a tabix TSV on GRCh37 (`gerp_tabix=`) and from a bigWig on
GRCh38 (`gerp_bigwig=`), the two forms its upstream bundles ship.

## Batch Query Strategy

`TabixAnnotator` issues one tabix query per merged region rather than one per
variant:

```
per variant:
  tabix_seek(chr, pos)            # one seek per variant

per buffer (TabixAnnotator):
  sort the buffer's variants by position
  merge overlapping regions (bridging gaps of up to 1000 bp)
  tabix_seek per merged region    # one seek per merged region
  index the records by position
  for each variant:
    hash lookup                   # no I/O
```

`BinaryAnnotator` answers the same merged-region queries by binary search over a
memory-mapped `.vpd` file, with no decompression.

## Usage

Plugins are loaded via the standard `--plugin` CLI flag:

```bash
vep -i input.vcf \
    --assembly GRCh38 \
    --plugin CADD,snv=/path/to/cadd_snv.tsv.gz,indels=/path/to/cadd_indels.tsv.gz \
    --plugin gnomADc,/path/to/gnomad_coverage.tsv.gz \
    --plugin REVEL,file=/path/to/revel.tsv.gz \
    --offline --cache
```

Pass `--assembly` whenever plugins are enabled. CADD, dbscSNV and AlphaMissense use it
to reject a data file whose own build marker contradicts the run; REVEL uses it to
pick which of its two coordinate columns (`hg19_pos`, `grch38_pos`) the index was
built on; LoFTEE uses it to select the END_TRUNC cutoff and the GERP reader. Without
it those guards are inert and REVEL falls back to the GRCh37 column with only a
warning. gnomADc carries no in-file build marker and is not checked.

The runner resolves `--plugin Name,...` arguments against the built-in registry
first. A name no built-in claims is loaded as a C-ABI dynamic library through
`vep-plugin`.
