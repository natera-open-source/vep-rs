# Compatibility with Ensembl VEP

This page gathers the Ensembl VEP behaviour vep-rs does not implement. The flags the group pages list under "Accepted for compatibility" are the whole of it; this page adds the three refused flags, the consequence terms never assigned, the input formats not read and the absent database mode.

## Consequence terms never assigned

Six regulatory-feature terms need a regulatory-feature cache that is not loaded: `TFBS_ablation`, `TFBS_amplification`, `TF_binding_site_variant`, `regulatory_region_ablation`, `regulatory_region_amplification` and `regulatory_region_variant`. `--regulatory`, alone or through `--everything`, is accepted and warns that regulatory-feature annotation is not implemented.

The root term `sequence_variant` has no predicate, and Ensembl VEP's consequence engine never assigns it either. The other 34 of Ensembl VEP's 41 Sequence Ontology terms are assigned.

## Input formats

`--format` takes `vcf`, `ensembl` (the tab-delimited default format), `hgvs` and `region` (`chr:start-end/allele` and `chr:start-end:strand/allele`), or `guess`. HGVS input is genomic (`:g.`) only; `:c.`, `:n.`, `:p.` and `:r.` notation is rejected as not implemented.

The `id`, `spdi` and `caid` formats are not read: Ensembl VEP resolves all three through its variation database and refuses them in `--offline` mode.

## Refused flags

Three Ensembl VEP 116 flags are declared so that the command line names them, and refused with an error that says what each would do. `--regulatory_gff` and `--extended_promoters` belong to regulatory annotation, `--custom_suppress_filter` to custom annotation; neither is supported.

| Flag | What Ensembl VEP does with it |
| --- | --- |
| [`--regulatory_gff`](#opt_regulatory_gff) | Annotates regulatory features and motifs from a GFF3 file in place of the cache |
| [`--extended_promoters`](#opt_extended_promoters) | Widens the promoters read with `--regulatory_gff` to their extended bounds |
| [`--custom_suppress_filter`](#opt_custom_suppress_filter) | Keeps a `--custom` VCF source's `FILTER` column out of the output unless the source's `fields` argument names it |

### `--regulatory_gff` {#opt_regulatory_gff}

**`--regulatory_gff <file>`**, string, no default. Ensembl VEP: refused.

Refused with an error: it would annotate regulatory features and motifs from a GFF3 file in place of the cache, and regulatory annotation is not supported.

### `--extended_promoters` {#opt_extended_promoters}

**`--extended_promoters`**, boolean, default off. Ensembl VEP: refused.

Refused with an error: it would widen the promoters read with `--regulatory_gff` to their extended bounds, and regulatory annotation is not supported.

### `--custom_suppress_filter` {#opt_custom_suppress_filter}

**`--custom_suppress_filter`**, boolean, default off. Ensembl VEP: refused.

Refused with an error: it would keep a `--custom` VCF source's `FILTER` column out of the output unless the `--custom` source's `fields` argument names it, and custom annotation is not supported.

## No database mode

vep-rs is offline-first and never contacts a database. `--database`, `--host`, `--port`, `--user` and `--password` do not exist and are rejected as unknown flags. The transcript set is the `--json_cache` contents ([Cache setup](../cache-setup.md)), so `--refseq`, `--merged` and `--gencode_basic` select nothing.

## Accepted without effect

Every flag below is parsed and changes no annotation. The group page of each carries its entry.

| Flag | Why it has no effect | Group page |
| --- | --- | --- |
| `--input_data` | Input is read from `--input_file` | [Input and output](input-output.md#opt_input_data) |
| `--compress_output` | The output is written uncompressed | [Input and output](input-output.md#opt_compress_output) |
| `--fields` | The columns follow the field flags | [Input and output](input-output.md#opt_fields) |
| `--cache`, `--offline` | The annotation cache is `--json_cache`; vep-rs always runs offline | [Cache, species and runtime](cache-species-runtime.md#opt_cache) |
| `--dir`, `--dir_cache`, `--cache_version` | They locate no cache; they name the cache path the output header prints | [Cache, species and runtime](cache-species-runtime.md#opt_dir) |
| `--species` | Names the species in the header's cache path; the data annotated is the `--json_cache` contents | [Cache, species and runtime](cache-species-runtime.md#opt_species) |
| `--hgvsg` | The `HGVSg` column is added and always empty | [Output control](output-control.md#opt_hgvsg) |
| `--regulatory` | Regulatory-feature annotation is not implemented; the flag warns | [Output control](output-control.md#opt_regulatory) |
| `--coding_only`, `--no_intergenic` | Not implemented | [Output control](output-control.md#opt_coding_only) |
| `--pick_order` | The pick order is fixed to Ensembl VEP's default | [Filtering, picking and transcript sets](filtering-picking.md#opt_pick_order) |
| `--total_length` | Not implemented | [Filtering, picking and transcript sets](filtering-picking.md#opt_total_length) |
| `--refseq`, `--merged`, `--gencode_basic` | The transcript set is the `--json_cache` contents | [Filtering, picking and transcript sets](filtering-picking.md#opt_refseq) |
| `--custom` | Custom annotation sources are not loaded | [Plugins](plugins.md#opt_custom) |
| `--no_stats`, `--stats_file`, `--stats_text` | No statistics file is written, with or without them | [Statistics, verbosity and miscellaneous](statistics-verbosity-misc.md#opt_no_stats) |
| `--mirna` | The `miRNA` column is added and always empty | [Statistics, verbosity and miscellaneous](statistics-verbosity-misc.md#opt_mirna) |

Protein domains are never mapped. The `DOMAINS` field exists and protein features are read from the cache, but `--domains` leaves the column empty, and the native cache builder writes no domain features: the Ensembl GFF3 it reads carries none.
