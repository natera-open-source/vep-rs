# Regulatory and motif features

`--regulatory`, which `--everything` switches on (`Config.pm` line 362), annotates each record against the regulatory features and transcription-factor binding motifs of the cache, as Ensembl VEP does from the `_reg.gz` files beside its transcript files (`AnnotationSource/Cache/RegFeat.pm`). The rows, terms and fields follow Ensembl VEP's `OutputFactory.pm` lines 1867 to 1942; the two fields vep-rs does not write and the one `--spdi` difference are on [the divergence page](divergences/regulatory-rows.md).

## The cache

The features come from an Ensembl VEP cache converted with `scripts/data/storable_to_json.pl` ([Cache setup](cache-setup.md#convert-an-existing-ensembl-vep-cache)). Beside `transcripts/<chr>/` the script writes `regulatory/<chr>/<start>-<end>.json`, one shard per region holding the region's regulatory features and motif features, and sets `regulatory: true` in `info.json`; `vep-cache-converter` copies the shards, as `.json` or `.json.gz`, and writes the marker only after the tree has been written. The Ensembl VEP caches that carry a regulatory build convert with it; the human caches do, on both assemblies.

A cache converted without the regulatory tree holds no regulatory features, and a cache built by `vep-cache-builder` holds none: convert the Ensembl VEP cache again. Every motif feature carries its binding matrix; the converter and the loader refuse a motif without one as a corrupt shard. The shards of a chromosome are read the first time a record of that chromosome is annotated, so a run without `--regulatory` never opens them.

On a cache with no regulatory data, `--regulatory` writes one warnings-file line: `--regulatory: this cache holds no regulatory data; no RegulatoryFeature or MotifFeature row is written`. The regulatory columns of the tab and VCF headers are declared and print empty. Ensembl VEP skips the source without a word (`CacheDir.pm` lines 168 to 176).

## The rows

A record inside a regulatory feature gets a `RegulatoryFeature` row with `regulatory_region_variant`; a record inside a motif gets a `MotifFeature` row with `TF_binding_site_variant`. A variant that covers the whole feature adds `regulatory_region_ablation` or `TFBS_ablation` when it is a deletion, or `regulatory_region_amplification` or `TFBS_amplification` when it is a duplication, a tandem duplication or a tandem repeat, a sequence variant whose allele repeats its reference included. Five of the six terms are `MODIFIER`; `TFBS_ablation` is `MODERATE`. An insertion overlaps a feature when both of its flanking bases lie inside it.

A record that overlaps no transcript keeps its `intergenic_variant` row beside its regulatory and motif rows, with `-` in the Gene, Feature and Feature_type cells. The rows of a record come in Ensembl VEP's order: transcript rows, regulatory-feature rows, motif-feature rows, then the intergenic row. In JSON the rows are the `regulatory_feature_consequences`, `motif_feature_consequences` and `intergenic_consequences` arrays, and `most_severe_consequence` ranks every row's terms.

`--pick`, `--pick_allele`, `--flag_pick` and `--flag_pick_allele` rank regulatory and motif rows with the transcript rows, and `--most_severe` and `--summary` read their terms. `--per_gene` and `--pick_allele_gene` keep every regulatory, motif and intergenic row beside the one transcript row they keep per gene (`OutputFactory.pm` lines 842 to 844). `--no_intergenic` drops the rows of a record with no transcript row, and `--coding_only` drops every regulatory and motif row.

## The fields

A `RegulatoryFeature` row carries `BIOTYPE` (the feature's type, such as `promoter`, `enhancer` or `CTCF_binding_site`) whether or not `--biotype` is given, and no `STRAND`. A `MotifFeature` row carries `STRAND`, `MOTIF_NAME` (the binding matrix's stable id), `MOTIF_POS` (the position of the variant's first base in the matrix, counted along the motif's own strand) and `TRANSCRIPTION_FACTORS` (the matrix's factors, comma-separated; empty when the matrix names none). `HIGH_INF_POS` and `MOTIF_SCORE_CHANGE` are empty on every row. The fields a transcript row alone carries print `-`, and the record-level fields print on every row.

A structural variant's regulatory and motif rows carry none of `BIOTYPE`, `STRAND` and the motif fields, and carry `OverlapBP` and `OverlapPC` as a transcript row does ([VEP default and tab](output/vep-and-tab.md)).

A plugin annotates the rows of the feature types it declares, transcript rows unless it says otherwise ([Writing a plugin](developers/writing-a-plugin.md)); a plugin that declares regulatory or motif features switches `--regulatory` on. `CELL_TYPE` is not written, and `--cell_type` is refused as an unknown argument. `--regulatory_gff` and `--extended_promoters` are refused with an error: regulatory features are read from the cache alone.
