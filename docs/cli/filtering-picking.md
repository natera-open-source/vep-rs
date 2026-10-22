# Filtering, picking and transcript sets

These flags reduce the rows reported per record, allele or gene, flag the ones a pick would drop, shape the allele columns and filter transcripts by their attributes. The transcript-set flags select the cache read and the GENCODE set annotated.

## Filtering and picking

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--pick` | bool | false | Pick one consequence per variant |
| `--pick_allele` | bool | false | Pick one consequence per variant allele |
| `--per_gene` | bool | false | Pick one consequence per gene |
| `--pick_allele_gene` | bool | false | Pick one consequence per allele per gene |
| `--most_severe` | bool | false | Report only the single most severe consequence |
| `--summary` | bool | false | Report a summary of all consequences |
| `--flag_pick` | bool | false | Flag (but don't remove) non-picked consequences |
| `--flag_pick_allele` | bool | false | Flag non-picked consequences per allele |
| `--flag_pick_allele_gene` | bool | false | Flag non-picked consequences per allele per gene |
| `--pick_order` | `mane_select`, `mane_plus_clinical`, `canonical`, `appris`, `tsl`, `biotype`, `ccds`, `rank`, `length`, `ensembl`, `refseq` (comma-separated) | none | The criteria the pick family compares, in order, as a comma-separated list; a shorter list compares its criteria and then keeps the first transcript in stable-id order, and an empty value keeps it at once. Default: mane\_select,mane\_plus\_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq |
| `--transcript_filter` | String (repeatable) | none | Annotate only the transcripts that pass a filter, in Ensembl VEP's filter grammar: units of `field operator value` joined by `and` and `or`, grouped by parentheses and inverted by `not`; repeated flags are joined by `and`. The fields are the cache's transcript keys, named exactly, in another case, without the leading `_`, or by a unique prefix: `stable_id`, `version`, `biotype`, `start`, `end`, `strand`, `is_canonical`, `source`, `dbID`, `cdna_coding_start`, `cdna_coding_end`, `_gene_stable_id`, `_gene_version`, `_gene_symbol`, `_gene_symbol_source`, `_gene_hgnc_id`, `_ccds`, `_refseq`, `_protein`, `_swissprot`, `_trembl`, `_uniparc`, `_uniprot_isoform`, `_gene_phenotype`; the object keys `translation` (coding transcripts), `slice`, `attributes`, `_gene`, `_trans_exon_array`, `_variation_effect_feature_cache` and `_vep_lazy_loaded`, which are defined and read meaningfully by `exists` and `nex` alone; and `description`, `coding_region_start` and `coding_region_end`, which the cache leaves undefined. `_source_cache` is absent, as it is from every cache but a merged one. The operators are `is` (`eq`, `=`), `ne` (`!=`), `>`, `<`, `>=`, `<=`, `exists` (`defined`, the default), `nex`, `match` (`matches`, `regex`, `re`), `nre`, `in` (a comma-separated list or a file with one value per line) and `is_child`. `match` and `nre` take the `regex` crate's dialect, case-insensitive and unanchored: a back-reference or an escape the crate lacks is an error naming the pattern. A record whose every transcript is filtered out is intergenic |
| `--exclude_predicted` | bool | false | Drop RefSeq's predicted transcript models, the `XM_` and `XR_` accessions, by adding the filter `not stable_id match ^X._` after the `--transcript_filter` values, as Ensembl VEP's option set does. An Ensembl cache carries no predicted model, so the flag changes nothing there; it takes effect on a converted RefSeq or merged cache (no RefSeq or merged cache is built natively) |
| `--allele_number` | bool | false | Add allele number from input VCF |
| `--show_ref_allele` | bool | false | Show reference allele in output |
| `--uploaded_allele` | bool | false | Add the UPLOADED\_ALLELE column: the record's allele string as uploaded, REF and every ALT joined by / |
| `--ambiguity` | bool | false | Add the AMBIGUITY column: the IUPAC ambiguity code of the record's alleles (R for A/G); absent on an indel or any other allele set without a code |
| `--minimal` | bool | false | Trim every bi-allelic record's alleles from both ends and split a multi-allelic record whose alleles trim into one allele per row, each on its own trimmed span, as Ensembl VEP minimises; switches on `--allele_number` as in Ensembl VEP |
| `--total_length` | bool | false | Append /\<total length\> to the cDNA, CDS and protein positions: the transcript's cDNA length, its coding sequence length and its peptide length. JSON keeps Ensembl VEP's shape: a single base prints the string "428/1200" as both cdna\_start and cdna\_end, and a range prints its start as both |

## Transcript sets

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--refseq` | bool | false | Reads the cache under the species directory suffixed \_refseq (the transcripts annotated are that cache's contents), and resolves an HGVS gene name among the RefSeq transcripts before the Ensembl ones, as Ensembl VEP's `--refseq` orders its core groups; under `--json_cache` only that order applies |
| `--merged` | bool | false | Reads the cache under the species directory suffixed \_merged (the transcripts annotated are that cache's contents); no effect under `--json_cache` |
| `--gencode_basic` | bool | false | Annotate against the GENCODE basic transcript set alone: a transcript without a gencode\_basic attribute is dropped as the cache loads, before the transcript index is built |
| `--gencode_primary` | bool | false | Annotate against the GENCODE primary transcript set alone: a transcript without a gencode\_primary attribute is dropped as the cache loads; the set exists for human GRCh38 only, so `--assembly` naming another assembly is warned about |
| `--flag_gencode_primary` | bool | false | Add the GENCODE\_PRIMARY column: 1 on a transcript of the GENCODE primary set, empty otherwise; the set exists for human GRCh38 only, so `--assembly` naming another assembly is warned about |
| `--use_transcript_ref` | bool | false | Read each row's reference allele from the transcript instead of the input, and add GIVEN\_REF (the input's reference allele) and USED\_REF (the row's) to the output. The row's allele is the transcript's spliced cDNA over the variant where the variant lies wholly inside the exons, `-` for an insertion, and the input's allele elsewhere; an alternate allele equal to the transcript's sequence takes the terms an allele equal to the reference takes, and the reference row of `--process_ref_homs` is the transcript's allele. Needs `--fasta`, which supplies a non-coding transcript's cDNA; `--use_given_ref` switches it off. Not read: the edited transcript sequences of a BAM-edited (RefSeq) cache, which Ensembl VEP reads with this flag switched on by the cache itself |
| `--use_given_ref` | bool | false | Accepted without effect: the default it names, the input's reference allele, stands without it, and given beside `--use_transcript_ref` it keeps that default |

