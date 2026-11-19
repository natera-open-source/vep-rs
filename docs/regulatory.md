# Regulatory and motif features

`--regulatory`, which `--everything` switches on (`Config.pm` line 362), annotates each record against the regulatory features and transcription-factor binding motifs of the cache, as Ensembl VEP does from the `_reg.gz` files beside its transcript files (`AnnotationSource/Cache/RegFeat.pm`).

## The rows

The rows, terms and fields follow Ensembl VEP's `OutputFactory.pm` lines 1867 to 1942.

A record inside a regulatory feature gets a `RegulatoryFeature` row; a record inside a motif gets a `MotifFeature` row. A record that covers the whole feature adds a term to its row.

| The record | `RegulatoryFeature` row | `MotifFeature` row |
| --- | --- | --- |
| lies inside the feature | `regulatory_region_variant` | `TF_binding_site_variant` |
| covers the whole feature as a deletion | adds `regulatory_region_ablation` | adds `TFBS_ablation` |
| covers the whole feature as a duplication, a tandem duplication or a tandem repeat, a sequence variant whose allele repeats its reference included | adds `regulatory_region_amplification` | adds `TFBS_amplification` |

Five of the six terms are `MODIFIER`; `TFBS_ablation` is `MODERATE`. An insertion overlaps a feature when both of its flanking bases lie inside it.

A record that overlaps no transcript keeps its `intergenic_variant` row beside its regulatory and motif rows, with `-` in the Gene, Feature and Feature_type cells. The rows of a record come in Ensembl VEP's order:

1. transcript rows
2. regulatory-feature rows
3. motif-feature rows
4. the intergenic row

A sequence variant's `--regulatory_gff` motif rows stand in stable-id order, since such a motif has no database id. Ensembl VEP orders them by the text of a read counter ([the motif row order class](divergences/regulatory-gff-motif-rows.md#1-two-motif-rows-of-one-record-in-the-order-of-a-read-counter)). Two motifs at one position can therefore come in the other order there, as the records before them in the input decide.

In JSON the rows are the `regulatory_feature_consequences`, `motif_feature_consequences` and `intergenic_consequences` arrays, and `most_severe_consequence` ranks every row's terms.

A structural variant's regulatory and motif rows come as their source lists them, the same for every input. The source's order is the cache's order region by region ascending, or the file order of `--regulatory_gff`. Ensembl VEP writes them in the order the input's records first reach each 1 Mb region ([the structural row order class](divergences/regulatory-structural-rows.md#1-a-structural-variants-rows-in-the-order-its-input-reaches-each-1-mb-region)). A record before the variant, a breakend's mate on either side of it or an earlier buffer can therefore put another region's rows first there.

The pick, summary and filter flags treat regulatory and motif rows as follows:

- `--pick`, `--pick_allele`, `--flag_pick` and `--flag_pick_allele` rank regulatory and motif rows with the transcript rows.
- `--most_severe` and `--summary` read the terms of regulatory and motif rows.
- `--per_gene` and `--pick_allele_gene` keep every regulatory, motif and intergenic row beside the one transcript row they keep per gene. They write those rows first, then the transcript picks, the order `pick_VariationFeatureOverlapAllele_per_gene` returns them in (`OutputFactory.pm` lines 838 to 856).
- `--flag_pick_allele_gene` flags every regulatory, motif and intergenic row beside the picks (`OutputFactory.pm` line 626).
- `--no_intergenic` drops the rows of a record with no transcript row.
- `--coding_only` drops every regulatory and motif row.

## The fields

A `RegulatoryFeature` row carries `BIOTYPE` whether or not `--biotype` is given, and no `STRAND`. `BIOTYPE` is the feature's type, such as `promoter`, `enhancer` or `CTCF_binding_site`. A `MotifFeature` row carries:

- `STRAND`
- `MOTIF_NAME`, the binding matrix's stable id
- `MOTIF_POS`, the position of the record's first base in the matrix, counted along the motif's own strand
- `TRANSCRIPTION_FACTORS`, the matrix's factors, comma-separated; empty when the matrix names none

The fields a transcript row alone carries print `-`, and the record-level fields print on every row.

A plugin annotates the rows of the feature types it declares, transcript rows unless it says otherwise ([Writing a plugin](developers/writing-a-plugin.md)). A plugin that declares regulatory or motif features switches `--regulatory` on.

### Motif scores

`HIGH_INF_POS` and `MOTIF_SCORE_CHANGE` are computed from the matrix's position counts as Ensembl VEP computes them (`MotifFeatureVariationAllele.pm` lines 153 to 265, funcgen's `BindingMatrix.pm` and `BindingMatrix/Converter.pm`). Each position's counts give:

- probabilities: each count plus a pseudocount of 0.1, over the position's total plus 0.4
- weights: the log2 of each probability over 0.25
- bits: each probability times `2 - H`, with `H` the position's entropy

A position is informative when its bits exceed 1.5. The relative similarity score of a sequence is the sum of its bases' weights. The sum is placed between the sums of the smallest and of the largest weight of every position. The motif's best site therefore scores 1 and its worst 0.

`HIGH_INF_POS` is `Y` when the record changes one base at an informative position and `N` for any other record, a deletion included. `MOTIF_SCORE_CHANGE` is the motif's score with the allele written over its sequence minus the motif's own score, to three decimals. It is written when both of these hold:

- the alternate allele is A, C, G and T bases of the reference allele's length
- the record starts inside the motif

An allele running past the motif's end is scored on the positions inside it.

The alleles are read on the forward strand of the genome whatever strand the record names. An Ensembl-format, region or HGVS record given on the reverse strand therefore scores as its forward-strand form does. On a reverse-strand motif the allele is reverse-complemented and written over the positions it covers. Ensembl VEP writes a multi-base allele there one or more positions too far along the motif ([the divergence page](divergences/motif-score-reverse-strand.md)).

A motif whose matrix carries no counts writes neither field, as Ensembl VEP omits both for a matrix without elements.

### Structural variants and breakends

A structural variant's regulatory and motif rows carry none of `BIOTYPE`, `STRAND` and the motif fields. They carry `OverlapBP` and `OverlapPC` as a transcript row does ([VEP default and tab](output/vep-and-tab.md)). Their terms follow the structural class:

- A `<DEL>`, a `<DEL:ME>` or a lone `<CN0>` covering the feature ablates it.
- A `<DUP>`, a `<DUP:TANDEM>`, a `<CNV:TR>` or a lone `<CN2>` amplifies it.
- A `<CNV>` or a list of copy numbers (`<CN0>,<CN2>`) carries the overlap term alone.

A paired breakend's bracket allele gets a row, with the overlap term, for each feature its mate lies in. It keeps its own position's features only when the mate lies within 5000 bases of them on the same chromosome. The own allele gains rows, with `intergenic_variant`, for the features a same-chromosome mate lies in within 5000 bases of it.

A feature's row is written once per allele. Ensembl VEP writes it once per listing of the record inside the feature, four times when both ends lie in it ([the breakend listing class](divergences/regulatory-rows.md#4-a-breakend-whose-two-ends-lie-in-one-regulatory-feature)).

A mate on another chromosome is looked up on that chromosome alone. Ensembl VEP's regulatory lookup names no chromosome. It therefore also writes the own allele's `intergenic_variant` row for a feature of the record's chromosome that the mate's position number falls in. It does so when the buffer's span there reaches the feature ([the other-chromosome mate class](divergences/regulatory-rows.md#5-a-breakend-whose-mate-on-another-chromosome-selects-a-feature-of-its-own-chromosome)).

### The `CELL_TYPE` field

`--cell_type <list>` adds `CELL_TYPE`, and switches `--regulatory` on. A `RegulatoryFeature` row prints `<epigenome>:<activity>` for each listed epigenome the feature has a record for, in the order listed. The activity is one of:

- `ACTIVE`
- `INACTIVE`
- `POISED`
- `REPRESSED`
- `NA`

An epigenome the feature's map lacks is skipped, and a feature with none of them prints an empty value.

A `MotifFeature` row prints an empty value, since a motif carries no activity map. The list joins with `,` in the default and tab formats and `&` in a VCF entry, and is a JSON array.

The names are the cache's `cell_types`, and a name may hold a comma (`A549_(m,_58_y)`). The value is read by longest match against them, so `--cell_type A549_(m,_58_y),BJ_(m)` names two epigenomes. Ensembl VEP splits the value on every comma, so a name with one prints an empty `CELL_TYPE` there ([`--cell_type` names with a comma](divergences/cell-type-comma-names.md)).

A name the cache lacks ends the run before the output file and statistics report exist. The error is `ERROR: [cli] Cell type "<name>" unavailable; available cell types are:` with the list on the next line. On a cache without regulatory data the names are not checked, as Ensembl VEP checks them only when it reads regulatory features. The flag then adds the column and no row fills it.

A value that names no epigenome, `''` or `','`, is the flag absent. Ensembl VEP prints an empty `CELL_TYPE=` on every regulatory and motif row of the default format and `"cell_type":[]` in JSON ([the empty value](divergences/cell-type-comma-names.md#an-empty---cell_type-value)).

## The cache

The features come from an Ensembl VEP cache converted with `scripts/data/storable_to_json.pl` ([Cache setup](cache-setup.md#convert-an-existing-ensembl-vep-cache)). Beside `transcripts/<chr>/` the script writes:

- `regulatory/<chr>/<start>-<end>.json`, one shard per region holding the region's regulatory features and motif features
- `regulatory/binding_matrices.json`, every binding matrix the motifs name stored once under its stable id

In `info.json` it sets `regulatory: true`, `regulatory_schema` (the shape of the tree) and `cell_types`, the epigenomes of the regulatory build. `vep-cache-converter` copies the tree, the shards as `.json` or `.json.gz`, and writes the three keys only after it has been written. The Ensembl VEP caches that carry a regulatory build convert with it; the human caches do, on both assemblies.

The entries of the tree carry these values.

| Entry | Carries |
| --- | --- |
| a regulatory feature | its stable id, span, feature type, database id and activity |
| a motif feature | its stable id, span, strand, database id, score, its sequence on its own strand and the stable id of its matrix |
| a matrix of `regulatory/binding_matrices.json` | its name, source, threshold, length, counts per position (A, C, G, T) and transcription factor complexes |

The activity is one code per epigenome of `cell_types`, in order.

| Code | Activity |
| --- | --- |
| `A` | active |
| `I` | inactive |
| `P` | poised |
| `R` | repressed |
| `N` | not assessed |
| `-` | no record for that epigenome |

A run refuses a regulatory tree without `regulatory_schema` with `ERROR: [cache] regulatory shards of <dir> predate the current converter: reconvert the cache`.

A cache converted without the regulatory tree holds no regulatory features. Convert the Ensembl VEP cache again, or build the tree from the regulatory GFF files ([below](#regulatory-features-from-gff-files)).

Every motif feature names a matrix of `regulatory/binding_matrices.json`; the converter and the loader refuse a motif naming none as a corrupt shard.

The shards of a chromosome are read together, the first time a record of that chromosome is annotated. A run without `--regulatory` therefore never opens them.

On a cache with no regulatory data, `--regulatory` writes one warnings-file line. The line reads `--regulatory: this cache holds no regulatory data; no RegulatoryFeature or MotifFeature row is written`. The regulatory columns of the tab and VCF headers are declared and print empty. Ensembl VEP skips the source without a word (`CacheDir.pm` lines 168 to 176).

## Regulatory features from GFF files

`--regulatory_gff` reads the regulatory annotation files the Ensembl site distributes under `regulation/<species>/<assembly>/annotation/` in place of the cache's regulatory data, as Ensembl VEP reads them (`AnnotationSource/File/RegFeat.pm`). The value is a bare file name, the regulatory-feature GFF3, or comma-separated `key=value` pairs in the style of `--custom`, any one or more of them.

| Key | File |
| --- | --- |
| `file` | the regulatory-feature GFF3 |
| `motifs` | the transcription factor binding motif GFF3 |
| `emars` | the EMAR GFF |

```bash
vep -i input.vcf --cache --regulatory_gff \
  file=regulatory_features.gff3.gz,motifs=motif_features.gff3.gz,\
emars=EMARs.gff.gz
```

The flag switches `--regulatory` on. The files are annotation sources of the run. An input record on a chromosome they list and the cache lacks is therefore admitted without a skip warning. It gains the files' rows at its position, as [Chromosomes](custom-annotation.md#chromosomes) describes.

### Files and records that are refused

Each file is bgzip-compressed and tabix-indexed; the files Ensembl distributes are bgzipped, so `tabix -p gff <file>` writes the index a file lacks. A file that does not exist, or has no index, stops the run with an error naming the key and the file. A value the grammar refuses stops it too. These errors carry Ensembl VEP's texts. Each refusal comes before the output file and the statistics report exist, as in Ensembl VEP, so a refused run leaves neither behind.

A record whose start or end is not a whole number, whose start is 0 or whose end is below its start stops the run when its chromosome is first read, after the output file is opened, naming the record and the file. So does any other line whose start or end no number leads or whose end is below zero. A start or end no number leads at the file's first line leaves an index vep-rs refuses as malformed; such a start on a later sequence leaves none. Ensembl VEP reads on past a line its `tabix` could not place and, when the end was the fault, loses the later records of its block, a documented intended divergence, [`gxf_malformed_coordinate_line_refused`](divergences/gxf-malformed-coordinate-line-refused.md).

### What a record becomes

A record of the regulatory or EMAR file is a regulatory feature whose type prints as `BIOTYPE` when that type is one of:

- `promoter`
- `enhancer`
- `CTCF_binding_site`
- `open_chromatin_region`
- `promoter_flanking_region`
- `TF_binding_site`
- `regulatory_region`
- `EMAR`

A record of the motif file of type `TF_binding_site` is a motif feature. A record of any other type is skipped, with one warnings-file line per type and file naming the type as decoded. In that line the replacement character stands for a byte that is not UTF-8 where Ensembl VEP writes the byte.

A feature's identifier is its `ID`, else its `Name`. Each has every `%XX` escape decoded where Ensembl VEP prints them as written ([the attribute escapes class](divergences/regulatory-rows.md#2-a---regulatory_gff-attribute-value-holding-a-gff3-escape)). Without either it is `<chr>:<start>-<end>:<type>`. `<chr>` is the chromosome as the file spells it, whatever the cache or the input calls it. `<type>` is the type column with every `%XX` escape decoded, as both engines read it.

Two distinct records of one file sharing an identifier are one feature. The run writes Ensembl VEP's warning once per identifier and file, and an input record inside both gets the later record's row.

A byte that is not UTF-8 in a record is read as U+FFFD by `--regulatory_gff` and printed so ([input bytes that are not UTF-8](divergences/non-utf8-input.md)), a type so read skipping its record with a warning, where `vep-cache-builder` refuses the file, naming the line.

A record of the regulatory file and one of the EMAR file sharing an identifier are two features with two rows. Ensembl VEP keeps the EMAR row alone ([the shared identifier class](divergences/regulatory-rows.md#3-one-identifier-in-both-the-regulatory-file-and-the-emar-file-of---regulatory_gff)). An EMAR is a container region that overlaps many regulatory features. An input record inside a promoter or an enhancer therefore usually gains an EMAR row beside it.

### Motif rows

A motif's matrix is known from its record alone:

- `MOTIF_NAME` is the record's `binding_matrix_id`
- `TRANSCRIPTION_FACTORS` is its `transcription_factor` list
- `MOTIF_POS` is counted along the motif's strand over its span, with or without a FASTA

Ensembl VEP omits `MOTIF_POS` without a FASTA, since the stand-in slice it then builds has no strand ([the `MOTIF_POS` class](divergences/regulatory-gff-motif-rows.md#2-motif_pos-omitted-from-a-motif-row-without-a-fasta)). The files carry no position counts, so `HIGH_INF_POS` and `MOTIF_SCORE_CHANGE` are not written. They carry no epigenome activity, so `CELL_TYPE` is empty and `--cell_type` is warned about once and ignored.

### Row order

A structural variant's regulatory and motif rows come file by file in the flag's source order (`file`, `emars`, `motifs`), each file's records in file order. A sequence variant's rows stand in stable-id order.

Ensembl VEP's order of a structural variant's rows can change the variant's `--pick` row ([the structural row order class](divergences/regulatory-structural-rows.md#1-a-structural-variants-rows-in-the-order-its-input-reaches-each-1-mb-region)). The rows are otherwise the rows of the cache's features.

### Extended promoters

`--extended_promoters` widens each promoter to its `extended_start` and `extended_end` attributes, the wider region the files record beside the core promoter. The extension never narrows a feature, affects promoters alone, and changes nothing without `--regulatory_gff`.

A bound that is not a whole number is ignored, the promoter keeping its own bound on that side. The run writes one warnings-file line naming the record and the attribute. Ensembl VEP reads such a bound as 0 ([Defect N, class 27](divergences/defect-n.md#27-a-promoter-whose-extended_start-is-not-a-number-is-widened-to-position-0)).

A widened promoter is found for every input record its extension covers, however far its core lies from the input record's 1 Mb region. Ensembl VEP reads each region with a 10,000-base margin ([the seek margin class](divergences/regulatory-rows.md#1-a-promoter-whose-extension-reaches-past-the-seek-margin-of-the-next-1-mb-region)).

### The warning and the headers

When the cache holds regulatory data too, the run reads the files alone; the cache's transcripts are read all the same. It writes `WARNING: Using regulatory features from --regulatory_gff; ignoring the regulatory data in this cache`. The output headers name the first of the `file`, `emars` and `motifs` files as `regulatory_gff version <file>` among the source versions, as Ensembl VEP's do. `--show_cache_info` prints that entry too, as `regulatory_gff<TAB><file>` beside the cache's source versions, and refuses a file that is not there before printing.

### Baking the files into a cache

`vep-cache-builder --regulatory-gff <file> --emars-gff <file> --motif-gff <file> [--extended-promoters]` bakes the same files into a cache's regulatory tree ([Cache setup](cache-setup.md#regulatory-data)). The tree holds one shard per 1 Mb region of each feature's start. In `binding_matrices.json` it holds a stub of each motif's matrix (its span and factors, no counts). `info.json` carries `regulatory: true` and `regulatory_schema` and no `cell_types`.

Each baked feature records the file it came from and its place in that file. A run with `--regulatory` over that cache therefore writes the rows a run with `--regulatory_gff` over the files writes, the shared-identifier rules included. The bake itself writes the warnings the run-time reader writes.

## Where next

- [The divergence page of `--regulatory`](divergences/regulatory-rows.md): the warning on a cache without regulatory data, the `--spdi` chromosome name, the GFF index refusal and five defects of Ensembl VEP's regulatory file reader and breakend lookup.
- [The `--regulatory_gff` motif rows](divergences/regulatory-gff-motif-rows.md): two more defects, in the order of those rows and in `MOTIF_POS` without a FASTA.
- [The motif score on a reverse strand](divergences/motif-score-reverse-strand.md): the scoring of a multi-base allele on a reverse-strand motif.
- [`--cell_type` names with a comma](divergences/cell-type-comma-names.md): how a `--cell_type` name holding a comma is read.
