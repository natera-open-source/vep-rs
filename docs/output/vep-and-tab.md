# VEP default and tab formats

This page describes the two tab-delimited text formats. The default format gathers the optional fields into one `Extra` column. The `tab` format gives every Extra key a named column of its own. The rules that decide what a row says, where a reader of Ensembl VEP's output might expect something else, close the page.

## VEP default format

One line per (allele, transcript) consequence, plus one intergenic line per allele that overlaps no feature.

### Header

```
## ENSEMBL VARIANT EFFECT PREDICTOR v116.2
## Output produced at [TIME]
## Using cache in [PATH]/116
## Using API version 116, DB version ?
## 1000genomes version phase3
## COSMIC version 98
## ClinVar version 202306
## ...one line per source version recorded in the cache's info.json...
## Column descriptions:
## Uploaded_variation : Identifier of uploaded variant
## Location : Location of variant in standard coordinate format (chr:start or chr:start-end)
## ...one line per column...
## Extra column keys:
## IMPACT : Subjective impact classification of consequence type
## DISTANCE : Shortest distance from variant to transcript
## STRAND : Strand of the feature (1/-1)
## FLAGS : Transcript quality flags
## VEP command-line: vep --assembly GRCh37 --force_overwrite --input_file [PATH]/input.vcf --json_cache [PATH]/116 --offline --output_file [PATH]/out.txt
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
```

The source-version lines come from the cache's `info.json` (`source_versions`). That file is written by `scripts/data/storable_to_json.pl` and `vep-cache-builder`, or by `scripts/golden/prune_json_cache.py --perl-info` from a VEP cache's `info.txt`. The `Using cache in` line prints the cache directory as given. Under `--mask_header_cache_path` it prints `[PATH]/<leaf>`, the last component alone; the same masking applies in the `--tab` header, in the VCF `cache=` token and in the Parquet footer.

The `Extra column keys` block lists every key the active flags can produce; plugin keys are not listed. `--no_headers` suppresses the whole block, column line included.

### Columns

| #   | Column             | Description                                                                                     |
| --- | ------------------ | ----------------------------------------------------------------------------------------------- |
| 1   | Uploaded_variation | Input ID, else the record's own line as `CHROM_POS_REF/ALT1/ALT2` (POS and alleles as written)  |
| 2   | Location           | `chr:pos` or `chr:start-end` (`chr:end-start` for an insertion)                                 |
| 3   | Allele             | The ALT allele of this row; a structural variant prints the SO term of its class (see below)   |
| 4   | Gene               | Ensembl gene stable ID                                                                          |
| 5   | Feature            | Ensembl transcript stable ID                                                                    |
| 6   | Feature_type       | `Transcript`, or `-` on an intergenic row; VEP's `RegulatoryFeature` and `MotifFeature` rows are not produced |
| 7   | Consequence        | SO consequence term(s), comma-separated                                                         |
| 8   | cDNA_position      | Position or range in the cDNA; `?` for an unmappable end                                        |
| 9   | CDS_position       | Position or range in the coding sequence                                                        |
| 10  | Protein_position   | Position or range in the protein                                                                |
| 11  | Amino_acids        | Reference/variant peptides (`G/A`), one value when unchanged                                    |
| 12  | Codons             | Reference/variant codon windows, allele bases in upper case (`gGa/gCa`)                         |
| 13  | Existing_variation | Co-located known variant IDs, comma-separated                                                   |
| 14  | Extra              | `KEY=VALUE` pairs joined by `;`, in VEP's key order (see below); `-` when empty                  |

The Allele column of a structural variant carries the SO term of its class (`deletion`, `duplication`, `Alu_insertion`). A breakend keeps its literal ALT instead (`N.`, `A[21:27348902[`, or `N[21:00030425700[` for a mate named by `CHR2`). `chromosome_breakpoint` is printed for a record without a parsed mate, a single breakend or a bare `<BND>`. Every ALT of a multi-allelic structural variant prints the record's class term, and an unsupported `<CPX>` prints `CPX`. A sequence allele whose reference base is ambiguous (`N G`) prints `G` with `coding_sequence_variant` and no amino acids, since no peptide is read from an ambiguous reference.

### Extra field

Keys appear in Ensembl VEP's order: `IMPACT`, `DISTANCE`, `STRAND` and `FLAGS` first, then the fields of the active flags in VEP's flag order (the table names the common ones), then every other key alphabetically.

| Key                                              | Type or values                                                             | Notes                                                                                            |
| ------------------------------------------------ | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `ALLELE_NUM`                                     | integer, from 1                                                            | Under `--allele_number`; precedes `IMPACT`                                                       |
| `REF_ALLELE`                                     | allele string                                                              | Under `--show_ref_allele`; precedes `IMPACT`                                                     |
| `IMPACT`                                         | `HIGH`, `MODERATE`, `LOW`, `MODIFIER`                                      | On every row; an intergenic row carries `IMPACT=MODIFIER` and no `STRAND`, `DISTANCE` or `FLAGS` |
| `DISTANCE`                                       | integer, bases                                                             | Up- and downstream rows only; the shortest gap between the variant and the transcript            |
| `STRAND`                                         | `1` or `-1`                                                                | The feature's strand                                                                             |
| `FLAGS`                                          | `cds_start_NF`, `cds_end_NF`                                               | Every `cds_` transcript attribute, comma-joined                                                  |
| Flag fields: class and gene                      | `VARIANT_CLASS`, `SYMBOL`, `SYMBOL_SOURCE`, `HGNC_ID`, `BIOTYPE`           | In VEP's flag order, each only under its own flag; `VARIANT_CLASS` also joins an intergenic row  |
| Flag fields: transcript                          | `CANONICAL`, `MANE`, `TSL`, `APPRIS`, `CCDS`, `ENSP`                       | Same order rule                                                                                  |
| Flag fields: predictions and structure           | `SIFT`, `PolyPhen`, `EXON`, `INTRON`, `DOMAINS`                            | Same order rule                                                                                  |
| Flag fields: HGVS                                | `HGVSc`, `HGVSp`, `HGVS_OFFSET`                                            | Same order rule; the `=` of an `HGVSp` value is escaped as `%3D`                                 |
| Flag fields: frequencies and co-located variants | the allele frequencies, `MAX_AF`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `PUBMED` | Same order rule; `CLIN_SIG` also joins an intergenic row                                         |
| `OverlapBP`, `OverlapPC`                         | integer; percentage                                                        | Structural-variant rows that overlap their feature; alphabetical after the flag fields           |
| Plugin keys                                      | per plugin                                                                 | Alphabetical after the flag fields; a variant-level plugin key also joins an intergenic row      |

Inside a value `;` is escaped as `%3B`, and the `=` of an HGVSp value as `%3D`; the JSON output, unescaped as under Ensembl VEP's `--no_escape`, keeps the `=`. Lists join with `,`. An `Extra` column with no key prints `-`.

### Example

The example shows two records of the GRCh37 golden corpus as Ensembl VEP 116.2 wrote them (`tests/golden/116/GRCh37/expected/default.txt.gz`). The first, `1 213061271 806345 G C`, is a missense change reported against four transcripts. The second, the ID-less multi-allelic site `21 10795228 . C G,T`, overlaps no transcript.

```
806345	1:213061271	C	ENSG00000162769	ENST00000366971	Transcript	missense_variant	1433	1235	412	G/A	gGa/gCa	-	IMPACT=MODERATE;STRAND=1
806345	1:213061271	C	ENSG00000162769	ENST00000419102	Transcript	missense_variant	631	632	211	G/A	gGa/gCa	-	IMPACT=MODERATE;STRAND=1;FLAGS=cds_start_NF
806345	1:213061271	C	ENSG00000162769	ENST00000474693	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2364;STRAND=1
806345	1:213061271	C	ENSG00000162769	ENST00000483790	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
21_10795228_C/G/T	21:10795228	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
21_10795228_C/G/T	21:10795228	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

## Tab format

Selected with `--tab` or `--output_format tab`. A row holds the 13 default columns, then one column per Extra key the active flags can produce, in the Extra key order above, then one column per built-in plugin field. An absent value prints `-`; lists join with `,`. The header repeats the default format's block with a `## <column> : <description>` line for every column and no `Extra column keys` block. Its column line names every column.

## Differences a reader will notice

Every formatter follows Ensembl VEP's own rules, cited to their Perl sources in the code: `crates/vep-io/src/output/fields.rs` holds the single field table every format shares, `crates/vep-effects/src/display.rs` the positions, codons and peptides, and the formatters sit beside them. The rules below are the ones a reader of these two formats most often needs. Structural variants have [a page of their own](structural-variants.md).

### Row order

Rows of one input record are transcript-major and allele-minor: sorted by transcript stable id, then by ALT order within a transcript. Intergenic rows come last. A record's rows are contiguous in the default and tab formats; the VCF output folds them into one line and the JSON output into one object.

### Uploaded_variation

The column takes the first `;`-element of a VCF record's ID; an Ensembl-format ID column is taken whole. The element is not used when it is `.`, empty or `0`, the values VEP reads as no ID, or when it is the `chr_start_alleles` name VEP generates for an empty ID column. So `rs2;.` names the record `rs2`, while `.;rs2`, `;rs2`, `0` and `0;rs2` do not name it.

A record with no ID is named from its own input line, as Ensembl VEP 116 names it: `CHROM_POS_REF/ALT1/ALT2`, with POS and every allele exactly as the VCF line spells them, case included. `21_10795228_C/G/T` is one such name. `21_30425650_a/g` names the lowercase `21 30425650 . a g`. `21_47612512_C/CTG` names the insertion `21 47612512 . C CTG`, whose variant starts at 47612513. `21_33867341_C/<CN2>` names the symbolic record `21 33867341 . C <CN2>`, whose variant starts at 33867342. An Ensembl-format line is named `chr_start_alleles`, start the lower of its two coordinates. VEP names a region or HGVS input by the input string itself; vep-rs names it `chr_start_alleles` from the parsed variant.

### Positions

`cDNA_position`, `CDS_position` and `Protein_position` are `start-end` ranges, printed as a single value when both ends coincide. An end that falls in an intron or outside the transcript prints `?`; for the CDS and protein columns an end in a UTR prints `?` as well (`?-335`, `493-?`). An insertion prints the two flanking bases (`1161-1162`). `cDNA_position` appears only on exon-overlapping rows; the CDS and protein columns, `Amino_acids` and `Codons` only on coding rows.

### Amino_acids and Codons

`Amino_acids` holds the `ref/alt` peptides, collapsed to one value when unchanged (`P` for a synonymous change). A partial codon appends `X` (`M/X`, `D/DX`). An allele whose codon window is empty prints `-` (`-/S` for an in-frame insertion, `M/-` for a whole-codon deletion). `Codons` holds the reference codon window in lower case with the allele's own bases in upper case, on both sides (`gTGCTGGgc/ggc` for a deletion, `-/TCT` for an insertion). Structural variants carry positions but no codons.

### DISTANCE

The value is the shortest gap between the variant and the transcript: the smallest of the four gaps between the record's start or end and the transcript's start or end. So the pair `POS+1..POS` of a zero-length structural variant measures from whichever flank is nearer. A `<CNV:TR>` contraction read under the cap measures from the run it is located on, not from the units it removes. The key is printed on up- and downstream rows only.

### FLAGS, OverlapBP and OverlapPC

`FLAGS` lists every `cds_` transcript attribute (`cds_start_NF,cds_end_NF`).

`OverlapBP` and `OverlapPC` appear on every structural variant row that overlaps its feature, breakends excepted, in the default Extra column and in JSON (`bp_overlap`, `percentage_overlap`). They never appear in `tab`, `vcf` or Parquet: VEP's `--overlaps` flag, which adds them as columns there, is not implemented.

### ENSP

Under `--protein` every transcript row of a transcript with a translation carries its ENSP, structural-variant rows (`transcript_ablation`, `feature_truncation`, breakends) included.
