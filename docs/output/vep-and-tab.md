# VEP default and tab formats

This page describes the two tab-delimited text formats. The default format gathers the optional fields into one `Extra` column. The `tab` format gives every Extra key a named column of its own. The rules that decide what a row says, where a reader of Ensembl VEP's output might expect something else, close the page.

## VEP default format

One line per (allele, transcript) consequence, plus one intergenic line per allele that overlaps no feature.

### Header

```
## ENSEMBL VARIANT EFFECT PREDICTOR v116.2
## vep-rs 0.5.0
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

The default and `tab` headers are Ensembl VEP's `##` block. It opens with the version line, then the engine line `## vep-rs <version>` ([Engine identity](../output-formats.md#engine-identity)), the time, the cache path and the API line. One `## <source> version <value>` line per source version of the cache follows, then the column descriptions, the command line and the column line.

Two kinds of line Ensembl VEP prints there have no counterpart. The Perl API component versions (`## ensembl version 116.0d85231`, `## ensembl-variation version ...`) are not written, since no Perl API runs here. The `## Connected to <database>` line of a database run is not written, since vep-rs reads a local cache only. The `## <short name> : [PATH]/<file>` lines of `--custom` sources are not written either, since custom annotation sources are not supported. The VCF header's `##INFO` definition and plugin lines follow the same rule ([VCF format](vcf.md#header)).

The `Extra column keys` block lists every key the active flags can produce. Then comes one `## <KEY> : <description>` line per field of each loaded built-in plugin, with the plugin's own text, as Ensembl VEP lists its plugins' header info. `--fields` changes nothing in the default format: it keeps its thirteen columns and its `Extra` column whatever the list names, as in Ensembl VEP.

Every text format opens with a header block. `--no_headers` suppresses it whole, the column line included, in the default, `tab` and `vcf` formats; `json` has none. The Parquet TSV intermediate keeps its column line, which DuckDB reads as the schema, whatever the flag.

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

Keys appear in Ensembl VEP's order: `IMPACT`, `DISTANCE`, `STRAND` and `FLAGS`, with `ALLELE_NUM` under `--allele_number` and `REF_ALLELE` under `--show_ref_allele` before them. Then come the fields of the active flags in VEP's flag order (the table names the common ones, in that order), then every other key alphabetically: the structural-variant `OverlapBP` and `OverlapPC`, and plugin keys. The co-located variant keys, the frequency keys, `FREQS` and `SV` have [a page of their own](known-variants.md).

| Key                                              | Type or values                                                             | Notes                                                                                            |
| ------------------------------------------------ | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `ALLELE_NUM`                                     | integer, from 1                                                            | Under `--allele_number`; precedes `IMPACT`                                                       |
| `REF_ALLELE`                                     | allele string                                                              | Under `--show_ref_allele`; precedes `IMPACT`                                                     |
| `IMPACT`                                         | `HIGH`, `MODERATE`, `LOW`, `MODIFIER`                                      | On every row; an intergenic row carries `IMPACT=MODIFIER` and no `STRAND`, `DISTANCE` or `FLAGS` |
| `DISTANCE`                                       | integer, bases                                                             | Up- and downstream rows only; the shortest gap between the variant and the transcript            |
| `STRAND`                                         | `1` or `-1`                                                                | The feature's strand                                                                             |
| `FLAGS`                                          | `cds_start_NF`, `cds_end_NF`                                               | Every `cds_` transcript attribute, comma-joined                                                  |
| Flag fields: class, notations and gene           | `VARIANT_CLASS`, `SPDI`, `GA4GH_VRS`, `CLINVAR_SOMATIC_CLASSIFICATION`, `SYMBOL`, `SYMBOL_SOURCE`, `HGNC_ID`, `BIOTYPE` | In VEP's flag order, each only under its own flag; `VARIANT_CLASS` also joins an intergenic row  |
| Flag fields: transcript                          | `CANONICAL`, `MANE`, `TSL`, `APPRIS`, `CCDS`, `ENSP`, `NEAREST`            | Same order rule; `NEAREST` also joins an intergenic row                                          |
| Flag fields: predictions and structure           | `SIFT`, `PolyPhen`, `EXON`, `INTRON`, `DOMAINS`                            | Same order rule                                                                                  |
| Flag fields: HGVS                                | `HGVSc`, `HGVSp`, `HGVS_OFFSET`                                            | Same order rule; the `=` of an `HGVSp` value is escaped as `%3D`                                 |
| Flag fields: frequencies and co-located variants | the allele frequencies, `MAX_AF`, `FREQS`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `PUBMED` | Same order rule; `CLIN_SIG` also joins an intergenic row; [Co-located variants](known-variants.md) |
| Flag fields: the last two                        | `SHIFT_LENGTH`, `VAR_SYNONYMS`                                             | Same order rule; after every other flag field                                                    |
| `OverlapBP`, `OverlapPC`                         | integer; percentage                                                        | Structural-variant rows that overlap their feature; alphabetical after the flag fields           |
| Plugin keys                                      | per plugin                                                                 | Alphabetical after the flag fields; a variant-level plugin key also joins an intergenic row      |

Inside a value `;` is escaped as `%3B` and, unless `--no_escape`, the `=` of an HGVSp value as `%3D`; the JSON output, unescaped as under Ensembl VEP's `--no_escape`, keeps the `=`. Lists join with `,`. An `Extra` column with no key prints `-`.

#### Consequence vocabulary and total lengths

`--terms display` or `--terms NCBI` renders every consequence term in that vocabulary, repeats included on a feature row (`SPLICE_SITE,SPLICE_SITE,INTRONIC`) and each label once on a `--summary` row. The ranks, the impacts and the most severe term stay on the Sequence Ontology term. `DISTANCE` is printed only when the rendered list names `upstream_gene_variant` or `downstream_gene_variant`, so the display vocabulary, whose labels are `UPSTREAM` and `DOWNSTREAM`, never prints it, as Ensembl VEP does not.

`--total_length` appends `/<total>` to `cDNA_position` (the transcript's exon length sum), `CDS_position` (the coding sequence length) and `Protein_position` (the peptide length) wherever the column applies, a `-` placeholder included (`-/1200`).

#### HGVS_OFFSET and SHIFT_LENGTH

`HGVS_OFFSET` is the number of bases the HGVS notations moved an insertion or deletion 3' along the transcript, negative on the reverse strand, present only when something moved. Under `--shift_hgvs 0` nothing moves and the key never appears.

`SHIFT_LENGTH` (`--shift_length`) is the number of bases `--shift_3prime 1` moved the allele before its terms and positions were computed. It is written on every transcript row, `0` included, and on no intergenic row; without `--shift_3prime 1` the column is listed and left empty. Under `--shift_3prime 1 --shift_genomic 1` the `Location` column prints the moved position (a transcript row along its transcript's strand, an intergenic row towards higher coordinates) and `Allele` the inserted bases rotated along the move; `Uploaded_variation` keeps the input's name.

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

Selected with `--tab` or `--output_format tab`. A row holds the 13 default columns, then one column per Extra key the active flags can produce, in the Extra key order above, then one column per built-in plugin field. An absent value prints `-`; lists join with `,`. The header repeats the default format's block with a `## <column> : <description>` line for every column (a plugin column takes the plugin's own description) and no `Extra column keys` block. Its column line names every column.

### Choosing the columns (--fields)

`--fields` names the columns, comma-separated, and the list stands in for the whole default set, defaults included, in the order given: `--fields Location,Allele,Feature,Consequence,SYMBOL` writes exactly those five columns, each with its description line. Any field a flag or a plugin produces is a valid name; a field whose flag is off prints `-`, so `SYMBOL` needs `--symbol`. A name nothing produces prints `-` in every row under a `## <name> : ?` description line, as Ensembl VEP writes it. Nothing is trimmed and a trailing comma opens no column.

An empty list writes the tab header block and no data line, as Ensembl VEP's output loop does (every row is the empty string, which ends it); the VCF output's entries are then empty (`CSQ=,,`) on both engines. That loop also ends at a row whose whole line is `0`, so `--tab --fields SIFT` makes Ensembl VEP stop at the first SIFT score of `0` and lose every later row; vep-rs writes every row.

The list applies to the `tab` and `vcf` outputs; the default and JSON formats ignore it, as Ensembl VEP's do. `--output_format parquet` refuses it, since the Parquet column set is fixed (`ERROR: [output] --fields applies to --tab and --vcf output; the Parquet column set is fixed`, one of the [rules of vep-rs's own](../cli-reference.md#rules-of-vep-rss-own)). The `GRCh38-writers` corpus runs a list that reorders the columns and names both a field whose flag is off and a field nothing produces, in every format.

## Differences a reader will notice

Every formatter follows Ensembl VEP's own rules, cited to their Perl sources in the code: `crates/vep-io/src/output/fields.rs` holds the single field table every format shares, `crates/vep-effects/src/display.rs` the positions, codons and peptides, and the formatters sit beside them. The rules below are the ones a reader of these two formats most often needs. Structural variants have [a page of their own](structural-variants.md).

### Row order

Rows of one input record are transcript-major and allele-minor: sorted by transcript stable id, then by ALT order within a transcript. Intergenic rows come last. A record's rows are contiguous in the default and tab formats; the VCF output folds them into one line and the JSON output into one object.

### Uploaded_variation

The column takes the first `;`-element of a VCF record's ID; an Ensembl-format ID column is taken whole. The element is not used when it is `.`, empty or `0`, the values VEP reads as no ID, or when it is the `chr_start_alleles` name VEP generates for an empty ID column. So `rs2;.` names the record `rs2`, while `.;rs2`, `;rs2`, `0` and `0;rs2` do not name it.

A record with no ID is named from its own input line, as Ensembl VEP 116 names it: `CHROM_POS_REF/ALT1/ALT2`, with POS and every allele exactly as the VCF line spells them, case included. `21_10795228_C/G/T` is one such name. `21_30425650_a/g` names the lowercase `21 30425650 . a g`. `21_47612512_C/CTG` names the insertion `21 47612512 . C CTG`, whose variant starts at 47612513. `21_33867341_C/<CN2>` names the symbolic record `21 33867341 . C <CN2>`, whose variant starts at 33867342.

An Ensembl-format line is named `chr_start_alleles`, start the lower of its two coordinates. An Ensembl-format line whose alleles are written in lower case is the exception: it is named at the start and from the alleles its record carries after Ensembl VEP's own minimisation, the name Ensembl VEP builds before it upper-cases them. `21_30348410_ta/-` names `21 30348409 30348411 cta/c`. An HGVS or SPDI input is named by its notation as written, whitespace removed (`21:g.30282169C>A`, `NC_000021.9:10538884:A:C`), and a region input by its line as written (`21:25585733-25585733:1/G`), in both engines.

### Positions

`cDNA_position`, `CDS_position` and `Protein_position` are `start-end` ranges, printed as a single value when both ends coincide. An end that falls in an intron or outside the transcript prints `?`; for the CDS and protein columns an end in a UTR prints `?` as well (`?-335`, `493-?`). An insertion prints the two flanking bases (`1161-1162`). `cDNA_position` appears only on exon-overlapping rows; the CDS and protein columns, `Amino_acids` and `Codons` only on coding rows.

### Amino_acids and Codons

`Amino_acids` holds the `ref/alt` peptides, collapsed to one value when unchanged (`P` for a synonymous change). A partial codon appends `X` (`M/X`, `D/DX`). An allele whose codon window is empty prints `-` (`-/S` for an in-frame insertion, `M/-` for a whole-codon deletion). `Codons` holds the reference codon window in lower case with the allele's own bases in upper case, on both sides (`gTGCTGGgc/ggc` for a deletion, `-/TCT` for an insertion). Structural variants carry positions but no codons.

### DISTANCE

The value is the shortest gap between the variant and the transcript: the smallest of the four gaps between the record's start or end and the transcript's start or end. So the pair `POS+1..POS` of a zero-length structural variant measures from whichever flank is nearer. A `<CNV:TR>` contraction read under the cap measures from the run it is located on, not from the units it removes. The key is printed on up- and downstream rows only.

### FLAGS, OverlapBP and OverlapPC

`FLAGS` lists every `cds_` transcript attribute (`cds_start_NF,cds_end_NF`).

**OverlapBP** and **OverlapPC** appear on every structural variant row that overlaps its feature, breakends excepted, in the default Extra column and in JSON (`bp_overlap`, `percentage_overlap`). `--overlaps` adds them as columns to the `tab` and `vcf` outputs, as in Ensembl VEP, and to Parquet. They are reported on transcript rows only: Ensembl VEP also reports them on regulatory-feature and motif-feature rows, and regulatory annotation is not implemented.

### ENSP

Under `--protein` every transcript row of a transcript with a translation carries its ENSP, structural-variant rows (`transcript_ablation`, `feature_truncation`, breakends) included.

### The row filters

`--no_intergenic` and `--coding_only` run before the pick family, where Ensembl VEP's `get_all_VariationFeatureOverlapAlleles` applies them. `--no_intergenic` removes the intergenic row of every allele that overlaps no transcript. `--coding_only` keeps the rows whose allele affects the transcript's coding sequence, as Ensembl VEP's `within_cds` judges it: a CDS position, a coding consequence term where the position columns print `-`, or a variant covering a translated transcript whole. A deletion whose ends both lie outside the coding sequence keeps its row by the start, stop or coding-sequence term its inner span earns. Every other transcript row and the intergenic rows go.

A record left with no row is written as an unannotated record: no row in the default and tab formats, its VCF line without a CSQ entry, and in JSON the record's identity, coordinates and allele string with `"most_severe_consequence":"?"` and no consequence list. A pick flag then ranks the rows the filters leave, so `--coding_only --pick` keeps one coding row per record and none for a record without one.

### The pick family

The pick family acts on a whole input record, as Ensembl VEP's `filter_VariationFeatureOverlapAlleles` ranks every overlap allele of a variation feature together. The table gives each flag's effect; the `PICK` column exists under `--flag_pick`, `--flag_pick_allele` and `--flag_pick_allele_gene` alone.

| Flag                      | Effect on the record                                              | On a structural record                          |
| ------------------------- | ----------------------------------------------------------------- | ----------------------------------------------- |
| `--pick`                  | keeps one row for the record, whatever its allele count           | the same                                        |
| `--per_gene`              | keeps one row per gene over the record                            | the same                                        |
| `--flag_pick`             | marks one row of the record with `PICK=1`                         | the same                                        |
| `--pick_allele`           | groups the rows by the Allele column first, then picks            | collapses to `--pick`                           |
| `--pick_allele_gene`      | groups the rows by the Allele column first, then picks per gene   | collapses to `--per_gene`                       |
| `--flag_pick_allele`      | groups the rows by the Allele column first, then flags            | collapses to `--flag_pick`                      |
| `--flag_pick_allele_gene` | groups the rows by allele, then gene, then flags                  | one flagged row per gene over the whole record  |

The ranking follows `--pick_order` (default `mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq`): each criterion in turn, a lower value better, the rows tied on it compared on the next. `rank` is the lowest consequence rank of the row. `length` is the translateable-sequence length of a translated transcript and the spliced length otherwise, longer better. `ensembl` and `refseq` read the source cache of a merged cache and tie on every other cache. A shorter list falls through to the first row in Ensembl VEP's order, transcript-major by stable id and allele-minor, so a full tie keeps the lexicographically first stable id.

A record overlapping no transcript has one intergenic row per allele, which `--pick` keeps for the first allele alone and the flagging variants mark like a transcript row. Ensembl VEP emits the gene buckets of `--per_gene`, the allele buckets of `--pick_allele` and both of `--pick_allele_gene` in Perl hash order, which differs from run to run. vep-rs keeps the row order above, and the golden runs of those three sets compare the rows without their order (`"unordered": true` in the corpus manifest).

### --most_severe and --summary {#most_severe-and-summary}

`--most_severe` and `--summary` squash the rows left after the pick family into one row per record. The row carries `Uploaded_variation`, `Location` and `Consequence` (the one lowest-ranked term, or every term once in rank order), and `-` in every other fixed column, `Existing_variation` included even under `--check_existing`. Its Extra or tab cells hold the variant-level keys alone (`VARIANT_CLASS`, `MINIMISED`, `NEAREST`, `SV`) or `-`; there is no IMPACT, STRAND, ALLELE_NUM or plugin key. `--most_severe --per_gene` summarises the per-gene picks, which can differ from the unpicked set when a gene's picked transcript is not its most severe.

The JSON object files the row under `intergenic_consequences` as `{"consequence_terms": [...]}` and keeps `most_severe_consequence`, the variant-level keys and `colocated_variants` at the top level. `--vcf` is refused with either flag, as Ensembl VEP refuses it. Each also refuses the per-feature flags of Ensembl VEP's incompatibility list (`biotype`, `symbol`, `canonical`, `pick`, `flag_pick`, ...) with its message ([incompatible flags](../cli-reference.md#incompatible-flags-configpm615-632)).

### NEAREST

Under `--nearest transcript|gene|symbol` every row of a record, intergenic rows and the `--summary` row included, carries the identifier(s) of the protein-coding transcription start site(s) nearest the record, found as Ensembl VEP's `TranscriptTree::nearest` finds them. The start site of a protein-coding transcript is its start on the forward strand and its end on the reverse.

The two candidates are the highest site at or below the record's end and the lowest site at or above its start plus two. A site exactly one base above the start is the interval containing the start in Ensembl VEP's tree and is never a candidate. Each candidate is measured to the nearer end of the record, and both are reported when they tie, the site at or below the end first. A site on an end of the record is at distance zero.

The identifier is the transcript stable id, the gene stable id, or the gene symbol falling back to the gene stable id; values repeat once. Where several transcripts start at the winning coordinate, Ensembl VEP reports the one its interval tree holds first, which depends on the order its cache directory was listed in. vep-rs reports every transcript at that coordinate in stable-id order; `gene` and `symbol` agree with Ensembl VEP whenever those transcripts belong to one gene.

The site table is built from the run's own transcript set, so `--gencode_basic` or `--gencode_primary` narrows it. Ensembl VEP writes its table (`transcript_gene_tss.txt`) under the filters of the first `--nearest` run on a cache directory and reuses that file afterwards. A chromosome with no protein-coding transcript prints `NEAREST=` with no value in the default format, an empty tab cell and `"nearest": []`, as Ensembl VEP does. The JSON key is the top-level `nearest` array; the VCF CSQ field joins values with `&`.

### HGVSg

Under `--hgvsg` (which needs a reference FASTA) every row of a sequence variant, intergenic rows included, carries the HGVS genomic notation of its allele, as Ensembl VEP's `hgvs_genomic` writes it. The notation names the chromosome as the input names it, uses `g.` numbering, or `m.` on the mitochondrial chromosome, and reads the reference bases from the FASTA. An insertion or deletion is described at its most 3' position along the chromosome, and an insertion repeating the bases beside it as a duplication (`21:g.46925315_46925323dup` for `21 46925297 C CCCCCCTGGG`). A same-length change that is the reverse complement is an inversion, and any other change is clipped to the bases that differ (`21:g.34775851T>C` for `21 34775850 ATG ACG`).

Each allele of a multi-allelic record is trimmed against the reference on its own to place and class it. An insertion or deletion is then written from the trimmed bases, any other allele with the bases the record's working allele string holds for it, at the trimmed span, as Ensembl VEP writes it. `21:g.34775850delinsGT` is the `CGT` allele of `21 34775849 CAT C,CGT`, whose working allele string is `AT/-/GT`. An allele outside `ACGT`, an allele equal to the reference bases and a structural variant carry no value. The JSON key is `hgvsg`.

### SPDI

Under `--spdi` every row of a sequence variant carries its allele in SPDI notation, as Ensembl VEP's `spdi_genomic` writes it from the record as the parser leaves it: the chromosome, the 0-based position, the deleted bases and the inserted bases joined by `:`. `21:10538884:A:C` is `21 10538885 A C`, `21:10578406::A` an insertion, `21:20000325:T:` a one-base deletion, and `21:30000000:25:` a deletion of more than 20 bases, written as their count. A multi-allelic record's alleles are the anchor-trimmed pair at the record's span, so `21 14144642 AT CT,A` prints `21:14144641:AT:CT` and `21:14144641:AT:A`. An allele given on the reverse strand is complemented; nothing is shifted.

The chromosome is named as the cache names it on every row overlapping a transcript, and on an intergenic row as the reference FASTA names it. Ensembl VEP needs a FASTA for those rows and stops without one; vep-rs then prints the cache's name.

Ensembl VEP keeps the record's notations by allele and looks them up by the row's `Allele`. So under `--shift_3prime 1 --shift_genomic 1` a row whose inserted bases the shift rotated carries no `SPDI`, unless the rotation lands on another allele of the record, in which case the row carries that allele's notation. For `21 13620283 G GTCA,GCAT` the `TCA` allele moved 10 bases prints `CAT` as its Allele and `21:13620283::CAT` as its `SPDI`. A structural variant and an allele outside `ACGTN` carry none. The JSON key is `spdi`.

### GA4GH_VRS

Under `--ga4gh_vrs` (refused with `--vcf`) every row that has an SPDI, whether or not `--spdi` prints it, also carries the GA4GH Variation Representation Specification `Allele` Ensembl VEP derives from that SPDI. The object is a `SequenceLocation` on `refseq:<NC_ accession>` with a `SimpleInterval` from the SPDI's position to that position plus the deleted length, and a `SequenceState` holding the inserted bases; it has no digest.

Ensembl VEP builds it only from an SPDI named by an `NC_` accession, which it has with a database, or offline on an intergenic row alone when the FASTA record carries that name. vep-rs takes the accession from the chromosome synonyms table (`--synonyms` or the cache's `chr_synonyms.txt`, the first `NC_` synonym of the cache's chromosome name) and writes the object Ensembl VEP writes with a database, or nothing for a chromosome with no such synonym.

The text formats carry the object JSON-encoded on one line, `start` a string and `end` a number as Ensembl VEP's encoder writes them, with the keys in a fixed order where Ensembl VEP's follow its hash order. The JSON format nests the object under `ga4gh_vrs`, both bounds numbers.
