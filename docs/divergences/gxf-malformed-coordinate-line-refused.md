# A GFF3 or GTF line whose coordinates are malformed

**The defect in one sentence.** Ensembl VEP 116.3's `--gff`, `--gtf`, `--regulatory_gff` and `--custom` GFF sources read the records the `tabix` index returns. The `tabix` the release ships, htslib 1.9, indexes a file past a line inside a sequence whose end no number leads or reads below zero, warning and exiting 0. It does the same past a start no number leads at the file's first line, and fails the build when that line sits on a later sequence. A start of 0 and an end below its start, not below zero, it indexes silently. The source then builds its model without the unparsed line, or on the span as written, and reports no fault in the file.

A variant in a CDS whose line was lost is a `5_prime_UTR_variant`; a transcript whose gene line was lost is intergenic; a GTF transcript whose first CDS line was lost has every codon renumbered.

**Why this is a defect.** A start or end that is not a number and an end below its start are not coordinates; a start of 0 is none in GFF3 (positive one-based integers), and GTF 2.2 numbers from 1 while accepting and discouraging values outside the sequence. Nothing in the route says so. Where `tabix` exits 0, with its message or without, Ensembl VEP builds the rest of the transcript and prints its rows as on a sound file, and the only sign of the fault is a consequence that is wrong. On some shapes the run dies in a Perl method with a message naming no line of the file.

**What vep-rs does instead.** `vep-cache-builder` reads the file itself, line by line, and refuses the build at the first kept line whose start or end is not a number, whose start is 0 or whose end is below its start. The error names the line: `<file> line <n>: GFF3 start "<value>" is not a number`, or `<file> line <n>: GFF3 coordinates <start>-<end> are not 1-based and ordered` (`GTF` for a GTF); no cache is written. A line it skips, for its type or for holding four to eight columns, is read for one thing more: a start or end no number leads, or an end below zero, refuses the build there with the same error.

`--regulatory_gff` refuses a record whose start or end is not a whole number, whose start is 0 or whose end is below its start, naming the record. The builder's bake refuses the same, and both refuse a start or end no number leads, or an end below zero, on a line of any type. A `--custom` GFF or `--regulatory_gff` source is read through its index as Ensembl VEP reads it. An index that decompresses to more than 268435456 bytes is refused as malformed, the size `tabix` writes past an unparsed line at the file's head.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-io` at release/116, htslib 1.9):

| Source | What it does |
| --- | --- |
| `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, line 211, `_get_records_by_coords` | The records of a region come from the tabix index (`seek`). Line 222 reads each record's raw start as the loop bound; line 224 tests the type and line 237 keeps the record. A line the index does not return is no record. |
| `Bio/EnsEMBL/IO/ColumnBasedParser.pm`, line 250, `read_record`; `Parser/GXF.pm`, lines 163-166, `get_start` | The columns are split at the tabs and the start is returned as written; nothing reads a coordinate as a number. |
| htslib `tbx.c`, lines 95 and 105, `tbx_parse1` | The start and end columns are read with `strtol` in base 0, which takes blanks, a sign and a `0x` or `0` prefix, and reads on past an empty or blank column into the next. Text no number so leads fails the line, and so does an end that reads below zero (line 142); a start that does is read as the first base (line 99). |
| `tbx.c`, lines 235-236, `tbx_index` | The 1.9 build prints `Failed to parse TBX_GENERIC` and files the line all the same, under the sequence id of the record before it, a failed start read as 0 and a failed end keeping its line's start. The build fails, `Unsorted positions` and exit 1, when the start so filed is below that record's, as a failed start is after a record starting beyond the first base. At the file's first line the id was never set, and the index names 67.7 million sequences, 541 MB uncompressed. A query that reaches a line it cannot parse stops there and returns no later line. |
| `Bio/EnsEMBL/VEP/AnnotationSource/File/RegFeat.pm`, lines 291-325, `_get_regfeats_by_coords` | The regulatory records of a region come from the same `seek`, with the same raw start as the loop bound, so a line the index drops is no feature and a later line of its block is lost with it. |
| `BaseGXF.pm`, lines 363-435, `_get_parent_child_structure`; line 719, `_add_translation` | A record whose parent the fetched records lack is an orphan, warned about and dropped with its children. The translation runs from the first CDS record in transcript order to the last, so a lost first CDS line moves the coding start to the next one. |
| `BaseGXF.pm`, line 544, `_create_transcript`; lines 740-756, `_add_translation` | A CDS record is given to the exon that holds its start, and one whose start no exon holds is left out. The translation's ends are the first and last given records' distances from their exons' edges, with no test that either lies inside its exon. |

**What correcting it would change.** With the source refusing a line whose coordinates are not a one-based ordered pair of numbers, every malformed file would fail before a row is written. Every file whose coordinates are plain decimal numbers would read as it reads in 116.3; a zero-padded coordinate is the exception, section 7. htslib's own index build refuses an end below its start from 1.10. A line no number leads it skips, with the same message and exit 0, from 1.10 to 1.21, and refuses from 1.22.

**The class it produces.** `gxf_malformed_coordinate_line_refused`: a run over a GTF or GFF3 transcript model, a `--regulatory_gff` file or a `--custom` GFF holding a line of the class. On a kept line that is a start or end that is not a number, a start of 0 or an end below its start; on any line, a start or end no number leads or an end below zero. Ensembl VEP builds without the line or on the span as written, or dies in a Perl method. vep-rs refuses the build or the run at the line, naming it and the value, or refuses the source's index, each reader as tabled. A zero-padded coordinate is the class's extension (section 7).

| vep-rs reader | A kept line of the class | A skipped line whose start or end no number leads, or whose end is below zero |
| --- | --- | --- |
| `vep-cache-builder`: the GFF3 and GTF models, the regulatory bake | refuses the line | refuses the line |
| `--regulatory_gff` at run time | refuses the record | refuses the line |
| `--custom` GFF overlay | refuses no line: reads each column as Perl reads it (section 8), and lists no window whose end the index build fails | a line of another type yields no window |
| `vep-cache-builder`: a CDS line that sets the translation's start or end in no exon of its transcript, [section 9](gxf-cds-line-outside-its-exons.md) | refuses the line | no such line: a CDS line is a kept line |

At run time `--regulatory_gff` refuses such a skipped line as it does a kept one, where Ensembl VEP loses the records after an end of the class in its block (section 6). A `--custom` start of `-5` is matched (section 8).

A skipped-type line the index build parses is matched (section 4). A `--custom` VCF row whose `END` lies below its position is [`tabix_row_with_an_end_below_its_start`](tabix-row-with-an-end-below-its-start.md)'s. No committed file carries a line of the class: it is documented only.

**How often the shape occurs.** Tier 3: the shape is a coordinate column no release file carries. None of the 28,494 feature lines of the committed GTF, promoter and regulatory GFF files carries it, and none of their starts or ends is zero-padded. 26 committed flag sets run a cache built from a GTF or GFF3, a `--custom` GFF or `--regulatory_gff`, and none runs Ensembl VEP over a file with the shape. htslib 1.22.1's `tabix` refuses a column no number leads and an end below its start with exit 1, so a file indexed with it meets the refusal before either engine runs.

## 1. A CDS line whose start is not a number

**In short.** The GRCh38 release-116 GFF3 lines of the gene TMPRSS3 and its transcript `ENST00001113783`, with the start of the transcript's first CDS line in transcript order, the line at `42395324-42395417`, rewritten `abc`, beside the lines of PWP2 and its transcript `ENST00000456705` as released. On the record `tmprss3_codon6`, `21 42395401 G>A`, Ensembl VEP 116.3 reading the file through `--gff` writes `5_prime_UTR_variant` and `HGVSc=ENST00001113783.1:c.-78C>T`, where the released file gives `missense_variant` and `c.17C>T`; vep-rs refuses to build a cache from the file at that CDS line. On `pwp2_codon3` both engines write `missense_variant`.

### 1.1 What a reader needs to know first

Ensembl VEP's `--gff` source does not read the file: it reads the records `tabix` returns for each region from the index built with `tabix -p gff`, so a line the index cannot place is a line the source never sees. The translation of a transcript runs from its first CDS line in transcript order to its last, so losing the first moves the coding start to the second.

### 1.2 The record

The transcript's two CDS lines, the first in transcript order (the higher coordinates, on the reverse strand) with its start rewritten `abc`; the gene, transcript and exon lines, and the PWP2 lines, are as released and omitted here:

```
21	havana_tagene	CDS	42394915	42395084	.	-	2	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1
21	havana_tagene	CDS	abc	42395417	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1
```

The input:

```
21	42395401	tmprss3_codon6	G	A	.	.	.
21	44107424	pwp2_codon3	T	C	.	.	.
```

Both engines run `--symbol --hgvs`. Ensembl VEP reads the file sorted, bgzipped and indexed as `cds_start_not_a_number.gff3.gz` through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 1.3 What each engine writes

The index build, `tabix -p gff` of htslib 1.9, prints the line and exits 0:

```
[E::get_intv] Failed to parse TBX_GENERIC, was wrong -p [type] used?
The offending line was: "21	havana_tagene	CDS	abc	42395417	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1"
```

Ensembl VEP 116.3, exit 0:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	5_prime_UTR_variant	249	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;SYMBOL=TMPRSS3;SOURCE=cds_start_not_a_number.gff3.gz;HGVSc=ENST00001113783.1:c.-78C>T
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;SOURCE=cds_start_not_a_number.gff3.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

vep-rs, building the cache:

```
ERROR: [cache-builder] Failed to parse GFF3: cds_start_not_a_number.gff3 line 6: GFF3 start "abc" is not a number
```

The build exits 1 and writes no cache, so vep-rs annotates nothing from the file. The same start on another line, or the same line elsewhere:

| The line with the start `abc` | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| the gene line | drops the gene, warning `Parent entries with the following IDs were not found or skipped due to invalid types: gene:ENSG00000160183`, drops the transcript as an orphan and writes `intergenic_variant` for `tmprss3_codon6` | refuses the file at line 2 |
| the exon line at `42395942-42396122` | builds the transcript without an exon under its CDS and dies with `Can't call method "seq" on an undefined value at .../Bio/EnsEMBL/Transcript.pm line 828`, exit 255 | refuses the file at line 8 |
| the CDS line, on the later sequence of a two-sequence file | has no file to read: `tabix` prints the message above, then `[E::hts_idx_push] Unsorted positions on sequence #1: 44115691 followed by 1` and `tbx_index_build failed`, exits 1 and writes no index | refuses the build at the line |
| a GTF's first CDS line in transcript order: line 618 of the committed `chr21_genes.gtf.gz`, of `ENST00000398405` | for the record `2736990` at `21:42375772` its `--gtf` source writes `missense_variant` with `CDS_position` 1188, `Protein_position` 396, `Amino_acids` `K/N` and `HGVSc=ENST00000398405.6:c.1188C>T`, where the committed file gives 1282, 428, `P/S` and `c.1282C>T` | refuses the file with `GTF start "abc" is not a number` at line 618 |

The file of this section holds one sequence, so the unparsed line is its first once sorted.

### 1.4 Why Ensembl VEP writes that

`tabix` could not place the line, so no query returns it and `_get_records_by_coords` never meets it. `_add_translation` then takes the remaining CDS line as the transcript's first, the coding start moves 94 bases into the transcript, and `tmprss3_codon6` falls before it, in what the model now calls the 5' untranslated region. Every value on the row follows from that model; no step reads the file's own line.

### 1.5 Why vep-rs's refusal is the correct outcome

The line states a coding segment the model needs and gives it no start; a model that leaves the segment out describes a transcript the file does not. Refusing the file names the line to correct, and the corrected file builds the transcript the file describes; building without the line prints a consequence that is false of the record.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 cds_start_not_a_number.gff3 --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
vep -i records.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 --symbol --hgvs
```

vep-rs's cache builder refuses the file at its sixth line with the error above and writes no cache. Ensembl VEP 116.3 with `--gff cds_start_not_a_number.gff3.gz`, prepared by its documented recipe (`grep -v "#"`, `sort -k1,1 -k4,4n -k5,5n`, `bgzip`, its own `tabix -p gff`), writes `5_prime_UTR_variant`. Reading that index takes it about 10 GB of memory; under 4 GB it is killed. Restoring the start makes vep-rs build the cache and both engines write `missense_variant`.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/gxf_malformed_coordinate_line_refused.toml` holds the file of this section, the records, Ensembl VEP's rows and vep-rs's refusal. On files its steps write it asserts the refusals of sections 2 to 6, the index build's two outcomes for a start no number leads, both engines' rows of section 7 and their listing of section 8's `-5` window. The class covers every field of every row of a run over a file holding a line the class names. No committed dataset reaches it, so the class is documented only.

## 2. A CDS line whose end is below its start, and a line whose start is 0

**In short.** The same lines, with the first CDS line's coordinates swapped, `42395417-42395324`. The 1.9 `tabix` indexes the file silently and Ensembl VEP 116.3 builds on the span, writing `5_prime_UTR_variant` and `HGVSc=ENST00001113783.1:c.-77C>T` for `tmprss3_codon6`, exit 0; vep-rs refuses the build at the line. With the same line's start rewritten `0` instead, `tabix` indexes the file silently and Ensembl VEP dies before a row. Its message is `You must specify a region in the format chr, chr:start or chr:start-end at .../Bio/DB/HTS/Tabix.pm line 104`, exit 255; vep-rs refuses the build at the line.

### 2.1 What a reader needs to know first

A GFF3 or GTF coordinate is one-based and a feature's end is never below its start; `tabix` of htslib 1.9 indexes both shapes without a word, a start of 0 placed at the sequence's first base.

### 2.2 The record

The first CDS line with its coordinates swapped, and the same line with its start rewritten `0`:

```
21	havana_tagene	CDS	42395417	42395324	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1
21	havana_tagene	CDS	0	42395417	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1
```

The input and the flags are those of 1.2, the files `cds_end_below_start.gff3` and `cds_start_zero.gff3`.

### 2.3 What each engine writes

Ensembl VEP 116.3 on the swapped line, exit 0:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	5_prime_UTR_variant	249	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;SYMBOL=TMPRSS3;SOURCE=cds_end_below_start.gff3.gz;HGVSc=ENST00001113783.1:c.-77C>T
```

vep-rs, building the cache from each file:

```
ERROR: [cache-builder] Failed to parse GFF3: cds_end_below_start.gff3 line 6: GFF3 coordinates 42395417-42395324 are not 1-based and ordered
ERROR: [cache-builder] Failed to parse GFF3: cds_start_zero.gff3 line 6: GFF3 coordinates 0-42395417 are not 1-based and ordered
```

On the zero start Ensembl VEP writes no row: the run dies with the message above. The swapped coordinates on the exon line at `42395942-42396122` instead end Ensembl VEP's run in `Bio::EnsEMBL::Feature::new` with `Start (42396122) must be less than or equal to end+1 (42395943)`, exit 255. A zero start on the gene line dies as the CDS line's does; vep-rs refuses each at its line. In a file of two sequences both shapes read the same on either, the zero start's run writing the rows of an earlier sequence's sound records, and no other row, before it dies.

### 2.4 Why Ensembl VEP writes that

The swapped CDS line reaches `_add_translation` with its end below its start, and the coding start the model takes from it leaves the record 77 bases before the coding sequence, where section 1's model left it 78. The zero start reaches `_get_records_by_coords` as the lowest start of the region, and the rescan seeks the index from base 0, a region `Bio::DB::HTS::Tabix` refuses as not of the form `chr:start-end`.

### 2.5 Why vep-rs's refusal is the correct outcome

A pair that is not a one-based ordered span describes no segment, and a model built from it carries a coding start the file never states. The refusal names the line; the corrected file builds the model the file describes.

### 2.6 Try it yourself

The commands of 1.6 on `cds_end_below_start.gff3` and `cds_start_zero.gff3`: vep-rs's cache builder refuses each at its sixth line with the errors above; Ensembl VEP 116.3 writes `5_prime_UTR_variant` on the first and dies on the second.

## 3. A line whose end is not a number, mid-file

**In short.** The same lines, with the end of the exon line at `42395942-42396122` rewritten `xyz`. The 1.9 `tabix` prints its message and exits 0, and a query returns neither that line nor the exon line after it in its block. Ensembl VEP 116.3 builds the transcript without those exons and dies under its CDS, `Can't call method "seq" on an undefined value at .../Bio/EnsEMBL/Transcript.pm line 828`, exit 255; vep-rs refuses the build at the line. A `five_prime_UTR` line with the same end, a type the model skips, kills Ensembl VEP's run the same way, and vep-rs refuses the build at that line too.

### 3.1 What a reader needs to know first

`tabix` reads a region's lines block by block, and the 1.9 reader stops a block at a line whose end it cannot parse, so the lines after it in the block are lost to the query whatever their own columns hold.

### 3.2 The record

The exon line with its end rewritten `xyz`, and the `five_prime_UTR` line, a type neither engine's model keeps, added with the same end:

```
21	havana_tagene	exon	42395942	xyz	.	-	.	Parent=transcript:ENST00001113783;Name=ENSE00004778870;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00004778870;rank=1;version=1
21	havana_tagene	five_prime_UTR	42395942	xyz	.	-	.	Parent=transcript:ENST00001113783
```

The input and the flags are those of 1.2, the files `exon_end_not_a_number.gff3` and `utr_end_not_a_number.gff3`.

### 3.3 What each engine writes

Ensembl VEP 116.3 on either file writes no row and exits 255 with the `Transcript.pm` message above. vep-rs, building the cache from each file:

```
ERROR: [cache-builder] Failed to parse GFF3: exon_end_not_a_number.gff3 line 8: GFF3 end "xyz" is not a number
ERROR: [cache-builder] Failed to parse GFF3: utr_end_not_a_number.gff3 line 8: GFF3 end "xyz" is not a number
```

In a file of two sequences the 1.9 `tabix` exits 0 with the `five_prime_UTR` line on either sequence. On the first Ensembl VEP writes `intergenic_variant` for a record under that sequence's transcript, exit 0; on the later it dies as above. vep-rs refuses the build at the line on either.

### 3.4 Why Ensembl VEP writes that

The query stops at the unparsed end and the exon line after it is never returned, so the transcript reaches `_create_transcript` with a CDS segment no exon covers, and `Transcript.pm` asks that exon for its sequence.

### 3.5 Why vep-rs's refusal is the correct outcome

Either line states an end no reader can place. No htslib parses it, and the release's index build passes it only by losing the lines after it to every query that reaches it. Refusing the file names the line to correct, whatever its type, and the corrected file builds the transcript as released.

### 3.6 Try it yourself

The commands of 1.6 on `exon_end_not_a_number.gff3` and `utr_end_not_a_number.gff3`: vep-rs's cache builder refuses each at its eighth line; Ensembl VEP 116.3 dies on both.

## 4. A line of a type the model skips

**In short.** A `five_prime_UTR` line added to the released lines, its columns rewritten. A start or end no number leads, or an end below zero, refuses vep-rs's build; Ensembl VEP writes the released rows, loses later lines or has no index. Other shapes the index build parses are matched.

| The line's columns | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| start `abc`, the line the file's first | the released file's rows, through the index of section 1 | refuses the build at the line |
| start `abc`, the line on a later sequence | no index to read: `tabix` exits 1 | refuses the build at the line |
| end `xyz` or `-9` | loses the lines after it in its block (section 3) | refuses the build at the line |
| a line of four to eight columns whose start or end no number leads, or whose end is below zero | the index build of the three rows above: `tabix` reads no column past the end | refuses the build at the line, which it skips whatever its type |
| a line that opens with a tab, or with white space and a tab, a type no number leads where a start is read | fails the index build's parse, as the start `abc` does | refuses the build at the line, the record's type kept or skipped |
| start `0`; end below start, not below zero; a start a number leads without being one (`42395942x`, `-5`) | the released file's rows on either sequence of a two-sequence file, the line skipped by its type with `Ignoring 'five_prime_UTR' feature_type` | the released file's rows, the line skipped at its type |
| an empty or blank start, the end a number below every later start of its sequence | the released file's rows on either sequence: `tabix` reads the end's number as the start, in silence | the released file's rows, the line skipped at its type |
| a seqid of any text, the line at the file's end | the released file's rows, the line indexed as a sequence of its own | the released file's rows, the line skipped at its type |
| a line whose columns NUL bytes end, numbers in its fourth and fifth: among the records, led by `>` before a placed line, or behind a placed line that begins with `>` | the released file's rows: `tabix` places the line, a column ending at a NUL byte to it, and the source reads one column in it, `Ignoring '' feature_type` | matched: the released file's rows, the line skipped as any placed line of fewer than nine columns |
| one NUL byte inside the sequence name, the source or the type of a line | no index to read, `tabix` exits 1; on the file's first line, the index of section 1 | refuses the build at the line: the build reads a start one column on, the type or its part behind the NUL byte |

A line the model keeps by its identifier is no skipped line; a `--custom` or `--regulatory_gff` source is read through its index:

| The line | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| an `unconfirmed_transcript` line whose `ID` opens with `gene:` or `transcript:`, the cache built with the release's protein FASTA or core dumps | skips it by its type | keeps it, as every line so named whatever its type, and holds it to the whole rule of section 2: a start of 0 refuses it |
| the same line, the cache built from the genome alone | skips it by its type | skips it by its type: only a start or end no number leads, or an end below zero, refuses it |
| a `--regulatory_gff` line of a type the reader skips, an `exon` line, whose end is `xyz` or `-9` | loses the records after it in its block, exit 0 (section 6) | refuses the file at run time, `the end of exon record <ID> of chromosome 21 is not a whole number: xyz`, and at the bake |
| a `--regulatory_gff` line that opens with a tab, or with white space and a tab, indexed as written behind a record at its sequence's first base | writes the file's other rows, exit 0: `tabix` fails the line's parse and files it there, and no query returns it | refuses the file at run time, `the start of a record of type Ensembl of chromosome 21 is not a whole number: enhancer`, the source column read as the type, and at the bake |
| a `--regulatory_gff` record whose start or end a NUL byte ends | builds the feature on the digits, with Perl's `isn't numeric` warning | refuses the record, as for `100x` (section 8) |
| a `--regulatory_gff` record whose start or end holds a NUL byte among its digits | writes no row for the record: `tabix` ends the column at the byte and files the record where no query returns it | refuses the record |
| a `--regulatory_gff` line of any type with a NUL byte in its sequence name or its source, or inside its type, indexed as written behind a record at its sequence's first base | writes the file's other rows and none for the line, exit 0: `tabix` ends a column at the byte, reads the type or its tail where a start stands and fails the line, and no query returns it | refuses the file at run time, `the start of enhancer record <ID> of chromosome 21 is not a whole number: enhancer`, and at the bake: both read the start where `tabix` reads one |
| a `--regulatory_gff` line of any type whose start or end a NUL byte leads | as for a start `abc` or an end `xyz`: `tabix` reads an empty column, and `strtol` stops at the byte | refuses the file at run time and at the bake, naming the empty column |
| a `--regulatory_gff` line whose every column a NUL byte ends, of any type, numbers in its fourth and fifth | skips it, `Ignoring '' feature type`, a record it would state lost | matched: skipped at run time and at the bake; with a start no number leads, refused by both |
| a `--regulatory_gff` record of nine columns whose sequence name ends in three NUL bytes before a source that is a number | builds the feature where a query reaches the base `tabix` files the line at, the source's number: `tabix` reads the name up to its first NUL byte, and the source reads the record the index returns whatever its first column (`RegFeat.pm`, lines 311-315) | passes over the line at run time and skips it at the bake: the sequence name as written is no chromosome's, as on the line whose every column a NUL byte ends |
| a `--custom` GFF promoter window with a start of 0 or an end below its start, not below zero | reads the index as `tabix` wrote it | matched: reads the same index |
| a `--custom` GFF window, BED record or VCF row whose start is `-5` | lists it from the sequence's first base to its end, under the number Perl reads | matched: lists the same record under the same number |
| a `--custom` window inside the file whose end no number leads, or a minus sign leads | lists no such window, and loses the window after it for a record that listed one before it | lists no such window, and every window after it; at the file's head the line leaves the index of section 5 |

**A line with no start column.** The release's `tabix` fails a line of fewer than four columns as it fails a start no number leads. The GFF3 specification 1.26 has parsers ignore a blank line and ends the annotation at `##FASTA` or a line beginning with `>`.

| The line | The release's `tabix` | `vep-cache-builder` |
| --- | --- | --- |
| a blank line | sorted to the file's head by the recipe, it leaves the index of section 1 | skips it |
| a `##FASTA` directive | a comment: `tabix` skips it and the recipe's `grep` drops it, so through the recipe the records behind a bare one are indexed and Ensembl VEP reads them; indexed as written, a query of the sequence behind it returns nothing | a comment: skips it and reads the records behind it |
| the first line beginning with `>` that a sequence follows: the next line that is neither blank nor a comment exists and has fewer than four columns, or a start or end `tabix` fails; a column ends at a tab or at a NUL byte, here and in every row, as `tabix` ends one; a line longer than 16777216 bytes is a sequence's only where no tab stands in it, its 16777217th byte aside: one in which a tab stands elsewhere opens no section and is refused for its length, as on any line of the file | fails the `>` line's parse or, where it can parse it, places it as a sequence of its own; indexed as written it then fails the line behind it, and by that line the build has exited 1 where the line has fewer than four columns, as a sequence line has; section 4's table gives the outcome of a line of four or more | opens the sequence section: reads nothing from there on, and logs the line |
| a line beginning with `>` that a line of four or more columns follows, its start and end ones `tabix` parses, or that nothing follows | where it can parse it, places it as a sequence of its own, the first table's last row, and indexes the lines behind it; fails its parse otherwise: the build exits 1 or, on the file's first data line, leaves the index of section 1 | reads it as any other line: a record of its type at nine columns, skipped at four to eight where `tabix` parses its start and end, and refused at fewer than four, a record cut short |
| a record inside the sequence section, nine or more columns whose fourth and fifth are whole numbers, neither a comment nor led by `>` | indexed as written, has exited 1 before the record where the line behind the section's `>` line has fewer than four columns, as a sequence line has; section 4's table gives the outcome of a line of four or more | refuses it, annotation the section has ended: `<file> line <n>: GFF3 annotation after the sequence section, which opens at line <m>` |
| any other line of fewer than four columns | fails its parse | refuses it, a record cut short: `<file> line <n>: GFF3 line ends at column <k>, before its start` (`GTF` in a GTF) |
| a comment behind a byte-order mark | no comment to `tabix`: the recipe's `grep` drops it; indexed as written it fails its parse and, on the file's first line, leaves the index of section 1, a query of the first sequence returning nothing where more comments follow | skips it: the mark is taken off the start of every line before the line is read, in the protein FASTA as in a GFF3, a GTF and a regulatory file, so a file with marks builds the cache of the same file without them |

A tab-delimited file whose every data line `tabix` places opens no section.

**A kept line `tabix` files apart.** A byte-order mark before a file's first record names a sequence of its own to `tabix`: Ensembl VEP loses the line, which vep-rs reads with the mark off. A kept line led by spaces and no tab is lost likewise; vep-rs reads it trimmed. So is a kept record whose seqid begins with `>`, which vep-rs reads unless a sequence follows.

## 5. A `--custom` GFF promoter window whose start is not a number

**In short.** A promoter file of three windows, the second with its start rewritten `abc`, read as a `--custom` source in `gencode_promoter` mode beside the released transcript lines. The 1.9 `tabix` prints its message and writes an index naming 67.7 million sequences, 541 MB uncompressed for a three-line file; Ensembl VEP 116.3 reads it and lists the two sound windows, `PW1` on `tmprss3_codon6` and `PW2` on `pwp2_codon3`. vep-rs refuses to open the source, since the index decompresses past its 268435456-byte bound, and writes no row.

### 5.1 What a reader needs to know first

A `--custom` GFF is read through its index on both engines; neither reads the file line by line. The window with the start `abc` sorts first, so the 1.9 `tabix` files it under a sequence id never set, and the index is written for 67.7 million sequences: a three-line file gets an index larger than most genomes'.

### 5.2 The record

The promoter file, its second window with the start `abc`:

```
21	.	promoter_window	42395390	42395410	.	+	.	ID=PW1;gene_id=ENSG00000160183
21	.	promoter_window	abc	44107430	.	+	.	ID=PWBAD
21	.	promoter_window	44107420	44107430	.	+	.	ID=PW2
```

The input is the two records of 1.2. Both engines read the released transcript lines (Ensembl VEP through `--gff`, vep-rs from a cache built from them) and the promoter file as `--custom file=promoters.gff.gz,format=gff,gff_type=gencode_promoter`.

### 5.3 What each engine writes

Ensembl VEP 116.3, exit 0:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SOURCE=released.gff3.gz;GENCODE_PROMOTER=21:42395390-42395410;GENCODE_PROMOTER_associated_gene=ENSG00000160183;GENCODE_PROMOTER_feature_id=PW1;GENCODE_PROMOTER_type=promoter_window
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SOURCE=released.gff3.gz;GENCODE_PROMOTER=21:44107420-44107430;GENCODE_PROMOTER_associated_gene=.;GENCODE_PROMOTER_feature_id=PW2;GENCODE_PROMOTER_type=promoter_window
```

vep-rs, exit 1 and no row:

```
ERROR: [custom] cannot open the --custom source promoters.gff.gz: tabix index promoters.gff.gz.tbi is malformed: it decompresses to more than 268435456 bytes
```

### 5.4 Why Ensembl VEP writes that

`tabix` filed the window it could not parse under a sequence id it never set, so no query of chromosome 21 meets it, and the two sound windows are returned as from a sound file. Nothing in the route reports the window or the size of the index.

### 5.5 Why vep-rs's refusal is the correct outcome

An index of 541 MB for three lines describes no file; it is the trace of a line the indexer could not read. Refusing it names the index to rebuild, and the file corrected and re-indexed is read by both engines with the same three windows; reading it lists two windows and loses the third without a word.

### 5.6 Try it yourself

```
vep -i records.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --custom file=promoters.gff.gz,format=gff,gff_type=gencode_promoter
```

With the file sorted, bgzipped and indexed by the 1.9 `tabix`, vep-rs refuses the source with the error above. Ensembl VEP 116.3 with the same `--custom` string and `--gff released.gff3.gz` lists `PW1` and `PW2`. Restoring the start makes both engines list the three windows; htslib 1.22.1's `tabix` writes no index for the file as it stands, so neither engine opens it. Nor does the 1.9 `tabix` when the window sits on the later sequence of a two-sequence file: it exits 1.

## 6. The same shapes through `--regulatory_gff`

**In short.** A regulatory file of four records on chromosome 21, the second an enhancer `ENSR21_BZTD` at `5013566-5014200`, read through `--regulatory_gff` over the chromosome 21 cache, with records `in_first` at 21:5011400 inside `ENSR21_85ZF5`, `in_lost` at 21:5014000 inside `ENSR21_BZTD` and `in_last` at 21:5015600 inside `ENSR21_93ZQB`. With the enhancer's end rewritten `xyz`, Ensembl VEP 116.3 loses it and the enhancer after it in its block, writing `intergenic_variant` alone for `in_lost` and `in_last`, exit 0; vep-rs refuses the file, naming the record.

With its start rewritten `abc`, Ensembl VEP loses the enhancer alone; vep-rs refuses the index. With its coordinates swapped, Ensembl VEP writes a `RegulatoryFeature` row whose consequence is `intergenic_variant`; vep-rs refuses the file, naming the record. With its start rewritten `0`, Ensembl VEP gives the enhancer to `in_first` as well, the feature running from the chromosome's first base; vep-rs refuses the file. In a file of two sequences each shape reads the same on either sequence, but for the start `abc`, which on the later one fails the index build.

### 6.1 What a reader needs to know first

`--regulatory_gff` reads a regulatory GFF3 through its index on both engines, as `--gff` does, and vep-rs reads every record of a chromosome's index chunks into its tree before the first variant, refusing a record whose coordinates are not a one-based ordered pair of numbers. The builder's `--regulatory-gff` bake reads the plain file line by line and refuses the same shapes at their line. Both refuse a start or end no number leads, or an end below zero, on a line of any type.

### 6.2 The record

The enhancer line in its four shapes, the other three lines as written:

```
21	Ensembl	enhancer	5011320	5011550	.	.	.	ID=ENSR21_85ZF5;color=#faca00
21	Ensembl	enhancer	5013566	xyz	.	.	.	ID=ENSR21_BZTD;color=#faca00
21	Ensembl	CTCF_binding_site	5014322	5014343	.	-	.	ID=ENSR21_536SRH;color=#40e0d0
21	Ensembl	enhancer	5015556	5015987	.	.	.	ID=ENSR21_93ZQB;color=#faca00
```

The input:

```
21	5011400	in_first	A	C	.	.	.
21	5014000	in_lost	A	C	.	.	.
21	5015600	in_last	A	C	.	.	.
```

Both engines run `--regulatory_gff file=<file>` over the GRCh38 release-116 cache of chromosome 21 with its FASTA.

### 6.3 What each engine writes

Ensembl VEP 116.3 on the end `xyz`, exit 0:

```
in_first	21:5011400	C	-	ENSR21_85ZF5	RegulatoryFeature	regulatory_region_variant	-	-	-	-	-	-	IMPACT=MODIFIER;BIOTYPE=enhancer
in_first	21:5011400	C	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
in_lost	21:5014000	C	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
in_last	21:5015600	C	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

vep-rs on the same file, exit 1 and no row:

```
ERROR: [cache] reg_end_not_a_number.gff3.gz: the end of enhancer record ENSR21_BZTD of chromosome 21 is not a whole number: xyz
```

On the start `abc`, Ensembl VEP writes the released rows for `in_first` and `in_last` and `intergenic_variant` alone for `in_lost`. On the swapped coordinates, `5014200-5013566`, Ensembl VEP writes for `in_lost`:

```
in_lost	21:5014000	C	-	ENSR21_BZTD	RegulatoryFeature	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;BIOTYPE=enhancer
```

and the released rows for the other two records. On the start `0`, Ensembl VEP writes the released rows and one more for `in_first`:

```
in_first	21:5011400	C	-	ENSR21_BZTD	RegulatoryFeature	regulatory_region_variant	-	-	-	-	-	-	IMPACT=MODIFIER;BIOTYPE=enhancer
```

The enhancer runs from the chromosome's first base there. vep-rs refuses the three files, the first at its index:

```
ERROR: [cli] --regulatory_gff file=reg_start_not_a_number.gff3.gz: tabix index reg_start_not_a_number.gff3.gz.tbi is malformed: it decompresses to more than 268435456 bytes
ERROR: [cache] reg_end_below_start.gff3.gz: the coordinates 5014200-5013566 of enhancer record ENSR21_BZTD of chromosome 21 are not 1-based and ordered
ERROR: [cache] reg_start_zero.gff3.gz: the coordinates 0-5014200 of enhancer record ENSR21_BZTD of chromosome 21 are not 1-based and ordered
```

The builder's bake of each of the four files refuses it, the line counted in the sorted file:

```
ERROR: [cache-builder] Failed to build the regulatory tree: Failed to read reg_end_not_a_number.gff3.gz: the end of line 2 of reg_end_not_a_number.gff3.gz is not a whole number: xyz: invalid digit found in string
ERROR: [cache-builder] Failed to build the regulatory tree: Failed to read reg_start_not_a_number.gff3.gz: the start of line 1 of reg_start_not_a_number.gff3.gz is not a whole number: abc: invalid digit found in string
ERROR: [cache-builder] Failed to build the regulatory tree: Failed to read reg_end_below_start.gff3.gz: GFF coordinates 5014200-5013566 at line 2 are not 1-based and ordered
ERROR: [cache-builder] Failed to build the regulatory tree: Failed to read reg_start_zero.gff3.gz: GFF coordinates 0-5014200 at line 1 are not 1-based and ordered
```

### 6.4 Why Ensembl VEP writes that

`_get_regfeats_by_coords` reads the records `tabix` returns for the region, so a line the index dropped is no feature, and a query that stops at an unparsed end returns no later line of the block. A record with its end below its start is built as a feature whose end precedes its start; the variant inside the written span overlaps the feature by the row's test and misses it by the consequence's, hence a regulatory row with an intergenic term. A start of 0 is built as a feature from the chromosome's first base.

### 6.5 Why vep-rs's refusal is the correct outcome

The file's record states no span the reader can place; refusing the run names the record to correct, where reading on loses the record and, through the block, its neighbours. A record whose end is below its start, or whose start is 0 in a one-based format, is no feature, and a row built on it describes a span the file never stated.

### 6.6 Try it yourself

```
vep -i records.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --regulatory_gff file=reg_end_not_a_number.gff3.gz
```

vep-rs refuses the file with the error above, and the swapped and zero-start files with theirs. Ensembl VEP 116.3 with the same `--regulatory_gff` string over its chromosome 21 cache writes `intergenic_variant` alone for `in_lost` and `in_last`. Restoring the end makes both engines write the three regulatory rows.

## 7. A zero-padded coordinate: the index reader's octal

**In short.** The released lines with the end of the transcript's first CDS line in transcript order, `42395417`, written `0042395417`. The 1.9 `tabix` reads a column in base 0, a leading `0` as octal up to the first digit above 7, so the end reads 275. The line is indexed silently with an end below its start, and no query over the CDS returns it. Ensembl VEP 116.3 builds the transcript without the line and writes `5_prime_UTR_variant` and `HGVSc=ENST00001113783.1:c.-78C>T` for `tmprss3_codon6`, exit 0, as in section 1. vep-rs reads the column in decimal, builds the released model and writes `missense_variant` and `c.17C>T`.

### 7.1 What a reader needs to know first

`strtol` in base 0 reads a leading `0` as an octal prefix and stops at the first character that is no digit of the base. Perl and vep-rs read the column in decimal. The index and the model so disagree on where a zero-padded line lies, and the index decides which lines the model is given.

### 7.2 The record

The transcript's first CDS line in transcript order with its end zero-padded; every other line is as released:

```
21	havana_tagene	CDS	42395324	0042395417	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1
```

The input and the flags are those of 1.2, the file `cds_end_zero_padded.gff3`.

### 7.3 What each engine writes

`tabix -p gff` of htslib 1.9 prints nothing and exits 0. Ensembl VEP 116.3, exit 0:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	5_prime_UTR_variant	249	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;SYMBOL=TMPRSS3;SOURCE=cds_end_zero_padded.gff3.gz;HGVSc=ENST00001113783.1:c.-78C>T
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;SOURCE=cds_end_zero_padded.gff3.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

vep-rs builds the cache and writes:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

In a file of two sequences the zero-padded end loses its line on either sequence. A zero-padded start fails the index build wherever a line before it on its sequence starts above the octal reading: `tabix` prints `Unsorted positions on sequence #2: 42394915 followed by 275` and exits 1. On the lowest line of its sequence the index is built, the line is returned, and both engines write the released rows.

### 7.4 Why Ensembl VEP writes that

The index holds the line from 42395324 to 275, an interval no query over the transcript overlaps, so `_add_translation` starts the translation at the next CDS line, as in section 1. Perl reads `0042395417` as 42395417, but is never given the line.

### 7.5 Why vep-rs's rows are the correct ones

A GFF3 or GTF coordinate is a decimal integer, and Ensembl VEP's own model reads it so. vep-rs builds the transcript the file describes, and refuses nothing: the line states a span in the format's own terms.

### 7.6 Try it yourself

The commands of 1.6 on `cds_end_zero_padded.gff3`: vep-rs builds the cache and writes `missense_variant` for both records. Ensembl VEP 116.3, the file sorted, bgzipped and indexed with its own `tabix -p gff`, writes `5_prime_UTR_variant` for `tmprss3_codon6`. htslib 1.10 and later refuse the file at index time: `Invalid record on sequence #1: end 275 < begin 42395324`.

### 7.7 How the comparator treats this extension

The extension is the same class. The reproduction's last step builds the file of 7.2, runs both engines and holds each to its rows above. No committed start or end is zero-padded, so the extension is documented only.

## 8. How each reader reads a column's spelling

A start or end a number leads without being a plain decimal number is read differently along the route. vep-rs's builder and its `--regulatory_gff` reader take a whole decimal number, a `+` allowed before it; its `--custom` sources take Perl's number, an `exact` VCF match the text. Each cell was measured on the pinned image over two-sequence files, the line on the first sequence and on a later one.

| Spelling | Ensembl VEP 116.3 | vep-rs | Class |
| --- | --- | --- | --- |
| `-5` as a start, a minus sign leading | index reader: -5, the line filed from the sequence's first base; Perl: -5, so a kept transcript line ends the run at `Tabix.pm` line 104, and a regulatory record or a promoter window runs from -5 | builder: refuses a kept line, leaves a skipped one unread; `--regulatory_gff`: refuses the record; `--custom`: -5, as Perl, the window listed from the first base | this entry; matched on `--custom` |
| `+7`, a plus sign leads | index reader: 7; Perl: 7 | 7 in every reader | matched |
| ` 9`, a blank leads | index reader: 9; Perl: 9 | builder: refuses a kept line, leaves a skipped one unread; `--regulatory_gff`: refuses the record; `--custom`: 9 | this entry |
| `100x`, text after the digits | index reader: 100; Perl: 100, with an `isn't numeric` warning in the warnings file | builder: refuses a kept line, leaves a skipped one unread; `--regulatory_gff`: refuses the record; `--custom`: 100, and no warning: Perl's `isn't numeric` line is its own diagnostic, not a row | this entry |
| `0x1A`, hexadecimal | index reader: 26; Perl: 0 | builder: refuses a kept line, leaves a skipped one unread; `--regulatory_gff`: refuses the record; `--custom`: 0, as Perl | this entry |
| `0100`, zero-padded | index reader: 64, octal; Perl: 100 | 100 in every reader | the octal extension |
| `1e2`, an exponent (`+1e2`, `1.5e2`, `2000e-1`) | index reader: the integer that leads, 1 (2000 for `2000e-1`); Perl: the number spelled, 100 (150, 200), so a promoter window runs from 100 for the records the index returns it to | builder and `--regulatory_gff`: as for `100x`; `--custom`: 100, as Perl | this entry; matched on `--custom` |
| `100.5`, a fraction | index reader: 100; Perl: 100.5, printed as 100: the walk to the record, `within`, `surrounding`, `distance` and the overlap percentage read the fraction, the last overlap test the whole number | builder and `--regulatory_gff`: as for `100x`; `--custom`: 100.5, as Perl | this entry; matched on `--custom` |
| `-9223372036854775809`, below the 64-bit range (`-1e19`) | index reader: the first base; Perl: a float, printed as -9223372036854775808 | builder and `--regulatory_gff`: as for `-5`; `--custom`: -9223372036854775808, as Perl | this entry; matched on `--custom` |
| `1e19` as a start, past the 64-bit range | index reader: 1; Perl: a number past every position, so the window is listed for no record | builder and `--regulatory_gff`: as for `100x`; `--custom`: listed for no record | this entry; matched on `--custom` |
| `-2147483649`, a start below -2^31 | index reader: the first base; Perl: the number written, but the compiled `overlap` the image swaps in (`Utils/VariationEffect.pm`, lines 73-74; ensembl-xs 2.3.2, `Bio-EnsEMBL.xs`, lines 468-472) takes 32-bit integers, so a window is listed from the low 32 bits of its start | `--custom`: the number written, the window listed from the first base | an Ensembl VEP defect, not replicated |
| `2147483648`, an end from 2^31 up (`9223372036854775808`, `1e19`) | index reader: an end cut to 32 bits, below the window's start (1 for `1e19`); Perl: the number written, of which the compiled `overlap` of the row above takes the low 32 bits, so the window is listed for no record | `--custom`: the number written, held at 9223372036854775807 from 2^63 up; the window is listed for the records the reader returns it to | an Ensembl VEP defect, not replicated |

The recipe's `sort -k4,4n` reads `+7` and `0x1A` as 0, so a start so spelled sorts first on its sequence, and the index build fails with `Unsorted positions` when a later line starts below the number the index reader reads. A row classed `this entry` is a column Ensembl VEP reads as a number and vep-rs refuses on a kept line.
