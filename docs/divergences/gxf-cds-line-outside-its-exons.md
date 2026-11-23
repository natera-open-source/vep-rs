# A GFF3 or GTF CDS line outside its transcript's exons

On three files whose CDS line puts the translation's start in no exon, Ensembl VEP 116.3 writes `5_prime_UTR_variant` for a base the line covers or moves a coding row 28 codons; vep-rs refuses to build a cache from each file at the line. This page holds one section of the divergence [`gxf_malformed_coordinate_line_refused`](gxf-malformed-coordinate-line-refused.md). That entry's page states the defect, the class and the sections before this one, and its row of the index covers this section.

## 9. A CDS line that sets the translation's start or end in no exon

**In short.** The released lines of section 1 with one CDS line moved off its exon, three ways, each putting the translation's start in no exon. The first CDS line of PWP2's `ENST00000456705` starts 13 bases before its exon. The first CDS line in transcript order of TMPRSS3's `ENST00001113783` ends 32 bases past its exon. That same line lies in the intron, in no exon. `tabix` indexes each file in silence and Ensembl VEP 116.3 reads it without a message, exit 0: it writes `5_prime_UTR_variant` for a base the moved line covers, or moves a coding row 28 codons. vep-rs refuses to build a cache from each file at the line.

### 9.1 What a reader needs to know first

A CDS line states the coding part of one exon, so its span lies inside that exon's. Ensembl VEP's `--gff` and `--gtf` sources do not hold a line to that. They attach each CDS line to the exon that holds the line's start coordinate, and take the translation's two ends from the first and the last attached line.

vep-rs holds the translation to the exons at its two ends. Its start is the 5' end of the first CDS line in transcript order and its end the 3' end of the last: the lowest and the highest coding base on the forward strand, the highest and the lowest on the reverse. Each lies inside an exon of the transcript, or the build is refused.

### 9.2 The record

The changed CDS line of each file above the exon it left; the other lines are those of section 1 as released:

```
21	havana	CDS	44107360	44107434	.	+	0	ID=CDS:ENSP00000411749;Parent=transcript:ENST00000456705;protein_id=ENSP00000411749;version=1
21	havana	exon	44107373	44107434	.	+	.	Parent=transcript:ENST00000456705;Name=ENSE00001605720;constitutive=0;ensembl_end_phase=0;ensembl_phase=-1;exon_id=ENSE00001605720;rank=1;version=1
```

```
21	havana_tagene	CDS	42395324	42395500	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1
21	havana_tagene	exon	42395324	42395468	.	-	.	Parent=transcript:ENST00001113783;Name=ENSE00003591562;constitutive=0;ensembl_end_phase=1;ensembl_phase=-1;exon_id=ENSE00003591562;rank=2;version=1
```

The third file's line lies between that exon and the one above it:

```
21	havana_tagene	CDS	42395500	42395900	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1
21	havana_tagene	exon	42395942	42396122	.	-	.	Parent=transcript:ENST00001113783;Name=ENSE00004778870;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00004778870;rank=1;version=1
```

The input and the flags are those of section 1.2; the files are `cds_starts_before_its_exon.gff3`, `cds_ends_past_its_exon.gff3` and `cds_in_no_exon.gff3`.

### 9.3 What each engine writes

Ensembl VEP 116.3 on the line that starts before its exon, the forward strand, exit 0:

```
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	5_prime_UTR_variant	52	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SYMBOL=PWP2;SOURCE=cds_starts_before_its_exon.gff3.gz;HGVSc=ENST00000456705.1:c.-11T>C
```

On the line that ends past its exon, the reverse strand, exit 0:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	missense_variant	249	100	34	R/C	Cgc/Tgc	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;SOURCE=cds_ends_past_its_exon.gff3.gz;HGVSc=ENST00001113783.1:c.100C>T;HGVSp=ENSP00000783588.1:p.Arg34Cys
```

On the line in no exon, exit 0:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	5_prime_UTR_variant	249	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;SYMBOL=TMPRSS3;SOURCE=cds_in_no_exon.gff3.gz;HGVSc=ENST00001113783.1:c.-78C>T
```

The released lines give `missense_variant` on both records, `c.8T>C` and `c.17C>T`. vep-rs, building the cache from each file:

```
ERROR: [cache-builder] Failed to build cache: cds_starts_before_its_exon.gff3 line 12: CDS 44107360-44107434 of transcript ENST00000456705 starts the translation at 44107360, in no exon of the transcript
ERROR: [cache-builder] Failed to build cache: cds_ends_past_its_exon.gff3 line 6: CDS 42395324-42395500 of transcript ENST00001113783 starts the translation at 42395500, in no exon of the transcript
ERROR: [cache-builder] Failed to build cache: cds_in_no_exon.gff3 line 6: CDS 42395500-42395900 of transcript ENST00001113783 starts the translation at 42395900, in no exon of the transcript
```

Each build exits 1 and writes no cache. Further shapes, each measured as a GFF3 and as a GTF, with equal rows unless the row says otherwise; the last two leave both ends of the translation inside an exon:

| The line | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| the last CDS line in transcript order of the reverse-strand transcript, starting before its exon | leaves the line out, so a coding base of it is a `3_prime_UTR_variant` | refuses the line: the translation ends in no exon |
| the last CDS line of the forward-strand transcript, ending past its exon | writes the released rows: the end it takes passes the exon, where the spliced sequence stops | refuses the line: the translation ends in no exon |
| the reverse-strand CDS line of this section ending 232 bases past its exon, more than the exon before it holds | prints `substr outside of string` and writes `intergenic_variant` with `c.300C>T` for `tmprss3_codon6` | refuses the line |
| a transcript's only CDS line, in no exon | writes `intergenic_variant` with `n.52T>C` for `pwp2_codon3` | refuses the line |
| exon lines: a GFF3's last two coding exons with their `rank` values exchanged, against position | places each exon at its `rank` (`BaseGXF.pm`, line 574), exit 0 without a message: on the forward strand a base of the last exon is `missense_variant`, `c.380C>T`, where the file ranked by position gives `synonymous_variant`, `c.516C>T`; on the reverse strand, in a transcript of three coding exons, 6 of 8 rows change and the coding bases of one exon read `intron_variant`. Neither engine reads a GTF's `exon_number`: both write the rows of the position order | refuses the GFF3: `<file> line <n>: CDS <start>-<end> of transcript <id> ends the translation at base <offset> of the <length>-base exon <start>-<end>, the last in transcript order that the coding sequence overlaps` |
| exon lines: one more exon nested in the last coding exon, ranked after it in a GFF3 | as a GFF3, places it at its `rank`: on the forward strand a base in both exons is `coding_sequence_variant` at cDNA `560-605`, and on the reverse strand the run dies, `First part of exon maps into gap`, exit 2; as a GTF, warns `Failed to add exon to transcript <id>` and writes `intergenic_variant` for the transcript's records, exit 0 | refuses the model with the same line where the nested exon comes last in transcript order: by its `rank` on either strand, by position on the forward strand |
| exon ranks: a GFF3's three exons ranked 1, 3, 2 against position, every CDS line in an exon and the translation ending inside the exon ranked last | places each exon at its `rank`, exit 0 without a message, on either strand: the bases of the exon ranked last read `intron_variant` beside their coding or UTR term, and a base of the exon ranked second is a `3_prime_UTR_variant` whose HGVSc is a coding position. As a GTF both engines write the rows of the position order | refuses the GFF3: `<file>: exon <start>-<end> of transcript <id>, ranked <n>, starts before exon <start>-<end>, ranked <m>: the exon ranks run against position` (`starts after` on the reverse strand), a transcript of no CDS line included |
| one more CDS line in an intron beside the released ones | writes the released rows: the line is given to no exon and sets neither end of the translation | matched: builds the cache of the released lines |
| one CDS line across the first intron in place of the two it spans | writes the released rows: the line's start, the released one, sets the translation's start, and its end is not read | matched: builds the cache of the released lines |

In the two rows of exon lines every CDS line lies inside an exon, and the exon the translation's end is counted in does not hold that end. No build writes a translation offset outside the exon it is counted in, so a further malformed shape either builds offsets inside their exons or is refused. In the row of exon ranks both offsets lie inside their exons, and the order itself is refused: Ensembl VEP reads the introns between neighbours in rank order, so an exon inside one reads as intron.

### 9.4 Why Ensembl VEP writes that

`_create_transcript` files a transcript's CDS records under their start coordinates and gives each exon the records whose start it holds (`BaseGXF.pm`, lines 511 and 544). A record whose start no exon holds is given to none and takes no part in the translation, which runs from the next record. On the first file the coding start moves to exon 2 and `pwp2_codon3` falls 11 bases before it; on the third the lost line is the one section 1 loses.

`_add_translation` takes the translation's start as the distance from its exon's edge to the start record's own (lines 740 to 745) and does not test the sign. On the second file that is the exon's end less the record's end plus one, -31. The coding start then lands 32 bases above the exon in the spliced transcript, on the last bases of the exon before it, and `tmprss3_codon6` is base 100 of a coding sequence the file never states.

### 9.5 Why vep-rs's refusal is the correct outcome

A translation that starts or ends in no exon states a coding base the transcript does not carry, and no reading of it gives the transcript the file describes. Leaving the line out drops a coding segment; taking its end as written moves the reading frame onto bases of another exon. Where Ensembl VEP still writes the released rows, the forward-strand line ending past its last exon, the spliced sequence happens to stop at the exon, and the file's coding end is a base the transcript lacks all the same. The refusal names the line, the transcript, the line's span and the base, and the corrected file builds the transcript the file describes.

A CDS line that leaves its exon between the translation's two ends is not refused, the last two rows of the table. The translation the file states starts and ends on bases the transcript carries, vep-rs builds the cache the released lines build, and both engines write the released rows.

### 9.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 cds_ends_past_its_exon.gff3 --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
```

vep-rs's cache builder refuses the file at its sixth line with the error above and writes no cache. Ensembl VEP 116.3 with `--gff cds_ends_past_its_exon.gff3.gz`, the file sorted, bgzipped and indexed as in section 1.6, writes `c.100C>T`. Restoring the line's end makes vep-rs build the cache and both engines write `c.17C>T`.

### 9.7 How the comparator treats this section

The entry's reproduction writes the three files of this section in a step of its own. It holds `tabix` to exit 0 in silence on each, and asserts Ensembl VEP's rows above and vep-rs's three refusals. No committed GTF or GFF3 model holds a translation that starts or ends in no exon, so the section is documented only.
