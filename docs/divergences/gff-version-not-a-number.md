# A GFF3 line whose `version` is not a number

**The defect in one sentence.** Ensembl VEP 116.2's `--gff` source takes the text of a `version` attribute as the version of the gene, the transcript or the translation without reading it as a number. A value such as `v2` is printed after the identifier wherever a version is printed: `ENSG00000160183.v2` in `Gene`, `ENST00001113783.v2` in `Feature` and in the HGVSc reference name, `ENSP00000783588.v2` in `ENSP` and in the HGVSp reference name. A value written `01` is printed as `01`. vep-rs refuses to build a cache from a file whose `version` is not a whole number, naming the line and the text, and prints `01` as `1`.

**Why this is a defect.** A version is a whole number: every Ensembl file that carries `version` writes one, a versioned identifier such as `ENST00001113783.1` is read everywhere as an accession and a release number, and an HGVS reference name needs the same form to resolve. `ENST00001113783.v2` names no release of the transcript, and `ENST00001113783.01` and `ENST00001113783.1` are two spellings of one. Ensembl VEP applies no test to the text: the only test on the way to the output is whether the identifier already ends in a dot and digits, which `v2` and `01` both pass through.

The HGVSc reference name takes any defined text, so a transcript line's `transcript_version` of `0`, or an empty one, names `ENST00001113783.0` or `ENST00001113783.` there, while `Feature` omits a value that is not true to Perl.

**What vep-rs does instead.** vep-rs holds every version as a number. A `version` whose text is not a whole number refuses the cache build with the error `<file> line <n>: version "<text>" is not a whole number`, on a gene, transcript or CDS line alike, so no cache is written from the file. A `transcript_version` beside such a `version` is read by neither engine. A whole number above 4294967295, the largest the cache holds, is refused the same way, naming the limit.

A whole number written with leading zeros is that number, so `01` prints as `1`. Every other field is read as Ensembl VEP reads it.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, line 481, `_create_transcript`: the transcript's version is the transcript line's `version`, else its `transcript_version`, as written. Line 639, `_add_identifiers`: the gene's version is the gene line's `version` as written. Line 727, `_add_translation`: the translation's version is the first ordered CDS record's `version`, `|| 1`, so any text but an empty or `0` value is kept as written.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1465: under `--transcript_version`, `Feature` gains `.` and the transcript's version unless the identifier already ends in `.` and digits. Line 1469: under `--gene_version`, `Gene` gains the gene's version the same way. Lines 1560-1565: under `--protein_version`, `ENSP` gains the translation's version the same way.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1499-1500, `hgvs_transcript`: the HGVSc reference name gains `.` and the transcript's version under the same test, whatever the flags. Lines 1759-1760, `hgvs_protein`: the HGVSp reference name gains the translation's version the same way.

**What correcting it would change.** With each attribute read as a whole number, as `_add_translation` already reads an absent one as 1, the version printed would be a number on every row, as vep-rs's is. The rows on files whose `version` values are whole numbers written plainly, every Ensembl release file among them, would not change.

**The class it produces.** `gff_version_not_a_number`: the version part of the `Gene`, `Feature`, `ENSP`, `HGVSc` and `HGVSp` values, and of the JSON `gene_id`, `transcript_id`, `protein_id`, `hgvsc` and `hgvsp`, of a transcript read from a GTF or GFF3. The transcript's gene line, transcript line or first CDS line in transcript order carries a `version` that is not a positive whole number written plainly; on a transcript line the attribute is `version`, else `transcript_version`, the GTF spelling. Ensembl VEP's `--gff` and `--gtf` sources print the text as written. vep-rs refuses to build a cache from the file, naming the line and the text, and reads text that is a whole number, leading zeros included, as the number.

The cache holds each version as a 32-bit number, so digits above 4294967295 on any of the three lines end the build with an error naming the line and the limit, where Ensembl VEP prints the digits ([section 4](#4-the-same-transcript-with-version4294967296-on-its-three-lines)).

The class is documented only: every `version` on chromosomes 21 and MT of the release-116 GRCh38 GFF3 is a plain whole number. The committed GTF's gene, transcript and CDS lines carry none, and the `transcript_version` of each of its 59 transcript lines is a plain whole number.

**How often the shape occurs.** Tier 3: the shape is a `version` whose text is not a plain positive whole number, which no release file carries. The 964 feature lines of the committed GTF carry none, 14 committed flag sets run `--protein_version`, `--transcript_version`, `--gene_version` or `--hgvs` on a cache built from a GTF or GFF3, and no committed flag set runs Ensembl VEP over a file with the shape.

## 1. A TMPRSS3 transcript whose first CDS line carries `version=v2`

**In short.** The GRCh38 release-116 GFF3 lines of the gene TMPRSS3 and its transcript `ENST00001113783`, with `version=1` rewritten `version=v2` on the transcript's first CDS line in transcript order, the line at `42395324-42395417`, beside the lines of PWP2 and its transcript `ENST00000456705` as released. On the record `tmprss3_codon6`, `21 42395401 G>A`, Ensembl VEP 116.2 reading the file through `--gff` writes `ENSP=ENSP00000783588.v2` and `HGVSp=ENSP00000783588.v2:p.Pro6Leu`; vep-rs refuses to build a cache from the file at that CDS line. On the PWP2 record `pwp2_codon3` Ensembl VEP writes `ENSP00000411749.1:p.Phe3Ser`; vep-rs, with no cache built, writes no row.

### 1.1 What a reader needs to know first

A GFF3 CDS line carries the protein's identifier as `protein_id` and its version as `version`, and Ensembl VEP's `--gff` source reads both from the first CDS line in transcript order. Under `--protein_version` the version follows the identifier in `ENSP`, and under `--hgvs` it follows the identifier in the HGVSp reference name, `ENSP00000783588.1:p.Pro6Leu`.

### 1.2 The record

The transcript's two CDS lines, the first in transcript order (the higher coordinates, on the reverse strand) with `version=v2`; the gene, transcript and exon lines, and the PWP2 lines, are as released and omitted here:

```
21	havana_tagene	CDS	42394915	42395084	.	-	2	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1
21	havana_tagene	CDS	42395324	42395417	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=v2
```

The input:

```
21	42395401	tmprss3_codon6	G	A	.	.	.
21	44107424	pwp2_codon3	T	C	.	.	.
```

Both engines run `--hgvs --protein --protein_version --symbol`. Ensembl VEP reads the file, bgzipped and indexed as `tmprss3_cds_version_v2.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.2:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;ENSP=ENSP00000783588.v2;SOURCE=tmprss3_cds_version_v2.gff3.gz;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.v2:p.Pro6Leu
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411749.1;SOURCE=tmprss3_cds_version_v2.gff3.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

vep-rs, building the cache:

```
ERROR: [cache-builder] Failed to parse GFF3: tmprss3_cds_version_v2.gff3 line 6: version "v2" is not a whole number
```

The build exits 1 and writes no cache, so vep-rs annotates nothing from the file. The same values stand in the `ENSP` and `HGVSp` fields of Ensembl VEP's VCF `CSQ` entry and in its JSON `protein_id` and `hgvsp`. With the line rewritten `version=01` instead, Ensembl VEP writes `ENSP00000783588.01` in both places and vep-rs builds the cache and writes `ENSP00000783588.1`. The same lines read a GTF: with `version "v2"` added to the committed GTF's copy of the line (a GTF CDS line carries `protein_version`, which neither engine reads as the version), Ensembl VEP's `--gtf` source writes `ENSP00000783588.v2` and vep-rs refuses the build at that line.

### 1.4 Why Ensembl VEP writes that

`_add_translation` takes the first ordered CDS record's `version`, the text `v2`, and builds the translation with it; `v2` is true to Perl, so the `|| 1` default does not apply. The `ENSP` writer and `hgvs_protein` each append `.` and that text, since `ENSP00000783588` does not already end in a dot and digits. No step between the attribute and the output reads the text as a number.

### 1.5 Why vep-rs's refusal is the correct outcome

The file's `version` is a version number by the convention of the files that carry it, and `v2` is not one: the transcript line of the same file carries `version=1`, a whole number, and so does the transcript's other CDS line. An identifier with the text appended, `ENSP00000783588.v2`, names no release of the protein. Refusing the file names the line to correct, and a corrected file prints a version the identifier can carry; printing the text prints a version that does not exist.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 tmprss3_cds_version_v2.gff3 --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
vep -i records.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --hgvs --protein --protein_version --symbol
```

vep-rs's cache builder refuses the file at its sixth line with the error above and writes no cache. Ensembl VEP 116.2 with `--gff tmprss3_cds_version_v2.gff3.gz`, the file bgzipped and indexed, writes `ENSP00000783588.v2`. Restoring `version=1` on the line makes vep-rs build the cache and both engines write `ENSP00000783588.1`.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/gff_version_not_a_number.toml` holds a file carrying the three shapes of this page at once, the records, Ensembl VEP's rows and vep-rs's refusal. The class covers the version part of the `Gene`, `Feature`, `ENSP`, `HGVSc` and `HGVSp` values, and of the JSON `gene_id`, `transcript_id`, `protein_id`, `hgvsc` and `hgvsp`. The transcript's gene line, transcript line or first CDS line in transcript order carries a `version` that is not a plain positive whole number. No committed dataset reaches the shape, so the class is documented only.

## 2. The same transcript with `version=v2` on its transcript line

**In short.** The same lines, with `version=1` rewritten `version=v2` on the transcript line of `ENST00001113783` instead, and every CDS line as released. On `tmprss3_codon6`, Ensembl VEP 116.2 writes `ENST00001113783.v2` in `Feature` and `HGVSc=ENST00001113783.v2:c.17C>T`; vep-rs refuses to build a cache from the file at the transcript line. On `pwp2_codon3` Ensembl VEP writes `ENST00000456705.1` in both places; vep-rs, with no cache built, writes no row.

### 2.1 What a reader needs to know first

A GFF3 transcript line carries the transcript's version as `version`; a GTF transcript line carries it as `transcript_version`, and Ensembl VEP reads `version` first, then `transcript_version`, from either. Under `--transcript_version` the version follows the identifier in `Feature`, and under `--hgvs` it follows the identifier in the HGVSc reference name whether or not `--transcript_version` is given. A transcript line without either attribute gives the transcript no version, and both engines print the identifier bare.

### 2.2 The record

The transcript line with `version=v2`; the gene, exon and CDS lines, and the PWP2 lines, are as released and omitted here:

```
21	havana_tagene	mRNA	42393803	42396122	.	-	.	ID=transcript:ENST00001113783;Parent=gene:ENSG00000160183;Name=TMPRSS3-224;biotype=protein_coding;tag=gencode_basic;transcript_id=ENST00001113783;version=v2
```

The input is the two records of 1.2. Both engines run `--hgvs --protein --protein_version --symbol --transcript_version --gene_version`. Ensembl VEP reads the file, bgzipped and indexed as `tmprss3_transcript_version_v2.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 2.3 What each engine writes

Ensembl VEP 116.2:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183.19	ENST00001113783.v2	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;ENSP=ENSP00000783588.1;SOURCE=tmprss3_transcript_version_v2.gff3.gz;HGVSc=ENST00001113783.v2:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu
pwp2_codon3	21:44107424	C	ENSG00000241945.10	ENST00000456705.1	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411749.1;SOURCE=tmprss3_transcript_version_v2.gff3.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

vep-rs, building the cache:

```
ERROR: [cache-builder] Failed to parse GFF3: tmprss3_transcript_version_v2.gff3 line 3: version "v2" is not a whole number
```

The build exits 1 and writes no cache. The same values stand in the `Feature` and `HGVSc` fields of Ensembl VEP's VCF `CSQ` entry and in its JSON `transcript_id` and `hgvsc`. Without `--transcript_version` Ensembl VEP writes `ENST00001113783` in `Feature` and `.v2` in the HGVSc reference name still.

With the line rewritten `version=01` instead, Ensembl VEP writes `ENST00001113783.01` in both places and vep-rs builds the cache and writes `ENST00001113783.1`. With `transcript_version=3` added beside `version=v2`, neither engine reads it: Ensembl VEP writes `.v2` and vep-rs refuses the `version`. The same transcript read from a GTF, its line's `transcript_version "1"` rewritten `transcript_version "v2"`: Ensembl VEP's `--gtf` source writes `ENST00001113783.v2` in both places and vep-rs refuses the build at that line.

A `transcript_version` of `0`, or an empty one, with no `version` beside it reaches one field alone. Ensembl VEP writes `HGVSc=ENST00001113783.0:c.17C>T`, or `ENST00001113783.:c.17C>T` for the empty value, and `ENST00001113783` in `Feature` even under `--transcript_version`; vep-rs writes the bare identifier in both. A GTF line's `transcript_version "0"` or `transcript_version ""` does the same through `--gtf`. A `version=0` with no `transcript_version` is no version to either engine, and the rows agree.

### 2.4 Why Ensembl VEP writes that

`_create_transcript` takes the transcript line's `version`, the text `v2`, and builds the transcript with it. The `Feature` writer appends `.` and that text under `--transcript_version`, and `hgvs_transcript` appends it to the HGVSc reference name, since `ENST00001113783` does not already end in a dot and digits and the version is defined. No step between the attribute and the output reads the text as a number.

### 2.5 Why vep-rs's refusal is the correct outcome

`v2` is not a version number: every other transcript line of the file carries a whole number, and `ENST00001113783.v2` names no release of the transcript; `ENST00001113783.v2:c.17C>T` is an HGVS description whose reference name resolves to nothing. Refusing the file names the line to correct, and a corrected file prints a name that exists; printing the text prints one that does not.

### 2.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 tmprss3_transcript_version_v2.gff3 \
  --genome-fasta <corpus>/reference.fa.gz --output-dir cache
vep -i records.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --hgvs --protein --protein_version --symbol --transcript_version --gene_version
```

vep-rs's cache builder refuses the file at its third line with the error above and writes no cache. Ensembl VEP 116.2 with `--gff tmprss3_transcript_version_v2.gff3.gz`, the file bgzipped and indexed, writes `ENST00001113783.v2` in both. Restoring `version=1` on the line makes vep-rs build the cache and both engines write `ENST00001113783.1`.

## 3. The same transcript with `version=v2` on its gene line

**In short.** The same lines, with `version=19` rewritten `version=v2` on the gene line of `ENSG00000160183` instead, and the transcript and CDS lines as released. On `tmprss3_codon6`, Ensembl VEP 116.2 writes `ENSG00000160183.v2` in `Gene`; vep-rs refuses to build a cache from the file at the gene line. On `pwp2_codon3` Ensembl VEP writes `ENSG00000241945.10`; vep-rs, with no cache built, writes no row.

### 3.1 What a reader needs to know first

A GFF3 gene line carries the gene's version as `version`, and Ensembl VEP's `--gff` and `--gtf` sources read it from that attribute alone: a GTF gene line's `gene_version` is not read, so a gene read from an Ensembl GTF has no version in either engine. Under `--gene_version` the version follows the identifier in `Gene`; no HGVS name carries it. A gene line without `version` gives the gene no version, and both engines print the identifier bare.

### 3.2 The record

The gene line with `version=v2`; the transcript, exon and CDS lines, and the PWP2 lines, are as released and omitted here:

```
21	ensembl_havana	gene	42371840	42396122	.	-	.	ID=gene:ENSG00000160183;Name=TMPRSS3;biotype=protein_coding;description=transmembrane serine protease 3 [Source:HGNC Symbol%3BAcc:HGNC:11877];gene_id=ENSG00000160183;logic_name=ensembl_havana_gene_homo_sapiens;version=v2
```

The input and the flags are those of 2.2. Ensembl VEP reads the file, bgzipped and indexed as `tmprss3_gene_version_v2.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 3.3 What each engine writes

Ensembl VEP 116.2:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183.v2	ENST00001113783.1	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;ENSP=ENSP00000783588.1;SOURCE=tmprss3_gene_version_v2.gff3.gz;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu
pwp2_codon3	21:44107424	C	ENSG00000241945.10	ENST00000456705.1	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411749.1;SOURCE=tmprss3_gene_version_v2.gff3.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

vep-rs, building the cache:

```
ERROR: [cache-builder] Failed to parse GFF3: tmprss3_gene_version_v2.gff3 line 2: version "v2" is not a whole number
```

The build exits 1 and writes no cache. The same value stands in the `Gene` field of Ensembl VEP's VCF `CSQ` entry and in its JSON `gene_id`. Without `--gene_version` Ensembl VEP writes `ENSG00000160183`. With the line rewritten `version=01` instead, Ensembl VEP writes `ENSG00000160183.01` and vep-rs builds the cache and writes `ENSG00000160183.1`. The same gene read from a GTF, with `version "v2"` added to the committed GTF's copy of its line: Ensembl VEP's `--gtf` source writes `ENSG00000160183.v2` and vep-rs refuses the build at that line.

### 3.4 Why Ensembl VEP writes that

`_add_identifiers` takes the gene line's `version`, the text `v2`, as the gene version of every transcript of the gene. The `Gene` writer appends `.` and that text under `--gene_version`, since `ENSG00000160183` does not already end in a dot and digits and the text is true to Perl. No step between the attribute and the output reads the text as a number.

### 3.5 Why vep-rs's refusal is the correct outcome

`v2` is not a version number: the file's other gene line carries `version=10`, a whole number, and `ENSG00000160183.v2` names no release of the gene. Refusing the file names the line to correct, and a corrected file prints a name that exists; printing the text prints one that does not.

### 3.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 tmprss3_gene_version_v2.gff3 \
  --genome-fasta <corpus>/reference.fa.gz --output-dir cache
vep -i records.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --hgvs --protein --protein_version --symbol --transcript_version --gene_version
```

vep-rs's cache builder refuses the file at its second line with the error above and writes no cache. Ensembl VEP 116.2 with `--gff tmprss3_gene_version_v2.gff3.gz`, the file bgzipped and indexed, writes `ENSG00000160183.v2`. Restoring `version=19` on the line makes vep-rs build the cache and both engines write `ENSG00000160183.19`.

## 4. The same transcript with `version=4294967296` on its three lines

**In short.** The same lines, with `version=4294967296`, one above the largest 32-bit number, on the gene line, the transcript line and the transcript's first CDS line in transcript order, and the PWP2 lines as released. Ensembl VEP 116.2 writes the digits after every identifier, `ENSG00000160183.4294967296`, `ENST00001113783.4294967296` and `ENSP00000783588.4294967296`, in `Gene`, `Feature`, `ENSP` and both HGVS reference names. vep-rs refuses to build the cache: the build ends at the first such line with an error naming the file, the line and the limit, and writes no row. With `version=4294967295` on the same lines both engines write `.4294967295` everywhere.

### 4.1 What a reader needs to know first

The cache holds a gene's, a transcript's and a translation's version as a 32-bit unsigned number, so 4294967295 is the largest version a built cache carries. Ensembl VEP holds the attribute's text and prints it at any width. No release file carries a version near the limit: the largest `version` on chromosomes 21 and MT of the release-116 GRCh38 GFF3 is 24.

### 4.2 The record

The three lines of 1.2, 2.2 and 3.2 with `version=4294967296` in place of `v2`; the exon lines and the PWP2 lines are as released and omitted here:

```
21	ensembl_havana	gene	42371840	42396122	.	-	.	ID=gene:ENSG00000160183;Name=TMPRSS3;biotype=protein_coding;description=transmembrane serine protease 3 [Source:HGNC Symbol%3BAcc:HGNC:11877];gene_id=ENSG00000160183;logic_name=ensembl_havana_gene_homo_sapiens;version=4294967296
21	havana_tagene	mRNA	42393803	42396122	.	-	.	ID=transcript:ENST00001113783;Parent=gene:ENSG00000160183;Name=TMPRSS3-224;biotype=protein_coding;tag=gencode_basic;transcript_id=ENST00001113783;version=4294967296
21	havana_tagene	CDS	42395324	42395417	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=4294967296
```

The input and the flags are those of 2.2. Ensembl VEP reads the file, bgzipped and indexed as `tmprss3_version_above_limit.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs is asked to build a cache from the file and that FASTA.

### 4.3 What each engine writes

Ensembl VEP 116.2:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183.4294967296	ENST00001113783.4294967296	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;ENSP=ENSP00000783588.4294967296;SOURCE=tmprss3_version_above_limit.gff3.gz;HGVSc=ENST00001113783.4294967296:c.17C>T;HGVSp=ENSP00000783588.4294967296:p.Pro6Leu
pwp2_codon3	21:44107424	C	ENSG00000241945.10	ENST00000456705.1	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411749.1;SOURCE=tmprss3_version_above_limit.gff3.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

vep-rs, from `vep-cache-builder`, with no cache written and no row:

```
ERROR: [cache-builder] Failed to parse GFF3: tmprss3_version_above_limit.gff3 line 2: version 4294967296 is above 4294967295, the largest version the cache holds
```

The same values stand in the JSON `gene_id`, `transcript_id`, `protein_id`, `hgvsc` and `hgvsp` and in the VCF `CSQ` entry of Ensembl VEP's output. With the digits on one line alone, Ensembl VEP writes them in that line's fields only and the build names that line. A GTF transcript line's `transcript_version "4294967296"` is read the same way: Ensembl VEP's `--gtf` source writes `ENST00001113783.4294967296` in `Feature` and HGVSc, and the build names the line in the GTF.

### 4.4 Why Ensembl VEP writes that

The three reads of section 1.4, 2.4 and 3.4 keep the attribute's text, and the writers append it; a string of digits passes the only test on the way to the output, that the identifier does not already end in a dot and digits. Nothing in Ensembl VEP bounds a version.

### 4.5 Why vep-rs's refusal is the correct outcome

The cache format stores each version in a 32-bit field. Reading a larger number as no version, or as the number's low bits, would print a version the file does not state. Printing the digits as Ensembl VEP does would need a wider version field in the cache, a format change for an input no release file carries. The error names the line, so the file can be corrected and the build repeated.

### 4.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 tmprss3_version_above_limit.gff3 \
  --genome-fasta <corpus>/reference.fa.gz --output-dir cache
```

The build ends with the error above and exit status 1. Ensembl VEP 116.2 with `--gff tmprss3_version_above_limit.gff3.gz`, the file bgzipped and indexed, writes `ENST00001113783.4294967296` for `tmprss3_codon6`. Rewriting the three values `4294967295` builds the cache and makes both engines write `.4294967295` in every version-bearing field.
