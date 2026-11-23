# Two CDS lines of one GFF3 transcript sharing a start

**The defect in one sentence.** Ensembl VEP 116.3's `--gff` source stores a transcript's CDS lines in a hash keyed by their start, so when two CDS lines of one transcript share a start one replaces the other. The later line of the file is kept when the two follow their transcript line, as a sorted file reads them. Everything the dropped line alone stated is lost without a warning: its protein identifier and version when the pair opens the coding sequence, and its end when the pair closes it.

**Why this is a defect.** The file states two coding segments at that start, and a reader that keeps one of them reads a coding sequence the file does not describe. Which segment survives depends on the order the two lines stand in the file, and on whether their transcript line precedes them, not on anything either line says. When the pair closes the coding sequence and the line it keeps ends first, the bases between the two ends are coding by the file and untranslated to Ensembl VEP. A variant there is `3_prime_UTR_variant` with no protein change, where the file's coding sequence gives a `missense_variant`. Nothing in the output or on standard error says a line was dropped.

**What vep-rs does instead.** vep-rs reads every CDS line of the transcript: the coding span is the union of the lines, so a base any line calls coding is coding. The translation's identifier and version come from the line that opens the coding sequence, the lowest start on the forward strand and the highest on the reverse. Where two lines share that start, the forward strand takes the earlier line in the file and the reverse strand the later, the same line for a file however many times it is read. The lines sharing a start are never both dropped and never read as two translations.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, lines 385-401, `_get_parent_child_structure`: a record's children are pushed in the order the records were read, the file's order within the region.
- `BaseGXF.pm`, lines 504-519, `_create_transcript`: line 511 stores each `CDS` or `stop_codon` child under its start, `$cdss{$child->{start}} = $child`, so a later child at the same start replaces the earlier one. Lines 544-556 collect the stored records exon by exon into the ordered list; a start is one key, so one record.
- `BaseGXF.pm`, lines 406-422, the orphan path for a record read before its parent: line 412 skips a record whose parent, identifier and start were already filed. Two CDS lines that share their `ID`, as Ensembl's files write a transcript's CDS lines, and are read before their transcript line therefore keep the earlier. A sorted file reads the transcript line first, since it starts no later than its CDS lines.
- `BaseGXF.pm`, lines 719-756, `_add_translation`: lines 723-734 take the protein identifier and version from the first ordered record; lines 749-756 set the translation's end from the last ordered record's end on the forward strand, its start on the reverse.

**What correcting it would change.** With the CDS children kept as a list, every line's span would reach the translation, as it does in vep-rs, and the pair's identifiers would come from the file rather than from the order of two lines. The rows on a file with one CDS line per start would not change.

**The class it produces.** `gff_cds_lines_sharing_a_start`: the `ENSP` and `HGVSp` values of a transcript whose GTF or GFF3 CDS lines include two at one start. Where the pair closes the coding sequence and the kept line ends first, it also covers the terms, positions, codons, amino acids and HGVSc of a variant between the two ends. Ensembl VEP writes the identifier, version and end of the one line it keeps, the later as a sorted file reads them. vep-rs writes the identifiers of the line that opens the coding sequence, the earlier of the pair on the forward strand and the later on the reverse, and reads the longer end. Documented only: no committed file carries such a pair.

**How often the shape occurs.** Tier 3: the shape is two CDS or stop-codon lines of one transcript at one start, which no release file carries. The 964 feature lines of the committed GTF carry none, 14 committed flag sets run a cache built from a GTF or GFF3, and no committed flag set runs Ensembl VEP over a file with the shape.

## 1. A PWP2 transcript whose first and last CDS lines are each written twice

**In short.** The GRCh38 release-116 GFF3 lines of the gene PWP2 and its transcript `ENST00000456705`, on the forward strand, beside the lines of TMPRSS3 and its transcript `ENST00001113783` as released. The transcript's first CDS line, `44107417-44107434`, is followed by a copy naming the protein `ENSP00000411750` version 2. Its last CDS line, `44115691-44115744`, is followed by a copy naming the same protein and ending at `44115738`, six bases earlier. On `pwp2_codon3`, `21 44107424 T>C`, Ensembl VEP 116.3 reading the file through `--gff` writes `ENSP00000411750.2` and `ENSP00000411750.2:p.Phe3Ser`, vep-rs `ENSP00000411749.1` and `ENSP00000411749.1:p.Phe3Ser`. On `pwp2_codon186`, `21 44115741 G>A`, between the two ends, Ensembl VEP writes `3_prime_UTR_variant` and `c.*3G>A` with no HGVSp, vep-rs `missense_variant`, `c.557G>A` and `p.Arg186His`.

### 1.1 What a reader needs to know first

A GFF3 transcript's coding sequence is the set of its CDS lines, each a span with a phase, a `protein_id` and a `version`. Ensembl VEP's `--gff` source builds one translation per transcript from them: its identifiers from the first line in transcript order, its start from that line and its end from the last. A file with two CDS lines at one start says two things about that start; the format gives it no meaning, and a reader has to decide what to keep.

### 1.2 The record

The transcript's eight CDS lines, the released six and the two copies; the gene, transcript and exon lines, and the TMPRSS3 lines, are as released and omitted here:

```
21	havana	CDS	44107417	44107434	.	+	0	ID=CDS:ENSP00000411749;Parent=transcript:ENST00000456705;protein_id=ENSP00000411749;version=1
21	havana	CDS	44107417	44107434	.	+	0	ID=CDS:ENSP00000411750;Parent=transcript:ENST00000456705;protein_id=ENSP00000411750;version=2
21	havana	CDS	44108984	44109096	.	+	0	ID=CDS:ENSP00000411749;Parent=transcript:ENST00000456705;protein_id=ENSP00000411749;version=1
21	havana	CDS	44113753	44113847	.	+	1	ID=CDS:ENSP00000411749;Parent=transcript:ENST00000456705;protein_id=ENSP00000411749;version=1
21	havana	CDS	44114625	44114768	.	+	2	ID=CDS:ENSP00000411749;Parent=transcript:ENST00000456705;protein_id=ENSP00000411749;version=1
21	havana	CDS	44115264	44115399	.	+	2	ID=CDS:ENSP00000411749;Parent=transcript:ENST00000456705;protein_id=ENSP00000411749;version=1
21	havana	CDS	44115691	44115744	.	+	1	ID=CDS:ENSP00000411749;Parent=transcript:ENST00000456705;protein_id=ENSP00000411749;version=1
21	havana	CDS	44115691	44115738	.	+	1	ID=CDS:ENSP00000411750;Parent=transcript:ENST00000456705;protein_id=ENSP00000411750;version=2
```

The input:

```
21	44107424	pwp2_codon3	T	C	.	.	.
21	44115741	pwp2_codon186	G	A	.	.	.
21	42395401	tmprss3_codon6	G	A	.	.	.
```

Both engines run `--hgvs --protein --protein_version --symbol`. Ensembl VEP reads the file, bgzipped and indexed as `pwp2_cds_shared_start.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.3:

```
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411750.2;SOURCE=pwp2_cds_shared_start.gff3.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411750.2:p.Phe3Ser
pwp2_codon186	21:44115741	A	ENSG00000241945	ENST00000456705	Transcript	3_prime_UTR_variant	601	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411750.2;SOURCE=pwp2_cds_shared_start.gff3.gz;HGVSc=ENST00000456705.1:c.*3G>A
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;ENSP=ENSP00000783588.1;SOURCE=pwp2_cds_shared_start.gff3.gz;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu
```

vep-rs:

```
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411749.1;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
pwp2_codon186	21:44115741	A	ENSG00000241945	ENST00000456705	Transcript	missense_variant	601	557	186	R/H	cGc/cAc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411749.1;HGVSc=ENST00000456705.1:c.557G>A;HGVSp=ENSP00000411749.1:p.Arg186His
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;ENSP=ENSP00000783588.1;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu
```

The same values stand in the `CSQ` entries of the VCF output and in the JSON `protein_id`, `hgvsp`, `consequence_terms`, `cds_start`, `protein_start`, `codons` and `amino_acids`. Ensembl VEP's JSON `most_severe_consequence` for `pwp2_codon186` is `3_prime_UTR_variant` and vep-rs's `missense_variant`. The TMPRSS3 row differs only in what the class [`file_cache_writes_no_source_field`](file-cache-writes-no-source-field.md) covers, its `SOURCE` item above, in every format. With the last line's copy ending where the released line ends, so that each pair differs in `protein_id` and `version` alone, the two engines differ in `ENSP` and `HGVSp` alone and both write `missense_variant` for `pwp2_codon186`. With the first CDS line of the reverse-strand TMPRSS3 transcript written twice the same way, both engines take the later line and agree.

The orphan path shows with the copy of the first PWP2 line given the first line's `ID`, as Ensembl's files write a transcript's CDS lines. With the transcript and first exon lines moved to the CDS start and written after that pair, Ensembl VEP reads the pair before its transcript line and keeps the earlier, so it writes `ENSP00000411749.1` as vep-rs does.

The same lines read a GTF: with the committed GTF's lines of the two transcripts and the two copies written in its form, Ensembl VEP's `--gtf` source writes `ENSP00000411750.1` on the PWP2 rows and `3_prime_UTR_variant` for `pwp2_codon186`, vep-rs `ENSP00000411749.1` and `missense_variant`.

### 1.4 Why Ensembl VEP writes that

`_create_transcript` stores the eight CDS children under their starts, and the two copies, each later in the file than the line it copies, replace them: six records remain, the first and the last the copies. `_add_translation` reads `ENSP00000411750` and `2` from the first record, so every row's `ENSP` and HGVSp name that protein, and sets the translation's end from the last record's end, `44115738`. The variant at `44115741` lies after that end, inside the exon, so it is untranslated: `3_prime_UTR_variant`, `c.*3G>A`, no HGVSp.

### 1.5 Why vep-rs's row is the correct one

Every base a CDS line of the file covers is coding by the file, and `44115741` lies inside the released last line, so the file places codon 186 there. A reading that ends the coding sequence six bases earlier because a second line at that start came later is not a reading of the file. The identifiers are a tie the file does not break; vep-rs breaks it the same way every time, and the protein it names is one of the two the file names at that start. Ensembl VEP's choice is one of the two as well, but it changes with the order of two lines the format gives no meaning to.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 pwp2_cds_shared_start.gff3 --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
vep -i records.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --hgvs --protein --protein_version --symbol
```

vep-rs writes `ENSP00000411749.1` on both PWP2 rows and `missense_variant` for `pwp2_codon186`. Ensembl VEP 116.3 with `--gff pwp2_cds_shared_start.gff3.gz`, the file bgzipped and indexed, writes `ENSP00000411750.2` and `3_prime_UTR_variant`. Removing the two copies makes the engines agree on every row.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/gff_cds_lines_sharing_a_start.toml` holds the file, the records and both engines' rows. The class covers the `ENSP` and `HGVSp` values, and the JSON `protein_id` and `hgvsp`, of a transcript two of whose CDS lines share a start. It also covers the consequence terms, positions, codons, amino acids and HGVSc of a variant between the two ends of a pair that closes the coding sequence. No committed dataset reaches the shape, so the class is documented only.
