# The HGVSp of a GFF3 translation whose first CDS line names no protein

**The defect in one sentence.** Ensembl VEP 116.3's `--gff` source takes a translation's protein identifier from the transcript's first CDS line in transcript order, and when that line carries no `protein_id` it builds the translation with none. Its HGVSp writer then names the reference sequence as the empty identifier followed by the version, so the protein notation is written `.1:p.Pro6Leu`, a change on a sequence the notation does not name.

**Why this is a defect.** An HGVS protein description is a reference sequence name, a colon, `p.` and the change, and the name is what makes it resolvable. `.1:p.Pro6Leu` has no name: nothing can be looked up, and the `.1` is the version of nothing. The `ENSP` field of the same row is empty, so one row says both that the transcript's protein has no identifier and that version 1 of it carries a change at residue 6. The writer appends the version to whatever the translation's display identifier is, and for a translation with neither a stable nor a database identifier that is the empty string.

**What vep-rs does instead.** A cache vep-rs builds from the file and the genome writes no `HGVSp` and no `ENSP` for a translation whose first CDS line names no protein, as both engines write none for a transcript without a translation. A build of the cache's content that reads the core table dumps takes their identifier and writes both. The rest of the row, the consequence term, the positions, the codons and the amino acids, is written as Ensembl VEP writes it. Later CDS lines of the transcript that do name the protein change nothing in either engine: both read the identifier from the first line in transcript order alone.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` and `ensembl` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, lines 719-734, `_add_translation`: line 723 reads `protein_id` from the first ordered CDS record, undefined when the attribute is absent, and line 726 copies it to the transcript's `_protein` only when it is set. Lines 730-734 build the translation with that identifier and the version of line 727, the line's `version` or 1.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1747-1760, `hgvs_protein`: line 1748 takes the reference name from the translation's `display_id`, and lines 1759-1760 append `.` and the translation's version unless the version is undefined or the name already ends in one. Line 1905, in `_get_hgvs_protein_format`, writes the name, `:p.` and the change.
- `Bio/EnsEMBL/Translation.pm`, lines 907-910, `display_id`: the stable identifier, else the database identifier, else the empty string.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1565-1574: `ENSP` is written from the transcript's `_protein` when that is set, with the version appended under `--protein_version`; for this translation nothing is written.

**What correcting it would change.** With `hgvs_protein` returning nothing for a translation without a name, as it returns nothing for a transcript without a translation, the row would carry its HGVSc and no HGVSp, as vep-rs's row does; the `ENSP` field is already empty in both engines.

**The class it produces.** `gff_cds_without_protein_id_hgvsp`: the `HGVSp` value, and the JSON `hgvsp`, of a row on a transcript whose first CDS line in transcript order, read from a GTF or GFF3, names no `protein_id`. Ensembl VEP's `--gff` and `--gtf` sources write the version-only name, `.1:p.Pro6Leu`; a cache vep-rs built from the file and the genome writes none, and one built with the core table dumps writes the dumps' identifier. The class is documented only: the release-116 GRCh38 GFF3 names a protein on every CDS line of chromosomes 21 and MT, and the committed GTF's 313 CDS lines all carry `protein_id`.

**How often the shape occurs.** Tier 3: the shape is a CDS line without `protein_id` at the 5' end of a transcript's coding sequence, which no release file carries. The 964 feature lines of the committed GTF carry none, no committed flag set runs `--hgvs` on a cache built from a GTF or GFF3, and none runs Ensembl VEP over a file with the shape.

## 1. A TMPRSS3 transcript whose first CDS line names no protein

**In short.** The GRCh38 release-116 GFF3 lines of the gene TMPRSS3 and its transcript `ENST00001113783`, with `protein_id=ENSP00000783588` removed from the transcript's first CDS line in transcript order, beside the lines of PWP2 and its transcript `ENST00000456705` as released. That line, at `42395324-42395417`, holds the 5' end of the coding sequence on the reverse strand. On the record `tmprss3_codon6`, `21 42395401 G>A`, a `missense_variant` at codon 6, Ensembl VEP 116.3 reading the file through `--gff` writes `HGVSp=.1:p.Pro6Leu` and no `ENSP`; vep-rs writes neither. On the PWP2 record `pwp2_codon3`, `21 44107424 T>C`, both engines write `ENSP00000411749.1:p.Phe3Ser`.

### 1.1 What a reader needs to know first

A GFF3 CDS line carries the protein's identifier as `protein_id` and its version as `version`; Ensembl's files write both on every CDS line of a transcript, and Ensembl VEP's `--gff` source reads them from the first CDS line in transcript order alone. An HGVSp notation names the protein sequence the change is described against, `ENSP00000783588.1:p.Pro6Leu`: the identifier, its version, then the change.

### 1.2 The record

The transcript's two CDS lines, the first in transcript order (the higher coordinates, on the reverse strand) without its `protein_id`; the gene, transcript and exon lines, and the PWP2 lines, are as released and omitted here:

```
21	havana_tagene	CDS	42394915	42395084	.	-	2	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;protein_id=ENSP00000783588;version=1
21	havana_tagene	CDS	42395324	42395417	.	-	0	ID=CDS:ENSP00000783588;Parent=transcript:ENST00001113783;version=1
```

The input:

```
21	42395401	tmprss3_codon6	G	A	.	.	.
21	44107424	pwp2_codon3	T	C	.	.	.
```

Both engines run `--hgvs --protein --protein_version --symbol`. Ensembl VEP reads the file, bgzipped and indexed as `tmprss3_cds_no_protein_id.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.3:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;SOURCE=tmprss3_cds_no_protein_id.gff3.gz;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=.1:p.Pro6Leu
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411749.1;SOURCE=tmprss3_cds_no_protein_id.gff3.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

vep-rs:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;HGVSc=ENST00001113783.1:c.17C>T
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;ENSP=ENSP00000411749.1;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

In the VCF output the `HGVSp` field of the TMPRSS3 `CSQ` entry is `.1:p.Pro6Leu` in Ensembl VEP's and empty in vep-rs's; in the JSON output Ensembl VEP's `hgvsp` is `".1:p.Pro6Leu"` and vep-rs's transcript object has no `hgvsp`; neither engine writes `ENSP` or `protein_id` for it. Beyond `HGVSp`, the TMPRSS3 rows differ only in what the class [`file_cache_writes_no_source_field`](file-cache-writes-no-source-field.md) covers, the `SOURCE` item above: the consequence, the positions, the codons and the amino acids agree. The PWP2 rows differ only in what that class covers, in every format. The same lines read a GTF: with `protein_id "ENSP00000783588"` removed from the committed GTF's copy of the line, Ensembl VEP's `--gtf` source writes `HGVSp=.1:p.Pro6Leu` and vep-rs none.

### 1.4 Why Ensembl VEP writes that

`_add_translation` reads `protein_id` from the first ordered CDS record and finds none, so the translation is built with no stable identifier and with the version 1 from the line's `version`, and the transcript's `_protein` stays unset, which leaves `ENSP` empty. `hgvs_protein` takes the reference name from the translation's `display_id`, the empty string for a translation with neither a stable nor a database identifier, appends `.1` because the name does not already end in a version, and writes `.1:p.Pro6Leu`.

### 1.5 Why vep-rs's row is the correct one

The file names no protein for this translation, and the row says so in its empty `ENSP`. An HGVS protein description without a reference name describes a change on no sequence: `.1` is the version of an identifier the file never gave, and a reader has nothing to resolve it against. Writing no HGVSp, as both engines write none for a transcript without a translation, states what the file states; `.1:p.Pro6Leu` states a version for a protein that was never named.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 tmprss3_cds_no_protein_id.gff3 --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
vep -i records.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --hgvs --protein --protein_version --symbol
```

vep-rs writes no `HGVSp` for `tmprss3_codon6`. Ensembl VEP 116.3 with `--gff tmprss3_cds_no_protein_id.gff3.gz`, the file bgzipped and indexed, writes `HGVSp=.1:p.Pro6Leu`. Restoring `protein_id=ENSP00000783588` on the line makes both engines write `ENSP=ENSP00000783588.1` and `HGVSp=ENSP00000783588.1:p.Pro6Leu`.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/gff_cds_without_protein_id_hgvsp.toml` holds the file, the records and both engines' rows. The class covers the `HGVSp` value, and the JSON `hgvsp`, of a row on a transcript whose first CDS line in transcript order names no protein. No committed dataset reaches the shape, so the class is documented only.
