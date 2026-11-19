# The `ncRNA_gene` lines of an Ensembl GFF3

**The defect in one sentence.** Ensembl VEP 116.2's `--gff` source reads a line only when its type is on a fixed list, and `ncRNA_gene`, the type an Ensembl GFF3 gives every non-coding gene, is not on it. The gene line is skipped, every transcript whose `Parent` names it is pruned, and a record inside the gene is written as `intergenic_variant`. vep-rs reads the gene line, keeps its transcripts in the cache it builds and annotates the record against them.

**Why this is a defect.** An Ensembl GFF3 types each gene line by the Sequence Ontology term of its biotype: `gene` for a protein-coding gene, `pseudogene` for a pseudogene and `ncRNA_gene` for every non-coding gene, lncRNA, miRNA, misc_RNA, snRNA, snoRNA and the mitochondrial tRNA and rRNA genes alike. The list the file source reads admits the transcript types those genes carry (`lnc_RNA`, `miRNA`, `snRNA`, `snoRNA`, `rRNA`, `tRNA`, `transcript`) and six older gene types, but not the one type the file writes for all of them. Of the 646,577 transcript lines of the GRCh38 release-116 GFF3, 199,386 belong to the 42,155 genes written as `ncRNA_gene`, and the file source annotates none of them. Release 116.3 adds `ncRNA_gene` to the list.

**What vep-rs does instead.** vep-rs reads an `ncRNA_gene` line as a gene, as it reads `gene` and `pseudogene`, and the cache it builds carries the gene's transcripts. A record on a non-coding gene is annotated against them, as Ensembl VEP annotates it from its own cache or from the same genes in a GTF, and as release 116.3 annotates it from the GFF3. Every other type outside the list is skipped as the file source skips it, the skipped types logged once with their counts.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `ensembl-io` and `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/GFF.pm`, lines 71-116, `%INCLUDE_FEATURE_TYPES`: the 35 types the `--gff` source reads, `ncRNA_gene` not among them. Release 116.3 adds it, at line 84 of that release's file.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, lines 207-260, `_get_records_by_coords`: a line is kept only when the list names its type (line 224); any other is dropped with the warning `Ignoring '<type>' feature_type from <file> GFF/GTF file` (line 242).
- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, lines 363-437, `_get_parent_child_structure`: a record whose `Parent` names a line that was not read is an orphan (lines 388-398 and 406-422), reported once as `Parent entries with the following IDs were not found or skipped due to invalid types` (lines 424-427). Every record with a parent is then removed from the top level (lines 430-435), so the transcripts of a skipped gene build nothing.
- `Bio/EnsEMBL/Variation/VariationFeature.pm`, lines 822-840, `get_IntergenicVariation`: a record that overlaps no transcript is written as `intergenic_variant`.
- `Bio/EnsEMBL/IO/Translator/Feature.pm`, lines 198-203 and 469-493: the GFF3 writer types a line by its biotype's Sequence Ontology term, `protein_coding_gene` written as `gene`; the term of every non-coding gene biotype is `ncRNA_gene`.

**What correcting it would change.** With `ncRNA_gene` on the list, as release 116.3 has it, the gene line is read and its transcripts keep their parent. A `--gff` run then writes the rows on them that a `--gtf` run over the same genes and Ensembl VEP's own cache write, and the rows below would agree.

**The class it produces.** `file_cache_reads_ncrna_gene_lines`: the rows a cache vep-rs built from an Ensembl GFF3 writes on the transcripts of an `ncRNA_gene` line, and the `intergenic_variant` row Ensembl VEP 116.2's `--gff` source writes where those transcripts were a record's only features. The class is documented only.

The committed caches hold 313 transcripts built from a file whose biotype an Ensembl GFF3 writes under an `ncRNA_gene` line. 23 are in the cache built from a GTF, 266 in the cache built from a GFF3 and 24 in a mitochondrial shard built from the GFF3 that a third corpus splices into a converted cache. The one committed corpus whose reference Ensembl VEP wrote from a file read a GTF, which types every gene `gene`, and the corpora built from a GFF3 compare with Ensembl VEP reading its own cache, where the engines agree.

**How often the shape occurs.** Tier 2: the shape is a non-coding transcript of a cache built from a GTF or GFF3, met where Ensembl VEP reads a GFF3 as its transcript source. The committed caches hold 77,601 transcripts, 313 of them non-coding transcripts built from a file, in three corpora; no committed flag set compares with Ensembl VEP reading a GFF3.

## 1. LINC00308 and a TMPRSS3 transcript, as the GFF3 writes them

**In short.** A two-gene file holds the GRCh38 release-116 GFF3 lines of LINC00308, a lncRNA gene typed `ncRNA_gene` with the transcripts `ENST00000419431` and `ENST00000657104`, and of TMPRSS3, a protein-coding gene typed `gene`, with its transcript `ENST00000398397`. On the record `linc00308_exon_snv`, `21 22098664 T>C`, in the first exon of both LINC00308 transcripts, vep-rs annotating the cache built from the file writes two `non_coding_transcript_exon_variant` rows; Ensembl VEP 116.2 reading the file through `--gff` writes one `intergenic_variant` row. On the TMPRSS3 record `1995820`, `21 42395417 T>C`, both engines write `start_lost`.

### 1.1 What a reader needs to know first

A GFF3 line's third column is its feature type, a Sequence Ontology term. An Ensembl GFF3 types a gene line `gene`, `ncRNA_gene` or `pseudogene` by its biotype, and a transcript line by its biotype's term, `mRNA` for a protein-coding transcript, `lnc_RNA` for a lncRNA; a transcript's `Parent` attribute names its gene and an exon's names its transcript. Ensembl VEP's `--gff` reads such a file in place of its cache; vep-rs builds a cache from the file first ([A GTF or GFF3 as the transcript model](../cache-from-gtf-or-gff3.md)). `non_coding_transcript_exon_variant` is the term for a variant inside an exon of a transcript without a coding sequence.

### 1.2 The record

The gene and transcript lines of the file, the gene lines less their `description` and `logic_name` attributes, each transcript's exon lines following as in the release file (four exons for `ENST00000657104`, six for `ENST00000419431`, nine exons and eight CDS lines for `ENST00000398397`):

```
21	havana	ncRNA_gene	22098616	22116605	.	+	.	ID=gene:ENSG00000184856;Name=LINC00308;biotype=lncRNA;gene_id=ENSG00000184856;version=7
21	havana_tagene	lnc_RNA	22098616	22116605	.	+	.	ID=transcript:ENST00000657104;Parent=gene:ENSG00000184856;Name=LINC00308-202;biotype=lncRNA;tag=gencode_primary,Ensembl_canonical;transcript_id=ENST00000657104;version=1
21	havana	lnc_RNA	22098617	22116528	.	+	.	ID=transcript:ENST00000419431;Parent=gene:ENSG00000184856;Name=LINC00308-201;biotype=lncRNA;tag=gencode_basic;transcript_id=ENST00000419431;transcript_support_level=5;version=2
21	ensembl_havana	gene	42371840	42396122	.	-	.	ID=gene:ENSG00000160183;Name=TMPRSS3;biotype=protein_coding;gene_id=ENSG00000160183;version=19
21	havana	mRNA	42381876	42396122	.	-	.	ID=transcript:ENST00000398397;Parent=gene:ENSG00000160183;Name=TMPRSS3-201;biotype=protein_coding;ccdsid=CCDS42939.1;tag=gencode_basic;transcript_id=ENST00000398397;transcript_support_level=1;version=4
```

The input:

```
21	22098664	linc00308_exon_snv	T	C	.	.	.
21	42395417	1995820	T	C	.	.	.
```

Both engines run `--symbol --biotype --canonical --mane --tsl --protein --protein_version --gene_version --transcript_version`. Ensembl VEP reads the file, bgzipped and indexed as `linc00308_tmprss3.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.2:

```
linc00308_exon_snv	21:22098664	C	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
1995820	21:42395417	C	ENSG00000160183.19	ENST00000398397.4	Transcript	start_lost	233	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=1;ENSP=ENSP00000381434.3;SOURCE=linc00308_tmprss3.gff3.gz
```

vep-rs:

```
linc00308_exon_snv	21:22098664	C	ENSG00000184856.7	ENST00000419431.2	Transcript	non_coding_transcript_exon_variant	48	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SYMBOL=LINC00308;BIOTYPE=lncRNA;TSL=5
linc00308_exon_snv	21:22098664	C	ENSG00000184856.7	ENST00000657104.1	Transcript	non_coding_transcript_exon_variant	49	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SYMBOL=LINC00308;BIOTYPE=lncRNA;CANONICAL=YES
1995820	21:42395417	C	ENSG00000160183.19	ENST00000398397.4	Transcript	start_lost	233	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=1;ENSP=ENSP00000381434.3
```

The TMPRSS3 rows differ only in the `SOURCE` item, the class [`file_cache_writes_no_source_field`](file-cache-writes-no-source-field.md). The tab output carries the same rows in its columns. In the VCF output Ensembl VEP's `CSQ` for the first record is one entry, `C|intergenic_variant|MODIFIER|` followed by empty fields, and vep-rs's two entries name the transcripts; in the JSON output Ensembl VEP's record carries an `intergenic_consequences` list and vep-rs's a `transcript_consequences` list of two. Ensembl VEP's warnings file reports the skipped line, `Ignoring 'ncRNA_gene' feature_type from linc00308_tmprss3.gff3.gz GFF/GTF file`, and the orphaned transcripts, `Parent entries with the following IDs were not found or skipped due to invalid types: gene:ENSG00000184856`.

### 1.4 Why Ensembl VEP writes that

`_get_records_by_coords` keeps a line only when `%INCLUDE_FEATURE_TYPES` names its type. `gene`, `mRNA`, `lnc_RNA`, `exon` and `CDS` are named, `ncRNA_gene` is not, so the LINC00308 gene line is dropped with the `Ignoring` warning and the two `lnc_RNA` lines and their exons are read. `_get_parent_child_structure` then finds no top-level record for `gene:ENSG00000184856`, lists the identifier in the orphan warning, and removes every record that names a parent from the top level, the two transcripts with the rest. Nothing remains at the record's position, so `get_IntergenicVariation` writes `intergenic_variant`. The TMPRSS3 gene line is typed `gene`, which the list names, so its transcript is built and annotated.

### 1.5 Why vep-rs's row is the correct one

The file states the gene, its two transcripts and their exons, and the record lies inside the first exon of both, at cDNA positions 48 and 49. Ensembl VEP writes the two `non_coding_transcript_exon_variant` rows itself when it reads the same two transcripts from a GTF, where the gene line is typed `gene`, and when it reads them from its own cache. The type `ncRNA_gene` is the Sequence Ontology term for a non-coding gene, written by Ensembl's own GFF3 writer for every such gene; it changes nothing about the transcripts under it. vep-rs's rows are the rows Ensembl VEP writes for the record from any source that reads the gene.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 linc00308_tmprss3.gff3 --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
vep -i records.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --symbol --biotype --canonical --mane --tsl --protein --protein_version \
  --gene_version --transcript_version
```

vep-rs writes the two LINC00308 rows. Ensembl VEP 116.2 with `--gff linc00308_tmprss3.gff3.gz`, the file bgzipped and indexed, writes `intergenic_variant` for the record. Retyping the gene line's third column from `ncRNA_gene` to `gene` makes Ensembl VEP 116.2 write vep-rs's two rows, and release 116.3 writes them from the file as it stands.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/file_cache_reads_ncrna_gene_lines.toml` holds the file, the two records and both engines' rows for `linc00308_exon_snv`. The class covers the rows a cache built from an Ensembl GFF3 writes on the transcripts of an `ncRNA_gene` line where Ensembl VEP 116.2's `--gff` source writes none, and the `intergenic_variant` row that source writes where those transcripts were the record's only features. No committed dataset compares a GFF3-built cache with Ensembl VEP reading the GFF3, so the class is documented only.
