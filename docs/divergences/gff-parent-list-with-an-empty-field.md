# A GFF3 transcript whose `Parent` list opens with an empty field

**The defect in one sentence.** Ensembl VEP 116.2's `--gff` source links a record to every set identifier in its `Parent` list but keeps the record at the top level whenever the list's first field is empty. A transcript line written `Parent=,gene:X` is therefore visited twice, as the gene's child and as a record of its own, and its gene identifier, symbol and gene version are those of whichever visit comes last in Perl's hash order.

**Why this is a defect.** The line names its gene, and the source reads the name: the link step splits the list on commas, skips the empty field and files the record under `gene:X`. The prune step that decides which records stand at the top level then reads the first field alone, so the same record is both a child and a top-level record.

Each top-level record is handed to the transcript list, and a child's gene is set on every visit, so the gene the transcript ends with depends on the order a Perl hash yields its values. On the record below that order gives the gene under seventeen of twenty-one hash seeds and no gene under the other four. Output that changes with a hash seed is not a reading of the file.

**What vep-rs does instead.** vep-rs reads a transcript line whose `Parent` list opens with an empty field as the transcript of the gene line a later field names, once, and gives it that gene's identifier, symbol and version. A list that opens with an empty field and names no gene line of the file leaves the record a record of its own, as it is in Ensembl VEP, and a list that opens with an empty field and names nothing set, `Parent=,0`, likewise. A GTF transcript line with `gene_id ",X"` prints `,X` as its gene in both engines, since the printed identifier is the line's own attribute, so only the GFF3 form shows the difference.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/GFF.pm`, lines 242-256, `_record_get_parent_ids`: `split(',', $parents)` keeps a leading empty field, so `Parent=,gene:X` reads as the two fields `''` and `gene:X`.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, lines 385-401, `_get_parent_child_structure`: line 389 skips an empty parent, line 391 files the record under `gene:X`, and line 400 enters every record into the top-level hash.
- `BaseGXF.pm`, lines 430-435: the prune deletes a record from the top-level hash only when its parent list is non-empty and its first field is not the empty string, so the record stays.
- `BaseGXF.pm`, lines 311-345, `_create_transcripts`: line 320 walks the top-level hash in hash order; a gene record hands over its children with that gene (lines 325-330 and 336), and a non-gene record hands over itself with no gene (lines 331-333 and 336). The same record hash is handed over twice, and line 336 overwrites its gene on the second visit.
- `BaseGXF.pm`, lines 863-875, `merge_features`: the transcript list is made unique by each record's MD5, so the record is annotated once, with the gene of the last visit.
- `BaseGXF.pm`, lines 840-847, `lazy_load_transcript`, and lines 632-703, `_add_identifiers`: the gene identifier comes from the gene record's `gene_id` attribute (lines 651-653), else from its `ID` with the `gene:` prefix removed (lines 655-658), and the symbol from its `Name` (lines 669-671). With no gene record, line 657 substitutes on an undefined value, which Perl reports on standard error, and both stay unset.

**What correcting it would change.** With the prune reading the same set fields as the link step, the record would leave the top-level hash and be visited once, as its gene's child. The row would then carry the gene's identifier, symbol and version under every hash seed, as vep-rs's row does.

**The class it produces.** `gff_parent_list_with_an_empty_field`: the `Gene`, `SYMBOL` and gene-version values, and the JSON `gene_id` and `gene_symbol`, of a transcript whose GFF3 line's `Parent` list opens with an empty field and names a gene line of the file. Ensembl VEP's `--gff` source writes the gene's values or leaves them empty by hash order; a cache vep-rs built from the file writes the gene's. The class is documented only: no Ensembl or RefSeq file writes such a list, and the committed GTF carries none.

**How often the shape occurs.** Tier 3: the shape is a GFF3 or GTF line whose parent list opens with an empty field, which no release file carries. The 964 feature lines of the committed GTF carry none, 14 committed flag sets run a cache built from a GTF or GFF3, and no committed flag set runs Ensembl VEP over a file with the shape.

## 1. A TMPRSS3 transcript whose `Parent` list opens with an empty field

**In short.** The GRCh38 release-116 GFF3 lines of the gene TMPRSS3 and its transcript `ENST00001113783`, with the transcript's `Parent=gene:ENSG00000160183` rewritten `Parent=,gene:ENSG00000160183`, make a one-transcript file. On the record `1995820`, `21 42395417 T>C`, a `start_lost` on that transcript, vep-rs writes `ENSG00000160183.19` and `SYMBOL=TMPRSS3`. Ensembl VEP 116.2 reading the file through `--gff` writes the same row under one hash seed and a row with `-` for the gene and no `SYMBOL` under another, with a Perl warning about an undefined value.

### 1.1 What a reader needs to know first

A GFF3 record's `Parent` attribute lists the identifiers of the records it belongs to, comma-separated; an exon shared by two transcripts lists both. Ensembl VEP's `--gff` source builds a tree from those lists when a record overlaps the lines, then makes a transcript from each tree node under a gene, or from each node with no parent whose type is not a gene's. A Perl hash yields its keys in an order that depends on the process's hash seed, which Perl draws at random unless `PERL_HASH_SEED` fixes it.

### 1.2 The record

The transcript line, with its gene line above it (its `description` and `logic_name` left out here) and its three exon and two CDS lines below (omitted here), as in the release file but for the empty first field of `Parent`:

```
21	ensembl_havana	gene	42371840	42396122	.	-	.	ID=gene:ENSG00000160183;Name=TMPRSS3;biotype=protein_coding;gene_id=ENSG00000160183;version=19
21	havana_tagene	mRNA	42393803	42396122	.	-	.	ID=transcript:ENST00001113783;Parent=,gene:ENSG00000160183;Name=TMPRSS3-224;biotype=protein_coding;tag=gencode_basic;transcript_id=ENST00001113783;version=1
```

The input:

```
21	42395417	1995820	T	C	.	.	.
```

Both engines run `--symbol --biotype --canonical --mane --tsl --protein --protein_version --gene_version --transcript_version`. Ensembl VEP reads the file, bgzipped and indexed as `tmprss3_parent.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.2 under `PERL_HASH_SEED=0`, and under sixteen of the twenty other seeds from 1 to 20:

```
1995820	21:42395417	C	ENSG00000160183.19	ENST00001113783.1	Transcript	start_lost	233	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;ENSP=ENSP00000783588.1;SOURCE=tmprss3_parent.gff3.gz
```

Ensembl VEP 116.2 under `PERL_HASH_SEED=11`, and under seeds 12, 17 and 19:

```
1995820	21:42395417	C	-	ENST00001113783.1	Transcript	start_lost	233	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;BIOTYPE=protein_coding;ENSP=ENSP00000783588.1;SOURCE=tmprss3_parent.gff3.gz
```

with, on its standard error, Perl's own warning for line 657 (under `--fork` the warnings file carries it instead, behind the child's process id):

```
Use of uninitialized value in substitution (s///) at /opt/vep/src/ensembl-vep/modules/Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm line 657, <__ANONIO__> line 1.
```

vep-rs:

```
1995820	21:42395417	C	ENSG00000160183.19	ENST00001113783.1	Transcript	start_lost	233	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;ENSP=ENSP00000783588.1
```

In the VCF output the second seed's `CSQ` entry has empty `SYMBOL` and `Gene` fields, and in the JSON output its `gene_id` and `gene_symbol` are `null`; the first seed's and vep-rs's carry the gene's values. Beyond the gene's values, the rows differ only in what the class [`file_cache_writes_no_source_field`](file-cache-writes-no-source-field.md) covers, the `SOURCE` item above: the consequence, the positions, the codons and the protein are the same in every run.

### 1.4 Why Ensembl VEP writes that

The link step files the transcript under `gene:ENSG00000160183` and enters it into the top-level hash; the prune keeps it there because its first `Parent` field is empty. `_create_transcripts` then meets the transcript twice: once among the gene's children, where its gene record is set to the gene, and once as a top-level record of a non-gene type, where its gene record is set to nothing. The hash order decides which visit is the last, and `merge_features` keeps one copy with that visit's gene record. With no gene record, `_add_identifiers` finds no gene identifier on the transcript line, substitutes on an undefined value at line 657 and leaves the identifier, the symbol and the version unset.

### 1.5 Why vep-rs's row is the correct one

The transcript line names one gene, `gene:ENSG00000160183`, and the file holds that gene's line with its identifier, symbol and version. The empty field before it names nothing and is skipped as no parent by the link step itself. A reading that files the record under the gene it names, once, and prints that gene's values is the reading of the file as written; a result that is the gene's values under one hash seed and empty under another is not.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 tmprss3_parent.gff3 --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
vep -i record.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --symbol --biotype --canonical --mane --tsl --protein --protein_version \
  --gene_version --transcript_version
```

vep-rs writes `ENSG00000160183.19` and `SYMBOL=TMPRSS3`. Ensembl VEP 116.2 with `--gff tmprss3_parent.gff3.gz`, the file bgzipped and indexed, writes the same under `PERL_HASH_SEED=0` and the gene-less row under `PERL_HASH_SEED=11`. Restoring `Parent=gene:ENSG00000160183` makes the engines agree under every seed.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/gff_parent_list_with_an_empty_field.toml` holds the file, the record and both engines' rows, Ensembl VEP's under the seed that leaves the gene empty. The class covers the `Gene`, `SYMBOL` and gene-version values, and the JSON `gene_id` and `gene_symbol`, of a transcript whose GFF3 `Parent` list opens with an empty field and names a gene line. No committed dataset reaches the shape, so the class is documented only.
