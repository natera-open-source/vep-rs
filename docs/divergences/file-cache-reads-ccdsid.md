# The CCDS identifier of an Ensembl GFF3 or GTF under `--ccds` and `--pick`

**The defect in one sentence.** Ensembl VEP 116.3 fills a transcript's CCDS identifier from the database's xrefs alone, when it reads the database or builds its cache, so its `--gff` and `--gtf` sources never set one. A transcript built from an Ensembl GFF3 line that carries `ccdsid=CCDS42939.1`, or from a GTF line that carries `ccds_id "CCDS42939"`, therefore prints no `CCDS` under `--ccds` and ranks as a transcript without a CCDS at the `ccds` step of `--pick`. vep-rs reads the identifier from the line into its cache, prints it under `--ccds` and ranks the transcript as one that has a CCDS.

**Why this is a defect.** The Ensembl GFF3 and GTF write each transcript's CCDS identifier on its line, as `ccdsid` and `ccds_id`, from the same xref the database holds, and the file source reads every attribute of the line into the record it builds. It copies the transcript identifier, the version, the biotype, the support level and the tags from those attributes and stops before the CCDS identifier, so a fact the file states is dropped on the way to the output. The `--pick` order then ranks every transcript of the file as lacking a CCDS, and a transcript that has one can lose the pick to a longer one without.

**What vep-rs does instead.** vep-rs reads an Ensembl GFF3 line's `ccdsid` or a GTF line's `ccds_id` into the cache as the transcript's CCDS identifier, the value as the line carries it. It prints the identifier under `--ccds` and ranks the transcript at the `ccds` step of `--pick` as Ensembl VEP ranks the same transcript from its cache. The two files spell the identifier differently: the GFF3 carries the xref's display label, `CCDS42939.1`, the form Ensembl VEP's own cache prints, and the GTF the xref's primary identifier, `CCDS42939`, without the version. A RefSeq GFF3's CCDS is not read into the cache, so a cache built from that file prints no `CCDS`, as the file source does.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl` core and `ensembl-io` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource/Database/Transcript.pm`, lines 616-624, `prefetch_transcript_ids`: the only writer of the transcript's `_ccds` key, the display identifier of the first `CCDS` xref (line 622), or `-` (line 623); the cache dumper stores the key with the transcript, so the cache source carries it.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, lines 632-703, `_add_identifiers`: the file sources copy the gene identifier, the symbol, the support level, the tags and the canonical flag from the record's attributes and read no `ccdsid` or `ccds_id`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1543-1546: `CCDS` is written under `--ccds` when `_ccds` is defined and not `-`, which a file-source transcript never satisfies.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 740: the `ccds` criterion of the pick order ranks a transcript with a set `_ccds` ahead of one without; a file-source transcript always ranks as without.
- `Bio/EnsEMBL/Transcript.pm`, line 3146, `summary_as_hash`: the core API writes `ccdsid` into the GFF3 transcript line from the transcript's CCDS xref, the xref `prefetch_transcript_ids` reads, as the xref's display label.
- `Bio/EnsEMBL/Utils/IO/GTFSerializer.pm`, lines 459-464: the GTF transcript line gets `tag "CCDS"` and one `ccds_id` per CCDS xref, the xref's primary identifier, which carries no version.

**What correcting it would change.** With `_add_identifiers` copying `ccdsid` or `ccds_id` into `_ccds`, a `--gff` or `--gtf` run would print the file's CCDS identifiers under `--ccds` and rank the transcripts that have one ahead at the `ccds` step, as they rank from Ensembl VEP's cache. The rows below would agree.

**The class it produces.** `file_cache_reads_ccdsid`: the `CCDS` value a cache vep-rs built from an Ensembl GFF3 or GTF writes under `--ccds` where Ensembl VEP's `--gff` or `--gtf` source writes none. It also covers the transcript the `ccds` step of `--pick` selects on such a cache where the file source selects by a later criterion. The class is documented only. Of the two committed caches built from a file, the one built from the GFF3 holds 223 transcripts with a CCDS identifier and the one built from the GTF 12. The one committed run that passes `--ccds` over a file-built cache compares it with Ensembl VEP reading its own cache, where the engines agree.

**How often the shape occurs.** Tier 2: the shape is a transcript of a cache built from an Ensembl GFF3 or GTF that carries a CCDS identifier, met under `--ccds` or a `--pick` family flag. The two committed caches built from a file hold 1,045 transcripts, 235 with a CCDS identifier. 1 committed flag set passes `--ccds` over such a cache, and no committed flag set passes `--ccds` or a pick flag over a file Ensembl VEP reads through `--gff` or `--gtf`.

## 1. TMPRSS3-201 under `--ccds`, and two TMPRSS3 transcripts under `--pick`

**In short.** The GRCh38 release-116 GFF3 lines of the gene TMPRSS3 and its transcript `ENST00000398397`, whose line carries `ccdsid=CCDS42939.1`, make a one-transcript file. On the record `1995820`, `21 42395417 T>C`, vep-rs annotating the cache built from that file under `--ccds` writes `CCDS=CCDS42939.1`; Ensembl VEP 116.3 reading the file through `--gff` writes no `CCDS`. With the two transcripts `ENST00000896630`, `ccdsid=CCDS13686.1`, and `ENST00000896631`, which has none, under `--ccds --pick`, vep-rs picks `ENST00000896630` and Ensembl VEP picks `ENST00000896631`.

### 1.1 What a reader needs to know first

`--ccds` adds the transcript's Consensus CDS identifier to the output. `--pick` keeps one row per record, chosen by an ordered list of criteria: MANE Select, MANE Plus Clinical, canonical, APPRIS, support level, biotype, CCDS, consequence rank, length, then source. Two transcripts tied on every criterion before CCDS are separated by whether each has a CCDS identifier; tied there too, the longer translation, then the longer transcript, wins.

### 1.2 The record

The transcript line of the first file, with its gene line above it less its `description` and `logic_name` attributes; its nine exon and eight CDS lines follow as in the release file:

```
21	ensembl_havana	gene	42371840	42396122	.	-	.	ID=gene:ENSG00000160183;Name=TMPRSS3;biotype=protein_coding;gene_id=ENSG00000160183;version=19
21	havana	mRNA	42381876	42396122	.	-	.	ID=transcript:ENST00000398397;Parent=gene:ENSG00000160183;Name=TMPRSS3-201;biotype=protein_coding;ccdsid=CCDS42939.1;tag=gencode_basic;transcript_id=ENST00000398397;transcript_support_level=1;version=4
```

The second file holds the same gene line and the lines of `ENST00000896630`, 2,387 transcript bases with a 1,365-base CDS and `ccdsid=CCDS13686.1`, and of `ENST00000896631`, 2,580 bases with a 1,368-base CDS and no `ccdsid`. Both are `havana_tagene` transcripts spanning `21:42371840-42395962` on the reverse strand with no support level or canonical tag.

The input:

```
21	42395417	1995820	T	C	.	.	.
```

Both engines run `--symbol --biotype --canonical --mane --tsl --protein --protein_version --gene_version --transcript_version --ccds`, and for the second file `--pick` as well. Ensembl VEP reads each file, bgzipped and indexed as `tmprss3_ccds.gff3.gz` and `tmprss3_pick.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 1.3 What each engine writes

On the first file, Ensembl VEP 116.3:

```
1995820	21:42395417	C	ENSG00000160183.19	ENST00000398397.4	Transcript	start_lost	233	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=1;ENSP=ENSP00000381434.3;SOURCE=tmprss3_ccds.gff3.gz
```

vep-rs:

```
1995820	21:42395417	C	ENSG00000160183.19	ENST00000398397.4	Transcript	start_lost	233	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=1;CCDS=CCDS42939.1;ENSP=ENSP00000381434.3
```

The tab output carries a `CCDS` column that is `-` in Ensembl VEP's and `CCDS42939.1` in vep-rs's, and the VCF output a `CCDS` field in the `CSQ` entry, empty in one and filled in the other. The JSON output carries a `ccds` key in vep-rs's consequence entry that Ensembl VEP's lacks.

On the second file under `--pick`, Ensembl VEP 116.3 keeps the longer transcript:

```
1995820	21:42395417	C	ENSG00000160183.19	ENST00000896631.2	Transcript	start_lost	291	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;ENSP=ENSP00000566690.1;SOURCE=tmprss3_pick.gff3.gz
```

vep-rs keeps the transcript with the CCDS identifier:

```
1995820	21:42395417	C	ENSG00000160183.19	ENST00000896630.2	Transcript	start_lost	101	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;CCDS=CCDS13686.1;ENSP=ENSP00000566689.1
```

Without `--pick` both engines write both transcripts' rows, `start_lost` on each, vep-rs's `ENST00000896630` row alone carrying `CCDS=CCDS13686.1`.

The same transcripts' lines of the GRCh38 release-116 GTF part the engines the same way, the identifier spelled as the GTF carries it. On the one-transcript file, whose transcript line carries `ccds_id "CCDS42939"`, Ensembl VEP 116.3 reading it through `--gtf` writes the record's one row with no `CCDS` in four formats, and vep-rs on the cache built from that file writes `CCDS=CCDS42939`. On the two-transcript file under `--ccds --pick`, Ensembl VEP picks `ENST00000896631` and vep-rs `ENST00000896630`, whose row carries `CCDS=CCDS13686`; without `--pick` both engines write both transcripts' rows.

### 1.4 Why Ensembl VEP writes that

The `--gff` source builds each transcript in `_create_transcript` and names it in `_add_identifiers`, which reads the gene identifier, the symbol, the support level and the tags from the line's attributes and nothing else, so the transcript has no `_ccds` key. `OutputFactory.pm` writes `CCDS` only from that key, and ranks the `ccds` criterion of `--pick` from it, so both transcripts of the second file rank as lacking a CCDS, the tie passes to the consequence rank, equal, and then to length, where `ENST00000896631` wins.

### 1.5 Why vep-rs's row is the correct one

The file states `ccdsid=CCDS42939.1` on the transcript's line, written by Ensembl's own API from the xref Ensembl VEP's cache carries as the transcript's CCDS. The same transcript from Ensembl VEP's cache prints that identifier under `--ccds` and ranks ahead of a transcript without one at the `ccds` step. vep-rs's rows are the rows Ensembl VEP writes for these records from its cache.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 tmprss3_ccds.gff3 --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
vep -i record.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --symbol --biotype --canonical --mane --tsl --protein --protein_version \
  --gene_version --transcript_version --ccds
```

vep-rs writes `CCDS=CCDS42939.1`. Ensembl VEP 116.3 with `--gff tmprss3_ccds.gff3.gz`, the file bgzipped and indexed, writes no `CCDS`. Removing `ccdsid=CCDS42939.1;` from the transcript line makes the engines agree, and running Ensembl VEP from its cache instead of the file makes it print the identifier. The transcript's GTF lines given as `--gtf` to the builder and to Ensembl VEP show the same pair, vep-rs writing `CCDS=CCDS42939`.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/file_cache_reads_ccdsid.toml` holds the first file, the record and both engines' `CCDS` values. The class covers the `CCDS` value a cache built from an Ensembl GFF3 or GTF writes under `--ccds` where the `--gff` or `--gtf` source writes none. It also covers the row the `ccds` step of `--pick` selects on such a cache where the file source selects by length or source. No committed dataset reaches the shape, so the class is documented only.
