# `--gencode_basic` and `--gencode_primary` on a cache built from a GTF or GFF3

**The defect in one sentence.** Ensembl VEP 116.3 applies `--gencode_basic` and `--gencode_primary` in one routine, `filter_transcript`, that its cache and database sources call and its `--gff` and `--gtf` sources never do, so under either file source both flags are accepted and change nothing: every transcript of the file is annotated. On a cache built from the same file, vep-rs leaves out every transcript that lacks the flag's attribute.

**Why this is a defect.** The two flags promise a transcript set, the GENCODE basic or primary set, whatever the source of the transcripts. The file source reads the `tag` attribute of each transcript line and stores the `gencode_basic` and `gencode_primary` tags as the same transcript attributes the cache source's filter tests, so the information the filter needs is present on that path. The filter is simply never called there, and no warning says so. The source does apply `--transcript_filter` to every transcript it builds, which shows the step where the set filter belongs.

**What vep-rs does instead.** vep-rs reads a GTF or GFF3 into its cache and annotates against that cache, so the two flags apply as they do on every cache: a transcript without the attribute is left out of the annotation. On an Ensembl GFF3 or GTF line the tag `basic`, the spelling the GRCh37 release-87 files and the GENCODE GTF carry, is read as `gencode_basic`, the attribute Ensembl VEP's own GRCh37 cache carries for those transcripts. Ensembl VEP's file sources store a tag only under the exact spelling `gencode_basic`, which those files never write.

A GTF or GFF3 read beside a cache in the same run is filtered the same way: a transcript of the file without the attribute is left out, where Ensembl VEP's file source annotates it ([A file beside a cache in one run](../cache-from-gtf-or-gff3.md#a-file-beside-a-cache-in-one-run)).

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, lines 215-251, `filter_transcript`: line 223 drops a transcript without a `gencode_basic` attribute under `--gencode_basic`, line 228 one without `gencode_primary` under `--gencode_primary`.
- `Bio/EnsEMBL/VEP/AnnotationSource/Cache/Transcript.pm`, line 241, and `Bio/EnsEMBL/VEP/AnnotationSource/Database/Transcript.pm`, line 251: the only two callers, one per source that reads Ensembl's own transcripts.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, lines 140-164, `get_features_by_regions_uncached`, and lines 311-345, `_create_transcripts`: the file source's path from the records of a region to its transcript list, with no filter call.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, lines 682-697, `_add_identifiers`: the `tag` list is read (lines 683-686) and a `gencode_primary` or `gencode_basic` attribute is added for a tag spelled exactly so (lines 690-696); `basic` matches neither.
- `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, lines 92-175, `annotate_InputBuffer`: line 116 applies the `--transcript_filter` set to each transcript of any source, the cache sources' `filter_transcript` having run earlier on theirs.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1613: the `gencode_primary` attribute the file source stored is read for `--flag_gencode_primary`, so that flag writes its column on the file source while `--gencode_primary` filters nothing.

**What correcting it would change.** With `filter_transcript` called where the file source lists a region's transcripts, a `--gff` or `--gtf` run under either flag would annotate the tagged transcripts alone, as the same run on Ensembl VEP's cache does, and the rows below would agree.

**The class it produces.** `file_cache_honours_gencode_basic`: the rows on a transcript without the `gencode_basic` attribute under `--gencode_basic`, or without `gencode_primary` under `--gencode_primary`, that Ensembl VEP's `--gff` or `--gtf` source writes and a cache vep-rs built from the same file does not. A record left with no transcript takes the intergenic row that every cache gives it. The class is documented only. Two committed corpora annotate a cache built from a file: the GTF corpus and the GFF3 corpus. The two committed runs that pass the flags over such a cache compare it with Ensembl VEP reading its own cache, where the engines agree, so no committed run reaches the shape.

**How often the shape occurs.** Tier 2: the shape is a transcript of a cache built from a GTF or GFF3 that lacks one of the two attributes, met under a non-default flag. The two committed caches built from a file hold 1,045 transcripts, 876 of which lack one attribute or the other: 424 without `gencode_basic` and 865 without `gencode_primary`. Two committed flag sets pass one of the flags over such a cache, and no committed flag set passes either over a file Ensembl VEP reads through `--gff` or `--gtf`.

## 1. A record on two LINC00308 transcripts, one outside the basic set

**In short.** The GTF corpus's transcript model holds two transcripts of the lncRNA gene LINC00308 over `21:22098664`: `ENST00000419431`, tagged `gencode_basic`, and `ENST00000657104`, tagged `gencode_primary` and `Ensembl_canonical` but not `gencode_basic`. On the record `linc00308_exon_snv`, `21 22098664 T>C`, both engines write two `non_coding_transcript_exon_variant` rows without the flag. Under `--gencode_basic`, Ensembl VEP 116.3 reading the GTF through `--gtf` writes the same two rows, and vep-rs annotating the cache built from that GTF writes the `ENST00000419431` row alone.

### 1.1 What a reader needs to know first

Ensembl's GTF marks each transcript's set memberships with repeated `tag` attributes: `gencode_basic` for the GENCODE basic set, `gencode_primary` for the primary set. Ensembl VEP's `--gencode_basic` keeps the transcripts tagged basic, and `--gencode_primary` the transcripts tagged primary; on its cache a record whose every transcript is dropped is written as `intergenic_variant`. The `--gtf` source builds a transcript from the file's lines when a record overlaps it and stores the tags as attributes of that transcript.

### 1.2 The record

The two transcript lines of the corpus's GTF, their exon and CDS lines omitted:

```
21	havana	transcript	22098617	22116528	.	+	.	gene_id "ENSG00000184856"; transcript_id "ENST00000419431"; transcript_version "2"; gene_name "LINC00308"; transcript_biotype "lncRNA"; tag "gencode_basic"; transcript_support_level "5";
21	havana_tagene	transcript	22098616	22116605	.	+	.	gene_id "ENSG00000184856"; transcript_id "ENST00000657104"; transcript_version "1"; gene_name "LINC00308"; transcript_biotype "lncRNA"; tag "gencode_primary"; tag "Ensembl_canonical";
```

The input, one record of the corpus:

```
21	22098664	linc00308_exon_snv	T	C	.	.	.
```

Both engines run the corpus's flags, `--symbol --biotype --canonical --mane --tsl --protein --protein_version --gene_version --transcript_version`, plus `--gencode_basic`. Ensembl VEP reads the GTF through `--gtf` with the corpus's FASTA; vep-rs reads the corpus's committed cache, which was built from that GTF and FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.3 writes two rows, the same two it writes without the flag:

```
linc00308_exon_snv	21:22098664	C	ENSG00000184856	ENST00000419431.2	Transcript	non_coding_transcript_exon_variant	48	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SYMBOL=LINC00308;BIOTYPE=lncRNA;TSL=5;SOURCE=chr21_genes.gtf.gz
linc00308_exon_snv	21:22098664	C	ENSG00000184856	ENST00000657104.1	Transcript	non_coding_transcript_exon_variant	49	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SYMBOL=LINC00308;BIOTYPE=lncRNA;CANONICAL=YES;SOURCE=chr21_genes.gtf.gz
```

vep-rs writes the first transcript's row alone, without the `SOURCE` item, the class [`file_cache_writes_no_source_field`](file-cache-writes-no-source-field.md):

```
linc00308_exon_snv	21:22098664	C	ENSG00000184856	ENST00000419431.2	Transcript	non_coding_transcript_exon_variant	48	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SYMBOL=LINC00308;BIOTYPE=lncRNA;TSL=5
```

The tab output carries the same two rows against one, the VCF output two `CSQ` entries against one, and the JSON output two `transcript_consequences` entries against one. Over the corpus's 74 records in four formats, Ensembl VEP writes 862 rows with the flag as without, its JSON output identical byte for byte. vep-rs writes 862 rows without the flag and 708 with it, dropping the 154 rows on 23 of the 26 transcripts the file leaves untagged. Under `--gencode_primary` Ensembl VEP writes the same 862 rows again and vep-rs 107, on the 9 transcripts tagged `gencode_primary`.

On an Ensembl GFF3 the engines part the same way. On the GRCh37 release-87 chromosome 21 file, whose 2,178 transcripts carry `tag=basic` on 1,217 and `gencode_basic` on none, Ensembl VEP writes 568 rows on 80 records with the flag as without, its JSON output identical. vep-rs on the cache built from that file writes 568 rows without the flag and 320 with it.

The GRCh38 release-116 chromosome 21 slice `crates/vep-cache-builder/tests/fixtures/grch38_chr21/slice.gff3` holds 65 transcript lines, which carry `gencode_basic` on 44 and `gencode_primary` on 30. On the 150 chromosome 21 records of the corpus `tests/golden/116/GRCh38-native-builder`, Ensembl VEP 116.3 writes 186 rows under each of the three flag states, and vep-rs on the cache built from the slice writes 186, 174 and 159. Ensembl VEP 116.2 wrote 173 under each, skipping the slice's `ncRNA_gene` lines, which hold 44 of the transcripts and which Ensembl VEP 116.3 reads.

The GRCh37 release-87 GTF of the same chromosome carries `tag "basic"` on the same 1,217 of its 2,178 transcript lines. On the same 80 records Ensembl VEP reading it through `--gtf` writes 568 rows with the flag as without, and vep-rs on the cache built from the GTF writes 568 and 320, the 320 rows the GFF3 build's.

### 1.4 Why Ensembl VEP writes that

`Runner.pm` builds the `--gtf` source with the flags in its configuration, and the source reads each transcript line's tags into transcript attributes. The cache and database sources pass every transcript through `filter_transcript`, which tests those attributes; the file source lists the transcripts of a region and hands them to `annotate_InputBuffer` without that call. `annotate_InputBuffer` applies the `--transcript_filter` set alone, so `ENST00000657104`, which carries no `gencode_basic` attribute, is annotated like any other transcript, and nothing in the output or the warnings file says that the flag was passed.

### 1.5 Why vep-rs's row is the correct one

The flag asks for the transcripts of the GENCODE basic set, and the file states which transcripts those are: `ENST00000419431` carries the tag and `ENST00000657104` does not. A cache built from the file keeps the file's tags, and the same flag on Ensembl VEP's own cache, where the attribute comes from the same release's annotation, drops the untagged transcript. vep-rs's row is the row Ensembl VEP writes for this record from its cache under the flag.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 \
  --symbol --biotype --canonical --mane --tsl --protein --protein_version \
  --gene_version --transcript_version --gencode_basic
```

vep-rs writes the `ENST00000419431.2` row alone. Ensembl VEP 116.3 with `--gtf <corpus>/chr21_genes.gtf.gz --fasta <corpus>/reference.fa.gz` and the same flags writes both rows. Dropping `--gencode_basic`, or running Ensembl VEP from its cache instead of the GTF, makes the engines agree.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/file_cache_honours_gencode_basic.toml` holds the record, the flags and both engines' rows on the corpus's GTF and cache. The class covers the rows on a transcript without the `gencode_basic` attribute under `--gencode_basic`, or without `gencode_primary` under `--gencode_primary`, that Ensembl VEP's `--gff` or `--gtf` source writes and the cache built from the file does not. No committed dataset reaches the shape, so the class is documented only and the format comparator counts no row under it.
