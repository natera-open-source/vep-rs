# The `SOURCE` field of a GTF or GFF3 read as the transcript source

**The defect in one sentence.** Ensembl VEP's `--gtf` and `--gff` register the file as a custom annotation source and stamp every transcript read from it with the file's base name. That name is printed as a `SOURCE` field on every transcript row, with a `SOURCE` column, an empty column named after the file and header lines describing both, none of which a run from a cache writes. A cache vep-rs built from the same file names no source, so vep-rs writes none of them.

**Why this is a defect.** The `SOURCE` field exists for a merged cache, where it tells an Ensembl transcript from a RefSeq one. A GTF or GFF3 read as the transcript source carries one set of transcripts, so the field repeats the file name on every row, in every format. The custom column the registration adds is empty on every row, because a transcript source annotates no record. The same genes read from Ensembl VEP's own cache print no `SOURCE`, no file column and no header line for either, so a user comparing the two routes finds a field, a column and two header lines the file's content never asked for.

**What vep-rs does instead.** A cache vep-rs builds from a GTF or GFF3 records no source name, and a run on it writes `SOURCE` only from a merged cache's `Ensembl` or `RefSeq` stamp, as Ensembl VEP's own cache run does. The default output's Extra column, the tab columns, the VCF `CSQ` fields, the JSON consequence keys and the header lines are those of a cache run.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/Config.pm`, lines 531-543: `--gff <file>` and `--gtf <file>` set `--custom file=<file>,format=gff` or `format=gtf`, so the file is a custom annotation source.
- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 154-158: a custom source's short name is the last component of its path, the file's base name.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, line 635, `_add_identifiers`: every transcript built from the file gets the source's short name as its `_source_cache`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1583: `SOURCE` is written from `_source_cache` whenever it is defined, flag or no flag.
- `Bio/EnsEMBL/VEP/Constants.pm`, line 101: a run with a custom source lists `SOURCE` among its output fields (line 98 lists it for `--merged`), described as `Source of transcript` (line 167).
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2435-2453, `get_custom_headers`: each custom source adds a header line `## <short name> : [PATH]/<file>` and a column of its own. `OutputFactory/Tab.pm`, lines 187-202, and `OutputFactory/VCF.pm`, lines 276 and 449, append that column to the tab columns, the `CSQ` fields and the VCF `##INFO` lines. A transcript source adds no custom annotation to any record (lines 1004-1012 read `_custom_annotations`, which `BaseGXF.pm` never fills), so the column is empty.
- `Bio/EnsEMBL/VEP/OutputFactory/JSON.pm`, line 345: the JSON output lowers the key, so each transcript consequence carries `"source": "<file>"`.

**What correcting it would change.** With the transcript source stamping no `_source_cache` and registering no output column, a `--gtf` or `--gff` run would write the rows, columns and header lines of a cache run, and the rows below would agree.

**The class it produces.** `file_cache_writes_no_source_field`: the `SOURCE` field Ensembl VEP's `--gtf` or `--gff` source writes on every transcript row, the empty column named after the file, and the header lines describing both, which a run on a cache vep-rs built from the file writes none of. The field is the default format's `SOURCE=<file>` item, the tab column, the `CSQ` field and the JSON `source` key. The class is coded: the golden harness sets the field, the column and the header lines aside before comparing a run on a file-built cache with Ensembl VEP reading the file.

The committed caches hold 1,082 transcripts built from a file (of 77,675): 1,045 in the two caches built whole from a GTF or a GFF3, and 37 in a mitochondrial shard built from the GFF3 that a third corpus splices into a converted cache. The one corpus whose reference read the file, through `--gtf`, reaches the class on every transcript row of its two runs, 862 and 691 rows in each format.

**How often the shape occurs.** Tier 2: the shape is a transcript of a cache built from a GTF or GFF3, met where Ensembl VEP reads the file as its transcript source. The committed caches hold 77,675 transcripts, 1,082 of them built from a file, in three corpora; 2 committed flag sets compare such a cache with Ensembl VEP reading the file.

## 1. The eight genes of a GTF, read from the file and from a cache built from it

**In short.** The record `1995820`, `21 42395417 T>C`, lies in the first codon of TMPRSS3's transcript `ENST00000398397`. Ensembl VEP 116.2 reading the GTF `chr21_genes.gtf.gz` through `--gtf` writes `start_lost` with `SOURCE=chr21_genes.gtf.gz` as the row's last item; vep-rs annotating the cache built from that GTF writes the same row without the item. The file is the committed corpus `tests/golden/116/GRCh38-gtf-builder`, eight genes of chromosome 21 with 59 transcripts, and the difference is on every one of the 862 transcript rows its 74 records produce.

### 1.1 What a reader needs to know first

`--gtf` and `--gff` make Ensembl VEP read its transcripts from a file in place of a cache; vep-rs builds a cache from the file first ([A GTF or GFF3 as the transcript model](../cache-from-gtf-or-gff3.md)). Ensembl VEP's custom annotation sources (`--custom`) each add a named column to the output and a header line naming the file; `SOURCE` is the field a merged cache fills with `Ensembl` or `RefSeq` to tell the two transcript sets apart. The default output's Extra column is a `;`-joined list of `key=value` items; the tab output gives each field a column; the VCF output lists the fields in the `CSQ` entry's `Format` and the JSON output lowers each field's name to a key.

### 1.2 The record

The GTF's gene and transcript lines for `ENST00000398397`, its exon and CDS lines following as in the release file:

```
21	ensembl_havana	gene	42371840	42396122	.	-	.	gene_id "ENSG00000160183"; gene_version "19"; gene_name "TMPRSS3"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
21	havana	transcript	42381876	42396122	.	-	.	gene_id "ENSG00000160183"; gene_version "19"; transcript_id "ENST00000398397"; transcript_version "4"; gene_name "TMPRSS3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TMPRSS3-201"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS42939"; tag "gencode_basic"; transcript_support_level "1";
```

The input:

```
21	42395417	1995820	T	C	.	.	.
```

Both engines run `--symbol --biotype --canonical --mane --tsl --protein --protein_version --gene_version --transcript_version`. Ensembl VEP reads `chr21_genes.gtf.gz`, the corpus's bgzipped and indexed GTF, through `--gtf` with the corpus FASTA; vep-rs reads the corpus's committed cache, built from that GTF and that FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.2:

```
1995820	21:42395417	C	ENSG00000160183	ENST00000398397.4	Transcript	start_lost	233	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=1;ENSP=ENSP00000381434.1;SOURCE=chr21_genes.gtf.gz
```

vep-rs:

```
1995820	21:42395417	C	ENSG00000160183	ENST00000398397.4	Transcript	start_lost	233	1	1	M/V	Atg/Gtg	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=1;ENSP=ENSP00000381434.1
```

Ensembl VEP's default output also carries the header lines `## SOURCE : Source of transcript` and `## chr21_genes.gtf.gz : [PATH]/chr21_genes.gtf.gz`; vep-rs's carries neither.

Its tab output has 29 columns, `SOURCE` the 28th, filled on every transcript row, and `chr21_genes.gtf.gz` the 29th, `-` on every row; vep-rs's has the other 27. Its VCF `CSQ` entry has 31 fields, `SOURCE` the 30th and `chr21_genes.gtf.gz` the 31st, empty, with a header line `##INFO=<ID=chr21_genes.gtf.gz,Number=.,Type=String,Description="[PATH]/chr21_genes.gtf.gz">`; vep-rs's entry has 29 fields and no such line. Its JSON transcript consequences each carry `"source": "chr21_genes.gtf.gz"`; vep-rs's carry no `source` key. On the corpus's 74 records the field is on all 862 transcript rows of the default run and on the 691 of the run under `--transcript_filter`, whose 23 `intergenic_variant` rows carry none; every other value agrees.

### 1.4 Why Ensembl VEP writes that

`--gtf` is an option set for `--custom file=chr21_genes.gtf.gz,format=gtf`, so the file is a custom source whose short name is its base name. `_add_identifiers` gives every transcript it builds that name as `_source_cache`, the key a merged cache sets to `Ensembl` or `RefSeq`, and `OutputFactory.pm` prints `SOURCE` whenever the key is defined. Because a custom source is loaded, `SOURCE` joins the output fields, and `get_custom_headers` adds the source's own column and header line, as it does for a `--custom` VCF or BED; a transcript source fills no custom annotation, so its column stays empty.

### 1.5 Why vep-rs's row is the correct one

The row's annotation is the transcript's, and both engines agree on every value of it. `SOURCE` names where a transcript came from when a run has two transcript sets to tell apart, and a GTF or GFF3 read alone has one. Ensembl VEP reading the same transcripts from its cache writes no `SOURCE`, no file column and no header line for either. vep-rs's row is the row Ensembl VEP writes for the record from its cache.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline \
  --json_cache <corpus>/json_cache --fasta <corpus>/reference.fa.gz \
  --assembly GRCh38 --symbol --biotype --canonical --mane --tsl --protein \
  --protein_version --gene_version --transcript_version
```

with `<corpus>` the committed `tests/golden/116/GRCh38-gtf-builder`. vep-rs writes the row without a `SOURCE` item. Ensembl VEP 116.2 with `--gtf <corpus>/chr21_genes.gtf.gz` in place of the cache writes `SOURCE=chr21_genes.gtf.gz` on the row and the two header lines above it; with `--tab`, `--vcf` or `--json` it writes the column, the `CSQ` field and the `source` key. A GFF3 read through `--gff` writes the same field with its own base name.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/file_cache_writes_no_source_field.toml` holds the record and both engines' `SOURCE` values on the committed corpus. The class covers the `SOURCE` field, the column named after the file and the header lines describing both, which Ensembl VEP's `--gtf` or `--gff` source writes and a run on a file-built cache does not. The golden harness sets them aside before comparing a run against a reference Ensembl VEP wrote from a file, in each format: the default format's `SOURCE=` item, the tab columns, the `CSQ` fields with their `Format` names and `##INFO` line, and the JSON `source` key.

A corpus with such a reference names the class in its manifest's `masks`, so the comparison states what it sets aside. The two runs of the corpus `GRCh38-gtf-builder` reach the class on every transcript row.
