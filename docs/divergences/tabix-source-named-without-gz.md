# A tabix-indexed source named without `.gz` or `.bgz`

**The defect in one sentence.** Ensembl VEP 116.2 opens every `--custom` source of the vcf, bed, gff and gtf formats, and every `--regulatory_gff` file, through a tabix reader that tests the file's name before it opens the file. A name that does not end in `.gz` or `.bgz` ends the run with `ERROR: Input file is not bgzipped, cannot use tabix`, whatever the file holds, so a bgzipped and indexed source named `sites.vcf`, `sites.vcf.GZ` or `regions.bed` is never read. vep-rs opens the same file by its content and annotates the run.

**Why this is a defect.** Whether a file is bgzipped is a property of its bytes: a BGZF file opens with the gzip magic bytes and the BGZF extra field, and its tabix index sits beside it under the file's own name with `.tbi` or `.csi` appended. The name is a convention. The reader's own library, htslib, opens a file by its content and finds the index by the name it is given, so the file the test refuses is one the library would read. The test is also case-sensitive, so `.GZ` fails it, and the message blames the file's compression, which is not at fault.

**What vep-rs does instead.** vep-rs opens a `--custom` or `--regulatory_gff` source by content: it reads the BGZF blocks and finds the index beside the file (`<file>.csi`, `<file>.tbi`, or either on the file's stem) whatever the file is called. A file that is not bgzipped, or has no index beside it, exits 1 with an error naming the file, as [Custom annotation](../custom-annotation.md) states.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `ensembl-io` at release/116):

- `Bio/EnsEMBL/IO/TabixParser.pm`, line 52, `open`: `confess "ERROR: Input file is not bgzipped, cannot use tabix\n" unless $filename =~ /\.b?gz$/;` runs before line 63 hands the file to `Bio::DB::HTS::Tabix->new`. The pattern matches a lower-case `.gz` or `.bgz` at the end of the name and nothing else.
- `Bio/EnsEMBL/IO/Parser/VCF4Tabix.pm`, line 53, `Parser/BedTabix.pm`, line 50, `Parser/GFF3Tabix.pm`, line 50, and `Parser/GTFTabix.pm`, line 50: each format's parser calls that `open` first, so every tabix format takes the test.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/VCF.pm`, lines 152-155, `File/BED.pm`, line 82, `File/GFF.pm`, line 136, and `File/GTF.pm`, line 101: a `--custom` source of each format opens through the matching parser. `File/RegFeat.pm`, lines 204-206: a `--regulatory_gff` file (`file=`, `motifs=` or `emars=`) opens through `GFF3Tabix`.
- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 369-372, `valid_chromosomes`, reached from `BaseRunner.pm`, line 453, while `Runner.pm`, line 823, builds the input parser: the first call to `parser`, so the run ends before the output opens.
- `File.pm`, lines 183-188: a `format=bigwig` source opens through `Bio::DB::BigFile`, not a tabix parser, so the test does not reach it.

**What correcting it would change.** With the test removed, or replaced by a read of the file's first bytes, Ensembl VEP would read the files vep-rs reads and both engines would write the same rows. A file that is not bgzipped would fail inside the library with its own message.

**The class it produces.** `tabix_source_named_without_gz`: a run whose `--custom` source of a tabix format, or whose `--regulatory_gff` file, is bgzipped and indexed under a name not ending in `.gz` or `.bgz`. vep-rs annotates the run; Ensembl VEP writes no output. The class is documented only: every committed source is named `.gz` or `.bgz`, so no manifest names the class and no comparator applies it.

## 1. A VCF source copied under a name without `.gz`

**In short.** The `sites.vcf.gz` source of the GRCh38 custom corpus and its `.tbi` index are copied to `sites_plain.vcf` and `sites_plain.vcf.tbi`, byte for byte. One record of the corpus is run under `--custom` with the copy as an `exact` source named `sites`. vep-rs writes the same two rows it writes for the original, `sites=multi_1` with the source's `AF`, `CLASS` and `FILTER` on the matched allele. Ensembl VEP 116.2 writes no output and ends with `ERROR: Input file is not bgzipped, cannot use tabix`.

### 1.1 What a reader needs to know first

A `--custom` VCF source annotates a variant with the rows that overlap it, or under `type=exact` the rows whose alleles match it. Each engine asks the source's tabix index for the rows over the variant's window, so the source must be bgzipped and indexed. Both engines refuse a plain-text file or a file without an index; the question here is what decides that a file is bgzipped. Ensembl VEP reads the name; vep-rs reads the bytes.

### 1.2 The record

The source is the corpus file `tests/golden/116/GRCh38-custom/custom/sites.vcf.gz` with its index, copied under the new name:

```
cp sites.vcf.gz sites_plain.vcf
cp sites.vcf.gz.tbi sites_plain.vcf.tbi
```

The input, one record of the same corpus, a two-allele insertion against the Ensembl 116 GRCh38 cache:

```
21	5227532	rs1979768284	T	TTTCCTCTCCTCTCTT,TTTCCTCTCTT	.	.	.
```

Both engines run `--custom file=sites_plain.vcf,short_name=sites,format=vcf,type=exact,fields=AF%CLASS`.

### 1.3 What each engine writes

vep-rs, in the default format:

```
rs1979768284	21:5227532-5227533	TTCCTCTCCTCTCTT	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;sites=multi_1;sites_AF=0.003;sites_CLASS=multi_allelic;sites_FILTER=PASS
rs1979768284	21:5227532-5227533	TTCCTCTCTT	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

Ensembl VEP 116.2 writes no output file in any format and prints, with exit status 255:

```
ERROR: Input file is not bgzipped, cannot use tabix
 at .../Bio/EnsEMBL/IO/TabixParser.pm line 52.
```

The tab, VCF and JSON outputs of vep-rs carry the same two rows, and every one is byte for byte the output of the same run on `sites.vcf.gz`, apart from the header lines that name the file. The same pattern holds for the copy named `sites.vcf.GZ` and for the corpus's `regions_chr.bed.gz` copied to `regions_plain.bed` with its index, where vep-rs writes `regions=exact_1,overlap_1` on the record at 21:5028999. It holds for the regulatory corpus's `regulatory_features.gff3.gz` copied to `regulatory_plain.gff3` under `--regulatory_gff file=` too: vep-rs writes the record's two regulatory rows, and Ensembl VEP ends the same way, exit status 2. A copy named `sites.vcf.bgz` is read by both engines.

### 1.4 Why Ensembl VEP writes that

Before the output opens, the runner asks each annotation source for its chromosomes, and a `--custom` VCF source answers by opening its parser. The parser's `open` tests the name against `/\.b?gz$/` and dies when it does not match. `sites_plain.vcf` ends in `.vcf`, so the run dies there; the file's bytes, the index beside it and the record are never read. `sites.vcf.GZ` fails the same test because the pattern is lower-case. A `--regulatory_gff` file opens through the GFF3 parser, which inherits the same `open`.

### 1.5 Why vep-rs's row is the correct one

The copy is the same file as `sites.vcf.gz`: the same BGZF blocks, the same index, the same rows. Nothing about the data changed, so the annotation cannot change. vep-rs reads the blocks and the index and writes the rows the original gives, the rows both engines write for `sites.vcf.gz`. A file that is not bgzipped is still refused, by its bytes.

### 1.6 Try it yourself

```
cp sites.vcf.gz sites_plain.vcf && cp sites.vcf.gz.tbi sites_plain.vcf.tbi
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --custom file=sites_plain.vcf,short_name=sites,format=vcf,type=exact,fields=AF%CLASS
```

vep-rs writes the two rows of 1.3. Ensembl VEP 116.2 ends with the message of 1.3 and writes no output. Renaming the copy and its index to end in `.gz` or `.bgz` makes the engines agree.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/tabix_source_named_without_gz.toml` holds the copy, the record and both engines' results. The class covers a run on a tabix-format `--custom` source or a `--regulatory_gff` file that is bgzipped and indexed under a name not ending in `.gz` or `.bgz`, a run Ensembl VEP does not complete. No committed source is named that way, so no manifest names the class and the class is documented only.
