# The order of a transcript's rows when a file is read beside a cache

**The defect in one sentence.** A transcript that a cache and a GTF or GFF3 read beside it both hold, under one identifier, has two rows, and Ensembl VEP orders them by a count of the transcripts it has built from the file. vep-rs writes the cache's row first, then each file's row in argument order.

**Why this is a defect.** The record and the transcript determine both rows, and the file's row carries the file's name in `SOURCE`. Their order in Ensembl VEP's output comes from neither. Ensembl VEP numbers the transcripts it builds from a file, counting from zero, and sorts a record's rows by the text of the identifier and that number. The number follows the records the input held before and, under `--fork`, the buffer size. Where the input reaches two genes in one 1 Mb window of the file, it also follows Perl's hash order. So an unrelated record moves the two rows, and two runs of one command can differ byte for byte.

**What vep-rs does instead.** vep-rs writes the cache's row of an identifier first, then each file's row in argument order: the `--custom` strings in the order given, then `--gff`'s file, then `--gtf`'s. That is the order of Ensembl VEP's own list of sources. It depends on the command line alone, so every run of an input writes the same bytes. The two rows themselves are the ones Ensembl VEP writes. [A file beside a cache in one run](../cache-from-gtf-or-gff3.md#a-file-beside-a-cache-in-one-run) describes such a run.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, line 489, `_create_transcript`: `dbID => $self->{_tr_dbID}++`, one counter per source, advanced per transcript built. Lines 840-847: a transcript is built when a record of the input first reaches it. Line 320: the genes of a window are walked as the `values` of a hash, each gene's transcripts in an order the file's lines fix (lines 329-338 and 385-417).
- `Bio/EnsEMBL/VEP/AnnotationSource.pm`, lines 219-263: the windows are gathered in the order of the input's records. `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, lines 101-113: the transcripts are built in the order of that list, and a transcript that no record reaches is skipped.
- `Bio/EnsEMBL/VEP/Runner.pm`, lines 387-388 and 482-528: under `--fork` a buffer is annotated in child processes, a chunk each, and the count starts again in every child.
- `Bio/EnsEMBL/Variation/BaseVariationFeature.pm`, lines 310-321, `_get_transcript_key`: the identifier, then `_` and the number when the number is true, so a count of zero leaves the bare identifier. `VariationFeature.pm`, line 899: the row is stored under that key. Line 666: the rows are returned in `sort keys` order, a text sort.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 486: the rows are written from that list, and no later step sorts them.

**What correcting it would change.** Keying a file's transcript by its source and its identifier, in place of a count of the transcripts built, would write a record's rows in one order for every input and every run, as vep-rs's are. Ensembl VEP's own order of two such rows would change for some inputs.

**The class it produces.** `file_beside_cache_same_id_row_order`: two `Transcript` rows of one record, allele and transcript identifier, one read from the cache and one from a file (JSON `transcript_consequences` entries), whose fields agree between the engines and whose order differs. The class is documented only: no comparison counts rows under it. The two committed flag sets that carry such rows are compared without their order ([1.7](#17-how-the-comparator-treats-this-class)). A row changed, missing or added is not the class and stays a difference.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md). The shape is a transcript of a gene-model file whose identifier the cache also holds, which a run reaches when it names the file with `--gtf`, `--gff` or `--custom file=<file>,format=gff|gtf` beside a cache. Of the 964 feature lines of the committed GTF `chr21_genes.gtf.gz`, 56 are such transcripts. Five committed flag sets read a gene-model file beside a cache, all in the `GRCh38-gene-model-beside-cache` corpus, the one corpus that carries the shape. Two of the five read that GTF, and in each the 74 records give 862 pairs of rows that share a record, an allele and an identifier.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. `ENST00000452969` read from the cache and from a GTF

**In short.** The record `21 43783348 tmem97p1_exon_snv` lies in the one exon of `ENST00000452969`, which the cache and the corpus's GTF both hold. Both engines write the transcript twice: a row from the cache and a row with `SOURCE=chr21_genes.gtf.gz`. Run alone, Ensembl VEP 116.3 writes the file's row first. Behind two earlier records of a sorted input, it writes the cache's row first. vep-rs writes the cache's row first in both runs.

### 1.1 What a reader needs to know first

`--gtf`, `--gff` and `--custom file=<file>,format=gff|gtf` read a gene-model file as a transcript source. Beside a cache, a run annotates the transcripts of both, and a row of a file's transcript carries the file's short name in `SOURCE`. A file that describes a gene of the cache under the cache's own identifiers gives each such transcript one row per source.

Ensembl VEP reads its input a buffer at a time, `--buffer_size` records, and reads a gene-model file in windows of 1 Mb. Under `--fork`, `--fork 1` included, it annotates a buffer in child processes.

A Perl hash holds its keys in no order of its own. `values` and `keys` return them in an order set by the hash function and its seed, and Perl draws a new seed each time it starts unless the environment variable `PERL_HASH_SEED` fixes one. Perl also perturbs the order as keys are stored, unless `PERL_PERTURB_KEYS=0` switches that off. So the order is the same from run to run only with both variables set.

### 1.2 The record

The record is the 55th of the `GRCh38-gene-model-beside-cache` golden corpus (`tests/golden/116/GRCh38-gene-model-beside-cache/variants.vcf`), a single-base substitution written from the reference:

```
21 43783348 tmem97p1_exon_snv A G
```

`ENST00000452969` is the one transcript of the processed pseudogene TMEM97P1 in the corpus's GTF, `chr21_genes.gtf.gz`. Its one exon spans `21:43783123-43783573` on the forward strand, in the file and in the cache.

Both engines run `--gtf <corpus>/chr21_genes.gtf.gz` beside the cache, with the corpus's FASTA, in the default format. Ensembl VEP runs with `--fork 1` and `--buffer_size 5000`.

### 1.3 What each engine writes

Ensembl VEP 116.3, the record alone:

```
tmem97p1_exon_snv	21:43783348	G	ENSG00000234030	ENST00000452969	Transcript	non_coding_transcript_exon_variant	226	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=chr21_genes.gtf.gz
tmem97p1_exon_snv	21:43783348	G	ENSG00000234030	ENST00000452969	Transcript	non_coding_transcript_exon_variant	226	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
tmem97p1_exon_snv	21:43783348	G	ENSG00000310147	ENST00000847659	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=3543;STRAND=1
```

vep-rs, the record alone:

```
tmem97p1_exon_snv	21:43783348	G	ENSG00000234030	ENST00000452969	Transcript	non_coding_transcript_exon_variant	226	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
tmem97p1_exon_snv	21:43783348	G	ENSG00000234030	ENST00000452969	Transcript	non_coding_transcript_exon_variant	226	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=chr21_genes.gtf.gz
tmem97p1_exon_snv	21:43783348	G	ENSG00000310147	ENST00000847659	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=3543;STRAND=1
```

The two outputs hold the same three rows. Ensembl VEP wrote the file's row first on each of six runs of this input, three with `PERL_HASH_SEED` and `PERL_PERTURB_KEYS` set and three with them unset.

Two earlier records of the corpus, `ankrd30bp1_down_snv` at `21:13381749` and `linc00308_up_snv` at `21:22097116`, lie 30.4 Mb and 21.7 Mb away in other genes. With both before the record, in a sorted input, Ensembl VEP writes the cache's row first, as vep-rs does in both runs. The record's own rows are unchanged. With `--buffer_size 1` under `--fork 1`, Ensembl VEP writes the file's row first again on that three-record input.

The order also changes between runs of one command. An input of `linc00308_up_snv`, `mir6070_exon_snv` and the record gives `ENST00000613967` (MIR6070) and `ENST00000452969` two rows each. With the two environment variables unset, ten runs of Ensembl VEP were made. Five wrote the file's row first on `ENST00000452969` and the cache's row first on `ENST00000613967`. The other five wrote each pair the other way. vep-rs wrote the cache's row first on both, every time.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP stores a record's transcript rows in a hash and writes them in the text order of its keys. A cache transcript's key is its identifier, `_` and its database identifier: `ENST00000452969_281439`. A file's transcript has no database identifier, so the source gives it a number, the count of the transcripts it built before this one. The key drops a number of zero.

Alone, the record reaches one transcript of the file, the first the source builds. Its key is the bare `ENST00000452969`, which sorts before `ENST00000452969_281439`, so the file's row comes first.

The two earlier records reach one transcript of ANKRD30BP1 and two of LINC00308 in the file, and the source builds them first. The record's file transcript is the fourth, with the key `ENST00000452969_3`. As text, `3` follows the `2` that opens `281439`, so that key sorts after `ENST00000452969_281439` and the cache's row comes first. Under `--fork` each chunk of a buffer is annotated in a child process, where the count starts again. With `--buffer_size 1` the record is alone in its process, and its file transcript is the first built once more.

MIR6070 and TMEM97P1 lie in one 1 Mb window of the file. Ensembl VEP walks the genes of a window as the values of a hash, so either gene's transcript can be built before the other's. After the two transcripts of LINC00308, the two take the numbers 2 and 3. `_2` sorts before `_281602` and `_281439`, the endings of the cache copies' keys, and `_3` after them. The gene Perl's hash order lists first has its file row first, and the other its cache row.

Source (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `BaseGXF.pm`, line 489: the counter. Line 320: the walk over the hash's values. `Runner.pm`, lines 482-528: the child processes of `--fork`.
- `BaseVariationFeature.pm`, lines 310-321: the key. `VariationFeature.pm`, line 666: the text sort.

### 1.5 Why vep-rs's row order is the correct one

The record's coordinates determine both rows. The exon's first base is 43783123, so the record is base 43783348 - 43783123 + 1 = 226 of the transcript. A pseudogene's transcript has no coding sequence, so each copy's row reads `non_coding_transcript_exon_variant` at cDNA position 226 with `STRAND=1`, and the file's copy carries the file's name. Both engines print those two rows.

Only the order differs, and Ensembl VEP's own runs disagree about it. The same record, transcript, cache, file and flags give one order alone and the other behind two unrelated records. A buffer size changes it under `--fork`, and two runs of one command order `ENST00000452969`'s rows differently. So neither order can be read from the record, the transcript or the command line.

vep-rs's order is the order in which Ensembl VEP lists a run's sources: the cache, then each file as the command line gives it (`AnnotationSourceAdaptor.pm`, lines 94-105; `Config.pm`, lines 531-543). It is the same for every input and on every run, so two runs of one input can be compared byte for byte.

### 1.6 Try it yourself

The command as vep-rs takes it, with the record of 1.2 in `record.vcf` after the two header lines of [Reproducing these examples](reproducing.md#writing-a-probe-input):

```
vep -i record.vcf -o out.txt --offline \
  --json_cache <corpus>/json_cache \
  --fasta <corpus>/reference.fa.gz \
  --species homo_sapiens --assembly GRCh38 \
  --gtf <corpus>/chr21_genes.gtf.gz
```

`<corpus>` is `tests/golden/116/GRCh38-gene-model-beside-cache`. vep-rs writes the row without `SOURCE` first. Ensembl VEP 116.3 takes the same `--gtf` flag with the cache options of [Reproducing these examples](reproducing.md), the corpus's FASTA decompressed, `--fork 1` and `--buffer_size 5000`. It writes the row with `SOURCE=chr21_genes.gtf.gz` first. With the lines `21 13381749 ankrd30bp1_down_snv` and `21 22097116 linc00308_up_snv` of the corpus's `variants.vcf` placed before the record, it writes that row second.

With the lines `21 22097116 linc00308_up_snv` and `21 43609938 mir6070_exon_snv` before the record, and `PERL_HASH_SEED` and `PERL_PERTURB_KEYS` unset, repeated runs write that row first on some runs and second on others.

### 1.7 How the comparator treats this class

No comparison counts rows under `file_beside_cache_same_id_row_order`. The golden harness and the format comparator check that the entries both engines share come in one order, outside a flag set marked `unordered`. The corpus's two flag sets whose file shares identifiers with the cache, `gtf` and `two-files`, are marked so. Their entries pair by what they hold, the cache's row with the cache's and the file's row with the file's, and every field compares exactly.

A row changed, missing or added stays a difference under the mark. The corpus's three other flag sets read a file whose identifiers the cache lacks, and they compare in the order written.
