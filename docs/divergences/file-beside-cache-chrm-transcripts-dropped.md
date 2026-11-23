# A gene-model file's `chrM` transcripts beside a cache

**The defect in one sentence.** Beside a cache, Ensembl VEP writes no row for a transcript of a GTF or GFF3 whose lines name the mitochondrial sequence `chrM`, and vep-rs writes the transcript's rows, with the file's name in `SOURCE`.

**Why this is a defect.** The synonyms table of each human cache of release 116 links `MT` and `chrM`. Ensembl VEP's documentation of `--synonyms` says synonyms let the input and any annotation source, a GFF among them, name a chromosome differently. Ensembl VEP follows the link to find the file's lines for an input record on `MT`, and reads them. It then looks each transcript up under the file's name less its `chr`, `M`, a name the synonyms table does not hold. The transcript is dropped without a warning. Ensembl VEP compares no sequence in doing so, and it annotates the same lines written `MT`.

**What vep-rs does instead.** vep-rs reads a file's `chrM` as the cache's `MT` and annotates the file's transcripts beside the cache's. In the runs of this page each row is the one Ensembl VEP writes when the file names the sequence `MT`, but for the file's name in `SOURCE`. On GRCh37 alone there are two bounds. Over a FASTA whose mitochondrial record has another length than the assembly's 16,569 bases, vep-rs leaves the file's `chrM` transcripts out with a warning, as [A file beside a cache in one run](../cache-from-gtf-or-gff3.md#a-file-beside-a-cache-in-one-run) describes. It does the same for a mitochondrial region that has tested coding transcripts, none of which reads as a coding sequence on the FASTA record it is cut from.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-io` at release/116):

- `Bio/EnsEMBL/IO/Parser/GXF.pm`, lines 88-94, `get_seqname`: a line's sequence name with a leading `chr` removed. `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, line 282: that name becomes the chromosome of the transcript read from the line.
- `Bio/EnsEMBL/VEP/AnnotationSource.pm`, lines 219-263: each source tables the lowest and highest position of a buffer's input records by chromosome. `InputBuffer.pm`, lines 463-498, `min_max`: the buffer keeps the first table it is handed. `AnnotationSourceAdaptor.pm`, lines 94-105: the cache's transcript source is listed before a file's, and the sort at line 98 orders on a key neither sets. So the table is keyed by the cache's name.
- `AnnotationSource.pm`, lines 288-311, `_check_overlap`, called through lines 143 and 331-341: a transcript is kept only where its chromosome is a key of that table or resolves to one. `BaseVEP.pm`, lines 689-737, `get_source_chr_name`: a name resolves in one step, to a synonym of it or to the name with `chr` added or removed.

**What correcting it would change.** Looking a file's transcript up under the name its line writes would, in Ensembl VEP, resolve `chrM` to `MT` through the cache's synonym. Input records on the mitochondrial sequence would gain the rows of such a file's transcripts in Ensembl VEP's output.

**The class it produces.** `file_beside_cache_chrm_transcripts_dropped`: the `Transcript` rows vep-rs alone writes for an input record on the mitochondrial sequence, each with a gene-model file's name in `SOURCE`, where the file's lines write `chrM`, the run reads it beside a cache and vep-rs reads its `chrM` as `MT`. The class is documented only. No comparison sets its rows aside, because a row does not say how its file names the sequence. A comparison that meets such a row counts it as a difference, as it counts any row one engine alone writes.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md). The shape is a line of a gene-model file, read beside a cache, that names its sequence `chrM`. Five committed flag sets read a gene-model file beside a cache, all in the `GRCh38-gene-model-beside-cache` corpus. Their two files, `chr21_genes.gtf.gz` and `novel_models.gff3.gz`, hold 1,808 lines outside their comments, each on chromosome 21, so no committed corpus carries the shape. The made-up file of section 1 is a test fixture that no corpus reads. GENCODE's release 19 GTF, the GRCh37 gene set, writes its 37 mitochondrial genes on `chrM`.

Every run of this page names its file with `--gtf` or `--gff`, and the runs of section 1 are on GRCh37. The same input record and file read beside the GRCh38 cache give the same outputs: the row of `MITOT1` in vep-rs, none in Ensembl VEP 116.3, and the row in both where the file's lines name `MT`.

## 1. A made-up gene on `chrM` read beside the GRCh37 cache

**In short.** The input record `MT 150 chrm_gene_exon_snv` lies in the one exon of the made-up transcript `MITOT1`, which a GTF writes on `chrM`. Read beside the GRCh37 cache, vep-rs writes the row of `MITOT1` and the ten rows of the cache's transcripts. Ensembl VEP 116.3 writes the ten rows alone. It writes the row of `MITOT1` when the file's lines name `MT`, or when a synonyms table links `MT` with `M`.

### 1.1 What a reader needs to know first

`--gtf`, `--gff` and `--custom file=<file>,format=gff|gtf` read a gene-model file, a GTF or GFF3 of genes and transcripts, as a transcript source. Named beside a cache, the file adds its transcripts to the cache's: a run annotates both, and a row of a file's transcript carries the file's short name in `SOURCE`. Both engines need a FASTA for such a file.

A cache directory holds `chr_synonyms.txt`, a table of pairs of names for one sequence, and a run on the cache reads it. The table of the GRCh37 human cache of release 116 links `MT` with `chrM` and with `NC_012920.1`. It holds no `M`.

Ensembl's own files name the human mitochondrial sequence `MT`. GENCODE's release 19 GTF writes `chrM`.

### 1.2 The record

A single-base substitution, written from the made-up reference:

```
MT 150 chrm_gene_exon_snv G A
```

The gene-model file is `chrm_gene_synthetic.gtf.gz` in `crates/vep-cli/tests/fixtures/gene-model-chrm-beside-cache`. It holds one gene, `MITOG1`, with one transcript, `MITOT1`, of the biotype `Mt_tRNA`. The one exon spans 101 to 200 on the forward strand, and the file's three lines write the sequence `chrM`.

vep-rs reads the committed GRCh37 cache slice `crates/vep-cli/tests/fixtures/bayesdel-mitochondrial-contig/json_cache`, which holds `MT` and the synonyms of 1.1. Ensembl VEP reads its GRCh37 cache of release 116. The FASTA is `reference_synthetic.fa` in `crates/vep-cli/tests/fixtures/mitochondrial-frame-check`: one made-up record named `MT`, of the assembly's 16,569 bases. Both engines run `--gtf <file>` beside the cache, in the default format.

### 1.3 What each engine writes

vep-rs writes eleven rows: ten of the cache's transcripts and this one, whose eighth cell is the cDNA position.

```
chrm_gene_exon_snv	MT:150	A	MITOG1	MITOT1	Transcript	non_coding_transcript_exon_variant	50	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=chrm_gene_synthetic.gtf.gz
```

Ensembl VEP 116.3 writes ten rows, each of a transcript of the cache and each equal to vep-rs's row of that transcript. It writes no row of `MITOT1`, exits 0 and writes no warning.

Three more runs of Ensembl VEP on the input record change one thing each:

| The file's lines name | Synonyms | Row of `MITOT1` |
| --- | --- | --- |
| `chrM` | the cache's table, named again with `--synonyms` | none |
| `chrM` | the cache's table and a line `MT M`, with `--synonyms` | written |
| `MT` | the cache's table | written |

Each written row is vep-rs's row, but for the name of the file read in `SOURCE`. A GFF3 of the same gene read with `--gff` behaves the same way in Ensembl VEP: no row where its lines name `chrM`, and the row where they name `MT`. vep-rs writes the row for the GFF3 whose lines name `chrM`.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP annotates a buffer of input records one source after another, the cache's first. Each source tables the lowest and highest input record position on each chromosome, under its own name for the chromosome, and hands the table to the buffer. The buffer keeps the first table, the cache's, whose one key here is `MT`.

The file's source resolves the input record's `MT` to the file's `chrM` through the synonym and reads the three lines. Its parser returns each line's sequence name with a leading `chr` removed, so the transcript built from the lines carries the chromosome `M`.

Before annotation, the source keeps a transcript only where its chromosome is a key of the buffer's table or resolves to one. A name that is no key is resolved in one step: to a synonym of it that is a key, or else to the name with `chr` added or removed where that is a key. `M` is not `MT`, the synonyms pair `MT` with `chrM` and not with `M`, and `chrM` is no key. The transcript is dropped, and no message is written.

With the line `MT M` in the synonyms, `M` resolves to `MT` and the transcript is kept. With the file's lines written `MT`, the name is the table's key. By the source, the lookup is the same for any name a file writes. Of the names run here, `chrM` alone has a form without `chr` that the cache neither holds nor lists as a synonym.

Source (`ensembl-vep` at release/116.3, `ensembl-io` at release/116):

- `AnnotationSourceAdaptor.pm`, lines 94-105: the order of the sources. `AnnotationSource.pm`, lines 243-247 and 260: the table and its hand-over. `InputBuffer.pm`, lines 466-467: the first table is kept.
- `AnnotationSource.pm`, line 179, and `BaseGXF.pm`, lines 207-260: the lines are found and read. `Parser/GXF.pm`, line 92, and `BaseGXF.pm`, line 282: the name less `chr`.
- `AnnotationSource.pm`, lines 296-306: the lookup that drops the transcript. `BaseVEP.pm`, lines 706-731: a synonym, or else `chr` added or removed.

### 1.5 Why vep-rs's row is the correct one

The input record's coordinates and the file's lines determine the row. Position 150 lies in the exon that spans 101 to 200 on the forward strand, so it is base 150 - 101 + 1 = 50 of the transcript. The biotype is non-coding, so the term is `non_coding_transcript_exon_variant`, with `STRAND=1` and the file's name in `SOURCE`. vep-rs prints that row.

Ensembl VEP's own run treats the file's `chrM` as the input record's `MT`. The cache's table pairs the two names, and Ensembl VEP follows that pair to read the file's lines for the input record. Its documentation of `--synonyms` names GFF files among the sources a synonym serves. The transcript is then lost on `M`, a name that neither the file, the cache nor the input writes.

Ensembl VEP's other runs agree with the row. It writes the row when the lines name `MT` and when `M` is linked to `MT`, over the same FASTA. So the missing row follows from how a name is spelt, not from the sequence or the transcript.

Ensembl VEP's reader of custom annotation files allows for the shortened name. It compares a line's name as written as well, and its comment names files that write `chr1` where the cache's synonyms hold no such name (`AnnotationSource/File.pm`, lines 326-334). The gene-model reader makes no second comparison.

On GRCh37 the name `chrM` alone does not settle the sequence, since it also names the 16,571-base mitochondrial sequence of the UCSC hg19 naming. A file's transcripts take their bases from the run's FASTA. So vep-rs does not read a file's `chrM` as `MT` there where that FASTA's mitochondrial record has another length than the assembly's 16,569 bases. This run's has that length.

A published gene set shows the same. Take the mitochondrial lines of GENCODE's release 19 GTF, a file the repository does not hold, read beside the GRCh37 cache over Ensembl's GRCh37 FASTA record `MT`. The input record `MT 3308` of the committed corpus `tests/golden/116/GRCh37` gets 20 rows of the file's transcripts in vep-rs, one for each transcript of the file that the input record lies in or within 5,000 bases of, and none in Ensembl VEP 116.3. With `M` linked to `MT`, Ensembl VEP writes the same 20 rows.

### 1.6 Try it yourself

The command as vep-rs takes it, with the record of 1.2 in `record.vcf` after the two header lines of [Reproducing these examples](reproducing.md#writing-a-probe-input):

```
vep -i record.vcf -o out.txt --offline \
  --json_cache <fixtures>/bayesdel-mitochondrial-contig/json_cache \
  --fasta <fixtures>/mitochondrial-frame-check/reference_synthetic.fa \
  --species homo_sapiens --assembly GRCh37 \
  --gtf <fixtures>/gene-model-chrm-beside-cache/chrm_gene_synthetic.gtf.gz
```

`<fixtures>` is `crates/vep-cli/tests/fixtures`. vep-rs writes the row of `MITOT1`. Ensembl VEP 116.3 takes the same `--gtf` and `--fasta` flags with the cache options of [Reproducing these examples](reproducing.md) and writes none. Its row appears with `--synonyms` naming a copy of the cache's `chr_synonyms.txt` that ends with the line `MT M`, the two names separated by a tab.

### 1.7 How the comparator treats this class

No comparison sets rows aside under `file_beside_cache_chrm_transcripts_dropped`. No committed corpus reads a gene-model file whose lines name `chrM`, so none holds a row of the class. Where a comparison meets one, the row is one vep-rs alone writes, and it is counted as a difference.
