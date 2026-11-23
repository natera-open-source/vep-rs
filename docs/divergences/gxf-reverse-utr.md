# The 3' UTR of a reverse-strand transcript read from a GTF or GFF file

**The defect in one sentence.** When Ensembl VEP reads its transcripts from a GTF or GFF file through `--gtf` or `--gff`, it reads a reverse-strand transcript's 3' UTR as the reverse complement of the span mirrored about the transcript. The stop codon a frameshift or a lost stop runs into then lies at the wrong distance, and the count after `fsTer` or `extTer` in the HGVSp is wrong.

**Why this is a defect.** The UTR is the transcript's own sequence past its stop codon, and the file names the exons that hold it. Ensembl VEP reads those exons correctly for every other purpose on the row: the consequence terms, the codons, the amino acids and the HGVSc agree with vep-rs. One routine reads them again for the HGVSp's stop distance, and that read maps the coordinates of a reverse-strand slice from the wrong end. A cache transcript carries its UTR from the cache and never reaches the routine, so the same record annotated against the Ensembl cache gets the right count.

**What vep-rs does instead.** vep-rs builds a cache from the file ([Building a cache from a GTF or GFF3 file](../cache-from-gtf-or-gff3.md)) and reads a transcript's UTR once, from the exons the file names on the transcript's own strand, as it reads an Ensembl cache's. The stop distance is counted along that sequence.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, `ensembl` at release/116):

- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 1101-1113, `_generic_utr`: when the transcript carries no cached UTR, the transcript is transferred onto its feature slice and the UTR is read from the copy.
- `Bio/EnsEMBL/Exon.pm`, line 985: `transfer` drops the exon's cached bases, so the copy's exons read their sequence from the slice again.
- `Bio/EnsEMBL/Variation/Utils/FastaSequence.pm`, lines 369-371: the FASTA-backed `subseq` maps a slice-relative coordinate to `slice start + coordinate - 1` whatever the slice's strand, and multiplies the strand. The core API counts a reverse-strand slice's coordinates from its genomic end, so on such a slice this read is mirrored about the slice.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, line 2492, `_stop_loss_extra_AA`, and line 2333, `_get_fs_peptides`: the extension and frameshift counts read the UTR through `_three_prime_utr`.

**What correcting it would change.** The class vanishes: the HGVSp of every reverse-strand frameshift and stop-loss row of a `--gtf` or `--gff` run carries the count the transcript's sequence gives, as the same transcript's row from the Ensembl cache does.

**The class it produces.** `gxf_reverse_strand_utr_sequence`, applied by the format comparator to the `HGVSp` field of a `STRAND` -1 row of a transcript Ensembl VEP read from the file, where the two values differ in the digits after `fsTer` or `extTer` alone. The concordance figures score consequence sets and never see the field, so nothing is set aside in them.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. The stop distance of a reverse-strand transcript's HGVSp

**In short.** Under `--gtf` or `--gff` with `--hgvs`, a frameshift or a lost stop on a reverse-strand transcript gets an HGVSp whose `fsTer` or `extTer` count Ensembl VEP reads from a UTR it has misread. vep-rs counts along the transcript's own UTR: on the record `1032152`, Ensembl VEP writes `p.Ter452ArgextTer71` and vep-rs `p.Ter452ArgextTer9`. On the GRCh38 chromosome 21 ClinVar records against the release 116 GTF, and against the RefSeq GFF3 (the record under it is in [1.6](#16-try-it-yourself)), every row that differs lies on a reverse-strand transcript and differs in that count alone. Nothing else on the row differs, and the same record against the Ensembl cache agrees in both engines.

### 1.1 What a reader needs to know first

A frameshift's HGVSp, `p.Lys34ArgfsTer41`, names the first changed residue and the distance, in codons, to the first stop of the new reading frame. A lost stop's, `p.Ter452ArgextTer9`, names the residue the stop becomes and the distance to the next stop. When the new stop lies past the end of the coding sequence, the engine reads on into the 3' UTR, so the count depends on the UTR's bases.

`--gtf` and `--gff` give Ensembl VEP its transcripts from a file instead of a cache: the exons' coordinates from the file, their bases from `--fasta`. vep-rs reads such a file through a cache built from it with the cache builder, so both engines annotate the same transcript model. The UTR of a cache transcript is stored in the cache by the dumper; the UTR of a file transcript is read at run time.

### 1.2 The record

The record is `21 42372764 1032152 A G` from `variants.vcf` of the `GRCh38-gtf-builder` golden corpus (`tests/golden/116/GRCh38-gtf-builder/`), a ClinVar GRCh38 substitution in the stop codon of TMPRSS3, a protein-coding gene on the reverse strand:

```
21	42372764	1032152	A	G	.	.	AF_ESP=0.00077;AF_EXAC=0.00024;AF_TGP=0.001;ALLELEID=1022044;CLNDISDB=MONDO:MONDO:0010987,MedGen:C1832827,OMIM:601072,Orphanet:90636|MedGen:CN169374|MedGen:C3661900;CLNDN=Autosomal_recessive_nonsyndromic_hearing_loss_8|not_specified|not_provided;CLNHGVS=NC_000021.9:g.42372764A>G;CLNREVSTAT=criteria_provided,_conflicting_classifications;CLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Uncertain_significance(3)|Likely_benign(1);CLNSIGSCV=SCV001756215|SCV004021199|SCV004267804|SCV004810384;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=ClinGen:CA10042220;GENEINFO=TMPRSS3:64699;MC=SO:0001578|stop_lost;ORIGIN=1;RS=149919890
```

The corpus's `chr21_genes.gtf.gz` carries the Ensembl 116 GTF lines of the gene. The transcript `ENST00000398405.6` ends its coding sequence at the stop codon `21:42372762-42372764` (reverse strand, `TGA` read along the transcript) and its last exon spans `21:42371840-42372779`, so its 3' UTR is `21:42371840-42372761`, 922 bases. Both engines run the corpus's flags with `--hgvs` added:

```
--gtf chr21_genes.gtf.gz --fasta reference.fa --symbol --biotype --canonical --mane --tsl --protein --protein_version --gene_version --transcript_version --hgvs
```

### 1.3 What each engine writes

The transcript's row in the default format. Ensembl VEP:

```
1032152	21:42372764	G	ENSG00000160183	ENST00000398405.6	Transcript	stop_lost	1944	1354	452	*/R	Tga/Cga	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000381442.1;SOURCE=chr21_genes.gtf.gz;HGVSc=ENST00000398405.6:c.1354T>C;HGVSp=ENSP00000381442.1:p.Ter452ArgextTer71
```

vep-rs:

```
1032152	21:42372764	G	ENSG00000160183	ENST00000398405.6	Transcript	stop_lost	1944	1354	452	*/R	Tga/Cga	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000381442.1;HGVSc=ENST00000398405.6:c.1354T>C;HGVSp=ENSP00000381442.1:p.Ter452ArgextTer9
```

Every column but the HGVSp's count agrees (`SOURCE` is the file source's own column, which the comparison of a `--gtf` run sets aside on every row). The record's 18 transcript rows all differ this way, with Ensembl VEP's count 18, 56, 71 or 86 by transcript and vep-rs's 9 on each, since every TMPRSS3 transcript shares the 3' UTR downstream of the stop.

The `fsTer` shape is the same read. The hand-written one-base deletion `21 42372766 fs_del_c1351 GT G`, `c.1351del` on the same transcript, in the codon before the stop, runs its new frame through the stop codon into the UTR; under the same flags Ensembl VEP writes:

```
fs_del_c1351	21:42372767	-	ENSG00000160183	ENST00000398405.6	Transcript	frameshift_variant	1941	1351	451	T/X	Acc/cc	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000381442.1;SOURCE=chr21_genes.gtf.gz;HGVSc=ENST00000398405.6:c.1351del;HGVSp=ENSP00000381442.1:p.Thr451ProfsTer19
```

and vep-rs:

```
fs_del_c1351	21:42372767	-	ENSG00000160183	ENST00000398405.6	Transcript	frameshift_variant	1941	1351	451	T/X	Acc/cc	-	IMPACT=HIGH;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;TSL=5;ENSP=ENSP00000381442.1;HGVSc=ENST00000398405.6:c.1351del;HGVSp=ENSP00000381442.1:p.Thr451ProfsTer8
```

From the transcript's UTR the new frame reads `CCT GAA GAG GAA GGG GAC AAG TAG`, a stop at the eighth codon. From the misread sequence it reads `CCT GAA CTA TCT CTT CCG TGG TGC CAT CTA CAT TTT TGG GAC TCG GGA ATT ATG TGA`, nineteen codons. A frameshift whose new frame meets a stop inside the coding sequence, `21 42372777 fs_del_c1340 TC T`, gets `p.Arg447LysfsTer3` from both engines.

### 1.4 Why Ensembl VEP writes that

`_stop_loss_extra_AA` translates the new reading frame from the lost stop into the 3' UTR until it meets a stop codon, and reads the UTR through `_three_prime_utr`. For a file transcript the UTR is not cached, so `_generic_utr` transfers the transcript onto its feature slice, the span `21:42371840-42396007` on the reverse strand, and calls the core API's `three_prime_utr` on the copy. The transfer drops each exon's cached bases, so the copy's exons read their bases from the slice.

That read goes through the FASTA-backed `subseq` of `FastaSequence.pm`. The core API numbers a reverse-strand slice from its genomic end: relative position 1 is `21:42396007`. The FASTA read instead adds the relative position to the slice's genomic start, `21:42371840`, and reverse-complements the result. The last exon, relative positions 23,229 to 24,168 on the reverse slice, is read from `21:42395068-42396007` in place of `21:42371840-42372779`: the exon's span mirrored about the transcript. The UTR Ensembl VEP translates therefore begins `ACTATCTCTTCCGTGGTGCCATCTACATTT`, and its first in-frame stop lies 71 codons on.

Source (`ensembl-variation` at release/116, `ensembl` at release/116):

- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 1101-1113: `_generic_utr`, the transfer and the read from the copy.
- `Bio/EnsEMBL/Exon.pm`, line 985: the cached bases dropped on transfer.
- `Bio/EnsEMBL/Variation/Utils/FastaSequence.pm`, lines 369-371: the slice-relative coordinate mapped from the slice's start.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2333 and 2492: the frameshift and extension counts that read the UTR.

### 1.5 Why vep-rs's row is the correct one

The transcript's 3' UTR is the reverse complement of `21:42371840-42372761`, and it begins `AGAGGAAGGGGACAAGTAGCCACCTGAGTT`. The record turns the stop codon `TGA` into `CGA`, arginine, so translation continues into that UTR: `AGA GGA AGG GGA CAA GTA GCC ACC TGA`, and the ninth codon is `TGA`, a stop. The extension is nine codons, `extTer9`.

Ensembl VEP's own annotation of the same record against the Ensembl cache gives `extTer9` too: the cache carries the UTR the dumper read on the transcript's strand. Its 71 comes from a sequence that is not the transcript's UTR at all, the reverse complement of `21:42395068-42396007`, the transcript's first exon and the start of its first intron, about 22 kb upstream of the stop codon.

### 1.6 Try it yourself

Run the record of 1.2 through both engines over the corpus's GTF and FASTA. vep-rs reads the cache the corpus built from the file; Ensembl VEP reads the file itself with the FASTA decompressed:

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --fasta <corpus>/reference.fa.gz --species homo_sapiens --assembly GRCh38 \
  --symbol --biotype --canonical --mane --tsl --protein --protein_version \
  --gene_version --transcript_version --hgvs
```

| HGVSp on `ENST00000398405.6` | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `--gtf <corpus>/chr21_genes.gtf.gz --fasta reference.fa` | `p.Ter452ArgextTer71` | `p.Ter452ArgextTer9` |
| the Ensembl 116 GRCh38 cache | `p.Ter452ArgextTer9` | `p.Ter452ArgextTer9` |
| `--gff` over the chromosome 21 lines of NCBI's GCF_000001405.40 GFF3, on the three RefSeq models the record lies on | `extTer86` on `NM_024022.4` and `NM_001256317.3`, `extTer15` on `NM_032404.3` | `extTer9` on each |

The bases behind both counts can be read from the FASTA: the reverse complement of `21:42371840-42372761` is the UTR vep-rs and the cache run translate, and the reverse complement of `21:42395068-42396007` is the sequence the `--gtf` run translates.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) applies `gxf_reverse_strand_utr_sequence` when a comparison case names it in `masks`. Under the class a differing `HGVSp` (the JSON `hgvsp`) is accepted on a row whose `STRAND` is -1 and whose `SOURCE` (the JSON `source`) in Ensembl VEP's output names the file. Both values must be present and equal once the digits, or `?`, after `fsTer` or `extTer` are set aside. A difference anywhere else in the notation, on a forward-strand row, or where one engine prints no HGVSp stays a difference.

The `SOURCE` says whether the row's transcript was read from the file or from a cache. Ensembl VEP stamps every transcript it builds from a GTF or GFF file with the file's base name and writes the stamp as the row's `SOURCE` ([The `SOURCE` field of a GTF or GFF3 read as the transcript source](file-cache-writes-no-source-field.md)). A row of a cache transcript carries no `SOURCE`, or the `Ensembl` or `RefSeq` of a merged cache, and its `HGVSp` stays a difference under the class: the cache stores that transcript's UTR, so the misread cannot be the cause of a difference there. A file whose own name is `Ensembl` or `RefSeq` reads as a merged cache's row and stays a difference too.

Each accepted row is counted once under the class, so a comparison case of the `--gtf` or `--gff` route under `--hgvs` over the chromosome 21 ClinVar records reports under the class the number of rows it accepts. No golden corpus runs those routes under `--hgvs`, so no corpus manifest documents a row of the class.
