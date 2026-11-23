# The reference codons of an insertion standing on the bases an edited RefSeq model holds and the genome lacks

**The defect in one sentence.** For an insertion written exactly where an edited RefSeq model holds bases the genome lacks, inside the coding sequence, Ensembl VEP cuts those bases out of the model's sequence before it reads the reference codons. The parts of `Codons` and `Amino_acids` before the slash, their reference halves, then leave out the bases the row's own `CDS_position` names and take in as many bases past the row's last residue. vep-rs reads the reference codons from the model's coding sequence as the cache holds it, over the residues the row prints, and the reference amino acids are their translation.

**Why this is a defect.** An edited model, one the cache corrected against an alignment, can hold bases the genome lacks. The row of an insertion written at such a point prints the span of those added bases as its `CDS_position`, in both engines. The reference codons are the model's bases over the residues the row prints in `Protein_position`, and the model's protein, which the reference amino acids name, is translated from them. Ensembl VEP's reference codons hold no base of the printed span and, in its place, the codons of later residues, so its reference amino acids are not the model's residues there. Its alternate codons end at the row's last residue, so its two halves disagree.

**What vep-rs does instead.** vep-rs reads the reference codons from the model's coding sequence as the cache holds it, over the residues the row prints. The reference amino acids are their translation, the model's residues at those positions. The alternate codons and amino acids and the positions are the ones Ensembl VEP writes. So are the terms, with one exception: where the inserted bases are the span's own with three left out, vep-rs writes `inframe_deletion` ([`refseq_edited_model_shorter_insertion_term`](refseq-edited-model-shorter-insertion-term.md)).

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `ensembl` at release/116, `ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/Mapper.pm`, lines 485-558, `map_insert`: an insertion is mapped through the two bases beside its point (lines 489-491). The first mapped coordinate is then moved to start one past the first base (lines 522-528), and the last to end one before the second (lines 536-542).
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 263-264, `cds_start`: the row's CDS span runs from the first coordinate's start to the last coordinate's end. Where the two bases are a model's added bases apart, the span is those bases.
- `Bio/EnsEMBL/Variation/TranscriptVariation.pm`, lines 441-468, `_get_ref_allele`: under `--use_transcript_ref`, which a cache of corrected models, a BAM-edited cache, switches on (`Bio/EnsEMBL/VEP/CacheDir.pm`, lines 411-416), a reference allele is cut from the model's sequence (lines 462-463). An insertion's is `-` whatever the model holds (line 451).
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 872-873, `codon`: on a model carrying an `_rna_edit` attribute (line 827), the reference allele is written over the CDS span, `substr($cds, $cds_start - 1, $vf_nt_len) = $seq`. `$seq` is empty for `-` (lines 821-822) and `$vf_nt_len` is the span's length (line 851).
- The same routine, lines 841-843 and 877: the reference codons are `$codon_len` bases from the first base of the first printed residue, a length fixed before the write. Line 868 and the second arm of line 877 read the alternate codons from a sequence the write does not touch.
- The same file, lines 724-736, `peptide`: the reference amino acids are the translation of that codon string.
- The write of lines 872-873 is also behind [`refseq_edited_model_deleted_span_codons`](refseq-edited-model-deleted-span-codons.md), the opposite shape. There a deletion spans bases a model's alignment deleted, the reference allele is longer than the CDS span, and the window fills with genomic bases the model lacks. Here the allele is shorter than the span, empty against the added bases, and the window loses the model's own bases.

**What correcting it would change.** A write made only when the reference allele has the span's length would leave the model's sequence whole on these rows. Every such row's reference codons would be the model's bases over the residues it prints, and its reference amino acids the model's residues there, the values vep-rs writes. The alternate codons and amino acids and the three position columns would not change. Neither would the terms, except on a row whose inserted bases are the span's own with three left out, where Ensembl VEP decides the term on these reference codons.

**The class it produces.** `refseq_edited_model_inserted_span_codons`: the reference half of `Codons`, and with it the reference half of `Amino_acids`, on the row of an insertion written exactly where a RefSeq model carrying an `_rna_edit` attribute holds bases the genome lacks, inside its coding sequence. The two engines' reference codons have one length and different bases, all in lower case. The alternate halves, `cDNA_position`, `CDS_position` and `Protein_position` are equal. So are the terms, unless the inserted bases are the span's own with three left out: that row carries `protein_altering_variant` in Ensembl VEP and `inframe_deletion` in vep-rs. The format comparator sets the two reference halves aside where a comparison names the class and gives it the cache vep-rs read (section 1.7).

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an insertion written exactly where a RefSeq model holds bases the genome lacks, inside its coding sequence. This page calls a place where a model holds such bases a point of added bases. The Ensembl caches of release 116 hold none, 0 in the 238,823 transcripts of GRCh37 and 0 in the 669,092 of GRCh38, so the divergence needs a cache that holds RefSeq models.

Of the golden corpora, 7 hold points of added bases on their models. Three hold them at 5 genomic positions of chromosome 21, none inside a coding sequence: `GRCh37-merged` 18, `GRCh37-refseq` 11 and `GRCh38-merged` 4. `GRCm39-mouse-merged` and `GRCm39-mouse-refseq` hold 3 each, `AaegL5-aedes` 6 and `AgamP4-anopheles` 4, all inside the coding sequence: the 16 models the tier counts. No insertion of the 7 corpora lies on a point, as written or after the 3' shift.

Against the RefSeq GRCh38 cache of release 116, the gnomAD v4.1 chromosome 21 dataset has 9 insertions written on a point of added bases and the 1000 Genomes high-coverage chromosome 21 dataset has 1. All ten insertions lie on points outside a coding sequence, so their rows carry no codons.

The genome-wide ClinVar GRCh38 file has 2 insertions written on a point inside a coding sequence of that cache, both at `2:73385942`. Both were run on both engines under `--refseq`. Section 2 works ClinVar variation 1167157, whose row is of this class. The other, an insertion as long as the model's bases there, shows the same two reference halves of `Codons` on a row whose term differs too, which makes it a row of [`refseq_edited_model_equal_span_insertion_term`](refseq-edited-model-equal-span-insertion-term.md). A further 29 reach such a point only after the 3' shift and do not show it, since the two columns are read at the written position.

Against the RefSeq GRCh37 cache of release 116, the genome-wide ClinVar GRCh37 file shows the difference on 3 rows. Two are the rows of the same two records, written there at `2:73613070`. The third is ClinVar variation 3035246 at `2:96780986`, 18 bases inserted on the nine bases the reverse-strand model `NM_000682.7` holds there, its coding bases 903 to 911. Its row differs in the reference half of `Codons` alone.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. One base inserted on the six bases `NM_002111.8` holds and the genome lacks

**In short.** `4 3074935 ins_on_model_bases A AT` on GRCh38 inserts one base between 4:3074935 and 4:3074936. The RefSeq model `NM_002111.8` (HTT, forward strand) holds six bases between those two positions that the genome lacks, its coding bases 111 to 116. Under `--refseq` both engines write `frameshift_variant`, `CDS_position` 111-116 and `Protein_position` 37-39 on the model's row, with the same alternate codons and amino acids. vep-rs writes the reference codons `cagcagcaa`, the model's coding bases 109 to 117, and the reference amino acids `QQQ`, the model's residues 37 to 39. Ensembl VEP 116.3 writes `caacagccg`, coding bases 109, 110 and 117 to 123, and `QQP`.

### 1.1 What a reader needs to know first

Ensembl's RefSeq cache holds RefSeq transcript models placed on the genome. Where a model's own sequence differs from the genome, the cached model is corrected against an alignment of the RefSeq sequences, and each difference is recorded as an `_rna_edit` attribute: a start, an end and the model's bases. An attribute whose end is one below its start adds bases at that point, bases the genome does not have. A model the correction succeeded on prints `BAM_EDIT=OK` on its rows. A run on such a cache switches on `--use_transcript_ref`, which reads a row's reference allele from the model in place of the input ([RefSeq and merged caches](../cache-conversion.md#refseq-and-merged-caches)).

A model's mapper pairs its cDNA positions with genomic positions, one pair per unbroken stretch. Added bases belong to no pair: the pair before them ends at the genomic base before the point, and the pair after them starts at the next genomic base, several cDNA positions later. A cDNA position counts from the model's first base and a CDS position from the first base of its coding sequence. Residue n is translated from CDS positions 3n - 2 to 3n.

The `Codons` column prints the codons of the residues in `Protein_position`, the reference before the slash and the alternate after it, with the changed bases in upper case. `Amino_acids` prints their translations the same way. This page calls the part before the slash the reference half and the part after it the alternate half.

### 1.2 The record

One base inserted between two neighbouring genomic positions, written for this example, against the RefSeq GRCh38 cache of release 116:

```
4	3074935	ins_on_model_bases	A	AT	.	.	.
```

Both engines drop the shared first base, so the record is one base, `T`, inserted between 4:3074935 and 4:3074936.

The cache holds, for `NM_002111.8`, the `_rna_edit` attribute `256 255` with six bases, and two mapper pairs that meet at the point. The first pairs cDNA 1-255 with 4:3074681-3074935 (3074935 - 3074681 + 1 = 255 bases). The second pairs cDNA 262-414 with 4:3074936-3075088 (153 bases on each side). The six bases at cDNA 256-261 lie between the two pairs, and the genome has no base for them.

The coding sequence starts at cDNA 146, so a CDS position is the cDNA position less 145. The six bases are CDS 256 - 145 = 111 to 261 - 145 = 116. Both engines run `--refseq` with `--fasta`, which a run under `--use_transcript_ref` requires.

### 1.3 What each engine writes

The row of `NM_002111.8` in the default format. Ensembl VEP 116.3:

```
ins_on_model_bases	4:3074935-3074936	T	3064	NM_002111.8	Transcript	frameshift_variant	256-261	111-116	37-39	QQP/HX	caacagccg/caTa	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK
```

vep-rs:

```
ins_on_model_bases	4:3074935-3074936	T	3064	NM_002111.8	Transcript	frameshift_variant	256-261	111-116	37-39	QQQ/HX	cagcagcaa/caTa	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK
```

The two rows differ in `Amino_acids` and `Codons`, before the slash, and in nothing else:

| Part of the row | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Codons`, reference half | `caacagccg` | `cagcagcaa` |
| `Amino_acids`, reference half | `QQP` | `QQQ` |
| `Codons`, alternate half | `caTa` | `caTa` |
| `Amino_acids`, alternate half | `HX` | `HX` |

Both rows print `GIVEN_REF=-` and `USED_REF=-`: the insertion's reference allele is empty. The record's row on `NM_001388492.1`, a model of the same gene without the six bases, is equal in both engines: `CDS_position` 110-111, `Protein_position` 37, `Q/HX` and `caa/caTa`.

### 1.4 Why Ensembl VEP writes that

An insertion has no reference bases, so its span on the model is read from the two bases beside its point. Through the two pairs of 1.2, `map_insert` maps 4:3074935 to cDNA 255 and 4:3074936 to cDNA 262. It then starts the first coordinate one past the first base, 255 + 1 = 256, and ends the last one before the second, 262 - 1 = 261. On an ordinary model the two bases are neighbours, so the start is one above the end and the span holds nothing.

Here they are seven positions apart and the span is real: cDNA 256-261, CDS 111-116, the model's six added bases. Both engines print that span and count the allele against it for the term. One base in place of six is a change of five, no multiple of three, hence `frameshift_variant`.

The residues follow from the span. CDS 111 is the last base of residue 37 (3 × 37 = 111), and CDS 116 lies in residue 39, which is CDS 115 to 117. The codon window therefore starts at 3 × 37 - 2 = 109 and ends at 3 × 39 = 117, which is 117 - 109 + 1 = 9 bases (lines 841-843).

For the reference codons of a model with an `_rna_edit` attribute, line 873 writes the reference allele over the span before the window is read. The line assumes that the allele has the span's length. On such a model the reference allele of a substitution or a deletion is the model's own bases, so the write puts back what the sequence already holds. An insertion keeps its empty reference allele (line 451). Here the write therefore replaces CDS 111-116 with nothing. The sequence gets six bases shorter, and the base that was CDS 117 stands at 111.

Line 877 then reads nine bases from position 109 of the shortened sequence: CDS 109 and 110, then the seven bases that were CDS 117 to 123, 2 + 7 = 9. They are `ca`, then `acagccg`, the last of them the base that was CDS 123. The six bases the row's own `CDS_position` names are not in the string, and CDS 118-123, the codons of residues 40 and 41, are. The reference amino acids are the translation of the string (lines 724-736): `caa`, `cag` and `ccg` give `QQP`.

The alternate codons are read without the write, with the window shortened by the span and lengthened by the allele: 9 + (1 - 6) = 4 bases. They are CDS 109 and 110, the inserted base and CDS 117, `caTa`, and both engines agree on them. Ensembl VEP's two halves thus disagree about what follows the variant. The alternate half ends at CDS 117, the last base of residue 39, and the reference half runs on to CDS 123.

The entry [`refseq_edited_model_deleted_span_codons`](refseq-edited-model-deleted-span-codons.md) is the same write on the opposite shape. Its reference allele is longer than the CDS span, so the window shows genomic bases the model lacks, in upper case. On this shape the allele is empty and the span is not. The window then shows bases of the model from beyond the row, in lower case, and the amino acids can differ too.

### 1.5 Why vep-rs's row is the correct one

The record inserts a base between two neighbouring genomic positions, and the model's own mapper pairs place those positions on either side of its six added bases. The row's `CDS_position` is therefore 111-116 in both engines, and its `Protein_position` the residues that span touches, 37 to 39. The reference codons of a row are the model's bases over the residues the row prints. Those are CDS 109-117, which the cache holds as `cagcagcaa`. Residues 37, 38 and 39 of the model's cached protein, `QQQ`, are translated from them. vep-rs writes both.

Ensembl VEP's string cannot be read off its own row. It leaves out the six bases `CDS_position` names and holds the codons of residues 40 and 41, which `Protein_position` 37-39 does not include. Its reference amino acids, `QQP`, are not residues 37 to 39 of the model's protein. On every row where the write changes nothing, the same routine returns the model's window.

Of 51 records written on and around 8 points of added bases inside coding sequences, 2 of them mitochondrial, 92 rows belong to records written on a nuclear point. The reference codons differ on 68 of the 92 and the reference amino acids on 66. On all 68 the model's bases over the printed residues are vep-rs's string, and cutting the span and reading the unchanged window gives Ensembl VEP's. On all 66 the cached protein's residues at the printed positions are vep-rs's reference amino acids and never Ensembl VEP's. The 68 rows are 5 distinct windows at 5 points. The alternate halves, the three position columns and the terms are equal on every one.

On the 6 rows of records written on the two mitochondrial points, the same two rules give both engines' values.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf` under the two VCF header lines, the command as vep-rs takes it:

```
vep -i record.vcf -o out.txt --offline --refseq \
  --json_cache /path/to/json_cache/homo_sapiens_refseq/116_GRCh38 \
  --fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --species homo_sapiens --assembly GRCh38
```

The cache is the conversion of Ensembl's RefSeq GRCh38 cache of release 116 ([RefSeq and merged caches](../cache-setup.md#refseq-and-merged-caches)). The committed GRCh38 RefSeq and merged caches hold no model with added bases inside its coding sequence, so the example runs on the whole converted cache. vep-rs writes `QQQ/HX` and `cagcagcaa/caTa` on the row of `NM_002111.8`. Ensembl VEP 116.3 takes `--refseq` with the cache options of [Reproducing these examples](reproducing.md), its cache root holding `homo_sapiens_refseq/116_GRCh38`. It writes `QQP/HX` and `caacagccg/caTa`.

The reproduction `scripts/divergences/reproductions/refseq_edited_model_inserted_span_codons.toml` runs this record and the record of section 2 on both engines. It copies the two region files it needs from the RefSeq cache it is given and converts them with `tools/vep-storable-to-json`. It passes `--use_given_ref`, under which neither engine asks for a FASTA. An insertion's reference allele is empty with the flag and without it, and the compared columns are the same.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) holds the class. A comparison that names it gives the comparator the cache vep-rs read, as `--json-cache`. From the model's coding sequence in that cache the comparator derives both strings of 1.5: the model's bases over the printed residues, and as many bases read with the span cut out.

It sets `Codons` aside where vep-rs's reference half is the first string, Ensembl VEP's is the second, the two differ and the alternate halves are equal. It sets `Amino_acids` aside beside such `Codons` where the alternate halves are equal and each engine's reference half is the translation of its own reference codons under the model's genetic code, as that engine reads the code. Where the model's cache record names a genetic code the comparator holds no table for, such `Amino_acids` end the comparison with an error that names the model.

The row is an insertion's on a model whose rows carry `BAM_EDIT=OK` in both engines, with an empty `USED_REF` and the same `CDS_position` and `Protein_position` in both. Everything else stays a difference and is charged to vep-rs: any other pair of strings, a term, a position, an alternate half, and every such row of a comparison that gives no cache.

Compared over the genome-wide ClinVar files under `--refseq`, the class sets aside the `Codons` of the 2 GRCh38 rows and the 3 GRCh37 rows counted above, and no other row. The golden harness compares `Codons` and `Amino_acids` on the rows of the RefSeq and merged corpora, and none of their records has the shape.

The nearest class, `refseq_edited_model_deleted_span_codons`, accepts no such row. Its test asks for a `USED_REF` longer than the CDS span and a difference in an upper-case base, and it takes the amino acids only beside codons of that shape. A row of this shape has an empty reference allele and differs in lower-case bases, and on 66 of the 68 rows above the amino acids differ as well.

Rows near the shape agree. Of the 92 rows of 1.5, 24 are equal in `Codons` and 26 in `Amino_acids`. They are the rows of models without the added bases there, the rows on chromosome 7, where cutting the model's one added base leaves the same string, and for the amino acids the 2 rows of section 2's model. A row of a record one base beside a point, or of one that only the `--hgvs` shift moves onto it, is equal in both columns, 162 of 162 on the nuclear points. The counts are the same with and without `--hgvs`.

## 2. A ClinVar insertion on the three bases `NM_015120.4` holds and the genome lacks

**In short.** `2 73385942 1167157 A AGG`, a ClinVar record, inserts two bases between 2:73385942 and 2:73385943. `NM_015120.4` (ALMS1, forward strand) holds three bases there that the genome lacks, its coding bases 75 to 77. Under `--refseq` both engines write `frameshift_variant`, `CDS_position` 75-77 and `Protein_position` 25-26. vep-rs writes the reference codons `gaggaa`, the model's coding bases 73 to 78. Ensembl VEP 116.3 writes `gaagag`, coding bases 73, 74 and 78 to 81. The two strings translate alike, so `Amino_acids` reads `EE/EX` in both engines and `Codons` alone differs.

### 2.2 The record

```
2	73385942	1167157	A	AGG	.	.	.
```

Both engines drop the shared first base: two bases, `GG`, are inserted between 2:73385942 and 2:73385943. The model's `_rna_edit` attribute `186 185` holds three bases, its cDNA 186 to 188. The coding sequence starts at cDNA 112, so they are CDS 186 - 111 = 75 to 188 - 111 = 77. CDS 75 is the last base of residue 25 (3 × 25 = 75), and CDS 77 lies in residue 26, which is CDS 76 to 78. The codon window is CDS 3 × 25 - 2 = 73 to 3 × 26 = 78, six bases.

### 2.3 What each engine writes

The row of `NM_015120.4` in the default format. Ensembl VEP 116.3:

```
1167157	2:73385942-73385943	GG	7840	NM_015120.4	Transcript	frameshift_variant	186-188	75-77	25-26	EE/EX	gaagag/gaGGa	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK
```

vep-rs:

```
1167157	2:73385942-73385943	GG	7840	NM_015120.4	Transcript	frameshift_variant	186-188	75-77	25-26	EE/EX	gaggaa/gaGGa	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK
```

The record's row on `NM_001378454.1`, a model of the same gene without the three bases, is equal in both engines: `CDS_position` 74-75, `Protein_position` 25, `E/EX` and `gaa/gaGGa`.

### 2.4 Why Ensembl VEP writes that

Line 873 replaces CDS 75-77 with the insertion's empty reference allele, and line 877 reads six bases from CDS 73 of the shortened sequence. They are CDS 73 and 74, then the four bases that were CDS 78 to 81, 2 + 4 = 6: `ga`, then `agag`. The alternate codons are 6 + (2 - 3) = 5 bases, CDS 73 and 74, the two inserted bases and CDS 78, `gaGGa` in both engines.

### 2.5 Why vep-rs's row is the correct one

The row prints `Protein_position` 25-26, and the model's bases over those residues are CDS 73 to 78, `gaggaa`. Ensembl VEP's string leaves out CDS 75-77, the span the row prints, and holds CDS 79-81, the codon of residue 27. Both strings happen to translate to `EE`, the model's residues 25 and 26, so the difference shows in `Codons` alone.

### 2.6 Try it yourself

The command of 1.6 on this record writes `gaggaa/gaGGa` on the row of `NM_015120.4` from vep-rs, and Ensembl VEP 116.3 writes `gaagag/gaGGa`. The reproduction of 1.6 holds both rows, and the two rows of the models without the added bases, to their recorded values.
