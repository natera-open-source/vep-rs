# The term of an insertion equal to the bases an edited RefSeq model holds and the genome lacks

**The defect in one sentence.** Where a RefSeq model edited to the genome holds bases the genome lacks, inside its coding sequence, an insertion of exactly those bases written at that point leaves the model's sequence as it is. Ensembl VEP cuts the bases out of the model before it reads the reference codons, and where the residues it then reads differ from the alternate's it writes its default `intergenic_variant` on a row inside a coding exon. vep-rs writes `synonymous_variant` with `IMPACT=LOW` on such a row, beside the positions Ensembl VEP prints.

**Why this is a defect.** The row names a transcript, a gene and a strand and prints a cDNA, a CDS and a protein position, and `intergenic_variant` is defined as a variant between genes. Both engines map the insertion onto the model's added bases and build the same alternate coding sequence, with the allele in their place. The allele is those bases, so the alternate sequence is the model's own and no residue changes, which is what Ensembl VEP's table of terms says of `synonymous_variant`. Ensembl VEP writes that term itself where the bases after the span happen to translate alike, so its term turns on bases the record does not touch.

**What vep-rs does instead.** vep-rs maps the insertion onto the same span and builds the same alternate coding sequence. It reads the reference codons from the model's coding sequence as the cache holds it, over the residues the row prints. The two windows are equal, and vep-rs writes `synonymous_variant` with `IMPACT=LOW` beside the positions Ensembl VEP prints.

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `ensembl` at release/116, `ensembl-vep` at release/116.2, the modules the Ensembl VEP 116.2 image holds):

- `Bio/EnsEMBL/Mapper.pm`, lines 485-580, `map_insert`: an insertion is mapped through the two bases beside its point (lines 489-491). Where they lie in two mapper pairs (line 504), the first mapped coordinate is moved to start one past the first base (lines 524-528) and the last to end one before the second (lines 538-542).
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, line 513, `cds_coords`, and lines 263-264, `cds_start`: the row's CDS span runs from the first coordinate's start to the last coordinate's end. Where the two bases are a model's added bases apart, the span is those bases.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, `codon`, lines 821-822 and 851: an insertion's reference allele, `-`, becomes the empty string `$seq`, and `$vf_nt_len` is the length of the CDS span, here the count of the model's added bases.
- The same routine, lines 872-874: on a model carrying an `_rna_edit` attribute (line 827), the reference allele is written over the CDS span, `substr($cds, $cds_start - 1, $vf_nt_len) = $seq`. For an insertion the model's own bases are taken out and nothing is put in.
- The same routine, lines 841-843 and 877: the reference codons are `$codon_len` bases from the first base of the first printed residue, read from the shortened sequence. Line 868 and the second arm of line 877 read the alternate codons from a sequence the write does not touch.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, line 809, `_get_ref_pep`: the reference peptide is the translation of that window. Line 1082, `synonymous_variant`: the term needs two equal peptides.
- `Bio/EnsEMBL/Variation/Utils/Constants.pm`, lines 886-890: the `include` rule of `missense_variant` asks for `increase_length => 0`, so the predicate is never run on an insertion. Lines 1051-1053: the rule of `synonymous_variant` asks for `coding => 1` alone, so that predicate is run. Line 1046 holds its description: "A sequence variant where there is no resulting change to the encoded amino acid".
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, three more predicates that are run: `frameshift` (lines 1556-1564) takes the allele's length less the CDS span's, modulo three. `inframe_insertion` (line 1111) needs alternate codons longer than the reference codons. `protein_altering_variant` (line 384) refuses two peptides of one length.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 321-325, `_skip_oc`: the `include` rules decide which predicates are run. Line 285: with no predicate true, `$DEFAULT_OVERLAP_CONSEQUENCE` is assigned. `Bio/EnsEMBL/Variation/Utils/Constants.pm`, lines 636-649: the default is `intergenic_variant`, with `IMPACT` `MODIFIER`.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 831-836: an alternate allele equal to the bases of the model's first `_rna_edit` attribute sets `invalid_alleles`. `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1713 and 1717-1718: the comment reads "the ref and alt alleles are the same", and the flag's one effect is a warning, written unless `--use_given_ref` is passed.
- The write of lines 872-874 is also behind [`refseq_edited_model_deleted_span_codons`](refseq-edited-model-deleted-span-codons.md), where a deletion's reference allele is longer than its CDS span.

**What correcting it would change.** A write made only when the reference allele has the span's length would leave the model's sequence whole on these rows. The reference window would be the model's own, equal to the alternate's, and `synonymous_variant` would accept. Each such row would carry that term and `IMPACT=LOW`, as vep-rs writes them, and its `Amino_acids` would print the unchanged residues once. No position column changes: both engines print the same three.

**The class it produces.** `refseq_edited_model_own_bases_insertion_term`: a `Transcript` row of an insertion whose `Consequence` is `intergenic_variant` in Ensembl VEP and `synonymous_variant` in vep-rs, with `IMPACT` `MODIFIER` against `LOW` and the same cDNA, CDS and protein positions in both. No comparison tool sets a row of the shape aside: a disagreement of this shape is reported as a difference and charged to vep-rs. The concordance figures are measured on the Ensembl caches, whose models hold no such bases, so they hold no row of the shape.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an insertion written exactly where a RefSeq model holds bases the genome lacks, inside its coding sequence, whose inserted bases are the model's own bases there. The two Ensembl caches the concordance datasets run on hold 0 `_rna_edit` attributes and 0 such places, in 727,139 transcript entries on GRCh38 and 253,803 on GRCh37, an entry being a transcript counted once for each region file of the cache that stores it. The concordance figures are measured without `--refseq` or `--merged`, so they hold 0 rows of the shape.

A set of 78 hand-written records, run with `--refseq --hgvs` on the RefSeq GRCh38 cache of release 116, holds a record of the shape at five places inside coding sequences. They give 21 rows of the shape: 15 at `1:16567177`, 3 at `7:94664512` and 1 each at `2:73385942`, `4:3074935` and `19:14090084`. Of the 21, 16 are on reverse-strand models and 5 on forward-strand ones. On all 21 both engines print the span of the model's own bases as `CDS_position`, the inserted bases equal the cached coding sequence over it, and no residue changes. The 15 models at `1:16567177` hold their own bases at `1:16567174` too, three bases apart, and this set's records stand at `1:16567177`.

vep-rs writes `synonymous_variant` on all 21 rows. Ensembl VEP writes `intergenic_variant` on 17 and `synonymous_variant` on 4, the row at `2:73385942` and the three at `7:94664512`. The 17 are the rows where the residues of its reference window differ from the alternate's, and the 4 the rows where they are equal.

Among the 4,398,127 records of the genome-wide ClinVar GRCh38 file, 3 are insertions written on one of the 613 places where a model of that RefSeq cache holds bases the genome lacks. Two lie on a place inside a coding sequence, and neither inserts the model's bases. The third inserts the model's own base on a non-coding model, where both engines write `non_coding_transcript_exon_variant`. The file therefore holds 0 records of the shape.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. The six bases `NM_002111.8` holds and the genome lacks, inserted where the model holds them

**In short.** `4 3074935 own_bases_4_3074935 A AGCAGCA` on GRCh38 inserts six bases between 4:3074935 and 4:3074936. The forward-strand RefSeq model `NM_002111.8` holds six bases at that point that the genome lacks, its coding positions 111 to 116, and the record's six bases are those six. Both engines print `CDS_position` 111-116 and `Protein_position` 37-39 on the model's row, with the same alternate codons. Ensembl VEP 116.2 writes `intergenic_variant` and `IMPACT=MODIFIER` there; vep-rs writes `synonymous_variant` and `IMPACT=LOW`.

### 1.1 What a reader needs to know first

A RefSeq cache carries each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). A model the alignment corrected prints `BAM_EDIT=OK` on its rows, and an `_rna_edit` attribute names each difference as `<start> <end> <bases>` in cDNA coordinates. An insertion edit has a start one past its end: `256 255` with six bases says the model holds six bases after cDNA position 255 that the genome lacks. This page calls them the model's added bases, or its own bases.

A model's mapper pairs its cDNA positions with genomic positions, one pair per unbroken stretch. Added bases belong to no pair: the pair before them ends at the genomic base before the point, and the pair after them starts at the next genomic base, several cDNA positions later. A cDNA position counts from the model's first base and a CDS position from the first base of its coding sequence. Residue n is translated from CDS positions 3n - 2 to 3n.

An insertion is written between two neighbouring genomic bases. Its span on a transcript runs from one past the base before it to one before the base after it, which on an ordinary model is empty. Where the model holds added bases between the two, the span is those bases. The row prints them as its `cDNA_position` and `CDS_position`, and the alternate coding sequence both engines build carries the allele in their place.

The `Codons` column prints the codons of the residues in `Protein_position`, the reference before the slash and the alternate after it, with the changed bases in upper case. `Amino_acids` prints their translations, once when the two are equal. This page calls the part before the slash the reference half. Ensembl VEP decides a row's term by running a list of tests, its predicates, each of which accepts or refuses the row. When none accepts, it writes a default term, `intergenic_variant`.

### 1.2 The record

Six bases inserted between two neighbouring genomic positions, written for this example, against the RefSeq GRCh38 cache of release 116:

```
4	3074935	own_bases_4_3074935	A	AGCAGCA	.	.	.
```

Both engines drop the shared first base, so the record is six bases, `GCAGCA`, inserted between 4:3074935 and 4:3074936. The cache holds, for the forward-strand model `NM_002111.8`, the `_rna_edit` attribute `256 255 GCAGCA`, a coding sequence that starts at cDNA 146, and two mapper pairs that meet at the insertion point:

| cDNA | Genome |
| --- | --- |
| 1-255 | 4:3074681-3074935 |
| 262-414 | 4:3074936-3075088 |

The six bases at cDNA 256-261 lie between the two pairs, and the genome has no base for them. A CDS position is the cDNA position less 145, so they are CDS 256 - 145 = 111 to 261 - 145 = 116. The cached coding sequence holds `gcagca` there, the record's six bases one for one. The model carries a second attribute, `13476 13475`, far from the record.

Both engines run `--refseq --use_given_ref`. The second flag lets a run read an edited cache without a reference FASTA. The 21 rows counted above were run with a FASTA and `--hgvs` in its place and carry the same terms.

### 1.3 What each engine writes

The model's row in the default format. Ensembl VEP 116.2:

```
own_bases_4_3074935	4:3074935-3074936	GCAGCA	3064	NM_002111.8	Transcript	intergenic_variant	256-261	111-116	37-39	QQP/QQQ	caacagccg/caGCAGCAa	-	IMPACT=MODIFIER;STRAND=1;BAM_EDIT=OK
```

vep-rs:

```
own_bases_4_3074935	4:3074935-3074936	GCAGCA	3064	NM_002111.8	Transcript	synonymous_variant	256-261	111-116	37-39	QQQ	cagcagcaa/caGCAGCAa	-	IMPACT=LOW;STRAND=1;BAM_EDIT=OK
```

The rows agree on the three position columns and on the alternate half of `Codons`. They differ in four places:

| Part of the row | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- |
| `Consequence` | `intergenic_variant` | `synonymous_variant` |
| `IMPACT` | `MODIFIER` | `LOW` |
| `Amino_acids` | `QQP/QQQ` | `QQQ` |
| `Codons`, reference half | `caacagccg` | `cagcagcaa` |

This entry is the `Consequence` and the `IMPACT` that follows it. The reference half of `Codons` and the `Amino_acids` read from it are a difference of their own, in the bases Ensembl VEP reads as the reference codons, and another entry's subject. 1.4 shows how the term follows from them.

The record's row on `NM_001388492.1`, a model of the same gene without the six bases, is equal in both engines: `inframe_insertion`, `CDS_position` 110-111, `Protein_position` 37 and `Q/QQQ`. On that model the six bases are new, and two residues are added.

A second record inserts, between 2:73385942 and 2:73385943, the three bases the model `NM_015120.4` holds there and the genome lacks: `2 73385942 own_bases_2_73385942 A AGGA`. Both engines write `synonymous_variant` and `IMPACT=LOW` on that model's row, with `CDS_position` 75-77, `Protein_position` 25-26 and `EE` in `Amino_acids`. It is the row at `2:73385942` among the 4 counted above.

### 1.4 Why Ensembl VEP writes that

An insertion has no reference bases, so its span on the model is read from the two bases beside its point. Through the two pairs of 1.2, `map_insert` maps 4:3074935 to cDNA 255 and 4:3074936 to cDNA 262. It starts the first coordinate one past the first base, 255 + 1 = 256, and ends the last one before the second, 262 - 1 = 261. The span is cDNA 256-261, CDS 111-116, the model's six added bases, and both engines print it.

The residues follow from the span. CDS 111 is the last base of residue 37 (3 × 37 = 111), and CDS 116 lies in residue 39, which is CDS 115 to 117. The codon window therefore starts at 3 × 37 - 2 = 109 and ends at 3 × 39 = 117, which is 117 - 109 + 1 = 9 bases (lines 841-843).

The alternate codons are read from the alternate coding sequence, the model's with the allele in place of the span. They are CDS 109 and 110, the six inserted bases and CDS 117, 2 + 6 + 1 = 9 bases: `caGCAGCAa`, three glutamines, `QQQ`. Both engines write that string.

For the reference codons of a model with an `_rna_edit` attribute, line 873 writes the reference allele over the span before the window is read. An insertion's reference allele is empty, so the write replaces CDS 111-116 with nothing. The sequence gets six bases shorter, and the base that was CDS 117 stands at 111. Line 877 then reads nine bases from position 109: CDS 109 and 110, then the seven bases that were CDS 117 to 123, 2 + 7 = 9. They are `caacagccg`, and their translation is `QQP`.

The two peptides, `QQP` and `QQQ`, have one length and differ in their third residue, residue 39. `synonymous_variant` needs equal peptides and refuses (line 1082). `missense_variant`, the term for changed residues of unchanged length, is not run. An insertion's genomic span holds no base, so Ensembl VEP counts the allele as one that increases the length, and the term's `include` rule leaves such alleles out.

The other predicates refuse as well. `frameshift` subtracts the CDS span from the allele's length: 6 - (116 - 111 + 1) = 0, a multiple of three. `inframe_insertion` needs alternate codons longer than the reference codons, and both strings hold nine bases. `protein_altering_variant` returns false for two peptides of one length.

On each of the 17 rows counted above, 21 of the 33 predicates a transcript row of a sequence variant can meet are left out by their `include` rules, the other 12 are run and none accepts. Line 285 then assigns the default, `intergenic_variant`.

On `NM_015120.4` the same write decides the other way. The row's window is CDS 73 to 78, and with CDS 75-77 cut out Ensembl VEP reads CDS 73, 74 and 78 to 81. That string, `gaagag`, is not the model's `gaggaa`, and it translates to the same two residues, `EE`. The peptides are equal, `synonymous_variant` accepts, and Ensembl VEP writes the term vep-rs writes. Which of the two terms a row gets is decided by the bases that follow the span.

Ensembl VEP's source has a test for this allele. Lines 831-836 compare each alternate allele with the bases of the model's first `_rna_edit` attribute and set a flag when the two are equal, as they are on 3 of the 21 rows, this one among them. The comment where the flag is read says the reference and alternate alleles are then the same. The flag's one effect is a warning, `Transcript-assembly mismatch in` and the record's name, which `--use_given_ref` silences. The row is built as before.

### 1.5 Why vep-rs's row is the correct one

The record's flanking bases are 4:3074935 and 4:3074936. The first mapper pair ends at 3074935, cDNA 255, and the second starts at 3074936, cDNA 262. Counted from the coding start at cDNA 146, those are coding positions 255 - 145 = 110 and 262 - 145 = 117. Between them lie cDNA 256 to 261, coding positions 111 to 116: the six bases the genome has no position for, and the row's `CDS_position` in both engines.

Coding positions 109 to 111 are the codon of residue 37, since 3 × 37 = 111, 112 to 114 the codon of residue 38 and 115 to 117 that of residue 39: the row's `Protein_position` 37-39. The model's nine bases there are `cagcagcaa`, the reference codons vep-rs prints, three glutamines. The alternate window keeps positions 109 and 110, holds the six inserted bases in place of 111 to 116, and keeps 117: `caGCAGCAa`. The six inserted bases equal the model's own 111 to 116 one for one, so the two windows are the same nine bases.

Six bases stand on the same six: 9 - 6 + 6 = 9 bases, the same length, the same frame and the same three residues. Read on the model, the changed coding sequence is the model's with those bases exchanged for themselves. That is a `synonymous_variant` by the description Ensembl VEP's table gives the term, no resulting change to the encoded amino acid. `intergenic_variant` says the variant lies between genes, and the same row places it on coding positions 111 to 116 of the model.

Ensembl VEP's reference string cannot be read off its own row. It leaves out CDS 111-116, the span `CDS_position` names, and holds CDS 118-123, the codons of residues 40 and 41, which `Protein_position` 37-39 does not include. Its `QQP` is not residues 37 to 39 of the model's protein, which are `QQQ`. On all 21 rows counted above the inserted bases equal the cached coding sequence over the span and no residue of the model changes, and Ensembl VEP's two terms follow its own reference window on all 21.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache <cache>/homo_sapiens_refseq/116_GRCh38 \
  --species homo_sapiens --assembly GRCh38 \
  --refseq --use_given_ref
grep NM_002111.8 out.txt
```

`<cache>` is the conversion of Ensembl VEP's RefSeq GRCh38 cache of release 116 ([Cache conversion](../cache-conversion.md)). vep-rs writes `synonymous_variant` on the model's row. Ensembl VEP 116.2 takes the same flags with the cache options of [Reproducing these examples](reproducing.md) and its RefSeq cache, and writes `intergenic_variant`. With the second record of 1.3 and `grep NM_015120.4`, both write `synonymous_variant`.

No committed corpus holds a model with such bases inside its coding sequence. The entry's reproduction, `scripts/divergences/reproductions/refseq_edited_model_own_bases_insertion_term.toml`, therefore copies the two region files of Ensembl VEP's cache that hold the two models, `4/3000001-4000000` and `2/73000001-74000000`. It converts them in the Ensembl VEP image and runs both records on both engines. It holds four rows to their recorded values: the two models' rows and the rows of the two models of the same genes without the bases.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`scripts/concordance/compare_vep_formats.py`) has no class of this shape, and no committed corpus holds a row of it. A disagreement of the shape is reported as a difference in `Consequence` and `IMPACT` and charged to vep-rs.

Other columns of such a row differ as well and stay differences: `Amino_acids` and the reference half of `Codons`, which come from the same write and are another entry's subject, as 1.3 says, and, under `--hgvs`, the `HGVSp`. On the row of 1.3 Ensembl VEP's `HGVSp` names a change of residue 39, and vep-rs's says that the three residues from 37 are unchanged.

Rows beside the shape are no part of this entry. Where Ensembl VEP's reference window translates to the alternate's residues, as on `NM_015120.4`, both engines write `synonymous_variant` and no term differs. An insertion of as many other bases at such a point also gets Ensembl VEP's default term, and is another entry's subject.
