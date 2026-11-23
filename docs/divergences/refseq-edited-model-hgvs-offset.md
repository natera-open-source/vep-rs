# The HGVS_OFFSET of an insertion or deletion on a RefSeq model whose alignment edited bases

**The defect in one sentence.** On a RefSeq model whose alignment edited bases, Ensembl VEP gives an insertion or deletion the 3' shift of whichever transcript it shifted first. A transcript on the other strand therefore carries that transcript's shift and, where it is zero, leaves `HGVS_OFFSET` out. vep-rs computes every row's shift on its own transcript's strand and prints the distance in `HGVS_OFFSET` whenever the notation was moved.

**Why this is a defect.** `HGVS_OFFSET` states how many bases the HGVS notations of a row were moved 3' along the transcript. The move is a property of the allele and of its own transcript: the bases following the insertion point on that transcript's strand. A transcript on the opposite strand reads the genome in the other direction, so its shift says nothing about this one.

A duplication's notation survives, because Ensembl VEP's dup test looks 3' from the unshifted position first. The row therefore carries the shifted notation beside the inherited offset, none here. A deletion's or a plain insertion's notation would move with the inherited shift and differ as well.

**What vep-rs does instead.** vep-rs computes the shift of every row on its own transcript's strand and prints the distance in `HGVS_OFFSET` whenever the notation was moved, on RefSeq and merged caches alike.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 150-157 and 164-199, `_return_3prime`: on an `NM_` or `XM_` model carrying an `_rna_edit` attribute, the method walks every shift hash already attached to the variant, those of other transcripts and other alleles included, and copies the first whose flanking sequences it cannot read.
- The same file, lines 175-176, 179 and 182-183: for an insertion the method sets `start` to the base before it and `end` to the base after it, so the sequence it cuts is two bases longer than the length the test at line 182 expects. For a deletion the test adds the length of the `REF/ALT` string, not of the deleted bases. The test fails for every insertion and deletion, the first hash is copied at line 186, and the method returns at line 187 before the transcript-level shift at lines 211-256 runs.
- The same file, lines 1446-1447, `hgvs_transcript`: the copied hash's `_hgvs_offset` becomes the row's offset; `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1761-1764: an offset of zero prints nothing.
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 570-589, `hgvs_variant_notation`: an insertion is tested against the bases 3' of the unshifted position first, so the duplication is named even though no shift was applied.

**What correcting it would change.** Every insertion and deletion on an edited model would carry its own transcript's shift, as vep-rs's rows do and as Ensembl VEP's rows do when the record overlaps no transcript on the other strand.

**The class it produces.** `refseq_edited_model_hgvs_offset`, applied by the format comparator to the `HGVS_OFFSET` of a coding row carrying `BAM_EDIT=OK` whose `HGVSc` names a duplication and agrees, when Ensembl VEP prints no offset and vep-rs prints one. No golden corpus carries a row of it: the GRCh37 RefSeq corpus is compared there without `--hgvs`, and the GRCh38 RefSeq corpus's `hgvs` flag set has none. The concordance figures score consequence sets and never see the field, so nothing is set aside in them.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A 21-base duplication in KRTAP10-7 beside two TSPEAR models on the other strand

**In short.** A 21-base insertion in the single exon of KRTAP10-7 duplicates the 21 bases that follow it, so both engines write `NM_198689.3:c.986_1006dup`. vep-rs adds `HGVS_OFFSET=21`, the distance the duplication moved 3'; Ensembl VEP adds nothing, because it reused the zero shift of the two TSPEAR models on the opposite strand, which it had computed first. Run on the KRTAP10-7 model alone, Ensembl VEP writes `HGVS_OFFSET=21` as well.

### 1.1 What a reader needs to know first

A converted RefSeq or merged cache carries the attributes Ensembl VEP writes about each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). A model the alignment corrected prints `BAM_EDIT=OK` on its rows, and its `_rna_edit` attributes name each edit as `<start> <end> <bases>`. Under `--hgvs` both engines move an insertion or deletion to its most 3' position along the transcript before writing its notation, and `HGVS_OFFSET` prints the distance moved. An insertion that moves is written as a duplication of the bases it passed over. A row whose notation did not move carries no offset.

### 1.2 The record

The record is `21 46021521 rs782070127 T TGCGCCCCCACCTCCTCCTGCC`, record 52 of `variants.vcf` in the `GRCh37-refseq` golden corpus (`tests/golden/116/GRCh37-refseq/`), a gnomAD v2.1.1 insertion of 21 bases after `21:46021521`:

```
21	46021521	rs782070127	T	TGCGCCCCCACCTCCTCCTGCC	.	.	.
```

It lies in the single exon of `NM_198689.3` (KRTAP10-7, forward strand), a model whose alignment carries 17 edits, among them the 15-base deletion `160 174` (`REFSEQ_OFFSET=-15`), and in an intron of `NM_001272037.2` and `NM_144991.3` (TSPEAR, reverse strand). Both engines run `--refseq --hgvs`.

### 1.3 What each engine writes

The KRTAP10-7 row, default format. Ensembl VEP:

```
rs782070127	21:46021521-46021522	GCGCCCCCACCTCCTCCTGCC	386675	NM_198689.3	Transcript	inframe_insertion	1010-1011	985-986	329	C/CAPTSSCR	tgc/tGCGCCCCCACCTCCTCCTGCCgc	-	IMPACT=MODERATE;STRAND=1;REFSEQ_OFFSET=-15;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK;HGVSc=NM_198689.3:c.986_1006dup;HGVSp=NP_941962.1:p.Cys335_Gln336insArgAlaProThrSerSerCys
```

vep-rs:

```
rs782070127	21:46021521-46021522	GCGCCCCCACCTCCTCCTGCC	386675	NM_198689.3	Transcript	inframe_insertion	1010-1011	985-986	329	C/CAPTSSCR	tgc/tGCGCCCCCACCTCCTCCTGCCgc	-	IMPACT=MODERATE;STRAND=1;REFSEQ_OFFSET=-15;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK;HGVSc=NM_198689.3:c.986_1006dup;HGVSp=NP_941962.1:p.Cys335_Gln336insArgAlaProThrSerSerCys;HGVS_OFFSET=21
```

Every column agrees, the HGVSc and HGVSp included; vep-rs's row ends with `HGVS_OFFSET=21`. The two TSPEAR rows are identical on both engines: `c.-122-33653_-122-33633dup` and `c.83-33653_83-33633dup`, with no offset.

### 1.4 Why Ensembl VEP writes that

`hgvs_transcript` calls `_return_3prime` to build a shift hash for the allele. The method first shifts the allele along the genome in the transcript's direction, 21 bases here, and attaches that hash to the variant. On an Ensembl model or an unedited RefSeq model it stops there and the 21 is the row's offset. On an edited `NM_` model it goes on to compare the allele's flanking sequences with every shift hash already attached to the variant, whatever transcript or allele made it, so that a matching hash can be reused.

For an insertion the method sets `start` to the base before the insertion point and `end` to the base after it, and cuts a sequence two bases longer than the length its sanity test expects. For a deletion the test adds the length of the `REF/ALT` string. The test fails on every insertion and deletion, so the first hash in the list is copied and the method returns before the transcript-level shift is computed.

Here the first hash is the one the two TSPEAR models made. On the reverse strand the 3' direction runs leftwards along the genome; the base before the insertion point is `T` and the inserted bases end in `C`, so that shift is zero. The KRTAP10-7 row inherits a zero offset, which Ensembl VEP does not print.

The notation survives because `hgvs_variant_notation` tests an insertion against the bases 3' of the unshifted position first. The 21 bases after c.985 equal the insertion, so the duplication `c.986_1006dup` is named from the unshifted position, the same text a 21-base shift would give.

Source (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 140-145: the genomic shift; lines 150-157: the edited-model branch entered; line 164: every shift hash of the variant; lines 175-176 and 179: the two extra bases; lines 182-188: the failing test and the copy; lines 211-256: the transcript-level shift that is never reached.
- The same file, lines 1446-1447 and 1835-1838: the copied hash's offset becomes the row's; `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1761-1764: a zero offset prints nothing.
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 570-589: the 3'-first duplication test.

### 1.5 Why vep-rs's row is the correct one

The reference bases `21:46021522-46021542` are `GCGCCCCCACCTCCTCCTGCC`, the inserted sequence, and the next base `21:46021543` is `A`. The edited model reads the same bases at cDNA 1011-1031 (CDS 986-1006), followed by `A` at cDNA 1032. An insertion of those 21 bases after c.985 can therefore be moved 3' by exactly 21 bases, no more, and its most 3' description is the duplication of c.986-1006. `HGVS_OFFSET` states that move, 21, on both the genome and the edited transcript. Zero would say the notation was not moved, which `c.986_1006dup` contradicts.

Ensembl VEP gives the same answer when nothing interferes: on the same record with `--transcript_filter "stable_id match NM_198689"`, so that no TSPEAR model is annotated first, its row ends with `HGVS_OFFSET=21`, byte for byte vep-rs's row. On an unedited model it prints the offset for the same shape, as `rs1242762494` in the `GRCh37-hgvs` corpus shows (`c.43_47dup;HGVS_OFFSET=5`).

### 1.6 Try it yourself

```
vep -i tests/golden/116/GRCh37-refseq/variants.vcf -o out.txt --offline --cache \
  --dir_cache tests/golden/116/GRCh37-refseq/cache_root --cache_version 116 \
  --species homo_sapiens --assembly GRCh37 --refseq --hgvs \
  --fasta tests/golden/116/GRCh37-refseq/reference.fa.gz --force_overwrite --no_stats --quiet
grep rs782070127 out.txt | grep NM_198689 | grep -o 'HGVS[^=]*=[^;]*'
```

| Run | Ensembl VEP | vep-rs |
| --- | --- | --- |
| record as given | `c.986_1006dup`, no `HGVS_OFFSET` | `c.986_1006dup;HGVS_OFFSET=21` |
| `--transcript_filter "stable_id match NM_198689"` | `c.986_1006dup;HGVS_OFFSET=21` | `c.986_1006dup;HGVS_OFFSET=21` |

With Ensembl VEP 116.3 the first run needs its RefSeq GRCh37 cache and the GRCh37 reference; the second run adds the filter and prints the offset.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) applies `refseq_edited_model_hgvs_offset` when a comparison case names it in `masks`. The shape is a differing `HGVS_OFFSET` (the JSON `hgvs_offset`) on a row carrying `BAM_EDIT=OK`, a CDS position and an `HGVSc` that names a duplication. The difference is accepted when Ensembl VEP's value is blank and vep-rs's is a non-zero integer, negative on a reverse-strand transcript. Every other field of the row, the `HGVSc` included, must agree or be accepted by its own class.

The shape is that narrow for two reasons. The edited-model branch is entered only for an allele with CDS coordinates, and only a duplication keeps its notation whatever the shift. A substitution is never shifted, and a deletion or plain insertion that inherits the foreign shift changes its `HGVSc` too, a difference of another shape. Each accepted row is counted once. The GRCh37 RefSeq comparison case under `--hgvs` expects 1 row in each of the four formats and the GRCh38 case 0.

No golden corpus documents a row of the class. The GRCh37 RefSeq corpus's committed comparison runs without `--hgvs`, and the GRCh38 RefSeq corpus's `hgvs` flag set has none. Its three insertion and deletion records on edited models carry the same `HGVS_OFFSET` on both engines, 6, 7 and 5 on each of two models. The class is distinct from `refseq_edited_model_hgvsc_offset`: that class accepts a renumbered `HGVSc` whose offsets agree, this one an absent offset beside an agreeing `HGVSc`. No row of either corpus carries both.
