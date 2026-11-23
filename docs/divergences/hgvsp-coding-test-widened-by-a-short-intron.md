# The HGVSp of an insertion before an intron's end on a transcript with a very short intron

**The defect in one sentence.** On a transcript that has an intron of 13 bases or fewer, Ensembl VEP 116.3 widens every exon by 12 bases in the test that decides whether an allele is coding. Under `--hgvs`, an insertion written up to 12 bases before the end of another, ordinary intron inside the coding region, and carried to that end by the 3' shift, then gets a definite `HGVSp` beside columns that place it in the intron. vep-rs writes no `HGVSp` for such an insertion.

**Why this is a defect.** The comment above the widening gives its purpose: a variant that falls in such an intron is to be called exonic. The test is asked of every exon of the transcript, so it also answers yes beside an exon thousands of bases from the short intron. The `HGVSp` that follows states one protein consequence, the inserted bases taken as the exon's first bases. The row's `HGVSc` names intron positions alone, and its CDS and protein positions are empty. The HGVS Nomenclature gives `p.?`, an unknown effect, where a change most probably affects splicing. On a transcript with the same junction and no short intron, Ensembl VEP writes no `HGVSp` for the same record.

**What vep-rs does instead.** For an insertion the record writes inside an intron and the 3' shift carries to the intron's end, vep-rs asks the test of the allele as written, without the widening. Both bases of these insertions are intron bases, so vep-rs writes no `HGVSp`, on a transcript with a short intron as on one without. Its `HGVSc`, its position columns and its consequence terms on these rows equal Ensembl VEP's.

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `ensembl` at release/116, the modules the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 1032-1035, `_create_intron_trees`: an intron whose end lies at most 12 above its start, 13 bases or fewer, sets `_has_frameshift_intron` on the transcript. Lines 912-915 set it where no interval tree is installed.
- The same file, lines 858-866, `_overlapped_exons`: the comment of lines 858-860 says such introns "can be up to 12 bases long and if a variant falls in one we actually want to call it exonic". Line 861 sets a stretch of 12 when the flag is set, and line 865 asks the exon tree for the allele's span widened by it on both sides.
- The same file, line 1062, `_exon_interval_tree`: the tree holds every exon of the transcript, so the widened span is tested against all of them and not against the exons beside the short intron alone.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 497-516, `_bvfo_preds`: the introns are looked up first, which builds the intron trees and sets the flag (`BaseTranscriptVariation.pm`, lines 813 and 970). The exons follow, and a hit sets `exon`.
- The same routine, lines 535-537: without an exon the allele is `non_coding` and the test ends. Lines 540 and 555-557: with one, and where lines 523-525 have not flagged the allele `utr` for lying wholly outside the coding region, the CDS coordinates of the allele as written are fetched, and a list that is neither one coordinate nor empty sets `coding`. An insertion inside an intron maps to one gap.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, `hgvs_protein`: line 1711 shifts the allele 3', and lines 1730-1738 write nothing unless `coding` is set and the shifted position has a first and a last residue. The flags are the ones computed for the allele as written and kept (`BaseVariationFeatureOverlapAllele.pm`, lines 340 and 363).
- `Bio/EnsEMBL/Mapper.pm`, lines 485-580, `map_insert`: an insertion is mapped through the two bases beside it. Where one is an intron base and the other an exon base, lines 522-528 and 536-542 place the insertion against the exon base. The shifted insertion therefore stands directly before the exon's first base.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2409-2421, `_get_alternate_cds`: the inserted bases are joined into the coding sequence at that point. Lines 2053-2054 and 2012-2014 then write a frameshift, and lines 2063-2065, 2464-2471 and 1947 a duplication.
- The same stretch, set on a transcript by an earlier record, is behind [`nearest_exon_jb_stretch_priming`](../plugins/transcript-model.md#3-nearestexonjb-a-value-on-an-upstream-or-downstream-row-that-depends-on-the-records-before-it).

**What correcting it would change.** A widening kept to the exons beside a short intron, the case its comment names, would leave `exon` unset on these rows. `coding` would stay unset, `hgvs_protein` would return nothing, and the `HGVSp` would be empty, as vep-rs writes it and as Ensembl VEP writes it on a transcript without such an intron. No other column of these rows differs between the engines. The same flag chooses which term tests are run, so a correction there would also reach the consequence terms of such rows, on which the engines agree.

**The class it produces.** `hgvsp_coding_test_widened_by_a_short_intron`: a `Transcript` row of an insertion whose `HGVSp` Ensembl VEP fills and vep-rs leaves empty, with the same intron-position `HGVSc`, the same splice or intron terms and an empty `CDS_position` and `Protein_position` in both. No comparison tool sets a row of the shape aside: a disagreement of this shape is reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an insertion written 1 to 12 bases before the end of an intron, which the 3' shift carries to that end. The intron stands before a coding exon, on a transcript that has an intron of 13 bases or fewer elsewhere, in a run with `--hgvs`.

The rows Ensembl VEP 116.3 writes for the genome-wide ClinVar files hold 265 rows on GRCh37 and 1,212 on GRCh38 whose `HGVSc` is a duplication ending on an intron's last base before a coding base, or an insertion between the two. On every one Ensembl VEP's test answers the same with the widening as without it, so the files hold 0 rows of the shape.

Two records written for this page show it, one on each strand. They give 22 transcript rows in each engine, 12 and 10, and the engines differ on 4 of them, in the `HGVSp` alone, in each of the default, tab, VCF and JSON formats. The records stand 1 and 3 bases before the intron's end; the reach of 12 bases is the stretch of line 861.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A copy of an intron's last base, written one base before the intron's end, on ENST00000288304

**In short.** `21 38610911 intron_end_dup_21_38610911 C CC` on GRCh37 inserts one base between 21:38610911 and 21:38610912, inside an intron of 1,886 bases of the reverse-strand transcript ENST00000288304. The base copies the intron's last base and is written one base before the intron's end. Both engines write `splice_acceptor_variant` and the `HGVSc` `c.76-1dup`, with no CDS or protein position. Ensembl VEP 116.3 adds `HGVSp=ENSP00000288304.5:p.Pro26AlafsTer27`, and vep-rs writes no `HGVSp`. The transcript has introns of 2 and 4 bases more than 10,000 bases away. On ENST00000309117, which shares the junction and has no such intron, neither engine writes an `HGVSp`.

### 1.1 What a reader needs to know first

`--hgvs` adds up to two notations to a transcript row. `HGVSc` describes the change on the transcript's coding sequence, and `HGVSp` the change predicted on its protein. A coding base is counted from the first base of the coding sequence, and residue n is translated from coding bases 3n - 2 to 3n. A base of an intron inside the coding region is named from the nearest coding base: `c.76-1` is the last intron base before coding base 76, and `c.76-2` the one before it.

An insertion that repeats the bases beside it can be written at more than one place, since the changed sequence is the same wherever the copy stands in the repeat. The HGVS Nomenclature fixes one place, the last in the transcript's direction, and calls this the 3' rule. Both engines shift an insertion there before they write a notation. On a reverse-strand transcript that direction runs toward lower genomic coordinates. A copy of the bases directly before it is written as a duplication, `dup`.

An intron ends at its acceptor site. Its last two bases, an `A` and a `G` on the transcript in most introns, mark where the next exon starts. A change at those bases, or at the exon's first base, can move or remove the splice. What the protein becomes then depends on the splice the cell makes, which the record does not show.

Ensembl VEP writes an `HGVSp` only for an allele it has flagged `coding`. An insertion is flagged from the two bases beside it as the record writes them. It is `coding` when one of them lies in an exon, the `exon` flag, and inside the coding region. This page calls that the coding test. Ensembl VEP's source calls an intron of 13 bases or fewer a frameshift intron, and this page calls it a short intron and any longer one an ordinary intron.

### 1.2 The record

One base inserted between two neighbouring positions, written for this example, against the Ensembl GRCh37 cache of release 116:

```
21	38610911	intron_end_dup_21_38610911	C	CC	.	.	.
```

Both engines drop the shared first base, so the record is one base, `C`, inserted between 21:38610911 and 21:38610912. ENST00000288304.5 is a reverse-strand transcript of 10 exons. Four of its parts matter here:

| Part | Genomic span | Bases |
| --- | --- | --- |
| Intron 2 | 21:38610911-38612796 | 1,886 |
| Exon 3 | 21:38610761-38610910 | 150 |
| Intron 6 | 21:38600643-38600644 | 2 |
| Intron 8 | 21:38600094-38600097 | 4 |

On this strand the transcript runs from high to low coordinates. Intron 2 ends with the genome's `T` at 38610912 and `C` at 38610911, which the transcript carries as `A` and `G`. Exon 3 opens at 38610910, the genome's `G` and the transcript's `C`. The exon opens at cDNA 869 and the coding sequence at cDNA 794, so that base is coding base 869 - 794 + 1 = 76. It is the first base of residue 26, since 3 × 26 - 2 = 76.

The inserted base is a second `C` beside the `C` at 38610911, a `G` on the transcript: a copy of the intron's last base. It can be written before that base or after it, with one changed sequence. The exon's first base is another base, so the 3' rule stops between 38610911 and 38610910, the intron's last base and the exon's first.

The coding sequence has 762 bases, 253 codons and the stop codon: 3 × 253 + 3 = 762. The two short introns lie far from the record. The nearer one ends at 38600644, and 38610910 - 38600644 = 10,266 positions separate it from the first base of exon 3.

Both engines run `--hgvs` with a reference FASTA.

### 1.3 What each engine writes

The row on ENST00000288304 in the default format. Ensembl VEP 116.3:

```
intron_end_dup_21_38610911	21:38610911-38610912	C	ENSG00000157538	ENST00000288304	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1;HGVSc=ENST00000288304.5:c.76-1dup;HGVSp=ENSP00000288304.5:p.Pro26AlafsTer27;HGVS_OFFSET=-1
```

vep-rs:

```
intron_end_dup_21_38610911	21:38610911-38610912	C	ENSG00000157538	ENST00000288304	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1;HGVSc=ENST00000288304.5:c.76-1dup;HGVS_OFFSET=-1
```

The rows differ in one item of the last column, the `HGVSp`. Ensembl VEP's value says that residue 26, a proline, becomes an alanine, that the frame shifts there, and that the new frame ends at a stop 27 codons on. `HGVS_OFFSET=-1` in both rows is the one-base shift, which runs toward lower coordinates on this strand.

The record writes 12 rows in each engine, and the other 11 are equal in every column. Four of the 12 are rows of coding transcripts with an exon at 21:38610761-38610910:

| Transcript | `HGVSc` | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| ENST00000288304 | `c.76-1dup` | `p.Pro26AlafsTer27` | no `HGVSp` |
| ENST00000309117 | `c.202-1dup` | no `HGVSp` | no `HGVSp` |
| ENST00000476950 | `c.202-1dup` | no `HGVSp` | no `HGVSp` |
| ENST00000398998 | `c.58-1dup` | no `HGVSp` | no `HGVSp` |

All four carry `splice_acceptor_variant`, and the insertion stands one base before that exon on each. Of the four, ENST00000288304 alone has a short intron.

A second record makes the same change on the forward strand, three bases at once: `7 100333339 intron_end_dup_7_100333339 C CAAG`. It inserts `AAG` between 7:100333339 and 7:100333340, a copy of the last three bases of intron 1 of ENST00000538115, 7:100331845-100333342, written three bases before the intron's end. Exon 2 opens at 100333343, coding base 54. On that transcript's row both engines write `splice_region_variant,intron_variant` and the `HGVSc` `c.54-3_54-1dup`. ENST00000538115 has a short intron of 2 bases, intron 30, at 7:100371476-100371477.

The second record writes 10 rows in each engine, all on transcripts whose exon opens at 100333343:

| Transcripts | Short intron | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| ENST00000538115, ENST00000542585, ENST00000546292 | 2 bases | `p.Arg18dup` | no `HGVSp` |
| ENST00000546213 | 2 bases | no `HGVSp` | no `HGVSp` |
| ENST00000348028 and five others | none | no `HGVSp` | no `HGVSp` |

The five others are ENST00000349350, ENST00000421100, ENST00000427578, ENST00000443370 and ENST00000449052. On ENST00000546213 the junction lies in the 5' UTR, and its `HGVSc` is `c.-4216-3_-4216-1dup`. The rows of the six transcripts without a short intron carry a third term, `splice_polypyrimidine_tract_variant`, in both engines; 1.7 says why.

### 1.4 Why Ensembl VEP writes that

The record lies inside the transcript, so Ensembl VEP looks up the introns it overlaps, and that builds the transcript's intron trees. Two introns have an end at most 12 above their start: 38600644 - 38600643 = 1 and 38600097 - 38600094 = 3. Each sets the frameshift-intron flag on the transcript (lines 1032-1035).

The exon lookup follows. The bases beside the insertion are 38610911 and 38610912. With the flag set the stretch is 12, and line 865 widens the span on both sides, from 38610911 - 12 = 38610899 to 38610912 + 12 = 38610924. Exon 3 ends at 38610910, inside that span, so `exon` is set. Without the stretch the span is 38610911 to 38610912, which holds no base of an exon.

With `exon` set, the CDS coordinates of the insertion as written are fetched. Both of its bases are intron bases, so the list is one gap, neither one coordinate nor empty, and lines 555-557 set `coding`. The flags are kept on the row.

`hgvs_protein` then shifts the insertion 3', one base toward the exon, to the gap between 38610911 and 38610910. At the shifted place `map_insert` finds an exon base beside the insertion, cDNA 869, and puts the insertion directly before it, between cDNA 868 and 869. That is before coding base 76, so the shifted position has a first residue, 26. `coding` is set, and lines 1730-1738 let the notation be written.

The alternate coding sequence is built with the inserted base at that point: the 75 bases before it, the base, and the 687 after it, 75 + 1 + 687 = 763 bases. One base more shifts the frame. The codon of residue 26 opens with the inserted `G` in place of its own `C`, and proline becomes alanine. The shifted frame ends at a stop 27 codons on: `p.Pro26AlafsTer27`.

On ENST00000538115 the 2-base intron sets the flag, since 100371477 - 100371476 = 1. The bases beside the second insertion are 100333339 and 100333340, and the widened span runs from 100333339 - 12 = 100333327 to 100333340 + 12 = 100333352. It holds 100333343, the first base of exon 2. The shift carries the insertion three bases on, to the gap before 100333343. `map_insert` places it before cDNA 202, coding base 202 - 149 + 1 = 54. Three bases joined in there, 53 + 3 + 7,819 = 7,875, keep the frame and repeat residue 18: `p.Arg18dup`.

The rows without an `HGVSp` stop earlier. ENST00000309117 has no short intron, so its stretch is 0. The span 38610911 to 38610912 holds no base of its exon 3, `exon` stays unset, lines 535-537 flag the allele `non_coding`, and `hgvs_protein` returns nothing. ENST00000546213 gets `exon` from the stretch, but the junction lies outside its coding region. Lines 523-525 flag the allele `utr`, and `coding` is not set.

### 1.5 Why vep-rs's row is the correct one

The insertion lies between 38610911 and 38610912. Intron 2 runs from 38610911 to 38612796, 38612796 - 38610911 + 1 = 1,886 bases, so both bases are intron bases. Exon 3 ends one base lower, at 38610910. The 3' rule moves the copy one base, to the gap between the intron's last base and the exon's first. The described change is a duplication of the intron's last base, `c.76-1dup`, and both engines write it. It names one intron position and no coding base.

The change touches the acceptor site whichever way the copy is placed. Taken as the intron's new last base, it leaves `G` and `G` before the exon where the acceptor's `A` and `G` stood, in an intron of 1,886 + 1 = 1,887 bases. Taken as the exon's first base, it leaves the acceptor whole and changes the base the exon opens with, in an exon of 150 + 1 = 151 bases and a coding sequence of 762 + 1 = 763. Ensembl VEP's `HGVSp` is the second outcome stated as the consequence. The record does not show which splice the cell makes.

The HGVS Nomenclature (`hgvs-nomenclature.org/stable/`) does not choose between the two. Its Splicing page, `recommendations/RNA/splicing/`, gives `r.spl` for "variants changing the +1, +2, -2 and -1 position of an intron", where "RNA was not analysed, but splicing is most probably affected". It gives `r.(spl?)` where "normal splicing might be affected as a consequence of variants in the first or last nucleotide of an exon". For the protein level of a variant that most probably affects splicing, the same page answers: "The best format seems to use p.?, meaning "I do not know what to expect on the protein level"."

The protein page for substitutions, `recommendations/protein/substitution/`, describes its own splicing example by `p.?` as well: "it can not be excluded that the variant affects splicing, having unknown consequences". The general recommendations, `recommendations/general/`, make the protein level optional: "Descriptions on the RNA and/or protein level may be given in addition."

Neither engine writes `p.?` on a row of this record, which writes the insertion inside the intron. Of the two values an engine does write, a definite notation or none, the empty field is the one that agrees with those pages: it states nothing the record does not show. vep-rs writes it on this record. Ensembl VEP writes it too on ENST00000309117, ENST00000476950 and ENST00000398998, where the same record stands one base before the same exon.

What separates the row with a notation from the three without is an intron 10,266 positions away, which the record does not touch. The comment on the widening names another case, a variant that falls in a short intron. Here the widening changes the answer for a record in an intron of 1,886 bases.

The second record repeats the argument on the forward strand. Its bases 100333339 and 100333340 lie in intron 1, which runs from 100331845 to 100333342: 100333342 - 100331845 + 1 = 1,498 bases. The described change, `c.54-3_54-1dup`, names three intron positions and no coding base. Placed in the intron, the copy keeps the acceptor's two bases and puts a second `A` and `G` three bases before them. Placed in the exon, it adds a residue. Ensembl VEP states the second on three of the four transcripts with the 2-base intron. It writes nothing on the fourth, where the junction lies in the 5' UTR, or on the six without such an intron.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --fasta tests/golden/116/GRCh37-hgvs/reference.fa.gz \
  --species homo_sapiens --assembly GRCh37 --hgvs
grep -E 'ENST00000288304|ENST00000309117' out.txt
```

The committed GRCh37 cache holds the twelve transcripts of the record, and the committed FASTA holds chromosome 21. vep-rs writes no `HGVSp` on either row. Ensembl VEP 116.3 takes `--hgvs` with the cache options of [Reproducing these examples](reproducing.md), its GRCh37 cache and the GRCh37 FASTA. It writes `p.Pro26AlafsTer27` on the first row and no `HGVSp` on the second.

The entry's reproduction, `scripts/divergences/reproductions/hgvsp_coding_test_widened_by_a_short_intron.toml`, runs both records on both engines. No committed cache holds the chromosome 7 transcripts, so it copies the two region files of Ensembl VEP's GRCh37 cache that hold the records, `21/38000001-39000000` and `7/100000001-101000000`, and converts them with `tools/vep-storable-to-json`. It holds seven rows to their recorded values: the four that differ, and the rows of ENST00000309117, ENST00000348028 and ENST00000546213, which do not.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`tools/vep-formats-compare`) has no class of this shape. A disagreement of the shape is reported as a difference in `HGVSp` and charged to vep-rs.

The terms of such a row stay equal and are no part of this entry. Both engines apply the stretch when they set `exon` for the consequence terms. The rule of `splice_polypyrimidine_tract_variant` asks for `exon => 0` (`Bio/EnsEMBL/Variation/Utils/Constants.pm`, lines 977-980), so the term is left off the four rows of the second record whose transcripts have the 2-base intron, in both engines, and written on the other six.

Rows beside the shape, of a change the record does not write inside the intron, are no part of this entry, whatever either engine writes on them. The same change written at the junction, `21 38610910 at_junction_21_38610910 G GC`, has the exon's first base beside it. Ensembl VEP writes `frameshift_variant`, a CDS and a protein position and an `HGVSp` for it on all four coding transcripts. A variant inside a short intron, the case the comment on the widening names, is another shape.
