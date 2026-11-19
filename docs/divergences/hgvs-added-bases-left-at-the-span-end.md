# The HGVSc of bases an allele adds to a reference it keeps, at a transcript's end

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.2 moves inserted bases to their most 3' position only on a record classed as an insertion, one whose alleles keep no reference base. The bases that an allele of a multi-allelic record adds to reference bases it keeps stay at the end of the record's span and are written as an insertion there. Where that end is a reverse-strand transcript's last base, the insertion is numbered one base past it, `c.*1871+1_*1871insA`, and the HGVSIntronOffset plugin writes the offsets 1 and 0 beside it. There vep-rs places the added bases by the 3' rule and writes the duplication of the last bases they repeat, `c.*1871dup`, with both offsets 0.

**Why this is a defect.** A record of several alternates can hold an allele that keeps the reference bases and adds bases to them. That allele makes the same sequence as an insertion given in a record of its own, and Ensembl VEP writes two notations for the one change. Given alone, the inserted bases are moved along the run of bases they repeat and named a duplication where they copy the bases before them. Given beside another allele, they are left at the end of the record's span, and named a duplication only where the whole alternate is the kept bases written twice.

The HGVS Nomenclature (`hgvs-nomenclature.org/stable/`) gives a change one description. Its general recommendations, `recommendations/general/`, state the rule: "3'rule: for all descriptions, the most 3' position possible of the reference sequence is arbitrarily assigned to have been changed." They add: "the 3'rule also applies for changes in single residue stretches and tandem repeats (nucleotide or amino acid)." Its Duplication page, `recommendations/DNA/duplication/`, says under Notes: "when a variant can be described as a duplication, it must be described as a duplication and not as, e.g., an insertion". None of these sentences turns on the other alleles of a record.

The coordinate Ensembl VEP writes past a last base, such as `c.*1871+1`, is an intron's form at a place where the transcript has no intron. The Numbering page, `background/numbering/`, says: "it is not allowed to describe variants in nucleotides beyond the boundaries of a transcript reference sequence, using that transcript reference sequence."

**What vep-rs does instead.** Added bases that repeat the bases beside them give the same sequence at every point of that run of bases. vep-rs places the added bases by the 3' rule where the run of such points reaches the transcript's first or last base:

- The run ends directly past the last base, and the added bases repeat the transcript's last bases. vep-rs writes the duplication of those last bases.
- The run starts at or before the first base and ends between two bases of the transcript. vep-rs writes the notation at that last point: a duplication where the copied bases all lie in the transcript, an insertion where the added bases copy nothing.
- The added bases stand directly past the last base and can be written at no point inside the transcript. vep-rs writes no HGVSc.

`HGVS_OFFSET` states the count of bases the rule moved the added bases, negative on a reverse-strand transcript, and is absent where nothing moved. The HGVSIntronOffset plugin's two fields are the intron parts of the coordinates vep-rs writes, and absent where it writes no HGVSc. Where the run reaches neither end of the transcript, both engines write the same notation, the added bases at the end of the record's span.

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `VEP_plugins` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 120-123, `_return_3prime`: the 3' shift returns at once unless the record's class is `insertion` or `deletion`. `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 362-380, gives a record that class only where no reference base is left and every alternate is sequence, or where every alternate is empty. A record with an allele that keeps a reference base has another class, and none of its alleles is shifted.
- `TranscriptVariationAllele.pm`, lines 1474-1483, `hgvs_transcript`: the whole alternate allele is passed with the record's span, and `Utils/Sequence.pm` line 526 reads the reference bases of that span.
- The same file, lines 2204-2242, `_clip_alleles`: the bases the two alleles share at their start are stripped and the start is moved past them, then the bases they share at their end. Lines 2271-2284 re-type an emptied reference as `ins`, its start one above its end. On the transcript's strand, the added bases so stand after as many leading bases of the span as the alternate begins with: past the span's last base where that is the whole span, before the span's first base where it is none.
- No duplication is looked for. `Utils/Sequence.pm` line 570 compares the neighbouring bases only for a span of no bases, the span of a record that is an insertion. The clip's own test, at lines 2274-2279, needs the notation's numbering, which lines 1548-1553 set after the clip.
- `TranscriptVariationAllele.pm`, lines 2800-2827, `_get_cDNA_position`: the position after a reverse-strand transcript's last base lies before the first exon in genomic order. The loop reaches its intron branch at index 0, line 2802 reads the exon at index -1, and lines 2820-2826 build the coordinate off the edge exon, its last base with `+1`.
- The same file, lines 1530-1533, 1543 and 1545: the intron parts of the two coordinates are stored, 1 and 0. `HGVSIntronOffset.pm`, lines 62-72, writes them beside whatever notation the call returns.

**What correcting it would change.** Ensembl VEP would describe the added bases of such an allele as it describes the same bases given alone: at the most 3' point of their run, and as a duplication where they copy the bases before them. No coordinate one base past a last base would be written, and the plugin's start offset 1 would go with it.

**The class it produces.** `hgvs_added_bases_left_at_the_span_end`, on a row whose allele holds more bases than its Location spans, where Ensembl VEP prints an HGVSc of the row's transcript and no `HGVS_OFFSET`. The format comparator holds the rule (`scripts/concordance/compare_vep_formats.py`, `added_bases_left_at_the_span_end`) and judges the two engines' rows together. Three forms are the class:

- Past a reverse-strand transcript's last base, against a duplication. Ensembl VEP's notation is an `ins` of the added bases between a coordinate and the same coordinate with `+1`, beside the plugin's pair 1 and 0. vep-rs's is a `dup` of as many bases ending at that coordinate, with no `HGVS_OFFSET` and the pair 0 and 0. The last added base repeats the kept base beside it. The HGVSc and the start offset are set aside.
- At the same place, against no notation. vep-rs writes no HGVSc, no `HGVS_OFFSET` and neither field, and the last added base does not repeat the kept base beside it. The HGVSc and both fields are set aside.
- Carried along a run. Ensembl VEP's notation is an `ins` between two neighbouring exonic bases, or a `dup`. vep-rs's is a `dup` of as many bases, or an `ins` of the same bases turned, further 3' by the count of bases its `HGVS_OFFSET` states. In every form, bases are counted (`exonic_run`) only between exonic coordinates of one region (5' UTR, coding sequence, 3' UTR): a duplication or a carry over an intron offset or two regions stays a difference. Both rows write the same plugin offsets. The HGVSc and the `HGVS_OFFSET` are set aside. Where nothing moves, vep-rs's `dup` names the bases that end at Ensembl VEP's `ins`, beside no `HGVS_OFFSET`, and the HGVSc alone is set aside.

The class does not take a row whose allele holds no more bases than its Location spans, a row with an `HGVS_OFFSET` on Ensembl VEP's side, or a pair other than 1 and 0 past the last base. It does not take Ensembl VEP's other form, an insertion numbered one base before a forward-strand transcript's first base beside the pair 0 and -1. On that row `hgvs_clipped_reference_against_whole_alternate` sets the end offset aside, the fourth shape of [how a comparison treats the offset fields](hgvs-clipped-reference-intron-offset-fields.md), and takes the HGVSc where its own coordinate test holds, and the `HGVS_OFFSET` stays a difference.

A notation of another transcript, numbering or kind stays a difference, and so does a count of bases that is not the allele's surplus over its span. Where neither row carries a plugin field (`without_plugin`), the HGVSc alone is the class past a last base; along a run, HGVSc and `HGVS_OFFSET`.

The class is built into the format comparator, and the golden corpus builder and the test harness name it. The plugin tuple scorer's class table declares it for the HGVSIntronOffset plugin's two fields. The concordance figures score consequence sets and never see these fields. No manifest of a committed corpus names the class.

**How often the shape occurs.** Tier 2, `--hgvs` being needed. The shape is a record of several alternates with an allele that keeps its reference bases and adds bases to them. The transcript decides whether the run of the added bases reaches one of its ends. 1,188 of the 15,646,481 records of the chromosome 21 cuts carry the shape, all in the 1000 Genomes Phase 3 cut, the one cut of the five whose records hold several alternates. 22 of the 111 committed corpora hold such a record.

Counted by row on chromosomes 2, 3, 5, 10, 13, 15, 17 and 19 of 1000 Genomes Phase 3, 27 rows of 19 records are rows of such an allele whose run reaches a transcript's first or last base with a placement inside. 10 lie at a last base and 17 at a first base, and chromosome 13 holds none. On 18 of the 27 the engines differ in one of the class's forms. 3 are an `ins` past a reverse-strand transcript's last base against a `dup`, 7 are carried along a run on a forward-strand transcript, and 8 are an `ins` against a `dup` with nothing moved, on a reverse-strand transcript.

On 8 more of the 27 both engines write the same `dup`, and on 1 Ensembl VEP writes its other form, whose HGVSc and `HGVS_OFFSET` stay differences. The row of section 3, where vep-rs writes no notation, is not among the 27. The same count finds no such row on chromosome 6 of 1000 Genomes high-coverage or in four files of gnomAD v4.1 genomes records, of chromosomes 6, 16, 19 and X. The genome-wide ClinVar files hold none: each of their 4,388,172 GRCh37 and 4,398,127 GRCh38 records has one alternate.

## 1. One base added at the last base of ENST00000337288

**In short.** A 1000 Genomes Phase 3 record on GRCh37 at 2:96850596 has three alternates. The first keeps the four reference bases at 96850597 to 96850600 and adds a fifth `T`. ENST00000337288, on the reverse strand, ends at 96850597. Its last eight bases are a run of `T` on the genome, `A` on the transcript's strand. Ensembl VEP 116.2 writes `HGVSc=ENST00000337288.5:c.*1871+1_*1871insA`, `HGVS_IntronStartOffset=1` and `HGVS_IntronEndOffset=0`. vep-rs writes `HGVSc=ENST00000337288.5:c.*1871dup` and both offsets 0.

### 1.1 What a reader needs to know first

Both engines drop the first base of a VCF record when every allele shares it and some allele differs in length from the reference. Each alternate is then held against the reference bases that remain. A record of one alternate is trimmed further, to its minimal pair, so an added base becomes an insertion with no reference base. Beside other alternates it stays as written: here five bases against four.

ENST00000337288 (ENSG00000084090) has eight exons on the reverse strand, 96850597 to 96874556 on the genome, and 3,368 spliced bases. Its coding sequence ends at base 1,497, so its last base, 96850597, is c.*1871. The Location 2:96850597-96850600 is c.*1868 to c.*1871, the transcript's last four bases.

The HGVSIntronOffset plugin writes two fields beside a transcript notation, the intron parts of its start and end coordinates: `45` for `c.123+45`, `0` for an exonic coordinate.

### 1.2 The record

```
2	96850596	.	CTTTT	CTTTTT,CTTT,C	.	.	.
```

Both engines run `--hgvs --plugin HGVSIntronOffset` on GRCh37 with a reference FASTA. The base at 96850596 is shared by all four alleles and dropped.

### 1.3 What each engine writes

Ensembl VEP 116.2, default format, the first alternate's row on ENST00000337288:

```
2_96850596_CTTTT/CTTTTT/CTTT/C	2:96850597-96850600	TTTTT	ENSG00000084090	ENST00000337288	Transcript	3_prime_UTR_variant	3365-3368	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000337288.5:c.*1871+1_*1871insA;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=1
```

vep-rs:

```
2_96850596_CTTTT/CTTTTT/CTTT/C	2:96850597-96850600	TTTTT	ENSG00000084090	ENST00000337288	Transcript	3_prime_UTR_variant	3365-3368	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000337288.5:c.*1871dup;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

Every other column is equal. The record's other two alternates delete one base and all four, and both engines write `c.*1871del` and `c.*1868_*1871del` for them, with the offsets 0 and 0.

### 1.4 Why Ensembl VEP writes that

Four reference bases are left after the shared base is dropped, so the record's class is neither `insertion` nor `deletion`, and line 123 returns before any shift. The alternate, five bases on the transcript's strand, is held against the transcript's last four. `_clip_alleles` strips the four leading bases the two share and moves the start past them, to the position after the last base. One alternate base is left against an empty reference, re-typed `ins` and held from that position to the last base.

`_get_cDNA_position` turns the position after the last base into the genomic base 96850596, which lies before the first exon in genomic order. The loop reaches its intron branch at index 0, line 2802 reads the exon at index -1, and lines 2820-2826 build the coordinate off the last exon: its last base with `+1`, printed `*1871+1` in the 3' UTR's numbering. Lines 1530-1533 split the two coordinates into exon and intron parts, lines 1543 and 1545 store 1 and 0, and the plugin writes the pair.

### 1.5 Why vep-rs's row is the correct one

The alternate keeps the bases at 96850597 to 96850600 and holds one more of the same base. 96850597 to 96850604 are one run of that base, c.*1864 to c.*1871 of the transcript, and the base at 96850596, outside the transcript, is another. The added base gives the same sequence at each of the nine points of that run. By the 3' rule quoted above, the change is assigned to the last of them on the transcript's strand, directly after c.*1871.

There the added base is a copy of c.*1871. The Duplication page defines a duplication as "a sequence change where, compared to a reference sequence, a copy of one or more nucleotides is inserted directly 3' of the original copy of that sequence", and requires that such a change be described as one. The notation is `c.*1871dup`. It names one base, which the transcript holds, and the coordinate is exonic, so both intron offsets are 0.

`c.*1871+1` names no base of the transcript. The Numbering page gives the `+` form to "nucleotides at the 5' end of an intron", numbered "relative to the last nucleotide of the directly upstream exon", and no intron follows the last exon. The start offset 1 is counted at 96850596, outside the transcript.

The same base given in a record of its own, directly past the last base, is a row outside the transcript. [`hgvs_added_bases_written_outside_the_transcript`](hgvs-added-bases-written-outside-the-transcript.md) is the entry for that row.

### 1.6 Try it yourself

Save the record of 1.2 as `record.vcf` and annotate it against a GRCh37 cache with the reference FASTA. vep-rs's command is below; Ensembl VEP takes the same record with its own cache and the FASTA under `--offline --cache`:

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh37 --hgvs --plugin HGVSIntronOffset
```

vep-rs's row for the five-base allele on ENST00000337288 carries `c.*1871dup` and both offsets 0. Ensembl VEP 116.2's carries the insertion at `c.*1871+1_*1871` and the offsets 1 and 0. The reproduction `scripts/divergences/reproductions/hgvs_added_bases_left_at_the_span_end.toml` runs this record and the record of section 4 on both engines.

### 1.7 How the comparator treats this class

On this row the rule accepts the HGVSc and `HGVS_IntronStartOffset`, `hgvs_intronstartoffset` in JSON. The allele holds one base more than the Location spans. Ensembl VEP's notation is an `ins` of that base between `*1871+1` and `*1871` beside the pair 1 and 0, on a row whose `STRAND` is -1. vep-rs's is a `dup` of one base ending at `*1871`, with no `HGVS_OFFSET` and the pair 0 and 0. The end offset is 0 on both rows and is compared like any other field.

## 2. Two and three bases added at the last base of ENST00000578576

**In short.** A 1000 Genomes Phase 3 record on GRCh37 at 17:13399486 has four alternates. Three keep the reference base at 13399487 and add one, two and three bases, each a `T`, and the fourth deletes it. 13399487 is the last base of ENST00000578576, on the reverse strand, and the last of a run of ten, `A` on the transcript's strand. For the alleles that add two and three bases, Ensembl VEP 116.2 writes `c.*27+1_*27insAA` and `c.*27+1_*27insAAA`, each with the offsets 1 and 0. vep-rs writes `c.*26_*27dup` and `c.*25_*27dup` with 0 and 0.

### 2.1 What a reader needs to know first

ENST00000578576 (ENSG00000153976) has two exons on the reverse strand, 13399487 to 13459360 on the genome, and 750 spliced bases. Its coding sequence ends at base 723, so its last base, 13399487, is c.*27. The run of ten bases is 13399487 to 13399496, c.*18 to c.*27.

### 2.2 The record

```
17	13399486	.	AT	ATT,ATTT,ATTTT,A	.	.	.
```

Both engines run `--hgvs` and the HGVSIntronOffset plugin on GRCh37 with a reference FASTA, in a run that loads further plugins. The base at 13399486 is shared by all five alleles and dropped, and the Location is 17:13399487.

### 2.3 What each engine writes

Ensembl VEP 116.2, default format, the alleles of three and four bases on ENST00000578576:

```
17_13399486_AT/ATT/ATTT/ATTTT/A	17:13399487	TTT	ENSG00000153976	ENST00000578576	Transcript	3_prime_UTR_variant	750	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000578576.1:c.*27+1_*27insAA;AA=T;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=1;LOEUF=1.1360e+00;NearestExonJB=ENSE00002235077|0|start|649;gnomAD_100x_cov=0.0000;gnomAD_10x_cov=0.9984;gnomAD_15x_cov=0.9642;gnomAD_1x_cov=1.0000;gnomAD_20x_cov=0.8698;gnomAD_25x_cov=0.7209;gnomAD_30x_cov=0.5167;gnomAD_50x_cov=0.0317;gnomAD_5x_cov=1.0000;gnomAD_mean_cov=30.4430;gnomAD_median_cov=30.0000;pHaplo=0.575305175361985;pTriplo=0.198646110735125
17_13399486_AT/ATT/ATTT/ATTTT/A	17:13399487	TTTT	ENSG00000153976	ENST00000578576	Transcript	3_prime_UTR_variant	750	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000578576.1:c.*27+1_*27insAAA;AA=T;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=1;LOEUF=1.1360e+00;NearestExonJB=ENSE00002235077|0|start|649;gnomAD_100x_cov=0.0000;gnomAD_10x_cov=0.9984;gnomAD_15x_cov=0.9642;gnomAD_1x_cov=1.0000;gnomAD_20x_cov=0.8698;gnomAD_25x_cov=0.7209;gnomAD_30x_cov=0.5167;gnomAD_50x_cov=0.0317;gnomAD_5x_cov=1.0000;gnomAD_mean_cov=30.4430;gnomAD_median_cov=30.0000;pHaplo=0.575305175361985;pTriplo=0.198646110735125
```

vep-rs:

```
17_13399486_AT/ATT/ATTT/ATTTT/A	17:13399487	TTT	ENSG00000153976	ENST00000578576	Transcript	3_prime_UTR_variant	750	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000578576.1:c.*26_*27dup;AA=T;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0;LOEUF=1.1360e+00;NearestExonJB=ENSE00002235077|0|start|649;gnomAD_100x_cov=0.0000;gnomAD_10x_cov=0.9984;gnomAD_15x_cov=0.9642;gnomAD_1x_cov=1.0000;gnomAD_20x_cov=0.8698;gnomAD_25x_cov=0.7209;gnomAD_30x_cov=0.5167;gnomAD_50x_cov=0.0317;gnomAD_5x_cov=1.0000;gnomAD_mean_cov=30.4430;gnomAD_median_cov=30.0000;pHaplo=0.575305175361985;pTriplo=0.198646110735125
17_13399486_AT/ATT/ATTT/ATTTT/A	17:13399487	TTTT	ENSG00000153976	ENST00000578576	Transcript	3_prime_UTR_variant	750	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000578576.1:c.*25_*27dup;AA=T;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0;LOEUF=1.1360e+00;NearestExonJB=ENSE00002235077|0|start|649;gnomAD_100x_cov=0.0000;gnomAD_10x_cov=0.9984;gnomAD_15x_cov=0.9642;gnomAD_1x_cov=1.0000;gnomAD_20x_cov=0.8698;gnomAD_25x_cov=0.7209;gnomAD_30x_cov=0.5167;gnomAD_50x_cov=0.0317;gnomAD_5x_cov=1.0000;gnomAD_mean_cov=30.4430;gnomAD_median_cov=30.0000;pHaplo=0.575305175361985;pTriplo=0.198646110735125
```

The two notations and the start offset are the only differences: the fields of the run's other plugins stand on both engines' rows with the same values. For the allele of two bases both engines write `c.*27dup` with 0 and 0. Its whole alternate is the kept base written twice, which `Utils/Sequence.pm` lines 602-606 name a duplication before the clip. For the deletion both write `c.*27del`.

### 2.4 Why the rows differ

As in section 1, no allele is shifted, the clip strips the kept base, and the two or three added bases are left after the last base, at `*27+1`. They are copies of the transcript's last two and three bases, c.*26 to c.*27 and c.*25 to c.*27, and stand directly 3' of them. The Duplication page's sentences of 1.5 make each a duplication, `c.*26_*27dup` and `c.*25_*27dup`. Ensembl VEP writes the duplication for one added base and an insertion for two and three, at one base of one transcript.

## 3. An added base that repeats nothing at the last base of ENST00000260442

**In short.** A 1000 Genomes Phase 3 record on GRCh37 at 15:52401459 has two alternates. The first keeps the base at 52401460, the last base of the reverse-strand transcript ENST00000260442, and adds one base past it. The added base is not the base at 52401460. Ensembl VEP 116.2 writes `HGVSc=ENST00000260442.3:c.*585+1_*585insT`, `HGVS_IntronStartOffset=1` and `HGVS_IntronEndOffset=0`. vep-rs writes no HGVSc and neither field.

### 3.1 What a reader needs to know first

ENST00000260442 (ENSG00000137875) has two exons on the reverse strand, 52401460 to 52404972 on the genome, and 1,249 spliced bases. Its coding sequence ends at base 664, so its last base, 52401460, is c.*585.

### 3.2 The record

```
15	52401459	.	GT	GAT,G	.	.	.
```

Both engines run `--hgvs` and the HGVSIntronOffset plugin on GRCh37 with a reference FASTA, in a run that loads further plugins. The base at 52401459 is shared by all three alleles and dropped, and the Location is 15:52401460.

### 3.3 What each engine writes

Ensembl VEP 116.2, default format, the first alternate's row on ENST00000260442:

```
15_52401459_GT/GAT/G	15:52401460	AT	ENSG00000137875	ENST00000260442	Transcript	3_prime_UTR_variant	1249	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000260442.3:c.*585+1_*585insT;AA=T;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=1;LOEUF=1.7350e+00;NearestExonJB=ENSE00001015560|0|start|711;gnomAD_100x_cov=0.0000;gnomAD_10x_cov=0.9992;gnomAD_15x_cov=0.9943;gnomAD_1x_cov=1.0000;gnomAD_20x_cov=0.9626;gnomAD_25x_cov=0.8609;gnomAD_30x_cov=0.6306;gnomAD_50x_cov=0.0407;gnomAD_5x_cov=0.9992;gnomAD_mean_cov=33.2150;gnomAD_median_cov=32.0000;pHaplo=0.416589592254294;pLI_gene_value=0.00;pTriplo=0.14887513770265
```

vep-rs:

```
15_52401459_GT/GAT/G	15:52401460	AT	ENSG00000137875	ENST00000260442	Transcript	3_prime_UTR_variant	1249	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;AA=T;LOEUF=1.7350e+00;NearestExonJB=ENSE00001015560|0|start|711;gnomAD_100x_cov=0.0000;gnomAD_10x_cov=0.9992;gnomAD_15x_cov=0.9943;gnomAD_1x_cov=1.0000;gnomAD_20x_cov=0.9626;gnomAD_25x_cov=0.8609;gnomAD_30x_cov=0.6306;gnomAD_50x_cov=0.0407;gnomAD_5x_cov=0.9992;gnomAD_mean_cov=33.2150;gnomAD_median_cov=32.0000;pHaplo=0.416589592254294;pLI_gene_value=0.00;pTriplo=0.14887513770265
```

The HGVSc and the two offsets are the only differences. The record's other alternate deletes the last base, and both engines write `c.*585del` for it with the offsets 0 and 0.

### 3.4 Why vep-rs writes no notation

Ensembl VEP's row comes about as in section 1: the clip leaves the added base after the last base, numbered `*585+1`. That base stands between 52401459 and 52401460. It differs from the base at 52401460, so it can be written at no point inside the transcript, and it is no copy of the transcript's last base.

An insertion is named by the two bases it lies between. The Insertion page, `recommendations/DNA/insertion/`, says under Notes: "`positions flanking` should contain two flanking nucleotides, e.g., `123_124`, not `123_125`." One of the two is 52401459, which the transcript does not hold, and the Numbering page's sentence above forbids describing a change there with this transcript. The change has no notation of ENST00000260442. vep-rs writes none, and the plugin's fields stand beside an HGVSc and nowhere else.

## 4. A base carried five bases to the last base of ENST00000318578

**In short.** A 1000 Genomes Phase 3 record on GRCh37 at 15:39547041 has two alternates: a substitution of the base there, and an allele that keeps it and adds a `A` after it. ENST00000318578, on the forward strand, ends five bases on, at 39547046, and those five bases are all `A`. Ensembl VEP 116.2 writes `HGVSc=ENST00000318578.3:c.*2154_*2153insA` and no `HGVS_OFFSET`. vep-rs writes `HGVSc=ENST00000318578.3:c.*2158dup` and `HGVS_OFFSET=5`. Both write the offsets 0 and 0.

### 4.1 What a reader needs to know first

ENST00000318578 (ENSG00000175746) has two exons on the forward strand, 39542885 to 39547046 on the genome, and 3,078 spliced bases. Its coding sequence ends at base 920, so 39547041 is c.*2153 and the last base, 39547046, is c.*2158.

`HGVS_OFFSET` is the count of bases by which the 3' shift moved a notation. Both engines print it beside the HGVSc of an insertion or deletion that the shift moved.

### 4.2 The record

```
15	39547041	.	T	A,TA	.	.	.
```

Both engines run `--hgvs --plugin HGVSIntronOffset` on GRCh37 with a reference FASTA. The two alternates do not share their first base, so no base is dropped and the Location is 15:39547041.

### 4.3 What each engine writes

Ensembl VEP 116.2, default format, the second alternate's row on ENST00000318578:

```
15_39547041_T/A/TA	15:39547041	TA	ENSG00000175746	ENST00000318578	Transcript	3_prime_UTR_variant	3073	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000318578.3:c.*2154_*2153insA;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

vep-rs:

```
15_39547041_T/A/TA	15:39547041	TA	ENSG00000175746	ENST00000318578	Transcript	3_prime_UTR_variant	3073	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000318578.3:c.*2158dup;HGVS_OFFSET=5;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

The HGVSc and the `HGVS_OFFSET` are the only differences. For the substitution both engines write `c.*2153T>A`.

### 4.4 Why Ensembl VEP writes that

One alternate is a substitution, so the record's class is neither `insertion` nor `deletion`, and line 123 returns before any shift. The clip strips the kept base and leaves an `ins` directly after it, between c.*2153 and c.*2154, the first point of the run. No duplication is looked for there. Lines 1535-1538 order a notation's two coordinates for printing unless the second carries a star. This one does, so the 3' flank is printed first: `*2154_*2153`.

### 4.5 Why vep-rs's row is the correct one

39547042 to 39547046, c.*2154 to c.*2158, are five of the added base, and 39547047, outside the transcript, is another base. The added base gives the same sequence directly after 39547041 and after each of the five: six points. The 3' rule assigns the change to the last of them, directly after c.*2158, where the added base is a copy of the transcript's last base: `c.*2158dup`. The rule carried it five bases from the point the record writes, which `HGVS_OFFSET=5` states.

For a record that is an insertion alone, Ensembl VEP makes the same move. Lines 1441-1449 shift the insertion along its run, and `Utils/Sequence.pm` lines 570-597 name the copy a duplication.

### 4.6 The same run with nothing to move

A record can write the added bases at the last point of their run. At 17:7835440, on GRCh37, an allele keeps the base there, c.-183 of the reverse-strand transcript ENST00000572656, and adds a copy of the base before it, c.-184, the transcript's first base. The run of that base starts before the transcript and ends between c.-184 and c.-183, where the record writes the added base. Ensembl VEP 116.2 writes `c.-184_-183insC`, and vep-rs `c.-184dup`. Neither prints an `HGVS_OFFSET`, and both write the offsets 0 and 0.
