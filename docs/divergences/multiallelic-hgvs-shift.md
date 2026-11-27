# The later alleles of a multi-allelic insertion record under `--hgvs`

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP annotates the second and later alleles of a multi-allelic insertion record at the coding and protein positions the preceding allele's HGVS shift left cached on the object the alleles share. Their consequence terms, CDS position, protein position, codons, amino acids and HGVSp can therefore differ from the annotation of the same alleles written one per record. A plugin field read from those amino acids, the Blosum62 score, differs with them ([section 3](#3-a-record-whose-first-allele-gives-the-second-a-blosum62-score)). vep-rs keeps each allele's shift with that allele, so every row is annotated at its own allele's coordinates.

**Why this is a defect.** An allele's consequence depends on the allele and the transcript alone. Ensembl VEP's own annotation of each allele written on its own record, and of the whole record under `--shift_hgvs 0`, reads the allele at the positions, codons and amino acids vep-rs reads it at. Only the multi-allelic record under the default shifting differs, and only on the later alleles. The HGVS shift is a notation convention that moves where an insertion is written; it is not meant to move the allele the consequence is computed for.

**What vep-rs does instead.** vep-rs keeps each allele's shift with that allele. Every row of a record is annotated at its own allele's coordinates, whatever the record's other alleles and whatever `--shift_hgvs` says.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1724-1733, `hgvs_protein`: the alleles of one record share one TranscriptVariation per transcript; when an allele's HGVSp is written, the shared object's `translation_start`, `translation_end` and `_translation_coords` are deleted and recomputed with that allele's 3' shift, and the recomputed values stay cached.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 365-380, `translation_start`: a shifted coordinate the transcript mapper returns as a gap is stored as an undefined position, and a stored value, defined or not, is returned to every later caller.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 385-404 and 1691-1779: the rows are written one allele at a time in the record's ALT order, each row complete before the next begins and HGVSp its last field. A later allele's row therefore reads what the preceding allele's HGVSp left on the shared object.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 798-830, `codon`: the codons and amino acids of a row are read at the positions the shared object holds.
- The same file, lines 844-858 and 2404-2410, `codon` and `_get_alternate_cds`: where the shifted position has protein positions, the shared object's `cds_start` and `cds_end` are deleted and read again with the allele's shift, and that pair stays cached too.

**What correcting it would change.** The class vanishes: every allele of a multi-allelic record is annotated at its own shifted position, as the same allele is on a record of its own.

**The class it produces.** `hgvs_shift_shared_across_alleles`, documented in a corpus's manifest for the rows of a later allele whose consequence set or fields differ for this reason. The golden harness and the format comparator apply it from the manifest; no committed corpus holds such a row ([section 1.7](#17-how-the-comparator-treats-this-class)). With the candidate cache the format comparator also decides four columns of a later allele's row without one: `CDS_position`, `Protein_position`, `Amino_acids` and `Codons` ([section 2.7](#27-how-the-comparator-treats-this-class)). A class of its own, `blosum62_hgvs_shift_shared_across_alleles`, takes the Blosum62 plugin's score on such a row ([section 3.7](#37-how-the-comparator-treats-this-class)).

A consequence-set difference of the class is a pair the concordance figures would charge to vep-rs. The dataset measurements behind those figures run both engines without plugins and without `--hgvs`, so no such pair arises in them and none is set aside. The runs that measure the plugins verified on real data, on the genome-wide ClinVar files and on the population datasets, load those plugins under `--hgvs`. They are scored on plugin fields alone, never on a consequence set.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A record with two insertion alleles at a stop codon

**In short.** A mouse record inserts `AAAGGG` or `AAAG` after the same base, at the stop codon of a transcript on the forward strand. Written on its own, or under `--shift_hgvs 0`, the second allele is an in-frame insertion at residue 69 with codons `taa/taAAAGa` on both engines: `inframe_insertion,coding_sequence_variant,NMD_transcript_variant` on Ensembl VEP's row and `inframe_insertion,stop_retained_variant,NMD_transcript_variant` on vep-rs's, a difference of one term that [`stop_codon_overlap_insertion_range`](stop-codon-overlap-insertion-range.md) documents. On one record under `--hgvs`, in a run without statistics, Ensembl VEP writes the second allele as `frameshift_variant,NMD_transcript_variant` with no protein position, codons or amino acids. vep-rs writes the second allele's own row in every case.

### 1.1 What a reader needs to know first

A multi-allelic VCF record lists several alternate alleles on one line, and both engines write one set of rows per allele. Under `--hgvs`, an insertion or deletion is shifted to the most 3' position along the transcript at which the altered sequence is the same before its HGVS notation is written, and `HGVS_OFFSET` reports the distance moved. Each allele of a record has its own shift, since the alleles differ in length and bases.

### 1.2 The record

The record is a public mouse variant from the European Variation Archive, release 8, against the Ensembl 116 GRCm39 cache (the `GRCm39-mouse` golden corpus, `tests/golden/116/GRCm39-mouse/`, carries the cache and the FASTA slice; this record is not among its records):

```
19	42047557	rs265482355	A	AAAAGGG,AAAAG	.	.	EVA_8;TSA=insertion
```

Both alleles insert after `19:42047557`, inside the stop codon of `ENSMUST00000405524`, a transcript of `ENSMUSG00000025176` on the forward strand, where the first allele inserts six bases and the second four. The same alleles written one per record are:

```
19	42047557	rs265482355_a1	A	AAAAGGG	.	.	EVA_8;TSA=insertion
19	42047557	rs265482355_a2	A	AAAAG	.	.	EVA_8;TSA=insertion
```

Both engines run with `--hgvs`.

### 1.3 What each engine writes

The second allele's row on `ENSMUST00000405524`, default format. Ensembl VEP, the two alleles on one record:

```
rs265482355	19:42047557-42047558	AAAG	ENSMUSG00000025176	ENSMUST00000405524	Transcript	frameshift_variant,NMD_transcript_variant	597-598	206-207	-	-	-	-	IMPACT=HIGH;STRAND=1;HGVSc=ENSMUST00000405524.1:c.207_*1insAAGA;HGVS_OFFSET=1
```

Ensembl VEP, the second allele on a record of its own:

```
rs265482355_a2	19:42047557-42047558	AAAG	ENSMUSG00000025176	ENSMUST00000405524	Transcript	inframe_insertion,coding_sequence_variant,NMD_transcript_variant	597-598	206-207	69	*/*KX	taa/taAAAGa	-	IMPACT=MODERATE;STRAND=1;HGVSc=ENSMUST00000405524.1:c.207_*1insAAGA;HGVS_OFFSET=1
```

vep-rs, the two alleles on one record:

```
rs265482355	19:42047557-42047558	AAAG	ENSMUSG00000025176	ENSMUST00000405524	Transcript	inframe_insertion,stop_retained_variant,NMD_transcript_variant	597-598	206-207	69	*/*KX	taa/taAAAGa	-	IMPACT=MODERATE;STRAND=1;HGVSc=ENSMUST00000405524.1:c.207_*1insAAGA;HGVS_OFFSET=1
```

The first allele's row on the transcript is the same in both engines and in both of Ensembl VEP's runs, with its own shift of one base, `HGVS_OFFSET=1`, and no HGVSp:

```
rs265482355	19:42047557-42047558	AAAGGG	ENSMUSG00000025176	ENSMUST00000405524	Transcript	inframe_insertion,stop_retained_variant,NMD_transcript_variant	597-598	206-207	69	*/*KG	taa/taAAAGGGa	-	IMPACT=MODERATE;STRAND=1;HGVSc=ENSMUST00000405524.1:c.207_*1insAAGGGA;HGVS_OFFSET=1
```

Both alleles shift one base, since the base after the insertion point, `c.207`, is an `A` like each allele's first inserted base and the next, `c.*1`, is not. On this record and two more multi-allelic insertion records of mouse chromosome 19, rs3389480233 and rs3408733889, Ensembl VEP's multi-allelic rows differ from its own one-per-record rows on 24 rows. vep-rs's rows equal its own one-per-record rows on every one of the 201 rows the three records produce, and Ensembl VEP's one-per-record rows on 199 of them. The other two are the second allele's rows of this record, on this transcript and on ENSMUST00000577435, which differ in the one term [`stop_codon_overlap_insertion_range`](stop-codon-overlap-insertion-range.md) documents.

### 1.4 Why Ensembl VEP writes that

The alleles of a record share one TranscriptVariation per transcript, and that object caches the protein positions a row reads. The rows are written one allele at a time, in the record's ALT order, each row complete before the next begins, with the HGVSp its last field. The first allele's row is computed from positions read fresh, residue 69.

Its HGVSp computation then deletes the shared object's `translation_start` and `translation_end` and recomputes them at the allele's shift of one base. Shifted, the insertion's first base is `c.*1`, past the stop codon, where the transcript mapper returns a gap, so the recomputed protein position is undefined and stays cached as such; finding no position, the computation writes no HGVSp on the row.

The second allele's row reads the shared object next. `codon` finds no protein position, so there is no reference codon to compare and the in-frame test fails. The frameshift test, which reads the unshifted CDS positions, fires on a four-base insertion over no deleted base, and the row prints no protein position, codons or amino acids. Written on its own record, the allele is the first and only one the object serves, and its positions are read fresh.

On this record the two engines' cDNA and CDS positions agree. The first allele's computation returned before it read the CDS pair again, because its shifted position has no protein position (`TranscriptVariationAllele.pm`, lines 1730-1738). Where it has one, the shared CDS pair moves too ([section 2](#2-a-record-whose-first-allele-moves-the-coding-positions-of-the-second)).

Source (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1724-1733: the shared protein positions deleted and recomputed with the calling allele's shift.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 365-380: a gap stored as an undefined position and returned to every later caller.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 385-404 and 1691-1779: the rows written one allele at a time, HGVSp last.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 798-830: the codons read at the shared positions.

### 1.5 Why vep-rs's row is the correct one

The allele decides. Inserting `AAAG` after `19:42047557` puts four bases inside the stop codon `taa`, giving `taAAAGa`: the codon `taA` is kept as a stop and the insertion runs into the 3' UTR. That is an in-frame insertion at the coding sequence's last codon, residue 69, which is what both engines write for the allele on its own record. There, and on the two-allele record under `--shift_hgvs 0`, which switches the shifting off, Ensembl VEP's row has the position, codons and amino acids of vep-rs's row and differs from it in one consequence term, `coding_sequence_variant` for vep-rs's `stop_retained_variant`, which [`stop_codon_overlap_insertion_range`](stop-codon-overlap-insertion-range.md) documents. The presence of a six-base sibling allele on the same line changes none of the second allele's bases.

### 1.6 Try it yourself

Annotate the record of 1.2 and its two one-allele twins with `--hgvs` against the GRCm39 cache and FASTA. vep-rs's command is below; Ensembl VEP takes the same records with its mouse cache under `--offline --cache --no_stats`, the statistics setting this entry's reproduction runs with:

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --species mus_musculus --assembly GRCm39 --hgvs
```

| Run | Ensembl VEP, second allele | vep-rs, second allele |
| --- | --- | --- |
| the two alleles on one record | `frameshift_variant,NMD_transcript_variant`, no protein position | `inframe_insertion,stop_retained_variant,NMD_transcript_variant`, residue 69 |
| the second allele alone | `inframe_insertion,coding_sequence_variant,NMD_transcript_variant`, residue 69 | `inframe_insertion,stop_retained_variant,NMD_transcript_variant`, residue 69; [`stop_codon_overlap_insertion_range`](stop-codon-overlap-insertion-range.md) documents the stop term |
| the two alleles, `--shift_hgvs 0` | the same terms and residue as the allele alone | the same terms and residue as the allele alone |

The statistics setting decides where the second allele's terms are read. Under `--no_stats` an allele's consequences are first evaluated while its own row is built, at the positions the shared object then holds, which gives the row of 1.3. A run that keeps statistics evaluates every allele's consequences before the record's rows, at the record's own positions ([section 2.4](#24-why-ensembl-vep-writes-that)). In both, the positions a row prints come from the shared object.

### 1.7 How the comparator treats this class

The golden harness (`crates/vep-cli/tests/common/mod.rs`) and the format comparator (`tools/vep-formats-compare`) apply `hgvs_shift_shared_across_alleles` from a corpus manifest. A `divergences` entry of the class names a later allele's row with the consequence set vep-rs writes, and a `field_divergences` entry a field with the value vep-rs prints; the comparison checks vep-rs's row against the documented reading in place of Ensembl VEP's.

No predicate decides a consequence set of the class from the two outputs, since the rows do not state which allele's shift set their positions. The four columns of section 2 and the plugin score of section 3 are decided from the candidate cache.

No committed corpus holds a row of the class, so no manifest names it. The entry's row of [Intended divergences from Ensembl VEP](../intended-divergences.md) counts the corpora that hold a record of two or more insertion or deletion alleles, the shape a record needs before the difference can occur. The GRCm39 and TAIR10 golden corpora hold such records and annotate them under `--hgvs`; they leave out the records on which a later allele's rows differ. The oracle's population cases under `--hgvs` carry one alternate allele per record, so none has the shape.

## 2. A record whose first allele moves the coding positions of the second

**In short.** `21 47545369 . A AC,ACC`, a 1000 Genomes Phase 3 record on GRCh37, inserts one `C`, or two, after coding base 478 of ENST00000413758 (COL6A2, forward strand), in front of seven more `C`. Under `--hgvs` both notations move seven bases: `HGVS_OFFSET=7`. On the second allele's row Ensembl VEP 116.3 prints `CDS_position` `485-486`, `Protein_position` `162` and `Amino_acids` `P/TX`. vep-rs prints `478-479`, `160` and `T/TX`, the record's own position, as both engines do on the first allele's row.

### 2.1 What a reader needs to know first

An insertion into a run of one base can be written after any base of the run, and HGVS notation writes it after the last. `HGVSc`, `HGVSp` and `HGVS_OFFSET` describe that moved position. `CDS_position`, `Protein_position`, `Amino_acids` and `Codons` describe the record where the input writes it.

The coding sequence of ENST00000413758 in the cache is 750 bases, incomplete at both ends (`FLAGS=cds_start_NF,cds_end_NF`). Its codons 160 to 162 read `ACC CCC CCA`: threonine, proline, proline. Coding bases 479 to 485 are the seven `C`.

### 2.2 The record

```
21	47545369	.	A	AC,ACC	100	PASS	.
```

The record is of the 1000 Genomes Phase 3 call set on GRCh37; its INFO column is left out here. 21:47545369 is coding base 478, the first base of codon 160. Both alleles insert between coding bases 478 and 479: the first one `C`, the second two. Both engines run with `--hgvs`.

### 2.3 What each engine writes

The rows on ENST00000413758, default format, with the `Extra` column cut to the fields this example reads. The first allele's row, the same in both engines:

```
21_47545369_A/AC/ACC	21:47545369-47545370	C	ENSG00000142173	ENST00000413758	Transcript	frameshift_variant	478-479	478-479	160	T/TX	acc/aCcc	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_start_NF,cds_end_NF;HGVSc=ENST00000413758.1:c.485dup;HGVSp=ENSP00000395751.1:p.Asp163ArgfsTer3;HGVS_OFFSET=7
```

The second allele's row, Ensembl VEP 116.3:

```
21_47545369_A/AC/ACC	21:47545369-47545370	CC	ENSG00000142173	ENST00000413758	Transcript	frameshift_variant	478-479	485-486	162	P/TX	acc/aCCcc	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_start_NF,cds_end_NF;HGVSc=ENST00000413758.1:c.484_485dup;HGVSp=ENSP00000395751.1:p.Asp163GlnfsTer?;HGVS_OFFSET=7
```

vep-rs:

```
21_47545369_A/AC/ACC	21:47545369-47545370	CC	ENSG00000142173	ENST00000413758	Transcript	frameshift_variant	478-479	478-479	160	T/TX	acc/aCCcc	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_start_NF,cds_end_NF;HGVSc=ENST00000413758.1:c.484_485dup;HGVSp=ENSP00000395751.1:p.Asp163GlnfsTer?;HGVS_OFFSET=7
```

On the second allele's row the two engines agree in `cDNA_position`, `Codons`, `HGVSc`, `HGVSp` and `HGVS_OFFSET`, and differ in the three columns named above.

### 2.4 Why Ensembl VEP writes that

The first allele's row prints the record's own positions, read fresh. Its HGVSp computation then deletes the shared object's protein positions and reads them again seven bases on, in codon 162. `codon` and `_get_alternate_cds` do the same with the CDS pair, which then reads `485-486`. The computation also deletes the peptide of the reference allele, one object for the record's alleles, and reads it again at codon 162: proline.

The second allele's row prints what the shared object holds: `CDS_position` `485-486` and `Protein_position` `162`. Its `Amino_acids` joins the reference allele's peptide, the proline of codon 162, with the allele's own peptide, `TX`.

That peptide stands at codon 160 because it was read before any row. The run shown kept statistics (`--stats_text`), and such a run logs each record's most severe consequence before its rows. That evaluates the consequences of every allele, whose tests read the allele's peptide at the record's own position, and the allele keeps it. The row's `Codons` reads codon 160 as well, since each allele keeps the codon it read first.

The rows of section 1 are of a run without statistics (`--no_stats`). There an allele's consequences are first evaluated on its own row, at the positions the shared object then holds.

Source (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1724-1733: the shared protein positions deleted and read again with the calling allele's shift.
- The same file, lines 844-858 and 2404-2410: the shared CDS pair deleted and read again with that shift, and left there.
- The same file, lines 1780-1783 and 1793: the reference allele's peptide deleted and read again at the shifted positions.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1742, 1746 and 1750: `Amino_acids`, `CDS_position` and `Protein_position` printed from the allele objects and the shared object.
- `OutputFactory.pm`, lines 347 and 1022, and `Bio/EnsEMBL/VEP/Stats.pm`, line 241: the record's most severe consequence logged before its rows.
- `Bio/EnsEMBL/Variation/VariationFeature.pm`, lines 1000-1002, and `Utils/VariationEffect.pm`, line 785: every allele's consequences evaluated, reading its peptide.

### 2.5 Why vep-rs's row is the correct one

The record's own coordinates decide. An insertion after 21:47545369 lies between coding bases 478 and 479, inside codon 160, whose residue is threonine: `CDS_position` `478-479`, `Protein_position` `160`, `Amino_acids` `T/TX`. Ensembl VEP prints exactly that on the first allele's row, which inserts at the same point. Its own `Codons` column on the second allele's row still reads codon 160, `acc/aCCcc`. `P/TX` pairs a reference residue of codon 162 with an alternate reading of codon 160, and no one position of the transcript gives both.

### 2.6 Try it yourself

Annotate the record of 2.2 with `--hgvs`, against the GRCh37 cache and FASTA of [Reproducing these examples](reproducing.md):

```
vep -i record.vcf -o out.txt --offline \
  --json_cache /path/to/json_cache/homo_sapiens/116_GRCh37 \
  --fasta /path/to/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
  --assembly GRCh37 --hgvs
```

vep-rs writes `478-479`, `160` and `T/TX` on the second allele's row, and the same when that allele stands on a record of its own. Ensembl VEP 116.3, run on the two-allele record against its own GRCh37 cache with `--stats_text`, writes `485-486`, `162` and `P/TX`.

### 2.7 How the comparator treats this class

With the candidate cache (`--json-cache`) the format comparator decides `CDS_position`, `Protein_position`, `Amino_acids` and `Codons` without a manifest. It must know the record's alternate alleles in order. A row name such as `21_47545369_A/AC/ACC` spells them, and a row under an identifier takes them from its record's line in the input VCF (`--input`). Without either, the row stays a difference. The row's allele must be the second or a later one, and Ensembl VEP's row of the first allele on the transcript must carry a non-zero `HGVS_OFFSET`.

vep-rs's `CDS_position` must be the record's span on the coding sequence, and Ensembl VEP's that span moved by the first allele's offset, both mapped through the cache's transcript model. Each `Protein_position` must be the codons of its `CDS_position`. The two `Amino_acids` must share their alternate part, each reference part the translation of its side's codons; a vep-rs insertion between two codons has the reference part `-`. That is the run with statistics shown above.

Without statistics, Ensembl VEP reads an in-frame insertion's codons at the moved pair too. It inserts the record's own allele into the window from the moved codon (`TranscriptVariationAllele.pm` lines 841-843 and 877). It upper-cases the allele from the own pair's codon position (lines 911-923), and it keeps the own pair's reference codons. On such a row each side's `Amino_acids` must be the translation of its window, each `Codons` its window, and Ensembl VEP's reference codons vep-rs's. `chr1 149390802 A AGGC,AGGCGGCGGC` on ENST00000650865 (forward strand) gives `G/GAAA` and `-/GGGGCGGCGgca` against `-/GGG` and `-/GGCGGCGGC`. Records inside a codon and on the reverse strand (ATXN1) were checked against Ensembl VEP 116.3 the same way.

A first allele's row, a first allele without an offset, a move of another size, a vep-rs value that is not the record's own, every other column and a transcript whose mapping the cache does not hold stay differences. The row's `DownstreamProtein` falls under [`plugin_hgvs_shifted_cds_coordinates`](plugin-hgvs-shifted-cds-coordinates.md#9-a-later-allele-of-a-multi-allelic-record). A `BLOSUM62` that Ensembl VEP alone writes on a later allele's row falls under the class of [section 3.7](#37-how-the-comparator-treats-this-class).

## 3. A record whose first allele gives the second a Blosum62 score

**In short.** `11 1092804 . C CAATA,CCCA`, a 1000 Genomes Phase 3 record on GRCh37, inserts four bases, or three, after coding base 4626 of ENST00000359061 (MUC2, forward strand), the last base of codon 1542. Under `--hgvs` the first allele's notation moves one base, `HGVS_OFFSET=1`; the second allele's does not move. On the second allele's row Ensembl VEP 116.3 prints `Amino_acids` `T/P`, the residue of codon 1543 over the inserted one, and its Blosum62 plugin scores that pair: `BLOSUM62=-1`. vep-rs prints `-/P`, an insertion between two codons, and writes no `BLOSUM62`.

### 3.1 What a reader needs to know first

The Blosum62 plugin writes `BLOSUM62`, a substitution score, on a row whose `Amino_acids` is one residue over another ([Blosum62](../plugins/transcript-model.md#blosum62)). Three bases inserted between two codons add a residue and replace none: the row's `Amino_acids` has `-` for its reference part, and the plugin writes nothing.

Codons 1542 and 1543 of ENST00000359061 both read `ACC`, threonine. Coding base 4626 is the last base of codon 1542, and coding base 4627 the first of codon 1543.

### 3.2 The record

```
11	1092804	.	C	CAATA,CCCA	100	PASS	.
```

The record is of the 1000 Genomes Phase 3 call set on GRCh37; its INFO column is left out here. 11:1092804 is coding base 4626. Both alleles insert between coding bases 4626 and 4627: the first `AATA`, the second `CCA`. Coding base 4627 is an `A`, like the first allele's first base and unlike the second's, so the first allele's notation alone moves. Both engines run with `--hgvs` and the Blosum62 plugin.

### 3.3 What each engine writes

The rows on ENST00000359061, default format, with the `Extra` column cut to the fields this example reads. The first allele's row, the same in both engines:

```
11_1092804_C/CAATA/CCCA	11:1092804-1092805	AATA	ENSG00000198788	ENST00000359061	Transcript	frameshift_variant	4626-4627	4626-4627	1542-1543	-/NX	-/AATA	-	IMPACT=HIGH;STRAND=1;HGVSc=ENST00000359061.5:c.4627_4628insATAA;HGVSp=ENSP00000351956.5:p.Thr1543AsnfsTer?;HGVS_OFFSET=1
```

The second allele's row, Ensembl VEP 116.3:

```
11_1092804_C/CAATA/CCCA	11:1092804-1092805	CCA	ENSG00000198788	ENST00000359061	Transcript	inframe_insertion	4626-4627	4627-4628	1543	T/P	-/CCA	-	IMPACT=MODERATE;STRAND=1;HGVSc=ENST00000359061.5:c.4626_4627insCCA;HGVSp=ENSP00000351956.5:p.Thr1543Pro;BLOSUM62=-1
```

vep-rs:

```
11_1092804_C/CAATA/CCCA	11:1092804-1092805	CCA	ENSG00000198788	ENST00000359061	Transcript	inframe_insertion	4626-4627	4626-4627	1542-1543	-/P	-/CCA	-	IMPACT=MODERATE;STRAND=1;HGVSc=ENST00000359061.5:c.4626_4627insCCA;HGVSp=ENSP00000351956.5:p.Thr1542_Thr1543insPro
```

The second allele's two rows agree in `cDNA_position`, `Codons` and `HGVSc`, and neither carries an `HGVS_OFFSET`. They differ in `BLOSUM62` and in `CDS_position`, `Protein_position`, `Amino_acids` and `HGVSp`, columns this page's lead names. The second allele's rows on ENST00000441003, where the point lies between codons 1541 and 1542, differ the same way.

### 3.4 Why Ensembl VEP writes that

The first allele's row prints the record's own positions, read fresh. Its HGVSp computation then deletes the shared object's protein positions and reads them again one base on, inside codon 1543, and the CDS pair moves with them to `4627-4628` ([section 2.4](#24-why-ensembl-vep-writes-that)). It also deletes the peptide of the reference allele, one object for the record's alleles, and reads it again at codon 1543: threonine.

The second allele's row prints what the shared object holds. Its `Amino_acids` joins the reference allele's peptide with the allele's own: threonine over proline. The allele's own peptide was read at the record's own position before any row, since the run shown kept statistics (`--stats_text`). The second allele's notation does not move, so its own HGVSp computation deletes nothing, and the pair stands.

The plugins run once a row is built. Blosum62 reads the same joined pair, finds one residue over one, and writes the matrix's score for threonine over proline, -1.

Source (`VEP_plugins` and `ensembl-variation` at release/116, `ensembl-vep` at release/116.3):

- `Blosum62.pm`, line 110: a value only where the peptide pair is one capital letter over one. Lines 112-118 return the pair's score from the matrix of lines 50-71, in the residue order of line 73, as lines 84-88 index it.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 618-630, `pep_allele_string`: the reference allele's peptide, a slash, the allele's own.
- The same file, lines 1724-1733: the shared protein positions deleted and read again where the calling allele's shift is not zero. Lines 1780-1783 and 1793: the reference allele's peptide deleted there and read again.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1742: `Amino_acids` printed from `pep_allele_string`; line 398: the plugins run once the row is built.
- The statistics lines of section 2.4: every allele's peptide read before the record's rows.

### 3.5 Why vep-rs's row is the correct one

The record's own coordinates decide. 11:1092804 is coding base 4626, the last base of codon 1542, so the three inserted bases stand between codons 1542 and 1543 as a codon of their own, `CCA`, proline. The allele replaces no residue: `Amino_acids` is `-/P`, and a substitution score has no pair to score. Ensembl VEP's row of the first allele, which inserts at the same point, has `-` for its reference part, and its own `Codons` column on the second allele's row reads `-/CCA`, no reference codon. `T/P` pairs the threonine of codon 1543, which the allele leaves in place, with the inserted proline.

### 3.6 Try it yourself

Annotate the record of 3.2 with `--hgvs` and the plugin, against the GRCh37 cache and FASTA of [Reproducing these examples](reproducing.md):

```
vep -i record.vcf -o out.txt --offline \
  --json_cache /path/to/json_cache/homo_sapiens/116_GRCh37 \
  --fasta /path/to/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
  --assembly GRCh37 --hgvs --plugin Blosum62
```

vep-rs writes `-/P` and no `BLOSUM62` on the second allele's rows of ENST00000359061 and ENST00000441003. Ensembl VEP 116.3, run with `--hgvs`, `--stats_text` and the plugin against its own GRCh37 cache, writes `T/P` and `BLOSUM62=-1` on both.

### 3.7 How the comparator treats this class

The comparator takes the score under a class of its own, `blosum62_hgvs_shift_shared_across_alleles`. The class reads `BLOSUM62` alone. It needs the candidate cache (`--json-cache`) and the record's alternate alleles in order, as [section 2.7](#27-how-the-comparator-treats-this-class) describes.

Ensembl VEP alone must write the value, on a later allele's row whose `CDS_position` and `Protein_position` section 2.7 takes. Ensembl VEP's `Amino_acids` must be one capital letter over one: the first the translation of the moved pair's codon, the second the alternate part of vep-rs's own `Amino_acids`. vep-rs's own `Amino_acids` must not have that form, and the value must be the matrix's score for Ensembl VEP's pair.

A first allele's row, a first allele without an offset, another value or residue, a value on vep-rs's side and a run without the cache stay differences. The class has no row of its own in [Intended divergences from Ensembl VEP](../intended-divergences.md): it is this page's mechanism met on a plugin's field, and the row of `hgvs_shift_shared_across_alleles` there counts the shape.
