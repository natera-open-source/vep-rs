# Intended divergences from Ensembl VEP

This page documents every disagreement between vep-rs and Ensembl VEP that the adjusted concordance figure sets aside, and shows for each one why it is set aside.

vep-rs 0.3.0 reports two concordance figures against Ensembl VEP release 115.2. **Raw concordance** keeps every tuple both engines write: a tuple is one output row reduced to its Location, Allele, Feature, Feature_type and Consequence set, and raw F1 counts the tuples the two engines share against the tuples each writes alone, with nothing set aside. **Adjusted concordance** sets aside the tuples of the divergence classes documented on this page and nothing else.

A class is one shape of disagreement. Each class on this page has been traced to a specific mechanism in one engine's source, reproduced on a public record with a command a reader can run, and judged. Eleven of the twelve classes are the output of five defects in Ensembl VEP's source that vep-rs does not reproduce (classes 8 and 11 are also changes in vep-rs's own behaviour between 0.2.0 and 0.3.0). The twelfth, `<CNV:TR>`, is a difference of representation in which neither engine is wrong. A disagreement that is not one of these classes is charged to vep-rs in both figures, whatever its cause.

The page is organised by defect. Each defect's section opens with the defect in one sentence, why it is a defect, what vep-rs does instead, where it is in Ensembl VEP's source and what correcting it would change. The classes it produces follow.

Each class is a walkthrough of one real record. It gives the background a reader needs, the record, both engines' rows, why Ensembl VEP writes what it writes, why vep-rs's row is the correct one, what a reader can reproduce, and how the comparator scripts recognise the class. The accompanying paper's supplement documents five classes at the paper's pinned build, four of them defects: they correspond to defects A to D and the representation difference. The two classes vep-rs 0.3.0 stopped reproducing, 8 and 11, are the two the paper does not document.

Grouped by the place in Ensembl VEP's source that produces them, the eleven defect classes are five defects. Two classes belong to one defect when the same wrong lines produce both, so that correcting those lines changes both.

| Defect | What goes wrong in Ensembl VEP | Where | Classes | What correcting it would change |
| --- | --- | --- | --- | --- |
| A | The splice-region verdict is stored by assignment, so a later differing region overwrites an earlier one | `ensembl-variation`, `BaseTranscriptVariationAllele.pm` line 215 | 1, 2 | both classes vanish: `splice_region_variant` is kept wherever any differing region earned it |
| B | `start_lost` and `start_retained_variant` are computed by two independent tests, and nothing keeps both from being true of one allele | `ensembl-variation`, `Utils/VariationEffect.pm` lines 851-963 and 1028-1075 | 3, 4 | both classes vanish: one term per allele, the one the record supports |
| C | A record marked skipped loads no transcript regions of its own but stays in its batch and is annotated against whatever its neighbours loaded | `ensembl-vep`, `AnnotationSource.pm` line 238 against `InputBuffer.pm` lines 284-330 | 5, 6, 7, 8 | the batch dependence vanishes: Ensembl VEP writes for a skipped record what it writes for it alone, the same from any file; class 8 then agrees with vep-rs, and classes 5 to 7 differ only in that vep-rs annotates an oversize or unsupported record in full where Ensembl VEP skips it |
| D | Transcripts are matched to records by coordinate alone, with no chromosome, and a record without a parsed mate is never chromosome-checked | `ensembl-vep`, `InputBuffer.pm` lines 345-366; `ensembl-variation`, `StructuralVariationOverlap.pm` lines 64-71 | 9, 10 | both classes vanish: a record is matched only to transcripts on its own chromosome, and its `intergenic_variant` row is written when there are none |
| E | A mate-side row is created from the mate's coordinate, but every predicate except `feature_truncation` is evaluated at the record's own coordinate | `ensembl-variation`, `BaseVariationFeatureOverlapAllele.pm` lines 257 and 273 | 11 | the class vanishes: a mate-side row carries the terms the mate breakpoint supports |
| none | A `<CNV:TR>` allele is expanded to literal sequence and annotated as an insertion or deletion, where vep-rs keeps it symbolic; neither reading is wrong | `ensembl-vep`, `Parser/VCF.pm` lines 382-431 | 12 | nothing to correct; the comparator pairs the two readings |

Source citations name the Ensembl repository, module and line at release 115 (the `release/115` branch of `ensembl-variation` and of `ensembl-vep`). The comparator scripts of this repository are `scripts/concordance/compare_vep_outputs.py` for the SNP/indel datasets and `scripts/validation/compare_sv_concordance.py` for the structural-variant sets. Every row on this page is in Ensembl VEP's default tab output, whose columns are:

```
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
```

Coordinates follow VEP's output conventions: a symbolic deletion's Location starts one base after the VCF POS (the padding base), a breakend's Location is the base after the breakpoint, and an insertion's Location names the two flanking bases.

The twelve classes, with the defect each belongs to:

| # | Class | Defect | Kind | Datasets affected | Side set aside | Count |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Splice-region term lost to an assignment order | A | Ensembl VEP defect | ClinVar GRCh37, ClinVar GRCh38 | both members of each pair | 440, 1,444 pairs |
| 2 | Class 1 on the structural-variant sets | A | Ensembl VEP defect | synthetic `02_mnp_complex` (GRCh38) | both members of each pair | 10 pairs |
| 3 | `start_lost` and `start_retained_variant` emitted together | B | Ensembl VEP defect | ClinVar GRCh37 and GRCh38, gnomAD v2.1.1 and v4.1 chr21, 1000 Genomes high-coverage chr21 | both members of each pair | 261, 928, 3, 77, 1 pairs |
| 4 | Class 3 on the structural-variant sets | B | Ensembl VEP defect | ClinVar SV, gnomAD-SV, 1000 Genomes SV and the synthetic `05_symbolic_del_ins`, `07_cnv_repeat`, `11_multi_allelic`, `12_complex_imprecise` | both members of each pair | 3 GRCh37 and 24 GRCh38 pairs |
| 5 | Structural variants above `--max_sv_size` annotated against whichever transcripts the batch loaded | C | Ensembl VEP defect | gnomAD-SV v2.1 chr21 (GRCh37), gnomAD-SV v4.1 chr21 (GRCh38) | vep-rs (transcripts Ensembl VEP never loaded); Ensembl VEP where it names a transcript vep-rs does not | 10,934 vep-rs tuples; 183,964 vep-rs and 1,823 Ensembl VEP tuples |
| 6 | The mate end of a breakend above `--max_sv_size` | C | Ensembl VEP defect | gnomAD-SV v2.1 chr21, gnomAD-SV v4.1 chr21 | vep-rs | 56 (GRCh37), 310 (GRCh38) tuples |
| 7 | `<CPX>`, class 5 reaching an unsupported allele type | C | Ensembl VEP defect | gnomAD-SV v4.1 chr21 | vep-rs | 218 tuples on 3 records |
| 8 | `<NON_REF>` is not a variant | C | Ensembl VEP defect; vep-rs 0.3.0 behaviour change | synthetic special-allele set (`10_special_alleles`), both assemblies | Ensembl VEP transcript rows and the vep-rs `intergenic_variant` row | 1,838 + 180 (GRCh37), 5,403 + 200 (GRCh38) rows |
| 9 | Transcripts from another chromosome | D | Ensembl VEP defect | gnomAD-SV v2.1 chr21, gnomAD-SV v4.1 chr21, synthetic breakend set (`09_breakends`) | Ensembl VEP | 1,838 (GRCh37), 22,481 (GRCh38) tuples |
| 10 | The orphaned `intergenic_variant` on the other side of class 9 | D | Ensembl VEP defect | gnomAD-SV v2.1 chr21, gnomAD-SV v4.1 chr21 | vep-rs | 155 (GRCh37), 518 (GRCh38) rows |
| 11 | Breakend mate rows read from the other chromosome's coordinate | E | Ensembl VEP defect; vep-rs 0.3.0 behaviour change | synthetic `09_breakends` and `13_vcf45_features`, gnomAD-SV v4.1 chr21 | both members of each pair | 1,407 (GRCh37), 7,258 (GRCh38) pairs |
| 12 | `<CNV:TR>` read symbolically, not expanded | none | representation difference | synthetic tandem-repeat set (`07_cnv_repeat`), both assemblies | both members of each pair | 1,368 (GRCh37), 2,644 (GRCh38) pairs |

How to read the counts:

- Setting a tuple aside is also called **masking**: a class is fully masked when its rule sets aside every tuple of the class, and a rule is shape-limited when it sets aside only pairs of one exact shape.
- The **six SNP/indel datasets** are ClinVar GRCh37 and GRCh38 (`clinvar_grch37.vcf.gz`, `clinvar_grch38.vcf.gz`), gnomAD v2.1.1 chr21 (GRCh37, `gnomad_genomes_v2.1.1_chr21.vcf.gz`), gnomAD v4.1 chr21 (GRCh38, `gnomad_genomes_v4.1_chr21.vcf.gz`), 1000 Genomes Phase 3 chr21 (GRCh37, `1kg_integrated_chr21.vcf.gz`) and 1000 Genomes high-coverage chr21 (GRCh38, `1kg_highcov_chr21.vcf.gz`), each scored as one file.
- The **two structural-variant sets** are one pool per assembly of sixteen files, scored file by file and pooled: the thirteen synthetic files under `tests/sv_validation/<assembly>/` plus `clinvar_sv_chr21.vcf.gz`, the gnomAD-SV file (`gnomad_sv_v2.1_chr21.vcf.gz` for GRCh37, `gnomad_sv_v4.1_chr21.vcf.gz` for GRCh38) and `1kg_sv_chr21.vcf.gz`.
- Counts are tuples or pairs measured with the comparators of this repository on vep-rs 0.3.0 against Ensembl VEP 115.2's output for the same inputs.
- A **pair** is one tuple from each engine that the comparator matches as two readings of one record on one transcript. For classes 1 to 4 and 11 the two tuples share Location, Allele, Feature and Feature_type and differ in the consequence set. For class 12 they share Feature and Feature_type only, because the two engines write the same record at a different Location and Allele (12.3).
- Counts are per input file. The comparator also reports a pooled count per assembly in which a tuple that recurs in two files is counted once (the gnomAD-SV and 1000 Genomes files share 1 GRCh37 and 49 GRCh38 records at identical coordinates). Where the two differ this page gives the per-file count and names the pooled one (classes 4 and 10).
- The classes overlap: a tuple two rules both match is set aside once. Class 7's 218 tuples are inside class 5's 183,964, and class 5's 1,823 Ensembl VEP tuples are class 9 tuples, so the Count column does not sum to the number of tuples set aside.
- On the two ClinVar datasets classes 1 and 3 together are every tuple the adjusted figure sets aside. On the other four SNP/indel datasets class 3 is. On the two structural-variant sets classes 2 and 4 to 12 together are, and the adjusted figure is 1.0 on both.

## Defect A: a later differing region erases the splice-region term

**The defect in one sentence.** Ensembl VEP decides the splice-region question once per changed place of a variant and keeps only the last answer, so a change inside the splice-region window loses its `splice_region_variant` term whenever a later change outside the window follows it.

**Why this is a defect.** A consequence term describes where a changed base falls. Whether base 3 of an exon is in the splice-region window does not depend on whether the same variant also changes base 5. Ensembl VEP agrees when the bases are changed one at a time: it writes `splice_region_variant` for the first change alone, and only the presence of the second change removes the term. An annotation that an unrelated second change can switch off is not an annotation of the first change.

**What vep-rs does instead.** vep-rs keeps the term when any differing region earns it. That is how Ensembl VEP already treats its other thirteen splice flags, which are only ever switched on and never overwritten; the splice-region slot is the one exception, and the exception is the defect.

**Where it is in Ensembl VEP's source** (`ensembl-variation`):
- `Bio/EnsEMBL/Variation/BaseTranscriptVariationAllele.pm`, `_intron_effects`, line 215: the verdict is stored by assignment, `$intron_effects->{splice_region} = _intron_overlap(...)`, inside the `foreach my $region` loop that starts at line 117.
- `Utils/VariationEffect.pm`, lines 87-110: `_intron_overlap` returns a definite 0 or 1 for each region, so a later 0 replaces an earlier 1.
- `VariationFeatureOverlapAllele.pm`, line 405: `_get_differing_regions` splits a multi-base change into the regions the loop visits.

**What correcting it would change.** Classes 1 and 2 would vanish: Ensembl VEP would keep `splice_region_variant` wherever any region earned it, which is what vep-rs writes.

**The classes it produces.**
- Class 1: the ClinVar GRCh37 and GRCh38 datasets, 440 and 1,444 pairs; both members of each pair are set aside.
- Class 2: the same shape on the structural-variant sets, scored by the second comparator, 10 pairs on the synthetic complex-substitution set; both members of each pair are set aside.

### 1. Splice-region term lost to an assignment order

**In short.** ClinVar variation 630437 changes two bases of BRCA2 exon 6 at once: base 3, which is inside the splice-region window, and base 5, which is not. Ensembl VEP writes `missense_variant` alone; vep-rs writes `missense_variant,splice_region_variant`. The changed base 3 is in the window by the term's own definition, and Ensembl VEP itself writes the term when that base is changed alone, so vep-rs's row is the correct one.

#### 1.1 What a reader needs to know first

Three Sequence Ontology terms appear in this class:

- `splice_region_variant` (SO:0001630): "a sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron". Ensembl VEP applies exactly that window on each side of every intron boundary.
- `missense_variant` (SO:0001583): "a sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved".
- `synonymous_variant` (SO:0001819): "a sequence variant where there is no resulting change to the encoded amino acid".

A variant that changes several bases at once (a multi-nucleotide variant or a deletion-insertion) may change bases that are not adjacent. Ensembl VEP compares the reference and alternate alleles position by position and evaluates each run of differing bases, called a differing region, on its own. A term that holds for any of those regions holds for the variant.

#### 1.2 The record

ClinVar variation 630437, from the ClinVar GRCh37 VCF (`clinvar_grch37.vcf.gz` in the layout of `scripts/data/download_real_world_vcfs.sh`, which pins the ClinVar release; the variation id and alleles are stable across releases):

```
13	32900381	630437	GTA	CTT	.	.	CLNHGVS=NC_000013.10:g.32900381_32900383delinsCTT;CLNVC=Indel;GENEINFO=BRCA2:675
```

The transcript is ENST00000380152 (BRCA2, forward strand). Its exon 6 of 27 spans 32900379-32900419 and is preceded by intron 5, 32900288-32900378, so exon bases 1 to 3 are 32900379, 32900380 and 32900381.

The reference GTA at 32900381-32900383 is exon bases 3, 4 and 5. The alternate CTT changes base 3 (G to C) and base 5 (A to T) and leaves base 4 (T) as it was, so the variant is two differing regions of one base each, 32900381 and 32900383. Base 3 lies inside the 1-3 exonic window of the splice-region definition; base 5 does not.

#### 1.3 What each engine writes

```
# Ensembl VEP
630437	13:32900381-32900383	CTT	ENSG00000139618	ENST00000380152	Transcript	missense_variant	711-713	478-480	160	V/L	GTA/CTT	-	IMPACT=MODERATE;STRAND=1
# vep-rs
630437	13:32900381-32900383	CTT	ENSG00000139618	ENST00000380152	Transcript	missense_variant,splice_region_variant	711-713	478-480	160	V/L	GTA/CTT	-	IMPACT=MODERATE;STRAND=1
```

The rows agree in every column but Consequence, where vep-rs adds `splice_region_variant`. The other two BRCA2 transcripts the record touches, ENST00000530893 and ENST00000544455, carry the same difference.

#### 1.4 Why Ensembl VEP writes that

Ensembl VEP evaluates the splice-region window once per differing region and stores each verdict by assignment into a single slot, so the last region's verdict replaces the earlier ones. For this record the first region (exon base 3) sets the slot to 1 and the second region (exon base 5) sets it back to 0. By the time the consequence is read, the term the first region earned is gone.

The other thirteen splice flags computed in the same loop are only ever set to 1, so they cannot be lost this way. `splice_region_variant` is the one term the order of regions can remove.

Source (`ensembl-variation`):

- `Bio/EnsEMBL/Variation/BaseTranscriptVariationAllele.pm`, `_intron_effects`, line 215: `$intron_effects->{splice_region} = _intron_overlap(...)`, inside the `foreach my $region` loop that starts at line 117.
- `Utils/VariationEffect.pm`, lines 87-110: `_intron_overlap` returns a defined 0 or 1.
- `VariationFeatureOverlapAllele.pm`, line 405: `_get_differing_regions` supplies the regions.

#### 1.5 Why vep-rs's row is the correct one

The record changes base 32900381, the third base of exon 6. That position is inside the window the term's own definition names ("within 1-3 bases of the exon").

A term defined on the position of a changed base cannot depend on what other bases the same record changes. Ensembl VEP agrees when the base is changed alone (1.6): the same nucleotide change at the same position is a `splice_region_variant` by Ensembl VEP's own reading. Only the presence of a second changed base two positions downstream removes it.

#### 1.6 Try it yourself

Annotate the two differing bases as separate records and then together, with the Ensembl VEP command of the closing section (GRCh37 cache and FASTA):

| Input | Ensembl VEP, ENST00000380152 | vep-rs, ENST00000380152 |
| --- | --- | --- |
| `13 32900381 . G C` | `missense_variant,splice_region_variant` | `missense_variant,splice_region_variant` |
| `13 32900383 . A T` | `synonymous_variant` | `synonymous_variant` |
| `13 32900381 . GTA CTT` | `missense_variant` | `missense_variant,splice_region_variant` |

Each part alone is annotated identically by both engines. The two parts as one record lose, in Ensembl VEP only, the term the first part earned. No flag changes the outcome: this class is deterministic, unlike classes 5 to 10.

#### 1.7 How the comparator treats this class

`_is_covered_splice_region_swap` in `scripts/concordance/compare_vep_outputs.py` matches a pair when Ensembl VEP's set has nothing vep-rs's lacks and vep-rs's set adds `splice_region_variant`, alone or together with `intron_variant`. Both members of a matching pair are set aside.

The `intron_variant` allowance exists because Ensembl VEP's per-region loop steps over a region that lies in an intron of twelve bases or fewer before it sets either the intronic flag or the splice-region slot (`BaseTranscriptVariationAllele.pm` lines 139-141 and 172-174). A region skipped that way withholds both terms at once. None of the ClinVar pairs uses the allowance; every one of them adds `splice_region_variant` alone.

The rule is deliberately narrower than the family of splice differences. A pair in which Ensembl VEP lacks a splice sub-term such as `splice_donor_variant`, or in which vep-rs adds `splice_polypyrimidine_tract_variant`, is never set aside, because the assignment at line 215 can remove no term but `splice_region_variant`.

Counts:

- ClinVar GRCh37: 440 pairs; ClinVar GRCh38: 1,444 pairs, every one of that shape.
- gnomAD and 1000 Genomes datasets: none.
- Structural-variant sets: 10 pairs (class 2).

The SNP/indel comparator's summary groups every splice-family pair into descriptive buckets for reporting. The buckets do not decide exclusion; this rule does, and the summary prints beside each bucket how many of its pairs the rule excluded.

### 2. Class 1 on the structural-variant sets

**In short.** The synthetic complex substitution `synth_cpx_sub_0780`, in a GRCh38 file that the structural-variant comparator scores, replaces `AAAG` with `CA` in SOD1 exon 2. That gives two differing regions: exon base 2, inside the splice-region window, and exon bases 4 and 5, outside it. Ensembl VEP writes `frameshift_variant` alone on ENST00000270142; vep-rs writes `frameshift_variant,splice_region_variant`. This is class 1 on a structural-variant set: exon base 2 is inside the window the term names, and Ensembl VEP writes the term when that base is changed alone, so vep-rs's row is the correct one.

#### 2.1 What a reader needs to know first

The structural-variant sets are scored by a different script from the SNP/indel datasets, but they contain sequence variants too (the synthetic complex-substitution set is nothing else). Nothing new happens in either engine; the same predicate is applied by the second comparator. The terms are those of 1.1.

#### 2.2 The record

`synth_cpx_sub_0780` from the synthetic complex-substitution set `tests/sv_validation/grch38/02_mnp_complex.vcf.gz` (GRCh38):

```
21	31663791	synth_cpx_sub_0780	AAAG	CA	.	PASS	.
```

The transcript is ENST00000270142 (SOD1, forward strand, canonical, GRCh38 31659693-31668931). Its exon 2 of 5 spans 31663790-31663886 and is preceded by intron 1 (31659842-31663789), so 31663791 is exon base 2 and 31663793-31663794 are exon bases 4 and 5.

Comparing `AAAG` with `CA` position by position gives two differing regions: 31663791 (A to C) and 31663793-31663794 (the two reference bases with no alternate counterpart). Position 31663792 (A) is unchanged. Base 2 of the exon lies inside the 1-3 window; bases 4 and 5 do not.

#### 2.3 What each engine writes

```
# Ensembl VEP
synth_cpx_sub_0780	21:31663791-31663794	CA	ENSG00000142168	ENST00000270142	Transcript	frameshift_variant	151-154	74-77	25-26	ES/AX	gAAAGt/gCAt	-	IMPACT=HIGH;STRAND=1
# vep-rs
synth_cpx_sub_0780	21:31663791-31663794	CA	ENSG00000142168	ENST00000270142	Transcript	frameshift_variant,splice_region_variant	151-154	74-77	25-26	ES/AX	gAAAGt/gCAt	-	IMPACT=HIGH;STRAND=1
```

The rows agree in every column but Consequence, where vep-rs adds `splice_region_variant`. The record touches ten SOD1 transcripts with the same difference: eight coding ones (`frameshift_variant`) and two non-coding ones (`non_coding_transcript_exon_variant`).

#### 2.4 Why Ensembl VEP writes that

As 1.4: the first differing region (exon base 2) sets `splice_region` to 1 and the second (exon bases 4 and 5) sets it back to 0, on each of the ten transcripts.

#### 2.5 Why vep-rs's row is the correct one

As 1.5: base 31663791 is the second base of exon 2, inside the window the term names, and Ensembl VEP calls it so when it is changed alone.

#### 2.6 Try it yourself

Ensembl VEP 115.2, GRCh38 cache and FASTA, the command of the closing section:

- `21 31663791 . A C` alone gives ENST00000270142 `missense_variant,splice_region_variant`;
- `21 31663792 . AAG A` alone (the second differing region, Location `21:31663793-31663794`, Allele `-`) gives `frameshift_variant`;
- the full record `21 31663791 . AAAG CA` gives `frameshift_variant`.

vep-rs gives the same two rows for the two parts and `frameshift_variant,splice_region_variant` for the full record.

#### 2.7 How the comparator treats this class

`filter_registry_swap_pairs` in `scripts/validation/compare_sv_concordance.py` builds the pairs of the structural-variant comparison, one tuple on each side per Location, Allele, Feature and Feature_type. It applies to them the same two rules the SNP/indel comparator exports (`EXCLUDING_RULES` of `scripts/concordance/compare_vep_outputs.py`, imported so the two scripts cannot drift), setting both members of a matching pair aside.

Counts: the class 1 rule matches the 10 transcripts of `synth_cpx_sub_0780` (GRCh38) and nothing else. Fully masked; the report counts them as `excluded_registry_swap_rust` and `excluded_registry_swap_perl`, split by rule in `excluded_registry_swap_by_bucket`.

## Defect B: `start_lost` and `start_retained_variant` are both written for one allele

**The defect in one sentence.** Ensembl VEP decides "the start codon is lost" and "the start codon is retained" with two separate tests that never consult each other, so one allele can be given both terms, and on the records of this class it is.

**Why this is a defect.** The two terms contradict each other: a start codon cannot be both lost and kept by the same change. Which of the two is wrong depends on the allele. On a sequence deletion the loss test reads the rebuilt sequence at the old UTR length, finds `TGA` where `ATG` stood, and reports the start altered; `start_lost` is the wrong member. On a structural allele the retention test returns "not altered" for every structural allele before reading anything, and the term is its negation; `start_retained_variant` is the wrong member. In both cases one term was computed from the record and the other was not.

**What vep-rs does instead.** vep-rs writes one of the two terms, the one the record supports: `start_retained_variant` on a sequence deletion whose rebuilt sequence still begins with `ATG`, `start_lost` on a structural deletion whose interval contains the codon. The two are never written together.

**Where it is in Ensembl VEP's source** (`ensembl-variation`, `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`):
- `start_lost`, lines 851-903. The sequence arm (lines 869-883) calls `_inv_start_altered` (lines 906-944), whose `substr($utr_and_translateable, $atg_start, 3)` at line 939 reads at the old UTR length and whose `ne 'ATG'` test is line 941. The structural arm (lines 886-896, at 891 and 894) is an overlap of the span with the codon's three coordinates.
- `start_retained_variant`, lines 947-963, returning `!_ins_del_start_altered(@_)` at line 961 without consulting the loss test.
- `_ins_del_start_altered`, lines 1028-1075, returning 0 for a `TranscriptStructuralVariationAllele` at line 1037 and comparing the rebuilt string's tail with the original coding sequence at line 1071.

**What correcting it would change.** Correcting the misread at line 939 removes the sequence pairs; correcting the structural shortcut at line 1037 removes the structural pairs. With both corrected, classes 3 and 4 vanish and Ensembl VEP writes one term per allele, as vep-rs does.

**The classes it produces.**
- Class 3: the SNP/indel datasets (ClinVar GRCh37 and GRCh38, gnomAD v2.1.1 and v4.1 chr21, 1000 Genomes high-coverage chr21: 261, 928, 3, 77 and 1 pairs); both members of each pair are set aside.
- Class 4: the structural-variant sets, 3 pairs on GRCh37 and 24 on GRCh38; both members of each pair are set aside.

### 3. `start_lost` and `start_retained_variant` emitted together

**In short.** Two deletions touch a start codon, one of each allele kind. ClinVar variation 4629713 deletes one base of a two-base run that continues into the start codon of a RUNX1 transcript, so the codon re-forms one base upstream; the symbolic `<DEL>` `gnomAD-SV_v3_DEL_chr21_88237371` removes the whole start codon of an ADAMTS5 transcript. On both, Ensembl VEP writes `start_lost` and `start_retained_variant` together; vep-rs writes only one of them, `start_retained_variant` on the sequence deletion and `start_lost` on the structural one. The two terms are alternatives by definition, and each record's own sequence or coordinates say which one holds, so vep-rs's rows are the correct ones.

#### 3.1 What a reader needs to know first

Five Sequence Ontology terms appear in this class:

- `start_lost` (SO:0002012): "a codon variant that changes at least one base of the canonical start codon".
- `start_retained_variant` (SO:0002019): "a sequence variant where at least one base in the start codon is changed, but the start remains".
- `frameshift_variant` (SO:0001589): "a sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three".
- `5_prime_UTR_variant` (SO:0001623): "a UTR variant of the 5' UTR".
- `feature_truncation` (SO:0001906): "a sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence".

`start_retained_variant` exists for the case `start_lost` excludes: the two are alternatives for one allele on one transcript, and Ensembl's IMPACT ranking treats them so (HIGH against LOW).

Ensembl VEP decides the two start terms with two independent tests that were written for different variant kinds, and on some alleles both return true.

On a sequence variant (an allele given as bases), Ensembl VEP rebuilds the 5' UTR plus coding sequence with the edit applied. `start_retained_variant` is true when the rebuilt sequence still ends with the original coding sequence, or still reads ATG at the coding start with the UTR unchanged. `start_lost` is true when the three bases read at the *old* UTR length in the edited sequence are not ATG.

When a deletion removes one base of a run that continues into the start codon, the codon re-forms one base upstream. The first test finds the coding sequence intact and says retained; the second reads the edited sequence at the old offset and says lost.

On a structural allele (a symbolic `<DEL>` and the like) the retention test is skipped altogether and returns 0 without reading any sequence. `start_retained_variant` is the negation of that test, so it is true for every structural allele that overlaps a start codon, including a deletion that removes it.

#### 3.2 The record

Two records show the class, one of each allele kind.

**A sequence deletion.** ClinVar variation 4629713, from the ClinVar GRCh37 VCF:

```
21	36265245	4629713	AT	A	.	.	CLNHGVS=NC_000021.8:g.36265247del;CLNVC=Deletion;GENEINFO=RUNX1:861
```

VEP trims the shared A and annotates the deletion of the T at 36265246 (Location `21:36265246`, Allele `-`). The transcript is ENST00000486278 (RUNX1, reverse strand). Its exon 3 of 7 spans 36265222-36265260. The 5' UTR is exons 1 and 2 plus the first 14 bases of exon 3, 36265260 down to 36265247, and the coding sequence begins at 36265246 (cDNA position 217).

On the forward strand the reference at 36265244-36265248 reads `CATTC`. Read on the transcript strand, from 36265248 down to 36265244, it is `GAATG`: the UTR ends `...GA` at 36265248-36265247 and the start codon ATG occupies 36265246, 36265245 and 36265244. The deleted base is the A of that ATG.

Removing it from `...GA|ATGAATCCT...` gives `...G|ATGAATCCT...`: the UTR is one base shorter and the coding sequence, ATG included, is byte for byte what it was. ClinVar's own name for the record, `g.36265247del`, describes the same edit as the deletion of the UTR's last base.

**A structural deletion.** `gnomAD-SV_v3_DEL_chr21_88237371`, from gnomAD-SV v4.1 chr21 (GRCh38, `gnomad_sv_v4.1_chr21.vcf.gz`; the `v3` in the record id is gnomAD's own naming of the call set inside the v4.1 release):

```
21	26965974	gnomAD-SV_v3_DEL_chr21_88237371	N	<DEL>	1	LOWQUAL_WHAM_SR_DEL;OUTLIER_SAMPLE_ENRICHED	END=26966647;SVLEN=673;SVTYPE=DEL
```

The deletion removes 26965975-26966647. The transcript is ENST00000284987 (ADAMTS5, reverse strand, canonical). Its exon 1 of 8 spans 26965288-26967088. The 5' UTR is 26967088 down to 26966392 (697 bases) and the start codon ATG occupies 26966391, 26966390 and 26966389 (forward-strand reference `CAT` at 26966389-26966391).

The deleted interval contains the last 256 UTR bases, the whole start codon and the first 417 coding bases (26966391 down to 26965975). That is what the row's CDS position `?-417` and protein position `?-139` say.

#### 3.3 What each engine writes

The sequence deletion, on ENST00000486278:

```
# Ensembl VEP
4629713	21:36265246	-	ENSG00000159216	ENST00000486278	Transcript	frameshift_variant,start_lost,start_retained_variant	217	1	1	M/X	Atg/tg	-	IMPACT=HIGH;STRAND=-1
# vep-rs
4629713	21:36265246	-	ENSG00000159216	ENST00000486278	Transcript	frameshift_variant,start_retained_variant	217	1	1	M/X	Atg/tg	-	IMPACT=HIGH;STRAND=-1
```

The structural deletion, on ENST00000284987:

```
# Ensembl VEP
gnomAD-SV_v3_DEL_chr21_88237371	21:26965975-26966647	deletion	ENSG00000154736	ENST00000284987	Transcript	frameshift_variant,start_lost,feature_truncation,start_retained_variant,5_prime_UTR_variant	442-1114	?-417	?-139	-	-	-	IMPACT=HIGH;STRAND=-1;OverlapBP=673;OverlapPC=1.37
# vep-rs
gnomAD-SV_v3_DEL_chr21_88237371	21:26965975-26966647	deletion	ENSG00000154736	ENST00000284987	Transcript	frameshift_variant,start_lost,feature_truncation,5_prime_UTR_variant	442-1114	?-417	?-139	-	-	-	IMPACT=HIGH;STRAND=-1;OverlapBP=673;OverlapPC=1.37
```

In each pair the rows agree in every column but Consequence: Ensembl VEP writes both start terms, vep-rs writes one. The second record's other coding transcript, ENST00000970346, carries the same difference.

#### 3.4 Why Ensembl VEP writes that

On the sequence deletion, the retention test rebuilds UTR plus coding sequence with the A removed, finds that the rebuilt string still ends with the original coding sequence, and reports the start not altered, so `start_retained_variant` is emitted. The loss test reads three bases at the old UTR length in the rebuilt string, finds `TGA` where `ATG` stood, and reports the start altered, so `start_lost` is emitted too.

On the structural deletion, `start_lost` is emitted because the deleted span overlaps the codon's three bases. `start_retained_variant` is emitted because the retention test returns 0 for every structural allele before reading anything, and the term is its negation.

Source (`ensembl-variation`):

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, four functions:
  - `start_lost`, lines 851-903. The sequence-variant arm is lines 869-883: line 871 calls `_inv_start_altered`. The structural arm is lines 886-896: lines 891 and 894 test an overlap of the span with the codon's three coordinates.
  - `_inv_start_altered`, lines 906-944: its `substr($utr_and_translateable, $atg_start, 3)` at line 939 reads at the old UTR length, and its `ne 'ATG'` test is line 941.
  - `start_retained_variant`, lines 947-963: returns `!_ins_del_start_altered(@_)` at line 961.
  - `_ins_del_start_altered`, lines 1028-1075: returns 0 for a `TranscriptStructuralVariationAllele` at line 1037 and compares the rebuilt string's tail with the original coding sequence at line 1071.

#### 3.5 Why vep-rs's row is the correct one

The two terms are mutually exclusive by definition: a start that "remains" is not lost.

On the sequence deletion the edited molecule reads `ATGAATCCT...` from one base upstream of where it read before, with the coding sequence identical to the last base. The start remains, so `start_retained_variant` is the term the sequence supports; ClinVar's `g.36265247del` describes the same molecule as an edit that does not touch the codon at all. vep-rs keeps `start_retained_variant` and drops `start_lost`.

On the structural deletion the codon's three bases, 26966389-26966391, lie inside the deleted interval 26965975-26966647. Nothing of the codon is left to be retained, so `start_lost` is the term the coordinates support, and vep-rs drops `start_retained_variant`.

In both cases vep-rs emits the member of the pair that the record's own sequence or coordinates support, and every other term of the row unchanged.

Both rows also carry `frameshift_variant`, which both engines derive from the deletion of one base at the annotated coding coordinates. Ensembl VEP's sequence rule counts deleted bases and does not ask whether the frame re-forms upstream; on the structural row its structural rule fires on the span. That term is identical on both sides, is not part of the class, and is not defended here.

#### 3.6 Try it yourself

The sequence deletion has two VCF spellings of one edit, because the deleted A is part of a two-base run. Annotate both with Ensembl VEP 115.2 (GRCh37 cache and FASTA, the command of the closing section):

- The left-aligned `21 36265245 . AT A` above gives ENST00000486278 `frameshift_variant,start_lost,start_retained_variant`.
- The right-aligned `21 36265246 . TT T` (ClinVar's own `g.36265247del`) gives `5_prime_UTR_variant` with no start term at all; vep-rs writes the same `5_prime_UTR_variant` for that spelling.

One molecule, and Ensembl VEP calls its start codon both lost and untouched depending on which base of the run the record names.

The structural deletion has an explicit spelling too: write the same 673 deleted bases out. REF is the reference sequence of 21:26965974-26966647, 674 bases (`samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 21:26965974-26966647`); ALT is its first base, `G`.

- The explicit spelling: Ensembl VEP evaluates the edited sequence and writes ENST00000284987 `start_lost,5_prime_UTR_variant`, with no `start_retained_variant`. Its sequence-variant rules also drop `feature_truncation`, a structural term, and `frameshift_variant`, because the 417 deleted coding bases are a multiple of three; the probe is about the start pair.
- The `<DEL>` spelling of the same interval receives both start terms.

#### 3.7 How the comparator treats this class

`_is_start_cooccurrence_swap` in `scripts/concordance/compare_vep_outputs.py` matches a pair when Ensembl VEP's set contains both `start_lost` and `start_retained_variant` and vep-rs's set is the same set minus exactly the erroneous member for the allele kind: `start_lost` on a sequence allele, `start_retained_variant` on a structural one. An allele written as bases or `-` is a sequence allele; a Sequence Ontology class name such as `deletion` is structural. A vep-rs set that adds any other term, or drops the other member, is not matched and stays charged. Both members of a matching pair are set aside.

The class is fully masked. Counts, from the SNP/indel summary's `start_cooccurrence_swap` line:

- ClinVar GRCh37: 261 pairs.
- ClinVar GRCh38: 928 pairs.
- gnomAD v2.1.1 chr21: 3 pairs.
- gnomAD v4.1 chr21: 77 pairs.
- 1000 Genomes high-coverage chr21: 1 pair.

Its reach into the structural-variant sets is class 4.

### 4. Class 3 on the structural-variant sets

**In short.** Structural deletions and gains over a start codon in the structural-variant files reach the structural arm of class 3. Ensembl VEP writes `start_lost` and `start_retained_variant` together on them; vep-rs keeps one of the two. For the structural deletions, whose interval contains the codon, the argument of 3.5 holds and vep-rs's row is the correct one. For the five structural gains the record decides only which member came from evidence, and the comparator sets both members aside rather than deciding it.

#### 4.1 What a reader needs to know first

The structural-variant sets are scored by a different script from the SNP/indel datasets, and structural deletions and gains over a start codon reach the structural arm of class 3. Nothing new happens in either engine; the same predicate is applied by the second comparator. The terms are those of 3.1.

#### 4.2 The record

The structural deletion of 3.2, `gnomAD-SV_v3_DEL_chr21_88237371` on ENST00000284987 (gnomAD-SV v4.1 chr21, GRCh38), and the RUNX1 deletion of 3.2, which the ClinVar structural-variant file `clinvar_sv_chr21.vcf.gz` (GRCh37) also carries.

#### 4.3 What each engine writes

The rows are those of 3.3.

#### 4.4 Why Ensembl VEP writes that

As 3.4 for the deletions. On a structural gain, the same two structural tests apply:

- `start_lost` from the overlap of the stated span with the codon;
- `start_retained_variant` from the retention test that returns 0 for every structural allele.

#### 4.5 Why vep-rs's row is the correct one for the deletions, and what the five structural gains show

As 3.5 for the 21 structural deletions of the class, whose interval contains the codon.

The five structural gains of the class, all over a start codon, are:

- an Alu insertion on two transcripts;
- an insertion;
- two duplications.

For them the deletion argument does not apply, and a symbolic gain does not say what becomes of the codon.

What the record does decide is which member came from evidence. Ensembl VEP's structural `start_lost` is the overlap of the stated span with the codon's three coordinates (3.4). Its `start_retained_variant` on a structural allele comes from a test that returns "retained" for every structural allele without reading the sequence. So the retained member is the one emitted from nothing and is the contradiction, and vep-rs keeps the member computed from the record.

That is a statement about which test read the record, not about the biology. Whether a gain over the codon leaves the start in place is not decidable from a symbolic record (a tandem duplication keeps its first copy of the codon). So for these five pairs the page claims only that vep-rs's row drops the term emitted without evidence, and the comparator sets both members aside rather than deciding it.

#### 4.6 Try it yourself

For the deletions, 3.6. The five structural gains are in the files named in 4.7; the `<CN2>` record is `21 37367279 synth_cn2_0016 T <CN2>` with `END=37521417`.

#### 4.7 How the comparator treats this class

The same `filter_registry_swap_pairs` as class 2 applies the class 3 rule to the pairs of the structural-variant comparison.

Counts: the class 3 rule matches 3 pairs on GRCh37 and 24 on GRCh38.

- GRCh37, 3 pairs: one each in `11_multi_allelic.vcf.gz`, `clinvar_sv_chr21.vcf.gz` (where it is the RUNX1 deletion of 3.2) and `gnomad_sv_v2.1_chr21.vcf.gz`.
- GRCh38, 24 pairs: `gnomad_sv_v4.1_chr21.vcf.gz` 20, and one each in `05_symbolic_del_ins.vcf.gz`, `07_cnv_repeat.vcf.gz`, `12_complex_imprecise.vcf.gz` and `1kg_sv_chr21.vcf.gz`; 23 distinct, one deletion recurring in two files.

The tandem-repeat file, `07_cnv_repeat.vcf.gz`, carries 700 `<CNV>` and `<CN0>` to `<CN9>` copy-number records beside its 200 `<CNV:TR>` records, and Ensembl VEP expands only the `<CNV:TR>` ones (12.4). Its class 3 pair is `synth_cn2_0016`, a `<CN2>` copy gain at 21:37367279-37521417 over the start codon of ENST00000647188.

Of the 27 pairs, 21 are structural deletions whose interval contains the codon and 5 are structural gains, where vep-rs keeps `start_lost`. The RUNX1 record is the one sequence allele, where it keeps `start_retained_variant`.

Fully masked; the report counts them as `excluded_registry_swap_rust` and `excluded_registry_swap_perl`, split by rule in `excluded_registry_swap_by_bucket`.

## Defect C: a record marked skipped is still annotated, against whatever its neighbours loaded

**The defect in one sentence.** A structural variant that Ensembl VEP marks skipped, for a span above `--max_sv_size` or an ALT type with no Sequence Ontology term, loads no transcript regions of its own but is still annotated, against whatever regions its batch neighbours loaded.

**Why this is a defect.** Ensembl VEP reports the skip in its warnings file, then annotates the record anyway. The output depends on things that have nothing to do with the record: which other records share its chunk, and so `--fork`, `--buffer_size` and the order of lines in the file. One record gets 51 transcripts, or 497, or 68, or none, depending on its neighbours, while the warnings file states that it was skipped. An annotation that changes with the rest of the file is not an annotation of the record, and a record declared skipped and then annotated in part is a self-contradiction.

**What vep-rs does instead.** vep-rs's output for a record depends on the record's coordinates alone. In the default and tab outputs it annotates an oversize record against every transcript it overlaps, the same set from any file; in the VCF and JSON outputs it does what Ensembl VEP does there, carrying the line without consequences or omitting it. For a `<NON_REF>` reference block, which asserts the absence of a variant, vep-rs writes one `intergenic_variant` row, the row Ensembl VEP itself writes when the block is annotated alone.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):
- `Bio/EnsEMBL/VEP/Parser.pm`, lines 491-498: the size test writes the warning and sets `vep_skip`; the record is dropped only when VEP runs as a REST server. `Config.pm` line 310 holds the 10,000,000 default.
- `Parser/VCF.pm`, lines 477-481 and 575: an ALT type with no Sequence Ontology term gets the same mark.
- `AnnotationSource.pm`, line 238: `next if $vf->{vep_skip}` when regions are collected, before the breakend-mate loop at lines 250-257; line 143, `filter_features_by_min_max`, then drops loaded transcripts outside the chunk's coordinate range.
- `InputBuffer.pm`, lines 284-330: `get_overlapping_vfs` hands each loaded transcript every record of the buffer that overlaps it, skipped or not.
- `Runner.pm`, line 465: a chunk is at most `int(buffer_size / (2 * fork))` records.

**What correcting it would change.** If a skipped record were left out of the annotation pass, Ensembl VEP would write for it what it writes for it alone, the same from any file. Class 8 would then agree with vep-rs. Classes 5 to 7 would differ only in that vep-rs annotates an oversize or unsupported record in full on the default and tab outputs where Ensembl VEP writes nothing or one `intergenic_variant` row (5.5).

**The classes it produces.**
- Class 5: a structural variant above `--max_sv_size`, annotated against the batch's transcripts (gnomAD-SV v2.1 and v4.1 chr21). Set aside: vep-rs's rows, and Ensembl VEP's rows that name a transcript vep-rs does not, each of which is also a class 9 tuple.
- Class 6: the mate end of such a breakend, which Ensembl VEP writes only when a neighbour loaded the region (gnomAD-SV v2.1 and v4.1 chr21). Set aside: vep-rs's rows.
- Class 7: `<CPX>`, an ALT type Ensembl VEP has no term for (gnomAD-SV v4.1 chr21). Set aside: vep-rs's rows.
- Class 8: `<NON_REF>`, a reference block annotated as a variant (the synthetic special-allele set). Set aside: Ensembl VEP's transcript rows and vep-rs's `intergenic_variant` row.

### 5. Structural variants above `--max_sv_size` annotated against whichever transcripts the batch loaded

**In short.** The record is a gnomAD-SV breakend whose span of 27 Mb on chromosome 21 exceeds the `--max_sv_size` cap (default 10,000,000). Ensembl VEP's warnings file declares the record skipped, yet its output carries 51 transcript rows for it, a count that changes with `--fork`, `--buffer_size` and the neighbouring lines of the file. vep-rs writes a row for each of the 3,763 transcripts the span overlaps or comes within 5,000 bases of, the same set from any file that contains the record. The record's coordinates alone determine which transcripts it reaches, so vep-rs's complete set is the correct one.

#### 5.1 What a reader needs to know first

A breakend record (`SVTYPE=BND`) describes one end of a rearrangement: the position where the reference is broken and the position it is joined to, its mate. gnomAD-SV writes the record with a symbolic `<BND>` allele and puts the mate in `INFO/CHR2` and `INFO/END2`; VCF's own notation writes the mate into the allele, `N[chr21:41770788[`. Ensembl VEP reads either form.

When the mate is on the same chromosome, Ensembl VEP represents the record as a feature spanning from the base after POS to the mate (Location `21:14604132-41770788`). It can then write rows under two alleles: the local end, `N.`, and the mate, `N[chr21:41770788[`. For the record below the reference output has rows under the local allele only; class 6 explains the missing mate allele.

A breakend record with a parsed mate never gets an `intergenic_variant` row. Ensembl VEP's intergenic object creates its alleles through the same distance gate as a transcript row (9.4), and with no transcript to measure against it creates none.

Such a record annotated against nothing therefore writes nothing at all (`ensembl-variation`, `StructuralVariationOverlap.pm` lines 72-87 with an undefined feature). A record without a parsed mate (a `<DEL>`, `<CPX>` or `<NON_REF>`) takes the unconditional allele of lines 64-71 and writes one `intergenic_variant` row instead. This is why 5.6 below shows 0 rows for the record alone and 7.6 shows 1.

`--max_sv_size` (default 10,000,000) is the documented size cap. A structural variant whose span exceeds it is reported as skipped in the warnings file that Ensembl VEP writes beside its output (`<output>_warnings.txt`). The warning for the record below reads `WARNING: line 9777 skipped (21 14604131 gnomAD-SV_v3_BND_chr21_3f3f16e7 N ...): variant size (27166656) is bigger than --max_sv_size (10000000)`.

Ensembl VEP annotates its input in batches of `--buffer_size` records (default 5,000) and splits each batch into chunks for the `--fork` workers. For each chunk it loads the transcript-cache regions (1 Mb bins) that the chunk's records overlap.

Four Sequence Ontology terms are involved in this class:

- `feature_truncation` (SO:0001906): "a sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence".
- `intron_variant` (SO:0001627): "a transcript variant occurring within an intron".
- `non_coding_transcript_exon_variant` (SO:0001792): "a sequence variant that changes non-coding exon sequence in a non-coding transcript".
- `coding_transcript_variant` (SO:0001968): "a transcript variant of a protein coding gene".

#### 5.2 The record

`gnomAD-SV_v3_BND_chr21_3f3f16e7`, gnomAD-SV v4.1 chr21 (GRCh38, `gnomad_sv_v4.1_chr21.vcf.gz`), a breakend at 21:14604131 whose mate lies 27 Mb away on the same chromosome:

```
21	14604131	gnomAD-SV_v3_BND_chr21_3f3f16e7	N	<BND>	251	UNRESOLVED	CHR2=chr21;END=14604131;END2=41770788;POS2=41770787;SVLEN=27166657;SVTYPE=BND
```

Two transcripts stand for the two kinds of row. ENST00000400562 (SAMSN1, reverse strand, a retained-intron transcript, 14593720-14658821) has the local breakpoint 14604132 in its intron 5 (14602100-14609481), and the feature's span covers its exons 1 to 5. Ensembl VEP names it. ENST00000270142 (SOD1, forward strand, 31659693-31668931) lies 17 Mb inside the span. Ensembl VEP does not name it.

Between 14604132 and 41770788 the span overlaps, or comes within 5,000 bases of, 3,763 transcripts of the GRCh38 chromosome 21 set.

#### 5.3 What each engine writes

Ensembl VEP's reference output (the whole file, `--fork 4 --buffer_size 5000`) carries 51 rows for the record, all under the local allele `N.`. vep-rs writes 3,783 rows: 3,763 under `N.` and 20 under the mate allele (class 6). All 51 of Ensembl VEP's rows are among them.

The SAMSN1 row is identical in both. The SOD1 row, one of the 3,712 local-allele rows vep-rs writes and Ensembl VEP does not, is the divergence:

```
# Ensembl VEP (1 of 51 rows)
gnomAD-SV_v3_BND_chr21_3f3f16e7	21:14604132-41770788	N.	ENSG00000155307	ENST00000400562	Transcript	feature_truncation,non_coding_transcript_exon_variant,intron_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1
# vep-rs (the same row, 1 of 3,783)
gnomAD-SV_v3_BND_chr21_3f3f16e7	21:14604132-41770788	N.	ENSG00000155307	ENST00000400562	Transcript	feature_truncation,non_coding_transcript_exon_variant,intron_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1
# Ensembl VEP
(no row for ENST00000270142)
# vep-rs (1 of the 3,712 rows Ensembl VEP lacks)
gnomAD-SV_v3_BND_chr21_3f3f16e7	21:14604132-41770788	N.	ENSG00000142168	ENST00000270142	Transcript	feature_truncation,coding_transcript_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
```

Every row Ensembl VEP writes for the record, vep-rs writes too. The 3,712 rows like the SOD1 row are the ones only vep-rs writes.

#### 5.4 Why Ensembl VEP writes that

The size check does not remove the record. It marks it skipped and writes the warning, but the record stays in its batch. The mark only stops the record from loading any transcript-cache region of its own. The annotation pass then runs it against whatever regions the other records of its chunk happened to load, further limited to the coordinate range those other records span.

Which transcripts are reported therefore depends on which records share the chunk, and so on `--fork`, `--buffer_size` and the order and content of the file.

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/Parser.pm` lines 491-498: the size test calls `skipped_variant_msg` and sets `vep_skip`; the record is dropped only when VEP runs as a REST server.
- `Config.pm` line 310: the 10,000,000 default.
- `AnnotationSource.pm` line 238: `next if $vf->{vep_skip}` in `get_all_regions_by_InputBuffer`, before any region is requested for the record; line 143: `filter_features_by_min_max`, which drops loaded transcripts outside the chunk's coordinate range.
- `InputBuffer.pm` lines 284-330: `get_overlapping_vfs`, which hands each loaded transcript every record of the buffer that overlaps it, skipped or not.
- `Runner.pm` line 465: a chunk is at most `int(buffer_size / (2 * fork))` records.

#### 5.5 Why vep-rs's row is the correct one

The record's coordinates alone determine which transcripts it reaches: every transcript the span 14604132-41770788 overlaps or comes within 5,000 bases of on chromosome 21, 3,763 of them. vep-rs reports that set, the same set from any file that contains the record. Ensembl VEP reports 51 of them, or 497, or 68, or none, depending on the neighbouring lines and the fork count (5.6), while its own warnings file says the record was skipped.

An annotation that changes with the rest of the file is not an annotation of the record. A record that is declared skipped and then annotated in part is a self-contradiction. Both point at the same defect.

vep-rs takes `--max_sv_size` too, with the same default of 10,000,000. In the default and tab outputs it annotates an oversize record against every transcript it overlaps. In the VCF output it carries the record's line without consequences, and in the JSON output it omits the record, which is what Ensembl VEP does in those two formats.

Ensembl VEP's cap is a memory guard (its source comments "to avoid memory problems"), not a statement about the record. Skipping is therefore a legitimate output, and vep-rs's VCF and JSON outputs give it. What neither engine should write is a partial set that depends on the neighbouring lines; on the default and tab outputs vep-rs writes the complete set instead.

#### 5.6 Try it yourself

Run Ensembl VEP 115.2 (GRCh38 cache and FASTA, the command of the closing section) on inputs built from the record:

| Input                                                                  | Rows for the record (local allele `N.` / mate allele)                                 |
| ---------------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| the record alone                                                       | 0 / 0 (the warnings file says skipped; a breakend record gets no intergenic row, 5.1) |
| the record plus the SNV `21 14604131 . A C`                            | 4 / 0                                                                                 |
| the record plus the SNV `21 41770788 . C A`                            | 20 / 20                                                                               |
| the whole `gnomad_sv_v4.1_chr21.vcf.gz`, `--fork 4 --buffer_size 5000` | 51 / 0                                                                                |
| the same file, `--fork 1`                                              | 497 / 0                                                                               |
| the same file, `--fork 16`                                             | 68 / 0                                                                                |
| the same file, `--buffer_size 1`                                       | 0 / 0                                                                                 |

vep-rs writes 3,763 / 20 for every one of those inputs. Ensembl VEP's count moves with the neighbouring lines and the fork count, from 0 to 497; vep-rs's does not.

#### 5.7 How the comparator treats this class

`filter_intended_divergences` in `scripts/validation/compare_sv_concordance.py` acts on a (Location, Allele) key that both engines annotated and whose span exceeds 10,000,000 bases (or whose allele type is one Ensembl VEP has no Sequence Ontology term for, class 7). A transcript named by vep-rs and not by Ensembl VEP is set aside on the vep-rs side when vep-rs named more transcripts for the key than Ensembl VEP did. A transcript named by Ensembl VEP and not by vep-rs is set aside on the Ensembl VEP side when vep-rs named at least as many. The direction guards (those two conditions on how many transcripts each engine named) keep a vep-rs under-annotation charged.

The rule needs Ensembl VEP to have written at least one row for the key. Every same-chromosome record above the cap in the two files has one (34 records in gnomAD-SV v2.1 chr21, 32 in v4.1), so none is outside the rule. The 18 GRCh38 breakends whose mate is on chromosome X are point records (11.1), have no measured span, and are not in this class.

Counts:

- gnomAD-SV v2.1 chr21 (GRCh37): 10,934 vep-rs tuples and no Ensembl VEP tuple.
- gnomAD-SV v4.1 chr21 (GRCh38): 183,964 vep-rs tuples (218 of them class 7) and 1,823 Ensembl VEP tuples. Each of those Ensembl VEP tuples names a transcript on another chromosome and is therefore also a class 9 tuple (a tuple in two classes is set aside once).

The class is fully masked on both sides (`excluded_transcript_selection_rust` and `excluded_transcript_selection_perl` in the report).

Its blind spot is stated here because it is inherent to the shape: a wrong vep-rs term on a transcript Ensembl VEP's batch never loaded would be set aside unseen. The direction guards limit that exposure to transcripts Ensembl VEP did not name at all. The mate allele of the same records is class 6.

The SNP/indel comparator, `scripts/concordance/compare_vep_outputs.py`, applies the same rule under the same two guards to a whole-genome input, which carries such records among its sequence variants. Its sorted merge collects the transcripts each engine named on every record in the rule's scope, and `filter_snp_indel_intended_divergences` sets aside the one-sided rows the rule admits, reported as `sv_transcript_selection` beside the pair classes. The scope test (`_batch_dependent_annotation_scope`) is defined once, in the SNP/indel comparator, and imported by the structural-variant one.

### 6. The mate end of a breakend above `--max_sv_size`

**In short.** The record is the same-chromosome breakend of class 5, whose span 14604132-41770788 exceeds `--max_sv_size`; its output has two alleles, the local end and the mate. Ensembl VEP writes no row under the mate allele. vep-rs writes 20, one for each transcript within 5,000 bases of the mate breakpoint 41770788. Ensembl VEP's own admission rule gives the mate allele exactly those 20 rows, and Ensembl VEP writes them itself once a neighbouring record loads the region, so vep-rs's rows are the correct ones.

#### 6.1 What a reader needs to know first

The record of class 5 has two alleles in the output: the local end, written `N.`, and the mate, written `N[chr21:41770788[`. Ensembl VEP creates a mate-allele row for a transcript only when the mate lies within 5,000 bases of that transcript on the same chromosome (the breakend distance gate of 9.4). So the mate allele's rows are the transcripts around 41770788. The terms are those of 5.1.

#### 6.2 The record

The same record as class 5, `gnomAD-SV_v3_BND_chr21_3f3f16e7` (gnomAD-SV v4.1 chr21, GRCh38).

The mate breakpoint 41770788 lies in intron 2 (41768018-41774115) of ENST00000700893, a long non-coding RNA (lncRNA) transcript of ENSG00000289513 on the forward strand at 41767163-41774864. It also lies within 5,000 bases of 19 other transcripts, the two RIPK4 transcripts among them.

Ensembl VEP represents this same-chromosome record as the span 14604132-41770788 (5.1), and both engines evaluate the mate allele of a span record over that span, so the row's terms are the span's. Exons 1 and 2 of ENST00000700893 (41767163-41767234 and 41767859-41768017) lie inside it, hence `non_coding_transcript_exon_variant` beside `intron_variant` below. The point-breakend rule of class 11 does not apply to a span record.

#### 6.3 What each engine writes

Ensembl VEP's reference output has no mate-allele row for the record. vep-rs writes 20, one of them:

```
# Ensembl VEP
(no row under the mate allele)
# vep-rs (1 of 20 mate-allele rows)
gnomAD-SV_v3_BND_chr21_3f3f16e7	21:14604132-41770788	N[chr21:41770788[	ENSG00000289513	ENST00000700893	Transcript	feature_truncation,non_coding_transcript_exon_variant,intron_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
```

The 20 rows are the 20 transcripts of 6.2: ENST00000700893 and the 19 others within 5,000 bases of the mate breakpoint.

#### 6.4 Why Ensembl VEP writes that

The same skip as class 5, one step earlier in the loop. The region-loading pass skips a `vep_skip` record before it reaches the code that requests the regions around a breakend's mate. So no region near 41770788 is loaded for this record, and the mate allele's transcripts exist for the annotation pass only if another record of the chunk loaded them. In the reference run none did.

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/AnnotationSource.pm`, line 238: `next if $vf->{vep_skip}`, the skip, against lines 250-257, the breakend-mate region loop it skips.

Source (`ensembl-variation`):

- `Bio/EnsEMBL/Variation/StructuralVariationOverlap.pm`, lines 72-87 and 130-146: the mate allele is created for a transcript within `MAX_DISTANCE_FROM_TRANSCRIPT`.
- `Utils/VariationEffect.pm`, line 60: `MAX_DISTANCE_FROM_TRANSCRIPT` is 5,000 bases.

#### 6.5 Why vep-rs's row is the correct one

The mate breakpoint lies inside ENST00000700893, in its intron 2, and within 5,000 bases of the other 19 transcripts. By Ensembl VEP's own admission rule, then, the mate allele has rows on exactly those 20 transcripts. Whether they are written cannot depend on a neighbouring line loading the region. When that region is loaded, Ensembl VEP writes the same 20 rows with the same consequence sets (6.6).

The reservation of 5.5 applies here too. The documented outcome for an oversize record is a skip; vep-rs gives that skip in its VCF and JSON outputs, and on the default and tab outputs it writes the complete set rather than the neighbour-dependent one.

#### 6.6 Try it yourself

As 5.6: annotate the record alone and beside one single-nucleotide variant (SNV), and count Ensembl VEP's mate-allele rows.

- The record plus the SNV `21 41770788 . C A`, at the mate breakpoint: Ensembl VEP writes 20 mate-allele rows, and they equal vep-rs's 20 in Feature and Consequence, row for row.
- The record alone: none.
- The record plus the SNV `21 14604131 . A C` at the local end: none.
- The whole file at every `--fork` and `--buffer_size` setting of 5.6: none.

The same procedure on each of the 37 records of the class in the two gnomAD-SV files (the record plus one SNV at its `END2` position) reproduces vep-rs's mate-allele rows on 366 of 366.

#### 6.7 How the comparator treats this class

`filter_giant_breakend_mate_divergences` in `scripts/validation/compare_sv_concordance.py` sets aside a vep-rs tuple when all four of these hold:

- its Allele is bracket notation naming a mate on the same chromosome as the Feature;
- the Location span exceeds 10,000,000 bases;
- Ensembl VEP wrote at least one tuple for the Location (it saw the record) and none under the mate allele;
- the mate lies within 5,000 bases of the transcript's span in the cache, so that Ensembl VEP's own gate would have created the allele.

A mate farther from the transcript than that is a row Ensembl VEP never writes and stays charged to vep-rs.

Counts:

- gnomAD-SV v2.1 chr21 (GRCh37): 56 tuples on 15 records.
- gnomAD-SV v4.1 chr21 (GRCh38): 310 tuples on 22 records.

The spans run from 10.7 to 39.2 Mb. Fully masked, vep-rs side, under the counter `excluded_transcript_selection_mate_rust`.

### 7. `<CPX>`, class 5 reaching an unsupported allele type

**In short.** `gnomAD-SV_v3_CPX_chr21_df09f6c2` is a gnomAD-SV v4.1 chr21 record whose symbolic allele `<CPX>` marks a complex rearrangement, a structural variant type Ensembl VEP does not support. Ensembl VEP writes 412 rows for the record; vep-rs writes 508, the same 412 plus rows for 96 more transcripts, all in the last 200 kb of the span. The record states the span, every one of the 508 transcripts lies inside it, and Ensembl VEP itself names all 508 once the regions are loaded, so vep-rs's rows are the correct ones.

#### 7.1 What a reader needs to know first

`<CPX>` is gnomAD-SV's symbolic allele for a complex rearrangement. Here the rearrangement is a deletion followed by an inversion, `CPX_TYPE=delINV`.

Ensembl VEP has no Sequence Ontology class for it and says so in the warnings file. For a file holding the record alone the warning reads: `WARNING: line 1 skipped (21 14923958 gnomAD-SV_v3_CPX_chr21_df09f6c2 N ...): CPX is not a supported structural variant type`.

The Allele column then carries the type string `CPX`, and the record is annotated for the span alone, without a direction term (`feature_elongation` or `feature_truncation`, the term that says what a longer or shorter feature results): the rows claim overlap, not disruption. The terms written are:

- `coding_transcript_variant` (SO:0001968) and `non_coding_transcript_variant` (SO:0001619, "a transcript variant of a non coding RNA gene") for transcripts inside the span;
- `intron_variant` and the exon terms for transcripts the span partly overlaps;
- `mature_miRNA_variant` for a microRNA it covers;
- `intergenic_variant` when nothing is loaded.

#### 7.2 The record

`gnomAD-SV_v3_CPX_chr21_df09f6c2`, from gnomAD-SV v4.1 chr21 (GRCh38, `gnomad_sv_v4.1_chr21.vcf.gz`), is a 1.66 Mb complex event, well below `--max_sv_size`:

```
21	14923958	gnomAD-SV_v3_CPX_chr21_df09f6c2	N	<CPX>	1	UNRESOLVED	CHR2=chr21;CPX_INTERVALS=DEL_chr21:14923958-15004545,INV_chr21:15004545-16588073;CPX_TYPE=delINV;END=16588080;SVLEN=1664122;SVTYPE=CPX
```

Two transcripts stand for the two kinds of row:

- ENST00000285679 (USP25, forward strand, 15730025-15878911) lies entirely inside the span 14923959-16588080, and Ensembl VEP names it.
- ENST00000362160 (MIRLET7C, a microRNA at 16539828-16539911) also lies inside the span, and Ensembl VEP does not name it.

#### 7.3 What each engine writes

Ensembl VEP's reference output carries 412 rows for the record; vep-rs writes 508. The 412 are a subset of the 508. The 96 transcripts only vep-rs names all start at or beyond 16,391,614, in the last 200 kb of the span; 87 of them are transcripts of one lncRNA gene, ENSG00000215386.

```
# Ensembl VEP (1 of 412 rows)
gnomAD-SV_v3_CPX_chr21_df09f6c2	21:14923959-16588080	CPX	ENSG00000155313	ENST00000285679	Transcript	coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;OverlapBP=148887;OverlapPC=100.00
# vep-rs (the same row, 1 of 508)
gnomAD-SV_v3_CPX_chr21_df09f6c2	21:14923959-16588080	CPX	ENSG00000155313	ENST00000285679	Transcript	coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;OverlapBP=148887;OverlapPC=100.00
# Ensembl VEP
(no row for ENST00000362160)
# vep-rs (1 of the 96 rows Ensembl VEP lacks)
gnomAD-SV_v3_CPX_chr21_df09f6c2	21:14923959-16588080	CPX	ENSG00000199030	ENST00000362160	Transcript	mature_miRNA_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;OverlapBP=84;OverlapPC=100.00
```

For ENST00000285679 the two rows are identical in every column. For ENST00000362160 vep-rs writes a row and Ensembl VEP writes none.

#### 7.4 Why Ensembl VEP writes that

An ALT whose type has no Sequence Ontology term is kept with the same skip mark that the `--max_sv_size` size test sets, and everything of class 5 follows from it. No region is loaded for the record, and the annotation pass runs it against whatever its chunk loaded, bounded by the chunk's coordinate range.

Here the coordinate range of the chunk's other records ended between 16,357,273 and 16,391,614: the highest transcript start Ensembl VEP named and the lowest it omitted. The 96 transcripts starting beyond it were therefore dropped, although their 1 Mb bin was loaded.

Alone in a batch the record gets one `intergenic_variant` row rather than none, because it has no parsed mate and so takes the unconditional allele (5.1).

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 477-481 (`get_SO_term` returns nothing for the type, `$skip_line` is set and the type string stands in for the term) and line 575 (`vep_skip` is set on the feature).
- Then `AnnotationSource.pm`, lines 238 and 143, as in 5.4.

#### 7.5 Why vep-rs's row is the correct one

The event's span is stated by the record, and every transcript inside 14923959-16588080 lies inside the interval the record describes. The rows claim overlap, not disruption: `coding_transcript_variant` and its non-coding counterparts are MODIFIER terms. The 96 transcripts Ensembl VEP omits lie inside the interval as surely as the 412 it names.

The boundary between the two groups, 16,391,614, is a property of the other records in the chunk, not of this one. When the regions are loaded Ensembl VEP names the same 508 transcripts with the same consequence sets (7.6).

#### 7.6 Try it yourself

Run Ensembl VEP 115.2 with the GRCh38 cache and FASTA, using the command of the closing section:

| Input | Rows for the record |
| --- | --- |
| the record alone | 1, `intergenic_variant` (the warnings file says the type is unsupported) |
| the record plus the four SNVs `21 14923958 . T A`, `21 15500000 . G A`, `21 16400000 . G A` and `21 16590000 . A C` (the last one past the record's end, so that the chunk's coordinate range reaches the whole span) | 508, equal to vep-rs's 508 in Feature and Consequence |
| the whole `gnomad_sv_v4.1_chr21.vcf.gz` at `--fork 1`, `--fork 4` or `--fork 16` | 412 |
| the same file, `--buffer_size 1` | 1, `intergenic_variant` |

Alone, or at `--buffer_size 1`, the record gets one `intergenic_variant` row. With four SNVs extending the chunk's coordinate range across the whole span, Ensembl VEP writes the 508 rows vep-rs writes.

At `--buffer_size 1` all 169 `CPX` records of the GRCh38 file and all 60 of the GRCh37 file write exactly `intergenic_variant`. At the reference settings 143 and 42 of them carry transcript rows.

The other two giant `CPX` records of the GRCh38 file, `21:26860425-29034168` and `21:36134639-36403117`, carry:

- 112 and 79 rows at `--fork 4`;
- 265 and 103 at `--fork 1`;
- 112 and 70 at `--fork 16`.

#### 7.7 How the comparator treats this class

This class uses the same `filter_intended_divergences` as class 5, whose scope test admits a (Location, Allele) key for either of two reasons:

- the span exceeds 10,000,000 bases;
- the Allele is one of the unsupported symbolic types present in the datasets, `CPX` and `CTX` (`_VEP_SKIP_UNSUPPORTED_ALLELES`).

Both direction guards are unchanged.

The rule reaches 218 vep-rs tuples on the three giant `CPX` records of gnomAD-SV v4.1 chr21 (GRCh38) and nothing else:

- The two `CTX` records are concordant.
- A pair on a `CPX` record where both engines name the transcript and differ on a term is outside the scope test and stays charged.
- The 26 GRCh38 and 18 GRCh37 `CPX` records whose reference rows are `intergenic_variant` alone have no transcript within 5,000 bases, so vep-rs writes the same row for each and they are concordant.

The class is fully masked on the vep-rs side and counted inside `excluded_transcript_selection_rust`.

### 8. `<NON_REF>` is not a variant

**In short.** A gVCF reference-confidence block, a record whose ALT is `<NON_REF>`, states that the sample matches the reference over an interval: it asserts the absence of a variant. Ensembl VEP writes transcript rows for the block when another record in the batch has loaded the region, and `intergenic_variant` when nothing is loaded. vep-rs 0.3.0 annotates the block against no transcript and writes one `intergenic_variant` row for it under every setting. A record that changes no base can affect no transcript, so vep-rs's row is the correct one.

#### 8.1 What a reader needs to know first

A gVCF (genomic VCF) reports every position of a sample, not only the variant ones. Positions where the sample matches the reference are compressed into reference-confidence blocks. A block is one record whose ALT is the symbolic `<NON_REF>` (or its VCF 4.3 spelling `<*>`, the same construct), spanning POS to `INFO/END`, meaning "no variant was called in this interval, with this confidence".

Such a record asserts the absence of a variant; there is no allele to predict a consequence for.

Ensembl VEP has no Sequence Ontology class for the allele and says so in the warnings file (`NON_REF is not a supported structural variant type`). Yet it writes rows for the record:

- `intergenic_variant` when nothing is loaded;
- transcript rows for the span (`intron_variant`, `coding_sequence_variant`, `upstream_gene_variant`, and so on) when a neighbour has loaded the region.

vep-rs 0.3.0 annotates the block against no transcript. It writes the record's one output line in Ensembl VEP's own convention for a record that overlaps nothing it annotated: Feature `-`, consequence `intergenic_variant` (SO:0001628), IMPACT MODIFIER.

#### 8.2 The record

This class is shown on two records. Both are from the synthetic special-allele set `tests/sv_validation/grch37/10_special_alleles.vcf.gz` (GRCh37; the GRCh38 set is `tests/sv_validation/grch38/10_special_alleles.vcf.gz`).

The first is a 279-base block over a coding exon, and one of the 180 GRCh37 blocks that carry transcript rows in the reference output:

```
21	27347258	synth_nonref_0033	C	<NON_REF>	.	PASS	END=27347537
```

The block 27347259-27347537 lies in APP (ENST00000346798, reverse strand, canonical, GRCh37 27252861-27542972). It covers 155 of the 159 bases of exon 11 (27347383-27347541, a coding exon) and 124 bases of intron 11 (27328070-27347382). Eleven other APP transcripts overlap it too.

The second is a 228-base block in an intron, one of the 20 GRCh37 blocks whose reference row is `intergenic_variant`. It is the record on which the flag experiment of 8.6 is run:

```
21	38589780	synth_nonref_0188	A	<NON_REF>	.	PASS	END=38590008
```

The block 38589781-38590008 lies in or near eight transcripts in all:

- it lies in intron 2 (38581058-38592538) of ENST00000454482 (DSCR9, forward strand, lincRNA, 38580804-38594037);
- it lies inside an intron of four other DSCR9 transcripts;
- it is 3,712 bases upstream of ENST00000440629;
- it is within 5,000 bases of two more transcripts.

#### 8.3 What each engine writes

For the first record the reference output (the whole file, `--fork 4 --buffer_size 5000`) carries 12 transcript rows, and vep-rs one row:

```
# Ensembl VEP (1 of 12 rows)
synth_nonref_0033	21:27347259-27347537	<NON_REF>	ENSG00000142192	ENST00000346798	Transcript	coding_sequence_variant,intron_variant	1338-?	1304-?	435-?	-	-	-	IMPACT=MODIFIER;STRAND=-1;OverlapBP=279;OverlapPC=0.10
# vep-rs
synth_nonref_0033	21:27347259-27347537	<NON_REF>	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

For the second record the reference output carries one row, identical to vep-rs's. Add one SNV beside the record and Ensembl VEP writes eight transcript rows for it instead:

```
# Ensembl VEP, reference output
synth_nonref_0188	21:38589781-38590008	<NON_REF>	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
# Ensembl VEP, the record plus the SNV 21 38589780 . A C (1 of 8 rows)
synth_nonref_0188	21:38589781-38590008	<NON_REF>	ENSG00000230366	ENST00000454482	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;OverlapBP=228;OverlapPC=1.72
# vep-rs
synth_nonref_0188	21:38589781-38590008	<NON_REF>	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

vep-rs also warns once per record, on standard error: `21:27347258 <NON_REF>: NON_REF is not a supported structural variant type; a reference-confidence block is written as one intergenic row`.

#### 8.4 Why Ensembl VEP writes that

The unsupported type is kept with the skip mark of class 7, and class 5's batch dependence follows. Alone, the record loads nothing and is written as `intergenic_variant` (the unconditional allele of 5.1, since the record has no mate). Beside a record that loads the region, it is annotated as a span over every transcript there, coding sequence included.

In the reference output:

- GRCh37 set: 180 of the 200 `<NON_REF>` records carry transcript rows (1,838 of them) and 20 carry `intergenic_variant`.
- GRCh38 set: all 200 carry transcript rows (5,403).

At `--buffer_size 1` all 200 of each set write `intergenic_variant`.

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 477-481 and 575: a type with no Sequence Ontology term is marked skipped, as in 7.4.
- `AnnotationSource.pm`, line 238: the skip mark stops the record's own region loading, as in 5.4.

#### 8.5 Why vep-rs's row is the correct one

A reference-confidence block states that the sample carries the reference over the interval. It changes no base, so no transcript's sequence, splicing or dosage is affected. Every transcript-level consequence Ensembl VEP writes for it describes the effect of a variant that the record says does not exist: a `coding_sequence_variant` of APP over 155 coding bases, an `upstream_gene_variant` 3,712 bases from a gene.

vep-rs 0.3.0 therefore annotates the block against no transcript, and writes the record's one output line as Ensembl VEP's fallback row, Feature `-` and `intergenic_variant`. That row is a convention, not a claim that the block lies between genes: the first block lies in a coding exon of APP, the second in a DSCR9 intron. It is the row Ensembl VEP itself writes for each of the 200 blocks at `--buffer_size 1`, and it keeps one output line per input record, which is what the default output promises.

vep-rs 0.2.0 annotated the block as a copy-number span over every transcript it covers. Those rows matched Ensembl VEP's transcript rows whenever Ensembl VEP's batch had loaded the region and differed otherwise, and both readings of a non-variant were wrong.

The warning names the record so that a caller who wants no line at all can drop it.

#### 8.6 Try it yourself

Ensembl VEP 115.2, GRCh37 cache and FASTA, the command of the closing section:

- `synth_nonref_0188` alone gives one row, `intergenic_variant`.
- The record plus the SNV `21 38589780 . A C` gives eight transcript rows for the block, the eight transcripts within 5,000 bases of it.
- The whole `10_special_alleles.vcf.gz` (GRCh37) at `--buffer_size 1` gives `intergenic_variant` for all 200 `<NON_REF>` records.
- The same file at `--fork 4 --buffer_size 5000`: 180 of them carry transcript rows, `synth_nonref_0033` with its 12 rows among them.
- The set of rows changes again between `--fork 4` and `--fork 1`.

vep-rs 0.3.0 writes one `intergenic_variant` row per block under every setting.

#### 8.7 How the comparator treats this class

`filter_non_ref_batch_divergences` in `scripts/validation/compare_sv_concordance.py` acts on a (Location, Allele) key whose Allele is `<NON_REF>` and whose vep-rs set is exactly the one `intergenic_variant` row. At such a key:

- every Ensembl VEP transcript tuple is set aside;
- the vep-rs intergenic row is set aside when Ensembl VEP wrote transcript tuples and no intergenic row there, so that neither engine is charged for a block neither should annotate.

A record on which vep-rs wrote anything else is never touched. A key where Ensembl VEP also wrote the intergenic row keeps it on both sides, matched.

Counts (`excluded_non_ref_batch_perl` and `excluded_non_ref_batch_rust`):

- GRCh37 set: 1,838 Ensembl VEP tuples and 180 vep-rs rows.
- GRCh38 set: 5,403 Ensembl VEP tuples and 200 vep-rs rows.

Fully masked, both sides.

Those transcript rows agreed with vep-rs 0.2.0's rows whenever Ensembl VEP's batch had loaded the region. Raw F1 on the two structural-variant sets is therefore lower for vep-rs 0.3.0 than for 0.2.0 by these tuples together with the pairs of class 11, and by nothing else, since raw F1 sets nothing aside.

`<*>` is the same gVCF construct written the other way, and the argument of 8.5 applies to it word for word; it is out of this class by measurement, not by principle. Ensembl VEP does not mark `<*>` skipped, both engines annotate it as a span, and their rows are identical (1,007 rows on the GRCh37 set, 2,595 on the GRCh38 set), so there is no divergence to set aside. The class covers the spelling on which the two engines differ.

## Defect D: transcripts matched to records by coordinate alone, without the chromosome

**The defect in one sentence.** Ensembl VEP matches loaded transcripts to the records of a batch through one interval tree keyed on coordinate numbers only, with no chromosome, so a chromosome 21 deletion can be handed a chromosome X transcript whose coordinates overlap its numbers.

**Why this is a defect.** A transcript on another chromosome cannot be affected by a deletion on chromosome 21; the row Ensembl VEP writes for it describes an overlap that does not exist. The row appears only when another record in the batch, a breakend whose mate is on chromosome X, loaded that chromosome's regions, so the output again depends on the record's neighbours. And because the record's transcript list is then not empty, the `intergenic_variant` row it should carry is never written.

**What vep-rs does instead.** vep-rs matches transcripts by chromosome and coordinate, so a record on chromosome 21 reaches chromosome 21 transcripts only, and a record with no transcript within range gets its `intergenic_variant` row.

**Where it is in Ensembl VEP's source**:
- `ensembl-vep`, `Bio/EnsEMBL/VEP/InputBuffer.pm`, lines 345-366: `interval_tree` inserts every record and every breakend mate by start and end, with no chromosome; lines 284-330: `get_overlapping_vfs` fetches by coordinates alone.
- `ensembl-vep`, `AnnotationSource.pm`, lines 250-257: the regions around each breakend's mate are loaded for the batch.
- `ensembl-variation`, `Bio/EnsEMBL/Variation/StructuralVariationOverlap.pm`, lines 64-71: a record without a parsed mate (a `<DEL>`, `<INS>`, `<CPX>`, `<NON_REF>` or a single breakend) gets its one allele unconditionally and is never chromosome-checked. Lines 72-87 and 130-146 are the contrast: a record with a parsed mate has its ends passed through `_close_to_feature`, which compares `seq_region_name`s and applies the 5,000-base window.
- `ensembl-variation`, `StructuralVariationFeature.pm`, `get_IntergenicStructuralVariation`, lines 674-689: the `intergenic_variant` row is created only when the transcript list is empty.

**What correcting it would change.** If the tree were keyed on chromosome and coordinate, or the unconditional allele were chromosome-checked, classes 9 and 10 would vanish: the off-chromosome rows would not be written and the record's `intergenic_variant` row would be.

**The classes it produces.**
- Class 9: the off-chromosome transcript rows Ensembl VEP writes (gnomAD-SV v2.1 and v4.1 chr21 and the synthetic breakend set: 1,838 and 22,481 tuples). Set aside: Ensembl VEP's rows.
- Class 10: the record's own `intergenic_variant` row, which those rows displace (155 and 518 rows). Set aside: vep-rs's rows.

### 9. Transcripts from another chromosome

**In short.** A 75-base deletion on chromosome 21 from gnomAD-SV v2.1, `gnomAD-SV_v2.1_DEL_21_178762`, has coordinates that fall, as numbers only, inside an intron of GRPR, a gene on chromosome X. Ensembl VEP writes `intron_variant` on that chromosome X transcript, ENST00000380289; vep-rs writes `intergenic_variant`, because the deletion has no chromosome 21 transcript within the 5,000-base window. A deletion on chromosome 21 changes no base of chromosome X, so vep-rs's row is the correct one, and Ensembl VEP itself writes `intergenic_variant` when the record is annotated alone.

#### 9.1 What a reader needs to know first

A transcript can only be affected by a variant on its own chromosome; a deletion on chromosome 21 cannot lie in an intron of a gene on chromosome X.

Two Sequence Ontology terms appear in this class:

- `intron_variant` (SO:0001627): "a transcript variant occurring within an intron".
- `intergenic_variant` (SO:0001628): "a sequence variant located in the intergenic region, between genes". This is the row Ensembl VEP writes for a variant with no transcript within 5,000 bases.

A breakend is a record for one side of a breakpoint; the other side, called its mate, is named in the allele. A breakend whose mate lies on a different chromosome (an inter-chromosomal breakend) is the one record type whose annotation legitimately reaches a second chromosome. In the allele `N[X:16162317[` the mate lies on chromosome X, and Ensembl VEP loads the transcript-cache regions around the mate so that the mate allele can be annotated there.

#### 9.2 The record

`gnomAD-SV_v2.1_DEL_21_178762`, gnomAD-SV v2.1 chr21 (GRCh37, `gnomad_sv_v2.1_chr21.vcf.gz`), a 75-base deletion:

```
21	16155573	gnomAD-SV_v2.1_DEL_21_178762	N	<DEL>	693	PASS	END=16155648;SVTYPE=DEL;SVLEN=75
```

The transcript Ensembl VEP names is ENST00000380289 (GRPR, forward strand, canonical), on chromosome X at X:16141679-16171144: exon 1 X:16141679-16142489, intron 1 X:16142490-16168427. The deletion's coordinates 16155574-16155648 fall, as numbers, inside that intron.

On chromosome 21 the nearest transcript to the deletion is ENST00000448152 (GAPDHP16, a processed pseudogene at 21:16146622-16147583), 7,991 bases away, beyond the 5,000-base window.

The neighbouring record that matters is a breakend 6.7 kb away whose mate lies on chromosome X at X:16162317, inside GRPR:

```
21	16162316	gnomAD-SV_v2.1_BND_21_56729	N	<BND>	473	UNRESOLVED	END=16162317;SVTYPE=BND;SVLEN=-1;CHR2=X;POS2=16162316;END2=16162317
```

#### 9.3 What each engine writes

```
# Ensembl VEP
gnomAD-SV_v2.1_DEL_21_178762	21:16155574-16155648	deletion	ENSG00000126010	ENST00000380289	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;OverlapBP=75;OverlapPC=0.25
# vep-rs
gnomAD-SV_v2.1_DEL_21_178762	21:16155574-16155648	deletion	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

Ensembl VEP names a transcript on chromosome X and writes `intron_variant`; vep-rs names no transcript and writes `intergenic_variant`.

Ensembl VEP's `OverlapPC=0.25` is the deletion's 75 bases as a share of GRPR's 29,466-base span: the overlap was computed as if the two lay on one chromosome.

#### 9.4 Why Ensembl VEP writes that

The breakend record in the same chunk loads the chromosome X cache regions around its mate. Ensembl VEP then matches loaded transcripts to the records of the batch through an interval tree. The tree holds every record and every mate by coordinate number only, with no chromosome. So GRPR, at X:16141679-16171144, is handed every record of the batch whose numbers fall in that range, the chromosome 21 deletion among them.

Whether a chromosome is then compared depends on the record's shape. A record with a parsed mate has both its ends checked: the local end and each mate pass through `_close_to_feature`, which compares chromosome names and applies the 5,000-base window (the loop over the record and its breakends). A record with no parsed mate gets its one allele unconditionally and is never chromosome-checked. That shape is a `<DEL>`, `<INS>`, `<CPX>` or `<NON_REF>`, and also a single breakend written `G.` or `.C` (a breakpoint whose partner is unknown).

Alone in a batch the deletion loads no chromosome X region and is written as `intergenic_variant`.

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/AnnotationSource.pm`, lines 250-257: `get_all_regions_by_InputBuffer` requests the regions around each breakend's mate.
- `InputBuffer.pm`, lines 345-366: `interval_tree` inserts every record and mate by start and end, no chromosome. Lines 284-330: `get_overlapping_vfs` fetches by coordinates alone.

Source (`ensembl-variation`):

- `Bio/EnsEMBL/Variation/StructuralVariationOverlap.pm`, lines 64-71: a record without breakends gets its allele unconditionally. Against that, lines 72-87 are the loop over the record's own end and its mates, and lines 130-146 are `_close_to_feature`, which compares `seq_region_name`s and applies the 5,000-base window.

#### 9.5 Why vep-rs's row is the correct one

A deletion of 75 bases on chromosome 21 changes no base of chromosome X, so it is not in GRPR's intron whatever the coordinate numbers say. The row is impossible on the record's own coordinates.

On its own chromosome the deletion has no transcript within 7,991 bases and the 5,000-base window is empty, which is `intergenic_variant`. That is the row vep-rs writes, and the row Ensembl VEP itself writes when no chromosome X region happens to be loaded (9.6).

#### 9.6 Try it yourself

Annotate the deletion with Ensembl VEP 115.2, GRCh37 cache and FASTA, the command of the closing section, under each of the four inputs below:

| Input | Ensembl VEP row for the deletion |
| --- | --- |
| the record alone | `intergenic_variant` |
| the record plus the `gnomAD-SV_v2.1_BND_21_56729` line of 9.2 only | ENST00000380289 (GRPR, chromosome X) `intron_variant` |
| the whole `gnomad_sv_v2.1_chr21.vcf.gz`, `--fork 1`, `--fork 4` or `--fork 16` | ENST00000380289 `intron_variant` |
| the same file, `--buffer_size 1` | `intergenic_variant` |

At `--buffer_size 1` every one of the 155 chromosome 21 records of this class in the file writes `intergenic_variant`, the vep-rs row.

#### 9.7 How the comparator treats this class

`filter_cross_chromosome_divergences` in `scripts/validation/compare_sv_concordance.py` sets aside, on Ensembl VEP's side, every tuple whose Feature the transcript cache does not hold for the chromosome of the tuple's Location. Bracket-notation breakend alleles are exempt, since their mate legitimately names another chromosome. A single breakend's `G.` or `.C` allele is not bracket notation and is not exempt.

The test is a fact about the cache, not a mechanism guess. It runs on vep-rs's side as an assertion: vep-rs annotates from a per-chromosome index, so if it ever produced the shape the comparison aborts rather than masking it.

Counts:

- GRCh37: 1,838 Ensembl VEP tuples (gnomAD-SV v2.1 chr21 1,790; `09_breakends` 48).
- GRCh38: 22,481 (gnomAD-SV v4.1 chr21 22,381; `09_breakends` 100).

The synthetic tuples are single breakends (`G.`, `.C`, printed with the Allele `chromosome_breakpoint`) at chromosome 21 positions near 33.03 Mb (GRCh37) and 26.0 Mb (GRCh38). Their numbers fall inside chromosome 1 transcripts that the mates of neighbouring inter-chromosomal records had loaded.

Ensembl VEP's side of the class is fully masked (`excluded_cross_chromosome_perl`); the vep-rs row it leaves standing is class 10.

### 10. The orphaned `intergenic_variant` on the other side of class 9

**In short.** `gnomAD-SV_v2.1_DEL_21_178762` is a 75-base deletion on chromosome 21 with no transcript within 5,000 bases on its own chromosome. Ensembl VEP's only row for it names a transcript on chromosome X, a row class 9 removes, so nothing of the record is left on Ensembl VEP's side; vep-rs writes one row, `intergenic_variant`. A variant with no transcript within 5,000 bases on its own chromosome is intergenic on either engine's terms, and Ensembl VEP itself writes that row whenever the chromosome X region is not loaded, so vep-rs's row is the correct one.

#### 10.1 What a reader needs to know first

Class 9 removes Ensembl VEP's off-chromosome rows for a record. For the records of this class the off-chromosome transcript was the only one Ensembl VEP named, so once class 9 has removed that row, nothing of the record is left on Ensembl VEP's side.

vep-rs's side has one row, `intergenic_variant`: the row for a variant with no transcript within 5,000 bases on its own chromosome. One Sequence Ontology term appears in this class:

- `intergenic_variant` (SO:0001628): "a sequence variant located in the intergenic region, between genes".

Without a rule for this shape, the comparator would count that row as a vep-rs tuple with no counterpart and charge vep-rs for it.

#### 10.2 The record

The same record as class 9, `gnomAD-SV_v2.1_DEL_21_178762` (gnomAD-SV v2.1 chr21, GRCh37):

```
21	16155573	gnomAD-SV_v2.1_DEL_21_178762	N	<DEL>	693	PASS	END=16155648;SVTYPE=DEL;SVLEN=75
```

On chromosome 21 the nearest transcript is ENST00000448152 (GAPDHP16, 21:16146622-16147583), 7,991 bases from the deletion.

#### 10.3 What each engine writes

```
# Ensembl VEP (the record's only row; removed by class 9)
gnomAD-SV_v2.1_DEL_21_178762	21:16155574-16155648	deletion	ENSG00000126010	ENST00000380289	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;OverlapBP=75;OverlapPC=0.25
# vep-rs (the record's only row)
gnomAD-SV_v2.1_DEL_21_178762	21:16155574-16155648	deletion	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

Each engine writes exactly one row for the record. Ensembl VEP's row names ENST00000380289, the off-chromosome transcript, and class 9 removes that row; vep-rs's row names no gene or transcript and carries `intergenic_variant`.

#### 10.4 Why Ensembl VEP writes that

The same as class 9: the off-chromosome transcript, matched by coordinate number, fills the record's transcript list, and Ensembl VEP writes an intergenic row only when that list is empty. With the chromosome X transcript in the list, the intergenic row is never written. Without it, as when the record is annotated alone, it is.

Source (`ensembl-variation`): as in 9.4, plus one module:

- `Bio/EnsEMBL/Variation/StructuralVariationFeature.pm`, `get_IntergenicStructuralVariation`, lines 674-689: the intergenic overlap is created if and only if the transcript overlap list is empty.

#### 10.5 Why vep-rs's row is the correct one

The deletion overlaps no chromosome 21 transcript, and none lies within 5,000 bases of it: the nearest is 7,991 bases away. That is the definition of the intergenic row on either engine's terms.

It is also the row Ensembl VEP writes for the record whenever no chromosome X region is loaded: for the record alone, and for all 155 records of the class at `--buffer_size 1` (9.6).

#### 10.6 Try it yourself

As 9.6: the outcome turns on whether a chromosome X region is loaded alongside the record.

- The record alone, or the whole `gnomad_sv_v2.1_chr21.vcf.gz` at `--buffer_size 1`: Ensembl VEP writes its `intergenic_variant` row, identical to vep-rs's.
- The whole file at `--fork 1`, `--fork 4` or `--fork 16`, or the record beside `gnomAD-SV_v2.1_BND_21_56729` alone: Ensembl VEP writes the GRPR row instead, the Ensembl VEP row of 10.3.

The 1000 Genomes structural-variant files carry records of the class at identical coordinates:

- GRCh37 (`1kg_sv_chr21.vcf.gz`, 1000 Genomes Phase 3, fetched by `scripts/data/download_real_world_vcfs.sh`): one record.
- GRCh38 (`1kg_sv_chr21.vcf.gz`, 1000 Genomes high-coverage): 49 records.

There, in a file with no inter-chromosomal breakend nearby, Ensembl VEP writes `intergenic_variant` for each of them inside the reference output itself.

#### 10.7 How the comparator treats this class

`filter_cross_chromosome_orphan_intergenic` in `scripts/validation/compare_sv_concordance.py` sets aside a vep-rs tuple when all four of these hold.

- Its Feature is `-` and its consequence `intergenic_variant`.
- Ensembl VEP wrote at least one tuple for the same Location and Allele.
- Every Ensembl VEP tuple there was removed by class 9.
- vep-rs's whole set for the record is that one row.

A record with any surviving Ensembl VEP tuple, a record Ensembl VEP wrote nothing for, or a vep-rs set with anything beside the intergenic row is left charged where it falls.

Counts:

- gnomAD-SV v2.1 chr21 (GRCh37): 155 rows.
- gnomAD-SV v4.1 chr21 (GRCh38): 518 rows.

The comparator's pooled per-assembly count (introduction) reads 154 and 469, because 1 and 49 of those records recur in the 1000 Genomes file, where both engines write the intergenic row.

Fully masked, vep-rs side (`excluded_cross_chromosome_rust`).

## Defect E: a breakend's mate-side row evaluated at the record's own coordinate

**The defect in one sentence.** For a breakend, Ensembl VEP chooses which transcripts get a mate-side row by the mate's chromosome and coordinate, then evaluates almost every consequence test on that row with the record's own coordinate, which may lie on another chromosome.

**Why this is a defect.** The row is written for the mate breakpoint: its transcript was chosen because the mate lands in or near it. Asking whether the record's chromosome 1 position lies inside a chromosome 21 transcript's intron answers a question about the wrong breakpoint. Only `feature_truncation` looks at the mate, so the row comes out as bare `feature_truncation` or, when even that does not fire, falls through to `intergenic_variant` on a row that names a transcript.

**What vep-rs does instead.** vep-rs evaluates the mate breakpoint where it is: inside the transcript it writes `feature_truncation` plus the exon, intron, UTR or coding region the breakpoint lands in; within 5,000 bases of the transcript, `upstream_gene_variant` or `downstream_gene_variant`.

**Where it is in Ensembl VEP's source** (`ensembl-variation`):
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 257 and 273: the record's own feature is passed to every predicate.
- `Utils/VariationEffect.pm`, lines 350-358: `feature_truncation` alone tests `within_feature` on the breakend; lines 444-458: `upstream` and `downstream` read the record; lines 629-634: `within_intron` reads the record's differing region, which `BaseTranscriptVariationAllele.pm` line 118 derives from the record's start.
- `StructuralVariationOverlap.pm`, lines 72-87 and 130-146: the mate allele is created from the mate's coordinate and chromosome.

**What correcting it would change.** If every predicate were evaluated at the breakpoint the row is for, class 11 would vanish.

**The classes it produces.**
- Class 11: mate-side rows of point breakends whose mate lies on another chromosome or far from the record (the synthetic `09_breakends` and `13_vcf45_features` sets and gnomAD-SV v4.1 chr21: 1,407 and 7,258 pairs). Set aside: both members of each pair. vep-rs 0.2.0 reproduced Ensembl VEP's rows here and 0.3.0 does not, which, together with class 8, is why its raw F1 on the structural-variant sets is lower for 0.3.0 than for 0.2.0.

### 11. Breakend mate rows read from the other chromosome's coordinate

**In short.** The synthetic breakend `synth_bnd_inter_0066_B` sits on chromosome 1 and names its mate at 21:33034355, inside intron 1 of SOD1 and 2,542 bases upstream of a lincRNA. Ensembl VEP picks the transcripts for the mate-allele rows by the mate's coordinate but evaluates the rows with the record's own chromosome 1 coordinate, so it writes `feature_truncation` alone for SOD1 and `intergenic_variant` for the lincRNA. vep-rs writes `feature_truncation,intron_variant` and `upstream_gene_variant` with `DISTANCE=2542`. vep-rs's rows state what the mate breakpoint does where it is, which is what Ensembl VEP itself writes for the same breakpoint when it is a record's own coordinate, so they are the correct ones.

#### 11.1 What a reader needs to know first

A breakend record has two coordinates: its own (POS on CHROM) and its mate's (the chromosome and position inside the allele's brackets).

Ensembl VEP represents a breakend in one of two ways:

- A `<BND>` whose mate is on the same chromosome and given in `INFO/CHR2` and `END2` becomes a span from POS+1 to the mate (classes 5 and 6). Every allele of a span record is evaluated over the span.
- A breakend written in bracket notation, or a `<BND>` whose mate is on another chromosome, becomes a point record whose Location is the single base after POS.

This class is about the mate-allele rows of point records. For a transcript within 5,000 bases of the mate, Ensembl VEP writes a mate-allele row (6.1). What that row should say is what the mate breakpoint does to that transcript:

- when the breakpoint falls inside it: `feature_truncation` (SO:0001906), plus the region the breakpoint lands in (`intron_variant`, SO:0001627, or an exon term);
- when it falls near it: `upstream_gene_variant` (SO:0001631, "a sequence variant located 5' of a gene") and `downstream_gene_variant` (SO:0001632, "a sequence variant located 3' of a gene") with a DISTANCE.

`intergenic_variant` (SO:0001628) is the row for a position near no transcript. On a row that names a transcript it is a contradiction in terms.

#### 11.2 The record

**An inter-chromosomal breakend** from the synthetic breakend set `tests/sv_validation/grch37/09_breakends.vcf.gz` (GRCh37; the GRCh38 set `tests/sv_validation/grch38/09_breakends.vcf.gz` has the same construction at other coordinates). It is a pair of records that are each other's mates:

```
1	14374333	synth_bnd_inter_0066_B	A	]21:33034355]A	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0066_A
21	33034355	synth_bnd_inter_0066_A	G	]1:14374333]G	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0066_B
```

The rows in question are the first record's mate-allele rows: the transcripts near 21:33034355.

- ENST00000270142 (SOD1, forward strand, canonical, GRCh37 33031935-33041244) contains the mate breakpoint in its intron 1 (33032155-33036102).
- ENST00000449339 (a lincRNA on the reverse strand at 33025906-33031813) ends 2,542 bases before it, which for a reverse-strand transcript is upstream.

The record's own coordinate, 1:14374334, is on another chromosome altogether.

**The coincidence** `gnomAD-SV_v2.1_BND_21_56716`, gnomAD-SV v2.1 chr21 (GRCh37, `gnomad_sv_v2.1_chr21.vcf.gz`):

```
21	15437970	gnomAD-SV_v2.1_BND_21_56716	N	<BND>	892	UNRESOLVED;UNSTABLE_AF_PCRMINUS;LOW_CALL_RATE	END=15437971;SVTYPE=BND;SVLEN=-1;CHR2=Y;POS2=15437970;END2=15437971
```

Its mate, Y:15437971, lies in intron 17 (Y:15436587-15438100) of ENST00000329134 (UTY, reverse strand, Y:15434914-15592550). Its own coordinate, 21:15437971, is the same number on chromosome 21.

#### 11.3 What each engine writes

```
# Ensembl VEP, record synth_bnd_inter_0066_B, mate allele (2 of 7 rows)
synth_bnd_inter_0066_B	1:14374334	]21:33034355]A	ENSG00000142168	ENST00000270142	Transcript	feature_truncation	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
synth_bnd_inter_0066_B	1:14374334	]21:33034355]A	ENSG00000234509	ENST00000449339	Transcript	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
# vep-rs, the same rows
synth_bnd_inter_0066_B	1:14374334	]21:33034355]A	ENSG00000142168	ENST00000270142	Transcript	feature_truncation,intron_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
synth_bnd_inter_0066_B	1:14374334	]21:33034355]A	ENSG00000234509	ENST00000449339	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2542;STRAND=-1
```

The record has seven mate-allele rows and every one differs the same way:

- The four transcripts that contain the breakpoint get `feature_truncation` alone from Ensembl VEP and `feature_truncation` plus the region from vep-rs.
- The three nearby ones get `intergenic_variant` from Ensembl VEP and `upstream_gene_variant` or `downstream_gene_variant` with a distance from vep-rs.

Ensembl VEP's own rows for the second record of the pair are what vep-rs writes for the first. The second record is the same breakpoint written as the record's own line (its local allele, printed `.G`):

```
# Ensembl VEP, record synth_bnd_inter_0066_A, local allele (2 of 7 rows)
synth_bnd_inter_0066_A	21:33034356	.G	ENSG00000142168	ENST00000270142	Transcript	feature_truncation,intron_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
synth_bnd_inter_0066_A	21:33034356	.G	ENSG00000234509	ENST00000449339	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2543;STRAND=-1
```

The distance differs by one because the two rows are the two sides of one junction. Ensembl VEP places a record's own breakpoint on the first base after POS (33034356 for the second record) and a mate on the base the bracket names (33034355 in the first record's allele). The two are one base apart by convention, in both engines.

On the coincidence record both engines write the same row:

```
# Ensembl VEP and vep-rs
gnomAD-SV_v2.1_BND_21_56716	21:15437971	N[Y:15437971[	ENSG00000183878	ENST00000329134	Transcript	feature_truncation,intron_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1
```

#### 11.4 Why Ensembl VEP writes that

Ensembl VEP decides which transcripts get a mate-allele row by the mate's coordinate and chromosome. It then evaluates almost every consequence predicate (the test that decides one consequence term) on that row with the record's own coordinate. Only `feature_truncation` looks at the mate breakpoint; `within_intron`, the exon tests, `upstream` and `downstream` all read the record's POS, on chromosome 1, against a transcript on chromosome 21.

For SOD1 the chromosome 1 number 14374334 is nowhere near the transcript's chromosome 21 coordinates, so the intron test finds nothing and only `feature_truncation` fires. For the lincRNA nothing fires at all, and the row falls through to `intergenic_variant`. When the record's own number happens to land inside or near the transcript, as on the coincidence record, the row comes out right.

Source (`ensembl-variation`):

- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 257 and 273: `$bvf` is the record's own feature and is passed to every predicate.
- `Utils/VariationEffect.pm`, lines 350-358: `feature_truncation` alone tests `within_feature` on `$bvfoa->breakend`; lines 444-458: `upstream` and `downstream` read `$bvf`; lines 629-634: `within_intron` reads the record's differing region, which `BaseTranscriptVariationAllele.pm` line 118 derives from `$vf->{start}`.
- `StructuralVariationOverlap.pm`, lines 72-87 and 130-146: the mate allele is created from the mate's coordinate and chromosome.

#### 11.5 Why vep-rs's row is the correct one

The mate breakpoint 21:33034355 lies in intron 1 of SOD1 (33032155-33036102) and 2,542 bases upstream of the lincRNA's end at 33031813. Those are facts about the record's own coordinates, and they are what vep-rs writes.

Ensembl VEP's row for SOD1 says the transcript is truncated but names no region. Its row for the lincRNA calls a position 2,542 bases from the transcript intergenic, on a line that names the transcript.

Ensembl VEP itself annotates the same breakpoint correctly when it is the record's own coordinate (11.3, the second record of the pair). Its mate row changes when the record's chromosome 1 position is moved while the mate stays fixed (11.6), which no correct reading of the mate row could do.

#### 11.6 Try it yourself

Annotate with Ensembl VEP 115.2, GRCh37 cache and FASTA, using the command of the closing section. Keep the mate at `21:33034355` and move the record's own position:

| Record | Ensembl VEP, mate-allele rows on ENST00000270142 and ENST00000449339 |
| --- | --- |
| `1 14374333 . A ]21:33034355]A` (as in the file) | `feature_truncation`; `intergenic_variant` |
| `1 33034300 . A ]21:33034355]A` | `feature_truncation,intron_variant`; `upstream_gene_variant` with `DISTANCE=2488` |

The second input places the chromosome 1 number inside SOD1's chromosome 21 span and 2,488 bases from the lincRNA, and the rows follow the number. The mate rows describe wherever the record's own coordinate happens to point. vep-rs writes `feature_truncation,intron_variant` and `upstream_gene_variant` with `DISTANCE=2542` for both inputs.

#### 11.7 How the comparator treats this class

`filter_breakend_mate_local_read_pairs` in `scripts/validation/compare_sv_concordance.py`, with the derivation in `scripts/validation/breakend_mate_context.py`, is computed, not pattern-matched. For a pair whose row is a mate-allele row of a point breakend (11.1) it derives both of these sets from the transcript cache:

- the set the mate coordinate supports for the transcript: inside it, `feature_truncation` plus the exon, intron, UTR or coding region the breakpoint lands in, with the non-coding and NMD context terms; within 5,000 bases of it, `upstream_gene_variant` or `downstream_gene_variant`;
- the set the record's own coordinate yields by Ensembl VEP's reading: `feature_truncation` alone when the mate is inside, `intergenic_variant` when it is not.

It sets a pair aside only when the two derivations differ, Ensembl VEP's set equals the second and vep-rs's set equals the first. A vep-rs over- or under-annotation matches neither and stays charged.

Counts, reported as `excluded_breakend_mate_local_read_rust` and `excluded_breakend_mate_local_read_perl`:

- GRCh37: 1,407 pairs (`09_breakends.vcf.gz` 1,369; `13_vcf45_features.vcf.gz` 38).
- GRCh38: 7,258 pairs (`09_breakends.vcf.gz` 4,055; `13_vcf45_features.vcf.gz` 22; `gnomad_sv_v4.1_chr21.vcf.gz` 3,181).

Rows like the coincidence record, where the two derivations agree, are concordant and not in the class. The class is fully masked, both sides.

vep-rs 0.2.0 evaluated no position on the mate side and so reproduced Ensembl VEP's fall-through rows. vep-rs 0.3.0 evaluates the mate breakpoint where it is. Raw F1 on the structural-variant sets is lower for 0.3.0 than for 0.2.0 by these pairs together with the tuples of class 8, and by nothing else, since raw F1 sets nothing aside.

## A representation difference, not a defect: `<CNV:TR>`

**In one sentence.** A `<CNV:TR>` record describes a change in the copy number of a tandem repeat without spelling the alternate allele out, and the two engines read it in two different, equally faithful ways.

**The two readings.** Ensembl VEP materialises the allele: it repeats the unit the stated number of times and annotates the difference as an ordinary insertion or deletion of bases. vep-rs keeps the allele symbolic, as a copy-number change of the run, and annotates the run's span for its direction. Ensembl VEP itself gives the symbolic reading when the expanded allele would exceed `--max_sv_size`.

**Why neither is wrong.** Ensembl VEP's row is a faithful annotation of a real molecule: a run of 121 units where the reference has nine is, base for base, an insertion. vep-rs's row keeps the record in the terms its author chose: a copy-number change of a repeat, located on the repeat, with the region term and the direction term. Nothing in either row is false, so the comparator pairs the two readings rather than ranking them.

**Where the readings part** (`ensembl-vep`): `Bio/EnsEMBL/VEP/Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`, lines 382-431, called from `create_StructuralVariationFeatures` at line 520.

**The class.**
- Class 12: the synthetic tandem-repeat set on both assemblies, 1,368 and 2,644 pairs. Set aside: both members of each pair.

### 12. `<CNV:TR>` read symbolically, not expanded

**In short.** A `<CNV:TR>` record (VCF 4.4) states a copy-number change of a tandem repeat without spelling the alternate allele out. Ensembl VEP rebuilds the alternate sequence and annotates its difference from the reference run as an ordinary insertion or deletion of bases: on a gain its Allele is the inserted bases and, in a coding exon, its term is `inframe_insertion`. vep-rs keeps the allele symbolic: its Location is the run, its Allele is `tandem_repeat`, and its consequence is the region the run lies in, plus `feature_elongation` where the run touches exonic sequence. Neither engine is wrong: the class is a difference of representation, not a defect, and the comparator pairs the two readings rather than ranking them.

#### 12.1 What a reader needs to know first

A `<CNV:TR>` record (VCF 4.4) describes a change in the copy number of a tandem repeat without spelling the alternate allele out. Its fields are:

- POS is the base before the repeat run;
- `INFO/END` is the run's last base;
- `SVLEN` is its length;
- `RUS` is the repeat unit;
- `RUL` is the unit's length;
- `RUC` is the number of units in the alternate allele;
- `RB` is the alternate allele's length in bases;
- `RN=1` says the alternate allele consists of one repeat sequence.

The record `T <CNV:TR>` with `END=25896436;SVLEN=18;RUS=CA;RUC=121;RB=242` says: the run of nine CA units at 25896419-25896436 carries 121 units in the alternate allele, a gain of 112 units (224 bases).

There are two ways to annotate such a record. One is to reconstruct the alternate sequence, compare it with the reference run and annotate the difference as an ordinary insertion or deletion of bases. The other is to keep the allele symbolic, as a copy-number change of the run itself, and annotate the run's span for its direction. Ensembl VEP does the first; vep-rs does the second.

The terms involved are:

- `intron_variant` (SO:0001627).
- `feature_elongation` (SO:0001907): "a sequence variant that causes the extension of a genomic feature, with regard to the reference sequence".
- `coding_sequence_variant` (SO:0001580): "a sequence variant that changes the coding sequence".
- `inframe_insertion` (SO:0001821): "an inframe non synonymous variant that inserts bases into in the coding sequence".

Both engines follow Ensembl's rule for a span allele (a deletion, duplication, insertion or repeat): `feature_elongation` and `feature_truncation` are written only when the span touches the transcript's exonic sequence. The rule is in `ensembl-variation`, `Utils/VariationEffect.pm` lines 337-348 and 350-373, the `within_cdna` tests at lines 344 and 362. So a run inside an intron gets the region term alone, and a run inside an exon gets the region term plus the direction term. A breakend's `feature_truncation` is a different test, the breakpoint falling anywhere inside the transcript, intron included (lines 357-358; classes 6 and 11).

#### 12.2 The record

Both records are from the synthetic tandem-repeat set `tests/sv_validation/grch38/07_cnv_repeat.vcf.gz` (GRCh38; the GRCh37 set is `tests/sv_validation/grch37/07_cnv_repeat.vcf.gz`). Every `<CNV:TR>` record of those sets sits on a tandem repeat the chromosome 21 reference actually carries. The first is a gain in an intron:

```
21	25896418	synth_cnv_tr_0001	T	<CNV:TR>	.	PASS	SVTYPE=CNV;END=25896436;SVLEN=18;RN=1;RUS=CA;RUL=2;RUC=121;RB=242
```

The reference at 21:25896410-25896445 reads `AGTAACACT` `CACACACACACACACACA` `AAACAAAAC` (9 + 18 + 9 bases). The run is the 18 bases 25896419-25896436 (nine CA units), and the alternate has 121 units, a gain of 112 units or 224 bases.

The transcript is ENST00000346798 (APP, reverse strand, canonical, 25880550-26170770); the run lies in its intron 16 (25891869-25897572). The record has 32 rows on each side, one per APP transcript the run lies in.

The second is a gain inside a coding exon:

```
21	26021865	synth_cnv_tr_0023	T	<CNV:TR>	.	PASS	SVTYPE=CNV;END=26021886;SVLEN=21;RN=1;RUS=GTG;RUL=3;RUC=46;RB=138
```

The run is the 21 bases 26021866-26021886 (seven GTG units, reference `GACTCTGTGGTGGTGGTGGTGGTGGTGGCAATGCTG` at 26021860-26021895). The alternate has 46 units, a gain of 39 units or 117 bases, a multiple of three. On ENST00000346798 the run lies in exon 6 (26021840-26022042), coding positions 819-839.

#### 12.3 What each engine writes

```
# Ensembl VEP, record 1
synth_cnv_tr_0001	21:25896436-25896437	CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA	ENSG00000142192	ENST00000346798	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
# vep-rs, record 1
synth_cnv_tr_0001	21:25896419-25896436	tandem_repeat	ENSG00000142192	ENST00000346798	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;OverlapBP=18;OverlapPC=0.01
```

```
# Ensembl VEP, record 2
synth_cnv_tr_0023	21:26021886-26021887	GTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTG	ENSG00000142192	ENST00000346798	Transcript	inframe_insertion	968-969	818-819	273	A/ATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	gcc/gcCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACc	-	IMPACT=MODERATE;STRAND=-1
# vep-rs, record 2
synth_cnv_tr_0023	21:26021866-26021886	tandem_repeat	ENSG00000142192	ENST00000346798	Transcript	feature_elongation,coding_sequence_variant	969-989	819-839	273-280	-	-	-	IMPACT=HIGH;STRAND=-1;OverlapBP=21;OverlapPC=0.01
```

On the first record the two engines agree on the consequence: the run is intronic, so neither writes a direction term. The whole difference is the Location and the Allele. Ensembl VEP's Location is an insertion between the run's last base and the next, and its Allele is the 224 inserted bases; vep-rs's Location is the run, and its Allele is the symbolic class `tandem_repeat`. The comparator keys on both columns, so it sees two tuples.

On the second record the two readings also name different terms: an in-frame insertion of 39 codons, against a coding-sequence change that lengthens the feature.

#### 12.4 Why Ensembl VEP writes that

Ensembl VEP builds the alternate allele as a literal sequence (it materialises the allele). It reads the reference bases of the run and builds the alternate as `RUS` repeated `RUC` times (or `RB / RUL` times when `RUC` is absent). It then hands the pair to its sequence-variant machinery, which trims the common prefix. What is left is a gain as an insertion of whole units after the run's last base, or a loss as a deletion of whole units ending at the run's last base.

The materialisation is skipped when the alternate allele is longer than `--max_sv_size` bases (the same 10,000,000 default), and the record is then annotated symbolically.

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`, lines 382-431: builds the literal alternate allele; called from `create_StructuralVariationFeatures` at line 520 once the reference sequence of the run has been fetched.

#### 12.5 Why vep-rs keeps the symbolic reading, and why Ensembl VEP's reading is not a defect

Ensembl VEP's reading is a faithful annotation of a real molecule. A run of 121 CA units where the reference has nine is, base for base, an insertion of 224 bases after the run, and an in-frame gain of 39 GTG units inside a coding exon does insert 39 codons. Nothing in its rows is false, and this class is not a defect.

When `RUC` is exact the literal reading also carries what the symbolic one does not: the reading frame and the inserted residues, and an IMPACT of MODERATE for an in-frame gain where the symbolic reading's `feature_elongation` is HIGH.

vep-rs keeps the record in the terms the VCF author chose: the variant is a copy-number change of a repeat, its location is the repeat, its allele is the repeat class, and its consequence is the region the run lies in. Where the run touches exonic sequence, the consequence adds what a longer or shorter run does to the feature (`feature_elongation` or `feature_truncation`).

That reading is stable under the information the record carries and does not carry:

- it does not depend on the reference sequence of the run being fetched;
- it does not turn a copy-number call into a literal sequence the caller never asserted (a `<CNV:TR>` record may carry confidence intervals on the unit count, `CIRUC` and `CIRB`, which Ensembl VEP's expansion ignores with a warning);
- it is the reading Ensembl VEP itself gives the record when it does not expand the allele (12.6).

Neither reading is wrong; the comparator pairs them rather than ranking them.

#### 12.6 Try it yourself

Ensembl VEP 115.2, GRCh38 cache and FASTA, the command of the closing section. Lower the expansion cap below the alternate allele's length and Ensembl VEP annotates the record symbolically. Two settings do that for the two records:

- `--max_sv_size 200` on the first record (alternate 242 bases);
- `--max_sv_size 100` on the second record (alternate 138 bases).

With `--max_sv_size 200` the Ensembl VEP row for ENST00000346798 on the first record becomes

```
synth_cnv_tr_0001	21:25896419-25896436	tandem_repeat	ENSG00000142192	ENST00000346798	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;OverlapBP=18;OverlapPC=0.01
```

With `--max_sv_size 100` the row on the second record becomes

```
synth_cnv_tr_0023	21:26021866-26021886	tandem_repeat	ENSG00000142192	ENST00000346798	Transcript	feature_elongation,coding_sequence_variant	969-989	819-839	273-280	-	-	-	IMPACT=HIGH;STRAND=-1;OverlapBP=21;OverlapPC=0.01
```

Both rows are byte for byte the vep-rs rows of 12.3. At the default `--max_sv_size` of 10,000,000 every `<CNV:TR>` record of the synthetic sets is expanded.

#### 12.7 How the comparator treats this class

`filter_cnv_tr_expansion_divergences` in `scripts/validation/compare_sv_concordance.py` pairs a vep-rs tuple whose Allele is `tandem_repeat` with the Ensembl VEP tuple on the same chromosome, Feature and Feature_type that is its literal reading (`_cnv_tr_pair_kind`):

- for a gain, an insertion (Allele all bases) located between the run's last base and the next;
- for a loss, a deletion (Allele `-`) ending at the run's last base.

The pair is set aside, both members, when the two consequence sets agree once the symbolic reading's direction term, where present (`feature_elongation` on a gain, `feature_truncation` on a loss), is disregarded. One further correspondence holds inside a coding sequence: the literal reading's `inframe_insertion`, `inframe_deletion` or `frameshift_variant` against the symbolic reading's `coding_sequence_variant` (`_is_cnv_tr_representation_pair`). Any other difference, a start or splice term one side carries and the other lacks, stays charged.

The class is fully masked:

- GRCh37 set: 1,368 of 1,368 pairs;
- GRCh38 set: 2,644 of 2,644 pairs.

The counters are `excluded_cnvtr_rust` and `excluded_cnvtr_perl`, against the class sizes `cnvtr_swap_pairs_total_rust` and `cnvtr_swap_pairs_total_perl`.

## Scoring against Ensembl VEP 116.2

The classes on this page are the ones the comparators apply when the reference is Ensembl VEP 115.2 (`--reference-release 115.2`), the release the paper was measured against. Ensembl VEP 116 drops a record it marks skipped before its input buffer, so defect C's classes 5 to 8 have no reference row to compare against, and 116.1's chromosome-aware overlap removes defect D's classes 9 and 10 while dropping rows vep-rs writes in three breakend shapes of its own; defects A, B and E and the `<CNV:TR>` representation are unchanged. The classes the comparators apply against 116.2, the default, are documented in [docs/concordance-provenance/116.2-divergence-classes.md](concordance-provenance/116.2-divergence-classes.md); under 115.2 they write exactly the reports they wrote for the paper.

## Reproducing these examples

Every example above is one or two VCF records from a public dataset run through both engines with their default consequence settings.

**The engines.** vep-rs is built from this repository with `cargo build --release`, which gives `target/release/vep`; the commands below write it that way and Ensembl VEP's binary as `vep`, the name it has inside the `ensemblorg/ensembl-vep:release_115.2` image (or a local install of release 115.2).

**The inputs**, all fetched and laid out by `scripts/data/download_real_world_vcfs.sh`, which pins the ClinVar release:

- ClinVar's GRCh37 and GRCh38 VCFs (`clinvar_grch37.vcf.gz`, `clinvar_grch38.vcf.gz`);
- gnomAD-SV v2.1 (GRCh37, `gnomad_sv_v2.1_chr21.vcf.gz`) and v4.1 (GRCh38, `gnomad_sv_v4.1_chr21.vcf.gz`), chromosome 21;
- the 1000 Genomes structural-variant releases, chromosome 21 (`1kg_sv_chr21.vcf.gz` per assembly: Phase 3 for GRCh37, high-coverage for GRCh38);
- ClinVar's structural variants, chromosome 21 (`clinvar_sv_chr21.vcf.gz`);
- the synthetic sets committed under `tests/sv_validation/<assembly>/`, of which this page uses `02_mnp_complex`, `05_symbolic_del_ins`, `07_cnv_repeat`, `09_breakends`, `10_special_alleles`, `11_multi_allelic`, `12_complex_imprecise` and `13_vcf45_features` (each `.vcf.gz`).

Contig names are Ensembl's (`21`, not `chr21`; `scripts/README.md` describes the canonicalisation).

**The transcripts.** The transcript sets are Ensembl release 115's: Ensembl VEP reads its own cache (`homo_sapiens/115_GRCh37` and `homo_sapiens/115_GRCh38`), and vep-rs reads the JSON conversion of the same cache (`scripts/data/storable_to_json.pl`, described in `docs/cache-setup.md`), so both engines see identical transcripts. Transcript coordinates quoted above are read from that cache.

**The reference genome.** Both engines take the assembly's Ensembl primary-assembly FASTA on every example, the SNP/indel ones included: `Homo_sapiens.GRCh38.dna.primary_assembly.fa` (Ensembl release 115) and `Homo_sapiens.GRCh37.75.dna.primary_assembly.fa` (Ensembl's GRCh37 primary assembly, first served under release 75 and unchanged since), each with its `.fai` index. The reference bases quoted above come from them (`samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 21:26965974-26966647` gives the 674 bases of class 3's explicit spelling).

**Writing a probe input.** Every probe input on this page is written as its five leading VCF columns, CHROM POS ID REF ALT, tab-separated. To run one, put it in a file after the two header lines `##fileformat=VCFv4.2` and `#CHROM POS ID REF ALT QUAL FILTER INFO` (tab-separated). Use `.` in QUAL and FILTER and, in INFO, `.` or the keys the class shows (`END=`, `SVTYPE=` and the like); `--format vcf` needs all eight columns. The `alone` inputs are such a file with one record, the `plus one SNV` inputs the same file with the SNV line added.

**The commands.** vep-rs:

```
target/release/vep -i record.vcf -o vep-rs.txt --offline \
  --json_cache /path/to/json_cache/homo_sapiens/115_GRCh38 \
  --fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --species homo_sapiens --assembly GRCh38 \
  --fork 16 --buffer_size 5000 --force --no_stats
```

Ensembl VEP 115.2, with the flags the reference outputs were generated with:

```
vep -i record.vcf -o vep.txt --offline --cache \
  --dir_cache /path/to/vep_cache --cache_version 115 \
  --species homo_sapiens --assembly GRCh38 --format vcf \
  --fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --buffer_size 5000 --fork 4 --no_stats --force_overwrite
```

Substitute `GRCh37` and its cache and FASTA for the GRCh37 examples.

**The batch-dependent classes.** The reference outputs behind the counts used `--fork 4` for the structural-variant sets and `--fork 16` for the SNP/indel datasets. Where a class depends on the batch (5 to 10), the tables above say which `--fork` and `--buffer_size` values give which rows.

The three `--fork` values those tables use are the two reference settings, 4 and 16, and 1, the serial run. The serial run's chunk is the largest the buffer allows (a chunk is at most `buffer_size / (2 × fork)` records, 5.4), so it loads the most regions and writes the most rows; `--buffer_size 1` is the setting under which every record is annotated alone. Both engines accept `--max_sv_size` (default 10,000,000; 5.5 says what each does with it), and class 12 lowers it on Ensembl VEP alone. Ensembl VEP's warnings are in `vep.txt_warnings.txt`; vep-rs's are on standard error.

**The comparators** are `scripts/concordance/compare_vep_outputs.py` (SNP/indel) and `scripts/validation/compare_sv_concordance.py` (structural variants); their per-class counters are named in each class's last section, and `scripts/README.md` describes the full measurement.
