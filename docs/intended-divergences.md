# Intended divergences from Ensembl VEP

This page documents every disagreement between vep-rs and Ensembl VEP that the adjusted concordance figure sets aside, and shows for each one why it is set aside.

vep-rs reports two concordance figures against Ensembl VEP release 116.2. **Raw concordance** keeps every tuple both engines write: a tuple is one output row reduced to its Location, Allele, Feature, Feature_type and Consequence set, and raw F1 counts the tuples the two engines share against the tuples each writes alone, with nothing set aside. **Adjusted concordance** sets aside the tuples of the divergence classes documented on this page and nothing else.

A class is one shape of disagreement. Each class on this page has been traced to a specific mechanism in one engine's source, reproduced with a command a reader can run, on a public record or on a hand-written record printed in full on the page, and judged. Seventeen of the nineteen classes are the output of fourteen defects in Ensembl VEP's source that vep-rs does not reproduce. Class 5 is the records Ensembl VEP 116 drops before annotation: it writes no row for them in any output format, vep-rs annotates them, and with no Ensembl VEP row to judge, vep-rs's rows are set aside. Class 12, `<CNV:TR>`, is a difference of representation in which neither engine is wrong. A disagreement that is not one of these classes is charged to vep-rs in both figures, whatever its cause.

The page is organised by defect. Each defect's section opens with the defect in one sentence, why it is a defect, what vep-rs does instead, where it is in Ensembl VEP's source and what correcting it would change. The classes it produces follow.

Each class is a walkthrough of one record, from a public dataset or hand-written and printed in full on the page. It gives the background a reader needs, the record, both engines' rows, why Ensembl VEP writes what it writes, why vep-rs's row is the correct one, what a reader can reproduce, and how the comparator scripts recognise the class.

**The paper's classes.** The accompanying paper measured vep-rs against Ensembl VEP release 115.2, and its supplement documents five classes at the paper's pinned build, four of them defects: defects A to D below and the representation difference. Against release 115.2 those five produced twelve classes, 1 to 12, and release 116 changed the behaviour behind two of the defects. Defect C, a record marked skipped but annotated against whatever its batch had loaded, is corrected in release 116 by dropping the record before annotation: its classes 5 to 8 compared the two engines' transcript sets on such a record, and against 116.2 there is no Ensembl VEP row to compare, so class 5 now names the dropped records themselves and classes 6 to 8 are retired. Defect D, transcripts matched to records by coordinate alone, is corrected in release 116.1 by a chromosome-aware overlap: its classes 9 and 10 have no rows against it and are retired, and the same overlap introduces two shapes of its own, classes 23 and 24. Defects A, B and E and the `<CNV:TR>` representation are unchanged. Class numbers are identifiers and are kept as the paper's and the comparators' reports use them: 6 to 10 are retired, and 23 and 24 are placed with the mechanism that produces them.

Grouped by the place in Ensembl VEP's source that produces them, the seventeen defect classes are fourteen defects. Two classes belong to one defect when the same wrong lines produce both, so that correcting those lines changes both.

| Defect | What goes wrong in Ensembl VEP | Where | Classes | What correcting it would change |
| --- | --- | --- | --- | --- |
| A | The splice-region verdict is stored by assignment, so a later differing region overwrites an earlier one | `ensembl-variation`, `BaseTranscriptVariationAllele.pm` line 215 | 1, 2 | both classes vanish: `splice_region_variant` is kept wherever any differing region earned it |
| B | `start_lost` and `start_retained_variant` are computed by two independent tests, and nothing keeps both from being true of one allele | `ensembl-variation`, `Utils/VariationEffect.pm` lines 851-963 and 1028-1075 | 3, 4 | both classes vanish: one term per allele, the one the record supports |
| C | A structural record marked skipped is dropped before annotation and gets no row in any output format: a span over `--max_sv_size`, measured as `ensembl-io` derives it (`SVLEN` first, else `INFO/END`, a breakend's `INFO/END` discarded), an ALT type with no Sequence Ontology term, a multi-allelic symbolic list with no type, a deletion with neither `END` nor `SVLEN` or whose start exceeds its end, or a record whose coordinates `validate_vf` refuses. Up to release 115.2 a record on one of the four marked routes stayed in its batch and was annotated against whatever its neighbours had loaded, while a refused coordinate was already dropped; release 116 corrected that by dropping the marked record too | `ensembl-vep`, `Parser.pm` lines 189-191, 475-488, 491-501 and 821-827, `Parser/VCF.pm` lines 477-481, 487-489, 555-558 and 575; `ensembl-io`, `BaseVCF4.pm` lines 297-317 | 5 | nothing is left to correct in Ensembl VEP: the batch dependence is gone, and what remains is a record it declines to annotate and vep-rs annotates |
| D | Retired. Up to release 115.2, transcripts were matched to records by coordinate alone, with no chromosome, so a record could be handed a transcript from another chromosome and lose its `intergenic_variant` row to it; release 116.1's chromosome-aware overlap (`InputBuffer.pm` lines 286-345, `_chr_matches_for_vf_overlap` at lines 348-357) hands a record to a transcript only when the record's chromosome, or a breakend mate's, is the transcript's | `ensembl-vep`, `InputBuffer.pm` | 9, 10, retired | nothing: neither class has a row against release 116.2 |
| P | The tests after the chromosome-aware overlap compare chromosome names the overlap did not resolve the same way: the slice test names a record that arrives without a slice after the first transcript that reaches it, and the mate-allele test compares the mate's chromosome as written, through a synonym table matched verbatim and a case-folded equality, where the overlap had resolved it | `ensembl-vep`, `AnnotationType/Transcript.pm` lines 122 and 125-143, `AnnotationSource.pm` lines 235-258, `InputBuffer.pm` lines 191-207, 223-234 and 348-357; `ensembl-variation`, `StructuralVariationOverlap.pm` lines 86-102 and 144-160, `Utils/VariationEffect.pm` lines 111-124 | 23, 24 | both classes vanish: a bracket breakend's own-end rows and mate rows are written whatever the order of the input and the spelling of the mate's chromosome |
| E | A mate-side row is created from the mate's coordinate, but every predicate except `feature_truncation` is evaluated at the record's own coordinate | `ensembl-variation`, `BaseVariationFeatureOverlapAllele.pm` lines 257 and 273 | 11 | the class vanishes: a mate-side row carries the terms the mate breakpoint supports |
| none | A `<CNV:TR>` allele is expanded to literal sequence and annotated as an insertion or deletion, where vep-rs keeps it symbolic; neither reading is wrong | `ensembl-vep`, `Parser/VCF.pm` lines 382-431 | 12 | nothing to correct; the comparator pairs the two readings |
| F | A breakpoint record without a bracket mate is split on `/` into one allele per part, the reference base included, so the reference allele gets a consequence row of its own | `ensembl-variation`, `StructuralVariationOverlap.pm` lines 64-77 | 13 | the reference-allele row disappears and every row is written once, as vep-rs writes it |
| G | Any coordinate of the input buffer that, with the `--distance` flank added, reaches 2^31 (a breakend mate position, a span's end or a record's own position) turns every record of the buffer intergenic, the records loaded beside it included | `ensembl-vep`, `AnnotationSource.pm` line 190 (the flank added to each coordinate) with lines 251-254 and 306-309 | 14 | the buffer's records keep their transcript rows and the breakend its own-end rows, as vep-rs writes them |
| H | A gVCF variant site, `A G,<NON_REF>`, is routed to the structural path by the `<NON_REF>` token and dropped whole as an unsupported type, the real allele with it | `ensembl-vep`, `Parser/VCF.pm` lines 237-245 and 469-481, `Parser.pm` lines 663-728 | 15 | the real allele is annotated, as vep-rs annotates it |
| I | A dot-notation breakend without `SVTYPE` takes the sequence rules for its coordinates, so its start is POS + 1 when the ALT's length differs from REF's and POS when it does not, its end is POS plus the length of REF minus one, and the row is a span rather than the point (`POS-(POS+1)` with the positions of an insertion for a REF of one or two bases) | `ensembl-io`, `BaseVCF4.pm` lines 255-275 and 297-316 | 16 | the record is the single base at POS + 1, as it is with `SVTYPE=BND` and as vep-rs writes it |
| J | The missing value `.` in the ALT column is annotated as an allele when `SVTYPE` routes the record to the structural path | `ensembl-io`, `BaseVCF4.pm` lines 486-490; `ensembl-vep`, `Parser/VCF.pm` lines 237-245, `Parser.pm` lines 710-711 | 17 | the record writes nothing, as it does without `SVTYPE` and as vep-rs writes |
| K | An `INFO/END` or `SVLEN` that is not an integer, or a POS in a non-integer spelling on the sequence path, is taken as written, passes `looks_like_number`, and the record is annotated with the text in its Location and a fractional distance | `ensembl-io`, `BaseVCF4.pm` lines 302-308; `ensembl-vep`, `Parser.pm` line 475 | 18 | the record is refused, as vep-rs refuses a non-integer `END`, or read at the number its text denotes, as vep-rs reads `SVLEN` and POS |
| L | With a FASTA, a `<CNV:TR>` whose expansion yields no alternate allele, or one equal to the reference run, is dropped without a row or a warning | `ensembl-vep`, `Parser/VCF.pm` lines 382-431 and 516-543 | 19 | the record is annotated symbolically over its run, as it is without a FASTA and as vep-rs annotates it |
| M | `*` as the only ALT names no allele, but a true `SVTYPE` routes the record to the structural path, where it is annotated from INFO alone: a deletion over `END`, or a breakpoint under two alleles | `ensembl-vep`, `Parser/VCF.pm` lines 237-245, 469 and 472-475; `ensembl-variation`, `StructuralVariationOverlap.pm` lines 64-77 | 20 | the record writes nothing, as it does without `SVTYPE` and as vep-rs writes |
| N | A POS that is not a number is numified by the structural path's `$start++` (`.` to 1, `30_425_949` to 31, `-5` to -4) and the record is annotated at a position its line never named | `ensembl-io`, `BaseVCF4.pm` lines 255-264; `ensembl-vep`, `Parser.pm` lines 475-488 | 21 | the record is refused, as Ensembl VEP's own sequence path refuses it and as vep-rs refuses every form |
| O | With a FASTA, a multi-allelic `<CNV:TR>` whose `RN` is 0 for one allele expands that allele to an empty part of the joined allele string, and the other allele is lost: an empty middle part hides the multi-allelic shape, an empty last part is dropped by `split` | `ensembl-vep`, `Parser/VCF.pm` lines 411 and 428, `Parser.pm` lines 881 and 978 | 22 | both alleles are annotated, as vep-rs annotates them |

Source citations name the Ensembl repository, module and line at release 116: the `release/116.2` tag of `ensembl-vep` and the `release/116` branch of `ensembl-variation` and of `ensembl-io`, as bundled in the `ensemblorg/ensembl-vep:release_116.2` image, on which every worked example runs. Where a defect's section recounts what an earlier release did, it says so and names that release. The comparator scripts of this repository are `scripts/concordance/compare_vep_outputs.py` for the SNP/indel datasets and `scripts/validation/compare_sv_concordance.py` for the structural-variant sets. Every row on this page is in Ensembl VEP's default tab output, whose columns are:

```
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
```

Coordinates follow VEP's output conventions: a symbolic deletion's Location starts one base after the VCF POS (the padding base), a breakend's Location is the base after the breakpoint, and an insertion's Location names the two flanking bases.

The nineteen classes, in page order, with the defect each belongs to:

| # | Class | Defect | Kind | Datasets affected | Side set aside | Count |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Splice-region term lost to an assignment order | A | Ensembl VEP defect | ClinVar GRCh37, ClinVar GRCh38 | both members of each pair | 440, 1,994 pairs |
| 2 | Class 1 on the structural-variant sets | A | Ensembl VEP defect | synthetic `02_mnp_complex` (GRCh38) | both members of each pair | 20 pairs |
| 3 | `start_lost` and `start_retained_variant` emitted together | B | Ensembl VEP defect | ClinVar GRCh37 and GRCh38, gnomAD v2.1.1 and v4.1 chr21, 1000 Genomes high-coverage chr21 | both members of each pair | 261, 1,466, 3, 127, 1 pairs |
| 4 | Class 3 on the structural-variant sets | B | Ensembl VEP defect | ClinVar SV, gnomAD-SV, 1000 Genomes SV and the synthetic `05_symbolic_del_ins`, `07_cnv_repeat`, `11_multi_allelic`, `12_complex_imprecise` | both members of each pair | 3 GRCh37 and 38 GRCh38 pairs |
| 5 | A record Ensembl VEP drops before annotation (`reference_skipped_record`) | C | no Ensembl VEP row: the record is dropped | gnomAD-SV v2.1 chr21 (GRCh37) and v4.1 chr21 (GRCh38), 1000 Genomes SV (GRCh38), synthetic `10_special_alleles` and `12_complex_imprecise` (both assemblies) | vep-rs | 20,524 tuples on 389 records (GRCh37); 278,470 on 594 (GRCh38), 278,309 pooled |
| 23 | A bracket breakend's own-end rows lost to the first transcript that reaches it (`bnd_own_chromosome_rows_lost`) | P | Ensembl VEP defect | synthetic breakend set (`09_breakends`), both assemblies | vep-rs, in the structural-variant comparator under the cache authority; counted and not set aside in the SNP/indel comparator | 202 (GRCh37), 1,012 (GRCh38) tuples |
| 24 | A bracket breakend's mate rows lost to an unresolved chromosome prefix (`bnd_synonym_mate_rows_lost`) | P | Ensembl VEP defect | none of the measured sets (every mate is named as the cache names its chromosome) | vep-rs, as class 23 | 0 |
| 11 | Breakend mate rows read from the other chromosome's coordinate | E | Ensembl VEP defect | synthetic `09_breakends` and `13_vcf45_features`, gnomAD-SV v4.1 chr21 | both members of each pair | 1,407 (GRCh37), 10,518 (GRCh38) pairs |
| 12 | `<CNV:TR>` read symbolically, not expanded | none | representation difference | synthetic tandem-repeat set (`07_cnv_repeat`), both assemblies | both members of each pair | 1,368 (GRCh37), 3,751 (GRCh38) pairs |
| 13 | The reference-allele row of a bare breakpoint | F | Ensembl VEP defect | none of the measured sets (every gnomAD-SV `<BND>` names its mate through `CHR2`) | the reference-allele rows are set aside: `bnd_bare_reference_allele_row` in both comparators | 0 |
| 14 | A coordinate at or past 2,147,478,648 (a breakend mate, a span end or the record's own position) turns the reference's buffer intergenic | G | Ensembl VEP defect | none of the measured sets (no chromosome reaches the position) | nothing: counted as `bnd_mate_32bit_buffer_intergenic` in both comparators and not set aside | 0 |
| 15 | A gVCF variant site dropped with its real allele | H | Ensembl VEP defect | none of the measured sets (no gVCF among them) | vep-rs's rows at the site: `gvcf_variant_site_dropped` in both comparators | 0 |
| 16 | A dot-notation breakend without `SVTYPE` written as a pair | I | Ensembl VEP defect | none of the measured sets (the synthetic breakends carry `SVTYPE`) | both rows of the pair: `dot_breakend_without_svtype_pair` in both comparators | 0 |
| 17 | The missing ALT value annotated as an allele | J | Ensembl VEP defect | none of the measured sets | the reference's rows: `missing_alt_annotated` in the structural-variant comparator | 0 |
| 18 | A non-integer `END`, `SVLEN` or POS text annotated and printed | K | Ensembl VEP defect | none of the measured sets | the reference's text rows: `non_integer_end_annotated` in both comparators; vep-rs's rows of the same record in the structural-variant comparator | 0 |
| 19 | An empty or unchanged `<CNV:TR>` expansion dropped without a row | L | Ensembl VEP defect | none of the measured sets (every `<CNV:TR>` of the synthetic sets is one allele with `RN=1` and a positive `RUC` that changes the run) | vep-rs's rows at the run: `cnv_tr_empty_expansion_dropped` in the structural-variant comparator | 0 |
| 20 | `*` alone under a true `SVTYPE` annotated as the type's class | M | Ensembl VEP defect | none of the measured sets | the reference's rows: `star_alone_svtype_annotated` in the structural-variant comparator | 0 |
| 21 | A non-numeric POS numified by the structural path's increment | N | Ensembl VEP defect | none of the measured sets | the reference's rows: `non_numeric_pos_numified` in both comparators (the SNP/indel comparator the negative Locations alone) | 0 |
| 22 | A real `<CNV:TR>` allele dropped beside an `RN=0` allele | O | Ensembl VEP defect | none of the measured sets (no `<CNV:TR>` of the synthetic sets is multi-allelic) | vep-rs's rows under the dropped allele's direction term: `cnv_tr_rn_zero_allele_dropped` in the structural-variant comparator | 0 |

Classes 6 to 10 are retired and their numbers are not reused. Against release 115.2 they were: 6, the mate end of a breakend above `--max_sv_size`, which Ensembl VEP wrote only when a neighbouring record had loaded the region; 7, `<CPX>`, an unsupported type annotated against the batch's transcripts; 8, `<NON_REF>`, a reference block annotated as a variant; 9, transcripts from another chromosome; 10, the `intergenic_variant` row those transcripts displaced. Release 116 drops the records of 6 to 8 before annotation, so they are class 5 now, and release 116.1's overlap never hands a record a transcript from another chromosome, so 9 and 10 have no rows. The paper's class 5, an oversize structural variant annotated against whichever transcripts the batch had loaded, keeps its number: the records are the same, and the difference is that Ensembl VEP now writes nothing for them.

How to read the counts:

- Setting a tuple aside is also called **masking**: a class is fully masked when its rule sets aside every tuple of the class, and a rule is shape-limited when it sets aside only pairs of one exact shape.
- The **six SNP/indel datasets** are ClinVar GRCh37 and GRCh38 (`clinvar_grch37.vcf.gz`, `clinvar_grch38.vcf.gz`), gnomAD v2.1.1 chr21 (GRCh37, `gnomad_genomes_v2.1.1_chr21.vcf.gz`), gnomAD v4.1 chr21 (GRCh38, `gnomad_genomes_v4.1_chr21.vcf.gz`), 1000 Genomes Phase 3 chr21 (GRCh37, `1kg_integrated_chr21.vcf.gz`) and 1000 Genomes high-coverage chr21 (GRCh38, `1kg_highcov_chr21.vcf.gz`), each scored as one file.
- The **two structural-variant sets** are one pool per assembly of sixteen files, scored file by file and pooled: the thirteen synthetic files under `tests/sv_validation/<assembly>/` plus `clinvar_sv_chr21.vcf.gz`, the gnomAD-SV file (`gnomad_sv_v2.1_chr21.vcf.gz` for GRCh37, `gnomad_sv_v4.1_chr21.vcf.gz` for GRCh38) and `1kg_sv_chr21.vcf.gz`.
- Counts are tuples or pairs measured with the comparators of this repository against Ensembl VEP 116.2's output for the same inputs.
- A **pair** is one tuple from each engine that the comparator matches as two readings of one record on one transcript. For classes 1 to 4 and 11 the two tuples share Location, Allele, Feature and Feature_type and differ in the consequence set. For class 12 they share Feature and Feature_type only, because the two engines write the same record at a different Location and Allele (12.3).
- Counts are per input file. The comparator also reports a pooled count per assembly in which a tuple that recurs in two files is counted once (35 tuples on each side recur on GRCh37; 16,731 Ensembl VEP and 16,892 vep-rs tuples on GRCh38, on records two files carry at identical coordinates). Where the two differ this page gives the per-file count and names the pooled one (classes 4 and 5).
- The classes do not overlap on the measured sets: on each structural-variant set the Count column sums to the tuples set aside per file, 23,504 vep-rs and 2,778 Ensembl VEP tuples on GRCh37, 293,809 and 14,327 on GRCh38; pooled, GRCh38 sets aside 293,647 and 14,326, the 161 class 5 tuples and the one class 4 pair that recur in two files counted once.
- On the two ClinVar datasets classes 1 and 3 together are every tuple the adjusted figure sets aside. On gnomAD v2.1.1 and v4.1 chr21 and 1000 Genomes high-coverage chr21 class 3 is, and on 1000 Genomes Phase 3 chr21 nothing is set aside. On the two structural-variant sets classes 2, 4, 5, 11, 12 and 23 together are. The adjusted figure is 1.0 on every dataset.

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
- Class 1: the ClinVar GRCh37 and GRCh38 datasets, 440 and 1,994 pairs; both members of each pair are set aside.
- Class 2: the same shape on the structural-variant sets, scored by the second comparator, 20 pairs on the synthetic complex-substitution set; both members of each pair are set aside.

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

Each part alone is annotated identically by both engines. The two parts as one record lose, in Ensembl VEP only, the term the first part earned. No flag changes the outcome: this class is deterministic.

#### 1.7 How the comparator treats this class

`_is_covered_splice_region_swap` in `scripts/concordance/compare_vep_outputs.py` matches a pair when Ensembl VEP's set has nothing vep-rs's lacks and vep-rs's set adds `splice_region_variant`, alone or together with `intron_variant`. Both members of a matching pair are set aside.

The `intron_variant` allowance exists because Ensembl VEP's per-region loop steps over a region that lies in an intron of twelve bases or fewer before it sets either the intronic flag or the splice-region slot (`BaseTranscriptVariationAllele.pm` lines 139-141 and 172-174). A region skipped that way withholds both terms at once. None of the ClinVar pairs uses the allowance; every one of them adds `splice_region_variant` alone.

The rule is deliberately narrower than the family of splice differences. A pair in which Ensembl VEP lacks a splice sub-term such as `splice_donor_variant`, or in which vep-rs adds `splice_polypyrimidine_tract_variant`, is never set aside, because the assignment at line 215 can remove no term but `splice_region_variant`.

Counts:

- ClinVar GRCh37: 440 pairs; ClinVar GRCh38: 1,994 pairs, every one of that shape.
- gnomAD and 1000 Genomes datasets: none.
- Structural-variant sets: 20 pairs (class 2).

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

The rows agree in every column but Consequence, where vep-rs adds `splice_region_variant`. The record has 47 rows on each side, and 19 of them, all SOD1 transcripts, differ the same way: eleven coding transcripts (`frameshift_variant`), five coding transcripts flagged as NMD targets (`frameshift_variant,NMD_transcript_variant`), one whose exon 2 is untranslated (`3_prime_UTR_variant,NMD_transcript_variant`) and two non-coding ones (`non_coding_transcript_exon_variant`); the other 28 rows are identical.

#### 2.4 Why Ensembl VEP writes that

As 1.4: the first differing region (exon base 2) sets `splice_region` to 1 and the second (exon bases 4 and 5) sets it back to 0, on each of the nineteen transcripts.

#### 2.5 Why vep-rs's row is the correct one

As 1.5: base 31663791 is the second base of exon 2, inside the window the term names, and Ensembl VEP calls it so when it is changed alone.

#### 2.6 Try it yourself

Ensembl VEP 116.2, GRCh38 cache and FASTA, the command of the closing section:

- `21 31663791 . A C` alone gives ENST00000270142 `missense_variant,splice_region_variant`;
- `21 31663792 . AAG A` alone (the second differing region, Location `21:31663793-31663794`, Allele `-`) gives `frameshift_variant`;
- the full record `21 31663791 . AAAG CA` gives `frameshift_variant`.

vep-rs gives the same two rows for the two parts and `frameshift_variant,splice_region_variant` for the full record.

#### 2.7 How the comparator treats this class

`filter_registry_swap_pairs` in `scripts/validation/compare_sv_concordance.py` builds the pairs of the structural-variant comparison, one tuple on each side per Location, Allele, Feature and Feature_type. It applies to them the same two rules the SNP/indel comparator exports (`EXCLUDING_RULES` of `scripts/concordance/compare_vep_outputs.py`, imported so the two scripts cannot drift), setting both members of a matching pair aside.

Counts: 20 pairs on the GRCh38 set and none on the GRCh37 set. Nineteen are the transcripts of `synth_cpx_sub_0780`; the twentieth is ENST00001077173 on `synth_mnp_0852` (`21 31667270 synth_mnp_0852 CTT TTA`, whose two differing bases, 31667270 and 31667272, stand either side of an unchanged base), `NMD_transcript_variant,synonymous_variant` against `NMD_transcript_variant,splice_region_variant,synonymous_variant`. Fully masked; the report counts them as `excluded_registry_swap_rust` and `excluded_registry_swap_perl`, split by rule in `excluded_registry_swap_by_bucket` (`covered_splice_region_swap`).

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
- Class 3: the SNP/indel datasets (ClinVar GRCh37 and GRCh38, gnomAD v2.1.1 and v4.1 chr21, 1000 Genomes high-coverage chr21: 261, 1,466, 3, 127 and 1 pairs); both members of each pair are set aside.
- Class 4: the structural-variant sets, 3 pairs on GRCh37 and 38 on GRCh38 (37 distinct); both members of each pair are set aside.

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

In each pair the rows agree in every column but Consequence: Ensembl VEP writes both start terms, vep-rs writes one. The second record's other four coding transcripts, ENST00000970346, ENST00001108182, ENST00001108183 and ENST00001131497, carry the same difference, and its sixth row, `upstream_gene_variant` on ENST00000652031, is identical on both sides.

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

The sequence deletion has two VCF spellings of one edit, because the deleted A is part of a two-base run. Annotate both with Ensembl VEP 116.2 (GRCh37 cache and FASTA, the command of the closing section):

- The left-aligned `21 36265245 . AT A` above gives ENST00000486278 `frameshift_variant,start_lost,start_retained_variant`.
- The right-aligned `21 36265246 . TT T` (ClinVar's own `g.36265247del`) gives `5_prime_UTR_variant` with no start term at all; vep-rs writes the same `5_prime_UTR_variant` for that spelling.

One molecule, and Ensembl VEP calls its start codon both lost and untouched depending on which base of the run the record names.

The structural deletion has an explicit spelling too: write the same 673 deleted bases out. REF is the reference sequence of 21:26965974-26966647, 674 bases (`samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 21:26965974-26966647`); ALT is its first base, `G`.

- The explicit spelling: Ensembl VEP evaluates the edited sequence and writes ENST00000284987 `start_lost,5_prime_UTR_variant`, with no `start_retained_variant`. Its sequence-variant rules also drop `feature_truncation`, a structural term, and `frameshift_variant`, because the 417 deleted coding bases are a multiple of three; the probe is about the start pair. vep-rs writes the same `start_lost,5_prime_UTR_variant` row for the explicit spelling.
- The `<DEL>` spelling of the same interval receives both start terms.

#### 3.7 How the comparator treats this class

`_is_start_cooccurrence_swap` in `scripts/concordance/compare_vep_outputs.py` matches a pair when Ensembl VEP's set contains both `start_lost` and `start_retained_variant` and vep-rs's set is the same set minus exactly the erroneous member for the allele kind: `start_lost` on a sequence allele, `start_retained_variant` on a structural one. An allele written as bases or `-` is a sequence allele; a Sequence Ontology class name such as `deletion` is structural. A vep-rs set that adds any other term, or drops the other member, is not matched and stays charged. Both members of a matching pair are set aside.

The class is fully masked. Counts, from the SNP/indel summary's `start_cooccurrence_swap` line:

- ClinVar GRCh37: 261 pairs.
- ClinVar GRCh38: 1,466 pairs.
- gnomAD v2.1.1 chr21: 3 pairs.
- gnomAD v4.1 chr21: 127 pairs.
- 1000 Genomes high-coverage chr21: 1 pair.

Its reach into the structural-variant sets is class 4.

### 4. Class 3 on the structural-variant sets

**In short.** Structural deletions and gains over a start codon in the structural-variant files reach the structural arm of class 3. Ensembl VEP writes `start_lost` and `start_retained_variant` together on them; vep-rs keeps one of the two. For the structural deletions, whose interval contains the codon, the argument of 3.5 holds and vep-rs's row is the correct one. For the structural gains, ten pairs on five records, the record decides only which member came from evidence, and the comparator sets both members aside rather than deciding it.

#### 4.1 What a reader needs to know first

The structural-variant sets are scored by a different script from the SNP/indel datasets, and structural deletions and gains over a start codon reach the structural arm of class 3. Nothing new happens in either engine; the same predicate is applied by the second comparator. The terms are those of 3.1.

#### 4.2 The record

The structural deletion of 3.2, `gnomAD-SV_v3_DEL_chr21_88237371` on ENST00000284987 (gnomAD-SV v4.1 chr21, GRCh38), and the RUNX1 deletion of 3.2, which the ClinVar structural-variant file `clinvar_sv_chr21.vcf.gz` (GRCh37) also carries.

#### 4.3 What each engine writes

The RUNX1 deletion's rows in `clinvar_sv_chr21.vcf.gz`, which carries the record under the same ID, are those of 3.3. The `<CN2>` gain `synth_cn2_0016` of 4.5, on ENST00000647188 (DYRK1A, forward strand):

```
# Ensembl VEP
synth_cn2_0016	21:37367280-37521417	duplication	ENSG00000157540	ENST00000647188	Transcript	start_lost,feature_elongation,start_retained_variant,coding_sequence_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant	508-12083	-	-	-	-	-	IMPACT=HIGH;STRAND=1;OverlapBP=154138;OverlapPC=96.59
# vep-rs
synth_cn2_0016	21:37367280-37521417	duplication	ENSG00000157540	ENST00000647188	Transcript	start_lost,feature_elongation,coding_sequence_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant	508-12083	-	-	-	-	-	IMPACT=HIGH;STRAND=1;OverlapBP=154138;OverlapPC=96.59
```

The record's other four pairs, on ENST00000643624, ENST00000644942, ENST00000647504 and ENST00001131494, differ the same way; its other 29 rows are identical on both sides.

#### 4.4 Why Ensembl VEP writes that

As 3.4 for the deletions. On a structural gain, the same two structural tests apply:

- `start_lost` from the overlap of the stated span with the codon;
- `start_retained_variant` from the retention test that returns 0 for every structural allele.

#### 4.5 Why vep-rs's row is the correct one for the deletions, and what the structural gains show

As 3.5 for the 30 structural deletions of the class (2 on GRCh37, 28 on GRCh38; 29 pooled, one deletion recurring in two GRCh38 files), whose interval contains the codon.

The ten structural-gain pairs of the class, on five GRCh38 records, all over a start codon, are:

- an Alu insertion on two transcripts (gnomAD-SV v4.1 chr21, `21:39660415-39660694`);
- an insertion on one transcript (`05_symbolic_del_ins`, `21:37420264-37526076`, ENST00000646548);
- the `<CN2>` copy gain `synth_cn2_0016` on five transcripts (`07_cnv_repeat`, `21:37367280-37521417`);
- two duplications on one transcript each (gnomAD-SV v4.1 chr21, `21:30487315-30487396` and `21:46142001-46156000`).

For them the deletion argument does not apply, and a symbolic gain does not say what becomes of the codon.

What the record does decide is which member came from evidence. Ensembl VEP's structural `start_lost` is the overlap of the stated span with the codon's three coordinates (3.4). Its `start_retained_variant` on a structural allele comes from a test that returns "retained" for every structural allele without reading the sequence. So the retained member is the one emitted from nothing and is the contradiction, and vep-rs keeps the member computed from the record.

That is a statement about which test read the record, not about the biology. Whether a gain over the codon leaves the start in place is not decidable from a symbolic record (a tandem duplication keeps its first copy of the codon). So for these ten pairs the page claims only that vep-rs's row drops the term emitted without evidence, and the comparator sets both members aside rather than deciding it.

#### 4.6 Try it yourself

For the deletions, 3.6. The structural gains are the records named in 4.5; the `<CN2>` record is `21 37367279 synth_cn2_0016 T <CN2>` with `END=37521417`.

#### 4.7 How the comparator treats this class

The same `filter_registry_swap_pairs` as class 2 applies the class 3 rule to the pairs of the structural-variant comparison.

Counts: the class 3 rule matches 3 pairs on GRCh37 and 38 on GRCh38 per file; pooled per assembly, 3 and 37.

- GRCh37, 3 pairs: one each in `11_multi_allelic.vcf.gz` (`synth_multi_mix_0476`, a deletion), `clinvar_sv_chr21.vcf.gz` (where it is the RUNX1 deletion of 3.2) and `gnomad_sv_v2.1_chr21.vcf.gz` (`gnomAD-SV_v2.1_DEL_21_179957`).
- GRCh38, 38 pairs per file: `gnomad_sv_v4.1_chr21.vcf.gz` 30, `07_cnv_repeat.vcf.gz` 5, and one each in `05_symbolic_del_ins.vcf.gz`, `12_complex_imprecise.vcf.gz` and `1kg_sv_chr21.vcf.gz`. Pooled, 37: the deletion `21:30755356-30755423` on ENST00000335093 is `HGSV_242982` in the 1000 Genomes file and `gnomAD-SV_v3_DEL_chr21_d911172a` in the gnomAD-SV file, one key counted once.

The tandem-repeat file, `07_cnv_repeat.vcf.gz`, carries 700 `<CNV>` and `<CN0>` to `<CN9>` copy-number records beside its 200 `<CNV:TR>` records, and Ensembl VEP expands only the `<CNV:TR>` ones (12.4). Its five class 3 pairs are the `<CN2>` copy gain `synth_cn2_0016`, Location 21:37367280-37521417, on the five transcripts whose start codon it covers, ENST00000647188 among them.

Of the 41 per-file pairs of both assemblies (40 pooled), 30 are structural deletions whose interval contains the codon (29 pooled) and 10 are structural gains, where vep-rs keeps `start_lost`. The RUNX1 record is the one sequence allele, where it keeps `start_retained_variant`.

Fully masked; the report counts them as `excluded_registry_swap_rust` and `excluded_registry_swap_perl`, split by rule in `excluded_registry_swap_by_bucket`.

## Defect C, restated for release 116: a record marked skipped is dropped before annotation

**The difference in one sentence.** Ensembl VEP 116 drops a structural record it marks skipped before the record reaches the annotation pass, so the record has no row in any output format, while vep-rs annotates it.

**What release 116 corrected.** Up to release 115.2 the skipped record stayed in its batch: it loaded no transcript regions of its own and was annotated against whatever regions its neighbours had loaded (`AnnotationSource.pm` line 238 against `InputBuffer.pm` lines 284-330 of that release). Its output changed with `--fork`, `--buffer_size` and the order of lines in the file while the warnings file reported the record as skipped; one record got 51 transcripts, or 497, or none, depending on its neighbours. That batch dependence was defect C, and it produced the paper's classes 5 to 8: an oversize structural variant annotated against whichever transcripts the batch had loaded, the mate end of such a breakend, `<CPX>`, and a `<NON_REF>` block annotated as a variant. From release 116 the parser drops the record instead (`validate_svf` returns 0 for a `vep_skip` record and `next` discards it unless `--dont_skip` is set), which is the correction those classes called for. What remains is not a defect but a difference of coverage: Ensembl VEP writes nothing for the record, vep-rs annotates it, and no row of either engine can be judged against the other's. This section documents that difference and the records it reaches.

**Which records are dropped.** A record takes one of four routes to the mark, or a fifth to the same drop:

- a structural variant whose span exceeds `--max_sv_size` (default 10,000,000; `Config.pm` line 315). The span is the one `ensembl-io` derives for the record (`BaseVCF4.pm` lines 297-317): from `SVLEN` first, else from `INFO/END`, else the record's start itself; and the size test compares the span's end minus its start with the cap (`Parser.pm` lines 491-501). A breakend's `INFO/END` is discarded before the span is derived (`Parser/VCF.pm` lines 487-489), so a breakend's span is its `SVLEN` where the record carries one. A gnomAD-SV breakend record carries an `SVLEN` to its mate, so a same-chromosome breakend whose mate lies more than 10,000,000 bases away is dropped whole, its mate end included;
- an ALT type with no Sequence Ontology term, `<CPX>`, `<CTX>` and `<NON_REF>` among them: `get_SO_term` (`Parser.pm` lines 663-728) finds no term for the type, and `Parser/VCF.pm` lines 477-481 mark the line, the type string standing in for the term;
- a multi-allelic symbolic list with no type of its own (`<INS>,<DEL>`, `<INV>,<DUP>`): `get_SO_term` reads one type from the ALTs joined by `/` and finds none;
- a `<DEL>` with neither `END` nor `SVLEN`, or whose start exceeds its end (`Parser/VCF.pm` lines 555-558, "deletion looks incomplete");
- a record whose start or end is not a number, or whose end lies more than one base before its start (`END=0`, `END=abc`, a breakpoint's `END=`): `validate_vf`'s own checks return 0 before any mark is read (`Parser.pm` lines 475-488), and the record is dropped the same way.

The mark is set on the feature at `Parser/VCF.pm` line 575. `validate_vf` then hands a structural feature to `validate_svf` (`Parser.pm` line 551), which returns 0 for a marked record (lines 821-827), and `next` at lines 189-191 discards it unless `--dont_skip` is set. The warning is still written to the warnings file, naming the line and the reason.

**What vep-rs does instead.** vep-rs annotates every record of the four marked routes, and refuses the fifth route's records at warning level as Ensembl VEP does, so neither engine writes a row for those. In the default, tab and Parquet outputs an oversize or unsupported record is annotated against every transcript it overlaps, and a list with no type gets one row per ALT under each ALT's own class; an unsupported ALT whose letters name a class is written under that class's term (`deletion` for `<del>`, `<Del>`, `<DEL:>` or `<DEL >`, `duplication` for `<Dup>`, `inversion` for `<Inv>`, `copy_number_variation` for `<CN=X>`, `insertion` for `<Ins>`, `chromosome_breakpoint` for `<Bnd>`) and one whose letters name none under its raw type (`CPX`); a `<NON_REF>` block, which asserts the absence of a variant, gets one `intergenic_variant` row and a warning on standard error; in the VCF output the line of an oversize or unsupported record is carried without consequences; in the JSON output an unsupported record is an `input`-only object and an oversize one is omitted; a `<DEL>` with neither `END` nor `SVLEN` is read as a deletion of the base after POS, and one whose start exceeds its end (`END` at POS, `SVLEN=0`, a `SVLEN` that is not a number) as the pair `POS-(POS+1)` with an insertion's positions (`21:30425650-30425651 deletion` for `21 30425650 x N <DEL> . . SVTYPE=DEL;END=30425650`), each written in every format. vep-rs takes `--max_sv_size` too, with the same default, and `--max_sv_size -1` lifts the cap on both engines.

**Why the rows are set aside.** The divergence is every row vep-rs writes for such a record, and there is no Ensembl VEP row at the record to judge it against: Ensembl VEP's output does not say what the record's transcripts are, only that the record was not annotated. The comparators therefore set vep-rs's rows aside rather than charge them, and the golden corpora record those rows in their manifests under the class `reference_skipped_record` (`tests/golden/116/GRCh37`: four `<BND>` records above `--max_sv_size` and two `<CPX>`; `tests/golden/116/GRCh38-release116`: `<INS>,<DEL>`, `<INV>,<DUP>` and a `<DEL>` without `END`).

**The class it produces.**
- Class 5: the dropped records of the measured sets. GRCh37: 389 records and 20,524 vep-rs tuples, in gnomAD-SV v2.1 chr21 (89 records), the synthetic special-allele set (its 200 `<NON_REF>` blocks) and the synthetic complex set (its 100 `<CPX>` records). GRCh38: 594 records and 278,470 tuples per file (278,309 pooled: 161 tuples of records the gnomAD-SV and 1000 Genomes files carry at identical coordinates are counted once), in gnomAD-SV v4.1 chr21 (245 records), 1000 Genomes SV (49), the special-allele set (200) and the complex set (100). Set aside: vep-rs's rows.

### 5. A record Ensembl VEP drops before annotation

**In short.** Three records of the measured sets take three of the routes to the drop. `gnomAD-SV_v3_BND_chr21_3f3f16e7` is a gnomAD-SV breakend whose `SVLEN` to its mate, 27 Mb on the same chromosome, exceeds `--max_sv_size`; `gnomAD-SV_v3_CPX_chr21_df09f6c2` is a 1.66 Mb complex rearrangement whose allele `<CPX>` has no Sequence Ontology term; `synth_nonref_0033` is a 279-base `<NON_REF>` block over a coding exon. Ensembl VEP 116.2 writes no row for any of the three, in any output format, and names each line in its warnings file. vep-rs annotates the breakend and the complex event against every transcript their spans overlap or come within 5,000 bases of, and writes one `intergenic_variant` row for the block. With no Ensembl VEP row to compare, vep-rs's rows are set aside.

#### 5.1 What a reader needs to know first

A breakend record (`SVTYPE=BND`) describes one end of a rearrangement: the position where the reference is broken and the position it is joined to, its mate. gnomAD-SV writes the record with a symbolic `<BND>` allele and puts the mate in `INFO/CHR2` and `INFO/END2`; VCF's own notation writes the mate into the allele, `N[chr21:41770788[`. Both engines read either form, and when the mate is on the same chromosome both represent the record as a feature spanning from the base after POS to the mate (Location `21:14604132-41770788`), with rows under two alleles: the local end, `N.`, and the mate, `N[chr21:41770788[`.

`--max_sv_size` (default 10,000,000) is the documented size cap of both engines. A structural variant whose span exceeds it is reported as skipped in the warnings file Ensembl VEP writes beside its output (`<output>_warnings.txt`); Ensembl VEP's source calls the cap a guard "to avoid memory problems". The span Ensembl VEP measures is derived from `SVLEN` first: for the breakend below, `SVLEN=27166657` gives a span of 27,166,656 bases, whatever `END` says.

`<CPX>` is gnomAD-SV's symbolic allele for a complex rearrangement; here the rearrangement is a deletion followed by an inversion, `CPX_TYPE=delINV`. Ensembl VEP has no Sequence Ontology class for it. vep-rs annotates the record for the span alone, under the type string `CPX` and without a direction term (`feature_elongation` or `feature_truncation`): its rows claim overlap, not disruption, with `coding_transcript_variant` (SO:0001968) and `non_coding_transcript_variant` (SO:0001619, "a transcript variant of a non coding RNA gene") for transcripts inside the span, `intron_variant` and the exon terms for transcripts the span partly overlaps, and `mature_miRNA_variant` for a microRNA it covers.

A gVCF (genomic VCF) reports every position of a sample, not only the variant ones, and compresses the positions where the sample matches the reference into reference-confidence blocks: one record whose ALT is the symbolic `<NON_REF>` (or its VCF 4.3 spelling `<*>`, the same construct), spanning POS to `INFO/END`, meaning "no variant was called in this interval, with this confidence". Such a record asserts the absence of a variant; there is no allele to predict a consequence for. Ensembl VEP has no Sequence Ontology class for `<NON_REF>` and drops the block. vep-rs writes the record's one output line in Ensembl VEP's own convention for a record that overlaps nothing it annotated: Feature `-`, consequence `intergenic_variant` (SO:0001628), IMPACT MODIFIER, and warns once per record on standard error. That row is a convention, not a claim that the block lies between genes, and it keeps one output line per input record, which is what the default output promises. `<*>` is not marked skipped by Ensembl VEP: both engines annotate it as a span and their rows agree, so it is outside the class by measurement, not by principle.

#### 5.2 The records

`gnomAD-SV_v3_BND_chr21_3f3f16e7`, gnomAD-SV v4.1 chr21 (GRCh38, `gnomad_sv_v4.1_chr21.vcf.gz`), a breakend at 21:14604131 whose mate lies 27 Mb away on the same chromosome:

```
21	14604131	gnomAD-SV_v3_BND_chr21_3f3f16e7	N	<BND>	251	UNRESOLVED	CHR2=chr21;END=14604131;END2=41770788;POS2=41770787;SVLEN=27166657;SVTYPE=BND
```

`gnomAD-SV_v3_CPX_chr21_df09f6c2`, from the same file, a 1.66 Mb complex event, well below `--max_sv_size`:

```
21	14923958	gnomAD-SV_v3_CPX_chr21_df09f6c2	N	<CPX>	1	UNRESOLVED	CHR2=chr21;CPX_INTERVALS=DEL_chr21:14923958-15004545,INV_chr21:15004545-16588073;CPX_TYPE=delINV;END=16588080;SVLEN=1664122;SVTYPE=CPX
```

`synth_nonref_0033`, from the synthetic special-allele set `tests/sv_validation/grch37/10_special_alleles.vcf.gz` (GRCh37; the GRCh38 set is `tests/sv_validation/grch38/10_special_alleles.vcf.gz`), a 279-base block whose interval 27347259-27347537 lies in APP (ENST00000346798, reverse strand, canonical, GRCh37 27252861-27542972) and covers 155 of the 159 bases of exon 11 (27347383-27347541, a coding exon) and 124 bases of intron 11 (27328070-27347382):

```
21	27347258	synth_nonref_0033	C	<NON_REF>	.	PASS	END=27347537
```

#### 5.3 What each engine writes

Ensembl VEP 116.2 writes no row for any of the three records, in any output format; its warnings file names each line as skipped, with the reason. vep-rs writes 4,652 rows for the breakend, 4,628 under the local allele `N.` for the transcripts its span 14604132-41770788 overlaps or comes within 5,000 bases of, and 24 under the mate allele `N[chr21:41770788[` for the transcripts within 5,000 bases of the mate breakpoint 41770788; 534 rows for the complex event, every transcript inside 14923959-16588080 under the allele `CPX`; and one `intergenic_variant` row for the block.

```
# Ensembl VEP, warnings file, each record run alone
WARNING: line 1 skipped (21 14604131 gnomAD-SV_v3_BND_chr21_3f3f16e7 N ...): variant size (27166656) is bigger than --max_sv_size (10000000)
WARNING: line 1 skipped (21 14923958 gnomAD-SV_v3_CPX_chr21_df09f6c2 N ...): CPX is not a supported structural variant type
WARNING: line 1 skipped (21 27347258 synth_nonref_0033 C <NON_REF> . PA...): NON_REF is not a supported structural variant type
# Ensembl VEP, output
(no row for any of the three records)
# vep-rs, the breakend (2 of 4,652 rows: the SAMSN1 row under the local allele and a mate-allele row)
gnomAD-SV_v3_BND_chr21_3f3f16e7	21:14604132-41770788	N.	ENSG00000155307	ENST00000400562	Transcript	feature_truncation,non_coding_transcript_exon_variant,intron_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1
gnomAD-SV_v3_BND_chr21_3f3f16e7	21:14604132-41770788	N[chr21:41770788[	ENSG00000289513	ENST00000700893	Transcript	feature_truncation,non_coding_transcript_exon_variant,intron_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
# vep-rs, the complex event (2 of 534 rows: USP25 and the microRNA MIRLET7C)
gnomAD-SV_v3_CPX_chr21_df09f6c2	21:14923959-16588080	CPX	ENSG00000155313	ENST00000285679	Transcript	coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;OverlapBP=150140;OverlapPC=100.00
gnomAD-SV_v3_CPX_chr21_df09f6c2	21:14923959-16588080	CPX	ENSG00000199030	ENST00000362160	Transcript	mature_miRNA_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;OverlapBP=84;OverlapPC=100.00
# vep-rs, the block
synth_nonref_0033	21:27347259-27347537	<NON_REF>	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

In the reference outputs of the whole files the same three warnings name line 9777 and line 9978 of the gnomAD-SV file and line 188 of the special-allele set, and the outputs carry no row for the records. vep-rs warns once for the block, on standard error: `21:27347258 <NON_REF>: NON_REF is not a supported structural variant type; a reference-confidence block is written as one intergenic row`.

#### 5.4 Why Ensembl VEP writes nothing

Each record is marked skipped on its route and dropped by the parser before it reaches the input buffer:

- The breakend: `Parser/VCF.pm` lines 487-489 discard its `INFO/END` and re-derive the end, which `BaseVCF4.pm` lines 302-306 take from `SVLEN`: the span is 27,166,656 bases, and the size test at `Parser.pm` lines 491-501 marks the record and writes `variant size (27166656) is bigger than --max_sv_size (10000000)`.
- The complex event: `get_SO_term` (`Parser.pm` lines 663-728) has no term for `CPX`, so `Parser/VCF.pm` lines 477-481 mark the line and write `CPX is not a supported structural variant type`.
- The block: the same route, `NON_REF is not a supported structural variant type`.

Then `validate_svf` (`Parser.pm` lines 821-827) returns 0 for each, and `next` at lines 189-191 discards it. Nothing of the record reaches the annotation pass or any output writer. With `--dont_skip` the marked record enters the buffer, loads no transcript regions of its own (`AnnotationSource.pm` line 238) and is annotated against whatever its batch has loaded, as releases up to 115.2 annotated it.

#### 5.5 Why vep-rs annotates the record

The record's coordinates alone determine which transcripts it reaches. The breakend's span 14604132-41770788 overlaps, or comes within 5,000 bases of, a definite set of chromosome 21 transcripts, the same from any file; the complex event's span 14923959-16588080 contains its transcripts as surely as any deletion contains its own; and a `<NON_REF>` block changes no base, so no transcript's sequence, splicing or dosage is affected by it, which is what the intergenic row says.

Ensembl VEP's cap is a memory guard, not a statement about the record, and its missing term for `<CPX>` is a gap in its vocabulary, not a reading of the event. Skipping is a legitimate output for both, and vep-rs gives that skip in its VCF and JSON outputs, where it carries the line without consequences or omits it. On the default, tab and Parquet outputs it writes the complete set, which is what a caller who kept the record in the file asked for. Neither engine's output here can be checked against the other's, which is why the rows are set aside rather than defended as the correct ones.

#### 5.6 Try it yourself

Run Ensembl VEP 116.2 and vep-rs (GRCh38 cache and FASTA for the gnomAD-SV records, GRCh37 for the block, the commands of the closing section) on inputs built from the records:

| Input | Ensembl VEP 116.2 rows for the record | vep-rs rows for the record |
| --- | --- | --- |
| the breakend alone | 0, the size warning | 4,652 (4,628 local, 24 mate) |
| the breakend plus the SNV `21 14604131 . A C` | 0 (the SNV gets its 4 rows) | 4,652 |
| the breakend plus the SNV `21 41770788 . C A`, at the mate breakpoint | 0 (the SNV gets its 24 rows) | 4,652 |
| the whole `gnomad_sv_v4.1_chr21.vcf.gz`, `--fork 4 --buffer_size 5000` | 0, line 9777 skipped | 4,652 |
| the complex event alone | 0, the type warning | 534 |
| the complex event plus the four SNVs `21 14923958 . T A`, `21 15500000 . G A`, `21 16400000 . G A` and `21 16590000 . A C` | 0 (the SNVs get their 368 rows) | 534 |
| the block alone | 0, the type warning | 1 |
| `synth_nonref_0188` (`21 38589780 synth_nonref_0188 A <NON_REF>` with `END=38590008`, a 228-base block in a DSCR9 intron) plus the SNV `21 38589780 . A C` | 0 (the SNV gets its 8 rows) | 1 |

No neighbouring record changes the outcome: under release 115.2 the breakend beside the SNV at its mate breakpoint got 20 mate-allele rows and the block beside its SNV eight transcript rows, and under 116.2 both get none. Other records of the measured sets on the same routes: the GRCh37 breakend `gnomAD-SV_v2.1_BND_21_56677` (`SVLEN=14688819`; the warning reads `variant size (14688818)`) gets 0 rows against vep-rs's 381; the GRCh37 complex event `gnomAD-SV_v2.1_CPX_21_5151` 0 against 4; the GRCh38 `gnomAD-SV_v3_CTX_chr21_d895c6b5` 0 (`CTX is not a supported structural variant type`) against vep-rs's one `intergenic_variant` row under the allele `CTX`; the GRCh38 breakend `gnomAD-SV_v3_BND_chr21_f6d2016e` (`SVLEN=39158730`; the warning reads `variant size (39158729)`) 0 against 6,967; and `gnomAD-SV_v3_CPX_chr21_896ebdf0`, a `<CPX>` spanning 39,354,624 bases, is on two routes at once and gets both warnings, the type's and the size's, against vep-rs's 7,123 rows.

#### 5.7 How the comparator treats this class

The two comparators recognise a dropped record from what each reads, and both set vep-rs's rows for it aside.

The structural-variant comparator, `filter_reference_skipped_records` in `scripts/validation/compare_sv_concordance.py`, reads the input VCF and names the records from it with a port of Ensembl VEP's own routes: `get_SO_term` over the ALT list and `SVTYPE`, the incomplete-deletion test, the size test over the span `get_end` derives (`SVLEN` first, else `END`, a breakend's `END` discarded), and the coordinate checks of `validate_vf`. Each such record is confirmed dropped when Ensembl VEP wrote no row under the record's own Location start and an allele string the record produces: the Sequence Ontology class of its ALT and of its ALT list, the class term vep-rs prints for an ALT whose letters name a class but whose spelling has no term, the raw type of an ALT whose letters name none (`CPX`), a breakend's bracket, own-end and `CHR2`/`END2` mate strings, a sequence ALT as written. Exactly that record's vep-rs rows, under those alleles, are set aside. A deletion whose start exceeds its end is written by vep-rs at the pair `POS-(POS+1)`, printed low to high and starting one base before the record's own start, and the comparator locates its rows there. A kept record sharing the start, or the whole Location string, with a dropped one neither hides it nor joins it (gnomAD-SV places a `<CPX>` at the Location of a kept `<DEL>`; a `<NON_REF>` block can start where a kept record starts): two records at one start producing the same allele are told apart by the end of the Location. The records and the outputs are the authority, not the warnings file, which is absent on a cached reference and truncates its record excerpt. A disagreement between them stops the comparison with an `ERROR:` naming the record: a record on a route whose own rows Ensembl VEP wrote (a reference run under `--dont_skip` or under another `--max_sv_size`, declared with `--reference-max-sv-size`, `-1` for no cap); vep-rs structural rows at a Location Ensembl VEP wrote nothing for that no record on a route explains; or a dropped and a kept record at one start whose alleles and ends coincide.

The SNP/indel comparator, `scripts/concordance/compare_vep_outputs.py`, reads the two outputs and no input VCF, so it applies an output-only rule: Ensembl VEP 116.2 writes at least one row for every record it keeps, an `intergenic_variant` row when no transcript is in range, so a Location at which vep-rs wrote rows and Ensembl VEP wrote nothing, under any allele, is a record Ensembl VEP dropped. At least one vep-rs row at the Location must carry a structural allele (every route to the mark is structural), and `tandem_repeat` does not admit a Location, because Ensembl VEP keeps a `<CNV:TR>` record and writes it at another Location (class 12); once admitted, every vep-rs row at the Location is set aside, a sequence allele beside the symbolic one of a mixed ALT list included, because Ensembl VEP skips the line whole. The SNP/indel datasets carry no such record.

Counts (`excluded_reference_skipped_rust`, per file in `excluded_reference_skipped_by_file`, the records in `reference_skipped_records`):

- GRCh37: 20,524 vep-rs tuples on 389 records. gnomAD-SV v2.1 chr21: 19,135 tuples on 89 records, its 60 `<CPX>` records and 29 `<BND>` records whose `SVLEN` exceeds the cap. `10_special_alleles`: 200 tuples on its 200 `<NON_REF>` blocks, one `intergenic_variant` row each. `12_complex_imprecise`: 1,189 tuples on its 100 `<CPX>` records.
- GRCh38: 278,470 vep-rs tuples on 594 records per file; pooled, 278,309, because 161 tuples of records the gnomAD-SV and 1000 Genomes files carry at identical coordinates are counted once. gnomAD-SV v4.1 chr21: 272,627 tuples on 245 records, its 169 `<CPX>` records, its two `<CTX>` records, and 72 `<BND>` records and two mobile-element deletions whose `SVLEN` exceeds the cap. `1kg_sv_chr21`: 662 tuples on its 49 `<CPX>` records. `10_special_alleles`: 200 tuples on 200 blocks. `12_complex_imprecise`: 4,981 tuples on its 100 `<CPX>` records.

The class is fully masked on the vep-rs side. No record of the measured sets takes the fifth route.

## Defect P: the tests after the chromosome-aware overlap compare chromosome names the overlap did not resolve

**Defect D, retired.** Up to release 115.2 Ensembl VEP matched loaded transcripts to the records of a batch through one interval tree keyed on coordinates alone, with no chromosome, so a chromosome 21 record could be handed a chromosome X transcript whose coordinates overlapped its numbers and lose its `intergenic_variant` row to it (the paper's defect D, classes 9 and 10); release 116.1 keys the overlap on the chromosome as well (`InputBuffer.pm` lines 286-345, `_chr_matches_for_vf_overlap` at lines 348-357: a record is handed to a transcript only when the record's chromosome, or one of its breakend mates', is the transcript's, names resolved through the synonym table), and neither class has a row against release 116.2.

**The defect in one sentence.** After the chromosome-aware overlap has handed a record to a transcript, two further tests compare chromosome names the overlap did not resolve the same way: the slice test names a record that arrives without a slice after the first transcript that reaches it, and the mate-allele test compares the mate's chromosome as written, so a bracket breakend loses the rows of its own end when a mate-chromosome transcript reaches it first, and the rows of its mate when the mate's chromosome is spelled in a way that test does not resolve, `Chr21` for `21`.

**Why this is a defect.** A record's chromosome is a fact of its own line, and a mate's chromosome is a fact of its allele; neither depends on which transcript the annotation pass happened to reach first, or on the spelling the file uses for a chromosome the cache knows under another name. Ensembl VEP's own overlap admits the transcript, and a later test then refuses the row for a reason that is not in the record. The result changes with the order of the input: the same record keeps its own-end rows in one file and loses them in another.

**What vep-rs does instead.** vep-rs annotates each breakend of a record on its own chromosome, resolved once from the record and the allele: the own end against the transcripts of CHROM within 5,000 bases of the base after POS, the mate against the transcripts of the mate's chromosome within 5,000 bases of the mate position, whatever the order of the input and whichever name the file uses for the chromosome.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):
- `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, line 122: `$vf->{slice} ||= $slice ||= $tr->{slice}`, inside the loop over the records a transcript reaches: a record without a slice takes the slice of the first transcript that reaches it, and in offline mode a record arrives without one.
- The same module, lines 125-143: the slice test. A `chromosome_breakpoint` record passes when the transcript's chromosome equals the record's slice name (line 128) or any of its breakends' `seq_region_name` with a leading `chr` stripped (lines 130-136); any other record passes only on its slice name (line 141); `next if $seq_mismatch` at line 143 drops the pair.
- `AnnotationSource.pm`, lines 235-258: the regions of a buffer are loaded in the order its records first ask for them, each record's own region before its breakend mates', a region listed once at its first request; `AnnotationType/Transcript.pm` lines 277-300 (`merge_features`) remove duplicates and keep that order, so the transcripts are visited region by region in it. A record whose mate lies within reach of a transcript in a region an earlier record asked for, before the record's own region was asked for, is therefore reached by that transcript before any of its own chromosome's.
- `InputBuffer.pm`, lines 191-207 and 223-234: a buffer holds the records of one chromosome, so the records that share a buffer with a breakend are the ones on its own chromosome before it in the input.
- `InputBuffer.pm`, lines 348-357: `_chr_matches_for_vf_overlap` resolves a chromosome name through `get_source_chr_name`, the synonym table, so a mate written `NC_000002.12` reaches the transcripts of chromosome 2; the slice test at `AnnotationType/Transcript.pm` line 134 compares the raw name, `chr` stripped, and `NC_000002.12` is not `2`. A mate written `Chr21` passes the slice test (the strip is case-insensitive) and fails the mate-allele test below, whose synonym lookup is verbatim and whose fallback is a case-folded equality of the whole name.

(`ensembl-variation`): `StructuralVariationOverlap.pm` lines 86-102 and 144-160: a mate allele is added to a row only when `_close_to_feature` passes, which compares the mate's `seq_region_name` with the transcript's through `_compare_seq_region_names` (`Utils/VariationEffect.pm` lines 111-124: the synonym table's keys verbatim, else a case-folded equality) and applies the 5,000-base window (`MAX_DISTANCE_FROM_TRANSCRIPT`, line 60). `StructuralVariationFeature.pm` lines 1300-1304: a breakend's slice is named after the chromosome exactly as the allele writes it.

**What correcting it would change.** If the slice test and the mate-allele test compared the record's own chromosome and each breakend's, resolved as the overlap resolves them, classes 23 and 24 would vanish: a bracket breakend's own-end rows and mate rows would be written whatever the order of the input and the spelling of the mate's chromosome.

**The classes it produces.**
- Class 23: the own-end rows of a bracket breakend whose mate is on another chromosome, lost when a mate-chromosome transcript names its slice first (the synthetic breakend set: 202 vep-rs tuples on 53 records in GRCh37, 1,012 on 71 in GRCh38). Set aside: vep-rs's rows, in the structural-variant comparator under the cache authority; counted and not set aside by the SNP/indel comparator.
- Class 24: the mate rows of a bracket breakend whose mate's chromosome is spelled with a prefix the mate-allele test does not resolve, `Chr21` (none on the measured sets). Set aside the same way.

A third shape of the same test is counted and never set aside, under the name `fasta_named_slice_intergenic`: under `--check_ref` or `--lookup_ref` with a FASTA whose sequence names differ from the cache's (`chr1` against `1`), every record's slice is named after the FASTA (`Parser.pm` line 749) and every transcript fails the test at line 141, so Ensembl VEP writes `intergenic_variant` where vep-rs writes transcript rows. The condition is the run's configuration, not a property of any record, and it reaches every record of the run; both comparators count it, and a run that shows the count is re-run with matching sequence names rather than scored. It is 0 on every measured set.

### 23. A bracket breakend's own-end rows lost to the first transcript that reaches it

**In short.** `synth_bnd_inter_0062_B` is a chromosome 1 breakend of the synthetic set whose mate lies in SOD1 on chromosome 21. Run alone, both engines write the same seven rows for it: one under its own-end allele, for the transcript its own breakpoint lands in, and six under the mate allele. In the whole file Ensembl VEP writes the six mate rows and nothing under the own-end allele, because the file's first record, a chromosome 1 breakend whose mate lies 1,016 bases from this one's, asked for the chromosome 21 region first: a SOD1 transcript reaches the record before any chromosome 1 transcript does, names its slice `21`, and the chromosome 1 transcript then fails the slice test. vep-rs writes the seven rows in every context. The transcripts the own end reaches are a fact of the record's own line, and Ensembl VEP writes the row itself when the record is alone, with its mate, under `--buffer_size 1`, or before the other record, so vep-rs's row is the correct one.

#### 23.1 What a reader needs to know first

A bracket breakend has two coordinates, its own (POS on CHROM) and its mate's (the chromosome and position inside the allele's brackets), and Ensembl VEP writes rows under two alleles for it: the own end, printed as the reference base and a dot (`A.` or `.A`, the dot on the side the mate joins), and the mate, printed as the bracket string. The own-end rows are the transcripts within 5,000 bases of the base after POS on CHROM; the mate rows are the transcripts within 5,000 bases of the mate position on the mate's chromosome (class 11 describes what those rows should say).

Ensembl VEP annotates in batches of at most `--buffer_size` records (default 5,000) that never span two chromosomes, split into chunks for the `--fork` workers, and for each chunk loads the transcript-cache regions (1 Mb bins) its records overlap, a record's own region and each mate's, in the order the records ask for them. The transcripts are then visited region by region in that order, and a transcript that reaches a record decides, by a test on chromosome names, whether to annotate it. The term involved is `feature_truncation` (SO:0001906) with the region the breakpoint lands in, or `upstream_gene_variant` and `downstream_gene_variant` for a breakpoint near a transcript (11.1).

#### 23.2 The record

`synth_bnd_inter_0062_B`, the second record of the synthetic breakend set `tests/sv_validation/grch37/09_breakends.vcf.gz` (GRCh37), a chromosome 1 breakend whose mate lies in intron 1 of SOD1 on chromosome 21, with the file's first record, `synth_bnd_inter_0066_B` (the record of class 11), whose mate lies 1,016 bases from it in the same intron:

```
1	14374333	synth_bnd_inter_0066_B	A	]21:33034355]A	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0066_A
1	16602616	synth_bnd_inter_0062_B	A	A]21:33033339]	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0062_A
```

The record's own breakpoint, 1:16602617, lies in an intron of ENST00000375592 (ENSG00000037637, reverse strand), the one chromosome 1 transcript within 5,000 bases of it; its mate, 21:33033339, lies within 5,000 bases of six chromosome 21 transcripts, the SOD1 transcripts and their neighbours of 11.2.

The GRCh38 set (`tests/sv_validation/grch38/09_breakends.vcf.gz`) shows the same shape on `synth_bnd_inter_0007_B`, whose own breakpoint 1:57985330 lies in an intron of the non-coding transcript ENST00000485760 (ENSG00000173406) and whose mate 21:25898426 lies in APP, after the file's first record, `synth_bnd_inter_0015_B`, whose mate lies 43 kb away in the same gene:

```
1	10305888	synth_bnd_inter_0015_B	A	]21:25941574]A	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0015_A
1	57985329	synth_bnd_inter_0007_B	A	[21:25898426[A	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0007_A
```

#### 23.3 What each engine writes

Alone, both engines write seven rows for the GRCh37 record: one under its own-end allele `A.` and six under the mate allele `A]21:33033339]` (the mate rows differ between the engines by class 11, as the rows of 11.3 do). In the whole file Ensembl VEP writes the six mate rows, unchanged, and nothing under `A.`; vep-rs writes the same seven rows as alone:

```
# Ensembl VEP, whole file (2 of its 6 rows, all under the mate allele)
synth_bnd_inter_0062_B	1:16602617	A]21:33033339]	ENSG00000142168	ENST00000270142	Transcript	feature_truncation	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
synth_bnd_inter_0062_B	1:16602617	A]21:33033339]	ENSG00000234509	ENST00000449339	Transcript	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
# Ensembl VEP, whole file
(no row under A.)
# vep-rs, whole file (1 of 7 rows: the own-end row, identical to the row both engines write for the record alone)
synth_bnd_inter_0062_B	1:16602617	A.	ENSG00000037637	ENST00000375592	Transcript	feature_truncation,intron_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1
```

The GRCh38 record has 42 rows on each side when run alone, one of them under its own-end allele `.A`; in the whole file Ensembl VEP writes the 41 mate rows and vep-rs all 42:

```
# vep-rs, whole file (1 of 42 rows; Ensembl VEP writes the other 41)
synth_bnd_inter_0007_B	1:57985330	.A	ENSG00000173406	ENST00000485760	Transcript	feature_truncation,intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1
```

#### 23.4 Why Ensembl VEP writes that

A buffer holds records of one chromosome (`InputBuffer.pm` lines 191-207 and 223-234): the GRCh37 file's seventeen chromosome 1 records form one buffer, with `synth_bnd_inter_0066_B` first and `synth_bnd_inter_0062_B` second. The regions of a buffer are loaded in the order its records first ask for them, each record's own region before its mate's, a region asked for twice listed once, at its first request (`AnnotationSource.pm` lines 235-258): the first record asks for its own region, chromosome 1 at 14 Mb, then its mate's, chromosome 21 at 33 Mb; the second record's own region, chromosome 1 at 16 Mb, comes third, and its mate's region is the one already listed second.

The transcripts are visited in that order. The chromosome-aware overlap hands the record to the SOD1 transcripts near its mate (through the mate, `InputBuffer.pm` lines 336-342) and to ENST00000375592 near its own breakpoint (through its own coordinates, lines 324-326). The SOD1 transcripts come first, and line 122 of `AnnotationType/Transcript.pm` names the record's slice after the first of them, `21`. From then on the slice test (lines 125-143) passes a transcript only if its chromosome is `21`, the slice's name, or the mate's, also `21`; ENST00000375592 fails it and is dropped at line 143, and the own-end row is never created. When the record's own region comes first in the load order, ENST00000375592 names the slice `1`, the SOD1 transcripts still pass through the mate, and both sets of rows are written. The whole GRCh38 file does the same to `synth_bnd_inter_0007_B`: the first chromosome 1 record's mate region, chromosome 21 at 25 Mb, holds the APP transcripts that reach the record's mate, and it is listed before the record's own region at 57 Mb.

#### 23.5 Why vep-rs's row is the correct one

The transcripts within 5,000 bases of the record's own breakpoint on chromosome 1 are a fact of the record's own line, and Ensembl VEP's own admission rule gives the own-end allele exactly those rows. Whether they are written cannot depend on which record of the file asked for the mate's region first. Ensembl VEP agrees whenever the order is the other way round: alone, with its mate, under `--buffer_size 1`, or placed before the other record, it writes the own-end row, byte for byte the row vep-rs writes in every context (23.6).

#### 23.6 Try it yourself

Annotate the GRCh37 record with Ensembl VEP 116.2 (GRCh37 cache and FASTA, the command of the closing section) in each of these inputs, and count its rows under the own-end allele `A.` and under the mate allele:

| Input | Ensembl VEP rows for `synth_bnd_inter_0062_B`, own end / mate |
| --- | --- |
| the record alone | 1 / 6 |
| the record with its mate record `synth_bnd_inter_0062_A` (`21 33033339 ... A [1:16602616[A`) | 1 / 6 |
| `synth_bnd_inter_0066_B` then the record (the file's first two lines) | 0 / 6 |
| the record then `synth_bnd_inter_0066_B` | 1 / 6 |
| a hand-written `1 14374333 x A ]21:33034355]A` with `SVTYPE=BND` and no MATEID, then the record | 0 / 6 |
| a hand-written chromosome 1 breakend at 14374333 whose mate is far from the record's, `]21:27367670]A`, then the record | 1 / 6 |
| a hand-written chromosome 1 breakend at a higher position, `1 16954684 x A ]21:33034355]A`, then the record (both own breakpoints lie in the same 1 Mb region, which is then listed first) | 1 / 6 |
| a hand-written chromosome 2 breakend, `2 1000000 x A ]21:33034355]A`, then the record (two buffers: a buffer holds one chromosome) | 1 / 6 |
| the SNV `21 33033339 . A C` then the record, or the SNV `1 14374333 . A C` then the record, or the SNV `1 16602616 . C A` then the record | 1 / 6 |
| the mate record `synth_bnd_inter_0066_A` (on chromosome 21) then the record | 1 / 6 |
| the whole `09_breakends.vcf.gz`, `--fork 4 --buffer_size 5000` (the reference setting), or `--fork 1 --buffer_size 5000` | 0 / 6 |
| the whole file at `--buffer_size 1`, `--fork 1` or `--fork 4` | 1 / 6 |

vep-rs writes 1 / 6 for every one of those inputs. The loss needs, in the record's own buffer and before it, a record that asks for the region holding the transcripts within reach of the record's mate before the record's own region is asked for; `--fork` does not change the result and `--buffer_size 1`, which annotates every record alone, removes it. The GRCh38 record behaves the same way: 1 / 41 alone, with its mate, after the SNV `21 25898426 . T A` at its mate's position, or after SNVs at the first record's two breakpoints; 0 / 41 after the file's first record `synth_bnd_inter_0015_B`, after `synth_bnd_inter_0037_B` (`1 11286046 ... A A[21:26072808[`, whose mate's region holds the APP transcripts that reach this record's mate), and in the whole file.

#### 23.7 How the comparator treats this class

`filter_bnd_own_rows_lost` in `scripts/validation/compare_sv_concordance.py` runs under the cache authority (`--vep-rs-cache`, the transcript cache that says which chromosome a transcript is on). It sets aside a vep-rs tuple when all of these hold: the Location is a bracket breakend's whose mate is on another chromosome; Ensembl VEP wrote a Transcript row at the Location under the bracket allele and none under any own-end allele (`N.`, `.N` or `chromosome_breakpoint`); and the tuple's transcript is one the cache files on the Location's chromosome within 5,000 bases of the breakpoint, Ensembl VEP's own admission rule (`StructuralVariationOverlap.pm` lines 144-160). A vep-rs own-end row on a transcript farther than that, or on another chromosome, stays charged.

The SNP/indel comparator, `scripts/concordance/compare_vep_outputs.py`, has no transcript-to-chromosome authority, so it counts the shape under the same name and sets nothing aside: a Location at which Ensembl VEP wrote rows under a bracket allele and none under an own-end allele, and vep-rs wrote own-end rows.

Counts (`excluded_bnd_own_rows_lost_rust`, per file in `excluded_bnd_own_rows_lost_by_file`):

- GRCh37: 202 vep-rs tuples on 53 records, all in `09_breakends.vcf.gz`: the records on chromosomes 1 to 5 whose mates are on chromosome 21 (12 on chromosome 1, 17 on 2, 8 on 3, 5 on 4, 11 on 5); none of the set's chromosome 21 records, and none in any other file.
- GRCh38: 1,012 vep-rs tuples on 71 records, all in `09_breakends.vcf.gz` (21 on chromosome 1, 18 on 2, 12 on 3, 11 on 4, 9 on 5).

The gnomAD-SV files, the 1000 Genomes files and the other synthetic files show none. Fully masked, vep-rs side.

### 24. A bracket breakend's mate rows lost to an unresolved chromosome prefix

**In short.** A hand-written chromosome 21 breakend whose mate, 50 bases away on the same chromosome, is written `N[Chr21:30425700[` reaches its transcripts through its own coordinates and passes the slice test; the test that adds the mate allele then compares `Chr21` with the transcript's `21` through the synonym table, whose keys hold `chr21` and not `Chr21`, and a case-folded equality, which `chr21` against `21` fails, so Ensembl VEP writes the own-end rows and nothing under the mate allele. vep-rs reads `Chr21` as chromosome 21 and writes the mate rows too. The mate's chromosome is a fact of the allele, and Ensembl VEP writes the same mate rows when the allele spells it `chr21` or `21`, so vep-rs's rows are the correct ones.

#### 24.1 What a reader needs to know first

A VCF may name a chromosome by any of the names its reference assembly report lists: `21`, `chr21`, `NC_000021.9`, `CM000683.2`. Ensembl VEP's cache carries a synonym table for them (`--synonyms` adds to it), and its chromosome-aware overlap resolves a mate's name through that table before matching it to a transcript's chromosome (`InputBuffer.pm` lines 348-357). Two later tests compare names again: the slice test of `AnnotationType/Transcript.pm` lines 125-143, which strips a leading `chr` from a breakend's name, case-insensitively, and compares the rest as written, and `_close_to_feature` (`StructuralVariationOverlap.pm` lines 144-160), which accepts a synonym table key verbatim or a case-folded equality of the whole name. The mate rows are written only when every test passes.

#### 24.2 The record

Hand-written, GRCh38: a breakend at 21:30425650 whose mate lies 50 bases away on the same chromosome, written with a capitalised prefix:

```
21	30425650	x	N	N[Chr21:30425700[	.	.	SVTYPE=BND
```

The own breakpoint, 21:30425651, lies in the coding sequence of ENST00000390690 (reverse strand) and 4,569 bases upstream of ENST00000334068; the mate, 50 bases away, reaches the same two transcripts. The GRCh37 counterpart is `21 33034300 x N N[Chr21:33034355[` with `SVTYPE=BND`, the breakpoint of class 11 as its mate, which reaches seven transcripts.

#### 24.3 What each engine writes

```
# Ensembl VEP (both rows)
x	21:30425651	N.	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=4569;STRAND=1
x	21:30425651	N.	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant	276	262	88	-	-	-	IMPACT=HIGH;STRAND=-1
# vep-rs (all four rows)
x	21:30425651	N.	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=4569;STRAND=1
x	21:30425651	N[Chr21:30425700[	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=4569;STRAND=1
x	21:30425651	N.	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant	276	262	88	-	-	-	IMPACT=HIGH;STRAND=-1
x	21:30425651	N[Chr21:30425700[	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant	276	262	88	-	-	-	IMPACT=HIGH;STRAND=-1
```

Ensembl VEP writes the two own-end rows and nothing under the mate allele; vep-rs writes the same two own-end rows and the two mate rows, identical to the own-end rows in every column but Allele because the two breakpoints are 50 bases apart. On GRCh37 the record writes 7 rows against 14 the same way.

#### 24.4 Why Ensembl VEP writes that

The chromosome 21 transcripts reach the record through its own coordinates and pass the slice test on the record's slice, `21`. For each of them the mate allele is added to the row only if `_close_to_feature` passes for the breakend (`StructuralVariationOverlap.pm` lines 90-102), and its name test, `_compare_seq_region_names` (`Utils/VariationEffect.pm` lines 111-124), finds `Chr21` neither among the synonym table's keys, which hold `chr21`, nor equal to `21` once both are lower-cased. The own-end rows are written; the mate allele is never added. With the mate written `chr21` or `21` the same test passes and the mate rows are written.

A mate written with a RefSeq accession takes the other route. On another chromosome, `N[NC_000002.12:321682[`, the overlap resolves the accession (`InputBuffer.pm` line 356) and hands the record to the chromosome 2 transcripts, and the slice test then compares the breakend's `seq_region_name`, the accession as the allele wrote it (`StructuralVariationFeature.pm` lines 1300-1304), `chr` stripped, with the transcript's `2` (`AnnotationType/Transcript.pm` lines 130-136), finds no match and drops the pair at line 143: Ensembl VEP writes no mate row. vep-rs does not resolve the accession and writes no mate row either, so the two outputs agree and nothing is set aside. On the record's own chromosome, `N[NC_000021.9:30425700[`, the transcripts reach the record through its own coordinates, the accession is a key of the synonym table, and Ensembl VEP writes the mate rows; vep-rs writes none, a vep-rs under-annotation that stays charged to vep-rs in both figures.

#### 24.5 Why vep-rs's row is the correct one

The mate's chromosome is a fact of the allele, and `Chr21` names the chromosome the cache holds as `21` and lists as `chr21`. The transcripts within 5,000 bases of the mate position on that chromosome are the mate rows under Ensembl VEP's own admission rule, and Ensembl VEP writes them itself when the allele spells the chromosome `chr21` or `21` (24.6). A prefix that differs from a listed synonym by its case cannot be a reason to withhold them.

#### 24.6 Try it yourself

Annotate the record with each spelling of the mate's chromosome, with Ensembl VEP 116.2 and vep-rs (GRCh38 cache and FASTA, the commands of the closing section; the GRCh37 counterpart of 24.2 in parentheses):

| Mate | Ensembl VEP 116.2 rows | vep-rs rows |
| --- | --- | --- |
| `N[21:30425700[` or `N[chr21:30425700[` | 4 (14): own end and mate | 4 (14) |
| `N[Chr21:30425700[` | 2 (7): own end only | 4 (14) |
| `N[NC_000021.9:30425700[` (GRCh37 `NC_000021.8`) | 4 (14) | 2 (7): own end only, a vep-rs under-annotation charged to vep-rs |
| `N[2:321682[` or `N[chr2:321682[` (GRCh37 `Y:15437971`, `chrY:15437971`) | 13 (15) | 13 (15); the mate rows differ by class 11 |
| `N[NC_000002.12:321682[` (GRCh37 `NC_000024.9:15437971`) | 2 (7): own end only | 2 (7): own end only, so the outputs agree and nothing is set aside |

Only the capitalised prefix produces the class: Ensembl VEP writes the own-end rows alone where vep-rs also writes the mate rows.

#### 24.7 How the comparator treats this class

`filter_bnd_synonym_mate_rows` in `scripts/validation/compare_sv_concordance.py` runs under the cache authority. It sets aside a vep-rs tuple under a bracket allele when the mate name, `chr` stripped, is not a chromosome the cache files, Ensembl VEP wrote rows at the Location and none under the bracket allele, and the tuple's transcript is not on the Location's own chromosome; or when the mate name is the record's own chromosome spelled with a prefix Ensembl VEP does not resolve (`Chr21`), in which case vep-rs's mate rows on the own chromosome are set aside too. The SNP/indel comparator counts the shape and sets nothing aside, a synonym there being any mate name outside the canonical contig set.

Counts (`excluded_bnd_synonym_mate_rows_rust`): 0 on every measured set. Every mate in the gnomAD-SV files is named through `CHR2` as `chr21`, `chrX` or `chrY`, synonyms the cache lists, and the synthetic breakends write their mates with the cache's own names.

## Defect E: a breakend's mate-side row evaluated at the record's own coordinate

**The defect in one sentence.** For a breakend, Ensembl VEP chooses which transcripts get a mate-side row by the mate's chromosome and coordinate, then evaluates almost every consequence test on that row with the record's own coordinate, which may lie on another chromosome.

**Why this is a defect.** The row is written for the mate breakpoint: its transcript was chosen because the mate lands in or near it. Asking whether the record's chromosome 1 position lies inside a chromosome 21 transcript's intron answers a question about the wrong breakpoint. Only `feature_truncation` looks at the mate, so the row comes out as bare `feature_truncation` or, when even that does not fire, falls through to `intergenic_variant` on a row that names a transcript.

**What vep-rs does instead.** vep-rs evaluates the mate breakpoint where it is: inside the transcript it writes `feature_truncation` plus the exon, intron, UTR or coding region the breakpoint lands in; within 5,000 bases of the transcript, `upstream_gene_variant` or `downstream_gene_variant`.

**Where it is in Ensembl VEP's source** (`ensembl-variation`):
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 257 and 273: the record's own feature is passed to every predicate.
- `Utils/VariationEffect.pm`, lines 350-358: `feature_truncation` alone tests `within_feature` on the breakend; lines 444-458: `upstream` and `downstream` read the record; lines 629-634: `within_intron` reads the record's differing region, which `BaseTranscriptVariationAllele.pm` line 118 derives from the record's start.
- `StructuralVariationOverlap.pm`, lines 86-102 and 144-160: the mate allele is created from the mate's coordinate and chromosome.

**What correcting it would change.** If every predicate were evaluated at the breakpoint the row is for, class 11 would vanish.

**The classes it produces.**
- Class 11: mate-side rows of point breakends whose mate lies on another chromosome or far from the record (the synthetic `09_breakends` and `13_vcf45_features` sets and gnomAD-SV v4.1 chr21: 1,407 and 10,518 pairs). Set aside: both members of each pair.

### 11. Breakend mate rows read from the other chromosome's coordinate

**In short.** The synthetic breakend `synth_bnd_inter_0066_B` sits on chromosome 1 and names its mate at 21:33034355, inside intron 1 of SOD1 and 2,542 bases upstream of a lincRNA. Ensembl VEP picks the transcripts for the mate-allele rows by the mate's coordinate but evaluates the rows with the record's own chromosome 1 coordinate, so it writes `feature_truncation` alone for SOD1 and `intergenic_variant` for the lincRNA. vep-rs writes `feature_truncation,intron_variant` and `upstream_gene_variant` with `DISTANCE=2542`. vep-rs's rows state what the mate breakpoint does where it is, which is what Ensembl VEP itself writes for the same breakpoint when it is a record's own coordinate, so they are the correct ones.

#### 11.1 What a reader needs to know first

A breakend record has two coordinates: its own (POS on CHROM) and its mate's (the chromosome and position inside the allele's brackets).

Ensembl VEP represents a breakend in one of two ways:

- A `<BND>` whose mate is on the same chromosome and given in `INFO/CHR2` and `END2` becomes a span from POS+1 to the mate (class 5 shows one). Every allele of a span record is evaluated over the span.
- A breakend written in bracket notation, or a `<BND>` whose mate is on another chromosome, becomes a point record whose Location is the single base after POS.

This class is about the mate-allele rows of point records. For a transcript within 5,000 bases of the mate, Ensembl VEP writes a mate-allele row (`StructuralVariationOverlap.pm` lines 86-102 and 144-160). What that row should say is what the mate breakpoint does to that transcript:

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
- `StructuralVariationOverlap.pm`, lines 86-102 and 144-160: the mate allele is created from the mate's coordinate and chromosome.

#### 11.5 Why vep-rs's row is the correct one

The mate breakpoint 21:33034355 lies in intron 1 of SOD1 (33032155-33036102) and 2,542 bases upstream of the lincRNA's end at 33031813. Those are facts about the record's own coordinates, and they are what vep-rs writes.

Ensembl VEP's row for SOD1 says the transcript is truncated but names no region. Its row for the lincRNA calls a position 2,542 bases from the transcript intergenic, on a line that names the transcript.

Ensembl VEP itself annotates the same breakpoint correctly when it is the record's own coordinate (11.3, the second record of the pair). Its mate row changes when the record's chromosome 1 position is moved while the mate stays fixed (11.6), which no correct reading of the mate row could do.

#### 11.6 Try it yourself

Annotate with Ensembl VEP 116.2, GRCh37 cache and FASTA, using the command of the closing section. Keep the mate at `21:33034355` and move the record's own position:

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
- GRCh38: 10,518 pairs (`09_breakends.vcf.gz` 6,333; `13_vcf45_features.vcf.gz` 27; `gnomad_sv_v4.1_chr21.vcf.gz` 4,158).

Rows like the coincidence record, where the two derivations agree, are concordant and not in the class. The class is fully masked, both sides.

vep-rs releases before 0.3.0 evaluated no position on the mate side and so reproduced Ensembl VEP's fall-through rows; from 0.3.0 vep-rs evaluates the mate breakpoint where it is, and these pairs are what that costs the raw figure, which sets nothing aside.

## A representation difference, not a defect: `<CNV:TR>`

**In one sentence.** A `<CNV:TR>` record describes a change in the copy number of a tandem repeat without spelling the alternate allele out, and the two engines read it in two different, equally faithful ways.

**The two readings.** Ensembl VEP materialises the allele: it repeats the unit the stated number of times and annotates the difference as an ordinary insertion or deletion of bases. vep-rs keeps the allele symbolic, as a copy-number change of the run, and annotates the run's span for its direction. Ensembl VEP itself gives the symbolic reading when the expanded allele would exceed `--max_sv_size`.

**Why neither is wrong.** Ensembl VEP's row is a faithful annotation of a real molecule: a run of 121 units where the reference has nine is, base for base, an insertion. vep-rs's row keeps the record in the terms its author chose: a copy-number change of a repeat, located on the repeat, with the region term and the direction term. Nothing in either row is false, so the comparator pairs the two readings rather than ranking them.

**Where the readings part** (`ensembl-vep`): `Bio/EnsEMBL/VEP/Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`, lines 382-431, called from `create_StructuralVariationFeatures` at line 520.

**The class.**
- Class 12: the synthetic tandem-repeat set on both assemblies, 1,368 and 3,751 pairs. Set aside: both members of each pair.

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

Both engines follow Ensembl's rule for a span allele (a deletion, duplication, insertion or repeat): `feature_elongation` and `feature_truncation` are written only when the span touches the transcript's exonic sequence. The rule is in `ensembl-variation`, `Utils/VariationEffect.pm` lines 337-348 and 350-373, the `within_cdna` tests at lines 344 and 362. So a run inside an intron gets the region term alone, and a run inside an exon gets the region term plus the direction term. A breakend's `feature_truncation` is a different test, the breakpoint falling anywhere inside the transcript, intron included (lines 357-358; class 11).

#### 12.2 The record

Both records are from the synthetic tandem-repeat set `tests/sv_validation/grch38/07_cnv_repeat.vcf.gz` (GRCh38; the GRCh37 set is `tests/sv_validation/grch37/07_cnv_repeat.vcf.gz`). Every `<CNV:TR>` record of those sets sits on a tandem repeat the chromosome 21 reference actually carries. The first is a gain in an intron:

```
21	25896418	synth_cnv_tr_0001	T	<CNV:TR>	.	PASS	SVTYPE=CNV;END=25896436;SVLEN=18;RN=1;RUS=CA;RUL=2;RUC=121;RB=242
```

The reference at 21:25896410-25896445 reads `AGTAACACT` `CACACACACACACACACA` `AAACAAAAC` (9 + 18 + 9 bases). The run is the 18 bases 25896419-25896436 (nine CA units), and the alternate has 121 units, a gain of 112 units or 224 bases.

The transcript is ENST00000346798 (APP, reverse strand, canonical, 25880550-26170770); the run lies in its intron 16 (25891869-25897572). The record has 41 rows on each side, one per APP transcript within reach of the run.

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

Ensembl VEP builds the alternate allele as a literal sequence (it materialises the allele). It reads the reference bases of the run and builds the alternate as `RUS` repeated `RUC` times (or `RB` divided by the unit's length, in whole units, when `RUC` is absent). It then hands the pair to its sequence-variant machinery, which trims the common prefix. What is left is a gain as an insertion of whole units after the run's last base, or a loss as a deletion of whole units ending at the run's last base.

The materialisation is skipped for the whole record when any repeat sequence of any of its alleles is longer than `--max_sv_size` bases (the same 10,000,000 default; under `-1` it is always skipped), and every allele of the record is then annotated symbolically. vep-rs reads the direction of change under the same limit, so above it the two engines agree on the symbolic gain. The materialisation also needs the run's reference bases: a run without `--fasta` (`--offline` with no FASTA, the protocol of the golden corpora; the measured structural-variant sets ran with the FASTA, so their reference rows are the expanded alleles) fetches no slice, builds no literal allele, and every `<CNV:TR>` record is annotated as a symbolic `tandem_repeat` over the run, which `copy_number_gain` reads as a gain, so a loss is labelled `feature_elongation` with the region terms (lines 516-527 of the same file; `Utils/VariationEffect.pm` lines 246-256). vep-rs reads the direction from the INFO fields with or without a FASTA, so against a FASTA-less reference a loss is `feature_truncation` with the deletion's terms where the reference writes `feature_elongation`; a gain keeps the run as its span on both sides and agrees.

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

Ensembl VEP 116.2, GRCh38 cache and FASTA, the command of the closing section. Lower the expansion cap below the alternate allele's length and Ensembl VEP annotates the record symbolically. Two settings do that for the two records:

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

Both rows are byte for byte the vep-rs rows of 12.3. Running Ensembl VEP without `--fasta` gives the same two symbolic rows, with the warning `could not fetch sequence for tandem repeats; consequence calculation for tandem repeats is less precise when using --offline without --fasta`. At the default `--max_sv_size` of 10,000,000 and with the FASTA every `<CNV:TR>` record of the synthetic sets is expanded.

#### 12.7 How the comparator treats this class

`filter_cnv_tr_expansion_divergences` in `scripts/validation/compare_sv_concordance.py` pairs a vep-rs tuple whose Allele is `tandem_repeat` with the Ensembl VEP tuple on the same chromosome, Feature and Feature_type that is its literal reading (`_cnv_tr_pair_kind`):

- for a gain, an insertion (Allele all bases) located between the run's last base and the next;
- for a loss, a deletion (Allele `-`) ending at the run's last base;
- for a literal reading the expansion does not trim to a unit indel, any other sequence allele (`N` bases included) whose span touches the run or the base after it: a multi-allelic `<CNV:TR>` record, whose allele string is never prefix-trimmed, a unit that is not the reference's repeat, a `.` unit read as `N`, or an `SVLEN` that disagrees with the unit count, each a substitution across the run (`AGCAGCAGCAGC protein_altering_variant`, `AGCAGC inframe_deletion`, `NNNN frameshift_variant`);
- for an unexpanded reading, the reference's own `tandem_repeat` tuple at the run, written without a FASTA (12.4).

The pair is set aside, both members, when the two consequence sets agree once the symbolic reading's direction term, where present (`feature_elongation` on a gain, `feature_truncation` on a loss, either on a literal reading), is disregarded. One further correspondence holds inside a coding sequence: the literal reading's coding effect (`inframe_insertion`, `inframe_deletion` or `frameshift_variant` for a unit indel; any coding term, `protein_altering_variant` or `missense_variant` included, for a literal substitution) against the symbolic reading's `coding_sequence_variant`; an unexpanded pair may differ in the direction and coding terms alone (`_is_cnv_tr_representation_pair`). Any other difference, a start or splice term one side carries and the other lacks, stays charged.

The pairing reads every shape the expansion can take, the ones below included; the run at 30425588-30425596 is three copies of `AGC`, the run at 30425756-30425763 four copies of `AG`, the run at 30430410-30430417 four copies of `CT`, each with `--fasta` unless stated:

- a unit that is not DNA, which the reference repeats verbatim across the run and vep-rs keeps symbolic: `RUS=123;RUC=4` gives `21:30425588-30425596 123123123123 ENSG00000240432 ENST00000390690 Transcript coding_sequence_variant`, `RUS=R;RUC=10` `RRRRRRRRRR` with `frameshift_variant`, `RUS=*;RUC=9` `*********` with `coding_sequence_variant`, `RUS=X;RUC=10` `XXXXXXXXXX` with `frameshift_variant`, each against vep-rs's `tandem_repeat feature_elongation,coding_sequence_variant` at the run;
- a loss whose literal reading is a substitution with a coding term other than a deletion's: `RN=3;RUS=AG,AG,C;RUC=1,1,1` gives `21:30425760-30425763 C ENSG00000240432 ENST00000390690 Transcript protein_altering_variant`, `RUS=.;RB=5` gives `21:30425756-30425763 NNNNN ENSG00000240432 ENST00000390690 Transcript coding_sequence_variant`, `SVLEN=30_425_949;RUS=AGC;RUC=4` (an `SVLEN` Perl reads as 30) gives `21:30425597-30425617 AGC ENSG00000240432 ENST00000390690 Transcript protein_altering_variant`, and `SVLEN=4000;RUS=CT;RUC=5` on the run at 30430410 gives `21:30430418-30434408 C ENSG00000186971 ENST00000334068 Transcript stop_retained_variant,3_prime_UTR_variant` against vep-rs's `21:30430410-30434409 tandem_repeat` with `stop_lost,feature_truncation,coding_sequence_variant,3_prime_UTR_variant`, the literal reading's stop term against the symbolic reading's;
- a multi-allelic record whose reference rows collapse to one tuple, without a FASTA: `<CNV:TR>,<CNV:TR>,<CNV:TR>` with `RN=1,1,1;RUS=AGC,AGC,AGC;RUC=4,2,3` is one tuple in the reference, `21:30425588-30425596 tandem_repeat ENSG00000240432 ENST00000390690 Transcript feature_elongation,coding_sequence_variant` for all three alleles (every allele a gain without a FASTA, 12.4), and a gain tuple and a loss tuple (`feature_truncation,inframe_deletion`) in vep-rs, both paired with it.

An expansion the reference drops without a row is not a pair but class 19 (19.7): an unchanged expansion over several repeat sequences (`RN=3;RUS=AGC,AGC,AGC;RUC=1,1,1`, `RN=2;RUS=AGC,AG;RUC=3,0`), a unit of no bases under `RUC` (`RUS=,AGC;RUC=4,4`) and a count in fullwidth digits (`RUC=４`, which Perl numifies to 0).

The class is fully masked:

- GRCh37 set: 1,368 of 1,368 pairs;
- GRCh38 set: 3,751 of 3,751 pairs.

The counters are `excluded_cnvtr_rust` and `excluded_cnvtr_perl`, against the class sizes `cnvtr_swap_pairs_total_rust` and `cnvtr_swap_pairs_total_perl`.

## Defect F: a breakpoint record without a mate gives its reference allele a row

**The defect in one sentence.** A breakpoint record that names no mate (`<BND>` with no `CHR2`, or with `CHR2` but neither `END2` nor `END`) is held by Ensembl VEP as the allele string `N/<BND>`, which parses to no breakend, so the overlap object builds one allele per slash-separated part of the string, the reference base included, and the reference allele gets a consequence row of its own.

**Why this is a defect.** A reference allele has no variant effect; a consequence row describes what an alternate allele does to a feature. The second row restates the first under an allele number that points at the reference base, and nothing in it is information about the variant.

**What vep-rs does instead.** vep-rs holds the record as one single breakend, prints the class term `chromosome_breakpoint` and the terms of the span or the point once, and takes no `feature_truncation` for it, as Ensembl VEP's own ALT row does not.

**Where it is in Ensembl VEP's source** (`ensembl-vep` and `ensembl-variation`):
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 491-497: a bracket mate is built only under `defined $info->{CHR2}` with `END2 || END`; line 499 prepends the reference base to the allele string unless it opens or closes with `.`, so the record's string is `N/<BND>`.
- `ensembl-variation`, `StructuralVariationFeature.pm`, `_parse_breakends`, lines 1274-1283: finds no bracket in the string, so the record has no breakends.
- `ensembl-variation`, `StructuralVariationOverlap.pm`, lines 64-77: with no breakends, `split(/\//, $vf->{allele_string})` gives one allele per part, numbered from 1: the reference base `N` is allele 1 and `<BND>` allele 2.
- `ensembl-vep`, `OutputFactory.pm`, line 2029: `$vfoa->{breakend}->{string} || $svf->class_SO_term` prints `chromosome_breakpoint` for both alleles, neither of which carries a breakend string.

**What correcting it would change.** The reference-allele row would disappear and every row of the record would be written once, as vep-rs writes it; class 13 would vanish.

### 13. The reference-allele row of a bare breakpoint

#### 13.1 What a reader needs to know first

A VCF breakend that names no mate (`<BND>` with no `CHR2`, or with `CHR2` but neither `END2` nor `END`, or with a `POS2` Ensembl VEP never reads) is held by Ensembl VEP as the allele string `N/<BND>`: the reference base, a slash, the symbolic ALT. The string parses to no breakend, so the overlap object builds its alleles by splitting the string on `/`, one allele per part, and the first part is the reference base. Both alleles print the class term `chromosome_breakpoint`, so every row of the record appears twice, distinguished only by `ALLELE_NUM`.

#### 13.2 The record

```
21	30425650	bare_bnd	N	<BND>	.	.	SVTYPE=BND
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

#### 13.3 What each engine writes

Ensembl VEP 116.2, four rows for two transcripts:

```
bare_bnd	21:30425651	chromosome_breakpoint	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4569;STRAND=1
bare_bnd	21:30425651	chromosome_breakpoint	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	ALLELE_NUM=2;IMPACT=MODIFIER;DISTANCE=4569;STRAND=1
bare_bnd	21:30425651	chromosome_breakpoint	ENSG00000240432	ENST00000390690	Transcript	coding_sequence_variant	ALLELE_NUM=1;IMPACT=MODIFIER;STRAND=-1
bare_bnd	21:30425651	chromosome_breakpoint	ENSG00000240432	ENST00000390690	Transcript	coding_sequence_variant	ALLELE_NUM=2;IMPACT=MODIFIER;STRAND=-1
```

vep-rs, two rows:

```
bare_bnd	21:30425651	chromosome_breakpoint	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
bare_bnd	21:30425651	chromosome_breakpoint	ENSG00000240432	ENST00000390690	Transcript	coding_sequence_variant
```

The same doubling appears for a `<BND>` with `SVLEN` and no `CHR2` (over its span), for `CHR2` with `POS2` alone, for `CHR2` with an `END2` that is not a coordinate, for a `<BND>` without `SVTYPE`, and for a `<BND>` with `MATEID`; a native single breakend (`.N`, `N.`) is held without the reference base and writes each row once, and a bracket record's own breakend and mate are two different alleles, written once each.

The rows carry no `feature_truncation` under `--no_stats`, the flag the reference protocol runs with: `OutputFactory.pm` line 2029 reads `$vfoa->{breakend}->{string}` before the terms are computed at line 2034, which autovivifies an empty breakend whose chromosome matches no feature's, so the predicate's live arm (`VariationEffect.pm` lines 357-358) fails. With statistics on, `Stats.pm` line 241 computes the terms first, `within_feature` falls back to the whole feature, and the same rows carry `feature_truncation`. vep-rs follows the reference protocol's rows.

#### 13.4 Why Ensembl VEP writes that

`Parser/VCF.pm` builds a bracket mate only under `defined $info->{CHR2}` with `END2 || END` (lines 491-497) and prepends the reference base to the allele string unless it opens or closes with `.` (line 499), so the record's string is `N/<BND>`. `StructuralVariationFeature.pm` `_parse_breakends` (lines 1274-1283) finds no bracket in it. `StructuralVariationOverlap.pm` lines 64-77 then run `if (!@$breakends) { foreach (split(/\//, $vf->{allele_string})) { add_StructuralVariationOverlapAllele(...) } }`, numbering the parts from 1: the reference base `N` is allele 1 and `<BND>` allele 2. Neither allele carries a breakend string, so `OutputFactory.pm` line 2029 (`$vfoa->{breakend}->{string} || $svf->class_SO_term`) prints `chromosome_breakpoint` for both.

#### 13.5 Why vep-rs's row is the correct one

A reference allele has no variant effect; a consequence row describes what an alternate allele does to a feature. The second row restates the first under a number that points at the reference base, and nothing in it is information about the variant. vep-rs holds the record as one single breakend, prints the class term and the terms of the span or the point once, and takes no `feature_truncation` for it, as Ensembl VEP's own ALT row does not.

#### 13.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tbare_bnd\tN\t<BND>\t.\t.\tSVTYPE=BND\n' > bare.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/bare.vcf -o /work/bare.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i bare.vcf -o bare.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

#### 13.7 How the comparator treats this class

The class is `bnd_bare_reference_allele_row` in the 116.2 registry of both comparators, and both set the reference-allele rows aside. Each comparator compares unique tuples of (Location, Allele, Feature, Feature_type, Consequence), and the reference-allele row is identical in every column to the ALT row, so the set-aside happens where the tuples are built: the duplicate row is dropped before the raw and the adjusted figure, neither of which ever carries a tuple for it, and each comparator reads the reference's raw rows to count the rows it set aside (`bnd_bare_reference_allele_rows` in the structural-variant report, `bnd_bare_reference_allele_row` with its `excluded` count in the SNP/indel report, both per file). The class therefore moves no tuple between the raw and the adjusted figure; it names and counts what the construction excluded. None of the measured sets carries the shape: every gnomAD-SV `<BND>` on chromosome 21 names its mate through `CHR2` with `END` or `END2`, and the synthetic sets write their breakends in bracket and dot notation.

## Defect G: a coordinate past the 32-bit bound turns the buffer intergenic

**The defect in one sentence.** Ensembl VEP adds the `--distance` flank to every coordinate of an input buffer when it chooses the cache regions and the buffer's bounds, and when any coordinate plus the flank reaches 2^31 (a breakend's mate position, a span's end or a record's own position) every record of the buffer is written as `intergenic_variant`.

**Why this is a defect.** A record's annotation depends on its own coordinates and the transcripts it overlaps, not on the other records read with it. Here the record that supplied the coordinate loses its transcript rows under all its alleles, and so does every other record loaded with it, whatever transcripts they overlap, while the output does not identify which records were damaged.

**What vep-rs does instead.** vep-rs holds every coordinate in 64 bits, writes the breakend's own-end rows against the transcripts at its own position, writes nothing under the mate allele when no transcript lies at the mate, annotates a span to 2,147,483,647 against every transcript it covers, and annotates the other records of the batch as if the record were absent, which is what Ensembl VEP itself writes for the same records at any coordinate below the bound.

**Where it is in Ensembl VEP's source** (`ensembl-vep`, `Bio/EnsEMBL/VEP/AnnotationSource.pm`):
- line 190: `get_regions_from_coords` adds the `--distance` flank to each coordinate it is given when it chooses the cache regions and the buffer's bounds.
- lines 251-254: `get_all_regions_by_InputBuffer` hands it every breakend's mate position, beside every record's start and end.
- lines 306-309: `_check_overlap` keeps a transcript only when it overlaps the bounds with the flank added again.

**What correcting it would change.** The buffer's records would keep their transcript rows and the breakend its own-end rows, as vep-rs writes them; class 14 would vanish.

### 14. A coordinate at or past 2,147,478,648 (a breakend mate, a span end or the record's own position) turns the reference's buffer intergenic

#### 14.1 What a reader needs to know first

Ensembl VEP annotates its input in buffers of 5,000 records. For each buffer it works out the span of coordinates the records cover, adds the `--distance` flank (5,000 bases by default) at both ends, and loads the transcripts of the cache regions inside that span; every record's start and end and every breakend's mate position are part of the span. When any coordinate of the buffer plus the flank reaches 2^31 (2,147,483,648), a breakend's mate position, a span's end (`END=2147483647`, or `SVLEN=2147483647` from POS 30425650, each annotated under `--max_sv_size -1`) or a record's own position (`21 2147483648 x N N[21:30425700[ . . SVTYPE=BND`), every record of the buffer is written as `intergenic_variant`: the record that supplied the coordinate, under all its alleles, and every other record loaded with it, whatever transcripts they overlap. The boundary follows the flank exactly: with the default flank the first position that does this is 2,147,478,648 and 2,147,478,647 annotates normally; under `--distance 1` the first such position is 2,147,483,647 and 2,147,483,646 annotates normally. No human chromosome reaches these positions, so the shape appears only in a malformed or synthetic record, and when it does the records it damages are its neighbours in the buffer, which the output does not identify.

#### 14.2 The record

```
21	30425650	control_dup	N	<DUP>	.	.	SVTYPE=DUP;END=30425899
21	30425650	end2_2147478648	N	<BND>	.	.	SVTYPE=BND;CHR2=21;END2=2147478648
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

#### 14.3 What each engine writes

Ensembl VEP 116.2, three rows, every one intergenic:

```
control_dup	21:30425651-30425899	duplication	-	-	-	intergenic_variant	ALLELE_NUM=1;IMPACT=MODIFIER
end2_2147478648	21:30425651	N.	-	-	-	intergenic_variant	ALLELE_NUM=0;IMPACT=MODIFIER
end2_2147478648	21:30425651	N[21:2147478648[	-	-	-	intergenic_variant	ALLELE_NUM=1;IMPACT=MODIFIER
```

The same file with `END2=2147478646`, four rows:

```
control_dup	21:30425651-30425899	duplication	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	ALLELE_NUM=1;IMPACT=MODIFIER;STRAND=1
control_dup	21:30425651-30425899	duplication	ENSG00000240432	ENST00000390690	Transcript	feature_elongation,coding_sequence_variant	ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1
end2_2147478646	21:30425651	N.	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	ALLELE_NUM=0;IMPACT=MODIFIER;STRAND=1
end2_2147478646	21:30425651	N.	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant	ALLELE_NUM=0;IMPACT=HIGH;STRAND=-1
```

vep-rs, for the first file, four rows: the same four as Ensembl VEP's second file, under `end2_2147478648`:

```
control_dup	21:30425651-30425899	duplication	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
control_dup	21:30425651-30425899	duplication	ENSG00000240432	ENST00000390690	Transcript	feature_elongation,coding_sequence_variant
end2_2147478648	21:30425651	N.	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
end2_2147478648	21:30425651	N.	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant
```

Both engines write no row under the mate allele for a mate this far from any transcript, at 2,147,478,646 as at 46,709,983 or 1,000,000,000; the difference is the transcript rows of the whole buffer.

A span's end or a record's own position past the bound turns the buffer the same way. Under `--max_sv_size -1`, `21 30425650 x N <DEL> . . SVTYPE=DEL;END=2147483647` is one row in Ensembl VEP 116.2 and 5,074 transcript rows in vep-rs, every transcript from 30425651 to the end of the chromosome (`SVLEN=2147483647` in place of the `END` gives `21:30425651-2177909297 deletion - intergenic_variant`):

```
# Ensembl VEP 116.2
del_end_2e31m1	21:30425651-2147483647	deletion	-	-	-	intergenic_variant	ALLELE_NUM=1;IMPACT=MODIFIER
# vep-rs (1 of 5,074 rows)
del_end_2e31m1	21:30425651-2147483647	deletion	ENSG00000141956	ENST00000269844	Transcript	transcript_ablation
```

The breakend `21 2147483648 x N N[21:30425700[ . . SVTYPE=BND`, whose own position is past the bound, is two intergenic rows in Ensembl VEP 116.2, one per allele, where vep-rs writes the mate allele's rows at the mate (its own end, at 2,147,483,649, reaches no transcript on either side):

```
# Ensembl VEP 116.2
pos_2e31_brk	21:2147483649	N.	-	-	-	intergenic_variant	ALLELE_NUM=0;IMPACT=MODIFIER
pos_2e31_brk	21:2147483649	N[21:30425700[	-	-	-	intergenic_variant	ALLELE_NUM=1;IMPACT=MODIFIER
# vep-rs
pos_2e31_brk	21:2147483649	N[21:30425700[	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
pos_2e31_brk	21:2147483649	N[21:30425700[	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant
```

At POS 2147483647 the breakend's own position plus the flank also passes 2^31 and the rows are the same; a mate at 2,147,478,647 annotates normally (14.1).

#### 14.4 Why Ensembl VEP writes that

`AnnotationSource.pm` `get_regions_from_coords` (line 190) adds the `--distance` flank to each coordinate it is given when it chooses the cache regions and the buffer's bounds, and `get_all_regions_by_InputBuffer` hands it every record's start and end and every breakend's mate position (lines 251-254 for the mates); `_check_overlap` (lines 306-309) then keeps a transcript only when it overlaps the bounds with the flank added again. The buffer's largest coordinate plus the flank is where the behaviour turns, whichever record supplies it and whether it is a mate position, a span's end or a start: the measurements of 14.1 place the turn at 2^31 for both flank sizes, the bound of a signed 32-bit integer, and past it every record of the buffer is written as intergenic.

#### 14.5 Why vep-rs's rows are the correct ones

A record's annotation depends on its own coordinates and the transcripts it overlaps, not on the other records read with it. vep-rs holds every coordinate in 64 bits, writes the breakend's own-end rows against the transcripts at its own position, writes nothing under the mate allele when no transcript lies at the mate, annotates a span to 2,147,483,647 against every transcript it covers, and annotates the other records of the batch as if the record were absent, which is what Ensembl VEP itself writes for the same records at any coordinate below the bound.

#### 14.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tcontrol_dup\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=30425899\n21\t30425650\tend2_2147478648\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=2147478648\n' > overflow.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/overflow.vcf -o /work/overflow.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i overflow.vcf -o overflow.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

Changing `END2` to `2147478646` restores the transcript rows in Ensembl VEP's output; adding `--distance 1` moves the turn to `END2=2147483647`. Replacing the breakend line with `21\t30425650\tdel_end_2e31m1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=2147483647` and adding `--max_sv_size -1` to both commands shows the same turn from a span's end.

#### 14.7 How the comparator treats this class

The class is `bnd_mate_32bit_buffer_intergenic` in the 116.2 registry of both comparators, counted and never set aside: the records it damages beyond the record's own are its neighbours in a buffer whose boundaries the outputs do not record, so no tuple can be attributed to it with certainty, and crediting vep-rs for an unattributable row is not done. The SNP/indel comparator counts the reference's rows at the Location of a record whose bracket allele names a mate at or past 2,147,478,648; the structural-variant comparator counts such records from the input (a bracket ALT, or `CHR2` with `END2`, else `END`) and the reference's tuples at their Locations (`bnd_mate_32bit_buffer_records`, `bnd_mate_32bit_buffer_reference_tuples`). A span's end (`SVLEN`, else `END`) or a record's own POS at or past the bound is counted as a mate position is, both comparators counting a reference Location whose start or end is at or past 2,147,478,648; the rows of such a record stay counted, not set aside. None of the measured sets carries the shape.

## Defect H: a gVCF variant site is dropped with its real allele

**The defect in one sentence.** Ensembl VEP's parser sends a line to its structural-variant path whenever an ALT carries an angle-bracket token, so a gVCF variant site such as `A G,<NON_REF>` is read as a structural variant of type `G/<NON_REF>`, which has no Sequence Ontology term, and the line is dropped whole as an unsupported type, the real allele with it.

**Why this is a defect.** `<NON_REF>` is not an allele; it is the gVCF's statement that the site may carry alleles the caller did not see. The site's real allele is a variant like any other, and dropping it loses every variant of a gVCF.

**What vep-rs does instead.** vep-rs reads `<NON_REF>` as the absence of a variant: it removes the token from the ALT list when another allele is present and annotates the real alleles, numbered as they stand in the ALT column. A record whose only ALT is `<NON_REF>`, a reference block, is one `intergenic_variant` row.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):
- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 237-245: a line is routed to `create_StructuralVariationFeatures` when REF and the ALTs together are not a run of `ACGT` and an ALT matches `[<\[\]][^\*]+[>\]\[]`, which `<NON_REF>` does.
- `Parser/VCF.pm`, lines 469-481: the ALTs are joined (`G/<NON_REF>`) and passed to `get_SO_term`; with no term for them, the line is marked `vep_skip`.
- `Parser.pm`, lines 712-717: the last branch of `get_SO_term` strips the angle brackets and finds no `%SO_TERMS` entry for `G/NON_REF`.
- `Parser.pm`, lines 821-827: the marked line is dropped before the input buffer.

**What correcting it would change.** The real allele would be annotated, as vep-rs annotates it; class 15 would vanish.

### 15. A gVCF variant site dropped with its real allele

#### 15.1 What a reader needs to know first

A GATK gVCF writes every variant site with `<NON_REF>` beside its real alleles (`A G,<NON_REF>`), the symbolic token standing for "any other allele" so that the site can be merged later. Ensembl VEP's parser sends a line to its structural-variant path whenever an ALT carries an angle-bracket token, reads the type of the whole ALT list (`G/<NON_REF>`), finds no Sequence Ontology term for it and drops the line as an unsupported structural variant: the `G` allele is never annotated. vep-rs reads `<NON_REF>` as the absence of a variant and annotates the real allele; a record whose only ALT is `<NON_REF>` (a reference block) is one `intergenic_variant` row.

#### 15.2 The record

```
21	30425650	gvcf_a_g	A	G,<NON_REF>	.	.	.
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

#### 15.3 What each engine writes

Ensembl VEP 116.2, no row, and in its warnings file:

```
WARNING: line 1 skipped (21 30425650 gvcf_a_g A G,<NON_REF> . . .): G/NON_REF is not a supported structural variant type
```

vep-rs, two rows, the SNV annotated as `A G` alone would be:

```
gvcf_a_g	21:30425650	G	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
gvcf_a_g	21:30425650	G	ENSG00000240432	ENST00000390690	Transcript	missense_variant
```

The same happens with a reference base of `N` (`N G,<NON_REF>`, where vep-rs writes `coding_sequence_variant`, the peptide being unreadable from an ambiguous reference) and with an indel (`N NA,<NON_REF>`, `frameshift_variant`).

#### 15.4 Why Ensembl VEP writes that

`Parser/VCF.pm` lines 237-245 route a line to `create_StructuralVariationFeatures` when REF and the ALTs together are not a run of `ACGT` and an ALT matches `[<\[\]][^\*]+[>\]\[]`, which `<NON_REF>` does; lines 469-481 join the ALTs (`G/<NON_REF>`) and call `get_SO_term`, whose last branch (`Parser.pm` lines 712-717) strips the angle brackets and finds no `%SO_TERMS` entry for `G/NON_REF`, so the line is marked `vep_skip` and dropped before the input buffer (`Parser.pm` lines 821-827).

#### 15.5 Why vep-rs's rows are the correct ones

`<NON_REF>` is not an allele; it is the gVCF's statement that the site may carry alleles the caller did not see. The site's real allele is a variant like any other, and dropping it loses every variant of a gVCF. vep-rs removes the token from the ALT list when another allele is present and parses the record as its real alleles, numbered as they stand in the ALT column.

#### 15.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tgvcf_a_g\tA\tG,<NON_REF>\t.\t.\t.\n' > gvcf.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/gvcf.vcf -o /work/gvcf.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i gvcf.vcf -o gvcf.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

#### 15.7 How the comparator treats this class

The class is `gvcf_variant_site_dropped` in the 116.2 registry of both comparators, set aside: vep-rs's rows at the site, recognised by the record's shape, `<NON_REF>` beside a sequence allele and no true `SVTYPE` (a list of `<NON_REF>` alone, or a block beside a symbolic allele, is the registry class `reference_skipped_record`; under a true `SVTYPE` both engines read the block as an allele of the structural record and agree). The SNP/indel comparator reads the shape from the input VCFs when `--input-dir` (paired with each output by stem) or `--input-vcf` names them, as the runners do, and for an ID-less record from the record's own name in vep-rs's output, `CHROM_POS_REF/ALT/<NON_REF>`, the name Ensembl VEP and vep-rs give such a record, at a Location the reference wrote nothing for; the structural-variant comparator reads the input record and sets aside vep-rs's rows at the Locations and under the alleles vep-rs's parser writes for the site: the anchor base trimmed from an indel's alleles, a lone allele minimised from both ends (`ATG AG,<NON_REF>` is the deletion of `T` at POS + 1, Allele `-`), an insertion at its pair. The skipped-record class names every other unsupported type and leaves this one to its own count.

## Defect I: a dot-notation breakend without `SVTYPE` is written as a pair

**The defect in one sentence.** A dot-notation single breakend (`N.`, `.N`) whose INFO carries no `SVTYPE` takes Ensembl VEP's sequence rules for its coordinates, so its start is POS + 1 when the ALT's length differs from REF's and POS when it does not, its end is POS plus the length of REF minus one, and the record is written as a span rather than the point: the pair `POS-(POS+1)` with the cDNA, CDS and protein positions of an insertion for a REF of one or two bases, a span beginning at POS + 1 for a REF of three or more bases beside an ALT of another length, and one beginning at POS when the two lengths are equal (`NAC` against `NA.`).

**Why this is a defect.** The dot is the breakend's notation, and a single breakend is a point, the base after POS. The same record with `SVTYPE=BND` is that point in Ensembl VEP itself, so the pair follows from which INFO key is present, not from what the record says.

**What vep-rs does instead.** vep-rs reads every breakpoint record's span by the breakpoint rule, `SVTYPE` or not: the single base at POS + 1 unless `SVLEN` gives a span. When INFO names a mate through `CHR2` and `END2`, the `N.` row and the bracket row are both written at the point.

**Where it is in Ensembl VEP's source** (`ensembl-io` and `ensembl-vep`):
- `ensembl-io`, `BaseVCF4.pm`, `get_start`, lines 255-275: adds one to POS for a symbolic ALT or an ALT whose length differs from REF's, and leaves POS alone when the lengths are equal, as a two-base REF's is against `N.`.
- `ensembl-io`, `BaseVCF4.pm`, `get_end`, lines 297-316: ends a record without `SVLEN` or `END` at the start only when the raw INFO carries `SVTYPE` or the ALT carries `<`, `[` or `]`; a dot breakend carries none, so the last branch gives POS plus the length of REF minus one: POS for a one-base REF, POS + 1 for two bases, POS + 2 for three.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 237-245: the line is still routed to the structural path (`^\.\w+|\w+\.$`), and `get_SO_term` reads `BND` from the dot, so the record is a `chromosome_breakpoint` with a sequence record's coordinates.

**What correcting it would change.** The record would be the single base at POS + 1, as it is with `SVTYPE=BND` and as vep-rs writes it; class 16 would vanish.

### 16. A dot-notation breakend without `SVTYPE` written as a pair

#### 16.1 What a reader needs to know first

A single breakend in VCF 4.4 dot notation (`N.`, `.N`) is a point: the base after POS, where one side of the break is known and the other is not. Ensembl VEP reads it as a point when INFO carries `SVTYPE=BND`. Without `SVTYPE` the same ALT falls through the parser's sequence rules for its coordinates, and the length of REF sets the span it gets: the start moves to POS + 1 when the ALT's length differs from REF's and stays at POS when the two are equal, while the end, with no `SVLEN`, no `END` and no `SVTYPE` in the raw INFO and no bracket in the ALT, is POS plus the length of REF minus one. With a one-base REF (`N` against `N.`) the start is POS + 1 and the end POS; with a two-base REF (`NA` against `N.`) the lengths are equal, the start stays at POS and the end is POS + 1; either way the record is written as the pair `POS-(POS+1)` with the cDNA, CDS and protein positions of an insertion. With a REF of three or more bases (`NAC` against `.N`) the start is POS + 1 and the end POS plus the length of REF minus one, `(POS+1)-(POS+2)` for three bases, a span that begins at the point; when the ALT has the REF's length (`NAC` against `NA.` or `.AC`, `NACG` against `NAC.`) the start stays at POS and the span is `POS-(POS+2)` or `POS-(POS+3)`, the point its second base. When INFO names a mate through `CHR2` and `END2`, still without `SVTYPE`, the mate is appended as a bracket allele and both the `N.` row and the `N[21:30425700[` row carry the pair.

#### 16.2 The record

```
21	30425650	ndot_nosvtype	N	N.	.	.	.
21	30425650	ndot_nosvtype_mate	N	N.	.	.	CHR2=21;END2=30425700
21	30425650	ndot_nosvtype_ref3	NAC	.N	.	.	.
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

#### 16.3 What each engine writes

Ensembl VEP 116.2, for the first record:

```
ndot_nosvtype	21:30425650-30425651	chromosome_breakpoint	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
ndot_nosvtype	21:30425650-30425651	chromosome_breakpoint	ENSG00000240432	ENST00000390690	Transcript	coding_sequence_variant	276-277	262-263	88
```

Ensembl VEP 116.2 for the same record with `SVTYPE=BND`, and vep-rs for both:

```
ndot_nosvtype	21:30425651	chromosome_breakpoint	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
ndot_nosvtype	21:30425651	chromosome_breakpoint	ENSG00000240432	ENST00000390690	Transcript	coding_sequence_variant	276	262	88
```

For the second record, whose INFO names a mate, Ensembl VEP 116.2 writes the pair under both alleles:

```
ndot_nosvtype_mate	21:30425650-30425651	N.	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
ndot_nosvtype_mate	21:30425650-30425651	N[21:30425700[	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
ndot_nosvtype_mate	21:30425650-30425651	N.	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant	276-277	262-263	88
ndot_nosvtype_mate	21:30425650-30425651	N[21:30425700[	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant	276-277	262-263	88
```

and vep-rs the same four alleles at the point, `21:30425651`, with the positions `276`, `262` and `88`.

For the third record, whose REF is three bases, Ensembl VEP 116.2 writes the span that begins at the point:

```
ndot_nosvtype_ref3	21:30425651-30425652	chromosome_breakpoint	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
ndot_nosvtype_ref3	21:30425651-30425652	chromosome_breakpoint	ENSG00000240432	ENST00000390690	Transcript	coding_sequence_variant
```

and vep-rs the point, `21:30425651`, under the same allele on both features.

#### 16.4 Why Ensembl VEP writes that

`BaseVCF4.pm` `get_start` (lines 255-275) adds one to POS for a symbolic ALT or an ALT whose length differs from REF's and leaves POS alone when the lengths are equal, and `get_end` (lines 297-316) ends a record without `SVLEN` or `END` at the start only when the raw INFO carries `SVTYPE` or the ALT carries `<`, `[` or `]`; a dot breakend carries none, so the last branch gives POS plus the length of REF minus one. The length of REF therefore sets the span: against the two-character ALT a one-base REF starts at POS + 1 and ends at POS, a two-base REF starts at POS and ends at POS + 1, and both are written as the pair `POS-(POS+1)`; a REF of three or more bases starts at POS + 1 and ends at POS plus the length of REF minus one, `21:30425651-30425652` for `NAC` against `.N` at 30425650, and against an ALT of its own length starts at POS, `21:30425650-30425652` for `NAC` against `NA.`. `Parser/VCF.pm` lines 237-245 still route the line to the structural path (`^\.\w+|\w+\.$`), and `get_SO_term` reads `BND` from the dot, so the record is a `chromosome_breakpoint` with a sequence record's coordinates, an insertion's for a REF of one or two bases.

#### 16.5 Why vep-rs's rows are the correct ones

The dot is the breakend's notation, and the same record with `SVTYPE=BND` is a point in Ensembl VEP itself. vep-rs reads every breakpoint record's span by the breakpoint rule, `SVTYPE` or not: the single base at POS + 1 unless `SVLEN` gives a span.

#### 16.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tndot_nosvtype\tN\tN.\t.\t.\t.\n' > dot.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/dot.vcf -o /work/dot.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i dot.vcf -o dot.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

#### 16.7 How the comparator treats this class

The class is `dot_breakend_without_svtype_pair` in the 116.2 registry of both comparators, set aside: the reference's row at a two-coordinate span and vep-rs's row at the point `chr:(P+1)`, the span's start or the base after it (the end of `chr:P-(P+1)` for a REF of one or two bases, the start of `21:30425651-30425652` for `NAC` against `.N`, the second base of `21:30425650-30425652` for `NAC` against `NA.`), on the same feature and under the same allele are a pair and both leave the adjusted figure, whatever their consequence sets; the allele is the bare `chromosome_breakpoint` of a record without a mate, or the `N.` and `N[21:30425700[` alleles of one whose INFO names a mate. A reference row on a feature vep-rs did not name stays charged.

## Defect J: the missing ALT value is annotated as an allele

**The defect in one sentence.** A record whose ALT column is the missing value `.` has no alternate allele, but when INFO carries `SVTYPE` Ensembl VEP routes the line to its structural path, reads the `.` as a type and annotates the record from INFO alone: a `chromosome_breakpoint` at POS + 1 for `SVTYPE=BND`, a deletion over `END` for `SVTYPE=DEL`.

**Why this is a defect.** `.` is the VCF's missing value, not an allele, and the rows describe an allele the record never states. Ensembl VEP's own sequence path agrees: without `SVTYPE` it writes nothing for the record.

**What vep-rs does instead.** vep-rs writes no row for the record in either case and reports the line at warning level (`No alternate allele`).

**Where it is in Ensembl VEP's source** (`ensembl-io` and `ensembl-vep`):
- `ensembl-io`, `BaseVCF4.pm`, `get_alternatives`, lines 486-490: returns `['.']` for the missing value.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 237-245: the line is routed to the structural path when `SVTYPE` is set, since `N.` is not a run of `ACGT`; lines 472-475 let `SVTYPE` name the type when the ALT is outside the VCF 4.4 forms.
- `ensembl-vep`, `Parser.pm`, lines 710-711: `get_SO_term` reads a breakend from `^\.`.

**What correcting it would change.** The record would write nothing, as it does without `SVTYPE` and as vep-rs writes; class 17 would vanish.

### 17. The missing ALT value annotated as an allele

#### 17.1 What a reader needs to know first

A VCF record whose ALT column is `.` has no alternate allele. Without `SVTYPE` Ensembl VEP writes nothing for it, as vep-rs does. With `SVTYPE` in INFO the line is routed to the structural path, where the `.` is read as a type: `get_SO_term` reads a breakend from a leading dot, and `SVTYPE` names the type when the ALT is outside the VCF 4.4 forms, so `SVTYPE=BND` yields a `chromosome_breakpoint` at POS + 1 and `SVTYPE=DEL;END` a deletion over END, each with rows.

#### 17.2 The record

```
21	30425650	dot_svtype_bnd	N	.	.	.	SVTYPE=BND
21	30425650	dot_svtype_del_end	N	.	.	.	SVTYPE=DEL;END=30425949
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

#### 17.3 What each engine writes

Ensembl VEP 116.2:

```
dot_svtype_bnd	21:30425651	chromosome_breakpoint	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
dot_svtype_bnd	21:30425651	chromosome_breakpoint	ENSG00000240432	ENST00000390690	Transcript	coding_sequence_variant
dot_svtype_del_end	21:30425651-30425949	deletion	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
dot_svtype_del_end	21:30425651-30425949	deletion	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant,5_prime_UTR_variant
```

vep-rs: no row for either, the line reported at warning level (`No alternate allele`).

#### 17.4 Why Ensembl VEP writes that

`BaseVCF4.pm` `get_alternatives` (lines 486-490) returns `['.']` for the missing value; `Parser/VCF.pm` lines 237-245 route the line to the structural path when `SVTYPE` is set, since `N.` is not a run of `ACGT`; `get_SO_term` (`Parser.pm` lines 710-711) reads a breakend from `^\.`, and lines 472-475 of `Parser/VCF.pm` let `SVTYPE` name the type when the ALT is outside the VCF 4.4 forms.

#### 17.5 Why vep-rs's output is the correct one

`.` is the VCF's missing value, not an allele, and Ensembl VEP's own sequence path writes nothing for it. vep-rs refuses the line in every case and says so.

#### 17.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tdot_svtype_bnd\tN\t.\t.\t.\tSVTYPE=BND\n' > missing.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/missing.vcf -o /work/missing.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i missing.vcf -o missing.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

#### 17.7 How the comparator treats this class

The class is `missing_alt_annotated` in the 116.2 registry. The structural-variant comparator reads the input and sets aside the reference's rows at the Location of a record whose ALT is `.` and whose `SVTYPE` Perl reads as true; the SNP/indel comparator reads no input and cannot tell such a row from a record with an allele, so there the class counts nothing and the rows stay charged.

## Defect K: a non-integer coordinate text is annotated and printed

**The defect in one sentence.** Ensembl VEP takes `INFO/END`, `SVLEN` and POS as the text they are written in, checks only that the start and end look like numbers, and annotates a record whose value is a fraction, an exponent or `NaN` with the text in its Location and a fractional distance.

**Why this is a defect.** `END`, `SVLEN` and POS are integers in the VCF specification, and a coordinate is a base position: a half base or a floating-point notation names none. A Location such as `21:30425651-30425949.5` or `21:30425651-NaN` describes no interval of the genome, and the `NaN` record's JSON object carries an end that is not a JSON value.

**What vep-rs does instead.** vep-rs refuses a non-integer `END` as a coordinate that is not a coordinate and reports the line. It reads `SVLEN` as Perl numifies a length (`300.7` is 300 bases, `NaN` is 0, the pair) and a POS by the integer its text denotes, so those records are annotated at the number where Ensembl VEP prints the text.

**Where it is in Ensembl VEP's source** (`ensembl-io` and `ensembl-vep`):
- `ensembl-io`, `BaseVCF4.pm`, `get_end`, lines 307-308: `INFO/END` is returned as the string it is; lines 302-305: a `SVLEN` end is `$start + abs($svlen) - 1`, so `300.7` gives the end `30425950.7` and `NaN` the end `NaN`.
- `ensembl-io`, `BaseVCF4.pm`, `get_start`, lines 255-275: a sequence record's POS is kept as its text, the increment applying only to a structural record's or an indel's start.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser.pm`, `validate_vf`, line 475: the values pass `looks_like_number`, which admits a fraction, an exponent and `NaN`; every comparison with `NaN` is false, so line 483's start-against-end test skips nothing.

**What correcting it would change.** The record would be refused, as vep-rs refuses a non-integer `END`, or read at the number its text denotes, as vep-rs reads `SVLEN` and POS; class 18 would vanish.

### 18. A non-integer `END`, `SVLEN` or POS text annotated and printed

#### 18.1 What a reader needs to know first

`INFO/END` and `SVLEN` are integers in the VCF specification, and so is POS. Ensembl VEP takes each value as written, checks only that the start and end look like numbers, and annotates the record with the text in its Location and a fractional distance: an `END` of `30425949.5`, a `SVLEN` of `300.7` (end `30425950.7`) or `NaN` (end `NaN`), and on the sequence path a POS of `30425949.0`, which denotes an integer but is printed as its text. vep-rs refuses a non-integer `END` as a coordinate that is not a coordinate and reports the line; it reads `SVLEN` as Perl numifies a length (`300.7` is 300 bases, `NaN` is 0, the pair) and a POS by the integer its text denotes, so the record is annotated at the number where Ensembl VEP prints the text.

#### 18.2 The record

```
21	30425650	del_end_float	N	<DEL>	.	.	SVTYPE=DEL;END=30425949.5
21	30425650	del_end_exp	N	<DEL>	.	.	SVTYPE=DEL;END=3.0425949e7
21	30425650	del_svlen_float	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=300.7
21	30425650	del_svlen_nan	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=NaN
21	30425949.0	pos_float_snv	A	G	.	.	.
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

#### 18.3 What each engine writes

Ensembl VEP 116.2:

```
del_end_float	21:30425651-30425949.5	deletion	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	DISTANCE=4270.5
del_end_float	21:30425651-30425949.5	deletion	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant,5_prime_UTR_variant
del_end_exp	21:30425651-3.0425949e7	deletion	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	DISTANCE=4271
del_end_exp	21:30425651-3.0425949e7	deletion	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant,5_prime_UTR_variant
```

For the `SVLEN` and POS records, Ensembl VEP 116.2 prints the text too (the `NaN` record's JSON object carries `"start":30425651,"end":nan`, which is not a JSON value):

```
del_svlen_float	21:30425651-30425950.7	deletion	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	DISTANCE=4269.30000000075
del_svlen_float	21:30425651-30425950.7	deletion	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant,5_prime_UTR_variant
del_svlen_nan	21:30425651-NaN	deletion	-	-	-	intergenic_variant
pos_float_snv	21:30425949.0	G	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	DISTANCE=4271
pos_float_snv	21:30425949.0	G	ENSG00000240432	ENST00000390690	Transcript	upstream_gene_variant	DISTANCE=23
```

vep-rs: no row for the two `END` records, each line reported at warning level (`Invalid END: '30425949.5'`); for the other three the rows of the same record at the number, `21:30425651-30425950` for `SVLEN=300.7`, the pair `21:30425650-30425651` for `SVLEN=NaN` (two transcript rows, `coding_sequence_variant` on ENST00000390690, where the reference wrote one intergenic row), and `21:30425949` for the POS.

#### 18.4 Why Ensembl VEP writes that

`BaseVCF4.pm` `get_end` (lines 307-308) returns `INFO/END` as the string it is, and `Parser.pm` `validate_vf` (line 475) accepts it through `looks_like_number`, which admits a fraction, an exponent and `NaN`; the end then takes part in every arithmetic as a number and is printed as the text. A `SVLEN` reaches the same place through lines 302-305, `$start + abs($svlen) - 1`, so `300.7` gives the end `30425950.7` and `NaN` the end `NaN` (every comparison with it is false, so line 483's start-against-end test skips nothing). On the sequence path POS is kept as its text by `get_start` (lines 255-275, which increment only a structural record's or an indel's start), passes `looks_like_number` and is printed; on the structural path `$start++` numifies it (Defect N). An `END` written with leading zeros or a `+` sign (`030425949`, `+30425949`) is the same shape without the arithmetic: the record is annotated correctly and the Location carries the text, which both comparators read as the number it names.

#### 18.5 Why vep-rs's output is the correct one

A coordinate is a base position; a half base or a floating-point notation names none. vep-rs reads `END` as an integer coordinate and refuses the record otherwise, as it refuses a negative one. `SVLEN` is a length, not a coordinate, and vep-rs reads it as Perl numifies a length, so `300.7` is 300 bases and `NaN` is 0, the pair, where Ensembl VEP prints an end that is not a position; a POS whose text denotes an integer (`30425949.0`, `3e7`) is read as that integer, which is what Ensembl VEP's own structural path makes of it.

#### 18.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tdel_end_float\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425949.5\n' > float.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/float.vcf -o /work/float.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i float.vcf -o float.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

#### 18.7 How the comparator treats this class

The class is `non_integer_end_annotated` in the 116.2 registry of both comparators, set aside: the reference's rows whose Location carries a coordinate that is not an integer. The structural-variant comparator, which reads the input, also sets aside vep-rs's rows of the same record where vep-rs read the value as a number (the rows at the numified end for a fractional `SVLEN`, `21:30425651-30425950`, at the pair for `SVLEN=NaN`, `21:30425650-30425651`, at the integer for a POS text, `21:30425949`), so that neither side is charged and the skipped-record cross-check does not stop on them; the SNP/indel comparator reads no input, locates no vep-rs row of the record, and sets aside the reference's text rows alone. A Location whose coordinates are integers written with a sign or leading zeros is not the class: both comparators read such a coordinate as its number before comparing, so those rows agree with vep-rs's.

## Defect L: an empty or unchanged `<CNV:TR>` expansion is dropped without a row

**The defect in one sentence.** With a FASTA, Ensembl VEP expands a `<CNV:TR>` record to a literal allele string, and when the expansion yields no alternate allele (`RN=.` on the only allele, a `RUC` that is not a positive number) or an allele identical to the reference run (`RUC` times the unit length equal to `SVLEN`), the record passes through the parser without a feature to annotate: no row in any format, and no warning.

**Why this is a defect.** A record the input states is a variant and the parser accepts should be annotated or refused with a message, never dropped in silence. Ensembl VEP itself annotates the same record as a symbolic `tandem_repeat` over the run when no FASTA is given.

**What vep-rs does instead.** vep-rs reads the record as the symbolic tandem repeat it is, with or without a FASTA, and writes the symbolic rows over the run.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):
- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`, lines 382-431: line 410 skips an `RN` of `.` (`next if $records eq '.'`); lines 414 and 419 numify `RUC` as Perl does, so `abc` and `-2` give zero or fewer repeats and `$seq x 0` the empty string; the function returns `REF/` with no alternate, or `AGCAGCAGC/AGCAGCAGC` when the count matches the run.
- `Parser/VCF.pm`, `create_StructuralVariationFeatures`, lines 532-543: builds a VariationFeature from that string.
- `Parser.pm`, lines 553-592: pass the feature, and no allele reaches a consequence.

**What correcting it would change.** The record would be annotated symbolically over its run, as it is without a FASTA and as vep-rs annotates it; class 19 would vanish.

### 19. An empty or unchanged `<CNV:TR>` expansion dropped without a row

#### 19.1 What a reader needs to know first

With a FASTA, Ensembl VEP expands a `<CNV:TR>` record to a literal allele string (12.4). When the expansion yields no alternate allele (`RN=.` on the only allele, a `RUC` that is not a positive number) or an allele identical to the reference run (`RUC` times the unit length equal to `SVLEN`), the allele string has no alternate and the record passes through the parser without a feature to annotate: no row in any format, and no warning. Without a FASTA the same record is a symbolic `tandem_repeat` with rows, and vep-rs writes those rows with or without a FASTA.

#### 19.2 The record

```
21	30425587	same_ruc	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=3
21	30425587	gain_rn_dot	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=.;RUS=AGC;RUC=4
21	30425587	gain_ruc_abc	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=abc
21	30425587	gain_ruc_neg	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=-2
21	30425587	gain1_ruc	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=4
```

The run at 30425588-30425596 is three copies of `AGC` in the GRCh38 sequence. Executed with `--fasta` (Ensembl's chromosome 21 of release 116) on Ensembl VEP 116.2.

#### 19.3 What each engine writes

Ensembl VEP 116.2 with `--fasta`, two rows for the fifth record and none for the first four, with no warning:

```
gain1_ruc	21:30425596-30425597	AGC	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant
gain1_ruc	21:30425596-30425597	AGC	ENSG00000240432	ENST00000390690	Transcript	inframe_insertion
```

vep-rs, for each of the five records (the gain's literal reading pairs with the fifth under class 12):

```
21:30425588-30425596	tandem_repeat	ENSG00000240432	ENST00000390690	Transcript	feature_elongation,coding_sequence_variant
```

#### 19.4 Why Ensembl VEP writes that

`_expand_tandem_repeat_allele_string` (`Parser/VCF.pm` lines 382-431) skips an `RN` of `.` (`next if $records eq '.'`, line 410), numifies `RUC` as Perl does (`abc` and `-2` give zero or fewer repeats, `$seq x 0` the empty string, lines 414 and 419) and returns `REF/` with no alternate, or `AGCAGCAGC/AGCAGCAGC` when the count matches the run; `create_StructuralVariationFeatures` (lines 532-543) builds a VariationFeature from that string, which `Parser.pm` lines 553-592 pass, and no allele reaches a consequence.

#### 19.5 Why vep-rs's rows are the correct ones

A record the input states is a variant and the parser accepts should be annotated or refused with a message, never dropped in silence; and the same record is annotated by Ensembl VEP itself when no FASTA is given. vep-rs reads the record as the symbolic tandem repeat it is in both cases.

#### 19.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425587\tsame_ruc\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=3\n' > empty.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache -v /path/to/fasta:/fasta ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/empty.vcf -o /work/empty.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats --fasta /fasta/Homo_sapiens.GRCh38.dna.chromosome.21.fa
vep -i empty.vcf -o empty.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

#### 19.7 How the comparator treats this class

The class is `cnv_tr_empty_expansion_dropped` in the 116.2 registry. The structural-variant comparator reads the input record (a single `<CNV:TR>` whose `RN` is `.`, whose repeat counts, `RUC` or `RB` in whole units, truncated as Perl's `x` operator truncates them, repeat nothing, or whose counts times the unit lengths equal `SVLEN`, `RUC=3.5` on a three-unit run included, as is an unchanged expansion over several repeat sequences, `RN=3;RUS=AGC,AGC,AGC;RUC=1,1,1` or `RN=2;RUS=AGC,AG;RUC=3,0` on the nine-base run, and a count in fullwidth digits, `RUC=４`, which Perl numifies to 0; a unit of no bases repeats nothing under `RUC`, `"" x n` being the empty expansion, and is set aside like `RN=.`, while under `RB` it divides by zero in Perl, which ends the reference run, so no reference holds such a record) and sets aside vep-rs's `tandem_repeat` rows at the run when the reference wrote nothing on the chromosome from the base before the run to the base after it; a reference that wrote the symbolic rows (no FASTA) or the literal rows (a real gain) leaves nothing to set aside. The SNP/indel comparator reads no input and does not judge the class.

## Defect M: `*` as the only ALT is annotated from `SVTYPE`

**The defect in one sentence.** The ALT `*` marks a base covered by an upstream deletion and names no allele of the record it stands in, but when INFO carries a true `SVTYPE` Ensembl VEP routes the line to its structural path, where `*` falls outside the VCF 4.4 forms, so `SVTYPE` names the type and the record is annotated from INFO alone: a deletion over `END`, or a breakpoint at POS + 1 written under two alleles.

**Why this is a defect.** The specification defines `*` as the mark of a base deleted by an overlapping record, so the record's alleles are the reference and nothing else, and a deletion of 299 bases, or a breakpoint, read from INFO alone describes an allele the record never states. Ensembl VEP's own sequence path agrees, writing no consequence for the record without `SVTYPE`.

**What vep-rs does instead.** vep-rs annotates nothing from INFO when the ALT column names no allele: no row in the default and tab formats, the VCF line carried without a CSQ value, and a JSON object with `"most_severe_consequence":"?"`, which is what Ensembl VEP writes for the same record without `SVTYPE`.

**Where it is in Ensembl VEP's source** (`ensembl-vep` and `ensembl-variation`):
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 237-245: the line is routed to the structural path, REF and ALT together (`N*`) not being a run of `ACGT` and `SVTYPE` being true.
- `ensembl-vep`, `Parser/VCF.pm`, line 469: the ALTs are joined into the allele string `*`; lines 472-475 test it against the VCF 4.4 forms, `^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)`, which `*` fails, so `$type = $info->{SVTYPE}` and `get_SO_term` names a deletion or a breakpoint.
- `ensembl-vep`, `Parser/VCF.pm`, line 499: for the breakpoint the reference base is prepended, giving `N/*`.
- `ensembl-variation`, `StructuralVariationOverlap.pm`, lines 64-77: with no bracket parsed from `N/*`, one allele is built per slash-separated part, `N` and `*`, numbered 1 and 2 (the mechanism of Defect F), so every row of the breakpoint is written twice.

**What correcting it would change.** The record would write nothing, as it does without `SVTYPE` and as vep-rs writes; class 20 would vanish.

### 20. `*` alone under a true `SVTYPE` annotated as the type's class

#### 20.1 What a reader needs to know first

The ALT `*` is the VCF 4.4 marker for a base covered by an upstream deletion (section 1.6.1.5 of the specification): it says that the base at POS is missing from the haplotype that carries the overlapping deletion, and it names no allele of the record it stands in. A record whose only ALT is `*` therefore has no alternate allele. Without `SVTYPE` Ensembl VEP writes nothing for it, the line carried through the VCF output without consequences and the JSON object written with `"most_severe_consequence":"?"` and no consequence array, and vep-rs writes the same. With a true `SVTYPE` in INFO Ensembl VEP routes the line to the structural path, where `*` falls outside the VCF 4.4 forms, so `SVTYPE` names the type and the record is annotated from INFO alone: `SVTYPE=DEL;END` becomes a deletion over `END`, `SVTYPE=BND` a breakpoint at POS + 1 written under two alleles.

#### 20.2 The record

```
21	30425650	star_svtype_del	N	*	.	.	SVTYPE=DEL;END=30425949
21	30425650	star_svtype_del_refa	A	*	.	.	SVTYPE=DEL;END=30425949
21	30425650	star_svtype_bnd	N	*	.	.	SVTYPE=BND
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

#### 20.3 What each engine writes

Ensembl VEP 116.2, for the first record (the second, with REF `A`, writes the same two rows under its own name):

```
star_svtype_del	21:30425651-30425949	deletion	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4271;STRAND=1
star_svtype_del	21:30425651-30425949	deletion	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,coding_sequence_variant,5_prime_UTR_variant	ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1;OverlapBP=276;OverlapPC=50.46
```

and for the third, four rows, every row twice under `ALLELE_NUM` 1 and 2:

```
star_svtype_bnd	21:30425651	chromosome_breakpoint	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4569;STRAND=1
star_svtype_bnd	21:30425651	chromosome_breakpoint	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	ALLELE_NUM=2;IMPACT=MODIFIER;DISTANCE=4569;STRAND=1
star_svtype_bnd	21:30425651	chromosome_breakpoint	ENSG00000240432	ENST00000390690	Transcript	coding_sequence_variant	ALLELE_NUM=1;IMPACT=MODIFIER;STRAND=-1
star_svtype_bnd	21:30425651	chromosome_breakpoint	ENSG00000240432	ENST00000390690	Transcript	coding_sequence_variant	ALLELE_NUM=2;IMPACT=MODIFIER;STRAND=-1
```

with `Use of uninitialized value $region1` warnings in its warnings file (`Utils/VariationEffect.pm` line 120). The VCF output carries each line with a CSQ value, and the JSON object carries the consequences with `"allele_string":"*"` for the deletions and `"N/*"` for the breakpoint.

vep-rs: no row for any of the three in the default and tab formats, no warning; the VCF output carries each line without a CSQ value, and the JSON output writes an object with the record's `id`, `seq_region_name`, `start` (POS), `end` (POS plus the length of REF minus one), `strand`, `allele_string` (REF and every ALT, `N/*`), `assembly_name`, `input` and `"most_severe_consequence":"?"`, which is what Ensembl VEP writes for the same record without `SVTYPE`.

#### 20.4 Why Ensembl VEP writes that

`Parser/VCF.pm` lines 237-245 route the line to the structural path: REF and ALT together (`N*`) are not a run of `ACGT`, and `SVTYPE` is true. Line 469 joins the ALTs into the allele string `*`; lines 472-475 test it against the VCF 4.4 forms, `^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)`, which `*` fails, so `$type = $info->{SVTYPE}` and `get_SO_term` names a deletion or a breakpoint. For the deletion the allele string stays `*` and one allele is built from it. For the breakpoint line 499 prepends the reference base (`N/*`), `_parse_breakends` finds no bracket in it, and `StructuralVariationOverlap.pm` lines 64-77 build one allele per slash-separated part, `N` and `*`, numbered 1 and 2 (the mechanism of Defect F), so every row is written twice.

#### 20.5 Why vep-rs's output is the correct one

`*` names no allele of this record: the specification defines it as the mark of a base deleted by an overlapping record, so the record's alleles are the reference and nothing else, and a deletion of 299 bases, or a breakpoint, read from INFO alone describes an allele the record never states. Ensembl VEP's own sequence path agrees, writing no consequence for the record without `SVTYPE`. vep-rs annotates nothing from INFO when the ALT column names no allele, and carries the line as Ensembl VEP carries the `SVTYPE`-less record.

#### 20.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tstar_svtype_del\tN\t*\t.\t.\tSVTYPE=DEL;END=30425949\n' > star.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/star.vcf -o /work/star.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i star.vcf -o star.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

Replacing `SVTYPE=DEL;END=30425949` with `.` gives no row from either engine.

#### 20.7 How the comparator treats this class

The class is `star_alone_svtype_annotated` in the 116.2 registry. The structural-variant comparator reads the input and sets aside the reference's rows at the Location of a record whose ALTs are all `*` and whose `SVTYPE` Perl reads as true; the SNP/indel comparator reads no input and cannot tell such a row from a record with an allele, so there the class counts nothing and the rows stay charged.

## Defect N: a POS that is not a number is numified on the structural path

**The defect in one sentence.** On its structural path Ensembl VEP derives a record's start by adding one to the raw POS text (`$start++`), and Perl's increment numifies whatever the text is (`.` to 1, `30_425_949` to 31, `-5` to -4), so the record is annotated at a position its line never named.

**Why this is a defect.** POS is a 1-based integer in the VCF specification; a position that is not a number names no base, and a negative one names none either. A `<DEL>` with `END=30425949` and POS `.` becomes a deletion of the chromosome from base 1, 3,254 transcript rows describing an event the line never asserted, while Ensembl VEP's own sequence path refuses the same POS (`Invalid start '.' or end '0' coordinate`).

**What vep-rs does instead.** vep-rs reads POS as an integer coordinate and refuses the line with a warning (`Invalid POS`), in every form.

**Where it is in Ensembl VEP's source** (`ensembl-io` and `ensembl-vep`):
- `ensembl-io`, `BaseVCF4.pm`, `get_start`, lines 255-264: for a record whose raw INFO carries `SVTYPE` or whose ALT carries `<`, `[` or `]`, line 263 runs `$start++` on the raw POS text; Perl's increment reads the text's leading digits and nothing else (`.`, an empty string, `0x1D04A1D` and fullwidth digits are 0, `30_425_949` is 30) and adds one. A sequence record takes the other branch, which increments only when an ALT differs in length from REF, so an SNV keeps the text.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser.pm`, `validate_vf`, lines 475-488: finds a number at `start` and at `end` with `start` at most `end + 1` and keeps the record; the sequence record's text fails `looks_like_number` at line 475.

**What correcting it would change.** The record would be refused, as Ensembl VEP's own sequence path refuses it and as vep-rs refuses every form; class 21 would vanish.

### 21. A non-numeric POS numified by the structural path's increment

#### 21.1 What a reader needs to know first

POS is a 1-based integer in the VCF specification. On the structural path Ensembl VEP derives a record's start by adding one to the raw POS text (`$start++`), and Perl's increment numifies whatever the text is: `.`, an empty column, `0x1D04A1D` and fullwidth digits numify to 0 and increment to 1; `30_425_949` reads its leading digits, 30, and increments to 31; `-5`, a number but no coordinate, increments to -4. The numified start passes `looks_like_number`, the start-against-end test passes, and the record is annotated at a position its line never named: a `<DEL>` with `END=30425949` and POS `.` is a deletion of the chromosome from base 1. A sequence record keeps the text and is refused (`Invalid start '.' or end '0' coordinate`), as vep-rs refuses every form.

#### 21.2 The record

```
21	.	pos_dot_del	N	<DEL>	.	.	SVTYPE=DEL;END=30425949
21	30_425_949	pos_underscore_del	N	<DEL>	.	.	SVTYPE=DEL;END=30425949
21	-5	pos_neg_del	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=-300
21	.	pos_dot_snv	A	G	.	.	.
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`; `--max_sv_size -1` for the first two records).

Once numified, the first two span more than 10,000,000 bases, so at the default `--max_sv_size` Ensembl VEP drops them as oversize (`variant size (30425948) is bigger than --max_sv_size (10000000)`) and annotates them under `--max_sv_size -1`.

#### 21.3 What each engine writes

Ensembl VEP 116.2 under `--max_sv_size -1`, 3,254 rows for the first record, one per transcript the span from base 1 to 30425949 reaches, 3,253 of them `transcript_ablation`; two of them:

```
pos_dot_del	21:1-30425949	deletion	ENSG00000228592	ENST00000262354	Transcript	transcript_ablation	ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1;OverlapBP=23932;OverlapPC=100.00
pos_dot_del	21:1-30425949	deletion	ENSG00000229336	ENST00000282964	Transcript	transcript_ablation	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;OverlapBP=442;OverlapPC=100.00
```

The same 3,254 rows for the second record at `21:31-30425949`, and for the third, at the default cap, one row:

```
pos_neg_del	21:-4-295	deletion	-	-	-	intergenic_variant	ALLELE_NUM=1;IMPACT=MODIFIER
```

An empty POS, `0x1D04A1D` and fullwidth digits write the first record's rows. The fourth record is skipped (the warnings file reads `Invalid start '.' or end '0' coordinate`), as is a sequence record whose POS is empty, `0x1D04A1D` or `30_425_949`; a negative POS on a sequence record passes `looks_like_number` and is written as `21:-5 G - intergenic_variant`.

vep-rs: no row for any of them, each line reported at warning level (`Invalid POS: '.'`, `Invalid POS: '30_425_949'`, `Invalid POS: '-5'`).

#### 21.4 Why Ensembl VEP writes that

`BaseVCF4.pm` `get_start` (lines 255-264) takes the raw POS text and, for a record whose raw INFO carries `SVTYPE` or whose ALT carries `<`, `[` or `]`, runs `$start++` on it (line 263); Perl's increment numifies the text, reading its leading digits and nothing else (`.`, an empty string, `0x1D04A1D` and fullwidth digits are 0, `30_425_949` is 30), and adds one. `Parser.pm` `validate_vf` (lines 475-488) then finds a number at `start` and at `end` and `start` at most `end + 1`, and keeps the record. A sequence record takes the other branch of `get_start`, which increments only when an ALT differs in length from REF, so an SNV keeps the text, and `looks_like_number('.')` fails at line 475.

#### 21.5 Why vep-rs's output is the correct one

A position that is not a number names no base, and a negative one names none either. The 3,254 rows of the first record describe the ablation of every transcript on the first 30 Mb of the chromosome, an event the line never asserted: its POS is the missing value. vep-rs reads POS as an integer coordinate and refuses the line with a warning, as Ensembl VEP's own sequence path refuses the same POS.

#### 21.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t-5\tpos_neg_del\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-300\n' > pos.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/pos.vcf -o /work/pos.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i pos.vcf -o pos.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

Replacing `-5` with `.` and `SVLEN=-300` with `END=30425949`, and adding `--max_sv_size -1` to both commands, gives the 3,254 rows from Ensembl VEP and the refusal from vep-rs.

#### 21.7 How the comparator treats this class

The class is `non_numeric_pos_numified` in the 116.2 registry of both comparators, set aside: the reference's rows of a record whose POS is not a non-negative number (not digits with an optional fraction, exponent, leading `+` or surrounding blanks: `.`, an empty column, `0x1D04A1D`, fullwidth digits, `30_425_949`, `-5`), at Locations vep-rs wrote nothing for; a fractional POS is class 18, and a POS above 2^62 - 1, or at it for a structural record whose first base is POS + 1, is a number and not this class: vep-rs refuses the line as past its coordinate bound and the reference's rows stay charged. The structural-variant comparator reads the input and matches the reference's rows by the start the reference prints for the record, POS + 1 on the structural path and POS on the sequence path (`21:-5 G - intergenic_variant`, the one form a sequence record reaches). The SNP/indel comparator reads no input and sets aside a reference-only Location whose first coordinate is negative, the one trace no written POS leaves; a POS numified to 0 prints a start of 1, which a record at POS 0 prints in both engines, and `30_425_949` prints 31, so those are recognised by the structural-variant comparator alone.

## Defect O: a `<CNV:TR>` allele with `RN=0` drops a real allele of the record

**The defect in one sentence.** With a FASTA, Ensembl VEP expands every allele of a multi-allelic `<CNV:TR>` record to literal sequence and joins them into one allele string, reference first, and an allele whose `RN` is 0 becomes an empty part of that string: an empty middle part hides the multi-allelic shape, so the record is minimised as a single-allele indel and only the empty allele survives, and an empty last part is dropped by `split`, so only the first allele remains.

**Why this is a defect.** The record states two alleles, and both are annotated or neither. Which allele is lost depends on the order the ALTs are listed in, not on what they say, and the lost allele vanishes without a row or a warning.

**What vep-rs does instead.** vep-rs reads `RN` per allele, gives the `RN=0` allele the length of no bases, a loss of the whole run, and annotates every allele of the record, the other as its own gain or loss.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):
- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`, line 408: walks `RN` per allele; line 411: for `RN=0` the loop `for ( 1 .. $records )` runs no iteration, so `$repeat` stays undefined and is pushed as the allele at line 424; line 428: the reference and the alleles are joined with `/`, `AGCAGCAGC//AGCAGCAGCAGC` for `RN=0,1` and `AGCAGCAGC/AGCAGCAGCAGC/` for `RN=1,0`.
- `Parser.pm`, line 978: a multi-allelic string is recognised by `/.+\/.+\/.+/`, which neither string matches, an empty part having no character, so the record is not split per allele.
- `Parser.pm`, line 881: the string is minimised as one indel, where Perl's `split` on `/` drops a trailing empty part and keeps an empty middle one, and the first minimised feature is kept: `AGCAGCAGC/-` for `RN=0,1`, `-/AGC` for `RN=1,0` and `RN=2,0`.

**What correcting it would change.** Both alleles would be annotated, as vep-rs annotates them; class 22 would vanish.

### 22. A real `<CNV:TR>` allele dropped beside an `RN=0` allele

#### 22.1 What a reader needs to know first

A multi-allelic `<CNV:TR>` record gives `RN` per ALT, the number of repeat sequences each alternate allele consists of (12.1); `RN=0` says the allele has no repeat sequence, the run being absent from it. With a FASTA, Ensembl VEP expands every allele to literal sequence (12.4) and joins the alleles into one allele string, reference first; an allele of no repeat sequence is an empty part of that string. When the empty part is in the middle (`RN=0,1`), the string `AGCAGCAGC//AGCAGCAGCAGC` is not recognised as multi-allelic, the record is minimised as a single-allele indel and only the first alternate, the empty one, survives: the real gain of one unit is gone. When the empty part is last (`RN=1,0`, `RN=2,0`), the split drops it and only the first allele remains: the loss of the run is gone. Either way one allele of the record is annotated and the other vanishes without a row or a warning. vep-rs annotates both, the `RN=0` allele as the loss of the whole run and the other as its own gain or loss.

#### 22.2 The record

```
21	30425587	rn_0_1	N	<CNV:TR>,<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=4
21	30425587	rn_1_0	N	<CNV:TR>,<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=1,0;RUS=AGC;RUC=4
21	30425587	rn_2_0	N	<CNV:TR>,<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=2,0;RUS=AGC,AGC;RUC=2,2
21	30425587	rn_0_1_ruc3	N	<CNV:TR>,<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=3
```

The run at 30425588-30425596 is three copies of `AGC` in the GRCh38 sequence (19.2), so the fourth record's real allele, three units, is the run unchanged. Executed with `--fasta` (Ensembl's chromosome 21 of release 116) on Ensembl VEP 116.2.

#### 22.3 What each engine writes

Ensembl VEP 116.2 with `--fasta`, for the first record two rows under the empty allele alone (JSON `"allele_string":"AGCAGCAGC/-"`), the four-unit allele gone, and for the fourth record the same two rows under its own ID, the unchanged allele gone:

```
rn_0_1	21:30425588-30425596	-	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4624;STRAND=1
rn_0_1	21:30425588-30425596	-	ENSG00000240432	ENST00000390690	Transcript	inframe_deletion	ALLELE_NUM=1;IMPACT=MODERATE;STRAND=-1
```

and for the second and third records two rows under the gain alone (JSON `"allele_string":"-/AGC"`), the `RN=0` allele gone:

```
rn_1_0	21:30425596-30425597	AGC	ENSG00000186971	ENST00000334068	Transcript	upstream_gene_variant	ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4623;STRAND=1
rn_1_0	21:30425596-30425597	AGC	ENSG00000240432	ENST00000390690	Transcript	inframe_insertion	ALLELE_NUM=1;IMPACT=MODERATE;STRAND=-1
```

vep-rs, for each record, both alleles at the run, the `RN=0` allele as the loss of the run and the other as its gain, an unchanged allele read as a gain; for the first record, and for the fourth under its own ID, on ENST00000390690:

```
rn_0_1	21:30425588-30425596	tandem_repeat	ENSG00000240432	ENST00000390690	Transcript	feature_truncation,inframe_deletion	ALLELE_NUM=1
rn_0_1	21:30425588-30425596	tandem_repeat	ENSG00000240432	ENST00000390690	Transcript	feature_elongation,coding_sequence_variant	ALLELE_NUM=2
```

The loss row is the symbolic reading of Ensembl VEP's `inframe_deletion` row and pairs with it under class 12 (12.7); the gain row has no partner, for the fourth record as for the first. For the second and third records the gain row pairs with Ensembl VEP's `inframe_insertion` row and the loss row has none. Without a FASTA both engines write the symbolic rows for both alleles and agree.

#### 22.4 Why Ensembl VEP writes that

`Parser/VCF.pm` `_expand_tandem_repeat_allele_string` walks `RN` per allele (line 408); for `RN=0` the loop `for ( 1 .. $records )` at line 411 runs no iteration, so `$repeat` stays undefined and is pushed as the allele (line 424), and line 428 joins the reference and the alleles with `/`: `AGCAGCAGC//AGCAGCAGCAGC` for `RN=0,1`, `AGCAGCAGC/AGCAGCAGCAGC/` for `RN=1,0`. `Parser.pm` line 978 recognises a multi-allelic string by `/.+\/.+\/.+/`, which neither string matches, an empty part having no character, so the record is not split per allele; line 881 minimises it as one indel, where Perl's `split` on `/` drops a trailing empty part and keeps an empty middle one, and the first minimised feature is kept: `AGCAGCAGC/-` for `RN=0,1`, `-/AGC` for `RN=1,0` and `RN=2,0`.

#### 22.5 Why vep-rs's rows are the correct ones

The record states two alleles, and both are annotated or neither. An allele of no repeat sequence is a real allele, the loss of the whole run, and Ensembl VEP itself annotates it as a deletion when it is the one that survives; the four-unit allele is a real allele too, and Ensembl VEP itself annotates it as an insertion when it comes first. Which allele is lost depends on the order the ALTs are listed in, not on what they say. vep-rs reads `RN` per allele, gives the `RN=0` allele the length of no bases, a loss of the run, and annotates every allele of the record.

#### 22.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425587\trn_0_1\tN\t<CNV:TR>,<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=4\n' > rnzero.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache -v /path/to/fasta:/fasta ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/rnzero.vcf -o /work/rnzero.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats --fasta /fasta/Homo_sapiens.GRCh38.dna.chromosome.21.fa
vep -i rnzero.vcf -o rnzero.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats --fasta /path/to/fasta/Homo_sapiens.GRCh38.dna.chromosome.21.fa
```

Swapping the two `RN` values (`RN=1,0`) loses the other allele.

#### 22.7 How the comparator treats this class

The class is `cnv_tr_rn_zero_allele_dropped` in the 116.2 registry. The structural-variant comparator reads the input record, a two-ALT `<CNV:TR>,<CNV:TR>` whose `RN` list has exactly one entry Perl's loop runs no iteration for (0 or negative; `.` is skipped at line 410 and leaves no part, so it is not the shape), when the reference ran with a FASTA and expanded the record (a row at or beside the run under an allele other than `tandem_repeat`). It sets aside vep-rs's `tandem_repeat` rows at the run whose tuple is the dropped allele's alone: the `feature_elongation` rows when the `RN=0` allele comes first and the real allele is a gain or unchanged, which vep-rs reads as a gain, the reference having kept the `RN=0` deletion (the first and fourth records); the `feature_truncation` rows when the `RN=0` allele comes last and the real allele is a gain, the reference having kept the gain (the second and third). A real allele that is a loss (`RUC=2` over the three-unit run) shares the `feature_truncation` tuple with the `RN=0` allele, so the reference's row pairs with that one tuple under class 12 and nothing is left for this class to judge; an unchanged real allele before a trailing `RN=0` (`RN=1,0;RUC=3`) leaves the reference an allele string with no alternate and no row, which is class 19, not this one. The allele the reference kept pairs under class 12 as every expanded `<CNV:TR>` does. Without a FASTA both engines write the symbolic rows for both alleles and nothing is set aside. The SNP/indel comparator reads no input and does not judge the class.

## Records Ensembl VEP 116.2 cannot process

Two shapes end an Ensembl VEP 116.2 run with no output for any record of the file. Nothing is set aside for them: the reference wrote no output for the run, so there is no row to compare against and no class. They are listed because vep-rs reads the same files to completion, with exit code 0 and every other record annotated.

**A coordinate at or above 2^62, 10^30 or `Inf` under `--max_sv_size -1`.** A `<DEL>` whose `END` is 2^62 (`21 30425650 x N <DEL> . . SVTYPE=DEL;END=4611686018427387904`), 2^63, 2^64, 10^30 or `Inf`, or whose `SVLEN` is 2^62 - 1 or more (`SVLEN=4611686018427387903`), a `<BND>` whose `END2` is such a number, a bracket mate at 10^30 (`N[21:1000000000000000000000000000000[`), and any record whose POS is 10^30 (this one in every run, the default cap included): the forked worker dies on the coordinate and the parent aborts the run with exit code 2 and no output file:

```
ERROR: Forked process(es) died: read-through of cross-process communication detected
```

(`Runner.pm` line 583, `_forked_buffer_to_output`). Under the default `--max_sv_size` the spans are oversize and dropped with the usual warning (`variant size (4611686018396962253) is bigger than --max_sv_size (10000000)`), and the run completes. vep-rs refuses the record whose `END` or POS is past 2^62 - 1 at warning level (`Invalid END: '4611686018427387904'`, `Invalid POS`), saturates a `SVLEN` span at 2^62 - 1 (its coordinate bound; under `--max_sv_size -1` that record is annotated against every transcript from POS to the end of the chromosome), and annotates the rest of the file.

**A `<CNV:TR>` repeat unit of no bases under `RB`, with a FASTA.** `21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RUS=,AGC;RB=3,12`: the expansion divides each allele's `RB` by its unit's length (`Parser/VCF.pm` line 414, `shift(@RB) / length($seq)`), and the length of an empty unit is 0, so Perl dies with exit code 255 and no output file:

```
Illegal division by zero at /opt/vep/src/ensembl-vep/modules/Bio/EnsEMBL/VEP/Parser/VCF.pm line 414.
```

Without a FASTA the record is not expanded and the run completes. vep-rs reads the record symbolically with or without a FASTA, the empty unit contributing no bases (`21:30425588-30425596 tandem_repeat` with `feature_truncation,inframe_deletion` on ENST00000390690), and completes the run.

## Scoring against Ensembl VEP 116.2

The comparators score against one Ensembl VEP release at a time, chosen by `--reference-release`, and `116.2` is the default: against it they apply the classes of this page, and every class the 116.2 registry sets aside or counts is listed with its treatment in [docs/concordance-provenance/116.2-divergence-classes.md](concordance-provenance/116.2-divergence-classes.md). Under `--reference-release 115.2`, the release the paper was measured against, they apply the paper's twelve classes instead, classes 1 to 12 as the paper's supplement and this page's previous version document them, and write exactly the reports they wrote for the paper; the five masks of the paper's classes 5 to 10 are listed as not applied in every 116.2 report, each with its reason.

## Reproducing these examples

Every example above is one or two VCF records, from a public dataset or hand-written and printed in full in its section, run through both engines with their default consequence settings.

**The engines.** Every command on this page writes both binaries as `vep`: vep-rs's, a release binary or a source build placed on PATH (`cargo build --release` leaves it at `target/release/vep`), and Ensembl VEP's, the name it has inside the `ensemblorg/ensembl-vep:release_116.2` image (or a local install of release 116.2); the text beside each command names the engine it runs.

**The inputs**, all fetched and laid out by `scripts/data/download_real_world_vcfs.sh`, which pins the ClinVar release:

- ClinVar's GRCh37 and GRCh38 VCFs (`clinvar_grch37.vcf.gz`, `clinvar_grch38.vcf.gz`);
- gnomAD-SV v2.1 (GRCh37, `gnomad_sv_v2.1_chr21.vcf.gz`) and v4.1 (GRCh38, `gnomad_sv_v4.1_chr21.vcf.gz`), chromosome 21;
- the 1000 Genomes structural-variant releases, chromosome 21 (`1kg_sv_chr21.vcf.gz` per assembly: Phase 3 for GRCh37, high-coverage for GRCh38);
- ClinVar's structural variants, chromosome 21 (`clinvar_sv_chr21.vcf.gz`);
- the synthetic sets committed under `tests/sv_validation/<assembly>/`, of which this page uses `02_mnp_complex`, `05_symbolic_del_ins`, `07_cnv_repeat`, `09_breakends`, `10_special_alleles`, `11_multi_allelic`, `12_complex_imprecise` and `13_vcf45_features` (each `.vcf.gz`).

Contig names are Ensembl's (`21`, not `chr21`; `scripts/README.md` describes the canonicalisation).

**The transcripts.** The transcript sets are Ensembl release 116's: Ensembl VEP reads its own cache (`homo_sapiens/116_GRCh37` and `homo_sapiens/116_GRCh38`), and vep-rs reads the JSON conversion of the same cache (`scripts/data/storable_to_json.pl`, described in `docs/cache-setup.md`, laid out as `json_cache/homo_sapiens/116_GRCh37` and `json_cache/homo_sapiens/116_GRCh38`), so both engines see identical transcripts. Transcript coordinates quoted above are read from that cache. The worked examples of classes 13 to 22 run vep-rs on the pruned slice of the GRCh38 cache committed with the golden corpus, `tests/golden/116/GRCh38-release116/json_cache`, which holds the transcripts the corpus's records reach.

**The reference genome.** Both engines take the assembly's Ensembl primary-assembly FASTA on every example, the SNP/indel ones included: `Homo_sapiens.GRCh38.dna.primary_assembly.fa` (Ensembl release 116) and `Homo_sapiens.GRCh37.75.dna.primary_assembly.fa` (Ensembl's GRCh37 primary assembly, first served under release 75 and unchanged since), each with its `.fai` index. The reference bases quoted above come from them (`samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 21:26965974-26966647` gives the 674 bases of class 3's explicit spelling).

**Writing a probe input.** Every probe input on this page is written as its five leading VCF columns, CHROM POS ID REF ALT, tab-separated. To run one, put it in a file after the two header lines `##fileformat=VCFv4.2` and `#CHROM POS ID REF ALT QUAL FILTER INFO` (tab-separated). Use `.` in QUAL and FILTER and, in INFO, `.` or the keys the class shows (`END=`, `SVTYPE=` and the like); `--format vcf` needs all eight columns. The `alone` inputs are such a file with one record, the `plus one SNV` inputs the same file with the SNV line added.

**The commands.** vep-rs:

```
vep -i record.vcf -o vep-rs.txt --offline \
  --json_cache /path/to/json_cache/homo_sapiens/116_GRCh38 \
  --fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --species homo_sapiens --assembly GRCh38 \
  --fork 16 --buffer_size 5000 --force --no_stats
```

Ensembl VEP 116.2, with the flags the reference outputs were generated with:

```
docker run --rm -v "$PWD":/work -v /path/to/vep_cache:/cache -v /path/to/fasta:/fasta ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/record.vcf -o /work/vep.txt --offline --cache \
  --dir_cache /cache --cache_version 116 \
  --species homo_sapiens --assembly GRCh38 --format vcf \
  --fasta /fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --buffer_size 5000 --fork 4 --no_stats --force_overwrite
```

Substitute `GRCh37` and its cache and FASTA for the GRCh37 examples.

**The order-dependent class.** The reference outputs behind the counts used `--fork 4` for the structural-variant sets and `--fork 16` for the SNP/indel datasets. Two classes depend on the other records of a buffer: class 23 on the order in which the buffer's records ask for their regions, and class 14 on the records loaded beside the one whose coordinate reaches 2,147,478,648, which the turn makes intergenic too; their sections say which inputs give which rows. No neighbouring record, `--fork` or `--buffer_size` value changes the outcome of any other class. Flags that do change an outcome are named in the sections that use them: `--max_sv_size` (default 10,000,000 on both engines; `-1` lifts the cap on both) decides class 5's oversize route, and class 12 lowers it on Ensembl VEP alone; `--dont_skip` lets Ensembl VEP annotate a record it would drop (5.4); `--fasta`, which the reference protocol passes, decides whether Ensembl VEP expands a `<CNV:TR>` record (classes 12, 19 and 22); `--distance` moves the coordinate at which class 14's turn happens (14.6). Ensembl VEP's warnings are in `vep.txt_warnings.txt`; vep-rs's are on standard error.

**The comparators** are `scripts/concordance/compare_vep_outputs.py` (SNP/indel) and `scripts/validation/compare_sv_concordance.py` (structural variants); their per-class counters are named in each class's last section, and `scripts/README.md` describes the full measurement.
