# Intended divergences from Ensembl VEP

This page indexes every disagreement between vep-rs and Ensembl VEP that the adjusted concordance figure sets aside. Each defect has a page of its own, which shows why its disagreement is set aside.

## The defects

The defects in page order, each with the classes it produces and its page:

| Defect | In one sentence | Classes | Page |
| --- | --- | --- | --- |
| [A](divergences/defect-a.md) | A later differing region erases the splice-region term | 1, 2 | [Defect A](divergences/defect-a.md) |
| [B](divergences/defect-b.md) | `start_lost` and `start_retained_variant` are both written for one allele | 3, 4 | [Defect B](divergences/defect-b.md) |
| [C](divergences/defect-c.md) | Restated for release 116: a record marked skipped is dropped before annotation | 5 | [Defect C](divergences/defect-c.md) |
| [P](divergences/defect-p.md) | The tests after the chromosome-aware overlap compare chromosome names the overlap did not resolve | 23, 24 | [Defect P](divergences/defect-p.md) |
| [E](divergences/defect-e.md) | A breakend's mate-side row evaluated at the record's own coordinate | 11 | [Defect E](divergences/defect-e.md) |
| [none](divergences/cnv-tr.md) | A representation difference, not a defect: `<CNV:TR>` | 12 | [`<CNV:TR>`](divergences/cnv-tr.md) |
| [F](divergences/defect-f.md) | A breakpoint record without a mate gives its reference allele a row | 13 | [Defect F](divergences/defect-f.md) |
| [G](divergences/defect-g.md) | A coordinate past the 32-bit bound turns the buffer intergenic | 14 | [Defect G](divergences/defect-g.md) |
| [H](divergences/defect-h.md) | A gVCF variant site is dropped with its real allele | 15 | [Defect H](divergences/defect-h.md) |
| [I](divergences/defect-i.md) | A dot-notation breakend without `SVTYPE` is written as a pair | 16 | [Defect I](divergences/defect-i.md) |
| [J](divergences/defect-j.md) | The missing ALT value is annotated as an allele | 17 | [Defect J](divergences/defect-j.md) |
| [K](divergences/defect-k.md) | A non-integer coordinate text is annotated and printed | 18 | [Defect K](divergences/defect-k.md) |
| [L](divergences/defect-l.md) | An empty or unchanged `<CNV:TR>` expansion is dropped without a row | 19 | [Defect L](divergences/defect-l.md) |
| [M](divergences/defect-m.md) | `*` as the only ALT is annotated from `SVTYPE` | 20 | [Defect M](divergences/defect-m.md) |
| [N](divergences/defect-n.md) | A POS that is not a number is numified on the structural path | 21 | [Defect N](divergences/defect-n.md) |
| [O](divergences/defect-o.md) | A `<CNV:TR>` allele with `RN=0` drops a real allele of the record | 22 | [Defect O](divergences/defect-o.md) |

## Concordance figures and classes

vep-rs reports two concordance figures against Ensembl VEP release 116.2. **Raw concordance** keeps every tuple both engines write. A tuple is one output row reduced to its Location, Allele, Feature, Feature_type and Consequence set. Raw F1 counts the tuples the two engines share against the tuples each writes alone, with nothing set aside. **Adjusted concordance** sets aside the tuples of the divergence classes documented on these pages and nothing else.

A class is one shape of disagreement. Each class has been traced to a specific mechanism in one engine's source, reproduced with a command a reader can run, and judged. The record it is reproduced on is public, or hand-written and printed in full on the class's page.

Seventeen of the nineteen classes are the output of fourteen defects in Ensembl VEP's source that vep-rs does not reproduce. Class 5 is the records Ensembl VEP 116 drops before annotation: it writes no row for them in any output format, vep-rs annotates them, and with no Ensembl VEP row to judge, vep-rs's rows are set aside. Class 12, `<CNV:TR>`, is a difference of representation in which neither engine is wrong. A disagreement that is not one of these classes is charged to vep-rs in both figures, whatever its cause.

## How the pages are organised

The pages are organised by defect. There is one page for each of the fifteen defects that produce a class (A to C and E to P; D is retired) and one for the `<CNV:TR>` representation difference. Each defect page opens with the defect in one sentence, why it is a defect, what vep-rs does instead, where it is in Ensembl VEP's source and what correcting it would change. The classes it produces follow.

Each class is a walkthrough of one record, from a public dataset or hand-written and printed in full on the page. It opens with an **In short** paragraph that states the class in a few sentences. It then gives the background a reader needs, the record, both engines' rows and why Ensembl VEP writes what it writes. It closes with why vep-rs's row is the correct one, what a reader can reproduce, and how the comparator scripts recognise the class. The engines, inputs, caches and commands every walkthrough uses are on [Reproducing these examples](divergences/reproducing.md).

## The paper's classes

The accompanying paper measured vep-rs against Ensembl VEP release 115.2, and its supplement documents five classes at the paper's pinned build, four of them defects: defects A to D below and the representation difference. Against release 115.2 those five produced twelve classes, 1 to 12, and release 116 changed the behaviour behind two of the defects.

Defect C, a record marked skipped but annotated against whatever its batch had loaded, is corrected in release 116 by dropping the record before annotation. Its classes 5 to 8 compared the two engines' transcript sets on such a record. Against 116.2 there is no Ensembl VEP row to compare, so class 5 names the dropped records themselves and classes 6 to 8 are retired.

Defect D, transcripts matched to records by coordinate alone, is corrected in release 116.1 by a chromosome-aware overlap. Its classes 9 and 10 have no rows against it and are retired, and the same overlap introduces two shapes of its own, classes 23 and 24. Defects A, B and E and the `<CNV:TR>` representation are unchanged. Class numbers are identifiers and are kept as the paper's and the comparators' reports use them: 6 to 10 are retired, and 23 and 24 are placed with the mechanism that produces them.

## Where each defect is in Ensembl VEP's source

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

Source citations name the Ensembl repository, module and line at release 116: the `release/116.2` tag of `ensembl-vep` and the `release/116` branch of `ensembl-variation` and of `ensembl-io`, as bundled in the `ensemblorg/ensembl-vep:release_116.2` image, on which every worked example runs. Where a defect's page recounts what an earlier release did, it says so and names that release. Every engine row on the defect pages is in Ensembl VEP's default tab output, whose columns are:

```
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
```

Coordinates follow VEP's output conventions: a symbolic deletion's Location starts one base after the VCF POS (the padding base), a breakend's Location is the base after the breakpoint, and an insertion's Location names the two flanking bases.

## The classes

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

Classes 6 to 10 are retired and their numbers are not reused. Against release 115.2 they were:

- 6, the mate end of a breakend above `--max_sv_size`, which Ensembl VEP wrote only when a neighbouring record had loaded the region;
- 7, `<CPX>`, an unsupported type annotated against the batch's transcripts;
- 8, `<NON_REF>`, a reference block annotated as a variant;
- 9, transcripts from another chromosome;
- 10, the `intergenic_variant` row those transcripts displaced.

Release 116 drops the records of 6 to 8 before annotation, so they are class 5. Release 116.1's overlap never hands a record a transcript from another chromosome, so 9 and 10 have no rows. The paper's class 5, an oversize structural variant annotated against whichever transcripts the batch had loaded, keeps its number: the records are the same, and the difference is that Ensembl VEP writes nothing for them.

## How to read the counts

The divergence pages share these scoring terms:

- **Tuple**: one output row reduced to its Location, Allele, Feature, Feature_type and Consequence set.
- **Raw and adjusted figure**: raw F1 counts the tuples the two engines share against the tuples each writes alone, with nothing set aside. Adjusted F1 sets aside the tuples of the classes on these pages and nothing else.
- **Set aside**: a tuple set aside is left out of the adjusted figure. Setting a tuple aside is also called **masking**: a class is fully masked when its rule sets aside every tuple of the class, and a rule is shape-limited when it sets aside only pairs of one exact shape.
- **Charged**: a disagreement that is not one of these classes is charged to vep-rs in both figures, whatever its cause.
- **Pair**: one tuple from each engine that the comparator matches as two readings of one record on one transcript. For classes 1 to 4 and 11 the two tuples share Location, Allele, Feature and Feature_type and differ in the consequence set. For class 12 they share Feature and Feature_type only, because the two engines write the same record at a different Location and Allele ([12.3](divergences/cnv-tr.md#123-what-each-engine-writes)).
- **Pooled**: counts are per input file, and the comparator also reports a pooled count per assembly in which a tuple that recurs in two files is counted once.
- **The reference**: Ensembl VEP's output for the same inputs, which the comparators score vep-rs against. `--reference-release` chooses the release, and `116.2` is the default.
- **The comparators**: the SNP/indel comparator is `scripts/concordance/compare_vep_outputs.py`, for the six SNP/indel datasets; the structural-variant comparator is `scripts/validation/compare_sv_concordance.py`, for the two structural-variant sets.
- **The registry**: the 116.2 registry is the list of classes the comparators apply against Ensembl VEP 116.2, named by identifiers such as `reference_skipped_record`. It sets each class aside or counts it, as [Scoring against Ensembl VEP 116.2](#scoring-against-ensembl-vep-1162) describes.

The **six SNP/indel datasets** are ClinVar GRCh37 and GRCh38 (`clinvar_grch37.vcf.gz`, `clinvar_grch38.vcf.gz`), gnomAD v2.1.1 chr21 (GRCh37, `gnomad_genomes_v2.1.1_chr21.vcf.gz`), gnomAD v4.1 chr21 (GRCh38, `gnomad_genomes_v4.1_chr21.vcf.gz`), 1000 Genomes Phase 3 chr21 (GRCh37, `1kg_integrated_chr21.vcf.gz`) and 1000 Genomes high-coverage chr21 (GRCh38, `1kg_highcov_chr21.vcf.gz`), each scored as one file.

The **two structural-variant sets** are one pool per assembly of sixteen files, scored file by file and pooled: the thirteen synthetic files under `tests/sv_validation/<assembly>/` plus `clinvar_sv_chr21.vcf.gz`, the gnomAD-SV file (`gnomad_sv_v2.1_chr21.vcf.gz` for GRCh37, `gnomad_sv_v4.1_chr21.vcf.gz` for GRCh38) and `1kg_sv_chr21.vcf.gz`.

Counts are tuples or pairs measured with the comparators of this repository against Ensembl VEP 116.2's output for the same inputs.

On GRCh37 35 tuples on each side recur in two files; on GRCh38 16,731 Ensembl VEP and 16,892 vep-rs tuples recur, on records two files carry at identical coordinates. Where the per-file and the pooled count differ, the per-file count is given and the pooled one is named with it (classes [4](divergences/defect-b.md#4-class-3-on-the-structural-variant-sets) and [5](divergences/defect-c.md#5-a-record-ensembl-vep-drops-before-annotation)).

The classes do not overlap on the measured sets. On each structural-variant set the Count column sums to the tuples set aside per file: 23,504 vep-rs and 2,778 Ensembl VEP tuples on GRCh37, 293,809 and 14,327 on GRCh38. Pooled, GRCh38 sets aside 293,647 and 14,326, the 161 class 5 tuples and the one class 4 pair that recur in two files counted once.

On the two ClinVar datasets classes [1](divergences/defect-a.md#1-splice-region-term-lost-to-an-assignment-order) and [3](divergences/defect-b.md#3-start_lost-and-start_retained_variant-emitted-together) together are every tuple the adjusted figure sets aside. On gnomAD v2.1.1 and v4.1 chr21 and 1000 Genomes high-coverage chr21 class 3 is, and on 1000 Genomes Phase 3 chr21 nothing is set aside. On the two structural-variant sets classes [2](divergences/defect-a.md#2-class-1-on-the-structural-variant-sets), 4, 5, [11](divergences/defect-e.md#11-breakend-mate-rows-read-from-the-other-chromosomes-coordinate), [12](divergences/cnv-tr.md#12-cnvtr-read-symbolically-not-expanded) and [23](divergences/defect-p.md#23-a-bracket-breakends-own-end-rows-lost-to-the-first-transcript-that-reaches-it) together are. The adjusted figure is 1.0 on every dataset.

## Records Ensembl VEP 116.2 cannot process

Two shapes end an Ensembl VEP 116.2 run with no output for any record of the file. Nothing is set aside for them: the reference wrote no output for the run, so there is no row to compare against and no class. They are listed because vep-rs reads the same files to completion, with exit code 0 and every other record annotated.

**A coordinate at or above 2^62, 10^30 or `Inf` under `--max_sv_size -1`.** The records that reach it are:

- a `<DEL>` whose `END` is 2^62, 2^63, 2^64, 10^30 or `Inf`, or whose `SVLEN` is 2^62 - 1 or more (`SVLEN=4611686018427387903`); the 2^62 record is:

  ```
  21 30425650 x N <DEL> . . SVTYPE=DEL;END=4611686018427387904
  ```

- a `<BND>` whose `END2` is such a number;
- a bracket mate at 10^30 (`N[21:1000000000000000000000000000000[`);
- any record whose POS is 10^30 (this one in every run, the default cap included).

The forked worker dies on the coordinate and the parent aborts the run with exit code 2 and no output file:

```
ERROR: Forked process(es) died: read-through of cross-process communication detected
```

(`Runner.pm` line 583, `_forked_buffer_to_output`). Under the default `--max_sv_size` the spans are oversize and the run completes, each such record dropped with the usual warning:

```
variant size (4611686018396962253) is bigger than --max_sv_size (10000000)
```

vep-rs refuses the record whose `END` or POS is past 2^62 - 1 at warning level (`Invalid END: '4611686018427387904'`, `Invalid POS`) and annotates the rest of the file. It saturates a `SVLEN` span at 2^62 - 1, its coordinate bound; under `--max_sv_size -1` that record is annotated against every transcript from POS to the end of the chromosome.

**A `<CNV:TR>` repeat unit of no bases under `RB`, with a FASTA.** The record is:

```
21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RUS=,AGC;RB=3,12
```

The expansion divides each allele's `RB` by its unit's length (`Parser/VCF.pm` line 414, `shift(@RB) / length($seq)`), and the length of an empty unit is 0, so Perl dies with exit code 255 and no output file:

```
Illegal division by zero at /opt/vep/src/ensembl-vep/modules/Bio/EnsEMBL/VEP/Parser/VCF.pm line 414.
```

Without a FASTA the record is not expanded and the run completes. vep-rs reads the record symbolically with or without a FASTA, the empty unit contributing no bases (`21:30425588-30425596 tandem_repeat` with `feature_truncation,inframe_deletion` on ENST00000390690), and completes the run.

## Scoring against Ensembl VEP 116.2

The comparators score against one Ensembl VEP release at a time, chosen by `--reference-release`, and `116.2` is the default. Against it they apply the classes of these pages, and every class the 116.2 registry sets aside or counts is listed with its treatment in [docs/concordance-provenance/116.2-divergence-classes.md](concordance-provenance/116.2-divergence-classes.md). Under `--reference-release 115.2`, the release the paper was measured against, they apply the paper's twelve classes instead, classes 1 to 12 as the paper's supplement and this page's previous version document them, and write exactly the reports they wrote for the paper. The five masks of the paper's classes 5 to 10 are listed as not applied in every 116.2 report, each with its reason.
