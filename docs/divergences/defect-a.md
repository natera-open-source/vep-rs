# Defect A: a later differing region erases the splice-region term

**The defect in one sentence.** Ensembl VEP decides the splice-region question once per changed place of a variant and keeps only the last answer, so a change inside the splice-region window loses its `splice_region_variant` term whenever a later change outside the window follows it.

**Why this is a defect.** A consequence term describes where a changed base falls. Whether base 3 of an exon is in the splice-region window does not depend on whether the same variant also changes base 5. Ensembl VEP agrees when the bases are changed one at a time: it writes `splice_region_variant` for the first change alone, and only the presence of the second change removes the term. An annotation that an unrelated second change can switch off is not an annotation of the first change.

**What vep-rs does instead.** vep-rs keeps the term when any differing region, a run of changed bases, earns it. That is how Ensembl VEP already treats its other thirteen splice flags, computed [in the same loop](#14-why-ensembl-vep-writes-that), which are only ever switched on and never overwritten; the splice-region slot is the one exception, and the exception is the defect.

**Where it is in Ensembl VEP's source** (`ensembl-variation`):
- `Bio/EnsEMBL/Variation/BaseTranscriptVariationAllele.pm`, `_intron_effects`, line 215: the verdict is stored by assignment, `$intron_effects->{splice_region} = _intron_overlap(...)`, inside the `foreach my $region` loop that starts at line 117.
- `Utils/VariationEffect.pm`, lines 87-110: `_intron_overlap` returns a definite 0 or 1 for each region, so a later 0 replaces an earlier 1.
- `VariationFeatureOverlapAllele.pm`, line 405: `_get_differing_regions` splits a multi-base change into the regions the loop visits.

**What correcting it would change.** Classes 1 and 2 would vanish: Ensembl VEP would keep `splice_region_variant` wherever any region earned it, which is what vep-rs writes.

**The classes it produces.**
- Class 1: the ClinVar GRCh37 and GRCh38 datasets, 440 and 1,994 pairs; both members of each pair are set aside.
- Class 2: the same shape on the structural-variant sets, scored by the structural-variant comparator, 20 pairs on the synthetic complex-substitution set; both members of each pair are set aside.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Splice-region term lost to an assignment order

**In short.** ClinVar variation 630437 changes base 3 of BRCA2 exon 6, inside the splice-region window, and after it base 5, outside the window. Ensembl VEP keeps only the splice-region verdict of the last changed place, so base 5's verdict erases the term base 3 earned and Ensembl VEP writes `missense_variant` alone. vep-rs keeps the term any changed place earns and writes `missense_variant,splice_region_variant`, which is correct: base 3 is in the window by the term's own definition, and Ensembl VEP itself writes the term when that base is changed alone.

### 1.1 What a reader needs to know first

Three Sequence Ontology terms appear in this class:

- `splice_region_variant` (SO:0001630): "a sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron". Ensembl VEP applies exactly that window on each side of every intron boundary.
- `missense_variant` (SO:0001583): "a sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved".
- `synonymous_variant` (SO:0001819): "a sequence variant where there is no resulting change to the encoded amino acid".

A variant that changes several bases at once (a multi-nucleotide variant or a deletion-insertion) may change bases that are not adjacent. Ensembl VEP compares the reference and alternate alleles position by position and evaluates each run of differing bases, called a differing region, on its own. A term that holds for any of those regions holds for the variant.

### 1.2 The record

ClinVar variation 630437, from the ClinVar GRCh37 VCF (`clinvar_grch37.vcf.gz` in the layout of `scripts/data/download_real_world_vcfs.sh`, which pins the ClinVar release; the variation id and alleles are stable across releases):

```
13	32900381	630437	GTA	CTT	.	.	CLNHGVS=NC_000013.10:g.32900381_32900383delinsCTT;CLNVC=Indel;GENEINFO=BRCA2:675
```

The transcript is ENST00000380152 (BRCA2, forward strand). Its exon 6 of 27 spans 32900379-32900419 and is preceded by intron 5, 32900288-32900378, so exon bases 1 to 3 are 32900379, 32900380 and 32900381.

The reference GTA at 32900381-32900383 is exon bases 3, 4 and 5. The alternate CTT changes base 3 (G to C) and base 5 (A to T) and leaves base 4 (T) as it was, so the variant is two differing regions of one base each, 32900381 and 32900383. Base 3 lies inside the 1-3 exonic window of the splice-region definition; base 5 does not.

### 1.3 What each engine writes

The two rows as a table: ClinVar variation 630437 on ENST00000380152, one row from each engine, the differing cells in bold.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `630437` | `630437` |
| Location | `13:32900381-32900383` | `13:32900381-32900383` |
| Allele | `CTT` | `CTT` |
| Gene | `ENSG00000139618` | `ENSG00000139618` |
| Feature | `ENST00000380152` | `ENST00000380152` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | **`missense_variant`** | **`missense_variant,splice_region_variant`** |
| cDNA_position | `711-713` | `711-713` |
| CDS_position | `478-480` | `478-480` |
| Protein_position | `160` | `160` |
| Amino_acids | `V/L` | `V/L` |
| Codons | `GTA/CTT` | `GTA/CTT` |
| Existing_variation | `-` | `-` |
| Extra | `IMPACT=MODERATE;STRAND=1` | `IMPACT=MODERATE;STRAND=1` |

The rows agree in every column but Consequence, where vep-rs adds `splice_region_variant`. The other two BRCA2 transcripts the record touches, ENST00000530893 and ENST00000544455, carry the same difference.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP evaluates the splice-region window once per differing region and stores each verdict by assignment into a single slot, so the last region's verdict replaces the earlier ones. For this record the first region (exon base 3) sets the slot to 1 and the second region (exon base 5) sets it back to 0. By the time the consequence is read, the term the first region earned is gone.

The other thirteen splice flags computed in the same loop are only ever set to 1, so they cannot be lost this way. `splice_region_variant` is the one term the order of regions can remove.

Source (`ensembl-variation`):

- `Bio/EnsEMBL/Variation/BaseTranscriptVariationAllele.pm`, `_intron_effects`, line 215: `$intron_effects->{splice_region} = _intron_overlap(...)`, inside the `foreach my $region` loop that starts at line 117.
- `Utils/VariationEffect.pm`, lines 87-110: `_intron_overlap` returns a defined 0 or 1.
- `VariationFeatureOverlapAllele.pm`, line 405: `_get_differing_regions` supplies the regions.

### 1.5 Why vep-rs's row is the correct one

The record changes base 32900381, the third base of exon 6. That position is inside the window the term's own definition names ("within 1-3 bases of the exon").

A term defined on the position of a changed base cannot depend on what other bases the same record changes. Ensembl VEP agrees when the base is changed alone ([1.6](#16-try-it-yourself)): the same nucleotide change at the same position is a `splice_region_variant` by Ensembl VEP's own reading. Only the presence of a second changed base two positions downstream removes it.

### 1.6 Try it yourself

Annotate the two differing bases as separate records and then together, with the Ensembl VEP command of [the closing section](reproducing.md) (GRCh37 cache and FASTA):

| Input | Ensembl VEP, ENST00000380152 | vep-rs, ENST00000380152 |
| --- | --- | --- |
| `13 32900381 . G C` | `missense_variant,splice_region_variant` | `missense_variant,splice_region_variant` |
| `13 32900383 . A T` | `synonymous_variant` | `synonymous_variant` |
| `13 32900381 . GTA CTT` | `missense_variant` | `missense_variant,splice_region_variant` |

Each part alone is annotated identically by both engines. The two parts as one record lose, in Ensembl VEP only, the term the first part earned. No flag changes the outcome: this class is deterministic.

### 1.7 How the comparator treats this class

`_is_covered_splice_region_swap` in `scripts/concordance/compare_vep_outputs.py` matches a pair when every term in Ensembl VEP's set is also in vep-rs's and vep-rs's set adds `splice_region_variant`, alone or together with `intron_variant`. Both members of a matching pair are set aside.

The rule accepts `intron_variant` beside `splice_region_variant` for one case, a region Ensembl VEP skips in a short intron. The `intron_variant` allowance exists because Ensembl VEP's per-region loop steps over a region that lies in an intron of twelve bases or fewer before it sets either the intronic flag or the splice-region slot (`BaseTranscriptVariationAllele.pm` lines 139-141 and 172-174). A region skipped that way withholds both terms at once. None of the ClinVar pairs uses the allowance; every one of them adds `splice_region_variant` alone.

The rule is deliberately narrower than the family of splice differences. A pair in which Ensembl VEP lacks a splice sub-term such as `splice_donor_variant`, or in which vep-rs adds `splice_polypyrimidine_tract_variant`, is never set aside, because the assignment at line 215 can remove no term but `splice_region_variant`.

Counts:

- ClinVar GRCh37: 440 pairs; ClinVar GRCh38: 1,994 pairs, every one of that shape.
- gnomAD and 1000 Genomes datasets: none.
- Structural-variant sets: 20 pairs (class 2).

The SNP/indel comparator's summary groups every splice-family pair into descriptive buckets for reporting. The buckets do not decide exclusion; this rule does, and the summary prints beside each bucket how many of its pairs the rule excluded.

## 2. Class 1 on the structural-variant sets

**In short.** The synthetic complex substitution `synth_cpx_sub_0780`, in a GRCh38 file that the structural-variant comparator scores ([2.1](#21-what-a-reader-needs-to-know-first)), replaces `AAAG` with `CA` in SOD1 exon 2. That gives two differing regions: exon base 2, inside the splice-region window, and exon bases 4 and 5, outside it. Ensembl VEP writes `frameshift_variant` alone on ENST00000270142; vep-rs writes `frameshift_variant,splice_region_variant`. This is class 2, the shape of class 1 on a structural-variant set: exon base 2 is inside the window the term names, and Ensembl VEP writes the term when that base is changed alone, so vep-rs's row is the correct one.

### 2.1 What a reader needs to know first

The structural-variant sets are scored by a different script from the SNP/indel datasets, but they contain sequence variants too (the synthetic complex-substitution set is nothing else). Nothing new happens in either engine; the same predicate is applied by the second comparator. The terms are those of [1.1](#11-what-a-reader-needs-to-know-first).

### 2.2 The record

`synth_cpx_sub_0780` from the synthetic complex-substitution set `tests/sv_validation/grch38/02_mnp_complex.vcf.gz` (GRCh38):

```
21	31663791	synth_cpx_sub_0780	AAAG	CA	.	PASS	.
```

The transcript is ENST00000270142 (SOD1, forward strand, canonical, GRCh38 31659693-31668931). Its exon 2 of 5 spans 31663790-31663886 and is preceded by intron 1 (31659842-31663789), so 31663791 is exon base 2 and 31663793-31663794 are exon bases 4 and 5.

Comparing `AAAG` with `CA` position by position gives two differing regions: 31663791 (A to C) and 31663793-31663794 (the two reference bases with no alternate counterpart). Position 31663792 (A) is unchanged. Base 2 of the exon lies inside the 1-3 window; bases 4 and 5 do not.

### 2.3 What each engine writes

The two rows as a table: `synth_cpx_sub_0780` on ENST00000270142, one row from each engine, the differing cells in bold.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `synth_cpx_sub_0780` | `synth_cpx_sub_0780` |
| Location | `21:31663791-31663794` | `21:31663791-31663794` |
| Allele | `CA` | `CA` |
| Gene | `ENSG00000142168` | `ENSG00000142168` |
| Feature | `ENST00000270142` | `ENST00000270142` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | **`frameshift_variant`** | **`frameshift_variant,splice_region_variant`** |
| cDNA_position | `151-154` | `151-154` |
| CDS_position | `74-77` | `74-77` |
| Protein_position | `25-26` | `25-26` |
| Amino_acids | `ES/AX` | `ES/AX` |
| Codons | `gAAAGt/gCAt` | `gAAAGt/gCAt` |
| Existing_variation | `-` | `-` |
| Extra | `IMPACT=HIGH;STRAND=1` | `IMPACT=HIGH;STRAND=1` |

The rows agree in every column but Consequence, where vep-rs adds `splice_region_variant`. The record has 47 rows on each side. The rows of 19 SOD1 transcripts differ the same way: eleven coding transcripts (`frameshift_variant`), five coding transcripts flagged as nonsense-mediated decay (NMD) targets (`frameshift_variant,NMD_transcript_variant`), one whose exon 2 is untranslated (`3_prime_UTR_variant,NMD_transcript_variant`) and two non-coding ones (`non_coding_transcript_exon_variant`). The other 28 rows are identical.

### 2.4 Why Ensembl VEP writes that

As [1.4](#14-why-ensembl-vep-writes-that): the first differing region (exon base 2) sets `splice_region` to 1 and the second (exon bases 4 and 5) sets it back to 0, on each of the nineteen transcripts.

### 2.5 Why vep-rs's row is the correct one

As [1.5](#15-why-vep-rss-row-is-the-correct-one): base 31663791 is the second base of exon 2, inside the window the term names, and Ensembl VEP calls it so when it is changed alone.

### 2.6 Try it yourself

Ensembl VEP 116.2, GRCh38 cache and FASTA, the command of [the closing section](reproducing.md):

- `21 31663791 . A C` alone gives ENST00000270142 `missense_variant,splice_region_variant`;
- `21 31663792 . AAG A` alone (the second differing region, Location `21:31663793-31663794`, Allele `-`) gives `frameshift_variant`;
- the full record `21 31663791 . AAAG CA` gives `frameshift_variant`.

vep-rs gives the same two rows for the two parts and `frameshift_variant,splice_region_variant` for the full record.

### 2.7 How the comparator treats this class

`filter_registry_swap_pairs` in `scripts/validation/compare_sv_concordance.py` builds the pairs of the structural-variant comparison, one tuple on each side per Location, Allele, Feature and Feature_type. It applies to them the same two rules the SNP/indel comparator exports (`EXCLUDING_RULES` of `scripts/concordance/compare_vep_outputs.py`, imported so the two scripts cannot drift), setting both members of a matching pair aside.

Counts: 20 pairs on the GRCh38 set and none on the GRCh37 set. Nineteen are the transcripts of `synth_cpx_sub_0780`. The twentieth is ENST00001077173 on `synth_mnp_0852` (`21 31667270 synth_mnp_0852 CTT TTA`, whose two differing bases, 31667270 and 31667272, stand either side of an unchanged base), with `NMD_transcript_variant,synonymous_variant` against `NMD_transcript_variant,splice_region_variant,synonymous_variant`. Fully masked; the report counts them as `excluded_registry_swap_rust` and `excluded_registry_swap_perl`, split by rule in `excluded_registry_swap_by_bucket` (`covered_splice_region_swap`).
