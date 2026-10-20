# Defect B: `start_lost` and `start_retained_variant` are both written for one allele

**The defect in one sentence.** Ensembl VEP decides "the start codon is lost" and "the start codon is retained" with two separate tests that never consult each other, so one allele can be given both terms, and on the records of this class it is.

**Why this is a defect.** The two terms contradict each other: a start codon cannot be both lost and kept by the same change. Which of the two is wrong depends on the allele. On a sequence deletion, an allele given as bases, Ensembl VEP rebuilds the leading untranslated region (UTR) and the coding sequence with the edit applied. The loss test behind `start_lost` reads that rebuilt sequence at the old UTR length, finds `TGA` where `ATG` stood, and reports the start altered; `start_lost` is the wrong member.

On a structural allele, a symbolic one such as `<DEL>`, the retention test behind `start_retained_variant` returns "not altered" for every structural allele before reading anything, and the term is its negation; `start_retained_variant` is the wrong member. In both cases one term was computed from the record and the other was not.

**What vep-rs does instead.** vep-rs writes one of the two terms, the one the record supports: `start_retained_variant` on a sequence deletion whose rebuilt sequence still begins with `ATG`, `start_lost` on a structural deletion whose interval contains the codon. The two are never written together.

**Where it is in Ensembl VEP's source** (`ensembl-variation`, `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`):
- `start_lost`, lines 851-903. The sequence branch (lines 869-883) calls `_inv_start_altered` (lines 906-944), whose `substr($utr_and_translateable, $atg_start, 3)` at line 939 reads at the old UTR length and whose `ne 'ATG'` test is line 941. The structural branch (lines 886-896, at 891 and 894) is an overlap of the span with the codon's three coordinates.
- `start_retained_variant`, lines 947-963, returning `!_ins_del_start_altered(@_)` at line 961 without consulting the loss test.
- `_ins_del_start_altered`, lines 1028-1075, returning 0 for a `TranscriptStructuralVariationAllele` at line 1037 and comparing the rebuilt string's tail with the original coding sequence at line 1071.

**What correcting it would change.** Correcting the misread at line 939 removes the sequence pairs; correcting the structural shortcut at line 1037 removes the structural pairs. With both corrected, classes 3 and 4 vanish and Ensembl VEP writes one term per allele, as vep-rs does.

**The classes it produces.**
- Class 3: the SNP/indel datasets (ClinVar GRCh37 and GRCh38, gnomAD v2.1.1 and v4.1 chr21, 1000 Genomes high-coverage chr21: 261, 1,466, 3, 127 and 1 pairs); both members of each pair are set aside.
- Class 4: the structural-variant sets, 3 pairs on GRCh37 and 38 on GRCh38 (37 distinct); both members of each pair are set aside.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 3. `start_lost` and `start_retained_variant` emitted together

**In short.** Two deletions touch a start codon, one of each allele kind, bases and symbolic. ClinVar variation 4629713 deletes one base of a two-base run that continues into the start codon of a RUNX1 transcript, so the codon re-forms one base upstream; the symbolic `<DEL>` `gnomAD-SV_v3_DEL_chr21_88237371` removes the whole start codon of an ADAMTS5 transcript. On both, Ensembl VEP writes `start_lost` and `start_retained_variant` together; vep-rs writes only one of them, `start_retained_variant` on the sequence deletion and `start_lost` on the structural one. The two terms are alternatives by definition, and each record's own sequence or coordinates say which one holds, so vep-rs's rows are the correct ones.

### 3.1 What a reader needs to know first

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

### 3.2 The record

Two records show the class, one of each allele kind.

**A sequence deletion.** ClinVar variation 4629713, from the ClinVar GRCh37 VCF:

```
21	36265245	4629713	AT	A	.	.	CLNHGVS=NC_000021.8:g.36265247del;CLNVC=Deletion;GENEINFO=RUNX1:861
```

VEP trims the shared A and annotates the deletion of the T at 36265246 (Location `21:36265246`, Allele `-`). The transcript is ENST00000486278 (RUNX1, reverse strand). Its exon 3 of 7 spans 36265222-36265260. The 5' UTR is exons 1 and 2 plus the first 14 bases of exon 3, 36265260 down to 36265247, and the coding sequence begins at 36265246 (cDNA position 217).

On the forward strand the reference at 36265244-36265248 reads `CATTC`. The transcript is on the reverse strand, so its positions run downward: read on the transcript strand, from 36265248 down to 36265244, it is `GAATG`. The UTR ends `...GA` at 36265248-36265247 and the start codon ATG occupies 36265246, 36265245 and 36265244. The deleted base is the A of that ATG.

Removing it from `...GA|ATGAATCCT...` gives `...G|ATGAATCCT...`: the UTR is one base shorter and the coding sequence, ATG included, is byte for byte what it was. ClinVar's own name for the record, `g.36265247del`, describes the same edit as the deletion of the UTR's last base.

**A structural deletion.** `gnomAD-SV_v3_DEL_chr21_88237371`, from gnomAD-SV v4.1 chr21 (GRCh38, `gnomad_sv_v4.1_chr21.vcf.gz`; the `v3` in the record id is gnomAD's own naming of the call set inside the v4.1 release):

```
21	26965974	gnomAD-SV_v3_DEL_chr21_88237371	N	<DEL>	1	LOWQUAL_WHAM_SR_DEL;OUTLIER_SAMPLE_ENRICHED	END=26966647;SVLEN=673;SVTYPE=DEL
```

The deletion removes 26965975-26966647. The transcript is ENST00000284987 (ADAMTS5, reverse strand, canonical). Its exon 1 of 8 spans 26965288-26967088. The 5' UTR is 26967088 down to 26966392 (697 bases) and the start codon ATG occupies 26966391, 26966390 and 26966389 (forward-strand reference `CAT` at 26966389-26966391).

The deleted interval contains the last 256 UTR bases, the whole start codon and the first 417 coding bases (26966391 down to 26965975). That is what the row's CDS position `?-417` and protein position `?-139` say.

### 3.3 What each engine writes

The sequence deletion, on ENST00000486278, one row from Ensembl VEP and one from vep-rs, the differing cells in bold:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | 4629713 | 4629713 |
| Location | `21:36265246` | `21:36265246` |
| Allele | `-` | `-` |
| Gene | `ENSG00000159216` | `ENSG00000159216` |
| Feature | `ENST00000486278` | `ENST00000486278` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | **`frameshift_variant,start_lost,start_retained_variant`** | **`frameshift_variant,start_retained_variant`** |
| cDNA_position | `217` | `217` |
| CDS_position | `1` | `1` |
| Protein_position | `1` | `1` |
| Amino_acids | `M/X` | `M/X` |
| Codons | `Atg/tg` | `Atg/tg` |
| Existing_variation | `-` | `-` |
| Extra | `IMPACT=HIGH;STRAND=-1` | `IMPACT=HIGH;STRAND=-1` |

The structural deletion, on ENST00000284987, one row from Ensembl VEP and one from vep-rs, the differing cells in bold:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `gnomAD-SV_v3_DEL_chr21_88237371` | `gnomAD-SV_v3_DEL_chr21_88237371` |
| Location | `21:26965975-26966647` | `21:26965975-26966647` |
| Allele | `deletion` | `deletion` |
| Gene | `ENSG00000154736` | `ENSG00000154736` |
| Feature | `ENST00000284987` | `ENST00000284987` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | **`frameshift_variant,start_lost,feature_truncation,start_retained_variant,5_prime_UTR_variant`** | **`frameshift_variant,start_lost,feature_truncation,5_prime_UTR_variant`** |
| cDNA_position | `442-1114` | `442-1114` |
| CDS_position | `?-417` | `?-417` |
| Protein_position | `?-139` | `?-139` |
| Amino_acids | `-` | `-` |
| Codons | `-` | `-` |
| Existing_variation | `-` | `-` |
| Extra | `IMPACT=HIGH;STRAND=-1;OverlapBP=673;OverlapPC=1.37` | `IMPACT=HIGH;STRAND=-1;OverlapBP=673;OverlapPC=1.37` |

In each pair the rows agree in every column but Consequence: Ensembl VEP writes both start terms, vep-rs writes one. The second record's other four coding transcripts, ENST00000970346, ENST00001108182, ENST00001108183 and ENST00001131497, carry the same difference, and its sixth row, `upstream_gene_variant` on ENST00000652031, is identical on both sides.

### 3.4 Why Ensembl VEP writes that

On the sequence deletion, the retention test rebuilds UTR plus coding sequence with the A removed, finds that the rebuilt string still ends with the original coding sequence, and reports the start not altered, so `start_retained_variant` is emitted. The loss test reads three bases at the old UTR length in the rebuilt string, finds `TGA` where `ATG` stood, and reports the start altered, so `start_lost` is emitted too.

On the structural deletion, `start_lost` is emitted because the deleted span overlaps the codon's three bases. `start_retained_variant` is emitted because the retention test returns 0 for every structural allele before reading anything, and the term is its negation.

Source (`ensembl-variation`):

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, four functions:
  - `start_lost`, lines 851-903. The sequence-variant branch is lines 869-883: line 871 calls `_inv_start_altered`. The structural branch is lines 886-896: lines 891 and 894 test an overlap of the span with the codon's three coordinates.
  - `_inv_start_altered`, lines 906-944: its `substr($utr_and_translateable, $atg_start, 3)` at line 939 reads at the old UTR length, and its `ne 'ATG'` test is line 941.
  - `start_retained_variant`, lines 947-963: returns `!_ins_del_start_altered(@_)` at line 961.
  - `_ins_del_start_altered`, lines 1028-1075: returns 0 for a `TranscriptStructuralVariationAllele` at line 1037 and compares the rebuilt string's tail with the original coding sequence at line 1071.

### 3.5 Why vep-rs's row is the correct one

The two terms are mutually exclusive by definition: a start that "remains" is not lost.

On the sequence deletion the edited molecule reads `ATGAATCCT...` from one base upstream of where it read before, with the coding sequence identical to the last base. The start remains, so `start_retained_variant` is the term the sequence supports; ClinVar's `g.36265247del` describes the same molecule as an edit that does not touch the codon at all. vep-rs keeps `start_retained_variant` and drops `start_lost`.

On the structural deletion the codon's three bases, 26966389-26966391, lie inside the deleted interval 26965975-26966647. Nothing of the codon is left to be retained, so `start_lost` is the term the coordinates support, and vep-rs drops `start_retained_variant`.

In both cases vep-rs emits the member of the pair that the record's own sequence or coordinates support, and every other term of the row unchanged.

Both rows also carry `frameshift_variant`, which both engines derive from the deletion of one base at the annotated coding coordinates. Ensembl VEP's sequence rule counts deleted bases and does not ask whether the frame re-forms upstream; on the structural row its structural rule fires on the span. That term is identical on both sides, is not part of the class, and is not defended here.

### 3.6 Try it yourself

The sequence deletion has two VCF spellings of one edit, because the deleted A is part of a two-base run. Annotate both with Ensembl VEP 116.2 (GRCh37 cache and FASTA, the command of [the closing section](reproducing.md)):

- The left-aligned `21 36265245 . AT A` above gives ENST00000486278 `frameshift_variant,start_lost,start_retained_variant`.
- The right-aligned `21 36265246 . TT T` (ClinVar's own `g.36265247del`) gives `5_prime_UTR_variant` with no start term at all; vep-rs writes the same `5_prime_UTR_variant` for that spelling.

One molecule, and Ensembl VEP calls its start codon both lost and untouched depending on which base of the run the record names.

The structural deletion has an explicit spelling too: write the same 673 deleted bases out. REF is the reference sequence of 21:26965974-26966647, 674 bases (`samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 21:26965974-26966647`); ALT is its first base, `G`.

- The explicit spelling: Ensembl VEP evaluates the edited sequence and writes ENST00000284987 `start_lost,5_prime_UTR_variant`, with no `start_retained_variant`. Its sequence-variant rules also drop `feature_truncation`, a structural term, and `frameshift_variant`, because the 417 deleted coding bases are a multiple of three; the probe is about the start pair. vep-rs writes the same `start_lost,5_prime_UTR_variant` row for the explicit spelling.
- The `<DEL>` spelling of the same interval receives both start terms.

### 3.7 How the comparator treats this class

A pair matches when Ensembl VEP's set contains both `start_lost` and `start_retained_variant` and vep-rs's set is the same set minus exactly the erroneous member for the allele kind: `start_lost` on a sequence allele, `start_retained_variant` on a structural one. The test is `_is_start_cooccurrence_swap` in `scripts/concordance/compare_vep_outputs.py`. An allele written as bases or `-` is a sequence allele; a Sequence Ontology class name such as `deletion` is structural. A vep-rs set that adds any other term, or drops the other member, is not matched and stays charged. Both members of a matching pair are set aside.

The class is fully masked. Counts, from the SNP/indel summary's `start_cooccurrence_swap` line:

- ClinVar GRCh37: 261 pairs.
- ClinVar GRCh38: 1,466 pairs.
- gnomAD v2.1.1 chr21: 3 pairs.
- gnomAD v4.1 chr21: 127 pairs.
- 1000 Genomes high-coverage chr21: 1 pair.

Its reach into the structural-variant sets is class 4.

## 4. Class 3 on the structural-variant sets

**In short.** Ensembl VEP writes `start_lost` and `start_retained_variant` together on the structural deletions and gains over a start codon in the structural-variant files, and on one sequence allele there, the RUNX1 deletion of class 3. vep-rs keeps one term, `start_lost` on the deletions and gains and `start_retained_variant` on the RUNX1 record, which is right for the deletions because their interval contains the codon ([3.5](#35-why-vep-rss-row-is-the-correct-one)). For the ten gains, Ensembl VEP's `start_retained_variant` comes from a test that answers "retained" for every structural allele without reading the record, so vep-rs keeps the member computed from the record, and the comparator sets both members aside rather than deciding the biology.

### 4.1 What a reader needs to know first

The structural-variant sets are scored by a different script from the SNP/indel datasets, and structural deletions and gains over a start codon reach the structural branch of class 3. Nothing new happens in either engine; the same predicate is applied by the structural-variant comparator. The terms are those of [3.1](#31-what-a-reader-needs-to-know-first).

### 4.2 The record

The structural deletion of [3.2](#32-the-record), `gnomAD-SV_v3_DEL_chr21_88237371` on ENST00000284987 (gnomAD-SV v4.1 chr21, GRCh38), and the RUNX1 deletion of [3.2](#32-the-record), which the ClinVar structural-variant file `clinvar_sv_chr21.vcf.gz` (GRCh37) also carries.

### 4.3 What each engine writes

The RUNX1 deletion's rows in `clinvar_sv_chr21.vcf.gz`, which carries the record under the same ID, are those of [3.3](#33-what-each-engine-writes). The `<CN2>` gain `synth_cn2_0016` of [4.5](#45-why-vep-rss-row-is-the-correct-one-for-the-deletions-and-what-the-structural-gains-show), on ENST00000647188 (DYRK1A, forward strand), one row from Ensembl VEP and one from vep-rs, the differing cells in bold:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `synth_cn2_0016` | `synth_cn2_0016` |
| Location | `21:37367280-37521417` | `21:37367280-37521417` |
| Allele | `duplication` | `duplication` |
| Gene | `ENSG00000157540` | `ENSG00000157540` |
| Feature | `ENST00000647188` | `ENST00000647188` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | **`start_lost,feature_elongation,start_retained_variant,coding_sequence_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant`** | **`start_lost,feature_elongation,coding_sequence_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant`** |
| cDNA_position | `508-12083` | `508-12083` |
| CDS_position | `-` | `-` |
| Protein_position | `-` | `-` |
| Amino_acids | `-` | `-` |
| Codons | `-` | `-` |
| Existing_variation | `-` | `-` |
| Extra | `IMPACT=HIGH;STRAND=1;OverlapBP=154138;OverlapPC=96.59` | `IMPACT=HIGH;STRAND=1;OverlapBP=154138;OverlapPC=96.59` |

The record's other four pairs, on ENST00000643624, ENST00000644942, ENST00000647504 and ENST00001131494, differ the same way; its other 29 rows are identical on both sides.

### 4.4 Why Ensembl VEP writes that

As [3.4](#34-why-ensembl-vep-writes-that) for the deletions. On a structural gain, the same two structural tests apply:

- `start_lost` from the overlap of the stated span with the codon;
- `start_retained_variant` from the retention test that returns 0 for every structural allele.

### 4.5 Why vep-rs's row is the correct one for the deletions, and what the structural gains show

As [3.5](#35-why-vep-rss-row-is-the-correct-one) for the 30 structural deletions of the class (2 on GRCh37, 28 on GRCh38; 29 pooled, one deletion recurring in two GRCh38 files), whose interval contains the codon.

The ten structural-gain pairs of the class, on five GRCh38 records, all over a start codon, are:

- an Alu insertion on two transcripts (gnomAD-SV v4.1 chr21, `21:39660415-39660694`);
- an insertion on one transcript (`05_symbolic_del_ins`, `21:37420264-37526076`, ENST00000646548);
- the `<CN2>` copy gain `synth_cn2_0016` on five transcripts (`07_cnv_repeat`, `21:37367280-37521417`);
- two duplications on one transcript each (gnomAD-SV v4.1 chr21, `21:30487315-30487396` and `21:46142001-46156000`).

For them the deletion argument does not apply, and a symbolic gain does not say what becomes of the codon.

Ensembl VEP's structural `start_lost` is the overlap of the stated span with the codon's three coordinates ([3.4](#34-why-ensembl-vep-writes-that)). Its `start_retained_variant` on a structural allele comes from a test that returns "retained" for every structural allele without reading the sequence. What the record does decide, then, is which member came from evidence: the retained member is the one emitted from nothing and is the contradiction, and vep-rs keeps the member computed from the record.

That is a statement about which test read the record, not about the biology. Whether a gain over the codon leaves the start in place is not decidable from a symbolic record (a tandem duplication keeps its first copy of the codon). So for these ten pairs the page claims only that vep-rs's row drops the term emitted without evidence, and the comparator sets both members aside rather than deciding it.

### 4.6 Try it yourself

For the deletions, [3.6](#36-try-it-yourself). The structural gains are the records named in [4.5](#45-why-vep-rss-row-is-the-correct-one-for-the-deletions-and-what-the-structural-gains-show); the `<CN2>` record is `21 37367279 synth_cn2_0016 T <CN2>` with `END=37521417`.

### 4.7 How the comparator treats this class

The same `filter_registry_swap_pairs` as class 2 applies the class 3 rule to the pairs of the structural-variant comparison.

Counts: the class 3 rule matches 3 pairs on GRCh37 and 38 on GRCh38 per file; pooled per assembly, 3 and 37.

- GRCh37 has 3 pairs: one each in `11_multi_allelic.vcf.gz` (`synth_multi_mix_0476`, a deletion), `clinvar_sv_chr21.vcf.gz` (where it is the RUNX1 deletion of [3.2](#32-the-record)) and `gnomad_sv_v2.1_chr21.vcf.gz` (`gnomAD-SV_v2.1_DEL_21_179957`).
- GRCh38 has 38 pairs per file: `gnomad_sv_v4.1_chr21.vcf.gz` 30, `07_cnv_repeat.vcf.gz` 5, and one each in `05_symbolic_del_ins.vcf.gz`, `12_complex_imprecise.vcf.gz` and `1kg_sv_chr21.vcf.gz`. Pooled, 37: the deletion `21:30755356-30755423` on ENST00000335093 is `HGSV_242982` in the 1000 Genomes file and `gnomAD-SV_v3_DEL_chr21_d911172a` in the gnomAD-SV file, one key counted once.

The tandem-repeat file, `07_cnv_repeat.vcf.gz`, carries 700 `<CNV>` and `<CN0>` to `<CN9>` copy-number records beside its 200 `<CNV:TR>` records, and Ensembl VEP expands only the `<CNV:TR>` ones ([12.4](cnv-tr.md#124-why-ensembl-vep-writes-that)). Its five class 3 pairs are the `<CN2>` copy gain `synth_cn2_0016`, Location 21:37367280-37521417, on the five transcripts whose start codon it covers, ENST00000647188 among them.

Of the 41 per-file pairs of both assemblies (40 pooled), 30 are structural deletions whose interval contains the codon (29 pooled) and 10 are structural gains, where vep-rs keeps `start_lost`. The RUNX1 record is the one sequence allele, where it keeps `start_retained_variant`.

Fully masked; the report counts them as `excluded_registry_swap_rust` and `excluded_registry_swap_perl`, split by rule in `excluded_registry_swap_by_bucket`.
