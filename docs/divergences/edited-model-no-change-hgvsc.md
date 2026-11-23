# An HGVSc naming a change of a base to itself on an edited RefSeq model

**The defect in one sentence.** On a RefSeq transcript model whose alignment to the genome substituted a base, Ensembl VEP 116.3 writes, under `--hgvs`, an HGVSc for an alternate allele equal to the model's own base at that position: `NR_048535.1:n.2023T>T`, a substitution of a base by itself. vep-rs writes no HGVSc for such an allele.

**Why this is a defect.** An HGVS description names a difference between a sequence and its reference. An allele equal to the reference is no difference: the nomenclature has no `T>T` form, and Ensembl VEP's own notation builder returns nothing when the two alleles it compares are equal. On an edited model the builder compares the alternate with the genome's base, which differs, and only afterwards are both alleles respelled against the model, where they agree. The row's other cells already say so: `USED_REF` is the allele itself.

**What vep-rs does instead.** vep-rs compares the alternate allele with the reference allele the row was analysed against, the model's edited base when the transcript's reference is in force, and writes no HGVSc when they are equal, as it writes none for any allele equal to its reference. Every other cell of the row agrees with Ensembl VEP's, and an alternate that differs from the edited base gets the same notation in both engines.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, `ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/CacheDir.pm`, lines 411-416: a cache built with a BAM switches the transcript's reference on unless `--use_given_ref` is given. `Bio/EnsEMBL/Variation/TranscriptVariation.pm`, lines 441-468, `_get_ref_allele`: with it on, the reference allele is cut from the model's edited cDNA, so at an edited position it is the edited base, the `USED_REF` cell, and not the genome's.
- `Bio/EnsEMBL/Variation/VariationFeatureOverlap.pm`, line 481: an alternate allele is dropped only when it equals the given reference, so an alternate equal to the transcript-derived reference keeps its row.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1474-1483, `hgvs_transcript`: `hgvs_variant_notation` is handed the alternate allele and the genomic slice's sequence. `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 526-532: the reference allele is cut from that sequence, the genome's base, and the notation is abandoned only when the two are equal, so a substitution is built here.
- `TranscriptVariationAllele.pm`, lines 1508-1514: on a transcript with an `_rna_edit` attribute, the notation's reference and alternate are replaced by `feature_seq` of the reference allele and of this allele, both the edited base. `Utils/Sequence.pm`, lines 654-657, `format_hgvs_string`: a substitution prints `<position><reference>><alternate>`, so the base twice.

**What correcting it would change.** Comparing the alleles after they are respelled, or comparing the alternate with the reference the row was analysed against, would give such an allele no HGVSc, as vep-rs writes, or the nomenclature's `n.2023=`. Every other allele's notation would be unchanged.

**The class it produces.** `edited_model_no_change_hgvsc`: under `--hgvs`, the `HGVSc` of a row on a model edited to the genome (`BAM_EDIT=OK`) by an edit that replaces bases, whose alternate allele equals the edited base. The edit is a substitution, or a span replaced by another number of bases, whose paired bases are read one to one. Ensembl VEP writes the base twice around `>`; vep-rs writes no notation. The class is documented only: the committed RefSeq and merged corpora hold no record whose alternate is an edited base, so no manifest names the class. The concordance figures score consequence sets and never see the field.

## 1. The alternate allele that is the model's edited base

**In short.** `NR_048535.1`, a non-coding RefSeq model on the reverse strand, carries the edit `2023 2023 T`: at cDNA position 2023 the model reads `T` where the genome reads `G`, so on the forward strand the model's base is `A` against the genome's `C`. Two substitutions at that base, `C` to `A` and `C` to `T`, are run on the GRCh38 RefSeq corpus's cache under `--refseq --hgvs`. Both engines write `GIVEN_REF=C;USED_REF=A;BAM_EDIT=OK` on both rows and `HGVSc=NR_048535.1:n.2023T>A` on the second. On the first, the allele that is the edited base, Ensembl VEP writes `HGVSc=NR_048535.1:n.2023T>T` and vep-rs writes no HGVSc.

### 1.1 What a reader needs to know first

A converted RefSeq cache carries the edits Ensembl VEP applied to each model from its alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). A substitution edit, `<position> <position> <base>`, makes the model's base differ from the genome's at one cDNA position. A run on such a cache uses the transcript's reference by default: each row's reference allele is read from the model, printed as `USED_REF` beside the input's `GIVEN_REF`, and an alternate allele is compared against it. Under `--hgvs` the HGVSc of a substitution on a non-coding model is `n.<position><reference>><alternate>`, both bases on the transcript's strand.

### 1.2 The record

Two substitutions at the edited base, against the GRCh38 RefSeq corpus's cache and FASTA:

```
21	42849245	edit_base	C	A	.	.	.
21	42849245	other_alt	C	T	.	.	.
```

The model's edit, from the cache record of `NR_048535.1`: `_rna_edit` `2023 2023 T`. Both engines run `--refseq --hgvs` with the corpus's FASTA.

### 1.3 What each engine writes

The `NR_048535.1` rows in the default format. Ensembl VEP 116.3:

```
edit_base	21:42849245	A	10785	NR_048535.1	Transcript	non_coding_transcript_exon_variant	2023	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;GIVEN_REF=C;USED_REF=A;BAM_EDIT=OK;HGVSc=NR_048535.1:n.2023T>T
other_alt	21:42849245	T	10785	NR_048535.1	Transcript	non_coding_transcript_exon_variant	2023	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;GIVEN_REF=C;USED_REF=A;BAM_EDIT=OK;HGVSc=NR_048535.1:n.2023T>A
```

vep-rs:

```
edit_base	21:42849245	A	10785	NR_048535.1	Transcript	non_coding_transcript_exon_variant	2023	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;GIVEN_REF=C;USED_REF=A;BAM_EDIT=OK
other_alt	21:42849245	T	10785	NR_048535.1	Transcript	non_coding_transcript_exon_variant	2023	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;GIVEN_REF=C;USED_REF=A;BAM_EDIT=OK;HGVSc=NR_048535.1:n.2023T>A
```

The tab output carries the same values in its `HGVSc` column, `-` on vep-rs's first row; the VCF output in the `HGVSc` field of `CSQ`, empty there; and the JSON output under `hgvsc`, which vep-rs's first entry lacks. Every other row of the two records, the eleven unedited models, is byte for byte the same in both engines. Under `--use_given_ref`, which keeps the input's reference, both engines write `n.2023G>T` and `n.2023G>A`, the genome's base first.

### 1.4 Why Ensembl VEP writes that

For `edit_base` the reference allele is read from the model: `A` on the forward strand, so `USED_REF=A`, and the alternate `A` equals it. The notation builder does not compare those two. It is handed the alternate, `T` on the transcript's strand, and the genomic slice, whose base at that position is `G`, and finds them different, so it builds a substitution `G>T`. Because the model carries an `_rna_edit` attribute, the reference and alternate of the notation are then replaced by the two alleles' sequences relative to the transcript, `T` and `T`, and the formatter prints them around `>`.

For `other_alt` the alternate is `A` on the transcript's strand, the respelled reference is `T`, and `n.2023T>A` is the same in both engines.

### 1.5 Why vep-rs's row is the correct one

The row says the reference used is `A` and the allele is `A`: there is no change on the transcript to describe. The nomenclature describes a difference, and Ensembl VEP's own builder writes none when the alleles it compares agree; here they are compared against the wrong sequence, the genome's, which the model does not carry at that base. vep-rs compares the allele with the reference the row was analysed against and writes no notation, the same treatment every allele equal to its reference gets. The alternate that does change the base, `other_alt`, gets `n.2023T>A` in both engines.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --cache --dir_cache <corpus>/cache_root \
  --cache_version 116 --species homo_sapiens --assembly GRCh38 --refseq --hgvs \
  --fasta <corpus>/reference.fa.gz
```

vep-rs writes the rows of 1.3 with no HGVSc on the `edit_base` row. Ensembl VEP 116.3 writes `HGVSc=NR_048535.1:n.2023T>T` there. Adding `--use_given_ref` makes the engines agree, `n.2023G>T` on both.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/edited_model_no_change_hgvsc.toml` holds the two records and both engines' `HGVSc` cells. The class covers the `HGVSc` of a row on a substitution-edited model whose alternate allele equals the edited base, under `--hgvs` with the transcript's reference in force. No committed dataset reaches the shape, so the class is documented only and the format comparator counts no cell under it.
