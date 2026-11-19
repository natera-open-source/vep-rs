# The HGVSc of a RefSeq model whose alignment inserted or deleted bases

**The defect in one sentence.** On a RefSeq transcript model whose alignment to the genome inserted or deleted bases, Ensembl VEP numbers the HGVSc from the unedited exons. An intronic distance counts from the unedited exon end, an exonic position on an `NM_` or `XM_` model is slid along the genome by the edit length before it is read, and an `NR_`, `XR_` or Ensembl model with such an edit is never corrected. vep-rs numbers every HGVSc position against the edited model, the transcript the identifier names.

**Why this is a defect.** HGVS coordinates on a transcript refer to the transcript's own sequence, which for a RefSeq model is the edited one. The cache records each difference from the genome as an `_rna_edit` attribute and stores the model's coordinate mapper with the edits folded in, so the cDNA numbering of every exon past an insertion or deletion is shifted by the edit's length. The cDNA, CDS and protein position columns of the same row use that mapper and agree between the engines; only the HGVSc is read through a mapper rebuilt without the edits.

**What vep-rs does instead.** vep-rs numbers every HGVSc position against the edited model, the transcript the identifier names, on RefSeq and merged caches alike.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `ensembl-variation` at release/116, `ensembl` at release/116):

- `Bio/EnsEMBL/VEP/Pipeline/DumpVEP/Dumper/Core.pm`, lines 164 and 169, inside `clean_transcript` (line 150): the dumper deletes the transcript's `edits_enabled` flag and its `transcript_mapper`.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1524-1525 and 2755, `_get_cDNA_position`: each end of the allele is turned into a cDNA coordinate by walking the exons and reading each exon's cDNA start and end through `Exon::cdna_start` (`Bio/EnsEMBL/Exon.pm`, lines 383-410). That method rebuilds a mapper at run time, and with the flag gone `TranscriptMapper.pm`'s `_load_mapper` folds no edit in.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1517 and 1594-1624, `get_misalignment_offset`: for an `NM_` or `XM_` model alone, the length changes of the indel edits at or before the allele's cDNA start are summed (the `REFSEQ_OFFSET` column) and added to the genomic position of both ends. An allele whose start has no cDNA position, an intronic one, counts no edit (lines 1608-1610).

**What correcting it would change.** The class vanishes: every HGVSc on an edited model carries the edited numbering, as the cDNA and CDS position columns of the same row already do.

**The class it produces.** `refseq_edited_model_hgvsc_offset`, applied by the format comparator to the `HGVSc` of a row carrying `BAM_EDIT=OK` whose two notations name the same reference and differ in their coordinate tokens alone. For the golden harness it is documented in the manifest of the one corpus run under `--refseq --hgvs`, the GRCh38 RefSeq corpus's `hgvs` flag set (two rows). The concordance figures score consequence sets and never see the field, so nothing is set aside in them.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion over a donor site downstream of a nine-base edit

**In short.** A 31-base deletion over the donor site of a COL18A1 exon lies on two RefSeq models whose alignment deleted nine bases upstream. Ensembl VEP writes `c.3546_3553+23del`, eight exonic bases none of which is deleted and an intronic distance nine short; vep-rs writes `c.3544+2_3544+32del`, the span the genome has, the same shape the unedited sibling model `NM_001379500.1` gets from both engines.

### 1.1 What a reader needs to know first

A converted RefSeq or merged cache carries the attributes Ensembl VEP writes about each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). A model the alignment corrected prints `BAM_EDIT=OK` on its rows, and its `_rna_edit` attributes name each insertion or deletion as `<start> <end> <bases>`. `REFSEQ_OFFSET` prints the summed length change of the edits upstream of the allele. Under `--hgvs` both engines shift an insertion or deletion 3' along the transcript before writing its notation, and `HGVS_OFFSET` prints the distance moved.

### 1.2 The record

The record is `21 45505272 rs2037126010 AGGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG A`, record 34 of `variants.vcf` in the `GRCh38-refseq` golden corpus (`tests/golden/116/GRCh38-refseq/`), a gnomAD v4.1 deletion of `21:45505273-45505303`:

```
21	45505272	rs2037126010	AGGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG	A	.	.	.
```

Two models carry a nine-base deletion edit upstream of it: `NM_030582.4` (`_rna_edit` `3387 3395`) and `NM_130444.3` (`_rna_edit` `4092 4100`). Both engines run the corpus's `hgvs` flag set, `--refseq --hgvs`, and shift the deletion seven bases 3' along the transcript, to `21:45505280-45505310` (`HGVS_OFFSET=7`).

### 1.3 What each engine writes

The two models' rows, default format. Ensembl VEP (the corpus's `expected_hgvs/default.txt.gz`):

```
rs2037126010	21:45505273-45505303	-	80781	NM_030582.4	Transcript	splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant	3548-?	3539-?	1180-?	-	-	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=-9;GIVEN_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;USED_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;BAM_EDIT=OK;HGVSc=NM_030582.4:c.3546_3553+23del;HGVS_OFFSET=7
rs2037126010	21:45505273-45505303	-	80781	NM_130444.3	Transcript	splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant	4253-?	4244-?	1415-?	-	-	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=-9;GIVEN_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;USED_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;BAM_EDIT=OK;HGVSc=NM_130444.3:c.4251_4258+23del;HGVS_OFFSET=7
```

vep-rs:

```
rs2037126010	21:45505273-45505303	-	80781	NM_030582.4	Transcript	splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant	3548-?	3539-?	1180-?	-	-	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=-9;GIVEN_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;USED_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;BAM_EDIT=OK;HGVSc=NM_030582.4:c.3544+2_3544+32del;HGVS_OFFSET=7
rs2037126010	21:45505273-45505303	-	80781	NM_130444.3	Transcript	splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant	4253-?	4244-?	1415-?	-	-	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=-9;GIVEN_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;USED_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;BAM_EDIT=OK;HGVSc=NM_130444.3:c.4249+2_4249+32del;HGVS_OFFSET=7
```

Every column but the HGVSc agrees, `REFSEQ_OFFSET` included: vep-rs keeps Ensembl VEP's value of the column and applies no offset to the notation.

### 1.4 Why Ensembl VEP writes that

`hgvs_transcript` turns each end of the shifted allele into a cDNA coordinate with `_get_cDNA_position`, which walks the transcript's exons and reads each exon's cDNA start and end through `Exon::cdna_start`. That method maps through the transcript's own mapper, built at run time, and the cache dumper deleted the transcript's edit flag and stored mapper, so the rebuilt mapper folds no `_rna_edit` in and numbers every exon by exon lengths alone. The edited numbering survives only in the cached mapper the cDNA, CDS and protein position columns read, which is why those columns agree while the HGVSc does not.

For an `NM_` or `XM_` model, `get_misalignment_offset` then sums the length changes of the indel edits at or before the allele's cDNA start, -9 here, and adds the sum to the genomic position of both ends, a move along the genome by the edit length. An allele whose start has no cDNA position, an intronic one, counts no edit, so an intronic distance is anchored on the unedited exon end. Here the -9 slides the start from `21:45505280` into the exon at `21:45505271`, numbers the exon end 3553 where the edited model's last exon base is 3544, and counts the intronic distance 23 for 32.

Source (`ensembl-vep` at release/116.2, `ensembl-variation` at release/116, `ensembl` at release/116):

- `Bio/EnsEMBL/VEP/Pipeline/DumpVEP/Dumper/Core.pm`, lines 164 and 169: the flag and the mapper deleted.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1524-1525, 2755 and 2881-2905: the cDNA positions read through the rebuilt mapper.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1517 and 1594-1624: the offset, and lines 1608-1610, the intronic start that counts no edit.
- `Bio/EnsEMBL/Exon.pm`, lines 383-410; `Bio/EnsEMBL/Transcript.pm`, lines 2268-2271: the run-time mapper.

### 1.5 Why vep-rs's row is the correct one

The model's exon is `21:45505134-45505278`, edited cDNA 3409-3553 and CDS 3400-3544; the nine bases the edit `3387 3395` removed lie upstream, so the unedited numbering of the exon runs nine higher. The shifted deletion removes `21:45505280-45505310`, every base of it past the exon end `21:45505278` by 2 to 32 bases: `c.3544+2_3544+32del`. Ensembl VEP's `c.3546_3553+23del` names eight exonic bases, c.3546 to c.3553, none of which is deleted, and an intronic run nine bases short. The unedited sibling model `NM_001379500.1` of the same gene, which shares the exon's genomic span, is written by both engines as `c.3013+2_3013+32del`, the shape vep-rs gives the edited models. Both notations describe 31 bases; only vep-rs's places them where the genome has them.

### 1.6 Try it yourself

```
vep -i tests/golden/116/GRCh38-refseq/variants.vcf -o out.txt --offline --cache \
  --dir_cache tests/golden/116/GRCh38-refseq/cache_root --cache_version 116 \
  --species homo_sapiens --assembly GRCh38 --refseq --hgvs \
  --fasta tests/golden/116/GRCh38-refseq/reference.fa.gz --force_overwrite --no_stats --quiet
grep rs2037126010 out.txt | grep -o 'HGVSc=[^;]*'
```

| Model | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `NM_030582.4` | `c.3546_3553+23del` | `c.3544+2_3544+32del` |
| `NM_130444.3` | `c.4251_4258+23del` | `c.4249+2_4249+32del` |
| `NM_001379500.1` (unedited) | `c.3013+2_3013+32del` | `c.3013+2_3013+32del` |

On the `GRCh37-refseq` corpus's cache under `--refseq --hgvs`, `rs1232890876` at `21:38639535` on `NM_001331018.1`, `NM_001331021.1` and `NM_001331022.1` reads `c.56+4A>G` from Ensembl VEP and `c.57+4A>G` from vep-rs, the models carrying the insertion edit `647 646 G`.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies `refseq_edited_model_hgvsc_offset` when a comparison case names it in `masks`. A differing `HGVSc` (the JSON `hgvsc`) on a row carrying `BAM_EDIT=OK` is accepted when both notations name the same reference, the text before the colon, and are equal once every cDNA coordinate token is replaced by one marker. So `c.3546_3553+23del` and `c.3544+2_3544+32del` compare equal, and a changed allele text or reference version does not. Each accepted row is counted once; the RefSeq comparison cases under `--hgvs` expect 2 rows on GRCh38 and 3 on GRCh37.

The golden harness documents the two rows of the GRCh38 RefSeq corpus's `hgvs` flag set in its manifest under the class; the other RefSeq and merged corpora run without `--hgvs` and document none. Every changed row carries `BAM_EDIT=OK` and only its `HGVSc` changes.
