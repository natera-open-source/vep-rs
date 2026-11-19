# The inserted bases of a shifted insertion's HGVSc on a RefSeq model edited to the genome

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.2 describes an insertion on a transcript model carrying an `_rna_edit` attribute (`BAM_EDIT=OK`) at its 3'-shifted coordinates but spells the inserted bases as the record wrote them, rotating them only on a coding row under `--shift_3prime 1`. `21 35467700 T TGA` on `NM_006933.7` is written `c.204_205insGA`, where the shifted insertion is `c.204_205insAG`, the notation vep-rs writes.

**Why this is a defect.** An insertion's notation names the two bases it falls between and the bases inserted there. The record inserts `GA` after c.203 of `CTGGC`, so the cDNA reads `CTG AG GC` under either description of it. Written as given it is `c.203_204insGA`; shifted one base 3', past the `G` at c.204 that it repeats, it is `c.204_205insAG`, that `G` moved from the front of the inserted sequence to its back. `c.204_205insGA` reads `CTG GA GC`, a sequence the variant does not produce. The 3' rule moves the description and never the bases, and Ensembl VEP writes the rotated bases itself on an unedited model and, on an edited one, on a coding row under `--shift_3prime 1`.

**What vep-rs does instead.** vep-rs spells a shifted insertion's bases as the shift rotates them on every model, with or without `--shift_3prime 1`, so an edited model's HGVSc reads as its unedited sibling's does.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, `ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1441-1444 and 1474-1483, `hgvs_transcript`: the 3' shift is computed, and `hgvs_variant_notation` receives the shift hash's `hgvs_allele_string`, the inserted bases rotated along the shift, with the shifted coordinates.
- The same file, lines 1508-1513: on a transcript with an `_rna_edit` attribute the notation's reference and alternate alleles are overwritten with the allele objects' `feature_seq`, the alleles as the record wrote them, read on the transcript's strand.
- The same file, lines 813-820 and 1330-1348, `codon` and `shift_feature_seqs`: `feature_seq` is rotated when the allele carries a shift hash as its codons are computed, which only a coding row does (line 815 returns before the rotation without a translation position). `TranscriptVariation.pm` line 134 attaches that hash under `--shift_3prime 1` alone (`Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm` line 162), and the hash `hgvs_transcript` builds for itself is deleted at line 1449 before the alleles are read.
- The same file, lines 1774-1778, `hgvs_protein`: the other caller of `shift_feature_seqs`, run after `hgvs_transcript` has read the alleles (`Bio/EnsEMBL/VEP/OutputFactory.pm` lines 1751 and 1760).

**What correcting it would change.** The notation would carry the rotated bases `hgvs_variant_notation` was given, as it does on an unedited model, and `c.204_205insGA` would read `c.204_205insAG`.

**The class it produces.** `refseq_edited_model_shifted_insertion_bases`: a differing `HGVSc` on an insertion's row carrying `BAM_EDIT=OK` whose `HGVS_OFFSET` is one non-zero integer on both sides. The two notations are insertions of one transcript at one coordinate pair, Ensembl VEP's bases the row's `Allele` read on the transcript's strand and vep-rs's those bases turned towards the transcript's 3' end by the offset modulo their length. Built into the format comparator and applied once per row. A duplication prints no bases and never differs. No committed corpus record reaches the shape, and the concordance figures score consequence sets and never see the field.

Under `--shift_3prime 1` both engines spell the rotated bases on a coding row, so the class has no coding row there. A non-coding row keeps the written bases on Ensembl VEP under the flag too, and the class applies: `21 35471002 T TGA`, a 3' UTR insertion on the same model, reads `c.*1349_*1350insGA` against `c.*1349_*1350insAG`. Where Ensembl VEP miscounts the coding row's rotation, its HGVSc follows its `Allele` column ([Allele on a coding row after a shift](../output/vep-and-tab.md#allele-on-a-coding-row-after-a-shift)).

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. `GA` inserted after c.203 of NM_006933.7

**In short.** `21 35467700 ins_cds T TGA` on the `GRCh37-refseq` golden corpus's cache under `--refseq --hgvs` inserts `GA` after c.203 of `NM_006933.7` (SLC5A3), whose three `_rna_edit` attributes are single-base substitutions far from the site, so the model equals the genome across the shift. Both engines shift the insertion one base (`HGVS_OFFSET=1`) and agree on every other column; Ensembl VEP writes `HGVSc=NM_006933.7:c.204_205insGA`, vep-rs `HGVSc=NM_006933.7:c.204_205insAG`. Under `--refseq --hgvs --shift_3prime 1` both write `c.204_205insAG`, and the Ensembl model of the same coding sequence, `ENST00000608209`, reads `c.204_205insAG` on both engines in both runs.

### 1.1 What a reader needs to know first

A converted RefSeq or merged cache carries the attributes Ensembl VEP writes about each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). A model the alignment corrected prints `BAM_EDIT=OK` on its rows, and its `_rna_edit` attributes name each difference as `<start> <end> <bases>`, a substitution when the span and the bases have one length. Under `--hgvs` both engines shift an insertion or deletion 3' along the transcript before writing its notation, `HGVS_OFFSET` prints the distance moved, and an insertion moved past a base it repeats carries that base from the front of its inserted sequence to the back.

### 1.2 The record

The record is hand-written on the `GRCh37-refseq` golden corpus's cache and reference (`tests/golden/116/GRCh37-refseq/`), an insertion of `GA` after `21:35467700`:

```
21	35467700	ins_cds	T	TGA	.	.	.
```

`NM_006933.7` is a forward-strand model whose cDNA 707-710 are `TGGC`, c.203-c.206, equal to the genome at `21:35467700-35467703`; its `_rna_edit` attributes are `652 652 G`, `2200 2200 A` and `3049 3049 T`, single-base substitutions none of which lies within 50 bases of the site. The inserted `G` repeats c.204, so both engines shift the insertion one base 3' to between c.204 and c.205, where the inserted sequence reads `AG`.

### 1.3 What each engine writes

The model's row under `--refseq --hgvs`, default format. Ensembl VEP:

```
ins_cds	21:35467700-35467701	GA	6526	NM_006933.7	Transcript	frameshift_variant	707-708	203-204	68	L/LX	ctg/ctGAg	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK;HGVSc=NM_006933.7:c.204_205insGA;HGVSp=NP_008864.4:p.Ala69ArgfsTer64;HGVS_OFFSET=1
```

vep-rs:

```
ins_cds	21:35467700-35467701	GA	6526	NM_006933.7	Transcript	frameshift_variant	707-708	203-204	68	L/LX	ctg/ctGAg	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK;HGVSc=NM_006933.7:c.204_205insAG;HGVSp=NP_008864.4:p.Ala69ArgfsTer64;HGVS_OFFSET=1
```

Every column but the HGVSc agrees, the HGVSp and the offset included. Under `--refseq --hgvs --shift_3prime 1` both engines write the row at the shifted span, `Allele` `AG`, `cDNA_position` `708-709`, `Codons` `-/AG`, and `HGVSc=NM_006933.7:c.204_205insAG`.

On a non-coding row of the same model both engines print the `Allele` column as written under the flag, Ensembl VEP's HGVSc keeps those bases, and the two notations differ there as in the plain run. `21 35471002 T TGA`, inserted in the 3' UTR, reads `HGVSc=NM_006933.7:c.*1349_*1350insGA;HGVS_OFFSET=1` on Ensembl VEP and `c.*1349_*1350insAG` on vep-rs in both runs.

### 1.4 Why Ensembl VEP writes that

`hgvs_transcript` shifts the insertion and hands `hgvs_variant_notation` the shift hash's `hgvs_allele_string`, the inserted bases rotated along the shift, so the notation it builds is `c.204_205insAG`. On a transcript carrying an `_rna_edit` attribute it then overwrites the notation's reference and alternate alleles with the allele objects' `feature_seq`. That overwrite serves the models whose sequence differs from the genome, where a substitution's reference base is the model's and not the genome's.

`feature_seq` is the allele as the record wrote it, read on the transcript's strand. It is rotated by `shift_feature_seqs`, which `codon` runs when the allele carries a shift hash and the row has a translation position: a coding row under `--shift_3prime 1`, where the hash is attached as the transcript variation is built. `hgvs_protein` runs the rotation too, but after `hgvs_transcript` has read the alleles. Without the flag `hgvs_transcript` builds its own hash for the shift and deletes it before the alleles are read, so the overwrite restores the unrotated `GA` at the shifted coordinates; on a non-coding row it does so under the flag as well.

Source (`ensembl-variation` at release/116, `ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1441-1449: `_return_3prime(1)`, the rotated `hgvs_allele_string` taken for an insertion, and the hash deleted when `hgvs_transcript` built it.
- The same file, lines 1474-1483: `hgvs_variant_notation` on the rotated allele at the shifted coordinates.
- The same file, lines 1508-1513: the alleles overwritten from `feature_seq` on a model with an `_rna_edit` attribute.
- The same file, lines 798-820 and 1330-1348: `codon` runs `shift_feature_seqs` when the allele carries a shift hash and the row has a translation position (line 815), rotating `feature_seq` by the shift. `Bio/EnsEMBL/Variation/TranscriptVariation.pm`, lines 132-135, attaches the hash under `--shift_3prime 1`, which `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, line 162, passes as `-no_shift`.
- The same file, lines 1774-1778: `hgvs_protein` runs `shift_feature_seqs` as well, after `hgvs_transcript` has read the alleles (`Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1751 and 1760).
- `Bio/EnsEMBL/Variation/VariationFeatureOverlapAllele.pm`, lines 246-264, `feature_seq`: the stored allele, complemented onto the transcript's strand when the record's strand differs.

### 1.5 Why vep-rs's row is the correct one

The model's cDNA at c.202-c.206 is `CTGGC`, the genome's bases at `21:35467699-35467703`, since no edit lies in the window. Inserting `GA` after c.203 gives `CT GA GGC`, and the inserted `G` equals c.204, so the insertion moves one base 3' and its sequence becomes `AG`: `CTG AG GC`, the same cDNA. `c.204_205insAG` describes that sequence at its most 3' position; `c.204_205insGA` describes `CTG GA GC`, a cDNA the variant does not produce, and no shift of the written insertion reaches it. The unedited Ensembl model of the same coding sequence, `ENST00000608209`, is written `c.204_205insAG` by both engines in both runs, and NM_006933.7 itself under `--shift_3prime 1`; vep-rs gives the edited model the notation both engines give everywhere else.

### 1.6 Try it yourself

```
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > ins.vcf
printf '21\t35467700\tins_cds\tT\tTGA\t.\t.\t.\n' >> ins.vcf
printf '21\t35471002\tutr3_ins\tT\tTGA\t.\t.\t.\n' >> ins.vcf
vep -i ins.vcf -o out.txt --offline --cache \
  --dir_cache tests/golden/116/GRCh37-refseq/cache_root --cache_version 116 \
  --species homo_sapiens --grch37 --refseq --hgvs \
  --fasta tests/golden/116/GRCh37-refseq/reference.fa.gz --force_overwrite --no_stats --quiet
grep NM_006933.7 out.txt | grep -o 'HGVSc=[^;]*'
```

| Run | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `--refseq --hgvs` | `c.204_205insGA` | `c.204_205insAG` |
| `--refseq --hgvs --shift_3prime 1` | `c.204_205insAG` | `c.204_205insAG` |
| `--hgvs` on `ENST00000608209`, either run | `c.204_205insAG` | `c.204_205insAG` |
| `21 35471002 T TGA` (3' UTR), either run | `c.*1349_*1350insGA` | `c.*1349_*1350insAG` |

The Ensembl rows come from the `GRCh37` golden corpus's cache (`tests/golden/116/GRCh37/json_cache`) with the same reference FASTA.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies `refseq_edited_model_shifted_insertion_bases` when a comparison case names it in `masks`. A differing `HGVSc` (the JSON `hgvsc`) is accepted on an insertion's row carrying `BAM_EDIT=OK` whose `HGVS_OFFSET` is one non-zero integer on both sides. Both notations must be insertions of the row's transcript, in one numbering, at one coordinate pair. Ensembl VEP's bases must be the row's `Allele` read on the transcript's strand: as written where `STRAND` is 1, reverse-complemented where it is -1, either reading on a row without `STRAND`. vep-rs's bases must be those bases turned towards the transcript's 3' end, the first bases moved to the end, by the offset modulo their length.

So `c.204_205insGA` against `c.204_205insAG` under `HGVS_OFFSET=1` compares equal, and other coordinates, other bases, a duplication, a longer insertion, a missing or differing offset, or the two sides reversed do not. Each accepted row is counted once in every format, under `--hgvs` alone and, on the non-coding rows that keep the written bases, under `--shift_3prime 1`. No committed record is an insertion the shift rotates on an edited model, so no comparison case or corpus manifest names the class.
