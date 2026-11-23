# The protein notation of an insertion the 3' shift moves onto the bases an edited RefSeq model holds and the genome lacks

**The defect in one sentence.** When the 3' shift of `--hgvs` carries an insertion onto bases that an edited RefSeq model holds and the genome lacks, Ensembl VEP cuts those bases out of the model's coding sequence before it reads the reference peptide of the protein notation. On the five rows of this page Ensembl VEP writes a `delins` of residues 44 to 46 that names the model's residues 47 and 49, and vep-rs writes the `dup` of residues 44 to 46, the model's own residues.

**Why this is a defect.** A protein notation, the `HGVSp` field, is made from a reference peptide, read from the model, and an alternate peptide, read from the model with the change. It names the reference protein's own residues at its positions. The model's protein is translated from its coding sequence as the cache holds it, with the bases the genome lacks in it. Read with those bases cut out, the reference peptide is not the model's protein at the positions the notation prints.

Ensembl VEP's notation on the five rows of this page is a `delins` of residues 44 to 46 by six residues. The two residues it names at 44 and 46 are the model's residues 47 and 49, three residues further along its protein. What the notation puts in their place holds: the six residues, put for positions 44 to 46 of the model's protein, give the protein of the genome that carries the record.

**What vep-rs does instead.** vep-rs reads the reference codons from the model's coding sequence as the cache holds it, at the span the shift moves the insertion to. Its reference peptide is the model's own residues there, residues 43 to 46 on the five rows, and its notation is the `dup` of residues 44 to 46. Every other cell of the five rows is equal in the two engines, except `HGVSc`, which differs by another entry's mechanism (section 1.3).

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `ensembl` at release/116; `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm` unless another file is named):

- In `hgvs_protein`, line 1711 calls `_return_3prime` unless the shifting flag is 0. Lines 1732-1733 take the translation start and end with the shifting offset. The protein notation is made at the shifted point.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, line 548, maps the record's start and end, each moved by the shifting offset, to protein positions. Line 513 maps them to coding positions in the same way.
- `Bio/EnsEMBL/Mapper.pm`, `map_insert`, lines 527 and 541: an insertion's span runs from one past the base before it to one before the base after it. `Bio/EnsEMBL/TranscriptMapper.pm`, lines 528-529, gives the residue positions of that span.
- Lines 1774-1778 run when the shift length is not 0. They delete the allele's stored `peptide`, `codon` and `feature_seq` and call `shift_feature_seqs`. Both peptides are read again at the shifted span, with the inserted bases as the shift leaves them (section 1.1).
- Line 1791 sets the notation's alternate half to the allele's peptide. Its codons come from `_get_alternate_cds`, lines 2409-2410 and 2421: the coding sequence before the span, the allele, and the coding sequence after the span.
- Line 1793 sets the reference half to the reference allele's peptide. In `codon`, line 868 takes the model's coding sequence for it. On a model with `_rna_edit` attributes, lines 872-874 write the reference allele over the span, and line 877 reads the codons after that write. The reference allele of an insertion is empty, so the span's bases are cut out.
- Line 1804 calls `_clip_alleles` on the two halves. On the five rows three reference residues are left against six alternate ones, and `_get_hgvs_protein_type`, lines 2072-2077, gives the notation the type `delins`. `_get_hgvs_protein_format`, line 1991, writes the first reference residue and the start position, an underscore, the last reference residue and the end position, the type and the alternate residues.
- The cut of lines 872-874 is the one the entry [`refseq_edited_model_inserted_span_codons`](refseq-edited-model-inserted-span-codons.md) describes for the reference halves of `Codons` and `Amino_acids` of an insertion written on such bases. Its page cites lines 868, 872, 873 and 877.

**What correcting it would change.** A correction is needed in the reference peptide alone. Read from the model's coding sequence as the cache holds it, without the write of lines 872-874, the reference peptide at the shifted span is the model's own residues there. On the five rows those are residues 43 to 46, the peptide vep-rs reads, and a notation made from it names the model's residues at its positions. The alternate peptide needs no correction: section 1.5 shows that the six residues Ensembl VEP writes are right.

**The class it produces.** None. A class is a named shape of difference that the comparator sets aside and counts under its name (section 1.7), and the comparator holds none for this entry. `refseq_edited_model_inserted_span_hgvsp` is the entry's name in the index, and no comparison counts rows under it. The entry is the `HGVSp` of the row of an insertion that the 3' shift moves onto bases a RefSeq model carrying an `_rna_edit` attribute holds and the genome lacks, inside its coding sequence.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an insertion that the 3' shift of a run under `--hgvs` moves onto bases a RefSeq model holds and the genome lacks, inside its coding sequence. The committed caches of `AaegL5-aedes`, `AgamP4-anopheles`, `GRCm39-mouse-merged` and `GRCm39-mouse-refseq` hold 6, 4, 3 and 3 models with such bases inside their coding sequence, the 16 the tier counts.

No row of the concordance datasets reaches the mechanism. Those datasets are annotated on two caches, one for each assembly, and neither holds a model with bases the genome lacks: none of the 238,823 models of the GRCh37 cache and none of the 669,092 of the GRCh38 cache.

The mechanism is reached on another cache, Ensembl's RefSeq GRCh37 cache, under `--refseq --hgvs`. No concordance figure is measured on that cache. One real record is shown to reach it there, on 5 rows: variation 768677 of the genome-wide ClinVar GRCh37 file, which section 1 works. Other insertions that only the shift moves onto such bases were not counted on GRCh37.

The page of `refseq_edited_model_inserted_span_codons` counts, on the RefSeq GRCh38 cache, 29 ClinVar insertions that only the 3' shift moves onto bases a model holds and the genome lacks, inside its coding sequence. None of the 29 is shown here.

All of this rests on one record. Its five rows are on five models of one gene, each on the forward strand and each the genome plus 9 bases at the same place, run under `--refseq --hgvs` and no other annotation flag. This page shows nothing for any other record, model or set of flags.

## 1. ClinVar variation 768677 on five RefSeq models of one gene

**In short.** On GRCh37, ClinVar variation 768677 is an insertion of 18 bases between 14:94582130 and 14:94582131. Five RefSeq models of one gene hold 9 bases the genome lacks three bases further on, at their coding positions 129 to 137. Under `--refseq --hgvs` the 3' shift moves the insertion onto those 9 bases, and both engines print `HGVS_OFFSET=3` on the five rows. Ensembl VEP 116.3 writes a `delins` of residues 44 to 46 that names there the models' residues 47 and 49; vep-rs writes the `dup` of residues 44 to 46, which names the models' own residues.

### 1.1 What a reader needs to know first

A cDNA position numbers a model's bases from its first base. A coding position numbers them from the first base of its coding sequence, and `CDS_position` prints coding positions. Three coding positions make one residue: residue n is coding positions 3n - 2 to 3n. A mapper pair ties one unbroken stretch of a model's cDNA positions to the genomic positions it lies on.

This page calls a RefSeq model that carries an `_rna_edit` attribute an edited model. On each of the five models here the attribute inserts 9 bases: the genome's exons joined, with that edit applied, are the model's cached sequence. The rows of the five models print `BAM_EDIT=OK`.

Under `--hgvs` both engines shift an insertion 3', as far as it can go without changing the sequence it describes, before they write its notations. They print the distance in bases as `HGVS_OFFSET`. Each step of the shift moves the insertion one base on and puts the first of the inserted bases at their end. This page calls the inserted bases after the shift the turned allele. The row's other cells are read at the point the record writes: the consequence term in `Consequence`, the three position columns, `Codons` and `Amino_acids`.

A protein notation, the `HGVSp` field, is made from a reference peptide and an alternate peptide read over the same residue positions. Before it is written, the residues the two share at their ends are clipped. A `delins` names the first and the last reference residue that are left, each with its position, and then the residues that replace them. A `dup` names a run of reference residues that the added residues repeat.

### 1.2 The record

The record is ClinVar variation 768677, as the genome-wide ClinVar GRCh37 file spells it. It inserts 18 bases between 14:94582130 and 14:94582131, and its rows print the `Location` `14:94582130-94582131`. The 18 bases are one string of 9, twice.

The five rows are the record's rows on `NM_001130080.3`, `NM_001288952.2`, `NM_001288956.2`, `NM_001366993.1` and `NM_001366994.1`, five forward-strand models of one gene. Each holds 9 bases the genome lacks after 14:94582133, three bases past the record's point, at its coding positions 129 to 137.

The worked row is the row of `NM_001130080.3`. The cache holds for that model:

- the `_rna_edit` attribute `249 248` with 9 bases: the edit puts those 9 bases after cDNA position 248;
- a coding start at cDNA 121: coding position n is cDNA position n + 120;
- two mapper pairs that meet three bases after the record's point, cDNA 242-248 on 14:94582127-94582133 and cDNA 258-403 on 14:94582134-94582279;
- between the two pairs, cDNA 249 to 257, coding positions 129 to 137: the 9 bases the genome lacks;
- a coding sequence of 369 bases, 123 codons, whose translation is the protein the cache holds for the model.

Turned by 3, the distance from the record's point to those 9 bases, the record's 18 bases are the model's own 9 added bases, twice. That holds on each of the five models.

Both engines are run under `--refseq --hgvs` and no other annotation flag.

### 1.3 What each engine writes

On the row of `NM_001130080.3` the engines differ in two cells, `HGVSp` and `HGVSc`, and agree in every other. The table gives `HGVSp` and seven of the cells they agree in, each by its positions, lengths and type. The `HGVSc` is another entry's subject, named below the table.

| Cell | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `inframe_insertion` | `inframe_insertion` |
| `cDNA_position` | 245-246 | 245-246 |
| `CDS_position` | 125-126 | 125-126 |
| `Protein_position` | 42 | 42 |
| `Codons` | 3 letters, a slash, 21 letters | the same cell |
| `Amino_acids` | 1 letter, a slash, 7 letters | the same cell |
| `HGVS_OFFSET` | 3 | 3 |
| `HGVSp` | a `delins` of 44 to 46 by 6 residues | a `dup` of 44 to 46 |

Ensembl VEP's `HGVSp` names a residue at position 44 and a residue at position 46, and six residues stand after the word `delins`. The two residues it names at 44 and 46 are not the model's residues 44 and 46. They are the model's residues 47 and 49. vep-rs's `HGVSp` names a residue at position 44 and a residue at position 46 before the word `dup`, and they are the model's own residues 44 and 46.

The record's rows on the other four models show the same. On 5 rows of 5 Ensembl VEP writes a `delins` of residues 44 to 46 by six residues, and vep-rs a `dup` of residues 44 to 46.

The `HGVSc` of the same five rows differs for a cause of its own: the mechanism of [`refseq_edited_model_seam_insertion_hgvsc`](refseq-edited-model-seam-insertion-hgvsc.md), which the rows reach through the shift. That difference is the other entry's, and this page does not work it.

### 1.4 Why Ensembl VEP writes that

The record writes its insertion between 14:94582130 and 14:94582131. On the first mapper pair those are cDNA 245 and 246, coding positions 125 and 126, both in residue 42, which is coding positions 124 to 126. The row's position columns, its consequence term, `Codons` and `Amino_acids` are read there, and the two engines agree on them.

`hgvs_protein` works at the shifted point, three bases on (lines 1711 and 1732-1733). There the insertion stands where the two mapper pairs meet: the base before it is 14:94582133, cDNA 248, and the base after it is 14:94582134, cDNA 258. The model numbers these two neighbours coding positions 128 and 138. An insertion's span runs from one past the base before it to one before the base after it, so the span is coding positions 129 to 137, the model's 9 added bases.

Coding position 129 is the last base of residue 43, and position 137 lies in residue 46, which is coding positions 136 to 138. The span therefore covers residue positions 43 to 46: four residues, coding positions 127 to 138, 12 bases. Both peptides are read again over those four positions (lines 1774-1778).

The alternate peptide is read from the model with the turned allele in place of the span's 9 bases (`_get_alternate_cds`). Of the 12 bases, 3 are left and 18 come in: 21 bases from coding position 127, seven residues. The turned allele is the 9 added bases twice (section 1.2), so its first base is the base of coding position 129. The alternate peptide's first codon is therefore coding positions 127 to 129, the codon of the model's residue 43.

The reference peptide is read over the same four residue positions, from the model's coding sequence with the 9 bases cut out and nothing in their place. The cut is the write of lines 872-874: the reference allele is written over the span, and an insertion's reference allele is empty. Line 877 reads 12 bases from position 127 of the shortened sequence: coding positions 127 and 128, and after them 138 to 147.

Those 12 bases are four codons. The first is coding positions 127, 128 and 138, read across the cut. The other three are coding positions 139 to 147, the codons of the model's residues 47, 48 and 49. After its first residue, the reference peptide holds the model's residues 47 to 49.

The codon read across the cut ends in the base of coding position 138, where the codon of residue 43 ends in that of 129. The two codons read the same residue, so the two peptides open with the same residue. `_clip_alleles` strips that one leading residue, which the two share, and no more.

Three reference residues are left against six alternate ones, and lines 2072-2077 give the notation the type `delins`. The stripped residue stood at position 43, so the notation starts at 44, and it ends at 46, the span's last residue position. Line 1991 writes the first and the last of the three reference residues beside those two positions: the model's residues 47 and 49.

### 1.5 Why vep-rs's row is the correct one

What follows is built from the reference FASTA, the record and the cached models alone, and each step holds on all five models:

- The genome's exons joined, with the model's edit applied, are the model's cached sequence: the model is the genome plus 9 bases after 14:94582133.
- The record's 18 bases, turned by 3, are those 9 bases twice (section 1.2).
- The exons of the genome that carries the record are the model's cached sequence plus 9 bases. The smallest change from the one to the other, placed 3'-most, is the duplication of coding positions 129 to 137.
- Translated from the model's coding start, the duplicated sequence gives the model's protein with residues 44 to 46 repeated.

The repeat starts at the third base of a codon, so the last step is read from the translation and not from the positions alone. One codon of the duplicated sequence is coding positions 136 and 137 followed by the copy of 129, where the model has 136 to 138. It reads the same residue as the model's residue 46.

A protein notation names the reference protein's own residues at its positions. Residues 44 and 46 of each model's protein are the ones vep-rs names, and not the ones Ensembl VEP names. The `dup` of residues 44 to 46 states the change the record makes to the model's protein.

vep-rs reads the reference codons from the model's coding sequence as the cache holds it, on the shifted span. Its reference peptide is the model's residues 43 to 46, the clip leaves three added residues that repeat the three before them, and the notation is the `dup` of 44 to 46.

Ensembl VEP's six residues, put for positions 44 to 46 of the model's protein, give the same protein: they are the model's residues 44 to 46, twice. The alternate half of its notation is right. Its reference half, which names the model's residues 47 and 49 at positions 44 and 46, is the part that does not hold.

### 1.6 Try it yourself

Save the line of variation 768677 from ClinVar's GRCh37 VCF as `record.vcf`, under the two VCF header lines; [Reproducing these examples](reproducing.md) names the file and gives the header lines. vep-rs takes this command:

```
vep -i record.vcf -o out.txt --offline --refseq --hgvs \
  --json_cache /path/to/json_cache/homo_sapiens_refseq/116_GRCh37 \
  --fasta /path/to/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
  --species homo_sapiens --assembly GRCh37
grep NM_001130080.3 out.txt | grep -o 'HGVSp=[^;]*'
```

The cache is the conversion of Ensembl's RefSeq GRCh37 cache of release 116 ([RefSeq and merged caches](../cache-setup.md#refseq-and-merged-caches)). The last line prints vep-rs's notation, a `dup` of residues 44 to 46. Ensembl VEP 116.3 takes `--refseq --hgvs` with the cache options of [Reproducing these examples](reproducing.md), its cache root holding `homo_sapiens_refseq/116_GRCh37`, and prints a `delins` of residues 44 to 46.

The reproduction `scripts/divergences/reproductions/refseq_edited_model_inserted_span_hgvsp.toml` runs the record on both engines. It holds the five rows to their recorded values in six columns: `Feature`, `Consequence`, `CDS_position`, `Protein_position`, `HGVS_OFFSET` and `HGVSp`. It is given Ensembl's RefSeq GRCh37 cache and the directory of the GRCh37 FASTA, and it writes the record's line itself, from the cached model.

### 1.7 How the comparator treats these rows

The format comparator (`tools/vep-formats-compare`) compares the two engines' outputs field by field. A class of the comparator is a named shape of difference that it sets aside: it counts the row under the class's name and does not report the difference. It holds no class for this entry, so nothing is set aside, and it reports the `HGVSp` of each of the five rows as a difference between the engines.

Eight classes the comparator holds for other shapes were tried on the rows, and none sets either of their two differences aside. Each of the five rows was compared alone ten times: with no class applied, with each of the eight applied in turn, and with all eight applied together. All 50 comparisons report both `HGVSp` and `HGVSc` as differences. That includes the comparisons under `refseq_edited_model_seam_insertion_hgvsc`, the class of the entry whose mechanism the rows' `HGVSc` reaches through the shift: that class does not set these rows aside. The eight classes are `start_cooccurrence_swap`, `refseq_edited_model_seam_insertion_hgvsc`, `refseq_edited_model_shifted_insertion_bases`, `refseq_edited_model_hgvsc_offset`, `refseq_edited_model_hgvs_offset`, `refseq_edited_model_tail_term`, `refseq_edited_model_deleted_span_term` and `refseq_edited_model_deleted_span_codons`.

`refseq_edited_model_inserted_span_codons`, the class of the entry of the same cut, takes `Codons` and `Amino_acids` alone. Its page says that a row of a record which only the `--hgvs` shift moves onto such bases is equal in both columns, and it does not name `HGVSp`. The five rows are of that kind. Their `Codons` and `Amino_acids` are read at the written point and are equal in the two engines. The cut shows in the protein notation, whose reference peptide is read again at the shifted span.
