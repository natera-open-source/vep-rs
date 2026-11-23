# The HGVSp of an insertion the 3' shift carries off the bases that open an edited RefSeq model's coding sequence

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP makes the protein notation of an insertion at the point the genome's 3' shift moves it to, and never reads the bases an edited RefSeq model holds where the record is written. Where those bases open the model's coding sequence, its `HGVSp` names the two residues beside the shifted point with a `?` between their positions, and the change the record makes to the protein goes unnamed. vep-rs reads such a row where its record is written, and its `HGVSp` names that change: on the record of section 1, a `dup` of residues 2 to 4.

**Why this is a defect.** An HGVSc is a statement on the sequence of the transcript it names, and the HGVSp follows from it. The record is written between the two genomic neighbours of the model's added bases, and both engines print those bases as the row's `cDNA_position`. Ensembl VEP's notation is read at another point, one the genome's shift reaches because the genome has no position for the added bases. There its `HGVSp` carries the `?` of a lost start, on a record that leaves the 5' UTR and the model's first coding bases as they are.

**What vep-rs does instead.** vep-rs reads such a row where its record is written, on the model's first coding bases, in the frame the cache gives the model. Its `HGVSp` names the change the record makes to the model's protein: on the record of section 1, a `dup` of residues 2 to 4. The row's `HGVSc` and its absent `HGVS_OFFSET` are read at the same place, as the entry [`refseq_edited_model_seam_shift_and_duplication_hgvs`](refseq-edited-model-seam-shift-and-duplication-hgvs.md) describes.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, as the Ensembl VEP 116.3 image holds it):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1441 and 1711: `hgvs_transcript` and `hgvs_protein` each call `_return_3prime` before a notation is made, and the insertion is moved 3' along the genome.
- Lines 1474-1478 of the same module: the coding notation is typed on the transcript's genomic slice at the shifted point, where the two positions are neighbours and the model's added bases do not exist.
- Lines 1608-1617 and 1524-1525 of the same module: the length of an edit that starts at or before the row's first cDNA position is added to both slice positions, which `_get_cDNA_position` then numbers along the genome's exon.
- Lines 2140-2144 of the same module, in `_get_hgvs_peptides`: an insertion's protein notation takes the two residues beside it as its reference.
- Lines 2161-2164 and 2030-2031 of the same module: where `start_lost` holds, the alternate becomes `?` and the type is emptied, and `_get_hgvs_protein_format` prints the reference, the start, `_`, the alternate and the end.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 870-871, 938-941, 961 and 1058-1071: the two start-codon tests that give one allele `start_lost` beside `start_retained_variant`, the subject of [Defect B](defect-b.md).

**What correcting it would change.** The protein notation would be read on the model's sequence at the row's own span. With the start-codon tests of [Defect B](defect-b.md) corrected as well, `start_lost` would lay no `?` over it. The row of section 1 would then carry an `HGVSp` that names its three duplicated residues. No other column of the row is part of this entry.

**The class it produces.** `refseq_edited_model_opening_bases_insertion_hgvsp`: the `HGVSp` of a `Transcript` row of an insertion written on bases an edited RefSeq model holds and the genome lacks, where those bases open the model's coding sequence and the genome's 3' shift moves the insertion off them. The class is documented only. No comparison rule is built for it, because no comparison case and no row the concordance figures score reaches it. Where no other class sets the row aside, a disagreement of the shape is reported as a difference in `HGVSp` and charged to vep-rs.

The counts behind that statement are three. With this reading of the protein notation and without it, vep-rs writes the same output on all 2,168 outputs of the 479 committed runs of 109 corpora. It does so on all 306 comparison cases, 199 of the fixture tier and 107 of the population tier. And every concordance dataset is annotated on a cache whose models hold no base the genome lacks: 0 of 238,823 models on GRCh37 and 0 of 669,092 on GRCh38.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an insertion written on the added bases that open a RefSeq model's coding sequence, in a run under `--hgvs` with the default 3' shift of the notation. The committed caches hold 49 models with bases the genome lacks, in 7 corpora, and on none of the 49 do those bases open the coding sequence. Seven committed flag sets run a RefSeq or merged cache under `--hgvs`.

Real data carries the shape once. The ClinVar files hold 12 records written between the two genomic neighbours of a RefSeq model's added bases, 9 on GRCh37 and 3 on GRCh38. Under `--refseq --hgvs` they give 39 and 10 rows. With this reading vep-rs writes one of the 49 otherwise than without it: the row of section 1, in its `HGVSp` alone.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Eighteen bases inserted on the nine bases that open the coding sequence of `NM_001288959.2`

**In short.** The record with ID `768677` of the ClinVar GRCh37 VCF inserts 18 bases between 14:94582130 and 14:94582131. The forward-strand RefSeq model `NM_001288959.2` holds nine bases at that point that the genome lacks, its coding positions 1 to 9, and the 18 bases are those nine twice. The genome's 3' shift moves the insertion three places. Both engines print `CDS_position` `?-9`, the start unknown because the record's left neighbour lies in the 5' UTR. Under `--hgvs`, Ensembl VEP 116.3 writes the `HGVSp` `p.GlyCys5_?4`: two residue names, the position 5, `_?` and the position 4. vep-rs writes `p.His2_Gly4dup`, a `dup` of residues 2 to 4.

### 1.1 What a reader needs to know first

A RefSeq cache carries each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). Where the alignment found bases in the model that the genome lacks, the model holds them between two bases that are neighbours on the genome. An insertion written between those two neighbours has the added bases as its span on the model, and the row prints them as its `cDNA_position`. `REFSEQ_OFFSET` prints the summed length change of the model's edits that start at or before the row's first cDNA position.

The position columns count that span from the record's two neighbours: it runs from one past the first to one before the second. A neighbour outside the coding sequence has no coding position, and the span's end on that side is printed `?` in `CDS_position` and `Protein_position`.

Under `--hgvs` an insertion is shifted 3' as far as its bases repeat before its notation is written, and `HGVS_OFFSET` prints the places moved. The shift runs along the genome, which has no position for a model's added bases.

HGVS describes a change on the transcript's own sequence, and the protein notation describes the change to its protein. An `ins` names two neighbouring positions and what is put between them. A `dup` names a range that is copied. Both keep every reference base or residue.

### 1.2 The record

The record is written at 14:94582130 in the ClinVar GRCh37 VCF, with the ID `768677`, a reference of one base and an alternate of 19: 18 bases inserted between 14:94582130 and 14:94582131. The 18 bases are one string of nine written twice. The genome lets the insertion move 3' by three bases and no further. On the genome it is no duplication: neither the 18 bases before the written point nor the nine before it equal the allele or its half.

On the RefSeq GRCh37 cache the record has six rows, each on a forward-strand model with one edit that adds nine bases. On `NM_001288959.2` the nine stand at the record's own point, between 14:94582130 and 14:94582131: cDNA positions 97 to 105. The cache places the model's coding start at cDNA 97, so they are its coding positions 1 to 9. The record's 18 bases are the model's nine added bases twice.

The cached coding start, cDNA 97 to 99, holds no initiation codon in the model itself. One stands four bases on, at cDNA 101 to 103, inside the nine added bases. On the five other models the nine bases stand three bases on, after 14:94582133, and their rows are no part of this entry.

Both engines run `--refseq --hgvs` with a reference FASTA.

### 1.3 What each engine writes

The record's row on `NM_001288959.2`. Both engines print `cDNA_position` 97-105, `CDS_position` `?-9`, `Protein_position` `?-3`, empty `Codons` and `Amino_acids`, `BAM_EDIT=OK` and `REFSEQ_OFFSET=9`.

The `?` comes from the record's left neighbour. 14:94582130 is cDNA 96, the last base of the 5' UTR, and has no coding position, so the span's start has none. The right neighbour, 14:94582131, is cDNA 106, coding position 10, and the span ends one before it, at coding position 9 and residue 3. Both engines print `5_prime_UTR_variant` on the row.

Five fields differ:

| Field | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `HGVSp` | `p.GlyCys5_?4`: two residue names, the position 5, `_?` and the position 4 | `p.His2_Gly4dup`, a `dup` of residues 2 to 4 |
| `HGVSc` | an `ins` of the 18 bases as written between coding positions 12 and 13 | `c.4_12dup`, a `dup` of coding positions 4 to 12 |
| `HGVS_OFFSET` | 3 | none |
| `Consequence` | `start_lost,start_retained_variant,5_prime_UTR_variant` | `coding_sequence_variant,5_prime_UTR_variant` |
| `IMPACT` | HIGH | MODIFIER |

The `HGVSp` is this entry. The `HGVSc` and `HGVS_OFFSET` differ as the entry [`refseq_edited_model_seam_shift_and_duplication_hgvs`](refseq-edited-model-seam-shift-and-duplication-hgvs.md) describes: vep-rs reads the row where its record is written, so it reports no places moved, and its `dup` is placed by the 3' rule on the model's sequence. That entry's rows print the span of the added bases in all three position columns. On this row the left neighbour lies outside the coding sequence, so the `CDS_position` starts with `?`. The terms and `IMPACT` differ as [`start_codon_overlap_without_start_codon`](start-codon-overlap-without-start-codon.md) describes: the cached coding start holds no start codon (1.2), so vep-rs writes no start term.

### 1.4 Why Ensembl VEP writes that

Both notations are made after the insertion is shifted, by line 1441 for the `HGVSc` and line 1711 for the `HGVSp`. The genome has no position for the model's nine added bases, so on it the insertion point is followed directly by 14:94582131. The shift moves the insertion three places, to between 14:94582133 and 14:94582134: `HGVS_OFFSET=3`.

Lines 1474-1478 type the coding notation on the genomic slice at that point, where the two positions are neighbours: an `ins`. The row's cDNA start, 97, is not before the edit's start, 97, so lines 1608-1617 add the edit's length, 9, to both positions; it is the row's `REFSEQ_OFFSET`. Lines 1524-1525 number them along the genome's exon, as coding positions 12 and 13. The model's nine bases are kept and the whole allele is added beside them.

The protein notation is read at the same point. `_get_hgvs_peptides` gives an insertion the two residues beside it as its reference, here the model's residues 4 and 5, and the notation starts at 5 and ends at 4. Ensembl VEP gives the row `start_lost`, so lines 2161-2164 make the alternate `?` and empty the type. Lines 2030-2031 then print the reference, the start, `_`, the alternate and the end: the two residue names, 5, `_?` and 4.

The row carries `start_lost` beside `start_retained_variant` because the two start-codon tests of `Utils/VariationEffect.pm` disagree on it. The retention test, lines 1058-1071, puts the allele in place of cDNA 97 to 105 and finds the tail of the rebuilt sequence equal to the coding sequence: not altered, so line 961 gives `start_retained_variant`. The loss test, lines 938-941, reads three bases of the rebuilt sequence at the old start. They are the allele's first three bases, equal to the model's own coding bases 1 to 3, and they are not the literal initiation codon: line 871 gives `start_lost`.

### 1.5 Why vep-rs's row is the correct one

The record's two neighbours, 14:94582130 and 14:94582131, are the model's cDNA positions 96 and 106. Between them lie cDNA 97 to 105, the nine bases the genome has no position for, and both engines print that span as the row's `cDNA_position`.

The model's exons are the genome's with one copy of those nine bases at that point. The record's 18 bases are the nine twice, so the genome that carries the record holds two copies there: the model's sequence with nine bases more. The smallest change from the model to that sequence, placed 3'-most, is a duplication of nine bases, coding positions 4 to 12. Read from the model's coding start, the record's protein is the model's with three residues duplicated, residues 2 to 4 in the frame the cache gives the model. vep-rs's `HGVSp` states that.

Ensembl VEP's `HGVSp` carries the `?` of a lost start. The record leaves the whole 5' UTR and the model's coding bases 1 to 12 as they are, so neither the cached coding start nor the initiation codon four bases on is changed. Its `HGVSc` states another sequence as well: the `ins` keeps the nine added bases and adds 18, the model plus 18 bases where the record gives the model plus nine. Read at the shifted point with the nine bases kept, the allele adds six residues between residues 4 and 5; the record adds three.

### 1.6 Try it yourself

Cut the record from the ClinVar GRCh37 VCF with the file's header, then run the command as vep-rs takes it:

```
zcat clinvar_grch37.vcf.gz | awk '/^#/ || $3 == "768677"' > record.vcf
vep -i record.vcf -o out.txt --offline \
  --json_cache <cache>/homo_sapiens_refseq/116_GRCh37 \
  --fasta Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
  --species homo_sapiens --assembly GRCh37 --refseq --hgvs
grep NM_001288959.2 out.txt | grep -o 'HGVS[^;]*'
```

`<cache>` is the conversion of Ensembl VEP's RefSeq GRCh37 cache ([Cache conversion](../cache-conversion.md)). The last command prints the row's HGVS items. vep-rs prints `c.4_12dup` and `p.His2_Gly4dup`, and no `HGVS_OFFSET`. Ensembl VEP 116.3 takes the same flags with the cache options of [Reproducing these examples](reproducing.md) and its RefSeq cache, and prints the three HGVS values of its column in 1.3.

### 1.7 How the comparator treats this class

Nothing is built for this class. The format comparator (`tools/vep-formats-compare`) has no class of this shape and no test for it, because no comparison case and no row the concordance figures score reaches a row of the shape. The counts stand under "The class it produces" above.

On the row of this section Ensembl VEP writes both start terms and vep-rs neither, the pair of [`start_codon_overlap_without_start_codon`](start-codon-overlap-without-start-codon.md). A comparison that names that class and is given the candidate's cache sets aside the `Consequence` and the `IMPACT` that follows it, and the `HGVSp` stays a difference charged to vep-rs. A comparison case that names no class reports the five fields of the table in 1.3 as differences, charged to vep-rs.
