# The coding reference of a deletion that spans a whole intron of an edited RefSeq model

**The defect in one sentence.** On an edited RefSeq model, for a record whose reference span has both ends on coding bases of the model and holds a whole intron, Ensembl VEP reads the genome's reference allele, intron included, as the coding reference. The reference half of `Codons`, the part before the slash, then holds that allele from the row's first coding position to the end of the row's codon window (section 1.1), intron bases among them. The reference half of `Amino_acids` is the translation of that half. vep-rs reads the reference codons from the model's coding sequence itself, so its reference halves are the model's own codons and residues.

**Why this is a defect.** An edited model is a RefSeq model that carries one or more `_rna_edit` attributes. On such a model Ensembl VEP writes the row's reference allele over the model's coding sequence, at the row's coding span, before it reads the reference codons. Such a write fits only where the allele is as long as the coding span it replaces. A span with a whole intron inside keeps the genome's allele, which is longer than the coding span by the intron. The codon window is then filled from the allele, intron bases among them, in place of the model's codons.

The reference peptide is the translation of that window. Every predicate that compares the reference peptide or the reference codons with the alternate's then starts from residues the protein does not hold, and so does the protein notation. Section 1 shows the outcome on a real record: `stop_lost` on a deletion that touches no base of the model's stop codon.

**What vep-rs does instead.** vep-rs reads the reference codons from the model's coding sequence itself, so its reference halves are the model's own codons and residues at the row's places. It runs ports of the same predicates and of the same notation routine on them. On all fifteen rows of section 2, `cDNA_position`, `CDS_position`, `Protein_position`, `HGVSc` and the alternate halves of `Codons` and `Amino_acids` are equal in the two engines.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at commit `db7394a02523a6e519cb0d89f84b2889bac12737`, the files under `modules/Bio/EnsEMBL/Variation/`):

- `TranscriptVariation.pm`, line 459, `_get_ref_allele`: a BAM-edited cache, a cache whose models were corrected against an alignment ([RefSeq and merged caches](../cache-conversion.md#refseq-and-merged-caches)), switches on `--use_transcript_ref`. Under that flag the model's own bases become the reference allele only on a span of at most two cDNA segments. A span with a whole intron inside has three or more, so the genome's allele stands.
- `TranscriptVariationAllele.pm`, lines 872 to 874, `codon`: for the reference allele of a model with one or more `_rna_edit` attributes, the allele is written over the coding sequence at the row's CDS span.
- The same routine, line 877: the reference codons are the first `$codon_len` bases, the length of the codon window, counted from the row's first codon.
- `Utils/VariationEffect.pm`, lines 1259 to 1261, `stop_lost`: line 1260 accepts where the reference peptide holds a stop and the alternate peptide holds none.
- `Utils/VariationEffect.pm`, lines 1178 and 1185, `inframe_deletion`: the two ways the predicate accepts. Line 1178 accepts where the reference codons begin or end with the alternate codons. Line 1185 reads the two strings after the routine trims them, as section 1.5 shows on the worked row. It accepts where the alternate is empty and the length of the reference is a multiple of three.
- `TranscriptVariationAllele.pm`, lines 1916 and 1926, `_get_hgvs_protein_format`: on a row that passes `stop_lost`, a notation of type `del` or `>` is written with a stop-lost extension.
- `Utils/Constants.pm`, line 752, the definition of `stop_lost`: "A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript".
- The entry [`refseq_edited_model_deleted_span_codons`](refseq-edited-model-deleted-span-codons.md) names the same two places, line 459 and the write over the coding span, for another shape: a deletion over the bases a model's alignment deleted. That entry covers the `Codons` cell and, beside it, the `Amino_acids` cell and one pair of terms.

**What correcting it would change.** This is shown in full on the worked row of section 1. Fed the model's own codon window, the same lines give the row vep-rs writes there. The reference peptide holds no stop, so there is no `stop_lost`. `inframe_deletion` stands in place of `protein_altering_variant`, and the protein notation is a plain deletion with no extension. For the coding terms alone, the check of section 2.1 gives the same result on eleven rows: fed the model's reference peptide, ports of the predicates return the coding terms vep-rs writes.

**The class it produces.** `refseq_edited_model_intron_span_coding_reference`. In [the index](../intended-divergences.md) a class is the name a comparison counts an entry's rows under. The golden harness and the format comparator set the worked row of section 1 aside under it, in the corpus `GRCh38-refseq-intron-span`; section 1.7 gives the shape they read. The datasets behind the concordance figures run on Ensembl caches, which hold no model with an `_rna_edit` attribute, so their comparison has no rule for the shape and sets no row aside. The fifteen rows of this page come from a pairing those datasets do not make: genome-wide ClinVar under `--refseq`.

**How often the shape occurs.** Tier 3 of [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data, as that page defines it, does not carry the shape. The shape is a record whose reference, as the engines parse it, has two or more bases, with both ends on coding bases of an edited model and one or more whole introns inside. A scan of genome-wide ClinVar against the whole RefSeq cache of each assembly finds fifteen rows of the shape, a row being one record on one model. There are 3 on GRCh38, from 3 records, and 12 on GRCh37, from 9 records. Section 2 lists them.

Of the 116,165 transcripts of the committed caches, 65, in 7 corpora, are coding models with an RNA edit and an intron between two bases of their coding sequence, the models that can carry the shape. That count is an upper bound on what can carry the shape, not a count of records that show it. One committed record shows it: the worked record of section 1, which the corpus `GRCh38-refseq-intron-span` holds for this entry. The cuts that page counts hold none, and a corpus record of a kind the cuts can carry is not ordinary data there, which is why the entry is tier 3.

The reference halves of `Codons` and `Amino_acids` differ between the engines on all fifteen. The coding terms of `Consequence`, the impact and the protein notation are equal on ten rows. On the other five a term or the notation differs: the term on 3 of them, the protein notation on all 5 and the impact on none.

**What is not shown.** Four limits bound what this page shows:

- Every one of the fifteen rows removes its span, and one of them, record 2850981, puts 3 bases in place of its 822. Each of the models carries one `_rna_edit`. Nothing is shown for another kind of record or for a model with more edits.
- The rows were written under `--refseq --hgvs` on the RefSeq caches of release 116. `--use_given_ref`, `--shift_3prime 1`, `--minimal`, other output formats and a merged cache were not run.
- The check of section 2.1, which holds each engine's coding terms to its own reference peptide, was made on eleven rows. It was not made on the four rows whose span holds the stop codon. Both engines write the same coding terms on those four.
- On the worked row that check predicts the kind of each engine's protein notation, a deletion with or without an extension. It does not rebuild the places the notation names.

## 1. ClinVar record 2530348: 693 bases deleted across intron 1 of `NM_000555.3`

**In short.** ClinVar record 2530348 on GRCh38 deletes the 693 bases `X:111410073-111410765`: the last 12 bases of exon 1 of `NM_000555.3`, the 333 bases of intron 1 and the first 348 bases of exon 2. On the model that is 360 coding bases, 120 codons in frame, and both engines place the deletion at residues 72 to 191. vep-rs writes `inframe_deletion` and a plain deletion in the protein notation. Ensembl VEP 116.3 translates intron bases inside its reference half and finds two stops there. It writes `stop_lost,protein_altering_variant` and adds a stop-lost extension to the notation. The model's stop codon is codon 442, and the deletion touches no base of it.

### 1.1 What a reader needs to know first

The RefSeq caches are BAM-edited caches and hold edited models, and the Ensembl caches hold none. The edit need not be near the record. On this example's model it is a run of bases appended after the last exon base.

`Codons` and `Amino_acids` each print a reference half before the slash and an alternate half after it. The reference half is read over the row's codon window, the codons of the residues `Protein_position` names. On this example the window is the 363 bases CDS 208 to 570, the codons of residues 70 to 190, around the coding span CDS 210 to 569.

Ensembl VEP's predicates for the coding terms compare the reference peptide or the reference codons with the alternate's. The reference peptide is the translation of the reference codons, and the protein notation starts from it as well.

### 1.2 The record

ClinVar record 2530348 is the VCF line at `X:111410072` on GRCh38 with a reference of 694 bases and an alternate of one. Both engines read it as a deletion of the 693 bases `X:111410073-111410765`. Section 1.6 writes the same deletion as a region line.

`NM_000555.3` lies on the reverse strand and has seven exons. Its coding region is cDNA 173 to 1498, which is 442 codons with the stop codon at codon 442, and its cached peptide has 441 residues. Its one `_rna_edit` appends 27 bases after cDNA 9406, the last exon base, far outside the coding region. In the model's direction the span has three parts:

| Part of the span | Genome | Model |
| --- | --- | --- |
| the last 12 bases of exon 1 | `X:111410754-111410765` | cDNA 382-393, CDS 210-221 |
| the whole of intron 1, 333 bases | `X:111410421-111410753` | no transcript position |
| the first 348 bases of exon 2 | `X:111410073-111410420` | cDNA 394-741, CDS 222-569 |

Both engines run `--refseq --hgvs` on the RefSeq cache of release 116. GRCh37 holds the same record on the same model at `X:110653301-110653993`, with intron 1 again 333 bases. Each engine writes there the coding terms and the protein notation it writes on GRCh38.

### 1.3 What each engine writes

The table gives the row of `NM_000555.3` cell by cell. The halves of `Codons` and `Amino_acids` are given by their lengths and places.

| Cell | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `cDNA_position`, `CDS_position`, `Protein_position` | 382-741, 210-569, 70-190 | the same |
| `USED_REF` | 693 bases | 693 bases |
| `Codons`, reference half | 363 bases: the model's CDS 208 and 209, then the first 361 bases of the allele | 363 bases: the model's CDS 208 to 570 |
| `Amino_acids`, reference half | 121 residues, with stops at places 6 and 53 | 121 residues, the model's 70 to 190, with no stop |
| `Codons` and `Amino_acids`, alternate halves | 3 bases (CDS 208, 209 and 570) and one residue, with no stop | the same |
| `Consequence` | `splice_acceptor_variant,splice_donor_variant,stop_lost,protein_altering_variant,splice_donor_5th_base_variant,intron_variant` | `splice_acceptor_variant,splice_donor_variant,inframe_deletion,splice_donor_5th_base_variant,intron_variant` |
| `IMPACT` | `HIGH` | `HIGH` |
| `HGVSp` | a deletion of residues 72 to 191 that names the model's residues at both places, followed by a stop-lost extension | the same deletion with the same two residues named, and no extension |
| `HGVSc` | equal in both engines | equal in both engines |

The first 361 bases of the allele are the 12 bases CDS 210-221, the 333 bases of intron 1 and the 16 bases CDS 222-237. Of the 121 places of Ensembl VEP's reference amino acids, 112 differ from the model's. The two `Consequence` cells differ in their coding terms alone: `stop_lost,protein_altering_variant` against `inframe_deletion`.

`USED_REF` is the reference allele an engine used, and it is 693 bases in both. What differs is where the reference codons are read from. vep-rs reads them from the model's coding sequence, and Ensembl VEP from that sequence with the allele written over it.

`Protein_position` reads 70-190, and the deletion in the protein notation runs from residue 72 to residue 191. Both engines write both, so neither is part of the difference, and this page takes them as written.

### 1.4 Why Ensembl VEP writes that

The line numbers are those of the source list at the top of this page. The span maps to three cDNA segments: exon, intron gap, exon. That is more than two, so the 693 genome bases stand as the reference allele (line 459). The model has an `_rna_edit` attribute, so those bases are written over the 360 coding positions 210 to 569 (lines 872 to 874). The 363-base window from CDS 208 then holds the model's two bases CDS 208 and 209, the 12 bases CDS 210-221, the 333 bases of the intron and the 16 bases CDS 222-237 (line 877).

The translation of that window has stops at places 6 and 53, and the alternate peptide has none, so `stop_lost` accepts (line 1260). The alternate codons are neither the start nor the end of that reference window, and trimming leaves a part of them. So `inframe_deletion` accepts by neither line 1178 nor line 1185, and `protein_altering_variant` takes its place. The notation is a deletion on a row that passes `stop_lost`, and it is written with the extension (lines 1916 and 1926).

### 1.5 Why vep-rs's row is the correct one

The deletion removes 693 genome bases, and 333 of them are the whole of intron 1. 693 less 333 leaves 360 coding bases, CDS 210 to 569, the length of 120 codons. On the model that is a deletion of 120 residues in frame, and both engines place it at residues 72 to 191.

The reference of the row is the model. Its residues 70 to 190 are the translation of its CDS 208 to 570, and its cached peptide of 441 residues holds no stop. The 333 bases that Ensembl VEP translates inside its reference half are intron 1, and no transcript position of the model holds them. The stops at places 6 and 53 of that half lie in those bases.

The model's stop codon is codon 442, cDNA 1496-1498, `X:111301687-111301689`, in exon 7. The deletion ends at cDNA 741 and touches no base of it. The definition of `stop_lost` asks for a change of at least one base of the terminator codon, and this row does not meet it.

Fed the model's own window, the same lines of Ensembl VEP's source give vep-rs's row. The reference peptide holds no stop, so `stop_lost` does not accept. The alternate codons are the window's first two bases and its last. Trimmed of those three bases, the alternate is empty and 360 reference bases remain, a multiple of three. That is what line 1185 asks for, so `inframe_deletion` accepts and the notation is a plain deletion. The two engines' notations name the same residues at the same places and differ only in the text after them, the extension.

### 1.6 Try it yourself

The reproduction `scripts/divergences/reproductions/refseq_edited_model_intron_span_coding_reference.toml` runs both engines on this deletion. It writes the deletion as a region line, which holds no base, and both engines read the 693 reference bases from the FASTA file:

```
X:111410073-111410765:1/-
```

Save that line as `record.txt` and run vep-rs on it:

```
vep -i record.txt -o out.txt --offline --format region \
  --json_cache <json_cache> --species homo_sapiens --assembly GRCh38 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --refseq --hgvs
```

`<json_cache>` is the conversion of the RefSeq GRCh38 cache of release 116. The reproduction copies the one region file of that cache that holds the model and converts it with `tools/vep-storable-to-json`. The reproduction lists `--refseq`, `--hgvs` and the FASTA file as its flags and gives the format, the cache, the species and the assembly as fields of their own. The runner that reads it passes `--offline`. Ensembl VEP 116.3 takes the same flags with the cache options of [Reproducing these examples](reproducing.md) and its `homo_sapiens_refseq` cache.

On the row of `NM_000555.3` the reproduction compares `Consequence`, the three position cells, `Amino_acids`, `Codons` and `HGVSp`, and it expects from each engine the values of section 1.3. The notation it expects from vep-rs is `NP_000546.2:p.His72_Ser191del`. From Ensembl VEP 116.3 it expects the same notation followed by `extTer?`. It holds each engine's `Amino_acids` and `Codons` cells whole.

### 1.7 How the comparator treats this row

The corpus `GRCh38-refseq-intron-span` holds this record twice, as a VCF line and as the region line of section 1.6, with Ensembl VEP's output of each run. Its manifest names the class and documents the row of `NM_000555.3` under it, with the coding terms vep-rs writes. The golden harness holds vep-rs to those terms.

The format comparator's built-in class of the same name reads the shape from the row alone. The row is a `Transcript` row that carries `BAM_EDIT=OK` in both engines, the same `CDS_position` range and the same `USED_REF`, which has more bases than the range. Both engines write `splice_acceptor_variant`, `splice_donor_variant` and `intron_variant` on it. On such a row the class sets four cells aside and compares every other field:

- `Consequence`, where Ensembl VEP alone writes `stop_lost` and `protein_altering_variant` and vep-rs alone writes `inframe_deletion`.
- `Codons`, where the alternate halves are equal, the reference halves have one length, and Ensembl VEP's upper-case bases spell, on the model's strand, the leading bases of `USED_REF` or of a `GIVEN_REF` of its length.
- `Amino_acids`, where the alternate halves are equal and each engine's reference half is the translation of its own reference codons.
- `HGVSp`, where Ensembl VEP's value is vep-rs's followed by `extTer?`.

A protein notation that differs in another way stays a difference, as on the three rows of section 2.2 beyond record 2530348. The comparison behind the concordance figures has no rule for the shape, and its datasets hold no row of it. The reproduction of section 1.6 holds both engines' rows and sets nothing aside.

## 2. The fifteen rows of genome-wide ClinVar

**In short.** Genome-wide ClinVar under `--refseq` has fifteen rows of the shape, 3 on GRCh38 and 12 on GRCh37. On every one the reference halves of `Codons` and `Amino_acids` differ between the engines, and the position cells and `HGVSc` are equal. On ten rows the coding terms, the impact and the protein notation are equal as well. On three rows both the coding terms and the protein notation differ: the worked record on each assembly, and record 2081311. Ensembl VEP reads the reference halves from the genome's allele, intron included, and writes `stop_lost,protein_altering_variant` on those three rows, where vep-rs reads the model's own codons and writes `inframe_deletion`. On the two rows of record 3780823 the protein notation alone differs.

### 2.1 The fifteen rows

Every row carries `splice_acceptor_variant,splice_donor_variant` before its coding terms and `splice_donor_5th_base_variant,intron_variant` after them in both engines, and `IMPACT=HIGH` in both. The table gives the coding terms between them. Its last column gives the protein notation by its kind and by the residue numbers it names. The kinds are a deletion, a deletion written with a stop-lost extension, a deletion-insertion and a frameshift form. The second column gives the span, its number of bases and the number of whole introns inside it.

| Assembly, record, model | Span; bases; introns inside | Coding terms | `HGVSp`: kind and places |
| --- | --- | --- | --- |
| GRCh38, 973541, `NM_022874.2` | `5:70946066-70951991`; 5926; 1 | `stop_lost,inframe_deletion` in both | a deletion written with a stop-lost extension, 210-263; equal |
| GRCh38, 1452802, `NM_001354801.1` | `7:44145167-44145585`; 419; 1 | `frameshift_variant` in both | a frameshift form at 52; equal |
| GRCh38, 2530348, `NM_000555.3` | `X:111410073-111410765`; 693; 1 | Ensembl VEP `stop_lost,protein_altering_variant`; vep-rs `inframe_deletion` | Ensembl VEP a deletion written with a stop-lost extension, 72-191; vep-rs a deletion, 72-191 |
| GRCh37, 3358865, `NM_133437.4` | `2:179421796-179425394`; 3599; 4 | `frameshift_variant` in both | a frameshift form at 19616; equal |
| GRCh37, 3358865, `NM_001267550.2` | `2:179421796-179425394`; 3599; 4 | `frameshift_variant` in both | a frameshift form at 28489; equal |
| GRCh37, 2081311, `NM_001267550.2` | `2:179517030-179517223`; 194; 1 | Ensembl VEP `stop_lost,protein_altering_variant`; vep-rs `inframe_deletion` | Ensembl VEP a deletion-insertion, 13030-13058; vep-rs a deletion, 13030-13057 |
| GRCh37, 973541, `NM_022874.2` | `5:70241893-70247818`; 5926; 1 | `stop_lost,inframe_deletion` in both | a deletion written with a stop-lost extension, 210-263; equal |
| GRCh37, 1452802, `NM_001354801.1` | `7:44184766-44185184`; 419; 1 | `frameshift_variant` in both | a frameshift form at 52; equal |
| GRCh37, 3780519, `NM_001351834.2` | `11:108198477-108199749`; 1273; 1 | `frameshift_variant` in both | a frameshift form at 2361; equal |
| GRCh37, 3780519, `NM_000051.4` | `11:108198477-108199749`; 1273; 1 | `frameshift_variant` in both | a frameshift form at 2361; equal |
| GRCh37, 3780823, `NM_001351834.2` | `11:108235809-108236235`; 427; 1 | `stop_lost,inframe_deletion` in both | a deletion written with a stop-lost extension, 2951-3057, in both; the last named residue differs |
| GRCh37, 3780823, `NM_000051.4` | `11:108235809-108236235`; 427; 1 | `stop_lost,inframe_deletion` in both | a deletion written with a stop-lost extension, 2951-3057, in both; the last named residue differs |
| GRCh37, 3780844, `NM_145862.2` | `22:29105994-29108005`; 2012; 1 | `frameshift_variant` in both | a frameshift form at 229; equal |
| GRCh37, 2850981, `NM_001123383.1` | `X:39922239-39923060`; 822, with 3 bases put in their place; 1 | `frameshift_variant` in both | a frameshift form at 1183; equal |
| GRCh37, 2530348, `NM_000555.3` | `X:110653301-110653993`; 693; 1 | Ensembl VEP `stop_lost,protein_altering_variant`; vep-rs `inframe_deletion` | Ensembl VEP a deletion written with a stop-lost extension, 72-191; vep-rs a deletion, 72-191 |

On GRCh38 the coding terms, the impact and the protein notation are equal on 2 of the 3 rows, and on GRCh37 on 8 of the 12. Together the three cells are equal on 10 of the 15 rows. On the other 5 the term differs on 3, the protein notation on all 5 and the impact on none.

Counting every transcript they fall on, the three records of GRCh38 have 35 rows, the rows of the shape among them, and the nine records of GRCh37 have 55. The three cells are equal on 34 of the 35 and on 51 of the 55.

On eleven of the fifteen rows the cause was checked against the row's own cells. On each of them Ensembl VEP's reference codons are the model's coding sequence with the reference allele written over the row's CDS span, cut to the window. vep-rs's reference codons are the model's window. Ports of Ensembl VEP's predicates, fed each engine's reference peptide, return the coding terms that engine wrote. The four rows left out are those of section 2.3.

### 2.2 The rows that differ

Record 2530348 on `NM_000555.3` is the worked row of section 1, where the coding terms and the protein notation differ. Its GRCh37 row differs in the same two cells, and each engine writes there the values it writes on GRCh38.

Record 2081311 on `NM_001267550.2`, on GRCh37, deletes 194 bases with one intron inside, and its CDS span is 84 bases. Ensembl VEP's codon window holds 47 intron bases. Ensembl VEP writes `stop_lost,protein_altering_variant` and a deletion-insertion notation over residues 13030 to 13058. Line 1916 names two notation types, `del` and `>`. This notation is a deletion-insertion, and it carries no extension. vep-rs writes `inframe_deletion` and a deletion of residues 13030 to 13057. The model's stop codon is codon 35992, and the span does not hold it.

Record 3780823, on GRCh37, has two rows, on `NM_000051.4` and on `NM_001351834.2`: 427 bases with one intron inside, and a CDS span of 8851 to 9171. The span holds the stop codon, codon 3057 of 3057. Both engines write `stop_lost,inframe_deletion` and a deletion written with a stop-lost extension, places 2951 to 3057. The two notations differ in the residue named at 3057 alone, where each engine names the last residue of its own reference half. vep-rs names the stop, the model's own residue at that place. Ensembl VEP's reference half is read from the allele's bases, and the residue it names is another.

On these three rows, that of record 2081311 and the two of record 3780823, Ensembl VEP's protein notation names a last residue that the model does not hold at that place. On the worked row it names the model's own residues at both places.

### 2.3 The four rows on the stop codon

On four rows the span holds the stop codon: record 973541 on each assembly and the two rows of record 3780823. Both engines write `stop_lost,inframe_deletion` on all four. The reference half of `Codons` differs between the engines there as on every row, and Ensembl VEP writes `inframe_deletion` all the same: the alternate codons are empty on these rows, and line 1178 accepts an empty alternate.

The check of section 2.1 leaves out the branches of the predicates for a span on the stop codon, so it was not made on these four rows.
