# `start_lost` and `start_retained_variant` on a transcript that shows no start codon

**The defect in one sentence.** Ensembl VEP 116.3 reads the first three coding positions of every transcript without the `cds_start_NF` flag (coding sequence start not found) as its start codon, whatever they hold. For a change there it writes `start_lost` or `start_retained_variant`, also on a transcript whose record shows no start codon. vep-rs writes the two terms only where the record shows a start codon, and elsewhere the row carries the terms the sequence supports.

**Why this is a defect.** Both terms state a fact about a start codon: that the change loses it, or that the change touches it and it remains. On these transcripts the record shows no start codon at the coding start: the first codon names another residue than Met and is no start codon of the genetic code, the recorded translation opens with that residue, and nothing in the record asserts a start. On others, gene segments of the RefSeq caches, the coding sequence begins inside a codon instead. Either way the record is that of a coding sequence cut off at its 5' end, without the flag that says so.

The tests behind the terms never ask whether the three positions held a start codon. They compare three bases of the edited sequence with the codon of Met, which this page calls the standard start codon, or they compare the residues at protein position 1. Behind a 5' UTR a change that keeps the first residue is then a lost start, and a change that makes the standard start codon there is a start kept.

**What vep-rs does instead.** vep-rs asks the cache record whether a start codon stands at the coding start. One stands where the record carries no `cds_start_NF`, the first codon is whole, and one of three holds:

1. The first codon is a start codon of the transcript's genetic code.
2. The translation the cache records opens with Met. It does where the cache records a start the annotation asserts on another codon: the first residue is then Met.
3. The record carries a consensus coding sequence (CCDS) id, which names a coding sequence annotated as complete.

On every other transcript no start codon is shown, and no change takes `start_lost` or `start_retained_variant` from the first three coding positions. The row keeps the terms the other tests give: `synonymous_variant`, `missense_variant`, `stop_gained`, `frameshift_variant`, the in-frame terms, and `coding_sequence_variant` for a change written in bases that runs from the 5' UTR into the coding sequence. A structural allele keeps its own terms. They are the terms both engines write for the same change on a transcript with the flag.

The rule reads the record and not the biology behind it. That sets three limits, and on the records of the second and the third vep-rs's row is not offered as the correct one:

- A record with a CCDS id is read as showing a start codon. On three records of the release 116 Ensembl caches the reference assembly differs from the annotated model at the first codon, so the cached bases are no start codon: ENST00000295641 on GRCh37, and ENST00000641572 and ENST00000641848 on GRCh38. vep-rs keeps the start terms there.
- A start that an annotation places on another codon is shown only where the cache records it, by a translation that opens with Met or by a CCDS id. In the GRCh37 Ensembl cache no translation opens with Met on a codon that is no start codon. Ten transcripts of one gene that GENCODE 19 tags as starting on another codon carry no CCDS id either, so they are read as showing no start codon, and a change in their first codon gets no start term.
- Both engines take the coding start from the cache. Where a run of bases that a RefSeq model has and the genome lacks begins at the coding start, the cache places the coding start earlier than the model does. NM_001288959.2 of the GRCh37 RefSeq cache is such a record: its cached coding sequence opens four bases early, out of the model's frame, on a codon that is no start codon. The model's own coding sequence opens on the standard start codon. vep-rs's row there is right for the cached record only.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the module the Ensembl VEP 116.3 image holds), all in `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`:

- Lines 965-990, `_overlaps_start_codon`: line 975 returns 0 for a transcript that carries `cds_start_NF`. Lines 983-986 overlap the span's cDNA coordinates with the first three coding positions. No base, residue or exon phase is read.
- Lines 851-904, `start_lost`: once that overlap is found (line 862), line 870 takes the answer of `_ins_del_start_altered` for a change of length that is not in frame, and line 871 that of `_inv_start_altered`. Lines 878-882 compare the reference and alternate residues at protein position 1, where the allele has both (lines 874-875).
- Lines 906-945, `_inv_start_altered`: line 926 returns 0 on a transcript with no 5' UTR. Lines 938-941 make the edit in the 5' UTR and coding sequence and compare the three bases at the old coding start with the standard start codon.
- Lines 947-963, `start_retained_variant`: after the same overlap test (line 957), the term is the negation of `_snp_start_altered` (line 959) or of `_ins_del_start_altered` (line 961).
- Lines 992-1026 and 1028-1075, those two routines: lines 1021-1022 and 1061-1065 compare three bases of the edited sequence with the standard start codon. Line 1069 answers that the start is altered where the edited sequence is shorter than the coding sequence, and line 1071 where it does not end with the unedited coding sequence. Line 1037 returns 0, not altered, for every structural allele.
- Lines 886-896, the structural part of `start_lost`: after the test of line 862, the span's genomic coordinates are overlapped with the three genomic positions of the coding start.
- Lines 1082, 1090 and 1109: `synonymous_variant` is withheld by `start_retained_variant` and not by `start_lost`, and `missense_variant` and `inframe_insertion` are withheld by `start_lost`.

**What correcting it would change.** A test beside line 975 that asked whether the record shows a start codon would take the two terms off these rows. Each row would then carry the terms Ensembl VEP writes for the same change on a transcript with the flag, which are the terms vep-rs writes. On a transcript that shows a start codon nothing would change.

**The class it produces.** `start_codon_overlap_without_start_codon`: a `Transcript` row on which Ensembl VEP writes `start_lost`, `start_retained_variant` or both and vep-rs writes neither, on a transcript whose cache record carries no `cds_start_NF` and shows no start codon. Section 1.7 says how a comparison treats it.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data carries the shape, and the divergence needs no flag or plugin. That page is the index of the documented differences, and it counts each shape over 15,646,481 chromosome 21 records. The shape the record alone shows is any change: 15,646,449 of the records write an alternate allele in bases or a structural variant, 999,998 per million. The rate is an upper bound: the transcript decides whether the divergence then occurs.

The Ensembl caches of release 116 hold 168 transcript records of the kind on GRCh38 and 479 on GRCh37. Each lacks the flag, opens on a whole codon that is no start codon, has a translation that does not open with Met and carries no CCDS id. The RefSeq caches add 143 gene segments on GRCh38 and 138 on GRCh37 whose first codon is padded, the kind of section 3, reached with `--refseq --all_refseq`.

The counts that follow are rows on which Ensembl VEP 116.3 writes a start term and vep-rs none. The genome-wide ClinVar dataset has 14 such rows on GRCh38 and 14 on GRCh37, from 12 input records on each assembly. Under `--refseq` it has one such row, on GRCh37, on the record of the third limit. The whole-genome population datasets were counted over every record that comes within 1,000 bases of such a coding start on a primary chromosome. They have 62 rows in gnomAD v2.1.1, 17 in 1000 Genomes Phase 3, 15 in gnomAD v4.1 and 1 in the 1000 Genomes high-coverage dataset.

Set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A change that keeps the first residue behind a 5' UTR, on ENST00000420031

**In short.** `1 16901031 residue_kept_1_16901031 A T` on GRCh37 is a made-up substitution of one base that, on the reverse-strand transcript ENST00000420031, changes the third base of the first codon and turns one codon of Ile into another. Ensembl VEP 116.3 writes `start_lost,synonymous_variant`, and vep-rs writes `synonymous_variant`. The record shows no start codon: the first codon names Ile, the recorded translation opens with Ile, and the record carries neither the flag nor a CCDS id. Ensembl VEP compares the three edited bases with the standard start codon, and the bases were not that codon before the change either.

### 1.1 What a reader needs to know first

A complete coding sequence opens with a start codon, the codon at which translation begins. In the standard genetic code that is the codon of Met, the standard start codon, and the code lists two codons of Leu as start codons beside it. For a change at the start codon Ensembl VEP writes one of two terms. The Sequence Ontology defines `start_lost` as "A codon variant that changes at least one base of the canonical start codon". It defines `start_retained_variant` as "A sequence variant where at least one base in the start codon is changed, but the start remains".

Not every transcript in a cache is complete at its 5' end. Where the annotation marks a coding sequence as cut off at its start, the record carries the attribute `cds_start_NF`, and Ensembl VEP writes neither term on that transcript. Other records are cut off and carry no flag. Their first codon is a codon from inside a reading frame, and the Sequence Ontology says of `start_lost`: "This term should not be applied to incomplete transcripts."

A cache record says three things about the coding start. Its coding sequence gives the first codon. Its translation gives the first residue, which is Met where the first codon is a start codon or the cache records an annotated start on another codon. And its consensus coding sequence (CCDS) id, where it has one, names a coding sequence that the CCDS project annotates as complete.

A cDNA coordinate is the position of a base in the spliced transcript, and a CDS coordinate is its position in the coding sequence. The 5' UTR is the part of the transcript before the coding sequence.

### 1.2 The record

A made-up substitution, run against the Ensembl GRCh37 cache of release 116:

```
1	16901031	residue_kept_1_16901031	A	T	.	.	.
```

ENST00000420031 is a reverse-strand transcript of six exons of the gene NBPF1. It runs from high to low coordinates and carries the complement of the genome's bases. Its record carries neither `cds_start_NF` nor a CCDS id. Three of its parts matter here:

| Part | Genomic span | cDNA |
| --- | --- | --- |
| Exon 1 | 1:16900982-16901051 | 1-70 |
| 5' UTR | 1:16901034-16901051 | 1-18 |
| First codon | 1:16901031-16901033 | 19-21 |

The coding sequence runs from cDNA 19 to cDNA 120, 120 - 19 + 1 = 102 bases, which is 34 codons. The recorded translation has 34 residues, so every codon names a residue. The first codon is cDNA 19 to 21, the bases 1:16901033, 16901032 and 16901031 in the transcript's direction, and names Ile. The translation opens with Ile.

The record changes 1:16901031, the third base of the first codon, from A to T on the genome and from T to A on the transcript. The codon before the change and the codon after it both name Ile. The committed GRCh37 cache holds the transcript, and both engines write the rows of 1.3 with no flag beyond their cache options.

### 1.3 What each engine writes

The row on ENST00000420031, the differing cells in bold:

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| Location | `1:16901031` | `1:16901031` |
| Allele | `T` | `T` |
| Gene | `ENSG00000219481` | `ENSG00000219481` |
| Feature | `ENST00000420031` | `ENST00000420031` |
| Consequence | **`start_lost,synonymous_variant`** | **`synonymous_variant`** |
| cDNA_position | `21` | `21` |
| CDS_position | `3` | `3` |
| Protein_position | `1` | `1` |
| Amino_acids | `I` | `I` |
| IMPACT | **`HIGH`** | **`LOW`** |

The rows differ in the `Consequence` and in the `IMPACT`, which follows the most severe term. The `Codons` column is equal in both rows and shows the codon's third base changed. `Amino_acids` holds one residue, Ile, since the change keeps it.

The record writes five rows in each engine, and the other four are equal in every column. On three of them, the rows of other transcripts of the gene, the same base is the third base of a codon inside the coding sequence: codon 148 of ENST00000287968, codon 241 of ENST00000432949 and codon 783 of ENST00000430580. Both engines write `synonymous_variant` there. The fourth is a `3_prime_UTR_variant` row of another transcript of the gene.

### 1.4 Why Ensembl VEP writes that

`start_lost` first asks `_overlaps_start_codon` (line 862). The transcript carries no `cds_start_NF`, so line 975 passes. The span is cDNA 21, which lies in cDNA 19 to 21, the first three coding positions, and lines 983-986 answer 1. No base has been read.

The allele keeps the length, so `_ins_del_start_altered` answers 0 at line 1042 and line 870 does not fire. Line 871 asks `_inv_start_altered`, and the transcript has a 5' UTR, so line 926 passes.

The routine joins the 5' UTR and the coding sequence, makes the edit at cDNA 21 (line 938) and reads the three bases that follow the 18 bases of the UTR (line 939). They are the edited first codon, a codon of Ile. Line 941 finds that they are not the standard start codon and answers 1, and line 871 sets `start_lost`.

The test of line 941 could not have come out otherwise. The change keeps the residue, the residue is Ile, and Met has one codon in the standard genetic code. The three bases were not that codon before the change either, and no line reads them as they stood.

`synonymous_variant` is written beside the term. The reference and alternate residues are equal, and line 1082 withholds the term for `start_retained_variant` but not for `start_lost`. `start_retained_variant` is not written: line 959 negates `_snp_start_altered`, which answers 1 unless the three edited bases are the standard start codon (lines 1021-1022).

The 5' UTR decides the term. On the same coding sequence with no 5' UTR, line 926 would return 0, lines 878-882 would find the residues equal, and the row would be `synonymous_variant` alone.

### 1.5 Why vep-rs's row is the correct one

The record's own coordinates say what the three bases are. cDNA 19 to 21 are the first of 34 codons, and the codon names Ile. The translation the cache records for the transcript opens with Ile. The standard genetic code does not start translation at a codon of Ile, and the record carries no CCDS id and no Met at residue 1 to assert a start there. No start codon is shown.

The substitution changes the codon's third base and leaves a codon of Ile. The first residue is the same before and after, so the change is a `synonymous_variant`, and that is all the sequence shows. `start_lost` would need a start codon at cDNA 19 to 21 before the change and none after it. None stood there before it.

Ensembl VEP's row contradicts itself. The Sequence Ontology places `start_lost` under `nonsynonymous_variant`, and the row carries `synonymous_variant` beside it. Its `IMPACT` is `HIGH` for a change that leaves the protein as it was.

The other rows of the record point the same way. Three transcripts of the gene with longer coding sequences read the same three bases as a codon of Ile inside the coding sequence, in the same frame, and both engines call the record `synonymous_variant` there. The first codon of ENST00000420031 is a codon of that reading frame.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --species homo_sapiens --assembly GRCh37
grep ENST00000420031 out.txt
```

The committed GRCh37 cache holds the transcript. vep-rs writes `synonymous_variant` on the row. Ensembl VEP 116.3 takes the record with the cache options of [Reproducing these examples](reproducing.md) and its GRCh37 cache, and writes `start_lost,synonymous_variant`.

The entry's reproduction, `scripts/divergences/reproductions/start_codon_overlap_without_start_codon.toml`, runs this record and twelve more on both engines and holds each row to both recorded values. Sections 2 and 4 give the other twelve.

### 1.7 How the comparison treats this class

A comparison sets a row of this class aside only when it is given the cache vep-rs annotated with, since the class reads the transcript's record there. It sets a row aside where the record carries no `cds_start_NF` and shows no start codon, Ensembl VEP's row holds `start_lost` or `start_retained_variant`, and vep-rs's row holds neither. A comparison that is given no cache sets nothing aside under this class, and the row stays a difference charged to vep-rs.

Under `--hgvs` the class reads the row's `HGVSp` with its terms, and sets aside these pairs:

| Change | Ensembl VEP's `HGVSp` | vep-rs's `HGVSp` |
| --- | --- | --- |
| A substitution of one residue at protein position 1 | The lost start: the residue, the number 1 and `?` | The two residues |
| A frameshift, under the standard genetic code | The lost start at the residue vep-rs's notation opens with: that residue, its number and `?` | The frameshift at its first changed residue, with its count |
| The same frameshift on a row without `HGVS_OFFSET`, where that residue is not the last of the row's protein position | A range: that residue and its number, `_`, `?` and the last number of the protein position | The frameshift at its first changed residue, with its count |

## 2. A change of the first residue of a transcript with no 5' UTR, on ENST00000600312

**In short.** `21 37402368 residue_changed_at_first_base_21_37402368 C G` on GRCh37 is a made-up substitution of the first base of the forward-strand transcript ENST00000600312. The transcript has no 5' UTR, its first codon names Leu and is no start codon, and the change turns Leu into Val. Ensembl VEP 116.3 writes `start_lost`, and vep-rs writes `missense_variant`. With no 5' UTR Ensembl VEP compares no base: the term comes from the residues at protein position 1, which differ.

### 2.1 What a reader needs to know first

Section 1.1 gives the two terms and what a cache record says about the coding start. Here the coding sequence begins on the transcript's first base, so Ensembl VEP does not make its test of the three bases at the old coding start.

### 2.2 The record

A made-up substitution, run against the same cache:

```
21	37402368	residue_changed_at_first_base_21_37402368	C	G	.	.	.
```

ENST00000600312 is a forward-strand transcript of two exons, 21:37402368-37402428 and 21:37417766-37417818, of 61 and 53 bases. All 114 bases are coding: the coding sequence runs from cDNA 1 to cDNA 114, and the transcript has no UTR at either end. The recorded translation has 38 residues, one for each of the 114 / 3 = 38 codons, and opens with Leu. The record carries neither `cds_start_NF` nor a CCDS id.

The first codon is 21:37402368-37402370 and names Leu. It is neither of the two codons of Leu that the standard genetic code lists as start codons. The record changes the codon's first base, and the codon then names Val.

### 2.3 What each engine writes

The row on ENST00000600312, the differing cells in bold:

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| Location | `21:37402368` | `21:37402368` |
| Allele | `G` | `G` |
| Feature | `ENST00000600312` | `ENST00000600312` |
| Consequence | **`start_lost`** | **`missense_variant`** |
| cDNA_position | `1` | `1` |
| CDS_position | `1` | `1` |
| Protein_position | `1` | `1` |
| Amino_acids | `L/V` | `L/V` |
| IMPACT | **`HIGH`** | **`MODERATE`** |

The record writes six rows in each engine. The other five are `downstream_gene_variant` rows of another gene's transcripts and are equal in every column.

### 2.4 Why Ensembl VEP writes that

The transcript carries no flag, and cDNA 1 lies in the first three coding positions, so `_overlaps_start_codon` answers 1. The allele keeps the length, and line 870 does not fire. `_inv_start_altered` returns 0 at line 926, since the transcript has no 5' UTR, and no base is compared.

Lines 878-882 then compare the residues. The span begins at protein position 1, the reference residue is Leu and the alternate is Val, and the alternate neither begins nor ends with the reference, so `start_lost` is set. Line 1090 withholds `missense_variant` because of it. On this path a first residue that changes into another is a lost start on any transcript without the flag, whatever the residue was.

### 2.5 Why vep-rs's row is the correct one

The transcript's own coordinates show a fragment. Its first base is its first coding base, and its last base is its last coding base. Each of the 38 codons names a residue. The model has neither a start codon nor a stop codon, and its translation opens with Leu on a codon at which the genetic code does not start.

The record changes residue 1 from Leu to Val. That is a `missense_variant`, the term both engines write for a changed residue at any other codon. `start_lost` would need a start codon at cDNA 1 to 3 before the change, and the record shows none.

### 2.6 Try it yourself

Run the command of 1.6 with this record as `record.vcf`, and read the row with `grep ENST00000600312 out.txt`. vep-rs writes `missense_variant`, and Ensembl VEP 116.3 writes `start_lost`.

## 3. A first codon that the start exon's phase pads, on the gene segment IGKC

**In short.** Some gene segments of the RefSeq caches begin inside a codon and carry no `cds_start_NF`. Their start exon has phase 1 or 2, and the cached coding sequence opens with as many unknown bases, which keep the frame. The transcript's first three coding bases are then no codon of the frame. For a deletion of one base there Ensembl VEP 116.3 writes `frameshift_variant,start_lost`, and vep-rs writes `frameshift_variant`.

### 3.1 What a reader needs to know first

An exon's phase says how many bases of its first codon lie before the exon. A coding sequence that begins in an exon of phase 1 begins on the second base of a codon. The cache then records it with one unknown base in front, so that every later codon reads in frame. In the Ensembl caches every such record carries `cds_start_NF`, and the RefSeq gene segments of this section carry no flag.

### 3.2 The record

A made-up deletion of one base, run against the GRCh38 RefSeq cache of release 116:

```
2	88857681	one_base_deleted_at_padded_codon_2_88857682	TT	T	.	.	.
```

IGKC is a gene segment of one exon on the reverse strand, 2:88857361-88857683, 323 bases, all of them coding. Its exon has phase 1. The cached coding sequence has 324 positions, one unknown base and the 323 bases, which is 108 codons, and the recorded translation opens with an unknown residue. The record carries neither `cds_start_NF` nor a CCDS id.

The record deletes 2:88857682, the transcript's second base, and its line is written at the base before it, as VCF writes a deletion. The deleted base is cDNA 2 and position 3 of the cached coding sequence: the last base of the codon whose first base lies outside the transcript.

### 3.3 What each engine writes

Both rows on IGKC carry cDNA position 2, CDS position 3, protein position 1 and `IMPACT=HIGH`. Ensembl VEP 116.3 writes `frameshift_variant,start_lost` as the `Consequence`, and vep-rs writes `frameshift_variant`. A substitution of the transcript's first base, 2:88857683, is `coding_sequence_variant` in both engines.

### 3.4 Why Ensembl VEP writes that

The transcript carries no flag, and cDNA 2 lies in cDNA 1 to 3, which line 985 takes for the first three coding positions. `_ins_del_start_altered` then makes the edit in the padded coding sequence at the position the cDNA coordinate names, which the unknown base puts one position before the deleted base (line 1058). The edited sequence is shorter than the coding sequence, and line 1069 answers that the start is altered. The deletion is not in frame, so line 870 sets `start_lost`.

### 3.5 Why vep-rs's row is the correct one

cDNA 1 to 3 are the second and third base of one codon and the first base of the next. They are no codon, so they are no start codon, and the first codon's first base is not part of the transcript at all. The deletion removes one base of a coding sequence and shifts the frame: `frameshift_variant`.

### 3.6 Try it yourself

No committed cache holds such a gene segment, so the entry's reproduction does not run this record. Ensembl VEP 116.3 takes it with its GRCh38 RefSeq cache, the reference file and `--refseq --all_refseq`. vep-rs takes it with the same flags on the conversion of that cache, which [Cache conversion](../cache-conversion.md) describes. Without `--all_refseq` neither engine writes a row on IGKC.

## 4. The other records of the reproduction

The reproduction runs eleven more made-up records on the cache of sections 1 and 2. Eight lie on ENST00000420031 and two on ENST00000600312. One lies on ENST00000601068, a reverse-strand transcript with no 5' UTR whose first codon, 21:47190003-47190005, names Cys.

| Record | The change | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `residue_changed_1_16901033` | First coding base, Ile to Val | `start_lost` | `missense_variant` |
| `standard_start_codon_made_1_16901031` | Third coding base, Ile to Met | `start_lost,start_retained_variant` | `missense_variant` |
| `one_base_deleted_1_16901033` | First coding base deleted | `frameshift_variant,start_lost` | `frameshift_variant` |
| `one_base_deleted_start_codon_made_1_16901032` | Second coding base deleted | `frameshift_variant,start_retained_variant` | `frameshift_variant` |
| `first_codon_deleted_1_16901031` | First codon deleted | `start_lost,inframe_deletion` | `inframe_deletion` |
| `three_bases_inserted_1_16901031` | Three bases added inside the first codon | `start_lost` | `inframe_insertion` |
| `deleted_from_5_prime_utr_1_16901033` | Two UTR bases and the first coding base deleted | `start_lost,5_prime_UTR_variant` | `coding_sequence_variant,5_prime_UTR_variant` |
| `symbolic_deletion_1_16901025` | `<DEL>` of 1:16901025-16901040 | `frameshift_variant,start_lost,feature_truncation,start_retained_variant,5_prime_UTR_variant` | `frameshift_variant,feature_truncation,5_prime_UTR_variant` |
| `one_base_deleted_at_first_codon_21_37402369` | Second coding base deleted, no 5' UTR | `frameshift_variant,start_lost` | `frameshift_variant` |
| `symbolic_deletion_at_first_base_21_37402368` | `<DEL>` of 21:37402368-37402373 | `start_lost,feature_truncation,inframe_deletion,start_retained_variant` | `feature_truncation,inframe_deletion` |
| `stop_codon_made_at_first_codon_21_47190003` | Third base of the first codon, Cys to a stop codon | `stop_gained,start_lost` | `stop_gained` |

Four of Ensembl VEP's rows carry `start_retained_variant`. Where a substitution makes the standard start codon, line 1022 answers that the start is kept, while lines 878-882 call it lost because Ile became Met. Where a deletion leaves the standard start codon at the old coding start after an unchanged 5' UTR, line 1065 answers that the start is not altered, on a frameshift. On a symbolic deletion `_ins_del_start_altered` returns 0 at line 1037 before any base is read, and `start_retained_variant` is the negation of that answer. The `start_lost` of those two rows comes from the genomic overlap of lines 886-896.

The two in-frame rows differ by one line. `inframe_insertion` is withheld where `start_lost` is set (line 1109). `inframe_deletion` has no such test, so the deletion's row carries both terms.

On a transcript that shows a start codon, a row with both start terms belongs to [`start_cooccurrence_swap`](defect-b.md), where vep-rs keeps the one term the record supports. On these transcripts the record shows no start codon for either term to describe, and vep-rs writes neither.

A deletion written in bases of the whole first codon of a transcript with no 5' UTR leaves no alternate residue for lines 878-882 to compare (line 875). It is `inframe_deletion` in both engines and no row of this class.
