# A start codon removed whole on a transcript with no 5' UTR

**The defect in one sentence.** A deletion can remove whole codons from the first base of a coding sequence, the start codon first. On a transcript that holds no 5' UTR, Ensembl VEP writes `inframe_deletion` and no start term for it. It writes `start_lost` for the same removal of the same codon on a transcript that holds one. vep-rs writes `start_lost` on both.

**Why this is a defect.** `start_lost` names a change of at least one base of the start codon, and a start codon removed whole is changed in all three. Whether a transcript model holds bases before its start codon changes nothing in that. Many models hold none: shorter isoforms of a gene, every mitochondrial transcript, and most gene models of a genome annotated from predictions. On such a row Ensembl VEP's answer is `IMPACT=MODERATE`, where no other term gives more, and, under `--hgvs`, a protein that lacks its first residues alone. The row of a transcript that shares the codon reads `HIGH` and `p.Met1?`. A filter on `HIGH` or on `start_lost` keeps one row and drops the other.

**What vep-rs does instead.** Where a deletion takes whole codons from the first coding base of a transcript without a 5' UTR, vep-rs edits the coding sequence as the record does and reads the three bases at its start. It writes `start_lost` unless they code methionine under the transcript's genetic code.

Under every genetic code `ATG` is such a triplet, and it is the one Ensembl VEP's test asks for behind a 5' UTR. Under five mitochondrial codes (2, 3, 5, 13 and 21) `ATA` is one as well. vep-rs makes the test where the coding sequence opens on a start codon that codes methionine. That is the coding sequence as the transcript holds it, before the removal. `IMPACT` and, under `--hgvs`, the HGVSp follow the terms as on any row that carries `start_lost`. Section 3 lists what stays as it is.

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `ensembl` at release/116, the modules the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 851-904, `start_lost`, has three ways to the term, tried in turn: the length test (line 870), `_inv_start_altered` (line 871) and the comparison of the two peptides (lines 872-882). Line 870 runs the length test of lines 1028-1075 and takes its answer unless the row is an in-frame insertion or deletion. That test reads the start as altered here, at line 1069, but the row is an in-frame deletion, so the line passes on.
- The same file, lines 906-945, `_inv_start_altered`, the second way: it edits the 5' UTR and the coding sequence and holds the three bases at the old start position against `ATG` (lines 936-941). Line 926 returns 0 before the edit when the transcript has no 5' UTR. Line 936 takes the length of the UTR's sequence, which such a transcript lacks, and line 926 keeps the routine from reaching it.
- The same file, lines 872-882, the third way, compare the two peptides. Line 875 returns 0 when the alternate peptide is empty, as it is when whole codons are removed.
- The same file, lines 947-963, `start_retained_variant`: the negation of the length test, false here, since that test reads the start as altered.
- `Bio/EnsEMBL/Transcript.pm`, lines 1753-1755: a transcript whose coding sequence starts at its first base has no 5' UTR sequence.

**What correcting it would change.** `_inv_start_altered` would make its edit with an empty UTR. A start codon removed whole would then read the same on every transcript that shares the codon.

**The class it produces.** `start_codon_removed_without_five_prime_utr`: a `Transcript` row of a deletion whose cDNA position and CDS position both start at 1 and whose `Amino_acids` holds no alternate residue. Both engines write `inframe_deletion`, and vep-rs adds `start_lost`. `IMPACT` and HGVSp follow each side's terms. Section 1.7 says what each comparison reads.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a record that removes a multiple of three bases. The five chromosome 21 record sets the index calls the cuts hold 229,923 such records among 15,646,481, the index's 14,695 per million. The rate is an upper bound: the removal must also start on the first coding base of a transcript without a 5' UTR. Section 2 counts the records with a single alternate allele that do in five whole-genome datasets.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. The start codon shared by ENST00000303645 and ENST00000339944

**In short.** A record written for this page removes 21:32298972-32298974 on GRCh38, the start codon of seven coding transcripts of one gene. Six hold a 5' UTR before it, and both engines write `start_lost,inframe_deletion` with `IMPACT=HIGH` on them. `ENST00000339944` starts on the codon. Ensembl VEP 116.3 writes `inframe_deletion` with `IMPACT=MODERATE` there, and vep-rs writes the row of the other six.

### 1.1 What a reader needs to know first

The transcripts of a gene are models of its messenger RNAs. Several can share a start codon, the three bases at which translation begins, and differ in how many bases they hold before it. Those bases are the 5' UTR. A model can hold none, so that its first base is its first coding base.

A deletion of three coding bases keeps the reading frame and is an `inframe_deletion`. When the three bases are the start codon, the Sequence Ontology term `start_lost` applies as well: a codon variant that changes at least one base of the canonical start codon.

### 1.2 The records

```
21	32298971	start_codon_removed	CATG	C	.	.	.
21	32298971	two_codons_removed	CATGGCC	C	.	.	.
```

The reference reads `C` at 21:32298971, `ATG` at 32298972-32298974 and `GCC` after it. `ENST00000303645` begins 117 bases before the `ATG`, at 21:32298855, so the codon is its `c.1_3` at cDNA positions 118 to 120. `ENST00000339944` begins at 21:32298972, so the same three bases are its `c.1_3` at cDNA positions 1 to 3. Neither transcript carries `cds_start_NF`.

The first record removes the start codon. The second removes the start codon and the alanine codon after it. After either record no `ATG` stands at 21:32298972.

Both engines run `--hgvs` on the release 116 cache and the reference of the corpus `tests/golden/116/GRCh38-plugin-loftee-fixture`.

### 1.3 What each engine writes

Ensembl VEP 116.3:

| Record | Feature | Consequence | HGVSp |
| --- | --- | --- | --- |
| `start_codon_removed` | `ENST00000303645` | `start_lost,inframe_deletion` | `p.Met1?` |
| `start_codon_removed` | `ENST00000339944` | `inframe_deletion` | `p.Met1del` |
| `two_codons_removed` | `ENST00000303645` | `start_lost,inframe_deletion` | `p.MetAla1_?2` |
| `two_codons_removed` | `ENST00000339944` | `inframe_deletion` | `p.Met1_Ala2del` |

vep-rs:

| Record | Feature | Consequence | HGVSp |
| --- | --- | --- | --- |
| `start_codon_removed` | `ENST00000303645` | `start_lost,inframe_deletion` | `p.Met1?` |
| `start_codon_removed` | `ENST00000339944` | `start_lost,inframe_deletion` | `p.Met1?` |
| `two_codons_removed` | `ENST00000303645` | `start_lost,inframe_deletion` | `p.MetAla1_?2` |
| `two_codons_removed` | `ENST00000339944` | `start_lost,inframe_deletion` | `p.MetAla1_?2` |

`IMPACT` is `HIGH` beside `start_lost` and `MODERATE` beside `inframe_deletion` alone. The HGVSp is on `ENSP00000306697.5` and `ENSP00000343661.4`. The engines agree on every other cell: `CDS_position` `1-3` and `1-6`, `Protein_position` `1` and `1-2`, `Amino_acids` `M/-` and `MA/-`, and the HGVSc, `c.1_3del` and `c.2_7del`. The second notation is one base 3' of the record, where the HGVS 3' rule places the removal, with `HGVS_OFFSET=1` in both engines.

### 1.4 Why Ensembl VEP writes that

On both transcripts the span overlaps the first three coding bases, so `start_lost` goes on to its three ways. The length test reads the start as altered on both: the edited sequence does not end with the whole coding sequence (line 1069 without a 5' UTR, line 1071 with one). Line 870 sets that aside, because the row is an in-frame deletion.

On `ENST00000303645` `_inv_start_altered` then edits the UTR and the coding sequence. The three bases at the old start position are `GCC` after the first record and `AAC` after the second. Neither is `ATG`, and the routine gives `start_lost`.

On `ENST00000339944` the routine returns at line 926: the transcript has no 5' UTR. The peptide comparison is left. The reference peptide is `M`, or `MA`, and the alternate peptide is empty, so line 875 returns 0. No line has found the start codon kept, and none names it lost.

### 1.5 Why vep-rs's row is the correct one

The two transcripts read the same three bases as their start codon, and the record removes those bases from both. After it, each coding sequence opens `GCC AAC`, or `AAC` for the second record. No start codon stands where translation began. The 117 bases `ENST00000303645` holds before the codon are not changed by the record and change nothing in what the record removes.

`start_lost` is defined by the start codon's bases, and all three are removed. `start_retained_variant` asks that the start remains, and it does not. So the term belongs on both rows, and with it `IMPACT=HIGH`.

The HGVS recommendations, in their text on the translation initiation codon, describe a variant of that codon as `p.(Met1?)`, its consequence unknown, and advise against a description that treats the codon as an ordinary residue. `p.Met1del` is such a description. `p.MetAla1_?2` is the form Ensembl VEP's formatter gives a lost start over two residues, and both engines write it on `ENST00000303645`.

### 1.6 Try it yourself

With the records of 1.2 saved as `records.vcf`, the command as vep-rs takes it:

```
vep -i records.vcf -o out.txt --offline --hgvs \
  --json_cache tests/golden/116/GRCh38-plugin-loftee-fixture/json_cache \
  --fasta tests/golden/116/GRCh38-plugin-loftee-fixture/reference.fa.gz
grep -e ENST00000303645 -e ENST00000339944 out.txt
```

Ensembl VEP 116.3 takes the same records and `--hgvs` with the cache options of [Reproducing these examples](reproducing.md). The entry's reproduction, `scripts/divergences/reproductions/start_codon_removed_without_five_prime_utr.toml`, runs both engines on these records.

### 1.7 How the comparisons treat this class

The format comparator (`tools/vep-formats-compare`) takes a `Transcript` row whose two lists of terms are equal but for vep-rs's `start_lost`, both holding `inframe_deletion`. The row's cDNA position and CDS position both start at 1, and its `Amino_acids` holds no alternate residue. On such a row it sets aside the `Consequence`, an `IMPACT` that is the one each side's own terms give, and an HGVSp that is a `del` notation from residue 1 on Ensembl VEP's side and the notation of a lost start on vep-rs's. Any other cell that differs stays a difference.

The tuple comparison (`tools/vep-outputs-compare`) scores consequence terms. For a pair of tuples whose terms differ by vep-rs's `start_lost` alone, it reads the two output rows behind the pair. It sets the pair aside where both rows hold a cDNA position and a CDS position that start at 1 and no alternate residue. The structural-variant comparator applies the same rule. The manifests of the committed corpora that hold such a record list its rows under the class, and the corpus comparison sets them aside. A row off the shape stays a difference charged to vep-rs.

The comparisons read neither the removed codons nor the triplet left at the start of the coding sequence, so they do not hold vep-rs's term to the bounds of section 3. Under `--most_severe` or `--summary` a record has one row and no `Transcript` row. Where every coding transcript that holds the removal has no 5' UTR, that row carries vep-rs's `start_lost`, which Ensembl VEP's row lacks, and the comparisons leave it a difference.

The plugin tuple scorer (`tools/vep-plugin-pairs`) holds the class with the same test of a pair and the same `IMPACT` and HGVSp cells. No plugin field follows these terms, so it sets no plugin value aside under the class.

## 2. The same removal in five whole-genome datasets

**In short.** gnomAD v4.1 and v2.1.1, 1000 Genomes high coverage and ClinVar on both assemblies hold 14 records, each with a single alternate allele, that remove whole codons from the first coding base of a transcript without a 5' UTR. On 17 of their 20 such rows vep-rs writes `start_lost`, which Ensembl VEP's rule gives on none.

| Dataset | Records | Rows | Rows with `start_lost` |
| --- | --- | --- | --- |
| gnomAD v4.1 (GRCh38) | 5 | 6 | 5 |
| gnomAD v2.1.1 (GRCh37) | 2 | 3 | 2 |
| 1000 Genomes high coverage (GRCh38) | 1 | 1 | 0 |
| ClinVar (GRCh38) | 4 | 7 | 7 |
| ClinVar (GRCh37) | 2 | 3 | 3 |

The rows are counted with vep-rs on the Ensembl release 116 caches under the default flags, and the last column is vep-rs's. The committed corpora hold five of the records. Ensembl VEP 116.3's stored rows for them, nine in all, read `inframe_deletion` without `start_lost`.

One gnomAD v4.1 record removes 16:11942043-11942045, the start codon of the MANE Select transcript `ENST00000399147`, which holds no 5' UTR. Another removes 21:44817966-44817968, the start codon of the reverse-strand `ENST00000397893`; the committed corpus `tests/golden/116/GRCh38` holds that record too.

The three rows without the term are one removal, the first two codons of `ENST00000536027`, which three of the datasets hold. Its coding sequence opens `ATG CCT ATG`, so `ATG` stands at its start afterwards. Ensembl VEP writes no start term there, and vep-rs writes `start_retained_variant`, the entry [`start_codon_kept_without_five_prime_utr`](start-codon-kept-without-five-prime-utr.md).

A ClinVar record removes 9,912 bases, 7:5999108-6009019 on GRCh38, from the first coding base of `ENST00000382321` through its first exons and introns. It takes the first 705 coding bases, 235 whole codons, and the row's `Amino_acids` holds no alternate residue. Both engines write `splice_acceptor_variant`, `splice_donor_variant`, `inframe_deletion`, `splice_donor_5th_base_variant` and `intron_variant`, and vep-rs adds `start_lost`. `IMPACT` is `HIGH` on both sides there.

## 3. What the rule does not read

- A removal after which `ATG` stands at the start of the coding sequence. Ensembl VEP's test behind a 5' UTR reads that as the start in place, and vep-rs reads it so on every transcript. Behind a 5' UTR the length test finds the UTR and `ATG` in place (lines 1061-1066), and both engines write `start_retained_variant` beside `inframe_deletion`. Without one Ensembl VEP writes `inframe_deletion` alone and vep-rs `start_retained_variant` beside it, the entry [`start_codon_kept_without_five_prime_utr`](start-codon-kept-without-five-prime-utr.md).
- A removal after which `ATA` stands at the start of the coding sequence, on a transcript of a mitochondrial genetic code under which `ATA` codes methionine. The first residue is kept. Both engines read a substitution of `ATG` by `ATA` on such a transcript as no lost start (`VariationEffect.pm`, lines 880-881), and vep-rs reads the removal the same way.
- A removal after which the triplet at the start of the coding sequence holds a base that is none of A, C, G and T. Its residue is not known, and the row stays as Ensembl VEP writes it.
- A coding sequence that opens on a codon that is no methionine start codon of its genetic code, on a transcript without `cds_start_NF`: `GCA`, or `CTG` under the standard code. A model's coding sequence can be cut off at its start without carrying the flag, and it then opens on any codon. The row stays as Ensembl VEP writes it.
- A transcript that carries `cds_start_NF`: neither engine writes a start term on it (`VariationEffect.pm`, line 975).
- A deletion that is not whole codons from the first coding base: one or two bases of the start codon, or three bases from its second. The frame shifts or the alternate peptide holds a residue, and Ensembl VEP's first and third ways decide, as in vep-rs.
- A removal that begins on the second or third base of the start codon and leaves fewer than three coding bases, on a coding sequence whose length is no multiple of three. Its row has no alternate residue either, and the rule of this page does not read it.
- A deletion that begins before the transcript's first base. Its span has no cDNA start, both engines write no start term for it, and the rule of this page does not read it.
- An insertion, and a deletion-insertion. Their alternate peptide holds residues, and Ensembl VEP's first and third ways decide, as in vep-rs.
