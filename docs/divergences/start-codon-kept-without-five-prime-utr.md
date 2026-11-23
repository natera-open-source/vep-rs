# A start codon kept on a transcript with no 5' UTR

How each engine reads an insertion or deletion at a start codon that leaves the codon in place, on a transcript that holds no 5' UTR.

**The difference in one sentence.** For an insertion or deletion at the start codon of a transcript with no 5' UTR that leaves a start codon there, Ensembl VEP 116.3 writes `start_retained_variant` only where the edited sequence ends with the whole coding sequence, and vep-rs always writes it. On a row it leaves without `start_retained_variant`, Ensembl VEP writes `start_lost` where the allele shifts the frame or where its peptide test holds, and vep-rs writes no `start_lost` there. Behind a 5' UTR Ensembl VEP writes `start_retained_variant` for every one of these alleles.

**What each engine does.** Ensembl VEP decides the start term of an allele that changes the length of the sequence in one routine. The routine joins the 5' UTR and the coding sequence and makes the record's edit in them. It reads the start as kept where the UTR is unchanged and the three bases after it read `ATG`. That test stands inside a condition that asks for a 5' UTR.

On a transcript with none, the routine reads a deletion as an altered start, because the edited sequence is shorter than the coding sequence. It reads an insertion as an altered start unless the edited sequence still ends with the whole coding sequence.

vep-rs makes the kept-start test on such a transcript as well, at the transcript's first base. It makes the edit in the coding sequence and reads the start as kept where the first three bases read `ATG`. Its other two tests, of the length and of the end of the edited sequence, are Ensembl VEP's. On a transcript with no 5' UTR vep-rs therefore writes the row it writes for the same allele behind a 5' UTR. That is Ensembl VEP's own row there, but for the `start_lost` that the entry [`start_cooccurrence_swap`](defect-b.md) sets aside.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the modules the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1028-1075, `_ins_del_start_altered`: the routine. Lines 1053 and 1058 join the 5' UTR and the coding sequence and make the edit.
- The same file, lines 1062-1065: the kept-start test, inside `if ($utr)` at line 1061. A transcript whose coding sequence opens on its first base has no 5' UTR, and the test is not made.
- The same file, line 1069: an edited sequence shorter than the coding sequence is answered as altered. Line 1071 compares the last bases of the edited sequence, as many as the coding sequence holds, with the coding sequence, and answers altered where they differ.
- The same file, lines 947-963, `start_retained_variant`: for an allele that changes the length, the negation of the routine's answer (line 961).
- The same file, lines 851-904, `start_lost`: line 870 writes the term where the routine answers altered and the row is no in-frame insertion or deletion. Lines 878-882 write it from the peptides, where the row's first residue is residue 1 and the alternate peptide neither opens nor ends with the reference peptide. Line 875 returns before that test where the alternate peptide is empty or a single unknown residue. Line 871 asks `_inv_start_altered` (lines 906-945), which returns at line 926 on a transcript with no 5' UTR.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2161-2164: a row that carries `start_lost` gets `?` for the new peptide in its HGVSp.

**The class it produces.** `start_codon_kept_without_five_prime_utr`: a `Transcript` row of an insertion or deletion whose cDNA position equals its CDS position and starts on one of the first three coding bases. vep-rs's terms are Ensembl VEP's without `start_lost` and with `start_retained_variant`. `IMPACT` and, under `--hgvs`, the HGVSp follow each side's terms. Section 1.8 gives what each comparison reads.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a record with an allele that adds or removes bases. The five chromosome 21 record sets the index calls the cuts hold 2,154,724 such records among 15,646,481, the index's 137,713 per million. The rate is an upper bound: the allele must also overlap the start codon of a transcript without a 5' UTR and leave a start codon at the start of its coding sequence. Section 2 counts the rows that do in six whole-genome datasets.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Four alleles at the start codon shared by ENST00000303645 and ENST00000339944

**In short.** Four records written for this page change the length of the coding sequence at 21:32298972-32298974 on GRCh38, the start codon two transcripts of MRAP share. Each leaves `ATG` as the first three coding bases. `ENST00000303645` holds a 5' UTR, and both engines write `start_retained_variant` on it for all four. `ENST00000339944` starts on the codon, and there Ensembl VEP 116.3 writes `start_lost` for three of the four and no start term for the fourth, where vep-rs writes `start_retained_variant` for all four.

### 1.1 What a reader needs to know first

A transcript is a model of a messenger RNA. Its coding sequence begins at a start codon, three bases that read `ATG`, and the bases the model holds before that codon are its 5' UTR. A model with no 5' UTR begins on the codon itself. Its coding sequence opens on the transcript's first base, so a base's cDNA position and its CDS position are the same number.

Ensembl VEP has two Sequence Ontology terms for a change at a start codon. `Constants.pm` of `ensembl-variation` words them so:

- `start_lost` (SO:0002012): "A codon variant that changes at least one base of the canonical start codon".
- `start_retained_variant` (SO:0002019): "A sequence variant where at least one base in the start codon is changed, but the start remains".

`start_lost` carries the impact `HIGH` and `start_retained_variant` the impact `LOW`. A row's `IMPACT` is that of its most severe term.

An allele that adds or removes a number of bases that is no multiple of three shifts the reading frame, and its term is `frameshift_variant`. One that adds or removes whole codons is an `inframe_insertion` or an `inframe_deletion`. On this page whole codons means a multiple of three bases, wherever in a codon the change starts.

### 1.2 The records

Both transcripts lie on the forward strand. `ENST00000303645` begins at 21:32298855, 117 bases before the codon, so the codon is its cDNA positions 118 to 120. `ENST00000339944` begins at 21:32298972, the codon's first base, so the codon is its cDNA positions 1 to 3. On both the codon is CDS positions 1 to 3, and the codon after it, `GCC`, CDS positions 4 to 6. Neither transcript carries `cds_start_NF`.

The four records, by their identifiers:

| Record | Change | CDS position |
| --- | --- | --- |
| `frameshift_insertion_inside_the_codon` | four bases written between 21:32298972 and 21:32298973 | `1-2` |
| `whole_codon_insertion_inside_the_codon` | three bases written between the same two positions | `1-2` |
| `third_base_removed_before_its_repeat` | 21:32298974 removed, one base | `3` |
| `three_codons_removed_from_the_third_base` | 21:32298974-32298982 removed, nine bases | `3-11` |

The record lines, which the entry's reproduction holds as the `text` of its `[input]` table (1.7 names it):

```
21 32298972 frameshift_insertion_inside_the_codon A ATGCC . . .
21 32298972 whole_codon_insertion_inside_the_codon A ATGC . . .
21 32298973 third_base_removed_before_its_repeat TG T . . .
21 32298973 three_codons_removed_from_the_third_base TGGCCAACGG T . . .
```

The columns are the chromosome, the position, the identifier, the reference bases and the alternate bases; the last three are unset. Both engines run `--hgvs` on the release 116 cache and the reference of the corpus `tests/golden/116/GRCh38-plugin-loftee-fixture`.

### 1.3 What each engine writes

Each table gives the rows of the two transcripts, Ensembl VEP 116.3's first.

**A frameshift insertion inside the codon.** The record writes four bases between the codon's first and second bases. The four open with `TG`, the codon's second and third bases. The first three coding bases therefore read `ATG`, and the frame shifts after them.

| Engine | Feature | Consequence | HGVSp |
| --- | --- | --- | --- |
| Ensembl VEP | `ENST00000303645` | `frameshift_variant,start_retained_variant` | `p.Ala2ProfsTer15` |
| Ensembl VEP | `ENST00000339944` | `frameshift_variant,start_lost` | `p.Ala2_?1` |
| vep-rs | `ENST00000303645` | `frameshift_variant,start_retained_variant` | `p.Ala2ProfsTer15` |
| vep-rs | `ENST00000339944` | `frameshift_variant,start_retained_variant` | `p.Ala2ProfsTer15` |

**An insertion of one whole codon inside the codon.** The record writes three bases at the same place. They open with `TG` as well, so `ATG` stands at the start of the coding sequence and a new codon, one of leucine, follows it.

| Engine | Feature | Consequence | HGVSp |
| --- | --- | --- | --- |
| Ensembl VEP | `ENST00000303645` | `inframe_insertion,start_retained_variant` | `p.Met1_Ala2insLeu` |
| Ensembl VEP | `ENST00000339944` | `inframe_insertion` | `p.Met1_Ala2insLeu` |
| vep-rs | `ENST00000303645` | `inframe_insertion,start_retained_variant` | `p.Met1_Ala2insLeu` |
| vep-rs | `ENST00000339944` | `inframe_insertion,start_retained_variant` | `p.Met1_Ala2insLeu` |

**The codon's third base removed before its repeat.** The record removes 21:32298974, the codon's third base, `G`. The base after it, 21:32298975, is `G` as well and takes its place. The first three coding bases therefore read `ATG`, and the frame shifts after them.

| Engine | Feature | Consequence | HGVSp |
| --- | --- | --- | --- |
| Ensembl VEP | `ENST00000303645` | `frameshift_variant,start_retained_variant` | `p.Ala2ProfsTer30` |
| Ensembl VEP | `ENST00000339944` | `frameshift_variant,start_lost` | `p.Ala2?` |
| vep-rs | `ENST00000303645` | `frameshift_variant,start_retained_variant` | `p.Ala2ProfsTer30` |
| vep-rs | `ENST00000339944` | `frameshift_variant,start_retained_variant` | `p.Ala2ProfsTer30` |

**Three codons removed from the codon's third base.** The record removes 21:32298974-32298982, nine bases from the codon's third base on. The base after them, 21:32298983, is `G` as well and takes the codon's third place. The first three coding bases therefore read `ATG`, and the protein lacks its residues 2 to 4.

| Engine | Feature | Consequence | HGVSp |
| --- | --- | --- | --- |
| Ensembl VEP | `ENST00000303645` | `start_lost,inframe_deletion,start_retained_variant` | `p.AlaAsnGly2_?4` |
| Ensembl VEP | `ENST00000339944` | `start_lost,inframe_deletion` | `p.AlaAsnGly2_?4` |
| vep-rs | `ENST00000303645` | `inframe_deletion,start_retained_variant` | `p.Ala2_Gly4del` |
| vep-rs | `ENST00000339944` | `inframe_deletion,start_retained_variant` | `p.Ala2_Gly4del` |

`IMPACT` is `HIGH` on the eight rows of the two frameshifts and `MODERATE` on the four of the insertion of a whole codon. On the deletion of three codons it is `HIGH` on Ensembl VEP's two rows and `MODERATE` on vep-rs's two. The HGVSp is on `ENSP00000306697.5` for `ENST00000303645` and on `ENSP00000343661.4` for `ENST00000339944`.

The engines agree on the positions of every row. The `cDNA_position` on `ENST00000303645` is `118-119` for the two insertions, `120` for the one base removed and `120-128` for the nine. On `ENST00000339944` it equals the `CDS_position` of 1.2.

On the transcript with a 5' UTR Ensembl VEP writes `start_retained_variant` for each of the four alleles. For the deletion of three codons it writes `start_lost` beside it. That pair of start terms on one row is the entry [`start_cooccurrence_swap`](defect-b.md), under which vep-rs keeps `start_retained_variant` alone.

### 1.4 Why Ensembl VEP writes that

On `ENST00000303645` the routine joins the 117 bases of the 5' UTR and the coding sequence and makes the edit. No record changes a base of the UTR, and the three bases after it read `ATG`. Line 1065 answers that the start is not altered, and `start_retained_variant` is the negation of that answer (line 961).

`start_lost` has another way to a row, the peptide test of lines 878-882. It holds where the row's first residue is residue 1 and the alternate peptide neither opens nor ends with the reference peptide. On the deletion of three codons it holds, so the row carries both start terms: the reference peptide is the four residues of the changed codons and the alternate peptide the first of them alone. On the other three alleles it does not. The alternate peptide of each insertion opens with the reference peptide (line 881). That of the one base removed is a single unknown residue, on which line 875 returns before the test.

On `ENST00000339944` the transcript has no 5' UTR, and the test of lines 1062-1065 is not made. For the two deletions the edited sequence is shorter than the coding sequence, and line 1069 answers altered. For the two insertions line 1071 compares the end of the edited sequence with the coding sequence. Neither record's inserted bases end with the codon's first base, `A`, so the coding sequence does not stand whole behind them, and the line answers altered.

An insertion whose last bases repeat the codon's bases before it leaves the whole coding sequence at the end of the edited sequence. Line 1071 answers that one as kept, whatever the first three bases then read. vep-rs makes the same comparison, and both engines write `start_retained_variant` for it on either transcript.

With that answer no row of the four carries `start_retained_variant`. Line 870 writes `start_lost` on the two frameshift rows. It excepts an in-frame insertion and an in-frame deletion. Line 871 asks `_inv_start_altered` (lines 906-945), a second test of the edited sequence, which returns at line 926 without an answer of altered on a transcript with no 5' UTR. That leaves the peptide test, which holds on the deletion of three codons and writes `start_lost`. On the insertion of a whole codon the alternate peptide opens with the reference peptide (line 881), and the row carries no start term.

### 1.5 Why vep-rs's row is the correct one

The two transcripts read the same three bases, 21:32298972-32298974, as their start codon. Each record leaves `ATG` as the first three bases of both coding sequences, as 1.3 derives from the record's positions. After each of the four changes, translation begins at the same codon on both transcripts.

`start_retained_variant` names a change of at least one base in the start codon after which the start remains. Each record changes the sequence at the codon, and the start remains. The term therefore belongs on both rows. `start_lost` and `start_retained_variant` exclude each other, since one change does not both lose and keep a start codon, so the row carries no `start_lost`.

No record changes any of the 117 bases that `ENST00000303645` holds before the codon. Those bases change nothing in what a record leaves at the start of the coding sequence. Behind them Ensembl VEP's routine reads the same edited bases and answers that the start is kept.

For a substitution the same module makes its kept-start test with or without a 5' UTR. `_snp_start_altered` joins an empty UTR where the transcript has none (line 1012) and reads the three bases at the start of the coding sequence (lines 1021-1022).

### 1.6 The cells that follow the terms

`IMPACT` is the impact of a row's most severe term. Beside `frameshift_variant`, a `HIGH` term, the row reads `HIGH` with either start term. Beside `inframe_deletion`, a `MODERATE` term, Ensembl VEP's row reads `HIGH` for its `start_lost`, and vep-rs's reads `MODERATE`.

Under `--hgvs` the HGVSp follows `start_lost`. On a row that carries the term, Ensembl VEP writes the notation of a lost start, with `?` for the new peptide (`TranscriptVariationAllele.pm`, lines 2161-2164). Those are `p.Ala2_?1`, `p.Ala2?` and `p.AlaAsnGly2_?4` in the tables of 1.3.

On a row without `start_lost` the same routine writes the notation of the change itself, and vep-rs's rows hold it. That is a frameshift notation, `p.Ala2ProfsTer15` or `p.Ala2ProfsTer30`, or the `del` notation `p.Ala2_Gly4del`. The two frameshift notations are the ones both engines write on `ENST00000303645`. The row of the insertion of a whole codon carries no `start_lost` on either side, and its HGVSp is the same on both.

### 1.7 Try it yourself

Put the four record lines of 1.2 in `records.vcf`, with a tab between the columns, after the two header lines [Reproducing these examples](reproducing.md) gives. The command as vep-rs takes it, from the repository's root:

```
vep -i records.vcf -o out.txt --offline --hgvs \
  --json_cache tests/golden/116/GRCh38-plugin-loftee-fixture/json_cache \
  --fasta tests/golden/116/GRCh38-plugin-loftee-fixture/reference.fa.gz \
  --species homo_sapiens --assembly GRCh38
grep -e ENST00000303645 -e ENST00000339944 out.txt
```

Ensembl VEP 116.3 takes the same records and `--hgvs` with the cache options of that page. The entry's reproduction is `scripts/divergences/reproductions/start_codon_kept_without_five_prime_utr.toml`. It holds the four record lines, runs both engines on them and holds the eight rows of the two transcripts to both recorded values.

### 1.8 How the comparisons treat this class

The format comparator (`tools/vep-formats-compare`) takes a `Transcript` row under the class `start_codon_kept_without_five_prime_utr` where all of these hold:

- vep-rs's terms are Ensembl VEP's without `start_lost`, where Ensembl VEP's hold it, and with `start_retained_variant`, every other term being equal;
- Ensembl VEP's terms hold no `start_retained_variant` and vep-rs's no `start_lost`;
- both rows hold `frameshift_variant`, `inframe_insertion` or `inframe_deletion`;
- on both rows the cDNA position equals the CDS position, which shows a coding sequence that opens on its transcript's first base;
- the CDS position starts at 1, 2 or 3, a base of the start codon.

An insertion written between coding bases 3 and 4 touches no base of the codon: neither engine tests the start for it, and no pair of this kind arises on its row.

On such a row it sets aside the `Consequence` and an `IMPACT` that is the one each side's own terms give. Where Ensembl VEP's row holds `start_lost`, it sets aside an HGVSp that holds `?` on Ensembl VEP's side and none on vep-rs's. In the JSON format the record's `most_severe_consequence` follows the rows' terms and is set aside with them. Any other cell that differs stays a difference.

The tuple comparison (`tools/vep-outputs-compare`) scores consequence terms. It reads the two rows behind a pair whose terms differ as the class's do, and sets the pair aside where both rows meet the conditions above. The structural-variant comparator (`tools/vep-sv-compare`) applies the same rule. The plugin tuple scorer (`tools/vep-plugin-pairs`) holds the class with the same test of a pair and sets aside the same cells.

The corpus builder (`tools/vep-golden-corpus`) records each row of the class under the class's name in the manifest of its corpus, and the corpus comparison sets those rows aside. A row off the shape stays a difference charged to vep-rs.

## 2. The same shape in six whole-genome datasets

**In short.** The genome-wide ClinVar dataset holds 24 rows of the class on GRCh38 and 20 on GRCh37, and four population datasets hold 106 more. One ClinVar record, 3780850, removes the first 250 coding bases of FANCC and leaves a start codon at the start of what remains.

### 2.1 A deletion of 1,860 bases in FANCC

ClinVar 3780850 removes 9:95247432-95249291 on GRCh38, 1,860 bases of FANCC. The deletion runs from the first base of the start codon to the last base of the next exon, and the coding sequence reads a start codon after it. The committed corpus `tests/golden/116/GRCh38` holds the record.

Nineteen transcripts of the gene hold a 5' UTR before the codon, and both engines write `start_retained_variant` beside `frameshift_variant` on them. `ENST00000649872` and `ENST00000649334` hold none: their rows print `1-250` as the `cDNA_position` and as the `CDS_position`. On `ENST00000649872` the engines write:

| Engine | Consequence |
| --- | --- |
| Ensembl VEP | `splice_acceptor_variant,splice_donor_variant,frameshift_variant,start_lost,splice_donor_5th_base_variant,intron_variant` |
| vep-rs | `splice_acceptor_variant,splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant,start_retained_variant,intron_variant` |

The row of `ENST00000649334` differs in the same way and ends with `NMD_transcript_variant` on both sides. vep-rs's terms are the ones both engines write on the nineteen transcripts with a 5' UTR. `IMPACT` is `HIGH` on every one of these rows.

### 2.2 The counts

| Dataset | Records | Rows |
| --- | --- | --- |
| ClinVar (GRCh38) | 16 | 24 |
| ClinVar (GRCh37) | 11 | 20 |
| gnomAD v4.1 genomes (GRCh38) | 38 | 53 |
| gnomAD v2.1.1 genomes (GRCh37) | 31 | 43 |
| 1000 Genomes high coverage (GRCh38) | 4 | 7 |
| 1000 Genomes Phase 3 (GRCh37) | 2 | 3 |

A counted row is one whose terms vep-rs's kept-start test changes on a transcript with no 5' UTR. The rows are found by running vep-rs with and without the test on the Ensembl release 116 caches, and the default flags and `--hgvs` give the same rows. On each of the 44 ClinVar rows the row without the test is the row of Ensembl VEP 116.3's output. The rows of the four population datasets are counted the same way and are not read against an Ensembl VEP output here.

The ClinVar rows are of four shapes, counted on GRCh38 first and on GRCh37 second:

- a frameshift deletion, 9 and 5 rows;
- a frameshift insertion, 12 and 6 rows;
- a deletion of whole codons, 1 and 0 rows;
- an insertion of whole codons, 2 and 9 rows.

Each frameshift insertion is written between two bases of the start codon, and its inserted bases open with the codon's remaining bases, so the codon stands whole before them.

The four population datasets hold 106 rows of 75 records. On 58 of the rows the test puts `start_retained_variant` in the place of `start_lost`. Of those, 37 are `start_lost,inframe_deletion` rows, whose `IMPACT` goes from `HIGH` to `MODERATE`. On the other 48 rows it adds `start_retained_variant` to a row with no start term.

### 2.3 A frameshift deletion with both start terms behind a 5' UTR

ClinVar 3336827 removes 11:57598273-57598280 on GRCh38, eight bases from the third base of a start codon on. On the transcripts that hold a 5' UTR Ensembl VEP writes both start terms beside `frameshift_variant`, as in `frameshift_variant,start_lost,start_retained_variant`. Line 1065 keeps the start, and the peptide test of lines 878-882 holds as well. The changed codons keep four bases, `ATG` and one more, so the alternate peptide is the first residue and an unknown one, which neither opens nor ends with the reference peptide of four residues. vep-rs writes `start_retained_variant` alone of the two there, under `start_cooccurrence_swap`.

`ENST00000531133` holds no 5' UTR, and the row prints `3-10` as both positions. Ensembl VEP writes `frameshift_variant,start_lost,NMD_transcript_variant` on it and vep-rs `frameshift_variant,start_retained_variant,NMD_transcript_variant`.
