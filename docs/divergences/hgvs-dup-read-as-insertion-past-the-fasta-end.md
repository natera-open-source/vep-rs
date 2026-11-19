# The HGVSc of a duplication on a reference FASTA that ends inside the transcript

**The defect in one sentence.** Under `--hgvs`, on a reference FASTA whose chromosome ends before the end of a transcript, Ensembl VEP 116.2 fetches the sequence past the file's last base as a single `N` and keeps the transcript's sequence that many bases short. An insertion that duplicates the bases beside it is then written as an insertion, `c.745_746insTGC`, where the complete chromosome gives `c.743_745dup`. vep-rs writes the duplication on either reference.

**Why this is a defect.** The notation of an insertion depends on the bases beside the insertion point and on nothing else. Those bases lie inside the file: here the reference ends 144 bases past the insertion, and only the transcript's far end, 21,170 bases on, lies beyond it. Ensembl VEP writes the duplication for the same record on the complete chromosome, and the shorter file changes no base the notation reads.

**What vep-rs does instead.** vep-rs reads the bases the notation needs at their genomic coordinates from the FASTA, never the transcript's whole span, so a reference that ends inside the transcript but past those bases gives the notation the complete chromosome gives.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/Utils/FastaSequence.pm`, lines 619-624, `_raw_seq`: a request whose first base lies past the FASTA's last base returns one `N`, whatever length was asked; lines 626-639 pad a request that merely runs over the end with one `N` per missing base.
- The same file, lines 409-520, the sequence cache behind a slice's `seq`: a request overlapping a cached region fetches only the bases the region lacks (lines 437-476). When they lie past the file's end the region gains that one `N` while its recorded end moves to the requested end (line 474), so its string is shorter than its span. Every later read cuts the string at an offset computed from the recorded coordinates (lines 491-498); on a reverse-strand region the offset counts from the end, where the bases are missing, so the cut comes back short, or displaced by the missing length.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1472-1483, `hgvs_transcript`: the check that the variant lies inside the slice reads the slice's coordinates, not its sequence (line 1472), and the sequence as returned is passed to `hgvs_variant_notation` with the shifted insertion's coordinates (lines 1476-1478).
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 575-596, `hgvs_variant_notation`: the duplication test reads the inserted length of bases after and before the insertion point from that string (lines 577-579). Outside the string `substr` yields nothing, and a displaced cut yields other bases, so neither equals the allele and the notation is an insertion (lines 592-594).

**What correcting it would change.** The row would carry the duplication Ensembl VEP writes for the same record on the complete chromosome, and the `substr outside of string` warnings its warnings file collects on such a run (`Utils/Sequence.pm` lines 526 and 577) would go with it.

**The class it produces.** `hgvs_dup_read_as_insertion_past_the_fasta_end`, documented: a differing `HGVSc` on an insertion's row, on a run whose reference FASTA ends inside the transcript. Ensembl VEP's notation is an `ins` between two positions; vep-rs's is a `dup` of as many bases as the allele, ending at the first of them. No committed corpus runs a FASTA shorter than its chromosome, so no manifest names the class and no comparator applies it. The concordance figures score consequence sets and never see the field.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. `CAG` inserted after 21:46302556 on a reference that ends at 46302700

**In short.** `21 46302556 . G GCAG` inserts `CAG` after 46302556, inside a `CAG` run, on the reverse-strand transcript ENST00000291691, which spans 46301130-46323870. Under `--hgvs` on the complete GRCh38 chromosome 21 both engines write `HGVSc=ENST00000291691.12:c.743_745dup` with `HGVS_OFFSET=-4`. On a copy of the chromosome cut after base 46302700, 144 bases past the insertion and 21,170 bases before the transcript's end, Ensembl VEP writes `c.745_746insTGC` and vep-rs `c.743_745dup` again.

### 1.1 What a reader needs to know first

For HGVS, an insertion is placed at the most 3' position along the transcript at which inserting the bases gives the same sequence, and `HGVS_OFFSET` reports how far that moved it from the input's position. An insertion whose bases repeat the bases just before them is written as a duplication of those bases. On a reverse-strand transcript the 3' direction is towards lower genomic coordinates, and every base is read as its complement.

Ensembl VEP reads a transcript's reference as one string over the transcript's whole genomic span, position 1 at the transcript's 5' end, and takes the notation's bases as substrings of it. The string is served from a cache of sequence regions already read, which grows a region by the bases a new request lacks.

A reference FASTA may end inside a chromosome: a file holding the first part of a chromosome, or a region extracted from one and renamed. Ensembl VEP takes the chromosome's length from the FASTA index.

### 1.2 The record

The record is a hand-written insertion on GRCh38, annotated with the release-116 cache and a chromosome 21 FASTA:

```
21	46302556	ins_cag_46302556	G	GCAG	.	.	.
```

The bases `21:46302550-46302565` read `TCTGCAGCAGCAGCAG`. ENST00000291691 lies on the reverse strand over 46301130-46323870, 22,741 bases, and the insertion falls in its coding sequence, between CDS positions 741 and 742 before the shift. Both engines run with `--hgvs` twice: against the complete chromosome 21, 46,709,983 bases, and against the same file cut after base 46302700, its index rewritten to that length.

### 1.3 What each engine writes

The transcript's row in the default format. Ensembl VEP on the cut reference:

```
ins_cag_46302556	21:46302556-46302557	CAG	ENSG00000160298	ENST00000291691	Transcript	inframe_insertion	1873-1874	741-742	247-248	-/L	-/CTG	-	IMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000291691.12:c.745_746insTGC;HGVSp=ENSP00000291691.8:p.Gln249delinsProGlu;HGVS_OFFSET=-4
```

vep-rs on either reference:

```
ins_cag_46302556	21:46302556-46302557	CAG	ENSG00000160298	ENST00000291691	Transcript	inframe_insertion	1873-1874	741-742	247-248	-/L	-/CTG	-	IMPACT=MODERATE;STRAND=-1;HGVSc=ENST00000291691.12:c.743_745dup;HGVSp=ENSP00000291691.8:p.Leu248dup;HGVS_OFFSET=-4
```

Ensembl VEP on the complete chromosome writes `HGVSc=ENST00000291691.12:c.743_745dup`, the rest of its row as above. The consequence, the positions, the codons and `HGVS_OFFSET` agree everywhere. The HGVSp differs on both references for another reason, the rotation of a shifted allele by the transcript's strand alone, the class [`shifted_allele_rotation_by_strand`](../output/vep-and-tab.md#allele-on-a-coding-row-after-a-shift).

### 1.4 Why Ensembl VEP writes that

The 3' shift reads a window of 1,000 bases on each side of the insertion first, 46301557-46303556, and the sequence cache keeps it as a region. The 856 bases past the file's end come back as `N`, one per base, so the region holds its 2,000 bases. The transcript's slice, 46301130-46323870, then overlaps that region. The 427 bases before it are fetched; the 20,314 bases after it, 46303557-46323870, begin past the file's end, so `_raw_seq` returns one `N`. The region's end moves to 46323870 while its string holds 2,428 bases for a span of 22,741.

The slice's sequence is cut from that string at the offset the coordinates give, 22,741 bases asked and 2,428 returned, and `hgvs_transcript` passes it on: its guard compares the shifted position with the slice's span, not with the string's length. The duplication test looks for `TGC` three bases after and three bases before slice position 21,318, the shifted insertion point (c.745, genomic 46302553). Both reads fall past the string's 2,428 bases, `substr` returns nothing, neither equals `TGC`, and the notation is the insertion between c.745 and c.746.

Of the 29 rows with an HGVSc Ensembl VEP writes for this record, 19 change between the two references, every one on a reverse-strand transcript. The nine forward-strand transcripts read the region's reverse complement, which is counted from the region's start, where no base is missing, and keep their `dup`. One reverse-strand transcript, ENST00000397679, asks for a cut that lies wholly past the short string, gets nothing back, fetches its slice afresh with the `N` padding and keeps `dup` too. The other reverse-strand transcripts read the string at offsets displaced by the missing bases, find other bases there, and write `ins`.

Source (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/Utils/FastaSequence.pm`, lines 619-624: the one `N` for a request that starts past the file's end; lines 457-476: the region grown by it, its end moved to the request's; lines 491-498: the cut at the coordinates' offset.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, line 1472: the guard on the slice's span; lines 1474-1483: the call with the slice's sequence.
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 575-596: the duplication test and the insertion it falls to.

### 1.5 Why vep-rs's row is the correct one

The sequence decides. Over `21:46302550-46302565`, `TCTGCAGCAGCAGCAG`, inserting `CAG` after 46302556 gives `TCTGCAGCAGCAGCAGCAG`, and so does inserting `GCA` after 46302552, four bases towards the transcript's 3' end and the furthest placement that does: `HGVS_OFFSET=-4` on both engines. Read on the transcript's strand the inserted bases are `TGC`, placed after c.745, genomic 46302553. The bases c.743_745 are 46302555-46302553, `GCA` on the genome and `TGC` on the transcript, the very bases the insertion repeats, so the notation is `c.743_745dup`. Every base this reading needs lies between 46302550 and 46302565, inside both references; the file's end at 46302700 changes none of them.

### 1.6 Try it yourself

Cut the chromosome after base 46302700, keeping the record's name `21`, and index the cut; then annotate the record of 1.2 with `--hgvs` against the GRCh38 cache and the cut file, and again against the complete chromosome. vep-rs's command is below; Ensembl VEP takes the same record with its own cache and the same FASTA under `--offline --cache`:

```
samtools faidx reference.fa 21:1-46302700 | sed '1s/.*/>21/' > cut.fa
samtools faidx cut.fa
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta cut.fa \
  --species homo_sapiens --assembly GRCh38 --hgvs
```

| Field | Ensembl VEP, cut reference | Ensembl VEP, complete chromosome | vep-rs, either reference |
| --- | --- | --- | --- |
| `HGVSc` | `ENST00000291691.12:c.745_746insTGC` | `ENST00000291691.12:c.743_745dup` | `ENST00000291691.12:c.743_745dup` |
| `HGVS_OFFSET` | `-4` | `-4` | `-4` |

On the cut reference Ensembl VEP's warnings file carries `substr outside of string` from `Utils/Sequence.pm` lines 526 and 577 and `Use of uninitialized value` from the lines that read the result; on the complete chromosome it writes no warnings file. The bases can be read from either FASTA: `samtools faidx <fasta> 21:46302550-46302565` prints `TCTGCAGCAGCAGCAG`.

### 1.7 How the comparator treats this class

No comparator applies the class. The committed corpora's references run their chromosomes to the assembly's full length, as their FASTA indexes record, so no committed flag set meets the shape and no manifest names the class. The concordance figures score consequence sets and never see the field.
