# Two transcripts of a merged cache that carry one `dbID`

**The defect in one sentence.** Under `--merged`, Ensembl VEP keeps one transcript per `dbID` among the transcripts a buffer loads. A merged cache numbers its Ensembl-set and its RefSeq-set transcripts apart, so one `dbID` can name a transcript of each set. The transcript read second loses every row, and a record on no other transcript is written as `intergenic_variant`. vep-rs keeps every transcript the cache distinguishes by set and `dbID`.

A buffer is the run of records Ensembl VEP annotates together: up to `--buffer_size` records, 5,000 by default, of one sequence. From the region files its records reach, it loads the transcripts that lie between its lowest and its highest record, 5,000 bases either side.

**Why this is a defect.** The two transcripts are different models: another set, another stable identifier, another span. The step that drops one exists for a transcript stored in two region files. It asks for the `dbID` alone, which is unique inside one set only. Ensembl VEP writes a row of the dropped transcript on the cache of that transcript's set alone, and on the merged cache in a buffer that does not load the pair's other transcript. The rows of a record therefore change with the other records of the run. Ensembl VEP's own step for regulatory features keys on the feature type and the `dbID`, because each type has its own `dbID` numbering.

**What vep-rs does instead.** vep-rs reads a chromosome's region files once and keeps one copy of each transcript, keyed on its `dbID` and on the set the cache records for it (`_source_cache` in Ensembl VEP's cache, `source_cache` in its conversion). A transcript stored in two region files is annotated once, and two transcripts of different sets under one `dbID` are both annotated. A record's rows are the same whatever else the input holds and at every `--buffer_size`.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `AnnotationType/Transcript.pm`, lines 266-300, `merge_features`: line 283 reads the transcript's `dbID`, lines 284-287 skip a transcript whose `dbID` an earlier one of the list carried, and line 297 records it. Neither the set (`_source_cache`) nor the stable identifier enters the key.
- `AnnotationType/Transcript.pm`, lines 236-248, `filter_transcript`: without `--all_refseq`, a RefSeq-set transcript of a merged cache is kept only where its identifier is a RefSeq accession, or on `MT` a mitochondrial gene identifier. Line 241 of `AnnotationSource/Cache/Transcript.pm` calls it as a region file is read, so a transcript it drops never reaches `merge_features`. Line 244 also keeps a transcript by the identifier of its `display_xref`, a field the cache dump deletes from every transcript it writes (`Pipeline/DumpVEP/Dumper/Core.pm`, lines 113 and 155-173).
- `AnnotationSource.pm`, lines 109-145, `get_all_features_by_InputBuffer`: the transcripts of every region file the buffer's records reach, less those outside the buffer's lowest to highest position with the flank, go through that method. Lines 219-263 of the same file list the region files in the order the records first ask for them.
- `InputBuffer.pm`, lines 195 and 224: a buffer ends where the sequence name changes, or at `--buffer_size` records.
- `AnnotationType/RegFeat.pm`, lines 193-195: the same step for regulatory features keys on the feature type and the `dbID`, under the comment that each type has its own `dbID` namespace.
- `Pipeline/DumpVEP/MergeVEP.pm`, lines 94-120: a merged region file is the two sets' transcript lists, each transcript flagged with its set (lines 103 and 110), joined and sorted by start (line 114). Each transcript keeps the `dbID` of its own set.

**What correcting it would change.** Keyed on the set and the `dbID`, the step would still drop the second copy of one transcript and would keep a transcript of the other set. Every record would get the rows of all its transcripts, the same in every run.

**The class it produces.** `merged_same_dbid_collapse`, the entry's name in the index: under `--merged`, the rows of a transcript whose `dbID` a transcript of the other set carries, which vep-rs writes and Ensembl VEP does not. Beside them stands the `intergenic_variant` row Ensembl VEP writes for a record left with no transcript. The SNP/indel comparator sets both aside in a comparison that names the converted cache with `--json-cache`, as [section 1.7](#17-how-the-comparator-treats-this-class) says. The committed corpora hold no such row, and no other comparison holds a rule for the class.

**How often the shape occurs.** The shape is a transcript of a merged cache whose `dbID` a transcript of the other set carries on the same sequence, where the run reads both. Without `--all_refseq`, both engines read a RefSeq-set transcript only where its identifier is a RefSeq accession, such as XM_015636078.1, or on `MT` a mitochondrial gene identifier. Tier 3: no committed cache carries the shape under the flags of a committed run. [How the entries are ordered](../intended-divergences.md#differences-outside-the-classes) defines the tiers and a flag set, one committed run of a corpus. Eight committed flag sets, of four corpora, pass `--merged`, and none of them passes `--all_refseq`.

The two committed human merged caches, GRCh38 and GRCh37 on chromosome 21, hold no such transcript, and neither does any other sequence of Ensembl's release 116 human merged caches. In each whole cache the two sets number their transcripts in ranges that do not meet. In the committed cattle and mouse merged caches, seven `dbID`s and one name a transcript of each set on one sequence. In each of those eight pairs, none of them on `MT`, the RefSeq-set transcript's identifier is no RefSeq accession, so a run reads that transcript under `--all_refseq` alone.

Ensembl's release 116 merged cache of the great tit (*Parus major*, `Parus_major1.1`) carries the shape. A region file that both sets' caches held flags each transcript with its set. Among those transcripts the Ensembl set's `dbID`s run from 38 to 35,314 and the RefSeq set's from 1 to 46,874. 34,586 `dbID`s name a transcript of each set. For 1,398 of them both transcripts lie on one sequence and the RefSeq-set transcript's identifier is a RefSeq accession. Any such pair collides in a buffer that loads both transcripts. For 54 pairs one region file holds both transcripts. For 7 of those a record can lie within 5,000 bases of both, and such a record loses a row annotated alone.

## 1. A great tit record that loses the row of one transcript

**In short.** `8 28659375 missense_forward_1 A G` lies in a coding exon of ENSPMJT00000005922 (Ensembl set) and of XM_015636078.1 (RefSeq set), which both carry `dbID` 5921 in one region file. Under `--merged`, Ensembl VEP 116.3 writes 19 rows for the record and none of ENSPMJT00000005922. vep-rs writes 20, the row of that transcript among them. On the Ensembl-set cache of the same release Ensembl VEP writes a row of that transcript for the record.

### 1.1 What a reader needs to know first

`--merged` reads a merged cache: the transcripts of Ensembl's own annotation, the Ensembl set, and of the RefSeq annotation, the RefSeq set, in one cache. Each row names its set in `SOURCE` ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). Ensembl builds the merged cache from the two sets' caches, region file by region file.

A `dbID` is the number a transcript has in the database its set was dumped from. The Ensembl set comes from one database and the RefSeq set from another, so the two sets of a merged cache are numbered apart.

Ensembl VEP's cache holds a chromosome's transcripts in region files of 1,000,000 bases, and the conversion vep-rs reads keeps them file for file. For each buffer Ensembl VEP reads the region files its records reach, counting 5,000 bases either side of a record.

### 1.2 The record

```
8	28659375	missense_forward_1	A	G	.	.	.
```

Both engines run `--merged` on `Parus_major1.1`: Ensembl VEP 116.3 on Ensembl's release 116 merged cache of the great tit, vep-rs on the conversion of that cache. The region file `8/28000001-29000000.gz` holds 114 transcripts, listed by start. Two of them carry `dbID` 5921:

| Place in the file | Transcript | Set | Span |
| --- | --- | --- | --- |
| 52 | XM_015636078.1 | RefSeq | 8:28,338,627-28,682,436 |
| 63 | ENSPMJT00000005922 | Ensembl | 8:28,444,521-28,682,465 |

### 1.3 What each engine writes

In the default format Ensembl VEP 116.3 writes 19 rows for the record: seven of gene ENSPMJG00000003614 with `SOURCE=Ensembl` and twelve of gene 107208101 with `SOURCE=RefSeq`. Among them is the row of XM_015636078.1:

```
missense_forward_1	8:28659375	G	107208101	XM_015636078.1	Transcript	missense_variant	1726	1498	500	T/A	Aca/Gca	-	IMPACT=MODERATE;STRAND=1;SOURCE=RefSeq
```

vep-rs writes the same 19 rows and one more, the row of ENSPMJT00000005922:

```
missense_forward_1	8:28659375	G	ENSPMJG00000003614	ENSPMJT00000005922	Transcript	missense_variant	1897	1432	478	T/A	Aca/Gca	-	IMPACT=MODERATE;STRAND=1;SOURCE=Ensembl
```

### 1.4 Why Ensembl VEP writes that

The record alone is a buffer that reads one region file, and both transcripts lie within its reach. The file lists its transcripts by start, and XM_015636078.1 starts 105,894 bases before ENSPMJT00000005922. `merge_features` meets XM_015636078.1 first and records `dbID` 5921. When it reaches ENSPMJT00000005922 the number is taken, and the transcript is skipped before any record is annotated against it.

A copy of the region file in which ENSPMJT00000005922 carries a `dbID` no other transcript has shows the cause. With that one value changed, Ensembl VEP 116.3 writes 20 rows for the record, byte for byte the rows vep-rs writes on the unchanged cache.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `AnnotationType/Transcript.pm`, lines 283-287: a transcript whose `dbID` was seen is skipped.
- `Pipeline/DumpVEP/MergeVEP.pm`, line 114: a merged region file lists the two sets' transcripts by start.

### 1.5 Why vep-rs's row is the correct one

ENSPMJT00000005922 is a transcript of the cache, on the forward strand, with ten exons. Position 28659375 is the ninth base of its ninth exon, and the eight exons before it hold 1,888 bases: cDNA position 1897. The coding sequence spans cDNA positions 466 to 1998, so the base is coding position 1432, the first base of codon 478. The record's reference base is the transcript's base there, and the changed codon gives another amino acid: `missense_variant`.

Ensembl VEP writes the same row where the RefSeq set is absent. On Ensembl's release 116 Ensembl-set cache of the great tit, without `--merged`, its row of ENSPMJT00000005922 for the record equals vep-rs's in every column but `SOURCE`, which only a merged cache's rows carry.

### 1.6 Try it yourself

Ensembl serves the cache as `parus_major_merged_vep_116_Parus_major1.1.tar.gz` under `release-116/variation/indexed_vep_cache/` of its FTP site. [Other species](../cache-species.md) gives the conversion command. Convert the directory under its own name, `parus_major_merged`: the conversion declares the cache merged from the ending of that name ([Cache conversion](../cache-conversion.md)). `record.vcf` holds the record of 1.2 under the two header lines of a VCF, `##fileformat=VCFv4.2` and the `#CHROM` column line. vep-rs's `vep`:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache /path/to/json_cache \
  --species parus_major --assembly Parus_major1.1 --merged
```

It writes 20 rows. Ensembl VEP 116.3's `vep` takes the same input and `--merged` with its own cache options:

```
vep -i record.vcf -o out.txt --offline --cache --dir_cache /path/to/cache \
  --species parus_major --assembly Parus_major1.1 --cache_version 116 --merged
```

It writes 19.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/merged_same_dbid_collapse.toml` holds the three records of [section 2](#2-a-pair-in-two-region-files-and-a-record-left-with-no-transcript), the flags and both engines' rows; the great tit record needs Ensembl's cache and is run as 1.6 says.

The class is `merged_same_dbid_collapse` in the 116.3 registry of the SNP/indel comparator, `tools/vep-outputs-compare`. The rows alone cannot tell a row of the class from a row vep-rs writes too many, so the comparator applies the class only in a comparison that names the converted cache with `--json-cache`. From Ensembl VEP's own output, the input VCF and that cache it replays which transcripts each of Ensembl VEP's buffers dropped.

It then sets aside vep-rs's rows of each dropped transcript at the records in reach of it, and Ensembl VEP's `intergenic_variant` row at a record left with no transcript. On the record of this section it sets aside the one row. A comparison that names no cache reports each such row as a row one engine alone writes, and so does one whose reference run the replay does not follow. [Class 18 of the registry's taxonomy](../concordance-provenance/116.2-divergence-classes.md#18-the-classes-that-read-the-converted-cache) gives the conditions.

## 2. A pair in two region files, and a record left with no transcript

**In short.** A synthetic merged cache holds two pairs. XR_000001.1 (RefSeq set) and ENSEXT00000000001 (Ensembl set) overlap in one region file under `dbID` 7. ENSEXT00000000002 (Ensembl set) and XR_000002.1 (RefSeq set) lie 2.7 million bases apart, in two region files, under `dbID` 9. Over three records in one buffer Ensembl VEP 116.3 writes three rows, one of them `intergenic_variant`; vep-rs writes four transcript rows. Under `--buffer_size 1` Ensembl VEP writes the row of XR_000002.1 in place of the `intergenic_variant` row.

### 2.1 The cache

`scripts/divergences/reproductions/inputs/merged_same_dbid_collapse/make_cache.pl` builds the cache with the Perl of an Ensembl VEP installation. It holds one sequence, `contig`, and four single-exon non-coding transcripts on the forward strand:

| Transcript | Set | `dbID` | Span |
| --- | --- | --- | --- |
| XR_000001.1 | RefSeq | 7 | contig:100,001-102,000 |
| ENSEXT00000000001 | Ensembl | 7 | contig:101,001-103,000 |
| ENSEXT00000000002 | Ensembl | 9 | contig:500,001-502,000 |
| XR_000002.1 | RefSeq | 9 | contig:3,200,001-3,202,000 |

The first three stand in the region file `1-1000000`, the fourth in `3000001-4000000`. vep-rs reads the committed conversion of the cache, the `json_cache` directory beside the script.

### 2.2 The records

```
contig	101500	shared_locus	A	G	.	.	.
contig	501000	first_of_far_pair	A	G	.	.	.
contig	3201000	second_of_far_pair	A	G	.	.	.
```

The first record lies on both transcripts of `dbID` 7, the second on ENSEXT00000000002 and the third on XR_000002.1. Both engines run `--merged`.

### 2.3 What each engine writes

Ensembl VEP 116.3, the three records in one buffer:

```
shared_locus	contig:101500	G	900001	XR_000001.1	Transcript	non_coding_transcript_exon_variant	1500	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq
first_of_far_pair	contig:501000	G	ENSEXG00000000002	ENSEXT00000000002	Transcript	non_coding_transcript_exon_variant	1000	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl
second_of_far_pair	contig:3201000	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

vep-rs:

```
shared_locus	contig:101500	G	ENSEXG00000000001	ENSEXT00000000001	Transcript	non_coding_transcript_exon_variant	500	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl
shared_locus	contig:101500	G	900001	XR_000001.1	Transcript	non_coding_transcript_exon_variant	1500	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq
first_of_far_pair	contig:501000	G	ENSEXG00000000002	ENSEXT00000000002	Transcript	non_coding_transcript_exon_variant	1000	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl
second_of_far_pair	contig:3201000	G	900002	XR_000002.1	Transcript	non_coding_transcript_exon_variant	1000	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq
```

Under `--buffer_size 1` Ensembl VEP writes the row of XR_000002.1 for the third record, as vep-rs does, and still no row of ENSEXT00000000001. vep-rs writes the same four rows at every buffer size. On the same cache built with four different `dbID`s, Ensembl VEP writes vep-rs's four rows.

### 2.4 Why Ensembl VEP writes that

The three records are one buffer, and its region files are read in the order of the records: `1-1000000`, then `3000001-4000000`. In the first file XR_000001.1 takes `dbID` 7 before ENSEXT00000000001, which is skipped, and ENSEXT00000000002 takes `dbID` 9. XR_000002.1, read from the second file, is skipped under `dbID` 9. The third record lies on no other transcript and is written as `intergenic_variant`.

With one record a buffer, the third record's buffer reads the second file alone and XR_000002.1 is kept. A pair that one record reaches, the first pair here and the pair of section 1, loses a row whenever that record is annotated: with the record alone, the row of the transcript that starts later. Any other pair loses a row only in a buffer that loads both transcripts. In each case the transcript the buffer reads first is kept.

Why vep-rs's rows are the correct ones is as in 1.5: each record lies in the one exon of its transcript, at the cDNA position the row states. This section does not repeat it, so its next part is 2.6.

### 2.6 Try it yourself

`records.vcf` holds the three records of 2.2 under the two header lines of a VCF, and `cache` is the directory the first command makes:

```
perl scripts/divergences/reproductions/inputs/merged_same_dbid_collapse/make_cache.pl cache
vep -i records.vcf -o out.txt --offline --cache --dir_cache cache \
  --species two_set_example --assembly EXAMPLE1 --cache_version 116 --merged
```

The first command runs under the Perl of an Ensembl VEP 116.3 installation, and the second is Ensembl VEP's `vep`: it writes the three rows of 2.3, and with `--buffer_size 1` the row of XR_000002.1 in place of the `intergenic_variant` row. vep-rs reads the committed conversion:

```
vep -i records.vcf -o out.txt --offline \
  --json_cache scripts/divergences/reproductions/inputs/merged_same_dbid_collapse/json_cache \
  --species two_set_example --assembly EXAMPLE1 --merged
```

It writes the four rows of 2.3.
