# The order of the `DOMAINS` list on a transcript stored in two region files

**The defect in one sentence.** Ensembl VEP's cache stores a transcript that reaches across a region boundary once in each region file, and two copies can hold the transcript's protein features in different orders. Ensembl VEP 116.2 writes `DOMAINS` in the order of whichever copy its buffer reads first, so one record gets two orders, decided by the other records of the run and by `--fork`. On a cache converted from Ensembl VEP's, vep-rs writes the order of the lowest region file's copy for every record, at every thread count.

**Why this is a defect.** The record, the flags and the cache fix the entries of a `DOMAINS` value. Their order in Ensembl VEP's output comes from the region file the run read the transcript from. One record of a ClinVar file is written with one order when the file is annotated without `--fork` or under `--fork 4`, and with the other under `--fork 8` or `--fork 16`. A record annotated alone and the same record after an earlier one differ the same way. Two runs of one input therefore differ byte for byte, and `--fork`, an option for speed, changes what is written.

**What vep-rs does instead.** A cache converted from Ensembl VEP's holds both copies. vep-rs reads a chromosome's region files in order of their start position and keeps the first copy of each transcript. It writes the order of the lowest region file's copy for every record, whatever else the input holds and at every thread count. A native cache stores each transcript once, its protein features in order of their identifiers ([Cache setup](../cache-setup.md)).

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`; `ensembl-variation` and `ensembl` at release/116):

- `OutputFactory.pm`, lines 1494-1512: the value is one label per feature of `get_overlapping_ProteinFeatures`, in that list's order. `BaseTranscriptVariation.pm` of `ensembl-variation`, lines 639-643, filters the transcript's stored `protein_features` array by overlap and keeps its order. Nothing sorts either list.
- `AnnotationSource/Database/Transcript.pm`, lines 459-480: the cache dumper stores the array `get_all_ProteinFeatures` returns. `Translation.pm` of `ensembl`, lines 772-819, groups a translation's features by analysis and joins the groups in the order of a hash's keys (lines 817-819), so two fetches of one transcript can return two orders.
- `AnnotationSource.pm`, lines 109-143 and 168-262: a buffer's region files are listed in the order of its records, each record reaching 5,000 bases either side. The features of files held in memory come first, then those read from disk (`AnnotationSource/Cache/BaseSerialized.pm`, lines 243-270).
- `AnnotationType/Transcript.pm`, lines 266-299, `merge_features`: the first record of each transcript identifier is kept and every later one dropped.
- `Runner.pm`, lines 470-528 and 643-692: under `--fork` the buffer is cut into slices of at most `int(buffer_size / (2 * fork))` records (line 486), and each worker lists the region files of its slice alone (lines 667-670).

**What correcting it would change.** The dumper could store a translation's features in a defined order, or the output could sort the labels. Either way every copy of a transcript would give one list, and a record's `DOMAINS` would be the same on every run and at every `--fork`. Only the order would move: no transcript of the release-116 caches has copies that hold different features.

**The class it produces.** `domains_order_region_file_copy`: the `DOMAINS` value (JSON `domains`) of a coding row on a transcript the cache stores in two or more region files with its protein features in different orders. Both engines list the same entries. The format comparator and the golden harness compare a `DOMAINS` value as a multiset of its items in a run that names their built-in class `order_only`. Such a run reports no row of this shape, and nothing is counted under the entry's own name. A run that does not name `order_only` reports the row as a difference in `DOMAINS`.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the divergence needs `--domains`, or `--everything`, which switches it on, and a transcript stored twice in two orders. The committed caches hold one such transcript, in one corpus, and 7 committed flag sets pass `--domains` or `--everything`.

The release-116 GRCh38 cache stores 57,293 of its 669,092 transcripts in more than one region file. 24 of them have copies in different orders, and 23 of the 24 show it at some residue. In the release-116 GRCh37 cache, read through its conversion, the figures are 14,677 of 238,823, then 4 and 4. On the chromosome 2 records of the ClinVar GRCh38 file, 40 rows of 10,835,497 differ under `--fork 8` or `--fork 16`, and none without `--fork` or under `--fork 4`.

## 1. One ClinVar record written in two orders

**In short.** `2 1209204 750974 G A`, a ClinVar record, is a synonymous change at codon 231 of ENST00000308624 (SNTG2, forward strand). Under `--domains` both engines list the same three entries on its row. Ensembl VEP 116.2 writes the PANTHER entry first when the record is annotated alone, and the Superfamily entry first when an earlier record of the buffer lies in the region file below. Over the chromosome's ClinVar records its output under `--fork 8` differs from its own output without `--fork` on 40 rows. vep-rs writes the Superfamily entry first on every run.

### 1.1 What a reader needs to know first

`--domains` adds the `DOMAINS` column: on a coding row, the protein features whose residues the variant overlaps, each as `label:identifier` ([Output control](../cli/output-control.md)).

Ensembl VEP's cache holds a chromosome's transcripts in region files of 1,000,000 bases. A transcript whose gene reaches across a boundary is stored whole in each file the gene overlaps, and every copy carries its own array of protein features. Ensembl VEP reads its input in buffers of `--buffer_size` records, 5,000 by default. For each buffer it fetches the transcripts of every region file a record reaches, counting 5,000 bases either side of the record.

ENST00000308624 spans 2:950,849-1,367,613 and is stored in `2/1-1000000.gz` and `2/1000001-2000000.gz`. Each copy holds the same 19 features over the same residues. The first lists the three Superfamily features before the PANTHER family, the second lists PANTHER before them, and both list the AFDB-ENSP mapping last. Residue 231 lies inside one Superfamily feature (147-277), the PANTHER family (23-523) and the AFDB-ENSP mapping (1-539).

### 1.2 The record

```
2	1209204	750974	G	A	.	.	.
```

A second input puts another ClinVar record, which lies in the first region file, before it:

```
2	951066	3446878	C	T	.	.	.
2	1209204	750974	G	A	.	.	.
```

The third input is the chromosome 2 records of the ClinVar GRCh38 file, 384,041 records sorted by position. Both engines run `--domains` on GRCh38, Ensembl VEP on its release-116 cache and vep-rs on the conversion of that cache.

### 1.3 What each engine writes

The record alone, in the default format. Ensembl VEP 116.2:

```
750974	2:1209204	A	ENSG00000172554	ENST00000308624	Transcript	synonymous_variant	841	693	231	R	agG/agA	-	IMPACT=LOW;STRAND=1;DOMAINS=PANTHER:PTHR10554,Superfamily:SSF50729,AFDB-ENSP_mappings:AF-Q9NY99-F1
```

vep-rs:

```
750974	2:1209204	A	ENSG00000172554	ENST00000308624	Transcript	synonymous_variant	841	693	231	R	agG/agA	-	IMPACT=LOW;STRAND=1;DOMAINS=Superfamily:SSF50729,PANTHER:PTHR10554,AFDB-ENSP_mappings:AF-Q9NY99-F1
```

The record's 13 other rows are equal, and the tab, VCF and JSON outputs carry the same two orders. On the two-record input Ensembl VEP writes the row of `2:1209204` with the Superfamily entry first, as vep-rs does, and no row differs. With `--buffer_size 1` on that input it writes the PANTHER entry first again.

Over the whole chromosome, `--domains` in the tab format gives 10,835,497 rows, 3,865,102 of them with a `DOMAINS` value:

| Ensembl VEP run | Rows whose `DOMAINS` differs from vep-rs's |
| --- | --- |
| no `--fork` | 0 |
| `--fork 4` | 0 |
| `--fork 8` | 40 |
| `--fork 16` | 40, the same rows |

The 40 rows are rows of ENST00000308624 on the records from `2:1209204` to `2:1267518`, each with the same entries in the other order. On the transcript's rows past them, the entries a row lists stand in one order in both copies, so the copy read does not show. Ensembl VEP's output without `--fork` and its output under `--fork 8` differ from each other on those 40 rows and on no other.

### 1.4 Why Ensembl VEP writes that

For each buffer Ensembl VEP lists the region files its records reach, in the order of the records, and gathers the transcripts of those files in that order. A transcript met twice is kept from the first file of the list that holds it. The record `2:1209204` alone reaches `1000001-2000000.gz` only, whose copy lists PANTHER first. After `2:951066` the list opens with `1-1000000.gz`, whose copy lists Superfamily first. Under `--buffer_size 1` each record is a buffer of its own.

Under `--fork` the buffer is cut into slices, and each worker lists the region files of its own slice. A slice holds at most `int(buffer_size / (2 * fork))` records: 625 at 4 forks, 312 at 8 and 156 at 16. The last record of the file that reaches the first region file is record 261: none after it lies at or below 1005000, the file's end plus the flank. The record `2:1209204` is record 313.

At 8 forks the second slice starts at record 313 and at 16 forks the third does. Neither holds a record that reaches the first region file, so its worker reads the second copy. Without `--fork` the buffer holds records 1 to 5,000, and at 4 forks the first slice holds records 1 to 625. Both hold record 261, so the first copy is read.

Source (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`):

- `AnnotationSource.pm`, lines 109-143: the region list, the features in its order and the call of `merge_features`.
- `AnnotationType/Transcript.pm`, lines 283-286: a transcript identifier already seen is skipped.
- `Runner.pm`, lines 486 and 513: the slice size and the cut.

### 1.5 Why vep-rs's row is the correct one

The row's three entries are the same in both engines and in both stored copies: the features that cover residue 231. Neither order can be read from the record. Each is the stored order of one copy, and the cache gives no rule for choosing between two copies of one transcript. Ensembl VEP's choice moves with the other records of the buffer and with `--fork`, so its own runs of one file disagree.

vep-rs's choice is fixed by the cache alone: the lowest region file's copy, for every record and at every thread count. It is also the order Ensembl VEP writes for a file sorted by position whenever an earlier record of the buffer or slice reaches the lower region file. Without `--fork` that is every row of this chromosome.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline \
  --json_cache /path/to/json_cache/homo_sapiens/116_GRCh38 \
  --species homo_sapiens --assembly GRCh38 --domains
```

vep-rs writes the Superfamily entry first for the one-record and the two-record input. Ensembl VEP 116.2 takes `--domains` with the cache options of [Reproducing these examples](reproducing.md). It writes PANTHER first for the one-record input, Superfamily first for the two-record input, and PANTHER first for the two-record input under `--buffer_size 1`.

For the fork counts, run Ensembl VEP over the chromosome 2 records of the ClinVar GRCh38 file with `--tab --domains`, once without `--fork` and once with `--fork 8`. The two outputs differ on 40 rows.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/domains_order_region_file_copy.toml` runs the record of section 2 on both engines and records the `DOMAINS` item of its row. vep-rs reads the committed cache of 2.2 there, and Ensembl VEP its release-116 cache.

The class covers the order of a `DOMAINS` list on a transcript whose stored copies differ in order. The golden harness and the format comparator compare the value as a multiset of its items where a run names their built-in class `order_only`. A list whose items differ is reported like any other field.

## 2. A record on a committed cache

**In short.** `21 26170588 740728 G A`, a ClinVar record, is a synonymous change at codon 11 of ENST00000348990 (APP, reverse strand). The transcript is stored in two region files of chromosome 21, and the two copies place one of its 49 protein features differently. Annotated alone, the record gets that feature before the AFDB-ENSP mapping from Ensembl VEP 116.2 and after it from vep-rs. The cache of the `GRCh38-everything` golden corpus holds both copies, so vep-rs's side of the reproduction runs on committed files.

### 2.1 The two copies

ENST00000348990 spans 21:25,880,550-26,170,830 and is stored in `21/25000001-26000000.gz` and `21/26000001-27000000.gz`. Each copy holds the same 49 features over the same residues. The first lists the `SignalP_GRAM_POSITIVE` feature (residues 1-17) last, in place 49, after the AFDB-ENSP mapping. The second lists it in place 29 and the AFDB-ENSP mapping last. Residue 11 lies inside seven features, those two among them.

### 2.2 The record

```
21	26170588	740728	G	A	.	.	.
```

The record lies 170,588 bases above the boundary at 26000000, so a buffer that holds it alone reaches the second region file only. A second input puts a ClinVar record of the first region file before it:

```
21	25975051	1600556	C	T	.	.	.
21	26170588	740728	G	A	.	.	.
```

Both engines run `--domains` on GRCh38. Ensembl VEP reads its release-116 cache, and vep-rs reads `tests/golden/116/GRCh38-everything/json_cache`, a slice of the conversion of that cache.

### 2.3 What each engine writes

The record alone, the row of ENST00000348990 in the default format. Ensembl VEP 116.2:

```
740728	21:26170588	A	ENSG00000142192	ENST00000348990	Transcript	synonymous_variant	243	33	11	A	gcC/gcT	-	IMPACT=LOW;STRAND=-1;DOMAINS=Low_complexity_(Seg):seg,Cleavage_site_(Signalp):SignalP-noTM,PANTHER:PTHR23103,Phobius:SIGNAL_PEPTIDE,Phobius:SIGNAL_PEPTIDE_H_REGION,SignalP_GRAM_POSITIVE:SignalP-TM,AFDB-ENSP_mappings:AF-P05067-F1
```

vep-rs:

```
740728	21:26170588	A	ENSG00000142192	ENST00000348990	Transcript	synonymous_variant	243	33	11	A	gcC/gcT	-	IMPACT=LOW;STRAND=-1;DOMAINS=Low_complexity_(Seg):seg,Cleavage_site_(Signalp):SignalP-noTM,PANTHER:PTHR23103,Phobius:SIGNAL_PEPTIDE,Phobius:SIGNAL_PEPTIDE_H_REGION,AFDB-ENSP_mappings:AF-P05067-F1,SignalP_GRAM_POSITIVE:SignalP-TM
```

The two rows hold the same seven entries, and the last two change places. On the two-record input Ensembl VEP writes the row with the AFDB-ENSP mapping before the `SignalP_GRAM_POSITIVE` entry, as vep-rs does. With `--buffer_size 1` on that input it writes the order of the record alone. vep-rs writes one order for all three inputs.

Why Ensembl VEP writes that and why vep-rs's row is the correct one are as in 1.4 and 1.5; this section does not repeat them, so its next part is 2.6.

### 2.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38-everything/json_cache \
  --species homo_sapiens --assembly GRCh38 --domains
```

vep-rs writes the AFDB-ENSP mapping before the `SignalP_GRAM_POSITIVE` entry on the row. Ensembl VEP 116.2 takes `--domains` with the cache options of [Reproducing these examples](reproducing.md) and writes the two entries the other way round for the record alone. The reproduction of section 1.7 is this record: `scripts/divergences/reevaluate.py` runs it on both engines and holds each value to the recorded one.
