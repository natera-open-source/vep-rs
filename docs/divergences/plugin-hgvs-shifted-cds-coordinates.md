# The Downstream and NMD plugins at the HGVS-shifted coordinates

Under `--hgvs`, Ensembl VEP's Downstream and NMD plugins answer for an insertion's or deletion's 3'-shifted CDS coordinates, vep-rs's for the record's own, as without the flag.

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP's protein notation leaves an insertion's or deletion's CDS coordinates at the 3'-shifted position on the row, and the Downstream and NMD plugins read them there. `DownstreamProtein`, `ProteinLengthChange` and `NMD` can then differ from the values the same plugins write for the same record without the flag.

**Why this is a defect.** `--hgvs` asks for notation columns and nothing else. The row's `CDS_position`, `Protein_position` and `Codons` columns print the record's own coordinates while the plugin answers for another span, and Ensembl VEP itself writes the own-span value once the flag is dropped. For Downstream the two spans spell one changed sequence when the shift stays inside an exon, so the two `DownstreamProtein` values read one changed sequence from two starting points. Across an intron, or onto a seam of an edited RefSeq model, the shifted reading builds another sequence, and on a seam `ProteinLengthChange` differs too. For NMD the rule on the first 100 coding bases is answered for a span the row does not print.

**What vep-rs does instead.** vep-rs's Downstream and NMD plugins read the record's own CDS coordinates, the pair the `CDS_position` column prints, and the allele as written, with and without `--hgvs`. A row's plugin values are the same in both runs. Ensembl VEP's plugins read that same pair when the flag is absent.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `ensembl-variation` and `VEP_plugins` at release/116):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1737, prints `CDS_position` from the pair the row's object holds: the record's own, except on a later allele of a multi-allelic record ([section 9](#9-a-later-allele-of-a-multi-allelic-record)). Line 1760, later in the same row's build, calls `hgvs_protein` under `--hgvs`.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 844-858, `codon`: when the shift is not zero, the row's kept `cds_start` and `cds_end` are deleted and read again at the shifted position. Lines 2406-2410, `_get_alternate_cds`, read them there once more, and that pair stays on the row.
- The same file, lines 1330-1345, `shift_feature_seqs`, called at line 1778, rotates the allele's bases for the shifted position. Lines 1730-1738 return from `hgvs_protein` before any of this when the shifted span has no translation coordinates, and the own pair stays.
- `OutputFactory.pm`, line 398, runs the plugins once the row is built. Nothing clears the pair in between.
- `Downstream.pm`, line 108, reads `($tv->cds_start, $tv->cds_end)`, the shifted pair, and line 110 splices the rotated `feature_seq`. `NMD.pm`, lines 88, 98 and 103, read the shifted pair, while line 137 reads the record's own genomic end.

**What correcting it would change.** The row's CDS pair and the allele's bases would be put back, once `hgvs_protein` returns, to the values the row prints. `OutputFactory.pm` already clears shifted coordinates before it builds a buffer's rows (lines 297-298), through `clear_shifting_variables` (`TranscriptVariationAllele.pm`, lines 494-528), whose description names the plugins among its readers. Each plugin would then give one value with and without `--hgvs`, the value of the row's own `CDS_position`.

**The class it produces.** `plugin_hgvs_shifted_cds_coordinates`: a `DownstreamProtein`, a `ProteinLengthChange` or an `NMD` (JSON `downstreamprotein`, `proteinlengthchange`, `nmd`) that differs on a row whose Ensembl VEP entry carries a non-zero `HGVS_OFFSET`. It is also either Downstream field that differs on a later allele's row without an `HGVS_OFFSET` whose earlier allele's row carries one ([section 10](#10-a-later-allele-whose-own-notation-moves-nothing)).

The class is built into the format comparator, which decides it from the candidate cache named by `--json-cache`: the transcript's strand, exons, coding bounds, mapper pairs, sequences and codon table. It maps the record's span, and that span moved by the offset, onto the coding sequence and recomputes each plugin's rule on both pairs. It accepts a row only where Ensembl VEP's value is its rule on the shifted pair, vep-rs's value is vep-rs's rule on the record's own pair, and Ensembl VEP's rule gives another value on the two pairs. On the row of section 10 the offset is the earlier allele's and the allele is read as written.

The same values on any other row without an `HGVS_OFFSET`, on a row whose shift does not reach the plugin, and on a row whose shifted pair gives the own pair's value stay differences. So do an `NMD` on any row without an `HGVS_OFFSET`, a wrong value on either side and a Downstream value one side lacks. A `DownstreamProtein` that Ensembl VEP's rule gives alike on both pairs belongs to [the start one codon late](downstream-start-past-own-codon.md) or to [the read past the coding end](downstream-utr-frame-past-coding-end.md), so the three classes share no row.

On some rows the start rule changes the value as well: vep-rs's value is not the one Ensembl VEP writes without the flag either, because the record's own span starts on a codon's third base. Such a row is this class's only where `downstream_start_past_own_codon` is named beside it. With either class named alone it stays a difference, since neither explains all three values.

**How often the shape occurs.** Tier 2, the flag and a plugin being needed. Under `--hgvs --plugin Downstream`, vep-rs writes the plugin's fields on 557,264 rows of the genome-wide ClinVar GRCh37 file. The `DownstreamProtein` of 122,020 of them, from 35,972 records, falls under the class, and no `ProteinLengthChange` does. On the GRCh38 file it is 594,975 of 2,495,728 rows, from 37,669 records, and one `ProteinLengthChange`. The start rule changes the value too on 11,573 of the GRCh37 rows and 56,235 of the GRCh38 rows, which need both classes named. With the NMD plugin loaded under `--hgvs`, 518 rows of the GRCh37 file differ in `NMD`: Ensembl VEP alone writes it on 173 and vep-rs alone on 345.

On the five chromosome 21 population cuts, run under `--hgvs --plugin Downstream --plugin NMD`, 1,954 of 15,646,481 records have a row of the class, 124.9 per million: `DownstreamProtein` on 1,948 of them, `NMD` on 9 and `ProteinLengthChange` on 1. Four committed flag sets carry rows of it, all in `DownstreamProtein`: 29 in each of the two runs of the `GRCh38-plugin-computed` corpus, 31 in the `--hgvs` run of `GRCh37-plugin-downstream-own-span` and ten in that of `GRCh37-plugin-downstream-coding-end`.

## 1. A four-base deletion in a repeat of NFIX

**In short.** `19 13136318 2948877 ACAAT A`, a ClinVar record, deletes coding bases 371 to 374 of ENST00000590027 (NFIX, forward strand). The same deletion can be written three bases on, and HGVS notation takes that position: `c.374_377del`, `HGVS_OFFSET=3`. vep-rs writes `DownstreamProtein=TKNW`, the translation from codon 124, where the record starts. Ensembl VEP 116.2 writes `KNW` under `--hgvs`, the translation from codon 125, where the shifted span starts, and `TKNW` without the flag.

### 1.1 What a reader needs to know first

An insertion or deletion inside a repeat can be written at several positions that give one changed sequence. HGVS notation names the most 3' of them on the transcript. Under `--hgvs` Ensembl VEP writes `HGVSc` and `HGVSp` there and the distance moved in `HGVS_OFFSET`. Every other column of the row describes the record as written: the consequence, `cDNA_position`, `CDS_position`, `Protein_position`, `Amino_acids` and `Codons`.

The Downstream plugin writes two fields on a frameshift row. `DownstreamProtein` is the translation of the changed transcript from the variant's codon on, 3' UTR included, up to the first stop; the plugin's description counts in "any amino acids overlapped by the variant itself". `ProteinLengthChange` is the length of the whole changed protein less the reference protein's.

The coding sequence of ENST00000590027 in the cache is 386 bases, 128 codons and two more bases: its 3' end is not annotated (`FLAGS=cds_end_NF`), and the cache spells no 3' UTR. Codons 124 to 128 read `ACA ATC AAA GAA CTG`, threonine, isoleucine, lysine, glutamate and leucine, and the last two bases are `GA`.

### 1.2 The record

```
19	13136318	2948877	ACAAT	A	.	.	.
```

The record deletes `CAAT`, 19:13136319-13136322: the last two bases of codon 124 and the first two of codon 125. The three bases after it are `CAA`, so deleting `TCAA`, coding bases 374 to 377, leaves the same sequence. It is a record of the `GRCh37-plugin-downstream-own-span` corpus, whose flags both engines run: `--hgvs --plugin Downstream`.

### 1.3 What each engine writes

Ensembl VEP 116.2, default format:

```
2948877	19:13136319-13136322	-	ENSG00000008441	ENST00000590027	Transcript	frameshift_variant	564-567	371-374	124-125	TI/X	aCAATc/ac	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000590027.1:c.374_377del;HGVSp=ENSP00000465616.1:p.Ile125LysfsTer?;HGVS_OFFSET=3;DownstreamProtein=KNW;ProteinLengthChange=-1
```

vep-rs:

```
2948877	19:13136319-13136322	-	ENSG00000008441	ENST00000590027	Transcript	frameshift_variant	564-567	371-374	124-125	TI/X	aCAATc/ac	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000590027.1:c.374_377del;HGVSp=ENSP00000465616.1:p.Ile125LysfsTer?;HGVS_OFFSET=3;DownstreamProtein=TKNW;ProteinLengthChange=-1
```

Every other column is equal. The tab, VCF and JSON outputs carry the same two values, the JSON under `downstreamprotein`. Without `--hgvs`, Ensembl VEP 116.2 writes `DownstreamProtein=TKNW` on the row, as vep-rs does in both runs. The record's rows on eight more transcripts of the gene differ the same way, 37 residues against 38.

### 1.4 Why Ensembl VEP writes that

Line 1737 of `OutputFactory.pm` has printed `371-374` when line 1760 computes the protein notation. `hgvs_protein` moves the deletion three bases, and `codon` and `_get_alternate_cds` leave the pair 374 and 377 on the row. Line 108 of `Downstream.pm` reads that pair, and line 110 removes those four bases. Line 113 computes `374 - (374 % 3)`, which is 372, the offset of codon 125. From codon 124 the changed sequence reads `ACA AAG AAC TGG A`, and the translation starts at `AAG`: lysine, asparagine, tryptophan. Without `--hgvs` line 108 reads 371 and 374, line 113 gives 369, and the translation starts at `ACA`.

### 1.5 Why vep-rs's row follows from the record's own coordinates

The record removes coding bases 371 to 374. Base 370, `A`, stays, and bases 375 and 376, `CA`, complete codon 124 as `ACA` again. The codons after it read `AAG AAC TGG`, and the sequence ends one base later. The row's own columns place the variant from codon 124: `CDS_position` 371-374, `Protein_position` 124-125, `Amino_acids` `TI/X`. Codon 124 is a codon the record overlaps, so the plugin's description counts its residue in, and `TKNW` is the sequence from it.

`KNW` is the same protein read from the codon where `c.374_377del` starts, a span no other column of the row describes. Ensembl VEP writes `TKNW` itself when `--hgvs` is absent.

### 1.6 Try it yourself

From `tests/golden/116/GRCh37-plugin-downstream-own-span/`:

```
vep -i variants.vcf -o out.txt --offline --json_cache json_cache \
  --assembly GRCh37 --fasta reference.fa.gz --hgvs --plugin Downstream
```

vep-rs writes `DownstreamProtein=TKNW` on the row, and the same value when `--hgvs` is dropped. The corpus's `expected/` directory holds Ensembl VEP 116.2's four outputs of the same records under `--hgvs`, with `KNW`. Its `expected_without-hgvs/` directory holds the four outputs of the run without the flag, with `TKNW`.

### 1.7 How the comparator treats this class

The format comparator accepts the row under `plugin_hgvs_shifted_cds_coordinates` when it is given the candidate cache. It maps 19:13136319-13136322 to coding bases 371 to 374, and the span three bases on to 374 to 377. Ensembl VEP's value is its rule's translation on the shifted pair, vep-rs's the translation from codon 124 on the own pair, and Ensembl VEP's rule gives `TKNW`, another value, on the own pair. The corpus manifest documents the `DownstreamProtein` of 31 rows of the `--hgvs` run under the class, and the golden harness reads vep-rs's value from it.

The same row with a wrong value on either side, with no value or with no `HGVS_OFFSET` stays a difference. Without an offset the class takes only a later allele's row whose earlier allele's row carries one ([section 10](#10-a-later-allele-whose-own-notation-moves-nothing)).

## 2. An insertion whose shift crosses an intron

**In short.** ClinVar record 2583096 inserts 133 bases between X:152960073 and X:152960074, 15 bases before the end of an exon of ENST00000253122 (SLC6A8, forward strand). The row prints `CDS_position` 1581-1582. The notation moves the insertion 143 bases along the genome, over the 85-base intron and 43 bases into the next exon: `HGVS_OFFSET=143`. On the coding sequence that is 58 bases, to between bases 1639 and 1640. vep-rs writes a `DownstreamProtein` of 74 residues, from codon 528. Ensembl VEP 116.2 writes 55 residues from codon 547 under `--hgvs`, and vep-rs's 74 without it.

Inside one exon the two pairs spell one changed sequence. Here they do not. Under `--hgvs` the plugin puts the 133 bases, rotated for the genomic shift, after coding base 1639, behind 58 exon bases that the record's own reading leaves after the insertion. The 55 residues are therefore not the last 55 of the 74. Both readings add 133 bases and meet the same stop, so `ProteinLengthChange=-34` on both engines. The record is in the `GRCh37-plugin-downstream-own-span` corpus, with its row on ENST00000430077 and with record 4730422, another insertion shifted over an intron.

## 3. An `NMD` value at the edge of the first 100 coding bases

**In short.** `3 10074546 1322879 CAAAG C`, a ClinVar record, deletes coding bases 96 to 99 of ENST00000287647 (FANCD2, forward strand). The NMD plugin marks a frameshift within the first 100 coding bases as escaping nonsense-mediated decay. vep-rs writes `NMD=NMD_escaping_variant`. Ensembl VEP 116.2 writes no `NMD` under `--hgvs`, where its plugin reads the shifted end, base 102, and writes the field without the flag.

### 3.1 What a reader needs to know first

The NMD plugin writes `NMD=NMD_escaping_variant` on a stop-gained, frameshift or splice-site row whose variant meets one of its positional rules. Two matter here. The variant falls in the first 100 coding bases, which the module tests as a CDS end of at most 101. Or the record's genomic end lies in the transcript's last exon, or within 51 bases of the 3' end of the exon before it. The plugin gives no value when either end of the CDS pair is undefined, as for a deletion that starts in an intron.

### 3.2 The record

```
3	10074546	1322879	CAAAG	C	.	.	.
```

The record deletes `AAAG`, 3:10074547-10074550, in the third of the transcript's 43 exons: the last base of codon 32 and the whole of codon 33. The three bases after it are `AAA`, so the deletion of `GAAA`, coding bases 99 to 102, is the same change: `c.99_102del`, `HGVS_OFFSET=3`. Both engines run `--hgvs --plugin NMD` on GRCh37 with a reference FASTA.

### 3.3 What each engine writes

Ensembl VEP 116.2, default format:

```
1322879	3:10074547-10074550	-	ENSG00000144554	ENST00000287647	Transcript	frameshift_variant	189-192	96-99	32-33	TK/X	acAAAG/ac	-	IMPACT=HIGH;STRAND=1;HGVSc=ENST00000287647.3:c.99_102del;HGVSp=ENSP00000287647.3:p.Lys33AsnfsTer16;HGVS_OFFSET=3
```

vep-rs:

```
1322879	3:10074547-10074550	-	ENSG00000144554	ENST00000287647	Transcript	frameshift_variant	189-192	96-99	32-33	TK/X	acAAAG/ac	-	IMPACT=HIGH;STRAND=1;HGVSc=ENST00000287647.3:c.99_102del;HGVSp=ENSP00000287647.3:p.Lys33AsnfsTer16;HGVS_OFFSET=3;NMD=NMD_escaping_variant
```

Every other column is equal. Run with `--plugin NMD` and no `--hgvs`, Ensembl VEP 116.2 writes `NMD=NMD_escaping_variant` on the row, as vep-rs does in both runs. The record's rows on four more transcripts of the gene differ the same way.

### 3.4 Why Ensembl VEP writes that

Under `--hgvs` the pair left on the row is 99 and 102. Line 88 of `NMD.pm` reads the end, 102, and line 98 tests it against 101, which fails. Line 137 reads the record's own genomic end, 3:10074550, in exon 3, so neither exon rule holds and the plugin writes nothing. Without the flag line 88 reads 99, and the rule holds.

### 3.5 Why vep-rs's row follows from the record's own coordinates

The record's last deleted base is coding base 99, as the row's `CDS_position`, `96-99`, prints. By the plugin's own rule a variant that ends there falls in the first 100 coding bases. Base 102 is the end of the span `c.99_102del` names. The plugin's exon rules do not use that span either: they read the record's own end. Under `--hgvs` the plugin therefore answers one rule for the shifted span and the others for the record's.

### 3.6 Try it yourself

No committed corpus holds the record. With the GRCh37 cache and FASTA of [Reproducing these examples](reproducing.md), and the record in `record.vcf`:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache /path/to/json_cache/homo_sapiens/116_GRCh37 \
  --fasta /path/to/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
  --assembly GRCh37 --hgvs --plugin NMD
```

vep-rs writes `NMD=NMD_escaping_variant` on the row, with and without `--hgvs`. Ensembl VEP 116.2, run on the same record against its own GRCh37 cache, writes the field only without `--hgvs`.

### 3.7 How the comparator treats this class

For `NMD` the comparator accepts a row where one entry carries `NMD_escaping_variant` and the other lacks it. Ensembl VEP's entry must carry it exactly when the plugin's rules hold on the shifted pair, and vep-rs's exactly when they hold on the record's own pair. Here the rule on the CDS end holds for 99 and fails for 102. A field both entries carry or both lack, and such a pair on a row without an `HGVS_OFFSET`, stay differences, a later allele's row whose earlier allele's row carries one included ([section 10](#10-a-later-allele-whose-own-notation-moves-nothing)).

## 4. An `NMD` value the shift adds

**In short.** `17 74729058 1184977 AGATTG A`, a ClinVar record, deletes five bases of ENST00000341249 (METTL23, forward strand). The first is the last intron base before the third exon, so the row prints `CDS_position` `?-88`, and the plugin gives no value for a pair with an undefined end. The notation moves the deletion two bases, onto coding bases 86 to 90: `c.86_90del`, `HGVS_OFFSET=2`. Under `--hgvs` Ensembl VEP 116.2 reads that pair and writes `NMD=NMD_escaping_variant`, the end being at most 101. Without the flag it writes no `NMD`, and vep-rs writes none in either run.

## 5. A deletion on which the flag and the start rule both change the value

**In short.** ClinVar record 1703957 deletes 25 bases, X:152960007-152960031: coding bases 1515 to 1539 of ENST00000253122 (SLC6A8), the last base of codon 505 and the eight codons after it. Both engines write `c.1519_1543del`, `HGVS_OFFSET=4`, `p.Ile507LeufsTer5` and `ProteinLengthChange=-125`. From codon 505 the changed sequence reads `GAC GAC CTT GCC CCT GGA` and a stop. vep-rs writes `DownstreamProtein=DDLAPG`. Ensembl VEP 116.2 writes `LAPG` under `--hgvs` and `DLAPG` without it.

Two mechanisms meet on the row. Base 1515 is a codon's third, so line 113 of `Downstream.pm` starts one codon late on the record's own span: `DLAPG`, from codon 506, where vep-rs starts at codon 505 ([the start one codon late](downstream-start-past-own-codon.md)). Under `--hgvs` the plugin reads the span four bases on, from the first base of codon 507, and writes `LAPG`. The comparator files the `--hgvs` row under this class, since Ensembl VEP's value changes with the flag, and the row of the run without it under `downstream_start_past_own_codon`. It sets the `--hgvs` row aside only where both classes are named, and the corpus's manifest names both on the row's entry. The record is in the `GRCh37-plugin-downstream-own-span` corpus.

On a reverse-strand transcript the shifted reading can also splice bases the variant does not produce: [the rotation of a shifted allele](../plugins/transcript-model.md#4-downstream-the-rotation-of-a-shifted-allele-on-a-reverse-strand-transcript).

## 6. A shift onto a codon's third base: residues past a stop

**In short.** `21 46430158 1338605 CTG C`, a ClinVar record, deletes coding bases 7840 and 7841 of ENST00000359568 (PCNT, forward strand, GRCh38), the first two bases of codon 2614, `TGT`. The base left and the two after it read `TGA`, a stop. Both engines write `p.Cys2614Ter` and `ProteinLengthChange=-723`, and vep-rs writes an empty `DownstreamProtein`. The notation moves the deletion two bases, onto the last base of codon 2614: `c.7842_7843del`, `HGVS_OFFSET=2`. Under `--hgvs` Ensembl VEP 116.2 writes 28 residues.

Line 113 of `Downstream.pm` computes `7842 - (7842 % 3)`, which is 7842, the offset of codon 2615. From codon 2613 the changed sequence reads `CTC TGA GAG CAG GCA`, and the plugin translates from `GAG`: the bases after the stop the variant makes. On the record's own span the same line gives the offset of codon 2614, and the value is empty, as vep-rs writes it. The shift alone puts the row in the shape of [the start one codon late](downstream-start-past-own-codon.md). The record is in the `GRCh38-plugin-computed` corpus.

## 7. A shift onto the last coding base: the 3' UTR read out of frame

**In short.** `19 50813042 932016 CAGTG C`, a ClinVar record, deletes coding bases 6107 to 6110 of ENST00000262269 (MYH14, forward strand), from the second base of the last codon, `CAG`, into the stop codon `TGA`. The notation moves it one base, `c.6108_6111del`, `HGVS_OFFSET=1`: a start on the third base of codon 2036 and an end on the last coding base. Both engines write `p.Gln2036HisfsTer32` and `ProteinLengthChange=30`. vep-rs writes the 31 residues before that stop. Ensembl VEP 116.2 writes 58 residues under `--hgvs`, and vep-rs's 31 without the flag.

The bases left are `CA`, and in frame the changed transcript reads on into the 3' UTR as `CAC CCT ACC CTG`, histidine first, as the `HGVSp` says. On the shifted span line 113 of `Downstream.pm` gives an offset one past the changed coding sequence, so the plugin translates the UTR from its first base, `CCC TAC CCT GTC`. Those 58 codons are of a frame the changed transcript does not have ([the read past the coding end](downstream-utr-frame-past-coding-end.md)). The record is in the `GRCh37-plugin-downstream-coding-end` corpus, whose `without-hgvs` flag set holds the run without the flag.

## 8. An insertion shifted onto a seam of an edited RefSeq model

**In short.** NM_015120.4 (ALMS1, GRCh38, forward strand) carries three bases the genome lacks, coding bases 75 to 77, between 2:73385942 and 2:73385943. The hand-written record `2 73385940 2_73385942_onto_d2 G GGA` inserts `GA` two bases before that seam, at `CDS_position` 72-73, and the notation moves it onto the seam: `c.73_74dup`, `HGVS_OFFSET=2`. Run with `--refseq`, vep-rs writes a `DownstreamProtein` of 17 residues and `ProteinLengthChange=-4128`, with and without `--hgvs`, and so does Ensembl VEP 116.2 without the flag. Under `--hgvs` Ensembl VEP writes 15 residues and `ProteinLengthChange=-4129`.

An insertion whose point lies on a seam takes the place of the bases the model carries there, in both engines. The shifted reading therefore puts the two bases in place of the model's three, where the record's own reading adds two and removes none: two changed sequences, three bases apart in length. Both fields differ, and the comparator accepts both under the class. No committed corpus holds the record, which needs Ensembl VEP's RefSeq cache.

## 9. A later allele of a multi-allelic record

**In short.** `21 47545369 . A AC,ACC`, a 1000 Genomes Phase 3 record on GRCh37, inserts one `C`, or two, between coding bases 478 and 479 of ENST00000413758 (COL6A2, forward strand). Both notations move seven bases: `HGVS_OFFSET=7`. On the second allele's row vep-rs writes a `DownstreamProtein` of 91 residues, from codon 160. Ensembl VEP 116.2 writes 89 residues, from codon 162, and prints `CDS_position` `485-486` where vep-rs prints `478-479`.

The alleles of a record share one object per transcript, and the protein notation of each allele it moves leaves the moved CDS pair on it (`TranscriptVariationAllele.pm`, lines 844-858 and 2406-2410). Line 1737 prints a row's `CDS_position` before line 1760 runs the row's own notation. A later allele's row therefore prints the record's span moved by the `HGVS_OFFSET` of the last earlier allele that has one ([the later alleles of a multi-allelic record](multiallelic-hgvs-shift.md#2-a-record-whose-first-allele-moves-the-coding-positions-of-the-second)). That is the row's own offset only where the alleles move alike, as here. The plugin still reads the pair at the row's own offset, with the allele turned by it.

`14 24646406 . G GGAAGAAGA,GGAAGAAGC`, of the same dataset, inserts eight bases between coding bases 681 and 682 of ENST00000311457 (REC8, forward strand), and its notations move seven bases and eight. On the second allele's row Ensembl VEP 116.2 prints `CDS_position` `688-689`, the span moved by seven, beside `HGVS_OFFSET=8`. Its plugin reads the pair eight bases on, 689-690, and writes 6 residues, from codon 230. vep-rs prints `681-682` and writes 8, from codon 228.

The comparator takes vep-rs's `CDS_position` as the own pair on such a row. The row must be a second or later alternate allele of its record, Ensembl VEP's pair the record's span moved by the row's `HGVS_OFFSET` or by the last earlier allele's, and vep-rs's the span itself. Ensembl VEP's value must be the plugin's rule at the row's own offset. The `NMD` reading is the same.

A pair moved by any other amount stays a difference, the first allele's offset included where an allele after it carries another. On any other row whose two engines print two `CDS_position` values, the plugin's values stay differences, the row of [section 10](#10-a-later-allele-whose-own-notation-moves-nothing) apart. The first allele's row of `21 47545369` prints one pair and differs as section 1 describes, three residues against five.

## 10. A later allele whose own notation moves nothing

**In short.** `14 29261309 . A AAAC,AC`, a record of 1000 Genomes Phase 3 (GRCh37), inserts three bases, or one, between coding bases 346 and 347 of ENST00000399387 (C14orf23, forward strand). The first allele's notation moves two bases, `HGVS_OFFSET=2`, and leaves the pair 348-349 on the object the alleles share. The second allele's notation moves nothing: neither engine's row carries an `HGVS_OFFSET`. Ensembl VEP 116.2 prints `CDS_position` `348-349` there, and its Downstream plugin reads that pair with the allele as written: 19 residues, from codon 117. vep-rs reads the record's own pair, `346-347`: 20 residues, from codon 116.

The second allele's row on each engine, `Extra` one field to a line, the differing cells in bold:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `14_29261309_A/AAAC/AC` | `14_29261309_A/AAAC/AC` |
| Location | `14:29261309-29261310` | `14:29261309-29261310` |
| Allele | `C` | `C` |
| Gene | `ENSG00000186960` | `ENSG00000186960` |
| Feature | `ENST00000399387` | `ENST00000399387` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `frameshift_variant` | `frameshift_variant` |
| cDNA_position | `450-451` | `450-451` |
| CDS_position | **`348-349`** | **`346-347`** |
| Protein_position | **`116-117`** | **`116`** |
| Amino_acids | **`-/TX`** | **`K/TX`** |
| Codons | `aaa/aCaa` | `aaa/aCaa` |
| Existing_variation | `-` | `-` |
| Extra: IMPACT | `HIGH` | `HIGH` |
| Extra: STRAND | `1` | `1` |
| Extra: HGVSc | `ENST00000399387.4:c.346_347insC` | `ENST00000399387.4:c.346_347insC` |
| Extra: HGVSp | **`ENSP00000382318.4:p.Thr117HisfsTer20`** | **`ENSP00000382318.4:p.Lys116ThrfsTer21`** |
| Extra: AA | `-` | `-` |
| Extra: DownstreamProtein | **`HKQTNNNNKQLSSTVSSRV`** | **`TNKQTNNNNKQLSSTVSSRV`** |
| Extra: HGVS_IntronEndOffset | `0` | `0` |
| Extra: HGVS_IntronStartOffset | `0` | `0` |
| Extra: LOEUF | `1.4690e+00` | `1.4690e+00` |
| Extra: NMD | `NMD_escaping_variant` | `NMD_escaping_variant` |
| Extra: NearestExonJB | `ENSE00001339307\|224\|start\|1349` | `ENSE00001339307\|224\|start\|1349` |
| Extra: ProteinLengthChange | `-32` | `-32` |
| Extra: gnomAD_100x_cov | `0.0000` | `0.0000` |
| Extra: gnomAD_10x_cov | `0.9976` | `0.9976` |
| Extra: gnomAD_15x_cov | `0.9902` | `0.9902` |
| Extra: gnomAD_1x_cov | `1.0000` | `1.0000` |
| Extra: gnomAD_20x_cov | `0.9268` | `0.9268` |
| Extra: gnomAD_25x_cov | `0.7559` | `0.7559` |
| Extra: gnomAD_30x_cov | `0.5224` | `0.5224` |
| Extra: gnomAD_50x_cov | `0.0220` | `0.0220` |
| Extra: gnomAD_5x_cov | `1.0000` | `1.0000` |
| Extra: gnomAD_mean_cov | `30.6790` | `30.6790` |
| Extra: gnomAD_median_cov | `30.0000` | `30.0000` |
| Extra: pHaplo | `0.432420029542035` | `0.432420029542035` |
| Extra: pTriplo | `0.270048787100055` | `0.270048787100055` |

`codon` deletes the kept pair only for an allele with a shift of its own (`TranscriptVariationAllele.pm`, lines 844-858), so the first allele's pair stays and `cds_start` returns it on the second allele's row (`BaseTranscriptVariation.pm`, lines 252-267). Lines 108 and 110 of `Downstream.pm` put the inserted base between coding bases 348 and 349. The record cannot be written there, since the base differs from the reference base after the insertion point, and the 19 residues are not the last 19 of the 20. The row's other bold cells differ for the same reason ([the later alleles of a multi-allelic record](multiallelic-hgvs-shift.md)).

The comparator takes a `DownstreamProtein` or a `ProteinLengthChange` on such a row under the class. The row must be a second or later alternate allele without an `HGVS_OFFSET`, vep-rs's pair the record's span, and Ensembl VEP's that span moved by the offset of the last earlier allele that has one. Ensembl VEP's value must be its rule on the printed pair with the allele as written and vep-rs's value vep-rs's rule on the own pair, and Ensembl VEP's rule must give another value there. An `NMD` on such a row stays a difference. So does a printed pair moved by any other amount.
