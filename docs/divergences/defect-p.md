# Defect P: the tests after the chromosome-aware overlap compare chromosome names the overlap did not resolve

**The defect in one sentence.** Depending on input order or on how the mate's chromosome is spelled, Ensembl VEP drops a bracket breakend's own-end rows, when a transcript of the mate's chromosome reaches it first, or its mate rows.

A bracket breakend names its mate inside the brackets of its allele, and Ensembl VEP writes its rows under two alleles, the own end and the mate. After the chromosome-aware overlap has handed a record to a transcript, two further tests compare chromosome names the overlap did not resolve the same way. A record's slice is the named chromosome it is placed on, and the slice test names a record that arrives without a slice after the first transcript that reaches it. The mate-allele test, which decides whether a row carries the mate allele, compares the mate's chromosome as written, and a mate spelled `Chr21` for `21` fails it.

**Why this is a defect.** A record's chromosome is a fact of its own line, and a mate's chromosome is a fact of its allele. Neither depends on which transcript the annotation pass happened to reach first, or on the spelling the file uses for a chromosome the cache knows under another name. Ensembl VEP's own overlap admits the transcript, and a later test then refuses the row for a reason that is not in the record. The result changes with the order of the input: the same record keeps its own-end rows in one file and loses them in another.

**What vep-rs does instead.** vep-rs annotates each breakend of a record on its own chromosome, resolved once from the record and the allele. The own end is annotated against the transcripts of CHROM within 5,000 bases of the base after POS, the mate against the transcripts of the mate's chromosome within 5,000 bases of the mate position. That holds whatever the order of the input and whichever name the file uses for the chromosome.

**Where it is in Ensembl VEP's source** (`ensembl-vep`; buffers, regions and chunks are described [below](#231-what-a-reader-needs-to-know-first)):
- `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, line 122: `$vf->{slice} ||= $slice ||= $tr->{slice}`, inside the loop over the records a transcript reaches: a record without a slice takes the slice of the first transcript that reaches it, and in offline mode (`--offline`) a record arrives without one.
- The same module, lines 125-143: the slice test. A `chromosome_breakpoint` record passes when the transcript's chromosome equals the record's slice name (line 128) or any of its breakends' `seq_region_name` with a leading `chr` stripped (lines 130-136); any other record passes only on its slice name (line 141). `next if $seq_mismatch` at line 143 drops the pair.
- `AnnotationSource.pm`, lines 235-258: the regions of a buffer are loaded in the order its records first ask for them, each record's own region before its breakend mates', a region listed once at its first request. `AnnotationType/Transcript.pm` lines 277-300 (`merge_features`) remove duplicates and keep that order, so the transcripts are visited region by region in it. So a transcript near a record's mate, in a region an earlier record asked for, can reach the record before any transcript of its own chromosome.
- `InputBuffer.pm`, lines 191-207 and 223-234: a buffer holds the records of one chromosome, so the records that share a buffer with a breakend are the ones on its own chromosome before it in the input.
- `InputBuffer.pm`, lines 348-357: `_chr_matches_for_vf_overlap` resolves a chromosome name through `get_source_chr_name`, the synonym table of other names the cache knows a chromosome by, so a mate written `NC_000002.12` reaches the transcripts of chromosome 2. The slice test at `AnnotationType/Transcript.pm` line 134 compares the raw name, `chr` stripped, and `NC_000002.12` is not `2`. A mate written `Chr21` passes the slice test (the strip is case-insensitive) and fails the mate-allele test below, whose synonym lookup is verbatim and whose fallback is a case-folded equality of the whole name.

(`ensembl-variation`): `StructuralVariationOverlap.pm` lines 86-102 and 144-160: a mate allele is added to a row only when `_close_to_feature` passes. That test compares the mate's `seq_region_name` with the transcript's through `_compare_seq_region_names` (`Utils/VariationEffect.pm` lines 111-124: the synonym table's keys verbatim, else a case-folded equality) and applies the 5,000-base window (`MAX_DISTANCE_FROM_TRANSCRIPT`, line 60). `StructuralVariationFeature.pm` lines 1300-1304: a breakend's slice is named after the chromosome exactly as the allele writes it.

**What correcting it would change.** If the slice test and the mate-allele test compared the record's own chromosome and each breakend's, resolved as the overlap resolves them, classes 23 and 24 would vanish. A bracket breakend's own-end rows and mate rows would be written whatever the order of the input and the spelling of the mate's chromosome.

**The classes it produces.**
- Class 23: the own-end rows of a bracket breakend whose mate is on another chromosome, lost when a mate-chromosome transcript names its slice first (the synthetic breakend set: 202 vep-rs tuples on 53 records in GRCh37, 1,012 on 71 in GRCh38). Set aside: vep-rs's rows, in the structural-variant comparator under the cache authority, the transcript cache that says which chromosome a transcript is on; counted and not set aside by the SNP/indel comparator.
- Class 24: the mate rows of a bracket breakend whose mate's chromosome is spelled in a way the later tests do not resolve, `Chr21` or an accession naming another chromosome (none on the measured sets). Set aside the same way.

A third shape of the same test, with no class number, is counted under the name `fasta_named_slice_intergenic` and never set aside. Under `--check_ref` or `--lookup_ref` with a FASTA whose sequence names differ from the cache's (`chr1` against `1`), every record's slice is named after the FASTA (`Parser.pm` line 749) and every transcript fails the slice test at line 141. So Ensembl VEP writes `intergenic_variant` where vep-rs writes transcript rows. The condition is the run's configuration, not a property of any record, and it reaches every record of the run; both comparators count it, and a run that shows the count is re-run with matching sequence names rather than scored. It is 0 on every measured set.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

**Defect D, retired.** Up to release 115.2 Ensembl VEP matched loaded transcripts to the records of a batch through one interval tree keyed on coordinates alone, with no chromosome. So a chromosome 21 record could be handed a chromosome X transcript whose coordinates overlapped its numbers and lose its `intergenic_variant` row to it (the paper's defect D, classes 9 and 10). Release 116.1 keys the overlap on the chromosome as well (`InputBuffer.pm` lines 286-345, `_chr_matches_for_vf_overlap` at lines 348-357). A record is handed to a transcript only when the record's chromosome, or one of its breakend mates', is the transcript's, names resolved through the synonym table. Neither class has a row against release 116.2.

## 23. A bracket breakend's own-end rows lost to the first transcript that reaches it

**In short.** `synth_bnd_inter_0062_B` is a chromosome 1 breakend whose mate lies in SOD1 on chromosome 21. Run alone, both engines write one row under its own-end allele and six under the mate allele. In the whole file Ensembl VEP drops the own-end row, because the file's first record, whose mate lies 1,016 bases from this record's mate, asked for the chromosome 21 region first. A SOD1 transcript therefore names the record's slice `21`, and the chromosome 1 transcript fails the slice test. vep-rs writes all seven rows in every context, and rightly: the own end's transcripts are a fact of the record's line, and Ensembl VEP writes that row whenever the order is the other way round.

### 23.1 What a reader needs to know first

A bracket breakend has two coordinates, its own (POS on CHROM) and its mate's (the chromosome and position inside the allele's brackets). Ensembl VEP writes rows under two alleles for it: the own end, printed as the reference base and a dot (`A.` or `.A`, the dot on the side the mate joins), and the mate, printed as the bracket string. The own-end rows are the transcripts within 5,000 bases of the base after POS on CHROM; the mate rows are the transcripts within 5,000 bases of the mate position on the mate's chromosome (class 11 describes what those rows should say).

Ensembl VEP annotates in batches of at most `--buffer_size` records (default 5,000) that never span two chromosomes, split into chunks for the `--fork` workers. For each chunk it loads the transcript-cache regions (1 Mb bins) its records overlap, a record's own region and each mate's, in the order the records ask for them. The transcripts are then visited region by region in that order, and a transcript that reaches a record decides, by a test on chromosome names, whether to annotate it. The term involved is `feature_truncation` (SO:0001906) with the region the breakpoint lands in, or `upstream_gene_variant` and `downstream_gene_variant` for a breakpoint near a transcript ([11.1](defect-e.md#111-what-a-reader-needs-to-know-first)).

### 23.2 The record

`synth_bnd_inter_0062_B`, the second record of the synthetic breakend set `tests/sv_validation/grch37/09_breakends.vcf.gz` (GRCh37), is a chromosome 1 breakend whose mate lies in intron 1 of SOD1 on chromosome 21. The block also carries the file's first record, `synth_bnd_inter_0066_B` (the record of class 11), whose mate lies 1,016 bases from it in the same intron:

```
1	14374333	synth_bnd_inter_0066_B	A	]21:33034355]A	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0066_A
1	16602616	synth_bnd_inter_0062_B	A	A]21:33033339]	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0062_A
```

The record's own breakpoint, 1:16602617, lies in an intron of ENST00000375592 (ENSG00000037637, reverse strand), the one chromosome 1 transcript within 5,000 bases of it. Its mate, 21:33033339, lies within 5,000 bases of six chromosome 21 transcripts, the SOD1 transcripts and their neighbours of [11.2](defect-e.md#112-the-record).

The GRCh38 set (`tests/sv_validation/grch38/09_breakends.vcf.gz`) shows the same shape on `synth_bnd_inter_0007_B`, whose own breakpoint 1:57985330 lies in an intron of the non-coding transcript ENST00000485760 (ENSG00000173406) and whose mate 21:25898426 lies in APP. It follows the file's first record, `synth_bnd_inter_0015_B`, whose mate lies 43 kb away in the same gene:

```
1	10305888	synth_bnd_inter_0015_B	A	]21:25941574]A	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0015_A
1	57985329	synth_bnd_inter_0007_B	A	[21:25898426[A	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0007_A
```

### 23.3 What each engine writes

Alone, both engines write seven rows for the GRCh37 record: one under its own-end allele `A.` and six under the mate allele `A]21:33033339]` (the mate rows differ between the engines by class 11, as the rows of [11.3](defect-e.md#113-what-each-engine-writes) do). In the whole file Ensembl VEP writes the six mate rows, unchanged, and nothing under `A.`; vep-rs writes the same seven rows as alone.

Ensembl VEP, whole file (2 of its 6 rows, all under the mate allele), one column per row:

| Column | Ensembl VEP, ENST00000270142 | Ensembl VEP, ENST00000449339 |
| --- | --- | --- |
| Uploaded_variation | `synth_bnd_inter_0062_B` | `synth_bnd_inter_0062_B` |
| Location | `1:16602617` | `1:16602617` |
| Allele | `A]21:33033339]` | `A]21:33033339]` |
| Gene | `ENSG00000142168` | `ENSG00000234509` |
| Feature | `ENST00000270142` | `ENST00000449339` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `feature_truncation` | `intergenic_variant` |
| cDNA_position | `-` | `-` |
| CDS_position | `-` | `-` |
| Protein_position | `-` | `-` |
| Amino_acids | `-` | `-` |
| Codons | `-` | `-` |
| Existing_variation | `-` | `-` |
| Extra | `IMPACT=HIGH;STRAND=1` | `IMPACT=MODIFIER;STRAND=-1` |

Ensembl VEP, whole file:

```
(no row under A.)
```

vep-rs, whole file (1 of 7 rows: the own-end row, identical to the row both engines write for the record alone):

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `synth_bnd_inter_0062_B` |
| Location | | `1:16602617` |
| Allele | | `A.` |
| Gene | | `ENSG00000037637` |
| Feature | | `ENST00000375592` |
| Feature_type | | `Transcript` |
| Consequence | | `feature_truncation,intron_variant` |
| cDNA_position | | `-` |
| CDS_position | | `-` |
| Protein_position | | `-` |
| Amino_acids | | `-` |
| Codons | | `-` |
| Existing_variation | | `-` |
| Extra | | `IMPACT=HIGH;STRAND=-1` |

The GRCh38 record has 42 rows on each side when run alone, one of them under its own-end allele `.A`; in the whole file Ensembl VEP writes the 41 mate rows and vep-rs all 42.

vep-rs, whole file (1 of 42 rows; Ensembl VEP writes the other 41):

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `synth_bnd_inter_0007_B` |
| Location | | `1:57985330` |
| Allele | | `.A` |
| Gene | | `ENSG00000173406` |
| Feature | | `ENST00000485760` |
| Feature_type | | `Transcript` |
| Consequence | | `feature_truncation,intron_variant,non_coding_transcript_variant` |
| cDNA_position | | `-` |
| CDS_position | | `-` |
| Protein_position | | `-` |
| Amino_acids | | `-` |
| Codons | | `-` |
| Existing_variation | | `-` |
| Extra | | `IMPACT=HIGH;STRAND=-1` |

### 23.4 Why Ensembl VEP writes that

A buffer holds records of one chromosome (`InputBuffer.pm` lines 191-207 and 223-234): the GRCh37 file's seventeen chromosome 1 records form one buffer, with `synth_bnd_inter_0066_B` first and `synth_bnd_inter_0062_B` second. The regions of a buffer are loaded in the order its records first ask for them, each record's own region before its mate's, a region asked for twice listed once, at its first request (`AnnotationSource.pm` lines 235-258). The first record asks for its own region, chromosome 1 at 14 Mb, then its mate's, chromosome 21 at 33 Mb. The second record's own region, chromosome 1 at 16 Mb, comes third, and its mate's region is the one already listed second.

The transcripts are visited in that order. The chromosome-aware overlap hands the record to the SOD1 transcripts near its mate (through the mate, `InputBuffer.pm` lines 336-342) and to ENST00000375592 near its own breakpoint (through its own coordinates, lines 324-326). The SOD1 transcripts come first, and line 122 of `AnnotationType/Transcript.pm` names the record's slice after the first of them, `21`. From then on the slice test (lines 125-143) passes a transcript only if its chromosome is `21`, the slice's name, or the mate's, also `21`; ENST00000375592 fails it and is dropped at line 143, and the own-end row is never created.

When the record's own region comes first in the load order, ENST00000375592 names the slice `1`, the SOD1 transcripts still pass through the mate, and both sets of rows are written. The whole GRCh38 file does the same to `synth_bnd_inter_0007_B`: the first chromosome 1 record's mate region, chromosome 21 at 25 Mb, holds the APP transcripts that reach the record's mate, and it is listed before the record's own region at 57 Mb.

### 23.5 Why vep-rs's row is the correct one

The transcripts within 5,000 bases of the record's own breakpoint on chromosome 1 are a fact of the record's own line, and Ensembl VEP's own admission rule gives the own-end allele exactly those rows. Whether they are written cannot depend on which record of the file asked for the mate's region first. Ensembl VEP agrees whenever the order is the other way round: alone, with its mate, under `--buffer_size 1`, or placed before the other record, it writes the own-end row, byte for byte the row vep-rs writes in every context ([23.6](#236-try-it-yourself)).

### 23.6 Try it yourself

Annotate the GRCh37 record with Ensembl VEP 116.2 (GRCh37 cache and FASTA, the command of [the closing section](reproducing.md)) in each of these inputs, and count its rows under the own-end allele `A.` and under the mate allele:

| Input | Ensembl VEP rows for `synth_bnd_inter_0062_B`, own end / mate |
| --- | --- |
| the record alone | 1 / 6 |
| the record with its mate record `synth_bnd_inter_0062_A` (`21 33033339 ... A [1:16602616[A`) | 1 / 6 |
| `synth_bnd_inter_0066_B` then the record (the file's first two lines) | 0 / 6 |
| the record then `synth_bnd_inter_0066_B` | 1 / 6 |
| a hand-written `1 14374333 x A ]21:33034355]A` with `SVTYPE=BND` and no MATEID, then the record | 0 / 6 |
| a hand-written chromosome 1 breakend at 14374333 whose mate is far from the record's, `]21:27367670]A`, then the record | 1 / 6 |
| a hand-written chromosome 1 breakend at a higher position, `1 16954684 x A ]21:33034355]A`, then the record (both own breakpoints lie in the same 1 Mb region, which is then listed first) | 1 / 6 |
| a hand-written chromosome 2 breakend, `2 1000000 x A ]21:33034355]A`, then the record (two buffers: a buffer holds one chromosome) | 1 / 6 |
| the SNV `21 33033339 . A C` then the record, or the SNV `1 14374333 . A C` then the record, or the SNV `1 16602616 . C A` then the record | 1 / 6 |
| the mate record `synth_bnd_inter_0066_A` (on chromosome 21) then the record | 1 / 6 |
| the whole `09_breakends.vcf.gz`, `--fork 4 --buffer_size 5000` (the reference setting), or `--fork 1 --buffer_size 5000` | 0 / 6 |
| the whole file at `--buffer_size 1`, `--fork 1` or `--fork 4` | 1 / 6 |

vep-rs writes 1 / 6 for every one of those inputs. The loss needs, in the record's own buffer and before it, a record that asks for the region holding the transcripts within reach of the record's mate before the record's own region is asked for. `--fork` does not change the result and `--buffer_size 1`, which annotates every record alone, removes it.

The GRCh38 record behaves the same way: 1 / 41 alone, with its mate, after the SNV `21 25898426 . T A` at its mate's position, or after SNVs at the first record's two breakpoints. It gets 0 / 41 after the file's first record `synth_bnd_inter_0015_B`, after `synth_bnd_inter_0037_B` (`1 11286046 ... A A[21:26072808[`, whose mate's region holds the APP transcripts that reach this record's mate), and in the whole file.

### 23.7 How the comparator treats this class

`filter_bnd_own_rows_lost` in `scripts/validation/compare_sv_concordance.py` runs under the cache authority (`--vep-rs-cache`, the transcript cache that says which chromosome a transcript is on). It sets aside a vep-rs tuple when all of these hold:

- the Location is a bracket breakend's whose mate is on another chromosome;
- Ensembl VEP wrote a Transcript row at the Location under the bracket allele and none under any own-end allele (`N.`, `.N` or `chromosome_breakpoint`);
- the tuple's transcript is one the cache files on the Location's chromosome within 5,000 bases of the breakpoint, Ensembl VEP's own admission rule (`StructuralVariationOverlap.pm` lines 144-160).

A vep-rs own-end row on a transcript farther than that, or on another chromosome, stays charged.

The SNP/indel comparator, `scripts/concordance/compare_vep_outputs.py`, has no transcript-to-chromosome authority, so it counts the shape under the same name and sets nothing aside. The shape is a Location at which Ensembl VEP wrote rows under a bracket allele and none under an own-end allele, and vep-rs wrote own-end rows.

Counts (`excluded_bnd_own_rows_lost_rust`, per file in `excluded_bnd_own_rows_lost_by_file`):

- GRCh37: 202 vep-rs tuples on 53 records, all in `09_breakends.vcf.gz`: the records on chromosomes 1 to 5 whose mates are on chromosome 21 (12 on chromosome 1, 17 on 2, 8 on 3, 5 on 4, 11 on 5). None of the set's chromosome 21 records, and none in any other file.
- GRCh38: 1,012 vep-rs tuples on 71 records, all in `09_breakends.vcf.gz` (21 on chromosome 1, 18 on 2, 12 on 3, 11 on 4, 9 on 5).

The gnomAD-SV files, the 1000 Genomes files and the other synthetic files show none. Fully masked, vep-rs side.

## 24. A bracket breakend's mate rows lost to an unresolved chromosome prefix

**In short.** A hand-written chromosome 21 breakend with its mate written `N[Chr21:30425700[`, 50 bases away, reaches its transcripts through its own coordinates and passes the slice test. The mate-allele test then compares `Chr21` with the transcript's `21` in two steps: through the synonym table, whose keys hold `chr21` and not `Chr21`, then by a case-folded equality, which `chr21` against `21` fails. So Ensembl VEP writes the own-end rows and nothing under the mate allele. vep-rs reads `Chr21` as chromosome 21 and writes the mate rows too. The mate's chromosome is a fact of the allele, and Ensembl VEP writes the same mate rows when the allele spells it `chr21` or `21`, so vep-rs's rows are the correct ones.

### 24.1 What a reader needs to know first

A VCF may name a chromosome by any of the names its reference assembly report lists: `21`, `chr21`, `NC_000021.9`, `CM000683.2`. Ensembl VEP's cache carries a synonym table for them (`--synonyms` names a file that replaces it, `CacheDir.pm` line 254), and its chromosome-aware overlap resolves a mate's name through that table before matching it to a transcript's chromosome (`InputBuffer.pm` lines 348-357). Two later tests compare names again. The slice test of `AnnotationType/Transcript.pm` lines 125-143 strips a leading `chr` from a breakend's name, case-insensitively, and compares the rest as written. `_close_to_feature` (`StructuralVariationOverlap.pm` lines 144-160) accepts a synonym table key verbatim or a case-folded equality of the whole name. The mate rows are written only when every test passes.

### 24.2 The record

Hand-written, GRCh38: a breakend at 21:30425650 whose mate lies 50 bases away on the same chromosome, written with a capitalised prefix:

```
21	30425650	x	N	N[Chr21:30425700[	.	.	SVTYPE=BND
```

The own breakpoint, 21:30425651, lies in the coding sequence of ENST00000390690 (reverse strand) and 4,569 bases upstream of ENST00000334068; the mate, 50 bases away, reaches the same two transcripts. The GRCh37 counterpart is `21 33034300 x N N[Chr21:33034355[` with `SVTYPE=BND`, the breakpoint of class 11 as its mate, which reaches seven transcripts.

### 24.3 What each engine writes

Ensembl VEP (both rows) and vep-rs (all four rows), one table per transcript, the differing cells in bold. ENST00000334068:

| Column | Ensembl VEP | vep-rs, own end | vep-rs, mate |
| --- | --- | --- | --- |
| Uploaded_variation | `x` | `x` | `x` |
| Location | `21:30425651` | `21:30425651` | `21:30425651` |
| Allele | `N.` | `N.` | **`N[Chr21:30425700[`** |
| Gene | `ENSG00000186971` | `ENSG00000186971` | `ENSG00000186971` |
| Feature | `ENST00000334068` | `ENST00000334068` | `ENST00000334068` |
| Feature_type | `Transcript` | `Transcript` | `Transcript` |
| Consequence | `upstream_gene_variant` | `upstream_gene_variant` | `upstream_gene_variant` |
| cDNA_position | `-` | `-` | `-` |
| CDS_position | `-` | `-` | `-` |
| Protein_position | `-` | `-` | `-` |
| Amino_acids | `-` | `-` | `-` |
| Codons | `-` | `-` | `-` |
| Existing_variation | `-` | `-` | `-` |
| Extra | `IMPACT=MODIFIER;DISTANCE=4569;STRAND=1` | `IMPACT=MODIFIER;DISTANCE=4569;STRAND=1` | `IMPACT=MODIFIER;DISTANCE=4569;STRAND=1` |

ENST00000390690:

| Column | Ensembl VEP | vep-rs, own end | vep-rs, mate |
| --- | --- | --- | --- |
| Uploaded_variation | `x` | `x` | `x` |
| Location | `21:30425651` | `21:30425651` | `21:30425651` |
| Allele | `N.` | `N.` | **`N[Chr21:30425700[`** |
| Gene | `ENSG00000240432` | `ENSG00000240432` | `ENSG00000240432` |
| Feature | `ENST00000390690` | `ENST00000390690` | `ENST00000390690` |
| Feature_type | `Transcript` | `Transcript` | `Transcript` |
| Consequence | `feature_truncation,coding_sequence_variant` | `feature_truncation,coding_sequence_variant` | `feature_truncation,coding_sequence_variant` |
| cDNA_position | `276` | `276` | `276` |
| CDS_position | `262` | `262` | `262` |
| Protein_position | `88` | `88` | `88` |
| Amino_acids | `-` | `-` | `-` |
| Codons | `-` | `-` | `-` |
| Existing_variation | `-` | `-` | `-` |
| Extra | `IMPACT=HIGH;STRAND=-1` | `IMPACT=HIGH;STRAND=-1` | `IMPACT=HIGH;STRAND=-1` |

Ensembl VEP writes the two own-end rows and nothing under the mate allele; vep-rs writes the same two own-end rows and the two mate rows, identical to the own-end rows in every column but Allele because the two breakpoints are 50 bases apart. On GRCh37 the record writes 7 rows against 14 the same way.

### 24.4 Why Ensembl VEP writes that

The chromosome 21 transcripts reach the record through its own coordinates and pass the slice test on the record's slice, `21`. For each of them the mate allele is added to the row only if `_close_to_feature` passes for the breakend (`StructuralVariationOverlap.pm` lines 90-102). Its name test, `_compare_seq_region_names` (`Utils/VariationEffect.pm` lines 111-124), finds `Chr21` neither among the synonym table's keys, which hold `chr21`, nor equal to `21` once both are lower-cased. The own-end rows are written; the mate allele is never added. With the mate written `chr21` or `21` the same test passes and the mate rows are written.

A mate written with a RefSeq accession takes the other route. On another chromosome, `N[NC_000002.12:321682[`, the overlap resolves the accession (`InputBuffer.pm` line 356) and hands the record to the chromosome 2 transcripts. The breakend's `seq_region_name` is the accession as the allele wrote it (`StructuralVariationFeature.pm` lines 1300-1304). The slice test then compares it, `chr` stripped, with the transcript's `2` (`AnnotationType/Transcript.pm` lines 130-136), finds no match and drops the pair at line 143. Ensembl VEP writes no mate row. vep-rs resolves the accession through the synonym table, as it resolves the record's own chromosome, and writes the mate rows it writes for `2`; the comparator sets them aside ([24.7](#247-how-the-comparator-treats-this-class)).

On the record's own chromosome, `N[NC_000021.9:30425700[`, the transcripts reach the record through its own coordinates, the accession is a key of the synonym table, and Ensembl VEP writes the mate rows. vep-rs resolves the accession the same way and writes the same mate rows, so the two outputs agree.

### 24.5 Why vep-rs's row is the correct one

The mate's chromosome is a fact of the allele, and `Chr21` names the chromosome the cache holds as `21` and lists as `chr21`. The transcripts within 5,000 bases of the mate position on that chromosome are the mate rows under Ensembl VEP's own admission rule, and Ensembl VEP writes them itself when the allele spells the chromosome `chr21` or `21` ([24.6](#246-try-it-yourself)). A prefix that differs from a listed synonym by its case cannot be a reason to withhold them, and an accession the synonym table lists names its chromosome as surely: Ensembl VEP's own overlap resolves it before the slice test refuses the pair on spelling.

### 24.6 Try it yourself

Annotate the record with each spelling of the mate's chromosome, with Ensembl VEP 116.2 and vep-rs (GRCh38 cache and FASTA, the commands of [the closing section](reproducing.md); the GRCh37 counterpart of [24.2](#242-the-record) in parentheses):

| Mate | Ensembl VEP 116.2 rows | vep-rs rows |
| --- | --- | --- |
| `N[21:30425700[` or `N[chr21:30425700[` | 4 (14): own end and mate | 4 (14) |
| `N[Chr21:30425700[` | 2 (7): own end only | 4 (14) |
| `N[NC_000021.9:30425700[` (GRCh37 `NC_000021.8`) | 4 (14) | 4 (14) |
| `N[2:321682[` or `N[chr2:321682[` (GRCh37 `Y:15437971`, `chrY:15437971`) | 13 (15) | 13 (15); the mate rows differ by class 11 |
| `N[NC_000002.12:321682[` (GRCh37 `NC_000024.9:15437971`) | 2 (7): own end only | 13 (15), as for `2` (`Y` on GRCh37), the chromosome the accession resolves to; the mate rows are set aside |

The capitalised prefix and an accession naming another chromosome produce the class: Ensembl VEP writes the own-end rows alone where vep-rs also writes the mate rows.

### 24.7 How the comparator treats this class

`filter_bnd_synonym_mate_rows` in `scripts/validation/compare_sv_concordance.py` runs under the cache authority. It sets aside a vep-rs tuple under a bracket allele when all of these hold:

- the mate name, `chr` stripped, is not a chromosome the cache files;
- Ensembl VEP wrote rows at the Location and none under the bracket allele;
- the tuple's transcript is not on the Location's own chromosome.

It also sets a tuple aside when the mate name is the record's own chromosome spelled with a prefix Ensembl VEP does not resolve (`Chr21`), in which case vep-rs's mate rows on the own chromosome are set aside too. The SNP/indel comparator counts the shape and sets nothing aside, a synonym there being any mate name outside the canonical contig set.

Counts (`excluded_bnd_synonym_mate_rows_rust`): 0 on every measured set. Every mate in the gnomAD-SV files is named through `CHR2` as `chr21`, `chrX` or `chrY`, synonyms the cache lists, and the synthetic breakends write their mates with the cache's own names.
