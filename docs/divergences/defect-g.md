# Defect G: a coordinate past the 32-bit bound turns the buffer intergenic

**The defect in one sentence.** When any coordinate of an input buffer plus the `--distance` flank reaches 2^31 (a breakend's mate position, a span's end or a record's own position), Ensembl VEP writes every record of the buffer as `intergenic_variant`.

Ensembl VEP adds the flank to every coordinate of the buffer when it chooses the cache regions and the buffer's bounds.

**Why this is a defect.** A record's annotation depends on its own coordinates and the transcripts it overlaps, not on the other records read with it. Here the record that supplied the coordinate loses its transcript rows under all its alleles, and so does every other record loaded with it, whatever transcripts they overlap, while the output does not identify which records were damaged.

**What vep-rs does instead.** vep-rs holds every coordinate in 64 bits. It writes the breakend's own-end rows against the transcripts at its own position, writes nothing under the mate allele when no transcript lies at the mate, and annotates a span to 2,147,483,647 against every transcript it covers. It annotates the other records of the batch as if the record were absent, which is what Ensembl VEP itself writes for the same records at any coordinate below the bound.

**Where it is in Ensembl VEP's source** (`ensembl-vep`, `Bio/EnsEMBL/VEP/AnnotationSource.pm`):
- line 190: `get_regions_from_coords` adds the `--distance` flank to each coordinate it is given when it chooses the cache regions and the buffer's bounds.
- lines 251-254: `get_all_regions_by_InputBuffer` hands it every breakend's mate position, beside every record's start and end.
- lines 306-309: `_check_overlap` keeps a transcript only when it overlaps the bounds with the flank added again.

**What correcting it would change.** The buffer's records would keep their transcript rows and the breakend its own-end rows, as vep-rs writes them; [class 14](#14-a-coordinate-at-or-past-2147478648-a-breakend-mate-a-span-end-or-the-records-own-position-turns-the-references-buffer-intergenic) would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 14. A coordinate at or past 2,147,478,648 (a breakend mate, a span end or the record's own position) turns the reference's buffer intergenic

**In short.** When any coordinate of an input buffer plus the `--distance` flank reaches 2^31, the bound of a signed 32-bit integer (with the default 5,000-base flank, any coordinate from 2,147,478,648), Ensembl VEP writes every record of the buffer as `intergenic_variant`. The record that supplied the coordinate loses its transcript rows under all its alleles, and so does every other record loaded with it. vep-rs holds every coordinate in 64 bits and annotates each record against the transcripts it overlaps, which is what Ensembl VEP itself writes for the same records below the bound.

### 14.1 What a reader needs to know first

Ensembl VEP annotates its input in buffers of 5,000 records. For each buffer it works out the span of coordinates the records cover, adds the `--distance` flank (5,000 bases by default) at both ends, and loads the transcripts of the cache regions inside that span. Every record's start and end and every breakend's mate position are part of the span.

When any coordinate of the buffer plus the flank reaches 2^31 (2,147,483,648), every record of the buffer is written as `intergenic_variant`: the record that supplied the coordinate, under all its alleles, and every other record loaded with it, whatever transcripts they overlap. The coordinate can be a breakend's mate position, a span's end (`END=2147483647`, or `SVLEN=2147483647` from POS 30425650, each annotated under `--max_sv_size -1`) or a record's own position (`21 2147483648 x N N[21:30425700[ . . SVTYPE=BND`).

The boundary follows the flank exactly: with the default flank the first position that does this is 2,147,478,648 and 2,147,478,647 annotates normally; under `--distance 1` the first such position is 2,147,483,647 and 2,147,483,646 annotates normally. No human chromosome reaches these positions, so the shape appears only in a malformed or synthetic record, and when it does the records it damages are its neighbours in the buffer, which the output does not identify.

### 14.2 The record

```
21	30425650	control_dup	N	<DUP>	.	.	SVTYPE=DUP;END=30425899
21	30425650	end2_2147478648	N	<BND>	.	.	SVTYPE=BND;CHR2=21;END2=2147478648
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

### 14.3 What each engine writes

Ensembl VEP 116.2, three rows, every one intergenic, one column per row:

| Column | `control_dup` | `end2_2147478648`, allele `N.` | `end2_2147478648`, mate allele |
| --- | --- | --- | --- |
| Uploaded_variation | `control_dup` | `end2_2147478648` | `end2_2147478648` |
| Location | `21:30425651-30425899` | `21:30425651` | `21:30425651` |
| Allele | `duplication` | `N.` | `N[21:2147478648[` |
| Gene | `-` | `-` | `-` |
| Feature | `-` | `-` | `-` |
| Feature_type | `-` | `-` | `-` |
| Consequence | `intergenic_variant` | `intergenic_variant` | `intergenic_variant` |
| Extra | `ALLELE_NUM=1;IMPACT=MODIFIER` | `ALLELE_NUM=0;IMPACT=MODIFIER` | `ALLELE_NUM=1;IMPACT=MODIFIER` |

The same file with `END2=2147478646`, four rows, beside vep-rs, for the first file, four rows: the same four as Ensembl VEP's second file, under `end2_2147478648`. One table per row, the differing cells in bold.

`control_dup` on ENST00000334068:

| Column | Ensembl VEP, `END2=2147478646` | vep-rs, `END2=2147478648` |
| --- | --- | --- |
| Uploaded_variation | `control_dup` | `control_dup` |
| Location | `21:30425651-30425899` | `21:30425651-30425899` |
| Allele | `duplication` | `duplication` |
| Gene | `ENSG00000186971` | `ENSG00000186971` |
| Feature | `ENST00000334068` | `ENST00000334068` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `upstream_gene_variant` | `upstream_gene_variant` |
| Extra | **`ALLELE_NUM=1;IMPACT=MODIFIER;STRAND=1`** |  |

`control_dup` on ENST00000390690:

| Column | Ensembl VEP, `END2=2147478646` | vep-rs, `END2=2147478648` |
| --- | --- | --- |
| Uploaded_variation | `control_dup` | `control_dup` |
| Location | `21:30425651-30425899` | `21:30425651-30425899` |
| Allele | `duplication` | `duplication` |
| Gene | `ENSG00000240432` | `ENSG00000240432` |
| Feature | `ENST00000390690` | `ENST00000390690` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `feature_elongation,coding_sequence_variant` | `feature_elongation,coding_sequence_variant` |
| Extra | **`ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1`** |  |

The breakend's own-end allele `N.` on ENST00000334068:

| Column | Ensembl VEP, `END2=2147478646` | vep-rs, `END2=2147478648` |
| --- | --- | --- |
| Uploaded_variation | **`end2_2147478646`** | **`end2_2147478648`** |
| Location | `21:30425651` | `21:30425651` |
| Allele | `N.` | `N.` |
| Gene | `ENSG00000186971` | `ENSG00000186971` |
| Feature | `ENST00000334068` | `ENST00000334068` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `upstream_gene_variant` | `upstream_gene_variant` |
| Extra | **`ALLELE_NUM=0;IMPACT=MODIFIER;STRAND=1`** |  |

The breakend's own-end allele `N.` on ENST00000390690:

| Column | Ensembl VEP, `END2=2147478646` | vep-rs, `END2=2147478648` |
| --- | --- | --- |
| Uploaded_variation | **`end2_2147478646`** | **`end2_2147478648`** |
| Location | `21:30425651` | `21:30425651` |
| Allele | `N.` | `N.` |
| Gene | `ENSG00000240432` | `ENSG00000240432` |
| Feature | `ENST00000390690` | `ENST00000390690` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `feature_truncation,coding_sequence_variant` | `feature_truncation,coding_sequence_variant` |
| Extra | **`ALLELE_NUM=0;IMPACT=HIGH;STRAND=-1`** |  |

Both engines write no row under the mate allele for a mate this far from any transcript, at 2,147,478,646 as at 46,709,983 or 1,000,000,000; the difference is the transcript rows of the whole buffer.

A span's end or a record's own position past the bound turns the buffer the same way. Under `--max_sv_size -1`, `21 30425650 x N <DEL> . . SVTYPE=DEL;END=2147483647` is one row in Ensembl VEP 116.2 (`SVLEN=2147483647` in place of the `END` gives `21:30425651-2177909297 deletion - intergenic_variant`) and 5,074 transcript rows in vep-rs. The vep-rs rows are every transcript from 30425651 to the end of the chromosome.

Ensembl VEP 116.2, one row, and vep-rs (1 of 5,074 rows), the differing cells in bold:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `del_end_2e31m1` | `del_end_2e31m1` |
| Location | `21:30425651-2147483647` | `21:30425651-2147483647` |
| Allele | `deletion` | `deletion` |
| Gene | **`-`** | **`ENSG00000141956`** |
| Feature | **`-`** | **`ENST00000269844`** |
| Feature_type | **`-`** | **`Transcript`** |
| Consequence | **`intergenic_variant`** | **`transcript_ablation`** |
| Extra | **`ALLELE_NUM=1;IMPACT=MODIFIER`** |  |

The breakend `21 2147483648 x N N[21:30425700[ . . SVTYPE=BND`, whose own position is past the bound, is two intergenic rows in Ensembl VEP 116.2, one per allele. vep-rs writes the mate allele's rows at the mate (its own end, at 2,147,483,649, reaches no transcript on either side).

Ensembl VEP 116.2, two rows, one per allele, and vep-rs, two rows, both under the mate allele. The own-end allele `N.`:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `pos_2e31_brk` | no row |
| Location | `21:2147483649` |  |
| Allele | `N.` |  |
| Gene | `-` |  |
| Feature | `-` |  |
| Feature_type | `-` |  |
| Consequence | `intergenic_variant` |  |
| Extra | `ALLELE_NUM=0;IMPACT=MODIFIER` |  |

The mate allele `N[21:30425700[`, Ensembl VEP's one row against vep-rs's two:

| Column | Ensembl VEP | vep-rs, ENST00000334068 | vep-rs, ENST00000390690 |
| --- | --- | --- | --- |
| Uploaded_variation | `pos_2e31_brk` | `pos_2e31_brk` | `pos_2e31_brk` |
| Location | `21:2147483649` | `21:2147483649` | `21:2147483649` |
| Allele | `N[21:30425700[` | `N[21:30425700[` | `N[21:30425700[` |
| Gene | **`-`** | **`ENSG00000186971`** | **`ENSG00000240432`** |
| Feature | **`-`** | **`ENST00000334068`** | **`ENST00000390690`** |
| Feature_type | **`-`** | **`Transcript`** | **`Transcript`** |
| Consequence | **`intergenic_variant`** | **`upstream_gene_variant`** | **`feature_truncation,coding_sequence_variant`** |
| Extra | **`ALLELE_NUM=1;IMPACT=MODIFIER`** |  |  |

At POS 2147483647 the breakend's own position plus the flank also passes 2^31 and the rows are the same; a mate at 2,147,478,647 annotates normally ([14.1](#141-what-a-reader-needs-to-know-first)).

### 14.4 Why Ensembl VEP writes that

`AnnotationSource.pm` `get_regions_from_coords` (line 190) adds the `--distance` flank to each coordinate it is given when it chooses the cache regions and the buffer's bounds. `get_all_regions_by_InputBuffer` hands it every record's start and end and every breakend's mate position (lines 251-254 for the mates). `_check_overlap` (lines 306-309) then keeps a transcript only when it overlaps the bounds with the flank added again.

The buffer's largest coordinate plus the flank is where the behaviour turns, whichever record supplies it and whether it is a mate position, a span's end or a start. The measurements of [14.1](#141-what-a-reader-needs-to-know-first) place the turn at 2^31 for both flank sizes, the bound of a signed 32-bit integer, and past it every record of the buffer is written as intergenic.

### 14.5 Why vep-rs's rows are the correct ones

A record's annotation depends on its own coordinates and the transcripts it overlaps, not on the other records read with it. vep-rs holds every coordinate in 64 bits. It writes the breakend's own-end rows against the transcripts at its own position, writes nothing under the mate allele when no transcript lies at the mate, and annotates a span to 2,147,483,647 against every transcript it covers. It annotates the other records of the batch as if the record were absent, which is what Ensembl VEP itself writes for the same records at any coordinate below the bound.

### 14.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tcontrol_dup\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=30425899\n21\t30425650\tend2_2147478648\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=2147478648\n' > overflow.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/overflow.vcf -o /work/overflow.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i overflow.vcf -o overflow.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

Changing `END2` to `2147478646` restores the transcript rows in Ensembl VEP's output; adding `--distance 1` moves the turn to `END2=2147483647`. Replacing the breakend line with `21\t30425650\tdel_end_2e31m1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=2147483647` and adding `--max_sv_size -1` to both commands shows the same turn from a span's end.

### 14.7 How the comparator treats this class

The class is `bnd_mate_32bit_buffer_intergenic` in the 116.2 registry of both comparators, counted and never set aside. The records it damages beyond the record's own are its neighbours in a buffer whose boundaries the outputs do not record, so no tuple can be attributed to it with certainty, and crediting vep-rs for an unattributable row is not done.

The SNP/indel comparator counts the reference's rows at the Location of a record whose bracket allele names a mate at or past 2,147,478,648.

The structural-variant comparator counts such records from the input (a bracket ALT, or `CHR2` with `END2`, else `END`) and the reference's tuples at their Locations (`bnd_mate_32bit_buffer_records`, `bnd_mate_32bit_buffer_reference_tuples`).

A span's end (`SVLEN`, else `END`) or a record's own POS at or past the bound is counted in the same way as a mate position: both comparators count a reference Location whose start or end is at or past 2,147,478,648. The rows of such a record stay counted, not set aside. None of the measured sets carries the shape.
