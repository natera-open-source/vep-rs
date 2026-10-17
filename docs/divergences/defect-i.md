# Defect I: a dot-notation breakend without `SVTYPE` is written as a pair

**The defect in one sentence.** A dot-notation single breakend (`N.`, `.N`) whose INFO carries no `SVTYPE` takes Ensembl VEP's sequence rules for its coordinates, so it is written as a span rather than as the point it names, the base after POS.

Its start is POS + 1 when the ALT's length differs from REF's and POS when it does not, and its end is POS plus the length of REF minus one. For a REF of one or two bases the span is the pair `POS-(POS+1)`, with the cDNA, CDS and protein positions of an insertion. For a REF of three or more bases beside an ALT of another length the span begins at POS + 1, and when the two lengths are equal (`NAC` against `NA.`) at POS.

**Why this is a defect.** The dot is the breakend's notation, and a single breakend is a point, the base after POS. The same record with `SVTYPE=BND` is that point in Ensembl VEP itself, so the pair follows from which INFO key is present, not from what the record says.

**What vep-rs does instead.** vep-rs reads every breakpoint record's span by the breakpoint rule, `SVTYPE` or not: the single base at POS + 1 unless `SVLEN` gives a span. When INFO names a mate through `CHR2` and `END2`, the `N.` row and the bracket row are both written at the point.

**Where it is in Ensembl VEP's source** (`ensembl-io` and `ensembl-vep`):
- `ensembl-io`, `BaseVCF4.pm`, `get_start`, lines 255-275: adds one to POS for a symbolic ALT or an ALT whose length differs from REF's, and leaves POS alone when the lengths are equal, as a two-base REF's is against `N.`.
- `ensembl-io`, `BaseVCF4.pm`, `get_end`, lines 297-316: ends a record without `SVLEN` or `END` at the start only when the raw INFO carries `SVTYPE` or the ALT carries `<`, `[` or `]`. A dot breakend carries none, so the last branch gives POS plus the length of REF minus one: POS for a one-base REF, POS + 1 for two bases, POS + 2 for three.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 237-245: the line is still routed to the structural path (`^\.\w+|\w+\.$`), and `get_SO_term` reads `BND` from the dot, so the record is a `chromosome_breakpoint` with a sequence record's coordinates.

**What correcting it would change.** The record would be the single base at POS + 1, as it is with `SVTYPE=BND` and as vep-rs writes it; [class 16](#16-a-dot-notation-breakend-without-svtype-written-as-a-pair) would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 16. A dot-notation breakend without `SVTYPE` written as a pair

**In short.** A dot-notation single breakend (`N.`, `.N`) is the base after POS. Ensembl VEP reads it as that point only when INFO carries `SVTYPE=BND`; without `SVTYPE` the length of REF sets a span, the coordinate pair `POS-(POS+1)` for a REF of one or two bases. vep-rs writes the point in both cases, and the comparator sets aside both engines' rows.

### 16.1 What a reader needs to know first

A single breakend in VCF 4.4 dot notation (`N.`, `.N`) is a point: the base after POS, where one side of the break is known and the other is not. Ensembl VEP reads it as a point when INFO carries `SVTYPE=BND`. Without `SVTYPE` the same ALT falls through the parser's sequence rules for its coordinates, and the length of REF sets the span it gets. The start moves to POS + 1 when the ALT's length differs from REF's and stays at POS when the two are equal. The end, with no `SVLEN`, no `END` and no `SVTYPE` in the raw INFO and no bracket in the ALT, is POS plus the length of REF minus one.

With a one-base REF (`N` against `N.`) the start is POS + 1 and the end POS, an insertion's coordinates. With a two-base REF (`NA` against `N.`) the lengths are equal, the start stays at POS and the end is POS + 1. Either way the record is written as the pair `POS-(POS+1)` with the cDNA, CDS and protein positions of an insertion.

With a REF of three or more bases (`NAC` against `.N`) the start is POS + 1 and the end POS plus the length of REF minus one, `(POS+1)-(POS+2)` for three bases, a span that begins at the point. When the ALT has the REF's length (`NAC` against `NA.` or `.AC`, `NACG` against `NAC.`) the start stays at POS and the span is `POS-(POS+2)` or `POS-(POS+3)`, the point its second base.

When INFO names a mate through `CHR2` and `END2`, still without `SVTYPE`, the mate is appended as a bracket allele and both the `N.` row and the `N[21:30425700[` row carry the pair.

### 16.2 The record

```
21	30425650	ndot_nosvtype	N	N.	.	.	.
21	30425650	ndot_nosvtype_mate	N	N.	.	.	CHR2=21;END2=30425700
21	30425650	ndot_nosvtype_ref3	NAC	.N	.	.	.
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

### 16.3 What each engine writes

Ensembl VEP 116.2, for the first record, writes two rows, one per feature, at the pair. Ensembl VEP 116.2 for the same record with `SVTYPE=BND`, and vep-rs for both, write the same two rows at the point. One table per row, the differing cells in bold.

The ENST00000334068 row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `ndot_nosvtype` | `ndot_nosvtype` |
| Location | **`21:30425650-30425651`** | **`21:30425651`** |
| Allele | `chromosome_breakpoint` | `chromosome_breakpoint` |
| Gene | `ENSG00000186971` | `ENSG00000186971` |
| Feature | `ENST00000334068` | `ENST00000334068` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `upstream_gene_variant` | `upstream_gene_variant` |

The ENST00000390690 row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `ndot_nosvtype` | `ndot_nosvtype` |
| Location | **`21:30425650-30425651`** | **`21:30425651`** |
| Allele | `chromosome_breakpoint` | `chromosome_breakpoint` |
| Gene | `ENSG00000240432` | `ENSG00000240432` |
| Feature | `ENST00000390690` | `ENST00000390690` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `coding_sequence_variant` | `coding_sequence_variant` |
| cDNA_position | **`276-277`** | **`276`** |
| CDS_position | **`262-263`** | **`262`** |
| Protein_position | `88` | `88` |

For the second record, whose INFO names a mate, Ensembl VEP 116.2 writes the pair under both alleles, four rows, and vep-rs the same four alleles at the point, `21:30425651`, with the positions `276`, `262` and `88`. One table per Ensembl VEP row; the vep-rs column points to that sentence.

The `N.` allele on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `ndot_nosvtype_mate` | see above |
| Location | `21:30425650-30425651` |  |
| Allele | `N.` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |

The `N[21:30425700[` allele on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `ndot_nosvtype_mate` | see above |
| Location | `21:30425650-30425651` |  |
| Allele | `N[21:30425700[` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |

The `N.` allele on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `ndot_nosvtype_mate` | see above |
| Location | `21:30425650-30425651` |  |
| Allele | `N.` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `feature_truncation,coding_sequence_variant` |  |
| cDNA_position | `276-277` |  |
| CDS_position | `262-263` |  |
| Protein_position | `88` |  |

The `N[21:30425700[` allele on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `ndot_nosvtype_mate` | see above |
| Location | `21:30425650-30425651` |  |
| Allele | `N[21:30425700[` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `feature_truncation,coding_sequence_variant` |  |
| cDNA_position | `276-277` |  |
| CDS_position | `262-263` |  |
| Protein_position | `88` |  |

For the third record, whose REF is three bases, Ensembl VEP 116.2 writes the span that begins at the point, two rows, and vep-rs the point, `21:30425651`, under the same allele on both features. One table per Ensembl VEP row; the vep-rs column points to that sentence.

The ENST00000334068 row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `ndot_nosvtype_ref3` | see above |
| Location | `21:30425651-30425652` |  |
| Allele | `chromosome_breakpoint` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |

The ENST00000390690 row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `ndot_nosvtype_ref3` | see above |
| Location | `21:30425651-30425652` |  |
| Allele | `chromosome_breakpoint` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `coding_sequence_variant` |  |

### 16.4 Why Ensembl VEP writes that

`BaseVCF4.pm` `get_start` (lines 255-275) adds one to POS for a symbolic ALT or an ALT whose length differs from REF's and leaves POS alone when the lengths are equal. `get_end` (lines 297-316) ends a record without `SVLEN` or `END` at the start only when the raw INFO carries `SVTYPE` or the ALT carries `<`, `[` or `]`. A dot breakend carries none, so the last branch gives POS plus the length of REF minus one.

The length of REF therefore sets the span. Against the two-character ALT a one-base REF starts at POS + 1 and ends at POS, a two-base REF starts at POS and ends at POS + 1, and both are written as the pair `POS-(POS+1)`. A REF of three or more bases starts at POS + 1 and ends at POS plus the length of REF minus one, `21:30425651-30425652` for `NAC` against `.N` at 30425650. Against an ALT of its own length it starts at POS, `21:30425650-30425652` for `NAC` against `NA.`.

`Parser/VCF.pm` lines 237-245 still route the line to the structural path (`^\.\w+|\w+\.$`), and `get_SO_term` reads `BND` from the dot. The record is therefore a `chromosome_breakpoint` with a sequence record's coordinates, an insertion's for a REF of one or two bases.

### 16.5 Why vep-rs's rows are the correct ones

The dot is the breakend's notation, and the same record with `SVTYPE=BND` is a point in Ensembl VEP itself. vep-rs reads every breakpoint record's span by the breakpoint rule, `SVTYPE` or not: the single base at POS + 1 unless `SVLEN` gives a span.

### 16.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tndot_nosvtype\tN\tN.\t.\t.\t.\n' > dot.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/dot.vcf -o /work/dot.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i dot.vcf -o dot.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

### 16.7 How the comparator treats this class

The class is `dot_breakend_without_svtype_pair` in the 116.2 registry of both comparators, set aside. The reference's row at a two-coordinate span and vep-rs's row at the point `chr:(P+1)`, on the same feature and under the same allele, are a pair, and both leave the adjusted figure, whatever their consequence sets. The point is the span's start or the base after it. It is the base after the start for a REF of one or two bases (the end of `chr:P-(P+1)`) and for `NAC` against `NA.` (the second base of `21:30425650-30425652`), and the start itself for `NAC` against `.N` (the start of `21:30425651-30425652`).

The allele is the bare `chromosome_breakpoint` of a record without a mate, or the `N.` and `N[21:30425700[` alleles of one whose INFO names a mate. A reference row on a feature vep-rs did not name stays charged.
