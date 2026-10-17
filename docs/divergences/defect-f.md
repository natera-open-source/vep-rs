# Defect F: a breakpoint record without a mate gives its reference allele a row

**The defect in one sentence.** For a breakpoint record that names no mate, Ensembl VEP gives the reference allele a consequence row of its own.

Such a record (`<BND>` with no `CHR2`, or with `CHR2` but neither `END2` nor `END`) is held by Ensembl VEP as the allele string `N/<BND>`, which parses to no breakend. The overlap object therefore builds one allele per slash-separated part of the string, the reference base included.

**Why this is a defect.** A reference allele has no variant effect; a consequence row describes what an alternate allele does to a feature. The second row restates the first under an allele number that points at the reference base, and nothing in it is information about the variant.

**What vep-rs does instead.** vep-rs holds the record as one single breakend, prints the class term `chromosome_breakpoint` and the terms of the span or the point once, and takes no `feature_truncation` for it, as Ensembl VEP's own ALT row does not.

**Where it is in Ensembl VEP's source** (`ensembl-vep` and `ensembl-variation`):
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 491-497: a bracket mate is built only under `defined $info->{CHR2}` with `END2 || END`; line 499 prepends the reference base to the allele string unless it opens or closes with `.`, so the record's string is `N/<BND>`.
- `ensembl-variation`, `StructuralVariationFeature.pm`, `_parse_breakends`, lines 1274-1283: finds no bracket in the string, so the record has no breakends.
- `ensembl-variation`, `StructuralVariationOverlap.pm`, lines 64-77: with no breakends, `split(/\//, $vf->{allele_string})` gives one allele per part, numbered from 1: the reference base `N` is allele 1 and `<BND>` allele 2.
- `ensembl-vep`, `OutputFactory.pm`, line 2029: `$vfoa->{breakend}->{string} || $svf->class_SO_term` prints `chromosome_breakpoint` for both alleles, neither of which carries a breakend string.

**What correcting it would change.** The reference-allele row would disappear and every row of the record would be written once, as vep-rs writes it; [class 13](#13-the-reference-allele-row-of-a-bare-breakpoint) would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 13. The reference-allele row of a bare breakpoint

**In short.** A `<BND>` record that names no mate is held by Ensembl VEP as the allele string `N/<BND>`, which parses to no breakend, so it builds one allele per slash-separated part of the string, the reference base included. Every row of the record is therefore written twice, once per allele, distinguished only by `ALLELE_NUM`. vep-rs holds the record as one single breakend and writes each row once, matching Ensembl VEP's own ALT row.

### 13.1 What a reader needs to know first

Ensembl VEP holds a VCF breakend that names no mate as the allele string `N/<BND>`: the reference base, a slash, the symbolic ALT. Such a breakend is a `<BND>` with no `CHR2`, or with `CHR2` but neither `END2` nor `END`, or with a `POS2` Ensembl VEP never reads. The string parses to no breakend, so the overlap object builds its alleles by splitting the string on `/`, one allele per part, and the first part is the reference base. Both alleles print the class term `chromosome_breakpoint`, so every row of the record appears twice, distinguished only by `ALLELE_NUM`.

### 13.2 The record

```
21	30425650	bare_bnd	N	<BND>	.	.	SVTYPE=BND
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

### 13.3 What each engine writes

Ensembl VEP 116.2, four rows for two transcripts, and vep-rs, two rows, one table per transcript, the differing cells in bold. First ENST00000334068:

| Column | Ensembl VEP, `ALLELE_NUM=1` | Ensembl VEP, `ALLELE_NUM=2` | vep-rs |
| --- | --- | --- | --- |
| Uploaded_variation | `bare_bnd` | `bare_bnd` | `bare_bnd` |
| Location | `21:30425651` | `21:30425651` | `21:30425651` |
| Allele | `chromosome_breakpoint` | `chromosome_breakpoint` | `chromosome_breakpoint` |
| Gene | `ENSG00000186971` | `ENSG00000186971` | `ENSG00000186971` |
| Feature | `ENST00000334068` | `ENST00000334068` | `ENST00000334068` |
| Feature_type | `Transcript` | `Transcript` | `Transcript` |
| Consequence | `upstream_gene_variant` | `upstream_gene_variant` | `upstream_gene_variant` |
| Extra | **`ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4569;STRAND=1`** | **`ALLELE_NUM=2;IMPACT=MODIFIER;DISTANCE=4569;STRAND=1`** |  |

Record `bare_bnd` on ENST00000390690, Ensembl VEP's two rows against vep-rs's one:

| Column | Ensembl VEP, `ALLELE_NUM=1` | Ensembl VEP, `ALLELE_NUM=2` | vep-rs |
| --- | --- | --- | --- |
| Uploaded_variation | `bare_bnd` | `bare_bnd` | `bare_bnd` |
| Location | `21:30425651` | `21:30425651` | `21:30425651` |
| Allele | `chromosome_breakpoint` | `chromosome_breakpoint` | `chromosome_breakpoint` |
| Gene | `ENSG00000240432` | `ENSG00000240432` | `ENSG00000240432` |
| Feature | `ENST00000390690` | `ENST00000390690` | `ENST00000390690` |
| Feature_type | `Transcript` | `Transcript` | `Transcript` |
| Consequence | `coding_sequence_variant` | `coding_sequence_variant` | `coding_sequence_variant` |
| Extra | **`ALLELE_NUM=1;IMPACT=MODIFIER;STRAND=-1`** | **`ALLELE_NUM=2;IMPACT=MODIFIER;STRAND=-1`** |  |

The same doubling appears for a `<BND>` with `SVLEN` and no `CHR2` (over its span), for `CHR2` with `POS2` alone, for `CHR2` with an `END2` that is not a coordinate, for a `<BND>` without `SVTYPE`, and for a `<BND>` with `MATEID`. A native single breakend (`.N`, `N.`) is held without the reference base and writes each row once, and a bracket record's own breakend and mate are two different alleles, written once each.

A second behaviour concerns `feature_truncation`: the rows carry none under `--no_stats`, the flag the reference protocol runs with. `OutputFactory.pm` line 2029 reads `$vfoa->{breakend}->{string}` before the terms are computed at line 2034. The read autovivifies an empty breakend (Perl creates the missing hash on first access), whose chromosome matches no feature's, so the breakend test of `feature_truncation` (`VariationEffect.pm` lines 357-358) fails. With statistics on, `Stats.pm` line 241 computes the terms first, `within_feature` falls back to the whole feature, and the same rows carry `feature_truncation`. vep-rs follows the reference protocol's rows.

### 13.4 Why Ensembl VEP writes that

`Parser/VCF.pm` builds a bracket mate only under `defined $info->{CHR2}` with `END2 || END` (lines 491-497) and prepends the reference base to the allele string unless it opens or closes with `.` (line 499), so the record's string is `N/<BND>`. `StructuralVariationFeature.pm` `_parse_breakends` (lines 1274-1283) finds no bracket in it. `StructuralVariationOverlap.pm` lines 64-77 then run `if (!@$breakends) { foreach (split(/\//, $vf->{allele_string})) { add_StructuralVariationOverlapAllele(...) } }`, numbering the parts from 1: the reference base `N` is allele 1 and `<BND>` allele 2. Neither allele carries a breakend string, so `OutputFactory.pm` line 2029 (`$vfoa->{breakend}->{string} || $svf->class_SO_term`) prints `chromosome_breakpoint` for both.

### 13.5 Why vep-rs's row is the correct one

A reference allele has no variant effect; a consequence row describes what an alternate allele does to a feature. The second row restates the first under a number that points at the reference base, and nothing in it is information about the variant. vep-rs holds the record as one single breakend, prints the class term and the terms of the span or the point once, and takes no `feature_truncation` for it, as Ensembl VEP's own ALT row does not.

### 13.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tbare_bnd\tN\t<BND>\t.\t.\tSVTYPE=BND\n' > bare.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/bare.vcf -o /work/bare.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i bare.vcf -o bare.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

### 13.7 How the comparator treats this class

The class is `bnd_bare_reference_allele_row` in the 116.2 registry of both comparators, and both set the reference-allele rows aside. Each comparator compares unique tuples of (Location, Allele, Feature, Feature_type, Consequence), and the reference-allele row is identical in every column to the ALT row, so the set-aside happens where the tuples are built. The duplicate row is dropped before the raw and the adjusted figure, neither of which ever carries a tuple for it. Each comparator reads the reference's raw rows to count the rows it set aside (`bnd_bare_reference_allele_rows` in the structural-variant report, `bnd_bare_reference_allele_row` with its `excluded` count in the SNP/indel report, both per file).

The class therefore moves no tuple between the raw and the adjusted figure; it names and counts what the construction excluded. None of the measured sets carries the shape: every gnomAD-SV `<BND>` on chromosome 21 names its mate through `CHR2` with `END` or `END2`, and the synthetic sets write their breakends in bracket and dot notation.
