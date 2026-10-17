# Defect M: `*` as the only ALT is annotated from `SVTYPE`

**The defect in one sentence.** The ALT `*` names no allele of the record it stands in, yet when INFO carries a true `SVTYPE` Ensembl VEP annotates the record from INFO alone.

The `*` marks a base covered by an upstream deletion. Ensembl VEP routes the line to its structural path, where `*` falls outside the VCF 4.4 forms, so `SVTYPE` names the type. The record becomes a deletion over `END`, or a breakpoint at POS + 1 written under two alleles.

**Why this is a defect.** The specification defines `*` as the mark of a base deleted by an overlapping record, so the record's alleles are the reference and nothing else. A deletion of 299 bases, or a breakpoint, read from INFO alone describes an allele the record never states. Ensembl VEP's own sequence path agrees, writing no consequence for the record without `SVTYPE`.

**What vep-rs does instead.** vep-rs annotates nothing from INFO when the ALT column names no allele: no row in the default and tab formats, the VCF line carried without a CSQ value, and a JSON object with `"most_severe_consequence":"?"`. That is what Ensembl VEP writes for the same record without `SVTYPE`.

**Where it is in Ensembl VEP's source** (`ensembl-vep` and `ensembl-variation`):
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 237-245: the line is routed to the structural path, REF and ALT together (`N*`) not being a run of `ACGT` and `SVTYPE` being true.
- `ensembl-vep`, `Parser/VCF.pm`, line 469: the ALTs are joined into the allele string `*`; lines 472-475 test it against the VCF 4.4 forms, `^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)`, which `*` fails, so `$type = $info->{SVTYPE}` and `get_SO_term` names a deletion or a breakpoint.
- `ensembl-vep`, `Parser/VCF.pm`, line 499: for the breakpoint the reference base is prepended, giving `N/*`.
- `ensembl-variation`, `StructuralVariationOverlap.pm`, lines 64-77: with no bracket parsed from `N/*`, one allele is built per slash-separated part, `N` and `*`, numbered 1 and 2 (the mechanism of [Defect F](defect-f.md)), so every row of the breakpoint is written twice.

**What correcting it would change.** The record would write nothing, as it does without `SVTYPE` and as vep-rs writes; [class 20](#20--alone-under-a-true-svtype-annotated-as-the-types-class) would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 20. `*` alone under a true `SVTYPE` annotated as the type's class

**In short.** The ALT `*` marks a base covered by an upstream deletion and names no allele of the record it stands in. Under a true `SVTYPE`, Ensembl VEP annotates the record from INFO alone, a deletion over `END` or a breakpoint at POS + 1 written under two alleles; without `SVTYPE` it writes nothing. vep-rs annotates nothing from INFO when the ALT column names no allele, and carries the line as Ensembl VEP carries the `SVTYPE`-less record.

### 20.1 What a reader needs to know first

The ALT `*` is the VCF 4.4 marker for a base covered by an upstream deletion (section 1.6.1.5 of the specification). It says that the base at POS is missing from the haplotype that carries the overlapping deletion, and it names no allele of the record it stands in. A record whose only ALT is `*` therefore has no alternate allele. Without `SVTYPE` Ensembl VEP writes nothing for it: the line is carried through the VCF output without consequences, and the JSON object is written with `"most_severe_consequence":"?"` and no consequence array. vep-rs writes the same.

With a true `SVTYPE` in INFO Ensembl VEP routes the line to the structural path, where `*` falls outside the VCF 4.4 forms. `SVTYPE` therefore names the type, and the record is annotated from INFO alone: `SVTYPE=DEL;END` becomes a deletion over `END`, `SVTYPE=BND` a breakpoint at POS + 1 written under two alleles.

### 20.2 The record

```
21	30425650	star_svtype_del	N	*	.	.	SVTYPE=DEL;END=30425949
21	30425650	star_svtype_del_refa	A	*	.	.	SVTYPE=DEL;END=30425949
21	30425650	star_svtype_bnd	N	*	.	.	SVTYPE=BND
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

### 20.3 What each engine writes

Ensembl VEP 116.2, for the first record (the second, with REF `A`, writes the same two rows under its own name), two rows, one per feature, and for the third, four rows, every row twice under `ALLELE_NUM` 1 and 2. vep-rs writes no row for any of the three; one table per row.

`star_svtype_del` on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `star_svtype_del` | no row |
| Location | `21:30425651-30425949` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4271;STRAND=1` |  |

`star_svtype_del` on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `star_svtype_del` | no row |
| Location | `21:30425651-30425949` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `feature_truncation,coding_sequence_variant,5_prime_UTR_variant` |  |
| Extra | `ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1;OverlapBP=276;OverlapPC=50.46` |  |

`star_svtype_bnd` on ENST00000334068 under `ALLELE_NUM=1`:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `star_svtype_bnd` | no row |
| Location | `21:30425651` |  |
| Allele | `chromosome_breakpoint` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4569;STRAND=1` |  |

`star_svtype_bnd` on ENST00000334068 under `ALLELE_NUM=2`:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `star_svtype_bnd` | no row |
| Location | `21:30425651` |  |
| Allele | `chromosome_breakpoint` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `ALLELE_NUM=2;IMPACT=MODIFIER;DISTANCE=4569;STRAND=1` |  |

`star_svtype_bnd` on ENST00000390690 under `ALLELE_NUM=1`:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `star_svtype_bnd` | no row |
| Location | `21:30425651` |  |
| Allele | `chromosome_breakpoint` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `coding_sequence_variant` |  |
| Extra | `ALLELE_NUM=1;IMPACT=MODIFIER;STRAND=-1` |  |

`star_svtype_bnd` on ENST00000390690 under `ALLELE_NUM=2`:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `star_svtype_bnd` | no row |
| Location | `21:30425651` |  |
| Allele | `chromosome_breakpoint` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `coding_sequence_variant` |  |
| Extra | `ALLELE_NUM=2;IMPACT=MODIFIER;STRAND=-1` |  |

Ensembl VEP's warnings file carries `Use of uninitialized value $region1` warnings (`Utils/VariationEffect.pm` line 120). The VCF output carries each line with a CSQ value, and the JSON object carries the consequences with `"allele_string":"*"` for the deletions and `"N/*"` for the breakpoint.

vep-rs: no row for any of the three in the default and tab formats, no warning. The VCF output carries each line without a CSQ value. The JSON output writes an object with the record's `id`, `seq_region_name`, `start` (POS), `end` (POS plus the length of REF minus one), `strand`, `allele_string` (REF and every ALT, `N/*`), `assembly_name`, `input` and `"most_severe_consequence":"?"`. That is what Ensembl VEP writes for the same record without `SVTYPE`.

### 20.4 Why Ensembl VEP writes that

`Parser/VCF.pm` lines 237-245 route the line to the structural path: REF and ALT together (`N*`) are not a run of `ACGT`, and `SVTYPE` is true. Line 469 joins the ALTs into the allele string `*`; lines 472-475 test it against the VCF 4.4 forms, `^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)`, which `*` fails, so `$type = $info->{SVTYPE}` and `get_SO_term` names a deletion or a breakpoint. For the deletion the allele string stays `*` and one allele is built from it. For the breakpoint line 499 prepends the reference base (`N/*`), and `_parse_breakends` finds no bracket in it. `StructuralVariationOverlap.pm` lines 64-77 then build one allele per slash-separated part, `N` and `*`, numbered 1 and 2 (the mechanism of [Defect F](defect-f.md)), so every row is written twice.

### 20.5 Why vep-rs's output is the correct one

`*` names no allele of this record: the specification defines it as the mark of a base deleted by an overlapping record, so the record's alleles are the reference and nothing else. A deletion of 299 bases, or a breakpoint, read from INFO alone describes an allele the record never states. Ensembl VEP's own sequence path agrees, writing no consequence for the record without `SVTYPE`. vep-rs annotates nothing from INFO when the ALT column names no allele, and carries the line as Ensembl VEP carries the `SVTYPE`-less record.

### 20.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tstar_svtype_del\tN\t*\t.\t.\tSVTYPE=DEL;END=30425949\n' > star.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/star.vcf -o /work/star.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i star.vcf -o star.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

Replacing `SVTYPE=DEL;END=30425949` with `.` gives no row from either engine.

### 20.7 How the comparator treats this class

The class is `star_alone_svtype_annotated` in the 116.2 registry. The structural-variant comparator reads the input and sets aside the reference's rows at the Location of a record whose ALTs are all `*` and whose `SVTYPE` Perl reads as true. The SNP/indel comparator reads no input and cannot tell such a row from a record with an allele, so there the class counts nothing and the rows stay charged.
