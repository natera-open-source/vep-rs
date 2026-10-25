# Defect J: the missing ALT value is annotated as an allele

**The defect in one sentence.** A record whose ALT column is the missing value `.` has no alternate allele, yet when INFO carries `SVTYPE` Ensembl VEP annotates it from INFO alone.

Ensembl VEP routes the line to its structural path, reads the `.` as a type and annotates the record: a `chromosome_breakpoint` at POS + 1 for `SVTYPE=BND`, a deletion over `END` for `SVTYPE=DEL`. A `.` listed after another ALT (`C G,.`) is read on the sequence path as an allele of its own, written as a row whose Allele is `.` with the next `ALLELE_NUM`.

**Why this is a defect.** `.` is the VCF's missing value, not an allele, and the rows describe an allele the record never states. Ensembl VEP's own sequence path agrees: without `SVTYPE` it writes nothing for the record.

**What vep-rs does instead.** vep-rs writes no row for the record in either case and no warning: a record whose ALT column is `.` names no allele and is passed over, as Ensembl VEP's own sequence path passes it over. On the sequence path a `.` after another ALT takes no row and no allele number, and the other ALTs are annotated as written, so `C G,.` without a true INFO/SVTYPE writes the `G` row alone. Under a true INFO/SVTYPE the `.` part writes a row of the record's class, as Ensembl VEP's does.

**Where it is in Ensembl VEP's source** (`ensembl-io` and `ensembl-vep`):
- `ensembl-io`, `BaseVCF4.pm`, `get_alternatives`, lines 486-490: returns `['.']` for the missing value.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 237-245: the line is routed to the structural path when `SVTYPE` is set, since REF and ALT together (`N.`) are not a run of `ACGT`. Lines 472-475 let `SVTYPE` name the type when the ALT is outside the VCF 4.4 forms.
- `ensembl-vep`, `Parser.pm`, lines 710-711: `get_SO_term` reads a breakend from `^\.`.

**What correcting it would change.** The record would write nothing, as it does without `SVTYPE` and as vep-rs writes; [class 17](#17-the-missing-alt-value-annotated-as-an-allele) would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 17. The missing ALT value annotated as an allele

**In short.** A record whose ALT column is `.` has no alternate allele. With `SVTYPE` in INFO, Ensembl VEP annotates the record from INFO alone, a `chromosome_breakpoint` for `SVTYPE=BND` or a deletion over `END` for `SVTYPE=DEL`; without `SVTYPE` it writes nothing. vep-rs writes no row and no warning in either case. On the sequence path a `.` listed after another ALT gets a row of its own from Ensembl VEP and none from vep-rs, which annotates the other ALTs as written.

### 17.1 What a reader needs to know first

A VCF record whose ALT column is `.` has no alternate allele. Without `SVTYPE` Ensembl VEP writes nothing for it, as vep-rs does. With `SVTYPE` in INFO the line is routed to the structural path, where the `.` is read as a type. `get_SO_term` reads a breakend from a leading dot, and `SVTYPE` names the type when the ALT is outside the VCF 4.4 forms. So `SVTYPE=BND` yields a `chromosome_breakpoint` at POS + 1 and `SVTYPE=DEL;END` a deletion over END, each with rows.

A `.` listed after another ALT (`21 30289183 x C G,.`) stays on the sequence path, where Ensembl VEP's overlap builder makes an allele of it. Beside the `G` row it writes a row whose Allele is `.`, `coding_sequence_variant` with `ALLELE_NUM=2` and no amino-acid or codon change, under the allele string `C/G/.`. vep-rs writes the `G` row alone, byte for byte Ensembl VEP's, and no row for the `.`.

### 17.2 The record

```
21	30425650	dot_svtype_bnd	N	.	.	.	SVTYPE=BND
21	30425650	dot_svtype_del_end	N	.	.	.	SVTYPE=DEL;END=30425949
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

### 17.3 What each engine writes

Ensembl VEP 116.2 writes two rows for each record, one per feature, and vep-rs no row for either, with no warning on standard error and no warnings file. One table per row.

The first record, `dot_svtype_bnd`, on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `dot_svtype_bnd` | no row |
| Location | `21:30425651` |  |
| Allele | `chromosome_breakpoint` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |

The first record, `dot_svtype_bnd`, on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `dot_svtype_bnd` | no row |
| Location | `21:30425651` |  |
| Allele | `chromosome_breakpoint` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `coding_sequence_variant` |  |

The second record, `dot_svtype_del_end`, on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `dot_svtype_del_end` | no row |
| Location | `21:30425651-30425949` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |

The second record, `dot_svtype_del_end`, on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `dot_svtype_del_end` | no row |
| Location | `21:30425651-30425949` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `feature_truncation,coding_sequence_variant,5_prime_UTR_variant` |  |

For `C G,.` both engines write the `G` row; Ensembl VEP adds the Allele `.` row and vep-rs does not.

### 17.4 Why Ensembl VEP writes that

`BaseVCF4.pm` `get_alternatives` (lines 486-490) returns `['.']` for the missing value. `Parser/VCF.pm` lines 237-245 route the line to the structural path when `SVTYPE` is set, since REF and ALT together (`N.`) are not a run of `ACGT`. `get_SO_term` (`Parser.pm` lines 710-711) reads a breakend from `^\.`, and lines 472-475 of `Parser/VCF.pm` let `SVTYPE` name the type when the ALT is outside the VCF 4.4 forms.

### 17.5 Why vep-rs's output is the correct one

`.` is the VCF's missing value, not an allele, and Ensembl VEP's own sequence path writes nothing for it when it stands alone. vep-rs passes the lone `.` over in every case, with no row and no warning, and on the sequence path gives a `.` after another ALT no row, since no allele is stated for it to annotate.

### 17.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tdot_svtype_bnd\tN\t.\t.\t.\tSVTYPE=BND\n' > missing.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/missing.vcf -o /work/missing.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i missing.vcf -o missing.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

### 17.7 How the comparator treats this class

The class is `missing_alt_annotated` in the 116.2 registry. The structural-variant comparator reads the input and sets aside the reference's rows at the Location of a record whose ALT is `.` and whose `SVTYPE` Perl reads as true. The SNP/indel comparator, which opens the input only to recognise gVCF variant sites under `--input-dir` or `--input-vcf`, judges this class from the output alone. There it cannot tell such a row from a record with an allele, so the class counts nothing and the rows stay charged.
