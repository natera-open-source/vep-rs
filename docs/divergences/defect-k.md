# Defect K: a non-integer coordinate text is annotated and printed

**The defect in one sentence.** Ensembl VEP annotates a record whose `INFO/END`, `SVLEN` or POS is a fraction, an exponent or `NaN` with the text in its Location and a fractional `DISTANCE`, checking only that the start and end look like numbers.

It takes `INFO/END`, `SVLEN` and POS as the text they are written in.

**Why this is a defect.** `END`, `SVLEN` and POS are integers in the VCF specification, and a coordinate is a base position: a half base or a floating-point notation names none. A Location such as `21:30425651-30425949.5` or `21:30425651-NaN` describes no interval of the genome, and the `NaN` record's JSON object carries an end that is not a JSON value.

**What vep-rs does instead.** vep-rs refuses a non-integer `END`, which names no base, and reports the line. It reads `SVLEN` as Perl numifies a length (`300.7` is 300 bases, `NaN` is 0 and gives a two-base interval) and a POS by the integer its text denotes, so those records are annotated at the number where Ensembl VEP prints the text.

**Where it is in Ensembl VEP's source** (`ensembl-io` and `ensembl-vep`):
- `ensembl-io`, `BaseVCF4.pm`, `get_end`, lines 307-308: `INFO/END` is returned as the string it is; lines 302-305: a `SVLEN` end is `$start + abs($svlen) - 1`, so `300.7` gives the end `30425950.7` and `NaN` the end `NaN`.
- `ensembl-io`, `BaseVCF4.pm`, `get_start`, lines 255-275: a sequence record's POS is kept as its text, the increment applying only to a structural record's or an indel's start.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser.pm`, `validate_vf`, line 475: the values pass `looks_like_number`, which admits a fraction, an exponent and `NaN`; every comparison with `NaN` is false, so line 483's start-against-end test skips nothing.

**What correcting it would change.** The record would be refused, as vep-rs refuses a non-integer `END`, or read at the number its text denotes, as vep-rs reads `SVLEN` and POS; [class 18](#18-a-non-integer-end-svlen-or-pos-text-annotated-and-printed) would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 18. A non-integer `END`, `SVLEN` or POS text annotated and printed

**In short.** Ensembl VEP annotates a record whose `END`, `SVLEN` or POS text is a fraction, an exponent or `NaN`, and prints that text in its Location. vep-rs refuses a non-integer `END`, because a coordinate is a base position and a half base or a floating-point notation names none. `SVLEN` is a length, not a coordinate, so vep-rs reads it as Perl numifies a length (`300.7` is 300 bases, `NaN` is 0); a POS whose text denotes an integer is read as that integer, as Ensembl VEP's own structural path reads it.

### 18.1 What a reader needs to know first

`INFO/END` and `SVLEN` are integers in the VCF specification, and so is POS. Ensembl VEP takes each value as written, checks only that the start and end look like numbers, and annotates the record with the text in its Location and a fractional `DISTANCE`. The text may be an `END` of `30425949.5`, a `SVLEN` of `300.7` (end `30425950.7`) or `NaN` (end `NaN`), or on the sequence path a POS of `30425949.0`, which denotes an integer but is printed as its text.

vep-rs refuses a non-integer `END`, which names no base, and reports the line. It reads `SVLEN` as Perl numifies a length (`300.7` is 300 bases, `NaN` is 0 and gives a two-base interval) and a POS by the integer its text denotes, so the record is annotated at the number where Ensembl VEP prints the text.

### 18.2 The record

```
21	30425650	del_end_float	N	<DEL>	.	.	SVTYPE=DEL;END=30425949.5
21	30425650	del_end_exp	N	<DEL>	.	.	SVTYPE=DEL;END=3.0425949e7
21	30425650	del_svlen_float	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=300.7
21	30425650	del_svlen_nan	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=NaN
21	30425949.0	pos_float_snv	A	G	.	.	.
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

### 18.3 What each engine writes

Ensembl VEP 116.2 writes two rows for each of the two `END` records, one per feature, and vep-rs no row for either, each line reported at warning level (`Invalid END: '30425949.5'`). One table per row.

`del_end_float` on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `del_end_float` | no row |
| Location | `21:30425651-30425949.5` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `DISTANCE=4270.5` |  |

`del_end_float` on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `del_end_float` | no row |
| Location | `21:30425651-30425949.5` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `feature_truncation,coding_sequence_variant,5_prime_UTR_variant` |  |

`del_end_exp` on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `del_end_exp` | no row |
| Location | `21:30425651-3.0425949e7` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `DISTANCE=4271` |  |

`del_end_exp` on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `del_end_exp` | no row |
| Location | `21:30425651-3.0425949e7` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `feature_truncation,coding_sequence_variant,5_prime_UTR_variant` |  |

For the `SVLEN` and POS records, Ensembl VEP 116.2 prints the text too (the `NaN` record's JSON object carries `"start":30425651,"end":nan`, which is not a JSON value), five rows. For these three records vep-rs writes the rows of the same record at the number: `21:30425651-30425950` for `SVLEN=300.7`, the two-base interval `21:30425650-30425651` for `SVLEN=NaN` (two transcript rows, `coding_sequence_variant` on ENST00000390690, where the reference wrote one intergenic row), and `21:30425949` for the POS. One table per Ensembl VEP row; the vep-rs column points to that sentence.

`del_svlen_float` on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `del_svlen_float` | see above |
| Location | `21:30425651-30425950.7` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `DISTANCE=4269.30000000075` |  |

`del_svlen_float` on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `del_svlen_float` | see above |
| Location | `21:30425651-30425950.7` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `feature_truncation,coding_sequence_variant,5_prime_UTR_variant` |  |

`del_svlen_nan`, the one intergenic row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `del_svlen_nan` | see above |
| Location | `21:30425651-NaN` |  |
| Allele | `deletion` |  |
| Gene | `-` |  |
| Feature | `-` |  |
| Feature_type | `-` |  |
| Consequence | `intergenic_variant` |  |

`pos_float_snv` on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `pos_float_snv` | see above |
| Location | `21:30425949.0` |  |
| Allele | `G` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `DISTANCE=4271` |  |

`pos_float_snv` on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `pos_float_snv` | see above |
| Location | `21:30425949.0` |  |
| Allele | `G` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `DISTANCE=23` |  |

### 18.4 Why Ensembl VEP writes that

`BaseVCF4.pm` `get_end` (lines 307-308) returns `INFO/END` as the string it is, and `Parser.pm` `validate_vf` (line 475) accepts it through `looks_like_number`, which admits a fraction, an exponent and `NaN`. The end then takes part in every arithmetic as a number and is printed as the text. A `SVLEN` reaches the same place through lines 302-305, `$start + abs($svlen) - 1`, so `300.7` gives the end `30425950.7` and `NaN` the end `NaN` (every comparison with it is false, so line 483's start-against-end test skips nothing).

On the sequence path POS is kept as its text by `get_start` (lines 255-275, which increment only a structural record's or an indel's start), passes `looks_like_number` and is printed. On the structural path `$start++` numifies it ([Defect N](defect-n.md)). An `END` written with leading zeros or a `+` sign (`030425949`, `+30425949`) is the same shape without the arithmetic: the record is annotated correctly and the Location carries the text, which both comparators read as the number it names.

### 18.5 Why vep-rs's output is the correct one

A coordinate is a base position; a half base or a floating-point notation names none. vep-rs reads `END` as an integer coordinate and refuses the record otherwise, as it refuses a negative one. `SVLEN` is a length, not a coordinate, and vep-rs reads it as Perl numifies a length, so `300.7` is 300 bases and `NaN` is 0, a two-base interval, where Ensembl VEP prints an end that is not a position. A POS whose text denotes an integer (`30425949.0`, `3e7`) is read as that integer, which is what Ensembl VEP's own structural path makes of it.

### 18.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tdel_end_float\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425949.5\n' > float.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/float.vcf -o /work/float.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i float.vcf -o float.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

### 18.7 How the comparator treats this class

The class is `non_integer_end_annotated` in the 116.2 registry of both comparators, set aside: the reference's rows whose Location carries a coordinate that is not an integer. The structural-variant comparator, which reads the input, also sets aside vep-rs's rows of the same record where vep-rs read the value as a number, so that neither side is charged and the skipped-record cross-check does not stop on them. Those are the rows at the numified end for a fractional `SVLEN`, `21:30425651-30425950`, at the two-base interval for `SVLEN=NaN`, `21:30425650-30425651`, and at the integer for a POS text, `21:30425949`. The SNP/indel comparator reads no input, locates no vep-rs row of the record, and sets aside the reference's text rows alone.

A Location whose coordinates are integers written with a sign or leading zeros is not the class: both comparators read such a coordinate as its number before comparing, so those rows agree with vep-rs's.
