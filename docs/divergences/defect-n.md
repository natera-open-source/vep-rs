# Defect N: a POS that is not a number is numified on the structural path

**The defect in one sentence.** On its structural path Ensembl VEP derives a record's start by adding one to the raw POS text (`$start++`), and Perl's increment numifies whatever the text is, so the record is annotated at a position its line never named.

Numifying is Perl's conversion of text to a number: the increment turns `.` into 1, `30_425_949` into 31 and `-5` into -4.

**Why this is a defect.** POS is a 1-based integer in the VCF specification; a position that is not a number names no base, and a negative one names none either. A `<DEL>` with `END=30425949` and POS `.` becomes a deletion of the chromosome from base 1, 3,254 transcript rows describing an event the line never asserted, while Ensembl VEP's own sequence path refuses the same POS (`Invalid start '.' or end '0' coordinate`).

**What vep-rs does instead.** vep-rs reads POS as an integer coordinate and refuses the line with a warning (`Invalid POS`), in every form.

**Where it is in Ensembl VEP's source** (`ensembl-io` and `ensembl-vep`):
- `ensembl-io`, `BaseVCF4.pm`, `get_start`, lines 255-264: for a record whose raw INFO carries `SVTYPE` or whose ALT carries `<`, `[` or `]`, line 263 runs `$start++` on the raw POS text. Perl's increment reads the text's leading digits and nothing else (`.`, an empty string, the hexadecimal `0x1D04A1D` and fullwidth digits are 0, `30_425_949` is 30) and adds one. A sequence record takes the other branch, which increments only when an ALT differs in length from REF, so an SNV keeps the text.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser.pm`, `validate_vf`, lines 475-488: finds a number at `start` and at `end` with `start` at most `end + 1` and keeps the record; the sequence record's text fails `looks_like_number` at line 475.

**What correcting it would change.** The record would be refused, as Ensembl VEP's own sequence path refuses it and as vep-rs refuses every form; class 21 would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 21. A non-numeric POS numified by the structural path's increment

**In short.** On the structural path Ensembl VEP adds one to the raw POS text, and Perl's increment turns any text into a number (`.` to 1, `30_425_949` to 31, `-5` to -4), so the record is annotated at a position its line never named. vep-rs refuses every such line with a warning and writes no row, because a position that is not a number names no base and a negative one names none either. Ensembl VEP's own sequence path also refuses `.`, an empty POS, `0x1D04A1D` and `30_425_949`, whereas it writes a negative POS as `21:-5`.

### 21.1 What a reader needs to know first

POS is a 1-based integer in the VCF specification. On the structural path Ensembl VEP derives a record's start by adding one to the raw POS text (`$start++`), and Perl's increment numifies whatever the text is: `.`, an empty column, `0x1D04A1D` and fullwidth digits numify to 0 and increment to 1. `30_425_949` reads its leading digits, 30, and increments to 31; `-5`, a number but no coordinate, increments to -4.

The numified start passes `looks_like_number`, the start-against-end test passes, and the record is annotated at a position its line never named: a `<DEL>` with `END=30425949` and POS `.` is a deletion of the chromosome from base 1. A sequence record keeps the text and is refused (`Invalid start '.' or end '0' coordinate`), as vep-rs refuses every form.

### 21.2 The record

```
21	.	pos_dot_del	N	<DEL>	.	.	SVTYPE=DEL;END=30425949
21	30_425_949	pos_underscore_del	N	<DEL>	.	.	SVTYPE=DEL;END=30425949
21	-5	pos_neg_del	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=-300
21	.	pos_dot_snv	A	G	.	.	.
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`; `--max_sv_size -1` for the first two records).

Once numified, the first two span more than 10,000,000 bases, so at the default `--max_sv_size` Ensembl VEP drops them as oversize (`variant size (30425948) is bigger than --max_sv_size (10000000)`) and annotates them under `--max_sv_size -1`.

### 21.3 What each engine writes

Ensembl VEP 116.2 under `--max_sv_size -1`, 3,254 rows for the first record, one per transcript the span from base 1 to 30425949 reaches, 3,253 of them `transcript_ablation`; two of them:

The first record, `pos_dot_del`, on ENST00000262354: one of Ensembl VEP's 3,254 rows, no vep-rs row.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `pos_dot_del` | no row |
| Location | `21:1-30425949` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000228592` |  |
| Feature | `ENST00000262354` |  |
| Feature_type | `Transcript` |  |
| Consequence | `transcript_ablation` |  |
| Extra | `ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1;OverlapBP=23932;OverlapPC=100.00` |  |

The first record on ENST00000282964: another of the 3,254 rows, no vep-rs row.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `pos_dot_del` | no row |
| Location | `21:1-30425949` |  |
| Allele | `deletion` |  |
| Gene | `ENSG00000229336` |  |
| Feature | `ENST00000282964` |  |
| Feature_type | `Transcript` |  |
| Consequence | `transcript_ablation` |  |
| Extra | `ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;OverlapBP=442;OverlapPC=100.00` |  |

The same 3,254 rows for the second record at `21:31-30425949`, and for the third, at the default cap, one row:

The third record, `pos_neg_del`: Ensembl VEP's one row, no vep-rs row.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `pos_neg_del` | no row |
| Location | `21:-4-295` |  |
| Allele | `deletion` |  |
| Gene | `-` |  |
| Feature | `-` |  |
| Feature_type | `-` |  |
| Consequence | `intergenic_variant` |  |
| Extra | `ALLELE_NUM=1;IMPACT=MODIFIER` |  |

An empty POS, `0x1D04A1D` and fullwidth digits write the first record's rows. The fourth record is skipped (the warnings file reads `Invalid start '.' or end '0' coordinate`), as is a sequence record whose POS is empty, `0x1D04A1D` or `30_425_949`. A negative POS on a sequence record passes `looks_like_number` and is written as `21:-5 G - intergenic_variant`.

vep-rs: no row for any of them, each line reported at warning level (`Invalid POS: '.'`, `Invalid POS: '30_425_949'`, `Invalid POS: '-5'`).

### 21.4 Why Ensembl VEP writes that

`BaseVCF4.pm` `get_start` (lines 255-264) takes the raw POS text and, for a record whose raw INFO carries `SVTYPE` or whose ALT carries `<`, `[` or `]`, runs `$start++` on it (line 263). Perl's increment numifies the text, reading its leading digits and nothing else (`.`, an empty string, `0x1D04A1D` and fullwidth digits are 0, `30_425_949` is 30), and adds one. `Parser.pm` `validate_vf` (lines 475-488) then finds a number at `start` and at `end` and `start` at most `end + 1`, and keeps the record. A sequence record takes the other branch of `get_start`, which increments only when an ALT differs in length from REF, so an SNV keeps the text, and `looks_like_number('.')` fails at line 475.

### 21.5 Why vep-rs's output is the correct one

A position that is not a number names no base, and a negative one names none either. The 3,254 rows of the first record describe the ablation of every transcript on the first 30 Mb of the chromosome, an event the line never asserted: its POS is the missing value. vep-rs reads POS as an integer coordinate and refuses the line with a warning, as Ensembl VEP's own sequence path refuses the same POS.

### 21.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t-5\tpos_neg_del\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-300\n' > pos.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/pos.vcf -o /work/pos.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i pos.vcf -o pos.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

Replacing `-5` with `.` and `SVLEN=-300` with `END=30425949`, and adding `--max_sv_size -1` to both commands, gives the 3,254 rows from Ensembl VEP and the refusal from vep-rs.

### 21.7 How the comparator treats this class

The class is `non_numeric_pos_numified` in the 116.2 registry of both comparators. A non-negative number here is digits with an optional fraction, exponent, leading `+` or surrounding blanks. The class sets aside the reference's rows of a record whose POS is anything else (`.`, an empty column, `0x1D04A1D`, fullwidth digits, `30_425_949`, `-5`), at Locations vep-rs wrote nothing for. A fractional POS is [class 18](defect-k.md#18-a-non-integer-end-svlen-or-pos-text-annotated-and-printed), and a POS above 2^62 - 1, or at it for a structural record whose first base is POS + 1, is a number and not this class. vep-rs refuses the line as past its coordinate bound and the reference's rows stay charged.

The structural-variant comparator reads the input and matches the reference's rows by the start the reference prints for the record, POS + 1 on the structural path and POS on the sequence path (`21:-5 G - intergenic_variant`, the one form a sequence record reaches). The SNP/indel comparator reads no input, so it sets aside only a reference-only Location whose first coordinate is negative: no POS that is a non-negative number prints one. A POS numified to 0 prints a start of 1, which a record at POS 0 also prints in both engines, and `30_425_949` prints 31, an ordinary start. Only the structural-variant comparator, which reads the input, recognises those forms.
