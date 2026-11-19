# Defect N: a POS that is not a number is numified on the structural path

**The defect in one sentence.** Ensembl VEP derives a structural record's start, and an anchor-chopped indel's, by adding one to the raw POS text (`$start++`), and Perl's increment numifies whatever the text is, so the record is annotated at a position its line never named. vep-rs refuses every such line with a warning and writes no row.

Numifying is Perl's conversion of text to a number: the increment turns `.` into 1, `30_425_949` into 31 and `-5` into -4. An anchor-chopped indel is a VCF indel whose REF and ALT share their first base, which the parser removes, moving the start one base on.

**Why this is a defect.** POS is a 1-based integer in the VCF specification; a position that is not a number names no base, and a negative one names none either. A `<DEL>` with `END=30425949` and POS `.` becomes a deletion of the chromosome from base 1, 3,254 transcript rows describing an event the line never asserted. Ensembl VEP's own sequence path refuses the same POS on a substitution (`Invalid start '.' or end '0' coordinate`) and numifies it on an anchor-chopped indel, `21 . x GA G` annotated at `21:1`.

**What vep-rs does instead.** vep-rs reads POS as an integer coordinate and refuses the line with a warning in every form. A POS that is not a number gets the sequence path's own text, `Invalid start '<POS>' or end '<end>' coordinate`; a negative or fractional one gets `Invalid POS: '<POS>'`.

**Where it is in Ensembl VEP's source** (`ensembl-io` and `ensembl-vep`):
- `ensembl-io`, `BaseVCF4.pm`, `get_start`, lines 255-264: for a record whose raw INFO carries `SVTYPE` or whose ALT carries `<`, `[` or `]`, line 263 runs `$start++` on the raw POS text. Perl's increment reads the text's leading digits and nothing else (`.`, an empty string, the hexadecimal `0x1D04A1D` and fullwidth digits are 0, `30_425_949` is 30) and adds one. A sequence record never reaches `get_start`: its start is the raw text itself (`Parser/VCF.pm` lines 248-254, `get_raw_start`).
- `ensembl-vep`, `Parser/VCF.pm`, lines 306-310 and 328-334: an indel whose REF and every ALT share their first base loses that base, and line 310 (several ALTs) or line 334 (one ALT) runs `$start++` on the raw POS text. The numification is the same: `.` to 1, `30_425_949` to 31, `-5` to -4, `12abc` to 13. A text of letters alone is incremented as a string instead (`abc` to `abd`), and the coordinate test refuses it quoting the incremented text.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/Parser.pm`, `validate_vf`, lines 475-488: finds a number at `start` and at `end` with `start` at most `end + 1` and keeps the record; a substitution's text fails `looks_like_number` at line 475, and a chopped indel's incremented number passes it.

**What correcting it would change.** The record would be refused, as Ensembl VEP's own sequence path refuses a substitution's and as vep-rs refuses every form; class 21 would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 21. A non-numeric POS numified by the structural path's increment

**In short.** On the structural path Ensembl VEP adds one to the raw POS text, and Perl's increment turns any text into a number (`.` to 1, `30_425_949` to 31, `-5` to -4), so the record is annotated at a position its line never named. vep-rs refuses every such line with a warning and writes no row, because a position that is not a number names no base and a negative one names none either. Ensembl VEP's own sequence path also refuses `.`, an empty POS, `0x1D04A1D` and `30_425_949` on a substitution, writes a negative POS as `21:-5`, and increments an anchor-chopped indel's POS as the structural path does.

### 21.1 What a reader needs to know first

POS is a 1-based integer in the VCF specification. On the structural path Ensembl VEP derives a record's start by adding one to the raw POS text (`$start++`), and Perl's increment numifies whatever the text is: `.`, an empty column, `0x1D04A1D` and fullwidth digits numify to 0 and increment to 1. `30_425_949` reads its leading digits, 30, and increments to 31; `-5`, a number but no coordinate, increments to -4.

The numified start passes `looks_like_number`, the start-against-end test passes, and the record is annotated at a position its line never named: a `<DEL>` with `END=30425949` and POS `.` is a deletion of the chromosome from base 1. A substitution keeps the text and is refused (`Invalid start '.' or end '0' coordinate`). An indel whose REF and ALT share their first base loses it, and the sequence path's own `$start++` (`Parser/VCF.pm` lines 310 and 334) numifies the text the same way, so `21 . x GA G` is annotated at `21:1`. vep-rs refuses every form.

### 21.2 The record

```
21	.	pos_dot_del	N	<DEL>	.	.	SVTYPE=DEL;END=30425949
21	30_425_949	pos_underscore_del	N	<DEL>	.	.	SVTYPE=DEL;END=30425949
21	-5	pos_neg_del	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=-300
21	.	pos_dot_snv	A	G	.	.	.
21	abc	bad_indel	GA	G	.	.	.
21	12abc	num_prefix	GA	G	.	.	.
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`; `--max_sv_size -1` for the first two records). Both engines take `--format vcf` as well: neither detects the format from a first data line whose POS is not a number.

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

An empty POS, `0x1D04A1D` and fullwidth digits write the first record's rows. The fourth record is skipped (the warnings file reads `Invalid start '.' or end '0' coordinate`), as is a substitution whose POS is empty, `0x1D04A1D` or `30_425_949`. A negative POS on a substitution passes `looks_like_number` and is written as `21:-5 G - intergenic_variant`.

The fifth record, `bad_indel`, is skipped with the incremented text quoted: `Invalid start 'abd' or end '1' coordinate`. The sixth, `num_prefix`, is annotated: its POS numifies to 12 and the chop moves it to 13, so Ensembl VEP writes the deletion `A/-` there.

The sixth record, `num_prefix`: Ensembl VEP's one row, no vep-rs row.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `num_prefix` | no row |
| Location | `21:13` |  |
| Allele | `-` |  |
| Gene | `-` |  |
| Feature | `-` |  |
| Feature_type | `-` |  |
| Consequence | `intergenic_variant` |  |
| Extra | `ALLELE_NUM=1;IMPACT=MODIFIER` |  |

`21 12abc x ACG A` is written the same way at `21:13-14` as `CG/-`. The chop reaches the page's other spellings too: `21 . x GA G` is written at `21:1`, `21 30_425_949 x GA G` at `21:31` and `21 -5 x GA G` at `21:-4`, each `- intergenic_variant`, at the starts the structural records take.

vep-rs: no row for any of them, each line reported at warning level with the POS as written. The six descriptions, in record order:

```
Invalid start '.' or end '0' coordinate
Invalid start '30_425_949' or end '30' coordinate
Invalid POS: '-5'
Invalid start '.' or end '0' coordinate
Invalid start 'abc' or end '1' coordinate
Invalid start '12abc' or end '13' coordinate
```

### 21.4 Why Ensembl VEP writes that

`BaseVCF4.pm` `get_start` (lines 255-264) takes the raw POS text and, for a record whose raw INFO carries `SVTYPE` or whose ALT carries `<`, `[` or `]`, runs `$start++` on it (line 263). Perl's increment numifies the text, reading its leading digits and nothing else (`.`, an empty string, `0x1D04A1D` and fullwidth digits are 0, `30_425_949` is 30), and adds one. `Parser.pm` `validate_vf` (lines 475-488) then finds a number at `start` and at `end` and `start` at most `end + 1`, and keeps the record.

A sequence record's start is the raw text (`Parser/VCF.pm` lines 248-254, `get_raw_start`), so for a substitution `looks_like_number('.')` fails at line 475. An indel whose REF and ALT share their first base loses it first, and `$start++` at `Parser/VCF.pm` line 334 (line 310 for several ALTs) numifies the text as the structural path does. The record then passes the test and is annotated at the incremented position. A text of letters alone, `abc`, is incremented as a string to `abd`, which the test refuses.

### 21.5 Why vep-rs's output is the correct one

A position that is not a number names no base, and a negative one names none either. The 3,254 rows of the first record describe the ablation of every transcript on the first 30 Mb of the chromosome, an event the line never asserted: its POS is the missing value. vep-rs reads POS as an integer coordinate and refuses the line with a warning, as Ensembl VEP's own sequence path refuses the same POS on a substitution.

### 21.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t-5\tpos_neg_del\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-300\n' > pos.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/pos.vcf -o /work/pos.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --format vcf --allele_number --no_stats
vep -i pos.vcf -o pos.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --format vcf --no_stats
```

`--format vcf` is needed because neither engine detects the format from a first data line whose POS is not a number. Replacing `-5` with `.` and `SVLEN=-300` with `END=30425949`, and adding `--max_sv_size -1` to both commands, gives the 3,254 rows from Ensembl VEP and the refusal from vep-rs. Replacing the record with `21 12abc num_prefix GA G . . .` gives the `21:13` row from Ensembl VEP and the refusal from vep-rs.

### 21.7 How the comparator treats this class

The class is `non_numeric_pos_numified` in the 116.2 registry of both comparators. A non-negative number here is digits with an optional fraction, exponent, leading `+` or surrounding blanks. The class sets aside the reference's rows of a record whose POS is anything else (`.`, an empty column, `0x1D04A1D`, fullwidth digits, `30_425_949`, `-5`), at Locations vep-rs wrote nothing for. A fractional POS is [class 18](defect-k.md#18-a-non-integer-end-svlen-or-pos-text-annotated-and-printed), and a POS above 2^62 - 1, or at it for a structural record whose first base is POS + 1, is a number and not this class. vep-rs refuses the line as past its coordinate bound and the reference's rows stay charged.

The structural-variant comparator reads the input and matches the reference's rows by the start the reference prints for the record, POS + 1 on the structural path and POS on the sequence path (`21:-5 G - intergenic_variant`, the form a substitution reaches). An anchor-chopped indel prints the incremented start instead, which that match does not reach, so its rows stay charged there.

The SNP/indel comparator, which opens the input only to recognise gVCF variant sites under `--input-dir` or `--input-vcf`, judges this class from the output alone. It therefore sets aside only a reference-only Location whose first coordinate is negative: no POS that is a non-negative number prints one. A POS numified to 0 prints a start of 1, which a record at POS 0 also prints in both engines, and `30_425_949` prints 31, an ordinary start. Only the structural-variant comparator, which reads the input, recognises those forms.

## Two more numified texts: a `--custom` record limit and a promoter bound

Defect N's mechanism, Perl reading a text that is not a number as one, reaches two more places outside the coordinate columns. Each is a documented class: the section names the rows a comparison sets aside, and no committed dataset has its shape, so the comparison tools carry no rule for it, as [Intended divergences from Ensembl VEP](../intended-divergences.md) lists it. Each carries its tier under that page's rule: tier 1 is a shape the chromosome 21 population cuts carry with no flag, tier 2 one that needs a flag or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts are over the cuts' 15,646,481 records, the 111 golden corpora and their 392 flag sets.

## 26. A `--custom` record limit that is not a whole number

**The defect in one sentence.** Ensembl VEP compares a `--custom` source's running record count with its `num_records` text as a number, so a limit of `1.5` lists two records and never writes the `...` marker, and a limit of `-1` lists nothing and drops the source's fields.

**Why this is a defect.** The option's own documentation calls the value the maximum number of records to show, and the code marks a list it cut with a final `...` entry. Under `1.5` the count passes `0 < 1.5` and `1 < 1.5`, two records are listed, the loop ends at the third, and the equality test that writes the marker never holds. The list is cut and unmarked, and it exceeds the stated maximum. Under `-1` the loop ends before the first record, so the row carries no field of the source although records matched.

**What vep-rs does instead.** vep-rs reads the value as the count it denotes, its integer part and never below 0, lists at most that many records and writes the `...` marker when more matched. So `1.5` lists one record and the marker, `-1` and `0` the marker alone, and the warnings file names the value and its reading once per source ([Custom annotation](../custom-annotation.md#records-and-statistics-per-input-record)).

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, line 127, the POD: `maximum number of records to show (50 by default)`; line 167 stores the text as given.
- The same file, line 337 ends the loop once the count exceeds the value; lines 436-439 keep a record while the count is below it and write the `...` entry when the count equals it.

**What correcting it would change.** A fractional or negative limit would be read as a count or refused, and a cut list would end in `...`, as vep-rs's does.

**The class it produces.** `custom_non_integer_num_records`: the `<short_name>` fields of a `--custom` source given a `num_records` that is not a non-negative integer, where the reference lists one record more than the candidate and no `...` entry, or lists nothing where the candidate writes `...`. No row on the committed corpora or their comparison cases: the only `num_records` any of them names is `num_records=1`, in one flag set of the GRCh38 custom corpus.

**How often the shape occurs.** Tier 3: the shape is a `num_records` value that is neither a non-negative integer nor `all`, which none of the 392 flag sets passes; 14 of them run a `--custom` source.

### 26.1 An insertion overlapping two sites under `num_records=1.5`

**In short.** The record `ins_padded`, `21 10542414 A AT` of the GRCh38 custom corpus, overlaps two records of `sites.vcf.gz`, at 21:10542414 and 21:10542415. Under `num_records=1` both engines write `sites=21:10542415-10542414,...` with the first record's `AF`, `CLASS` and `FILTER`. Under `num_records=1.5` vep-rs writes the same, one record and the marker; Ensembl VEP writes both records and no marker, `sites=21:10542415-10542414,21:10542415-10542415` with both records' values, the output `num_records=2` gives. Under `-1` vep-rs writes `sites=...` and Ensembl VEP no `sites` field; under `0` both write `sites=...`.

## 27. A promoter whose `extended_start` is not a number is widened to position 0

**The defect in one sentence.** Under `--extended_promoters`, Ensembl VEP compares a promoter's `extended_start` attribute with the record's start as numbers, and Perl reads text that is not a number as 0, so a value such as `x` passes the test and becomes the promoter's start.

**Why this is a defect.** The widening's own comment calls the extension additive and the attributes the promoter's extended bounds; a bound that names no base can widen nothing. Perl numifies `x` to 0 under a runtime warning and `0 < start` holds, so the record's start becomes `x`, which every later numeric use reads as 0. Every variant from the chromosome's first base to the promoter's end then gains the promoter's `regulatory_region_variant` row, as far as the buffer's regions reach. A non-numeric `extended_end` changes nothing, since `0 > end` never holds.

**What vep-rs does instead.** vep-rs reads each bound as a whole number and leaves the promoter's own coordinate in place when a bound is not one, so the promoter keeps its core start while a valid other bound still widens it. The warnings file names the record and the attribute once ([Regulatory features from GFF files](../regulatory.md#regulatory-features-from-gff-files)).

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/RegFeat.pm`, lines 360-364, `_record_to_regfeat`: `$start = $xs if defined($xs) && $xs < $start` on the attribute's raw text; lines 372 and 583 build the fallback identifier from the widened pair with `%d`, printing 0.
- The same file, lines 401-403, stores the widened span, and `Bio/EnsEMBL/VEP/AnnotationType/RegFeat.pm`, line 125, overlaps variants against it.

**What correcting it would change.** The attribute would be ignored or refused and the promoter would keep its core start, as vep-rs's rows have it.

**The class it produces.** `promoter_window_non_numeric_bound`: under `--extended_promoters`, a `RegulatoryFeature` row with `BIOTYPE=promoter` the reference alone writes, at a position before the promoter's core start, for a regulatory GFF record whose `extended_start` is not a whole number. The class is documented only: the outputs do not carry the attribute, no manifest names the class and no comparator applies it. No row on the committed corpora or their comparison cases: the 984 `extended_start` and `extended_end` values of the GRCh38 regulatory GFF corpus are all integers, so the corpus's `--extended_promoters` comparison case has 0 rows in each of its four formats.

**How often the shape occurs.** Tier 3: the shape is a promoter whose `extended_start` is not an unsigned integer. None of the 27,530 records of the committed regulatory and promoter GFF files has one, and 1 flag set runs `--extended_promoters`.

### 27.1 The promoter ENSR21_B34QG with its `extended_start` written as `x`

**In short.** In the GRCh38 regulatory GFF corpus the promoter `ENSR21_B34QG` spans `21:5553169-5553746` with `extended_start=5553154`, and under `--extended_promoters` both engines write its row for the variants at 5553154, 5553168 and 5553169 and not for the one at 5553153. With the attribute changed to `x` in a copy of the file, Ensembl VEP's promoter starts at 0, so all four variants gain the row, with every other record of the buffer up to 5553746. Its warnings file carries `Argument "x" isn't numeric in numeric lt (<)`. vep-rs keeps the core start and writes the row for 5553169 alone among the four, since `x` names no base, warning `Ignoring extended_start 'x' of promoter ENSR21_B34QG`.
