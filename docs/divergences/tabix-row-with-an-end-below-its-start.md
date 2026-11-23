# A VCF source row whose `END` lies below its position

**The defect in one sentence.** Ensembl VEP 116.3 reads a `--custom` VCF source, and a plugin's VCF table, through the tabix reader of htslib 1.9, which spans a row from its position to the value of its INFO `END` key. A row whose `END` lies below its position, or whose `END` is not a number (`END=x`, `END=.`, `END=`), gets a span that ends before it starts, and no query over the row's own bases returns it. A variant inside the row's reference allele carries none of the source's fields. vep-rs matches such a variant against the row as against any other row and fills the source's field.

**Why this is a defect.** The VCF specification defines `END` as the end of the record's reference allele, so an `END` below the position or one that is no number is a malformed tag, a dot the missing value, and the row still spans the allele in its `REF` column. The reader keeps such a value in place of the allele's end and hands every query an interval that ends before it starts. htslib 1.11 and later ignore such a tag with a warning and span the row by its reference allele. Ensembl VEP's own overlap test spans a row by its reference allele too, so the reader drops rows the engine would have matched.

**What vep-rs does instead.** vep-rs's `--custom` reader and its gnomADMt reader read an INFO `END` to its leading digits and keep it when it lies at or past the row's position. Otherwise the row spans its reference allele, as htslib 1.11 and later span it. Its SpliceAI reader reads no `END` and places every row at its position. Either way a variant inside the row is matched against it as against any other.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-io` at release/116, `ensembl-variation` at release/116, Bio::DB::HTS 2.11, htslib 1.9):

- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, line 321: a `--custom` source is asked for each variant's window, `seek(chr, start - 1, end + 1)`. `Bio/EnsEMBL/IO/TabixParser.pm`, line 89: the seek is `Bio::DB::HTS::Tabix->query("chr:start-end")`. A plugin's table takes the same route through `Bio/EnsEMBL/Variation/Utils/BaseVepTabixPlugin.pm`, lines 596-602, `_get_data_hts`.
- `Bio/DB/HTS/Tabix/Iterator.pm`, line 57: each row of the answer is read with htslib's `tbx_iter_next`.
- htslib `tbx.c`, lines 121-135, `tbx_parse1`: under the VCF preset a row's end is the last base of its `REF` allele, the zero-based start plus the column's length (line 123). An `END=` in the INFO column replaces it with `strtol` of the text after the `=` (line 133): the value as written for `END=5029000`, `0` for `x`, `.` or an empty value. Line 142 refuses only an end below `0`, so the row is indexed with its inverted span.
- htslib `hts.c`, lines 2641-2650, `hts_itr_next`: in htslib's zero-based coordinates a row is returned when its end lies past the query's start and its start before the query's end (line 2645), that is, when the row's last base reaches the query's first. An end below the start passes the first test only for a query that starts before that end, so a query over the row's own bases never returns the row, and a row whose end is `0` is returned by no query at all.
- htslib `tbx.c`, lines 167-178, `tbx_readrec`: the reader parses each row again as it answers a query, so the index decides nothing. An index built by htslib 1.11 or later, which spans the row by its reference allele, changes nothing in the reader's answer.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/VCF.pm`, lines 389-421, `_record_overlaps_VF`: the overlap test spans a returned row from its position to the position plus the length of its reference allele less one (lines 389-404), and the `exact` test compares the alleles (lines 405-421). Neither reads `END`, so a returned row would have matched.

**What correcting it would change.** With the row spanned by its reference allele, the reader would return it for a variant inside that allele, Ensembl VEP's overlap and exact tests would match the variants vep-rs matches, and the rows would carry the source's fields.

**The class it produces.** `tabix_row_with_an_end_below_its_start`: a `--custom` VCF source's field, or the field of a plugin reading a VCF table (gnomADMt, SpliceAI), that vep-rs alone fills on a variant inside a row whose `END` lies below its position or is not a number. The class is documented only: no committed VCF source or VCF plugin table carries an `END` key, so no manifest names the class and no comparator applies it.

## 1. Six deletion rows of a `--custom` VCF source, four with an `END` that is no end

**In short.** A VCF source holds six deletion rows, each with a four-base reference allele, between 21:5029100 and 21:5029223. Two are well formed, one without `END` and one with `END=5029123`, the end of its reference allele; the other four carry `END=5029000`, below their positions, `END=x`, `END=.`, the missing value, and `END=`. One substitution inside each row's reference allele is run under `--custom` with the source as an `overlap` source named `SPAN`. Ensembl VEP writes `SPAN=ctrl_no_end` and `SPAN=ctrl_end` on the two well-formed rows' variants and no `SPAN` on the other four. vep-rs writes each variant's row: `SPAN=end_below`, `SPAN=end_x`, `SPAN=end_dot` and `SPAN=end_empty` beside the two.

### 1.1 What a reader needs to know first

A `--custom` VCF source annotates a variant with the rows that overlap it, the row's `ID` column being the value written under the source's short name. Each engine asks the source's tabix index for the rows over the variant's window and tests each returned row against the variant: under `overlap`, the row's reference allele must share a base with the variant; under `exact`, the alleles must match. The six variants here are intergenic substitutions at the third base of each row's reference allele, so each writes one row and the only cell that can differ is the `SPAN` item.

### 1.2 The record

The source, six deletion rows against the GRCh38 reference, compressed with bgzip and indexed with `tabix -p vcf`:

```
21	5029100	ctrl_no_end	TCTG	T	.	.	DP=3
21	5029120	ctrl_end	ATGA	A	.	.	END=5029123
21	5029140	end_below	CAGC	C	.	.	END=5029000
21	5029160	end_x	CCAG	C	.	.	END=x
21	5029200	end_dot	AGAA	A	.	.	END=.
21	5029220	end_empty	TCGT	T	.	.	END=
```

The input, one substitution at the third base of each row, against the Ensembl 116 GRCh38 cache:

```
21	5029102	in_ctrl_no_end	T	G	.	.	.
21	5029122	in_ctrl_end	G	T	.	.	.
21	5029142	in_end_below	G	T	.	.	.
21	5029162	in_end_x	A	C	.	.	.
21	5029202	in_end_dot	A	C	.	.	.
21	5029222	in_end_empty	G	T	.	.	.
```

Both engines run `--custom file=span.vcf.gz,short_name=SPAN,format=vcf,type=overlap,coords=0`.

### 1.3 What each engine writes

The rows in the default format. Ensembl VEP 116.3:

```
in_ctrl_no_end	21:5029102	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;SPAN=ctrl_no_end
in_ctrl_end	21:5029122	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;SPAN=ctrl_end
in_end_below	21:5029142	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
in_end_x	21:5029162	C	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
in_end_dot	21:5029202	C	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
in_end_empty	21:5029222	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

vep-rs:

```
in_ctrl_no_end	21:5029102	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;SPAN=ctrl_no_end
in_ctrl_end	21:5029122	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;SPAN=ctrl_end
in_end_below	21:5029142	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;SPAN=end_below
in_end_x	21:5029162	C	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;SPAN=end_x
in_end_dot	21:5029202	C	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;SPAN=end_dot
in_end_empty	21:5029222	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;SPAN=end_empty
```

The tab output carries the same values in its `SPAN` column, the VCF output in the `SPAN` field of `CSQ`, and the JSON output under `custom_annotations`, where Ensembl VEP has no `SPAN` list on the four variants and vep-rs has one with the row's name. Neither engine writes a warning, and both exit 0. Running the six deletion rows themselves as the input under `type=exact` gives the same pattern: Ensembl VEP writes `SPAN` on the two well-formed rows alone, vep-rs on all six.

### 1.4 Why Ensembl VEP writes that

For the variant at 5029142 the source is asked for the window 5029141 to 5029143. The reader parses the row at 5029140 with its end at 5029000, the `END` value, where the reference allele ends at 5029143. Its test for a returned row asks whether the row's end reaches the window's start: 5029000 does not reach 5029141, so the row is skipped and the variant carries no `SPAN`. For the rows with `END=x`, `END=.` and `END=` the parsed end is `0`, which reaches no window's start, so no query returns them.

The two well-formed rows end at 5029103 and 5029123, past their variants' window starts, and are returned and matched. The row at 5029140 is returned only by a query that starts at or before 5029000 and reaches 5029140, a window no variant inside the row asks for.

### 1.5 Why vep-rs's row is the correct one

A row's reference allele is the bases it describes, and the row at 5029140 describes `CAGC` at 5029140 to 5029143 whatever its `END` says. An `END` below the position cannot be the end of that allele, and `x`, `.` and an empty value are not positions. Such a tag is dropped, as htslib 1.11 and later drop it, and the row spans its allele. The variant at 5029142 lies inside it, and Ensembl VEP's own overlap test, which spans the row by its reference allele, would have matched it had the reader returned the row, as the two well-formed rows show.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --custom file=span.vcf.gz,short_name=SPAN,format=vcf,type=overlap,coords=0
```

vep-rs writes a `SPAN` item on all six rows. Ensembl VEP 116.3 writes one on the rows for `in_ctrl_no_end` and `in_ctrl_end` alone. Setting those four rows' `END` to the end of their reference allele, or removing it, makes the engines agree.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/tabix_row_with_an_end_below_its_start.toml` holds the six rows, the six variants and both engines' rows. The class covers a `--custom` field, or the field of a plugin reading a VCF table, that vep-rs alone fills on a variant inside a row whose `END` lies below its position or is not a number. No committed dataset reaches the shape, so no manifest names the class and the class is documented only.
