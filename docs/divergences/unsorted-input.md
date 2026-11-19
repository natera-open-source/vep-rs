# An input that is not sorted by position

**The defect in one sentence.** Ensembl VEP 116.2 counts, over the whole run, the records that start more than 100 bases before the record admitted before them in the same fill of its buffer, and ends the run once the count passes 100. The first record of each fill is compared with nothing, so the count follows the buffer and not the file: one input is refused at the default buffer size, written in part at 50 and written whole at 1. vep-rs annotates each record wherever it stands in the input and writes the rows in input order in every text format.

**Why this is a defect.** A record's annotation does not depend on the records around it. Ensembl VEP shows as much: under `--no_check_variants_order` it writes every record of the same file, with the rows it writes for the file sorted. The check that ends the run guards speed and memory, and its verdict is no property of the input. It compares afresh at every fill, so a smaller buffer counts fewer records out of order and a buffer of one counts none. Its message asks for an input sorted by chromosome, an order the check never reads. A run ended this way has written the rows of every earlier fill, so its output can hold part of the input.

**What vep-rs does instead.** vep-rs annotates each record against the cache wherever the record stands in the input, and writes the rows in input order in every text format. It accepts `--no_check_variants_order`, and its text outputs hold the same rows with and without it ([Input and output](../cli/input-output.md)). Its Parquet writer, an output Ensembl VEP does not have, declares a sort order and so checks the input's order itself ([Parquet format](../output/parquet.md)). There the first record out of order ends the run with exit 1 and leaves the output path as it was, unless `--no_check_variants_order` is passed.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`):

- `Constants.pm`, lines 348-349: the two thresholds, 100 records and 100 bases.
- `InputBuffer.pm`, lines 187-188, `next`: every fill opens with no previous chromosome and a previous start of 0.
- `InputBuffer.pm`, lines 240-250: a record admitted to the buffer is counted when its start lies more than 100 bases before the previous record's start, and its start becomes the previous one. The count is a field of the buffer object and lasts the run.
- `InputBuffer.pm`, lines 216-220 and 258-263: the fill dies with the message of line 189 once the count is over 100, before it admits the next record read and again when the fill is complete.
- `InputBuffer.pm`, lines 224-233: a record on another chromosome ends the fill and resets the previous start. No line tests the order of the chromosomes.
- `Config.pm`, lines 640-645 and 764: the identifier, HGVS, SPDI and region formats are exempt. The VCF and default formats are checked. Line 83 declares `no_check_variants_order`, which switches both the count and the exit off.

**What correcting it would change.** Without the exit, Ensembl VEP would write for every input what it writes under `--no_check_variants_order`, which is what vep-rs writes. A check made on the file alone would at least give one verdict per input, whatever the buffer size.

**The class it produces.** `unsorted_input_ends_the_run`: every row of a run whose VCF or default-format input holds more than 100 records the check counts, each starting over 100 bases before the record admitted before it in the same buffer fill. The reference, Ensembl VEP's output, ends with a non-zero exit status, holding the rows of the fills before the one that passed the count, where the candidate, vep-rs's output, holds every record's rows. The class is documented only: no committed input holds more than 100 such records, and no comparator codes it.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an input in which the check counts more than 100 records. The five cuts are sorted by position. Of the 392 committed flag sets, 384 read a VCF or default-format input without `--no_check_variants_order`, and the check can count at most 6 records of any of their inputs.

## 1. A file read from its last record to its first

**In short.** The 237 single-base chromosome 21 records of the `GRCh38-statistics` golden corpus are written in Ensembl VEP's default format, last record first. Ensembl VEP 116.2 ends with `The input file appears to be unsorted` and writes no row. vep-rs writes 3,136 rows. Ensembl VEP writes the same 3,136 rows for the same file under `--no_check_variants_order`, and under `--buffer_size 1` without that flag.

### 1.1 What a reader needs to know first

Ensembl VEP reads its input into a buffer of `--buffer_size` records, 5000 by default, and annotates one buffer at a time. It loads the cache for the regions a buffer covers, which is quickest when neighbouring records lie near each other. Its documentation says that an unsorted input slows a run down and needs more memory, and that `--no_check_variants_order` permits one. The order has no part in what a record's rows say.

### 1.2 The records

The input is built from the corpus's `variants.vcf` (`tests/golden/116/GRCh38-statistics/`): its chromosome 21 records with a one-base reference and a one-base alternate allele, as default-format lines, in reverse order. The file has 237 lines, from position 46253873 down to 8035105, and opens:

```
21	46253873	46253873	G/C	+
21	46219578	46219578	C/G	+
21	46136261	46136261	G/A	+
```

Of the 236 lines below the first, 108 start more than 100 bases before the line above them. Both engines read the file against the corpus's cache with the default output format. Beyond the cache options, the one flag of the run names the input's format, `--format ensembl`, as in the command of 1.6.

### 1.3 What each engine writes

Ensembl VEP 116.2 writes its 43 header lines and no row, and this line to standard error:

```
Exiting the program. The input file appears to be unsorted. Please sort by chromosome and by location and re-submit.
```

Its exit status is 25, with `--fork 1` and without it.

vep-rs writes 3,136 rows, the rows of all 237 records in the order of the file, and exits 0. The first two of the 24 rows of the file's first line:

```
21_46253873_G/C	21:46253873	C	ENSG00000160294	ENST00000291688	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
21_46253873_G/C	21:46253873	C	ENSG00000160294	ENST00000397708	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
```

Three more runs of Ensembl VEP on the same file show what the exit depends on. With `--no_check_variants_order` it writes 3,136 rows, equal to vep-rs's line for line. With `--buffer_size 1` and no other flag it writes the same 3,136 rows and exits 0. With `--buffer_size 50` it writes 2,879 rows, then the line above, and exits with a non-zero status.

The count of 100 is exact. The first 225 lines of the file hold 100 lines that start more than 100 bases before the line above, and Ensembl VEP writes their 3,064 rows. The first 226 lines hold 101, and it writes none. Read in ascending order, the 237 lines give the same 3,136 rows as the run under `--no_check_variants_order`, in another order.

### 1.4 Why Ensembl VEP writes that

While it fills a buffer, Ensembl VEP remembers the start of the record it admitted last. A record that starts more than 100 bases before that start adds one to a count kept for the whole run. Once the count is over 100, the fill ends the run.

The remembered start is set back to 0 when a fill begins. The first record of every fill is therefore compared with nothing. At the default size the 237 lines belong to one fill, the count passes 100 at the 226th line, and the run ends before its first row. At a size of 50 the count passes 100 in a later fill, after the earlier fills are written. At a size of 1 every record opens a fill, so none is counted.

Source (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`):

- `InputBuffer.pm`, lines 187-188: the previous start is 0 when a fill begins.
- `InputBuffer.pm`, lines 240-250: the comparison and the count.
- `InputBuffer.pm`, lines 216-220 and 258-263: the two exits.
- `Constants.pm`, lines 348-349: 100 records, 100 bases.

### 1.5 Why vep-rs's rows are the correct ones

Each of the 237 lines names a position and an allele on chromosome 21, and the cache holds the transcripts there. Nothing in a line refers to the line above it. Ensembl VEP's own output agrees: for this file it writes the 3,136 rows under `--no_check_variants_order`, under `--buffer_size 1`, and with the lines sorted. vep-rs writes those rows. The exit is a guard on the run's speed, and at the default buffer size it withholds an annotation that Ensembl VEP computes for the same file under either of two other settings.

### 1.6 Try it yourself

```
grep -v '^#' tests/golden/116/GRCh38-statistics/variants.vcf \
  | awk -F'\t' -v OFS='\t' '$1 == "21" && length($4) == 1 &&
      $5 ~ /^[ACGT]$/ { print $1, $2, $2, $4 "/" $5, "+" }' \
  | tac > descending.txt
vep -i descending.txt -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --format ensembl
```

vep-rs writes 3,136 rows. Ensembl VEP 116.2 takes the same file with the cache options of [Reproducing these examples](reproducing.md) and ends with the line of 1.3. Add `--no_check_variants_order` or `--buffer_size 1` to see it write the rows.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/unsorted_input_ends_the_run.toml` builds the file of 1.2 and runs it on both engines. It records the number of rows each writes, 3,136 for vep-rs and 0 for Ensembl VEP, and Ensembl VEP's exit status, 25. The class covers every row of such a run. Ensembl VEP's output ends before or inside the input, so a comparison has no whole reference to read. No committed input holds more than 100 records out of order, and no comparator codes the class.
