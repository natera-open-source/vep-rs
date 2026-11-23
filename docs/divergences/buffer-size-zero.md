# A `--buffer_size` of zero

**The defect in one sentence.** Ensembl VEP 116.3 admits a record to its buffer only while the buffer holds fewer records than `--buffer_size`. At `--buffer_size 0` the first fill admits none, the run takes the empty buffer for the end of the input, and it exits 0 with its header written and no row. vep-rs annotates every record at zero, and its output is the one any other value gives.

**Why this is a defect.** The buffer size is a memory setting. Ensembl VEP's documentation describes it as the number of variants read into memory at once, lower for less memory and higher for speed, and excludes no value. At a size of 1 and at the default of 5000 the run of the example below writes the same rows. At zero it writes none, prints no error and no warning, and returns the exit status of a finished run. A caller that computes the value and arrives at zero loses every record of its input without a sign of it.

**What vep-rs does instead.** vep-rs puts every record it reads into a batch, and closes a batch once it holds at least `--buffer_size` records. At zero each record is a batch of its own, as at `--buffer_size 1`, and every record is annotated. The output is the one any other value gives. The command-line reference states the difference under [Compatibility with Ensembl VEP](../cli/compatibility.md).

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `Config.pm`, line 118: `buffer_size=i`, any integer. Line 300 sets the default, 5000. No other line of the file names the option, so no bound is tested.
- `InputBuffer.pm`, lines 191 and 212, `next`: both loops that fill the buffer run while `@$buffer < $buffer_size`. With a size of 0 neither admits a record, and the second never asks the parser for one.
- `Runner.pm`, lines 423-424 and 476-477: an empty buffer gives an empty list of output lines, without and with `--fork`.
- `Runner.pm`, lines 375-397, `next_output_line`: an empty list gives no line, which ends the loop of line 208 that writes the output. Lines 204-206 wrote the header before it.

**What correcting it would change.** A buffer that admits at least one record per fill, or a refusal of a size below 1, would end the silent loss. With the first, Ensembl VEP would write at zero the rows it writes at any other size, which are the rows vep-rs writes.

**The class it produces.** `buffer_size_zero_writes_the_header_alone`: every row of a run given `--buffer_size 0`. The reference, Ensembl VEP's output, holds its header and no row where the candidate, vep-rs's output, holds every record's rows. The class is documented only: no committed flag set passes the value, and no comparator codes it.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a `--buffer_size` of zero. Every one of the 511 committed flag sets is run with `--buffer_size 5000`, and none with zero.

## 1. One record read with a buffer of zero

**In short.** The record `21 45573402 C T` overlaps five transcripts. Under `--buffer_size 0`, Ensembl VEP 116.3 writes the header of its output and no row, and exits 0. vep-rs writes the record's five rows. They are the five rows both engines write at `--buffer_size 1` and at the default of 5000.

### 1.1 What a reader needs to know first

Both engines read the input in portions and annotate one portion at a time. `--buffer_size` sets how many records a portion holds, 5000 unless the run says otherwise. A larger portion uses more memory and fewer reads of the cache. The setting is meant to change how a run uses its machine and never what the run writes.

### 1.2 The record

The record has the position and alleles of the record `ps_45573402` of the `GRCh38-pick-summary` golden corpus, under an identifier of its own, and is read against that corpus's cache (`tests/golden/116/GRCh38-pick-summary/json_cache`):

```
21	45573402	utr_and_intron	C	T	.	.	.
```

Both engines take the default output format and one flag:

```
--buffer_size 0
```

### 1.3 What each engine writes

Ensembl VEP 116.3 writes the 43 header lines of its default output, the last of them the column names, and no row below them. Its exit status is 0. It writes nothing to standard error and no warnings file.

vep-rs writes the record's five rows below its header:

```
utr_and_intron	21:45573402	T	ENSG00000173638	ENST00000651099	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=37;STRAND=-1
utr_and_intron	21:45573402	T	ENSG00000299414	ENST00000763267	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=350;STRAND=1
utr_and_intron	21:45573402	T	ENSG00000299414	ENST00000763268	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
utr_and_intron	21:45573402	T	ENSG00000299414	ENST00000763269	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=312;STRAND=1
utr_and_intron	21:45573402	T	ENSG00000173638	ENST00000859469	Transcript	5_prime_UTR_variant	163	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
```

At `--buffer_size 1` and at `--buffer_size 5000` both engines write those five rows. The other formats follow the default one. Under `--vcf` Ensembl VEP writes its header lines and no record line, and under `--json` an empty file. Without `--fork` it writes the same as with it.

The whole input is lost, whatever its length. On the 37 records of the corpus's own `variants.vcf`, Ensembl VEP writes no row at `--buffer_size 0`. vep-rs writes the same rows at 0, at 1 and at 5000, and its three outputs differ only in the header lines that record the run's time and its command.

A size below zero does the same in Ensembl VEP: `--buffer_size -1` writes the header alone and exits 0. vep-rs refuses a value with a minus sign as an argument syntax error, exit 2, as [Compatibility with Ensembl VEP](../cli/compatibility.md) states, so the two engines differ in rows at the value zero alone.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP fills its buffer in a loop that goes on while the buffer holds fewer records than the size. The test comes before the read. With a size of zero the buffer, which holds no record, does not hold fewer than zero, so the loop ends before it reads a line of the input.

The run then asks whether the buffer holds anything. An empty buffer is how the end of the input shows itself, so the run stops asking for output and closes its file. The header was written before the first fill. Nothing along the way tests the size, so nothing reports it.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `InputBuffer.pm`, lines 191 and 212: the two loops, each under `@$buffer < $buffer_size`.
- `Runner.pm`, lines 423-424 and 476-477: an empty buffer returns no output line.
- `Runner.pm`, line 208: the loop that writes lines ends on the first fill.

### 1.5 Why vep-rs's rows are the correct ones

The input holds one record, and the record overlaps five transcripts of the cache. Both engines write those five rows whenever the buffer size is 1 or more, so the rows are what the record's annotation is. A buffer size states how much to hold in memory. It is no filter on the input, and no value of it is documented to drop a record. A run that exits 0 has reported success, and a successful run over one record carries that record's rows. vep-rs's output at zero is byte for byte its output at any other size, the header lines that record the run's time and its command apart.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --buffer_size 0
```

vep-rs writes the five rows of 1.3. Ensembl VEP 116.3 takes the same flag with the cache options of [Reproducing these examples](reproducing.md) and writes its header alone.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/buffer_size_zero_writes_the_header_alone.toml` runs the record of 1.2 on both engines with a buffer size of 0. It records vep-rs's five rows and, for Ensembl VEP, no row. The class covers every row of such a run. A comparison of the two outputs would count each of the candidate's rows as one the reference lacks. No committed flag set passes the value, and no comparator codes the class.
