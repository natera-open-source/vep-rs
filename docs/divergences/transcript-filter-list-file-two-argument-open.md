# An `in` list file whose name carries a mode character

**The defect in one sentence.** Ensembl VEP 116.3 opens the file an `in` list names with Perl's two-argument `open`, which reads a mode out of the name. A leading `<` is stripped and the rest of the name is opened; a leading `>` or `>>` opens the rest of the name for writing, truncating or extending that file. A trailing `|` runs the rest of the name as a command and reads its output as the list; the bare name `-` reads the standard input. vep-rs opens the file of the written name and reads its lines, and nothing else is read, written or run.

**Why this is a defect.** The `in` operator takes a comma-separated list or a file with one value per line, and the run first tests that a file of the written name exists. Perl's two-argument `open` then reads the same string again, as a mode followed by a file name or a command. So a file named `<ids.txt` passes the existence test and the list is read from `ids.txt`, and a file named `true|` passes the test and the list is the empty output of the command `true`. A list name can therefore truncate a file or run a command, and the name the user wrote is tested and then not opened.

**What vep-rs does instead.** vep-rs opens the file of the written name and reads its lines, whatever characters the name holds. The mode reading is never replicated: a name beginning with `<` or `>`, or ending with `|`, is a file name and nothing else, so no other file is read, truncated or written and no program runs. A name that is no path is refused as `Could not find/parse list <name>`, and a file that cannot be read as `Could not read from file <name>`, the two messages Ensembl VEP prints for those cases. The same rule governs `vep filter`, whose `filter_vep` counterpart shares the grammar.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/FilterSet.pm`):

- Lines 385-393, `finish_filter_node`: an `in` value with text on both sides of a comma is a comma-separated list.
- Line 396: any other value is tested as a path with `-e`, which reads the name as written.
- Line 397: `open IN, $list or throw("ERROR: Could not read from file $list\n")`, the two-argument form, in which Perl reads the mode from the string's ends: a leading `<`, `>`, `>>`, `+<`, `+>` or `|`, or a trailing `|`, with the blanks around them dropped. The bare name `-` is the standard input.
- Lines 398-410: each line the handle yields, its `\r` and vertical whitespace removed, becomes a key of the list.
- Line 415: a value that is neither a comma-separated list nor a path is refused as `Could not find/parse list`.

**What correcting it would change.** With the three-argument form, `open IN, '<', $list`, the file tested would be the file read, and a name carrying such a character would read as every other name. A list file whose name begins with a letter, a digit, `.`, `/` or `_`, does not end in `|` and is not `-` reads the same on both engines. The mode is read at the ends of the string alone, so `data/<ids.txt` names the file of that name in both.

**The class it produces.** `transcript_filter_list_file_two_argument_open`: the rows a `--transcript_filter` keeps, and the lines a `vep filter` run keeps, under an `in` list naming a file whose name begins with `<`, `>`, `>>`, `+<`, `+>` or `|`, ends with `|` or is `-`. Ensembl VEP reads another file or the standard input, writes a file or runs a program; vep-rs reads the named file. The class is documented only: no committed flag set and no committed tool run names such a file.

**How often the shape occurs.** Tier 3: the shape is a value of `--transcript_filter`, or of `vep filter`'s `--filter`, that no ordinary command line carries. Of the 511 committed flag sets, 10 pass `--transcript_filter`, one of them names a list file, `transcript_ids.txt`, and none names a file whose name carries a mode character; of the 28 committed `vep filter` runs, the one that passes an `in` list writes it inline.

## 1. A list file named `<ids.txt` beside `ids.txt`, and one named `true|`

**In short.** The working directory holds `ids.txt`, one line `ENST00001020056`, and a file named `<ids.txt`, one line `ENST00000624153`. Under `--transcript_filter 'stable_id in <ids.txt'`, vep-rs keeps the one row of `ENST00000624153` on the record `tf_7130369`; Ensembl VEP 116.3 reads `ids.txt` instead and keeps the two rows of `ENST00001020056` on the record `tf_32493707`. With a file named `true|` holding `ENST00000624153`, `stable_id in true|` keeps that row on vep-rs, where Ensembl VEP runs the command `true`, reads an empty list and drops every transcript.

### 1.1 What a reader needs to know first

`--transcript_filter` keeps the transcripts a filter passes before any consequence is computed; a record whose every transcript is dropped is written as `intergenic_variant`, once per allele. The `in` operator takes a comma-separated list or the name of a file with one value per line, relative to the working directory or absolute. A filter string is split at blanks, so the name holds none. Perl's `open` has a two-argument form, in which the mode is part of the file name string, and a three-argument form, in which the mode is its own argument and the name is taken as written.

### 1.2 The records

The input is two records of the `tests/golden/116/GRCh38-transcript-filter` corpus:

```
21	7130369	tf_7130369	A	G	.	.	.
21	32493707	tf_32493707	C	T,G	.	.	.
```

The first lies within reach of 3 transcripts, `ENST00000624153` among them; the second within reach of 48, `ENST00001020056` among them, each written twice for the record's two alleles. Both engines run over the GRCh38 release 116 cache in the working directory described above, with the filter quoted below added.

### 1.3 What each engine writes

Ensembl VEP 116.3 with `--transcript_filter 'stable_id in <ids.txt'` reads `ids.txt` and writes:

```
tf_7130369	21:7130369	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
tf_32493707	21:32493707	T	ENSG00000166979	ENST00001020056	Transcript	intron_variant,NMD_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
tf_32493707	21:32493707	G	ENSG00000166979	ENST00001020056	Transcript	intron_variant,NMD_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
```

vep-rs with the same filter reads the file named `<ids.txt` and writes:

```
tf_7130369	21:7130369	G	ENSG00000280330	ENST00000624153	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=82;STRAND=1
tf_32493707	21:32493707	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
tf_32493707	21:32493707	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

With `--transcript_filter 'stable_id in true|'`, the file `true|` holding `ENST00000624153`, Ensembl VEP writes the three intergenic rows of the two records, and vep-rs the second block above. A file of the same content whose name ends in `|` and whose command part writes a file gives Ensembl VEP the same three rows: the command ran and a new file appeared in its working directory. vep-rs writes the second block and no file.

The counts below hold over the corpus's 32 records in each of the four output formats. Under `<ids.txt`, Ensembl VEP writes 2 transcript rows and 33 intergenic, vep-rs 1 and 34. Under `true|` and the file-writing name, Ensembl VEP writes 0 and 35, vep-rs 1 and 34. Under `stable_id in ids.txt`, both write 2 and 33.

### 1.4 Why Ensembl VEP writes that

The value `<ids.txt` has no comma, so the run tests `-e '<ids.txt'`, true for the file of that name. The two-argument `open IN, '<ids.txt'` then reads the leading `<` as the mode and opens `ids.txt`, whose one line fills the list, so `ENST00001020056` is in it and `ENST00000624153` is not.

For `true|`, the existence test passes on the file of that name and `open IN, 'true|'` reads the trailing `|` as a pipe from a command: Perl runs `true`, the handle yields nothing, and the list is empty, so no transcript is in it. A command part holding a shell character runs through the shell, so a name of that shape writes or truncates a file; the handle yields nothing. A name ending in `|` whose rest is no program that can run fails the `open`, and the run dies with `Could not read from file`.

### 1.5 Why vep-rs's row is the correct one

The filter names a file, and the run has just tested that a file of exactly that name exists. Reading a different file, writing one or running a program, none of them named by the user, is not a reading of the list; the name that passed the test is the name to open. vep-rs opens that name and reads its lines, so the transcripts it keeps are the ones the named file lists. A list written as `ids.txt` reads the same on both engines.

### 1.6 Try it yourself

```
printf 'ENST00001020056\n' > ids.txt
printf 'ENST00000624153\n' > '<ids.txt'
vep -i records.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --transcript_filter 'stable_id in <ids.txt'
```

vep-rs keeps `ENST00000624153`. Ensembl VEP 116.3, run in the same directory over the GRCh38 release 116 cache, keeps the two `ENST00001020056` rows. A name ending in `|` runs the rest of the name as a command on Ensembl VEP, so it is shown here in words alone.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/transcript_filter_list_file_two_argument_open.toml` stages a file named `true|` holding `ENST00000624153` beside a program named `true` that prints nothing, names the file as the list, and holds both engines' `Feature` on `tf_7130369`: `-` on Ensembl VEP, the transcript on vep-rs. The class covers the rows a `--transcript_filter` keeps, and the lines a `vep filter` run keeps, under an `in` list whose file name carries a mode character. No committed flag set or tool run reaches the shape, so no manifest names the class and the class is documented only.
