# An `in` list that names a directory

**The defect in one sentence.** Ensembl VEP 116.2 reads a directory named as an `in` list file as an empty list: the open succeeds, every read of the handle fails, the loop ends with nothing read, and the run goes on, dropping every transcript without a message. vep-rs refuses the run with `ERROR: [filter] Could not read from file <name>`, exit status 1 and no output.

**Why this is a defect.** The `in` operator takes a file with one value per line, and a directory is not such a file. Yet `-e` is true for a directory, `open` on a directory for reading succeeds on Linux, and the first `readline` on the handle fails with "Is a directory", a failure the `while` loop reads as the end of the file. The list is then empty, no transcript is in it, and every record is written as `intergenic_variant`: the output a list of unknown identifiers would give, with no sign that nothing was read. A list path mistyped as its directory silently empties the annotation.

**What vep-rs does instead.** vep-rs refuses the run with `ERROR: [filter] Could not read from file <name>`, the message Ensembl VEP prints when its own open fails, exit status 1 and no output. The empty reading is never replicated: a directory is a list that could not be read, not a list of no transcripts. The same rule governs `vep filter`, whose `filter_vep` counterpart shares the grammar and keeps no line under such a list.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/FilterSet.pm`):

- Line 396, `finish_filter_node`: an `in` value that is no comma-separated list is tested as a path with `-e`, true for a directory.
- Line 397: `open IN, $list or throw("ERROR: Could not read from file $list\n")`; `open(2)` on a directory for reading succeeds, so the `or` branch is not taken.
- Lines 398-410: `while(<IN>)` reads the handle; `read(2)` on a directory fails with `EISDIR`, `readline` returns nothing, the body never runs and the list stays empty.
- Line 718, `filter_in`: `exists` on the empty list, false for every value, so a unit with a value fails on every transcript.

**What correcting it would change.** With `-f $list` in place of `-e`, or the handle's read error checked, a directory would be refused as a list that cannot be read, as a path that cannot be opened is; a regular file would read as before.

**The class it produces.** `transcript_filter_list_directory_read_as_empty`: the outcome of a `--transcript_filter` or `vep filter` run whose `in` list names a directory. Ensembl VEP writes every record as intergenic and `filter_vep` keeps no line; vep-rs refuses the run. The class is documented only: no committed flag set and no committed tool run names a directory as a list.

**How often the shape occurs.** Tier 3: the shape is a value of `--transcript_filter`, or of `vep filter`'s `--filter`, that no ordinary command line carries. Of the 392 committed flag sets, 9 pass `--transcript_filter`, one of them names a list file, and none names a directory; of the 28 committed `vep filter` runs, the one that passes an `in` list writes it inline.

## 1. A directory named as the list

**In short.** Under `--transcript_filter 'stable_id in <dir>'`, with `<dir>` an existing directory, Ensembl VEP 116.2 writes the record `tf_7130369`, within reach of 3 transcripts, as one `intergenic_variant` row and prints nothing on its error stream. vep-rs ends with `ERROR: [filter] Could not read from file <dir>` and writes no output. Over the corpus's 32 records, Ensembl VEP writes 35 intergenic rows and none of the 952 transcript rows.

### 1.1 What a reader needs to know first

`--transcript_filter` keeps the transcripts a filter passes before any consequence is computed; a record whose every transcript is dropped is written as `intergenic_variant`, once per allele. The `in` operator takes a comma-separated list or the name of a file with one value per line. Both engines refuse a name that is no path, `Could not find/parse list <name>`, and a path that cannot be opened for reading, such as a socket, `Could not read from file <name>`.

### 1.2 The record

The input is the `tests/golden/116/GRCh38-transcript-filter` corpus's own record:

```
21	7130369	tf_7130369	A	G	.	.	.
```

Both engines run over the GRCh38 release 116 cache, with `--transcript_filter 'stable_id in <dir>'` added, `<dir>` an empty directory.

### 1.3 What each engine writes

Ensembl VEP 116.2 writes one row and nothing on its error stream:

```
tf_7130369	21:7130369	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

vep-rs writes no output file and ends with exit status 1 and one line:

```
ERROR: [filter] Could not read from file <dir>
```

Over the corpus's 32 records in each of the four output formats, Ensembl VEP writes the 35 intergenic rows, one per allele, and vep-rs the one error line. With the record's three transcripts listed in a regular file instead, both engines write its three transcript rows.

### 1.4 Why Ensembl VEP writes that

The value has no comma, so the run tests `-e <dir>`, true. The two-argument `open IN, <dir>` succeeds, since Linux opens a directory for reading. The first `readline` fails with `EISDIR`; the `while` loop takes the undefined return as the end of the file, its body never runs, and the list is an empty hash. `filter_in` tests `exists` on that hash, false for every stable identifier, so every transcript is dropped and the record is written as intergenic.

### 1.5 Why vep-rs's row is the correct one

A directory holds no lines, so no list was read, and a filter whose list was never read has no meaning: the run cannot tell the user's list from an empty one. Ensembl VEP itself refuses the two neighbouring cases, a path that does not exist and a path that cannot be opened, with a message; a path that opens but cannot be read is the same failure one call later. vep-rs reads the file whole, meets the read error and refuses with the message Ensembl VEP prints for an open that fails, so a mistyped path ends the run with its name rather than emptying the annotation.

### 1.6 Try it yourself

```
mkdir ids
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --transcript_filter 'stable_id in ids'
```

vep-rs ends with `ERROR: [filter] Could not read from file ids`. Ensembl VEP 116.2 with the same filter over the GRCh38 release 116 cache writes the record's intergenic row.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/transcript_filter_list_directory_read_as_empty.toml` stages an empty directory, names it as the list, and holds Ensembl VEP's `Feature` on `tf_7130369`, `-`, beside vep-rs's exit status and empty output. The class covers the outcome of a `--transcript_filter` or `vep filter` run whose `in` list names a directory. No committed flag set or tool run reaches the shape, so no manifest names the class and the class is documented only.
