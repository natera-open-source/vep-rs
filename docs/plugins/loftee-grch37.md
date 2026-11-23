# LoFTEE on GRCh37

What a GRCh37 run of LoFTEE passes, what vep-rs refuses there, and where the
`master` branch of Ensembl VEP's LoF plugin reads its data differently.

On GRCh37 Ensembl VEP runs the `master` branch of the LoF plugin, and vep-rs
reads the keys that branch reads; [LoFTEE](loftee.md) lists the parameters and
the assembly each is read on.

## The key a GRCh37 run passes

The `master` branch runs its splice-prediction layer unless
`run_splice_predictions` holds `0` (`LoF.pm` lines 98 and 186); the `grch38`
branch leaves the layer off by default. vep-rs does not implement that layer,
so on GRCh37 it refuses a run that leaves the layer on: the key's last value
must be `0`. A run without the key, or with an empty value, is refused at
plugin initialization under one line, wrapped here:

```text
LoFTEE: this GRCh37 run asks for the splice-prediction layer, which
vep-rs does not implement: Ensembl VEP's LoF plugin runs it on that
assembly unless run_splice_predictions is 0; pass
run_splice_predictions:0
```

A GRCh37 run therefore names the key. `$LOFTEE` is the unpacked bundle for the
run's assembly, `/path/to/loftee/grch37`:

```bash
--plugin LoFTEE,human_ancestor_fa=$LOFTEE/human_ancestor.fa,\
gerp_tabix=$LOFTEE/GERP_scores.final.sorted.txt.gz,\
run_splice_predictions=0
```

An Ensembl VEP `LoF` argument string carrying `run_splice_predictions:0` is
read the same way. The assembly is the one `--assembly` names or, without the
flag, the cache's. On GRCh38 a command line without the key starts: the
`grch38` branch's default is the layer off.

LoFTEE reads the names `GRCh37`, `hg19` and `37`, in any letter case, as
GRCh37. Under any other name for that assembly, such as a cache that calls
itself `GRCh37.p13`, the run warns that the name is not one LoFTEE recognises,
reads the keys and defaults of the `grch38` branch, and is not refused.

## What a refused plugin does to the run

LoFTEE refuses, at plugin initialization, what it does not compute: a GRCh37
run that leaves the splice-prediction layer on, a true `run_splice_predictions`
on either assembly, and a database connection asked for under the word `mysql`.
It refuses a data file it cannot open the same way. Without `--safe` the refused plugin is
dropped with a warning, the run exits 0, and no row carries a `LoF` field.
Under `--safe` the run stops with exit status 1 and writes no output. Any
plugin that fails to start is handled the same way.

## What the `master` branch does besides

Each row is a behaviour of the `master` branch beside what vep-rs does with it.

| Subject | The `master` branch | vep-rs |
| --- | --- | --- |
| The default of `human_ancestor_fa` | `human_ancestor.fa.rz`, a razip-compressed file in the working directory, read with `samtools faidx` (`LoF.pm` lines 77 and 553). The `ANC_ALLELE` test runs wherever that file and a program that reads it are found. | The same file, read when it and its `.fai` are in the working directory; without them the test is off, as under `human_ancestor_fa=false`. vep-rs reads the razip container itself. The `ensemblorg/ensembl-vep` image reads no base from such a file: it ships no `samtools`, and `samtools` 1.21 refuses the format. |
| A `gerp_bigwig` | Never read. With no `gerp_file`, the branch reads its default, `/vep/loftee/GERP_scores.final.sorted.txt.gz` (`LoF.pm` line 132). | The stand-in rule holds on both assemblies: with no `gerp_file`, the file named under `gerp_bigwig`, `gerp_tabix`, `gerp_scores` or `gerp` is read. |
| `gerp_file:mysql` | Opens a database connection and sums GERP from its `gerp_bases` table (`LoF.pm` lines 134-136; `gerp_dist.pl` line 82). | Refused at plugin initialization: `LoFTEE: gerp_file=mysql asks for a database connection, which vep-rs does not open: pass the GERP table's path under gerp_file or gerp_tabix`. `conservation_file:mysql` (`LoF.pm` line 122) is refused under the same reason. |
| Whole-exon GERP totals | Each whole exon the GERP walk crosses takes its total from the `gerp_exons` table of the `conservation_file` database (`gerp_dist.pl` lines 59-60 and 129); an exon the table lacks counts 0. | The same, under a `conservation_file` naming the database. Without one, each such exon's total is the sum of the GERP file over the exon. |
| [A failed GERP read](#a-failed-gerp-read-counted-as-zero) | A `tabix` error line counts 0 toward `GERP_DIST`, its third field being summed with the scores (`gerp_dist.pl` lines 92-96), and the row is written: a defect, the message being no score. | The row fails: it keeps no `LoF` field, under the plugin's `went wrong` warning. |

A GRCh37 run therefore passes its GERP file under `gerp_tabix` or `gerp_file`.
Under a `conservation_file` naming the database, the whole-exon totals and the
PhyloCSF items come from it. With no `conservation_file`, or
`conservation_file:false`, the whole-exon totals come from the GERP file.

The sixth is the row the `master` branch leaves unwritten. Without a
`conservation_file` it writes no `LoF` field on a stop-gained or frameshift row
whose GERP walk crosses a whole exon, and vep-rs writes the row: the class
`loftee_gerp_walk_without_conservation_database`, whose
[page](../divergences/loftee-gerp-walk-without-the-conservation-database.md)
gives the record and both engines' rows.

The table's last row is a class as well,
`loftee_gerp_read_error_counted_as_zero`. The section
[A failed GERP read counted as zero](#a-failed-gerp-read-counted-as-zero)
gives its record and both engines' rows.

## How GRCh37 is compared with Ensembl VEP

The comparison runs on three datasets: genome-wide ClinVar, gnomAD v2.1.1 and
1000 Genomes Phase 3. The variants are the real records of those datasets. The
conservation data are synthetic: the bundle's GERP table is not redistributed
here, so both engines read a made-up GERP table and a made-up conservation
database. A release's plugin record and figure name the row
`LoFTEE (synthetic conservation data)`.

`tools/vep-loftee-synthetic-conservation` writes the two files from the GRCh37
cache. Every value is a function of a seed and the row's key, and the
comparison reads the files of seed 1. `<json_cache>` is the converted GRCh37
cache and `<dir>` the directory the files are written to:

```bash
vep-loftee-synthetic-conservation --cache <json_cache> --seed 1 \
  --out-dir <dir>
```

| File | What it holds | How the `master` branch reads it |
| --- | --- | --- |
| `GERP_scores.synthetic.txt.gz` with its `.tbi` | A made-up score for every base that lies in an exon of the cache's transcripts, each base once. | With `tabix`, the third tab-separated field of a row as one base's score (`LoF.pm` lines 132 and 138-141; `gerp_dist.pl` lines 92-96 and 119-121). |
| `phylocsf_gerp.synthetic.sql`, table `gerp_exons` | A made-up whole-exon total for every exon of every transcript. | By transcript and exon number, for each whole exon the GERP walk crosses (`gerp_dist.pl` lines 58-62 and 129-133). |
| The same file, table `phylocsf_data` | A made-up PhyloCSF row for every exon of every transcript. | By transcript and exon number, for the PhyloCSF items (`LoF.pm` lines 290-297 and 567-570). |

The run's ancestral FASTA is Ensembl's GRCh37 ancestral FASTA, the file
AncestralAllele reads, written under bare sequence names. The plugin asks
`samtools` for a region under the name of the variant's sequence (`LoF.pm`
lines 552-553). The file is block-gzipped and indexed with `samtools faidx`.

Both engines take one plugin string. `$RUN` is the directory
`grch37/loftee_synthetic_conservation` of the plugin data, and `<loftee>` is a
checkout of the plugin's `master` branch at commit
`c7a87dffe3b96a729d299ca3ccb1c25db7b44e30`:

```bash
--plugin LoF,loftee_path:<loftee>,\
human_ancestor_fa:$RUN/human_ancestor.fa.gz,\
gerp_file:$RUN/GERP_scores.synthetic.txt.gz,\
conservation_file:$RUN/phylocsf_gerp.synthetic.sql,\
run_splice_predictions:0
```

[Reproducing the plugin concordance](reproducing-concordance.md) gives the
commands that write the files and run both engines.

## A failed GERP read counted as zero

**The defect in one sentence.** On GRCh37 the LoF plugin that Ensembl VEP runs
for LoFTEE reads its GERP file through a `tabix` command and sums the third
field of every line the command prints, on standard output or standard error.
When a block of the file cannot be read, the error line has no third field and
counts 0, and the row is written from the total that results. vep-rs fails the
row, which carries none of the plugin's fields, under a `went wrong` warning.

**Why this is a defect.** An error message is not a score. The plugin writes
`GERP_DIST:0`, the text it writes for bases whose scores sum to 0, and neither
its output nor a warning says that a read failed. The `END_TRUNC` filter
applies at a total of 180 or less, so the bases the plugin could not read can
turn an `HC` call into `LC`. The module treats an unreadable file as fatal
elsewhere: at start-up it stops when `tabix` reports a failure, and its
`grch38` branch dies on a GERP file it cannot open ([what a failed read
does](loftee.md#which-gerp-file-is-read-and-what-a-failed-read-does)).

**What vep-rs does instead.** vep-rs reads the file itself. A read that fails
ends the plugin's work on the row: the row is written without `LoF`,
`LoF_filter`, `LoF_flags` and `LoF_info`, and the warnings file names the file
and the failure. The run goes on, every other row is annotated, and the exit
status is 0. A file whose first block cannot be read is refused when the plugin
loads ([What a refused plugin does to the
run](#what-a-refused-plugin-does-to-the-run)).

**Where it is in Ensembl VEP's source** (`LoF.pm` and `gerp_dist.pl` of the
LOFTEE repository's `master` branch, the plugin for GRCh37):

- `LoF.pm`, line 130: the program is `tabix` unless `tabix_path` names another.
  Lines 138-139: at start-up the plugin runs `tabix -l` on the file and stops
  only when the answer holds the word `fail`. Line 141: the GERP source is the
  string `tabix <file>`.
- `LoF.pm`, line 258: a stop-gained or frameshift row with a coding position
  calls `get_gerp_weighted_dist` with that string.
- `gerp_dist.pl`, lines 46-57: the stop codon's exon and the variant's exon go
  to `get_interval_gerp`, which hands a `tabix` source to
  `get_interval_gerp_tabix` (lines 70-77).
- `gerp_dist.pl`, line 92: the command runs with `2>&1`, standard error joined
  to standard output, and its answer is split into lines. Lines 93-97 split
  each line on tabs and add its third field to the total, and line 98 returns
  the total. Nothing reads the command's exit status.
- `LoF.pm`, lines 259 and 268: the total is printed as `GERP_DIST`, and
  `END_TRUNC` applies when the distance from the last exon is at most 50 and
  the total at most 180.

**What correcting it would change.** Reading the command's standard error apart
from its rows, and ending the plugin's work on the row when it holds a line,
would leave the row without the plugin's fields under Ensembl VEP's
`Plugin 'LoF' went wrong` warning. That is the row vep-rs writes.

**The class it produces.** `loftee_gerp_read_error_counted_as_zero`: on a
GRCh37 run of the plugin, the stop-gained and frameshift rows whose walk reads
a part of the GERP file that cannot be read. The reference row carries `LoF`
and a `LoF_info` whose `GERP_DIST` counts 0 for that part, with
`LoF_filter=END_TRUNC` where the total is at most 180 and the 50-base rule
fails. The candidate row carries none of the four fields. The class is
documented only: no committed run reads a damaged GERP file, so no comparator
codes it.

**How often the shape occurs.** Tier 3 of [Intended divergences from Ensembl
VEP](../intended-divergences.md): the shape is a GERP file holding a compressed
block that cannot be inflated, which breaks the file's format. No committed
flag set passes such a file.

### 1. A stop gain whose GERP block cannot be inflated

**In short.** `21 37443370 A>T` on GRCh37 is a stop gain in the last exon of
CBR1's two-exon transcript ENST00000439427. A made-up GERP table holds the
record's base in a compressed block whose first data byte is overwritten.
Ensembl VEP 116.3, running the plugin's `master` branch, writes `LoF=LC` with
`LoF_filter=END_TRUNC` and `GERP_DIST:0`, and vep-rs writes the row without a
LoF field and warns that the read failed. On the table before the overwrite
both engines write `LoF=HC` with `GERP_DIST:181`.

#### 1.1 What a reader needs to know first

For a stop-gained or frameshift variant the plugin walks the transcript's exons
from the variant to the stop codon and totals the per-base GERP scores of each
stretch, as [the walk's own
page](../divergences/loftee-gerp-walk-without-the-conservation-database.md#11-what-a-reader-needs-to-know-first)
describes. `LoF_info` prints the total as `GERP_DIST`. The `END_TRUNC` filter
marks a variant near the end of the coding sequence. On GRCh37 it applies when
`DIST_FROM_LAST_EXON` is at most 50 and `GERP_DIST` at most 180, and a filtered
call is `LC`.

The plugin does not open the GERP file. It runs `tabix` on it once for each
stretch and reads what the program prints. `tabix` reads a text table
compressed in blocks, the BGZF format, through an index that maps a position to
the block holding its row. Each block is inflated on its own. The plugin takes
the third tab-separated field of a row as the score (`gerp_dist.pl` line 96).

#### 1.2 The record and the damaged table

The record is hand-written and read against the cache of the `GRCh37` golden
corpus (`tests/golden/116/GRCh37/json_cache`):

```
21	37443370	last_exon_stop	A	T	.	.	.
```

It turns codon 138 of CBR1 into a stop on ENST00000439427, whose exons are
`21:37442313-37442702` and `21:37443248-37443881` and whose stop codon ends at
`37443495`. The gene's three-exon coding transcripts end their second exon at
`37443355`, so on them the record is an intron variant. The walk on
ENST00000439427 is one stretch, `21:37443370-37443495`, inside the exon that
holds the stop codon. It crosses no whole exon, so [the difference of the
conservation
database](../divergences/loftee-gerp-walk-without-the-conservation-database.md)
does not arise.

The GERP table is made up. It holds two rows of three tab-separated columns, a
chromosome, a position and a score:

```
21	37443369	1
21	37443370	181
```

Each row is compressed as a BGZF block of its own, and the index is built on
that table. The second row scores 181, one above the 180 at which `END_TRUNC`
stops applying, so the table read whole leaves the call `HC`.

One byte is then overwritten: the first byte of the second block's compressed
data becomes `0xff`. Its low bits open a deflate block of type 3, which the
deflate format reserves as an error (RFC 1951, section 3.2.3), so no inflater
reads the block. The block's header, checksum and length field, the first
block and the index are untouched. The first block stays whole because vep-rs
reads a table's first lines when the plugin loads.

Both engines take the GRCh37 assembly and one plugin flag.
`human_ancestor_fa:false` switches the ancestral-allele filter off, and
`run_splice_predictions:0` the splice-prediction layer, which vep-rs does not
run:

```
--plugin LoF,loftee_path:<loftee>,gerp_file:GERP_scores.synthetic.txt.gz,human_ancestor_fa:false,run_splice_predictions:0
```

#### 1.3 What each engine writes

The record has nine rows, one of them `stop_gained`. That row in the default
format, from Ensembl VEP 116.3:

```
last_exon_stop	21:37443370	T	ENSG00000159228	ENST00000439427	Transcript	stop_gained	513	412	138	K/*	Aaa/Taa	-	IMPACT=HIGH;STRAND=1;LoF=LC;LoF_filter=END_TRUNC;LoF_info=GERP_DIST:0,BP_DIST:125,PERCENTILE:0.767225325884544,DIST_FROM_LAST_EXON:-122,50_BP_RULE:FAIL
```

and from vep-rs:

```
last_exon_stop	21:37443370	T	ENSG00000159228	ENST00000439427	Transcript	stop_gained	513	412	138	K/*	Aaa/Taa	-	IMPACT=HIGH;STRAND=1
```

The record's eight other rows are the same from both engines. Ensembl VEP exits
with status 0 and writes no warnings file. vep-rs exits with status 0 and
writes one line to its warnings file:

```
WARNING: Plugin 'LoFTEE' went wrong: LoFTEE: failed to read GERP_scores.synthetic.txt.gz: the data provided to a libdeflater *_decompress function call was invalid in some way (e.g. bad magic numbers, bad checksum); the plugin annotated every row but the ones it failed on, the first of them on the variant at 21:37443370; 1 variant affected
```

On the table before the overwrite both engines write the same nine rows, the
stop-gained one as:

```
last_exon_stop	21:37443370	T	ENSG00000159228	ENST00000439427	Transcript	stop_gained	513	412	138	K/*	Aaa/Taa	-	IMPACT=HIGH;STRAND=1;LoF=HC;LoF_info=GERP_DIST:181,BP_DIST:125,PERCENTILE:0.767225325884544,DIST_FROM_LAST_EXON:-122,50_BP_RULE:FAIL
```

Asked for `21:37443370-37443495` on the damaged table, the `tabix` that the
image ships, htslib 1.9, prints no row. It prints one line on standard error
and exits with status 0:

```
[E::bgzf_uncompress] Inflate operation failed: invalid block type
```

Other damage leaves that `tabix` silent. A second block cut short at the end of
the file, or one whose checksum does not match its data, gives no row, no
message and exit status 0. Ensembl VEP writes the same `GERP_DIST:0` row, and
vep-rs fails the row under a warning that names the failure.

#### 1.4 Why Ensembl VEP writes that

At start-up the plugin asks `tabix` for the file's chromosome names and stops
only if the answer holds the word `fail`. The index is whole, so the answer is
`21` and the plugin starts.

For the stop-gained row the walk has one stretch. The plugin runs `tabix` on
`21:37443370-37443495` with standard error joined to standard output. `tabix`
finds the block through the index, fails to inflate it, prints its error line
and exits. The plugin splits the answer into lines, splits each line on tabs
and adds the third field to the total. The error line holds no tab, so it has
no third field, and Perl adds an undefined value as 0.

The total is printed as `GERP_DIST:0`. The record lies in the last exon, so
`DIST_FROM_LAST_EXON` is negative and the 50-base rule fails. With a total of
at most 180, `END_TRUNC` applies and the call becomes `LC`. The plugin never
reads the command's exit status, which is 0 here in any case.

Source (LOFTEE `master` branch):

- `LoF.pm`, lines 138-141: the start-up test and the command string; line 258:
  the walk's call; lines 259 and 268: `GERP_DIST` and the `END_TRUNC` test.
- `gerp_dist.pl`, line 92: the command and `2>&1`; lines 93-97: the sum of
  every line's third field.

#### 1.5 Why vep-rs's row is the correct one

The table gives the walk's stretch one score, 181 at the record's own base.
Read whole, the total is 181. That is above 180, so `END_TRUNC` does not apply
and the call is `HC`, the row both engines write on the table before the
overwrite (1.3).

With the block unreadable the file gives the stretch no total, neither 181 nor
0. Ensembl VEP's row states a total of 0 and filters the variant on it, a
verdict the table's own score contradicts. Nothing in its output separates the
row from one whose bases do sum to 0.

vep-rs's row states nothing the file did not give, and its warning names the
file and the failure. The other items need no GERP score: `BP_DIST` is the
stretch's end minus its start, 37443495 - 37443370 = 125, on either table. The
GERP total alone is what the failed read leaves unknown, and the `END_TRUNC`
verdict rests on it.

#### 1.6 Try it yourself

```
printf '21\t37443369\t1\n' | bgzip > first.gz
printf '21\t37443370\t181\n' | bgzip > second.gz
head -c "$(( $(wc -c < first.gz) - 28 ))" first.gz > first_block
{ cat first_block; head -c 18 second.gz; printf '\xff'; tail -c +20 second.gz; } \
  > GERP_scores.synthetic.txt.gz
cat first_block second.gz > whole.synthetic.txt.gz
tabix -s 1 -b 2 -e 2 whole.synthetic.txt.gz
mv whole.synthetic.txt.gz.tbi GERP_scores.synthetic.txt.gz.tbi

vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh37 \
  --plugin LoF,gerp_file:GERP_scores.synthetic.txt.gz,human_ancestor_fa:false,run_splice_predictions:0
```

`<corpus>` is `tests/golden/116/GRCh37`, and `record.vcf` holds the record of
1.2. The last 28 bytes of a `bgzip` stream are its end marker, which the first
block of a longer file does not carry, and a block's compressed data starts at
its byte 19. vep-rs writes the row and the warning of 1.3.

Ensembl VEP 116.3 takes the same `--plugin` flag with `loftee_path:` added,
`--dir_plugins` naming a checkout of the LOFTEE `master` branch and the cache
options of [Reproducing these examples](../divergences/reproducing.md). It
writes the `LoF=LC` row.

#### 1.7 How the comparator treats this class

The reproduction
`scripts/divergences/reproductions/loftee_gerp_read_error_counted_as_zero.toml`
builds the table of 1.2, runs the record on both engines and records the
stop-gained row of each. The class covers the stop-gained and frameshift rows
of a GRCh37 LoFTEE run whose walk reads a part of the GERP file that cannot be
read. No committed run reads a damaged GERP file, so no comparator codes the
class.
