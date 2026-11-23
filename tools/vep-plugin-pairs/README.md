# vep-plugin-pairs

The plugin comparison of the concordance harness: one binary with two commands, each reading its
inputs on every core.

- `score` is the plugin tuple scorer. It scores each engine's output against one Ensembl VEP
  default-format output and writes, per source, engine and projection, the integer counts a raw
  and an adjusted F1 are derived from, with the documented classes of difference set aside and
  every tuple left over listed.
- `compare` is the per-field plugin comparison. It compares the fields plugins write, record by
  record, between an Ensembl VEP output and a vep-rs output in the default format or VCF, and
  passes or fails each plugin.

Both resolve field names through one registry, `data/plugin_fields.tsv`: which plugin writes a
field name.

```text
vep-plugin-pairs score --reference OUTPUT --engine ID=PATH [--engine ID=PATH ...] --report JSON \
    [--plugins A,B] [--custom-source NAME=PREFIX ...] [--run-command COMMAND] \
    [--engine-binary ID=SHA256 ...] [--engine-commit ID=COMMIT ...] \
    [--class-table TSV] [--classes A,B] [--class-population CLASS=N ...] \
    [--input FILE] [--gwas-catalog FILE] [--revel-table FILE] [--revel-no-match] \
    [--utrannotator-sweep PATH ...] [--json-cache DIR] [--spliceai-table FILE ...] \
    [--residuals DIR] [--work-dir DIR] [--pairing-threads N]

vep-plugin-pairs compare (--perl-output PATH --rust-output PATH | --perl-dir DIR --rust-dir DIR) \
    --plugins A,B [--fields A,B] [--json-output PATH] [--max-mismatches N] [--tolerance X] \
    [--min-annotated N] [--threads N] [--block-bytes N] [--shard-bytes N] [--work-dir DIR]
```

Both command lines are read as CPython's `argparse` reads one: `--name value` or `--name=value`,
a prefix that only one option has, the last value of an option given twice and every value of one
that repeats. `--help` after a command lists its options. `score` also takes `--pairing-tool
BINARY`, which must name a file, `--sort-memory SIZE` and `--sort-threads N`, and reads none of
them, so a command line that carries them runs as written.

## Contract of `score`

**Tuples and sources.** A plugin tuple is one populated plugin field of one output row:
(Location, Allele, Feature, field name, value). A field is populated when its value is neither
empty nor `-`, so a row where neither side writes the field contributes nothing. `--plugins A,B`
names the plugins the scored run loaded, by their names in the registry or by the names the
registry lists as their aliases. A field belongs to the one plugin among them that the registry
gives it, never to a plugin the run did not load, so a field name two plugins write is counted
once, for the run's own. `--run-command` hands the registry the words of the scored run's command
line, and a run is refused whose flags make the name of a loaded plugin's column a key of Ensembl
VEP's own as well, where both engines would score that key as the plugin's. `--custom-source
NAME=PREFIX` groups the fields of a `--custom` source: those starting with the prefix that no
plugin of the run claims. A run names its plugins, its custom sources or both. Every populated
Extra key that no source claims and that is not a column of Ensembl VEP's own is listed per input
under `unclaimed_keys` with the number of times the input writes it, up to 100,000 distinct keys
an input. So is a field two plugins of the run both write, where no custom source claims it.

**Inputs.** The reference is an Ensembl VEP default-format output, plain or gzip. A scored
engine's path is a default-format output or a tuple file, told apart by the column count of the
first data line (14 or 5). A tuple file is what a comparator's normaliser writes: tab-separated,
no header, five columns (Location, Allele and Feature as Ensembl VEP writes them, the Ensembl VEP
field name, the decoded value), one line per row that writes the value. Every input is opened and
its column header and first data line checked before any is scored. A gzip input is unpacked
whole, every member of it, zero padding after a member skipped, and digested as stored.

**Extraction.** A line ends at its line feed and loses every carriage return before it; an empty
line is skipped. In a default-format output a line that opens with `#` is no row; one that opens
with a single `#` names 14 columns, and the last such line is the input's column header. A row
has 14 columns and a tuple line five, or the run is refused with the line's number; a tuple file
has no header, so a line of it that opens with `#` is a tuple. The Location is read with each
coordinate written as its number; the Allele and the Feature lose their surrounding ASCII
whitespace, the vertical tab included. Of the Extra column, each `;`-separated item with an `=`
and a value that is neither empty nor `-` is a field, and its value has the percent-encoding of
`=` and `;` undone (`%3D`, `%3B`), the form a tuple file carries. A tuple line's value is taken
as written. Every line of a tuple file is counted, and those whose populated value belongs to no
source of the run are counted apart (`tuples_outside_every_source`).

**Pairing.** Under each key, (Location, Allele, Feature) in the feature projection and (Location,
Allele) in the allele projection, and each field, the two sides' values pair one to one, a value
written more than once as often as each side writes it. Rows that share a (Location, Allele,
Feature) are never collapsed: a tuple counts as often as the rows of its side write it and pairs
with as many of the other side's. The values equal as text pair first; what neither side pairs
that way is the remainder, which the value rule and the classes settle.

**Value rule**, the same for every engine. Two values are equal when their text is equal, or when
both are decimal numbers that differ, in exact decimal arithmetic, by at most 0.000001 of the
larger of the two and by at most 0.000001 outright: both, so that two small numbers are never
equal by their size alone and two whole numbers are equal only when they are the same number. A
pair equal as numbers alone is counted under `text_only_differences`.

**Classes.** `--class-table TSV` declares the documented classes, one per line, tab-separated:
the class (`order_only` or a class of `data/field_classes.tsv`), the plugins it belongs to joined
by `;`, and for each plugin, joined by `;` likewise, the field names it reads there joined by
`,`, a `*` standing for any run of characters. Further columns are carried into the report and
onto the class's lines of the summary, unread but for the fourth, and `#` lines are skipped. A
run that names no table reads the tree's `scripts/concordance/plugin_class_table.tsv`, which is
built into the binary, and its report names that table as `plugin_class_table.tsv`; under a table
with no row no class is declared and the adjusted counts equal the raw ones. `--classes A,B`
names the declared classes the run applies: none for an empty value, and where a run names none
every class but those whose fourth column reads `not met`. A class is applied to the fields of
its own plugins alone. The rules are in `src/rules/`, one per class the tree's table names; a
class of `data/field_classes.tsv` that the tree's table never names has no rule here, and a run
whose table declares one ends in the error line when the class is consulted
(`<class>: no rule of this build`).

What a class reads beside the two rows comes from the run's own files and command line, each named
by an option and listed in the report. Without its file a class accepts no row it reads the file to
accept, and a line of the summary says so where the class is applied on a run that loaded its
plugin; `hgvs_insertion_left_where_the_shift_stops` and
`plugin_value_follows_stop_codon_term_pair` read on without `--json-cache`, as that option's item
says.

- `--gwas-catalog FILE`: the associations file the scored run's GWAS plugin read, for the
  classes of `GWAS`.
- `--revel-table FILE`: the table the run's REVEL plugin read, its tabix index beside it, read at
  a position as a class asks, for the classes of `REVEL`. The index is the first that exists of
  `FILE.csi`, `FILE` with `.csi` for its last extension, `FILE.tbi` and `FILE` with `.tbi` for its
  last extension, the order `tabix` looks in, and is read in the CSI or the TBI layout, whichever
  its bytes name. `--revel-no-match` says the run set the plugin's `no_match`, under which the
  table's transcript ids are not read.
- `--utrannotator-sweep PATH`, repeatable: Ensembl VEP's default-format outputs of the input
  under fixed Perl hash seeds, one run a file, or a directory, of which every file that opens as
  such an output is taken, for the classes of `UTRAnnotator`. One of them reads the runs to
  refuse as well: `utrannotator_lost_stop_written_before_stop_codon` leaves a lost stop the
  engine alone reports where a run writes the lost stop on the row.
- `--json-cache DIR`: the candidate JSON cache of the scored run. The classes of `Downstream` and
  of `Blosum62` read its coding sequences, 3' UTRs and codon tables, and are handed the applied
  classes as well: a value two causes change is set aside only where the second cause's class is
  applied too. The same read keeps what stands at the two coding ends of every transcript: whether
  a stop codon closes its coding sequence and a start codon opens it, as vep-rs reads the record
  (`src/rules/downstream/ends.rs`, held to `tests/fixtures/coding_ends/verdicts.tsv` of the
  repository). By it the two Downstream start classes take a value on a row whose two Consequence
  sets differ as a class of a coding end takes them, and `plugin_value_follows_stop_codon_term_pair`
  takes a value beside such a pair. With the cache that class takes a value beside a pair of a
  stop codon class only on a transcript that shows that class's stop codon, and none on a
  transcript the cache holds no reading of, but for a pair of
  `stop_codon_settled_before_unknown_bases` or of `frameshift_inside_a_held_stop_codon`, which it
  reads by the two rows alone; without it the
  class takes nothing beside a pair of a coding end and reads the pairs of the stop codon classes
  by the two rows alone. The summary's line on a run without `--json-cache` is written for a run
  that loaded `Downstream`.
  The class of `MaxEntScan` reads the cache's introns, and `spliceai_multi_entry_cell`
  the gene of a row that writes no `SYMBOL`.
  `hgvs_insertion_left_where_the_shift_stops`, where the table declares it for
  `HGVSIntronOffset`, reads the exon spans of a row's transcript and by them tells whether the
  two coordinates of a `dup` are those of the first and last base of a run of as many bases as
  the allele holds beside the notation's point; two exonic coordinates of such a run count the
  allele's bases less those of the introns between them. On the row of an insertion of 1,001
  bases that the reference does not walk, the reference's `dup` is the class's only as such a
  run: a row whose transcript the cache does not hold on the row's strand, and a `dup` with one
  coordinate in the 3' UTR and the other before it, stay differences. An engine `dup` of two
  exonic coordinates of one region is the class's where it counts as many bases as the allele
  holds or is such a run. Without the cache the class reads a `dup` by its count alone: that
  reference `dup` is taken from 2 bases up to the allele's own, but for two fewer, whatever
  intron would explain the count, and an engine `dup` around an intron stays a difference.
- `--spliceai-table FILE`, repeatable: a score file the run's SpliceAI plugin read, its `snv=`
  file and its `indel=` file, each with its tabix index beside it, found as the index of
  `--revel-table` is, read through the index at the variant of a row a class asks about, never a
  file whole. Under the files the last lines of the summary count, per engine, the rows
  `spliceai_multi_entry_cell` was asked about and the rows it refused, by the ground of the
  refusal: counts alone, a row once under the first ground met for it. A row whose gene neither
  its `SYMBOL` nor the cache names is not set aside.
- `--input FILE`: the scored run's input VCF, read for the records the rows name (below).
- `--run-command COMMAND`: the scored run's command line. `dbnsfp_mitochondrial_record_named_m`
  reads the assembly its words name, GRCh37 under `--grch37` and else the value of the last
  `--assembly` or `-a`, and sets nothing aside in a run that gives no command line or names no
  assembly. The class takes a field the registry gives dbNSFP among the run's plugins, and a
  rule is handed that plugin with the field. The scorer reads the two outputs and never a
  plugin's data file: of a mitochondrial record it tells the spelling and which output fills
  the field, and neither the name the plugin's table gives the contig, nor the sequence its
  rows count on, nor whether the run gave the engine a reference sequence.

**Adjusted counts.** A tuple left unpaired is set aside on both sides when an applied class of
its own plugin reads its field and the class's rule accepts it, the first to accept in table
order credited under `masked`. A rule reads the row, so it is consulted only with a reference row
under the (Location, Allele, Feature) that wrote the reference's value, or wrote no value of the
field where the reference has none to pair, and an engine row that did the same on its side; a
tuple file states no rows, so there the engine's side of the row is the tuples it wrote under
that key. A side that wrote no value of the field is handed to a rule as absent, with one
exception: the classes that read a transcript's sequence (`#transcript_reading` of
`data/field_classes.tsv`: the two Downstream classes and `plugin_hgvs_shifted_cds_coordinates`)
read an empty value as a value, the plugin's own where the first codon it reads is a stop, so
they are handed the empty text where the side's row writes the field with an empty value. Such a
value is still no tuple and moves no raw count.

Rows that share a key each stand for one unpaired tuple of the field, in the order their input
writes them: the nth unpaired copy of a value reads the nth row that wrote it, and the nth value
left with none to pair reads the nth row of the other side that wrote no value of the field. A
value left over where the other side has no such row left, as under a key two input records
share, is put to no class and stays a residual. Under a key that rows of the reference or of a
default-format engine share, the values left pair by input record, which a row names in its
`Uploaded_variation`, before they pair in their order: a value pairs with the other side's value
whose row carries the same `Uploaded_variation`, else stands against that record's row of the
other side that wrote no value of the field, and the values no record places this way pair in
their order with the rows still free. The pairs equal as text or as numbers are counted before
any of this and move with none of it.

A rule is handed the reference's rows of the same input record on the Feature: those at the
Location under the same `Uploaded_variation`, where it spells the record's alleles. Under an
identifier the rows state neither that they are one record's nor its alleles in order, and
`--input` settles it: of the input's lines that carry the rows' identifier, the one line, or of
several the one that writes the row at its Location with its allele, gives the record's rows
with an entry of the input line. A line carries the first `;`-joined element of its ID column,
the name both engines print for its rows, and none where that element is a dot, empty or `0`.
Without `--input`, where no line carries the identifier, where none or more than one of several
lines writes the row, and where an allele repeats among the rows, which rows of more than one
record do, nothing is handed, and a class that needs the record's rows sets nothing aside. A
sibling-allele class (`#sibling_allele`) is handed the engine's rows of the same input record on
the Feature as well, found among a default-format engine's rows at the Location by the same
rule; a tuple file states no rows, so for it such a class sets nothing aside. A class that reads
the reference bases a row's allele stands in for (`#record_reference`, and the classes of
`#transcript_reading`, which ask on a row whose two Consequence sets differ whether a class of a
coding end takes the pair) takes them from a row
whose `Uploaded_variation` spells the record's alleles, and for a row under an identifier from
the lines of `--input` that carry the identifier, read by position; without `--input`, and where
no such line or two that differ give the row at its Location, the class reads the row without
its record.

One class is decided without a reference row, since its shape is that the reference wrote none:
`plugin_value_on_reference_skipped_record` takes a value a default-format engine alone wrote at
a Location of a record the reference dropped, told from the two sides' rows at the Location and
the reference's first column, with the gVCF variant sites its rows' names and the lines of
`--input` show left to their own class. A difference no applied class accepts stays in the raw
and the adjusted counts alike. A declared class that is not applied is consulted after the
applied ones: an entry in which it would set a tuple aside names it under `not_applied` and
carries no `adj_` count.

**Projections**, both reported for every engine. `feature` is the tuple above. `allele` drops
the Feature on both sides and keeps a value once: (Location, Allele, field name, value). An
allele tuple left unpaired is set aside when every unpaired feature tuple it projects from was,
at least one, and is credited to the first of their classes in table order. The projected values
pair anew, so an allele tuple can stay unpaired although its feature tuples all paired; such a
tuple stays in the raw and the adjusted counts alike.

**Report.** `--report JSON` is written whole or not at all: staged beside its path and moved
onto it. `scores` holds one entry per source, engine and projection, integers only, under the
count key names of the concordance provenance records: `perl` and `rust` (the reference's and
the engine's tuples), `intersection`, `adj_perl`, `adj_rust`, `adj_intersection`, `masked` (per
applied class of the source, the keys under which it set a field's tuple aside, zero included),
`not_applied`, `text_only_differences`, `residual_rows` and `shared_keys`, the keys of the
projection under which the reference or the engine writes one field of the source more than
once: in `feature` the keys that rows share, in `allele` those with more than one value of a
field. No F1 is stored: raw F1 is 2 x intersection / (perl + rust), adjusted F1 the same over
the `adj_` counts. `class_table` names the table and its `sha256`, `gwas_catalog` the catalog
and its `sha256`, `revel_table` the REVEL table, its `sha256` and whether the run set
`no_match`, `utrannotator_sweep` each output read with its `sha256`, `run_input` the input of
`--input` with its `sha256` and its number of records, `json_cache` the directory `--json-cache`
gave and `spliceai_table` the files `--spliceai-table` gave, by path alone. `spliceai_refusals`
lists, under those files, the engines whose rows `spliceai_multi_entry_cell` was asked about,
each with the rows `asked` about, the rows `refused` and the refused rows of each of the
`grounds`. `classes` lists the declared classes with their plugins, field names, further columns
and whether the run applied them, and `not_populated` each source the reference wrote no field
of. Each input is listed with its `sha256` and its `unclaimed_keys`, a tuple file with
`tuples_outside_every_source`, and an engine with the `binary_sha256` and the `commit` that
`--engine-binary` and `--engine-commit` give it; `run_command` holds the words `--run-command`
gave.

The summary on standard output prints each entry's counts and one line per class of it with the
class's further columns, marking a class that set nothing aside and one whose count is above the
population `--class-population CLASS=N` gives it, the keys the class can reach in the scored
input; it names the unclaimed keys and the sources not populated. A character the stream's
encoding lacks is written as its backslash escape; the encoding is the locale's, UTF-8 under the
`C` and `POSIX` locales.

**Residual files.** `residual_rows` counts the tuples of either side that are outside the
intersection and that no applied class set aside. `--residuals DIR` writes them, never
truncated, to `DIR/<source>.<engine>.<projection>.tsv`: one line per tuple, the side
(`reference` or `engine`) and then the tuple's columns, five under `feature` and four under
`allele`, tab-separated, in the order `sort` gives the tuples under C collation. Every scored
source has its file, empty when nothing is left, and `residual_rows` is the file's line count.
A run is refused whose `DIR` holds a `.tsv` file of an earlier run, or in which two scored pairs
of a source and an engine would share a file name.

**Refusals and exit codes.** Exit 0 with the report written and the summary printed. A refusal
is one `ERROR: [score_plugin_tuples] <message>` line on standard error and exit 2: a wrong
value of an option, an input that cannot be read, a line with the wrong number of columns (the
first such line of the file is the one named), a column header that names other than 14 columns,
a reference without a data row, too many distinct unclaimed keys, a class table that cannot be
read, a registry refusal, a class with no rule. A command line `argparse` would refuse is exit 2
with the usage line and `vep-plugin-pairs score: error: <message>`. An interrupt or a
termination signal ends the run in the error line with exit 130; a repeat of either is ignored
while the run stops. A run that ends in the error line, its report unwritten, keeps none of its
residual files and leaves no scratch. A summary that cannot be written to standard output ends
in exit 120.

**Threads, scratch and memory.** The run uses `--pairing-threads` threads, every core where none
is named, and every thread count writes the same bytes. Each input is extracted in blocks on
every thread into shards by Location, so that equal Locations stay together, one shard per 128
MiB of the longest input; each thread sorts and pairs one shard of each side at a time. The
remainder is settled on the same threads in ranges of whole Locations, about 1 MiB of each of
its sorted files a range: a thread holds one range's rows, verdict lines and residual lines at a
time, and what the ranges leave is joined in their order, so a run that a class refuses ends in
the refusal of the first Location that refuses. The shards of a chromosome of `--json-cache` are
read when a class first asks for it, on the run's threads and at most 16 at a time. The scratch
is a directory of the run's own under `--work-dir`, or under the system's temporary directory:
it holds one record per row with a claimed field, the unpacked copy of a gzip input, the
remainder as sorted tuple lines, when a class is declared both sides' rows at the Locations of
the remainder, and the sorted verdict lines of each range while it is settled, and it is removed
when the run ends. The named lines of `--input` are sorted by identifier there and searched in
place. The Locations of the feature remainder are held in memory, so a run is sized for outputs
that mostly agree.

`VEP_PLUGIN_PAIRS_TIMINGS` in the environment prints the time each pass ended on standard error.

**Bounds.** A whole number a rule reads (a coordinate, an offset, a count) is read up to 128
bits; a longer one is no number to the rule that reads it, and the HGVSIntronOffset rules, which
add such numbers, refuse one above 10^37. A Location coordinate of more than 4,300 digits stays
as written. A `--class-population` count is at most 2^64-1, and `--pairing-threads` at least 1.
A JSON object in a value that is nested deeper than 497 levels is refused. The command line is
read as UTF-8 under every locale. A damaged gzip stream is refused in the words of CPython's
`gzip` for a stream cut short and a wrong magic, trailer or method, and as `Error -3 while
decompressing data: invalid deflate data` for corrupt deflate data. Deflate data that copies
from before its own start is seen at the stream's check sum alone and refused in the check
sum's words or, inside the column header or the first data line, by that line's column count;
an output of 64 KiB or more damaged so in a `--utrannotator-sweep` directory ends the run.

## Contract of `compare`

**Inputs.** One file pair (`--perl-output`, `--rust-output`) or every pair of two directories
(`--perl-dir`, `--rust-dir`). A file whose name ends in `.vcf` or `.vcf.gz`, in any case, is
read as VCF and any other as the default format; a name ending in `.gz` is read as gzip, every
member of the file, by the name alone and never by the bytes. A directory lists its `*.txt`,
`*.vep.vcf`, `*.vep.vcf.gz`, `*.vcf` and `*.vcf.gz` files, and a file pairs with the one whose
name is the same without its ending and without a leading `perl_` or `rust_`; of two files of
one directory with one such name the one listed last within the last matching ending is taken.
Text is read as UTF-8 with replacement, and a line ends at a line feed or a carriage return.

**Records.** `--plugins A,B` names registry plugins, `all` every plugin of the registry; an
unknown name is refused. A record is the fields of one plugin under one key, (Location, Allele,
Feature). In a default-format output a data line is stripped of white space at both ends and
split at tabs, a line of fewer than 14 columns is skipped, and the fields are the `key=value`
items of the Extra column that an entry of the plugin's registry list claims; an entry with a
placeholder (`<word>`, `<N>`) claims every name of that shape. In a VCF output the CSQ field
names come from the header's `##INFO=<ID=CSQ` line, the key is `CHROM:POS` with the entry's
Allele and Feature, each value has its percent-encoding undone, and a value that is empty or
`.` is no field. A later row under a key replaces an earlier one. A field name two requested
plugins write counts for the first of them in registry order and is reported for both under
`shared_fields`. `--fields A,B` replaces the registry's entries where the run names exactly one
plugin.

**Comparison.** Per plugin and file pair, over the keys either output annotates in sorted order:
a key both annotate is a `match` when every field of either side has the same text on both, or
two numbers less than `--tolerance` apart (default 1e-06), and a `mismatch` otherwise; a key one
side alone annotates is `perl_only` or `rust_only`. Two texts of `5UTR_annotation` count as the
same text where they hold the same groups of the same `key=value` pairs with the groups, or the
pairs inside one, in another order and vep-rs's pairs stand in the order of their keys: Ensembl
VEP's UTRAnnotator writes both in the order of a hash, and the rule is the one of the scorer's
class `utrannotator_annotation_pair_order`. The concordance is match over annotated keys, and
none where neither engine annotated a record: a plugin that annotated nothing is not measured.

**Report and verdict.** Standard output carries, per plugin, each file pair's fields, record
counts, the four counts, the concordance and the first three mismatches, then the aggregate and
`STATUS: PASS` or `STATUS: FAIL` with its reasons: a value mismatch, a record on one side alone,
or an engine that annotated fewer than `--min-annotated` records (default 1). `--json-output
PATH` writes the plugins, the tolerance, the floor, the file count, per file pair and plugin the
fields compared, the counts and the first `--max-mismatches` mismatches (default 10; 0 keeps
none; a negative value is refused) with `total_mismatches`, and per plugin the aggregate with
its `pass`. Exit 0 when every plugin passes, 1 when one fails, when the directories hold no
matching pair, when neither input mode is named and when an input cannot be read
(`ERROR: [vep-plugin-pairs] <path>: <reason>`), 2 for an unknown plugin and for a command line
`argparse` would refuse (the usage line and `vep-plugin-pairs compare: error: <message>`), 130
when an interrupt or a termination signal ended the run. The summary and the error lines are
written as UTF-8 under every locale, and a standard stream that takes no line changes neither
the exit status nor the JSON.

**Threads, scratch and memory.** The run uses `--threads N` threads, every core where none is
named, and every thread count writes the same report. Of several file pairs as many are
compared at a time as there are threads, which are shared out among them. A plain output is
read in ranges of `--block-bytes` bytes (default 64 MiB), a thread a range; a gzip output and a
pipe are read as one stream on one thread, the other output on the rest. What a row gives the
run's plugins goes to the scratch as one record (the place of its line in the output, the key
and, per plugin, the fields' numbers and values) in the shard that the key's hash names, the
same for both sides; a pair has one shard per `--shard-bytes` bytes (default 128 MiB) of its
longer output, 4,096 at most. Each thread then reads one shard of both sides at a time, sorts
its records by key and place, takes per key and plugin the fields of the last row and walks
the two sides together in the keys' order.

A reading thread holds one line, its buffers for the shards (about 16 MiB together at most) and
the plugins of up to 100,000 field names. A comparing thread holds one shard of each side:
about `--shard-bytes` of each output at most where the keys spread evenly, whatever the outputs'
length. The pair keeps, per plugin, its counts, the names of the fields met and each shard's
first `--max-mismatches` mismatches, of which the first in the keys' order are reported. The
scratch is a directory of the run's own under `--work-dir`, or under the system's temporary
directory: it holds the two sides' shards of the pairs in work and is removed when the run
ends, on an error and on a stop signal as well; such a signal ends the run with
`ERROR: [vep-plugin-pairs] interrupted before the report was written`. `--threads`,
`--block-bytes`, `--shard-bytes` and `--work-dir` are taken by their whole names alone.
`VEP_PLUGIN_PAIRS_TIMINGS` in the environment prints on standard error when a pair was read and
when it was compared.

## Layout

- `src/main.rs`: the first word of the command line, which names the command.
- `data/plugin_fields.tsv` and `src/registry.rs`: the plugin registry, with its aliases and the
  plugin columns that carry the name of one of Ensembl VEP's own keys.
- `data/field_classes.tsv` and `src/rules/`: the fields each class reads and the class rules, a
  file per plugin or family of classes; `entry.rs` and `record.rs` hold an output row and the
  rows of an input record as the rules read them, `shared.rs`, `cache.rs` and `tabix.rs` what
  several plugins' rules read, `exon_spans.rs` the exon spans a class of `hgvs_offsets.rs` reads
  from the cache. `tools/vep-formats-compare` holds the rules of the same classes,
  with the fields of each in its `data/classes.tsv`; `tests/twin_class_data.rs` holds the two
  data files to each other for the classes both carry, so a change to such a class edits both
  tools.
- `../vep-formats-compare/src/coding_ends.rs` with `coding_end_terms.rs`,
  `../vep-formats-compare/src/rules/coding_end.rs` and
  `../vep-formats-compare/src/genetic_codes.rs`: source files of that tool, compiled into this
  crate as `rules::coding_ends`, `coding_end` and `genetic_codes`, so that the coding ends of a
  cache transcript, the terms derived on a row there and the shapes of
  `alternate_codon_past_the_coding_end` have one definition for both tools.
- `src/score/`: the `score` command. `cli.rs` the command line; `mod.rs` the run in its order;
  `pass.rs` the passes; `sources.rs` and `class_table.rs` the run's sources and classes;
  `reads.rs` and `input_index.rs` the files a class reads; `number.rs` the value rule;
  `settle.rs`, `merge.rs`, `projection.rs`, `context.rs` and `tally.rs` the settling of the
  remainder; `residuals.rs`, `report.rs` and `summary.rs` the outputs; `interrupt.rs` the stop
  signals.
- `src/input.rs`, `src/extract.rs`, `src/pair.rs`, `src/rows.rs`, `src/resolve.rs` and
  `src/order.rs`: the passes over the inputs. An input's digest, plain copy, kind and blocks;
  the extraction rules and the shard records; the pairing of one shard, the counts and the
  remainder; the rows of the remainder; the fields met and what settles a new one; the byte
  rules every pass shares.
- `src/slots.rs`: the threads of the passes that work at one time, the caller's and those that
  take a plain input's digest among them.
- `src/compare/`: the `compare` command. `cli.rs` the command line, `files.rs` the file pairs,
  `read.rs` the readers of the two formats, `fields.rs` the plugins, their entries and the
  field names met, `shards.rs` the scratch and its records, `tally.rs` the comparison of a
  shard and the joining of the shards, `report.rs` the printed and the JSON report.
- `src/pytext.rs`, `src/score/text.rs`, `src/score/paths.rs`, `src/score/shlex.rs` and
  `src/compare/text.rs`: the CPython text, path and number rules the reports and the messages
  are defined in.
- `tests/`: `pairs.rs` runs the passes over small outputs, one test per extraction rule;
  `commands.rs` the first word of the command line; `score_*.rs`, `settle_threads.rs` and
  `compare.rs` the two commands; the other files hold the class rules to made-up rows, tables
  and caches, a file per family.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-plugin-pairs/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

CI runs the test, clippy and fmt commands above in its Test, Clippy and Format jobs, builds the
package on the declared `rust-version` in MSRV and audits its lockfile. After the Test job's
build it runs the plugin concordance driver's tests that score and the driver itself on one
committed corpus through the binary, the printed rows held to the committed ones.

## Harness

- `scripts/concordance/reproduce_plugin_concordance.py`, the driver of `plugin_gate_runs.toml`,
  scores each run with `vep-plugin-pairs score` under the tree's class table. It runs the file
  `--pairing-tool` names, else this tree's release build
  (`tools/vep-plugin-pairs/target/release/vep-plugin-pairs`), else its debug build, and refuses
  a run that scores with none built, naming the build command.
- `scripts/concordance/run_concordance.sh` runs `vep-plugin-pairs compare` over its two output
  directories under `--plugins`. Unset or empty, `VEP_PLUGIN_COMPARATOR` builds this tool through
  cargo before any engine runs; any other value names a built binary, as a path or a command on
  `PATH`. A build that fails, or a value that names no executable, stops the script before its
  first engine run.
