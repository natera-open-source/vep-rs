# vep-formats-compare

The format comparator of the comparison oracle: one Ensembl VEP output against a candidate's in
one output format, every row, CSQ entry and JSON object, every field verbatim, with the
documented divergence classes a comparison names set aside and counted. Two further modes hold
one engine's formats against each other and two engines' warnings files against each other.

```text
vep-formats-compare --format {default,tab,vcf,json,parquet,recoder} --expected FILE \
    --actual FILE [--input VCF] [--manifest JSON] [--masks CLASSES] \
    [--reference-release {115.2,116.2,116.3}] [--report JSON] [options]
vep-formats-compare --self-parity --default FILE --tab FILE --vcf FILE --json FILE \
    [--parquet PATH] [--input VCF] [--report JSON]
vep-formats-compare --warnings --expected FILE --actual FILE [--expected-skipped FILE] \
    [--actual-skipped FILE] [--masks CLASSES] [--reference-release {115.2,116.2,116.3}] \
    [--report JSON]
vep-formats-compare field-class --class CLASS --questions FILE --answers FILE \
    [--spliceai-table FILE]... [--reference-fasta FILE] [--json-cache DIR]
```

`--help` lists every option with what it reads: the manifest's run (`--run`,
`--reference-sequence`), `--keep-only`, `--floors`, `--field-pairs`, `--unordered`, `--examples`
and the files a class reads beside the two outputs (`--nearest-sites`, `--json-cache`,
`--reference-fasta`, `--gwas-catalog`, `--revel-table` with `--revel-no-match`,
`--spliceai-table`, `--utrannotator-sweep`). An option is accepted as `--flag value`,
`--flag=value` and an unambiguous prefix.

## Contract

A run writes one JSON report to `--report`, when the option names a file, and then its verdict on
the standard output: `<format or mode>: green` or `red`, followed by one indented line a reason.
Nothing is written, to the report or to a standard stream, before the whole report exists. The
report is ASCII with one space of indentation a level, its members in the order they were set.
Exit codes: 0 for a green verdict; 1 for a red one; 2 for a refused command line, an unreadable
or malformed input and a report that cannot be written, each with one line on the standard error
(`ERROR: [compare_vep_formats] <detail>`, or the usage text and the parser's message). A
malformed input that no reading of the comparison expects ends the run with the same line and
exit code 1, and no report.

An entry is one row of the default or tab format, one CSQ entry or record line of the VCF format,
or one consequence object or record object of the JSON format. Its key is the Location (default,
tab) or the input record's ordinal (VCF, JSON), the allele, the feature and the feature type, with
a counter on exact duplicates, and the entries of one key are compared field by field. A VCF line
names its input record by its first five columns, and a JSON object by the first five columns of
its `input`. A line or object that names no line of `--input` is keyed by its own ordinal in its
output, counted over the data lines of a VCF output and over every line of a JSON output. The VCF
format of an input that is not a VCF writes such lines. A manifest item is looked up by an
entry's Location, then by its input record's ordinal. An item whose candidate sets hold the
reference's own set beside another documents a key several records share: a row of that key
on which the candidate writes the reference's set is compared like any other. A
`field_divergences` item that states `rows_alike` documents a field on such a key and consequence
set: a row of it on which both outputs print one value of the field is compared like any other,
a row whose two values differ is held to the item's `vep_rs`, and an item without the key holds
every row of its key. The key's number is the count of the rows alike, taken over the rows of the
item's Location in the default and tab formats and over the entries of its records in the VCF and
JSON formats: an item whose count is another number is a stale line,
`<key>: field <name> documents <n> rows alike, <m> are`, its key the item's Location, Allele
and Feature between `|`. An item whose key no compared row reaches with the field is held to no
count. Where `--input`
is an input no run of the
manifest names, the lists number another input's records and no item is looked up by ordinal.
Where the outputs hold an entry and none has a Location or an ordinal an item is looked up by, no
item meets a row and the items are left out of `masks.applied`. A built-in class then takes a row
from the row's own cells and counts it once, and a class the manifest alone carries sets nothing
aside. A reason,
and any reason is a red verdict, is an entry on one side only, a field that differs after the
masks and field pairs, a divergence the manifest documents that no entry shows, the shared
entries standing in another order (unless `--unordered`), and a floor that is not met. A class
named under `--masks` is applied from the corpus manifest (`--manifest`, whose run `--input` and
`--run` select) or, for a built-in class, from the rows alone and the files the class reads; one
built-in class, `dbnsfp_mitochondrial_record_named_m`, reads the manifest as well, for its assembly
and the columns its run hands dbNSFP. Every entry or field set aside is counted in the report
under `masks.applied`, and a requested class that set nothing aside is listed under
`masks.unused`.

**The reference release.** `--reference-release` names the Ensembl VEP release that wrote the
reference output, `115.2`, `116.2` or `116.3`, and `116.3` without the option. A class whose
defect Ensembl VEP corrects is applied under a release before the one that corrects it and under
no other. Under the correcting release and every later one the class is taken out of the classes
`--masks` names before the manifest or a rule reads them: a row or field it would set aside is a
difference, a manifest's items under it are those of a class `--masks` does not name, and the
report names the class under `masks.corrected` with the release that corrects its defect, a
member the report of a run that names no such class lacks. `data/classes.tsv` states that release
in its column `CORRECTED_IN`: `116.3` for `frequency_filter_zero_value` and
`frequency_filter_allele_key_order`. The warnings mode takes the option beside `--masks`, and its
report lists the classes it applied alone.

**The two outputs are streams.** In the default, tab, VCF and JSON formats each output is read
twice, in blocks cut at line ends, and never held whole. A record group is every data line of one
side that shares a key: a Location in the default and tab formats, the first five columns of a
VCF line, the first five columns of a JSON object's `input` text. The first reading keeps one
small row a run of lines of one group and one 128-bit digest a group. A group whose lines give
the same entries in the same order on both sides, and that nothing forces, is kept back: it adds
its entries to the counts, its fields to the populated columns and its place to the order of the
shared keys, and its lines are never parsed into entries. The same entries means the same text
in every format but JSON, whose objects are read member by member, so that two lines whose
members stand in another order are one form, and whose numbers are read as values, so that a
fraction or an exponent spelt two ways is one form while a whole number against the same value
with a fraction is two. The second reading counts the populated fields of the kept-back lines
and hands every other group, in file order, to the comparison of entries.

Three readings force a group whatever its text: a key the manifest documents (its staleness is
checked on equal values too), `regulatory_feature_rows` (reference rows are dropped before the
comparison) and `conservation_short_name_spelling` (the reference text is rewritten). In the VCF
and JSON formats the lines of one key name the input's lines of that key in order; a group of
more lines than the input holds of its key is handed over. A digest only routes: where the
kept-back lines stand in one order on both sides their text is held against the other side's as
it streams by, a group whose text differs after all is handed over, and one line on the standard
error says how often the outputs were read again for it. Under `--unordered` the kept-back lines
of two outputs whose groups stand in different orders pass each other, and the digests alone
settle those groups.

A Parquet comparison is read through the same stream. Each side is written as its canonical
form (**A Parquet output**, below) with the rows in the side's own order, a few thousand rows
at a time, and the two texts are read as two tab outputs whose order is not read: a record group
that both sides hold line for line is never held or parsed, whatever its place in either side.
The groups handed over are put in the byte order of their lines before they are parsed, the
order the sorted canonical form has them in, so the report is the one the whole outputs give.
Each side is written once a reading, by one thread, and no file is written.

The outputs are compared whole in the recoder format, under `--keep-only`, under
`individual_allele_order`, under a class that reads across records
(`reference_skipped_record`, `breakend_mate_local_read`, `bnd_own_chromosome_rows_lost`) in the
default and tab formats, under a JSON floor on a field whose
name holds a dot or on a column no single line holds, and in the two other modes. A Parquet
comparison compares them whole as well where the reference is refused, where a tab reference's
or a tab candidate's column line follows a data line, where a row of either side holds a line
end, and where a row of a Parquet side opens with `#`.
In the VCF and JSON formats those three classes hand over the groups of the input records
they can touch and stream the rest: a record on a skip route, a breakend (a bracket, `<BND>`,
a single breakend or `SVTYPE=BND`), and every record whose span lies within two bases of one
of those.
`VEP_FORMATS_COMPARE_WHOLE=1` compares the whole outputs in every format. Both paths write the
same report, the same standard output and the same exit code, but for the one input **Bounds**
names.

**The input.** `--input` is read as a stream, plain or gzip. Of each data line the comparison
keeps the part it reads: the line up to its fifth column, lengthened to the end of its fifth
white-space field where that field ends later, and the whole line where it has fewer than five
columns. The structural-variant classes read the file themselves, INFO column included, and a
gzip file as GNU gzip's `zcat` writes it: up to the first zero byte after a member, and no
member whose header checksum is not the header's or whose header carries the encryption flag or
a reserved flag. A comparison under such a class whose two readings of the input count
different numbers of records is refused with exit code 2.

**Floors.** `--floors COL=MIN,...` asks that the populated count of column `COL`, its entries
with a value other than `-`, reach `MIN` on both sides. `COL` is the default-format name in
every format. In the JSON format the count is taken under the column's JSON key: the name in
lower case, renamed as the JSON output renames it, SIFT and PolyPhen under either of their two
keys and a position under either of its start and end keys. A column outside Ensembl VEP's own
that no entry of either side carries under that key is counted under the key of the same name
inside a plugin's object. A column with no row count in the JSON format has no floor there and
is listed under `floors_not_applicable`: a record-level column, `Feature_type`, `Extra` and a
custom source's columns. In the tab and VCF formats a floor on a column that neither header
declares is not applicable in the same way, and a column one header declares and the other
does not stays a difference of the field comparison. Under `conservation_short_name_spelling`
a floor spelt with the misspelt key is read under the renamed one.

**A Parquet output**, a file, a directory of part files or a directory of `chrom=<c>/`
directories of them, is read column by column and written as tab rows: the key columns lead
every file and the columns after them are the tab columns, one the writer stored behind `tab_`
(a name an earlier column carries in any letter case) under the name the file's footer states
under `tab_columns`, a NULL cell `-`, a string as stored,
an integer in decimal, a double as the shortest decimal that reads back, a list `,`-joined, a
struct of positions as `start`, `start-end` or `-`, then `/total`. Both sides are brought to the
canonical form (every number in its shortest text) and compared as a multiset of rows; the
comparison of the whole outputs sorts the rows of each side.
A column the Parquet output stores as strings, as it stores a field a plugin or a
`--custom` source declares, is no number column whatever its name: its cells are compared as
stored and as the tab output writes them. A cell of a number column that is no number of its
type stays as written, so a NULL in its place is a difference. Pages are read uncompressed and under snappy, gzip, brotli, LZ4 and ZSTD; a file under
LZO, a stored type outside those cells (binary, date, time, timestamp, decimal, a list of lists
or of structs, a map) and a damaged file are the comparator's error.

**Text.** An output, an input or a side table whose name ends in `.gz` is read through gzip,
every member of the file, zero bytes after a member being padding. Text is UTF-8, and a line ends at `\n`, `\r\n` or a lone `\r`. An
argument that is not UTF-8 names its file by its bytes and is printed with each such byte
escaped.

**Bounds.** A whole number is held in 64 bits, in 128 by the structural-variant classes: by
where it is read, a longer one is the comparator's error, the end of the run, no number, or the
nearest value held, and a function that reads only a number's digits, its sign or its being
zero stays exact at any length. The lines of a kept-back group are never parsed into entries, so
a longer whole number in a group that both outputs write alike can be passed over as the outputs
stream by, where the comparison of the whole outputs is refused; a VCF line whose position is
longer than 18 characters is handed over, and both paths refuse it. A JSON
string that holds half a surrogate pair is refused, as is a JSON value nested deeper than 81,578
levels. The usage and the help are laid out at 80 columns under the program's own name, whatever
the terminal's width and colours, and a refused command line is worded as CPython 3.14's
`argparse` words it.

**The class data file.** `data/classes.tsv` is the tool's class tables as one tab-separated file,
which a program outside the tool reads by its path to hold its own lists of class names to the
tool's. The first line names the columns; every other line is one name and has as many cells as
the first line. The first column, `class`, is the name, a class's but for one row, and a name
stands on one line.
Every other column but the last is one table of the comparison under the table's name:
`FIELD_CLASSES`,
`REGISTRY_CLASSES`, `TERM_PAIR_CLASSES`, `RECODER_MASKS`, `PLUGIN_WARNING_CLASSES`,
`WARNING_TEXT_CLASSES`, `INTRON_GAP_CLASSES`, `DOWNSTREAM_CLASSES`, `TRANSCRIPT_READING_CLASSES`,
`STRUCTURAL_VARIANT_CLASSES`, `CACHE_READING_CLASSES`, `MANIFEST_CLASSES`,
`RECORD_DROPPING_CLASSES`, `RECORD_READING_CLASSES`, `ALLELE_READING_CLASSES`,
`SIBLING_ALLELE_CLASSES`, `TERM_PAIR_FOLLOWING_CLASSES`, `RECORD_REFERENCE_CLASSES`,
`FRONT_END_WHOLE_CLASSES` and `RECODER_CORPUS_MASKS`. A cell of a table's column is empty where
the class is not in the table and otherwise holds the class's place in it, counted from 0 in the
table's own order. A table that gives each of its classes a value has that value in the column
after its own: `FIELD_CLASSES.values`, the field-name patterns the class reads, and
`RECODER_MASKS.values`, the prefixes the class drops, each comma-joined, and
`PLUGIN_WARNING_CLASSES.value`, the warning's text. The classes of the corrected behaviours are
the rows with a place in `REGISTRY_CLASSES`, `FIELD_CLASSES`, `PLUGIN_WARNING_CLASSES` or
`MANIFEST_CLASSES`; the classes the `masks` of a recoder corpus may carry are the rows with a
place in `RECODER_CORPUS_MASKS`, the one column of places that is no table of the comparison. A
row with no place in any table is a class the comparison names outside its tables, but for one:
`second_expected_divergence` is no class, it is the manifest member under which an item names
its second class.
The last column, `CORRECTED_IN`, holds no place: its cell is the first Ensembl VEP release that
corrects the defect the row's class sets aside, and is empty for a class no release corrects.
`scripts/concordance/feature_oracle.sh` lists the built-in classes and the classes that read the
input, and its test holds both lists to this file. Its list of built-in classes leaves out
`dbnsfp_mitochondrial_record_named_m`, which accepts nothing without a manifest, so a case that
names the class needs a manifest.json beside its input.
`tools/vep-plugin-pairs` holds the rules of the classes
`scripts/concordance/plugin_class_table.tsv` names, with the fields of each class in its
`data/field_classes.tsv`; `tools/vep-plugin-pairs/tests/twin_class_data.rs` holds that file to
this one for the classes both carry, so a change to such a class edits both tools.

**The cache class file.** `data/json_cache_classes.tsv` names the classes whose rules read the
candidate's JSON cache: one line a class, with the side table a comparison builds from
`--json-cache` where `--masks` names the class (`structural_variants`, `sequence_edges`,
`partial_codon_models`, `length_edits`, `codon_tables`, `downstream_sequences`, `transcript_edges`,
`intron_ends` or `splice_windows`). The command line decides from this file which of the last
eight it builds, `sequence_edges` beside `--reference-fasta` alone; for the classes of
`TRANSCRIPT_READING_CLASSES`, which the file lists under `downstream_sequences`, it builds
`partial_codon_models` as well. The structural-variant
comparator is built where the masks name a class of `STRUCTURAL_VARIANT_CLASSES`, and it
refuses a run without a cache for the two classes this file lists under `structural_variants`,
which are the `CACHE_READING_CLASSES` column of the class data file. A comparison of the
recoder format builds `sequence_edges`, `transcript_edges` and `intron_ends` alone.
`scripts/concordance/feature_oracle.sh` lists the classes it hands the cache, and its test
holds those lists to this file. Two classes for which a comparison opens the cache have no
line: `spliceai_multi_entry_cell`, whose score table reads the cache on `--spliceai-table`
whatever the masks name, and `reference_skipped_record`, for which the structural-variant
comparator opens a cache it is given and asks nothing of it. A class whose rule reads the
cache gets its line here and its name in the runner's list in one change;
`tests/json_cache_classes.rs` fails on a class the command line opens the cache for and the
file lacks.

**The `field-class` command.** A command line whose first argument is `field-class` asks one field
class about differing fields, for a program that classifies the differences of two default-format
outputs itself: `vep-formats-compare field-class --class CLASS --questions IN --answers OUT`, with
the options the class reads. Anywhere but first, `field-class` is an argument of the comparison.
`IN` holds one question a line, the columns separated by tabs: the row's Location, Allele and
Feature; the name of the differing field; the reference's value and the candidate's, `-` for a
row that lacks the field; then the cells of the reference's row the class reads and the same
cells of the candidate's row, each field with its text as the default format prints it and `-`
for a field the row lacks. The row's consequence terms are no column: neither class reads them.
`OUT` takes one line a question, in the questions' order: `1` where the class accepts the field,
`0` where it does not, the answer of the comparison's own predicate of the class over the two
rows. The command answers for two classes.

- `--class spliceai_multi_entry_cell --spliceai-table FILE [--spliceai-table FILE] [--json-cache
  DIR]` reads questions of 30 columns, twelve fields a row: `SpliceAI_pred`,
  `SpliceAI_pred_SYMBOL`, `SpliceAI_pred_DS_AG`, `SpliceAI_pred_DS_AL`, `SpliceAI_pred_DS_DG`,
  `SpliceAI_pred_DS_DL`, `SpliceAI_pred_DP_AG`, `SpliceAI_pred_DP_AL`, `SpliceAI_pred_DP_DG`,
  `SpliceAI_pred_DP_DL`, `SpliceAI_cutoff`, `SYMBOL`. Each `--spliceai-table` names a score file
  the run's SpliceAI plugin read (its `snv=` file, its `indel=` file), a bgzipped VCF with its
  tabix index beside it at `FILE.tbi`, the index the plugin demands of a file it reads (a `.csi`
  index alone is refused); `--json-cache` names the candidate JSON cache, whose transcripts give
  the gene of a row that states no `SYMBOL`, and is given where the cache holds a `transcripts`
  directory.
- `--class hgvs_shift_sequence_edge_window --reference-fasta FILE --json-cache DIR` reads
  questions of 24 columns, nine cells a row: the row's `Feature_type`, then its `STRAND`, `HGVSc`,
  `HGVSp`, `HGVS_OFFSET`, `SHIFT_LENGTH`, `cDNA_position`, `CDS_position` and
  `Protein_position`. `--reference-fasta` names the run's reference sequence, an
  uncompressed FASTA beside its `.fai` index, and `--json-cache` the candidate JSON cache, whose
  transcripts the class reads; a command line of the class names both. The answer is the one a
  comparison under `--masks hgvs_shift_sequence_edge_window` gives a default-format row beside
  the same two files, the row's span read from its Location.

Every question is answered before `OUT` is opened, and nothing is written to the standard output.
Exit codes: 0 once `OUT` is written; 2, with `OUT` not written, on a usage error (the usage block
and `vep-formats-compare field-class: error: ...` on stderr; another class than the two the
command answers for is an invalid choice of `--class`, and an option the class reads that the
command line lacks is a required argument), on a score file, a reference sequence, a cache or a
questions file that cannot be read (`ERROR: [compare_vep_formats] --spliceai-table: ...`,
`--reference-fasta: ...`, `--json-cache: ...`, `cannot read ...`), on a question of another number
of columns (`--questions: IN line N has M columns, a question of CLASS has 30`, or `has 16`) and
on an `OUT` that cannot be written; 1, with `OUT` not written, where a comparison would end on the
same file with `ERROR: [compare_vep_formats] <what>` and exit code 1.

**The case table.** `tests/data/builder_parity/cases.tsv` states the answers of the items of the
comparison that the golden corpus builder keeps copies of, so that each side's tests hold its own
items to the same lines. One line is one case of three tab-separated cells: the item, its
arguments as a JSON list and its answer as a JSON value, in which a tuple is a list, a set a
sorted list and an absent answer `null`. An item is a function by its name, a constant by its
name, a pattern's text as `<NAME>.pattern` and its flags as `<NAME>.flags`.
`notation_inserted_bases` takes its pattern by the pattern's name; `lone_insertion_gaps` takes
four arguments where its last two are 0; `SpliceWindows.holds` takes a transcript's stable id, its
strand, its introns in the transcript's order as `[first, last]` pairs and the first and last of a
run of positions of chromosome 1, and answers the positions of the run that lie in two windows of
the transcript.

## Threads

The two readings of a default, tab, VCF, JSON or Parquet comparison, the keys of the input
beside a VCF or JSON comparison, the row groups of a Parquet output that is compared whole and
the transcript shards of a candidate JSON cache are read across the host's available cores. `VEP_FORMATS_COMPARE_THREADS=<n>` bounds the
threads at work; it is an environment variable and not an option, and a value that is not a
whole number of 1 or more is refused with exit code 2. A gzip output is inflated by one thread,
and each side of a Parquet comparison is written by one thread, while the workers scan the
blocks of both sides. The comparison of entries and the report stay on
one thread, and every thread count writes the same bytes to the report and to the standard
output.

## Layout

- `src/main.rs`, `src/handover.rs`, `src/cli.rs`, `src/argparse.rs`: the binary, how a run ends,
  the three modes and the command line; `data/usage.txt` and `data/help.txt` are its usage and
  help texts.
- `src/classes.rs`, `data/classes.tsv`, `src/consts.rs`: the class names, the class tables and the
  constants several modules read.
- `src/input.rs`, `src/parse/`: reading a file as text and the input records; the parsers of the
  default and tab formats (`text.rs`), VCF (`vcf.rs`), JSON (`json.rs`), Parquet (`parquet.rs`)
  and the recoder (`recoder.rs`), and the reference side's text classes (`mod.rs`).
- `src/stream/`, `src/remainder/`: the record groups of two outputs read as streams (blocks cut
  at line ends, the gzip members, the UTF-8 check, the digests, the two order conditions, the
  workers; `columns.rs` the default and tab formats, `vcf.rs` the VCF format and the input's
  keys, `json.rs` the JSON format and its canonical form) and what a comparison keeps back of
  them.
- `src/compare_files.rs`, `src/compare.rs`, `src/entry.rs`, `src/counter.rs`, `src/fields.rs`,
  `src/lists.rs`, `src/floors.rs`, `src/regulatory.rs`, `src/record_rows.rs`: the comparison-mode
  report of one format, the comparison of two entry lists, and what it counts and reads of an
  entry.
- `src/manifest.rs`: the corpus manifest and the divergences it documents.
- `src/genetic_codes.rs`: the genetic codes a transcript is translated under, read from
  `crates/vep-core/tests/data/genetic_codes.tsv` (the 64 codons under every table of
  `Bio::Tools::CodonTable`) with genetic code 26 beside them, each engine's reading of a code
  (Ensembl VEP reads code 26 under its module's last table, vep-rs as the code itself), and the
  codes that read every codon as the standard code does, found by comparing the tables. A row that
  asks for a transcript whose `codon_table` names no held code ends the run with
  `ERROR: [compare_vep_formats] transcript <id>: ...` and exit code 2. `tools/vep-golden-corpus`
  compiles this file as well.
- `src/rules/`: the built-in classes, one module a family (`mod.rs` holds what a rule reads and
  the dispatch by class name), with the readers they share: a record's alleles, a transcript
  model, transcript edges and splice windows.
- `src/cache.rs`, `src/nearest.rs`, `src/fasta.rs`: the transcript shards of a candidate JSON
  cache, its nearest transcription start sites, and a reference sequence read through its index.
- `src/shared.rs`, `src/shared/`: the row, span and stop-codon readings the comparison shares
  with the tuple comparison.
- `src/self_parity.rs`, `src/warnings.rs`: the two other modes.
- `src/field_class.rs`: the `field-class` command; `tests/field_class.rs` and
  `tests/data/field_class/` hold it to planted questions and answers, and
  `tests/data/builder_parity/` and `tests/data/builder_cases/` hold the items the golden corpus
  builder keeps copies of.
- `src/report.rs`, `src/error.rs`: the report and its writer; the comparator's error.
- `src/pyfmt.rs`, `src/pyfmt/`, `src/pyjson.rs`, `src/pattern.rs`: the text forms, number
  readings, Unicode character classes, JSON reading and writing, and pattern matches the reports
  and error lines are written in.
- `tests/`: one file a family of cases, each naming what it shows; `tests/cli_suite.rs`,
  `tests/cli_captures.rs` and `tests/cli_parquet.rs` run the binary on planted files and hold its
  report, standard output, standard error and exit code; `tests/compare_both_paths.rs` holds a
  streamed comparison to the comparison of the whole outputs, `tests/parquet_stream.rs` a
  streamed Parquet comparison, and `tests/parquet_memory.rs` holds the bytes one holds at once
  under the bytes of the reference's rows; `tests/data/` holds the planted
  cases and `tests/parquet_writer/` writes the Parquet files the tests read.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-formats-compare/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

CI runs the test, clippy and fmt commands above in its Test, Clippy and Format jobs, builds the
package on the declared `rust-version` in MSRV and audits its lockfile. Its Harness Tests job
builds the release binary and runs `scripts/concordance/test_feature_oracle.sh` through it.

## Harness

`scripts/concordance/feature_oracle.sh` compares every output format of a case, and the two
default runs' warnings files, with this binary, and reads back its report. One environment
variable, `VEP_FORMATS_COMPARATOR`, names the build:

- Unset or empty: this tree's release build,
  `tools/vep-formats-compare/target/release/vep-formats-compare`. The runner builds nothing.
- Any other value: a build of the tool, as a path or a command on `PATH`.

A comparator that is not an executable stops the runner, a dry run as well, before its first
engine run with one line that names the build command and the variable. The runner reads a
comparison as failed when the binary exits 2 or leaves no report, and puts the end of its
standard error into the case's record; exit codes 0 and 1 are both a result, judged from the
report. `scripts/concordance/test_feature_oracle.sh` runs through the build
`VEP_FORMATS_COMPARATOR` names, else a release build under `CARGO_TARGET_DIR`, else the tree's
own.
