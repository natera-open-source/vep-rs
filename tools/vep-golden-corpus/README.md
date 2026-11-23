# vep-golden-corpus

The builder of the golden corpora under `tests/golden/`. A corpus is a small VCF of exemplar
records, one or more for every consequence set an Ensembl VEP reference output holds, beside a
`manifest.json` that records how the records were selected, the runs the golden tests make over
them (the manifest's own and its `flag_sets`) and, for each run, the documented divergences between
Ensembl VEP's output and vep-rs's. The tool selects a corpus, classifies a run's divergences, checks
a corpus against its references, replaces the identifiers under restricted redistribution terms in
stored files, keeps the manifest's lists and records the files a corpus does not commit; every byte
it writes is fixed.

```text
vep-golden-corpus select --release R --assembly A --dataset SID=REFERENCE=INPUT --out DIR [options]
vep-golden-corpus classify --corpus DIR --vep-default FILE --vep-rs-output FILE [options]
vep-golden-corpus check --corpus DIR --dataset SID=REFERENCE=INPUT [options]
vep-golden-corpus identifiers PATH [PATH ...] [--check]
vep-golden-corpus covers DIR [--list] [--add ROW] [--flag-set NAME] [--matrix TSV]
vep-golden-corpus flag-set DIR --name NAME [--flag=FLAG] [options]
vep-golden-corpus flag-set DIR --list
vep-golden-corpus generated --corpus DIR --command WORD [WORD ...]
vep-golden-corpus generated --corpus DIR --check
```

`vep-golden-corpus --help` describes the seven subcommands and the manifest keys each reads and
writes; `vep-golden-corpus <subcommand> --help` lists a subcommand's options with their defaults.

## Contract

Wherever a byte could be written in more than one way it follows a CPython rule, named here by its
CPython name; those rules are part of the contract.

**Command line.** The arguments are read as CPython's `argparse` reads them. A long option is named
whole or by any prefix that fits it alone (`--li` for `--list`) and takes its value from the next
argument or after `=`; a repeatable option keeps every value in order and any other option its last
one. An argument that starts with `-` is a value only when it is `-` alone, starts like a negative
number (`-` or `-.` and a digit) or holds a space; after `--` every argument is a value. The corpus
of `covers` and `flag-set` stands anywhere among their options. `generated` reads no abbreviated
option, and where `generated` is the first argument every argument after the first `--command` is
a word of the command, a `--` and an option of the subcommand among them, so `--command` is the
last option of a run. Arguments are acted on from the
left, so the first refusal or `-h` ends the reading, a required option that is missing is named
once every argument is read, and a subcommand's refusal comes before the arguments nothing took. A
refused command line writes the usage of the parser that refuses it and `<prog>: error: <message>`
to standard error and ends with exit code 2. `-h` and `--help` write the help to standard output
and end with exit code 0. The help texts under `src/help/` are the ones `argparse` prints at 80
columns, the width it takes where it finds none; the width of a terminal changes no byte. A count
(`--k`, `--seed` and the other integer options of `select`) takes what `int()` takes within 64
bits: blanks around the number, a sign, single underscores between ASCII digits.

**The manifest.** `<corpus>/manifest.json` is read as
`json.loads(path.read_text(encoding="utf-8"))` reads it: `\r\n` and a lone `\r` are newlines, any
JSON value stands between blanks, `NaN`, `Infinity` and `-Infinity` are numbers, a control
character inside a string is refused, and a key an object repeats keeps its first place and its
last value. Every object keeps its keys in the order written. An integer keeps its digits at any
length; a number with a fraction or an exponent is a double. A manifest is written whole as
`json.dumps(manifest, indent=1)` and a newline: one space an indent level, `[]` and `{}` for an
empty list and object, every character outside ASCII as a `\u` escape (`ensure_ascii`), a double as
`float.__repr__` writes it. That text passes through the replacement of **identifiers** before it
is written, so that no key a subcommand fills from an input or an output carries an identifier
under restricted redistribution terms. A manifest a subcommand does not change is therefore written
back to its bytes when it was written that way, which holds for every committed corpus that has a
`variants.vcf`.

Before a subcommand reads on, the manifest is held to its shapes, and one that is refused ends the
run on `<path> <reason>`: the file is a JSON object; `covers` and `reference_warnings` are lists of
non-blank strings; `masks` is a list of mask classes (`regulatory_feature_rows`, `order_only`, the
four classes the golden harness alone gives a meaning and the corrected-behaviour classes of
`src/constants.rs`); `perl_mode` is `offline`, `database`, `gtf` or `gff`; `plugin_files` maps
plugin names to lists of paths inside the corpus without whitespace; `generated` holds a `command`,
a non-empty list of non-empty words, and a `sha256`, a non-empty map of relative paths without a
`..` part to 64 lower-case hexadecimal digits; `keep_samples` is a count and
`samples` names that many; `required_genotype_classes` lists genotype classes; `record_class_key`
is `regulatory_class`; `flag_sets` is a list of entries under distinct names, each with the shapes
`flag-set` writes; `species` is an Ensembl species name; `cache_layout` is `ensembl` beside a
`release` and an `assembly`; `included_records` maps dataset names to non-empty lists of record
names; `fasta_in_cache` is a boolean that is true only beside a `fasta`; and a
`vep_rs_consequence_sets_with_reference_sequence` key stands on a `divergences` entry alone, as a
sorted list of distinct, sorted consequence sets that neither repeats the entry's
`vep_rs_consequence_sets` nor names Ensembl VEP's set beside another.

**Reference outputs and inputs.** A reference output and an input VCF are read as CPython reads
`open(path, "rt", encoding="utf-8", errors="replace")`, through `gzip.open` for a name that ends
in `.gz` or `.bgz`: a line ends at `\n`, `\r` or `\r\n` (universal newlines), a last line needs no
line end, each run of bytes outside UTF-8 reads as U+FFFD, and a byte order mark stays in the first
line. A gzip file is every member in turn, the zero bytes between members skipped; an empty file
has no line. A row of a default-format output is a line of seven tab-separated columns or more
that does not start with `#`; a column the line lacks reads as `-`, columns past the fourteenth
are not read, and the fields of the Extra column are its `;`-separated `KEY=VALUE` items, a key
given twice holding its last value. The rows of one input record are the consecutive rows of one
Uploaded_variation and Location. In an input VCF a line that starts with `##` is skipped, another
line that starts with `#` is a header line, and a record is a line of eight columns or more. A
file that cannot be opened or read ends the run on `<path> is not a readable text file (<reason>)`
with the text of the CPython exception as the reason: `[Errno <n>] <strerror>: '<path>'`,
`Not a gzipped file (<repr() of the two bytes where a member starts>)`, `Unknown compression method`,
`Compressed file ended before the end-of-stream marker was reached` for a gzip file cut short,
`Error -3 while decompressing data: <zlib's reason>` for deflate data that does not inflate,
`CRC check failed <stored> != <computed>` in `hex()` digits and `Incorrect length of data produced`
for a member that fails its trailer. The whole lines before the place a file fails are read first.

**Records, names and classes.** A record of an input is found in a reference output under the
first `;`-separated element of its ID. Where that element is `.`, empty or `0` it is found under
the names Ensembl VEP writes for a record without one: `CHR_POS_REF/ALTS` from the raw line
(release 116), the same alleles at the start the parser leaves (release 115: POS+1 for an indel
whose alleles share their first base, moved further over the leading bases a bi-allelic indel
shares), `CHR_POS_ALTS` and `CHR_POS+1_ALTS` for a symbolic or breakend allele, and for an empty or
`0` element the name built from the alleles as the parser leaves them, in the case written. POS is
read as `int()` reads it. The variant class of a row is read from its Location and Allele: `snv`,
`deletion`, `insertion` and `mnv_or_complex`; `breakend`, `breakend_insertion` and
`single_breakend`; `symbolic_<allele>` for any other allele that is not made of `ACGTNacgtn*-`. A
span of two adjacent positions is `mnv_or_complex` when the input record of the row's name holds
the allele as a substitution at those positions, and `insertion` otherwise. The classes of
`--family-columns`, of the genotypes of kept samples and of `--record-classes` are the ones
`vep-golden-corpus --help` defines; a genotype is read as Ensembl VEP reads one for
`--individual`, its bits ASCII digits that index the record's alleles, and a record Ensembl VEP
reads as a structural variant shows no genotype class.

**select.** `select` reads the datasets `--dataset SID=REFERENCE=INPUT` in the order given, each
input for its substitutions and then each reference for its rows, and keeps three things of the
rows: for every consequence set a uniform sample of at most 64 rows and the count of its rows; with
`--family-columns`, for every family stratum (the family classes and the variant class) a sample of
at most 8; and every row of the records `--include SID=NAME` names and, under `--keep-records`, of
the records and `included_records` of the manifest being regenerated. A sample takes each row
until it is full; after that the row of count `n` draws `randrange(n)` and takes the place drawn
when the place is inside the sample. A named record without a row is one line and exit code 1.

A regeneration writes a manifest without the keys `classify` wrote into the one it replaces, and
no `classify` writes a `vep_rs_consequence_sets_with_reference_sequence`, a key written by hand. A
`divergences` entry of that manifest, its own or a flag set's, that carries the key therefore ends
the run before a reference is read, with nothing written, on `select would discard
vep_rs_consequence_sets_with_reference_sequence, a key written by hand, on <whose> divergences
entry on <Location> <Allele> <Feature>[; ...]: pass --drop-reference-statements to discard it`,
where `<whose>` is `the manifest's` or `flag set <name>`. Under `--drop-reference-statements` the
run goes on and, as it writes the manifest, writes one line an entry to standard error:
`WARNING: [build_golden_corpus] select discards <key> on <whose> divergences entry on <Location>
<Allele> <Feature>`. A run that is refused for another reason, or whose manifest cannot be
written, names no entry.

From the samples `select` chooses the exemplars. The consequence sets are taken in sorted order;
the sample of each is shuffled, then sorted by whether the row is wider than `--max-span`, by its
cache cost under `--cache`, and by its preference (the sum of 4 for a row with FLAGS, 3 for a
position that holds `?` and 2 for a position that is a range), highest first. Up to `--k` rows are taken, or
`--focus-k` for a set that holds a `--focus-term`: first the rows of distinct (variant class,
strand, family classes), then the others, and every row for a count of zero or less. Every variant
class no chosen row shows then takes the first row in that order of the whole sample, every family
stratum the first of its own sample, and the named rows are added. With `--record-classes` the
references are read once more, a record at a time, and each class short of
`--record-class-exemplars` records takes records of its own sample of 16 in the same order.

Each input is then scanned once. A record whose first name among those it can carry is the
Uploaded_variation of a chosen row on the row's chromosome joins the corpus, as its first eight
columns, or with `--keep-samples N` those, `FORMAT` and the first N samples. The line is kept with
every identifier under restricted redistribution terms replaced by the fabricated one of its shape,
read as **identifiers** reads one text. A record the replacement would rename, its ID being such an
identifier, is passed over by every sample, and when a reference row names one the run ends, once
its input is read, on `<input> names <n> record[s] the references chose by an identifier under
restricted redistribution terms in the ID column; a corpus carries no such record, since its
reference rows are named by the identifier`, with nothing written. `--multi-allelic-per-dataset`
sequence records a dataset and multi-allelic stratum are sampled beside them, and with kept samples
`--genotype-exemplars` records a genotype class the samples show; a class takes records of its
sample until it has that many. More records than `--max-records` are cut to the first record of
every consequence set, family stratum, genotype class and record class, the named records, and then
the others in sorted order up to the bound. `variants.vcf` is `##fileformat=VCFv4.2`, the first
header line of the inputs cut to the kept columns, and the records ordered by chromosome (a name of
digits by its number, then `X`, `Y`, `MT`, then the others as they came), POS, REF and ALT. The
manifest holds, in a fixed order, the options of the run, the keys a regeneration carries from the
manifest it replaces (`covers`, `masks`, `perl_mode`, `reference_warnings`, `plugin_files`,
`generated`, the `flag_sets` without the seven keys `classify` writes), the datasets with the file
names, row
counts and set counts of their references, every consequence set with its rows and exemplar
records, the counts of each class among the chosen rows and records, the records, and what stayed
uncovered or unresolved. One line of standard output counts the same, and the run ends with exit
code 1 when a set, a genotype class or a record class stays uncovered or an exemplar has no record.

**The draws.** Every draw is CPython's `random.Random(seed)` for the integer `--seed`: the Mersenne
Twister MT19937 seeded by `init_by_array` from the 32-bit digits of `abs(seed)`, least significant
first. `randrange(n)` draws as many bits as `n` has (`getrandbits`: the top bits of one 32-bit
output, or a whole output as the low half and the top bits of a second above it) until they spell a
number below `n`; `shuffle` walks from the last place down and changes each item with one drawn
among the places up to its own. Four generators start from the seed: the samples of the rows, the
shuffles of the candidates, the samples of `--record-classes`, and the samples of the input scan.
`tests/fixtures/random/draws.tsv` holds CPython's draws for 22 seeds: outputs, `randrange` over 29
ranges, shuffles, samples and mixed runs.

**The cache cost.** Under `--cache DIR` a row costs the transcripts of that cache it pulls into
the corpus cache: every record of the shards `DIR/transcripts/<chromosome>/<first>-<last>.json[.gz]`
whose `start` to `end` overlaps the row's Location widened by `--flank`, a record counted once in
each shard that holds it. A plain file under `transcripts/` names no sequence and is passed over. A
bracket breakend allele adds its mate: on the row's chromosome the
whole span to the mate is costed in place of the row's own, on another chromosome the mate's
position is added (`chr` before its name dropped). A record costs
`len(zlib.compress(json.dumps(record, separators=(",", ":"), ensure_ascii=False).encode("utf-8")))`:
its JSON text as CPython writes it without blanks (keys in the file's order, a repeated key at its
first place with its last value, a double as `float.__repr__` writes it, characters outside ASCII
as they are), in the zlib stream of zlib's deflate at level 6 with a 32 KiB window and memory level
8, header and checksum counted. The deflate library is therefore part of the contract: the tool
builds stock zlib from the source its `libz-sys` crate carries (zlib 1.4 under the lockfile) and
links it statically, so the host's library is not read.
`tests/fixtures/cache_cost/lengths.tsv` holds, for each of 113 caches of the committed corpora, the
count and the sum of the costs of its records and a hash of every cost in order, 88,528 records in
all, as CPython gives them on zlib 1.3.2; zlib 1.2.11 and zlib 1.4 give the same lengths. The test
that reads the table checks every cache it names; a corpus whose cache the table does not name owes
no row.

**check.** `check --corpus DIR` first reads every file under the corpus as `identifiers --check`
reads a directory, and refuses a corpus that holds an identifier under restricted redistribution
terms in another form than the fabricated one: one line a file on standard error,
`<file below the corpus> holds <n> identifier[s] under restricted redistribution terms in another
form than the fabricated one; the identifiers subcommand replaces them`, exit code 1, and no
reference read. Otherwise it reads the references again and prints
`observed N, covered by corpus M, missing K` for the consequence sets, then one `  MISSING <set>`
line for each of the first 50 missing in sorted order. With `--family-columns` it prints
`family strata observed N, covered by corpus M, missing K` and up to 50
`  MISSING family stratum <stratum>` lines, a stratum covered when a row that shows it is of a
corpus record of its dataset and chromosome; for a corpus with kept samples
`genotype classes observed N, covered by corpus M, missing K` and a `  MISSING genotype class`
line for each class the inputs or `--require-genotype-class` name and no record shows; for a
manifest with a `record_class_key` the same for `record classes`. Its error lines go to standard
error after those: a `covers` name that is no row of `--matrix`, a `variants.vcf` header that does
not name `FORMAT` and the manifest's samples, a record of another column count, a file without a
header line. The run ends with exit code 1 when anything is missing or any error was written; a
refusal that ends it earlier comes after the error lines found before it.

**identifiers.** `identifiers PATH...` replaces every identifier under restricted redistribution
terms by the fabricated one of the same shape, its letters and then a zero for each of its digits,
in each named file and in every file under each named directory: the step a stored Ensembl VEP
output or an annotated input passes through before it is committed. Two shapes are read. The first
is `COS`, a capital letter and digits, wherever it stands as a word: with no letter, digit or `_`
of ASCII on either side, the two characters `\t` or `\n` of a JSON string counting as the end of
the word before them. The second is a letter pair and a count of digits that the variation files of
Ensembl VEP's release 116 human caches give a name of that source, on GRCh37 and GRCh38: `C` and
one of `D`, `G`, `I`, `M`, `R`, `S` and `X` with six or seven digits, `CP` with five or six, `BM`
with seven, `H` and one of `D`, `I`, `M`, `R` and `S` with four or six, and `HX` with six. A public
sequence accession can share that shape, and so it is read in four positions alone: as a member of
a list, after `&`, `,`, `|` or `;` and before one of those, a tab, a carriage return, `"`, `\` or
the end of the line; as the member that opens a line or a field, after `=`, a tab, `"`, `\t` or
`\n` and before `&`, `,` or `|`; as a JSON `id`, `"id":"<identifier>"` with or without a space
after the colon; and as the whole `Existing_variation` cell of a row of an output table, a column
that names known variants and never a sequence. A table is the rows below a header line, one `#`
and tab-separated names, that names the column: a header line that does not name it ends the table,
and a line that opens with `##` is no header. A second-shape identifier read in one of those
positions in one file of a directory is replaced wherever it stands as a word in every file of that
directory, though not after `:` and not before `.` and a digit, so that the outputs of one record
agree; one that stands alone in another table cell in every file is not read. Two identifiers of
one prefix and length become the same one, and a text that holds fabricated ones alone stays as it
is.

The files of a directory are every file under it, a link to a file among them and a link to a
directory not followed, in the order CPython sorts paths: part by part, so that `a/b` comes before
`a.txt`. A file is read as text when its bytes, inflated as `gzip.decompress` inflates them where
they open as a gzip member does whatever the file's name, hold the first shape's letters before a
digit or the second shape with no letter, digit, `_` or version (`.` and a digit) after it, and
hold no zero byte and UTF-8 alone; any other file holds nothing to replace. `--check` writes no
file: it writes `ERROR: [build_golden_corpus] <file> holds <n> identifier[s] under restricted
redistribution terms in another form than the fabricated one; the identifiers subcommand replaces
them` for each file that holds one, then `<path>: <n> files read, <m> hold an identifier in another
form than the fabricated one` on standard output for each named path, and ends with exit code 1
when a file held one. Without `--check` every named path is read before a file is written, so a
refusal leaves each tree as it was. A plain file is written back as its text. A gzip file is
written back as one member, the bytes of `zlib.compress(data, level=9, wbits=31)`: zlib's deflate at
level 9 under the header zlib writes, with no name and no time stamp. A block-compressed file (a
gzip file whose first member carries the extra subfield `BC`, as `bgzip` writes) that holds an
identifier ends the run on `<file> holds <n> identifier[s] under restricted redistribution terms and
is block-compressed: decompress it, pass the text through this subcommand, then compress and index
it again; nothing was written`, since the index beside it describes the bytes as they are stored.
Under a named directory the nearest `provenance.json`, in a rewritten file's own directory or above
it and no higher than the named directory, has every SHA-256 digest of that file recorded again:
the digest of its bytes and the digest of its uncompressed bytes are replaced by the ones after the
rewrite wherever their 64 hexadecimal digits stand in the record, which is read as text with
universal newlines. Each named path ends on one line of standard output,
`<path>: <n> identifiers replaced in <m> files, <k> digests recorded again`. A path that is neither
a file nor a directory ends the run on `<path> is neither a file nor a directory`, and a file that
cannot be read on `<file> is not a readable file (<reason>)`, the reason the text of the CPython
exception: `[Errno <n>] <strerror>: '<path>'`, or for a gzip file `Not a gzipped file (<repr() of
the two bytes where a member starts>)`, `Unknown compression method`, `Compressed file ended before
the end-of-stream marker was reached`, `Error -3 while decompressing data: <zlib's reason>`,
`CRC check failed` and `Incorrect length of data produced`.

**Threads.** The rows of a reference are read on the host's cores, a block of whole lines a
thread, and taken in file order. The blocks of a BGZF file, the gzip form `bgzip` writes, are
inflated on as many threads and read in the file's order; a member that is no such block, and
every member behind it, is read as any gzip member is. The files `check` and `identifiers` read for
identifiers are read on as many threads, a file a thread, and counted in the order of the walk. No
byte depends on the count of threads.
`VEP_GOLDEN_CORPUS_THREADS` bounds it; a value that is no whole number of 1 or more is refused.
The format comparator `classify` asks reads the shards of the corpus cache on threads of its own,
which that tool's `VEP_FORMATS_COMPARE_THREADS` bounds.

**covers.** `covers DIR` prints the `covers` list of the manifest, or of the flag set
`--flag-set NAME`, one name a line. `--add ROW` (repeatable) adds rows of the parity matrix
`--matrix` (default `tests/parity/matrix.tsv`, read only when a row is added) to the list, each
name once and after the names the list holds, and writes the manifest; the list is printed
afterwards only with `--list`. A name that is not a row is refused with one line a name and exit
code 1, and nothing is written. The rows of a matrix are the first tab-separated column of its
lines after the header, whose first column is `item`; lines end where `str.splitlines()` ends them
and a blank line names no row.

**flag-set.** `flag-set DIR --name NAME` writes one entry of the manifest's `flag_sets`: `name`,
`flags` (every `--flag`), `input` (`--input`, default `variants.vcf`, a file inside the corpus),
`stdin`, `unordered`, `formats` (`--format`, a comma-separated subset of `default`, `tab`, `vcf` and
`json`; all four by default), `expected_dir` (`--expected-dir`, default `expected_<name>`) and
`covers` (`--covers`, rows of `--matrix`), in this order, then `perl_mode` when `--perl-mode` names
a mode other than `offline`. A name is a letter or digit followed by letters, digits, `_` and `-`,
and is not `default`. An entry of that name has those keys rewritten where they sit and keeps every
other key but the seven `classify` writes and `perl_mode`; another name goes last. A rewrite of an
entry whose `divergences` hold a hand-written `vep_rs_consequence_sets_with_reference_sequence` is
refused as `select` refuses a regeneration, the run named `a flag-set rewrite` and the entry
`flag set <name>`'s, and under `--drop-reference-statements` names each such entry on standard
error and rewrites; the manifest's own list and another set's are not the rewrite's and stop nothing. The
run ends on one line of standard output that states what was written. `--list` prints the entries instead, one
tab-separated line each: name, flags joined by spaces, input, `true` or `false` for stdin and for
unordered, formats joined by commas, expected directory, covers joined by commas, a key an entry
lacks read as its default. A `--stdin` run whose flags name no `--format`, an input that is no file
of the corpus, a path that leaves the corpus and a covers name that is no row are refused before
anything is written.

**generated.** A corpus that does not commit a set of files a seeded generator reproduces records
them in its manifest's `generated` block: the `command` that writes them and the `sha256` of each
file, by its path under the directory written. `generated --corpus DIR --command WORD...` runs the
words as one command, without a shell, from the repository root with standard input closed,
`{generated}` in a word replaced by a scratch directory and `{corpus}` by the corpus directory with
its links resolved, as the golden harness runs the block. It then writes the manifest with the
block, where one stood or else last: the words as given and the digest of every file under the
scratch directory in path order, a link to a file counted under its own name and a directory that
is a link not entered. The run ends on `<corpus>: <n> generated files recorded`.
`generated --corpus DIR --check` runs the manifest's own command the same way and writes nothing:
one line of standard error names each path, in path order, that `is recorded and not written`,
`is written and not recorded` or `hashes as <digest>, the manifest records <digest>`, then
`<corpus>: <n> generated files written, <m> differ from the manifest's record` goes to standard
output, and the exit code is 1 when a path differs. The repository root is the tree the executable
stands in: the nearest directory above it that holds `tools/vep-golden-corpus`. A command is
refused by name, its words joined by spaces between backquotes, when it `does not start (<reason>)`
(the text of the `OSError` that names the program, or `embedded null byte` for a word that holds
one), when it `failed with status <status>: <line>` (a signal's number below zero for a run a
signal ended; the last line of its standard error that holds anything, read as UTF-8 with a
replacement character for bytes outside it, or `it wrote nothing on standard error`), when it
`wrote a file that cannot be read (<reason>)` and when it `wrote no file under {generated}`. A
run with neither option, with `--command` and no word or an empty word, with `--check` beside
`--command`, or with `--check` on a manifest without a block is refused before a command runs.

**classify.** `classify --corpus DIR --vep-default FILE --vep-rs-output FILE` compares the
default-format outputs of the two engines for one run, the manifest's own or with `--flag-set NAME`
that set's, and writes the run's lists into the manifest or into the set's entry, after its other
keys and in this order. `divergences` holds one entry for each consequence set Ensembl VEP writes on
a key (Location, Allele, Feature, Feature_type) and vep-rs does not write there: the key, the set,
the sets vep-rs writes on the key and the class that explains the difference. `vep_rs_only_tuples`
holds one entry for each key vep-rs alone writes. `divergence_summary` and `vep_rs_only_summary`
count the entries under each class, the classes in the order they first come. Under `--field NAME`
(repeatable) `fields_compared` names the fields, `field_divergences` holds one entry for each named
field of the Extra column whose two values differ on a key and consequence set both outputs hold,
and `field_divergence_summary` counts them. A named field is read on every row of the key and set,
the rows of the two outputs paired in their order: the entry states the first pair whose two values
differ, each side's own row, which the classes are asked beside, and where other pairs hold equal
values, the rows of other records that share the key, the entry's last key `rows_alike` counts
them, the rows the golden harness compares like any other. A key whose two outputs hold other
numbers of rows on the set, or whose pairs that differ hold other values than each other, is read
on the first row of each output alone and writes no `rows_alike`. The keys come in the
order of their four texts, compared by code point; every entry lists `record_indices`, the ordinals
of the run's records that its rows name; a difference no class explains is an
`unexplained_residual`. A field difference takes the class of the first class that takes it, of
twenty-seven asked in one order (`src/classify/fields.rs`). A field named `Amino_acids` or `Codons`
is the row's column of that name, which the class `alternate_codon_past_the_coding_end` reads; it
reads the residues beside both rows' codons, so a run that names `Amino_acids` names `Codons` too.

**The further formats.** `--vep-tab`, `--vep-vcf` and `--vep-json` name Ensembl VEP's output of
the run in the tab, VCF and JSON formats, and `--vep-rs-tab`, `--vep-rs-vcf` and `--vep-rs-json`
vep-rs's. A format is read where both of its outputs are named; one of the two alone, or a format
without `--field`, is refused. The golden harness compares every field of a documented row in every
format and finds a documented field under the name the format gives it, in any letter case. So on
a key with an entry under `--field` the row is read in each format named, in the order tab, VCF,
JSON, and each field of it whose two values differ and that no entry of the key names, in any
letter case, is one more entry behind the key's, the fields of one format in the byte order of
their names. A key without an entry is read in no further format.

A tab row is found by its key and consequence set; a tab output whose columns leave out
`Feature_type` is read without it, and a `GA4GH_VRS` column that holds JSON is compared as the
object it holds. A VCF entry and a JSON consequence object are
found by the record they stand on, their Allele, their Feature and their consequence set: the
record is the next of the run's records with the first five columns of the line, or of the
object's `input`, so identical lines are matched in their order, and a line or an object of no
record is not read. For a flag set over `variants.vcf` re-encoded in another format, a JSON
object names its record by the input's own line, and the lines of the VCF output, which carry no
input line, are the records in turn; a VCF output of another number of lines is refused. The VCF
entries are those under `CSQ`, or under the key the run's
`--vcf_info_field` names, with the fields its last header line names. A value is written in the
default format's spelling: `-` for an empty VCF value and for a member a JSON object lacks, `,`
for the `&` the VCF format writes in its place. A JSON text is written as it is, a list as its
items joined by `,`, and any other value as the harness prints it: without blanks, the members of
an object in the byte order of their names, a number that is no 64-bit integer in the shortest
digits that read back as it. The class of such a field is asked as that of a compared one; a
value that holds a tab or a line end is no question for the format comparator and takes no class
of its. Under `order_only` the three multiset columns compare as multisets under their names in
any letter case, and so do the two further columns the golden harness compares so: the items of
`NearestExonJB` and the groups of `5UTR_annotation`, each group by its pairs.

A run whose flags hold `--total_length` and that compares the JSON format (its entry names no
`formats`, or names `json` among them) holds further entries, with or without `--field`, after
those of the named fields and of the further formats: one for each position column that both outputs print as one range with a
total, `484-487/2674` or `?-435/2674`. The entry's `field` is the JSON end key (`cdna_end`,
`cds_end` or `protein_end`), its `vep` the range's first position, `-` where that is `?`, its
`vep_rs` the last position with the total, and its class `total_length_json_range_end`
(`src/classify/total_length.rs`). An end key that a further format already states on the key is not
written a second time: its entry takes the class where it holds the column's two values and keeps
its own class where it holds others. `fields_compared` is
written under `--field` alone.

The records of a run are the manifest's `records`, or for a flag set with an `input` of its own that
file's lines that are no header. A line that shares its first five columns with a manifest record is
that record; another line of a file that holds a VCF record (five columns or more, a POS of decimal
digits as `str.isdecimal()` reads them, a REF made of `ACGTNacgtn.`) is a record of its own; a file
that holds no VCF record is the manifest's records in another encoding, line for line. A region line
of such a file (`chr:start-end`, a strand where one is written, then `/` or `:` and the allele)
names its record in both outputs, as Ensembl VEP names a region record by its line, and a refusal
names that record by the line's position, without its allele. A record no
Ensembl VEP row names is skipped: it carries a `reference_warning`, takes `reference_rows` 0 before
it as its last two keys, and the keys vep-rs alone writes for skipped records are
`reference_skipped_record`. A record without a row needs no warning and is not skipped when its ALT
is `.`; when the reference's JSON output (`json.jsonl.gz` in the run's expected directory, else in
`expected/`) holds an object whose `most_severe_consequence` is `?` and whose `input` opens with the
record's first five columns; and when the run's flags carry `--coding_only`, `--no_intergenic`,
`--check_frequency` or `--filter_common` and vep-rs writes no row for it either. The `masks` of a
run, a flag set's own where its entry has the key,
leave the reference's RegulatoryFeature and MotifFeature rows out under `regulatory_feature_rows`
and compare `DOMAINS`, `FREQS` and `ZYG` as multisets of their comma-separated items under
`order_only`. A `vep_rs_consequence_sets_with_reference_sequence` key on an entry of the run's
`divergences` is carried onto the entry written for the same Location, Allele, Feature, records and
Ensembl VEP set, as its last key. A flag set whose entry holds no `divergences` is held to the
manifest's own list by the golden harness, so `classify --flag-set` carries the keys that list
states into the set's first list, or refuses the run as it refuses one of the manifest's own.
`--flag` and `--fasta-name` describe the manifest's own run and set its `flags` and `fasta`.

The class of a `divergences` entry is asked of each set vep-rs writes on the key, in the order of
the sets, and the first of three readings that names one decides; the rules of the first and the
third are the tuple comparator's own (`tools/vep-outputs-compare/src/rules/`). A key whose sets the
three name under no class is asked a fourth, the last of the list.

- The set rules of the tuple comparator, in that tool's order, on the two sets and the Allele. The
  class is the `taxonomy_class` of the first rule that takes them, as that tool's class file
  states it (`data/divergence_classes.tsv`): `splice_family_swap` or `start_cooccurrence_swap`.
  A pair whose allele opens with its whole reference allele is put to the
  row-aware rules also where a set rule classifies it without setting it aside, and is counted
  under the row-aware rule that takes it.
- `breakend_mate_local_read`, on a key whose Location is one position and whose Allele names a mate
  in brackets. The corpus cache (`json_cache/`, or `cache_root/<species>/<release>_<assembly>`
  under `cache_layout` `ensembl`) files the Feature under the mate's chromosome, with the mate
  inside the transcript or within 5,000 positions of it and the Location not so on that
  chromosome, a `chr` that leads the mate's or the Location's chromosome set aside in any case.
  Ensembl VEP's set is what a breakend at the Location's position supports on the transcript
  (`src/mate_context.rs`), `intergenic_variant` where that is nothing, vep-rs's set is what one at
  the mate's position supports, and the two differ.
- The row-aware rules of the tuple comparator, in that tool's order: the five classes of a coding
  end, the five stop codon classes, `stop_codon_settled_before_unknown_bases`,
  `frameshift_inside_a_held_stop_codon`, `splice_term_past_the_record_span`,
  `start_codon_removed_without_five_prime_utr` and `start_codon_kept_without_five_prime_utr`.
  Every row
  of Ensembl
  VEP's output that carries the key and its set is put beside every row of vep-rs's that carries
  the key and its set, each row read by its Location, its Allele, its set, its `cDNA_position`
  and its `CDS_position`, by its `Protein_position`, `Amino_acids` and `Codons` where the row has
  all twelve columns, and Ensembl VEP's by its first column as well. A pair of rows takes the class of the first rule that takes
  it, and the entry the class every pair takes.
  `stop_codon_settled_before_unknown_bases` holds a pair to its `Amino_acids` and its `Codons`:
  the `Codons` cell is the same on both rows and settles a stop codon before its unknown bases,
  Ensembl VEP's row writes no residues, and vep-rs's residues hold their first stop at that
  codon's place.
  `frameshift_inside_a_held_stop_codon` holds a pair to the two cells as well: both rows hold the
  same residues, a stop first and the unknown residue after the `/`, and the same `Codons` cell,
  whose halves differ by a number of bases that is no multiple of three. A row
  that ends before its `Codons` column takes neither class anywhere.
  `start_codon_removed_without_five_prime_utr` holds a pair to its positions and its residues:
  the `cDNA_position` and the `CDS_position` each open at 1 and end on one number on both rows,
  both rows hold `inframe_deletion`, neither row's `Amino_acids` writes an alternate residue, and
  vep-rs's set is Ensembl VEP's with `start_lost` beside it. A row that ends before its
  `Amino_acids` column takes the class nowhere.
  `start_codon_kept_without_five_prime_utr` holds a pair to its two positions: the
  `cDNA_position` and the `CDS_position` are one span of one number or two on both rows whose
  first number is 1, 2 or 3, both rows hold `frameshift_variant`, `inframe_insertion` or
  `inframe_deletion`, and vep-rs's set is Ensembl VEP's with `start_retained_variant` added and
  without its `start_lost`, where it holds one. A row that ends at its `CDS_position` column
  takes the class as a whole row does. Where the first column spells no alleles, both
  rows are read with the reference bases the Allele stands in for at the Location, as the records
  of the run give them by position. Where two or more records read the row differently, the
  record whose identifier is the first column decides; the rows are read without the bases where
  no record gives them, where no record carries that identifier and where more than one does and
  they read the row differently. A row whose Feature_type is `-`, the one row of a record under `--summary` or
  `--most_severe`, is read by its terms alone: in a flag set's run it takes
  `stop_retained_without_stop_codon_test` or `stop_codon_removed_without_stop_codon_test` where
  the manifest's own `divergences` name its Location and Allele under one of the two, and
  `stop_codon_overlap_insertion_range` where they name its Location under that class, when its two
  sets differ as the class's do; it takes `start_codon_removed_without_five_prime_utr` and
  `start_codon_kept_without_five_prime_utr` nowhere,
  since it carries no position. A row of any other Feature_type than `Transcript` takes none.
  Ten of these rules ask the coding ends of the key's transcript, as the tuple comparator asks
  them (`tools/vep-outputs-compare/README.md`, Coding ends), and `classify` hands them the
  corpus's own cache and reference sequence: the record of the key's Feature in the corpus cache
  whose span shares a base with the key's Location, the transcripts of a sequence read when a key
  on it first asks, and the reference file the run's `fasta` names (`--fasta-name` names it for
  the manifest's own run), read beside the manifest or, under `fasta_in_cache`, in the cache
  directory where none stands beside the manifest; a corpus without the file hands none. Beside
  the cache a class
  of a coding end takes a key where vep-rs's coding terms are the ones derived for the row, a
  stop codon class takes one only on a transcript whose record shows the stop codon at the
  allele's place, and `start_cooccurrence_swap` among the set rules only on one whose record
  shows a start codon. A corpus whose cache holds no `transcripts` directory hands no cache: no
  class of a coding end takes a key there, and the others read the two rows alone. A corpus
  without the reference file leaves a deletion that runs past the last base of a transcript
  whose stop codon polyadenylation completes an `unexplained_residual`. A row whose
  Feature_type is `-` names no transcript and asks no record, so no class of a coding end takes
  the record-level row of `--summary` or `--most_severe`.
  That tool's fourteenth such rule, `pachysolen_nuclear_code_read_as_standard`, is the rule of no
  class a corpus documents rows under (`CORRECTED_BEHAVIOUR_CLASSES` of `src/constants.rs`); it
  is asked of no pair here, so it takes no entry.
- `refseq_edited_model_intron_span_coding_reference`, asked of each set vep-rs writes on the key,
  in the order of the sets, once the three readings above name none of them
  (`src/classify/intron_span.rs`). The key is a Transcript's, its Allele is a run of bases or `-`,
  and it is read from the first row of each output that carries the key and its set. Both sets hold
  `splice_acceptor_variant`, `splice_donor_variant` and `intron_variant`; Ensembl VEP's set less
  vep-rs's is `stop_lost` and `protein_altering_variant`, and vep-rs's less Ensembl VEP's is
  `inframe_deletion`. Both rows carry `BAM_EDIT` `OK`, the same `CDS_position` range with a number
  at both ends and the same `USED_REF`, which spells more bases than the range holds.

A key several records share can hold every set Ensembl VEP writes on it while the row of one
record writes another: every `<CNV:TR>` record of one position writes the Allele `tandem_repeat`
at one Location. Such a key takes one entry for each set Ensembl VEP writes on it that a record
writes beside a set vep-rs alone writes there, a row named by its `Uploaded_variation`. The
entry's `record_indices` are the records that write a row of both sets, and its
`vep_rs_consequence_sets` every set vep-rs writes on the key, Ensembl VEP's own among them. Its
class is the one every such pair of sets takes: `cnv_tr_symbolic_loss_named_a_gain` on a
Transcript's key whose Allele is `tandem_repeat`, where Ensembl VEP's set less vep-rs's is
`feature_elongation` and vep-rs's less Ensembl VEP's is `feature_truncation` or, in a set that
holds no `feature_truncation`, no term (`src/classify/tandem_repeat_loss.rs`).

The run ends on one line of standard output, `divergent VEP rows N (SUMMARY); vep-rs-only tuples N
(SUMMARY); reference-skipped records [ORDINALS]` and, under `--field` or where the run holds
entries of the JSON end keys, `; field divergences N (SUMMARY)`, each summary as `str()` writes a
`dict`. A run that ends with exit code 1 writes
nothing: a record without a row or a warning and a warned record with rows, one line each; a
with-reference key that no written entry carries, or that repeats the sets of the entry it is
carried onto; `--flag` or `--fasta-name` beside `--flag-set`; an input that is neither a VCF nor the
manifest's records line for line, and a line of a VCF input that is neither a twin nor a VCF record,
named by its ordinal and the `repr()` of its first 40 characters; an output, an input or a JSON
output that cannot be read. A line of the JSON output that is no JSON object is named by its number,
and the reason a line is no JSON is the text of the `json.JSONDecodeError` for the line with its
line end.

**The HGVS field classes.** Six classes take a difference in `HGVSc`, `HGVS_OFFSET` or a field of
the HGVSIntronOffset plugin (`HGVS_IntronStartOffset`, `HGVS_IntronEndOffset`), and one of them a
`REFSEQ_OFFSET`. Each reads the differing field beside the other compared fields of the two rows,
their `cDNA_position` and `BAM_EDIT` and the reference row's `STRAND` and `DISTANCE`, and its
`REFSEQ_OFFSET` where the run compares that field; a field the run does not compare reads as
empty on both rows. A notation is of
the row's transcript when the name before its `:` is the row's Feature, each without its version,
and its body is read whole by one of the two patterns of `src/classify/hgvs/notation.rs`, as
CPython's `re.fullmatch` reads a pattern.
`refseq_edited_model_hgvsc_offset` is an `HGVSc` on a transcript whose record in the corpus cache
carries an `_rna_edit` attribute that inserts or deletes bases: two notations of one name that
differ in their cDNA coordinates alone, some integer of the one an edit distance of the transcript
from some integer of the other. The edit distances are the running sums of the edits' length
changes in start order, as absolute values other than 0; an attribute's value is `<start> <end>
<bases>` split on single spaces with its trailing empty fields dropped, where three fields whose
bases count `end - start + 1` substitute and state no edit, three fields of another count insert
their bases, and any other count of fields deletes `end - start + 1` bases.
`refseq_edited_model_hgvs_offset` is an `HGVS_OFFSET` on such a transcript that Ensembl VEP leaves
empty and vep-rs fills with a non-zero integer. `hgvs_added_bases_written_outside_the_transcript`
is a field of vep-rs's `dup` or `ins` notation on a row whose one term is `upstream_gene_variant`
or `downstream_gene_variant`, where Ensembl VEP prints no `HGVSc` and no `HGVS_OFFSET`: the
notation is of as many bases as the allele holds or adds to its span, beside no offset on a
downstream row and beside a non-zero offset on an upstream row, negative on `STRAND` -1 and held
to the reference row's `DISTANCE`. `hgvs_shifted_indel_written_span_gate` is a field Ensembl VEP
leaves empty, on a row with no upstream or downstream term, where vep-rs writes a `del`, `dup` or
`ins` of the row's transcript beside a non-zero `HGVS_OFFSET`; on such a row whose `BAM_EDIT` is
`OK` on both sides it is also a `REFSEQ_OFFSET` Ensembl VEP leaves empty where vep-rs writes a
non-zero integer.
`hgvs_shifted_deletion_clamped_span` is an `HGVSc` on a deletion's row to which Ensembl VEP gives
a non-zero `HGVS_OFFSET`: two `del` notations of the row's transcript that end at one coordinate,
vep-rs's starting earlier. `hgvs_clipped_reference_against_whole_alternate` is an `HGVSc` on a row
of a sequence allele whose reference `cDNA_position` is open on a side (`?`): two notations of the
row's transcript and one numbering whose coordinates lie within the span's length plus one of each
other, a reference coordinate whose intron distance is the `REFSEQ_OFFSET` of a row with `BAM_EDIT`
`OK` and a position open at its end standing for every base from its exon base to itself and a
reference substitution that names two or more reference bases on a row with `BAM_EDIT` `OK` read
as one of one base, a vep-rs `del` against no reference notation, or a reference `ins` or `dup`
(on the row of an allele shorter than its span also a `del`, a `delins` or a substitution) against
no vep-rs notation; and a field of the plugin beside such a notation, where a side's two values
are the intron parts of its own notation's coordinates. The doc comments of `src/classify/hgvs/`
state every condition of each class.

`shifted_allele_rotation_by_strand` takes an `HGVSp` on a transcript row whose reference `STRAND`
is -1 and whose Allele is two or more of the bases `ACGT`: two values that name one protein before
their last `:p.`, where both rows carry the same `cDNA_position`, `CDS_position` and `HGVSc` and
one non-zero `HGVS_OFFSET` larger in size than the Allele's length, and the Allele as written
differs from the Allele with its last base moved to its start once per base of the offset.
Ensembl VEP turns a reverse-strand allele by its length less the shift, so past that length it
reads the protein notation from the allele as written (`TranscriptVariationAllele.pm` lines
1337-1339). A run names the class only where it compares `HGVSc` and `HGVS_OFFSET` beside
`HGVSp`, and asks it before the organelle class.

`hgvsp_shifted_indel_written_span_start_test` takes an `HGVSp` on the row of an insertion or
deletion to which both engines give one non-zero `HGVS_OFFSET` and whose terms hold
`5_prime_UTR_variant` and neither `start_lost` nor `start_retained_variant`: two values that name
one protein before their first `:p.`, an escaped `=` (`%3D`) read as `=` and a pair of prediction
parentheses that both carry set aside. On a row with no `CDS_position`, Ensembl VEP's value is a
frameshift at residue 1 with the unknown count and vep-rs's that residue unchanged. On a row that
holds `coding_sequence_variant` and whose `CDS_position` has an end and no start, Ensembl VEP's
value is that frameshift and vep-rs's the lost start named at residue 1, to the residue of the last
coding base of vep-rs's `c.` deletion; or Ensembl VEP's is a deletion from residue 1 to residue n
and vep-rs's the lost start over those n residues, each named. Ensembl VEP reads the residues on
the allele moved 3' and the start codon test on the allele as the record writes it
(`TranscriptVariationAllele.pm` lines 1796 and 2161). The row is an insertion's or a deletion's by
its Allele and Location, and two or more bases over as many positions by vep-rs's `dup` or `ins`
notation. The pairs are those of `tools/vep-formats-compare/src/rules/shifted_start_test.rs`, so
the class has one definition for both tools. A run names the class only where it compares
`HGVS_OFFSET` beside `HGVSp`, and holds a lost start to its last residue where it compares `HGVSc`
too.

`hgvsp_start_codon_test_past_the_shifted_allele` takes an `HGVSp` on the row of an insertion or
deletion to which both engines give one positive `SHIFT_LENGTH`, whose terms hold
`start_retained_variant` and no `start_lost` and whose `CDS_position` has a start and an end: two
values that name one protein before their first `:p.`, an escaped `=` (`%3D`) read as `=` and a pair
of prediction parentheses that both carry set aside. vep-rs's value is the first residue unchanged,
and Ensembl VEP's is a frameshift at residue 1 with the unknown count and the same first residue, or
an insertion between a residue 0 and residue 1 whose two flanks carry one name. Under
`--shift_3prime 1` Ensembl VEP makes the second start codon test of `start_lost` the shift length
past the allele (`Utils/VariationEffect.pm` lines 920 and 928-930) and writes the lost start only
under that term (`TranscriptVariationAllele.pm` line 2161). The row is an insertion's or a
deletion's by its Allele and Location, and two or more bases over as many positions by vep-rs's
`dup` or `ins` notation. The pairs are those of
`tools/vep-formats-compare/src/rules/start_test_past_the_allele.rs`, so the class has one definition
for both tools. A run names the class only where it compares `SHIFT_LENGTH` beside `HGVSp`, and
reads two or more bases over as many positions as an insertion's where it compares `HGVSc` too.

The stop codon classes take an `HGVSp` whose two values name one protein before their first `:p.`,
an escaped `=` (`%3D`) read as `=`. `hgvsp_stop_codon_protein_comparison` is Ensembl VEP's
description read from the codon window against the one the two proteins give, on a row whose `FLAGS`
hold no `cds_end_NF`: the stop replaced by a stop against the stop unchanged at the same residue
above 1; on a `stop_retained_variant` row, one residue and a stop inserted between two neighbouring
residues, the residue the second neighbour's, against the stop unchanged at the residue after that
neighbour; on a row of any terms, two or more residues and a stop inserted between two neighbouring
residues against the stop unchanged, where the coding sequence of the row's transcript in the
corpus cache, under its own genetic code, spells the first neighbour at its number, the inserted
residues from the next number to the last and a stop codon after them at the unchanged stop's
number, both values alone or both inside prediction parentheses; on a `stop_lost` row, a range from a residue to the stop replaced by residues against a
frameshift named at that residue; and, against an extension named at the stop, the stop replaced by
the residues that follow the protein, residues inserted before the stop, or, on a `stop_lost` row,
the stop unchanged, the stop deleted before the extension's count, or the same extension with an
unknown count. On a `stop_lost` row of any `FLAGS` the class also takes one deletion with an
extension's count after it against the same deletion with another count, where the row's
`CDS_position` and its transcript in the corpus cache give both counts: the transcript's peptide
holds a stop and is the translation of its coding sequence; the coding bases are three for each
residue the deletion names; with them out, and the 3' UTR read after the coding sequence,
the first residue that differs is the deletion's first; the first stop of the reference at or after
it is a stop the peptide holds or the stop codon that closes the coding sequence; and the first
stop of the alternate at or after it stands vep-rs's count of places on. Ensembl VEP's count is
the one `_stop_loss_extra_AA` gives on those bases (`TranscriptVariationAllele.pm` lines 2503, 2513,
2519 and 2483). Where Ensembl VEP writes no value the class takes the stop unchanged, on a
`stop_retained_variant` row whose `Protein_position` names one codon, whose `CDS_position` lies
inside it and whose `Amino_acids` are its residue against that residue and then a stop, at the
residue after that codon; a value inside prediction parentheses is not read there.
`hgvsp_stop_named_on_incomplete_coding_sequence` is a description that names the stop
alone (unchanged, deleted, replaced by one residue or extended), inside one pair of prediction
parentheses or without them, where vep-rs writes none, on a row whose `FLAGS` hold `cds_end_NF`,
whose terms hold `incomplete_terminal_codon_variant`, whose description is `Ter1=` beside
`cds_start_NF`, or whose transcript's `translateable_seq` in the corpus cache is no whole number of
codons. A residue number is a decimal number of any length.

`hgvsp_frameshift_with_no_changed_residue` takes an `HGVSp` that Ensembl VEP writes and vep-rs
leaves empty, on a row whose terms hold `frameshift_variant`: Ensembl VEP's description is a
frameshift with the unknown count that names one residue, no stop, on both sides of its residue
number (`Leu11LeufsTer?`), inside one pair of prediction parentheses or without them. The class
reads the row alone and is asked after every class but the last, so a row an earlier class takes
keeps that class.

`hgvsp_frameshift_count_at_the_first_residue` takes an `HGVSp` whose two values name one protein
and a frameshift at residue 1, Ensembl VEP's with the unknown count (`fsTer?`) and vep-rs's with
a number, on a Transcript row of a transcript on which no start codon stands. The rule is
`frameshift_counted_at_the_first_residue` of `tools/vep-formats-compare/src/coding_end_terms.rs`,
a source this tool compiles, and this tool keeps no copy of it; that function's doc comment
states every condition. `classify` hands it the record of the row's Feature in the corpus cache
whose span shares a base with the row's Location, as the tuple comparator reads a record: the
rule takes a row where that record's start end reads as not shown (`StartEnd::NotShown` of
`tools/vep-formats-compare/src/coding_ends.rs`), with `cds_start_NF` or without it. It hands it
the key's Allele and consequence set, which holds no `start_lost`, the `CDS_position` both rows
print, and the count of reference bases the Allele stands in for, read by position from the
run's records; a row that two records read with different reference bases is not read. A
deletion whose `CDS_position` is open at its start is read at the span its `HGVS_OFFSET` moves
it to, so a run names the class on such a row only where it compares `HGVS_OFFSET` beside
`HGVSp`. The class takes the pair where the two values are the two texts derived from the cached
coding sequence: vep-rs's with the count of residues from residue 1 to the first stop codon of
the new reading frame, Ensembl VEP's the same text with the unknown count. The reading of the
corpus cache hands the rule no 3' UTR, so the class takes a row of a transcript that ends at its
coding end, and the same pair on a transcript with a 3' UTR is an `unexplained_residual`.

The intron classes take an `HGVSp` Ensembl VEP writes and vep-rs leaves empty on a row whose two
`HGVSc` are one notation of the row's transcript, so a run names them only where it compares
`HGVSc` and `HGVS_OFFSET` beside `HGVSp`. `hgvsp_intron_end_insertion` is a notation `c.N-1dup`,
`c.N-K_N-1dup` (`K` at least 2), `c.M+J_N-1dup` (`M` one below `N`) or `c.N-1_Nins` with its bases,
on rows without an `HGVS_OFFSET`, at a gap in which an intron of more than 13 bases ends before a
coding exon base. `hgvsp_intron_start_insertion` is a notation `c.N_N+1ins` with its bases, on rows
without an `HGVS_OFFSET` or with one non-zero integer on both, at a gap in which a coding exon ends
before an intron, the Location read as printed and as moved by the offset. A gap is the lower of
the two positions either side of it. A Location names it as two adjacent positions, or as the span
of a record that carries another allele beside the insertion: there the Allele, in the upper case
`str.upper()` gives, is the span's bases with the inserted ones on either strand, which names the
gaps its clippings leave, or the insertion alone on either strand (the inserted bases as they stand
or turned by the offset, or as many bases as a `dup` notation duplicates), which names every gap of
the span.
The exons and the coding region are a cache transcript's `exons`, `strand`, `cdna_coding_start` and
`cdna_coding_end`, each number read as `int()` reads a number, a boolean, a double or a text; the
chromosome is looked up as the Location writes it and without a leading `chr`, the Feature as
written and without a trailing version. `maxentscan_overlapping_windows` is a `MaxEntScan_` field
on a row whose Location is one position in the acceptor window of one intron of the row's
transcript (from the intron's twentieth base from its end to three bases after it) and the donor
window of the next (from three bases before it to its sixth base), where the two values are finite
numbers as `float()` reads them, or for `MaxEntScan_ref_seq` and `MaxEntScan_alt_seq` Ensembl VEP's
nine bases against vep-rs's twenty-three; the positions are read only where a `MaxEntScan_` field
is compared.

`loftee_gerp_walk_without_conservation_database` takes a field of the LoF plugin (`LoF`,
`LoF_filter`, `LoF_flags`, `LoF_info`) that Ensembl VEP leaves unwritten and vep-rs writes, on a
key whose consequence set holds `stop_gained` or `frameshift_variant`, where Ensembl VEP's row
carries none of the four and vep-rs's row carries `LoF` beside a `LoF_info` that holds a
`GERP_DIST` item, or beside no `LoF_info`, as a run under `skip_lof_info` writes the row. The
plugin Ensembl VEP runs on GRCh37 reads the GERP total of each whole exon its walk crosses from a
conservation database (`gerp_dist.pl` lines 59-60 and 129 on the `master` branch of the LOFTEE
repository) and dies on the row in a run without one. The four fields are read of both rows
whichever of them the run compares, and the class is asked before
`hgvsp_frameshift_with_no_changed_residue` and before the sequence-edge class, which is asked
last, so that it takes a field no other class takes.
The rule is the format comparator's class of the same name (`src/rules/loftee.rs` and
`src/rules/mod.rs` of `tools/vep-formats-compare`), which this tool keeps a copy of, of one
meaning (`src/classify/loftee.rs`).

`allele_number_repeated_or_reference_equal_alt` takes an `ALLELE_NUM`, or the `allele_num` the
JSON format writes, on a row of a record that writes an alternate twice, or an alternate equal to
the reference allele, before the row's allele. Ensembl VEP drops such an entry from the record's
allele list and numbers the entries left (`VariationFeatureOverlap.pm` lines 459-467 and 477-491
at release/116.3), where vep-rs numbers an allele by its place among the written alternates:
the `T` of `C/G/G/T` is 2 in Ensembl VEP's row and 3 in vep-rs's. The class takes the pair where
both default-format rows of the key carry the written list under `UPLOADED_ALLELE` with one text
of three or more parts, the key's Allele is a written alternate or, on an indel record, one
without the first base Ensembl VEP's VCF parser chops from the record's alleles (`Parser/VCF.pm`
lines 294-335: the `AT` of `CA/C/C/CAT`), vep-rs's number is the place of the first such
alternate and Ensembl VEP's is the count of the alternates up to that place that equal neither
the reference allele nor an alternate written before them. A row without the list, an Allele the
list holds under neither spelling, an Allele that is a written alternate at one place and the
parser's spelling of an alternate at another, an Allele whose alternate equals the reference
allele, and any other pair of numbers stay unexplained. The
list is read of both rows whichever fields the run compares, and the class is asked directly
before the LOFTEE class. The rule is the format comparator's own file,
`src/rules/allele_number.rs` of `tools/vep-formats-compare`, which `src/lib.rs` compiles by its
path, so this tool keeps no copy of it (`src/classify/allele_number.rs` asks it).

`spliceai_multi_entry_cell` is asked under `--spliceai-table FILE` alone (repeatable: the `snv=`
and the `indel=` score file of the run's SpliceAI plugin, each a bgzipped VCF with its tabix index
at `FILE.tbi`), and its answers are the format comparator's: the class is a rule of
`tools/vep-formats-compare`, whose `src/rules/spliceai.rs` states every condition, and this tool
keeps no copy of it. A run that names a score file and compares a field writes every differing
field that the eighteen classes asked before this one pass as one question, a line of the 30
tab-separated columns that tool's README states for its `field-class` command: the row's Location,
Allele and Feature, the field, Ensembl VEP's value and vep-rs's, then `SpliceAI_pred`, the nine
`SpliceAI_pred_` fields, `SpliceAI_cutoff` and `SYMBOL` as Ensembl VEP's row writes them and as
vep-rs's row does, `-` for a field a row lacks. The run starts `vep-formats-compare field-class
--class=spliceai_multi_entry_cell` once over all of them, each score file given as
`--spliceai-table` and the corpus cache as `--json-cache` where it holds a `transcripts`
directory, and the class takes the fields the command answers `1`. The command is started for a
run without a question as well, since it reads the index and the header of every score file
whatever it is asked. The command accepts a field whose name starts with `SpliceAI_`, in any
case, where a record line that names the row's gene holds several entries, vep-rs's row writes
every entry of the last line that names the gene which names it, in that line's order (as
`SpliceAI_pred`, or as the nine `SpliceAI_pred_` fields, each its values joined by commas) and
Ensembl VEP's row writes the whole value of the last line whose first entry names the gene,
beside a `SpliceAI_cutoff` that is the same or that one cutoff gives that line's first entry
where it gives the gene's entries vep-rs's, or nothing where no line's first entry names it; the
row's gene is the `SYMBOL` either row writes, else the `gene_symbol` of the
row's transcript in the corpus cache. In a run that names no output of a further format, after
the entries of a key that holds one of the class, `field_divergences` holds the entries the other
formats need, each a copy of that entry
under another `field`: one a `SpliceAI_pred_` member field whose two parts differ and the run does
not compare, then one for `spliceai`, each side's JSON object as
`json.dumps(object, sort_keys=True, separators=(",", ":"))` writes it, a member other than `SYMBOL`
a number where its text is `-?\d+(\.\d+)?` (an integer where its value is whole), the row's
`SpliceAI_cutoff` among the members, and `-` for a row without the nine members. A run that names
the VCF and JSON outputs reads those fields from them and derives none: it writes the same
entries, the member fields of a row in the byte order of their names where the derivation writes
them in the order of the value's members. Ensembl VEP's nine members are the first nine
`|`-joined parts of its `SpliceAI_pred`; vep-rs's, where its `SpliceAI_pred` holds several
comma-joined entries, are each the entries' values of that member joined by commas, the nine
fields vep-rs writes on a row whose gene several entries of one record name. A key on which both
engines write one `SpliceAI_pred` text and no field asked so far is a `SpliceAI_` field, while
those two sets of nine members differ, is one more question of the run: the first member field
that differs, with the two members as its values. Where the command accepts it, the key holds
that member field's entry and, in a run that names no output of a further format, behind it the
entries of the other member fields and of `spliceai`; a run that names one reads no further
format on such a key and writes that entry alone. Where the command does not accept it, the
member field is an `unexplained_residual`.

`hgvs_shift_sequence_edge_window` is asked under `--reference-fasta FILE` alone (the run's
reference sequence as an uncompressed FASTA, its index at `FILE.fai`), and its answers are the
format comparator's as well: the class is a rule of `tools/vep-formats-compare`, whose
`src/rules/sequence_edges.rs` states every condition, and this tool keeps no copy of it. The class
is asked last. A run that names a reference sequence and compares a field writes every differing
field that the twenty-six classes before it pass, under the answers of the SpliceAI class where the
run names score files too, as one question, a line of the 24 tab-separated columns that tool's
README states for the class: the row's Location, Allele and Feature, the field, Ensembl VEP's
value and vep-rs's, then the row's Feature_type, `STRAND`, `HGVSc`, `HGVSp`, `HGVS_OFFSET`,
`SHIFT_LENGTH`, `cDNA_position`, `CDS_position` and `Protein_position` as Ensembl VEP's row
writes them and as vep-rs's row does. The `STRAND` is Ensembl VEP's row's, each of `HGVSc`,
`HGVSp`, `HGVS_OFFSET` and `SHIFT_LENGTH` is `-` on both rows where the run does not compare it,
and the three position cells are the row's columns. The run
starts `vep-formats-compare field-class --class=hgvs_shift_sequence_edge_window` once over all of
them, the file given as `--reference-fasta` and the corpus cache as `--json-cache`, and the class
takes the fields the command answers `1`. The command is started for a run without a question as
well, since it reads the index and the cache's transcripts whatever it is asked, and a corpus
cache without a `transcripts` directory is refused before it starts. The command accepts an
`HGVSc` or an `HGVS_OFFSET` on a Transcript row of an insertion or a deletion less than 1,000
bases from the edge of its sequence that the transcript's strand runs towards, where each
engine's notation and offset are the ones its own shift leaves on the transcript: Ensembl VEP
shifts along a flank it cut from a slice clipped at that edge, vep-rs as far as the sequence
admits. It accepts the row's `HGVSp` as well where the allele shifts the reading frame and each
engine's value is the frameshift its own shift gives on the coding sequence the corpus cache
states. Each side's offset is read beside its notation or its `HGVSp`, so a run names the class
only where it compares `HGVSc` and `HGVS_OFFSET`, and `HGVSp` where a row without a notation
writes one. Under `--shift_3prime 1` the command accepts a `SHIFT_LENGTH` as well, where Ensembl
VEP's value is the length of the shift along its cut flank and vep-rs's the length the sequence
admits, and a position column where each engine's cell is the one the transcript gives the span
its own shift moved. A corpus's own reference is stored bgzipped: the class reads a copy inflated beside a
copy of its `.fai` index, the file the format comparison of the same run reads.

**The format comparator.** `classify` asks the build of `tools/vep-formats-compare` that
`--formats-comparator PROGRAM` names; else the one the environment variable
`VEP_FORMATS_COMPARATOR` names, the variable `scripts/concordance/feature_oracle.sh` reads for the
same build; else the release build of the tree its own executable stands in,
`tools/vep-formats-compare/target/release/vep-formats-compare` under the nearest directory above
the executable that holds `tools/vep-formats-compare/`. A name is a path where it holds a `/` and
else a command of `PATH`, and an empty name is none. The build is looked for by a run that names a
score file or a reference sequence and compares a field, before the fields are read, and by no
other run. A name that is
no executable file is refused with no other build tried in its place, and so is a tree without
the release build: one line states what was looked for and the command that builds the tool,
`cargo build --release --locked --manifest-path tools/vep-formats-compare/Cargo.toml`, and nothing
is written. The tool builds nothing. The questions and the answers are files of a directory of the
run's own under the system's temporary directory, which leaves with the run.

Two of the class modules read the transcripts of the corpus cache as translated sequences. The
organelle class takes an `HGVSp` on a transcript whose `codon_table` in the cache, a number or a
text of digits, names a genetic code under which Ensembl VEP reads some codon unlike the standard
code (0 and 1 name the standard code; tables 11 and 24 read every codon as it does, and so does
code 26 as Ensembl VEP reads it; a table no genetic code is held for names no such transcript),
where the two values name one protein before their last `:p.` and the difference lies on a path
Ensembl VEP translates with the standard code: any row with
`frameshift_variant`; a `stop_lost` row whose two descriptions hold `ext`; an `inframe_insertion`
or `protein_altering_variant` row on which one description ends in `dup` and the other holds an
`ins` that no `del` stands before. The Downstream classes take a `DownstreamProtein` or a
`ProteinLengthChange` on a row with `frameshift_variant` and no term that names a splice site, on
which both engines wrote a value and print one `CDS_position` without a gap, where the transcript's
`translateable_seq` and `three_prime_utr` give exactly the two values under its `codon_table` (1
where the cache names none or 0): vep-rs's from the codon that holds the first changed base of the
row's own span, Ensembl VEP's from the plugin's own codon boundary at the span its HGVS shift
leaves, the record's span moved by the row's `HGVS_OFFSET` through the transcript's `mapper` pairs
between `cdna_coding_start` and `cdna_coding_end`. `plugin_hgvs_shifted_cds_coordinates` is a value
the shift changes, with `second_expected_divergence` naming the class of the row's own span where
that changes the value as well; `downstream_utr_frame_past_coding_end` and
`downstream_start_past_own_codon` are a `DownstreamProtein` on a row whose own span opens on a
codon's third base. A `ProteinLengthChange` is read against the length of the cached `peptide`. A
field of the plugin that takes none of the three is asked of no further class. Sequences are
translated as `Bio::Tools::CodonTable` translates one, in the upper case `str.upper()` gives, their
bases counted as `len()` counts characters; a number of the cache is a number, a boolean, a double
cut to its whole part or a text `int()` reads, and a transcript with another value for a coordinate
is left out.

**Error lines and exit codes.** A run ends with exit code 0 after its work. A run the tool refuses
writes `ERROR: [build_golden_corpus] <message>` to standard error, the name in brackets a fixed part
of the line, and ends with exit code 1. A run that goes on after it names what it discards writes
`WARNING: [build_golden_corpus] <message>` there. Where the message quotes a value it is the value's `repr()`,
and where it carries the reason a file could not be read it is the text of the CPython exception:
`[Errno <n>] <strerror>: '<path>'` for an operating-system error, the `UnicodeDecodeError` of the
`utf-8` codec for bytes outside UTF-8, the `json.JSONDecodeError` with its line, column and
character for a text that is not JSON. A refused command line ends with exit code 2.

**Outside the contract.** The tool refuses, with one line `ERROR: [vep-golden-corpus] <detail>` and
exit code 1, what CPython reads or ends on a traceback for: an argument that is not UTF-8; a
manifest that holds a surrogate escape without its partner (`\ud83d` alone), nesting past 1,000
levels or an integer of more than 4,300 digits; a manifest file that cannot be written; standard
output that cannot be written; a record line of fewer than five columns; a record whose POS is not
an integer in ASCII digits that 64 bits hold; a multi-allelic record with an empty allele among
alleles that do not share their first base; a `--cache` whose `transcripts/` cannot be listed or
holds a directory that cannot be listed, and a shard a cost reaches that is no JSON list of records
with a
`start` and an `end` that are whole numbers of 64 bits; under `--keep-records`, `records` that are
not a list of objects with a `dataset` name and a `vcf` line; under `check`, a manifest without such a
list, a record without `dataset` that a row or a record class asks for, `combinations` that are no
list of names, and a `seed` that is a list or an object beside a `record_class_key`; a corpus
directory, `variants.vcf` or manifest that cannot be written; under `identifiers` and the first
read of `check`, a directory that holds a file whose name is not UTF-8, and under `identifiers` a
file or a `provenance.json` that cannot be written and a `provenance.json` that is not UTF-8 text;
under `generated`, an executable that stands in no tree that holds `tools/vep-golden-corpus`, a
scratch directory that cannot be made, a corpus directory that cannot be resolved and a written
file whose name is not UTF-8.
`classify` refuses in the same way a
manifest without a list of `records` or with a record that holds no `vcf` line; a manifest's `flags`
that is a number; a run's `input` that is no text; a `divergences` list of the run that holds a
number, a boolean or `null`, or of the manifest, for a flag set's run, that holds anything but
entries, and an entry of either whose Location, Allele, Feature, records or Ensembl VEP set is a
list or an object; a line of the reference's JSON output and a shard of the corpus cache with a lone
surrogate escape or nesting past 1,000 levels, and such a line with an integer of more than 4,300
digits; a shard whose deflate data does not inflate or that keeps a number where the records of a
chromosome stand; a directory or shard of the cache whose name is not UTF-8; an output that is a
pipe or a socket, which can be read once and not again. Eleven readings differ from CPython's
without a refusal: a path named to `identifiers` that the system does not let the tool examine is
neither a file nor a directory; a count with a digit outside ASCII, or one that 64 bits do not hold,
is an invalid count;
an argument that starts with `-` and a digit outside ASCII is read as an option; `repr()` writes a
code point Unicode 16 leaves unassigned as it is; the coordinates of a Location are numbers when
they are ASCII digits, and within 64 bits where a span's width, a substitution's place or a cache
cost is read from them; letters are compared in the upper case of the Unicode version the tool is
built with; a shard name and the position of a breakend's mate are read in ASCII digits, a number
past 64 bits as the largest; a chromosome whose name holds a digit outside ASCII is ordered among
the names that are no number; a frequency written in digits outside ASCII is no number; a residue
number of an `HGVSp` description is a number when its digits are ASCII; and such a number is read at
any length, where CPython ends on an error past 4,300 digits.

The transcripts of the corpus cache that `classify` reads for `breakend_mate_local_read` are read
as CPython reads them: a start, an end, a strand, an exon's ends, the mapper's coding bounds and
its pairs as `int()` reads the value (an integer, a float cut at its point, `true` as 1, a text
between blanks with a sign, underscores between digits and the decimal digits of any script), a
record without one of the first three or with another kind of value there holding no model. The
tool refuses, on a line of its own, a record on which CPython ends on a traceback (an infinite
number; `exons`, `attributes` or the mapper's `pairs` that are no list of objects; a
`variation_effect_feature_cache` or a `mapper` that is no object; a `miRNA` attribute whose range
has more than 4,300 digits; a `stable_id` that is a list or an object) and a number 64 bits do not
hold, which CPython reads. Two readings differ without a refusal, where CPython ends on a
traceback: a mate position or a Location of more than 4,300 digits, and a Location written in
digits `int()` does not read, name no mate-side row. The decimal digits are those of the Unicode
version the tool is built with.

The HGVS field classes read a digit as one of the ten ASCII digits, where CPython's `\d`,
`str.isdigit()` and `int()` read every decimal digit of Unicode, and a number within 64 bits: a
notation, an offset, a `DISTANCE`, a `cDNA_position` or a Location with a longer number names no
class, and neither does a Location whose end is a number and whose start is not, on which CPython
ends on a traceback. `classify --field` refuses in the tool's own words a record of the corpus
cache whose `attributes` is `true` or a number other than zero, which CPython cannot walk, and an
`_rna_edit` value with a coordinate that 64 bits do not hold.

The stop codon and intron classes refuse in the same way, where CPython reads on or ends on a
traceback: a transcript record of the cache whose `variation_effect_feature_cache` is neither an
object nor empty, or whose `translateable_seq` is a number. Three readings differ without a
refusal: every number these classes read (of a Location, an `HGVS_OFFSET`, a notation, a
description, a Feature's version), and a member of a SpliceAI value that the JSON object writes as
a number, is one in ASCII digits; a number of a Location or an `HGVS_OFFSET` of more than 38
digits takes no intron class; and a coordinate of a cache transcript that is infinite or lies past
64 bits is no integer, so that the transcript has no gap.

`spliceai_multi_entry_cell` refuses what the format comparator's command refuses, with exit code 1
and nothing written. A run on which the command ends with exit code 2 on an error line of its own
(`ERROR: [compare_vep_formats] <detail>`: a score file, a cache shard or a question it does not
read) ends on `ERROR: [vep-golden-corpus] <detail>`. The detail is the command's but for two this
tool words: an index the command cannot read is `--spliceai-table: <index> is not a readable tabix
index (<reason>)`, the reason the one this tool's own reading of the index as a gzip file gives,
or `not a tabix index` for a file that reads; and a block of a score file that does not inflate is
`--spliceai-table: <file> holds a block that does not inflate`. A run on which the command ends
in another way ends on `the format comparator <build> ended with exit code <n>: <line>`, the line
the last the command wrote, an error line of its own without its opening, with `was ended by a
signal` for a run without an exit code and `and wrote no error line` for one that wrote none. The
tool refuses in the same way a build the system does not start, a command that ends with exit
code 0 and leaves no answers, answers that are not one `0` or `1` a question, and questions it
cannot write under the temporary directory.

The Downstream and organelle classes refuse in the same way, where CPython reads on or ends on a
traceback: a transcript record of the cache whose `variation_effect_feature_cache` is neither an
object nor empty; a `stable_id` that is a list or an object on a record with a coding sequence or
an organelle's codon table; a `codon_table` text of digits that holds one outside ASCII or more
than 4,300 of them; a `cdna_coding_start`, a `cdna_coding_end` or a number of a mapper pair that is
infinite, lies past 64 bits or is written in digits outside ASCII; and a record whose span lies
past 64 bits. One reading differs without a refusal: a `CDS_position` or an `HGVS_OFFSET` in digits
outside ASCII, or of more than 4,300 digits, takes no class.

## Layout

- `src/main.rs`, `src/lib.rs`: the binary, and the run from the arguments to the exit code.
- `src/cli.rs` and `src/help/`: the command line as `argparse` reads it, the options of the seven
  subcommands and their help texts.
- `src/error.rs`: the failure that ends a run with one line and an exit code, and the lines of
  standard output and standard error.
- `src/json.rs`: the JSON value (`Value`, `Object`, `Int`), its reader and its writer.
- `src/pytext.rs`: the CPython rules of the error lines and of the files read as text, and the
  SHA-256 digest a `provenance.json` records a file under, computed by the `sha2` crate;
  `src/lib.rs` takes `tools/vep-outputs-compare/src/pyfmt.rs` as its `pyfmt` module for the rules
  both tools follow.
- `input` and `rules`: `src/lib.rs` takes `tools/vep-outputs-compare/src/input.rs` and
  `src/rules/` of that tool as these two modules, under the names those files use for each other,
  so the rules of a divergence class have one definition for both tools.
- `genetic_codes`: `src/lib.rs` takes `tools/vep-formats-compare/src/genetic_codes.rs` as this
  module, the genetic codes by table id and the codes that read every codon as the standard code
  does, which the organelle class asks.
- `allele_number`: `src/lib.rs` takes `tools/vep-formats-compare/src/rules/allele_number.rs` as
  this module, the rule of the class `allele_number_repeated_or_reference_equal_alt`, which the
  field class of that name asks.
- `shifted_start_test`: `src/lib.rs` takes
  `tools/vep-formats-compare/src/rules/shifted_start_test.rs` as this module, the protein notation
  pairs of an insertion or deletion shifted into the first codon, which the class
  `hgvsp_shifted_indel_written_span_start_test` asks.
- `start_test_past_the_allele`: `src/lib.rs` takes
  `tools/vep-formats-compare/src/rules/start_test_past_the_allele.rs` as this module, the protein
  notation pairs of an insertion or deletion shifted onto a start codon it leaves standing, which
  the class `hgvsp_start_codon_test_past_the_shifted_allele` asks.
- `src/constants.rs`: the manifest's fixed values, the mask classes, the keys `classify` writes, the
  columns of the two file formats, the genotype and record classes.
- `src/manifest.rs`: the shape check, the manifest's reader and writer, the parity matrix and the
  readings of a manifest the subcommands share.
- `src/covers.rs`, `src/flag_set.rs`: the two subcommands that keep the manifest's lists.
- `src/generated.rs`: the `generated` subcommand, with the run of a command from the executable's
  tree and the digests of the files it writes.
- `src/rows.rs`, `src/names.rs`, `src/classes.rs`: the rows of a default-format output and the
  lines of an input VCF; the names a record carries in an output; the variant, family, genotype and
  record classes. `src/rows/text.rs`: a plain or gzip file read as CPython reads it, and the reason
  one that cannot be read gives. `src/rows/text/bgzf.rs`: the blocks of a BGZF file inflated on
  several threads.
- `src/select.rs`, `src/check.rs`: the two subcommands that read reference outputs, with the scan
  of the inputs for the records of a corpus.
- `src/identifiers.rs`: the two shapes of identifier under restricted redistribution terms, their
  read in a text and the fabricated form written in their place. `src/identifiers/stored.rs`: the
  files of a corpus read for them, and the counts across the files of one directory.
  `src/identifiers/gzip.rs`: a gzip file held in memory, read as `gzip.decompress` reads one and
  written as zlib writes one member at level 9. `src/identifiers/command.rs`: the `identifiers`
  subcommand.
- `src/sampling.rs`, `src/sampling/random.rs`: the pass that samples the rows of the references,
  the reservoir step and the choice of exemplars; the generator of `random.Random(seed)`.
- `src/cache_cost.rs`: the cost of a row's transcripts in the fixture cache.
- `src/rows/blocks.rs`: a reference read a block of whole lines a thread, the blocks taken in file
  order.
- `src/classify/`: `mod.rs` the subcommand; `outputs.rs` and `records.rs` the two outputs and the
  records of a run; `formats.rs` the rows of the tab, VCF and JSON outputs; `cache.rs` the corpus
  cache; `fields.rs` the question a field class is asked and
  the dispatch that asks the classes in order; `keys.rs` the classes of a pair of consequence sets;
  `intron_span.rs` the class of a whole intron inside the coding span of an edited model;
  `tandem_repeat_loss.rs` the class of a tandem repeat's stated loss or stated run length and the
  records of a key several records share;
  `hgvs.rs`, `introns.rs`, `protein.rs`, `spliceai.rs`, `downstream.rs`, `organelle.rs`,
  `allele_number.rs` and `loftee.rs` the field classes; `total_length.rs` the JSON end keys of a
  run under
  `--total_length`; `comparator.rs` the build of the format comparator a run asks and one run of
  its `field-class` command. `hgvs/` holds what the HGVS classes read and decide: `notation.rs` the
  HGVSc patterns
  and the place of a cDNA coordinate, `edited.rs` the edited transcript models of the cache,
  `shifted.rs`, `clipped.rs`, `outside.rs` and `unrotated.rs` the classes of a shifted insertion
  or deletion, of a clipped reference, of added bases outside the transcript and of an inserted
  allele left unrotated.
- `src/mate_context.rs`: the mate-side derivation of a breakend's row.
- `src/scratch.rs`: a scratch directory, of a test, of the questions `classify` asks the format
  comparator and of the files a `generated` command writes.
- `src/*_tests.rs`, `src/rows/*_tests.rs`, `src/rows/text/bgzf_tests.rs`,
  `src/sampling/random_tests.rs`, `src/identifiers/*_tests.rs`, `src/classify/*_tests.rs`,
  `src/classify/hgvs/*_tests.rs`: unit tests of the module they are named after.
- `tests/fixtures/two_record/`: a reference output of two rows, its input VCF, and the corpus
  `select` writes for them with `--covers symbol`. `tests/fixtures/random/draws.tsv`: CPython's
  draws, one row a seed, kind of draw and argument, a number of ten digits or more written with
  `_` between its groups of three. `tests/fixtures/cache_cost/lengths.tsv`:
  CPython's costs of the transcript records of the committed corpus caches, one row a cache.
- `tests/command_line.rs`, `tests/covers.rs`, `tests/flag_set.rs`, `tests/manifest_shapes.rs`: the
  binary's exit code and bytes on the fixture corpus. `tests/generated.rs`: `generated` on that
  corpus, run from a copy of the binary inside a tree of the test's own, with shell scripts for
  commands. `tests/select.rs`, `tests/genotypes.rs`,
  `tests/strata.rs`, `tests/check.rs`: `select` and `check` on references and inputs each test
  writes, with their options, kept samples, family strata and record classes.
  `tests/stored_manifests.rs`: the committed corpora under `tests/golden/`, each manifest written
  back to its bytes, held to the shape check and to the mask classes of the golden harness.
  `tests/stored_corpora.rs`: the same corpora, the expected outputs of every run held to the
  digests of the corpus's `provenance.json`, the corpus to its assembly and release, and every
  coding transcript of its cache to a peptide; the shift-length set of `GRCh38-regulatory`, the
  reference slice of `GRCh38-everything-fasta` and the edited-model row of
  `GRCh38-refseq-intron-span` to what their runs are for; every table a
  provenance record lists under `synthetic_tables` to a file of its corpus; every plugin corpus
  that holds a file named `synthetic` to a directory name that ends in `-fixture`; the manifest and
  input of `GRCh38-identifiers` to what `select` writes over its committed reference and input; and
  six corpora whose selection the tree cannot run again to manifests without counts of the rows it
  chose. The same file holds every transcript of a committed cache that carries a UniProt
  identifier key to `uniprot_isoform`, and the `info.json` of every committed cache to the
  regulatory data, the variation columns and the sequences its slice holds, one corpus excepted
  whose record says what its file keeps. `tests/cache_cost_lengths.rs`: the caches
  `tests/fixtures/cache_cost/lengths.tsv` names, every record costed and held to the table. `tests/identifiers.rs`: the read
  and the replacement of the two identifier shapes on made-up identifiers, `select`, `check` and
  `identifiers` on inputs and corpora that carry them, and every committed corpus read as `check`
  reads it. `tests/common/` holds what the test files share.
- `tests/fixture_manifests.rs`: every JSON file under `tests/fixtures/`, none of which holds
  `suite`, `suites` or `multi_allelic_per_suite`, keys no subcommand reads: a record of a made-up
  corpus names its dataset under `dataset`.
- `tests/classify_keys.rs`, `tests/classify_records.rs`, `tests/classify_fields.rs`: the binary's
  `classify` on made-up outputs of a run: the keys that differ, the records without a row, and the
  fields that differ. `tests/classify_formats.rs`: the same on made-up outputs of every format, the
  fields a further format gives another name. `tests/classify_total_length.rs`: the JSON end keys
  of a run under `--total_length`. `tests/classify_allele_number.rs`: the allele numbers of a run
  under `--allele_number`, in the default format and under the name each further format gives the
  field. `tests/classify_support/` holds what the test files of `classify` share.
- `tests/classify_key_classes.rs`: the binary's `classify` on rows of the key classes and of their
  nearest shapes, and the classes it writes held to the tuple comparator's pass over the same rows.
- `tests/classify_intron_span.rs`: the binary's `classify` on a made-up row of
  `refseq_edited_model_intron_span_coding_reference`, on the nearest rows it leaves unexplained and
  on a key a set rule names, which keeps that rule's class.
- `tests/classify_tandem_repeat_loss.rs`: the binary's `classify` on made-up rows of four
  `<CNV:TR>` records at one position, whose rows share one key: the entry of
  `cnv_tr_symbolic_loss_named_a_gain` with the records whose rows differ, the nearest keys it
  leaves unexplained and a key whose rows agree, which writes no entry.
- `tests/classify_loftee.rs`: the binary's `classify` over the four fields of the LoF plugin on
  made-up rows of `loftee_gerp_walk_without_conservation_database`, with `LoF_info` and without
  it, on the nearest rows it leaves unexplained, and in a run that compares one of the four.
- `tests/classify_hgvs.rs`: the binary's `classify` on the runs of `tests/fixtures/hgvs_classes/`,
  one directory a run: the corpus `c/` with its manifest and, where a class reads it, its cache,
  the two outputs `vep.txt` and `rs.txt`, and `expected.json`, the manifest the run writes.
  `tests/fixtures/hgvs_notations/` holds the tables the unit tests of `src/classify/hgvs/notation.rs`
  read: cDNA coordinates with their places, pairs of coordinates with the bases between them, and
  HGVSc bodies with what each pattern reads of them.
- `tests/classify_introns.rs`, `tests/classify_protein.rs`, `tests/classify_spliceai.rs`: the
  binary's `classify` on the made-up runs of `tests/fixtures/intron_classes/`,
  `tests/fixtures/maxentscan/`, `tests/fixtures/stop_codon/`, `tests/fixtures/frameshift/` and
  `tests/fixtures/spliceai/`, each a
  corpus before its first call, the outputs the calls name and the manifest call `n` writes
  (`expected.<n>.json`). `tests/recorded_runs/` holds what the three share; `tests/score_files/`
  writes the score files of the SpliceAI runs, which run through a build of the format
  comparator. `tests/fixtures/cache_shards/` holds the cache shards the intron and MaxEntScan runs
  place in a corpus, `tests/fixtures/stop_codon/*.tsv` the pair tables of the two stop codon
  classes, and `tests/fixtures/frameshift/*.tsv` the pair table of
  `hgvsp_frameshift_with_no_changed_residue`.
- `tests/classify_downstream.rs`: the binary's `classify` on made-up outputs of runs under the
  Downstream plugin and on rows of organelle transcripts, each over a corpus cache of its own.
- `tests/formats_comparator.rs`: the binary's `classify` under score files with the format
  comparator named by the option, by the variable and by the executable's tree; the refusal of a
  run that finds no build; the one run of the command a run makes, with its arguments and its
  questions; and the line of a run whose command leaves no answers.
- `tests/comparison_tables.rs`: two files of `tools/vep-formats-compare` read by their paths. One
  test holds the tool's lists of the comparison's class names to the columns of
  `data/classes.tsv`: the corrected-behaviour classes of `src/constants.rs` are the rows with a
  place in `REGISTRY_CLASSES`, `FIELD_CLASSES`, `PLUGIN_WARNING_CLASSES` or `MANIFEST_CLASSES`
  whose `CORRECTED_IN` cell names no release up to Ensembl VEP 116.3,
  the recoder classes of the tool corpora are the column `RECODER_CORPUS_MASKS`, and the classes
  the field dispatch names have a place in `FIELD_CLASSES` but for four that table lacks. The
  other holds the tool's own copy of each item of `tests/data/builder_parity/cases.tsv` to the
  table's answer over every case, 24 items; the one item the tool keeps no copy of is named there.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-golden-corpus/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

The root `.cargo/config.toml` applies, so a release binary carries the same `target-cpu` floor as
`vep`. The build compiles zlib from the source of the `libz-sys` crate, so it needs a C compiler
and nothing of the host's zlib. `cargo test` reads `tests/golden/` and
`crates/vep-cli/tests/common/mod.rs` of the repository beside the package's own fixtures, and
`data/classes.tsv` and `tests/data/builder_parity/cases.tsv` of `tools/vep-formats-compare`. The
tests of `classify` under score files (`tests/classify_spliceai.rs`, `tests/formats_comparator.rs`)
ask a build of that tool, a release build and else a debug build, under `CARGO_TARGET_DIR` where
the variable is set and else in the tool's own `target/`. They build nothing, and without a build
they fail on a line that names the build command; CI's Test job runs that tool's tests, which
build it, before this package's. The unit
tests of the tuple comparator's `input` and `rules` run with this package's as well, and read what
they read there: that tool's `data/` and `tests/fixtures/`,
`docs/concordance-provenance/116.2-divergence-classes.md` and
`scripts/parity/gvcf_site_rows_fixture.tsv`. The files under `src/help/` and `tests/fixtures/` are
the contract in bytes: nothing regenerates them, and a change of the contract changes them in the
same commit.
