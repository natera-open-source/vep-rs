# vep-sv-compare

The structural-variant comparator of the concordance harness: tuple concordance between Ensembl
VEP and a scored engine over a directory of VCF files and the two engines' default-format outputs,
with per-file, per-type and pooled counts, the divergence masks of a reference release and four
report files.

```text
vep-sv-compare [-h] [--assembly {grch37,grch38}] [--input-dir DIR] [--perl-dir DIR]
    [--rust-dir DIR] [--output-dir DIR] [--vep-rs-cache DIR] [--scored-engine NAME]
    [--reference-release {115.2,116.2,116.3}] [--reference-max-sv-size N]
```

| Option | Default | What it names |
| --- | --- | --- |
| `--assembly` | `grch37` | `grch37` or `grch38`: the last component of the four default directories, and the assembly the markdown report names |
| `--input-dir` | `tests/sv_validation/<assembly>` | The directory of the input VCF files |
| `--perl-dir` | `tmp/sv_perl/<assembly>` | Ensembl VEP's outputs: `<name>.txt` for the input `<name>.vcf` or `<name>.vcf.gz` |
| `--rust-dir` | `tmp/sv_rust/<assembly>` | The scored engine's outputs, under the same names |
| `--output-dir` | `tmp/sv_concordance/<assembly>` | Where the reports are written; made with its parents |
| `--vep-rs-cache` | none | A vep-rs JSON cache: the authority for the chromosome, the span and the structure of a transcript |
| `--scored-engine` | `vep-rs` | The engine that wrote `--rust-dir` |
| `--reference-release` | `116.3` | `115.2`, `116.2` or `116.3`: the Ensembl VEP release that wrote `--perl-dir`, which selects the masks |
| `--reference-max-sv-size` | `10000000` | The `--max_sv_size` the reference ran with, read under `116.2`; `-1` is no cap |

An option is written `--option value` or `--option=value`, whole or as a prefix only one option
begins with, and the last of a repeated option stands. An empty `--vep-rs-cache` is no cache.
`-h` or `--help` prints the help, a paragraph on each option, and compares nothing.

A run under `116.3` is scored by every rule of `116.2`: what this page says of a run under
`116.2` holds for it, and its line on stdout and its reports name `116.3`. A run that names no
release is a run under `116.3`.

A run writes `concordance_report.json`, `concordance_report.md`, `discordant.tsv` (every one-sided
tuple of every file) and `discordant_adjusted.tsv` (the one-sided tuples the masks leave in place)
into the output directory. It prints to stdout the size of the cache when one is given, two or
more lines for each input file, the pooled rows, the row of each variant type, the files whose
adjusted F1 is above their raw F1, under `116.2` one line of every class's count, and the names
of the reports. In the reports and on stdout the labels `Perl` and `Rust` name the reference and
the scored engine.

## Contract

### Command line

The arguments are read in their order, as Python's `argparse` reads them. An option is named by
its whole name or by a prefix that only one option begins with, `--help` being one of them, and
takes its value from the next argument or from the text after `=` in the same argument. The next
argument is a value when it does not begin with `-`, is `-` alone, or names no option and either
begins as a negative number (`-1`, `-.5`) or holds a blank; after `=` any text is a value, an
empty one included.

A usage error is the first defect in that order: a prefix several options begin with (the
message lists them; `--=x` begins every long option and `-=x` every option), an option without
its value, an assembly or a release that is not offered, a size cap that is no integer. An
argument that names no option is unrecognised, as is every argument from a bare `--` on; they
are listed in one usage error once every option is read. A refused value is quoted as Python's
`repr()` writes a text, and a byte of an argument that is not UTF-8 is written `\udcNN`.

`-h`, `--help` or a prefix of `--help` prints the help on stdout where it is reached and ends
the run with exit code 0, whatever follows it: an option after it is not read, and an
unrecognised argument before it is not reported. The help option takes no value (`--help=x` and
`-h=x` are usage errors), and letters directly after `-h` are read as further one-letter
options, so `-hh` and `-hx` print the help as well. The help is the text `argparse` writes for
the options at a width of 80 columns and without colour, whatever the terminal and `COLUMNS`:
the usage block, a line on the program and a paragraph on each option.

### Inputs

The input files are the entries of `--input-dir` whose names end in `.vcf.gz` or, when there is
none, those whose names end in `.vcf`, a name that begins with a dot included. They are compared
in the order `sorted()` gives `pathlib` paths, the bytewise order of the names. An input's name
is its file name without its last suffix and without a `.vcf` left after that, and its two
outputs are `<name>.txt` under `--perl-dir` and `--rust-dir`. An output that does not exist is an
engine that wrote nothing: a `WARNING` line on stdout and no tuple. The comparison asks once
whether it exists; under `116.2` a mask that reads the output again asks again. Two inputs that
give one name (`x.vcf` and `x.vcf.vcf`) share every per-file entry of the reports, which stands
at the first file's place with the second file's data, while the pooled sets, the record count,
the raw file sums and `discordant.tsv` take both files and the adjusted file sums take one.

Every file is read as text with universal newlines (`\n`, `\r\n` and a lone `\r` end a line) and
as strict UTF-8. An input whose name ends in `.gz` is read as `zcat` writes it to a reader that
ignores how it ends: every member of the stream, and nothing from the point where the stream
cannot be opened or decoded, without a message. Such a file that is empty, is not gzip or is in
another compression format holds no record. A file or directory name that is not UTF-8 is
printed with replacement characters on stdout, in a refusal and in a failure that names a file;
a usage error and an error of the operating system write each such byte as `\udcNN`.

Under `116.2` the row-aware pair rules read an input a second time, where a row put to them
does not spell its alleles, and that reading takes a name that ends in `.gz` or `.bgz` as
Python's `gzip` module reads one: every member in turn, zero bytes after a member skipped, a
file of no bytes holding no record. A file that ends inside a member, bytes that begin no
member, an unknown compression method and a checksum or a length the data does not have end the
run in the module's words (`Compressed file ended before the end-of-stream marker was reached`,
`Not a gzipped file (b'..')`, `CRC check failed 0x.. != 0x..`), and deflate data that cannot be
decoded ends it with `<path>: corrupt deflate data`. A truncated `.vcf.gz` is therefore compared
up to its cut in a run where no rule reads it again, and ends a run where one does.

An input line that begins with `#` or has fewer than five columns holds no record, and a record
of fewer than eight columns has an empty INFO. POS is read as Ensembl VEP numifies it, so a
record whose POS is not an integer is kept and counted. In an output, a line that begins with
`#` and an empty line hold no row, and a row of fewer than seven columns ends the run, because
leaving it out would turn a truncated output into a smaller tuple set. Of a file's short rows
and its integers of more than 4,300 digits the one on the earliest line is reported. A byte
sequence that is not UTF-8 is met as Python's text reader meets it, 8,192 bytes of the file at a
time (of the decompressed bytes of a `.gz` input): it is reported before a short row or a long
integer on any line that ends within the same 8,192 bytes or after them, and after one on a line
that ends before them. Its detail is the text of Python's `UnicodeDecodeError` (`'utf-8' codec
can't decode byte 0xff in position 12: invalid start byte`), which names no file and counts the
position within the bytes being decoded. Of several files that cannot be read, the failure is
that of the first in the order of the comparison, an input before the reference's output of it
and that before the scored engine's; stdout then holds the lines of the files before it and its
own `Processing` line.

### Tuples

A row is compared as the tuple of its Location, Allele, Feature, Feature_type and Consequence
columns (the second, third, fifth, sixth and seventh), a tuple once however many rows carry it.

- The Location is compared with a leading `chr` of any case dropped from the chromosome, `M`
  named `MT` and each coordinate written as its number: `030425949` and `+30425949` are
  `30425949`. A coordinate that is no integer stays as written, as does a text without `:`.
- The Allele is compared with the `chr` after a breakend bracket dropped.
- The Consequence is compared as the set of its terms: each stripped of blanks, once, in
  code-point order.
- Under `--scored-engine fastvep`, the name in any case and with blanks around it, the scored
  side's symbolic ALT tokens are compared as the class names the other engines write: `<DEL>` as
  `deletion`, `<INS>` as `insertion`, `<INS:ME:ALU>` as `Alu_insertion` (the subtypes `L1`,
  `LINE1`, `SVA` and `HERV` likewise, any other as `mobile_element_insertion`, and a `<DEL:ME>`
  token the same with `_deletion`), `<DUP>` as `duplication`, `<DUP:TANDEM>` as
  `tandem_duplication`, `<INV>` as `inversion`, `<CNV:TR>` as `tandem_repeat`, `<CN0>` as
  `deletion`, `<CN2>` as `duplication`, any other `<CN...>` as `copy_number_variation` and `<CPX>`
  as `CPX`. `<NON_REF>`, `<*>`, `<BND>` and any other token stay as written, and the reference
  side is never mapped.

### Counts

Each file has a row of counts: its records, how many of their `chrom:POS` Locations each side
printed, each side's tuples, the tuples both sides hold, those one side holds alone, and the
precision, recall and F1 of the scored engine against the reference. A ratio over no tuples is
0, so two sides without tuples score an F1 of 0 and not 1. A ratio is the quotient of two counts
as doubles, exact up to 2^53 tuples, and a count is held in 64 bits.

The pooled row (`overall`) holds the sums of the files' counts, so that tuples are pooled the
way records are, and beside them the counts over the union of all files' tuples (`*_union`) and,
for each side, the file sum less the union (`cross_file_duplicate_*`). The adjusted pooled row
holds the sums of each file's counts after that file's masks and, beside them, the union's
counts after the union's masks; its `cross_file_duplicate_*` is below zero when the files' masks
set aside more than the union's masks and the tuples several files hold together.

The per-type rows group the records by a type read off REF and ALT: `SNV`, `MNP`, `Small_INS`,
`Small_DEL`, `Complex`, `Symbolic_DEL`, `Symbolic_INS`, `Mobile_Element`, `Symbolic_DUP`,
`Symbolic_INV`, `Symbolic_CNV`, `Complex_SV`, `BND`, `NON_REF`, `Spanning`, `RefOnly`,
`Multi_allelic` or `Other`. An ALT of more than one field is `Multi_allelic`, a trailing empty
field being none: `T,` is one allele and an ALT of commas alone is `Other`. A tuple has the type
of the record whose `chrom:POS` its Location is, the last such record where several share one,
else of the record whose ID is the first column of the first row that prints its Location, the
reference's rows before the scored engine's; a one-sided tuple with neither is `Unknown`. A
record's ID is the first of the `;`-separated identifiers of its ID column, the one an engine
writes in the first column: an empty one or `.` names no record, and where several records share
an ID the first gives its type. Letter and case tests follow the Unicode tables of the Rust
toolchain and of the `regex` crate.

### Masks

`--reference-release` selects the masks. They run, in the order below, on the union of all
files' tuples and on each file's tuples, and the adjusted counts are the counts of what they
leave. A mask marked (cache) sets nothing aside without the view of the cache it reads: the
transcripts of a chromosome, their spans or their structures.

Under `115.2`:

1. `filter_intended_divergences`: at a record Ensembl VEP 115.2 marks `vep_skip` (a span over
   10,000,000 bases, a symbolic type without a term), the tuples on transcripts the other side
   lacks: the scored side's only where it names more transcripts there, the reference's unless
   the scored side names fewer.
2. `filter_cnv_tr_expansion_divergences`: a tandem repeat's literal reading on the reference
   side and its symbolic reading on the scored side, where the two differ only in representation.
   It reads the cache for one term: a `splice_polypyrimidine_tract_variant` on the symbolic
   reading alone is representation where a transcript's structure shows that the run overlaps
   the tract of an intron and no exon, and that the literal reading's span overlaps no tract;
   the term stays charged without a cache.
3. `filter_cross_chromosome_divergences` (cache): the reference's tuples on a transcript the
   cache does not file under the Location's chromosome. The same shape on the scored side ends
   the run for vep-rs (`vep-rs named ...`, the first such tuple in tuple order) and is counted
   for another engine (`scored_engine_transcripts_absent_from_cache`). A Location without `:` is
   looked up under the empty chromosome name, and the failure then reads `chromosome None`.
4. `filter_cross_chromosome_orphan_intergenic`: the scored side's intergenic row that stands
   where every reference tuple fell to the mask before it.
5. `filter_giant_breakend_mate_divergences` (cache): the scored side's mate-allele tuples of a
   breakend whose span is over 10,000,000 bases, on transcripts within 5,000 bases of the mate.
6. `filter_registry_swap_pairs`: both members of a pair of consequence sets on one transcript
   that a pair rule of the registry of divergence classes takes.
7. `filter_non_ref_batch_divergences`: the reference's transcript tuples on a `<NON_REF>` block
   for which the scored side wrote the intergenic row alone, and that row with them.
8. `filter_breakend_mate_local_read_pairs` (cache): both members of a pair on a point breakend's
   mate-side row, where the reference's set is what the local coordinate derives and the scored
   set what the mate coordinate derives. A pair is one one-sided tuple a side at a Location,
   Allele, Feature and Feature_type: where a side holds two there, nothing is set aside.

Under `116.2`, masks 2, 6 and 8 of that list, in that order, and then the one-sided classes of
the 116.2 registry. Mask 6 puts a pair none of its rules takes to the registry's row-aware pair
rules, which read the rows behind the pair and, for a row under an identifier, the reference
bases of its input record, the one that carries the row's identifier where two or more records
read the row differently; a tuple one of whose rows has fewer than nine columns in any output
goes to no row-aware rule, and a row is read to its Codons column, the twelfth, where it holds
one.

Ten of the row-aware pair rules, and `start_cooccurrence_swap` among the rules of mask 6, ask
what stands at the two ends of the coding sequence of the pair's transcript, which they read
from the transcript's record in the `--vep-rs-cache` directory (`tools/vep-outputs-compare`,
Coding ends). With a cache, `structural_deletion_stop_lost_without_stop_codon` and the four
other classes of a coding end take a pair where vep-rs's coding terms are the ones derived for
the row, and a stop codon class or `start_cooccurrence_swap` takes a pair only on a transcript
whose record shows the codon the class stands on. A pair whose transcript the cache answers with
no one record at the pair's Location is taken by none of them. Without a cache the classes of a
coding end take no pair, and the classes that stood before take a pair by its two rows alone.
The transcripts of a sequence are read for this when a pair on it first asks, by the reader of
`tools/vep-outputs-compare/src/rules/cache_transcripts.rs`. That reader passes over a file under
a shard's name that cannot be read whole, is not UTF-8 or is no JSON array or object, and a JSON
array or object that is neither an array of transcript objects nor an object of such arrays ends
the run with a failure that names it. No reference sequence is handed to these rules, so the one
shape of row that reads past a transcript's last base stays a difference.

- `filter_reference_skipped_records`: the scored side's tuples of a record the reference
  dropped before its input buffer (see the cross-check below).
- `filter_bnd_own_rows_lost` (cache): the scored side's own-end tuples of a bracket breakend
  whose mate is on another chromosome, where the reference kept the mate rows alone, on
  transcripts within 5,000 bases of the breakpoint.
- `filter_bnd_synonym_mate_rows` (cache): the scored side's mate tuples of a bracket breakend
  whose mate chromosome is spelled in a way the reference does not resolve (`Chr21`, an
  accession). Under a mate name the cache files nothing for, every transcript that is not on
  the Location's own chromosome counts as lost, an own chromosome the cache does not file
  included.
- `filter_bnd_bare_reference_allele_rows`: the reference's second and later rows of one
  `chromosome_breakpoint` tuple in an output file. They are rows, not tuples, counted for each
  file and outside both figures.
- `filter_gvcf_variant_sites`: the scored side's tuples of a gVCF variant site (`G,<NON_REF>`),
  at the alleles its parser writes.
- `filter_dot_breakend_pairs`: the reference's two-coordinate tuple and the scored side's point
  of a dot-notation breakend without SVTYPE, as pairs on one feature under one allele. The points
  go to the first reference tuple that reaches them, the reference tuples taken in tuple order.
- `filter_missing_alt_annotated`: the reference's tuples of a record whose ALT is `.` under a
  true SVTYPE.
- `filter_non_integer_end_rows`: the reference's tuples at a Location with a coordinate that is
  no integer, and the scored side's tuples of the same record at the numified value. The start
  of such a Location is the text before its first `-`, so a negative start reads as 0.
- `filter_cnv_tr_empty_expansion`: the scored side's `tandem_repeat` tuples of a record whose
  expansion is empty or unchanged, where the reference wrote nothing beside the run. An `RN`
  that is true and splits to no field (`RN=,`) is an empty expansion, a repeat sequence whose
  count list has run out takes its unit and adds nothing, and a unit of no bases repeats to
  nothing whatever its `RUC` count, while one counted by `RB` leaves the record outside the
  class. The tuples are taken from all the scored side's tuples, so where the record's Location
  is no span a tuple both sides hold is among them.
- `filter_star_alone_svtype_annotated`: the reference's tuples of a record whose every ALT is
  `*` under a true SVTYPE.
- `filter_non_numeric_pos_numified`: the reference's tuples of a record whose POS is not a
  non-negative number, at the start the reference numifies it to.
- `filter_cnv_tr_rn_zero_allele_dropped`: the scored side's tuples of the allele the reference
  drops from a two-allele tandem repeat of which one allele repeats nothing. The size cap bounds
  the repeat sequence the reference expands: a cap of 0 reads as 5,000, and under a cap of -1 no
  repeat sequence is expanded.

Two shapes are counted under `116.2` and set nothing aside: the records whose only reference
row is intergenic while the scored side wrote transcript rows and no intergenic row
(`fasta_named_slice_intergenic_records`), and the records with a coordinate at or past
2,147,478,648 with the reference's tuples at their Locations (`bnd_mate_32bit_buffer_records`
and `bnd_mate_32bit_buffer_reference_tuples`; two records at one Location count its tuples
once). The coordinate is a bracket mate, a mate named by `CHR2` with `END2` (an empty `END2`
falls back to `END`, and a bare `CHR2` flag counts as present), or the start or end of a record
the reference keeps.

The classes read a record as Ensembl VEP's parser does, the spellings no dataset holds included:

- `NaN`, `Inf` and `Infinity` are read in any case, with the blanks of any script around them,
  and a dotted capital I (U+0130) or a dotless i (U+0131) reads as their `i`.
- A POS that is an integer once the blanks of any script are stripped is that integer; any
  other POS is numified as written, where only ASCII blanks are skipped: a no-break space before
  `30425949` leaves 30425949 and before `3e7` leaves 0. An END is tested for being a number with
  the blanks of any script stripped and is then numified as written: a no-break space before
  `2000` leaves the end 0.
- Where the slash-separated parts of a type name one type, the type is the first part as
  written: `<DEL> / DEL` gives `DEL ` with its blank, which has no term, and `<DEL>/DEL` gives
  `DEL`.
- A mobile-element token is read without regard to case under Unicode folding: a dotless-i
  spelling of `<ins:me:l1>` is a LINE1 insertion, and a dotted-capital spelling of `<INS:ME>`
  has no term. Upper-casing is the full Unicode mapping, so a dotless-i spelling of `<ins>` is
  an insertion to the scored engine, while fullwidth and accented letters stay themselves.
- In the test for an ALT written in a structural form, a word character is a letter, a number
  of any kind or `_`: a dot before a superscript two is such a form.
- A literal allele of a tandem repeat is a text of characters other than `a` to `z`, `_`, `.`,
  `<`, `>`, `[` and `]`.
- A bracket mate's position and the position of a point Location are read in the digits of any
  script and at any length, while the mate name the reference resolves is read with a position
  in `0` to `9`: a mate whose position is in another script is under no chromosome of the cache.
  A `CHR2` mate position is read in the digits of any script without a sign, and one written in
  characters that are digits without a decimal value (a superscript two) is no mate; a Location
  whose position is so written is no point.

### Skipped-record cross-check

Under `116.2` the input records and the outputs are held to each other, and a disagreement
stops the run before a report is written. A record is on a skip route by its own ALT and INFO:
a symbolic type without a term, a gVCF variant site, an incomplete deletion, a span that is not
a number or runs backwards, or a span over `--reference-max-sv-size`. The run stops with exit
code 1 and one line a problem on stderr, sorted and each once
(`ERROR: [compare_sv_concordance] skipped-record cross-check: ...`), when:

- the reference wrote a row under an own allele of a record on a skip route, which is a
  reference run under `--dont_skip`, of an earlier release or with another size cap;
- a skipped record and a kept one share a start, an end and an allele, so that neither's rows
  can be told from the other's;
- vep-rs wrote a structural row at a Location the reference wrote nothing for, which no record
  on a skip route and no pair class explains.

For a scored engine other than vep-rs the last kind is the engine's own disagreement: the rows
are counted (`scored_engine_skip_route_residue_rows`), stay in both figures, and the run goes on.

### Reports

`concordance_report.json` is one object of `overall`, `per_file`, `per_variant_type`, `adjusted`
and `per_file_adjusted`, laid out as CPython's `json.dump(report, fh, indent=2)` writes it and
ended by one newline: keys in a fixed order, strings in ASCII with any other character as a
`\u` escape of four hexadecimal digits, a ratio rounded to six decimals and written as
`float.__repr__` writes it. The name `--scored-engine` gives is held as text, a byte of it that
is not UTF-8 as U+FFFD, and `scored_engine` holds that text (`\ufffd` for such a byte).
`per_file` and `per_file_adjusted` are in the order the files were compared;
`per_variant_type` and every
by-file map are in the code-point order of their names, and a by-file map names only the files
with a count. The `adjusted` block states, beside the adjusted pooled row, every count the mask
accounting is rebuilt from: the two unions the adjusted counts subtract (`excluded_rust_extra`,
which `excluded_divergences` repeats, and `excluded_perl_extra`), each layer's own size, the
overlaps between layers, the size of the tandem-repeat class before its mask and the per-file
maps, an entry of a pair mask holding `rust`, `perl` and the count of each bucket or kind. Under
`116.2` the block also holds `reference_release`, `masks_applied`, `masks_not_applied` (the five
115.2 masks with the reason each is retired), the counts of each class,
`breakend_masks_check_ran` (whether a cache was given) and `divergence_classes`: every one-sided
class of the registry with its definition, its citation, whether this comparator sets its rows
aside (`excludes`) and how many it set aside (`excluded`: 0 for the two classes that are counted,
the sum of both sides for the dot-breakend and the non-integer-coordinate classes). Under `115.2`
none of those keys is written.

`concordance_report.md` holds the headline counts, a per-file table and a per-type table, both
in the code-point order of their names, every ratio as `f"{x:.4f}"` writes it. The adjusted
headline is written when a mask set a tuple aside, and it calls the two unions transcript
divergences whatever layers fed them. The `Adj F1` cell of a file is blank when its adjusted F1
equals its raw F1. A type whose F1 is above 0 and below 0.99 before rounding is marked `**`, so
a row can read `0.9900 **`; on stdout the mark `<--` stands at every F1 below 0.99, 0 included,
and the per-file summary lists a file only when its adjusted F1 is above its raw F1 by more than
0.0001. The line of the tandem-repeat class is written when its mask set a pair aside and names
the assembly in upper case (`GRCH38`). Under `116.2` two more lines name the release with the
masks applied and every class with its counts.

`discordant.tsv` and `discordant_adjusted.tsv` have a header and one tab-separated row per
one-sided tuple: the file's name, the variant type, `only_perl` or `only_rust` and the five
columns of the tuple as they are compared. The files stand in the order compared, and within a
file the reference's tuples before the scored engine's, each side in the order of its tuples as
texts. The adjusted table is the raw one without the rows a mask of their file took.

### Numbers

A coordinate or a count read from an input is held in 128 bits. A POS, an END or an SVLEN past
that is held at the bound (2^127 - 1, or -2^127) and a sum stays there; a number Perl reads as
infinite (`1e400`) is 2^63 - 1 with its sign. A Location with a coordinate past 128 bits is no
span and belongs to no class, a bracket mate position above 2^127 - 1 is past every bound and
out of every transcript's reach, and a span whose coordinates differ by more than 128 bits hold
is outside the tool's range. A Location coordinate, a POS or a fastVEP copy number of more than
4,300 digits ends the run (`Exceeds the limit (4300 digits) for integer string conversion ...`).
`--reference-max-sv-size` takes an integer of a magnitude of at most 2^127 - 1, written as
Python's `int()` reads one: blanks around it, a sign, single underscores between digits, the
digits of any script.

### Cache

A cache is the entries named `*.json` or, gzipped, `*.json.gz` of each directory under
`<cache>/transcripts`, the directories and the entries of each read in name order; a file
directly under `transcripts` is not read, and a directory name that is not UTF-8 names its
chromosome with replacement characters. A shard is an array of transcript records or an object
of such arrays, where a repeated key keeps its first place with its last array. Its text is read
as Python's `json` module reads one: `NaN`, `Infinity` and `-Infinity` are numbers, an integer
of more than 4,300 digits is none, and a lone surrogate escape is read as U+FFFD; a value may be
nested to any depth. A shard that cannot be read, is not UTF-8 or is not such a text holds
nothing, as does a gzipped one that is cut short, fails its checksum or its length, or is
followed by bytes other than zeros. A cache that yields no chromosome is refused.

A record is kept when its `stable_id` is true to Python: one whose id is a number or `true`
files its chromosome and no transcript. The transcript total on stdout counts the distinct ids
of each chromosome, such ids among them, numbers that are equal as one id (`7` and `7.0`, `1`
and `true`). A later record of one chromosome and id replaces the span only when it gives
`start` and `end` as integers, and the structure only when it also gives `strand`. An integer is
a JSON number or a text Python's `int()` reads (`12.9` is 12, `true` is 1, `false` is 0,
`" +1_2"` is 12, the digits of any script), of a magnitude of at most 2^127 - 1 and, as a JSON
integer, also -2^127; a larger one is absent. Where a list or an object is read (the exons, the
feature cache, its mapper, the pairs, the attributes), a value Python reads as false stands for
an absent one. A coding bound of zero or below gives no coding region, and of two mapper pairs
that hold a bound the later stands. The mature-miRNA attribute is read with the blanks and the
digits of any script.

A cache on which Python's reading ends in an exception ends the run with the exception's text: a
record whose id is an array or an object that holds something (in the wording of Python 3.14,
`cannot use 'list' as a set element (unhashable type: 'list')`), an infinite number where an
integer is read, a number or `true` where a list is read, a value that is no object where one is
asked for a key (`'NoneType' object has no attribute 'get'`), a mature-miRNA range of more than
4,300 digits, a gzipped shard whose deflate stream is corrupt (`Error -3 while decompressing
data`, without the detail zlib adds), a `transcripts` directory that cannot be listed. The first
such record in shard order decides; within a record the id comes first, then `start`, `end` and
`strand`, then, only when those three are integers, the exons, the mapper (its pairs under two
coding bounds above zero) and the attributes.

### Stdout and stderr

Stdout that is no terminal is held as CPython 3.14 holds the stdout of a script. The text of the
lines is gathered until 8,192 bytes of it are held, and those chunks are gathered in a block of
128 KiB, or of the size the file states for itself where that is larger, up to 8 MiB. The block
is written when the next chunk does not fit in it and when the run ends, so a pipe or a file
receives nothing before about 123 KiB of lines are printed or the run is over. On a terminal
every line is written as it is printed, and stderr is written line by line everywhere.

Where one file or pipe takes both streams (`2>&1`), the lines of a refusal therefore stand
before the stdout text that was held when they were written, and the line of a failure stands
after all of stdout, which is flushed before it.

A write that stdout or stderr does not take (a pipe whose reader has gone, a device with no
space left, a file at the file-size limit) ends no run by itself, and the reports are written
where the run gets that far. A
run that ends with a stream holding what it could not write exits 120 in place of its own exit
code, and a stdout in that state is named on stderr in the two lines CPython writes: `Exception
ignored while flushing sys.stdout:` and the error as Python names it (`BrokenPipeError: [Errno
32] Broken pipe`, `OSError: [Errno 28] No space left on device`). A write of stdout that fails
in mid-run, which takes more than a block of lines, is a failure of the run: stderr holds
`ERROR: [vep-sv-compare] [Errno 32] Broken pipe` and then the two lines, and the exit code is
120. No report is written then, unless the write that fails is of the closing lines, which
follow the reports. The help and a usage error are written whether their stream takes them or
not. A stdout or stderr that is closed is the null device to the tool.

The signal of a write past the file-size limit (`ulimit -f`) is ignored from the start of the
run, as CPython ignores it, so such a write fails with `[Errno 27] File too large` and the file
holds its bytes up to the limit. A stdout or stderr file at the limit is a stream that does not
take what the run wrote (`OSError: [Errno 27] File too large` in the two lines, exit code 120).
A report that passes the limit is a failure of the run with exit code 1 and `ERROR:
[vep-sv-compare] [Errno 27] File too large` on stderr: the reports are written in the order
JSON, markdown, `discordant.tsv`, `discordant_adjusted.tsv`, the one that passes the limit is
cut at it, and those after it are not written.

### Exit codes

- 0: the reports are written, or the help is printed.
- 2: a usage error. Stderr holds the usage block that opens the help and then
  `vep-sv-compare: error: <detail>`, the detail worded as `argparse` words it.
  `--reference-max-sv-size` refuses an integer of a magnitude above 2^127 - 1.
- 1: the comparison refuses its inputs, each refusal one line on stderr that begins
  `ERROR: [compare_sv_concordance]`: an input directory that does not exist, a cache that yields
  no transcripts, no VCF file in the input directory, or the skipped-record cross-check. The
  output directory is made before the cache is read and the inputs are listed, so each of the
  last three leaves it behind.
- 1: any other failure, as one line `ERROR: [vep-sv-compare] <detail>` on stderr: a short row, a
  byte sequence that is not UTF-8, a file that cannot be read or written, a cache on which
  Python's reading ends in an exception (the detail is that exception's text), a scored tuple
  the cache rules out, or a `VEP_SV_COMPARE_THREADS` that is not a whole number of 1 or more. An
  error of the operating system is worded as Python words an `OSError`: its number, the system's
  description and the path the failed call named, quoted as `repr()` writes a text (`[Errno 13]
  Permission denied: 'out/discordant.tsv'`). An error met while an input or an output is being
  read names no path, and one met in a second reading by a mask names the file.
- 120: in place of 0, 1 or 2, when the run ends with stdout or stderr holding text the stream
  did not take (Stdout and stderr above).

A panic is a defect of the tool and no answer to an input. The release build aborts on one: the
process ends on the abort signal after the panic's line on stderr, with no exit code of its own.
A build that unwinds, as the test profile does, carries a panic of a worker thread to the
calling thread and ends with exit code 101.

## Threads

The comparison runs across the host's available cores. `VEP_SV_COMPARE_THREADS=<n>` bounds the
threads at work, the reading of a cache excepted (below), `1` keeping the whole comparison on
the calling thread; unset or empty, the bound is the number of cores the host makes available.
It is an environment variable and not an option, so that the command line stays the nine
options above, and a value that is not a whole number of 1 or more is refused. Every thread
count writes the same bytes to every report file, to stdout and to stderr and ends with the
same exit code: each fan-out returns its results in the order of its indices and reports the
failure of the lowest index that fails.

Across threads:

- the shards of the cache, read on up to the bound's threads while the calling thread files
  their records in the order of the shards, at most 256 shards ahead of the filing: under a
  bound of 2 or more, up to one thread more than the bound is at work;
- the input files and the two outputs of each, read as streams in runs of whole lines, a run
  split and normalised on whichever thread is free, and every column text held once in a table
  all the files share;
- the tuple sets and the raw counts of each file, and the unions of the two sides;
- the masks of the union of all files and of every file, side by side, and the masks of one pair
  of tuple sets side by side;
- inside the skipped-record mask, the gVCF-site mask and the count of mates past the 32-bit
  bound, the input records in runs of 1,024;
- the raw outputs behind the row-aware pair rules and behind the count of bare-breakpoint rows,
  one output a thread.

The fan-outs nested in one another (the files, the masks of a file, the records a mask walks,
the outputs it reads) share one crew, so at most the bound's threads are at work at a time. A
thread starts when work waits for one: a small input on a host of many cores uses fewer threads
than the bound.

On one thread: the filing of the cache's records, the ranking of the table of texts, the pooled
and the per-type counts, the lines of stdout, which are written in the order of the files, and
the four report files. One output's rows behind the pair rules are read by one thread, and the
input records such a rule asks for are read one file after the other.

## Layout

- `src/main.rs`, `src/cli.rs`, `src/help.txt`: the binary; the nine options with their
  defaults, the help and its text, the usage errors and the exit codes.
- `src/cli/streams.rs`: stdout and stderr as CPython holds and ends them: the blocks stdout is
  written in, and the exit code of a run whose stream does not take what it wrote.
- `src/compare.rs`: the comparison itself: the input files found and read, the raw and adjusted
  counts, the lines of stdout and the refusals.
- `src/text.rs`: the reading of a text file (universal newlines, strict UTF-8, runs of whole
  lines) and the path and number spellings the reports and messages depend on.
- `src/threads.rs`: the thread bound, the crew that nested fan-outs share, the fan-outs in
  index order and the reader that hands runs of lines to free threads.
- `src/vcf.rs`: the input VCF, plain or gzipped: its records, each with its variant type.
- `src/output.rs`: the reader of an engine's output: its tuples, its Locations and the variant
  type of each printed Location, over one table of column texts.
- `src/normalize.rs`: the Location, the Allele and the consequence set as they are compared,
  and the class names of fastVEP's symbolic tokens.
- `src/strings.rs`: the tuple as five ranks in the sorted table of texts, the sets of tuples
  and the two sides of a comparison.
- `src/metrics.rs`: the counts of a file, the pooled row, the per-type rows and the one-sided
  tuples behind them.
- `src/perl.rs`: Ensembl VEP's readings of a text: truth, a number, a POS column, the INFO
  column.
- `src/record.rs`: what both engines make of an input record: its skip route, its Location and
  the alleles each engine prints for it.
- `src/locus.rs`: the span and the chromosome of a Location and the mate of a breakend allele.
- `src/cache.rs`, `src/mate_context.rs`: the cache as three views by chromosome, and the terms a
  point breakend's mate-side row supports, derived from a transcript's structure.
- `src/registry.rs`: what this comparator adds to the registry of divergence classes: the
  verdict of each one-sided class and the names the masks report under.
- `src/masks/mod.rs`: the masks of each release in their order, their result and the function
  that runs them.
- `src/masks/release115.rs`: the masks of the 115.2 reference alone.
- `src/masks/cnv_tr.rs`: the tandem-repeat classes: the literal and the symbolic reading, the
  size of the class, the empty expansion and the dropped allele.
- `src/masks/pairs.rs`: the registry's pair rules on this comparator's pairs, with the rows and
  the input records the row-aware rules read.
- `src/masks/mate_read.rs`: the mate-side pair the local read explains.
- `src/masks/skipped.rs`: the records the 116.2 reference dropped, their rows, the cross-check
  and the gVCF variant site.
- `src/masks/one_sided.rs`: the 116.2 classes read off a record's own text: the dot breakend,
  the missing ALT, the star-only ALT, the non-integer coordinate and the numified POS.
- `src/masks/breakend116.rs`: the breakend classes of the 116.2 reference and the two counted
  shapes.
- `src/report.rs`: the JSON emitter, the markdown, the two tables and the release block.
- `src/failure.rs`, `src/tempdir.rs`, `src/lib.rs`: the one error type, the scratch directory
  of the tests and the crate root.
- `../vep-outputs-compare/src/pyfmt.rs`, `../vep-outputs-compare/src/input.rs` and
  `../vep-outputs-compare/src/rules/`: two source files and a directory of source files of the
  tuple comparator, compiled into this crate as `pyfmt`, `input` and `rules`, so that the CPython
  formatting rules, the splitting of a file into records, the writing of a Location's
  coordinates as numbers, the normalising of a consequence set and the registry of divergence
  classes with its pair rules have one definition for both comparators.
- `../vep-formats-compare/src/rules/coding_end.rs` and
  `../vep-formats-compare/src/genetic_codes.rs`: two source files of the per-column comparator,
  compiled into this crate as `coding_end` and `genetic_codes`, the names by which the coding
  end rules under `rules` ask the shapes of `alternate_codon_past_the_coding_end`.
- `tests/reports_pairs.rs`: the four report files of hand-built comparisons, byte for byte, and
  the per-type table the binary writes for a record of several identifiers.
- `tests/help.rs`: the binary's help and usage errors: stdout, stderr and the exit code.
- `tests/streams.rs`: the binary on a stdout or stderr that takes nothing and on one file for
  both: the exit code, what the other stream carries, whether the reports are written and where
  the lines of stderr stand among the text of stdout.
- `tests/cli.rs` and `tests/fixtures/`: the binary end to end. The tests write directories of
  inputs, outputs and caches, run the binary and read its reports, stderr and exit code. Two
  committed fixtures carry goldens: `reference_release`, one input whose reference output exists
  in the 115.2 shape (`perl/`) and in the 116.2 shape (`perl_116/`), and `release_classes`, four
  inputs with a record for every class of the 116.2 reference the first does not reach. Each
  `expected_*` directory holds the four report files of one run, `stdout.txt` with the output
  directory written `{{OUTPUT_DIR}}`, `exit_code.txt` and, for a run that is refused,
  `stderr.txt` and no report; `cargo test` runs the binary over the fixtures on one thread, on
  two and on the host's default and compares bytes. On the comparisons of the 116.2 goldens a
  run naming `116.3` is held to the run naming `116.2`, the release each states apart, and a run
  naming no release to the run naming `116.3`, byte for byte.
- `tests/splice_term_pairs.rs`: the binary on made-up pairs of rows of
  `splice_term_past_the_record_span` and of its nearest shapes, in three readings of a pair's
  record: the pairs the class takes and the counts the report gives them.
- `tests/tandem_repeat_loss_read_without_a_fasta.rs`: the binary on two records of the corpus
  `tests/golden/116/GRCh38-vcf-tandem-repeats`, a gain and a loss of one unit, the reference's
  rows from the run that read a FASTA and the scored side's from the run that read none: both
  readings of each record are set aside.
- `tests/tandem_repeat_tract_term_of_the_run.rs`: the binary on one record of the GRCh37
  synthetic tandem-repeat set whose run overlaps the polypyrimidine tract of an intron and whose
  deleted units do not: the tract term of the symbolic reading is set aside beside a cache that
  shows it and stays charged without one.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-sv-compare/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

CI runs the test, clippy and fmt commands above in its Test, Clippy and Format jobs, holds the
doc line of every source file in Module Docs, builds the package on the declared `rust-version`
in MSRV and audits its lockfile.

The goldens are regenerated from the repository root with the release build. Each call of
`golden` writes one `expected_*` directory; its arguments are the fixture, the directory, the
release, the reference's directory and any further options.

```bash
bin=$PWD/tools/vep-sv-compare/target/release/vep-sv-compare
golden() {
  f=$PWD/tools/vep-sv-compare/tests/fixtures/$1
  e=$f/$2 release=$3 perl=$4 work=$(mktemp -d)
  shift 4
  rm -rf "$e" && mkdir "$e"
  "$bin" --reference-release "$release" --assembly grch37 --input-dir "$f/in" \
    --perl-dir "$f/$perl" --rust-dir "$f/rust" --vep-rs-cache "$f/cache" \
    --output-dir "$work/out" "$@" > "$work/stdout" 2> "$work/stderr"
  echo $? > "$e/exit_code.txt"
  sed "s|$work/out|{{OUTPUT_DIR}}|g" "$work/stdout" > "$e/stdout.txt"
  [ ! -s "$work/stderr" ] || sed "s|$work/out|{{OUTPUT_DIR}}|g" "$work/stderr" > "$e/stderr.txt"
  for report in concordance_report.json concordance_report.md discordant.tsv discordant_adjusted.tsv; do
    [ ! -f "$work/out/$report" ] || cp "$work/out/$report" "$e/$report"
  done
  rm -rf "$work"
}
golden reference_release expected_115_2 115.2 perl
golden reference_release expected_116_2 116.2 perl_116
golden reference_release expected_116_2_refused 116.2 perl
golden release_classes expected_115_2 115.2 perl
golden release_classes expected_116_2 116.2 perl
golden release_classes expected_116_2_no_size_cap 116.2 perl --reference-max-sv-size -1
golden release_classes expected_116_2_other_engine 116.2 perl \
  --scored-engine fastvep --reference-max-sv-size -1
```
