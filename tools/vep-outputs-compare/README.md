# vep-outputs-compare

The tuple comparator of the concordance harness. It compares the consequence tuples of two
directories of VEP default-format outputs, the reference's and the scored engine's, in three
modes, and writes the report files (`summary.json`, `summary.md`, `discordant.tsv`,
`discordant_open.tsv`, the per-class and pooled CSVs, the per-term totals cache) and its stdout
lines. The comparison is the longest step of a whole-genome concordance run after the annotation
itself, so it streams its inputs and runs across the host's cores.

```text
vep-outputs-compare --perl-dir DIR --rust-dir DIR --report-dir DIR [options]
vep-outputs-compare --by-consequence-class --perl-dir DIR --discordant-in TSV \
    --dataset NAME --class-csv CSV --sweep-id ID --binary-md5 DIGEST [options]
vep-outputs-compare --pool-consequence-class --class-csv CSV --pooled-csv CSV \
    [--per-term TSV] [options]
```

Every option is accepted as `--flag=value` and by an unambiguous prefix as well. The others,
each with a default: `--glob`, `--assembly`, `--canonical-contigs` and `--no-canonical-contigs`,
`--json-cache`, `--reference-fasta`, `--ensembl-cache`, `--reference-release {115.2,116.2,116.3}`,
`--input-dir` and `--input-vcf`, `--summary-json`,
`--summary-md`, `--discordant-tsv`, `--fail-below-f1`, `--engine {vep-rs,fastvep}`,
`--assembly-label`, `--append`, `--min-tuples`, `--vep-totals-cache`, `--per-term` and
`--per-term-engine`.

## Contract

The bytes of every report follow the CPython rules this section names, which are part of the
contract. Records end at `\n`, `\r\n` or a lone `\r`, as universal newlines read them, and are
decoded as UTF-8 with `errors="replace"`. Fields are stripped with `str.strip()`'s whitespace set.
Every path option goes through `str(pathlib.Path(arg))`, and the directories the summary names
through `Path.resolve()`. The JSON is `json.dumps(indent=2, ensure_ascii=True)` with the insertion
order of every object, and its precision, recall and F1 values are `float.__repr__`.

The markdown and stdout lines format with `f"{x:.6f}"`; the per-class CSVs round with
`Decimal.quantize` under `ROUND_HALF_UP` and are written as `csv.DictWriter` writes them, `\r\n`
rows and minimal quoting. The top consequence sets keep `Counter.most_common`'s tie order.

The intended-divergence rules are the six pair rules in classification order, the sixteen row-aware
pair rules of the 116.2 registry after them, the 115.2 transcript-selection class over records
Ensembl VEP marks `vep_skip`, and the fourteen 116.2 one-sided classes. The summary of a run
against the 116.2 reference lists its classes in `divergence_classes`, each with its definition,
its taxonomy class, its citation and whether the adjusted F1 sets it aside. A row-aware pair rule
reads, beside the pair's two consequence sets and its allele, the first column, the Location,
the cDNA_position and the CDS_position of the rows that carry the pair's two tuples, and the
Protein_position, Amino_acids and Codons of a row that has them, fetched from the two outputs
for the pairs no other rule takes; the eleventh,
`stop_codon_settled_before_unknown_bases`, the twelfth, `frameshift_inside_a_held_stop_codon`,
and the fourteenth, `pachysolen_nuclear_code_read_as_standard`, read their Amino_acids and Codons,
the eleventh and twelfth columns, and take no pair with a row that ends before them. The eleventh
rule holds a pair
to the stop codon its Codons cell settles: the cell is the same on both rows, a stop codon of the
standard genetic code stands among the whole codons of its alternate half before the first codon
with a letter other than A, C, G and T, the reference's row writes no residues, and the alternate
half of vep-rs's residues holds its first stop at that codon's place. The twelfth reads the terms
and those two cells and no other column: vep-rs's row holds `frameshift_variant` where the
reference's holds `coding_sequence_variant` alone, `inframe_deletion` alone, or
`inframe_insertion` with `coding_sequence_variant`, both rows hold one Amino_acids cell, whose
half before the `/` opens with `*` and whose half after it holds `X`, and one Codons cell, whose
two halves differ by a number of bases that is no multiple of three. The fifteenth,
`start_codon_removed_without_five_prime_utr`, reads the Amino_acids alone of the two and, as they
do, takes no pair with a row that ends before the twelfth column. It holds a pair to a removal from
the first coding
base of a transcript with no 5' UTR: the cDNA_position and the CDS_position each open at 1 and end
on one number on both rows, both rows hold `inframe_deletion`, the Amino_acids of neither writes
an alternate residue, and the terms differ by the `start_lost` vep-rs's row alone holds. The
sixteenth, `start_codon_kept_without_five_prime_utr`, reads the two position columns alone and
takes a pair on rows that end at the ninth column. It holds a pair to a span on the start codon
of a
transcript with no 5' UTR: the cDNA_position and the CDS_position write one span of one number or
two on both rows, read up to a `/total` suffix, whose first number is 1, 2 or 3; a span with an
end that maps to no transcript base is none. The terms of vep-rs's row are those of the
reference's with `start_retained_variant` added and without the reference's `start_lost`, where
it holds one, and both rows hold `frameshift_variant`, `inframe_insertion` or
`inframe_deletion`. A pair whose allele opens with
its whole
reference allele is put to the row-aware rules also where a set rule
classifies it without setting it aside, and is counted under the row-aware rule that takes it. A
row among them whose first column spells no alleles takes the reference bases its allele stands in
for from the record of the input VCF behind its output file (`--input-vcf`, or the file of
`--input-dir` with the output's stem), which is read once more where such a row asks; these
rules read no other program. Where two or more records read such a row differently, the record
whose identifier is the row's first column decides; where no record carries that identifier, or
more than one does and they read the row differently, the row takes none. Of a `--json-cache`
directory these rules read the `info.json` for the names of the run's sequences (Sequence names,
below) and the transcripts of the sequences on which a pair asks for the coding ends of its
transcript (Coding ends, below). Three classes outside the pair rules read a cache as well:
`merged_same_dbid_collapse` and `merged_same_stable_id_collapse` the region files of that
directory (the paragraphs that follow), and `synonym_of_several_sequences_by_hash_order` the
sequence directories and the `chr_synonyms.txt` of an `--ensembl-cache` directory (the section
after Sequence names).

A comparison against the 116.2 or the 116.3 reference that names the converted cache both
engines read, with `--json-cache DIR`, also applies the classes `merged_same_dbid_collapse` and
`merged_same_stable_id_collapse`. They take their rows after the pair rules and before the
one-sided classes, the first before the second, and the reports list them in that order after
every other class in `excluded_categories` and `divergence_classes` but
`synonym_of_several_sequences_by_hash_order`, which stands last in a run that lists all three;
without the option neither class is applied and no report names one. Under `--merged` the
reference keeps one transcript per `dbID` among those an input buffer loads, and under `--refseq`
or `--merged` it then keeps one per stable id among those left, the `ensembl`-sourced one or
else the one of the lowest `dbID`. The rows alone cannot tell a row lost either way from a row
vep-rs writes too many, so the rules replay the two passes (`src/dbid_replay.rs`). The replay
reads the buffers from the reference's own output, in which every record the reference annotated
has a row in input order and whose header gives `--buffer_size` and `--fork`, and the transcripts
of each buffer from the cache's region files in the order the reference reads them.

The first class sets aside vep-rs's rows of a transcript the first pass dropped for a transcript
of the other set, at the records within 5,000 bases of it where the reference wrote no row of
that transcript. The second judges together vep-rs's rows of one stable id at one Location and
allele, since two transcripts of an id share the name their rows carry. It sets them aside where
the second pass dropped transcripts of the id within 5,000 bases of the record, at least as many
as the rows, the reference wrote no row of the id there that vep-rs lacks and, where it wrote
one, the transcript the pass kept lies within 5,000 bases of the record too. Either class sets
aside the reference's `intergenic_variant` row where every row vep-rs wrote at its Location and
allele is set aside so, the second where one of those rows is its own. For the rows the replay
names, the reference's output is read once more at their Locations, since a row both engines
wrote is not among the discordant rows.

A reference whose header holds no command line, or one that names neither `--refseq` nor
`--merged`, is no run of either class: none of its rows is read for them and nothing is written
to standard error. A run under `--refseq` is one of the second class alone. Under the 115.2
reference neither class is applied, and the cache gives the sequence names and, to
`start_cooccurrence_swap`, the coding start of the pair's transcript (Coding ends, below).

The replay follows such a run only where it can follow the reference line for line. For a file
of another kind it sets nothing aside and writes to standard error one line with the reason for
each class the run concerns:

- an option on the reference's command line outside `FOLLOWED_OPTIONS` of `src/dbid_replay.rs`,
  among them `--pick`, `--distance` and a `--fork` count above one, whose chunks depend on which
  child reports back first;
- an allele that is no run of bases, a record over a span the reference does not annotate, or
  one at a coordinate past which it writes a whole buffer as intergenic;
- no input VCF behind the output (`--input-vcf` or `--input-dir`), or one that holds a record
  the reference wrote no row for: such a record still sits in its buffer and asks for its
  regions, so the rows give the buffers only where the input holds no record beyond theirs;
- an empty line of the input VCF before a record, at which the reference ends a buffer;
- a region the cache holds no file for that one process keeps in memory from the buffer before,
  under which the reference holds the transcripts of the file it read next;
- a sequence name of which the cache files two synonyms.

Where no transcript of a stable id is `ensembl`-sourced and one of them carries a `dbID` that is
no plain number, the reference orders them by Perl's reading of that text: the replay names no
transcript of that id dropped and writes no line.

A region file the cache holds that is not the converter's JSON ends the run with an error.

A divergence class that reaches the whole-genome tuples is added in one change to the rule
tables (`EXCLUSION_REGISTRY` in `src/rules/mod.rs` for a pair rule, `DIVERGENCE_CLASSES_116_2` in
`src/rules/table.rs` for a class of the 116.2 registry, `MERGED_SAME_DBID_COLLAPSE` and
`SYNONYM_OF_SEVERAL_SEQUENCES_BY_HASH_ORDER` beside it for a class that reads a cache) and to
`data/divergence_classes.tsv`,
which states the classes for the programs that read them (Class table, below). The test of
`src/rules/table.rs` reads that file and fails while the two state a class, its place, a flag,
a citation or a definition differently.

A ratio over no tuple is 0, so the adjusted precision, recall and F1, taken over the tuples the
rules leave, read 0 where every tuple is set aside or neither output holds one.

The keys of each side are compared in bytewise order, the order `sort -u` gives them under
`LC_ALL=C`. On one thread a side's keys are sorted in memory, spilled as sorted runs under the
report directory past one GiB of keys and merged on finish. On more threads the keys of both
sides are routed by Location into partitions that cover that order end to end, written under the
report directory, and each partition is sorted and merged on its own; read in partition order the
keys, the counts and the discordant rows are those of the one sort. The merge holds every key
against the one before it on its side, and the first key of each partition against the last key
before it, and at the first that does not sort after it the run ends with exit code 1, writes no
summary and keeps the key files it stopped on in their directory under the report directory:
past such a key the merge would count a key both sides hold as missing and as extra at once. Its
error names that directory, the side's output file and the lines of that file that first write
the key and the key before it, counted from 1 over every line of the file; it names no key and no
partition, so it reads the same at any thread count. A run that stops while it compares the
paired files, on such a key or on a file it cannot read, leaves the discordant rows written so
far under the `--discordant-tsv` name with a `.partial` suffix, which its error names, and no
file under that name itself. An output or an input VCF that cannot be read again for the
row-aware pair rules ends the run with an error that names that file, as an error of the pass
over `discordant.tsv` names the report it is about. The `--glob` pattern names
entries of the output directory itself (`*`, `?` and `[...]`); a pattern with a `/` is refused.
The one-sided classes of the 116.2 reference read a Location's coordinates and a breakend mate's
position, written in the digits `0` to `9`, whole at any length; every other coordinate and
position is read as an integer up to 128 bits. The row-aware pair rules read within that bound
the two ends of a row's Location, its two cDNA positions and the position a generated record
name spells, and the input records are asked for a row's reference bases by the start of its
Location, read within the bound as well. A number past 2^127 - 1 reads there as no number: a
rule that needs the number refuses the pair, no record is asked for the row, and the pair's two
rows stay in `discordant_open.tsv`. Exit
codes: 2 on a usage error or an aggregate F1 below `--fail-below-f1`, 1 on a missing mode argument
(`error: mode <mode> requires ...` on stderr), on a `--canonical-contigs` value that names no
contig (`error: --canonical-contigs '<value>' names no contig` on stderr) or on any other failure
(`ERROR: [vep-outputs-compare] <detail>` on stderr), 0 otherwise.

## Sequence names

Two rules read the names of the run's sequences. The contig filter keys the rows whose Location
lies on a sequence of its set and no other. The class `bnd_synonym_mate_rows_lost` tells a
breakend mate written under the name the cache files its sequence under from a mate written under
a synonym.

- Under `--assembly GRCh37` or `--assembly GRCh38` both read the tool's own human names: `1` to
  `22`, `X`, `Y`, `MT` and `M` for the filter, the same without `M` for the cache. The bare
  `--canonical-contigs`, the default, is that filter, and the `info.json` of a `--json-cache`
  directory is not read for them.
- Under any other assembly, or none, the names come from the run. `--json-cache DIR` names a
  converted cache, whose `info.json` lists its sequences under `contig_order`: the bare filter
  takes that list, and so does the class. `--canonical-contigs <file or list>` names a set (a
  file of names separated by commas or white space with `#` lines ignored, or a comma-separated
  list), which the class takes where the run names no cache.
- With neither, the bare filter is off, which a `WARNING:` line on stdout says, and the class
  makes no test. The rows that reached the test are counted as `not_evaluated` in the class's
  statistics of `summary.json` and in a `Not evaluated` line of `summary.md`; a report with no
  such row carries neither.
- `--canonical-contigs <file or list>` names the filter's set under every assembly, and
  `--no-canonical-contigs` turns the filter off without the warning.

The summary's `canonical_contigs_filter_active` and `canonical_contigs_set` state the filter that
ran. Under an assembly that is not a human one, a `--json-cache` directory without an `info.json`
that lists at least one name ends the run with exit code 1; under every assembly a path that is
no directory is a usage error.

## A name the synonyms give to several sequences

A cache's `chr_synonyms.txt` can give one name to two or more of the sequences the cache holds:
Ensembl VEP's dump links names with one level of indirection (`BaseVEP.pm:596-621`), so a synonym
several sequences share links their other names to each other. The file then does not say which
sequence such a name stands for. Ensembl VEP reads a record under it on whichever of them comes
first in the order `keys` returns the name's synonyms (`BaseVEP.pm:711-716`), which changes
between runs of one command; vep-rs skips the record with a warning. That is the documented
difference `synonym_of_several_sequences_by_hash_order`.

- `--ensembl-cache DIR` names the Ensembl VEP cache directory the reference run read, the one
  that holds a directory a sequence. The tool takes the held sequences from those directories,
  as Ensembl VEP lists them (`CacheDir.pm:431`), and reads the `chr_synonyms.txt` beside them as
  both engines read the file: the first field of a line and each later field name each other,
  one hop. A name is of the class where the cache holds no sequence by it and two or more of
  the names it is paired with are held. The names come from the reference's own cache, never
  from vep-rs's converted one, so rows vep-rs lost on a sequence its conversion dropped are no
  rows of the class.
- Against the 116.2 or the 116.3 reference the reference's rows at a Location under such a name
  are set aside where vep-rs wrote no row at that Location: none among the discordant rows and
  none in vep-rs's output, which is read once more for those Locations. A record vep-rs wrote a
  row for stays open.
- Only a reference run that read the cache's own synonyms file is of the class: the header of
  the reference output must state a command line (`## VEP command-line:`) that names no
  `--synonyms`. A run given that flag, and an output that states no command line, keep their
  rows open. The rule reads no warnings file and no flag of vep-rs's run.
- The class stands in the registry of the classes that read a cache, after
  `merged_same_dbid_collapse` and `merged_same_stable_id_collapse` (`cache_116.2` of
  `data/divergence_classes.tsv`, the
  constant `SYNONYM_OF_SEVERAL_SEQUENCES_BY_HASH_ORDER` of `src/rules/table.rs`, the rule in
  `src/rules/shared_synonym.rs`), outside the one-sided classes. A run whose cache gives such a
  name lists the class last in `excluded_categories` and `divergence_classes`; every other run,
  one that names no such cache or a cache without such a name among them, writes its reports as
  before. A path that is no directory is a usage error.

## Coding ends

Ten of the sixteen row-aware pair rules ask what stands at the two ends of the coding sequence of
the pair's transcript, the thirteenth, `splice_term_past_the_record_span`, asks it for the rest of
a difference it shares with one of them, and the eleventh,
`stop_codon_settled_before_unknown_bases`, the twelfth, `frameshift_inside_a_held_stop_codon`,
the fourteenth, `pachysolen_nuclear_code_read_as_standard`, the fifteenth,
`start_codon_removed_without_five_prime_utr`, and the sixteenth,
`start_codon_kept_without_five_prime_utr`, ask nothing of a cache. vep-rs writes a stop codon's
terms and a start codon's terms only where that codon stands, so the eleven that ask read the
candidate's cache, `--json-cache DIR`, under every assembly.

- The five classes of a coding end (`structural_deletion_stop_lost_without_stop_codon`,
  `stop_codon_completed_by_polyadenylation`, `stop_codon_overlap_without_stop_codon`,
  `stop_codon_behind_last_residue_without_stop_codon`, `start_codon_overlap_without_start_codon`)
  take a pair only where vep-rs's coding terms are the ones derived for the row. The derivation
  evaluates Ensembl VEP's own term predicates with the stop codon and start codon tests answered
  by the transcript's cache record. It reads the five cells from cDNA_position to Codons, which
  both rows must print equal but for an Amino_acids and a Codons that differ as
  `alternate_codon_past_the_coding_end` documents, where it reads the reference's, and on a stop
  codon that polyadenylation completes the cached coding sequence.
- The five stop codon classes and `start_cooccurrence_swap` stood before the cache was read.
  Beside a cache each takes a pair only on a transcript whose record shows the codon the class
  stands on, at the allele's place.

The record is looked up by the pair's Feature, a version after the stable id dropped where the
cache's id has none, on the sequence its Location names. That sequence is the directory under
`transcripts/` of the Location's own name, of the name with or without a leading `chr`, or of
the first name the cache's `chr_synonyms.txt` pairs with one of those spellings that the cache
holds, as vep-rs finds a record's sequence in a cache; `M` and `MT` are spellings of each
other. A row whose first column spells a record with an empty reference allele is an
insertion's, whatever an input record reads at its Location. Where one stable id stands on several
records of a sequence, the record is the one whose span shares a base with the Location; a pair
that no record or more than one answers has no record. The transcripts of a sequence are read
when the first pair on it asks, every shard of its directory under `transcripts/` on the run's
threads, and no other sequence is read: a run whose pairs lie on every sequence of a cache reads
the whole of its transcripts once. A record is kept as the
reading of its two ends, its span and its coding length; the model's coordinates are kept for a
transcript with an end that shows no codon, and the coding sequence for one whose stop codon
polyadenylation completes. `tests/fixtures/coding_ends/verdicts.tsv` at the repository root
holds the reading to the one vep-rs's loader makes on the committed corpus caches.

What a comparison reads decides what it sets aside:

- Without `--json-cache`, no class of a coding end takes a pair, and the classes that stood
  before take a pair by its two rows alone, on any transcript. A row vep-rs writes with a stop
  term or a start term where no such codon stands can then be set aside under one of them.
- Beside `--json-cache`, a pair whose transcript the cache answers with no record, or with a
  record vep-rs does not read, is taken by none of the eleven rules.
- Without the input records (`--input-vcf` or `--input-dir`), a row of a sequence allele under
  an identifier carries no reference bases, and the four classes that read a sequence allele
  leave it a difference. The structural class reads no record.
- Without `--reference-fasta FILE`, a deletion that begins in a stop codon polyadenylation
  completes and runs past the transcript's last base is read at the record's own place, as
  vep-rs reads it in a run without a reference sequence. Beside the file, the bases after the
  transcript tell where else the same deletion can be written, and whether the stop codon is
  kept is read there. The file is a FASTA read
  through its `.fai` index where it has one and is not gzipped, and held whole otherwise.
  Beside a file that does not name the row's sequence no base is read, and such a deletion is
  read as it is without the file.

A `--json-cache` directory without a `transcripts/` directory ends the run with exit code 1 when
a pair first asks for a record, and so does one with a JSON array or object under a shard's name
that is neither an array of transcript objects nor an object of such arrays. A file under a
shard's name that cannot be read whole, is not UTF-8 or is no JSON array or object is passed
over. A `--reference-fasta` path that is no file is a usage error, and an index that states a
sequence longer than its FASTA holds ends the run with an error that names the index.

## Threads

The comparison runs across the host's available cores. `VEP_OUTPUTS_COMPARE_THREADS=<n>` bounds
the threads at work, `1` keeping the whole comparison on the calling thread; it is an environment
variable and not an option, since the command line names the comparison and the thread count
changes no byte of it, and a value that is not a whole number of 1 or more is refused. Every
thread count writes the same bytes to every report file and to stdout.

One thread reads an output file in runs of whole records while the others key them, and every
thread then sorts and merges partitions. The pass over `discordant.tsv` that applies the
intended-divergence rules, the read of an input VCF, the read of the record names behind the
skipped-record class, the read of the rows behind the row-aware pair rules, the replay over the
transcripts of a converted cache and the `--pool-consequence-class` mode stay on one thread; the
region files of that cache are read across the threads.

## Class table

`data/divergence_classes.tsv` states every divergence class the rules hold, one line a class
under a header line of column names: tab-separated, a field the bytes between two tabs with no
quoting, no comment lines. The lines of one reference release, read top to bottom, are that
release's classes in its registry's order: the pair rules, the row-aware pair rules, the
one-sided classes, then the classes that read a cache. Its columns:

- `registry` and `order`: the rule list the class stands in and its place there, counted from 1.
  `pair` holds the six rules over the two consequence sets of a pair and `row_pair_116.2` the
  sixteen that read the pair's two output rows, each list in classification order, the first rule
  that
  matches taking the pair; `one_sided_115.2` and `one_sided_116.2` hold the classes over rows one
  engine wrote and the other did not, and `cache_116.2` the three classes over such rows that a
  comparison applies only where it names a cache: the two classes of the converted cache
  (`--json-cache`), then the class told from the files of an `--ensembl-cache` directory.
- `bucket`: the name the reports count the class under.
- `reference_releases`: the `--reference-release` values whose registry holds the class,
  comma-separated.
- `reads`: what the rule is asked on: `sets` (the two consequence sets of a pair),
  `sets_and_allele` (those and the allele of the pair's tuples), `rows` (the pair's two output
  rows), `rows_and_coding_ends` (those rows and the coding ends of the pair's transcript in the
  candidate's cache), `record` (the rows one engine wrote for a record) or `cache` (such rows,
  told by a replay of the reference's buffers over the converted cache or from the files of the
  reference's cache).
- `classifies`: `true` or `false`, whether a pair the rule matches is counted under the bucket;
  empty for a one-sided class and for a cache class, whose rules carry no such flag.
- `excludes`: `true` or `false`, whether the adjusted F1 sets the class's rows aside.
- `excludes_with_cache`: `true` or `false`, whether the structural-variant comparator, which reads
  the input records and a cache, sets them aside; empty for a pair rule and for a cache class,
  which that comparator does not hold.
- `taxonomy_class`, `perl_citation` and `definition`: the class's taxonomy class, the Ensembl VEP
  source lines it cites and its definition, in the words of `divergence_classes`.

## Layout

- `src/input.rs`: records, lossy decoding, the five-column tuple key with its Location
  normalised and the canonical-contig filter, the record names of the rows it keys, and the
  reader of an input VCF, plain or gzipped.
- `src/sort.rs`: the in-memory sort-unique with its spilled runs, and the sorted keys of one
  partition.
- `src/partition.rs`: the bounds between partitions sampled from an output file, the scan that
  routes a side's keys into partition files across threads, and the files themselves.
- `src/threads.rs`: the thread bound and the indexed fan-out over partitions.
- `src/merge.rs`: the streaming merge of two sorted key streams, the discordant rows and the
  per-set counters, feeding the selected release's one-sided collector.
- `src/dbid_replay.rs`: the replay of the reference's `dbID` pass and of its pass per stable id
  over a RefSeq or a merged cache: the run read from the reference's header, its buffers read
  from its rows, the reader of a converted region file and the transcripts each buffer dropped
  for one of the other set or for another of their stable id.
- `src/rules/`: the pair rules (`pair.rs`), what a row-aware pair rule reads of an output row
  and of its record's alleles (`row.rs`), the transcripts of the candidate's cache
  (`cache_transcripts.rs`), the reference sequence (`reference_bases.rs`) and what a rule of a
  coding end is asked with (`coding_end.rs`), the Location and breakend readers (`breakend.rs`),
  the two collectors (`scope.rs`), the gVCF variant-site recognisers (`gvcf.rs`), the 116.2
  one-sided pass (`one_sided.rs`), the two classes of a comparison that names the converted
  cache (`shared_dbid.rs`, `shared_stable_id.rs`), the class read from a cache's synonyms file
  (`shared_synonym.rs`), the class table (`table.rs`) and the pass over `discordant.tsv`
  (`filter.rs`).
- `../vep-formats-compare/src/coding_ends.rs` with `coding_end_terms.rs`,
  `../vep-formats-compare/src/rules/coding_end.rs` and
  `../vep-formats-compare/src/genetic_codes.rs`: source files of the per-column comparator,
  compiled into this crate as `rules::coding_ends`, `coding_end` and `genetic_codes`, so that
  the coding ends of a cache transcript, the terms derived on a row there and the shapes of
  `alternate_codon_past_the_coding_end` have one definition for the comparison tools.
- `src/classes.rs`: the per-consequence-class counts, their pooling, the floor, the totals cache
  and the two CSV schemas.
- `src/report.rs`: the JSON emitter, the summary markdown and the per-file metrics.
- `src/compare.rs`, `src/per_class.rs`, `src/main.rs`: the three modes and the command line.
- `src/pyfmt.rs` and `src/csv.rs`: the CPython formatting and `csv` module rules listed above.
- `src/*_tests.rs` and `src/rules/*_tests.rs`: the unit tests, each file a test module of the
  source file its name opens with:
  - `src/classes_per_class_counts_tests.rs`: the per-class counts of a key file and of a
    discordance report, their pooling, the floor, the two CSVs and the totals cache.
  - `src/compare_canonical_contigs_tests.rs`: the contig set a `--canonical-contigs` value names.
  - `src/compare_one_sided_readers_tests.rs`: what the one-sided classes of the 116.2 reference
    read of an output pair beside its keys.
  - `src/compare_stop_tests.rs`: how a comparison that stops leaves the report directory, and
    the lines of an output file its error names.
  - `src/dbid_replay_tests.rs`: the run a header states, the buffers and fork chunks, the
    reader of a region file, the filter at load, and replays over the converted cache of the
    reproduction `merged_same_dbid_collapse` and over caches written for a test.
  - `src/dbid_replay_stable_id_tests.rs`: the pass per stable id: the transcript kept of those
    that share an id, what each record sees of the ones dropped, and replays over caches written
    for a test and over the merged cache of the corpus `GRCh38-merged`.
  - `src/input_location_numbers_tests.rs` and `src/input_ten_digits_tests.rs`: a Location keyed
    by the numbers of its coordinates, read from the digits `0` to `9` alone.
  - `src/main_mode_arguments_tests.rs` and `src/main_reference_release_tests.rs`: the mode an
    argument list selects, the reference release option and the two options that name an
    input VCF.
  - `src/merge_collectors_tests.rs` and `src/merge_key_order_tests.rs`: what each release's
    collector holds after a merge, and the stop on a key below the key before it.
  - `src/pyfmt_half_up_rounding_tests.rs`: `py_round_half_up` on the values a reader of a
    per-class CSV rounds.
  - `src/report_ratio_tests.rs`: precision, recall and F1 where a side holds no tuple.
  - `src/sort_collation_tests.rs`: the sorts order keys by byte.
  - `src/rules/breakend_coordinate_readers_tests.rs` and
    `src/rules/breakend_ten_digits_tests.rs`: Locations, positions and mate names as the
    one-sided classes read them.
  - `src/rules/filter_exclusion_bounds_tests.rs`, `src/rules/filter_one_sided_passes_tests.rs`
    and `src/rules/filter_pair_rule_rows_tests.rs`: the pass over `discordant.tsv` under each
    reference release and with a fabricated row-aware pair rule.
  - `src/rules/gvcf_variant_sites_tests.rs`: the records that are gVCF variant sites and the
    Locations an output or an input VCF names them at.
  - `src/rules/mod_registry_tests.rs`: the invariants of the registry.
  - `src/rules/shared_synonym_tests.rs`: the names a synonyms file gives to two or more held
    sequences, the sequences a cache directory holds, and the cache class's pass with and
    without vep-rs's output directory.
  - `src/rules/pair_substitution_tests.rs` and `src/rules/pair_swap_shapes_tests.rs`: the swaps
    each pair rule takes and refuses.
  - `src/rules/row_pair_rule_rows_tests.rs`: the row-aware pair rules on the fabricated rows of
    `tests/fixtures/row_aware_pair_rules/tables/`.
  - `src/rules/row_number_bound_tests.rs`: the numbers a row-aware pair rule reads within 128
    bits, and the pairs that stay differences where a row writes one past them.
  - `src/rules/row_coding_end_tests.rs`: the rules that ask a transcript's coding ends on the
    made-up cache of `tests/fixtures/coding_ends/cache`.
  - `src/rules/cache_transcripts_tests.rs`: the reader of a cache's transcripts on made-up
    shards, and held to the committed table of coding ends on the corpus caches.
  - `src/rules/row_settled_stop_codon_tests.rs`: the stop codon a Codons cell settles, the rows
    the rule of `stop_codon_settled_before_unknown_bases` takes and refuses, its place in the
    registry, the pairs it leaves to the other row-aware pair rules and the pair of its own
    that an earlier rule takes first.
  - `src/rules/row_frameshift_inside_a_held_stop_codon_tests.rs`: the bases a Codons cell holds
    in its two halves, the rows the rule of `frameshift_inside_a_held_stop_codon` takes and
    refuses, its place in the registry and the pairs it leaves to the other row-aware pair rules.
  - `src/rules/row_removed_start_codon_tests.rs`: the rows the rule of
    `start_codon_removed_without_five_prime_utr` takes and refuses, its place in the registry,
    the columns a row is read to and the pairs it leaves to the other row-aware pair rules.
  - `src/rules/row_kept_start_codon_tests.rs`: the two forms of the term difference the rule of
    `start_codon_kept_without_five_prime_utr` takes on each kind of insertion's or deletion's
    row, the rows it refuses, its place in the registry, the columns a row is read to and the
    pairs it leaves to the other row-aware pair rules.
  - `src/rules/shared_dbid_tests.rs`: the class `merged_same_dbid_collapse` over fabricated
    discordant rows and a collector filled by hand, and what the class states of itself.
  - `src/rules/shared_stable_id_tests.rs`: the class `merged_same_stable_id_collapse` over
    fabricated discordant rows and a collector filled by hand, and what the class states of
    itself.
  - `src/rules/table_definitions_tests.rs` and `src/rules/table_reference_classes_tests.rs`: what
    the class table states of two row-aware pair rules and of the 116.2 reference's classes.
- `data/divergence_classes.tsv`: the class table, held to the rule tables by the test of
  `src/rules/table.rs`.
- `tests/golden.rs` and `tests/fixtures/`: synthetic output directories (`perl/`, `rust/`), an
  input VCF and the goldens of the three modes, the tool's contract on them. The goldens carry
  `{{PERL_DIR}}`, `{{RUST_DIR}}`, `{{INPUT_DIR}}`, `{{REPORT_DIR}}`, `{{SCRATCH}}` and `{{GENERATED_AT}}` for
  the values that depend on the checkout and the run; `cargo test` runs the binary over the
  fixtures on one thread, on two and on the host's default, substitutes them and compares bytes,
  and compares the reports of a generated pair of directories across thread counts. Under the
  same thread bounds it holds a stop on keys out of order to one error, a stop while the files
  are compared to the `.partial` report, an output that is gone when it is read again to its own
  name in the error, and a pair whose rows write a number past 128 bits to stay a difference.
- `tests/both_cache_classes.rs`: the binary over one output pair that holds rows of two of the
  classes that read a cache, in a comparison that names the converted cache, the Ensembl VEP
  cache, both or neither: each class's own statistics, and the order the reports state the three
  classes in.
- `tests/exclusion_bounds_and_per_class_counts.rs`: the binary over output pairs and saved
  reports written for each test: the open rows beside the discordant rows, the contig filter of
  an assembly that is not a human one from a named set and from a cache, a breakend mate told
  from a synonym by the names such a run gives, the totals cache of the per-class mode and the
  usage errors.
- `tests/merged_cache_shared_dbid.rs`: the binary over the rows both engines write for the
  three records of the reproduction `merged_same_dbid_collapse`, with and without that
  reproduction's converted cache
  (`scripts/divergences/reproductions/inputs/merged_same_dbid_collapse/json_cache`): the rows
  set aside, the rows that stay open and what the reports state of the class.
- `tests/pair_rules_and_key_order.rs`: the binary over one small output pair: a side's keys in
  byte order under any locale, and what a run that stops on keys out of order keeps.
- `tests/reference_registries_and_one_sided_classes.rs` and
  `tests/fixtures/reference_registries_and_one_sided_classes/`: the binary under each reference
  release over a committed output pair, the 115.2 report held to `expected_115_2/` byte for byte
  and the 116.2 report to the classes it applies and defines.
- `tests/coding_end_classes.rs` and `tests/fixtures/coding_ends/`: the binary over one output
  pair beside the made-up cache, without it, without the input records and beside a cache that
  holds no transcripts, and over a pair of a deletion past a stop codon polyadenylation
  completes beside the reference sequence (`--reference-fasta`).
- `tests/refseq_cache_shared_stable_id.rs`: the binary over the rows both engines write under
  `--refseq` for four records of a converted cache written for the test, in which four stable
  ids each name two transcripts: the rows set aside under `merged_same_stable_id_collapse`, the
  rows that stay open and what the reports state of the class.
- `tests/row_aware_pair_rules.rs` and `tests/fixtures/row_aware_pair_rules/`: the compare mode
  over fabricated output pairs for the row-aware pair rules, and its reports over the pairs of
  `reports/` held to the committed ones byte for byte; `tables/` holds the fabricated rows of
  the rules.
- `tests/synonym_of_several_sequences.rs`: the compare mode over fabricated output pairs beside
  a fabricated Ensembl VEP cache directory whose synonyms file gives one name to two sequences:
  the rows the cache class sets aside and the runs that leave them open.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-outputs-compare/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

CI runs the test, clippy and fmt commands above in its Test, Clippy and Format jobs, builds the
package on the declared `rust-version` in MSRV and audits its lockfile. Nothing regenerates the
goldens: the files under `tests/fixtures/` are the contract, `cargo test` holds the binary to
them, and a change of the contract changes them in the same commit.

## Harness

`scripts/concordance/run_clone_measurement.sh`, `scripts/concordance/run_concordance.sh` and
`scripts/validation/compare_cache_outputs.sh` score tuples with this tree's build of the tool,
and no environment variable selects another comparator. `scripts/genomes/check_genome.sh`
compares a genome's two outputs with it as well, and names the converted cache for every
genome. `run_clone_measurement.sh` and
`compare_cache_outputs.sh` run the release build,
`tools/vep-outputs-compare/target/release/vep-outputs-compare`, and build nothing;
`run_concordance.sh` builds it through cargo, as it builds the engine it runs.

A release build that is not an executable, or a cargo build that fails, stops the script before
its first engine run with one line that names the build command.
`run_clone_measurement.sh` asks for a comparator only where it scores F1: not under
`--walltime-only`, and not for the Ensembl VEP engine.
