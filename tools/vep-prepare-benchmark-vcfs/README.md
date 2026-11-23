# vep-prepare-benchmark-vcfs

The input preparation of the concordance harness. It reads every `.vcf` and `.vcf.gz` file of a
directory and writes each as plain VCF text that every engine of a run annotates alike: header
lines as they are, the first column of every variant line under one contig naming, and, under the
canonical-contig filter, the variant lines on the canonical contigs alone. It writes
`manifest.json` beside the prepared files, with the run's options and, for each file, its counts
and the SHA-256 of what was written, and one line a file on stdout. The prepared files, the
manifest and stdout are fixed to the byte.

```text
vep-prepare-benchmark-vcfs --input-dir DIR --output-dir DIR [--chrom-mode {auto,keep,strip,add}]
    [--cache-dir DIR] [--max-variants N] [--assembly NAME]
    [--canonical-contigs [FILE_OR_LIST]] [--no-canonical-contigs]
```

- `--input-dir`: the directory whose `.vcf` and `.vcf.gz` files are prepared.
- `--output-dir`: where the prepared files and `manifest.json` are written. It is created with its
  parents.
- `--chrom-mode`: the contig naming of the variant lines. `keep` leaves the first column as it is,
  `strip` takes a leading `chr` off it and `add` puts `chr` before a name that has none. `auto`,
  the default, is decided file by file: `strip` where the file's first variant line names a `chr`
  contig and the cache names plain ones, `add` the other way round, `keep` in every other case.
- `--cache-dir`: a converted JSON cache directory, read for its contig naming alone and only by
  `auto`.
- `--max-variants`: writes at most that many variant lines a file.
- `--assembly`: recorded in the manifest. It also decides what the bare `--canonical-contigs`
  means.
- `--canonical-contigs`, on by default: keeps the variant lines on canonical contigs and drops the
  others. Bare, under no `--assembly` or under `GRCh37` or `GRCh38`, the set is the human one:
  `1` to `22`, `X`, `Y`, `MT` and `M`. With a value the set is the names in the file the value
  names, or the comma-separated names of the value itself. Bare under any other assembly it names
  no set: the filter is off for the run and one `WARNING:` line on stdout says so.
- `--no-canonical-contigs`: switches the filter off.

The command line is read as CPython's `argparse` reads it. An option may be shortened to any prefix
that fits it alone (`--max 5`), takes its value as the next argument or after `=`, and a later
value replaces an earlier one. A value that starts with `-` is taken only when it is `-` alone,
starts like a negative number or holds a space. `--canonical-contigs` takes the next argument as
its value only where that argument is not an option. `--max-variants` takes what `int()` takes:
blanks around the number, a sign, single underscores between digits, the decimal digits of any
script as Unicode 16.0 lists them, up to 4,300 digits. `-h` and `--help` print the usage and the
options.

## Contract

**Inputs.** The inputs are the entries of `--input-dir` whose names end in `.vcf`, then those whose
names end in `.vcf.gz`, each group in code point order; no other entry is read, and names are
matched in the case they are written in. A name that ends in `.gz` is read as CPython's `gzip`
module reads a file: member after member, the header checksum of a member not checked, the zero
bytes that pad the file after a member skipped, a zero-byte file an empty input, and any other
byte where a member would start a corrupt stream. A prepared file carries its input's name less
the `.gz`. Two inputs that are written under one name are prepared one after the other in input
order, so the file holds the later one, and the manifest lists both with the digest each had when
it was written. An input's naming is read first; then the input is opened, its prepared file is
created and its lines are read: where `--input-dir` is `--output-dir`, a `.vcf` input is its own
prepared file and is cut to nothing before a line of it is read.

**Lines.** A line ends at `\n`, `\r\n` or a lone `\r`, as universal newlines read it, and is
decoded as UTF-8 with `errors="replace"`, so a byte outside UTF-8 is written as U+FFFD. A line that
starts with `#` is a header line, wherever it stands: it is written as it is, its line end as
`\n`, and a last line without a line end stays without one. Every other line is a variant line,
an empty line and a line that starts with a byte order mark among them. Its first column, the
bytes before its first tab, goes through the naming mode. Under the filter, a name that still
starts with `chr`, in any letter case, and is not in the set then loses those three characters,
and a line whose name is not in the set is dropped; a set that names a prefixed contig keeps that
spelling. A variant line that is kept is written with its first column as it came out and a
`\n`, whatever its line end was.

**The cap.** Under `--max-variants N` a file is read up to the first variant line that follows
the `N`th written one: header lines met until then are written, dropped lines count as dropped and
not towards `N`, and nothing behind that line is read. `N` at or below zero writes the header
lines that precede the first variant line.

**The naming of a file and of a cache.** A file's naming is that of the first variant line whose
first column is not empty: `chr` where it starts with `chr` in any letter case, `plain` otherwise,
`unknown` where the file holds no such line. A cache's naming is that of the first directory, in
code point order, under its `transcripts` directory, else under its `variations` directory, and
`unknown` where neither holds a directory or no `--cache-dir` is given.

**A named contig set.** The value of `--canonical-contigs` names a file where one exists at that
path; otherwise it is the list itself. A file is read as UTF-8 text, strictly, in the lines of
`str.splitlines()`; a line whose first character after its leading blanks is `#` is left out.
Names are separated by commas and by the blanks of `str.isspace()`. A value that names no contig
ends the run before anything is created.

**The manifest.** `manifest.json` holds, in this order, `generated_at_utc`, `input_dir`,
`output_dir`, `requested_chrom_mode`, `cache_dir`, `max_variants`, `assembly`,
`canonical_contigs_policy` (`universal_ensembl_style`, `named_set` or `disabled`),
`canonical_contigs_filter_active`, `canonical_contigs_set` (the names in code point order, or
`null`) and `files`. Each entry of `files` holds `input_file`, `output_file`,
`detected_input_chrom_style`, `detected_cache_chrom_style`, `normalization_mode`,
`headers_written`, `variants_written`, `variants_dropped_non_canonical`, `output_sha256` and
`chr_stripped`, which is true where the filter took a `chr` off at least one line of the file.
The three directories are written as `pathlib.Path.resolve()` names them.

**Stdout.** One line a file in input order,
`prepared <input> -> <prepared file> (mode=<mode>, variants=<n>, headers=<n>)`, with
`, canonical_dropped=<n>` after the header count under the filter and ` chr_stripped=true` after
that where it holds; then `wrote manifest: <path>`. The `WARNING:` line of a bare filter that names
no set comes first.

**Bytes.** Wherever a byte could be written in more than one way it follows a CPython rule, named
here by its CPython name. The manifest is `json.dumps(indent=2, ensure_ascii=True)` with no newline
after it. `generated_at_utc` is `datetime.isoformat()` of the UTC time the manifest is written at:
the date, `T`, the time, six decimals of a second where the microseconds are not zero, and
`+00:00`. Every
directory option is read through `str(pathlib.Path(arg))` and resolved as `os.path.realpath`
resolves it without `strict`: against the working directory, through every symbolic link, a
component that does not exist kept as it is named. `max_variants` is written as `str(int)` writes
the value, at any length. An argument or a file name that is not UTF-8 is read as
`surrogateescape` reads it: each such byte is `\udcXX` in the manifest and sorts as that code
point, and stdout carries the name in its own bytes.

**Exit codes.** 0 when every file, the manifest and stdout are written. 2 on a command line that
is refused: the usage, then `vep-prepare-benchmark-vcfs: error: <message>`, on stderr. 1 on a
`--canonical-contigs` value that names no contig, with the one line
`error: --canonical-contigs <value> names no contig` on stderr, the value as `repr()` writes it.
1 on any other failure: an input directory that holds no input (the output directory is created
first), an input that cannot be read or is a corrupt gzip stream, an output that cannot be written,
a contig list that is not UTF-8, a thread bound that is no count. Such a run writes one line,
`ERROR: [vep-prepare-benchmark-vcfs] <detail>`, on stderr and no manifest; stdout holds the lines
of the inputs before the first one that failed, which are prepared whole.

**Where the bytes are not CPython's.** `repr()` escapes every character Unicode does not class as
printable; the refusal lines escape the ASCII ones, U+0080 to U+00A0, U+00AD, the blanks of
`str.isspace()` and a byte outside UTF-8, and write any other character as it is. After a
failure, inputs behind the failed one that were already being prepared are finished, so the output
directory may hold their prepared files. What a failed input leaves of its own prepared file is
not fixed.

## Threads

The preparation runs across the host's available cores. `VEP_PREPARE_BENCHMARK_VCFS_THREADS=<n>`
bounds the worker threads, an environment variable and not an option, since the command line names
the preparation and the thread count changes no byte of it; a value that is not a whole number of
1 or more fails the run. Every thread count writes the same bytes to every file and to stdout.

Up to eight inputs are prepared at once. An input of up to a mebibyte on disk is prepared on one
thread. For any other, one thread reads the input and cuts it into pieces, the worker threads
inflate the pieces of a BGZF file and process the whole lines of every piece, one thread joins the
lines that cross from a piece to the next and writes the prepared file, and one digests what is
written. A gzip file that is not BGZF is inflated on its reading thread. Under `--max-variants`
the lines are counted in order on one thread and the input is read no further than the cap needs.
No input is held in memory: at most twelve pieces of about two mebibytes each are in flight for
one input, and a single line longer than that is held whole.

## Layout

- `src/main.rs`: the run: the options resolved, the inputs prepared across threads and reported in
  input order, the manifest written.
- `src/cli.rs`: the command line as `argparse` reads it, the usage and the help text.
- `src/paths.rs`: a directory as `pathlib` resolves it, the inputs of a directory in their order,
  the naming of a cache.
- `src/contigs.rs`: the human contig set and the set a `--canonical-contigs` value names.
- `src/gzip.rs`: the members of a gzip file as the `gzip` module reads them, and the blocks of a
  BGZF file.
- `src/lines.rs`: what a run of whole lines is written as, and the join of an input's pieces at
  the lines that cross them.
- `src/pipeline.rs`: one input from its file to its prepared file and digest: the reading thread,
  the worker threads, the writing and the digest.
- `src/manifest.rs`: the manifest's JSON and its time.
- `src/pyfmt.rs`: the CPython rules named above, each with a unit test that holds its bytes.
- `src/lines_tests.rs`, `src/pipeline_tests.rs` and `src/cli_parse_tests.rs`: unit tests of the
  line rules, of an input cut at every byte, of both ways an input is prepared (on one thread and
  across threads) against the goldens, and of the command line.
- `tests/fixtures/`: made-up inputs in three directories, `contigs/` (the two namings and contigs
  that are dropped), `shapes/` (line ends, bytes outside UTF-8, empty inputs, two inputs under one
  name) and `gzip/` (BGZF blocks, members with every header field, zero padding), two cache
  layouts and two contig lists; `expected/<case>/` holds the prepared files, `manifest.json` and
  `stdout.txt` of one run each, and `expected_help.txt` the help. A golden carries
  `{{INPUT_DIR}}`, `{{OUTPUT_DIR}}`, `{{FIXTURES}}` and `{{GENERATED_AT}}` for the values that
  depend on the checkout and the run. A `.gitattributes` beside this page keeps the line endings of
  every fixture as committed.
- `tests/golden.rs`: the binary over the three fixture directories under the options of every
  case, and again under several thread bounds, against the goldens.
- `tests/input_shapes.rs`: the binary over inputs made in the test: files larger than a piece as
  plain text, one gzip member and BGZF blocks, under several thread bounds, with and without a cap;
  names outside UTF-8; an input directory that is the output directory; directories named through
  links.
- `tests/command_line.rs`: the binary's exit code and bytes for help, refused command lines, a
  value that names no contig, a thread bound that is no count and runs that fail. `tests/common/`
  holds what the test files share.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-prepare-benchmark-vcfs/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

The root `.cargo/config.toml` applies, so a release binary carries the same `target-cpu` floor as
`vep`. Nothing regenerates the goldens: the files under `tests/fixtures/expected/` are the
contract, `cargo test` holds the binary to them, and a change of the contract changes them in the
same commit.

## Harness

`scripts/concordance/run_concordance.sh` prepares its inputs with this tree's build of the tool,
built through cargo before the first engine run as the tuple comparator is: a build that fails
stops the script with one line that names the build command. It passes the benchmark directory,
`<work-dir>/inputs`, `--chrom-mode auto`, the JSON cache directory, the assembly and the contig
option of the run, and `--max-variants` in smoke mode, and reads the prepared `*.vcf` files and
`manifest.json` from `<work-dir>/inputs`.
