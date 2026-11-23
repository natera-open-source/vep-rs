# vep-file-census

The passes the concordance and validation launchers make over whole files: counting the records
of an input VCF and the rows of an output, digesting inputs, outputs and every file of a cache
tree, and the row passes of an `awk` program. Each subcommand writes on standard
output, byte for byte, what the shell pipeline it stands for writes for the same files, and leaves
with that pipeline's status where a launcher tests it. `line-sum` stands for no pipeline: it
reduces a file to the count of its lines and a sum that does not depend on their order, so that
two outputs that hold the same rows in another order are compared without a sort. The files are
read as streams across the host's cores: a release-scale input or output is tens of gigabytes,
and a pass over it on one core is minutes of a run that annotates nothing meanwhile.

```text
vep-file-census count [--gunzip] RULE [RULE ...] FILE
vep-file-census digest --algorithm {sha256,md5,xxh128} [--algorithm ...]
    [--names {lines,read,nul,words}] [--directory DIR] [--gunzip]
    [--without-lines-starting TEXT | --after-leading-lines-starting TEXT]
    [--chunk-bytes N] [FILE ...]
vep-file-census strip-chr FILE
vep-file-census metrics OUTPUT REPORT SECONDS
vep-file-census contigs [--gunzip] FILE
vep-file-census slice --records N --tally FILE [--label TEXT] [--table FILE] [--gunzip] FILE
vep-file-census fetched --list TSV --root DIR [--jobs N] [--report TSV]
vep-file-census line-sum [--gunzip] [--without-lines-starting TEXT] FILE [FILE ...]
```

An option takes its value as the next word or after `=`; `--` ends the options, so that a file
whose name starts with `-` can be named. A command line the tool cannot take is refused with
exit code 2, one `ERROR: [vep-file-census] <detail>` line and the usage on stderr. Every other
failure is one such line, without the usage.

## Contract

### count

One count a rule, one a line, in the order the rules are given, from one pass over FILE.

| Rule | The pipeline it stands for | A line counts where |
| --- | --- | --- |
| `--starting TEXT` | `grep -c '^TEXT' FILE` | it starts with TEXT |
| `--not-starting TEXT` | `grep -vc '^TEXT' FILE` | it does not start with TEXT |
| `--containing TEXT` | `grep -c 'TEXT' FILE` | it holds TEXT |
| `--info-field NAME` | `awk -F'\t' '!/^#/ && index(";" $8, ";NAME=") {c++} END {print c + 0}' FILE` | it does not start with `#` and its eighth tab-separated column holds the key NAME: `NAME=` at its start or after a `;` |
| `--extra-key NAME` | `awk -F'\t' '!/^#/ && index(";" $NF ";", ";NAME=") {c++} END {print c + 0}' FILE` | it does not start with `#` and its last tab-separated column holds the key NAME |

TEXT is a fixed text of one line, not a pattern. For the three `grep` rules a line ends at a
newline and at a zero byte, as `grep` reads a file that holds a zero byte, and a last line without
either is a line; for the two `awk` rules a line ends at a newline alone.

`--gunzip` counts what `gzip -dc FILE` writes (`gzip -dc FILE | grep -vc '^#'`, `zcat FILE | ...`).

Exit codes, those of the pipeline under `set -o pipefail`: 0; 1 where a `grep` rule counts no
line (the count `0` is printed), an `awk` rule never changing the status; 2 for a file that cannot
be opened, of which no count is printed, and for one that cannot be read to its end, of which the
lines read are counted. Under `--gunzip` a file that cannot be opened, is no gzip file or ends
early prints the count of what was inflated and leaves with 1, and bytes after the last member
that start no member (`gzip` warns of trailing garbage) leave with 2 after the count.

### digest

One line a file and algorithm, as `sha256sum`, `md5sum` and `xxh128sum` print it: the digest in
hexadecimal, two spaces, the name. The first two escape a name that holds a backslash or a
newline (`\\`, `\n`) and start its line with a backslash; `xxh128sum` writes the name as it is.
With several `--algorithm` the file is read once and its lines stand in the options' order. The
files of a run are read at once on as many threads as the run has, and their lines are written in
the order the names were given. A file that cannot be read has no line.

The names are the operands (`-` is standard input, named `-`), or with `--names` the names on
standard input, in one of the four ways a pipeline hands names to a digest program:

| `--names` | The pipeline it stands for | A name is | Exit code |
| --- | --- | --- | --- |
| `lines` | `... \| while IFS= read -r f; do sha256sum "$f"; done` | a line, whole; a last line without a newline is not read | that of the last name: 0, or 1 where it has no line |
| `read` | `... \| while read -r f; do md5sum "$f"; done` | a line without the blanks at either end | as `lines` |
| `nul` | `... \| xargs -0 sha256sum` | what a zero byte ends | 123 where a name has no line |
| `words` | `... \| xargs sha256sum` | what blanks and newlines separate; a quote of either kind and a backslash are read as `xargs` reads them | as `nul` |

Under `nul` and `words` no name at all prints the one line the digest program writes of an empty
standard input, since `xargs` starts its command once on an empty list. Over operands the exit
code is 1 where a file has no line. `--directory DIR` reads the names inside DIR and prints them
as given: the names of `(cd DIR && find . -type f -print0 | LC_ALL=C sort -z)`.

A body is digested as the digest program reads it from a pipe, so its line names standard input,
`-` or, for `xxh128`, `stdin`:

| Option | The pipeline it stands for |
| --- | --- |
| `--gunzip` | `gzip -dc FILE \| sha256sum` |
| `--without-lines-starting TEXT` | `grep -a -v '^TEXT' FILE \| sha256sum`: the lines that do not start with TEXT, a last line without a newline gaining one. A body without a line leaves with 1, as `grep` does. With `--gunzip`: `gzip -dc FILE \| grep -a -v '^TEXT' \| sha256sum` |
| `--after-leading-lines-starting TEXT` | `tail -c +$((n + 1)) FILE \| sha256sum`, `n` the bytes of the file's leading lines that start with TEXT |
| `--chunk-bytes N` | `... \| split -b N --filter=sha256sum \| awk '{print $1}' \| sha256sum`: the digest of the digests of the body's chunks of N bytes, one hexadecimal digest a line in chunk order; an empty body is the digest of nothing. The chunks are read and digested across the threads |

A body of a file that cannot be opened is the digest of nothing, with exit code 1; a gzip file
that ends early is digested to where it ends, with exit code 1, and trailing garbage after its
last member leaves with 2.

### strip-chr

The lines of a VCF without the `chr` prefix of its contig names, as this program writes them:

```text
zcat -f FILE | awk 'BEGIN{OFS="\t"} /^##contig=/ {sub(/ID=chr/,"ID=",$0); print; next}
                    /^#/ {print; next} {sub(/^chr/,"",$1); print}'
```

A `##contig=` line loses the `chr` of its first `ID=chr`; any other line that starts with `#` is
written as it is. A record whose first column starts with `chr` loses that prefix and is written
with one tab between its columns, the columns being what blanks (spaces and tabs) separate, as
`awk` builds a record again from its fields once one of them changes; any other record is written
as it is. Every line written ends in a newline. FILE is inflated where it starts as a gzip file
and read as it is otherwise, and bytes after a gzip member that start no member are read as they
are, as `zcat -f` reads them. Exit code 1 for a file that cannot be read.

### metrics

The annotation metrics of one default-format output: `REPORT` is written, then seven summary
lines on standard output. The bytes follow the CPython rules this paragraph names. OUTPUT is
UTF-8 text read under universal newlines: a row ends at `\n`, `\r\n` or a lone `\r`. A line that
starts with `##` is skipped; a line that starts with `#` names the columns, split at tabs after
`lstrip("#")`, for the rows after it; an empty line is skipped; every other line is a row. A row
under seven columns is counted and read no further. The Consequence column is split at `,` and
each term stripped with `str.strip()`; the Existing_variation column counts where it is neither
empty after `str.strip()` nor `-`; the Extra column is the first column named `Extra` after
`str.strip()` where a line named the columns and the row reaches it, else the fourteenth, and is
split at `;` into `KEY=VALUE` parts, a key that stands twice keeping its last value. The report
is `json.dump(indent=2)` of these keys in this order, then a newline: `total_annotation_lines`,
`consequence_counts` (the terms in `Counter.most_common()` order: the most frequent first, equal
counts in the order first met), `existing_variation_count`, `clin_sig_count`, `hgvsc_count`,
`hgvsp_count` (the rows whose `CLIN_SIG`, `HGVSc`, `HGVSp` value is not empty),
`frequency_field_counts` (the same for `AF`, `gnomADe_AF`, `gnomADg_AF`, `MAX_AF`, `EUR_AF`,
`AFR_AF`, `AMR_AF`, `EAS_AF`, `SAS_AF`) and `elapsed_seconds`, SECONDS as `int()` reads it in
ASCII digits. Exit code 1, with no report and no summary, for a SECONDS that is no whole number,
an OUTPUT that cannot be read or is not UTF-8, and a REPORT that cannot be written.

The rows are read across the threads. A line that names the columns after the first row changes
what the rows after it mean, so such a file is read again on one thread, in order.

### contigs and slice

`contigs` is the tally of a VCF's records by contig, in the order first met:

```text
awk -F'\t' '!/^#/ {if (!($1 in c)) o[++k] = $1; c[$1]++}
            END {for (i = 1; i <= k; i++) print o[i] "\t" c[o[i]]}'
```

`slice` is the pass that follows it over the same file: every line that starts with `#`, and of
each contig the records at which `floor(i * q / c)` moves, `i` the record's number on its contig
from 1, `c` the contig's count in the tally and `q` its share of the `--records N` asked for over
the `k` lines of `--tally`: `floor(N / k)`, one more for the first `N mod k` contigs of the
tally, and no more than `c`. The products and quotients are `awk`'s, in double precision.
`--table FILE` gets one line appended for each contig the records name, in that order:
`--label`, the contig, its count and its share; a file without a record leaves no table. The
slice is selected in file order on one thread, the file inflated across the others.

`--gunzip` reads what `gzip -dc FILE` writes. Both write what they read before a failure and
leave with 1 after it, as a pipeline whose reader failed.

### fetched

Fetched files held to a list of key, bytes and ETag: every file under `--root` against the size
and the ETag an object store states for its object. The list is tab-separated lines
`key<TAB>bytes<TAB>etag` without a header, cut as `str.splitlines()` cuts them; a blank line and a
line that starts with `#` are skipped. An ETag is 32 hexadecimal digits, the md5 of the file, or
`<32 hex>-<parts>` for an object uploaded in parts: the md5 of the parts' binary md5s joined,
then the count of parts. A third cell of 64 hexadecimal digits is the file's sha256.

The verdicts, one a key: `md5`, `sha256`, `multipart-md5 <MiB>` (the ETag is that of an upload in
parts of that size), `bytes-alone` (the size holds and no part size tried gave the ETag) and
`DIFFERS` (missing, unreadable, of another size or of another digest, or an ETag of neither
shape). The part sizes tried are the whole MiB for which the size makes the stated count of
parts: 8, 16, 32, 64, 128, 256, 512 and 1024 MiB first, then every other from 5 to 5120 MiB
upwards, twelve at most, and one alone where the object is one part.

`--report` gets one line a key in the list's order: key, bytes, verdict. Standard output names
every key that `DIFFERS` or is held by `bytes-alone`, key and verdict, and ends
`held: <a> by md5, <b> by multipart md5, <c> by bytes alone, <d> differ of <n>`, with
`<e> by sha256, ` before the differing where a key is held so. `--jobs` files are read at once
(default: the threads of the run, 32 at most). Exit codes: 0 when no key differs, 3 when one does,
2 for words the tool cannot take (a line without three columns, bytes that are no number, a root
that is no directory, fewer than one job), with nothing checked and no report.

### line-sum

The lines of a file in any order: one line a FILE, in the order the files are given, that holds
the count of the file's lines, a tab, and the sum of the lines' SHA256 digests modulo 2**256.
Two files that hold the same lines, each as often, print the same line wherever the lines stand
in them; a line that stands twice is counted and summed twice.

The bytes follow the CPython rules this paragraph names. A line ends as universal newlines end
one, at `\n`, at `\r\n` and at a lone `\r`, and a last line needs no line end; an empty line
is no line. The bytes of a line without its line end are digested as they stand. Each digest is
the integer `int.from_bytes(hashlib.sha256(line).digest(), "big")`, and the sum prints as
`format(total % 2**256, "064x")`: 64 hexadecimal digits in lower case. A file without a line
prints `0`, a tab and 64 zeros.

`--without-lines-starting TEXT` leaves out the lines that start with TEXT, a fixed text of one
line: with `#`, what is counted and summed is the rows of a tab or default-format output
without its header lines, or the records of a VCF. `--gunzip` sums what `gzip -dc FILE` writes.

Exit codes: 0; 1 for a file that cannot be opened or read to its end, a gzip file that ends
early among them, which has no line, since the sum of a part of a file would pass for the
file's, and after which no file is read; 2 where bytes after the last member of a gzip file
start no member, after the line of what the members inflate to.

### Files and threads

A plain file is cut into ranges that the threads read themselves. A block-gzipped file (BGZF) is
cut into runs of whole blocks, which inflate independently on the threads; any other gzip file
is one stream and is inflated member by member on one thread, its lines still processed on the
others. A member is held to its header checksum, its CRC32 checksum and its length. After a
member, zero bytes to the end of the file are padding; any other bytes that start no member are
trailing garbage and end the reading. A stream that ends early is read to where it ends; what it
inflates to there may be a few bytes more than `gzip` writes of it. No more of a file is held
than the pieces in flight and the line being joined across two of them; a digest holds its state,
and a line sum its count and its 32 bytes, to which the threads' sums are added as they arrive.

`VEP_FILE_CENSUS_THREADS=<n>` bounds the threads at work; unset or empty, the host's available
cores are the bound. It is an environment variable and not an option, since the command line
names the pass and the thread count changes no byte of it. A value that is not a whole number
of 1 or more is refused.

## Layout

- `src/main.rs`: the command line, the subcommands and their exit codes.
- `src/pipeline.rs`: the source thread, the worker threads and the ordered hand-over; a file as
  pieces, and a file as runs of whole lines joined across pieces; the thread bound.
- `src/gzip.rs`: the members of a gzip file as `gzip -dc` and `zcat -f` read them, and the blocks
  of a BGZF file.
- `src/count.rs`: the five rules.
- `src/digest.rs`: the three algorithms, the digest line, the bodies and the digest in chunks.
- `src/names.rs`: the four forms of a list of names.
- `src/strip.rs`: the contig-prefix pass.
- `src/metrics.rs`: the metrics of an output and the CPython rules its report follows.
- `src/contigs.rs`: the tally and the slice.
- `src/fetched.rs`: the list, the part sizes and the verdicts.
- `src/line_sum.rs`: the lines of a run, their digests and the sum modulo 2**256.
- `src/*_tests.rs`: unit tests of each of them.
- `tests/fixtures/`: made-up inputs (VCFs plain, gzipped, block-gzipped, cut short, with bytes
  after the member, a file with zero bytes, outputs in three formats, a small tree, lists of
  names, the same lines in two orders) and `expected/`, the goldens: each is the standard output
  of the shell pipeline its case names, written by that pipeline, and each `line_sum_*.out` is
  written by the CPython rule of the `line-sum` section. A `.gitattributes` beside this page
  keeps their bytes as committed.
- `tests/count.rs`, `tests/digest.rs`, `tests/rows.rs`, `tests/metrics.rs`, `tests/fetched.rs`,
  `tests/line_sum.rs`: the binary held to the goldens and to the exit codes, on one thread, on
  two and on the host's default.
- `tests/command_line.rs`: refused command lines, the two ways an option takes its value, names
  the digest programs escape, a gzip stream that ends early.
- `tests/threads.rs`: generated files of several pieces, plain, gzipped in members and
  block-gzipped, across thread counts.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-file-census/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

The root `.cargo/config.toml` applies, so a release binary carries the same `target-cpu` floor
as `vep`. Nothing regenerates the goldens: the files under `tests/fixtures/expected/` are the
contract, `cargo test` holds the binary to them, and a change of the contract changes them in
the same commit, each written by the pipeline its case names or, for `line-sum`, by the CPython
rule its section states.

## Harness

Four launchers read whole files through this tree's build of the tool, and no environment
variable selects another program:

- `scripts/concordance/generate_reference_output.sh` counts each input's records and each
  output's rows and annotated units (`count`) and digests the input and the output of every
  provenance sidecar (`digest`).
- `scripts/concordance/run_clone_measurement.sh`, in its Ensembl VEP cells, takes the `chr`
  prefix off a raw input (`strip-chr`) and counts the rows and records of its output-volume
  floor (`count`).
- `scripts/concordance/run_concordance.sh` digests both caches, the FASTA and every prepared
  input for its provenance record and its output cache key (`digest`).
- `scripts/validation/capture_baseline.sh` computes the metrics of the output (`metrics`) and
  the checksums of the cache's JSON files (`digest`).

The first, the second and the fourth run the release build,
`tools/vep-file-census/target/release/vep-file-census`, and build nothing: a build that is
missing or is not an executable stops the script with one line that names the build command,
before its first engine run. A dry run of `generate_reference_output.sh` and the vep-rs and
fastVEP engines of `run_clone_measurement.sh` read no whole file and ask for no build.
`run_concordance.sh` builds the tool through cargo, as it builds the engine it runs, and a build
that fails stops it the same way.
