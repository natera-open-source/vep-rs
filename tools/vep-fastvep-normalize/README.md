# vep-fastvep-normalize

Rewrites a fastVEP tab output so that its intergenic rows carry the Feature_type that
Ensembl VEP writes. fastVEP through v0.4.0 writes `Transcript` in the Feature_type cell of a row
whose Gene and Feature cells are `-`; Ensembl VEP writes `-` there. Feature_type is one of the five
fields the tuple comparison keys a row on (Location, Allele, Feature, Feature_type, Consequence),
so without the rewrite every intergenic row of such an output scores as a discordant tuple. From
v0.5.0 fastVEP writes `-` there itself, and the tool rewrites no row of such an output.

```text
vep-fastvep-normalize <fastVEP output> <rewritten output>
```

The command line is the two file names and nothing else: there is no option, and a name that
starts with `-` is a file name. The rewritten output is created, or emptied where it exists. One
line on stderr counts the rows, and nothing goes to stdout:

```text
normalize: rewrote <rows rewritten>/<rows> intergenic Feature_type=Transcript → -
```

## Contract

**Lines.** A line that starts with `#` is written as it stands and is not counted. Every other
line is a row and is counted, an empty line among them. Cells are separated by tabs alone. A row is
rewritten when it holds six cells or more, its fourth and fifth cells (Gene, Feature) are `-` and
its sixth (Feature_type) is `Transcript`, each cell compared whole and byte for byte: the sixth
cell is written `-`, and no other byte of the row changes.

**Bytes.** A line ends at `\n`, at `\r\n` or at a lone `\r`, as universal newlines read it, and
every line ending is written `\n`; a `\r` inside a row therefore splits it into two rows. No other
character ends a line. A last row without a line ending gains one; a last `#` line without one
stays without. The input is read as UTF-8 and refused where it is not. A byte order mark is not
skipped, so a first line that starts with one is a row. The arrow of the stderr line is U+2192,
written as UTF-8.

**Order of work.** The input is opened before the output is made. A run whose input cannot be
opened leaves no output file, and a run given one name for both files reads an emptied file,
writes nothing and counts `0/0`.

**Exit codes.** 0 when the output and the count line are written. 1 otherwise, with one line on
stderr and no count line:

- `usage: vep-fastvep-normalize <in> <out>` for any number of arguments but two.
- `ERROR: [vep-fastvep-normalize] <detail>` for an input that cannot be opened or read or is a
  directory, an output that cannot be made or written, a thread bound that is no count, and an
  input that is not UTF-8 (`<input>: the byte at offset <n> is not UTF-8`, the offset counted
  from 0).

A run stopped by its input or its output leaves in the output file the rows written until then. A
run stopped at a byte outside UTF-8 leaves the rewrite of every line before the line of that
byte.

## Threads and memory

The input is read as a stream and never held whole. One thread reads it in runs of whole lines of
about 1 MiB, rewriting threads take the runs as they come, and one thread writes the rewritten
runs in input order. A rewriting thread is started when a run is ready and none waits for one, up
to the host's available cores, so a run whose reading or writing is the slower part keeps few
threads. `VEP_FASTVEP_NORMALIZE_THREADS=<n>` bounds the rewriting threads, `1` keeping the whole
run on the calling thread, and a value that is not a whole number of 1 or more is refused; it is
an environment variable and not an option, so that the command line stays the two file names. An
input of one run stays on the calling thread under any bound. Every bound writes the same bytes.

Two runs a rewriting thread and two more are held at most, each as read and as rewritten, so
memory follows the threads at work and not the input; a line longer than a run is held whole. A
pipe is read as a file is.

## What is not rewritten

The Feature_type cell is the one difference the tool removes, in an output of fastVEP through
v0.4.0. Two more are left as fastVEP writes them, the first in a tab output of every fastVEP
release through v0.6.2 and the second in an output of fastVEP through v0.5.0:

- The column count. The harness's own command asks fastVEP for its tab format, whose header
  names 17 columns where Ensembl VEP's default format has 14. After `Existing_variation` it
  names `IMPACT`, `DISTANCE`, `STRAND` and `FLAGS` as columns of their own, which Ensembl VEP
  writes as keys of its one `Extra` column. The 13 columns before them stand in the same order
  in both, the five fields of the key among them, so the difference lies outside the key and
  does not affect concordance.
- The consequence vocabulary does. fastVEP through v0.5.0 emits VCF 4.4 copy-number terms that
  Ensembl VEP release 115 never emits (`copy_number_increase`, `copy_number_decrease`,
  `copy_number_change`, `transcript_variant`; no Ensembl VEP tuple carries any of them), and it
  emits them beside the term Ensembl VEP uses, so `{copy_number_increase, feature_elongation}`
  cannot match `{feature_elongation}` on a key that includes the whole consequence set. Per dataset
  its `feature_elongation` and `copy_number_increase` counts are equal (174 on the GRCh37 1000
  Genomes Phase 3 chromosome 21 dataset, 1,774 on the GRCh38 1000 Genomes high-coverage dataset;
  `manuscript/data/f1_by_consequence_class_by_dataset.csv`), which identifies the added term, and
  not a missing annotation, as the cause. From v0.6.0 fastVEP emits none of the four terms.

The vocabulary is left as fastVEP writes it, and the reason is the difference between the two
cases. Writing `-` for `Transcript` on an intergenic row corrects the formatting of a field both
engines agree about. Rewriting another engine's consequence vocabulary would change what is
measured: the reported quantity is concordance with Ensembl VEP's output on that key, because the
claim under test is drop-in replacement, and a pipeline whose filters are written against Ensembl
VEP's consequence sets does get a different answer from fastVEP's. That discordance is real and
belongs in the result.

## Layout

- `src/main.rs`: the command line, the order in which the two files are opened, the stderr lines
  and the exit codes.
- `src/rewrite.rs`: the rewrite of one run of whole lines and its row counts, with unit tests of
  every row shape and line ending.
- `src/stream.rs`: the reader that cuts the input into runs of whole lines, the rewriting threads,
  the writer that keeps input order and the bound on the buffers in work; its unit tests run
  inputs at every run length under several thread counts.
- `src/threads.rs`: the thread bound.
- `tests/fixtures/`: made-up fastVEP outputs (`<name>.tab`) and the goldens that hold the
  rewritten output and the stderr line of each to the byte (`expected_<name>.tab`,
  `expected_<name>.stderr`). `not_utf8.tab` and `not_utf8_at_the_end.tab` are `raw.tab` followed
  by a row that is not UTF-8, so the output they leave is `expected_raw.tab`. A `.gitattributes`
  beside this page keeps the line endings of the fixtures as committed.
- `tests/golden.rs`: the binary over every fixture under four thread bounds, the refused command
  lines and files, a pipe as input, and a generated input of several runs held to a rewrite
  written out cell by cell.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-fastvep-normalize/Cargo.toml
cargo test --locked --manifest-path tools/vep-fastvep-normalize/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-fastvep-normalize/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-fastvep-normalize/Cargo.toml -- --check
```

The root `.cargo/config.toml` applies, so a release binary carries the same `target-cpu` floor as
`vep` (`neoverse-v1` on aarch64 Linux, `x86-64-v3` on x86_64 Linux).

The goldens are the contract in bytes, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the goldens it moves with the binary's own output, so that the
diff of the goldens shows the change:

```bash
cd tools/vep-fastvep-normalize
target/release/vep-fastvep-normalize tests/fixtures/raw.tab tests/fixtures/expected_raw.tab \
    2> tests/fixtures/expected_raw.stderr
```

## Harness

`scripts/concordance/run_clone_measurement.sh` runs this tree's release build of the tool,
`tools/vep-fastvep-normalize/target/release/vep-fastvep-normalize`, on each SNP/indel cell of the
fastVEP engine: the raw output of the timed run goes in, and the cell's `rust/output.txt`, which
the tuple comparator scores, comes out. The script builds nothing: a build that is missing or is
not an executable stops it before its first engine run, with one line that names the build command
above. The tool is not asked for under `--walltime-only`, nor for another engine. A rewrite that
fails is recorded as a `normalizer_failed` row of the run's CSV and the cell is not scored; the
raw output is never compared in its place.
