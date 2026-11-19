# vep-plugin-pairs

The heavy path of `scripts/concordance/score_plugin_tuples.py`, the plugin tuple scorer, on every
core. The scorer writes a line per plugin value into `sort` from one process; this binary reads the
reference and each scored engine's output on every thread, extracts the populated plugin fields by
the scorer's extraction rules, splits them by Location so that equal Locations stay together, sorts
them and pairs the values that are equal as text.

The scorer drives it and stays the one place where the plugin registry, the number rule, the
documented classes, the residual files and the report are defined. Everything the binary does not
pair as equal text is the remainder, which it hands to the scorer in the order and under the keys
the scorer's own merge hands them to its settling step. A run through the binary writes the same
report and the same residual files, byte for byte, as a run on the scorer's own extraction, sort and
merge, which stays the reference implementation and the check.

```text
python3 scripts/concordance/score_plugin_tuples.py <the scorer's options> \
    --pairing-tool <vep-plugin-pairs> [--pairing-threads N]

vep-plugin-pairs --reference <output> --engine <output or tuple file> [--engine ...] --out <dir> \
    [--rows] [--threads N] [--shards N] [--distinct-keys N] [--block-bytes N]
```

## Contract

**The dialogue.** The registry lives in `compare_plugin_outputs.py`, so the binary resolves no field
name itself. The first time it meets a name it writes `field <hex of the name>` on standard output
and reads the answer on standard input: `source <name>` for a field of a plugin or custom source of
the run, `core` for a column of Ensembl VEP's own, `unclaimed` for a key no source claims, `refused`
where the registry refuses the run. A Location with a byte outside ASCII, or longer than 255 bytes,
is asked the same way (`location <hex>`) and answered with its canonical form, so that Python's
reading of a digit decides it. Nothing else is written on standard output.

**Extraction**, as `extract` of the scorer. A line ends at its line feed and loses every carriage
return before it; an empty line is skipped. The reference is a default-format output; an engine is
one too, or a tuple file when its first data line has five columns and no `#` line comes before it.
In a default-format output a line that opens with `#` is no row; one that opens with a single `#`
names 14 columns and the last of them is handed over as the column header. A row has 14 columns and
a tuple line five, or the run is refused with the line's number. The Location is written with each
coordinate as its number; the Allele and the Feature lose their surrounding ASCII whitespace, the
vertical tab included. Of the Extra column, each `;`-separated item with an `=` and a value that is
neither empty nor `-` is a field; its value has `%3D` and `%3B` decoded. A tuple line's value is
taken as written. Every line of a tuple file is counted, and those whose populated value belongs to
no source of the run are counted apart. A key no source claims is counted per input, up to
`--distinct-keys` distinct keys. A gzip input is unpacked whole, every member of it, zero padding
after a member skipped, and digested as stored.

**Pairing.** Under each key, (Location, Allele, Feature) in the feature projection and (Location,
Allele) in the allele projection, and each field, the two sides' values equal as text pair one to
one. The feature projection keeps a value as often as its rows write it; the allele projection keeps
it once. A key under which a side writes one field more than once is counted once for the field's
source.

**The hand-over**, in `--out`:

- `pairs.json`: per input its kind, the sha256 of the file as stored, its rows or tuple lines, the
  tuples outside every source, the column header and the unclaimed keys with their counts; per
  engine and projection, each field's values on the reference's side and on the engine's with the
  pairs equal as text, the shared keys per source and the remainder's lines per side. Field names
  and the header are hexadecimal, since an input can carry bytes no text holds.
- `input<n>.<projection>.reference` and `input<n>.<projection>.engine`, for engine `n` counted from
  1: the remainder as tuple lines, in the order `sort` gives them under C collation.
- With `--rows`, `input<n>.rows.reference` and `input<n>.rows.engine`: the reference's rows and the
  engine's rows at every Location of the feature remainder, each as its key, its ordinal among the
  input's data rows on twelve digits and the row as written, sorted likewise. A class reads the row
  of the pair in hand under its key and, for a record's other alleles, each side's rows of the same
  input record on the Feature. A tuple file states no rows, so its tuples under the remainder's
  keys stand for them.

The files are the same bytes whatever `--threads`, `--shards` and `--block-bytes` are.

**Refusals.** A refusal is one `ERROR: [vep-plugin-pairs] <message>` line on standard error and exit
status 2, and the scorer ends in the same message under its own name. The messages of a line with
the wrong number of columns, of a column header that names other than 14 columns, of a reference
without a data row and of too many distinct unclaimed keys are the scorer's own, word for word; a
gzip stream that cannot be read is named with the decoder's reason, which is not Python's. Where an
input gives more than one cause, the first line at fault is the one named.

**Scratch and memory.** The scratch, `--out/shards`, holds one record per row with a claimed field
(its key once, then its fields as numbers and values) and the unpacked copy of a gzip input, and is
removed on every exit the binary takes itself. Each thread holds one shard of each side at a time;
`--shards` defaults to one per 128 MiB of the longest input. The Locations of the feature
remainder, and with `--rows` the rows written for them, are held in memory, so the binary is sized
for outputs that mostly agree.

`VEP_PLUGIN_PAIRS_TIMINGS` in the environment prints the time each pass ended on standard error,
which the scorer passes on.

## Layout

- `src/main.rs`: the command line, the shards' bounds and the order of the passes.
- `src/input.rs`: an input's digest, its plain copy, its kind and its blocks.
- `src/resolve.rs`: the dialogue and the fields met.
- `src/extract.rs`: the extraction rules and the shard records.
- `src/pair.rs`: the pairing of one shard, the counts and the remainder.
- `src/rows.rs`: the rows of the remainder.
- `src/handover.rs`: `pairs.json`.
- `src/order.rs`: the byte rules every pass shares, each with a unit test.
- `tests/pairs.rs`: the binary over small outputs, one test per extraction rule, the planted
  missing, extra and changed values, the rows, and the same bytes on any number of threads.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-plugin-pairs/Cargo.toml
cargo test --locked --manifest-path tools/vep-plugin-pairs/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-plugin-pairs/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-plugin-pairs/Cargo.toml -- --check
```

The scorer's test module holds the two paths to each other. Its tests named for the pairing tool
score random outputs on both paths and compare the reports and the residual files byte for byte;
they find the binary under the cargo target directory, or where `VEP_PLUGIN_PAIRS` names it, and
skip without one. With `SCORE_PLUGIN_TUPLES_PATH=pairing-tool` every test of the module that scores
runs through the binary:

```bash
python3 -m pytest -q scripts/concordance/test_score_plugin_tuples.py -k pairing_tool
SCORE_PLUGIN_TUPLES_PATH=pairing-tool python3 -m pytest -q scripts/concordance/test_score_plugin_tuples.py
```

CI runs the test, clippy and fmt commands above in its Test, Clippy and Format jobs, runs both
pytest commands after the Test job's build, and builds the package on the declared `rust-version` in
MSRV.
