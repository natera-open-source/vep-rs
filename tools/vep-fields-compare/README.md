# vep-fields-compare

A Rust build of `scripts/concordance/compare_vep_fields.py`, the per-column field comparator of the
concordance harness: the same command line, the same `fields.json` byte for byte, and the same stdout
summary. The Python stays the reference implementation; this binary exists because the comparator is
the longest single step of a concordance run on the genome-wide cells.

```text
vep-fields-compare --perl <reference output> --rust <scored output> --out <fields.json> [--examples 20] [--tmp-dir DIR]
```

## Contract

Every behaviour of the Python is reproduced on the bytes it produces, including the parts that are
Python-specific: universal newlines (`\n`, `\r\n` and a lone `\r` end a record), UTF-8 decoding with
`errors="replace"`, `str.strip()`'s whitespace set (U+001C to U+001F included), `str(pathlib.Path(arg))`
for every path option and the `inputs` block, `json.dumps(indent=2, ensure_ascii=True)` escapes,
`float.__repr__` and `round(x, 6)` for the agreements, `f"{n:,}"` digit grouping in the stdout lines,
and `Counter.most_common` tie order for `top_mismatches`. An input whose name ends in `.gz` is
gunzipped on read, as the Python's `gzip.open` reads it, a zero-byte file included; a gzip member
followed by zero padding, which the Python skips, is read as a corrupt stream. Both sides sort
through a child `sort -S 20%` under `LC_ALL=C` with `LC_COLLATE` removed, as the Python does at
`-S 25%`; the two sorts run in parallel and inherit `TMPDIR`, and a host without `sort` sorts each
side in memory, as the Python does. The Python encodes the rows it pipes to `sort` with the locale's
codec, so its output holds under a UTF-8 or C locale; the binary writes UTF-8 under any locale. Exit
codes: 2 on a usage error or when an input is not a file (`ERROR: [compare_vep_fields] not a file: <p>`
on stderr), 1 on any other failure (a `sort` failure, an unreadable input, an unwritable `--out`),
reported as `ERROR: [vep-fields-compare] <detail>` on stderr, 0 otherwise.

## Layout

- `src/input.rs`: records, lossy decoding, the prepared row (five key columns, nine value columns), the `sort` child.
- `src/merge.rs`: the sorted merge, duplicate-key collapse, per-field and per-Extra-key counts, mismatch shapes.
- `src/report.rs`: the JSON emitter and the stdout lines.
- `src/pyfmt.rs`: the CPython formatting rules listed above, each with a unit test pinning Python's output.
- `tests/golden.rs` and `tests/fixtures/`: the synthetic fixture pairs and the goldens the Python wrote
  over them (`expected_fields.json`, `expected_fields_examples1.json`, `expected_stdout.txt`, the
  `expected_disjoint_*` pair). `cargo test` runs the binary over the fixtures and compares bytes.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-fields-compare/Cargo.toml
cargo test --locked --manifest-path tools/vep-fields-compare/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-fields-compare/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-fields-compare/Cargo.toml -- --check
```

The root `.cargo/config.toml` applies, so a release binary carries the same `target-cpu` floor as `vep`
(`neoverse-v1` on aarch64 Linux, `x86-64-v3` on x86_64 Linux). CI runs the test, clippy and fmt
commands above in its Test, Clippy and Format jobs, builds the package on the declared `rust-version`
in MSRV, and regenerates the goldens with the Python in Comparator Tests and compares them with the
committed files, so a change to either side that moves the contract fails there first.

Regenerating the goldens after a deliberate contract change, from this directory:

```bash
python3 ../../scripts/concordance/compare_vep_fields.py --perl tests/fixtures/reference.txt --rust tests/fixtures/scored.txt --out tests/fixtures/expected_fields.json --tmp-dir "${TMPDIR:-/tmp}" > tests/fixtures/expected_stdout.txt
python3 ../../scripts/concordance/compare_vep_fields.py --perl tests/fixtures/reference.txt --rust tests/fixtures/scored.txt --out tests/fixtures/expected_fields_examples1.json --tmp-dir "${TMPDIR:-/tmp}" --examples 1 > /dev/null
python3 ../../scripts/concordance/compare_vep_fields.py --perl tests/fixtures/disjoint_reference.txt --rust tests/fixtures/disjoint_scored.txt --out tests/fixtures/expected_disjoint_fields.json --tmp-dir "${TMPDIR:-/tmp}" > tests/fixtures/expected_disjoint_stdout.txt
```

## Harness

`scripts/concordance/run_clone_measurement.sh` runs this binary in place of the Python when
`VEP_FIELDS_COMPARATOR` names an executable (a path, or a command on `PATH`), and the Python otherwise;
a value that names no executable file is warned about on stderr and the Python runs. Nothing exports
the variable by default, so a measurement runs the Python unless the caller builds the binary and
exports it.
