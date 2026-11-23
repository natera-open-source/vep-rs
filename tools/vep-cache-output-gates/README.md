# vep-cache-output-gates

The gate scorer of the cache comparison. `scripts/validation/compare_cache_outputs.sh` annotates one
input with two JSON caches, a baseline and a builder's, and this tool decides from the reports of
that run whether the builder's cache annotates as the baseline does: it scores five gates, prints a
line a gate, writes the verdict as a JSON report and leaves with the verdict as its exit code. The
lines and the report are fixed to the byte.

```text
vep-cache-output-gates BASELINE_METRICS BUILDER_METRICS CACHE_DIFF CONCORDANCE VERDICT \
    F1_FLOOR TOTAL_DELTA COLOCATED_DELTA FREQUENCY_DELTA CONSEQUENCE_DELTA
```

The ten operands are all given, in this order, and there is no option; `-h` or `--help` alone on
the command line prints the operands.

- `BASELINE_METRICS`, `BUILDER_METRICS`: the metrics reports of the two outputs, as
  `vep-file-census metrics` writes them.
- `CACHE_DIFF`: the report of `vep-cache-diff` over the two caches.
- `CONCORDANCE`: the summary of `vep-outputs-compare` over the two outputs (its `--summary-json`).
  An empty operand, or a path that names nothing, is a summary that is missing.
- `VERDICT`: the report to write. Its directory is not created.
- `F1_FLOOR` and the four ceilings `TOTAL_DELTA`, `COLOCATED_DELTA`, `FREQUENCY_DELTA` and
  `CONSEQUENCE_DELTA`: the thresholds of the five gates.

## Contract

**The gates.** Five, scored in this order:

| Gate | What it reads | It passes where |
| --- | --- | --- |
| `consequence_f1` | `aggregate.f1` of the concordance summary | the F1 is at or above `F1_FLOOR` |
| `total_annotations_delta` | `total_annotation_lines` of the two metrics reports | the delta is below `TOTAL_DELTA` |
| `colocated_variant_delta` | `existing_variation_count` | the delta is below `COLOCATED_DELTA` |
| `frequency_field_delta` | `frequency_field_counts`, a count a field | the largest delta of a field is below `FREQUENCY_DELTA` |
| `per_consequence_type_delta` | `consequence_counts`, a count a term | the largest delta of a term is below `CONSEQUENCE_DELTA` |

A delta is `abs(builder - baseline) / baseline` over a baseline count above zero. Over any other
baseline count it is 0.0 where the builder's count equals zero and 1.0 where it does not. A count
that a report does not hold is 0, and the last two gates score every name that either report
counts. The first gate fails, with `N/A` on its line, where the summary is missing, is empty or
holds no `aggregate`; an `aggregate` without `f1` scores 0.0. The verdict is `PASS` where all five
gates pass and `FAIL` otherwise.

**The outputs.** Standard output holds one line a gate, an empty line, `  Verdict: <verdict>` and
`  Report:  <the verdict's path>`. The report holds `verdict`, `gates` and `cache_diff_summary`.
Each of the five entries of `gates` holds `gate`, `threshold`, `actual` and `passed`; the two
gates over one count hold `baseline` and `builder` as well, the two gates over many hold
`details`, an object of `baseline`, `builder` and `delta` for each name in the code point order of
the names, and a first gate that found no summary holds `actual` as `null` and a `note`.
`cache_diff_summary` carries three values over from the cache comparison report:
`transcript_match_rate` and `variation_match_rate`, the `match_rate` of its `transcripts` and
`variations` blocks, and `info_json_matches`, the `matches` of its `info_json` block; each is `null`
where the report holds no such block or key.

**Exit codes.** 0 for a `PASS` verdict, 1 for a `FAIL` verdict, 2 where no verdict was reached. A
count of operands other than ten is refused with one `ERROR: [vep-cache-output-gates] <detail>`
line and the usage on standard error; any other run without a verdict writes one such line. What
ends a run is looked for in this order: a threshold that is no number; a report that cannot be
read or is not JSON, in the order of the operands; then, gate by gate, a report whose shape the
gate cannot read; then a cache comparison report that is no object; then a verdict that cannot be
written. The lines of the gates scored before that point stay on standard output, so a
concordance summary of the wrong shape leaves none and an unwritable verdict leaves all five, and
no report is written.

**The bytes.** Every byte read or written follows the CPython rule that defines it:

- A threshold is read as `float()` reads a word in ASCII: blanks at both ends, a sign, digits with
  single underscores between digits, a fraction, an exponent, and `inf`, `infinity` and `nan` in any
  case. A digit or a blank outside ASCII, which `float()` reads, is refused.
- A report is read as `json.load` reads it: UTF-8 without a byte order mark, one value, the
  constants `NaN`, `Infinity` and `-Infinity` beside strict JSON, a `\u` escape of a lone
  surrogate kept, a key that stands twice keeping its first place and its last value, an integer
  of up to 4,300 digits. An array or object inside more than 1,000 others is refused.
- A count and the F1 are numbers as Python holds them: an `int`, `true` and `false` among them as
  1 and 0, or a `float`. Two integers subtract and divide exactly, the quotient rounded once to
  the nearest double; an integer beside a float is that double first; the F1 is held to its floor
  by value. An integer is read up to 128 bits: a longer one is written to the lines and the report
  as it stands, and ends the run where a delta or the F1 would be computed with it.
- A value of another type in a number's place ends the run where Python raises `TypeError`: as a
  baseline count, as a builder's count beside a baseline count above zero, as the F1. As a
  builder's count beside any other baseline count it scores 1.0 and is written as `str()` writes
  it on the line and as `json.dump` writes it in the report; a string that is a lone surrogate or
  holds one has no UTF-8 for the line and ends the run there.
- On the lines a delta and the F1 are `f"{x:.6f}"`, a threshold is `float.__repr__` and a count is
  `str()`. In the report `actual` and `delta` are `round(x, 6)`, an F1 that is an integer staying
  one.
- The report is `json.dump(indent=2)` under `ensure_ascii=True`, followed by one newline, with
  `NaN`, `Infinity` and `-Infinity` for the floats JSON has no number for.
- A path is `str(pathlib.Path(operand))`: the reports and the verdict are opened there, and the
  `Report:` line prints that form in the path's own bytes. The concordance summary is looked for
  at that form and opened at the operand as it is written.

**Memory and cores.** The four reports are read whole, each on a thread of its own. They are
reports of kilobytes; the tool reads no annotation output and no cache.

## Layout

- `src/main.rs`: the run: the operands, the four reports read at once, the order of the errors,
  the verdict report, the exit code.
- `src/gates.rs`: the five gates, their lines, the report's entries and the values carried over
  from the cache comparison.
- `src/number.rs`: the arithmetic of a delta and of the F1 test over integers and floats.
- `src/json.rs`: a report as `json.load` reads it.
- `src/pyfmt.rs`: the CPython text rules listed above, each with a unit test pinning Python's
  output.
- `src/unicode.rs`: the code points `str.isprintable()` refuses, which `str()` of a list or an
  object writes as escapes (Unicode 16.0.0, the database of CPython 3.14).
- `tests/fixtures/`: `reports/` holds made-up metrics reports, cache comparison reports and
  concordance summaries, whole, of every shape the gates tell apart, and files that are no JSON;
  `cases.tsv` names for each case its four reports (`-` for the empty operand, `absent` for a path
  that names nothing), its five thresholds (a line feed written `\n`) and its exit code;
  `expected/` holds each case's standard output and, where one is written, its verdict report. A
  `.gitattributes` beside this page keeps their bytes as committed.
- `tests/golden.rs`: the binary over every case, its standard output, its report and its exit code
  compared with `expected/` byte for byte.
- `tests/command_line.rs`: help, a count of operands other than ten, paths with repeated and
  trailing separators, a verdict that cannot be written, a verdict path outside UTF-8, the order
  of the errors, an integer of more than 128 bits, a threshold outside ASCII. `tests/common/` holds
  what the two test files share.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-cache-output-gates/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

The package has no dependency. The root `.cargo/config.toml` applies, so a release binary carries
the same `target-cpu` floor as `vep`. Nothing regenerates the goldens: the files under
`tests/fixtures/expected/` are the contract, `cargo test` holds the binary to them, and a change of
the contract changes them in the same commit.

## Harness

`scripts/validation/compare_cache_outputs.sh` scores its gates with this tree's release build of
the tool, `tools/vep-cache-output-gates/target/release/vep-cache-output-gates`, as its last step:
the two metrics reports of its work directory, `cache_diff.json`,
`concordance_report/concordance.json`, `verdict.json` and its five `--gate-*` values as the ten
operands. It removes a `concordance.json` an earlier run left in the work directory before its
tuple comparison, so a comparison that fails is scored as a summary that is missing. It takes the
verdict from the exit code, prints it and leaves with 1 on `FAIL`; any exit code but 0 and 1
stops the script with 1. It builds nothing: a build that is missing or is not an executable
stops the script before its first engine run, with one line that names the build command above. `scripts/validation/test_cache_output_gates.sh` holds the script to this, and runs
its last step on the tree's own release build over three cases of `tests/fixtures/` where there
is one.
