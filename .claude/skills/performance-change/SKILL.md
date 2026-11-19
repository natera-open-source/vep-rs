---
name: performance-change
description: 'Land a change made for speed or memory: prove every output byte-identical first, then measure it with a protocol a maintainer can repeat. Use for a change whose purpose is wall time, memory, threading or allocation, with no change to what vep-rs writes.'
---

# A performance change

Output never changes for speed. A change made for speed or memory leaves every byte of every output
format as it was, over every committed corpus, and that identity is shown before any timing is
quoted. The measurement then states what it ran, where and how often, so a maintainer re-measures it
on another machine. Correctness outranks speed: a defect met on the way is fixed first, in its own
pull request through [fix-a-divergence](../fix-a-divergence/SKILL.md), and a speed-up that would
need a wrong row is not made. Branching, signing and the pull request follow
[AGENTS.md](../../../AGENTS.md).

## 1. Ask first

Open a report with the Performance report form, or comment on the issue that exists, before writing
code: a release branch may carry a change in the same code, and in-flight maintainer work closes a
duplicate with credit rather than merging it. State the hot path and the measurement that shows it
(a profile, a timing of one phase), not a guess.

Done when a maintainer has answered on the issue.

## 2. Keep the public API

A change that moves a `pub` item of a crate's surface, a plugin trait or the C ABI is an API change,
not a performance change: it needs its own entry under Changed, and the minor version is the
maintainers' call. Keep the signatures; change the bodies.

Done when `git diff origin/main -- crates | grep -E '^[-+]\s*pub '` prints nothing, or each line it
prints is a private module's item and is named under Upgrade impact as such.

## 3. Prove the output identical

Build the branch and `main` as release binaries, run both over every committed corpus in every
format, and compare the files. The default and tab headers carry an `## Output produced at` line and
the VCF header a `time="..."` attribute, both set at run time; mask those two before comparing, and
nothing else. The headers also echo the command line with the output file's name, so both sides
write to one file name in two directories. A Parquet output carries a run identifier in every part,
so it is compared through its `--tab` projection.

```bash
cargo build --release --locked                                   # the branch: target/release/vep
mkdir -p /tmp/vep-main && git archive origin/main | tar -x -C /tmp/vep-main    # main's tree, built apart
cargo build --release --locked --manifest-path /tmp/vep-main/Cargo.toml --target-dir /tmp/vep-main/target
mkdir -p /tmp/identity/main /tmp/identity/branch
mask() { sed -E 's/^## Output produced at .*/## Output produced at <time>/; s/ time="[^"]*"/ time="<time>"/'; }
field() { python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2], sys.argv[3]))' "$@"; }
for corpus in tests/golden/116/*/; do
  [ -f "$corpus/variants.vcf" ] && [ -d "$corpus/json_cache" ] || continue   # a tool, RefSeq or cache-layout corpus: see below
  species=$(field "$corpus/manifest.json" species homo_sapiens); assembly=$(field "$corpus/manifest.json" assembly "")
  for format in vep tab vcf json; do
    for side in main branch; do
      bin=target/release/vep; [ "$side" = main ] && bin=/tmp/vep-main/target/release/vep
      "$bin" -i "$corpus/variants.vcf" -o "/tmp/identity/$side/out.$format" --offline \
        --json_cache "$corpus/json_cache" --species "$species" --assembly "$assembly" --format vcf \
        --output_format "$format" --buffer_size 5000 --force_overwrite --quiet --no_stats
    done
    cmp <(mask < "/tmp/identity/main/out.$format") <(mask < "/tmp/identity/branch/out.$format") \
      || echo "DIFFERS: $corpus $format"
  done
done
```

That loop is a flagless run over each corpus's records and cache, not the corpus's own run: the
golden harness runs a corpus under its manifest's `flags` list, with `{corpus}` resolved to the
corpus directory, `{plugins}` to its `plugin_data/` and `{generated}` to the directory the
manifest's `generated` block writes, and under `--fasta` when the manifest names a `fasta` outside
the cache. Repeat the loop with those flags, and with the flags of every `flag_sets` entry whose run
the change can touch (a plugin change with the plugin's flags, a structural-variant change under
`--max_sv_size`), and over the corpora the loop skips with the input and cache their manifests name.
Compare Parquet with `python3 scripts/adapters/parquet_to_vep_tab.py --canonical` on both sides.
Then run the golden harness: `cargo test -p vep-cli --test golden --test format_parity`, which holds
every committed run to Ensembl VEP's output.

A row met on the way that differs from Ensembl VEP's is not classed as a defect before the index is
read. Search [docs/intended-divergences.md](../../../docs/intended-divergences.md), its sections and
the pages they link, for the flag, field or shape. An entry there is intended and stays as it is,
cited as the index's first column names it (the text of its row's first cell, inner backticks
dropped) or by the path of its page. Removing an entry needs the newer Ensembl VEP behaviour stated
in plain words, with the record that shows it. A difference outside the index is a vep-rs defect for
its own pull request through [fix-a-divergence](../fix-a-divergence/SKILL.md). A change made for
speed edits no entry and nothing in the reproduction registry under
`scripts/divergences/reproductions/`; a page the index links that it must edit is named at the gate
(step 5).

Done when no corpus prints `DIFFERS` in any format and the harness is green; the pull request's
Verification states "output byte-identical to `main` in every format over the committed corpora"
with the loop's command.

## 4. Measure it

A wall-time claim is a median, never a single run: at least five runs per machine of each binary,
the runs of the two binaries interleaved, on one input and cache. The spread (the slowest and
fastest run, or the 5th and 95th percentile) stands beside the median. Measure on both architectures
when you have them (x86_64 and arm64). Each machine is stated with its CPU model and core count,
memory, operating system and the storage the input and cache were read from, and the `--fork` count
is part of the protocol. Peak memory is the maximum resident set size of the run (`/usr/bin/time -v`
on Linux).

The input is public and named with its version and record count, or it is a committed corpus;
nothing that could identify a person is ever used as an input. The two binaries are built from the
same toolchain with the same flags, and the commit of each is named. `cargo bench -p vep-cli` and
the `bench_plugins` tool (`cargo build --release -p vep-cli --features bench --bin bench_plugins`)
measure one phase; they support a claim about that phase, and the whole-run figure above is the one
the pull request states.

Done when the pull request's Verification lists, per machine: the input, the run count, both
medians, both spreads, peak memory, the `--fork` count and the commands.

## 5. Gate

Commit with `-s` and the model declaration, then run the full set in
[AGENTS.md](../../../AGENTS.md#before-pushing) (the admission check included, on the saved
description). The `CHANGELOG.md` line under Unreleased carries the measured effect inline, as
CONTRIBUTING.md's grammar asks; a shape such as `<N> percent faster on <the input>, <k> threads`
states it from the medians of step 4. A maintainer re-measures on another machine before merge, so
the figure is repeatable, not final.

Open the pull request with the template's four sections and nothing else. Its Upgrade impact is
"None" unless a default changed. Its Intended divergences checks the first statement, which a change
made for speed states truthfully: it touches no entry of the index, and one that would is not a
performance change. The admission check refuses that statement to a change that edits a page an
index row links (`docs/cache-setup.md` and several pages under `docs/output/` among them): a change that
must edit such a page checks the second statement instead and names the page's entries one line
each, as the index names them, each kept. Done when every command passes and its output is quoted
under Verification.
