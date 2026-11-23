# vep-storable-transcript-classes

Classes how two builds of an Ensembl VEP Storable cache differ on the transcripts whose content
differs. It reads two dumps of `tools/vep-storable-dump` for the same transcript ids and names, per
transcript, the groups of fields that differ once number quoting, `dbID` and `_gene_version` are
neutralised: `peptide` (the stored translation), `predictions` (SIFT and PolyPhen matrices),
`protein_features` (domain annotations), `other` (anything else, with its keys). With `--json-cache`
it also says which build the stored translation of a vep-rs JSON cache agrees with.

```text
vep-storable-transcript-classes --version-a 113 --dump-a a.txt --version-b 115 --dump-b b.txt \
    --out-dir <dir> [--json-cache <cache dir> [--json-cache-label <name>]]
```

- `--version-a`, `--version-b`: the names of the two builds, used in file names and messages.
- `--dump-a`, `--dump-b`: the two dumps.
- `--out-dir`: where the table and the summary go. It is made with its parents.
- `--json-cache`: a vep-rs JSON cache, the directory that holds `transcripts/`.
- `--json-cache-label`: the name the summary gives the JSON cache (default: the directory's name).

The command line is read as CPython's `argparse` reads it: an option may be shortened to any prefix
that fits it alone, takes its value as the next argument or after `=`, and a later value replaces an
earlier one. `-h` and `--help` print the usage and the options.
`VEP_STORABLE_TRANSCRIPT_CLASSES_THREADS` sets the number of worker threads (default: the host's
cores).

## Contract

Every output is fixed to the byte: the two files, standard output, the lines on standard error and
the exit code. `tests/fixtures/expected/` holds the reference bytes of 26 runs made by CPython 3.14,
and `cargo test` holds the binary to them. Wherever a byte could be written in more than one way it
follows a CPython rule, named here by its CPython name.

**A dump** is read as UTF-8 with `errors="replace"` (a byte sequence that is not UTF-8 reads as
U+FFFD) and universal newlines (`\r\n` and a lone `\r` end a line as `\n` does). The text is cut at
every `### `, wherever it stands; what precedes the first one is dropped. In each piece the first
line is the transcript id, without the white space `str.strip()` takes off, and the rest are the
block's lines. A later block of an id replaces an earlier one. In every line an integer in quotes
at the line's end loses its quotes (`=> '12',` becomes `=> 12,`; the digits are decimal digits of
any script, Unicode 16.0.0), and a line that holds `'dbID' =>` or `'_gene_version' =>` is then left
out.

**The ids** are those both dumps hold, in code point order. An id in one dump only ends the run with
`ERROR: [classify_storable_transcript_differences] <n> ids are in one dump only: <the first five as
repr() writes a list>` and exit code 1, before anything is written.

**The lines that differ** are those `difflib.unified_diff(a, b, n=0)` marks with `-` and `+`: every
line outside the matching blocks of `difflib.SequenceMatcher` with its default `autojunk`, under
which a line that fills more than one place in a hundred of a block of 200 lines or more takes no
part in the search of the blocks. A removed line that opens with `--` and an added one that opens
with `++` are passed over, as lines the diff itself would write.

**The classes.** In each line that differs the first `'<word>' =>` names a key, a word of ASCII
letters, digits and `_`. The key gives the class:

| Class | Keys |
|---|---|
| `peptide` | `peptide`, `translateable_seq` |
| `predictions` | `sift`, `polyphen_humdiv`, `polyphen_humvar`, `protein_function_predictions`, `matrix`, `matrix_compressed`, `peptide_length`, `sub_analysis`, `translation_md5` |
| `protein_features` | `hseqname`, `_display_label`, `analysis`, `protein_features`, `hstart`, `hend`, `score`, `percent_id`, `p_value`, `logic_name`, `db`, `db_version`, `interpro_ac`, `idesc`, `ilabel`, `start`, `end` |
| `other` | any other key |

A line without a key gives no class. A transcript with no line that differs is
`identical_after_normalisation`.

**The JSON cache.** The shards are `transcripts/<directory>/<name>.json` and `<name>.json.gz`, any
directory and any name (those that open with a dot included, a directory that is a symbolic link
followed), in the order `sorted()` gives their `pathlib` paths: by directory name, then by file
name. Each is read as `json.load` reads it and is a list of transcripts or an object of such lists.
A transcript is an object; its id is its `stable_id` up to the first `.` (a false value is the
empty id); the first transcript of an id in that order gives the id's stored translation, the
`peptide` of its `variation_effect_feature_cache` (a cache that is a false value has none). The
match of a transcript is `absent` when the JSON cache has no translation for it; else, where both
dumps hold the same translation (the first `'peptide' => '<capital letters and *>'` of the block, or
none), `both` or `neither`; else the version whose translation it equals, or `neither`.

**The table** is `differing_<a>_<b>_classes.tsv`: `transcript_id`, `classes` (joined by `+` in name
order), `other_keys` (joined by `,` in name order) and `json_cache_peptide_matches` (empty without
`--json-cache`), one line an id.

**The summary** is `differing_<a>_<b>_classes.json`, `json.dumps(summary, indent=2)` and a newline:
`version_a`, `version_b`, `transcripts`, `per_class`, `class_combinations` and, with `--json-cache`,
`json_cache` (the label, or the directory's `Path.name`) and `json_cache_peptide_matches`. The keys
of `class_combinations` and of `json_cache_peptide_matches` stand in the order of their first
transcript, as a `collections.Counter`'s do.

**Standard output** is one line:

```text
[classify_storable_transcript_differences] 310 transcripts: peptide 236, predictions 235, protein_features 265; combinations {"peptide+predictions+protein_features": 191, ...}; JSON cache peptide matches {"113": 236, "both": 74}
```

**Cores and memory.** The two dumps are read on two threads and both are held in memory: a block of
one is compared with the same id's block of the other. The transcripts are classed and the shards of
the JSON cache are read on the worker threads, one shard a thread, and what the shards say is taken
in their order. A shard is held in memory while it is read, and never the whole JSON cache.

**Outside the contract.**

- The order of the keys of `per_class`. CPython counts the classes of one transcript in the order of
  a `set`, which changes from run to run; the tool counts them in the order of their names. A class
  still stands where its first transcript puts it.
- `argparse`'s own usage and error text is not reproduced beyond the exit code 2 and the message
  after `error: `; the usage line names this tool.
- Where CPython ends on a traceback the tool prints one line `ERROR:
  [vep-storable-transcript-classes] <detail>` and exits with 1 as CPython does: a dump that cannot be
  read, a block whose first line has no end, an output directory that cannot be made, a file that
  cannot be written, a shard that is not UTF-8 or not JSON as `json.load` takes it (a nesting past
  1,000 levels and an integer of more than 4,300 digits included) or that is neither a list nor an
  object, a `stable_id` that is a true value and no string, a `variation_effect_feature_cache` that
  is a true value and no object under an id without a translation yet, a shard that cannot be
  opened or gunzipped, a version or a label that is not UTF-8.

## Layout

- `src/lib.rs`: the modules, and the two compiled from `tools/vep-storable-digest-join` by `#[path]`:
  the command line as `argparse` reads it (`argparse`) and text as CPython reads and writes it
  (`py`).
- `src/dumps.rs`: a dump as its blocks, a line's key and a block's translation.
- `src/difflib.rs`: the lines `unified_diff(a, b, n=0)` marks.
- `src/jsoncache.rs`: the shards of a JSON cache and what each says of the wanted ids.
- `src/classify.rs`: the classification, the table, the summary and the line; `src/main.rs`: the
  command line.
- `tests/golden.rs` with `tests/fixtures/inputs/` and `tests/fixtures/expected/`: the binary against
  the reference bytes of 26 runs, each on one thread and on four.
- `tests/fixtures/difflib/cases.txt`: 154 pairs of sequences with the lines CPython's `difflib`
  marks, which the unit test of `src/difflib.rs` holds the matcher to.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest. The package
compiles `src/argparse.rs` and `src/py.rs` of `tools/vep-storable-digest-join` by `#[path]`, so it
builds with that directory beside it, as in a checkout of the repository; a copy of this directory
alone does not build.

```bash
cargo build --release --locked --manifest-path tools/vep-storable-transcript-classes/Cargo.toml
cargo test --locked --manifest-path tools/vep-storable-transcript-classes/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-storable-transcript-classes/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-storable-transcript-classes/Cargo.toml -- --check
```

The reference bytes are the contract, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the files it moves under `tests/fixtures/expected/<case>/` with
the binary's own output, so that the diff of the fixtures shows the change. An output whose name
holds a character outside printable ASCII, a double quote or a backslash is kept as
`<case>/files/<plain name>`, and `<case>/files.tsv` spells the path the run writes it to.
