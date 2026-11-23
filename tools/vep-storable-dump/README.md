# vep-storable-dump

Writes named transcripts of an Ensembl VEP Storable cache as `Data::Dumper` text, one block per
transcript, for a field-level comparison of two cache builds. Each prediction matrix is summarised
as its length and its translation digest, so that the text stays readable.
`tools/vep-storable-transcript-classes` reads two such dumps.

```bash
vep-storable-dump <cache_dir> <ids.txt> <out.txt>
```

`cache_dir` is the directory of one cache (`.../homo_sapiens/116_GRCh38`) and `ids.txt` holds one
transcript id a line, without a version. The tool reads the Storable images itself and needs no
Perl. `VEP_STORABLE_DUMP_THREADS` sets the number of worker threads (default: the host's cores).

## Contract

The output is fixed to the byte: the file, every line on standard error and the exit code.
`tests/fixtures/expected/` holds the reference bytes of 32 runs, made by perl 5.34.0 with Storable
3.23, Data::Dumper 2.179 and IO::Uncompress::Gunzip 2.102, and `cargo test` holds the binary to
them. Where a byte could be written in more than one way it follows the rule of the Perl module
that defines it, named here by its Perl name.

**Command line.** Three arguments are needed; with fewer the tool prints `usage: <program>
<cache_dir> <ids.txt> <out.txt>` on standard error and exits 255. Arguments past the third are not
read.

**The ids.** The list is read first; one that cannot be opened is `cannot read <ids.txt>: <the
system's reason>` and the system's error number as the exit code. A line names an id when it opens
with `ENS`; its line feed and a final `.<digits>` are taken off and the rest is the id, whatever it
holds: a carriage return or further columns stay part of it and name no transcript. Every other
line is passed over.

**The walk.** The output file is opened next (`cannot write <out>: <reason>` otherwise). The
chromosome directories `1` to `22`, `X`, `Y` and `MT` are read in that order; in each the files
named `<digits>-<digits>.gz` are the transcript shards, in the byte order of their names. A shard
file is read as `IO::Uncompress::Gunzip::gunzip` and Storable's `fd_retrieve` read it, and its
transcripts are listed, as `tools/vep-storable-digest` does (its README's "A shard file").

**A transcript** is dumped when its `stable_id` without a final `.<digits>` is an id of the list
and no earlier transcript of the walk had that id. An item that is undef, one without a
`stable_id` and one whose `stable_id` is a reference are passed over; an item that is a string or
a reference to anything but a hash ends the run with Perl's message for it.

**A block** is `### <id>`, a line feed, the text, and an empty line. The text is what
`Data::Dumper::Dumper` returns under `$Data::Dumper::Sortkeys = 1`, `$Data::Dumper::Indent = 1`
and `$Data::Dumper::Terse = 1` for a copy of the transcript:

- the copy is an unblessed hash of the transcript's own values, so the text opens with `{`;
- where `_variation_effect_feature_cache` is an unblessed hash it is copied the same way, and
  where that copy's `protein_function_predictions` is an unblessed hash, every value of it is
  replaced by a string: `matrix:<length> md5:<translation_md5>` for a reference (the length of
  `matrix` in characters, `undef` where it is undef; the digest, empty where it is undef) and
  `none` for any other value. A matrix that is a reference to anything but a hash ends the run;
- every other value is the transcript's own, at any depth: a hash is `{`, its keys in the order of
  Perl's `cmp` one a line, `}`; an array `[`, its items one a line, `]`; each level two spaces
  deeper; a blessed value `bless( <value>, '<class>' )`; integers bare, strings in quotes; a value
  met a second time is the path of its first place. `tools/vep-storable-digest`'s README gives
  the rules in full under `md5_full`.

A transcript's translation holds a reference to the transcript itself. The copy is another hash,
so that reference is dumped as the whole transcript once more, matrices included, under
`'translation' => ... 'transcript' => bless( {`.

**The end.** One line on standard error, `[dump_storable_transcripts] <cache_dir>: <n> of <m>
requested transcripts dumped`, and exit code 0.

**Order and cores.** Shards are read and dumped on the worker threads, one shard a thread, and
their blocks are written by one thread in the order of the walk, so the file is the same at any
thread count. A worker starts no shard more than twice the thread count ahead of the writer: a
whole shard is held in memory for each shard in flight, and never a whole cache.

**Where the tool does not write what Perl writes.** A reference run that meets one of these cases
holds the tool's bytes for it.

- A line that Perl ends with the place of a `die` (` at <file> line <n>.`, also chained) is printed
  without it.
- A run that dies exits 255. Perl exits with the error number its last failed system call left,
  which depends on what the run did before (2 once a chromosome directory was absent).
- A matrix or a translation digest that is itself a reference: Perl writes its address in memory,
  or the length of that address, which no other run repeats; the tool stops and says so.
- A shard whose hash holds several lists of transcripts: Perl walks a hash in an order of its own
  that changes from run to run; the tool walks it in the order the image stores it.
- `\d` in an id's version is an ASCII digit, also on a string with Perl's UTF-8 flag, where Perl
  takes any Unicode digit.
- A cache directory whose path holds white space or one of `*?[]{}~\` is read as the directory it
  names; Perl's `glob` reads such a path as a pattern.
- A write to the output file that fails ends the run with `cannot write <out>: <reason>`, where
  Perl's `print` returns false and ends nothing.
- Images Storable 3.23 refuses are refused with its message, and what the reader does not decode
  with the reader's own, as in `tools/vep-storable-digest`.

## Layout

- `src/lib.rs`: the modules, all but one compiled from other tools by `#[path]`: the reader of
  Storable images, Perl's handling of a retrieved value, Storable's outcomes and the worker threads
  from `tools/vep-storable-to-json`; the shards of a cache, `Data::Dumper`'s text and Perl strings
  from `tools/vep-storable-digest`.
- `src/run.rs`: the list of ids, the copy of a transcript that is dumped, and the run;
  `src/main.rs`: the command line.
- `tests/golden.rs` with `tests/fixtures/`: the binary against the reference bytes of 32 runs, each
  on one thread and on four, over the caches of `tools/vep-storable-digest/tests/fixtures/caches/`
  and the three of `tests/fixtures/caches/` (`README.md` there says what they hold).

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest. The package
compiles files of `tools/vep-storable-to-json` and of `tools/vep-storable-digest` by `#[path]`
(`src/lib.rs` names them) and its tests read the caches of the second, so it builds and tests with
those two directories beside it, as in a checkout of the repository; a copy of this directory alone
does not build.

```bash
cargo build --release --locked --manifest-path tools/vep-storable-dump/Cargo.toml
cargo test --locked --manifest-path tools/vep-storable-dump/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-storable-dump/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-storable-dump/Cargo.toml -- --check
```

The reference bytes are the contract, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the files it moves under `tests/fixtures/expected/<case>/` with
the binary's own output, so that the diff of the fixtures shows the change.
