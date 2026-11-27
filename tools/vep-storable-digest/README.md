# vep-storable-digest

Writes one line per transcript of an Ensembl VEP Storable cache, so that two builds of a cache can be
compared transcript by transcript:

```text
stable_id  contig  start  end  md5_full  md5_content
```

`md5_full` is the MD5 of the transcript's `Data::Dumper` text. `md5_content` is the MD5 of the same
text after what a cache build changes without changing what Ensembl VEP annotates is taken out: the
keys `dbID`, `adaptor`, `created_date`, `modified_date` and `_gene_version` at every depth, the
difference between an integer stored as a number and one stored as a string, and the classes values
are blessed into. `tools/vep-storable-digest-join` joins two such files.

```bash
vep-storable-digest <cache_dir> <out.tsv>
```

`cache_dir` is the directory of one cache (`.../homo_sapiens/116_GRCh38`). The tool reads the
Storable images itself and needs no Perl. `VEP_STORABLE_DIGEST_THREADS` sets the number of worker
threads (default: the host's cores).

## Contract

The output is fixed to the byte: the file, every line on standard error and the exit code.
`tests/fixtures/expected/` holds the reference bytes of 24 runs, made by perl 5.34.0 with Storable
3.23, Data::Dumper 2.179, Digest::MD5 2.58 and IO::Uncompress::Gunzip 2.102, and `cargo test` holds
the binary to them. Where a byte could be written in more than one way it follows the rule of the
Perl module that defines it, named here by its Perl name.

**Command line.** Two arguments are needed; with fewer the tool prints `usage: <program> <cache_dir>
<out.tsv>` on standard error and exits 255. Arguments past the second are not read.

**The walk.** The output file is opened first; a file that cannot be opened is `cannot write <out>:
<the system's reason>` and the system's error number as the exit code. Only the chromosome
directories `1` to `22`, `X`, `Y` and `MT` are read, in the order their names sort as strings (`1`,
`10`, `11`, ..., `2`, `20`, ..., `MT`, `X`, `Y`); one that is absent is passed over. In a directory
the files named `<digits>-<digits>.gz` are the transcript shards, read in the byte order of their
names; `_var.gz` and `_reg.gz` files and every other name are not read.

**A shard file** is read as `IO::Uncompress::Gunzip::gunzip($file => \$buffer)` reads it, then as
Storable's `fd_retrieve` reads the result: the first gzip member alone, whatever follows it; no
check of the member's checksum and length; a file that does not open with the gzip magic as its own
bytes. The image is a hash of lists of transcripts (one list a chromosome; a list that is undef or
another false value is empty) or a list of transcripts. An image that ends early holds no
transcript and still counts as a shard file.

**A transcript** is an item that is a reference (one blessed into the class `0` is not taken). Its
id is `stable_id`, or `_stable_id` where that is undef, without a final `.<digits>`; an item with
neither is passed over. `start` and `end` are written as Perl writes the scalar (an integer as its
digits, a double with C's `%.15g`, undef as the empty string). The line is written by Perl's
`print`: a string with Perl's UTF-8 flag goes out one byte a character when every character of the
line is below 256, else as UTF-8.

**A transcript met again.** A transcript that spans a shard boundary is in two shards. The first
copy is written; a later copy of the same id is counted as a duplicate and, when its `md5_content`
is not the first copy's, as a conflict with the line `WARN: <id> differs between two shards of
<cache_dir>` on standard error.

**The end.** One line on standard error:

```text
[compare_storable_cache_versions] <cache_dir>: <n> shard files, <n> transcript records, <n> distinct transcripts, <n> shard-boundary duplicates, <n> conflicting duplicates
```

The exit code is 1 when a conflict was met, else 0.

**`md5_full`** is `Digest::MD5::md5_hex` of what `Data::Dumper::Dumper` returns for the transcript
under `$Data::Dumper::Sortkeys = 1`, `$Data::Dumper::Indent = 0` and `$Data::Dumper::Terse = 1`, by
the module's compiled implementation:

- a hash is `{'key' => value,...}` with its keys in the order of Perl's `cmp`, an array
  `[item,...]`, with no white space but the ` => `;
- a blessed value is `bless( <value>, '<class>' )`; a compiled pattern blessed into `Regexp` is
  `qr/<source>/<flags>`;
- a scalar the cache stores as an integer is its digits, in quotes when they are more than ten
  characters; one it stores as a string is in single quotes with `'` and `\` escaped, also when it
  is all digits; a double is the quoted text Perl gives it; undef is `undef`; a string with Perl's
  UTF-8 flag that holds a character above 127 is in double quotes with those characters as
  `\x{<hex>}` and `"`, `\`, `$` and `@` escaped; a version string is its literal;
- a reference to a scalar is `\<scalar>`, a blessed one `do{\(my $o = <scalar>)}`;
- a value met a second time is the path of its first place (`$VAR1->{'_trans_exon_array'}[0]`,
  `$VAR1` for the transcript itself): a hash or an array by the path of the reference that first
  led to it, a scalar a reference also points at as `\<path>` or `${<path>}`;
- a weak reference whose target nothing else holds is `undef`, as it is to Perl once the image is
  read.

**`md5_content`** is the same digest of a copy of the transcript. In the copy a hash or an array,
blessed or not, is an unblessed one; the five keys are left out of every hash; a scalar whose text
is an optional `-` and ASCII digits is the number `<text> + 0` (an integer while it fits 64 bits,
signed or unsigned; past that a double, which is written as a quoted string); a reference to a
scalar is a reference to a copy of that scalar as it stands; a hash or an array reached twice is
copied once, so the copy shares what the transcript shares; and a reference to anything else (a
reference, a compiled pattern, a version string, a value blessed into the class `0`) stays the
reference it is, to the value it points at, class and keys included.

**Order and cores.** Shards are read and digested on the worker threads, one shard a thread, and
their lines are written by one thread in the order of the walk, so the file and the lines of a run
are the same at any thread count. A worker starts no shard more than twice the thread count ahead
of the writer: a whole shard is held in memory for each shard in flight, and never a whole cache.

**Where the tool does not write what Perl writes.** A reference run that meets one of these cases
holds the tool's bytes for it.

- A line that Perl ends with the place of a `die` (` at <file> line <n>.`, also chained) is printed
  without it: `Not a HASH reference`, `File is not a perl storable`.
- A run that dies exits 255. Perl exits with the error number its last failed system call left,
  which depends on what the run did before (2 once a chromosome directory was absent).
- The reason after `gunzip <file> failed: ` is the tool's own wording (`the deflate stream is
  damaged or short`), not the module's.
- Perl's own warning lines (`Wide character in print`) are not printed.
- A reference where an id, a start or an end is expected: Perl writes its address in memory, which
  no other run repeats; the tool stops and says so.
- A shard whose hash holds several lists of transcripts: Perl walks a hash in an order of its own
  that changes from run to run; the tool walks it in the order the image stores it.
- `\d` in the id's version and in the test for an integer is an ASCII digit, also on a string with
  Perl's UTF-8 flag, where Perl takes any Unicode digit.
- An integer past 64 bits becomes the nearest double; Perl's own conversion can differ from it in
  the last place.
- A cache directory whose path holds white space or one of `*?[]{}~\` is read as the directory it
  names; Perl's `glob` reads such a path as a pattern.
- A write to the output file that fails ends the run with `cannot write <out>: <reason>`, where
  Perl's `print` returns false and ends nothing.
- Images Storable 3.23 refuses (another native layout, a later binary major, the boolean opcodes of
  binary minor 12) are refused with its message; tied values, hooks, code references and an opcode
  the format does not define are refused with the reader's own message, which names the opcode and
  its offset.

## Layout

- `src/lib.rs`: the modules, and the ones compiled from `tools/vep-storable-to-json` by `#[path]`:
  the reader of Storable images (`storable`), what Perl does with a retrieved value (`perl`,
  `number`, `json`), Storable's outcomes for an image (`shard`) and the worker threads (`pool`).
- `src/shards.rs`: the shards of a chromosome directory, `gunzip` and `fd_retrieve` on one, the
  list of its transcripts, and the weak references Perl clears.
- `src/dumper.rs`: `Data::Dumper`'s text for a value of an image or of an overlay beside it.
  `tools/vep-storable-dump` compiles this file and `src/shards.rs`, `src/text.rs` by `#[path]`.
- `src/digest.rs`: the copy behind `md5_content` and the line of a transcript; `src/md5.rs`: MD5;
  `src/text.rs`: Perl strings joined and printed; `src/run.rs`: the run; `src/main.rs`: the command
  line.
- `tests/dumper.rs` with `tests/fixtures/values/`: the dumper against `Data::Dumper`'s own text for
  16 images, at both indentations.
- `tests/crafted_images.rs`: the dumper and the content copy on images built in the test, nested to
  `Data::Dumper`'s recursion limit of 1,000 levels and past it, on a stack of a worker thread's
  size.
- `tests/golden.rs` with `tests/fixtures/caches/` and `tests/fixtures/expected/`: the binary against
  the reference bytes of 24 runs, each on one thread and on four.
  `tests/fixtures/caches/README.md` says what each cache holds.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest. The package
compiles files of `tools/vep-storable-to-json` by `#[path]` (`src/lib.rs` names them), so it builds
with that directory beside it, as in a checkout of the repository; a copy of this directory alone
does not build.

```bash
cargo build --release --locked --manifest-path tools/vep-storable-digest/Cargo.toml
cargo test --locked --manifest-path tools/vep-storable-digest/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-storable-digest/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-storable-digest/Cargo.toml -- --check
```

The reference bytes are the contract, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the files it moves under `tests/fixtures/expected/<case>/` with
the binary's own output, so that the diff of the fixtures shows the change.
