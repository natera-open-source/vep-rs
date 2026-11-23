# Reproductions of the intended divergences

One file per entry of the category tables of `docs/intended-divergences.md`, named `<entry>.toml`. Each records how to run the entry's example on Ensembl VEP and on vep-rs, where the two outputs differ and what each engine printed when the entry was last verified. `tools/vep-divergence-reevaluate` runs them all and says, per entry, whether the divergence is `still-present`, `closed` (the engines agree) or `changed` (they differ in another way), or `error` when the reproduction could not run, so that every entry is re-judged whenever vep-rs moves to a newer Ensembl VEP reference. `vep-divergence-reevaluate --check` confirms that every entry has a file, that every file is verified against the release's reference, that each names its entry's input class and that each names, in its `implementation` table, code and tests the repository holds.

## Shape

```toml
entry = "minimal_no_intergenic"                      # the file's own name
page = "docs/divergences/minimal-no-intergenic.md"   # where the entry is explained
verified_against = "116.3"                           # the Ensembl VEP release last run
image = "ensemblorg/ensembl-vep@sha256:..."          # that release's image, pinned by digest
input_class = "observed"                             # observed, unobserved or malformed (The input class)
input_class_basis = "2 · 121.0 per million · ..."    # one line: what gives the entry that class
notes = "..."                                        # optional: what the example shows

prepare = [ "<shell command>", ... ]                 # optional: build a synthetic input

[input]
format = "vcf"              # vcf (default), ensembl, hgvs, spdi, region, id, or guess
text = """                  # inline records; a VCF gets the two header lines unless it brings its own
21	8885975	str_ENST00000624291_8885976	GT	GA,GTC	.	.	.
"""
# or: path = "tests/golden/116/<corpus>/variants.vcf"   (a committed file, or "{stage}/<file>")
# header = ["##INFO=<ID=END,...>"]   extra header lines for an inline VCF
# samples = ["S1", "S2"]             sample columns of an inline VCF (FORMAT is added)
# name = "record.vcf"                the inline file's name

[cache]
species = "homo_sapiens"    # default homo_sapiens
assembly = "GRCh38"
source = "ensembl"          # ensembl (default), merged or refseq: the Ensembl VEP species directory
version = "116"             # default 116
json_cache = "tests/golden/116/GRCh38-pick-summary/json_cache"   # vep-rs's cache, committed
# fasta = "tests/golden/116/<corpus>/reference.fa.gz"   a bgzipped FASTA with .fai and .gzi
# gxf = "{corpus}/chr21_genes.gtf.gz"; gxf_kind = "gtf"   Ensembl VEP reads the file; vep-rs the json_cache
# corpus = "tests/golden/116/<corpus>"     the directory {corpus} names (default: json_cache's parent)
# perl_cache = "{stage}/perl-cache"         a prepared Ensembl VEP cache root instead of --perl-cache-dir
# chromosomes = ["MT"]                      the chromosomes the Ensembl VEP cache must hold (default: the input's)
# writable = true                           mount the Ensembl VEP cache writable (--nearest writes into it)

[run]
flags = ["--minimal", "--no_intergenic"]   # given to both engines after the fixed arguments
# ensembl_vep_flags = [...]; vep_rs_flags = [...]   one engine's extra flags
# format = "default"        default, tab, vcf or json (inferred from compare.kind when absent)
# format_flag = true        pass --tab, --vcf or --json; false when the flags select the format
# assembly_flag = true      pass --assembly <cache.assembly>; false for an entry about a run that names none
# tools = ["samtools"]      programs a plugin shells out to that the image lacks, found in the run's --tool-dir directories
# stats = true              --stats_text in place of --no_stats
# buffer_size = 5000
# stdin = true              pipe the input as -i STDIN
# perl_hash_seed = "0"      PERL_HASH_SEED for Ensembl VEP ("" leaves Perl's own)
# perl_mode = "offline"     database: Ensembl VEP runs --database --db_version <version> against Ensembl's public
#                           database and no cache, the route of the input formats its offline parsers refuse
# tool = "filter"           filter or recode: filter_vep and variant_recoder against vep filter and vep recode,
# args = ["-i", "{input}", "-o", "{out}", ...]   with these arguments in place of the fixed ones

[compare]
kind = "rows"
where = { record = "str_ENST00000624291_8885976" }
columns = ["Location", "Allele", "Feature", "Consequence", "DISTANCE"]

[expected]
ensembl_vep = []
vep_rs = ["21:8885976 | A | ENST00000624291 | downstream_gene_variant | 5000"]
# ensembl_vep_exit = 255; vep_rs_exit = 1   the exit status an engine's refusal of the input produces

# [[also]]                  a further run of the entry (Further runs), its tables written as
# [also.input], [also.run], [also.compare] and [also.expected]

[implementation]
code = ["crates/vep-cli/src/<file>.rs::<function>"]        # the items that hold the entry's behaviour
tests = ["crates/vep-cli/tests/<file>.rs::<test>"]         # the tests that pin vep-rs's row: at least one
# rules = ["tools/<tool>/src/<file>.rs::<function>"]       the comparison rules that recognise the difference
# none = "<one sentence>"   in place of code: why no item of vep-rs is specific to the entry
```

## What implements the entry

The `implementation` table names the smallest items that hold the behaviour the entry records, each as `<path>::<item>`: under `code` the functions, methods or constants that make vep-rs's value, under `tests` the tests that pin it, under `rules` the functions of the comparison tools that recognise the difference. An item is the chain of names from its file's top level, joined by `::`: `function`, `Type::method` for a method of any `impl` block of the type, `Trait::method`, `tests::function`, a constant or a static by its name; a function written inside another function is named by the outer one. An entry whose row vep-rs's ordinary path writes, the difference being Ensembl VEP's alone, gives that reason as one sentence under `none` and lists no `code`. `vep-divergence-reevaluate --check` refuses a record without the table, a path the repository lacks and an item its file does not define, and the admission check refuses a pull request that states it touches no intended divergence while a line it changes lies in a listed item.

## The input class

Every reproduction names the class of the input that reaches its entry and gives one line of basis for it. The class is judged on the triggering input as the entry's own page states it.

| `input_class` | The entry's input | `input_class_basis` |
| --- | --- | --- |
| `observed` | Measured real data carries the shape: the chromosome 21 cuts, a concordance dataset, a release file or a committed corpus built from real data. Every tier 1 or tier 2 entry of the index is `observed`, and so is a tier 3 entry whose own page records a measured figure of the shape in such data | At tier 1 or tier 2 the entry's `Tier and count` cell of the index, as printed; at tier 3 `seen in <the real data its page names>: <the figure its page measures there>` |
| `unobserved` | The entry is tier 3, its page measures no figure of the shape in real data, and the format and the documentation allow the input. An input that is unusual, deprecated, or undocumented but accepted is `unobserved` | `well-formed under <format or option contract>; not seen in the measured datasets` |
| `malformed` | The entry is tier 3, its page measures no figure of the shape in real data, and the input breaks the syntax or a stated constraint of its file format's specification, or the documented contract of the option, plugin data file or cache that reads it: a coordinate that is not a number, bytes that are not UTF-8 where the format requires it, an index that disagrees with its table, a value the documentation excludes | The rule the input breaks, with its source: the clause of the VCF, GFF3, GTF, BED, bigWig, tabix, BGZF or gzip specification, the UTF-8 requirement, or the line of Ensembl VEP's documentation or of a plugin's header that states the contract |

A classification in doubt is `unobserved`: `malformed` is for clear cases alone. A format's version is no part of the shape: a shape that the version the input declares, or a later one, admits is not `malformed`, and its basis names the version that admits it; a basis otherwise cites the declared version where it states the clause. A malformed shape that measured real data carries is `observed`, whatever rule it breaks. An entry that covers several shapes takes the class of its most reachable one, `observed` before `unobserved` before `malformed`, and its basis says so. The measured datasets an `unobserved` basis names are the measured real data of the table above, as the index and the entry's own page measure it; a corpus line written to make a shape is not among them.

The tier and the class answer two questions. The tier counts the index's own data alone, the chromosome 21 cuts and the committed corpora: at tier 3 their ordinary data does not carry the shape. The class says whether any measured real data carries the shape, so an entry whose page records a measured figure of it in a concordance dataset or a release file, or in a corpus record the tier does not count as ordinary data, is `observed` and keeps its tier 3. A page that names such data and gives no figure leaves its entry `unobserved`.

The index lists the `malformed` entries in its last section of entries, `Differences only malformed input reaches`, and every other entry in the category tables above it. `vep-divergence-reevaluate --check` refuses a reproduction that lacks either field, names a class outside the three or gives a basis of more than one line, an `observed` basis that is neither the entry's cell at tier 1 or tier 2 nor a `seen in` line at tier 3, an `unobserved` basis in another form and a `malformed` basis in the form of another class, and holds each class to the entry's tier and to the section its row stands in, as `scripts/divergences/check_divergence_order.py` does: a tier 1 or tier 2 entry is `observed`, a `malformed` entry is tier 3 and stands in the last section, and no other entry stands there.

## The fixed arguments

Ensembl VEP: `vep -i <input> -o <out> --offline --cache --dir_cache <root> --cache_version <version> --species <species> --assembly <assembly> --format <input.format> --buffer_size <n> --fork 1 --force_overwrite --quiet --no_stats [--fasta <plain FASTA>] <flags> [--tab|--vcf|--json]`, in the image, with `PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0`; when the flags name a plugin and no `--dir_plugins`, `--dir_plugins <plugins>` is added from `--plugins-dir`, which such a run needs. A `gxf` cache replaces the cache arguments with `--gtf <file>` or `--gff <file>`; `perl_mode = "database"` replaces them with `--database --db_version <version>`.

vep-rs: `vep -i <input> -o <out> --offline --json_cache <json_cache> --species <species> --assembly <assembly> --format <input.format> --buffer_size <n> --force_overwrite --quiet --no_stats [--fasta <FASTA>] <flags> [--tab|--vcf|--json]`.

The placeholders a flag or argument may carry: `{corpus}` (the committed corpus directory), `{stage}` (the directory the `prepare` steps write), `{plugins}` (the plugin directory: the corpus's `plugin_data` for vep-rs, the plugin modules for Ensembl VEP), `{repo}` (the repository root), `{input}`, `{out}`, and `{<name>}` for each `--extra-dir <name>=<path>` the run was given (a plugin repository outside the plugin modules, such as `{loftee}`). A `run.tools` list names programs the image lacks (LoFTEE's `samtools`); each is found in a `--tool-dir` directory, which Ensembl VEP then has first on its PATH. Ensembl VEP reads each at a read-only mount, the stage writable; vep-rs reads the host path. A reproduction fetches no file at run time: everything it reads is committed, prepared from committed files, or named by one of these options; a `perl_mode = "database"` entry reaches Ensembl's public database, and so does a `tool = "recode"` run, whose `variant_recoder` resolves identifiers through it. A `prepare` step may call `bgzip` and `tabix`, which the host then needs, and the one that rebuilds a bigWig with zoom levels calls `bigtools`. A step that converts a cut of an Ensembl VEP cache builds and runs `tools/vep-storable-to-json` through `cargo run`, so the host needs a Rust toolchain for it.

`prepare` steps run under `bash -euo pipefail` in the stage directory with `STAGE`, `REPO`, `CORPUS`, `JSON_CACHE`, `PERL_CACHE` (the resolved Ensembl VEP cache), `VEP` (the vep-rs binary), `VEP_CACHE_BUILDER` (its sibling) and `IMAGE` (the Ensembl VEP image of the run, empty without `--image`) in the environment. A step whose divergence lies in a tool of the release itself, the `tabix` an image ships, or that needs a run of Ensembl VEP beside the entry's own, runs it through `docker run` on `$IMAGE` and fails with a message when `IMAGE` is empty, so such an entry needs the image. A step writes only under `$STAGE`, but for the build of `tools/vep-storable-to-json`, which `cargo` keeps in that tool's `target` directory; one that hard-links a cache into the stage (`cp -al`) needs `--work-dir` on the cache root's filesystem, and one that rewrites a linked file removes the link first, since a write through it would reach the shared cache. Where a step's status is a pipeline's, the reader reads its pipe to the end: under `pipefail` a reader that leaves early (`head`, an `awk` that exits, a quiet `grep`) lets the closed pipe end the writer, and the step fails by timing with exit code 141. A pipeline that needs such a reader runs in a subshell of its own after `set +o pipefail`, and a test of `tools/vep-divergence-reevaluate` holds every committed step to this.

## What is compared

`compare.kind` names the value read from each engine's run; `expected.ensembl_vep` and `expected.vep_rs` record it as a list of strings (a single string stands for a list of one). The two lists are sorted unless `sorted = false`, trailing whitespace is removed from each value unless `strip = false`, absolute paths are replaced by `<path>` unless `mask_paths = false`, and each regular expression in `mask` is replaced by `<masked>`.

| kind | reads | keys |
| --- | --- | --- |
| `rows` | the default or tab rows matching `where`, each rendered as its `columns` joined by ` \| ` (default Location, Allele, Feature, Consequence) | `where`, `columns` |
| `field` | one column or `Extra` item of the matching rows, `-` when a row lacks it; `with_key = true` prefixes each value with its Feature | `where`, `field`, `with_key` |
| `count` | the number of matching rows | `where` |
| `header` | the output's `#` lines matching `pattern` | `pattern` |
| `lines` | the output's lines matching `pattern` | `pattern` |
| `warnings` | the lines of `<out>_warnings.txt` matching `pattern` | `pattern` |
| `stderr`, `stdout` | the engine's standard streams' lines matching `pattern` | `pattern` |
| `exit` | the exit status | |
| `stats` | the lines of one `[section]` of the `--stats_text` report (`run.stats = true`) | `section` |
| `json` | the values at `path` in the JSON records matching `where`: `a.b[].c` iterates a list, `a[transcript_id=ENST...].c` keeps the members whose key equals the value | `where`, `path` |
| `vcf` | the `field` of the `info` entries (default `CSQ`) of the VCF records matching `where`, or the whole data line with `line = true` | `where`, `field`, `info`, `line` |

`where` keys for `rows`, `field` and `count`: `record` (the first column), `location`, `allele`, `gene`, `feature`, `feature_type`, `consequence` (the set of terms, in any order) and `pattern` (a regular expression over the line); a value may be a list of alternatives. For `vcf` and `json`: `record` or `id`, `chrom`, `pos`, `location` (`chrom:pos`) and `pattern`; `json` also takes `input` (the record's own input line); a `vcf` `field` is further narrowed by `feature` and `allele`. A key the kind does not read, or a `where` on a kind that reads none, is refused.

## Further runs

An entry whose difference shows under more than one flag, where no one run can hold every shape, records each further shape as a table of the `also` list. A table opens with `[[also]]` and gives `[also.input]`, `[also.run]`, `[also.compare]` and `[also.expected]` with the keys of the record's own tables, and may carry a `notes` line. A table it leaves out is the record's own, `expected` excepted, which every further run records for itself. Every run reads the record's cache and runs the record's `prepare` steps in a stage of its own.

`vep-divergence-reevaluate --check` holds each further run to the schema as it holds the record's own run, and names a violation by the run's place in the list, `also[1]` for the first.

## The verdict

`closed` when the two engines' values are equal; `still-present` when they differ and each equals its recorded value; `changed` otherwise. An entry of several runs is `still-present` or `closed` when every run reads so and `changed` when its runs read differently, and a run that cannot be made is the entry's `error`. A reproduction that cannot run (a cache the roots lack, a prepare step that fails, an unknown key) is an `error`, and so is an engine run that fails: for a kind that reads the output, an engine that exits non-zero or writes no output file is a failed run, never a recorded empty value, unless `expected.ensembl_vep_exit` or `expected.vep_rs_exit` records the exit status the entry's divergence produces (an engine that refuses the input), in which case any other exit status is the error; the kinds `exit`, `stderr` and `stdout` read the failure itself. Every entry's file records the release it was last verified against; a re-evaluation against a newer release that finds an entry `closed` retires the entry with that release's number, and one that finds it `changed` re-traces the mechanism before the page and the file are updated.
