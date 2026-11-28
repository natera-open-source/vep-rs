# Preparing data files {#data-file-preparation}

A tabix plugin reads data files you stage first: the script that stages them
per assembly and tier, their format, and the `.vpd` store that replaces tabix.

## The `setup_plugin_data.sh` script

Stage the open plugin data files with `scripts/data/setup_plugin_data.sh`. A
plugin whose data file is not redistributed here reads the file you pass to it, and the
script stages no file for it. It does not
download from the upstream databases itself, and neither does Ensembl VEP's
installer: `INSTALL.pl` fetches a plugin's module and prints that the plugin
requires data, leaving the download to the user, so no plugin data is downloaded
for the user by either engine.

`--source upstream`, the default,
prints, per plugin, the upstream URL and the preparation each file needs (header
fixes, per-assembly tabix indexes, decompression) for you to run. `--source s3`
syncs an already-prepared copy from an S3 bucket of your own, named by
`--bucket` or `$VEP_BENCHMARK_BUCKET`, that you seeded from those steps:

```bash
# print sources for chr21 subsets
scripts/data/setup_plugin_data.sh --assembly GRCh38 --tier quick \
    --source upstream

# sync genome-wide files from your bucket
scripts/data/setup_plugin_data.sh --assembly GRCh38 --tier full \
    --source s3 --bucket my-bucket
```

`--assembly` is required: it selects which files are staged and which plugins
apply. An unsupported combination is skipped with its reason rather than
substituted from the wrong assembly.

Files land under `~/.vep/plugin_data/<assembly>/<tier>/` by default, in
per-plugin subdirectories (`alphamissense/alphamissense.tsv.gz`,
`revel/revel.tsv.gz`, ...). That is the layout
`scripts/data/plugin_data_layout.sh` defines for the data-preparation and
concordance scripts alike.
[Reproducing the plugin concordance](reproducing-concordance.md) lists the
file each plugin's concordance run reads from that layout.

## Tabix indexes and `.vpd` stores

All tabix plugins require bgzipped, tabix-indexed data files (`.tsv.gz` +
`.tsv.gz.tbi` or `.vcf.gz` + `.vcf.gz.tbi`). The one exception is a `.vpd`
binary store built from the data file by `vep-cache-converter
convert-plugin-data`, with its `.vpdi` index beside it.

Fifteen plugins read a store: AlphaMissense, AVI, BayesDel, CADD, ClinPred,
dbNSFP, dbscSNV, gnomADc, gnomADMt, OpenTargets, PrimateAI, REVEL, SpliceAI,
SpliceVault and TranscriptAnnotator. LOEUF reads its data file through tabix
alone, since it takes its coordinate columns from the file's tabix index.

Such a plugin reads its rows from the store instead of through tabix on two
conditions:

- the store is named after the data file with `.vpd` in place of `.gz` or
  `.bgz` and sits beside it (`whole_genome_SNVs.tsv.vpd` beside
  `whole_genome_SNVs.tsv.gz`)
- the store was built from that file's content

CADD takes the columns it reads from the data file's tabix index. A data file
with no index is read from its store alone, and CADD then takes the columns
from the names in the file's header. It reads a `Pos` column as one-based
positions and a `Start` column as zero-based starts, as the plugin's own
comparison reads them (`CADD.pm` lines 209 and 220). An `End` column, where the
header names one, ends each row.

A store built on other columns is refused there, and the refusal names the
flags that rebuild it. The converter checks a store through the data file's
index, so a data file with no index is converted with `--no-verify`.

### The layout flags of a store

A store spans each row as the plugin's tabix reader does. It answers a query
with every row reaching into it, the rows at one start in the data file's
order, as `tabix` prints them. The converter's layout flags name the columns
the plugin reads.

Every data file but two is read as single positions. The default flags fit a
data file whose first two columns the plugin reads. A data file holding both
assemblies' coordinates takes the columns of the run's assembly:

| Data file | Assembly | Flags |
| --- | --- | --- |
| REVEL | GRCh38 | `--start-col 2` |
| dbscSNV (`dbscSNV.pm` lines 53-55) | GRCh38 | `--chr-col 4 --start-col 5` |
| dbNSFP 4 (`dbNSFP.pm` lines 72-74) | GRCh37 | `--chr-col 7 --start-col 8` |

gnomADMt's VCF is read with `--end-col 7 --ref-col 3`, its INFO `END` and
reference allele. A CADD structural-variant data file whose index declares
zero-based starts and an end column is read with `--end-col 2 --zero-based`.

The converter holds its layout flags to the data file's tabix index before the
data file is read. It refuses `--chr-col` and `--start-col` naming other columns
than the index was built on. Beyond those two, the flags name either the end
column and the starts the index declares or single one-based positions: any
other `--end-col` given without `--ref-col` is refused. So is any other
`--zero-based`, unless it is given with `--end-col` and `--ref-col` under an
index that declares zero-based starts. An
index of tabix's `vcf` or `sam` preset declares one-based starts and no end
column. The refusal names the flags that match the index.

Where the flags name single one-based positions and the index declares an end
column or zero-based starts, the converter prints a warning naming the flags
that match the index and converts: CADD refuses that store. A data file without
an index is not compared.

### What the converter checks

The converter checks its store against the data file through the plugin's tabix
reader and leaves no store that fails the check. It checks every row start
under `--verify-all`, or when the data file holds at most 1,000 of them.
Otherwise it checks every n-th start, the same starts on every run.

At each start it compares the two readers at the start and at each row's end,
its midpoint and the positions on either side of its span. Each position is
asked through a single query and through batch queries: alone, then with the
positions beside it, in order and reversed.

A cell, column name or chromosome name over 65,535 bytes, or a header or row of
more than 65,535 cells, is refused naming the row, since a store's length
fields are 16 bits wide. The data file is then read through tabix.

The converter refuses an `--output` that is its input data file or a path the
data file's index is read from, or whose `.vpdi` path is, before it reads the
data file. The index paths are `<file>.tbi`, `<file>.csi` and the same two
names on the file's stem. A column index of 65,535 or more is refused at the
same point. An `--output` that is a symbolic link is written at the file the
link's chain names, and the link stays. A link that names a directory, and a
data path and `.vpdi` path that name one file through a link, are refused.

The converter reads the data file once and writes each row as it reads it, so
it keeps no row in memory while it reads: for each chromosome it holds the
name and a few counts. While it runs it needs room beside the output for the
store and the rows of one chromosome more, in scratch files it removes.

The converter writes each file whole under a scratch name beside its own and
renames it into place, the data first and the index last. A run that starts
while the files are written therefore finds no part-written store. A scratch
name ends in the converter's process id, a role and `.partial`. The two files
have the roles `data` and `index`. What the converter keeps of each chromosome
while it reads has the roles `rows`, `starts` and `spans`, each followed by
the chromosome's number.

An earlier store's index is removed once both scratch files are written
whole, before they are renamed. A conversion that is refused or fails before
that point, while it reads the data file included, leaves an earlier store
whole. One that fails at a rename, or whose store fails the check, leaves no
store.

The converter holds a lock on the `data` and `index` scratch files while it
writes, and the numbered scratch files are held with them. A conversion that
is killed can leave its scratch files beside the output. The next conversion
of the same `--output` removes those that no running conversion holds. On a
file system that grants no file lock they stay.

### When a plugin refuses a store

The store records its flags. A plugin refuses a store whose columns differ from
its data file's, naming the flags that rebuild it. It also refuses a store
whose index is of any version but the one it reads, an earlier version as a
later one, naming the conversion that rebuilds it. Every store built by the
converter of release 0.6.0 or earlier is one.

A refused store leaves the data file read through tabix, which needs the data
file's `.tbi` index beside it. Without that index the plugin is not loaded.
For a store refused for its index version, the plugin's refusal names the
missing `.tbi`, the store, the version its index carries and the conversion
that rebuilds it. Beside a data file with no index, a refusal that names the
conversion ends with `--no-verify`, the flag the conversion needs on such a
file.

The store's index also records the length, a chunked SHA256 and a sampled
digest of the data file it was built from and of its own `.vpd` file. The
chunked SHA256 is the SHA256 of the SHA256 of each 64 MiB of the file, taken
side by side, so it does not equal the file's `sha256sum`. The sampled digest
is the SHA256 of a file's first and last MiB. Every open compares the lengths
and the sampled digests. A plugin refuses a store whose data file differs from
that record, as a data file replaced under the same name does. It refuses a
store whose `.vpd` differs from it, as an index beside another store's data
does.

An open reads only those two MiB of each file. A change that keeps a file's
length and lies wholly between its first and last MiB is not seen by an open,
and the plugin reads the store as it stands. The converter's check reads both
files whole and refuses such a change by their chunked SHA256:

```bash
vep-cache-converter verify --store <file>.vpd --table <data file>
```

Run it after copying or restoring a store, or whenever a data file may have
been edited in place. A store whose index is of version 5, which records no
sampled digest, is held to the whole files' SHA256 on every open.

The index ends with its own length, the chunked SHA256 of its bytes and a sampled
digest of its first and last MiB. Every open checks the length and the sampled
digest, so a plugin refuses an index cut short, or altered in those MiB. An
index altered at the same length between them opens; `verify` checks the
chunked SHA256 of the whole index too and refuses it. A version 5 index
records the plain SHA256 of each file and of itself, and is held to them on
every open.

The converter refuses to write positions out of ascending order, and `verify`
refuses an index that holds them. An open of a version 6 index does not read
its position entries to check their order. An index altered between its
sampled MiB so that its positions no longer ascend therefore opens. Its
queries read only within the store's files, but they can return the wrong
rows or none until `verify` refuses the store. A version 5 index is
checked for its order on every open. It also
refuses an index that names one chromosome twice, or that holds a chromosome's
positions out of ascending order.

The warning names the store, the data file and what differs. The
data file is read through tabix, so the answer is the data file's either way.
An `INFO` line names what each plugin reads: the store, with the data file's
path and both whole-file digests, or the data file.

The data file must stay beside its store: each open reads the first and last
MiB of the data file, the store's `.vpd` and its index to check them, which
takes milliseconds whatever the files' sizes. `verify`, and an open of a
version 5 index, read all three files whole, which takes seconds for files of
several gigabytes where they are in memory.

## How a row's cells are read

A row is read with its line end alone removed, as Ensembl VEP's tabix iterator
hands it to a plugin. A plugin that prints a cell prints it as the row holds
it, the whitespace before the line end included. A whitespace-only last cell is
a cell.

A plugin whose cells are numbers reads them with the whitespace removed. The
JSON format writes a numeric cell as the number without its surrounding
whitespace ([JSON output](../output/json.md)). A byte that is not UTF-8 in a
row or a header prints as U+FFFD, through tabix and through a store alike
([Input bytes that are not UTF-8](../divergences/non-utf8-input.md)).

An empty last cell is a cell the row lacks for a plugin whose Ensembl VEP module
splits a row into a list of cells. That is how `split /\t/` leaves it.
PrimateAI reports such a score row as missing its values (`PrimateAI.pm` line
198). dbscSNV writes an empty or absent `ada_score` or `rf_score` as the empty
value, `rf_score=` in the default format. It skips `.` alone (`dbscSNV.pm`
lines 179 to 184).

PrimateAI writes an empty score cell that has a filled cell after it as the
empty string, an empty tab cell and `""` in JSON, since the module's split
keeps that cell (`PrimateAI.pm` lines 188, 196 and 204 to 207).

ClinPred's and BayesDel's modules split a row into five variables (`ClinPred.pm`
line 127, `BayesDel.pm` line 151), which keeps a trailing empty cell. A row
lacking its score cell is written as the empty value: `ClinPred=` in the
default format, `-` in the tab format, no JSON key. A present empty score cell
is written as the empty string: an empty tab cell, `""` in JSON. BayesDel's
score prints on every output row of the record, so its intergenic row carries
the same cell.

AVI (`AVI.pm` line 161) and SpliceVault (`SpliceVault.pm` line 201) split the
same way. A present empty cell is the empty string, and a cell the row lacks is
the empty value (`null` inside SpliceVault's JSON object). AVI skips a row
lacking its sixth cell. SpliceVault's events column is a list whether its cell
is empty or missing, so `SpliceVault_top_events` is then an empty tab cell and
`{}` in JSON.

## The mitochondrial contig of a table on GRCh37 {#mitochondrial-contig-grch37}

On GRCh37 the two pairs of mitochondrial names are two sequences: `MT` and
`chrMT` name the assembly's NC_012920.1, `chrM` and `M` the NC_001807.4 of the
UCSC hg19 naming, whose coding positions stand one base after the assembly's.
A plugin's table is read for a mitochondrial record under the record's own
spelling of the contig and its `chr` variant: a record on `MT` reads the rows a
table names `MT` or `chrMT`, a record written `chrM` those named `chrM` or `M`.

A contig of the other pair is read only where the table's rows state the bases
the run's reference sequence, `--fasta`, holds at their own positions. vep-rs
compares a bounded sample of the rows once a table and logs the verdict
([the check](../divergences/plugin-contig-mitochondrial-spelling.md)). A table
prepared from files of the UCSC naming therefore annotates no record on `MT`.
Without `--fasta`, or for a table with no reference-allele column, the
comparison cannot be made: the record's own pair alone is read and one warning
names the table. On GRCh38 the four names are one sequence and every spelling
is read.

## The cells of a TranscriptAnnotator table {#transcriptannotator-cells}

A cell a row has and leaves empty, or that `trim=1` reduces to nothing, is
written as an empty value: an empty cell in the tab format and `""` in JSON. A
column the row lacks is written as an absent one: `-` in the tab format and no
JSON key. The default and VCF formats write both as the empty value.

Both follow Ensembl VEP, whose `split /\t/` drops a row's trailing empty cells
before trimming (`TranscriptAnnotator.pm` lines 245 to 250). It drops a header
line's trailing empty name the same way.

`trim=1` removes the six ASCII whitespace characters, as Ensembl VEP's
`_trim_whitespaces` does (`TranscriptAnnotator.pm` line 91), so a no-break space
stays part of its cell. A `.vpd` store is built from rows read the same way,
the line end alone removed. A whitespace-only trailing cell is therefore a cell
through a store as through tabix.

## Where next

The download and indexing steps for one plugin are on its own page:

- [annotation-store plugins](annotation-store.md): CADD, REVEL and the other
  plugins that read a tabix-indexed data file or its store
- [LoFTEE](loftee.md): the two upstream data bundles
- [flat-file plugins](flat-file.md): LoFtool, pLI, GWAS and DosageSensitivity
- [SpliceAI](spliceai.md): the SpliceAI score VCFs
