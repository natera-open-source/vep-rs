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

The plugin reads its rows from the store instead of through tabix on two
conditions:

- the store is named after the data file with `.vpd` in place of `.gz` and sits
  beside it (`whole_genome_SNVs.tsv.vpd` beside `whole_genome_SNVs.tsv.gz`)
- the store was built from that file's content

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

The converter refuses an `--output` that is its input data file or the data
file's `.tbi`, or whose `.vpdi` path is, before it writes anything. It writes
each file whole under another name and renames it into place, the data first
and the index last. A run that starts while the files are written therefore
finds no part-written store.

### When a plugin refuses a store

The store records its flags. A plugin refuses a store whose columns differ from
its data file's, naming the flags that rebuild it. It also refuses a store of
an earlier index version, which carries neither spans nor flags. Every store
built by the converter of release 0.5.0 or earlier is one. A refused store
leaves the data file read through tabix.

The store's index also records the SHA256 and the length of the data file it
was built from and of its own `.vpd` file. A plugin refuses a store whose data
file differs from that record, as a data file replaced under the same name
does. It refuses a store whose `.vpd` differs from it, as an index beside
another store's data does.

The warning names the store, the data file and which of the two differs. The
data file is read through tabix, so the answer is the data file's either way.
An `INFO` line names what each plugin reads: the store, with the data file's
path and both SHA256 values, or the data file.

The data file must stay beside its store: each open reads both whole, once, to
check them. The time therefore grows with the two files' sizes: seconds for
files of several gigabytes where both are in memory.

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

On GRCh37 a plugin's table is read for a mitochondrial record under the
record's own spelling of the contig and its `chr` variant alone. A record on
`MT` reads the rows a table names `MT` or `chrMT`; a record written `chrM`
reads the rows named `chrM` or `M`. The two pairs name different sequences
there: the assembly's NC_012920.1, and NC_001807.4 of the UCSC hg19 naming. The
coding positions of the second stand one base after the assembly's. A table
prepared from files of the UCSC naming therefore annotates no record on `MT`,
in vep-rs as in Ensembl VEP.

AlphaMissense's published table states its `chrM` rows at the assembly's
positions and is read under every spelling. So are the dbNSFP and gnomADMt
tables, whose modules name the contig `M` themselves (`dbNSFP.pm` line 313,
`gnomADMt.pm` lines 259-260). On GRCh38 the four names are one sequence and
every spelling is read
([the documented difference](../divergences/plugin-contig-mitochondrial-spelling.md)).

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
