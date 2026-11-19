# LoFTEE on GRCh37

What a GRCh37 run of LoFTEE passes, what vep-rs refuses there, and where the
`master` branch of Ensembl VEP's LoF plugin reads its data differently.

On GRCh37 Ensembl VEP runs the `master` branch of the LoF plugin, and vep-rs
reads the keys that branch reads; [LoFTEE](loftee.md) lists the parameters and
the assembly each is read on.

## The key a GRCh37 run passes

The `master` branch runs its splice-prediction layer unless
`run_splice_predictions` holds `0` (`LoF.pm` lines 98 and 186); the `grch38`
branch leaves the layer off by default. vep-rs does not implement that layer,
so on GRCh37 it refuses a run that leaves the layer on: the key's last value
must be `0`. A run without the key, or with an empty value, is refused at
plugin initialization under one line, wrapped here:

```text
LoFTEE: this GRCh37 run asks for the splice-prediction layer, which
vep-rs does not implement: Ensembl VEP's LoF plugin runs it on that
assembly unless run_splice_predictions is 0; pass
run_splice_predictions:0
```

A GRCh37 run therefore names the key. `$LOFTEE` is the unpacked bundle for the
run's assembly, `/path/to/loftee/grch37`:

```bash
--plugin LoFTEE,human_ancestor_fa=$LOFTEE/human_ancestor.fa,\
gerp_tabix=$LOFTEE/GERP_scores.final.sorted.txt.gz,\
run_splice_predictions=0
```

An Ensembl VEP `LoF` argument string carrying `run_splice_predictions:0` is
read the same way. The assembly is the one `--assembly` names or, without the
flag, the cache's. On GRCh38 a command line without the key starts: the
`grch38` branch's default is the layer off.

LoFTEE reads the names `GRCh37`, `hg19` and `37`, in any letter case, as
GRCh37. Under any other name for that assembly, such as a cache that calls
itself `GRCh37.p13`, the run warns that the name is not one LoFTEE recognises,
reads the keys and defaults of the `grch38` branch, and is not refused.

## What a refused plugin does to the run

LoFTEE refuses, at plugin initialization, what it does not compute: a GRCh37
run that leaves the splice-prediction layer on, a true `run_splice_predictions`
on either assembly, and a `conservation_file` naming a database. It refuses a
data file it cannot open the same way. Without `--safe` the refused plugin is
dropped with a warning, the run exits 0, and no row carries a `LoF` field.
Under `--safe` the run stops with exit status 1 and writes no output. Any
plugin that fails to start is handled the same way.

## What the `master` branch does besides

Each row is a behaviour of the `master` branch that vep-rs does not have.

| Subject | The `master` branch | vep-rs |
| --- | --- | --- |
| The default of `human_ancestor_fa` | `human_ancestor.fa.rz`, a file in the working directory, read with `samtools faidx` (`LoF.pm` lines 77 and 553). The `ANC_ALLELE` test runs wherever that file and the program are found. | Off. The test runs only under a `human_ancestor_fa` naming an indexed FASTA. |
| A `gerp_bigwig` | Never read. With no `gerp_file`, the branch reads its default, `/vep/loftee/GERP_scores.final.sorted.txt.gz` (`LoF.pm` line 132). | The stand-in rule holds on both assemblies: with no `gerp_file`, the file named under `gerp_bigwig`, `gerp_tabix`, `gerp_scores` or `gerp` is read. |
| `gerp_file:mysql` | Opens a database connection and sums GERP from its `gerp_bases` table (`LoF.pm` lines 134-136; `gerp_dist.pl` line 82). | Refused at plugin initialization: the word is read as a path, and no such file exists. |
| Whole-exon GERP totals | Each whole exon the GERP walk crosses takes its total from the `gerp_exons` table of the `conservation_file` database (`gerp_dist.pl` lines 59-60 and 129). | Each such exon's total is the sum of the GERP file over the exon. A `conservation_file` naming a database is refused, since the same database asks for the PhyloCSF items. |
| A failed GERP read | A `tabix` error line counts 0 toward `GERP_DIST`, its third field being summed with the scores (`gerp_dist.pl` lines 92-96), and the row is written: a defect, the message being no score. | The row fails: it keeps no `LoF` field, under the plugin's `went wrong` warning. |

A GRCh37 run therefore passes its GERP file under `gerp_tabix` or `gerp_file`
and no `conservation_file`, or `conservation_file:false`: the whole-exon totals
come from the GERP file.

The sixth is the row the `master` branch leaves unwritten. Without a
`conservation_file` it writes no `LoF` field on a stop-gained or frameshift row
whose GERP walk crosses a whole exon, and vep-rs writes the row: the class
`loftee_gerp_walk_without_conservation_database`, whose
[page](../divergences/loftee-gerp-walk-without-the-conservation-database.md)
gives the record and both engines' rows.
