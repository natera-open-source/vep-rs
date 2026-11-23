# Mitochondrial frame check fixture

Every file here is made up. On GRCh37 `MT` and `chrMT` name the assembly's
mitochondrial sequence and `chrM` and `M` another one, so a plugin's table that
names its contig `M` is read for a record on `MT` only where its rows state the
bases the run's reference sequence holds at their own positions. The fixture
holds a reference, two records and three tables that differ in the name of the
contig and in the bases their rows state.

`reference_synthetic.fa` is one sequence named `MT` of 16,569 bases, each drawn
with `random.Random(37).choice("ACGT")` of Python's standard library, 60 to a
line, with the index `samtools faidx` writes. It copies no published sequence.

`variants.vcf` holds two single-base records on `MT` inside a protein-coding
transcript of the cache slice `../bayesdel-mitochondrial-contig/json_cache`,
each with the reference's base as REF and the next base in the order `A`, `C`,
`G`, `T` as ALT:

- `same_next` at `MT:3556`, where the reference holds the same base one
  position on;
- `other_next` at `MT:3557`, where it holds another base one position on.

Three BayesDel-shaped score files of five columns under the header
`#Chr Pos Ref Alt BayesDel_addAF`, bgzip-compressed and indexed with
`tabix -s 1 -b 2 -e 2`, the preparation `BayesDel.pm` describes (lines 54 to
62). Each holds 690 rows: three, one per alternate base, at every 97th position
from 100 to 16,493 and at every position from 3540 to 3600. A row's score is
`(4 x q + i) / 100000` to five decimals, `q` the position whose base the row
states and `i` the place of its alternate base in `A`, `C`, `G`, `T`, so a
score tells which position's row was read.

| File | Contig | The row at position p states |
| --- | --- | --- |
| `scores_stating_the_reference_synthetic.txt.gz` | `M` | the reference's base at p |
| `scores_one_position_on_synthetic.txt.gz` | `M` | the reference's base at p + 1 |
| `scores_own_name_synthetic.txt.gz` | `MT` | the reference's base at p |

In the second file the row at 3556 states the base of 3557, which is the base
of 3556 as well, so that row carries the alleles of `same_next` and the score
derived for position 3557.

`ensembl_vep.stating_the_reference.txt`, `ensembl_vep.one_position_on.txt` and
`ensembl_vep.own_name.txt` are the outputs of Ensembl VEP 116.3 in the default
format, one per score file, from the image
`ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3`
with this directory mounted read-only at `/work/corpus` as the working
directory:

```bash
vep -i variants.vcf -o <out> --offline --cache --dir_cache <cache> \
  --cache_version 116 --species homo_sapiens --assembly GRCh37 --format vcf \
  --buffer_size 5000 --fork 1 --no_stats --quiet --force_overwrite \
  --dir_plugins /work/corpus \
  --plugin BayesDel,file=/work/corpus/scores_<name>_synthetic.txt.gz --no_headers
```

`BayesDel.pm` line 144 hands the record's own name to
`BaseVepTabixPlugin.pm` lines 693 to 708, which read the file under that name
with or without `chr`. Ensembl VEP 116.3 wrote, on this fixture, 22 rows for
each record against every file. Against the two files naming the contig `M` no
row carries a `BayesDel` item. Against `scores_own_name_synthetic.txt.gz` every
row of both records carries `BayesDel=` with the score derived for the record's
own position and alternate base, 44 rows in all.

vep-rs, given `--fasta reference_synthetic.fa`, writes the rows of
`ensembl_vep.own_name.txt` against `scores_stating_the_reference_synthetic.txt.gz`
and no `BayesDel` item against `scores_one_position_on_synthetic.txt.gz`.
Without `--fasta` it writes no item against either and warns once, naming the
table.

## The dbNSFP-shaped table

`dbNSFP4.9a_one_position_on_synthetic.tsv.gz` is a made-up table of eleven
columns in the layout `dbNSFP.pm` reads on GRCh37: `#chr`, `pos(1-based)`,
`ref`, `alt`, `aaref`, `aaalt`, `rs_dbSNP`, `hg19_chr`, `hg19_pos(1-based)`,
`Ensembl_transcriptid` and one value column, `synthetic_score`. It is sorted
and indexed on the `hg19` pair, `tabix -s 8 -b 9 -e 9`, the preparation the
module documents (lines 72 to 74), and its name carries a version because the
module refuses a name without one (lines 197 to 211). It names the contig `M`
and holds the 690 rows of the score files above, each stating the reference's
base one position on, as `scores_one_position_on_synthetic.txt.gz` does. The
amino acids of a row are those of a record at the row's own position with the
row's alternate base, so the module's amino-acid comparison (lines 340 to 344)
passes wherever its position and alternate-allele comparison (lines 335 to 338)
does.

The row at 3577 with the alternate base `C` states `G` as its reference, the
reference's base at 3578, and carries the score derived for 3578, `0.14313`.
A record `MT 3577 A C` has the row's position and alternate base and another
reference base. `dbNSFP.pm` line 313 names the record `M` to the table and
lines 335 to 338 compare neither the name's sequence nor the reference base,
so Ensembl VEP 116.3 writes `synthetic_score=0.14313` on the record's
`missense_variant` row under `--plugin dbNSFP,<table>,synthetic_score`. vep-rs,
given `--fasta reference_synthetic.fa`, finds the table's rows to state another
sequence and writes no `synthetic_score`.
