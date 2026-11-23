# BayesDel mitochondrial contig fixture

Two fabricated single-base records on the mitochondrial contig of GRCh37,
spelt `MT`, inside a protein-coding transcript: `scored` at `MT:3555` and
`no_row` at `MT:3574`. `variants.vcf` holds the two records. `json_cache`
holds the 22 transcripts within 5,000 bases of them, cut from the release-116
GRCh37 cache by `scripts/golden/prune_json_cache.py --gzip --drop sorted_exons
--drop protein_features --drop protein_function_predictions --perl-info
<cache>/info.txt --cache-version 116`.

Two BayesDel-shaped score files of five columns, bgzip-compressed and indexed
with `tabix -s 1 -b 2 -e 2`, the preparation `BayesDel.pm` describes (lines 54 to 62). Each
holds three fabricated rows at the site of `scored`, one per alternate base,
and no row at the site of `no_row`. `scripts/fixtures/plugin_score_fixtures.py`
wrote both from `variants.vcf` with `--shape bayesdel --seed 37 --miss
MT:3574`, so their rows differ in the first cell alone:

- `scores_other_pair_synthetic.txt.gz` names the contig `M` (`--mitochondrial-name M`)
  under the header `#Chr pOs Ref Alt bAyEsDeL_aDdAf` (`--header-names
  Chr,pOs,Ref,Alt,bAyEsDeL_aDdAf`): the generator's header for the shape, the
  second and the fifth name in another letter case.
- `scores_own_pair_synthetic.txt.gz` names the contig `MT`, as the records do, under
  the header `#Chr second_place Ref Alt fifth_place` (`--header-names
  Chr,second_place,Ref,Alt,fifth_place`): other words in the second and the
  fifth place.

The generator writes under `<corpus>/data/`, so it wrote the two files over a
directory holding a copy of `variants.vcf`, and the files and their indexes
are copies of what it wrote there.

On GRCh37 `MT` and `chrMT` name the assembly's mitochondrial sequence and
`chrM` and `M` another one. `BayesDel.pm` line 144 hands the record's own name
to `BaseVepTabixPlugin.pm` lines 693 to 708, which read the file under that
name with or without `chr`, so a record on `MT` reads no row of a contig the
file names `M`. `BayesDel.pm` line 151 splits a row into five variables by
place, so the header's names decide nothing.

`ensembl_vep.other_pair.txt` and `ensembl_vep.own_pair.txt` are the outputs
of Ensembl VEP 116.3 in the default format, one per score file, from the image
`ensemblorg/ensembl-vep:release_116.3` with this directory mounted read-only
at `/work/corpus` as the working directory:

```bash
vep -i variants.vcf -o <out> --offline --cache --dir_cache <cache> \
  --cache_version 116 --species homo_sapiens --assembly GRCh37 --format vcf \
  --buffer_size 5000 --fork 1 --no_stats --quiet --force_overwrite \
  --dir_plugins /work/corpus \
  --plugin BayesDel,file=/work/corpus/scores_<pair>_synthetic.txt.gz --no_headers
```

Ensembl VEP 116.3 wrote, on this fixture, 22 rows for each record against
either file. Against `scores_other_pair_synthetic.txt.gz` no row of either record
carries a `BayesDel` item. Against `scores_own_pair_synthetic.txt.gz` every row of
`scored` carries `BayesDel=` with the score of the row for the record's own
alternate base, and no row of `no_row` carries the item.
