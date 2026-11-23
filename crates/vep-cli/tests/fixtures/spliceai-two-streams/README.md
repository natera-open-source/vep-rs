# SpliceAI score file of two joined streams

Seven records of the GRCh38 SpliceAI golden corpus
(`tests/golden/116/GRCh38-plugin-spliceai-fixture`), each a single-base change at one
of the corpus's own sites under a name of its own, in `variants.vcf`. The test
annotates them against the corpus's `json_cache` and reference, with the
corpus's made-up indel score file beside the made-up score file of this
directory.

`scores_two_streams.synthetic.vcf.gz` is a made-up score file, a VCF as the Ensembl VEP
plugin's documentation indexes one (`SpliceAI.pm:114-122`), with fabricated
scores, indexed with `tabix -p vcf`.
`scripts/fixtures/plugin_score_fixtures.py` wrote it from `variants.vcf` and
`symbols.tsv`, the gene symbols of each site, with `--shape spliceai_snv
--seed 41` and four options:

- `--second-stream-at 21:18279000`: the file is two bgzip streams joined. The
  second opens with the first row of `opens_second_stream`, so the first
  stream's end marker, a block that holds no bytes, stands inside the file
  before that row, as in a file made by joining two compressed files. `tabix`
  gives that row's chunk the empty block's offset.
- `--no-information 21:21324500`: the records at the site of `no_information`
  carry the missing-value `.` as their INFO column.
- `--join-alternates 21:21418500`: the site of `alternates_joined` has one
  record whose ALT column lists the three alternate bases.
- `--any-base-reference 21:26387000`: the records at the site of
  `any_base_reference` have `N` as their reference allele.

The generator writes under `<corpus>/data/`, so it wrote the file over a
directory holding a copy of `variants.vcf` and `symbols.tsv`, and the file and
its index are copies of what it wrote there.

`ensembl_vep.txt` is the output of Ensembl VEP 116.3 in the default format,
from the image `ensemblorg/ensembl-vep:release_116.3` with this directory
mounted read-only at `/work/corpus` as the working directory:

```bash
vep -i variants.vcf -o <out> --offline --cache --dir_cache <cache> \
  --cache_version 116 --species homo_sapiens --assembly GRCh38 --format vcf \
  --buffer_size 5000 --fork 1 --fasta <reference> --no_stats --quiet \
  --force_overwrite --dir_plugins <plugins> --no_headers \
  --plugin SpliceAI,snv=/work/corpus/scores_two_streams.synthetic.vcf.gz,indel=<indel>
```

`<reference>` is the corpus's reference, `<indel>` its indel score file, a
made-up file with fabricated scores, and `<plugins>` a `VEP_plugins` checkout
at release/116.

Ensembl VEP 116.3 wrote 33 rows. `SpliceAI_pred` stands on 18: the row of
`first_stream_one_gene` whose transcript has a gene symbol, the four rows of
`first_stream_two_genes`, the row of `opens_second_stream` and the twelve rows
of `second_stream_two_alternates`. No row of `no_information`,
`alternates_joined` or `any_base_reference` carries the field:
`SpliceAI.pm` line 272 takes a record whose REF and ALT columns equal the
variant's alleles, and lines 328 to 332 read no gene from an INFO column of
`.`.
