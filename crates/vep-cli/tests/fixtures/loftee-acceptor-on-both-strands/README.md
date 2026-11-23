# LoFTEE acceptor on both strands fixture

A fabricated sequence and a fabricated gene model in which two protein-coding
transcripts on opposite strands share their splice-acceptor bases, with two
fabricated single-base records, one on each shared acceptor.

`two_strand_acceptor.fa` is one sequence, `toy`, of 1,000 fabricated bases,
with its index. `two_strand_acceptor.gff3.gz` is the model, bgzip-compressed
and indexed with `tabix -p gff`:

| Transcript | Strand | Exons | Coding sequence |
| --- | --- | --- | --- |
| `FWDT1` of gene `FWDG1` | forward | 101-200, 401-500, 701-800 | 150-200, 401-500, 701-750 |
| `REVT1` of gene `REVG1` | reverse | 850-950, 620-698, 320-398 | 850-900, 620-698, 349-398 |

Bases 399-400 end the first intron of `FWDT1` and the second of `REVT1`, and
bases 699-700 end the second intron of `FWDT1` and the first of `REVT1`. Each
pair is the canonical acceptor of one of the two transcripts alone.
`variants.vcf` holds `acceptor_forward_motif` at `toy:400` and
`acceptor_reverse_motif` at `toy:699`, a `splice_acceptor_variant` of both
transcripts each. The nine reference bases around 400 hold the motif of
LoFTEE's `NAGNAG_SITE` flag read on the forward strand and not on the reverse
strand; the nine around 699 hold it on the reverse strand alone.

`ensembl_vep.txt` and the four files of `expected/` are the outputs of Ensembl
VEP 116.3 with the LoF plugin of the LOFTEE repository's `grch38` branch (tag
`v1.0.4_GRCh38`, its `maxEntScan/` files replaced by the stand-ins of
`scripts/data/write_loftee_stand_ins.py`), from the image
`ensemblorg/ensembl-vep:release_116.3` with this directory mounted read-only
at `/work/corpus` as the working directory and the plugin's checkout at
`/work/loftee`:

```bash
vep -i variants.vcf -o <out> --gff two_strand_acceptor.gff3.gz \
  --fasta two_strand_acceptor.fa --species homo_sapiens --assembly GRCh38 \
  --format vcf --buffer_size 5000 --fork 1 --no_stats --quiet \
  --force_overwrite --dir_plugins /work/loftee \
  --plugin LoF,loftee_path:/work/loftee,use_gerp_end_trunc:0,run_splice_predictions:0
```

`ensembl_vep.txt` is the default format under `--no_headers`. `expected/`
holds the default, `--tab`, `--vcf` and `--json` outputs with their headers,
gzip-compressed, under the file names of a golden corpus's `expected`
directory. `use_gerp_end_trunc:0` selects the plugin's rule that reads no GERP
file, in both engines.

Ensembl VEP 116.3 wrote two rows for each record. On the rows of `FWDT1` it
wrote `NAGNAG_SITE` at `toy:400` and not at `toy:699`, and on the rows of
`REVT1` the same: the plugin reads the nine bases once per record, on the
strand of the first transcript whose acceptor row it tests, and keeps them
for the record's other transcripts (`LoF.pm` lines 489-494). vep-rs reads
them on each row's own strand, so its rows of `REVT1` carry the flag at
`toy:699` and not at `toy:400`
([`loftee_nagnag_site_first_transcript_strand`](../../../../../docs/divergences/loftee-nagnag-site-first-transcript-strand.md)).
