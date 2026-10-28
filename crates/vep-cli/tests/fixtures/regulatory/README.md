# Regulatory test shards

`1-1000000.json` is hand-built: an enhancer at 100-120 of chromosome 21 holding a
forward motif at 105-108 with a two-factor matrix and a reverse motif at 115-118
with no factors, the shapes the tests of insertion edges, strand positions and
empty factor lists need.

`chr21-slice.json` is cut from the chromosome 21 regulatory shards of the
converted Ensembl VEP release 116 GRCh38 cache, unchanged: the enhancer
`ENSR21_C6M55` with the four `NRF1` motifs inside it at 46699248-46699265, and
the two motifs `ENSM00000574814` (`ENSPFM0290`, three factors) and
`ENSM00000574815` (`ENSPFM0484`) at 5090718-5090738, which lie in no regulatory
feature. Regenerate it from a converted cache with

```bash
jq -c '[.[] | select(.stable_id | IN("ENSR21_C6M55","ENSM00000588223","ENSM00000588224","ENSM00000588225","ENSM00000588226","ENSM00000574814","ENSM00000574815"))]' \
  <cache>/regulatory/21/*.json | jq -s 'add | sort_by(.start)'
```

and write one feature per line.
