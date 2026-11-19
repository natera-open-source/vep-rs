# Regulatory test shards

The shards hold the converter's shape: a regulatory feature carries its stable
id, span, feature type, database id and `activity`, one code per name of
`cell_types.json` in order (`A` active, `I` inactive, `P` poised, `R` repressed,
`N` not assessed, `-` not in the feature's map); a motif feature carries its
stable id, span, strand, database id, score, stranded sequence and the stable id
of its matrix in `binding_matrices.json`.

`cell_types.json` lists the 112 epigenomes of the release 116 GRCh38 regulatory
build, the cache's `info.json` `cell_types` list; both shards' activity strings
are aligned to it.

`1-1000000.json` is hand-built: an enhancer at 100-120 of chromosome 21, active
in the first forty epigenomes, inactive in the next forty, then one poised, one
repressed and one not assessed, holding a forward motif at 105-108 with a
two-factor matrix and a reverse motif at 115-118 with no factors, the shapes the
tests of insertion edges, strand positions, activity codes and empty factor
lists need.

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

`binding_matrices.json` holds the three matrices the slice names, copied from
the converted cache's table, and the two hand-built ones (`ENSPFM_FORWARD`,
`ENSPFM_REVERSE`: four positions, one count each).
