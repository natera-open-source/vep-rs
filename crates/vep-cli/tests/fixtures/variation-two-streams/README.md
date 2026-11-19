# Variation file of two joined streams

`all_vars.gz` holds the records of the GRCh38 co-located corpus's variation
file (`tests/golden/116/GRCh38-colocated/json_cache/variations/21/all_vars.gz`)
in that file's order, written as two bgzip streams joined and indexed with
`tabix -s 1 -b 5 -e 5`.

`write_fixture` of `scripts/fixtures/plugin_score_fixtures.py` wrote the file
and its index from those rows. The second stream opens with `rs1555830109` at
21:5131441, a site of the corpus's input, so the first stream's end marker, a
block that holds no bytes, stands inside the file before that row, as in a
file made by joining two compressed files. `tabix` gives that row's chunk the
empty block's offset.

The test puts the file in a copy of the corpus's cache in place of the
corpus's own file and holds the output to the one the corpus's cache gives.
