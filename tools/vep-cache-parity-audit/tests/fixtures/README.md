# Fixture caches of the audit

`storable/` holds small Storable caches written by perl 5.34.0 with Storable 3.23 (`nstore_fd`
into `gzip -n -c`), in the layout of an Ensembl VEP cache. `json/` holds JSON caches:
`tools/vep-storable-to-json` wrote `audit`, `tiny` and `audit-walk` from the Storable caches of
the same names, and `gaps` and `junk` are copies of `audit` edited by hand. The tests also read
two caches of `tools/vep-storable-digest/tests/fixtures/caches/`: `files` (a list whose second
item is a reference to an array) and `values` (an array that holds itself).

| Cache | What it holds |
|---|---|
| `storable/audit` | One shard of five transcripts in the shape a cache stores one in. Two carry a gene object, the flag of the cache loader, both UTR objects and a sequence edit; one has no translation. The fixture's transcripts also carry fields no real cache has under exons, introns and features (`dbID`, `adaptor`, dates, a slice), which the audit reports as gaps. |
| `storable/tiny` | One shard of five transcripts of nine fields each; one has an undef `_ccds`, two a field of their own. |
| `storable/audit-walk` | A cache under `homo_sapiens/116_GRCh38/`: chromosomes `1` (two shards), `2`, `MT`, `GL000191.1`; `3`, whose transcript stands alone under its chromosome; `4`, whose root is a list; `5`, whose root is blessed; `6`, whose list holds a string, an undef, a transcript without a `stable_id`, three whose `stable_id` is undef, `0` and empty, one whose `stable_id` is the number 12345 and one more; and in `1` a `_var.gz` and a `_reg.gz` file, a gzip file that is no shard and a file that is no gzip. |
| `json/audit`, `json/tiny`, `json/audit-walk` | The converter's output for the three. `audit-walk` lacks the shards of `3`, `4` and `5`, which the converter does not write. |
| `json/gaps` | `audit` with one transcript taken out; in the first, `biotype` and the peptide taken out, two fields added, `start` made a list, `strand` a boolean and `exons` a string; in the third, the translation made `null` and the attributes `false`; in the fifth, the cache made `null` and a field added whose keys are empty or hold a dot. |
| `json/junk` | `audit` with, beside its shard: a file that sorts last and holds the first transcript again with another `biotype`; a file that is no JSON, an empty one, a string and `true` at the top of a file; a file whose top is a list; a file of transcripts whose `stable_id` is every kind of value; a shard outside `transcripts/`; a text file and a file named `.JSON`. |
