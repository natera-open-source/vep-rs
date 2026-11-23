# Fixture caches of the dump

Three small Storable caches in the layout of an Ensembl VEP cache, written by perl 5.34.0 with
Storable 3.23 as the caches of `tools/vep-storable-digest/tests/fixtures/caches/` are, which the
tests of this tool read as well.

| Cache | What it holds |
|---|---|
| `walk` | One transcript in each of the chromosomes `1`, `2`, `9`, `10`, `22`, `X`, `Y` and `MT`, named after its chromosome; a second shard of `2` whose root is a list, with a versioned id and an id of the first shard's again; and a chromosome `GL000191.1`, which is not read. |
| `predictions` | Chromosome `1`: eight transcripts whose `_variation_effect_feature_cache` is an unblessed hash with every kind of prediction (a blessed matrix, undef, `0`, a string, a hash without a matrix, without a digest, with both undef, with a number for a matrix, with a matrix of characters above 255, a hash blessed into the class `0`, one matrix under two names and again outside the predictions), a blessed hash, a hash whose predictions are blessed or a list, undef, a string, and a transcript without a cache. Chromosome `2`: a transcript one of whose matrices is a list. Chromosome `3`: a transcript whose matrix is a reference to a list. |
| `items` | Chromosome `4`: a list of undef, a transcript, a string and a transcript. |
