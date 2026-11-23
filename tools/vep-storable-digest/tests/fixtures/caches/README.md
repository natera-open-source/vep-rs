# Fixture caches

Small Storable caches in the layout of an Ensembl VEP cache: a directory a chromosome, and in it
gzipped Storable images named `<start>-<end>.gz`. Every image is written by perl 5.34.0 with
Storable 3.23 (`nstore_fd` into `gzip -n -c`, unless the table names another writer). `tools/vep-storable-digest`
and `tools/vep-storable-dump` test on them.

A transcript of these caches has the shape a cache stores one in: a hash blessed into
`Bio::EnsEMBL::Transcript` with its exons under `_trans_exon_array`, the same exon objects again
under the translation's `start_exon` and `end_exon` and under the `sorted_exons` and the introns of
`_variation_effect_feature_cache`, a translation whose `transcript` is a weak reference back to the
transcript, attributes, protein features with their analysis, a prediction matrix of 256 byte
values, and `dbID`, `adaptor`, `created_date` and `modified_date` keys at several depths. Integers
are stored both as numbers and as strings of digits.

| Cache | What it holds |
|---|---|
| `plain` | Chromosomes `1`, `2`, `3`, `10`, `X`, `MT` and `GL000191.1` (not read). `1` has three shards: seven items in the first (three transcripts, a string, an undef, a transcript without an id, one whose id is under `_stable_id`), and in the second a transcript of the first again with other `dbID` values and another `created_date`. `2` has a shard whose list is empty and one whose list is undef. `10` holds a shard whose root is a list. `1` and `3` hold files whose names are not `<digits>-<digits>.gz`. |
| `conflict` | Chromosome `21`, two shards, a transcript in both whose `biotype` differs and one in both that is equal. |
| `values` | Chromosome `1`: one transcript that carries, under `scalars`, `flagged` and `refs`, every kind of value: integers and strings of digits at the limits of 64 bits, doubles, strings with quotes, backslashes and bytes above 127, strings and keys with Perl's UTF-8 flag, references to scalars and to references, a blessed scalar, a hash reached four times, an array that holds itself, an array with holes, a blessed array, a hash blessed into the class `0`, a compiled pattern, version strings, Perl's immortal scalars, and the five keys of the content rule nested in several of them. Chromosome `2`: one transcript stored in the machine's own byte order (`store_fd`), where integers and doubles keep their kind. Storable stores a whole double under two to the 53rd as an integer, so the image holds those as integers, one of them the smallest integer of eleven characters; one integer of seventeen digits is a value no double holds. |
| `ids` | Chromosome `5`, one shard of eighteen items: ids with and without a version, a version of several parts, an id that is an integer or a double, an empty id, an id that ends in a newline, ids with a byte above 127 with and without the UTF-8 flag, a character above 255, undef and non-numeric coordinates, an unblessed hash, a hash blessed into the class `0`, a reference to a scalar. |
| `files/plain-image` | A shard that is a Storable image without compression. |
| `files/two-members` | A shard followed by a second gzip member. |
| `files/trailing` | A shard followed by bytes that are no gzip member. |
| `files/no-trailer` | A shard whose gzip member lacks its last eight bytes. |
| `files/bad-checksum` | A shard whose gzip member carries a wrong checksum. |
| `files/short-image` | A shard whose image lacks its last twenty bytes, and a whole shard after it. |
| `files/no-image` | A gzipped line of text. |
| `files/empty-file` | A file of no bytes. |
| `files/empty-member` | A gzip member of no bytes. |
| `files/short-gzip` | The first half of a gzip member, and a whole shard after it. |
| `files/array-item` | A list whose second item is a reference to an array. |
| `files/hash-list` | A hash whose value is a hash in the place of a list. |
| `files/string-list` | A hash whose value is a string in the place of a list. |
| `files/blessed-root` | A blessed hash at the root. |
| `files/scalar-root` | A reference to a string at the root. |
| `files/id-reference` | A transcript whose `stable_id` is a reference to an array. |
