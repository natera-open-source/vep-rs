# The `FREQS` list order under `--check_frequency`

**The defect in one sentence.** Under `--check_frequency`, Ensembl VEP 116.3 keeps the frequency it judged each alternate allele by in a hash keyed by allele, and writes the `FREQS` list from that hash's keys. The items of a record with two or more alternate alleles come in Perl's hash order, which changes from one run to the next. vep-rs lists the items in the order of the record's alternate alleles, the same on every run.

**Why this is a defect.** The record and the cache fix the list's items, one `POPULATION:ALLELE:frequency` per alternate allele. Their order in Ensembl VEP's output comes from neither. The same input, flags and cache give `AF:G:1,AF:A:NA` on one run and `AF:A:NA,AF:G:1` on another, so two runs of one command differ byte for byte. The filter visits the alleles in the record's order, and the hash it stores them in discards that order.

**What vep-rs does instead.** vep-rs lists the items in the order of the record's alternate alleles, so every run of an input writes the same list.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `AnnotationSource/Cache/BaseCacheVariation.pm`, lines 405-441, `_add_check_freq_data_to_vf`: the alternate alleles are walked in the record's order (line 411). Each allele's frequency, or `NA`, is stored in the hash `%used_freqs` under the allele (line 435), and the hash is kept on the variant under the population's name (line 441).
- `OutputFactory.pm`, lines 1149-1163: the `FREQS` list is built by walking `keys` of the population hash and, inside it, `keys` of the allele hash (lines 1152-1153).

**What correcting it would change.** Walking the record's alternate alleles in place of the hash's keys would write the items in one order on every run, the order vep-rs writes.

**The class it produces.** `check_frequency_freqs_list_order`: the `FREQS` value of a record with two or more alternate alleles and frequency data, under `--check_frequency` or `--filter_common`. The reference's items are the candidate's, in an order the run's hash seed decided. The golden harness and the format comparator compare a `FREQS` value as a multiset of its items, so no comparison reports the order and nothing is counted under the class.

**How often the shape occurs.** Tier 2: the shape is a record with two or more distinct alternate alleles, 405.6 per million records of the cuts. 5 committed flag sets pass `--check_frequency` or `--filter_common`, and 80 corpora hold a record of the shape. The rate is an upper bound, since the cache decides whether a known variant gives a record frequency data. Every `FREQS` value in the committed reference outputs lists one item.

## 1. Two alleles listed in two orders

**In short.** The record `21 7759128 T A,G` lies on the known variant rs1555870136, whose `AF` frequency for `G` is 1. Under `--check_frequency --freq_pop AF --freq_freq 0.01 --freq_gt_lt gt --freq_filter exclude`, allele `G` passes and is dropped, and allele `A`'s row carries the list. Ensembl VEP 116.3 prints `FREQS=AF:G:1,AF:A:NA` under a hash seed fixed at 0, and either order under Perl's own seed. vep-rs prints `FREQS=AF:A:NA,AF:G:1`, the record's order, on every run.

### 1.1 What a reader needs to know first

`--check_frequency` judges each alternate allele of a record by the frequency a co-located known variant gives it in the population `--freq_pop` names. With `--freq_gt_lt gt` and `--freq_filter exclude`, an allele whose frequency is at or above `--freq_freq` is dropped and the others are kept. Every row of a record that had frequency data carries `FREQS`, which lists `POPULATION:ALLELE:frequency` for each alternate allele of the record, dropped alleles included, with `NA` for an allele the variant gives no frequency ([Filtering by frequency](../output/known-variants.md#filtering-by-frequency)).

A Perl hash holds its keys in no order of its own. `keys` returns them in an order set by the hash function and its seed, and Perl draws a new seed for each process unless the environment variable `PERL_HASH_SEED` fixes one.

### 1.2 The record

The cache is the `GRCh38-colocated` golden corpus's (`tests/golden/116/GRCh38-colocated/json_cache`). Its variation file holds rs1555870136 at `21:7759128` with the alleles `T/A/G` and the `AF` value `G:1`. The input is one VCF record with that variant's two alternate alleles, in its order:

```
21	7759128	two_alts	T	A,G	.	.	.
```

Both engines take the default format and these flags:

```
--check_frequency --freq_pop AF --freq_freq 0.01 --freq_gt_lt gt --freq_filter exclude
```

### 1.3 What each engine writes

The record's one row. Ensembl VEP 116.3, run with `PERL_HASH_SEED=0` and `PERL_PERTURB_KEYS=0`:

```
two_alts	21:7759128	A	-	-	-	intergenic_variant	-	-	-	-	-	rs1555870136	IMPACT=MODIFIER;FREQS=AF:G:1,AF:A:NA
```

vep-rs:

```
two_alts	21:7759128	A	-	-	-	intergenic_variant	-	-	-	-	-	rs1555870136	IMPACT=MODIFIER;FREQS=AF:A:NA,AF:G:1
```

The rows agree in every other column. With the environment variables unset, repeated runs of Ensembl VEP printed `AF:G:1,AF:A:NA` on some and `AF:A:NA,AF:G:1` on others; vep-rs printed `AF:A:NA,AF:G:1` every time. The tab format prints the same two values, and the VCF format joins the items with `&`: `AF:G:1&AF:A:NA` against `AF:A:NA&AF:G:1`. `--filter_common --freq_pop AF` gives the same two lists.

Ensembl VEP's order is the hash's: the alleles, their order in the record and the seed each move it. With the alleles written `G,A`, Ensembl VEP under seed 0 prints `AF:A:NA,AF:G:1` and vep-rs `AF:G:1,AF:A:NA`. For `21 44195140 G A,C`, on rs2838519 (`AF` `A:0.5625`), Ensembl VEP prints `AF:A:0.5625,AF:C:NA` under seed 0, as vep-rs does, and `AF:C:NA,AF:A:0.5625` under seed 11. A third allele, `T A,G,C`, gives `AF:G:1,AF:A:NA,AF:C:NA` against vep-rs's `AF:A:NA,AF:G:1,AF:C:NA`. A record with one alternate allele, or with one allele written twice, lists one item in both engines.

### 1.4 Why Ensembl VEP writes that

The frequency filter takes the record's alternate alleles in order and, for each, notes the frequency it judged the allele by in a hash keyed by the allele. It stores that hash on the variant under the population's name. When a row is written, the output factory asks the hash for its keys and prints one item per key, in the order `keys` returns them.

That order is the hash's own. It depends on where the seeded hash function places each key and on the order the keys were stored in, and by default Perl also perturbs it at random as keys are stored. Under that default a seed is drawn per process, so the order moves from run to run. Nothing sorts the keys, and the order in which the filter visited the alleles is not kept.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `AnnotationSource/Cache/BaseCacheVariation.pm`, lines 411 and 435: the walk over the alternate alleles and the store into the hash.
- `OutputFactory.pm`, lines 1152-1153: the two `keys` walks that build the list.

### 1.5 Why vep-rs's row is the correct one

The record names its alternate alleles `A` then `G`, and the list holds one item for each. Both engines print the same two items, `AF:A:NA` and `AF:G:1`. Only their order differs, and Ensembl VEP's own runs of this one input disagree with each other about it, so neither of its orders can be read from the record or the cache. vep-rs's order is the record's own, the order in which the filter judged the alleles in both engines. It is also the same on every run, which a file compared against an earlier run of itself needs.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 \
  --check_frequency --freq_pop AF --freq_freq 0.01 \
  --freq_gt_lt gt --freq_filter exclude
```

vep-rs writes `FREQS=AF:A:NA,AF:G:1`. Ensembl VEP 116.3 takes the same flags with the cache options of [Reproducing these examples](reproducing.md). With `PERL_HASH_SEED=0` and `PERL_PERTURB_KEYS=0` in its environment it writes `FREQS=AF:G:1,AF:A:NA`; without them, repeated runs write both orders.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/check_frequency_freqs_list_order.toml` runs the record of 1.2 on both engines and records the `FREQS` item of its row. The reproduction harness fixes Perl's hash seed at 0 for every Ensembl VEP run, so the recorded order is the one that seed gives.

The class covers the order of a `FREQS` list of two or more items. The golden harness and the format comparator compare the value as a multiset of its items, under their built-in class `order_only`. The same class compares the `ZYG` list of `--individual_zyg`, whose order is a hash's as well, and `DOMAINS`, whose order follows the region file Ensembl VEP reads a transcript's copy from ([the `DOMAINS` order of a transcript stored in two region files](domains-order-region-file-copy.md)). A list whose items differ is reported like any other field.
