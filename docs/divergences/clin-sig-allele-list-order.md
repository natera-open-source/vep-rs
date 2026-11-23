# The order of the `CLIN_SIG` lists of two known variants

**The defect in one sentence.** Under the default `--clin_sig_allele 1`, Ensembl VEP 116.3 stores, as a key of a hash, the list of ClinVar terms that each known variant gives a row's allele, and writes `CLIN_SIG` from that hash's keys. Two or more different lists come in Perl's hash order, which changes from one row to the next and from one run to the next. vep-rs writes the lists in the order of the row's known variants, the same on every row and every run.

**Why this is a defect.** The record and the cache fix the lists: the different lists that the row's known variants give its allele, each once. Their order in Ensembl VEP's output comes from neither. The same input, flags and cache give the tab-format cell `likely_pathogenic;likely_pathogenic,pathogenic` on one run and `likely_pathogenic,pathogenic;likely_pathogenic` on another, and the rows of one record in one output file can carry both. The walk that fills the hash takes the known variants in the order `Existing_variation` names them, and the hash discards that order.

**What vep-rs does instead.** vep-rs writes each different list once, at the place of the first known variant that carries it. The lists therefore stand in the order `Existing_variation` names their variants.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `OutputFactory.pm`, `add_colocated_variant_info`: the known variants are walked in the order `Existing_variation` names them (lines 1069-1078). Each variant's entries are grouped by allele, the terms of one allele joined with `,` (lines 1105-1117), and the list of the row's allele is stored as a key of the hash `%clin_sigs` (line 1120).
- Lines 1143-1144 of the same function: the cell is `keys(%clin_sigs)` joined with `;`.

**What correcting it would change.** Keeping the lists in an array, each list once, in place of the hash's keys would write them in one order on every row, the order vep-rs writes.

**The class it produces.** `clin_sig_allele_list_order`: the `CLIN_SIG` value of a row whose known variants give its allele two or more different lists, under `--check_existing` without `--clin_sig_allele 0`. Ensembl VEP's lists are vep-rs's, in an order a hash decided. In a run that names their built-in class `order_only`, the golden harness and the format comparator read the value as a multiset of its lists, so such a row agrees, and the format comparator counts it under `order_only`. Nothing is counted under the entry's own name, and a run that does not name `order_only` reports the row as a difference.

**How often the shape occurs.** Tier 3 of [the index's order](../intended-divergences.md#differences-outside-the-classes): no variation file of a committed corpus cache holds two known variants with one start and end that give one allele different lists. 52 committed flag sets, the committed runs of the corpora, switch `--check_existing` on and leave `--clin_sig_allele` at 1. The table counts the spans of that shape over every variation file of the human Ensembl caches of release 116. A known variant listed twice under one name is read once, and a known variant is counted for an allele only where its allele string names that allele.

| Cache | Known variant rows | Spans of the shape |
| --- | --- | --- |
| Human, GRCh38 | 1,487,415,922 | 1 |
| Human, GRCh37 | 1,129,087,916 | 0 |

One record of the genome-wide ClinVar dataset on GRCh38, of 4,398,127, stands on that span with the allele the lists are for. Under `--everything`, which switches `--check_existing` on, its 47 rows carry the two lists in the default, tab and VCF formats. In the default-format output of one Ensembl VEP 116.3 run, 25 of the 47 rows carry the lists in the order vep-rs writes and 22 in the other.

## 1. Two lists written in two orders

**In short.** Two known variants at `21:500` carry ClinVar entries for the allele `G`: the first `likely_pathogenic`, the second `likely_pathogenic` and `pathogenic`. For the record `21 500 A G` under `--check_existing`, vep-rs prints `CLIN_SIG=likely_pathogenic%3Blikely_pathogenic,pathogenic`, the lists in the variants' order, on every run. Ensembl VEP 116.3 prints `CLIN_SIG=likely_pathogenic,pathogenic%3Blikely_pathogenic` with Perl's hash order fixed (`PERL_HASH_SEED=0`, `PERL_PERTURB_KEYS=0`), and either order without.

### 1.1 What a reader needs to know first

`--check_existing` names the known variants at a record's position in `Existing_variation`. A known variant of the cache can carry `clin_sig_allele`, ClinVar's terms by allele, written `allele:term` and joined by `;`. Under the default `--clin_sig_allele 1`, a row's `CLIN_SIG` reports the terms recorded for the row's own allele: the terms of one known variant joined by `,`, the lists of different variants joined by `;` ([Co-located variant keys](../output/known-variants.md#co-located-variant-keys)). The default and VCF formats write that `;` as `%3B`, the tab format (`--tab`) writes it as it is, and the VCF format writes the `,` as `&`.

A Perl hash holds its keys in no order of its own. `keys` returns them in an order set by the hash function and its seed, and by default Perl also perturbs that order hash by hash. Two hashes filled with the same keys in one process can therefore return them in different orders, unless the environment variables `PERL_HASH_SEED` and `PERL_PERTURB_KEYS` fix the order.

### 1.2 The record

Every identifier is invented. The cache's variation file holds two known variants at `21:500`, in this order:

| Cell | First known variant | Second known variant |
| --- | --- | --- |
| `variation_name` | `rs9920000001` | `rs9920000002` |
| `allele_string` | `A/G` | `A/G/T` |
| `clin_sig` | `likely_pathogenic` | `likely_pathogenic,pathogenic` |
| `clin_sig_allele` | `G:likely_pathogenic` | `G:likely_pathogenic;G:pathogenic` |
| `clin_sig_ref_allele` | `A` | `A` |
| `phenotype_or_disease` | `1` | `1` |

`clin_sig_ref_allele` is the reference allele the entries are written against. It equals the record's, so the entry alleles are read as written; where it differs, both engines read them reverse-complemented first. The input is one VCF record on that position, with the allele both variants carry entries for:

```
21	500	.	A	G	.	.	.
```

Both engines take the flag `--check_existing`. Section 1.3 gives the row's `Existing_variation` cell, and its `CLIN_SIG` cell in the default, tab and VCF formats.

### 1.3 What each engine writes

The record has one row, an intergenic one, since the cache holds no transcript there. Ensembl VEP 116.3 was run with `PERL_HASH_SEED=0` and `PERL_PERTURB_KEYS=0`.

| Field of the row | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Existing_variation` | `rs9920000001,rs9920000002` | `rs9920000001,rs9920000002` |
| `CLIN_SIG`, default format | `likely_pathogenic,pathogenic%3Blikely_pathogenic` | `likely_pathogenic%3Blikely_pathogenic,pathogenic` |
| `CLIN_SIG`, tab format | `likely_pathogenic,pathogenic;likely_pathogenic` | `likely_pathogenic;likely_pathogenic,pathogenic` |
| `CLIN_SIG`, VCF format | `likely_pathogenic&pathogenic%3Blikely_pathogenic` | `likely_pathogenic%3Blikely_pathogenic&pathogenic` |

The rows agree in every other column, among them `PHENO=1,1`, one flag for each known variant. With the environment variables unset, twelve runs of Ensembl VEP printed the one order six times and the other six times; vep-rs printed its order every time. The JSON format has no joined cell: each object of `colocated_variants` carries its own `clin_sig` array, the same in both engines.

With the two known variants in the other order in the cache, vep-rs writes `likely_pathogenic,pathogenic%3Blikely_pathogenic`, since `Existing_variation` then names the second variant first. A third known variant whose list for the allele equals the first variant's adds no third list. vep-rs writes each different list once, and Ensembl VEP's hash holds each different list as one key. A row whose known variants give its allele one list, and every row under `--clin_sig_allele 0`, is the same in both engines.

### 1.4 Why Ensembl VEP writes that

When a row is written, the output factory walks the row's known variants in the order `Existing_variation` names them: germline before somatic, then dbSNP, HGMD, COSMIC and other identifiers, the cache's order among equals. For each variant with ClinVar entries it joins the terms of the row's allele with `,` and stores the joined list as a key of a hash. A list two variants share is one key. After the walk it asks the hash for its keys and joins them with `;`.

That order is the hash's own. It depends on where the hash function places each key, and by default Perl perturbs it at random for each hash, so it moves from row to row inside one run. Nothing sorts the keys, and the order of the walk is not kept.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `OutputFactory.pm`, lines 1069-1078: the walk over the known variants, sorted.
- `OutputFactory.pm`, lines 1105-1120: the terms of each entry allele joined with `,`, and the row's allele's list stored as a key.
- `OutputFactory.pm`, lines 1143-1144: the keys joined with `;`.
- `OutputFactory/VEP_output.pm`, line 118, and `OutputFactory/VCF.pm`, lines 401-402: the `;` written as `%3B`, and in VCF the `,` as `&`.

### 1.5 Why vep-rs's row is the correct one

The cache gives the first known variant one term for `G` and the second two, so the cell holds two lists: `likely_pathogenic`, and `likely_pathogenic,pathogenic`. Both engines print those two lists. Only their order differs, and Ensembl VEP's own runs of this one input disagree with each other about it, so neither of its orders can be read from the record or the cache.

vep-rs's order is the row's own: its lists follow the known variants in the order `Existing_variation` names them, and both engines write the same `Existing_variation`. It is also the same on every row and every run, so two runs of one input can be compared byte for byte.

### 1.6 Try it yourself

The test file `crates/vep-cli/tests/colocated_list_fields.rs` builds the two known variants in a cache of its own and holds vep-rs to the rows above:

```
cargo test -p vep-cli --test colocated_list_fields
```

For Ensembl VEP 116.3, copy a GRCh38 cache of release 116 and replace `21/all_vars.gz` with a file of two rows. Each row is tab-separated in the column order that `variation_cols` gives in `info.txt`. It holds `21` under `chr`, `500` under `start`, the six cells of 1.2 and `.` in every other cell, `end` among them, which a cache leaves unset where it equals the start. Compress the file with `bgzip` and index it with `tabix -C -s 1 -b 5 -e 5`, the chromosome and start columns of that order. Then run the record with the cache options of [Reproducing these examples](reproducing.md):

```
vep -i record.vcf -o out.txt --offline --cache --dir_cache <copy> \
  --cache_version 116 --species homo_sapiens --assembly GRCh38 \
  --format vcf --check_existing
```

With `PERL_HASH_SEED=0` and `PERL_PERTURB_KEYS=0` in its environment Ensembl VEP writes the longer list first; without them, repeated runs write both orders.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/clin_sig_allele_list_order.toml` builds a cache for each engine with the two known variants, runs the record of 1.2 on both engines and records the `CLIN_SIG` value of its row. The reproduction harness sets `PERL_HASH_SEED=0` and `PERL_PERTURB_KEYS=0` for every Ensembl VEP run, so the recorded order is the one that setting gives.

The entry covers the order of the lists of a `CLIN_SIG` value that holds two or more. The golden harness and the format comparator split the value at the `;` between lists, `%3B` in the default and VCF formats, and compare the lists as a multiset under their built-in class `order_only`. Each list is compared whole: a changed term, a list written twice and the terms of one list in another order are reported like any other field. The JSON format is compared as it stands.

No committed corpus holds the shape, so the reading changes no comparison of a committed `CLIN_SIG` value.
