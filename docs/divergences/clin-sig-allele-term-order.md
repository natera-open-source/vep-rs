# The order of one allele's ClinVar terms

**The defect in one sentence.** When Ensembl VEP 116.3 reads its variation database, it collects a known variant's ClinVar entries, each an `allele:term` pair, as the keys of a hash and writes `clin_sig_allele` from those keys. The pairs therefore stand in Perl's hash order, in a cache it dumps and in a run under `--database`. A cache that `vep-cache-builder` builds from the same database's tables holds the pairs in the order of the ClinVar records' identifiers. `CLIN_SIG` lists an allele's terms in the cell's order, so a row can list the same terms in another order on the two caches.

**Why this is a defect.** The database settles which pairs a variant has: each different pair of an allele and a term among its ClinVar records. It gives them no order, and Ensembl VEP's comes from the hash alone: two `--database` runs of one record, under two hash seeds, print the two orders. Its own GRCh38 cache of release 116 holds 20,864 lists of exactly the terms `likely_pathogenic` and `pathogenic`, 10,416 in that order and 10,448 in the other. An order by the terms' names would be the same for every such list, and the lists follow the record identifiers no more often than not, as the chromosome 21 counts below show.

**What vep-rs does instead.** `vep` writes a known variant's terms in the order of the cache's cell, as Ensembl VEP does, so on one cache the two engines write the same list for it. `vep-cache-builder` takes a variant's ClinVar records in the order of their identifiers and writes each different pair at the place of the first record that carries it. The same tables give the same cell on every build.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`; `ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/`):

- `AnnotationSource/Database/Variation.pm`, `get_features_by_regions_uncached`: each ClinVar record filed under the variant's name sets one key of the hash `%clin_sigs`, its allele and its term joined by `:` (lines 200-211). The cell is `keys(%clin_sigs)` joined with `;` (lines 212-213).
- `Pipeline/DumpVEP/Dumper/Variation.pm`, lines 207-209 and 247: the cache dump calls that routine and writes the cell as it stands.
- `OutputFactory.pm`, lines 1105-1117: for each known variant, a row's `CLIN_SIG` joins the terms of the row's allele in the cell's order.
- `DBSQL/PhenotypeFeatureAdaptor.pm` of `ensembl-variation`, `get_clinsig_alleles_by_location`, lines 1230-1258: the query orders the records by sequence region, start and end alone (line 1257).

**What correcting it would change.** Keeping the pairs in an array, each once, and ordering the query by the record identifier as well would write one order in every cache and every run: the order `vep-cache-builder` writes.

**The class it produces.** The entry's name is `clin_sig_allele_term_order`. It stands for the order of the `allele:term` pairs of one known variant's `clin_sig_allele`, between a cache Ensembl VEP's dump wrote, or a `--database` run, and a cache `vep-cache-builder` built from the same release's tables. It shows in `CLIN_SIG` on a row whose allele has two or more terms, under `--check_existing` without `--clin_sig_allele 0`, and in the `clin_sig_allele` string of the JSON format. No comparison names a class for it. Two runs on one cache read one cell and must write one list, so the golden harness and the format comparator report a list whose terms stand in another order as a difference.

**How often the shape occurs.** Tier 2 of [the index's order](../intended-divergences.md#differences-outside-the-classes). The shape is a known variant whose ClinVar records give one allele two or more terms, in a cache `vep-cache-builder` built with the ClinVar tables. Nine rows of the committed native cache, `tests/golden/116/GRCh38-native-builder/json_cache`, hold it, and one committed run of that corpus, its flag set `variation`, passes `--check_existing`. The table counts such lists in the `clin_sig_allele` column of every variation file of the human Ensembl caches of release 116, a list being the terms one known variant gives one allele.

| Cache | Known variant rows | Lists of two or more terms | Sets of terms written in more than one order |
| --- | --- | --- | --- |
| Human, GRCh38 | 1,487,415,922 | 306,349 | 387 of 571 |
| Human, GRCh37 | 1,129,087,916 | 159,646 | 327 of 542 |

A cache that `vep-cache-builder` built from the release's chromosome 21 variation VCF and ClinVar tables shares 2,901 such lists, term for term, with Ensembl VEP's GRCh38 cache: 1,313 stand in the same order in both and 1,588 do not. Of the 2,532 that hold two terms, 1,259 stand in the same order and 1,273 in the other.

## 1. Two terms of one allele in two orders

**In short.** The known variant `rs1569061768` at `21:34859477` has ClinVar records for two of its alleles. For the record `21 34859477 G A` under `--check_existing`, vep-rs on the committed native cache prints `CLIN_SIG=pathogenic,likely_pathogenic`, the order of the ClinVar records' identifiers. Ensembl VEP 116.3 on its own cache prints `CLIN_SIG=likely_pathogenic,pathogenic`, the order a hash had when the cache was dumped.

### 1.1 What a reader needs to know first

`--check_existing` names the known variants at a record's position. A known variant can carry `clin_sig_allele`, ClinVar's terms by allele, written `allele:term` and joined by `;`. Under the default `--clin_sig_allele 1` a row's `CLIN_SIG` reports the terms a known variant records for the row's own allele, joined by `,` in the order the cell holds them ([Co-located variant keys](../output/known-variants.md#co-located-variant-keys)). Where two known variants give the allele different lists, `CLIN_SIG` joins the lists by `;`, and their order is [an entry of its own](clin-sig-allele-list-order.md).

Ensembl VEP's cache is a dump of its variation database. A native cache is one `vep-cache-builder` builds. It reads the known variants from Ensembl's export of that database, the variation VCFs. Given the release's `phenotype_feature` and `phenotype_feature_attrib` tables, it fills the same cell from the same ClinVar records ([Variation data in a native cache](../cache-variation-data.md)).

A Perl hash holds its keys in no order of its own. `keys` returns them in an order set by the hash function and its seed. By default Perl draws the seed at random for each process and also perturbs the order hash by hash, unless the environment variables `PERL_HASH_SEED` and `PERL_PERTURB_KEYS` fix both.

### 1.2 The record

The variation database of release 116 holds six germline ClinVar records for `rs1569061768`, four for one alternate allele and two for another. They are the rows of its `phenotype_feature` table filed under the variant's name, with the `risk_allele` and `clinvar_clin_sig` values of `phenotype_feature_attrib`. By identifier:

| `phenotype_feature_id` | Allele | `clinvar_clin_sig` |
| --- | --- | --- |
| 298618465 | `A` | pathogenic |
| 298643442 | `A` | pathogenic |
| 298692902 | `A` | likely pathogenic |
| 298808066 | `A` | pathogenic |
| 299107564 | `C` | uncertain significance |
| 299198667 | `C` | likely pathogenic |

The two caches hold the same four pairs in the variant's `clin_sig_allele` cell, a term's blanks written as underscores:

| Cache | `clin_sig_allele` |
| --- | --- |
| Ensembl VEP's own, release 116, GRCh38 | `A:likely_pathogenic;C:likely_pathogenic;A:pathogenic;C:uncertain_significance` |
| The committed native cache | `A:pathogenic;A:likely_pathogenic;C:uncertain_significance;C:likely_pathogenic` |

The input is one VCF record with the allele of the first four records:

```
21	34859477	runx1_two_terms	G	A	.	.	.
```

Both engines take the flag `--check_existing`.

### 1.3 What each engine writes

The record has thirteen rows, one for each transcript in reach of the position. Each engine writes one `CLIN_SIG` value on all thirteen:

| Run | `CLIN_SIG` |
| --- | --- |
| Ensembl VEP 116.3 on its own cache, `homo_sapiens/116_GRCh38` | `likely_pathogenic,pathogenic` |
| vep-rs on the committed native cache | `pathogenic,likely_pathogenic` |

In the JSON format the known variant's `clin_sig_allele` is the cell as each cache holds it, the two strings of 1.2. Its `clin_sig` array, a cell the two caches hold alike, is the same in both. On a conversion of Ensembl VEP's cache vep-rs prints Ensembl VEP's order, since it then reads Ensembl VEP's cell.

The order also moves between two runs of Ensembl VEP on one record. The record at `21:33345322`, line 52 of `tests/golden/116/GRCh38-everything/variants.vcf`, has three alternate alleles. The ClinVar records of its known variant `rs34849215` give one of them the terms `pathogenic` and `likely pathogenic`. Ensembl VEP 116.3 was run on it under `--database --db_version 116 --check_existing`, where it builds the cell during the run:

| `PERL_HASH_SEED` | `CLIN_SIG` on the rows of that allele |
| --- | --- |
| `0` | `likely_pathogenic,pathogenic` |
| `11` | `pathogenic,likely_pathogenic` |

Both runs set `PERL_PERTURB_KEYS=0`. Their rows are equal but for that `CLIN_SIG` value.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP reads a region's known variants from the variation database in one routine, which its cache dump and a `--database` run both call. The routine fetches the region's ClinVar records, keyed by their span. For each variant it takes the records of the variant's span that are filed under the variant's own name, and stores each record's allele and term as one key of a hash. Two records with one allele and one term are one key. It then asks the hash for its keys and joins them with `;`.

The query that returns the records orders them by position alone, so the records of one variant arrive in no stated order, and the hash keeps none. A cache dump therefore freezes whichever order each variant's hash had in the process that wrote it.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`; `ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/`):

- `AnnotationSource/Database/Variation.pm`, line 161: the region's ClinVar records; line 197: the records of the variant's span.
- `AnnotationSource/Database/Variation.pm`, lines 200-211: one key for each pair; lines 212-213: the keys joined with `;`.
- `DBSQL/PhenotypeFeatureAdaptor.pm`, lines 1230-1258: the query, ordered by sequence region, start and end (line 1257).
- `Pipeline/DumpVEP/Dumper/Variation.pm`, lines 207-209 and 247: the dump's call and the cell written.

### 1.5 Why vep-rs's row is the correct one

The six records give the input's allele two terms, `pathogenic` and `likely pathogenic`, and both engines print those two, so both rows hold the right terms. Only their order differs. The database gives the two terms no order: each stands in records of its own, and the records carry no rank. Ensembl VEP's own cache orders the same two terms one way for some variants and the other way for others, and its own runs disagree with each other, so its order cannot be read from the record or the database.

vep-rs's order adds no rule of its own to the data: the pairs follow the ClinVar records as the database numbers them, each pair at its first record. Here the allele's first record reads `pathogenic` and its third `likely pathogenic`. A build from the same tables writes the same cell every time, so two native caches of one release can be compared byte for byte.

### 1.6 Try it yourself

Two tests hold vep-rs's side. The first builds a native cache from ClinVar tables whose rows stand in another order than their identifiers, and holds the cell. The second runs the record of 1.2 on the committed native cache and holds its rows:

```
cargo test -p vep-cache-builder --test variation_cache_from_dumps
cargo test -p vep-cli --test clinvar_terms_of_one_allele_on_a_native_cache
```

To see vep-rs's row, save the record of 1.2 as `record.vcf` and run vep-rs's `vep` on the committed cache:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38-native-builder/json_cache \
  --species homo_sapiens --assembly GRCh38 --format vcf --check_existing
```

For Ensembl VEP 116.3, run its `vep` on the same record and the release's own GRCh38 cache, with the cache options of [Reproducing these examples](reproducing.md):

```
vep -i record.vcf -o out.txt --offline --cache --dir_cache <cache> \
  --cache_version 116 --species homo_sapiens --assembly GRCh38 \
  --format vcf --check_existing
```

For the two runs of 1.3 under `--database`, which reach Ensembl's public database, save line 52 of `tests/golden/116/GRCh38-everything/variants.vcf` as `record.vcf`. Set the seed in Ensembl VEP's environment, `11` here and `0` for the other order, and read `CLIN_SIG` on the rows whose value holds both terms:

```
PERL_HASH_SEED=11 PERL_PERTURB_KEYS=0 vep -i record.vcf -o out.txt \
  --database --db_version 116 --species homo_sapiens --assembly GRCh38 \
  --format vcf --check_existing
```

### 1.7 How a comparison treats such a row

The reproduction `scripts/divergences/reproductions/clin_sig_allele_term_order.toml` runs the record of 1.2 on Ensembl VEP with its own cache and on vep-rs with the committed native cache. It records the `CLIN_SIG` value on the row of `ENST00000300305`.

No comparison sets the difference aside. The golden harness and the format comparator compare each list of a `CLIN_SIG` value whole, as [the entry on the order of two variants' lists](clin-sig-allele-list-order.md#17-how-the-comparator-treats-this-class) describes, because two runs on one cache read one cell and must write one list. A comparison of a run on a native cache with Ensembl VEP's run on its own cache therefore reports such a row as a difference.

`tests/golden/116/GRCh38-native-builder` is the corpus that makes that comparison. One of its records takes an allele with two terms, which the two caches hold in the same order. None takes an allele on which they differ, so no committed comparison meets the difference; the reproduction above is the committed run that shows it.
