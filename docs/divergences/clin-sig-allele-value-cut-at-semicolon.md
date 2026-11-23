# A ClinVar significance that holds several classifications

**The defect in one sentence.** ClinVar can give one record of a variant several classifications, which Ensembl's variation database stores as one value, such as `uncertain significance; risk factor`. When Ensembl VEP 116.3 reads that database, it joins each record's attributes with a semicolon and cuts the joined string at every semicolon. Such a value keeps its first classification alone, so the known variant's `clin_sig_allele` has no `allele:term` pair for a later one, in a cache Ensembl VEP dumps and in a run under `--database`. A cache that `vep-cache-builder` builds from the same database's tables holds a pair for each classification, so a row's `CLIN_SIG` can list a term there that it lacks on Ensembl VEP's cache.

**Why this is a defect.** The cut follows no rule about the classification. A record whose whole value is `risk factor` keeps its pair, and the same classification after a semicolon loses it. Ensembl's own import reads such a value as two classifications when it fills the variant's `clinical_significance`, so the variant-level `clin_sig` of Ensembl VEP's cache holds the term that the allele's list lacks. For a variant with one alternate allele, two runs of one record then disagree: `--clin_sig_allele 0` prints the term and the default does not. ClinVar's own VCF gives such a variant several `CLNSIG` values.

**What vep-rs does instead.** `vep` writes a row's `CLIN_SIG` from the cache's cell, as Ensembl VEP does, so on one cache the two engines write the same list. `vep-cache-builder` cuts a record's significance at its semicolons and writes one `allele:term` pair for each classification, each different pair once.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`; `ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/` and `scripts/import/`):

- `DBSQL/PhenotypeFeatureAdaptor.pm`, `get_clinsig_alleles_by_location`, lines 1233-1236: the query joins a record's `code=value` attributes with a semicolon and a blank.
- The same routine, lines 1273-1276: the joined string is cut at every `;` and each piece at `=`, so the text after a value's semicolon becomes a key with no value, which nothing reads.
- `AnnotationSource/Database/Variation.pm`, lines 203-213: each record's `risk_allele` and `clinvar_clin_sig` make one pair, and the cell holds each different pair once.
- `Pipeline/DumpVEP/Dumper/Variation.pm`, lines 207-209 and 247: the cache dump calls that routine and writes the cell as it stands.
- `OutputFactory.pm`, lines 1105-1117: a row's `CLIN_SIG` joins the terms the cell gives the row's allele. Line 1125 writes `clin_sig` under `--clin_sig_allele 0`.
- `scripts/import/clinvar_post_process.pl`, line 233: the import turns the semicolon into a comma for the variant's `clinical_significance`.

**What correcting it would change.** Cutting each record's significance into its classifications before the pairs are made would write a pair for each: the cell `vep-cache-builder` writes. Joining the attributes with a separator that no value holds would not be enough. The value would then reach the cell whole, and the cell joins its pairs with `;`, where the output cuts it: the later classification would be a piece that names no allele.

**The class it produces.** The entry's name is `clin_sig_allele_value_cut_at_semicolon`. It stands for the `allele:term` pairs of one known variant's `clin_sig_allele`, between a cache Ensembl VEP's dump wrote, or a `--database` run, and a cache `vep-cache-builder` built from the same release's tables. It shows where a ClinVar record's significance holds two or more classifications, and a later one stands in no record of the allele as the whole value or as its first classification. It shows in `CLIN_SIG` under `--check_existing` without `--clin_sig_allele 0`, and in the `clin_sig_allele` string of the JSON format.

No comparison names a class for it. Two runs on one cache read one cell and must write one list, and a rule that set a longer list aside would hide a difference between them. So a comparison reports a list with a term more as a difference.

**How often the shape occurs.** Tier 3 of [the index's order](../intended-divergences.md#differences-outside-the-classes): no committed table holds the shape. The shape is a known variant, in a cache `vep-cache-builder` built with the ClinVar tables, one of whose ClinVar records holds a significance of two or more classifications. No row of the committed native cache, `tests/golden/116/GRCh38-native-builder/json_cache`, holds it, and one committed run of that corpus, its flag set `variation`, runs under `--check_existing`. The table counts the shape in the ClinVar records of Ensembl's human variation database of release 116 on GRCh38, over chromosomes 1 to 22, X, Y and MT.

| Counted | Number |
| --- | --- |
| ClinVar records whose significance holds a semicolon | 127 |
| Of them, germline records that name an allele | 120 |
| Variants those 120 records belong to | 107 |
| Variants with a classification that no record gives the allele alone or first | 94 |
| Such pairs of an allele and a classification | 98 |

For each of the 107 variants, the cell of Ensembl VEP's GRCh38 cache of release 116 holds exactly the pairs of the values cut at their first semicolon. It lacks all 98 pairs of the 94 variants, and the variant-level `clin_sig` of the same row holds the term of each.

The later classification is `other` 65 times, `risk factor` 19, `association` 10 and `affects` 4. Of the 107 variants that carry the shape, the difference follows for the 94 of the table's fourth row. For 96 of the 98 pairs the allele is one of the variant's own alleles, so a row can show the term. The other two records name an allele the variant does not have, and no row shows theirs.

## 1. One of two classifications is lost

**In short.** The known variant `rs74315447` at `21:34370639` has eight ClinVar records for its one alternate allele. One reads `uncertain significance; risk factor`, and no other gives `risk factor`. For the record `21 34370639 T C` under `--check_existing`, vep-rs on a native cache prints `CLIN_SIG=uncertain_significance,risk_factor,likely_pathogenic,not_provided,likely_benign`. Ensembl VEP 116.3 on its own cache prints `CLIN_SIG=likely_benign,uncertain_significance,not_provided,likely_pathogenic`, four terms and no `risk_factor`.

### 1.1 What a reader needs to know first

`--check_existing` names the known variants at a record's position. A known variant can carry `clin_sig_allele`, ClinVar's terms by allele, written `allele:term` and joined by `;`. Under the default `--clin_sig_allele 1` a row's `CLIN_SIG` reports the terms a known variant records for the row's own allele ([Co-located variant keys](../output/known-variants.md#co-located-variant-keys)). Under `--clin_sig_allele 0` it reports the variant's `clin_sig`, the terms of all its alleles together.

ClinVar classifies a variant for a condition, and one such record can hold more than one classification: a pathogenicity classification such as `uncertain significance`, and beside it a category such as `risk factor`, `association`, `affects`, `drug response` or `other`. ClinVar's own VCF, which gives each variant one aggregate, writes such classifications as several values of `CLNSIG`, separated by a vertical bar. Ensembl's variation database stores them in one `clinvar_clin_sig` value, joined by a semicolon and a blank. A slash or a comma inside a `clinvar_clin_sig` value belongs to one classification (`benign/likely benign`, `pathogenic, low penetrance`), and both engines keep such a value whole.

Ensembl VEP's cache is a dump of its variation database. A native cache is one `vep-cache-builder` builds from Ensembl's variation VCFs and, for this cell, from the release's `phenotype_feature` and `phenotype_feature_attrib` tables ([Variation data in a native cache](../cache-variation-data.md)).

### 1.2 The record

The variation database of release 116 holds eight germline ClinVar records for `rs74315447`, all for its alternate allele `C`. By identifier:

| `phenotype_feature_id` | Allele | `clinvar_clin_sig` |
| --- | --- | --- |
| 299470652 | `C` | uncertain significance; risk factor |
| 299470960 | `C` | likely pathogenic |
| 299471037 | `C` | not provided |
| 299471368 | `C` | uncertain significance |
| 299474509 | `C` | likely pathogenic |
| 299476968 | `C` | likely benign |
| 299479825 | `C` | uncertain significance |
| 299480734 | `C` | uncertain significance |

The two caches hold these cells for the variant, a term's blanks written as underscores:

| Cache and cell | Value |
| --- | --- |
| Ensembl VEP's own, `clin_sig_allele` | `C:likely_benign;C:uncertain_significance;C:not_provided;C:likely_pathogenic` |
| A native cache, `clin_sig_allele` | `C:uncertain_significance;C:risk_factor;C:likely_pathogenic;C:not_provided;C:likely_benign` |
| Both caches, `clin_sig` | `uncertain_significance,not_provided,likely_benign,likely_pathogenic,risk_factor` |

The input is one VCF record with that allele:

```
21	34370639	kcne2_two_classifications	T	C	.	.	.
```

Both engines take the flag `--check_existing`.

### 1.3 What each engine writes

Each engine writes one `CLIN_SIG` value on every row of the record. Ensembl VEP writes 34 rows, one for each transcript its cache holds in reach of the position. The native cache of this example is built with the gene model of the builder's tests, a slice of chromosome 21 that holds none there, so vep-rs writes one row.

| Run | `CLIN_SIG` |
| --- | --- |
| Ensembl VEP 116.3 on its own cache | `likely_benign,uncertain_significance,not_provided,likely_pathogenic` |
| vep-rs on the native cache | `uncertain_significance,risk_factor,likely_pathogenic,not_provided,likely_benign` |
| Ensembl VEP 116.3 on its own cache, `--clin_sig_allele 0` | `uncertain_significance,not_provided,likely_benign,likely_pathogenic,risk_factor` |
| Ensembl VEP 116.3 under `--database` | `likely_pathogenic,likely_benign,not_provided,uncertain_significance` |

In the JSON format the known variant's `clin_sig_allele` is the cell as each cache holds it, the first two strings of the second table of 1.2. On Ensembl VEP's cache [converted for vep-rs](../cache-conversion.md), vep-rs prints Ensembl VEP's four terms, since it then reads Ensembl VEP's cell. The order of the terms differs between the runs for [a reason of its own](clin-sig-allele-term-order.md).

A second variant shows that the cut drops a classification only where it shares a record. `rs1805128` at `21:34449382` has ten records for its alternate allele. One reads `risk factor` alone, one `benign/likely benign; other` and two `likely benign; risk factor`. Ensembl VEP's cell holds `T:risk_factor`, from the first of these, and no `T:other`, while its `clin_sig` holds `other`. The native cell holds both.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP reads a region's ClinVar records with one query. For each record the query returns one string: the variant's name, the record's identifier and the record's attributes, each written `code=value`, all joined by a semicolon and a blank. The routine then cuts the string at every semicolon and each piece at its `=`, and files the value under the code.

A value that itself holds a semicolon is cut with the rest. The record's `clinvar_clin_sig` becomes the text before the semicolon, `uncertain significance`. The text after it, `risk factor`, holds no `=`, so it becomes a code with no value, which no later step reads. The cell is then built from each record's allele and its `clinvar_clin_sig`, so the second classification never reaches it. The cache dump and a `--database` run both call the routine.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`; `ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/`):

- `DBSQL/PhenotypeFeatureAdaptor.pm`, lines 1233-1236: the query's string, joined with `'; '`.
- `DBSQL/PhenotypeFeatureAdaptor.pm`, lines 1273-1276: the cut at every `;` and at `=`.
- `AnnotationSource/Database/Variation.pm`, line 161: the region's ClinVar records; lines 203-213: one pair for each record, joined with `;`.
- `Pipeline/DumpVEP/Dumper/Variation.pm`, lines 207-209 and 247: the dump's call and the cell written.

### 1.5 Why vep-rs's row is the correct one

The eight records give the allele five classifications: `uncertain significance`, `risk factor`, `likely pathogenic`, `not provided` and `likely benign`. vep-rs's row lists the five, each once. Ensembl VEP's row lists four. The one it leaves out is in the database, in a record of this allele, and the variant has no other alternate allele it could belong to.

Ensembl's own data and code read the value as two classifications. The import that fills the variant's `clinical_significance` replaces the semicolon with a comma (`scripts/import/clinvar_post_process.pl`, line 233, whose comment gives `uncertain significance; risk factor` as its example). So the `clin_sig` cell of Ensembl VEP's own cache holds `risk_factor` for the variant, and Ensembl VEP prints it under `--clin_sig_allele 0`. ClinVar's own VCF gives the same variant the `CLNSIG` values `Conflicting_classifications_of_pathogenicity` and `risk_factor`.

The cut is not a choice to keep pathogenicity classifications alone. Ensembl VEP's cell keeps `risk_factor`, `other` and the like wherever a record holds one as its whole value, as the second variant of 1.3 shows. It loses the classification only where the database joined it to another in one value.

### 1.6 Try it yourself

Tests of two files hold vep-rs's side. One test of the first file builds a native cache from ClinVar tables whose values hold one, two and three classifications, and holds each cell. The three tests of the second file build the variation file of the two variants of 1.3 from the release's own rows, which the repository keeps under `crates/vep-cli/tests/fixtures/clinvar-several-classifications`, and hold the row of 1.3 and both cells:

```
cargo test -p vep-cache-builder --test variation_cache_from_dumps
cargo test -p vep-cli \
  --test clinvar_significance_of_several_classifications_on_a_native_cache
```

To see vep-rs's row, build a native cache from those rows. `bgzip` compresses the variation VCF, which the builder reads compressed. One directory holds that file and the tables, so both variation options name it. The gene model is a slice of chromosome 21 that the repository keeps for the builder's tests, with the protein sequences of its transcripts:

```
F=crates/vep-cli/tests/fixtures/clinvar-several-classifications
G=crates/vep-cache-builder/tests/fixtures/grch38_chr21
cp -r $F variation && bgzip variation/homo_sapiens-chr21.vcf
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 $G/slice.gff3 --fasta $G/peptides.fa --include-variations \
  --variation-vcf-dir variation --variation-dumps variation \
  --output-dir native_cache
```

Save the record of 1.2 as `record.vcf` and run vep-rs's `vep` on that cache:

```
vep -i record.vcf -o out.txt --offline --json_cache native_cache \
  --species homo_sapiens --assembly GRCh38 --format vcf --check_existing
```

For Ensembl VEP 116.3, run its `vep` on the same record and the release's own GRCh38 cache, with the cache options of [Reproducing these examples](reproducing.md). Add `--clin_sig_allele 0` for the third row of 1.3:

```
vep -i record.vcf -o out.txt --offline --cache --dir_cache <cache> \
  --cache_version 116 --species homo_sapiens --assembly GRCh38 \
  --format vcf --check_existing
```

The run of the fourth row of 1.3 reaches Ensembl's public database:

```
vep -i record.vcf -o out.txt --database --db_version 116 \
  --species homo_sapiens --assembly GRCh38 --format vcf --check_existing
```

### 1.7 How a comparison treats such a row

The reproduction `scripts/divergences/reproductions/clin_sig_allele_value_cut_at_semicolon.toml` runs the record of 1.2 on Ensembl VEP with its own cache and on vep-rs with a native cache that its prepare steps build as 1.6 does. It records the known variant's `clin_sig_allele` string of the JSON format, four pairs against five.

No comparison sets the difference aside. The golden harness and the format comparator compare each list of a `CLIN_SIG` value whole, because two runs on one cache read one cell and must write one list. A comparison of a run on a native cache with Ensembl VEP's run on its own cache therefore reports such a row as a difference.

`tests/golden/116/GRCh38-native-builder` is the corpus that makes that comparison. Its native cache holds no variant with such a ClinVar record, so no committed comparison meets the difference; the reproduction above is the committed run that shows it.
