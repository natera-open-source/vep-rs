# A `--custom` key named `MAX_AF` under `--max_af`

**The defect in one sentence.** Under `--max_af`, Ensembl VEP 116.3 writes neither the maximum frequency nor `MAX_AF_POPS` on a row where a `--custom` source's key is spelled `MAX_AF` and holds a value. The key keeps the source's value alone, as in a run without the flag. vep-rs writes the source's values and then the maximum under `MAX_AF`, and the maximum's populations under `MAX_AF_POPS`.

**Why this is a defect.** The command asks for a source's field and for the flag's output, and the row holds the first alone. A source named `MAX` that is asked for its `AF` field declares the key `MAX_AF`: its short name, an underscore and the field. The flag declares `MAX_AF` and `MAX_AF_POPS`. Ensembl VEP's header of the run describes the source's key and both keys of the flag ([1.3](#13-what-each-engine-writes)).

Ensembl VEP's row then carries nothing the flag asked for: neither the maximum nor `MAX_AF_POPS`, a key the source never names. No line of its frequency routine tests for a source's key, so nothing shows that the loss is intended. It follows from one comparison, which finds the source's list where it expects a number ([1.4](#14-why-ensembl-vep-writes-that)).

**What vep-rs does instead.** vep-rs lists the values the source's key carries on the row, then the row's maximum, joined as the output format joins a list: by a comma in the default and tab outputs and by `&` in a VCF. `MAX_AF_POPS` holds the maximum's populations. The maximum follows the source's values on every such row, also where one of them is the same text. So the last value of the cell is the maximum wherever the flag writes `MAX_AF_POPS` on the row, and the values before it are the source's.

On two kinds of row the engines agree. A row with a source value and no frequency holds the source's values alone and no `MAX_AF_POPS`. A row with a frequency and no record of the source holds the maximum alone, with its populations. The JSON output is the same in both engines: the frequencies stand under `colocated_variants` and the source's record under `custom_annotations`.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/OutputFactory.pm`):

- Lines 1297-1312, `_add_custom_annotations_to_hash`: a source's values are a list under its key. Line 1301 pushes a record's name under the short name, and line 1303 each field's value under the short name, an underscore and the field.
- Lines 1003-1014 and 1426-1435: those lists are filled before the frequency routine runs. The first range pushes a source record that belongs to the whole input line, as an `overlap` source's does. The second pushes a source record that names the row's allele, as an `exact` source's does. Lines 1441-1443 call `add_colocated_frequency_data` afterwards.
- Lines 1265-1275: the routine takes `$hash->{MAX_AF}` as the running maximum (`my $current_max = $hash->{MAX_AF} ||= 0`, line 1267). It assigns the new maximum when `$max_af > $current_max` holds (lines 1269-1272) and pushes the populations when `$max_af >= $current_max` holds (line 1274). With the source's list under the key, `$current_max` is a reference, which Perl reads in a comparison as a number far above 1, so both tests fail.
- Line 1247: `merge_arrays($hash->{$out_key} ||= [], [$f])` adds each frequency column's value to the list its key holds, where the list lacks it (`Utils.pm`, lines 303-310). A frequency column is a key such as `gnomADg_AF`, which `--af_gnomadg` prints.

**What correcting it would change.** A correction that pushes the maximum onto the key's list after the source's values, and the populations under `MAX_AF_POPS`, would have Ensembl VEP write vep-rs's row.

**The class it produces.** `custom_key_named_max_af`: under `--max_af`, a row on which a `--custom` key spelled `MAX_AF` holds a value and a co-located variant gives the row's allele a population frequency. Ensembl VEP writes the source's values under `MAX_AF` and no `MAX_AF_POPS`; vep-rs writes the source's values, the maximum and the populations. The class is documented only: no committed flag set names such a source, so no manifest names the class and no comparator applies it.

## 1. A source named `MAX` that is asked for its `AF` field

**In short.** One record of the GRCh38 custom corpus is run under `--max_af` with the corpus's `sites.vcf.gz` as an `exact` source named `MAX` that asks for `AF` and `CLASS`. The source's record gives the first alternate allele `AF=0.004`, and the cache gives that allele `1.993e-05` in one gnomAD population. Ensembl VEP 116.3 writes `MAX_AF=0.004` and no `MAX_AF_POPS`. vep-rs writes `MAX_AF=0.004,1.993e-05` and `MAX_AF_POPS=gnomADg_NFE`.

### 1.1 What a reader needs to know first

A `--custom` source prints each field it is asked for under a key made of its short name, an underscore and the field's name, and the names of its records under the short name alone. A source of `type=exact` gives a record to the row whose allele the record names. A source of `type=overlap` gives it to every row of an input line the record overlaps.

`--max_af` prints two keys. `MAX_AF` is the largest frequency a co-located known variant holds for the row's allele over the 1000 Genomes, ESP and gnomAD populations, and `MAX_AF_POPS` lists the populations that reach it. The combined gnomAD columns never count. A source named `MAX` that asks for `AF` therefore shares the key `MAX_AF` with the flag.

### 1.2 The record

The input is line 9 of `tests/golden/116/GRCh38-custom/variants.vcf`, a substitution with two alternate alleles:

```
21	5380280	rs1444031311	A	C,G	.	.	.
```

The source is `tests/golden/116/GRCh38-custom/custom/sites.vcf.gz`. Its one record at the site names the first alternate allele:

```
21	5380280	exact_2	A	C	.	PASS	AF=0.004;CLASS=exact
```

The known variant at the site, `rs1444031311`, is read from the cache of `tests/golden/116/GRCh38-everything`. It holds gnomAD genomes frequencies for the first alternate allele: `7.994e-06` in the combined column, `1.993e-05` in `gnomADg_NFE` and `0` in the nine other populations. It holds none for the second alternate allele. Both engines run:

```
--max_af --custom file=sites.vcf.gz,short_name=MAX,format=vcf,type=exact,fields=AF%CLASS
```

### 1.3 What each engine writes

Ensembl VEP 116.3 describes the three keys in the header of its default output:

```
## MAX_AF : Maximum observed allele frequency in 1000 Genomes, ESP and ExAC/gnomAD
## MAX_AF_POPS : Populations in which maximum allele frequency was observed
## MAX_AF : Allele frequency, one value per alternate allele. AF field from [PATH]/sites.vcf.gz
```

The `Extra` column of the first alternate allele's row in that output, from Ensembl VEP 116.3:

```
IMPACT=MODIFIER;MAX_AF=0.004;MAX=exact_2;MAX_CLASS=exact;MAX_FILTER=PASS
```

From vep-rs:

```
IMPACT=MODIFIER;MAX_AF=0.004,1.993e-05;MAX_AF_POPS=gnomADg_NFE;MAX=exact_2;MAX_CLASS=exact;MAX_FILTER=PASS
```

The second alternate allele's row carries `IMPACT=MODIFIER` alone in both engines. In the tab output vep-rs's two cells are `0.004,1.993e-05` and `gnomADg_NFE`, where Ensembl VEP's are `0.004` and `-`. In the VCF output vep-rs's `MAX_AF` field is `0.004&1.993e-05`. The JSON records of the two engines are equal.

Ensembl VEP lists the key twice in the tab and VCF outputs: two `MAX_AF` columns, and two `MAX_AF` fields of the `CSQ` format, each holding `0.004`. That is [`duplicate_output_column_collapsed`](duplicate-output-column.md), and vep-rs lists the key once.

Two controls on the same record show that the difference needs both the key and the flag. With the source named `sites`, both engines write `MAX_AF=1.993e-05;MAX_AF_POPS=gnomADg_NFE` beside `sites_AF=0.004`. Without `--max_af`, both write `MAX_AF=0.004` for the source named `MAX`.

Two more sources reach the same key. The tab cells `MAX_AF` and `MAX_AF_POPS` of each row:

| Source | Row | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `short_name=MAX_AF`, whose key lists the record's name | first alternate allele | `exact_2` and `-` | `exact_2,1.993e-05` and `gnomADg_NFE` |
| `short_name=MAX` read as `type=overlap` | first alternate allele | `0.004` and `-` | `0.004,1.993e-05` and `gnomADg_NFE` |
| the same source | second alternate allele | `0.004` and `-` | `0.004` and `-` |

The `overlap` record reaches both rows. The second alternate allele has no frequency, and the engines agree on its row.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP builds each row as a Perl hash, a table from an output key to its value. The source's record is added first: `MAX_AF` becomes a list holding `0.004`. The frequency routine then works out the maximum, `1.993e-05` in `gnomADg_NFE`. It compares that number with what the hash holds under `MAX_AF`, expecting the running maximum of the known variants read so far, and finds the list.

Perl holds a list in a hash as a reference, and a reference used as a number is its address in memory (Perl's `perlref` documentation), a number far above any frequency. The maximum is therefore never larger and never equal, so it is not assigned and its population is not added.

The source named `MAX_AF` shows that the list is compared and not the source's own number. Its key holds the text `exact_2`. Compared as a plain value, a text that is no number counts as 0, and the maximum would replace it. Ensembl VEP keeps `exact_2` and writes no population.

### 1.5 Why vep-rs's row is the correct one

Each value follows from the two files alone. The source's record gives the first alternate allele `AF=0.004` under a source named `MAX`, so the row owes `MAX_AF=0.004`. The cache line gives that allele `1.993e-05` in `gnomADg_NFE` and less in every other population, so `--max_af` owes that maximum and that population. Ensembl VEP writes exactly that maximum and that population when the source has another name. vep-rs's row holds the source's value, the maximum and the population, and no other value.

The order inside the cell is the one Ensembl VEP's own hash gives a key that a source and the frequency routine share: the source's values first, as they are added first, then the routine's value. Both engines write a frequency column in that order. Under `--af_gnomadg`, the same record with the source named `gnomADg` prints `gnomADg_AF=0.004,7.994e-06` in each: the source's value, then the cache's frequency.

A frequency column lists a frequency once where a source value is the same text (line 1247). That merge guards the routine's own repeats, a column that two of the frequency flags report or a second co-located variant, and states no rule for a source's value. vep-rs writes the maximum after the source's values even then, for example `1.993e-05,1.993e-05`. Listed once, the cell would hold `1.993e-05` alone, the same cell as on a row whose source record lacks the field. Whether it is listed once would also depend on how the source spells the number. By the comparison of 1.4, Ensembl VEP keeps the source's value alone on that row too.

A reader of the cell therefore needs no other rule. Where the flag writes `MAX_AF_POPS` on the row, the last value is the maximum; where it writes none, every value is the source's.

### 1.6 Try it yourself

With vep-rs, where `record.vcf` holds the record of 1.2 under a VCF header, `<cache>` is the `json_cache` directory of `tests/golden/116/GRCh38-everything` and `sites.vcf.gz` is the source of 1.2:

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 --max_af \
  --custom file=sites.vcf.gz,short_name=MAX,format=vcf,type=exact,fields=AF%CLASS
```

vep-rs writes the row of 1.3 with both values and the population. Ensembl VEP 116.3 takes the same `--max_af` and `--custom` options with `--offline --cache` and the release's GRCh38 cache, and writes `MAX_AF=0.004` and no `MAX_AF_POPS`. Naming the source so that no key of it is `MAX_AF` makes the engines agree.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/custom_key_named_max_af.toml` holds the record, the command and both engines' values. The class covers a row under `--max_af` on which a `--custom` key spelled `MAX_AF` holds a value beside a population frequency. A flag set is one committed run of a test corpus and a manifest lists that corpus's documented rows; [How to read the counts](../intended-divergences.md#how-to-read-the-counts) defines both and the comparators. No committed flag set names such a source, so no manifest names the class and the class is documented only.
