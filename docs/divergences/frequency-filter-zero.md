# The frequency filter's reading of a zero and of a matched allele

Two defects of Ensembl VEP 116.2's `--check_frequency` and `--filter_common` filter, both in the lines that read a known allele's frequency for the filter's test, and both corrected in Ensembl VEP 116.3. vep-rs applies the corrected readings against the 116.2 reference, so the records the filter keeps can differ from 116.2's. Under `--freq_gt_lt lt --freq_filter include`, Ensembl VEP 116.2 writes no row for a record whose known allele is stored with the frequency `0`, and vep-rs writes the record's rows. Under `--freq_gt_lt gt --freq_filter exclude`, the same holds for a minimised insertion matched to a rare known allele stored with its anchor base.

**The defects in one sentence each.** A population frequency of exactly `0` is read as no frequency, so the allele fails the test instead of being judged. An input allele's frequency is looked up under the allele's own spelling before the spelling of the known allele it matched, so a minimised insertion can read the known reference allele's filled-in frequency in place of its own.

**Why they are defects.** The filter judges an allele by the frequency of the known allele it matched in the chosen population. A cache value of `0` is an observed frequency, the allele absent from the population's samples, and the rarest allele a rarity filter exists to keep; reading it as missing inverts the filter for exactly those alleles. The matched allele is what `compare_existing` established the input allele to be, and the annotation side of the same run prints the matched allele's frequency in its `gnomADe_AF` column; a filter that reads another key judges a frequency the run's own output contradicts.

Ensembl VEP's 116.3 release corrects both in `BaseCacheVariation.pm`, the first as "Zero frequency should not be NA" and the second as "Use matched alleles to look up frequencies for --check_frequency".

**What vep-rs does instead.** vep-rs reads the matched known allele's frequency, falls back to the input allele's own key only when no match was recorded, and judges any frequency that is present, `0` included; a value that is absent stays `NA`. It also reads the matched allele on the input's strand, so a match on a reverse-strand known variant finds its frequency.

**Where they are in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm`, line 419: `my $f = $freq_data->{$alt} || 'NA';`. The `||` tests Perl truth, and the string `0` is false.
- The same file, lines 419-429: the input allele's own key is read first; the matched allele (`%matched_alleles`, lines 206-210) is consulted only when the own key read `NA`, and its value is tested for truth too.
- The same file, lines 249-252: the one allele of a known variant without a value is filled in as one minus the others' total, so a known reference allele carries a frequency near 1.
- Lines 325-334 and 341: under `include` a record whose alleles all fail is left out; under `exclude` it is kept.

**What correcting them would change.** Both classes vanish, as they do against 116.3: the records the filter keeps are the ones vep-rs keeps, and `FREQS` prints `0` where 116.2 prints `NA`.

**The classes they produce.** `frequency_filter_zero_value` and `frequency_filter_allele_key_order`, applied by the golden harness from a corpus manifest and by the format comparator from the manifest or, for the first, built in. The concordance figures run without the filter flags, so nothing is set aside in them.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A frequency of exactly zero

**In short.** `--check_frequency` with `--freq_filter include` keeps a record only when every allele's frequency passes the test. A known allele absent from a population is stored in the cache with the frequency `0`. Ensembl VEP 116.2 reads that `0` as `NA`, fails the allele and drops the record; vep-rs judges `0` against `--freq_freq` and keeps it. Under `--freq_gt_lt gt --freq_filter exclude` both engines keep the record and differ in `FREQS` alone, `NA` against `0`. Under `lt exclude` the readings swap: vep-rs's `0` passes the at-or-below test and the record is dropped, while Ensembl VEP's `NA` fails it and the record is kept.

### 1.1 What a reader needs to know first

`--check_frequency` compares each alternate allele's frequency in the population `--freq_pop` names with `--freq_freq`, under `--freq_gt_lt`. `--freq_filter exclude`, the default, drops the alleles that pass the test; `include` keeps only the records whose every allele passes. `--filter_common` is the exclude form at 1 percent. The frequency read is the matched known allele's, and `FREQS` prints the population, the allele and the value read. [Filtering by frequency](../output/known-variants.md#filtering-by-frequency) describes the flags.

### 1.2 The record

The record is `21 29025626 . A G`, record 4 of `variants.vcf` in the `GRCh38-freqfilter-include` golden corpus (`tests/golden/116/GRCh38-freqfilter-include/`):

```
21	29025626	.	A	G	.	.	.
```

The cache's known variant at the site is `rs1445461377`, `A/G`, with `gnomADe` `G:0` and `gnomADe_AFR` `G:0`. The corpus flag set `freqfilter-include` runs:

```
--check_frequency --freq_pop gnomADe_AFR --freq_freq 0.05 --freq_gt_lt lt --freq_filter include
```

### 1.3 What each engine writes

Under the flag set, Ensembl VEP 116.2 writes no row for the record in any format: the committed `expected_freqfilter-include/` outputs carry nothing for `29025626`. vep-rs writes the record's 47 rows, the first of them:

```
21_29025626_A/G	21:29025626	G	ENSG00000156256	ENST00000334352	Transcript	intron_variant	-	-	-	-	-	rs1445461377	IMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE_AFR:G:0
```

Record 17 of the corpus, `21 45914006 . T C` (`rs2096456303`, `gnomADe_AFR` `C:0`), behaves the same, 40 rows. Under `--freq_gt_lt gt --freq_filter exclude` with the same population and threshold both engines keep the record and print, with `--pick`:

```
21_29025626_A/G	21:29025626	G	ENSG00000156256	ENST00000399976	Transcript	intron_variant	-	-	-	-	-	rs1445461377	IMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE_AFR:G:NA
21_29025626_A/G	21:29025626	G	ENSG00000156256	ENST00000399976	Transcript	intron_variant	-	-	-	-	-	rs1445461377	IMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE_AFR:G:0
```

The first row is Ensembl VEP's and the second vep-rs's.

### 1.4 Why Ensembl VEP writes that

The filter reads the allele's frequency as `$freq_data->{$alt} || 'NA'`. Perl's `||` returns its right side when the left is false, and the string `0` is false, so a frequency of exactly `0` reads `NA`. The `NA` branch sets the allele's pass to 0 and looks for the matched allele's value, which is the same `0` and fails the same truth test, so the allele fails. Under `include` a record none of whose alleles passed is left out of the buffer before annotation; under `exclude` it is kept and `FREQS` prints the `NA` that was read. Ensembl VEP's 116.3 release replaces the `||` with `//`, defined-or, so that only an undefined value reads `NA`.

Source (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm`):

- line 419: the truth test on the frequency.
- lines 421-429: the `NA` branch and the matched allele's truth test.
- lines 325-334: the `include` drop; line 341: the `exclude` keep.

### 1.5 Why vep-rs's row is the correct one

The cache states the allele's frequency in gnomAD exomes' African population: `0`, an observed value. The allele was not seen among the population's alleles, which is the strongest rarity the data can state, and `0 < 0.05` holds, so the allele passes an at-or-below test and the record stays. Reading `0` as "no frequency" makes the rarest alleles the ones the filter cannot judge. Ensembl VEP 116.3 reads the record as vep-rs does.

### 1.6 Try it yourself

Run the corpus input through both engines under the flag set. vep-rs's command is below; Ensembl VEP 116.2 takes the same flags with the cache options of [Reproducing these examples](reproducing.md):

```
vep -i variants.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --check_frequency \
  --freq_pop gnomADe_AFR --freq_freq 0.05 --freq_gt_lt lt --freq_filter include
```

| Run | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- |
| records written under `include` | 10 | 12 (records 4 and 17 added) |
| `FREQS` of record 4 under `gt exclude` | `GNOMADE_AFR:G:NA` | `GNOMADE_AFR:G:0` |

### 1.7 How the comparator treats this class

The golden harness (`crates/vep-cli/tests/common/mod.rs`) applies `frequency_filter_zero_value` from the corpus manifest. The flag set's `vep_rs_only_tuples` document the 87 row keys of records 4 and 17 under the class, and the records' VCF line and JSON object are counted, not compared, as a record Ensembl VEP dropped before annotation is.

The format comparator (`scripts/concordance/compare_vep_formats.py`) reads the same manifest and also recognises the class without one. It accepts a `FREQS` whose items agree except where the reference reads `NA` and the candidate `0`, and the rows of a record the candidate alone writes whose every `FREQS` item is `0`, with their record entry. Each accepted row is counted once under the class; a record kept on a frequency other than `0`, or a `FREQS` that differs in any other way, stays a difference. The record the reference alone keeps under `lt exclude` has no built-in shape: its rows count as missing from the candidate unless a corpus manifest documents them under the class.

## 2. The input allele's own key before the matched allele's

**In short.** An input allele and the known allele it matched can be spelled differently: a minimised insertion `-/A` matches a known variant stored anchored as `A/AA`. Ensembl VEP 116.2 looks the frequency up under the input's spelling first, and `A` is also the known variant's reference allele, whose frequency the cache filled in as nearly 1; so the filter judges a common allele where the input is a rare one. vep-rs judges the matched allele's frequency. Under `--freq_gt_lt gt --freq_filter exclude` Ensembl VEP writes no row for the record and vep-rs writes its row; under `include` the reverse.

### 2.1 What a reader needs to know first

The cache stores a known variant's alleles as the variant was reported, an insertion with its anchor base: `A/AA`. A VCF insertion `A AA` is minimised to `-/A` by the parser, and `compare_existing` records that the input's `A` is the known `AA`. The cache's frequency table lists a value per allele, and the one allele without a value is filled in as one minus the others, so the reference allele of a biallelic variant carries a frequency near 1. The same run's `--af_gnomade` column resolves the input allele through the match and prints the matched allele's frequency.

### 2.2 The record

A public GRCh38 record on chromosome 21, run against the cache of the `GRCh38-recoder` golden corpus (`tests/golden/116/GRCh38-recoder/json_cache`), which carries the variation rows of its region:

```
21	33543210	.	A	AA	.	.	.
```

The known variant at the site is `rs572248367`, `A/AA`, with `gnomADe` `AA:4.959e-05`. Both engines run:

```
--check_frequency --freq_pop gnomADe --freq_freq 0.01 --freq_gt_lt gt --freq_filter exclude --pick
```

and the same with `--freq_filter include`.

### 2.3 What each engine writes

Under `exclude`, Ensembl VEP 116.2 writes no row for the record; vep-rs writes:

```
21_33543210_A/AA	21:33543210-33543211	A	ENSG00000159140	ENST00000356577	Transcript	intron_variant	-	-	-	-	-	rs572248367	IMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE:A:4.959e-05
```

Under `include` the readings swap: Ensembl VEP writes

```
21_33543210_A/AA	21:33543210-33543211	A	ENSG00000159140	ENST00000356577	Transcript	intron_variant	-	-	-	-	-	rs572248367	IMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE:A:0.99995041
```

and vep-rs writes no row. The record `21 33543544 . C CC` (`rs11375188`, `C/CC`, `gnomADe` `CC:1`) run alongside is kept by both engines under `include` with `FREQS=GNOMADE:C:1` and dropped by both under `exclude`.

### 2.4 Why Ensembl VEP writes that

The filter reads `$freq_data->{$alt}` with `$alt` the input's allele, `A`. The known variant's frequency table has an entry under `A`, the reference allele, filled in at lines 249-252 as `1 - 4.959e-05`. The own-key read therefore succeeds with `0.99995041` and the matched allele is never consulted; the match is read only when the own key gives `NA`. The frequency judged is the reference allele's, far above 1 percent, so `exclude` drops the record and `include` keeps it.

For the control record the own key `C` reads the filled-in `0`, which the truth test turns into `NA`, and the filter falls back to the matched `CC` and reads `1`, so both engines agree on it. Ensembl VEP's 116.3 release reads the matched allele first and the own key only where no match was recorded.

Source (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm`):

- lines 419-429: the own key first, the match on a miss.
- lines 206-210: the matched alleles recorded by `compare_existing`.
- lines 249-252: the filled-in frequency of the allele without a value.

### 2.5 Why vep-rs's row is the correct one

`compare_existing` recorded that the input's `A` is the known variant's `AA`, and the cache gives that allele's frequency as `4.959e-05`, about one exome in twenty thousand. The own-key hit is a coincidence of spelling between the input's minimised allele and the known variant's reference allele, and the reference allele's filled-in frequency describes the alleles the input does not carry. Ensembl VEP's own `gnomADe_AF` column prints `4.959e-05` for the record in both engines, so 116.2's filter judges a value its own output denies. Ensembl VEP 116.3 reads the record as vep-rs does.

### 2.6 Try it yourself

Annotate the two records of 2.2 and 2.3 under both filter modes:

```
vep -i records.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --check_frequency --freq_pop gnomADe \
  --freq_freq 0.01 --freq_gt_lt gt --freq_filter exclude --pick
```

| Record and mode | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- |
| `21 33543210 A AA`, `exclude` | no row | kept, `FREQS=GNOMADE:A:4.959e-05` |
| `21 33543210 A AA`, `include` | kept, `FREQS=GNOMADE:A:0.99995041` | no row |
| `21 33543544 C CC`, either mode | as vep-rs | kept under `include`, dropped under `exclude` |

### 2.7 How the comparator treats this class

The golden harness and the format comparator apply `frequency_filter_allele_key_order` from a corpus manifest. The rows of a record one engine keeps and the other drops are documented as that engine's under the class, and a record the reference dropped is counted, not compared, at the record level as under class 1. No predicate decides the class from the two outputs, since the rows do not state which key the frequency was read under. No committed corpus carries an insertion input against an anchored known allele, so the class applies to none of their rows.
