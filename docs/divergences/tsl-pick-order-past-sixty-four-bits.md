# Two support levels past 64 bits in the pick order

**The defect in one sentence.** Ensembl VEP 116.3's pick order compares two transcript support levels with Perl's `<=>`, which holds a level past 64 bits as a double and rounds a 64-bit integer within one step of it to the same double. The two levels tie, and the stable sort keeps whichever row was listed first. vep-rs orders the levels by their exact values and keeps the lower level whatever the row order.

**Why this is a defect.** The `tsl` category of the pick order ranks the candidates by their level as a number, the lower the better. A level that fits 64 bits is an integer to Perl and two such levels compare exactly. A level past 64 bits is a double, and a comparison between a double and an integer converts the integer to a double as well. Two integers one apart therefore compare unequal against each other yet equal against the same double, so the comparison is not an order. The row a sort over it leaves first depends on the order the rows entered the sort, not on the levels alone.

**What vep-rs does instead.** vep-rs orders two levels by their exact values: two levels within 64 bits as integers, a level past 64 bits after every integer, and two levels past 64 bits as doubles. The lower level wins whatever the row order, and the order is the same for every pair, so the pick is well defined. Wherever Ensembl VEP's comparison is consistent, every level within 64 bits or every level past them, the two engines pick the same row.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/OutputFactory.pm`; `ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/DBSQL/VariationFeatureAdaptor.pm`):

- `OutputFactory.pm` line 726: the `tsl` category defaults to the number 100.
- Lines 751-753, `pick_worst_VariationFeatureOverlapAllele`: the level is `$1` of `m/tsl(\d+)/` over the `TSL` attribute, kept as the matched string when it is true.
- Line 789: the candidates are sorted by `$a->{$cat} <=> $b->{$cat}`, a numeric comparison of those strings. Perl reads a digit run as a 64-bit integer when it fits and as a double otherwise, and compares an integer with a double by converting the integer to a double (perlop, `<=>`; perlnumber, numeric operators).
- Line 796: the rows that tie with the first are those whose value is `eq` to it, a comparison of the text.
- Lines 799 and 808: one tied row is the pick; after the last category the first row left is.
- `VariationFeatureAdaptor.pm` lines 2163, 2183, 2207 and 2214, `_pick_likely_transcript`: the same default, reading, sort and tie decide the transcript a gene name resolves to in an HGVS input.

**What correcting it would change.** Comparing the levels by exact value would make the lower level the pick in every row order. Every level within 64 bits, every level an Ensembl gene set carries among them, is unaffected either way.

**The class it produces.** `tsl_pick_order_past_sixty_four_bits`: the row kept under `--pick`, `--pick_allele`, `--per_gene`, `--pick_allele_gene` or their flagging variants with `tsl` in the pick order, or the transcript a gene name resolves to in an HGVS input. The shape is two candidates whose levels have twenty or more digits and lie within one rounding step of each other across 2^64, the double's row listed first. Ensembl VEP keeps the row listed first; vep-rs keeps the lower level. The class is documented only: no comparator codes it, since no committed flag set reaches it.

**How often the shape occurs.** Tier 3: the shape is a transcript support level of twenty or more digits, which no Ensembl gene set carries (`tsl1` to `tsl5` and `tslNA`) and which reaches a cache only through a GTF or GFF3 given to `vep-cache-builder` or through an edited cache. Of the 511 committed flag sets, 14 pass a pick-family flag with `tsl` in the order and none passes such a level; no committed HGVS input names a gene.

## 1. Two levels that straddle 2^64 on one gene

**In short.** The record `tsl_seam`, `21 42395417 T>C`, lies in two TMPRSS3 transcripts, ENST00001110244 and ENST00001113783, whose GTF lines carry the levels `18446744073709551616` (2^64) and `18446744073709551615` (one below it). Both engines print the two rows with each level as written. Under `--pick --pick_order tsl`, Ensembl VEP 116.3 keeps ENST00001110244, the row listed first; vep-rs keeps ENST00001113783, the lower level.

### 1.1 What a reader needs to know first

The pick family keeps one row per record, chosen by the categories of `--pick_order` in turn; the default order places `tsl` fifth, and `--pick_order tsl` makes it the only category. In the `tsl` category the lower level wins, a transcript without a level ranking at 100. The level is the digit run after `tsl` in the transcript's `TSL` attribute, kept as text ([the default and tab output](../output/vep-and-tab.md#extra-field)); a GTF's `transcript_support_level` attribute becomes that `TSL` attribute in both engines. The candidates enter the sort in Ensembl VEP's row order, the transcripts' stable ids with their database ids.

### 1.2 The record

The GTF holds the gene line of TMPRSS3 and the lines of ENST00001110244 and ENST00001113783 from the committed `tests/golden/116/GRCh38-gtf-builder/chr21_genes.gtf.gz`, each transcript line given a `transcript_support_level` attribute: `18446744073709551616` on ENST00001110244 and `18446744073709551615` on ENST00001113783. The input is one VCF record:

```
21	42395417	tsl_seam	T	C	.	.	.
```

Both engines run `--tsl`; the second run adds `--pick --pick_order tsl`. Ensembl VEP reads the GTF, bgzipped and indexed, through `--gtf` with the corpus's chromosome 21 FASTA; vep-rs reads a cache built from the GTF and that FASTA.

### 1.3 What each engine writes

Without the pick, both engines write the same two rows, the level in each as the GTF spells it; the first seven columns and the `Extra` column of vep-rs's rows:

```
tsl_seam	21:42395417	C	ENSG00000160183	ENST00001110244	Transcript	start_lost	IMPACT=HIGH;STRAND=-1;TSL=18446744073709551616
tsl_seam	21:42395417	C	ENSG00000160183	ENST00001113783	Transcript	start_lost	IMPACT=HIGH;STRAND=-1;TSL=18446744073709551615
```

Ensembl VEP's rows end in `SOURCE=` with the GTF's file name as well, the documented [`file_cache_writes_no_source_field`](file-cache-writes-no-source-field.md) difference; the position, amino-acid and codon columns left out above agree. Under `--pick --pick_order tsl`, Ensembl VEP 116.3 keeps the first row, ENST00001110244 with `TSL=18446744073709551616`; vep-rs keeps the second, ENST00001113783 with `TSL=18446744073709551615`.

### 1.4 Why Ensembl VEP writes that

`18446744073709551616` is one more than the largest 64-bit unsigned integer, so Perl holds it as the double 2^64. `18446744073709551615` fits, so Perl holds it as an integer. `<=>` between the two converts the integer to a double, and the nearest double to 2^64 minus 1 is 2^64 itself, so the two compare equal. The stable sort of line 789 keeps the equal rows in their listed order, ENST00001110244 first, and that row becomes the pick. The `eq` test of line 796 finds the two texts unequal, so the tied group is that row alone and the loop returns it. Listed the other way round, with the integer's transcript first, Ensembl VEP would keep the integer's row.

### 1.5 Why vep-rs's row is the correct one

The category asks for the lower level, and 18446744073709551615 is the lower of the two numbers written, so the row carrying it is the pick the category defines. vep-rs's order is the same for every pair of levels, an integer below every double, so the pick does not depend on which row the sort met first. Ensembl VEP's comparison gives one answer for these two rows in one order and the other answer in the other order, and places two integers one apart on both sides of the same double, which no ranking can do. Where no level passes 64 bits, both engines compare exactly and agree.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gtf tmprss3_tsl_seam.gtf --genome-fasta <corpus>/reference.fa.gz \
  --output-dir json_cache
vep -i record.vcf -o out.txt --offline --json_cache json_cache \
  --species homo_sapiens --assembly GRCh38 --tsl --pick --pick_order tsl
```

vep-rs writes the ENST00001113783 row. Ensembl VEP 116.3 over the same GTF, `--gtf tmprss3_tsl_seam.gtf.gz --fasta <corpus>/reference.fa.gz --tsl --pick --pick_order tsl`, writes the ENST00001110244 row. Without `--pick` both write both rows.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/tsl_pick_order_past_sixty_four_bits.toml` builds the GTF from the committed corpus, gives the two transcripts their levels, builds the cache and records both engines' picked rows. The class covers the row the pick family keeps, and the transcript a gene name resolves to, when two candidates carry levels past 64 bits within one rounding step of each other. No committed flag set carries such a level, so no comparator codes the class.
