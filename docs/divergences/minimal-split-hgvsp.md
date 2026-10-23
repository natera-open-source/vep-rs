# The HGVSp, amino acids and positions of a `--minimal` split allele

A difference in the `HGVSp`, `Amino_acids`, `CDS_position` and `Protein_position` of an allele that `--minimal` splits off a multi-allelic record when its HGVS notations shift. The row's terms and positions differ under `--shift_3prime 1` as well, as do the `HGVSc` of a row `--coding_only` restores and the `miRNA` terms of a split allele on a miRNA transcript. So do the `SPDI` and `GA4GH_VRS` of a split allele and the `PICK` marks of the flagging variants; nothing is set aside, since no measured dataset runs the flags.

Under `--minimal`, Ensembl VEP splits a multi-allelic record whose alleles trim into one feature per allele (`InputBuffer.pm` lines 638-718) and merges them back at output (`OutputFactory.pm` lines 2272-2354). The rejoin builds each feature's rows first (lines 2289 and 2303). That first pass shifts a deletion's or insertion's HGVS notations as it does a bi-allelic record's, keeping the HGVSc and its offset (`TranscriptVariationAllele.pm` lines 1558 and 1838) and the peptides read at the shifted position (lines 1774-1778); the HGVSp itself is not kept. Both engines print that HGVSc and offset for an allele whose first pass built its row.

The default and tab factories then clear the cached coordinates (`OutputFactory.pm` lines 306-323) and build the rows again over the record's restored span. The HGVSp is recomputed from that span's translation coordinates (`TranscriptVariationAllele.pm` lines 1724-1737), so a record whose span reaches into an intron has none, and `Amino_acids` prints the peptides the first pass left (lines 618-630). The VCF and JSON factories clear the coordinates before the rejoin, so their rows keep the first pass's HGVSp, its peptides and the `CDS_position` and `Protein_position` the shift cached, beside an unshifted `cDNA_position`. vep-rs builds a split allele's row once, over the allele's own span, in every format.

The `miRNA` column reads the same cleared coordinates (`OutputFactory.pm` line 1630 takes `cdna_start` and `cdna_end`, which `clear_shifting_variables` deleted, `TranscriptVariationAllele.pm` lines 515-517). So the default and tab factories compute a split allele's hairpin terms over the record's restored span, and the VCF and JSON factories over the allele's own. For `21 8987003 multi ACT A,ACTTTA` under `--mirna --minimal`, the `TTA` allele's row on the miRNA transcript prints `miRNA=miRNA_loop` in Ensembl VEP's default output and `miRNA=miRNA_loop,miRNA_stem` here; the VCF and JSON of both engines print `miRNA_loop&miRNA_stem`.

The records below sit at the acceptor of exon 6 of CLIC6, the deletion allele `GATT/-` straddling the boundary and shifting one base into the exon for its notations:

```
21	34716319	multi_del_two	AGATT	AGATC,A	.	.	.
21	34716319	bi_del	AGATT	A	.	.	.
```

Under `--minimal --hgvs`, the two engines write for the deletion allele of `multi_del_two` on ENST00000349499, which both describe as `c.1900_1903del` with `HGVS_OFFSET=1`:

| Format | Engine | cDNA_position, CDS_position, Protein_position | Amino_acids and HGVSp |
| --- | --- | --- | --- |
| default, tab | Ensembl VEP | `?-2297`, `?-1902`, `?-634` | `IV/X`, none |
| default, tab | vep-rs | `?-2297`, `?-1902`, `?-634` | none, `p.Ile634TrpfsTer14` |
| VCF, JSON | Ensembl VEP | `?-2297`, `1900-1903`, `634-635` | `IV/X`, `p.Ile634TrpfsTer14` |
| VCF, JSON | vep-rs | `?-2297`, `?-1902`, `?-634` | none, `p.Ile634TrpfsTer14` |

The record's `C` allele, a synonymous substitution once trimmed, prints `p.Ile634=` in vep-rs and in Ensembl VEP's VCF and JSON, and no HGVSp in Ensembl VEP's default and tab outputs, whose second pass reads the record's span. The bi-allelic `bi_del` agrees in every format and prints no amino acids, since its one pass computes them before the HGVS shift.

Under `--shift_3prime 1` the second pass also rebuilds the deletion allele's terms and positions. Ensembl VEP's default output prints the shifted span's positions, `2295-2298`, `1900-1903` and `634-635`, no codons and, without `--hgvs`, no amino acids. Its terms are `splice_region_variant,start_retained_variant` with `IMPACT=LOW`: the `frameshift_variant` of the bi-allelic `bi_del` is lost, and `start_retained_variant` is [the start-codon test's](shift-start-codon-test.md). vep-rs prints the positions its default output gives a split allele, the record span's `?-2297`, `?-1902` and `?-634`, with the amino acids `IV/X`, the shifted codons `ATTGtg/tg` and `frameshift_variant,splice_region_variant`; the HGVS notations and their offset agree.

A row `--coding_only` restores has its HGVSc rebuilt as well. The `--coding_only` test runs over the record's span after the rejoin ([`--minimal` under `--coding_only`](../output-formats.md#parity-with-ensembl-vep)), so a row whose feature wrote nothing before the rejoin has no cached HGVSc. The second pass computes it over the restored span (`TranscriptVariationAllele.pm` line 1558 caches the first pass's; `OutputFactory.pm` lines 2291-2293 restore the record). The `G` allele of `21 34716318 CAGATT GAGATT,C` is an intronic substitution three bases before the exon's first base. Under `--minimal --coding_only --hgvs` Ensembl VEP prints the whole record as a deletion-insertion, `c.1900-3_1902delinsG` on ENST00000349499 and the same form in each other coding transcript's numbering; vep-rs prints the allele's own `c.1900-3C>G`, which both engines print without `--coding_only`.

Under `--coding_only` with `--minimal`, a trimmed allele whose own span lies in the UTR keeps its own consequence terms and codons. Ensembl VEP rebuilds the restored row's terms over the record's span: `21 31655249 AT ATG,A` on GRCh37 against ENST00000340345 prints `5_prime_UTR_variant` here and `start_lost M/L atg/Ctg` there.

## SPDI and GA4GH VRS of a split allele

`SPDI` and `GA4GH_VRS` are looked up by allele after the rejoin. The rejoin restores the record's own allele string (`OutputFactory.pm` lines 2291-2293), and the SPDI table is built from that string, one notation per allele of the record (line 1398). Each row's notation is then looked up by the row's `Allele`, the split allele's trimmed string (lines 1347 and 1400), and the VRS object is built from whatever that lookup returns (lines 1403-1410). A trimmed allele that is no allele of the record finds nothing, and one that equals another allele of the record finds that allele's notation. vep-rs builds each split allele's notation from the allele itself.

The records below are two multi-allelic records whose alleles trim, the second anchor-chopped to `TT/-/T`:

```
21	14144642	.	AT	CT,A	.	.	.
21	38218615	.	GTT	G,GT	.	.	.
```

Under `--spdi --minimal` the two engines write, on every transcript row of the allele:

| Allele | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `CT`, trimmed `C` | no `SPDI` | `21:14144641:A:C` |
| `A`, trimmed `-` | no `SPDI` | `21:14144642:T:` |
| `G`, the deletion `-` | `21:38218615:TT:` | `21:38218615:TT:` |
| `GT`, trimmed `-` | `21:38218615:TT:`, the `G` allele's | `21:38218616:T:` |

Ensembl VEP writes `GA4GH_VRS` with a database connection alone, since the `NC_` accession the object needs comes from the slice adaptor. Under `--ga4gh_vrs` the object follows the notation: none where Ensembl VEP finds no `SPDI`, the `G` allele's object on the `GT` rows, and the allele's own object in vep-rs. Without `--minimal` the two engines agree on every row.

## PICK marks under the flagging variants

The first pass marks rows too. Each of its calls (`OutputFactory.pm` lines 2289 and 2303) hands one split feature's alleles to `filter_VariationFeatureOverlapAlleles` (line 514). A split feature holds one alternate allele, so the flagging branches (lines 609-626) set `PICK` on one row of it: one row per feature under `--flag_pick` and `--flag_pick_allele`, one per gene of the feature under `--flag_pick_allele_gene`. The merge copies the marked rows under new keys (line 2316) and clears the consequence caches alone (line 2342), so the marks stay.

The output pass then filters the rejoined record once more (line 301), grouping its rows by allele string under the per-allele flags (lines 618 and 625). It marks a row of each group, a row the first pass had marked on every record below.

So Ensembl VEP marks one row per split feature. Under `--flag_pick` a record `--minimal` splits carries one mark per split allele, where a record it does not split carries one. Under `--flag_pick_allele` and `--flag_pick_allele_gene` the marks are doubled only where two alternate alleles trim to the same allele string: such a record carries two marks per allele string, or two per allele string and gene, in place of one. vep-rs runs the pick family once over the rejoined record, grouping the rows by the Allele column, so it marks one row per record under `--flag_pick` and one per allele string, or per allele string and gene, under the other two.

The records below are two insertions after one base that both trim to `T`, an SNV beside an insertion that both trim to `A`, and two insertions that trim to the distinct strings `T` and `A`:

```
21	43100486	min_gt_tg	G	GT,TG	.	.	.
21	17604928	min_a_ga	G	A,GA	.	.	.
21	43100486	min_gt_ag	G	GT,AG	.	.	.
```

Under `--minimal --allele_number`, which numbers each row's original allele, the two engines mark these rows; every other row of each run agrees, 76 rows for the records at 43100486 and 72 for the record at 17604928:

| Record and flag | Engine | Rows marked `PICK=1` (Allele, Feature, ALLELE_NUM) |
| --- | --- | --- |
| `min_gt_tg`, `--flag_pick_allele` | Ensembl VEP | `T` ENST00000291552 at 1 and at 2 |
| `min_gt_tg`, `--flag_pick_allele` | vep-rs | `T` ENST00000291552 at 2 |
| `min_a_ga`, `--flag_pick_allele_gene` | Ensembl VEP | `A` ENST00000348354 and ENST00001128578, each at 1 and at 2 |
| `min_a_ga`, `--flag_pick_allele_gene` | vep-rs | `A` ENST00000348354 and ENST00001128578, each at 2 |
| `min_gt_ag`, `--flag_pick` | Ensembl VEP | `A` ENST00000291552 at 2 and `T` ENST00000291552 at 1 |
| `min_gt_ag`, `--flag_pick` | vep-rs | `A` ENST00000291552 at 2 |

`min_gt_ag` under `--flag_pick_allele` and `21 17604928 G A,GC` under `--flag_pick_allele_gene`, whose alleles trim to distinct strings, mark the same rows in both engines: each allele group is one split feature, so the two passes mark one row. Without `--minimal` both engines mark the `GT` and `TG` rows of ENST00000291552 for `min_gt_tg` under `--flag_pick_allele`.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--minimal --hgvs`, and `--tab`, `--vcf` or `--json` for the other formats. `--shift_3prime 1` gives the shifted rows and `--coding_only` the restored row; the PICK marks take `--minimal --allele_number` with one of the three flagging variants and need no FASTA. The SPDI rows take `--minimal --spdi`, and the VRS objects `--minimal --ga4gh_vrs` with `--cache --db_version 116` in place of `--offline`. vep-rs takes the same flags over the same cache and FASTA, `--offline` throughout.
