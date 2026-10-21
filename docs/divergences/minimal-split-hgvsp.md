# The HGVSp, amino acids and positions of a `--minimal` split allele

A difference in the `HGVSp`, `Amino_acids`, `CDS_position` and `Protein_position` of an allele that `--minimal` splits off a multi-allelic record when its HGVS notations shift. The row's terms and positions differ under `--shift_3prime 1` as well, and so does the `HGVSc` of a row `--coding_only` restores; nothing is set aside, since no measured dataset runs the flags.

Under `--minimal`, Ensembl VEP splits a multi-allelic record whose alleles trim into one feature per allele (`InputBuffer.pm` lines 638-718) and merges them back at output (`OutputFactory.pm` lines 2272-2354). The rejoin builds each feature's rows first (lines 2289 and 2303). That first pass shifts a deletion's or insertion's HGVS notations as it does a bi-allelic record's, keeping the HGVSc and its offset (`TranscriptVariationAllele.pm` lines 1558 and 1838) and the peptides read at the shifted position (lines 1774-1778); the HGVSp itself is not kept. Both engines print that HGVSc and offset for an allele whose first pass built its row.

The default and tab factories then clear the cached coordinates (`OutputFactory.pm` lines 306-323) and build the rows again over the record's restored span. The HGVSp is recomputed from that span's translation coordinates (`TranscriptVariationAllele.pm` lines 1724-1737), so a record whose span reaches into an intron has none, and `Amino_acids` prints the peptides the first pass left (lines 618-630). The VCF and JSON factories clear the coordinates before the rejoin, so their rows keep the first pass's HGVSp, its peptides and the `CDS_position` and `Protein_position` the shift cached, beside an unshifted `cDNA_position`. vep-rs builds a split allele's row once, over the allele's own span, in every format.

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

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--minimal --hgvs`, and `--tab`, `--vcf` or `--json` for the other formats. `--shift_3prime 1` gives the shifted rows and `--coding_only` the restored row; vep-rs takes the same flags over the same cache and FASTA.
