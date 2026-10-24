# A VCF input without its `#CHROM` line

A difference in the header of the VCF output, not in any row, on an input shape no measured dataset has: a VCF whose header carries meta lines but not the mandatory `#CHROM` column line; nothing is set aside.

Ensembl VEP's VCF output copies the input's `#` lines into its header and takes the last of them as the column line (`OutputFactory/VCF.pm` line 224), which it prints after the `##VEP-command-line` line. A header without a `#CHROM` line hands its last meta line to that place. For an input whose header is `##fileformat=VCFv4.2` then `##source=probe`, the output's header is:

```
##fileformat=VCFv4.2
##VEP="v116.2" API="v116" time="..." cache="..." ...
##INFO=<ID=CSQ,...>
##VEP-command-line='...'
##source=probe
```

The output has no column line, and the meta line stands where the column line belongs.

vep-rs takes the input's `#CHROM` line as the column line when the header has one and writes the eight-column line `#CHROM POS ID REF ALT QUAL FILTER INFO` otherwise; every meta line stays in the meta block. For the same input it writes `##source=probe` after `##fileformat=VCFv4.2`, then the `##VEP`, `##vep-rs`, `##INFO` and `##VEP-command-line` lines, then the column line. Both engines annotate the record, and the `##INFO` definition and the data line agree byte for byte.

A VCF requires the column line, and a meta line in its place leaves the output outside the format, so the popped line is a defect of the copy and not a reading vep-rs reproduces. A single-`#` line after the `#CHROM` line, which the same reading takes for the column line and the sample list, is [its own page](vcf-hash-line-after-chrom.md).

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--vcf`, over a VCF of the two meta lines above and one record; vep-rs with the same flags over the same cache.
