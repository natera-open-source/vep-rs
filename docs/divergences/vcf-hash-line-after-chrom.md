# A `#` line after the `#CHROM` line

A difference in which samples a VCF header names under `--individual` and `--individual_zyg`, and in the header of the VCF output, on an input shape no measured dataset has: a single-`#` line after the mandatory `#CHROM` column line; nothing is set aside.

Ensembl VEP reads a VCF header through ensembl-io, which takes every single-`#` line as the column line in turn, so the last one stands (`BaseVCF4.pm` lines 96-97: `elsif ($line =~ /^#\s*(.+)$/) { $self->{'metadata'}->{'header'} = [split("\t",$1)];`). The sample names are that line's columns after its `FORMAT` column (lines 805-808), so a header whose `#CHROM` line is followed by `#note` names none.

Under `--individual all` the run writes no row and exits 0; under a list it stops with Ensembl VEP's unknown-sample error, which names the list against an empty header (`Parser/VCF.pm` lines 612-617); under `--individual_zyg all` every row is written without `ZYG`. The VCF output copies the `#CHROM` line among the meta lines and prints `#note` where the column line belongs (`OutputFactory/VCF.pm` line 224).

vep-rs reads the sample names from the `#CHROM` line alone and writes a single-`#` line after it nowhere. Take the eight-sample header of the `GRCh38-individual` corpus, then `#note`, then its `hb_multi` record (`21 34716340 G A,T`, genotypes `0/1 1/2 0/0 2/2 0/2 1/1 ./2 ./.`). Under `--individual all` vep-rs writes the 28 rows the header without the stray line gets, under `--individual S1` S1's four. Under `--individual_zyg all` it writes the record's eight rows, each carrying the seven called samples' `ZYG` entries, and the VCF output's column line is the `#CHROM` line.

A VCF 4.2 header is its `##` meta lines and one column line, the `#CHROM` line, so a single-`#` line after it is outside the format, and a run that names no sample without a word is a reading vep-rs does not reproduce. A header with no `#CHROM` line at all is [its own page](vcf-without-chrom-line.md).

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--individual all`, over the three-line input above; vep-rs with the same flags over the same cache.
