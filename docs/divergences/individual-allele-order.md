# An individual carrying two alternate alleles

A difference in the order of rows and alleles under one flag no measured dataset is run with; nothing is set aside.

Under `--individual`, Ensembl VEP builds each individual's copy of a record from the alleles its genotype names, collected in a hash (`Parser/VCF.pm` lines 635-639). An individual whose genotype names two alternate alleles (`1/2`) gets a copy whose allele string lists them in Perl's hash order, which changes from one run to the next. The copy's rows within each transcript, its JSON `allele_string` and its `ALLELE_NUM` values follow that order.

vep-rs lists the copy's alleles in the record's ALT order, an empty allele (a genotype bit past the record's alleles) last; Ensembl VEP's hash may put the empty allele first, where the copy's second allele reads as `-`. Under `--safe --individual --minimal` a copy of two alleles that trim is split, and the first allele of that order also decides its `REF_ALLELE`, `VARIANT_CLASS` and `UPLOADED_ALLELE` ([an individual's copy of a `--minimal` record under `--safe`](safe-individual-minimal.md)). The `GRCh38-individual` golden corpus carries two such records. Its reference outputs were written under a fixed Perl hash seed, and the two records list their ALTs in the order that seed gives, so a run of the corpus compares the row order too.
