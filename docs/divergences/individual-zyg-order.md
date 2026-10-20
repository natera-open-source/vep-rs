# The ZYG list order under `--individual_zyg`

A difference in the order of one list under one flag no measured dataset is run with; nothing is set aside.

Under `--individual_zyg`, Ensembl VEP keeps each named sample's genotype in a hash keyed by sample name (`Parser/VCF.pm` line 737) and builds the record's `ZYG` list from that hash's keys (`OutputFactory.pm` line 983). Perl's hash order changes from one process to the next, so the same input lists the samples in another order on another run, in the default, tab and VCF formats and in the JSON `zyg` array alike. The entries themselves, each sample's `HET`, `HOM`, `HETREF` or `HOMREF`, are the same on every run.

vep-rs lists the samples in the order the run names them, the header's column order under `all` and the list's order under a list, so two runs print the same list. The `GRCh38-individual-zyg` golden corpus compares the `ZYG` list as a multiset of its entries, and its reference outputs were written under a fixed Perl hash seed, so one Ensembl VEP order is on record.
