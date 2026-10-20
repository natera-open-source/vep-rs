# An empty REF column at POS 0

A difference in whether one malformed record is written, on an input shape no measured dataset has; nothing is set aside.

VCF places a telomere at POS 0 and requires a REF of at least one base. Ensembl VEP reads a record whose REF column is empty at POS 0 (two tabs where the base would be, the ALT `A`) as the allele string `/A` over the coordinates 0 to -1 (`Parser/VCF.pm` line 254, `$end += length($ref) - 1`). `validate_vf` admits those coordinates, and the record is written as `intergenic_variant` at `21:-1-0`.

vep-rs holds no coordinate before 0. It reads the record as the insertion `-/A` at 0-0 and skips it with Ensembl VEP's own warning for an insertion whose coordinates are not start = end + 1. The warning reads `Alleles look like an insertion (-/A) but coordinates are not start = end + 1 (START=0, END=0)`. Under `--dont_skip` the record is kept as it stands.

A record at POS 0 whose REF has a base is read alike by both engines: `21 0 x A G` is annotated at `21:0` and `21 0 x AT A` at `21:1`. `21 0 x A AT` is skipped by both with the insertion warning, since `trim_sequences` reads the chopped record's end of 0 as no end and resets it to the start (ensembl-variation `Utils/Sequence.pm` line 971).
