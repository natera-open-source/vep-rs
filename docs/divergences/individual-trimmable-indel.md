# An individual's copy of an indel the anchor chop leaves trimmable

A difference in the alleles annotated for an individual, on an input shape no measured dataset has (its indels carry one anchor base); nothing is set aside.

Ensembl VEP minimises a bi-allelic indel after the anchor chop (`Parser.pm` lines 860-881: `ACAC/AC` becomes `AC/-` two bases on). It then builds each individual's allele string from the chopped alleles, `CAC` and `C`, against the minimised reference `AC` (`Parser/VCF.pm` lines 603-639), which neither equals. For that record a reference-homozygous individual (`0/0`) is annotated as `AC/CAC` at the minimised span with `ZYG` `HOM`, a heterozygous carrier (`0/1`) as `AC/CAC/C` with both in Perl's hash order, and a homozygous carrier (`1/1`) as `AC/C`.

vep-rs compares the chopped alleles with the chopped reference, so the `0/0` individual is reference-homozygous and a carrier's copy is annotated with the record's minimised alleles, `AC/-`. A record whose chop leaves nothing further to trim, an indel whose alleles share the anchor base and nothing else, reads alike in both engines.
