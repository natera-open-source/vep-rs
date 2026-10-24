# Two region notations: `UPLOADED_ALLELE` under `--minimal`, and a `BND` word

Two differences on region input (`--format region`), a format no measured dataset uses; nothing is set aside. The first is in the text of one field and reaches every input format but VCF, the second in which records are written.

## `UPLOADED_ALLELE` of a region record `--minimal` trims

A region record's reference is a placeholder of one `N` per base of the span (`Parser/Region.pm` line 189), so `21:14221105-14221107/NA` is the record `NNN/NA`. Under `--minimal`, Ensembl VEP trims it to `NN/A` over 14221106 to 14221107 and keeps `NNN/NA` as the uploaded string (`Parser.pm` line 852, then line 1008). The alleles of the trimmed record still differ in length, so the test that follows (lines 862-881) minimises the record a second time. A region record carries no untrimmed string of its own, so that pass reads the trimmed `NN/A` as the uploaded string (line 1008 again). `UPLOADED_ALLELE` then prints `NN/A` on every row, a string the input never named (`OutputFactory.pm` line 1360).

vep-rs keeps the string as parsed, so under `--minimal --uploaded_allele` it prints `UPLOADED_ALLELE=NNN/NA`, and every other column of the 17 rows agrees with Ensembl VEP's, `MINIMISED=1` included. Without `--minimal` both engines print `NNN/NA`. The placeholder is the record as uploaded, and `--minimal`'s trim is not an upload, so the string the input named is the one vep-rs prints.

An Ensembl-format, HGVS, SPDI or identifier record behaves the same way: its parser hands the record to the same pass (`Parser/HGVS.pm` line 248, `Parser/SPDI.pm` line 204, `Parser/ID.pm` line 196, `Parser/VEP_input.pm` line 167) and sets no untrimmed string. For `21 10542414 10542417 ATGC/AG +` under `--minimal --uploaded_allele`, Ensembl VEP prints `TGC/G` and vep-rs `ATGC/AG`; for `21:g.30282168_30282170TCA>TG`, the same feature built from a notation, `CA/G` and `TCA/TG`, both `MINIMISED=1`, and without `--minimal` both engines print `TCA/TG`. VCF input is immune in both engines: the VCF parser sets the untrimmed string first (`Parser/VCF.pm` line 343), which the second pass keeps.

## A `BND` word

A region whose allele is a structural word (`DEL`, `DUP`, `INS:ME:ALU`, `BND`) is a structural variant over the span (`Parser/Region.pm` lines 208-216), built with no allele string. For `BND` Ensembl VEP then reads the record's breakends while collecting the cache regions, finds none to read, and dies as the buffer holding the record is annotated (`AnnotationSource.pm` line 251):

```
Can't use an undefined value as an ARRAY reference at
/opt/vep/src/ensembl-vep/modules/Bio/EnsEMBL/VEP/AnnotationSource.pm line 251.
```

The exit code is 255, and nothing is written for the record, nor for any record of the file the run had not written when it died. vep-rs annotates `21:14221000-14222000/BND` as a `chromosome_breakpoint` over the span, one row per overlapping transcript, 17 rows on the GRCh38 cache, and the run goes on. A `DEL`, `DUP` or `INS:ME:ALU` word, which names no breakend, is annotated alike by both engines.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), with `--format region --fasta <file>` and, for the first difference, `--minimal --uploaded_allele`; vep-rs with the same flags over the same cache and FASTA.
