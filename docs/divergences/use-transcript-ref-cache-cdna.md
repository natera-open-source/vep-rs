# A coding transcript's cDNA under `--use_transcript_ref`

A difference in one field's text, and in the terms of a row, under a flag no measured dataset is run with and a FASTA that disagrees with the cache; nothing is set aside.

Under `--use_transcript_ref`, Ensembl VEP reads every transcript's spliced cDNA from the FASTA, exon by exon (`Transcript::spliced_seq`, cached per transcript at `TranscriptVariation.pm` line 462), and cuts the reference allele from it. vep-rs assembles a coding transcript's cDNA from the cache's 5' UTR, coding sequence and 3' UTR, and reads the FASTA where that assembly gives none: for a non-coding transcript, or a coding transcript whose cached pieces do not fit together. On a FASTA of the genome the cache was built from, the two sequences are one.

They differ where the FASTA carries other bases at a coding transcript's exonic base: a chromosome masked to `N` outside a region of interest, or a FASTA of another build of the assembly. There Ensembl VEP's `USED_REF` is the FASTA's bases and the row's terms follow them; vep-rs's are the cache's. With a chromosome 21 FASTA masked to `N` at the record, the deletion `21 30541662 AG A` on ENST00000334046 writes `USED_REF=N` and no `Amino_acids` in Ensembl VEP, and `USED_REF=G` with `Amino_acids` `F/X` in vep-rs.
