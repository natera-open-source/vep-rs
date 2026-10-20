# HGVS transcript notation resolved from the cache

A difference in which inputs resolve, on an input format no measured dataset uses; nothing is set aside.

Ensembl VEP refuses HGVS input offline (`Parser/HGVS.pm` line 104) and resolves a `c.` or `n.` notation through its core database. The transcript is found by stable id, or by external name for an `NM_` or `XM_` reference, else as the likely transcript of the gene the reference names (`VariationFeatureAdaptor.pm` lines 1829-1849 and 2146-2226). The position is resolved through the transcript's exon mapper (lines 1608-1751), and the alleles against the genome on the transcript's strand (lines 2054-2142).

vep-rs resolves the same notations from its transcript cache and the reference FASTA, with the same positions, alleles, strand, warnings and refusals, and the `tests/golden/116/GRCh38-input-hgvs` corpus compares its rows with Ensembl VEP 116.2's database-mode output. The cache has no index by identifier, so a run with transcript notations in its input reads the cache's assembled chromosomes and scaffolds once before the first lookup. The cache's alternate loci, patches and LRG regions are read only for an id none of the assembled chromosomes and scaffolds carries.

Four readings differ:

- An `NM_` or `XM_` reference resolves against the cache's RefSeq transcripts, present under `--refseq` or `--merged`, and is otherwise refused as a reference no transcript carries (`Could not get a Transcript object`). Ensembl VEP reads such a reference from its RefSeq database whatever cache the run annotates with.
- The lookup reads the run's transcript set, so under `--gencode_basic` or `--gencode_primary` a transcript outside the set is not found, where Ensembl VEP's resolver reads its whole database.
- A gene name resolves among the transcripts on the assembled chromosomes and the unplaced scaffolds, in database-id order, where Ensembl VEP's `is_reference` test reads the database.
- A description whose position its regular expression leaves undefined (`c.-A>G`, a sign with no digits) is refused as one it could not parse, where Ensembl VEP carries the undefined position on.

Beyond those four, a notation naming a protein (`p.`) is refused as unimplemented.
