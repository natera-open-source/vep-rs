# HGVS transcript and protein notation resolved from the cache

A difference in which inputs resolve, on an input format no measured dataset uses; nothing is set aside.

Ensembl VEP refuses HGVS input offline (`Parser/HGVS.pm` line 104) and resolves a `c.` or `n.` notation through its core database. The transcript is found by stable id, or by external name for an `NM_` or `XM_` reference, else as the likely transcript of the gene the reference names (`VariationFeatureAdaptor.pm` lines 1829-1849 and 2146-2226). The position is resolved through the transcript's exon mapper (lines 1608-1751), and the alleles against the genome on the transcript's strand (lines 2054-2142).

vep-rs resolves the same notations from its transcript cache and the reference FASTA, with the same positions, alleles, strand, warnings and refusals, and the `tests/golden/116/GRCh38-input-hgvs` corpus compares its rows with Ensembl VEP 116.2's database-mode output. The cache has no index by identifier, so a run with transcript or protein notations in its input reads the cache's assembled chromosomes and scaffolds once before the first lookup. The cache's alternate loci, patches and LRG regions are read only for an id none of the assembled chromosomes and scaffolds carries.

Five readings differ:

- An `NM_` or `XM_` reference resolves against the cache's RefSeq transcripts, present under `--refseq` or `--merged`, and is otherwise refused as a reference no transcript carries (`Could not get a Transcript object`). Ensembl VEP reads such a reference from its RefSeq database whatever cache the run annotates with.
- The lookup reads the run's transcript set, so under `--gencode_basic` or `--gencode_primary` a transcript outside the set is not found, where Ensembl VEP's resolver reads its whole database.
- A gene name resolves among the transcripts on the assembled chromosomes and the unplaced scaffolds, in database-id order, where Ensembl VEP's `is_reference` test reads the database.
- A description whose position its regular expression leaves undefined (`c.-A>G`, a sign with no digits) is refused as one it could not parse, where Ensembl VEP carries the undefined position on.
- A notation every core group refuses is skipped with the first group's message alone, where Ensembl VEP joins every group's distinct exception text under its `Unable to parse HGVS notation` warning (`Parser/HGVS.pm` lines 199, 208 and 214). So in its database run a refused `ENST` notation also carries the `otherfeatures` group's `Could not get a Transcript object` message, and a refused gene name that group's own refusal.

## Protein notation

A `p.` notation resolves the way Ensembl VEP's protein branch resolves it (`VariationFeatureAdaptor.pm` lines 1900-1986). The transcript is found by the translation's stable id, by the transcript's stable id for an `ENST` reference written with `p.`, else as every transcript of the gene the reference names, tried in turn, the first to yield a variant taken. A substitution or `delins` is the one nucleotide path with the fewest changes from the genome's codon to a codon of the new residue (lines 2232-2319), the codon placed through `pep2genomic` (`TranscriptMapper.pm` lines 372-388). A deletion is the bases of its codons (lines 2323-2354). A frameshift is refused (line 1903).

Ensembl VEP lists every equally short path only under `--ambiguous_hgvs`, which vep-rs refuses. A residue change with several shortest paths (`Asp2Glu` on `GAC`, a base away from `GAA` and from `GAG`) is therefore skipped with the message Ensembl VEP writes for it without that flag, `Could not uniquely determine nucleotide change`. The corpus's twenty protein notations, on the forward-strand CRYAA and reverse-strand TFF3 transcripts, compare with Ensembl VEP 116.2's database-mode output.

Six readings of the protein branch differ:

- Ensembl VEP reads an `NP_`, `XP_` or `NM_` reference and a UniProt accession through its external-name tables (lines 1915-1916 and 1929), so in database mode `NP_000385.1:p.Thr4Ala`, `NM_000394.4:p.Thr4Ala` and `P02489:p.Thr4Ala` all resolve to CRYAA's `ENSP00000291554:p.Thr4Ala`. vep-rs reads no external name: a RefSeq protein or transcript id resolves where the cache carries the RefSeq model (under `--refseq` or `--merged`) and is otherwise refused, and a UniProt accession is read as a gene name and refused.
- The positions of a codon are ordered as numbers, where Ensembl VEP sorts them as strings (line 2304) and so misplaces a one- or two-base change in a codon whose positions straddle a power of ten.
- A residue position the exon mapper places in no exon is skipped as a sequence that cannot be found (`Unable to find the genomic reference sequence for protein`). Ensembl VEP's resolver dies on the mapper's gap there, and its parser writes the text of that `die` as the warning.
- A deletion range written end first (`Ile5_Thr4del`) is skipped with the same message, where Ensembl VEP builds a feature with no reference allele and writes neither a row nor a warning for it.
- When every transcript of a gene refuses a notation, the distinct messages are joined on one line in the order met, where Ensembl VEP joins them in hash order.
- A notation every core group refuses is skipped with the first group's message alone, as in the transcript branch, where Ensembl VEP joins every group's distinct exception text under its warning.

Two readings agree on a codon that crosses an exon boundary. A deletion whose codons cross one is refused in both engines, Ensembl VEP's resolver reading the bases of every exon piece against the span of the first (line 2351). A substitution whose shortest path changes the first and third base of such a codon writes, in both, a variant spanning the intron between them, since the path's positions are read as a genomic span (lines 2312-2314).
