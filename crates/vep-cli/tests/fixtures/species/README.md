# A cache of another species

`json_cache` is a hand-built cache of a species Ensembl ships no human-specific
data for, in the layout the converter writes: `info.json` names the species
`synthetic_species` on the assembly `SYN1`, no SIFT or PolyPhen version, and
variation columns holding no frequency column (the column set of an Ensembl
release 116 cache of a non-human species). Chromosome `1` holds two transcripts
with no `gencode_basic`, `gencode_primary`, `MANE_Select`, `MANE_Plus_Clinical`,
`appris` or `TSL` attribute and no CCDS, RefSeq, UniProt, phenotype or protein
domain key:

- `SYNT0000001` (gene `SYNG0000001`, `genA`): forward, protein coding, two exons
  at 1001-1060 and 1201-1260; a 9-base 5' UTR, a 96-base coding sequence of 32
  codons (`ATG` at 1010-1012, the stop `TGA` at 1243-1245) and a 15-base 3' UTR
  at 1246-1260, so 1014 is the middle base of the second codon (`GCT`, Ala).
- `SYNT0000002` (gene `SYNG0000002`, `genB`): forward, lncRNA, one exon at
  900-1080.

The one known variant, `rs1`, is a `C/T` at 1014 with no frequency, citation or
synonym. `CHECKSUMS.sha256` lists every file, so a run reads the cache with no
manifest warning.
