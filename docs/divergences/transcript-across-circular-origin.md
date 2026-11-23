# A transcript stored across the origin of a circular sequence

**The defect in one sentence.** An Ensembl cache stores a transcript that lies across the origin of a circular sequence under two sets of coordinates that disagree, and Ensembl VEP annotates it on neither side of the origin; vep-rs annotates it on both. A record among its bases at the start of the sequence gets the transcript's row with `intergenic_variant`, and a record among its bases at the end gets no row of the transcript.

**Why this is a defect.** The bases are the gene's own. In the cache examined below, the stored coding sequence equals the reference read from the gene's first base to the last base of the sequence and on from base 1. A row that names a transcript, a gene and a strand cannot also place the record between genes. On the other side of the origin a change inside the coding sequence gets no row of the gene at all. The two halves of one gene are treated differently only because the sequence's numbering starts inside it.

**What vep-rs does instead.** vep-rs reads such a transcript as one unbroken span that runs on past the last base of the sequence, and reads a record among the first bases one sequence length further on. Each record then gets the row the gene would get on a sequence whose numbering did not cut it. That is a coding, start codon or stop codon term with its cDNA, CDS and protein positions, or the upstream or downstream term with its distance. Location, and every field that belongs to the record, keeps the record's own coordinates.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116, `ensembl` at release/116):

- `Bio/EnsEMBL/DBSQL/TranscriptAdaptor.pm`, lines 1938-1947: fetched on a circular slice, a transcript whose start is above its end has the sequence length taken from its start. The cache dump stores that start, zero or negative, beside the unchanged end.
- `Bio/EnsEMBL/TranscriptMapper.pm`, lines 159-161: the mapper is built from each exon's `seq_region_start` and `seq_region_end`, which on a circular slice are read back from the database (`Bio/EnsEMBL/Feature.pm`, lines 1080-1086). The stored mapper pair therefore runs from the gene's first base, near the end of the sequence, to its last base, near the start.
- `Bio/EnsEMBL/Slice.pm`, lines 635-638: `is_circular` answers 0 for a slice without an adaptor. Offline no slice has one, so no later step treats the sequence as circular, whatever the cache's `circular` flag says.
- `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, lines 103-109: a transcript takes the records between its stored start less the flank and its stored end plus the flank. A record near the end of the sequence lies outside that window, and no row is built for it.
- `Bio/EnsEMBL/Mapper.pm`, lines 311-316 and 1232-1256: a genomic span is looked up in an interval tree over the mapper pairs, and an interval whose end is below its start is built with the two swapped (lines 1241-1243). A record among the first bases lies below that interval, no pair is found, and lines 431-435 return a gap.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 636-668, `within_cds`: true only for a `Bio::EnsEMBL::Mapper::Coordinate`. The coding predicates read the same coordinates, so each of them fails on the gap.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`: line 480 sets `within_feature` from the stored span, lines 540-557 set `coding` on a one-gap result, and line 285 assigns `$DEFAULT_OVERLAP_CONSEQUENCE` when no predicate is true. `Utils/Constants.pm`, lines 636-650, defines that default as `intergenic_variant`.

**What correcting it would change.** Every record within reach of such a transcript would get the transcript's row with its term and positions, as a record on any other gene does. No row of another transcript changes.

**The class it produces.** `transcript_across_circular_origin`, on a transcript whose cache record has a start below 1. Among the first bases of the sequence it is a `Transcript` row whose Consequence is `intergenic_variant` in Ensembl VEP and a term of the transcript in vep-rs. Among the last bases, and within 5,000 bases before the gene on that side, it is a row of the transcript that vep-rs alone writes. Nothing in either engine's output marks such a transcript, so no comparison sets these rows aside. A record that reaches one counts as a difference, charged to vep-rs in the concordance figures.

**Where the shape is found.** The Ensembl Genomes release 63 cache of *Synechocystis* sp. PCC 6803 stores five such transcripts among its 3,614: one on the chromosome and one on each of its four plasmids. Four of the five carry their coding sequence; the fifth is stored as unknown bases throughout, and vep-rs writes `coding_sequence_variant` inside it, the term that sequence gives. The committed caches hold no transcript of the shape.

## 1. Four substitutions on a synthetic circular sequence

**In short.** The sequence `ring` is circular and 12,000 bases long. The forward gene `RINGG1` has 150 coding bases and no untranslated region: its first base is base 11,911 and its last is base 60. Its transcript `RINGT1` is stored as an Ensembl cache dump stores one of this shape, with start -89, end 60, a mapper pair that runs from 11,911 to 60, and a slice flagged circular. Ensembl VEP 116.3 annotates none of the four records below on the gene; vep-rs annotates all four.

**The records.** Each is a substitution of one base, with default flags and no FASTA:

| Record | Base of `ring` | Where it lies |
| --- | --- | --- |
| `low_side_31` | 31 | coding base 121, the first base of codon 41 |
| `last_gene_base_60` | 60 | coding base 150, the third base of the stop codon |
| `before_the_gene_11000` | 11,000 | 911 bases before the gene's first base |
| `high_side_11950` | 11,950 | coding base 40, the first base of codon 14 |

**Ensembl VEP 116.3.** The columns are Location, Feature, Feature_type, Consequence, cDNA_position, CDS_position, Protein_position and Amino_acids, the records in the table's order:

```text
ring:31 | RINGT1 | Transcript | intergenic_variant | - | - | - | -
ring:60 | RINGT1 | Transcript | intergenic_variant | ?--11850 | - | - | -
ring:11000 | - | - | intergenic_variant | - | - | - | -
ring:11950 | - | - | intergenic_variant | - | - | - | -
```

**vep-rs.** The same columns:

```text
ring:31 | RINGT1 | Transcript | missense_variant | 121 | 121 | 41 | T/A
ring:60 | RINGT1 | Transcript | stop_retained_variant | 150 | 150 | 50 | *
ring:11000 | RINGT1 | Transcript | upstream_gene_variant | - | - | - | -
ring:11950 | RINGT1 | Transcript | missense_variant | 40 | 40 | 14 | I/V
```

**Why vep-rs's rows are the correct ones.** Base 12,000 is followed by base 1, so the gene's bases 1 to 90 are bases 11,911 to 12,000 of the sequence and its bases 91 to 150 are bases 1 to 60. Base 31 is therefore coding base 90 + 31 = 121, the first base of codon 41, and base 11,950 is coding base 11,950 - 11,911 + 1 = 40, the first base of codon 14.

Each of those two changes alters its codon's residue under the bacterial genetic code the cache names, a missense change. Base 60 is the third base of the stop codon, and the change leaves a stop codon. Base 11,000 lies 11,911 - 11,000 = 911 bases before the first base of a forward gene, inside the 5,000-base flank: an upstream row with `DISTANCE=911`.

**The row at the gene's last base.** Base 60 is the one base of the gene that meets the swapped interval, at its lower end. The mapper then returns a cDNA span from 1 to 1 + (60 - 11,911) = -11,850 (`Bio/EnsEMBL/Mapper.pm`, lines 345-352 and 380-407), which Ensembl VEP prints as the cDNA position `?--11850` beside the default term.

**Reproducing it.** The entry's reproduction, `scripts/divergences/reproductions/transcript_across_circular_origin.toml`, holds the four records and both engines' rows. Its prepare steps build the synthetic cache inside the Ensembl VEP image with `make_cache.pl`, under `scripts/divergences/reproductions/inputs/transcript_across_circular_origin/`, which lists the gene's fifty codons. vep-rs reads the `json_cache` beside that script, and a prepare step holds it to the conversion of the built cache by `tools/vep-storable-to-json`. `tools/vep-divergence-reevaluate`, given `--entry transcript_across_circular_origin`, runs it on both engines and holds each row to the recorded one.
