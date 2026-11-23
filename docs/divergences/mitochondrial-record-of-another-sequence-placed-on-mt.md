# An input record named `chrM` on another sequence, placed on `MT`

**The defect in one sentence.** On GRCh37, Ensembl VEP places an input record written `chrM` on the cache's `MT` through the cache's synonym table, whatever sequence the run's reference FASTA holds under that name. Where the FASTA record `chrM` is not the assembly's sequence, it writes the rows of the transcripts of `MT` at the input record's number, or one `intergenic_variant` row where a reference check reads the input record's base; vep-rs skips the input record with a warning. An input record written `M` is the same case in vep-rs; Ensembl VEP writes it as one `intergenic_variant` row wherever it admits it.

**Why this is a defect.** On GRCh37 `MT` names the assembly's mitochondrial sequence, NC_012920.1, of 16,569 bases, and the cache's transcripts stand at its positions. `chrM` is the UCSC hg19 spelling of NC_001807.4, of 16,571 bases, and `M` is that name without its prefix. The two sequences differ by four insertions or deletions, and a base of a coding position stands in NC_001807.4 one position after its place in the assembly ([the two sequences](plugin-contig-mitochondrial-spelling.md)).

A reference FASTA that holds a record of another length under `chrM` shows which sequence the name stands for in that run. An input record written `chrM` then states a position of that sequence. The two sequences run in step for about their first 300 positions and one or two positions apart after that. So at the input record's number the transcripts of `MT` hold the bases of a neighbouring position.

Ensembl VEP reads such an input record two ways. Without a reference check it resolves `chrM` to `MT` and writes the rows of the transcripts of `MT`. Under `--check_ref` or `--lookup_ref` it reads the reference base of an input record that states one, and for that names the input record's slice after the FASTA record, with that record's length. It then finds no transcript on a slice of that name and writes `intergenic_variant`.

Neither row is right for an input record that is a position of that FASTA record, where the two sequences are out of step. The first states a transcript position the input record does not change. The second is the row both engines write, under `--dont_skip` alone, for an input record whose chromosome name no annotation source or synonym holds. Without that flag both skip such an input record with a warning (`Parser.pm`, lines 189 and 542-546), as they do the input record on `Z` of the corpus `tests/golden/116/GRCh38-vcf-check-ref`. Ensembl VEP writes that row here without the flag and with no warning, for a sequence on which the cache holds no feature.

The GRCh37 corpus `tests/golden/116/GRCh37` holds one mitochondrial input record, `MT 3308 T G`. It changes the second base of the first codon of `ENST00000361390`, and both engines write `stop_gained,start_lost` for it. The file UCSC serves as hg19's `chrM` is not part of this repository. It holds that gene's coding sequence one position on, so the same substitution is there `chrM 3309 T G`, whose REF is the file's base. Over that file Ensembl VEP 116.3 writes it as `start_retained_variant` with Codons `atA/atG`, the codon's third base. Under `--check_ref` and under `--lookup_ref` it writes one `intergenic_variant` row.

**What vep-rs does instead.** On GRCh37 vep-rs reads `chrM` and `M` as the names of one sequence, and a FASTA record under either name whose length is not 16,569 as another sequence than the assembly's ([a FASTA record named `chrM` that is another sequence](fasta-synonym-reads-another-mitochondrial-sequence.md)). An input record written `chrM` or `M` whose name reads such a FASTA record is taken as on a chromosome the cache does not hold. vep-rs skips it.

The warning names the FASTA record and both lengths: `Chromosome chrM is the FASTA record chrM of 16571 bases, another sequence than the cache's MT of 16569 bases, and is read as no chromosome of the cache`. For an input record written `M` it opens `Chromosome M`, and it names the FASTA record as the file names it. It is one line of the warnings file for each such input record.

The chromosome test comes before the reference checks, so the warning is the same with and without `--check_ref` or `--lookup_ref`, and no notation flag changes it. Under `--dont_skip` the input record is kept as every input record on such a chromosome is, as one `intergenic_variant` row beside the warning.

Another annotation source can still hold the input record. Beside a `--custom` file that lists the chromosome, under the cache's name or the input record's, it is kept as any input record on a chromosome such a file lists is ([chromosomes](../custom-annotation.md#chromosomes)). It is then written on no transcript, with no warning, and takes the file's records under either name. A `--regulatory_gff` file was not measured.

The placement on `MT` stands wherever nothing tells the two sequences apart: in a run with no reference FASTA, over a FASTA that holds no record under `chrM` or `M`, and over a FASTA record under either name of 16,569 bases. The rule holds on GRCh37 alone: on GRCh38 `chrM` is the assembly's sequence. Where the cache itself holds a chromosome named `chrM` or `M`, the input record is annotated on it.

The length is the one test vep-rs makes, and its own bound. A FASTA record of another length that held the assembly's bases at their positions would be skipped as well, where Ensembl VEP's rows on `MT` would be right. A FASTA record of 16,569 bases that held other bases would not be skipped. Neither was measured, nor was an input record in the stretch where hg19's `chrM` runs in step with the assembly, which vep-rs skips too. Which FASTA record a name reads where the file holds two mitochondrial records is the FASTA reader's rule, and no such file was measured here.

Two names are outside the rule. Under a `--synonyms` table that links no name to `MT`, both engines rename a bare `M` to `MT` (`Parser.pm`, lines 794-796), and it is an input record on `MT` from then on. A breakend's mate is named through the table as in any other run. The rule was measured on VCF input.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Parser.pm`, lines 770-802 (`_have_chr`): an input record is admitted where its name is a chromosome of an annotation source or a synonym of one (lines 782-785), or reaches one with `chr` added or removed (lines 788-793), as `M` reaches `chrM`. The test reads names alone.
- `AnnotationSource.pm`, line 179, and `BaseVEP.pm`, lines 708-716: the cache's regions are read under the first synonym of the input record's name that the cache holds, `MT` for `chrM`. No line reads the FASTA.
- `AnnotationType/Transcript.pm`, line 122: an input record that has no slice takes the transcript's, so without a reference check its rows are written on the transcripts of `MT`.
- `Parser.pm`, line 749, and `BaseVEP.pm`, lines 413 and 425-433: under a reference check the slice of an input record that states a reference base is named after the FASTA record and takes its length, 16,571.
- `AnnotationType/Transcript.pm`, lines 141 and 143: a transcript whose slice has another name than the input record's is passed over. The input record is left with `intergenic_variant` ([defect P](defect-p.md)).
- `BaseVEP.pm`, lines 719-730: adding or removing the prefix `chr` gives an input record written `M` no name the cache holds, so it overlaps nothing with or without a reference check ([a bare `M`](bare-m-chromosome.md)).
- The GRCh37 cache's `chr_synonyms.txt` links `MT` to `chrM`, and no name to `M`.

**What correcting it would change.** An input record written `chrM` would be placed on `MT` only where the reference FASTA does not show `chrM` to be another sequence. Where it does, the input record would be skipped as on a chromosome no annotation source holds.

**The class it produces.** `mitochondrial_record_of_another_sequence_placed_on_mt`: on GRCh37, in a run whose reference FASTA holds under `chrM` or `M` a record of another length than 16,569, the rows of an input record written `chrM` or `M`. No comparison sets these rows aside: the two outputs do not tell which FASTA a run read, so such an input record counts as a difference in the concordance figures.

| Flag | Input record | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| none | written `chrM` | the rows of the transcripts of `MT` at its number | skipped with the warning |
| none | written `M` | one `intergenic_variant` row | skipped with the warning |
| `--lookup_ref` | written `chrM` or `M` | one `intergenic_variant` row, its reference the FASTA record's base | skipped with the warning |
| `--check_ref` | written `chrM` or `M`, REF the FASTA record's base | one `intergenic_variant` row | skipped with the warning |
| `--check_ref` | written `chrM` or `M`, REF another base | refused with the mismatch warning | skipped with the warning |

The table holds the cells this page measured, on one-base substitutions over a FASTA record of 16,571 bases named `chrM`: the reproduction's four input records and the real example above. Under `--hgvs` and `--hgvsg` Ensembl VEP's rows for the input record written `chrM` carry notations whose reference base is the FASTA record's ([a FASTA record named `chrM` that is another sequence](fasta-synonym-reads-another-mitochondrial-sequence.md)). No deletion and no other input format was measured. An insertion states no reference base: by the code Ensembl VEP reads no slice for it under either check (`Parser.pm`, lines 600-603 and 636) and so writes it as it does without one. None was measured.

**How often the shape occurs.** Tier 3 under the rule of [Intended divergences from Ensembl VEP](../intended-divergences.md#differences-outside-the-classes). The shape is a whole file beside a name: a reference FASTA whose record named `chrM` or `M` is not the assembly's sequence, and an input record written under one of those names. No committed run passes such a file, so none of the 1,198 mitochondrial input records the corpora hold between them is annotated over one. The entry's row of the index counts all 511 committed flag sets, since the placement needs no flag.

## 1. Four input records over a made-up reference two bases longer, named `chrM`

**In short.** `crates/vep-cli/tests/fixtures/mitochondrial-frame-check` holds a made-up reference of 16,569 bases named `MT`. The reproduction builds from it a file of one FASTA record, `chrM`, of 16,571 bases: the made-up reference with its first two bases put in front once more, so its base at a position is the made-up reference's base two positions back. It shifts every base by two positions, where hg19's `chrM` shifts the gene above by one. On the cache's `MT`, position 3349 is the first base of the codon `Atc` of `ENST00000361390`.

The input record `chrM 3349 named_chrm_other_base A C`, the third cell being its ID, states the made-up reference's base as its REF, where the longer file holds `G`. With no further flag Ensembl VEP 116.3 writes 20 rows for it on the transcripts of `MT`, among them `missense_variant` with Codons `Atc/Ctc` on `ENST00000361390`. vep-rs writes no row and the warning. The same input record written `M` gets one `intergenic_variant` row in Ensembl VEP.

The reproduction's other two input records, at 3350, state `T`, which the made-up reference and the longer file both hold there. They are the input records Ensembl VEP admits under `--check_ref`. The made-up file shows the mechanism alone: which row is right for an input record rests on the real example above.

| Input record | Flag | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `chrM 3349 A C` | none | the rows of the transcripts of `MT`, Codons `Atc/Ctc` | skipped |
| `M 3349 A C` | none | one `intergenic_variant` row | skipped |
| `chrM 3349 A C`, `M 3349 A C` | `--lookup_ref` | one `intergenic_variant` row each, `REF_ALLELE=G` | skipped |
| `chrM 3350 T A`, `M 3350 T A` | `--check_ref` | one `intergenic_variant` row each | skipped |
| `chrM 3349 A C`, `M 3349 A C` | `--check_ref` | refused: the mismatch warning | skipped |

Over the made-up reference itself, named `MT`, `chrM` or `M`, vep-rs writes all four input records on the transcripts of `MT`, as it does in a run with no FASTA. For the two written `M` that placement is the entry [a bare `M`](bare-m-chromosome.md).

**Try it yourself.** Build the longer file as the reproduction's `prepare` steps do (`scripts/divergences/reproductions/mitochondrial_record_of_another_sequence_placed_on_mt.toml`) and run both engines on its input records with `--fasta` on that file, under the caches of [Reproducing these examples](reproducing.md#reproducing-these-examples). The cache slice `crates/vep-cli/tests/fixtures/bayesdel-mitochondrial-contig/json_cache`, whose synonym table links `MT` to `chrM`, serves vep-rs.
