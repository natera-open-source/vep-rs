# A FASTA record named `chrM` that is another sequence, read for `MT` through a synonym

**The defect in one sentence.** On GRCh37, where a run's reference FASTA holds its only mitochondrial record under the name `chrM`, Ensembl VEP serves every read for the cache's `MT` from that FASTA record, because the cache's synonym table links the two names. It tests only that the FASTA record exists, so where it is another sequence Ensembl VEP builds its notations and its reference checks from that sequence's bases; vep-rs reads no base there.

**Why this is a defect.** On GRCh37 `MT` names the assembly's mitochondrial sequence, NC_012920.1, of 16,569 bases, and the cache's transcripts stand at its positions. `chrM` is the UCSC hg19 spelling of NC_001807.4, of 16,571 bases. The two differ by four insertions or deletions, and a base of a coding position stands in NC_001807.4 one position after its place in the assembly ([the two sequences](plugin-contig-mitochondrial-spelling.md)). Both engines place an input record written `MT`, `chrMT` or `NC_012920.1` on the cache's `MT`, so its rows are on the assembly's transcripts and the reads for them are reads for `MT`. Ensembl VEP places one written `chrM` there too.

A base read from that FASTA record at a position of `MT` is therefore the base of a neighbouring position. Ensembl VEP writes it as the reference base of HGVSc and HGVSg and shifts an insertion or a deletion along it. For an input record written `MT` it compares REF with it under `--check_ref` and gives it to the input record under `--lookup_ref`. It writes no warning that the base is another sequence's.

The GRCh37 corpus `tests/golden/116/GRCh37` holds one mitochondrial input record, `MT 3308 T G`. Its REF is the assembly's base and the second base of the codon `aTa` of `ENST00000361390`. Over the file UCSC serves as hg19's `chrM`, which the tree does not hold, Ensembl VEP 116.3 refuses that input record under `--check_ref`. Under `--hgvs` it writes the row's Codons, which it cuts from the cache's coding sequence, as `aTa/aGa` and its HGVSc as `c.2A>G`: two reference bases for one position.

**What vep-rs does instead.** On GRCh37 vep-rs takes a mitochondrial FASTA record whose length is not 16,569 for another sequence than the assembly's. Such a FASTA record is read under its own pair of names alone, `chrM` and `M`, or `MT` and `chrMT`. Neither a spelling vep-rs builds in nor a synonym, of the cache's table or of `--synonyms`, hands it to a name of the other pair or to an accession.

vep-rs reads the reference of an input record it places on a cache chromosome under the cache's name for it. A read for `MT` reaches no FASTA record in such a run, so `MT` is [a contig the FASTA lacks](../cache-setup.md#a-contig-the-fasta-lacks): a substitution or a deletion gets no HGVSc, an insertion the HGVSc it has unshifted, and no input record an HGVSg. Under the two reference checks an input record's own reference is read under the cache's name and, where that reaches no FASTA record, under the input record's own name. An input record written `MT` or `NC_012920.1` is skipped.

An input record written `chrM` or `M` is skipped before any check: over such a file it is on [a sequence the cache does not hold](mitochondrial-record-of-another-sequence-placed-on-mt.md).

A FASTA record named `chrM` of 16,569 bases is read for `MT` by both engines. The length is how vep-rs tells the other sequence, and the one test it makes. A FASTA record of 16,569 bases that held other bases would be read, and one of another length that held the assembly's bases at their positions would not; neither was measured. The rule holds on GRCh37 alone: on GRCh38 `chrM` is the assembly's sequence, and both engines read it for `MT`.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/Utils/FastaSequence.pm`, lines 378-397: the reader asks the FASTA for a slice's own name and, where the FASTA holds none, takes the first of the name's `chr` spellings and synonyms whose length is above zero. With its debug lines on it prints `USING SYNOYM chrM FOR MT` (line 394).
- `AnnotationSource/Cache/Transcript.pm`, lines 249-251: a transcript keeps the slice the cache stores for it, named `MT`. `AnnotationType/Transcript.pm`, line 122: an input record that has no slice takes the transcript's.
- `TranscriptVariationAllele.pm`, lines 1474-1476: HGVSc takes its reference allele from the transcript's slice. `VariationFeature.pm`, line 1927: HGVSg takes it from the flank read on the input record's slice.
- `BaseVEP.pm`, line 413: under a reference check the input record's slice is looked up among the FASTA's record names by the input record's name. Lines 708-716 give a name the FASTA lacks the first of its synonyms the FASTA holds, and lines 425-433 give the slice that FASTA record's name and length.
- `Parser.pm`, lines 749-751: the input record's reference is read from that slice, for the comparison of `--check_ref` (line 614) and the replacement of `--lookup_ref` (lines 636 and 645).
- The GRCh37 cache's `chr_synonyms.txt` links `MT` and `NC_012920.1` to `chrM`, and no name to `M` or to `chrMT`.

**What correcting it would change.** A read for `MT` would be served from a FASTA record named `chrM` only where that FASTA record is the assembly's sequence.

**The class it produces.** `fasta_synonym_reads_another_mitochondrial_sequence`: on GRCh37, in a run whose FASTA holds its only mitochondrial record under `chrM` at another length than 16,569, the cells of this table for an input record the cache places on `MT`. The cells of [section 2](#2-a-fasta-that-holds-chrmt-beside-chrm), over a FASTA of two mitochondrial records, are of the class as well. No comparison sets these rows aside: the two outputs do not tell which FASTA a run read, so an input record with such a row counts as a difference in the concordance figures.

| Flag | Input record | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `--hgvs`, `--hgvsg` | a substitution written `MT`, `chrMT` or `NC_012920.1` | the notation's reference base is the FASTA record's base | no HGVSc, no HGVSg |
| `--hgvs`, `--hgvsg` | a deletion or an insertion written `MT`, in a non-coding transcript | shifted and tested for a duplication along the FASTA record's bases | no HGVSc for a deletion, the unshifted HGVSc for an insertion, no HGVSg |
| `--check_ref` | written `MT` or `NC_012920.1` | REF compared with the FASTA record's base: refused with the mismatch warning alone, or admitted | skipped with the sub-slice warning, then the mismatch warning |
| `--lookup_ref` | written `MT` or `NC_012920.1` | the FASTA record's base given as the reference | skipped with the could-not-fetch warning |

The table holds the cells this page measured, each flag alone and the first two together: on forward-strand transcripts, over a FASTA record of 16,571 bases, at positions that record holds. A substitution was measured in a coding transcript, and a deletion and an insertion under the name `MT` alone, in a non-coding one. In a coding transcript Ensembl VEP's protein notation follows the same shift (`TranscriptVariationAllele.pm`, lines 411-466); no such row was measured.

An input record that is refused or skipped gets no row. vep-rs's warnings are those of a contig the FASTA lacks. Under `--check_ref` they are `Could not fetch sub-slice from MT:3349-3349(1)`, since no base could be read, and then the mismatch warning, which follows every check that did not pass. Under `--lookup_ref` the one warning is `Could not fetch reference allele for MT:3349-3349(1)`. An input record Ensembl VEP admits, or looks up, is written as one `intergenic_variant` row whatever its name: its slice is named after the FASTA record and no transcript's slice has that name ([defect P](defect-p.md), the shape counted as `fasta_named_slice_intergenic`).

Outside the class:

- Under either reference check an input record written `chrMT` reaches no FASTA record in Ensembl VEP, which skips it as vep-rs does.
- An insertion names no reference base and passes either reference check in both engines, which write the same rows for it.
- An input record written `chrM` or `M` is a position of the FASTA record `chrM`, and vep-rs skips it in every run of the table, with and without a reference check. Ensembl VEP writes one written `chrM` on the transcripts of `MT`, with the notations of the table's first row, and compares it with the FASTA record under a reference check. Its rows are the entry [`mitochondrial_record_of_another_sequence_placed_on_mt`](mitochondrial-record-of-another-sequence-placed-on-mt.md).
- A FASTA record named `M` of another length is reached from `MT` by neither engine.

**How often the shape occurs.** Tier 3 under the rule of [Intended divergences from Ensembl VEP](../intended-divergences.md#differences-outside-the-classes). The shape is a whole file: a reference FASTA whose only mitochondrial record is named `chrM` and is not the assembly's sequence. No committed run passes one, so none of the 1,198 mitochondrial input records the corpora hold between them is annotated over such a file. The entry's row counts 97 flag sets: those that reach [a contig the FASTA lacks](../cache-setup.md#a-contig-the-fasta-lacks), an entry of the same flags.

## 1. A made-up reference two bases longer, named `chrM`

**In short.** `crates/vep-cli/tests/fixtures/mitochondrial-frame-check` holds a made-up reference of 16,569 bases named `MT`. The reproduction builds from it a file of one FASTA record, `chrM`, of 16,571 bases: the reference with its first two bases put in front once more, so its base at a position is the reference's base two positions back. The input record `MT 3349 other_base_two_back A C`, the third cell being its ID, states the reference's base as its REF, which is also the first base of the codon `Atc` of `ENST00000361390` in the cache. The longer file holds `G` at 3349.

Under `--hgvs --hgvsg` over the longer file, Ensembl VEP 116.3 writes the transcript's row with Codons `Atc/Ctc`, `HGVSc=ENST00000361390.2:c.43G>C` and `HGVSg=MT:m.3349G>C`. vep-rs writes the same row without the two notations and keeps `HGVSp=ENSP00000354687.2:p.Ile15Leu`. Over the reference itself both engines write `c.43A>C` and `MT:m.3349A>C`. The same input record written `chrMT` or `NC_012920.1` gets the same notations in Ensembl VEP, its HGVSg under the name as written, and none in vep-rs. Written `chrM` it gets them in Ensembl VEP and no row in vep-rs. Each flag alone gives its own notation the same way.

The reproduction's other input records fill the table: `MT 3350 T A`, where the longer file holds the reference's base, and a deletion and an insertion in the non-coding `ENST00000387347`. Under the notation flags every row of the input records written `MT`, `chrMT` or `NC_012920.1` is equal in the two engines once HGVSc, HGVSg and `HGVS_OFFSET` are set apart.

| Input record, over the longer file | Flag | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `MT 2953 TT T`, a deletion | `--hgvs --hgvsg` | `n.1286del`, `MT:m.2956del`, `HGVS_OFFSET=2`; over the reference `n.1284del`, `MT:m.2954del` | no HGVSc, no HGVSg |
| `MT 2959 G GA`, an insertion | `--hgvs --hgvsg` | `n.1289dup`, `MT:m.2959dup`; over the reference `n.1289_1290insA`, `MT:m.2959_2960insA` | `n.1289_1290insA`, no HGVSg |
| `MT 3349 A C`, where the file holds `G` | `--check_ref` | refused: the mismatch warning alone | skipped: the sub-slice warning, then the mismatch warning |
| `MT 3350 T A`, where the file holds `T` | `--check_ref` | admitted: one `intergenic_variant` row | skipped the same way: no base of `MT` can be read to check it |
| `MT 3349 A C` | `--lookup_ref` | one `intergenic_variant` row, `REF_ALLELE=G` | skipped: the could-not-fetch warning |

**Try it yourself.** Build the longer file as the reproduction's `prepare` steps do (`scripts/divergences/reproductions/fasta_synonym_reads_another_mitochondrial_sequence.toml`) and run both engines on its input records with `--fasta` on that file, under the caches of [Reproducing these examples](reproducing.md#reproducing-these-examples). The cache slice `crates/vep-cli/tests/fixtures/bayesdel-mitochondrial-contig/json_cache`, whose synonym table links `MT` to `chrM`, serves vep-rs.

## 2. A FASTA that holds `chrMT` beside `chrM`

**In short.** Everything above this section is of a FASTA whose only mitochondrial record is `chrM`; this section is of a FASTA that holds a second one. The UCSC hg19 file set that adds the assembly's sequence, `hg19.p13.plusMT`, lists `chrM` of 16,571 bases as the 91st of its 298 sequences and `chrMT` of 16,569 bases as the last. Over a file of that layout vep-rs reads `chrMT` for the cache's `MT` under every flag, whichever of the two the file lists first. Ensembl VEP 116.3 reads `chrMT` for the notations. Under the two reference checks it reads `chrM` for an input record written `MT` or `NC_012920.1`.

**Where it is in Ensembl VEP's source.** Under a reference check a run with no database takes the name of an input record's slice from the FASTA's record names (`BaseVEP.pm`, line 413). Lines 708-716 give a name the FASTA lacks the first of its synonyms the FASTA holds, and lines 719-731 try the name with `chr` before it only where no synonym was taken. `MT` so becomes `chrM`, the synonym of the cache's table, and `chrMT` is never tried. For a notation the reader tries the name with `chr` before it first (`FastaSequence.pm`, line 384) and reads `chrMT`. Neither route compares a length or reads a base to choose.

**What vep-rs does.** A FASTA record's own name, in any letter case, reads that record. A name that no FASTA record bears goes to one record by the steps below, tried in turn, and the order of the file's records does not enter:

1. On GRCh37, for `MT` and `chrMT`: a mitochondrial record of 16,569 bases, one named `MT` or `chrMT` before one named `chrM` or `M`.
2. The record whose name is the asked name with `chr` put before it or taken off: `chrMT` for `MT`, `chrM` for `M`.
3. For a mitochondrial name, a record of the other pair of names, on GRCh37 only one of 16,569 bases.

Where two records stand alike at a step, as `M` and `chrM` do for `MT`, the one whose name sorts first by byte is read, `M` there. A synonym of the cache's table or of `--synonyms` is tried after the three steps, under the rule this page gives above for a record of another length.

The input record of every row below is the substitution of section 1 at position 3349, with ALT `C`: the made-up reference holds `A` there and the record two bases longer `G`. The first three rows are over `chrM` two bases longer beside `chrMT`, the made-up reference. The last is over another layout, the made-up reference named `chrM` beside a `chrMT` two bases longer, where step 1 has vep-rs read `chrM`.

| Flag | Input record | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `--check_ref` | written `MT` or `NC_012920.1`, REF the base of `chrMT` | refused: REF compared with the base of `chrM` | annotated on the transcripts of `MT` |
| `--check_ref` | written `MT` or `NC_012920.1`, REF the base of `chrM` | admitted: one `intergenic_variant` row | refused with the mismatch warning alone |
| `--lookup_ref` | written `MT` or `NC_012920.1` | `REF_ALLELE=G`, the base of `chrM`, on one `intergenic_variant` row | `REF_ALLELE=A`, the base of `chrMT`, on the transcripts of `MT` |
| `--hgvs`, `--hgvsg` | written `MT`, over the `chrMT` two bases longer | `c.43G>C` and `MT:m.3349G>C`, the base of `chrMT` | `c.43A>C` and `MT:m.3349A>C`, the base of `chrM` |

Each cell was measured with its flag alone, the two notation flags also together. Over the first file both engines write the same rows under `--hgvs` and `--hgvsg`, for an input record written `MT`, `chrMT` or `NC_012920.1`, and under either check both compare an input record written `chrMT` with `chrMT`. Ensembl VEP writes each input record it admits as one `intergenic_variant` row, which is the difference of [defect P](defect-p.md) and not of this class. The last row is charged to Ensembl VEP on the length alone: its reader takes the name with `chr` before it as it takes a synonym, testing only that the FASTA record exists. This section says nothing of an input record written `chrM` or `M`.

Of the real pair, `chrMT` is the sequence the cache's transcripts stand on: the cached coding sequence of `ENST00000361390` equals its bases from the transcript's start. A FASTA of those two sequences, UCSC's `chrM` and the assembly's sequence named `chrMT`, was run in both orders; the tree does not hold it. Ensembl VEP 116.3 refuses the corpus's input record `MT 3308 T G` under `--check_ref` there, and vep-rs annotates it.

The last row's file is made up: no file is known that names a sequence of another length `chrMT`. Two made-up records of 16,569 bases each, with different bases at the position, were measured as well under the four flags: vep-rs reads `chrMT` for `MT`, and Ensembl VEP reads `chrM` under the two reference checks. No length tells two such records apart, and this page charges neither engine there. No committed run passes a FASTA of two mitochondrial records.

**Try it yourself.** The reproduction's `prepare` steps build both files. Its first further run holds the two `--check_ref` rows for the input record written `MT`, its second the `--lookup_ref` row under both names and its third the last row, each with `--fasta` on the file the row names.
