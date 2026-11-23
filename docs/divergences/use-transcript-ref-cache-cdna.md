# A coding transcript's cDNA under `--use_transcript_ref`

A difference in `USED_REF`, and in the codons, residues and consequence terms that follow it, under `--use_transcript_ref`. It needs a FASTA that disagrees with the cache or, on the cache's own genome, one of ten RefSeq models whose alignment correction failed. Nothing is set aside in the concordance figures; a format comparison sets aside seven cells of the rows of [section 1](#1-a-refseq-model-whose-alignment-correction-failed). Where the FASTA carries other bases at a coding transcript's exonic base under `--use_transcript_ref`, Ensembl VEP's `USED_REF` and the row's terms follow the FASTA; vep-rs's follow the cache.

Under `--use_transcript_ref`, Ensembl VEP reads every transcript's spliced cDNA from the FASTA, exon by exon (`Transcript::spliced_seq`, cached per transcript at `TranscriptVariation.pm` line 462), and cuts the reference allele from it. vep-rs assembles a coding transcript's cDNA from the cache's 5' UTR, coding sequence and 3' UTR, and reads the FASTA where that assembly gives none: for a non-coding transcript, or a coding transcript whose cached pieces do not fit together. On a FASTA of the genome the cache was built from, the two sequences are one.

They differ where the FASTA carries other bases at a coding transcript's exonic base: a chromosome masked to `N` outside a region of interest, or a FASTA of another build of the assembly. There Ensembl VEP's `USED_REF` is the FASTA's bases and the row's terms follow them; vep-rs's are the cache's. With a chromosome 21 FASTA masked to `N` at the record, the deletion `21 30541662 AG A` on ENST00000334046 writes `USED_REF=N` and no `Amino_acids` in Ensembl VEP, and `USED_REF=G` with `Amino_acids` `F/X` in vep-rs.

One exception stands on the cache's own genome: a RefSeq model with `BAM_EDIT=FAILED` that carries an edit of its own. Past that edit Ensembl VEP reads the reference allele one position away from the record's, and vep-rs reads it at the record's position. Where the model holds a base the genome lacks, Ensembl VEP reads one position further along the transcript. Where the model lacks a base the genome holds, it reads one position back. Section 1 gives the records, the rows and the source lines.

## 1. A RefSeq model whose alignment correction failed

**In short.** `NM_015068.3` (gene `PEG10`, `BAM_EDIT=FAILED`) holds one cDNA position the genome lacks. For the ClinVar substitution `7 94665299 G T` on GRCh38, past that position, Ensembl VEP 116.3 prints `USED_REF=T`, the base one position on, and `synonymous_variant`. vep-rs prints `USED_REF=G`, the record's own reference base, and `missense_variant`.

### 1.1 What a reader needs to know first

A RefSeq or merged cache whose models were corrected against an alignment switches `--use_transcript_ref` on (`CacheDir.pm` lines 411-416). A model whose corrected sequence equals the RefSeq sequence prints `BAM_EDIT=OK`. A model whose sequence does not equal it prints `BAM_EDIT=FAILED`: the alignment's edits are removed, and every base of the model is a base of the genome's exons.

A model can carry an edit of its own, an `_rna_edit` attribute the alignment did not write, and a failed model keeps it. Ten failed models of the release 116 RefSeq caches, the same ten on GRCh38 and GRCh37, carry such an edit, and on each it changes the model's length. Three models of `PEG10` write one cDNA position as two bases, so each holds one position the genome lacks. Seven models of `OAZ1`, `OAZ2` and `OAZ3` write one cDNA position as no base, so each lacks one position the genome holds.

The cached map between cDNA and genomic positions counts the edit. Join a model's exons and count along them: past the edit, a base's cDNA number is one more than that count on the three models and one less on the seven.

### 1.2 The records and both engines' rows

Two ClinVar records on GRCh38, run under `--refseq`, each on a forward-strand model. The first, `7 94665299 G T`, lies past the edit of `NM_015068.3`, which adds a position:

| Cell | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `synonymous_variant` | `missense_variant` |
| `CDS_position` | 1744 | 1744 |
| `Amino_acids` | `F` | `V/F` |
| `Codons` | `Ttt/Ttt` | `Gtt/Ttt` |
| `IMPACT` | `LOW` | `MODERATE` |
| `GIVEN_REF` | `G` | `G` |
| `USED_REF` | `T` | `G` |

The second, `1 151770257 G A`, lies past the edit of `NM_016178.2` (gene `OAZ3`), which removes a position:

| Cell | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `splice_region_variant,synonymous_variant` | `missense_variant,splice_region_variant` |
| `CDS_position` | 565 | 565 |
| `Amino_acids` | `S` | `G/S` |
| `Codons` | `Agt/Agt` | `Ggt/Agt` |
| `IMPACT` | `LOW` | `MODERATE` |
| `GIVEN_REF` | `G` | `G` |
| `USED_REF` | `A` | `G` |

The genome holds `G` at both records' positions. It holds `T` at `7:94665300`, one position on, and `A` at `1:151770256`, one position back. In both rows above Ensembl VEP's `USED_REF` equals the alternate allele, so its two codons are one and the row reads as no change.

### 1.3 Where it is in Ensembl VEP's source

`ensembl-vep` at `release/116.3`, `ensembl-variation` at commit `db7394a02523a6e519cb0d89f84b2889bac12737`, `ensembl` at commit `0d852313c420a1d5128b6774635125ebbac03350`:

- `Pipeline/DumpVEP/Dumper/Core.pm`, lines 111-113 and 164: the cache writer corrects each model, caches its map and coding sequence, and then deletes the key `edits_enabled`, which tells a transcript to apply its edits. The `_rna_edit` attributes stay (line 198).
- `AnnotationType/Transcript.pm`, lines 617 and 621-625: a corrected model is cached with its edited cDNA; a failed one is cached with no cDNA.
- `TranscriptVariation.pm`, lines 457 and 462-463: the reference allele is cut from the cached sequence or, without one, from `spliced_seq`, at the row's cDNA start.
- `BaseTranscriptVariation.pm`, lines 478-492 and 1074-1082: the cDNA start comes from the cached map, which counts the model's own edit (`TranscriptMapper.pm`, lines 170-223).
- `Transcript.pm` of the `ensembl` repository, lines 823-879: `spliced_seq` joins the exons and applies the edits only when `edits_enabled` is set (line 866). On a cached model it is not, so the join holds no edit.
- `TranscriptVariationAllele.pm`, lines 872-877: the reference allele is written over the coding sequence before the reference codon is read.
- `TranscriptVariationAllele.pm`, lines 1510-1513 and 1793: under `--hgvs`, the `HGVSc` of a model that carries an `_rna_edit` names the reference allele so cut as the base changed, and the `HGVSp` takes its reference residue from the reference codon.
- `OutputFactory/JSON.pm`, line 388: a JSON record's `most_severe_consequence` is the most severe term among those of all its rows, the one of lowest rank.

A number that counts the edit therefore indexes a sequence without it.

### 1.4 Why vep-rs's row is the correct one

A failed model's sequence is the genome's exons, and its own edit changes no base after it. The transcript base a substitution changes is the genome's base at the record's position, which is the record's `REF`. The codon at the row's `CDS_position`, read from the genome through the cached map, is vep-rs's reference codon on each of the 18 coding rows that section 1.5 counts, and its residue is vep-rs's.

vep-rs reads the base at the row's cDNA number from a sequence numbered as the map numbers it: the cache's cDNA, which holds the model's own edit, or the genome through the map. Ensembl VEP's row names a base the record does not change. Six of its rows call a substitution that changes a residue `synonymous_variant`.

### 1.5 What reaches the shape

Of the models of the release 116 RefSeq caches that print `BAM_EDIT=FAILED`, 325 on GRCh38 and 340 on GRCh37, ten alone carry an edit of their own, the same ten on both assemblies:

| Gene | Models | Strand | The edit |
| --- | --- | --- | --- |
| `PEG10` | `NM_015068.3`, `NM_001172437.2`, `NM_001184961.1` | forward | adds a position |
| `OAZ1` | `NM_004152.3`, `NM_001301020.1` | forward | removes a position |
| `OAZ3` | `NM_016178.2`, `NM_001134939.1`, `NM_001301371.1` | forward | removes a position |
| `OAZ2` | `NM_002537.3`, `NM_001301302.1` | reverse | removes a position |

Of the 4,398,127 records of genome-wide ClinVar on GRCh38, nine differ on those models, on 24 rows. The 4,388,172 records of the GRCh37 file give the same nine records and 24 rows. All nine are one-base substitutions. On each assembly `USED_REF` differs on all 24 rows, `Codons` on 18, `Amino_acids` on 15, and `Consequence` and `IMPACT` on 6.

Under `--refseq` alone every other cell is equal, the position cells among them. Six of the 24 rows are 3' UTR rows, where `USED_REF` alone differs. On three coding rows the two reference codons give one residue. On six the base Ensembl VEP read is the alternate allele, which changes the term. A record past the edit whose neighbouring base equals its own differs nowhere, and no deletion or longer substitution of either file differs on these models.

Under `--refseq --hgvs`, two more cells differ on those rows, with the same counts on GRCh38 and on GRCh37. `HGVSc` differs on all 24: Ensembl VEP's names the base it read, its `USED_REF` on the transcript's strand, as the base changed.

The position differs on all 24 as well: vep-rs counts it along the genome's exons, and Ensembl VEP writes it with the model's edit counted. The first record of section 1.2 is `c.1744` against vep-rs's `c.1743` on `NM_015068.3`, and `7 94667457 T A`, in the 3' UTR, is `c.*1776` against `c.*1775`. That is a difference in the numbering, not in the base read. [The HGVSc position past the own edit of a failed-edit RefSeq model](refseq-failed-edit-model-hgvsc-offset.md) documents it, and a comparison that names that entry's class beside this one sets the cell aside.

`HGVSp` differs on the 15 rows whose `Amino_acids` differ: each output names the residue of its own reference codon. On the first record of section 1.2 Ensembl VEP writes `NM_015068.3:c.1744T>T` and `NP_055883.2:p.Phe582%3D`, a base changed to itself and a residue kept. vep-rs writes `NM_015068.3:c.1743G>T` and `NP_055883.2:p.Val582Phe`. In JSON the record's `most_severe_consequence` differs on two of the nine records, the two whose rows differ in `Consequence`: `synonymous_variant` against vep-rs's `missense_variant` on the first record of section 1.2, and `splice_region_variant` against `missense_variant` on the second.

No committed corpus's cache holds a model with `BAM_EDIT=FAILED`, and no record of the chromosome 21 datasets reaches one: the ten models lie on chromosomes 1, 7, 15 and 19.

### 1.6 How a comparison treats these rows

The format comparator `tools/vep-formats-compare` holds the class `use_transcript_ref_cache_cdna`. It reads a `Transcript` row of a one-base allele that carries `BAM_EDIT=FAILED` in both engines' outputs. The `GIVEN_REF` is one base and the same in both, vep-rs's `USED_REF` is that base and Ensembl VEP's is another.

On such a row the class sets aside seven cells, each only where the base Ensembl VEP read explains it:

- the `USED_REF`;
- the `Codons`, one codon in each engine's output, whose reference halves differ in one upper-case base: Ensembl VEP's is its `USED_REF` and vep-rs's is the `GIVEN_REF`, on the transcript's strand;
- the `Amino_acids` each engine's `Codons` translate to under the standard genetic code;
- a `Consequence` whose two sets differ in one term, in each output the term its own two residues give, every other term shared;
- the `IMPACT`, in each output that of its own terms, beside such a `Consequence`;
- an `HGVSc` of one substitution in each output, both at one position, whose base changed is Ensembl VEP's `USED_REF` in its output and the `GIVEN_REF` in vep-rs's, on the transcript's strand;
- an `HGVSp` of one residue in each output, both at one position, each naming the residue of its own reference codon as the residue changed, where the two reference codons give different residues and neither is a stop codon.

In the JSON format the class leaves uncompared the `most_severe_consequence` of a record that holds a row whose `Consequence` it sets aside.

Every other cell of the row is compared, and so is a row of a longer reference allele. An `HGVSc` whose two positions differ stays a difference, as does an `HGVSp` beside a reference codon that is a stop codon or beside two reference codons that give one residue. The rows of a FASTA that disagrees with the cache are not told from the rows alone, and nothing is set aside for them. One such row has the shape above and is set aside with the others: a row of a model that prints `BAM_EDIT=FAILED` whose record gives the cache's base as its reference.

### 1.7 Reproducing these rows

The entry's reproduction, `scripts/divergences/reproductions/use_transcript_ref_cache_cdna.toml`, runs the masked FASTA of the third paragraph. The rows of this section need the RefSeq cache, which the tree does not hold. The commands of [Reproducing these examples](reproducing.md#reproducing-these-examples), with `--refseq` and the RefSeq GRCh38 cache of release 116, print them for the two records of section 1.2; `grep -E 'NM_015068.3|NM_016178.2'` selects them. With `--hgvs` added they print the two notations of section 1.5, and with `--json` the record's `most_severe_consequence`.
