# Frameshift, extension, duplication and deletion peptides on a transcript whose genetic code is not the standard one

**The defect in one sentence.** On a transcript whose sequence region names a genetic code other than the standard one, Ensembl VEP translates the reference peptide with that code but the alternate sequence of a frameshift, a lost stop, a duplication test or an empty-window deletion with the standard code. Its HGVSp for those four shapes therefore names residues and stop distances the record's codons do not encode. vep-rs translates both with the transcript's table, so its HGVSp for the four shapes names the residues the transcript's own code gives. An empty-window deletion is an in-frame deletion whose deleted bases name no whole residue of the reference peptide (section 1.4).

**Which transcripts it reaches.** Every transcript whose sequence region names a genetic code other than the standard one, a nuclear transcript as much as an organelle's. The transcripts of the human `MT`, under table 2, are such transcripts. So is a nuclear transcript whose sequence region names another code: every transcript of the committed Spironucleus salmonicida cache is under table 6, the ciliate, dasycladacean and hexamita nuclear code. The class's name says `organelle`, but its rule reads the transcript's codon table and never the kind of sequence region.

**Why this is a defect.** The transcript's codon table is a property of the transcript, stated by the cache's `codon_table` attribute of its sequence region. Ensembl VEP honours it on the same row everywhere else: the reference peptide it compares against, the `peptide` method behind every substitution's HGVSp, the consequence terms and the codons. Four routines translate an alternate sequence with no table argument, so they fall back to the standard code. Under the vertebrate mitochondrial code `AGA` and `AGG` are stop codons, `ATA` is methionine and `TGA` is tryptophan; a frameshift peptide that reads those codons under the standard code describes a protein the mitochondrial ribosome does not make.

**What vep-rs does instead.** vep-rs translates both sides of every HGVSp path with the transcript's table, so a frameshift, extension, duplication or deletion string on such a transcript names the residues the transcript's own code gives. A missense, synonymous or nonsense HGVSp is identical on both engines, and a transcript under the standard code is unaffected.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`):

- line 2333, `_get_fs_peptides`: the frameshift peptide, `$alt_cds->translate()->seq()` with no codon table.
- line 2450, `_check_for_peptide_duplication`: the upstream reference translation an inserted peptide is tested against, translated bare.
- line 2492, `_stop_loss_extra_AA`: the `fsTerN` and `extTerN` count, translated bare.
- line 2555, `_get_del_peptides`: the deletion peptide, translated bare.
- lines 712-736, `peptide`, and `_codon_table`: the substitution path, which passes the transcript's table; `_peptide` translates the reference through `Transcript::translate`, which honours it too.

**What correcting it would change.** The class vanishes: every HGVSp on such a transcript follows the transcript's code, as the row's codons, terms and reference peptide already do.

**The class it produces.** `organelle_peptide_paths_translation_table`, documented in the manifests of the corpora section 1.7 lists, under `field_divergences` for the `HGVSp` of the rows that change, and applied by the golden harness and the format comparator from the manifest. For an input without a manifest the format comparator decides the class from the transcript's codon table in the candidate cache and the shape of the two values, as section 1.7 states. The concordance figures score consequence sets and never see the field, so nothing is set aside in them.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A frameshift in MT-ND1 whose new frame opens with a mitochondrial stop

**In short.** Deleting the eighth base of MT-ND1 turns codon 3 into `AGG`, a stop under the vertebrate mitochondrial code. Ensembl VEP reads `AGG` as arginine under the standard code and runs twelve codons on to a `TAA`: `p.Met3ArgfsTer12`. vep-rs translates the new frame with the transcript's table and writes `p.Met3Ter`. The consequence terms, codons, positions and every other column agree; only the HGVSp differs. The same difference reaches every transcript whose sequence region names a genetic code other than the standard one, nuclear transcripts among them; section 1.7 lists the committed corpora that carry it.

### 1.1 What a reader needs to know first

A transcript's cache record carries its codon table ([The codon table of each sequence region](../cache-setup.md#the-codon-table-of-each-sequence-region)). The human and mouse `MT` are table 2, the fly's `mitochondrion_genome` table 5 and yeast `Mito` table 3. HGVS describes a frameshift as the first changed residue and the distance to the first stop of the new frame (`p.Met3ArgfsTer12`), and an immediate stop as `Ter` (`p.Met3Ter`). A lost stop is `extTerN`, and an insertion of residues equal to those before it a duplication (`p.Met31dup`).

### 1.2 The record

A hand-written GRCh38 record against the Ensembl 116 human cache, with `--hgvs` (Ensembl VEP needs `--fasta`; vep-rs needs a cache carrying `MT`, such as the `GRCh38-native-builder` golden corpus's):

```
MT	3313	.	AT	A	.	.	.
```

MT-ND1 (`ENST00000361390`, forward strand, `MT:3307-4262`, codon table 2) opens `ATA CCC ATG GCC AAC CTC CTA CTC CTC ATT GTA CCC ATT CTA`. The record deletes c.8, the `T` of codon 3 (`ATG`, Met). The new frame reads `ATA CCC AGG CCA ACC TCC TAC TCC TCA TTG TAC CCA TTC TAA`: codon 3 is `AGG`.

### 1.3 What each engine writes

The transcript's row, default format. Ensembl VEP:

```
MT_3313_AT/A	MT:3314	-	ENSG00000198888	ENST00000361390	Transcript	frameshift_variant	8	8	3	M/X	aTg/ag	-	IMPACT=HIGH;STRAND=1;HGVSc=ENST00000361390.2:c.8del;HGVSp=ENSP00000354687.2:p.Met3ArgfsTer12
```

vep-rs:

```
MT_3313_AT/A	MT:3314	-	ENSG00000198888	ENST00000361390	Transcript	frameshift_variant	8	8	3	M/X	aTg/ag	-	IMPACT=HIGH;STRAND=1;HGVSc=ENST00000361390.2:c.8del;HGVSp=ENSP00000354687.2:p.Met3Ter
```

Both engines write 21 rows for the record; the other 20, the neighbouring tRNA, rRNA and protein-coding transcripts' rows, are identical.

### 1.4 Why Ensembl VEP writes that

`_get_fs_peptides` builds the alternate coding sequence and calls `translate()` on it with no codon table, so BioPerl translates it under table 1, while the reference it walks against is `_peptide`, the transcript's own translation under table 2. The walk from residue 3 compares the table-2 reference `M` with the table-1 alternate `R`, `AGG` being arginine in the standard code, so the first changed residue is `Arg`. The count then comes from `_stop_loss_extra_AA`, which translates bare again: the first `*` of the table-1 translation is the `TAA` at codon 14, so the count is 14 - 2 = 12.

On most rows of real data the count is not a number at all. `_stop_loss_extra_AA` counts from the first stop of the whole standard-code translation, wherever it lies, and drops a count that is not positive (lines 2503-2508 and 2519-2525). A codon before the changed residue that the standard code reads as a stop, and the transcript's table as a residue, supplies that first stop. Ensembl VEP then writes the stop distance of a frameshift after that codon as unknown, `?`, whatever the new frame holds. Over the mitochondrial records of the ClinVar dataset, on either assembly, 21 of the 24 differing rows take this path; on the other three no such codon precedes the variant.

The same bare `translate()` sits in `_check_for_peptide_duplication` and in `_get_del_peptides`. In the first, the reference sequence translated under table 1 is what an inserted peptide is tested against, so an inserted `ATA` Met after an `ATA` codon is compared with a table-1 `I` and reported as an insertion rather than a duplication.

The second serves the empty-window deletion: an in-frame deletion whose deleted bases name no whole residue of the reference peptide, the shape of a deletion that removes an intron with parts of both exons around it. Lines 2167-2173 send a deletion there when its reference peptide is empty, and it rebuilds both peptides from the variant's residue to the protein's end, the altered one bare (line 2555) and the reference through the transcript's own translation (line 2559). The four are the only translations in the module that pass no table.

Source (`ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`): lines 2333, 2450, 2492 and 2555, against `peptide` at lines 712-736.

### 1.5 Why vep-rs's row is the correct one

The record's own codons decide: after the deletion codon 3 of MT-ND1 is `AGG`, and MT-ND1 is translated by the mitochondrial ribosome under table 2, where `AGG` is a termination codon. There is no arginine at residue 3 of the altered peptide and no ten residues after it; `p.Met3Ter` is the description HGVS gives a frame whose first changed codon is a stop.

Ensembl VEP agrees that MT-ND1 reads under table 2 everywhere else on the row: its reference peptide for the comparison is the table-2 one, and a substitution in the same codon, `MT 3314 T G` (`ATG` to `AGG`), it reports as `stop_gained`, `p.Met3Ter`, through `peptide`. Only the four bare translations disagree with the rest of its own annotation.

### 1.6 Try it yourself

Annotate these GRCh38 records with `--hgvs` and the human 116 cache. HGVSp on `ENST00000361390` (MT-ND1), `ENST00000361453` (MT-ND2), `ENST00000361624` (MT-CO1) and `ENST00000361681` (MT-ND6):

```
vep -i records.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --species homo_sapiens --assembly GRCh38 --hgvs
```

| Input | Shape | Ensembl VEP | vep-rs |
| --- | --- | --- | --- |
| `MT 3313 . AT A` | frameshift, `AGG` stop at codon 3 | `p.Met3ArgfsTer12` | `p.Met3Ter` |
| `MT 3392 . GC G` | frameshift, `ATA` at codon 30 | `p.Tyr30IlefsTer11` | `p.Tyr30MetfsTer11` |
| `MT 3563 . GA G` | frameshift, `TGA` (Trp) before codon 87 | `p.Thr87ProfsTer?` | `p.Thr87ProfsTer14` |
| `MT 3313 . A AT` | frameshift, `AGG` stop at codon 29 | `p.Met3IlefsTer49` | `p.Met3IlefsTer27` |
| `MT 3400 . C CG` | frameshift | `p.Gln32ArgfsTer20` | `p.Gln32ArgfsTer5` |
| `MT 3399 . A AATA` | duplicated `ATA` (Met31) | `p.Met31_Gln32insMet` | `p.Met31dup` |
| `MT 3564 . A ATGA` | duplicated `TGA` (Trp86) | `p.Trp86_Thr87insTrp` | `p.Trp86dup` |
| `MT 6000 . CT C` | frameshift in MT-CO1 | `p.Leu33GlnfsTer?` | `p.Leu33GlnfsTer16` |
| `MT 14400 . C CC` | frameshift in MT-ND6 (reverse strand) | `p.Val92GlyfsTer19` | `p.Val92GlyfsTer10` |
| `MT 3464 . TA T` | frameshift into a `TAA`-reading frame | `p.Lys54AsnfsTer8` | `p.Lys54AsnfsTer8` |
| `MT 7443 . A G` | lost stop, no 3' sequence | `p.Ter514GlyextTer?` | `p.Ter514GlyextTer?` |
| `MT 4500 . TC T` | frameshift to a `TAA` under both codes | `p.Ser11LeufsTer12` | `p.Ser11LeufsTer12` |

The last three agree: where the new frame's first difference and first stop are codons the two codes read alike, the four paths give the same string under either table. A human `MT` transcript has no 3' UTR, so a lost stop runs off the transcript under both codes (`extTer?`). The extension shape changes only on a transcript under another code that has downstream sequence, and the empty-window deletion shape only on one with an intron, which no human `MT` transcript has.

### 1.7 How the comparator treats this class

The golden harness (`crates/vep-cli/tests/common/mod.rs`) and the format comparator (`tools/vep-formats-compare`) apply `organelle_peptide_paths_translation_table` from a corpus manifest. A `field_divergences` entry of the class names the row and its `HGVSp` with the value vep-rs prints, and the comparison checks vep-rs's value against it. The two outputs alone do not decide the class, since the rows do not state the transcript's codon table.

For an input without a manifest the format comparator reads that table from the candidate cache, named by `--json-cache`, and sets aside a differing `HGVSp` when all of these hold:

- The row's transcript is filed in the cache under a codon table that reads some codon unlike the standard one; tables 1 and 11 read every codon alike.
- Both values name the same protein.
- The row's consequence set and the two descriptions place the difference in a frameshift's description, a lost stop's extension count, or a duplication that one engine writes as an insertion.
- The two rows agree in everything else that states the coding change: the consequence set, the cDNA, CDS and protein positions, the amino acids, the codons, `HGVSc` and `HGVS_OFFSET`.

A row that fails one of them stays a difference: the same pair on a transcript under the standard table, a substitution on a transcript under another table, a value one engine leaves empty. So does a row of the fourth shape, the empty-window deletion, which the third condition does not name. Without a manifest the comparison checks the shape of the pair, not the value vep-rs prints, which only a manifest row records.

One more row is refused by a rule of its own: an insertion of two or more bases on a reverse-strand transcript that is shifted past its own length by a length that is no multiple of it. Its `HGVSp` differs for another documented reason ([Allele on a coding row after a shift](../output/vep-and-tab.md#allele-on-a-coding-row-after-a-shift)).

The committed corpora whose manifests document the class, and the rows each run of theirs documents:

| Corpus | Codon table | Sequence regions of the cache under it | Rows by run |
| --- | --- | --- | --- |
| `ARS-UCD2.0-cattle` | 2 | `MT` | `everything` 2, `hgvs-sift` 2 |
| `ARS-UI_Ramb_v3.0-sheep` | 2 | `MT` | `everything` 1, `hgvs-sift` 1 |
| `ASM294v2-fission-yeast` | 4 | `MT` | `everything` 1, `hgvs` 1 |
| `AgamP4-anopheles` | 5 | `Mt` | `everything` 2, `hgvs-sift` 2 |
| `BDGP6.54-fly` | 5 | `mitochondrion_genome` | `everything` 4, `hgvs` 4 |
| `GCA000497125v2-spironucleus` | 6 | all 3 | `everything` 20, `hgvs` 20 |
| `GRCz11-zebrafish` | 2 | `MT` | `everything` 2, `hgvs-sift` 2 |
| `R64-1-1-yeast` | 3 | `Mito` | `everything` 14, `hgvs` 14 |

Every documented row is the `HGVSp` of a frameshift: no corpus documents a row of the lost-stop, duplication or empty-window deletion shape. The `GCA000497125v2-spironucleus` rows stand on nuclear transcripts. Its cache names two chromosomes and one unplaced sequence, no organelle among them, and files every transcript under table 6. That table reads `TAA` and `TAG` as glutamine, where the standard code reads a stop. The other corpora's rows stand on mitochondrial transcripts.

The `GRCm39-mouse` corpus's `MT` rows are substitutions or shapes the two codes read alike, and the Arabidopsis `Mt` is table 1 and `Pt` table 11, whose readings the standard code shares on those rows, so no row of theirs changes.

On real data the class is counted over the 3,122 mitochondrial records of the ClinVar dataset (table 2). Under `--hgvs` 24 rows differ on GRCh37 and the same 24 on GRCh38, every one a frameshift: 21 with a `?` count on Ensembl VEP's side and 3 with two numbers. The 6 lost-stop rows, the 16 insertion or duplication rows and the one in-frame deletion row of those records agree on both assemblies.
