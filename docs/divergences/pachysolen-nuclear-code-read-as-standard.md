# A transcript on genetic code 26, which Ensembl VEP reads under the standard code

**The defect in one sentence.** On a transcript whose cache record names genetic code 26, the nuclear code of *Pachysolen tannophilus*, Ensembl VEP translates every codon under the standard code and reads `CTG` as leucine, where the code gives alanine. A substitution whose reference or alternate codon is `CTG` therefore carries leucine in its row, and `synonymous_variant` where the other codon is a leucine codon. vep-rs reads the code as NCBI's genetic code table states it.

**Why this is a defect.** A transcript's genetic code is stated by the `codon_table` attribute of its sequence region (a chromosome, scaffold or organelle genome) and carried on the transcript's cache record. Ensembl VEP hands that id to the BioPerl module it translates with, `Bio::Tools::CodonTable`, and the module in the Ensembl VEP 116.3 container image holds 24 tables. For an id it lacks, the module warns `Not a valid codon table ID [26]` and translates under its last table, whose residues are the standard code's. NCBI's genetic code table gives code 26 one codon unlike the standard code: `CTG` is alanine. A row that reads it as leucine names a residue the transcript's own code does not give.

**What vep-rs does instead.** vep-rs holds genetic code 26 as NCBI's table gives it and translates the codons of such a transcript under it. A row then names alanine for `CTG`, and its term follows from the two residues. A transcript on any other code is unaffected.

The peptide an Ensembl cache stores for the transcript was written through the same module and carries leucine at each `CTG`. vep-rs takes a row's reference residue from the stored peptide where that differs from the codon's own translation, so it reads the stored peptide of such a transcript under code 26 first. Where a stored residue is the one the module's last table gives its codon, vep-rs translates that codon again under code 26, and it keeps any other stored residue.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, and the BioPerl files of the Ensembl VEP 116.3 image):

- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 1146-1164, `_codon_table`: returns the transcript's cached `codon_table` with no check of the value.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 712 and 736, `peptide`: translates the allele's codon under that id. The reference residue comes from the same routine (`Utils/VariationEffect.pm`, lines 786 and 800-810), not from the stored peptide.
- `Bio/PrimarySeqI.pm`, line 619: builds the table with `Bio::Tools::CodonTable->new( -id => $codonTableId)`.
- `Bio/Tools/CodonTable.pm`, lines 203-249: the module's 24 names and tables. The last is `Strict`, whose residues are the standard code's and whose one start codon is `ATG`.
- `Bio/Tools/CodonTable.pm`, lines 328-338, `id`: an id with no table is warned about and set to 0.
- `Bio/Tools/CodonTable.pm`, lines 412-460, `translate`: reads `$TABLES[$id-1]`, which for id 0 is `$TABLES[-1]`, the last table.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1077-1083 and 1085-1098: `synonymous_variant` holds where the two residues are equal, `missense_variant` where they differ.

**What correcting it would change.** Corrected, Ensembl VEP would name alanine on every row on a `CTG` codon of such a transcript, and the row's term would follow from the residues the code gives. A row whose two codons hold no `CTG` stays as it is, and so does every row of a transcript on another code.

**The class it produces.** `pachysolen_nuclear_code_read_as_standard`, on a transcript whose cache record names genetic code 26 and a row whose `Codons` cell holds `CTG` as its reference or its alternate codon. A substitution of one base links `CTG` to a codon of leucine, methionine, valine, proline, glutamine or arginine, and to no alanine or stop codon. Where the other codon is a leucine codon, Ensembl VEP writes `synonymous_variant` with one residue and vep-rs `missense_variant` with two. Where it codes another residue, both write `missense_variant` and the `Amino_acids` cell alone differs. These are the rows of a codon inside the coding sequence; the entry records no row of a transcript's first codon.

The SNP/indel comparator (`tools/vep-outputs-compare`) compares tuples: a row's Location, Allele, Feature, Feature_type and set of `Consequence` terms. Two tuples, one from each engine, that differ in the set of terms alone are a pair. The comparator's 116.3 registry holds a pair rule under the class's name, and the structural-variant comparator applies it too. The rule sets aside both tuples of a pair whose sets differ by Ensembl VEP's `synonymous_variant` against vep-rs's `missense_variant` alone. It reads the `Amino_acids` and the `Codons` of the two rows to tell the pair from any other; [the class taxonomy](../concordance-provenance/116.2-divergence-classes.md#how-a-row-on-genetic-code-26-is-recognised) lists its conditions.

A row on which both engines write `missense_variant` has one tuple in both outputs, so it is no pair: nothing is set aside for it and nothing is charged. The format comparator (`tools/vep-formats-compare`) holds no class of this name and reports each differing field of such a row. Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

**Where the shape is found.** The Ensembl Genomes release 63 cache of *Pachysolen tannophilus* NRRL Y-2460 names genetic code 26 on 5,986 of its 5,991 transcripts; the other five, RNA genes, name the standard code. From the second codon on, the 5,986 coding sequences hold 15,128 `CTG` codons, and the cache stores leucine for every one. The public coding-sequence record of one of those proteins, `ODV97712`, declares `/transl_table=26` and holds alanine at residue 38, the place of its one `CTG` codon. No cache committed under `tests/golden/` holds a transcript on the code.

## 1. Three substitutions on a synthetic transcript of six codons

**In short.** The sequence `contig` holds one forward gene, `CODEG1`, of six codons from base 1,001 with no untranslated region. Its transcript `CODET1` names genetic code 26. Its third codon is `CTG` and its fourth `CTA`. Ensembl VEP 116.3 reads `CTG` as leucine on all three records below; vep-rs reads it as alanine.

**The records.** Each is a substitution of one base, run with default flags. Neither engine is given a reference FASTA: the cache record carries the coding sequence.

| Record | Base of `contig` | Alleles | The codon it changes |
| --- | --- | --- | --- |
| `alanine_codon_to_proline` | 1,008 | `T` to `C` | codon 3, `CTG` to `CCG` |
| `alanine_codon_to_leucine` | 1,009 | `G` to `A` | codon 3, `CTG` to `CTA` |
| `leucine_codon_to_alanine` | 1,012 | `A` to `G` | codon 4, `CTA` to `CTG` |

**Ensembl VEP 116.3.** The columns are Location, Consequence, CDS_position, Protein_position, Amino_acids, Codons and IMPACT, the records in the table's order:

```text
contig:1008 | missense_variant | 8 | 3 | L/P | cTg/cCg | MODERATE
contig:1009 | synonymous_variant | 9 | 3 | L | ctG/ctA | LOW
contig:1012 | synonymous_variant | 12 | 4 | L | ctA/ctG | LOW
```

**vep-rs.** The same columns:

```text
contig:1008 | missense_variant | 8 | 3 | A/P | cTg/cCg | MODERATE
contig:1009 | missense_variant | 9 | 3 | A/L | ctG/ctA | MODERATE
contig:1012 | missense_variant | 12 | 4 | L/A | ctA/ctG | MODERATE
```

**Why vep-rs's rows are the correct ones.** The gene's first base is base 1,001, so codon n starts at base 1,001 + 3 x (n - 1). Base 1,008 is the second base of codon 3 and base 1,009 its third; base 1,012 is the third base of codon 4. NCBI's table for genetic code 26 reads `CTG` as alanine, `CTA` as leucine and `CCG` as proline. The three records therefore change alanine to proline, alanine to leucine and leucine to alanine. Each changes the residue and keeps the peptide's length, which is a missense change.

**The warning Ensembl VEP writes.** For substitutions, Ensembl VEP writes the module's warning, `Not a valid codon table ID [26]`, once for each distinct codon a run translates under the id. It goes to standard error, or to the warnings file under `--quiet`. vep-rs writes no warning for the transcript.

**Reproducing it.** The entry's reproduction, `scripts/divergences/reproductions/pachysolen_nuclear_code_read_as_standard.toml`, holds the three records and both engines' rows. Beside it, `inputs/pachysolen_nuclear_code_read_as_standard/` holds `make_cache.pl`, which lists the gene's six codons, and the `json_cache` vep-rs reads. The prepare steps build the synthetic cache with that script inside the Ensembl VEP image, convert it with `tools/vep-storable-to-json` and check that the `json_cache` equals the conversion. `tools/vep-divergence-reevaluate`, given `--entry pachysolen_nuclear_code_read_as_standard`, runs both engines and holds each row to the recorded one. It needs Docker, `cargo` and the Ensembl VEP image, named with `--image`.
