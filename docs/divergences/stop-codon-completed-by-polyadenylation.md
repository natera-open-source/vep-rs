# A stop codon that polyadenylation completes

**The defect in one sentence.** Some coding sequences end on the transcript's last base, one or two bases into their stop codon, and the adenines of the poly(A) tail complete the codon. Ensembl VEP 116.3 reads those bases as an incomplete terminal codon: for an allele that starts in them it writes `incomplete_terminal_codon_variant` and neither `stop_lost` nor `stop_retained_variant`. vep-rs reads the completed stop codon, and the poly(A) tail after it, on a transcript whose record shows the model complete.

**Why this is a defect.** The bases are a stop codon, and the row calls them a codon whose end is unknown. The Sequence Ontology's term for a stop codon, SO:0000319, covers the case: "In genomic sequence, a stop codon may be only one or two nucleotides where the additional nucleotides are provided in the RNA via polyadenylation."

The human mitochondrial reference record, NC_012920.1, marks six of its thirteen coding sequences this way. Each carries a translation exception, the record's mark for a codon that is not read as its bases alone spell, which names the terminator at its last one or two positions. The note on each exception says that the stop codon "is completed by the addition of 3' A residues to the mRNA". Ensembl VEP's RefSeq cache of release 116 holds the same six genes with the stop codon's three bases written into the cached sequence by an edit.

**What vep-rs does instead.** vep-rs asks the transcript's record in the cache whether it shows the model complete: whether the annotation gives the coding sequence whole, up to its stop codon. Either of two things shows it:

- The codon table the record names, where that is table 2, 5, 9 or 13: the vertebrate, invertebrate, echinoderm and ascidian mitochondrial codes. Reference sequence records of mitochondrial genomes under these codes mark such an end as a stop codon that polyadenylation completes.
- A consensus coding sequence (CCDS) id on the record. The CCDS project lists coding regions annotated as full-length, and it marks one human id with the attribute `polyA required for stop codon`.

On such a transcript vep-rs looks at the coding end. The coding sequence must end one or two bases into a codon, on the transcript's last base, with no 3' UTR and no `cds_end_NF` flag. Those bases followed by adenines must spell a stop codon of the transcript's codon table. Then every allele is read on the completed message: the cached coding sequence, the adenines that complete its stop codon, and the adenines of the poly(A) tail after them.

The stop codon's terms follow as on any transcript. A change that leaves no stop codon at its place is `stop_lost`, and one that leaves a stop codon there is `stop_retained_variant`. An insertion before the codon is a frameshift or an in-frame insertion by its length. An insertion or deletion that starts inside the codon gets its stop term without `frameshift_variant`: Ensembl VEP's frameshift test stands back where the first residue touched is the stop (line 1554), and vep-rs follows it.

The `Amino_acids` cell holds the stop against the new residue, and the `Codons` cell shows the completing adenines as unchanged bases. The position cells are those of the cached coding sequence, since the adenines have no position in the genome.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the module the Ensembl VEP 116.3 image holds), in `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm` unless another file is named:

- Lines 1588-1614, `partial_codon`: lines 1601-1610 take the codon the allele starts in and count the coding bases from that codon's first base to the end of the coding sequence. One or two is an incomplete codon. The routine reads no base, no codon table and no transcript end.
- Line 1240, `stop_lost`, line 1296, `stop_retained`, and line 1548, `frameshift`: each returns 0 on that answer before it tests anything.
- Lines 1100-1128, `inframe_insertion`: the routine asks neither `partial_codon` nor the inserted length, so an insertion of one base directly before the stop codon is written `inframe_insertion`.
- Lines 1259-1260, `stop_lost`: for an allele that starts before those bases the routine reads the peptide alleles. The reference peptide ends on the unknown residue of the incomplete codon and holds no stop.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 774-776: a codon string that ends in an incomplete codon is translated with the unknown residue, `X`, at its end.

**What correcting it would change.** If `partial_codon` read such bases as the stop codon they begin, the stop codon tests would run. The rows would carry the stop codon's terms in place of `incomplete_terminal_codon_variant`, with the residue and codon cells of a stop codon, and `IMPACT` would follow the terms.

**The class it produces.** `stop_codon_completed_by_polyadenylation`: a `Transcript` row of an allele that holds a base of such a stop codon, or of an insertion directly before the codon or inside it, on a transcript whose record shows the model complete. Ensembl VEP's terms are the ones it gives a last codon held in part, and vep-rs's are the ones the allele gets on the completed message. Section 4 lists both for 21 alleles. `IMPACT`, `Amino_acids` and `Codons` follow the terms, and the position cells are equal on each of those rows.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape the record alone shows is a substitution, an insertion or a deletion, and no flag is needed. That page is the index of the documented differences, and it counts each shape over 15,646,481 chromosome 21 records. Of them 15,645,693 write a reference allele and an alternate allele in bases, 999,950 per million. The rate is an upper bound: the transcript decides, and the divergence occurs only on a transcript that ends inside its stop codon and whose record shows the model complete.

Twenty-two committed corpora hold 95 such transcripts under the four codon tables. In five of the corpora, 29 of their own input records in all hold a base of the codon, or are insertions directly before it or inside it. No committed cache holds a transcript that shows the completion by a consensus coding sequence id. The genome-wide ClinVar dataset holds no row of the class on either assembly.

The records of this page are made up for it: each is placed by hand on a transcript of a committed cache, and none comes from a dataset.

## 1. A substitution of the first of two stop codon bases, on ENST00000361390

**In short.** `MT 4261 T C` on GRCh38 changes the first base of the stop codon of ENST00000361390, a transcript of the gene MT-ND1. Its coding sequence ends at MT:4262, the transcript's last base, two bases into the codon, and the poly(A) tail supplies the third. Ensembl VEP 116.3 writes `incomplete_terminal_codon_variant,coding_sequence_variant` with the unknown residue. vep-rs writes `stop_lost` with the stop against glutamine, the residue the changed codon names.

### 1.1 What a reader needs to know first

A coding transcript ends its coding sequence with a stop codon, three bases that end translation. For a change at that codon Ensembl VEP writes one of two terms. `stop_lost` says that no stop codon stands at the stop codon's place after the change, and `stop_retained_variant` says that one does.

A messenger RNA ends in a poly(A) tail, a run of adenines added after its last transcribed base. Some messages end inside their stop codon: the genome holds the codon's first one or two bases, and the tail's adenines are the rest. In the human mitochondrial genome six of the thirteen protein-coding sequences end this way, each directly before a transfer RNA gene.

A codon table says which residue each codon names and which codons stop translation. A mitochondrial genome is read under a table of its own, and the cache record of each transcript names its table. Table 2 is the vertebrate mitochondrial code.

`incomplete_terminal_codon_variant` is Ensembl VEP's term for a change in a last codon that the coding sequence holds only in part. Such a codon is translated as the unknown residue, which the `Amino_acids` cell prints as `X`.

### 1.2 The record

A made-up substitution of one base, against the Ensembl GRCh38 cache of release 116:

```
MT	4261	substitution_of_the_first_stop_codon_base_MT_4261	T	C	.	.	.
```

ENST00000361390 is a forward-strand transcript of one exon with no 5' UTR and no 3' UTR, so a cDNA position and a coding position are the same number. Its cache record names codon table 2 and carries no `cds_end_NF` flag, the mark of a coding sequence cut short at its end.

| Part | Genomic span | Bases |
| --- | --- | --- |
| Coding sequence | MT:3307-4262 | 956 |
| Last whole codon, codon 318 | MT:4258-4260 | 3 |
| Stop codon, the bases the genome holds | MT:4261-4262 | 2 |

The coding sequence holds 956 bases, and 956 = 3 × 318 + 2. The first 954 are 318 codons, and the cached protein has 318 residues. The two bases left are coding positions 955 and 956, the first and second base of codon 319. The next base of the genome, MT:4263, is the first base of a transfer RNA gene.

The input record's reference allele gives the base at MT:4261, `T`, and the base at MT:4262 is an adenine. Followed by one more adenine they spell a stop codon of table 2. The reference record NC_012920.1 gives the coding sequence `3307..4262` the translation exception `(pos:4261..4262,aa:TERM)`, with the note quoted above.

### 1.3 What each engine writes

The row on ENST00000361390 in the default format. Ensembl VEP 116.3:

```
substitution_of_the_first_stop_codon_base_MT_4261	MT:4261	C	ENSG00000198888	ENST00000361390	Transcript	incomplete_terminal_codon_variant,coding_sequence_variant	955	955	319	X	Ta/Ca	-	IMPACT=LOW;STRAND=1
```

vep-rs:

```
substitution_of_the_first_stop_codon_base_MT_4261	MT:4261	C	ENSG00000198888	ENST00000361390	Transcript	stop_lost	955	955	319	*/Q	Taa/Caa	-	IMPACT=HIGH;STRAND=1
```

The rows differ in four cells. The `Consequence` is `incomplete_terminal_codon_variant,coding_sequence_variant` against `stop_lost`. The `Amino_acids` cell is the unknown residue against the stop and glutamine. The `Codons` cell holds two bases a side in Ensembl VEP's row and three in vep-rs's, the third an adenine in lower case. `IMPACT` is `LOW` against `HIGH`.

Both rows print the same positions: 955 for the cDNA and the coding position, and 319 for the protein position. The record writes 23 rows in each engine, and the other 22 are equal in every column.

### 1.4 Why Ensembl VEP writes that

`partial_codon` takes the codon the allele starts in, codon 319, whose first coding base is 3 × 319 - 2 = 955. The coding sequence holds 956 bases, so the codon holds 956 - 954 = 2 of them. That is fewer than three, and line 1610 answers that the codon is incomplete.

`stop_lost` returns 0 on that answer at line 1240 and `stop_retained` at line 1296. Neither reads a base. `partial_codon` is also the test behind the term `incomplete_terminal_codon_variant`, so the row carries it.

The two bases are translated as the unknown residue on both sides of the change, so the `Amino_acids` cell is `X`. `coding_sequence_variant` is the term of a coding change whose peptide holds the unknown residue (lines 1616-1649).

Nothing in these lines asks whether the transcript ends at its coding end, which table translates it, or what the two bases spell with adenines after them.

### 1.5 Why vep-rs's row is the correct one

The transcript ends at MT:4262, and its cache record gives it no 3' UTR. The message therefore reads codon 319 as the base of MT:4261, the base of MT:4262 and the first adenine of the poly(A) tail. Those three bases are a stop codon of table 2, the table the cache record names. The model is complete: 318 residues and a stop codon. The reference record says the same of the same two positions.

The input record changes the codon's first base. The new codon is the alternate base, the adenine of MT:4262 and the tail's adenine, and in table 2 it names glutamine. No stop codon stands at codon 319 after the change. The codons after it are adenines of the tail, which name lysine, so translation meets no stop codon further on. That is `stop_lost`.

The other cells say the same. The `Amino_acids` cell is the stop against glutamine. The `Codons` cell shows the changed base in upper case, and after it the adenine of MT:4262 and the tail's adenine in lower case, unchanged. The positions are those of the cached coding sequence, 955 and 319.

`incomplete_terminal_codon_variant` would be right for a codon whose third base the model does not give. Here the cache record shows the model complete, and the third base is the first adenine of the poly(A) tail.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38/json_cache \
  --species homo_sapiens --assembly GRCh38
grep ENST00000361390 out.txt
```

The committed GRCh38 cache holds the transcript. vep-rs writes `stop_lost` on the row. Ensembl VEP 116.3 takes the record with the cache options of [Reproducing these examples](reproducing.md) and its GRCh38 cache, and writes `incomplete_terminal_codon_variant,coding_sequence_variant`.

The entry's reproduction, `scripts/divergences/reproductions/stop_codon_completed_by_polyadenylation.toml`, holds this record and the twenty others of sections 2 and 4. `tools/vep-divergence-reevaluate` runs them on both engines and holds each row to both recorded values.

### 1.7 How the comparator treats this class

A comparison of the two engines' outputs that is asked to apply the class `stop_codon_completed_by_polyadenylation` reads the candidate's cache, vep-rs's own. It looks there for the transcript of the row and sees whether its record shows a stop codon that polyadenylation completes. Where it does, the comparison sets aside a row whose two term sets are the ones this page derives for such an allele, with the cells that follow the terms. It counts the row under the class and not as a difference.

A comparison that is given no cache sets nothing aside under this class, and such a row stays a difference.

## 2. A substitution of the one stop codon base, on ENST00000361381

**In short.** `MT 12137 T C` on GRCh38 changes the last base of ENST00000361381, a transcript of the gene MT-ND4. That base is the first base of its stop codon, and the poly(A) tail supplies the other two. Ensembl VEP 116.3 writes `incomplete_terminal_codon_variant,coding_sequence_variant`, and vep-rs writes `stop_lost`. The mechanism and the argument are those of section 1.

### 2.2 The record

```
MT	12137	substitution_of_the_stop_codon_base_MT_12137	T	C	.	.	.
```

ENST00000361381 is a forward-strand transcript of one exon, MT:10760-12137, with no UTR on either side, codon table 2 and no `cds_end_NF` flag. Its coding sequence holds 1,378 bases, and 1,378 = 3 × 459 + 1: 459 codons, a cached protein of 459 residues, and one base more. That base, MT:12137, is coding position 1378 and the first base of codon 460. MT:12138 is the first base of a transfer RNA gene.

Followed by two adenines the base spells a stop codon of table 2. NC_012920.1 gives the coding sequence `10760..12137` the translation exception `(pos:12137,aa:TERM)`.

### 2.3 What each engine writes

The row on ENST00000361381. Ensembl VEP 116.3:

```
substitution_of_the_stop_codon_base_MT_12137	MT:12137	C	ENSG00000198886	ENST00000361381	Transcript	incomplete_terminal_codon_variant,coding_sequence_variant	1378	1378	460	X	T/C	-	IMPACT=LOW;STRAND=1
```

vep-rs:

```
substitution_of_the_stop_codon_base_MT_12137	MT:12137	C	ENSG00000198886	ENST00000361381	Transcript	stop_lost	1378	1378	460	*/Q	Taa/Caa	-	IMPACT=HIGH;STRAND=1
```

Ensembl VEP's `partial_codon` finds one base in codon 460, since 1378 - 3 × 459 = 1, and the stop tests return 0. vep-rs reads the codon as the base of MT:12137 and two adenines of the tail. The changed codon names glutamine, and the tail after it holds no stop codon: `stop_lost`. The record writes 22 rows in each engine, and the other 21 are equal in every column.

## 3. The same shape under the invertebrate mitochondrial code, on FBtr0433501

**In short.** `mitochondrion_genome 6409 A G` on the fruit fly assembly BDGP6.54 changes the last base of FBtr0433501, a reverse-strand transcript of the gene mt:ND5 under codon table 5. The base is the first base of its stop codon. Ensembl VEP 116.3 writes `incomplete_terminal_codon_variant,coding_sequence_variant`, and vep-rs writes `stop_lost`. The mechanism and the argument are those of section 1.

### 3.2 The record

A made-up substitution of one base, against the Ensembl cache of *Drosophila melanogaster*, release 116:

```
mitochondrion_genome	6409	substitution_of_the_stop_codon_base_mitochondrion_genome_6409	A	G	.	.	.
```

FBtr0433501 is a reverse-strand transcript of one exon, 6409-8125 on the mitochondrial genome, with no UTR on either side. On this strand the transcript runs from high to low coordinates and carries the complement of the genome's bases, so its last base is 6409. Its coding sequence holds 1,717 bases, and 1,717 = 3 × 572 + 1: 572 codons and one base more, coding position 1717, the first base of codon 573.

The transcript's cache record names codon table 5, the invertebrate mitochondrial code. The transcript's base at 6409 followed by two adenines spells a stop codon of that table. The fly's mitochondrial reference record, NC_024511.2, gives this coding sequence the translation exception `(pos:complement(6409),aa:TERM)` with the same note as the human record.

### 3.3 What each engine writes

The row on FBtr0433501. Ensembl VEP 116.3:

```
substitution_of_the_stop_codon_base_mitochondrion_genome_6409	mitochondrion_genome:6409	G	FBgn0013684	FBtr0433501	Transcript	incomplete_terminal_codon_variant,coding_sequence_variant	1717	1717	573	X	T/C	-	IMPACT=LOW;STRAND=-1
```

vep-rs:

```
substitution_of_the_stop_codon_base_mitochondrion_genome_6409	mitochondrion_genome:6409	G	FBgn0013684	FBtr0433501	Transcript	stop_lost	1717	1717	573	*/Q	Taa/Caa	-	IMPACT=HIGH;STRAND=-1
```

The rows differ in the same four cells as in section 1, and both print 1717 and 573 for the positions. On the transcript the input record changes the codon's first base, and the changed codon names glutamine in table 5. The record writes 25 rows in each engine, and the other 24 are equal in every column.

### 3.6 Try it yourself

With the record of 3.2 saved as `record.vcf`:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/BDGP6.54-fly/json_cache \
  --species drosophila_melanogaster --assembly BDGP6.54
grep FBtr0433501 out.txt
```

The committed cache of the corpus `BDGP6.54-fly` holds the transcript. Ensembl VEP 116.3 takes the record with its own cache of the species and the cache options of [Reproducing these examples](reproducing.md). The entry's reproduction holds the human records alone, and the commands above run this one.

## 4. The other alleles at such a coding end

The reproduction holds 21 made-up records on the two human transcripts, those of sections 1 and 2 among them. The table gives both engines' term sets for each. Its first thirteen rows are the alleles on ENST00000361390, where the genome holds two bases of the stop codon, MT:4261 and MT:4262. Its last eight are those on ENST00000361381, where the genome holds one, MT:12137.

| Allele | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| The first base substituted | `incomplete_terminal_codon_variant,coding_sequence_variant` | `stop_lost` |
| The second base substituted | `incomplete_terminal_codon_variant,coding_sequence_variant` | `stop_lost` |
| Both bases substituted, to another stop codon of table 2 | `incomplete_terminal_codon_variant,coding_sequence_variant` | `stop_retained_variant` |
| The first base deleted | `incomplete_terminal_codon_variant,coding_sequence_variant` | `stop_lost` |
| The second base, an adenine, deleted | `incomplete_terminal_codon_variant,coding_sequence_variant` | `stop_retained_variant` |
| Both bases deleted | `protein_altering_variant,incomplete_terminal_codon_variant` | `stop_lost` |
| One base inserted before the codon | `inframe_insertion,incomplete_terminal_codon_variant,coding_sequence_variant` | `frameshift_variant` |
| Three bases inserted before the codon | `inframe_insertion,incomplete_terminal_codon_variant` | `inframe_insertion` |
| An adenine inserted between the two bases | `stop_gained,incomplete_terminal_codon_variant,coding_sequence_variant` | `stop_retained_variant` |
| A cytosine inserted between the two bases | `incomplete_terminal_codon_variant,coding_sequence_variant` | `stop_lost` |
| The last base of codon 318 and the first stop codon base deleted | `frameshift_variant` | `frameshift_variant,stop_lost` |
| Both bases and the two bases after the transcript deleted | `incomplete_terminal_codon_variant,coding_sequence_variant,3_prime_UTR_variant` | `stop_lost,3_prime_UTR_variant` |
| Both bases and the base after the transcript substituted | `incomplete_terminal_codon_variant,coding_sequence_variant,3_prime_UTR_variant` | `coding_sequence_variant,3_prime_UTR_variant` |
| The one base substituted | `incomplete_terminal_codon_variant,coding_sequence_variant` | `stop_lost` |
| The one base deleted | `protein_altering_variant,incomplete_terminal_codon_variant` | `stop_lost` |
| One base inserted before the codon | `inframe_insertion,incomplete_terminal_codon_variant,coding_sequence_variant` | `frameshift_variant` |
| Three bases inserted before the codon | `inframe_insertion,incomplete_terminal_codon_variant` | `inframe_insertion` |
| The last two bases of codon 459 and the one base deleted | `stop_lost,inframe_deletion` | `inframe_deletion` |
| Codon 459 and the one base deleted | `frameshift_variant` | `frameshift_variant,stop_lost` |
| The one base and the two bases after the transcript deleted | `incomplete_terminal_codon_variant,coding_sequence_variant,3_prime_UTR_variant` | `stop_lost,3_prime_UTR_variant` |
| The one base and the two bases after the transcript substituted | `incomplete_terminal_codon_variant,coding_sequence_variant,3_prime_UTR_variant` | `coding_sequence_variant,3_prime_UTR_variant` |

vep-rs's rows follow from the completed message:

- **A deleted adenine.** The second stop codon base of ENST00000361390 is an adenine, and the tail's adenines follow it. With it deleted the message reads the same bases as before, so the stop codon is kept: `stop_retained_variant`.
- **An insertion inside the codon.** With an adenine inserted between its two bases the codon's place still holds a stop codon of table 2. With a cytosine there it holds a codon of serine, and no stop codon follows in the tail.
- **A deletion that starts in the last whole codon.** Three bases deleted on ENST00000361381 leave the first base of codon 459 before the tail's adenines. That base and two adenines are a stop codon, so the protein loses its last residue and keeps its stop: `inframe_deletion` alone. A deletion of two bases on ENST00000361390, or of four on ENST00000361381, shifts the frame and takes a stop codon base: `frameshift_variant,stop_lost`.
- **An allele that runs past the transcript's end.** A deletion from the codon's first base takes every base the genome holds of it: `stop_lost`. A substitution that puts as many bases in place of those it takes gets `coding_sequence_variant,3_prime_UTR_variant`. Across a stop codon the genome holds whole, both engines write that row for such an allele: an allele that keeps its length and has no residues to compare gets no stop codon test (line 1447).

Ensembl VEP's rows follow from the same module. An allele that starts in the stop codon's bases is in an incomplete codon, and no stop test runs. `stop_gained` (lines 1211-1232) has no such guard, and the insertion of an adenine is read as a new stop codon in a codon of unknown residue. An insertion before the codon is `inframe_insertion` at any length.

A deletion that starts in the last whole codon is not in the incomplete codon by `partial_codon`'s count. Its reference peptide ends on the unknown residue, which is no stop. Where the alternate peptide holds no unknown residue, `stop_lost` reads the two peptides and finds no stop lost: `frameshift_variant` alone. Where the alternate peptide holds the unknown residue, `stop_lost` edits the coding sequence instead and calls any sequence that came out shorter a lost stop codon (line 1469). It does so even where the bases left spell a stop codon with the tail's adenines.

Both engines write `3_prime_UTR_variant` for an allele that crosses the last base of a transcript without a 3' UTR. That term is no part of this entry.
