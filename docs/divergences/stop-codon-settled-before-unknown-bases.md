# A stop codon that an allele's known bases settle before its unknown bases

**The defect in one sentence.** An alternate allele can hold unknown bases, written `N`, after known bases that spell a stop codon in the reading frame; the protein then ends at that stop whatever the unknown bases are. Ensembl VEP translates no allele that holds a letter other than A, C, G and T, so it writes no residues, no `stop_gained` and no HGVSp for it. vep-rs writes the row that every completion of the unknown bases shares. Where the stop codon is the transcript's own, kept at its place by an insertion, the term Ensembl VEP leaves out is `inframe_insertion`.

**Why this is a defect.** A record can spell an inserted sequence known at its ends and unknown in between, as ClinVar does for long insertions. Where the codons before the first unknown base hold a stop codon, translation ends there: no reading of the later bases reaches the protein. Ensembl VEP's own rules give one answer for every completion of such an allele, `stop_gained` beside the term its length gives, with `IMPACT=HIGH`. For the allele as written they give `frameshift_variant` alone, which is `IMPACT=HIGH` as well, or `coding_sequence_variant` with `IMPACT=MODIFIER` where the length is whole codons. A row of `coding_sequence_variant` files a premature stop with the changes whose effect is unknown.

**What vep-rs does instead.** vep-rs reads the allele's codons as it builds them for `Codons`. Where a stop codon stands among the whole codons before the first codon that holds a letter other than A, C, G and T, it translates the allele as any other. A codon that holds an unknown base reads `X`, or the residue its readings share. `Amino_acids`, `stop_gained`, the in-frame terms and `IMPACT` follow as for an allele of known bases. An insertion whose first bases keep the transcript's stop codon at its place is read the same way and takes `inframe_insertion` beside `stop_retained_variant` (section 3).

Under `--hgvs` the HGVSp is written where a stop codon is settled in the same way at the allele's most 3' placement, with its `HGVS_OFFSET`. An allele the record writes in an intron is read there too, where the 3' rule places it in the coding sequence (section 2.1). The rule reads codons, so it holds on either strand, and a gained stop is read for a deletion-insertion as for an insertion. Section 5 lists what it leaves as it is.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the modules the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/Variation/VariationFeatureOverlapAllele.pm`, line 74: the unambiguous letters are A, C, G, T and `-`. Lines 293-302, `seq_is_unambiguous_dna`, test an allele against them.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, line 701, in `peptide`: an allele that fails that test has no peptide, wherever its other letter stands. Line 815, in `codon`, admits the ambiguity letters, so the codons are built.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1211-1232, `stop_gained`: line 1228 looks for a stop in the alternate peptide, and is not reached without that peptide. Lines 1100-1159, `inframe_insertion`, need both peptides.
- The same file, lines 1541-1586, `frameshift`: line 1564 reads the lengths alone, so the term stands. Lines 1616-1649, `coding_unknown`: line 1627 holds where the allele has no peptide.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1683-1684, in `hgvs_protein`: no notation for an allele with a letter other than A, C, G, T and `-`. Lines 1396-1397 are the same test in `hgvs_transcript`.
- `Bio/EnsEMBL/Variation/VariationFeatureOverlap.pm`, lines 427-447, and the text of `Bio/EnsEMBL/Variation/TranscriptVariation.pm`, lines 85-90: the constructor option `-DISAMBIGUATE_SINGLE_NUCLEOTIDE_ALLELES` reads an ambiguity letter as each base it stands for, on an allele of one nucleotide alone. The source's reason for stopping there is the number of readings a longer allele has. A stop codon before the first unknown codon needs no such count: it is read once.

**What correcting it would change.** `peptide` would translate an allele whose codons hold a stop codon before the first codon with another letter, and `hgvs_protein` would describe it. The rows of such an allele would then read as the rows of its completions do. An allele whose known bases settle no stop codon would stay untranslated.

**The class it produces.** `stop_codon_settled_before_unknown_bases`: a `Transcript` row with the same `Codons` cell in both engines, whose alternate half holds a letter other than A, C, G and T after a stop codon. Ensembl VEP writes no `Amino_acids`; vep-rs writes the translation of that half. vep-rs adds `stop_gained` to a row both engines call `frameshift_variant`. Where Ensembl VEP writes `coding_sequence_variant`, vep-rs writes `stop_gained` in its place, with an in-frame term where one applies. Where Ensembl VEP writes `stop_retained_variant` without making its stop codon test, vep-rs writes the same terms in that term's place: `stop_gained` with `frameshift_variant`, with an in-frame term or alone. Section 4 gives such a row. Section 1.7 says what each comparison reads.

On a row of an insertion that keeps the stop codon, vep-rs writes `inframe_insertion` beside `stop_retained_variant`, with the residues, where Ensembl VEP writes `stop_retained_variant` alone, `frameshift_variant` or `coding_sequence_variant`.

On a row of an allele written in an intron (section 2.1) only the HGVSp and the `HGVS_OFFSET` differ: vep-rs alone writes them.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a record whose alternate allele holds a letter other than A, C, G and T. The five chromosome 21 record sets the index calls the cuts hold four such records, two in each ClinVar set, and on none of them is a stop codon settled. The index's rate, 0.26 per million records, is therefore an upper bound.

The genome-wide ClinVar files hold 668 such records on GRCh38 and 665 on GRCh37, every one with `N` as its only other letter. A stop codon is settled so on 143 and 139 of them, on 2,416 and 551 transcript rows; the GRCh38 cache holds several times as many transcripts of a gene. On 9 further rows, of 8 of the 668 GRCh38 records, and on 3 rows of 3 of the 665 GRCh37 records, the stop codon is the transcript's own, kept as section 3 describes. The committed corpus `tests/golden/116/GRCh38` holds one record with a gained stop, line 783 of its `variants.vcf`, with three rows.

On 14 rows of 2 more GRCh38 records and on 2 rows of 1 more GRCh37 record the allele is written in an intron, and vep-rs writes the notation read in the coding sequence (section 2.1).

Section 4 describes a part of the 2,416 and 551 rows. On 36 of those rows of GRCh38 and 8 of GRCh37, on 6 and 4 records, Ensembl VEP writes `stop_retained_variant`.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Twelve inserted bases, four of them unknown, after a codon of ENST00000622113

**In short.** A record written for this page inserts twelve bases after 21:10567700 on GRCh38: a glutamate codon, a stop codon, four unknown bases and two known ones. Ensembl VEP 116.3 writes `coding_sequence_variant` with `IMPACT=MODIFIER` and no residues. vep-rs writes `stop_gained,inframe_insertion` with `IMPACT=HIGH`, the residues `E*XX` and, under `--hgvs`, `p.Asn141_Ile142insGluTer`. Ensembl VEP writes vep-rs's terms, impact and notation itself once the four bases are known.

### 1.1 What a reader needs to know first

A VCF alternate allele is letters. A, C, G and T name a base; `N` and the other ambiguity letters say that the base is one of several. Both engines accept such an allele, place it on the transcript and print its codons in `Codons`. The reading frame runs from the transcript's first coding base, so the codons before the first unknown base are known codons.

A stop codon among them ends the protein. What follows it is not translated, so its bases, known or unknown, change nothing in the protein. A completion of the allele is the same record with each unknown base replaced by a known one.

### 1.2 The records

```
21	10567700	stop_then_unknown_inframe	C	CGAATAGNNNNCA	.	.	.
21	10567700	stop_then_completed_inframe	C	CGAATAGAAAACA	.	.	.
21	10567700	unknown_before_any_stop	C	CGAAGCTNNNNCA	.	.	.
```

Position 21:10567700 is the last base of codon 141, an asparagine, of the forward-strand ENST00000622113, between `c.423` and `c.424`. The first record inserts `GAA`, `TAG`, four `N` and `CA` after it. The second is a completion of the first, each `N` an `A`. The third holds two sense codons, `GAA` and `GCT`, before its unknown bases, so no stop codon is settled.

Both engines run `--hgvs` on the release 116 cache and the reference of the corpus `tests/golden/116/GRCh38-identifiers`.

### 1.3 What each engine writes

The rows of ENST00000622113. Ensembl VEP 116.3:

| Record | Consequence | Amino_acids | HGVSp |
| --- | --- | --- | --- |
| `stop_then_unknown_inframe` | `coding_sequence_variant` | `-` | `-` |
| `stop_then_completed_inframe` | `stop_gained,inframe_insertion` | `-/E*KT` | `p.Asn141_Ile142insGluTer` |
| `unknown_before_any_stop` | `coding_sequence_variant` | `-` | `-` |

vep-rs:

| Record | Consequence | Amino_acids | HGVSp |
| --- | --- | --- | --- |
| `stop_then_unknown_inframe` | `stop_gained,inframe_insertion` | `-/E*XX` | `p.Asn141_Ile142insGluTer` |
| `stop_then_completed_inframe` | `stop_gained,inframe_insertion` | `-/E*KT` | `p.Asn141_Ile142insGluTer` |
| `unknown_before_any_stop` | `coding_sequence_variant` | `-` | `-` |

`IMPACT` is `MODIFIER` beside `coding_sequence_variant` and `HIGH` beside `stop_gained`. The HGVSp is on `ENSP00000482040.1`. Both engines print `-/GAATAGNNNNCA` in `Codons` for the first record, and neither writes an HGVSc for an allele with an unknown base. The engines agree on every cell of the second and third records and on the position columns of all three.

### 1.4 Why Ensembl VEP writes that

`peptide` returns nothing for the first record's allele at line 701: it holds `N`. `stop_gained` and `inframe_insertion` need the alternate peptide and are false. `frameshift` compares lengths: twelve bases on an empty span are whole codons, so it is false too. `coding_unknown` then holds, since the allele has no peptide, and the row is `coding_sequence_variant`. `hgvs_protein` returns at lines 1683-1684.

The second record's allele passes the test of line 701. Its peptide is `E*KT`, which holds a stop where the reference peptide holds none, so line 1228 gives `stop_gained`. The third record fails line 701 as the first does.

### 1.5 Why vep-rs's row is the correct one

The insertion stands between two codons, so its first three bases are a codon of their own. `GAA` is a glutamate, residue 142 of the changed protein. `TAG` is a stop codon: the protein ends after 142 residues. The four unknown bases and the two after them lie past that stop. Whatever they are, the protein is the reference protein up to Asn141, a glutamate and a stop.

That is what `stop_gained` says, and what `p.Asn141_Ile142insGluTer` describes. Ensembl VEP writes both for the completed record, and writes them for every other completion as well: each differs from this one only in bases past the stop. `X` in `E*XX` marks the two codons that hold an unknown base. Such a codon reads a residue where every reading of it gives that one, and `B` or `Z` where its readings are aspartate and asparagine, or glutamate and glutamine.

The third record has no such answer. A completion of its unknown bases can spell a stop codon there or not, so its protein depends on them, and both engines leave the change unknown.

The frame of the bases after the unknown run, and the term the length gives, rest on the written count of unknown bases. Both engines read every allele's length as written. The gained stop, its place and the residues before it do not rest on that count.

### 1.6 Try it yourself

With the records of 1.2 saved as `records.vcf`, the command as vep-rs takes it:

```
vep -i records.vcf -o out.txt --offline --hgvs \
  --json_cache tests/golden/116/GRCh38-identifiers/json_cache \
  --fasta tests/golden/116/GRCh38-identifiers/reference.fa.gz
grep ENST00000622113 out.txt
```

Ensembl VEP 116.3 takes the same records and `--hgvs` with the cache options of [Reproducing these examples](reproducing.md). The entry's reproduction, `scripts/divergences/reproductions/stop_codon_settled_before_unknown_bases.toml`, runs both engines on these records, on those of section 2, on those of section 3 and on the two of section 4.

### 1.7 How the comparisons treat this class

The format comparator (`tools/vep-formats-compare`) holds a row of the class to its `Codons` cell. The cell is the same on both sides, a stop codon of the standard genetic code stands before the first codon with another letter, and vep-rs's residues are the translation of the alternate half. On such a row it sets aside the `Consequence`, the `Amino_acids`, an `IMPACT` that is the one each side's own terms give (another value stays a difference), an HGVSp that names a stop where Ensembl VEP writes none, and the `HGVS_OFFSET` beside it. It sets aside a `DownstreamProtein`, `ProteinLengthChange`, `NMD` or LoFTEE cell as well, where each side's terms decide whether the plugin writes it.

A row of section 2.1 has neither a `Codons` cell nor an HGVSc and is read so:

| A row of section 2.1 | How it is read |
| --- | --- |
| The format comparator | From the transcript's exons and coding sequence in the cache of the vep-rs run under comparison, which a comparison names with `--json-cache`; without it the HGVSp and the `HGVS_OFFSET` stay differences. The row is an insertion's, with the same terms on both sides and no cDNA, CDS or protein position. Ensembl VEP's row holds no HGVSc, no HGVSp and no offset; vep-rs's holds an HGVSp that names a stop, an offset and no HGVSc. The offset carries the transcript's strand as its sign. Its size is smaller than the count of known bases the allele opens with on the transcript, so the shift stops at a known base, and it places the allele, written outside every exon, inside one exon, in the coding sequence, where that exon's coding bases before the point are the last bases the shift passed and the coding base behind the point is not the allele's next base. Read from that point, the allele's codons hold a stop codon before the first codon with another letter, and a frameshift notation names and numbers the first changed residue, names the residue it becomes and counts to that stop; any other notation is held to the tests before it. The HGVSp and the `HGVS_OFFSET` are then set aside. |
| The tuple comparison and the plugin tuple scorer | The terms are equal on such a row, so neither reads a pair there. |

The tuple comparison (`tools/vep-outputs-compare`) scores consequence terms, and reads the two rows behind a pair whose terms differ as the class's do. It sets the pair aside where the two rows hold the same CDS position and the same `Codons` cell with such a stop codon, Ensembl VEP's row holds no residues and vep-rs's residues hold the stop at that codon's place. A pair whose rows end before the `Codons` column stays a difference. The manifest of `tests/golden/116/GRCh38` lists the three rows of its record under the class, and the corpus comparison sets them aside. A row off the shape stays a difference charged to vep-rs.

Both comparisons set a row of section 3 aside as well, held to the same `Codons` cell and residues as a row of section 1:

| Ensembl VEP's terms on a row of section 3 | How the row is read |
| --- | --- |
| `stop_retained_variant`, to which vep-rs adds `inframe_insertion` | Under this class. |
| `coding_sequence_variant` where vep-rs writes the two terms | Under this class, on a row of an insertion between two neighbouring positions that carries both cDNA positions and a CDS start. |
| `frameshift_variant` where vep-rs writes the two terms | On such an insertion's row too. The pair is also that of `stop_codon_overlap_insertion_range`: the tuple comparison counts it under that class alone, and the format comparator sets its `Consequence` and `IMPACT` aside under that class and its residues, notation and offset under this one. |

The plugin tuple scorer (`tools/vep-plugin-pairs`) takes a row by the format comparator's test of the `Codons` cell, the residues and the terms. On such a row it sets aside an `NMD` or LoFTEE field that vep-rs alone writes, where Ensembl VEP's terms keep the plugin out, under its class `plugin_value_follows_stop_codon_term_pair`. A plugin field both engines write stays a difference.

Each of the three takes a row of section 4 under this class as well. There Ensembl VEP's row holds `stop_retained_variant` where a row of sections 1 and 2 holds `frameshift_variant` or `coding_sequence_variant`. The two engines' rows of the record both show that Ensembl VEP's stop codon test is not made, by an allele with a letter other than A, C, G and T. vep-rs's row holds `stop_gained` with `frameshift_variant`, with one in-frame term or alone. The class `stop_retained_without_stop_codon_test` takes no row on which vep-rs writes `stop_gained`. It takes the rows on which vep-rs writes `coding_sequence_variant`, `frameshift_variant` or nothing in that term's place.

## 2. Two more forms: a frameshift and an allele the 3' rule moves

**In short.** Ten inserted bases with the same stop codon are a `frameshift_variant` in both engines, and vep-rs adds `stop_gained`. An insertion that begins with the four reference bases after its point is described four bases 3', where its first two bases complete a stop codon.

```
21	10567700	stop_then_unknown_frameshift	C	CGAATAGNNNN	.	.	.
21	10567700	stop_then_unknown_shifted	C	CATTTAGNNNNCA	.	.	.
```

The rows of ENST00000622113 under `--hgvs`:

| Record | Engine | Consequence | HGVSp |
| --- | --- | --- | --- |
| frameshift | Ensembl VEP | `frameshift_variant` | `-` |
| frameshift | vep-rs | `stop_gained,frameshift_variant` | `p.Ile142GlufsTer2` |
| shifted | Ensembl VEP | `coding_sequence_variant` | `-` |
| shifted | vep-rs | `stop_gained,inframe_insertion` | `p.Ile142_Leu143insTer` |

vep-rs writes `-/E*XX` and `-/I*XX` in `Amino_acids`, and `HGVS_OFFSET=4` on the second record's row; Ensembl VEP writes no offset beside its empty notation. On the first record the second `X` stands for the base left after three whole codons, as on any frameshift row. `Consequence` and `Amino_acids` are read where the record writes the allele, the HGVSp at its most 3' placement. Ensembl VEP writes vep-rs's terms, notation and offset for each completion of either record.

The second record's bases `ATTT` equal 21:10567701-10567704, so the same sequence results from inserting the rotated allele after 21:10567704. That point lies inside the codon of Leu143, after its first base `T`. The first two inserted bases there are `AG`, and `TAG` is a stop codon: Ile142 is kept and the stop follows it. The shift stops where a known base of the allele differs from the reference, so the placement does not depend on the unknown bases either.

### 2.1 An allele the record writes in an intron

**In short.** An insertion written three bases before a coding exon opens with the four reference bases after its point, so the 3' rule places it inside the exon. There its known bases read a glutamate and a stop codon before its unknown bases. Ensembl VEP 116.3 writes no notation. vep-rs writes the HGVSp of that placement with its `HGVS_OFFSET`, beside the terms of the intron.

```
21	10542391	intron_stop_then_unknown_inframe	C	CTAGTGAATAGNNNNC	.	.	.
21	10542391	intron_stop_then_unknown_frameshift	C	CTAGTGAATAGNNNN	.	.	.
21	10542391	intron_stop_then_completed_frameshift	C	CTAGTGAATAGAAAA	.	.	.
21	10542391	intron_unknown_before_any_stop	C	CTAGTGAANNNNTAGC	.	.	.
21	14166454	reverse_intron_stop_then_unknown	G	GGNNNNCTATTCCCCTG	.	.	.
21	14166454	reverse_intron_stop_then_completed	G	GGTTTTCTATTCCCCTG	.	.	.
```

| The records | Where they stand, on the cache of section 1 |
| --- | --- |
| Forward strand | 21:10542395 is the first base of an exon of ENST00000622113 and the last base of codon 22. The records write their bases after 21:10542391, in the intron before it. Their first four bases, `TAGT`, equal 21:10542392-10542395, so the same sequence results from inserting the rotated allele after 21:10542395, between codons 22 and 23. |
| Reverse strand | 21:14166451 is the first base of an exon of ENST00000681601 on its strand and the second base of codon 215, a glycine. The records write their bases, reverse-complemented, between 21:14166454 and 21:14166455, in the intron before that exon on the transcript. As the transcript reads them, their first four bases equal the four bases after the point, so the same reading places them after 21:14166451, inside codon 215. |
| What the transcript reads at that placement | `GAA`, `TAG`, then four unknown bases in the first, second and fifth records, the reverse-strand one behind a `G` that completes codon 215 as the glycine the reference reads there. The fourth record reads `GAA` and then the unknown bases; the `TAG` after them stands out of frame. The third and sixth are completions, each `N` an `A` on the transcript. |
| What every row carries | `splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant` in both engines, and no cDNA, CDS or protein position, no `Amino_acids` and no `Codons`. Both engines write an HGVSc for the completed records alone, and an offset only beside a notation, so the fourth record's row has none. |

The rows under `--hgvs`:

| Record | Engine | HGVSp | `HGVS_OFFSET` |
| --- | --- | --- | --- |
| first, whole codons | Ensembl VEP | `-` | `-` |
| first, whole codons | vep-rs | `p.Ser22_Pro23insGluTer` | 4 |
| second, a frameshift | Ensembl VEP | `-` | `-` |
| second, a frameshift | vep-rs | `p.Pro23GlufsTer2` | 4 |
| third, the second completed | Ensembl VEP | `-` | 4 |
| third, the second completed | vep-rs | `p.Pro23GlufsTer2` | 4 |
| fourth, no stop before the unknown bases | both | `-` | `-` |
| fifth, reverse strand | Ensembl VEP | `-` | `-` |
| fifth, reverse strand | vep-rs | `p.Leu216GlufsTer2` | -4 |
| sixth, the fifth completed | Ensembl VEP | `-` | -4 |
| sixth, the fifth completed | vep-rs | `p.Leu216GlufsTer2` | -4 |

| Why Ensembl VEP writes that | Lines of `TranscriptVariationAllele.pm` |
| --- | --- |
| An allele with a letter other than A, C, G and T has no notation; the test stands before the shift. | 1683-1684 |
| A completed record is shifted, which gives its HGVSc and its offset. `hgvs_protein` then asks whether the allele is coding on the span as the record writes it, in the intron, and returns: the defect of [The HGVSp of an insertion or deletion shifted from an intron into the coding sequence](hgvsp-shifted-indel-written-span-gate.md). | 1731 |

**Why vep-rs's row is the correct one.** The 3' rule describes the insertion inside the exon, where `GAA` and `TAG` follow the last unchanged codon. The protein ends at that stop whatever the unknown bases are, and the shift stops at a known base, so every completion is one change with one notation and offset, those vep-rs writes. In the fourth record a completion can spell a stop codon among the unknown bases, which follow `GAA` directly, so its completions share no notation and neither engine writes one.

## 3. A stop codon the insertion keeps

**In short.** Seven or nine bases inserted between a transcript's last sense codon and its stop codon open with `TAG` and go on with four unknown bases. The protein does not change. Ensembl VEP 116.3 writes `stop_retained_variant` alone on a reverse-strand transcript, and `frameshift_variant` or `coding_sequence_variant` on a forward-strand one. vep-rs writes the row it writes for every completion, `inframe_insertion,stop_retained_variant`.

```
21	10605549	kept_stop_then_unknown_frameshift	T	TTAGNNNN	.	.	.
21	10605549	kept_stop_then_unknown_inframe	T	TTAGNNNNCA	.	.	.
21	14108995	reverse_kept_stop_then_unknown_frameshift	A	ANNNNCTA	.	.	.
21	14108995	reverse_kept_stop_then_unknown_inframe	A	ATGNNNNCTA	.	.	.
```

| The records | Where they stand, on the cache of section 1 |
| --- | --- |
| Forward strand | 21:10605549 is the last base of codon 533, the last sense codon of ENST00000622113, and its stop codon `TAA` is 21:10605550-10605552. The records insert `TAG` and four `N` after that base, the second with `CA` after them. |
| Reverse strand | ENST00000681601 reads its stop codon `TAA` from 21:14108995 down to 21:14108993. The records insert the reverse complement of the same bases after 21:14108995, so the transcript reads them between codon 460, its last sense codon, and its stop codon. |

The rows of the two transcripts under `--hgvs`:

| Records | Engine | Consequence | HGVSp |
| --- | --- | --- | --- |
| forward, seven bases | Ensembl VEP | `frameshift_variant` | `-` |
| forward, nine bases | Ensembl VEP | `coding_sequence_variant` | `-` |
| reverse, both | Ensembl VEP | `stop_retained_variant` | `-` |
| forward, both | vep-rs | `inframe_insertion,stop_retained_variant` | `p.Ter534=` |
| reverse, both | vep-rs | `inframe_insertion,stop_retained_variant` | `p.Ter461=` |

| Cell | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `Amino_acids` | none | `-/*XX` |
| `IMPACT` | `HIGH` beside `frameshift_variant`, `MODIFIER` beside `coding_sequence_variant`, `LOW` beside `stop_retained_variant` | `MODERATE`, the impact both engines give `inframe_insertion`; the protein itself does not change |
| `HGVS_OFFSET` | none | 2 on the forward strand and -2 on the reverse, beside the HGVSp, the offsets Ensembl VEP writes for the completions; the output writes the notation's `=` as `%3D` |

| Why Ensembl VEP writes that | Lines of `VariationEffect.pm` |
| --- | --- |
| The allele has no peptide, so `inframe_insertion` returns and `stop_retained` takes its second arm, an overlap with the stop codon and an edit test of its bases. | 1115, 1318 |
| Reverse strand: the overlap holds for an allele of `A`, `C`, `G`, `T` and `N`, and the edit test returns before it reads a base. The stop codon counts as kept, untested. | 1408-1415, 1498 |
| Forward strand: the overlap fails, the defect of [An insertion that keeps the stop codon on a forward-strand transcript](stop-codon-overlap-insertion-range.md). The row follows the length: `frameshift_variant` for seven bases, the term of an unknown change for nine. | 1408-1415, 1564, 1616-1649 |
| A completion, each `N` an `A`, has a peptide, `*KX` for seven bases and `*KT` for nine. Ensembl VEP writes `inframe_insertion,stop_retained_variant` and those residues for nine bases on both strands and for seven on the reverse strand: the peptide is cut after its first stop, and `frameshift` returns where the stop codon is kept, before it reads a length. | 1126, 1549 |
| The seven-base completion on the forward strand is the overlap entry's row, `stop_gained,frameshift_variant`: its peptide ends in the `X` of a partial codon, so `stop_retained` takes its second arm again, whose overlap fails. Where Ensembl VEP calls a completion `inframe_insertion`, its notation is `p.Ter534delinsTer` or `p.Ter461delinsTer` and vep-rs's is the stop unchanged, the difference of [The HGVSp of a change at the stop codon](hgvsp-stop-codon-protein-comparison.md). | 1313-1318 |

**Why vep-rs's row is the correct one.** `TAG` stands where `TAA` stood, behind the same residues, so the protein is the reference protein under every completion. vep-rs's stop codon test reads those three bases and the codons before them, and gives the row `stop_retained_variant` without a peptide. Where the kept stop codon stands before the first codon with an unknown base and the two proteins, compared, are the same, vep-rs reads the peptide too. The row then takes `inframe_insertion` and the residues, as vep-rs's row of each completion does.

| A kept stop codon | What vep-rs writes |
| --- | --- |
| An insertion inside the stop codon, or inside a codon before it, where every codon up to the old stop position keeps its residue | The terms, the residues and the `IMPACT` of the rows above. The genome-wide ClinVar rows of this shape are such insertions, all on forward-strand transcripts: Ensembl VEP writes `frameshift_variant` on 8 of the 12 and `coding_sequence_variant` on 4. Their most 3' placement lies past the stop codon, so neither they nor their completions have an HGVSp, and the row writes no `HGVS_OFFSET`, which a completion carries beside its HGVSc. |
| A kept stop codon behind a codon that holds an unknown base | `stop_retained_variant` alone, also where every reading of that codon is the reference's residue: the peptide is read only where no codon before the stop holds an unknown base. |
| A kept stop codon where vep-rs does not compare the two proteins: an insertion at an exon's edge, or a transcript whose protein is not the translation of a coding sequence that one stop codon closes, as with a selenocysteine | `stop_retained_variant` alone. On such a transcript the completions share no row: Ensembl VEP's `stop_retained` compares the two peptides where the alternate peptide holds no `X` and reads the edited sequence where it holds one (line 1313), and the two readings can differ. At an exon's edge the row stays behind its completions, which share `inframe_insertion` and a notation in vep-rs; a gained stop is read there (section 5). |
| A kept stop codon whose peptide is the stop alone: one or two bases inserted inside the stop codon | The row of the stop codon test, without residues. No completion has `inframe_insertion`, which line 1121 withholds where both peptides are the stop alone. |
| A kept stop codon that an intron splits | The terms and the residues, and no HGVSp: the completions do not share one there. |
| A span that replaces reference bases and keeps the stop codon | The row of the stop codon test (section 5). |

## 4. A row Ensembl VEP calls `stop_retained_variant`

**In short.** On a reverse-strand transcript Ensembl VEP 116.3 reads an insertion that stands before the stop codon on the transcript, by fewer bases than it inserts, as one that overlaps the codon. For an allele with an unknown base it then writes `stop_retained_variant` without testing the codon. Where the allele's known bases settle a stop codon, vep-rs writes the gained stop, and the row is a row of this class.

```
21	30165951	stop_then_unknown_before_the_stop_codon	A	ANNNNGTGCTATTC	.	.	.
21	30165951	stop_then_completed_before_the_stop_codon	A	ATTTTGTGCTATTC	.	.	.
```

ENST00000286808 is a one-exon transcript of CLDN17 on the reverse strand, in the cache of section 1. Its coding sequence has 675 bases: 224 residues and the stop codon, 3 × 224 + 3 = 675. The stop codon is 21:30165943-30165945. Both records insert thirteen bases between 21:30165951 and 21:30165952, after the last base of codon 222, six coding bases before the stop codon. On the transcript the thirteen bases read `GAA`, `TAG`, `CAC` and four more, unknown in the first record and `A` in the second.

The rows of ENST00000286808 under `--hgvs`:

| Record | Engine | Consequence | HGVSp |
| --- | --- | --- | --- |
| unknown | Ensembl VEP | `stop_retained_variant` | `-` |
| unknown | vep-rs | `stop_gained,frameshift_variant` | `p.Tyr223GlufsTer2` |
| completed | Ensembl VEP | `stop_gained,frameshift_variant` | `p.Tyr223GlufsTer2` |
| completed | vep-rs | `stop_gained,frameshift_variant` | `p.Tyr223GlufsTer2` |

`IMPACT` is `LOW` beside `stop_retained_variant` and `HIGH` on the other three rows. Ensembl VEP writes no `Amino_acids` on the first record's row, where vep-rs writes `-/E*HXX`; both write `-/E*HKX` on the second's. On each record the two engines print the same `Codons` cell.

**Why Ensembl VEP writes that.** The first record's allele has no peptide (`TranscriptVariationAllele.pm`, line 701). `stop_retained` then asks two routines (`VariationEffect.pm`, line 1318). Ensembl VEP holds an insertion with its start one above its end, here 30165952 and 30165951. `_overlaps_stop_codon_cil` moves that start back by the inserted length on a reverse-strand transcript (line 1413): 30165952 - 13 = 30165939. The span 30165939 to 30165951 holds the stop codon, so the overlap holds.

`_ins_del_stop_altered_cil` returns 0 at line 1498 for an allele with a letter other than A, C, G and T, before it reads a base. Line 1318 reads that 0 as a stop codon kept. `frameshift` and `coding_unknown` return 0 beside `stop_retained` (lines 1549 and 1635), so the row holds that term alone.

Each of the two routines has an entry of its own. [`stop_retained_variant` on an insertion or deletion whose stop codon test is not made](stop-retained-without-stop-codon-test.md) describes the test. [An insertion that keeps the stop codon on a forward-strand transcript](stop-codon-overlap-insertion-range.md) describes the moved start, on both strands.

**Why vep-rs's row is the correct one.** The thirteen bases stand between codons 222 and 223. `GAA` is residue 223 and `TAG` ends the protein after it, where the reference protein has 224 residues. A stop codon is gained one residue early, whatever the four unknown bases are. The three bases at the old stop position, coding bases 673 to 675 of the changed sequence, are the seventh to ninth inserted bases, `CAC`. They are no stop codon, so no stop codon is kept there. Ensembl VEP writes `stop_gained,frameshift_variant` itself for the second record, whose bases it can test.

**How often.** The genome-wide ClinVar files hold 36 such rows on GRCh38 and 8 on GRCh37, every one on a reverse-strand transcript. vep-rs writes `stop_gained,frameshift_variant` on 22 and 7 of them, `stop_gained,inframe_insertion` on 13 and 1, and `stop_gained` alone on 1 row of GRCh38, an insertion of whole codons that stands inside a codon. vep-rs gives that row no in-frame term when its unknown bases are known either. Section 1.7 says how the comparisons take them.

## 5. What the rule does not read

- A stop codon that stands only after the unknown bases. Its frame rests on the written count of unknown bases, and a completion can end the protein earlier.
- A stop codon that itself holds an ambiguity letter, even one every reading of which is a stop codon. The protein ends there under every reading too, but the rule reads a stop codon of known bases alone, and such a row stays as Ensembl VEP writes it.
- An allele whose reference codons hold a stop codon: no stop is gained there. An insertion that keeps that stop codon is section 3's.
- An allele with an unknown base in its last codons, where reference residues are compared with them. A deletion-insertion, a substitution of several bases or an insertion inside a codon replaces reference residues. Where those residues do not begin the alternate peptide, Ensembl VEP's rules compare them with its end: `protein_altering_variant` asks whether the peptide ends with them, and the notation leaves out last residues the two share. Whether an unknown codon there equals a reference residue rests on the unknown bases. Such a row stays as Ensembl VEP writes it. An insertion between two codons replaces no residue.
- A stop codon that restates the protein's own end. A deletion-insertion can keep the residues before it, spell the protein's last residues and a stop codon again, and hold its unknown bases after them. The protein is then the reference protein under every completion, which Ensembl VEP's comparison of the two protein sequences finds (`VariationEffect.pm`, lines 1326-1361). No stop is gained, and the row stays as Ensembl VEP writes it.
- A record whose alleles share end bases up to an unknown base. Both engines cut the bases a record's two alleles share at their ends and annotate what is left. Where that cut stops at an unknown base of the alternate allele beside a reference base, the completion that spells the reference base there is cut further and read at another span, so the completions share no row. A record that writes an insertion after one shared base is not such a record: the cut uses up its reference allele.
- A row to which a test of the start codon or of the stop codon gives its term: `start_lost`, `start_retained_variant` or `stop_lost`, and `stop_retained_variant` on a span that replaces reference bases. Those tests decide such a row without a peptide, vep-rs keeps the row they give, and their own entries say what each engine writes there; see [`stop_retained_variant` on an insertion or deletion whose stop codon test is not made](stop-retained-without-stop-codon-test.md). The tests meant are vep-rs's: a row on which Ensembl VEP alone writes `stop_retained_variant`, for a test it does not make, is the row of section 4.
- An HGVSp for an allele that is not inside one exon at the placement its notation is read at: an insertion between an exon's first or last base and the intron beside it, or bases that cross an exon's end. Ensembl VEP's notation there turns on the inserted bases beside the edge, so the completions do not share one. The same holds on a record of two or more alternate alleles whose alleles share an end base with the reference allele: the engines leave those bases uncut and read the notation with them clipped.

  Nor is one written where the most 3' placement settles no stop codon. Such a row keeps its terms and residues and has no notation.
- An HGVSp beside an unknown residue of the reference protein, the `X` of a coding sequence that opens or ends inside a codon, on a row that does not shift the frame. The notation writes its own stop as `X` and tests whether the inserted residues repeat the residues beside them. An `X` there can match, and one completion then reads as a duplication where the others read as an insertion. Such a row keeps its terms and residues too.
- The HGVSc. The completions share its positions and not its inserted bases, so no exact description exists. The HGVS recommendations permit an incomplete one and do not require it.
