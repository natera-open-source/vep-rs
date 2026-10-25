# Defect C: a record marked skipped is dropped before annotation

**The difference in one sentence.** Ensembl VEP 116 drops a structural record it marks skipped, by setting its `vep_skip` flag, before the record reaches the annotation pass, so the record has no row in any output format, while vep-rs annotates it, an incomplete deletion excepted.

Defect C itself is corrected; what remains is not a defect but a difference of coverage: Ensembl VEP writes nothing for the record, vep-rs annotates it, and no row of either engine can be judged against the other's.

**What release 116 corrected.** Up to release 115.2 the skipped record stayed in its [batch](defect-p.md#231-what-a-reader-needs-to-know-first): it loaded no transcript regions of its own and was annotated against whatever regions its neighbours had loaded (`AnnotationSource.pm` line 238 against `InputBuffer.pm` lines 284-330 of that release). Its output changed with `--fork`, `--buffer_size` and the order of lines in the file while the warnings file reported the record as skipped; one record got 51 transcripts, or 497, or none, depending on its neighbours.

That batch dependence was defect C, and it produced the [paper's](../published-figures.md) classes 5 to 8: an oversize structural variant annotated against whichever transcripts the batch had loaded, the mate end of such a breakend, `<CPX>`, and a `<NON_REF>` block annotated as a variant. From release 116 the parser drops the record instead (`validate_svf` returns 0 for a `vep_skip` record and `next` discards it unless `--dont_skip` is set), which is the correction those classes called for. This section documents the coverage difference and the records it reaches.

**Which records are dropped.** A record takes one of four routes to the mark, or a fifth to the same drop:

- a structural variant whose span exceeds `--max_sv_size` (default 10,000,000; `Config.pm` line 315). The span is the one `ensembl-io` derives for the record (`BaseVCF4.pm` lines 297-317): from `SVLEN` first, else from `INFO/END`, else the record's start itself. The size test compares the span's end minus its start with the cap (`Parser.pm` lines 491-501).

  A breakend, one end of a rearrangement joined to its mate, has its `INFO/END` discarded before the span is derived (`Parser/VCF.pm` lines 487-489), so a breakend's span is its `SVLEN` where the record carries one. A gnomAD-SV breakend record carries an `SVLEN` to its mate, so a same-chromosome breakend whose mate lies more than 10,000,000 bases away is dropped whole, its mate end included;
- an ALT type with no Sequence Ontology term, `<CPX>`, `<CTX>` and `<NON_REF>` (a reference-confidence block of a genomic VCF) among them. `get_SO_term` (`Parser.pm` lines 663-728) finds no term for the type, and `Parser/VCF.pm` lines 477-481 mark the line, the type string standing in for the term;
- a multi-allelic symbolic list with no type of its own (`<INS>,<DEL>`, `<INV>,<DUP>`): `get_SO_term` reads one type from the ALTs joined by `/` and finds none;
- a deletion of any class (`<DEL>`, `<DEL:ME:ALU>`, or an unsupported type whose string carries `del`, `<FOO:DEL>` among them) with neither `END` nor `SVLEN`, or whose start exceeds its end (`Parser/VCF.pm` lines 555-558, "deletion looks incomplete"). The test reads the term or, for an unsupported type, the type string of lines 469-481;
- a record whose start or end is not a number, or whose end lies more than one base before its start (`END=0`, `END=abc`, a breakpoint's `END=`). `validate_vf`'s own checks return 0 before any mark is read (`Parser.pm` lines 475-488), and the record is dropped the same way.

For an oversize record `validate_vf` itself sets the mark and returns 0 (`Parser.pm` lines 498-499). On the other three marked routes the mark is set on the feature at `Parser/VCF.pm` line 575, and `validate_vf` hands the structural feature to `validate_svf` (`Parser.pm` line 551), which returns 0 for a marked record (lines 821-827). `next` at lines 189-191 discards the record unless `--dont_skip` is set. The warning is still written to the warnings file, naming the line and the reason.

**What vep-rs does instead.** vep-rs annotates every record of the first three marked routes, drops the fourth route's records as Ensembl VEP drops them, and refuses the fifth route's records at warning level as Ensembl VEP does, so neither engine writes a row for those. In the default, tab and Parquet outputs an oversize or unsupported record is annotated against every transcript it overlaps, and a list with no type gets one row per ALT under each ALT's own class.

An unsupported ALT whose letters name a class is written under that class's term: `deletion` for `<del>`, `<Del>`, `<DEL:>` or `<DEL >`, `duplication` for `<Dup>`, `inversion` for `<Inv>`, `copy_number_variation` for `<CN=X>`, `insertion` for `<Ins>`, `chromosome_breakpoint` for `<Bnd>`. One whose letters name none is written under the type string Ensembl VEP keeps for it: the ALT as written (`<CPX>`), or `INFO/SVTYPE` when that names the type (`CPX` for `<CPX>` with `SVTYPE=CPX`, the gnomAD-SV spelling). A `<NON_REF>` block, which asserts the absence of a variant, gets one `intergenic_variant` row, and its skip warning, in Ensembl VEP's form, goes to the warnings file and to standard error unless `--quiet`.

In the VCF output the line of an oversize or unsupported record is carried without consequences; in the JSON output an unsupported record is an `input`-only object and an oversize one is omitted.

A `<DEL>` with neither `END` nor `SVLEN`, and one whose start exceeds its end (`END` at POS, `SVLEN=0`, a `SVLEN` that is not a number), is warned `deletion looks incomplete` and written in no format, as Ensembl VEP writes it in none. `--dont_skip` keeps it on both engines, read as a deletion of the base after POS or as the pair `POS-(POS+1)` with an insertion's positions (`21:30425650-30425651 deletion` for `21 30425650 x N <DEL> . . SVTYPE=DEL;END=30425650`). Its VCF line is carried without consequences and its JSON object is `input`-only.

A `SVLEN` of `NaN` is the exception: Ensembl VEP's end is `NaN`, which no `start > end` test catches, so it writes the record's intergenic row at `21:30425651-NaN` ([Defect K](defect-k.md)); vep-rs reads the length as 0 and drops the record as incomplete.

Under the flag the default and tab outputs can differ. Ensembl VEP loads no transcript region for a marked record (`AnnotationSource.pm` line 238) and annotates it against whatever regions its batch-mates loaded: one `intergenic_variant` row for a record alone in its region, and transcript rows beside a neighbour. vep-rs annotates the kept record against the transcripts it overlaps whatever the batch holds. The pair `POS-(POS+1)` then carries `feature_truncation` on Ensembl VEP's side and not on vep-rs's, as the same pair written in the Ensembl format (`21 30165851 30165850 DEL +`), which neither engine marks, does without the flag.

vep-rs takes `--max_sv_size` too, with the same default, and `--max_sv_size -1` lifts the cap on both engines.

**Why the rows are set aside.** The divergence is every row vep-rs writes for such a record, and there is no Ensembl VEP row at the record to judge it against. Ensembl VEP's output does not say what the record's transcripts are, only that the record was not annotated. The comparators therefore set vep-rs's rows aside rather than charge them, and the golden corpora record those rows in their manifests under the class `reference_skipped_record` (`tests/golden/116/GRCh37`: four `<BND>` records above `--max_sv_size` and two `<CPX>`; `tests/golden/116/GRCh38-release116`: `<INS>,<DEL>` and `<INV>,<DUP>`; `tests/golden/116/GRCh38-vcf-edge-cases`: a `<NON_REF>` block).

**The class it produces.**
- Class 5: the dropped records of the measured sets. GRCh37: 389 records and 20,524 vep-rs tuples, in gnomAD-SV v2.1 chr21 (89 records), the synthetic special-allele set (its 200 `<NON_REF>` blocks) and the synthetic complex set (its 100 `<CPX>` records). GRCh38: 594 records and 278,470 tuples per file, in gnomAD-SV v4.1 chr21 (245 records), 1000 Genomes SV (49), the special-allele set (200) and the complex set (100). Pooled, 278,309: 161 tuples of records the gnomAD-SV and 1000 Genomes files carry at identical coordinates are counted once. Set aside: vep-rs's rows.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 5. A record Ensembl VEP drops before annotation

**In short.** Ensembl VEP drops every structural record it marks skipped, writing only a warning: a span over `--max_sv_size`, an ALT type with no Sequence Ontology term (`<CPX>`, `<CTX>`, `<NON_REF>`), a typeless multi-allelic list or an incomplete `<DEL>`. vep-rs drops an incomplete `<DEL>` too, with the same warning, and under the size check a record whose `SVLEN` or `END` overflowed the coordinate range, whose span is the ceiling's and not its own (`--dont_skip` keeps it; `--max_sv_size -1` lifts the check in both engines). It annotates every other such record against every transcript its span overlaps or comes within 5,000 bases of, since a record's coordinates alone decide which transcripts it reaches, and gives a `<NON_REF>` block one `intergenic_variant` row.

The cap is a memory guard and the missing term a gap in Ensembl VEP's vocabulary, so this is a coverage difference: with no Ensembl VEP row to compare, vep-rs's rows are set aside, 389 records on GRCh37 and 594 on GRCh38.

### 5.1 What a reader needs to know first

A breakend record (`SVTYPE=BND`) describes one end of a rearrangement: the position where the reference is broken and the position it is joined to, its mate. gnomAD-SV writes the record with a symbolic `<BND>` allele and puts the mate in `INFO/CHR2` and `INFO/END2`; VCF's own notation writes the mate into the allele, `N[chr21:41770788[`. Both engines read either form. When the mate is on the same chromosome both represent the record as a feature spanning from the base after POS to the mate (Location `21:14604132-41770788`), with rows under two alleles: the local end, `N.`, and the mate, `N[chr21:41770788[`.

`--max_sv_size` (default 10,000,000) is the documented size cap of both engines. A structural variant whose span exceeds it is reported as skipped in the warnings file Ensembl VEP writes beside its output (`<output>_warnings.txt`); Ensembl VEP's source calls the cap a guard "to avoid memory problems". The span Ensembl VEP measures is derived from `SVLEN` first: for the breakend below, `SVLEN=27166657` gives a span of 27,166,656 bases, whatever `END` says.

`<CPX>` is gnomAD-SV's symbolic allele for a complex rearrangement; here the rearrangement is a deletion followed by an inversion, `CPX_TYPE=delINV`. Ensembl VEP has no Sequence Ontology class for it. vep-rs annotates the record for the span alone, under the type string `CPX` and without a direction term (`feature_elongation` or `feature_truncation`). Its rows claim overlap, not disruption: `coding_transcript_variant` (SO:0001968) and `non_coding_transcript_variant` (SO:0001619, "a transcript variant of a non coding RNA gene") for transcripts inside the span, `intron_variant` and the exon terms for transcripts the span partly overlaps, and `mature_miRNA_variant` for a microRNA it covers.

A gVCF (genomic VCF) reports every position of a sample, not only the variant ones, and compresses the positions where the sample matches the reference into reference-confidence blocks. A block is one record whose ALT is the symbolic `<NON_REF>` (or its VCF 4.3 spelling `<*>`, the same construct), spanning POS to `INFO/END`, meaning "no variant was called in this interval, with this confidence". Such a record asserts the absence of a variant; there is no allele to predict a consequence for. Ensembl VEP has no Sequence Ontology class for `<NON_REF>` and drops the block.

vep-rs writes the record's one output line in Ensembl VEP's own convention for a record that overlaps nothing it annotated: Feature `-`, consequence `intergenic_variant` (SO:0001628), IMPACT MODIFIER. It also writes the record's skip warning in Ensembl VEP's form, to the warnings file and to standard error unless `--quiet`. That row is a convention, not a claim that the block lies between genes, and it keeps one output line per input record, which is what the default output promises. `<*>` is not marked skipped by Ensembl VEP: both engines annotate it as a span and their rows agree, so it is outside the class by measurement, not by principle.

### 5.2 The records

`gnomAD-SV_v3_BND_chr21_3f3f16e7`, gnomAD-SV v4.1 chr21 (GRCh38, `gnomad_sv_v4.1_chr21.vcf.gz`), a breakend at 21:14604131 whose mate lies 27 Mb away on the same chromosome:

```
21	14604131	gnomAD-SV_v3_BND_chr21_3f3f16e7	N	<BND>	251	UNRESOLVED	CHR2=chr21;END=14604131;END2=41770788;POS2=41770787;SVLEN=27166657;SVTYPE=BND
```

`gnomAD-SV_v3_CPX_chr21_df09f6c2`, from the same file, a 1.66 Mb complex event, well below `--max_sv_size`:

```
21	14923958	gnomAD-SV_v3_CPX_chr21_df09f6c2	N	<CPX>	1	UNRESOLVED	CHR2=chr21;CPX_INTERVALS=DEL_chr21:14923958-15004545,INV_chr21:15004545-16588073;CPX_TYPE=delINV;END=16588080;SVLEN=1664122;SVTYPE=CPX
```

`synth_nonref_0033`, from the synthetic special-allele set `tests/sv_validation/grch37/10_special_alleles.vcf.gz` (GRCh37; the GRCh38 set is `tests/sv_validation/grch38/10_special_alleles.vcf.gz`), is a 279-base block whose interval 27347259-27347537 lies in APP (ENST00000346798, reverse strand, canonical, GRCh37 27252861-27542972). It covers 155 of the 159 bases of exon 11 (27347383-27347541, a coding exon) and 124 bases of intron 11 (27328070-27347382):

```
21	27347258	synth_nonref_0033	C	<NON_REF>	.	PASS	END=27347537
```

### 5.3 What each engine writes

Ensembl VEP 116.2 writes no row for any of the three records, in any output format; its warnings file names each line as skipped, with the reason. vep-rs writes 4,652 rows for the breakend: 4,628 under the local allele `N.`, for the transcripts its span 14604132-41770788 overlaps or comes within 5,000 bases of. The other 24 are under the mate allele `N[chr21:41770788[`, for the transcripts within 5,000 bases of the mate breakpoint 41770788. It writes 534 rows for the complex event, every transcript inside 14923959-16588080 under the allele `CPX`, and one `intergenic_variant` row for the block.

Ensembl VEP, warnings file, each record run alone:

```
WARNING: line 1 skipped (21 14604131 gnomAD-SV_v3_BND_chr21_3f3f16e7 N ...): variant size (27166656) is bigger than --max_sv_size (10000000)
WARNING: line 1 skipped (21 14923958 gnomAD-SV_v3_CPX_chr21_df09f6c2 N ...): CPX is not a supported structural variant type
WARNING: line 1 skipped (21 27347258 synth_nonref_0033 C <NON_REF> . PA...): NON_REF is not a supported structural variant type
```

Ensembl VEP, output:

```
(no row for any of the three records)
```

Five of vep-rs's rows follow, one table each.

vep-rs, the breakend (2 of 4,652 rows: the SAMSN1 row under the local allele and a mate-allele row). The SAMSN1 row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `gnomAD-SV_v3_BND_chr21_3f3f16e7` |
| Location | | `21:14604132-41770788` |
| Allele | | `N.` |
| Gene | | `ENSG00000155307` |
| Feature | | `ENST00000400562` |
| Feature_type | | `Transcript` |
| Consequence | | `feature_truncation,non_coding_transcript_exon_variant,intron_variant` |
| cDNA_position | | `-` |
| CDS_position | | `-` |
| Protein_position | | `-` |
| Amino_acids | | `-` |
| Codons | | `-` |
| Existing_variation | | `-` |
| Extra | | `IMPACT=HIGH;STRAND=-1` |

The breakend's mate-allele row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `gnomAD-SV_v3_BND_chr21_3f3f16e7` |
| Location | | `21:14604132-41770788` |
| Allele | | `N[chr21:41770788[` |
| Gene | | `ENSG00000289513` |
| Feature | | `ENST00000700893` |
| Feature_type | | `Transcript` |
| Consequence | | `feature_truncation,non_coding_transcript_exon_variant,intron_variant` |
| cDNA_position | | `-` |
| CDS_position | | `-` |
| Protein_position | | `-` |
| Amino_acids | | `-` |
| Codons | | `-` |
| Existing_variation | | `-` |
| Extra | | `IMPACT=HIGH;STRAND=1` |

vep-rs, the complex event (2 of 534 rows: USP25 and the microRNA MIRLET7C). The USP25 row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `gnomAD-SV_v3_CPX_chr21_df09f6c2` |
| Location | | `21:14923959-16588080` |
| Allele | | `CPX` |
| Gene | | `ENSG00000155313` |
| Feature | | `ENST00000285679` |
| Feature_type | | `Transcript` |
| Consequence | | `coding_transcript_variant` |
| cDNA_position | | `-` |
| CDS_position | | `-` |
| Protein_position | | `-` |
| Amino_acids | | `-` |
| Codons | | `-` |
| Existing_variation | | `-` |
| Extra | | `IMPACT=MODIFIER;STRAND=1;OverlapBP=150140;OverlapPC=100.00` |

The complex event's MIRLET7C row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `gnomAD-SV_v3_CPX_chr21_df09f6c2` |
| Location | | `21:14923959-16588080` |
| Allele | | `CPX` |
| Gene | | `ENSG00000199030` |
| Feature | | `ENST00000362160` |
| Feature_type | | `Transcript` |
| Consequence | | `mature_miRNA_variant` |
| cDNA_position | | `-` |
| CDS_position | | `-` |
| Protein_position | | `-` |
| Amino_acids | | `-` |
| Codons | | `-` |
| Existing_variation | | `-` |
| Extra | | `IMPACT=MODIFIER;STRAND=1;OverlapBP=84;OverlapPC=100.00` |

vep-rs, the block:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `synth_nonref_0033` |
| Location | | `21:27347259-27347537` |
| Allele | | `<NON_REF>` |
| Gene | | `-` |
| Feature | | `-` |
| Feature_type | | `-` |
| Consequence | | `intergenic_variant` |
| cDNA_position | | `-` |
| CDS_position | | `-` |
| Protein_position | | `-` |
| Amino_acids | | `-` |
| Codons | | `-` |
| Existing_variation | | `-` |
| Extra | | `IMPACT=MODIFIER` |

In the reference outputs of the whole files the same three warnings name line 9777 and line 9978 of the gnomAD-SV file and line 188 of the special-allele set, and the outputs carry no row for the records. vep-rs writes the same skip warning for the block, naming its record as Ensembl VEP's does, and then writes the intergenic row.

### 5.4 Why Ensembl VEP writes nothing

Each record is marked skipped on its route and dropped by the parser before it reaches the [input buffer](defect-p.md#231-what-a-reader-needs-to-know-first):

- The breakend: `Parser/VCF.pm` lines 487-489 discard its `INFO/END` and re-derive the end, which `BaseVCF4.pm` lines 302-306 take from `SVLEN`. The span is 27,166,656 bases, and the size test at `Parser.pm` lines 491-501 marks the record and writes `variant size (27166656) is bigger than --max_sv_size (10000000)`.
- The complex event: `get_SO_term` (`Parser.pm` lines 663-728) has no term for `CPX`, so `Parser/VCF.pm` lines 477-481 mark the line and write `CPX is not a supported structural variant type`.
- The block: the same route, `NON_REF is not a supported structural variant type`.

Then `validate_vf` returns 0 for the breakend (`Parser.pm` line 499) and `validate_svf` (lines 821-827) for the complex event and the block, and `next` at lines 189-191 discards each. Nothing of the record reaches the annotation pass or any output writer. With `--dont_skip` the marked record enters the buffer, loads no transcript regions of its own (`AnnotationSource.pm` line 238) and is annotated against whatever its batch has loaded, as releases up to 115.2 annotated it.

### 5.5 Why vep-rs annotates the record

The record's coordinates alone determine which transcripts it reaches. The breakend's span 14604132-41770788 overlaps, or comes within 5,000 bases of, a definite set of chromosome 21 transcripts, the same from any file. The complex event's span 14923959-16588080 contains its transcripts as surely as any deletion contains its own. A `<NON_REF>` block changes no base, so no transcript's sequence, splicing or dosage is affected by it, which is what the intergenic row says.

Ensembl VEP's cap is a memory guard, not a statement about the record, and its missing term for `<CPX>` is a gap in its vocabulary, not a reading of the event. Skipping is a legitimate output for both, and vep-rs gives that skip in its VCF and JSON outputs, where it carries the line without consequences or omits it. On the default, tab and Parquet outputs it writes the complete set, which is what a caller who kept the record in the file asked for. Neither engine's output here can be checked against the other's, which is why the rows are set aside rather than defended as the correct ones.

### 5.6 Try it yourself

Run Ensembl VEP 116.2 and vep-rs (GRCh38 cache and FASTA for the gnomAD-SV records, GRCh37 for the block, the commands of [the closing section](reproducing.md)) on inputs built from the records:

| Input | Ensembl VEP 116.2 rows for the record | vep-rs rows for the record |
| --- | --- | --- |
| the breakend alone | 0, the size warning | 4,652 (4,628 local, 24 mate) |
| the breakend plus the SNV `21 14604131 . A C` | 0 (the SNV gets its 4 rows) | 4,652 |
| the breakend plus the SNV `21 41770788 . C A`, at the mate breakpoint | 0 (the SNV gets its 24 rows) | 4,652 |
| the whole `gnomad_sv_v4.1_chr21.vcf.gz`, `--fork 4 --buffer_size 5000` | 0, line 9777 skipped | 4,652 |
| the complex event alone | 0, the type warning | 534 |
| the complex event plus the four SNVs `21 14923958 . T A`, `21 15500000 . G A`, `21 16400000 . G A` and `21 16590000 . A C` | 0 (the SNVs get their 368 rows) | 534 |
| the block alone | 0, the type warning | 1 |
| `synth_nonref_0188` (`21 38589780 synth_nonref_0188 A <NON_REF>` with `END=38590008`, a 228-base block in a DSCR9 intron) plus the SNV `21 38589780 . A C` | 0 (the SNV gets its 8 rows) | 1 |

No neighbouring record changes the outcome: under release 115.2 the breakend beside the SNV at its mate breakpoint got 20 mate-allele rows and the block beside its SNV eight transcript rows, and under 116.2 both get none.

Other records of the measured sets on the same routes: the GRCh37 breakend `gnomAD-SV_v2.1_BND_21_56677` (`SVLEN=14688819`; the warning reads `variant size (14688818)`) gets 0 rows against vep-rs's 381, and the GRCh37 complex event `gnomAD-SV_v2.1_CPX_21_5151` 0 against 4. The GRCh38 `gnomAD-SV_v3_CTX_chr21_d895c6b5` gets 0 (`CTX is not a supported structural variant type`) against vep-rs's one `intergenic_variant` row under the allele `CTX`. The GRCh38 breakend `gnomAD-SV_v3_BND_chr21_f6d2016e` (`SVLEN=39158730`; the warning reads `variant size (39158729)`) gets 0 against 6,967. `gnomAD-SV_v3_CPX_chr21_896ebdf0`, a `<CPX>` spanning 39,354,624 bases, is on two routes at once and gets both warnings, the type's and the size's, against vep-rs's 7,123 rows.

### 5.7 How the comparator treats this class

The two comparators recognise a dropped record from what each reads, and both set vep-rs's rows for it aside.

The structural-variant comparator, `filter_reference_skipped_records` in `scripts/validation/compare_sv_concordance.py`, reads the input VCF and names the records from it with a port of Ensembl VEP's own routes. The routes are `get_SO_term` over the ALT list and `SVTYPE`, the incomplete-deletion test, the size test over the span `get_end` derives (`SVLEN` first, else `END`, a breakend's `END` discarded), and the coordinate checks of `validate_vf`.

Each such record is confirmed dropped when Ensembl VEP wrote no row under the record's own Location start and an allele string the record produces. Those strings are the Sequence Ontology class of its ALT and of its ALT list, and the class term vep-rs prints for an ALT whose letters name a class but whose spelling has no term. They include the type string of an ALT whose letters name none, the ALT as written (`<CPX>`) or `INFO/SVTYPE` when that names it (`CPX`), a breakend's bracket, own-end and `CHR2`/`END2` mate strings, and a sequence ALT as written.

For the breakend of the records above those strings are its own-end allele `N.` and its mate allele, and for the complex event, which carries `SVTYPE=CPX`, the type string `CPX`. Exactly that record's vep-rs rows, under those alleles, are set aside.

A deletion whose start exceeds its end, when `--dont_skip` keeps it, is written by vep-rs at the pair `POS-(POS+1)`, printed low to high and starting one base before the record's own start, and the comparator locates its rows there. A kept record that shares the start, or the whole Location string, with a dropped one neither hides it nor joins it. gnomAD-SV places a `<CPX>` at the Location of a kept `<DEL>`, and a `<NON_REF>` block can start where a kept record starts. Two records at one start producing the same allele are told apart by the end of the Location.

The records and the outputs are the authority, not the warnings file, which is absent on a cached reference and truncates its record excerpt. A disagreement between them stops the comparison with an `ERROR:` naming the record:

- a record on a route whose own rows Ensembl VEP wrote (a reference run under `--dont_skip` or under another `--max_sv_size`, declared with `--reference-max-sv-size`, `-1` for no cap);
- vep-rs structural rows at a Location Ensembl VEP wrote nothing for that no record on a route explains;
- or a dropped and a kept record at one start whose alleles and ends coincide.

The SNP/indel comparator, `scripts/concordance/compare_vep_outputs.py`, opens the input only to recognise gVCF variant sites under `--input-dir` or `--input-vcf`, so for this class it applies an output-only rule. Ensembl VEP 116.2 writes at least one row for every record it keeps, an `intergenic_variant` row when no transcript is in range. So a Location at which vep-rs wrote rows and Ensembl VEP wrote nothing, under any allele, is a record Ensembl VEP dropped.

At least one vep-rs row at the Location must carry a structural allele (every route to the mark is structural), and `tandem_repeat` does not admit a Location, because Ensembl VEP keeps a `<CNV:TR>` record and writes it at another Location (class 12). Once admitted, every vep-rs row at the Location is set aside, a sequence allele beside the symbolic one of a mixed ALT list included, because Ensembl VEP skips the line whole. The SNP/indel datasets carry no such record.

Counts (`excluded_reference_skipped_rust`, per file in `excluded_reference_skipped_by_file`, the records in `reference_skipped_records`):

- GRCh37: 20,524 vep-rs tuples on 389 records. gnomAD-SV v2.1 chr21: 19,135 tuples on 89 records, its 60 `<CPX>` records and 29 `<BND>` records whose `SVLEN` exceeds the cap. `10_special_alleles`: 200 tuples on its 200 `<NON_REF>` blocks, one `intergenic_variant` row each. `12_complex_imprecise`: 1,189 tuples on its 100 `<CPX>` records.
- GRCh38: 278,470 vep-rs tuples on 594 records per file; pooled, 278,309, because 161 tuples of records the gnomAD-SV and 1000 Genomes files carry at identical coordinates are counted once. gnomAD-SV v4.1 chr21: 272,627 tuples on 245 records, its 169 `<CPX>` records, its two `<CTX>` records, and 72 `<BND>` records and two mobile-element deletions whose `SVLEN` exceeds the cap. `1kg_sv_chr21`: 662 tuples on its 49 `<CPX>` records. `10_special_alleles`: 200 tuples on 200 blocks. `12_complex_imprecise`: 4,981 tuples on its 100 `<CPX>` records.

The class is fully masked on the vep-rs side. No record of the measured sets takes the fifth route.
