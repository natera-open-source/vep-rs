# Defect O: a `<CNV:TR>` allele with `RN=0` drops a real allele of the record

**The defect in one sentence.** With a FASTA, Ensembl VEP joins every allele of a multi-allelic `<CNV:TR>` record into one allele string, and an allele whose `RN` is 0 becomes an empty part of it, so one of the record's alleles is lost.

Each allele is expanded to literal sequence before the join, and the reference comes first in the string. An empty middle part hides the multi-allelic shape, so the record is minimised as a single-allele indel and only the empty allele survives. An empty last part is dropped by `split`, so only the first allele remains.

**Why this is a defect.** The record states two alleles, and both are annotated or neither. Which allele is lost depends on the order the ALTs are listed in, not on what they say, and the lost allele vanishes without a row or a warning.

**What vep-rs does instead.** With a FASTA vep-rs expands the alleles as Ensembl VEP expands them and reads the joined string as a multi-allelic one. An empty middle part is the `RN=0` allele, the deletion of the whole run that Ensembl VEP writes alone, and the real allele keeps its own literal reading beside it. An empty last part is dropped by `split` in both engines, so a trailing `RN=0` allele is lost in both. Without a FASTA both engines write the symbolic rows for both alleles.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):
- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`, line 408: walks `RN` per allele. Line 411: for `RN=0` the loop `for ( 1 .. $records )` runs no iteration, so `$repeat` stays undefined and is pushed as the allele at line 424. Line 428: the reference and the alleles are joined with `/`, `AGCAGCAGC//AGCAGCAGCAGC` for `RN=0,1` and `AGCAGCAGC/AGCAGCAGCAGC/` for `RN=1,0`.
- `Parser.pm`, line 978: a multi-allelic string is recognised by `/.+\/.+\/.+/`, which neither string matches, an empty part having no character, so the record is not split per allele.
- `Parser.pm`, line 881: the string is minimised as one indel, where Perl's `split` on `/` drops a trailing empty part and keeps an empty middle one, and the first minimised feature is kept: `AGCAGCAGC/-` for `RN=0,1`, `-/AGC` for `RN=1,0` and `RN=2,0`.

**What correcting it would change.** Both alleles would be annotated whatever their order; the class would vanish for the records where vep-rs keeps the real allele, and both engines would gain the trailing `RN=0` allele.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 22. A real `<CNV:TR>` allele dropped beside an `RN=0` allele

**In short.** A multi-allelic `<CNV:TR>` record lists `RN` per allele, and `RN=0` is a real allele, the loss of the whole run. With a FASTA, Ensembl VEP expands the alleles into one allele string, and the `RN=0` allele becomes an empty part of it: an empty middle part hides the multi-allelic shape, an empty last part is dropped by `split`. Either way one allele is annotated and the other vanishes without a row or a warning; vep-rs annotates both when the `RN=0` allele leads, while a trailing `RN=0` allele is lost to `split` in both engines.

### 22.1 What a reader needs to know first

A multi-allelic `<CNV:TR>` record gives `RN` per ALT, the number of repeat sequences each alternate allele consists of ([12.1](cnv-tr.md#121-what-a-reader-needs-to-know-first)); `RN=0` says the allele has no repeat sequence, the run being absent from it. With a FASTA, Ensembl VEP expands every allele to literal sequence ([12.4](cnv-tr.md#124-why-ensembl-vep-writes-that)) and joins the alleles into one allele string, reference first; an allele of no repeat sequence is an empty part of that string.

When the empty part is in the middle (`RN=0,1`), the string `AGCAGCAGC//AGCAGCAGCAGC` is not recognised as multi-allelic, the record is minimised as a single-allele indel and only the first alternate, the empty one, survives: the real gain of one unit is gone. When the empty part is last (`RN=1,0`, `RN=2,0`), the split drops it and only the first allele remains: the loss of the run is gone. Either way one allele of the record is annotated and the other vanishes without a row or a warning. Behind a leading `RN=0` vep-rs annotates both, that allele as the loss of the whole run and the other as its own gain or loss; both engines lose a trailing `RN=0`.

### 22.2 The record

```
21	30425587	rn_0_1	N	<CNV:TR>,<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=4
21	30425587	rn_1_0	N	<CNV:TR>,<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=1,0;RUS=AGC;RUC=4
21	30425587	rn_2_0	N	<CNV:TR>,<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=2,0;RUS=AGC,AGC;RUC=2,2
21	30425587	rn_0_1_ruc3	N	<CNV:TR>,<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=3
```

The run at 30425588-30425596 is three copies of `AGC` in the GRCh38 sequence ([19.2](defect-l.md#192-the-record)), so the fourth record's real allele, three units, is the run unchanged. Executed with `--fasta` (Ensembl's chromosome 21 of release 116) on Ensembl VEP 116.2.

### 22.3 What each engine writes

For the first record Ensembl VEP 116.2 with `--fasta` writes two rows under the empty allele alone (JSON `"allele_string":"AGCAGCAGC/-"`), and the four-unit allele is gone. For the fourth record it writes the same two rows under its own ID, and the unchanged allele is gone.

The first record, `rn_0_1`, on ENST00000334068: the first of Ensembl VEP's two rows; vep-rs's rows on this transcript are not shown.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `rn_0_1` | not shown |
| Location | `21:30425588-30425596` |  |
| Allele | `-` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4624;STRAND=1` |  |

For the second and third records it writes two rows under the gain alone (JSON `"allele_string":"-/AGC"`), and the `RN=0` allele is gone.

The second record, `rn_1_0`, on ENST00000334068: the first of Ensembl VEP's two rows; vep-rs's rows for this record are not shown.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `rn_1_0` | not shown |
| Location | `21:30425596-30425597` |  |
| Allele | `AGC` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |
| Extra | `ALLELE_NUM=1;IMPACT=MODIFIER;DISTANCE=4623;STRAND=1` |  |

The second record on ENST00000390690: Ensembl VEP's `inframe_insertion` row, the one vep-rs's gain row pairs with; vep-rs's rows for this record are not shown.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `rn_1_0` | not shown |
| Location | `21:30425596-30425597` |  |
| Allele | `AGC` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `inframe_insertion` |  |
| Extra | `ALLELE_NUM=1;IMPACT=MODERATE;STRAND=-1` |  |

With a FASTA vep-rs writes Ensembl VEP's row for the `RN=0` allele of the first and fourth records, byte for byte, and beside it a row for the real allele of the first record, its literal reading. For the second and third records it writes Ensembl VEP's `inframe_insertion` row alone, the trailing `RN=0` allele dropped by `split` in both engines. Without a FASTA both engines write the symbolic gain over the run for both alleles of every record, so nothing differs.

With `--variant_class` added, the `-` rows of the first and fourth records (`rn_0_1` and `rn_0_1_ruc3`, on both transcripts) differ in one cell: Ensembl VEP writes `VARIANT_CLASS=deletion` for its single-allele `AGCAGCAGC/-`, vep-rs `VARIANT_CLASS=sequence_alteration` for `AGCAGCAGC/-/AGCAGCAGCAGC` and `AGCAGCAGC/-/AGCAGCAGC`. `SO_variation_class` (ensembl-variation `Utils/Sequence.pm` lines 377-379) names a reference of bases whose every alternate is `-` a `deletion`; a string mixing `-` with a sequence allele keeps the default `sequence_alteration` (line 345), so the kept allele changes the record's class.

The first record on ENST00000390690, the real allele: the second of vep-rs's two rows, with no Ensembl VEP row.

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `rn_0_1` |
| Location |  | `21:30425588-30425596` |
| Allele |  | `AGCAGCAGCAGC` |
| Gene |  | `ENSG00000240432` |
| Feature |  | `ENST00000390690` |
| Feature_type |  | `Transcript` |
| Consequence |  | `protein_altering_variant` |
| cDNA_position |  | `331-339` |
| CDS_position |  | `317-325` |
| Protein_position |  | `106-109` |
| Amino_acids |  | `SCCS/SCCCS` |
| Codons |  | `aGCTGCTGCTcc/aGCTGCTGCTGCTcc` |
| Extra |  | `ALLELE_NUM=2;IMPACT=MODERATE;STRAND=-1` |

The allele is the run with one more unit, read against the run as a substitution of nine bases by twelve, so its term is `protein_altering_variant`. Ensembl VEP's reading of the same allele first in the list (`RN=1,0`, the second record) is an `inframe_insertion` of one unit after the run.

### 22.4 Why Ensembl VEP writes that

`Parser/VCF.pm` `_expand_tandem_repeat_allele_string` walks `RN` per allele (line 408). For `RN=0` the loop `for ( 1 .. $records )` at line 411 runs no iteration, so `$repeat` stays undefined and is pushed as the allele (line 424). Line 428 joins the reference and the alleles with `/`: `AGCAGCAGC//AGCAGCAGCAGC` for `RN=0,1`, `AGCAGCAGC/AGCAGCAGCAGC/` for `RN=1,0`.

`Parser.pm` line 978 recognises a multi-allelic string by `/.+\/.+\/.+/`, which neither string matches, an empty part having no character, so the record is not split per allele. Line 881 minimises it as one indel, where Perl's `split` on `/` drops a trailing empty part and keeps an empty middle one, and the first minimised feature is kept: `AGCAGCAGC/-` for `RN=0,1`, `-/AGC` for `RN=1,0` and `RN=2,0`.

### 22.5 Why vep-rs's rows are the correct ones

The record states two alleles, and both are annotated or neither. An allele of no repeat sequence is a real allele, the loss of the whole run, and Ensembl VEP itself annotates it as a deletion when it is the one that survives. The four-unit allele is a real allele too, and Ensembl VEP itself annotates it as an insertion when it comes first. Which allele is lost depends on the order the ALTs are listed in, not on what they say. vep-rs reads an empty middle part as the `RN=0` allele's sequence of no bases and keeps the real allele beside it.

### 22.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425587\trn_0_1\tN\t<CNV:TR>,<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=4\n' > rnzero.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache -v /path/to/fasta:/fasta ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/rnzero.vcf -o /work/rnzero.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats --fasta /fasta/Homo_sapiens.GRCh38.dna.chromosome.21.fa
vep -i rnzero.vcf -o rnzero.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats --fasta /path/to/fasta/Homo_sapiens.GRCh38.dna.chromosome.21.fa
```

Swapping the two `RN` values (`RN=1,0`) loses the other allele.

### 22.7 How the comparator treats this class

The class is `cnv_tr_rn_zero_allele_dropped` in the 116.2 registry. The structural-variant comparator reads the input record when the reference ran with a FASTA and expanded the record (a row at or beside the run under an allele other than `tandem_repeat`). The shape is a two-ALT `<CNV:TR>,<CNV:TR>` whose `RN` list has exactly one entry Perl's loop runs no iteration for (0 or negative; `.` is skipped at line 410 and leaves no part, so it is not the shape).

From a vep-rs run without a FASTA the comparator sets aside the `tandem_repeat` rows at the run that carry the dropped allele's direction term. These are the `feature_elongation` rows when the `RN=0` allele comes first and the real allele is a gain or unchanged, because the reference kept the `RN=0` deletion (the first and fourth records). vep-rs reads both alleles as a gain, so their `feature_elongation` rows are one tuple per transcript, and that tuple is set aside.

When the `RN=0` allele comes last and the real allele is a gain (the second and third), the reference kept the gain and the comparator looks for a `feature_truncation` row. vep-rs reads every symbolic `<CNV:TR>` allele as a gain, so no row carries it: nothing is set aside here, and the record's one symbolic tuple per transcript is left to [class 12](cnv-tr.md#12-cnvtr-read-symbolically-not-expanded).

A vep-rs run that read a FASTA writes the dropped allele's literal reading instead for the first record, the whole allele over the run (`AGCAGCAGCAGC`, the run with one more unit). The comparator reads that shape too: it rebuilds the dropped allele's sequence from `RUS` and `RUC` (or `RB`) as the expansion builds it and sets aside the vep-rs rows at the run under that allele. For the other three records such a run writes what the reference writes; the fourth record's unchanged allele gets no row.

A real allele that is a loss (`RUC=2` over the three-unit run) shares the symbolic tuple with the `RN=0` allele in a vep-rs run without a FASTA. Nothing is then left for this class to judge, and the reference's row is left to class 12. Behind a leading `RN=0`, a vep-rs run with a FASTA writes that allele's literal reading, `AGCAGC` as an `inframe_deletion` over the run, and the comparator sets those rows aside as it sets aside the gain's.

An unchanged real allele before a trailing `RN=0` (`RN=1,0;RUC=3`) leaves the reference an allele string with no alternate and no row, which is [class 19](defect-l.md#19-an-empty-or-unchanged-cnvtr-expansion-dropped-without-a-row), not this one. The allele the reference kept is left to class 12, as every expanded `<CNV:TR>` is. Without a FASTA both engines write the symbolic rows for both alleles and nothing is set aside. The SNP/indel comparator, which opens the input only to recognise gVCF variant sites under `--input-dir` or `--input-vcf`, does not judge the class.
