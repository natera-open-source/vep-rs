# Defect L: an empty or unchanged `<CNV:TR>` expansion is dropped without a row

**The defect in one sentence.** With a FASTA, Ensembl VEP drops a `<CNV:TR>` record without a row or a warning when its expansion to a literal allele string yields no alternate allele, or one identical to the reference run.

The expansion may yield no alternate allele (`RN=.` on the only allele, a `RUC` that is not a positive number) or an allele identical to the reference run (`RUC` times the unit length equal to `SVLEN`). The record then passes through the parser without a feature to annotate: no row in any format, and no warning.

**Why this is a defect.** A record the input states is a variant and the parser accepts should be annotated or refused with a message, never dropped in silence. Ensembl VEP itself annotates the same record as a symbolic `tandem_repeat` over the run when no FASTA is given.

**What vep-rs does instead.** vep-rs reads the record as the symbolic tandem repeat it is, with or without a FASTA, and writes the symbolic rows over the run.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):
- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`, lines 382-431: line 410 skips an `RN` of `.` (`next if $records eq '.'`). Lines 414 and 419 numify `RUC` as Perl does, so `abc` and `-2` give zero or fewer repeats and `$seq x 0` the empty string. The function returns `REF/` with no alternate, or `AGCAGCAGC/AGCAGCAGC` when the count matches the run.
- `Parser/VCF.pm`, `create_StructuralVariationFeatures`, lines 532-543: builds a VariationFeature from that string.
- `Parser.pm`, lines 553-592: pass the feature, and no allele reaches a consequence.

**What correcting it would change.** The record would be annotated symbolically over its run, as it is without a FASTA and as vep-rs annotates it; [class 19](#19-an-empty-or-unchanged-cnvtr-expansion-dropped-without-a-row) would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 19. An empty or unchanged `<CNV:TR>` expansion dropped without a row

**In short.** With a FASTA, Ensembl VEP expands a `<CNV:TR>` record to a literal allele string. When the expansion yields no alternate allele or an allele identical to the reference run, the record passes through the parser without a feature to annotate: no row and no warning. vep-rs writes the symbolic `tandem_repeat` rows over the run, with or without a FASTA, as Ensembl VEP itself does without one.

### 19.1 What a reader needs to know first

With a FASTA, Ensembl VEP expands a `<CNV:TR>` record to a literal allele string ([12.4](cnv-tr.md#124-why-ensembl-vep-writes-that)). The record's `RN`, `RUS`, `RUC` and `RB` fields are explained on [the `<CNV:TR>` page](cnv-tr.md#121-what-a-reader-needs-to-know-first). The expansion may yield no alternate allele (`RN=.` on the only allele, a `RUC` that is not a positive number) or an allele identical to the reference run (`RUC` times the unit length equal to `SVLEN`). The allele string then has no alternate, and the record passes through the parser without a feature to annotate: no row in any format, and no warning. Without a FASTA the same record is a symbolic `tandem_repeat` with rows, and vep-rs writes those rows with or without a FASTA.

### 19.2 The record

```
21	30425587	same_ruc	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=3
21	30425587	gain_rn_dot	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RN=.;RUS=AGC;RUC=4
21	30425587	gain_ruc_abc	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=abc
21	30425587	gain_ruc_neg	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=-2
21	30425587	gain1_ruc	N	<CNV:TR>	.	.	SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=4
```

The run at 30425588-30425596 is three copies of `AGC` in the GRCh38 sequence. Executed with `--fasta` (Ensembl's chromosome 21 of release 116) on Ensembl VEP 116.2.

### 19.3 What each engine writes

Ensembl VEP 116.2 with `--fasta`, two rows for the fifth record and none for the first four, with no warning; one table per row, the differing cells in bold. vep-rs, for each of the five records (the gain's literal reading pairs with the fifth under [class 12](cnv-tr.md#12-cnvtr-read-symbolically-not-expanded)), one row, its ID column left out in the last table.

The fifth record, `gain1_ruc`, on ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `gain1_ruc` | no row |
| Location | `21:30425596-30425597` |  |
| Allele | `AGC` |  |
| Gene | `ENSG00000186971` |  |
| Feature | `ENST00000334068` |  |
| Feature_type | `Transcript` |  |
| Consequence | `upstream_gene_variant` |  |

The fifth record, `gain1_ruc`, on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `gain1_ruc` | see the last table |
| Location | `21:30425596-30425597` |  |
| Allele | `AGC` |  |
| Gene | `ENSG00000240432` |  |
| Feature | `ENST00000390690` |  |
| Feature_type | `Transcript` |  |
| Consequence | `inframe_insertion` |  |

The vep-rs row, for each of the five records:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Location | no row for the first four | `21:30425588-30425596` |
| Allele |  | `tandem_repeat` |
| Gene |  | `ENSG00000240432` |
| Feature |  | `ENST00000390690` |
| Feature_type |  | `Transcript` |
| Consequence |  | `feature_elongation,coding_sequence_variant` |

### 19.4 Why Ensembl VEP writes that

`_expand_tandem_repeat_allele_string` (`Parser/VCF.pm` lines 382-431) skips an `RN` of `.` (`next if $records eq '.'`, line 410) and numifies `RUC` as Perl does (`abc` and `-2` give zero or fewer repeats, `$seq x 0` the empty string, lines 414 and 419). It returns `REF/` with no alternate, or `AGCAGCAGC/AGCAGCAGC` when the count matches the run. `create_StructuralVariationFeatures` (lines 532-543) builds a VariationFeature from that string, which `Parser.pm` lines 553-592 pass, and no allele reaches a consequence.

### 19.5 Why vep-rs's rows are the correct ones

A record the input states is a variant and the parser accepts should be annotated or refused with a message, never dropped in silence; and the same record is annotated by Ensembl VEP itself when no FASTA is given. vep-rs reads the record as the symbolic tandem repeat it is in both cases.

### 19.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425587\tsame_ruc\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=3\n' > empty.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache -v /path/to/fasta:/fasta ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/empty.vcf -o /work/empty.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats --fasta /fasta/Homo_sapiens.GRCh38.dna.chromosome.21.fa
vep -i empty.vcf -o empty.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

### 19.7 How the comparator treats this class

The class is `cnv_tr_empty_expansion_dropped` in the 116.2 registry. The structural-variant comparator reads the input record and sets aside vep-rs's `tandem_repeat` rows at the run when the reference wrote nothing on the chromosome from the base before the run to the base after it. A reference that wrote the symbolic rows (no FASTA) or the literal rows (a real gain) leaves nothing to set aside.

The input record is a single `<CNV:TR>` of one of three shapes: its `RN` is `.`; its repeat counts, `RUC` or `RB` in whole units, truncated as Perl's `x` operator truncates them, repeat nothing; or its counts times the unit lengths equal `SVLEN`.

`RUC=3.5` on a three-unit run is included, as is an unchanged expansion over several repeat sequences, `RN=3;RUS=AGC,AGC,AGC;RUC=1,1,1` or `RN=2;RUS=AGC,AG;RUC=3,0` on the nine-base run, and a count in fullwidth digits, `RUC=４`, which Perl numifies to 0. A unit of no bases repeats nothing under `RUC`, `"" x n` being the empty expansion, and is set aside like `RN=.`. Under `RB` it divides by zero in Perl, which stops Ensembl VEP with no output, so no reference holds such a record.

The SNP/indel comparator reads no input and does not judge the class.
