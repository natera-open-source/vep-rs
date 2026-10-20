# Defect H: a gVCF variant site is dropped with its real allele

**The defect in one sentence.** Ensembl VEP drops a gVCF variant site such as `A G,<NON_REF>` whole as an unsupported type, the real allele with it.

A gVCF writes `<NON_REF>` beside the real alleles of every variant site. Ensembl VEP's parser sends a line to its structural-variant path whenever an ALT carries an angle-bracket token, so the site is read as a structural variant of type `G/<NON_REF>`, which has no Sequence Ontology term.

**Why this is a defect.** `<NON_REF>` is not an allele; it is the gVCF's statement that the site may carry alleles the caller did not see. The site's real allele is a variant like any other, and dropping it loses every variant of a gVCF.

**What vep-rs does instead.** vep-rs reads `<NON_REF>` as the absence of a variant: it removes the token from the ALT list when another allele is present and annotates the real alleles, numbered as they stand in the ALT column. A record whose only ALT is `<NON_REF>`, a reference block, is one `intergenic_variant` row.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):
- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 237-245: a line is routed to `create_StructuralVariationFeatures` when REF and the ALTs together are not a run of `ACGT` and an ALT matches `[<\[\]][^\*]+[>\]\[]`, which `<NON_REF>` does.
- `Parser/VCF.pm`, lines 469-481: the ALTs are joined (`G/<NON_REF>`) and passed to `get_SO_term`; with no term for them, the line is marked `vep_skip`.
- `Parser.pm`, lines 712-717: the last branch of `get_SO_term` strips the angle brackets and finds no `%SO_TERMS` entry for `G/NON_REF`.
- `Parser.pm`, lines 821-827: the marked line is dropped before the input buffer.

**What correcting it would change.** The real allele would be annotated, as vep-rs annotates it; [class 15](#15-a-gvcf-variant-site-dropped-with-its-real-allele) would vanish.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 15. A gVCF variant site dropped with its real allele

**In short.** Ensembl VEP sends a line to its structural-variant path whenever an ALT carries an angle-bracket token. A gVCF variant site such as `A G,<NON_REF>` is therefore read as a structural variant of type `G/<NON_REF>`, which has no Sequence Ontology term, and the line is dropped whole with its real allele. vep-rs reads `<NON_REF>` as the absence of a variant and annotates the real alleles, numbered as they stand in the ALT column.

### 15.1 What a reader needs to know first

A GATK gVCF writes every variant site with `<NON_REF>` beside its real alleles (`A G,<NON_REF>`), the symbolic token standing for "any other allele" so that the site can be merged later. Ensembl VEP's parser sends a line to its structural-variant path whenever an ALT carries an angle-bracket token, reads the type of the whole ALT list (`G/<NON_REF>`), finds no Sequence Ontology term for it and drops the line as an unsupported structural variant. The `G` allele is never annotated. vep-rs reads `<NON_REF>` as the absence of a variant and annotates the real allele; a record whose only ALT is `<NON_REF>` (a reference block) is one `intergenic_variant` row.

### 15.2 The record

```
21	30425650	gvcf_a_g	A	G,<NON_REF>	.	.	.
```

Executed with Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`, the 116 GRCh38 cache, `--offline --cache --allele_number`).

### 15.3 What each engine writes

Ensembl VEP 116.2, no row, and in its warnings file:

```
WARNING: line 1 skipped (21 30425650 gvcf_a_g A G,<NON_REF> . . .): G/NON_REF is not a supported structural variant type
```

vep-rs, two rows, the SNV annotated as `A G` alone would be, one table per transcript against Ensembl VEP's missing row. First ENST00000334068:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `gvcf_a_g` |
| Location |  | `21:30425650` |
| Allele |  | `G` |
| Gene |  | `ENSG00000186971` |
| Feature |  | `ENST00000334068` |
| Feature_type |  | `Transcript` |
| Consequence |  | `upstream_gene_variant` |

Record `gvcf_a_g` on ENST00000390690:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | no row | `gvcf_a_g` |
| Location |  | `21:30425650` |
| Allele |  | `G` |
| Gene |  | `ENSG00000240432` |
| Feature |  | `ENST00000390690` |
| Feature_type |  | `Transcript` |
| Consequence |  | `missense_variant` |

The same happens with a reference base of `N` (`N G,<NON_REF>`, where vep-rs writes `coding_sequence_variant`, the peptide being unreadable from an ambiguous reference) and with an indel (`N NA,<NON_REF>`, `frameshift_variant`).

### 15.4 Why Ensembl VEP writes that

`Parser/VCF.pm` lines 237-245 route a line to `create_StructuralVariationFeatures` when REF and the ALTs together are not a run of `ACGT` and an ALT matches `[<\[\]][^\*]+[>\]\[]`, which `<NON_REF>` does. Lines 469-481 join the ALTs (`G/<NON_REF>`) and call `get_SO_term`, whose last branch (`Parser.pm` lines 712-717) strips the angle brackets and finds no `%SO_TERMS` entry for `G/NON_REF`. The line is then marked `vep_skip` and dropped before the input buffer (`Parser.pm` lines 821-827).

### 15.5 Why vep-rs's rows are the correct ones

`<NON_REF>` is not an allele; it is the gVCF's statement that the site may carry alleles the caller did not see. The site's real allele is a variant like any other, and dropping it loses every variant of a gVCF. vep-rs removes the token from the ALT list when another allele is present and parses the record as its real alleles, numbered as they stand in the ALT column.

### 15.6 Try it yourself

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\tgvcf_a_g\tA\tG,<NON_REF>\t.\t.\t.\n' > gvcf.vcf
docker run --rm -v "$PWD":/work -v /path/to/vep-cache:/cache ensemblorg/ensembl-vep:release_116.2 \
  vep -i /work/gvcf.vcf -o /work/gvcf.perl.txt --offline --cache --dir_cache /cache --cache_version 116 \
  --assembly GRCh38 --species homo_sapiens --allele_number --no_stats
vep -i gvcf.vcf -o gvcf.vep-rs.txt --offline --json_cache tests/golden/116/GRCh38-release116/json_cache --assembly GRCh38 --no_stats
```

### 15.7 How the comparator treats this class

The class is `gvcf_variant_site_dropped` in the 116.2 registry of both comparators, and both set aside vep-rs's rows at the site. They recognise the site by the record's shape: `<NON_REF>` beside a sequence allele and no true `SVTYPE`. (A list of `<NON_REF>` alone, or a block beside a symbolic allele, is the registry class `reference_skipped_record`; under a true `SVTYPE` both engines read the block as an allele of the structural record and agree.)

The SNP/indel comparator reads the shape from the input VCFs when `--input-dir` (paired with each output by stem) or `--input-vcf` names them, as the runners do. For an ID-less record it reads the shape from the record's own name in vep-rs's output, `CHROM_POS_REF/ALT/<NON_REF>`, the name Ensembl VEP and vep-rs give such a record, at a Location the reference wrote nothing for.

The structural-variant comparator reads the input record and sets aside vep-rs's rows at the Locations and under the alleles vep-rs's parser writes for the site. Those follow vep-rs's own trimming. The anchor base is trimmed from an indel's alleles, a lone allele is minimised from both ends (`ATG AG,<NON_REF>` is the deletion of `T` at POS + 1, Allele `-`), and an insertion is written at its pair of flanking bases. The skipped-record class names every other unsupported type and leaves this one to its own count.
