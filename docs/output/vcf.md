# VCF format

This page describes the VCF output: the input records written back with their consequences in a CSQ INFO field. Selected with `--vcf` or `--output_format vcf`, either of which also switches on `--symbol`, `--biotype` and `--numbers` (in VEP only the `--vcf` flag does). One line per input record: the record's own columns, with every allele's consequences comma-joined into the CSQ INFO field. An earlier CSQ entry is replaced.

## Header

The input's meta lines are kept, minus any earlier CSQ definition or `##VEP` line. A `##fileformat` line is added when the input lacks one. Three lines then precede the column line; the `##INFO` definition is one line, however wide it renders here.

```
##VEP="v116.2" API="v116" time="[TIME]" cache="[PATH]/116" 1000genomes="phase3" COSMIC="98" ...
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">
##VEP-command-line='vep --assembly GRCh37 ...'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
```

## CSQ field structure

The Format list is VEP's fixed prefix, `Allele` to `Existing_variation`, followed by the active flag fields and the built-in plugin fields. With the flags `--vcf` implies, the list is the one in the header above. Its 23 fields, numbered by position:

1. `Allele`: the ALT allele of the entry; the one field that is never empty.
2. `Consequence`: the SO consequence terms, joined by `&`.
3. `IMPACT`: `HIGH`, `MODERATE`, `LOW` or `MODIFIER`.
4. `SYMBOL`: the gene symbol (`--symbol`).
5. `Gene`: the Ensembl gene stable ID.
6. `Feature_type`: `Transcript`.
7. `Feature`: the Ensembl transcript stable ID.
8. `BIOTYPE`: the transcript biotype (`--biotype`).
9. `EXON`: `n/total` or `a-b/total` over the exons the variant strictly overlaps (`--numbers`).
10. `INTRON`: the same over the introns (`--numbers`).
11. `HGVSc`: empty unless `--hgvs`.
12. `HGVSp`: empty unless `--hgvs`.
13. `cDNA_position`: position or range in the cDNA.
14. `CDS_position`: position or range in the coding sequence.
15. `Protein_position`: position or range in the protein.
16. `Amino_acids`: reference/variant peptides.
17. `Codons`: reference/variant codon windows.
18. `Existing_variation`: co-located known variant IDs, joined by `&`.
19. `DISTANCE`: up- and downstream entries only.
20. `STRAND`: `1` or `-1`.
21. `FLAGS`: `cds_start_NF`, `cds_end_NF`.
22. `SYMBOL_SOURCE`: the source of the symbol (`--symbol`).
23. `HGNC_ID`: the HGNC identifier (`--symbol`).

Positions, peptides and codons follow the rules described for the [default format](vep-and-tab.md#differences-a-reader-will-notice). An absent value is empty, not `-`, except `Allele`. An intergenic entry carries `Allele`, `Consequence` and `IMPACT` alone, every later field empty:

```
G|intergenic_variant|MODIFIER||||...
```

The CSQ value of record 806345 from the [default format example](vep-and-tab.md#example), taken from the corpus's `expected/vcf.vcf.gz`. The file joins its four entries with `,` on one line:

```
CSQ=C|missense_variant|MODERATE|FLVCR1|ENSG00000162769|Transcript|ENST00000366971|protein_coding|6/10||||1433|1235|412|G/A|gGa/gCa|||1||HGNC|24682,
    C|missense_variant|MODERATE|FLVCR1|ENSG00000162769|Transcript|ENST00000419102|protein_coding|5/9||||631|632|211|G/A|gGa/gCa|||1|cds_start_NF|HGNC|24682,
    C|downstream_gene_variant|MODIFIER|FLVCR1|ENSG00000162769|Transcript|ENST00000474693|processed_transcript|||||||||||2364|1||HGNC|24682,
    C|intron_variant&non_coding_transcript_variant|MODIFIER|FLVCR1|ENSG00000162769|Transcript|ENST00000483790|processed_transcript||2/5||||||||||1||HGNC|24682
```

## Encoding

Inside a CSQ value: `,` and `|` become `&`, `;` becomes `%3B`, whitespace becomes `_`. Lists join with `&`.
