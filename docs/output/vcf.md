# VCF format

With `--vcf` or `--output_format vcf`, vep-rs writes the input records back with their consequences in a CSQ INFO field. Either flag also switches on `--symbol`, `--biotype` and `--numbers` (in Ensembl VEP only the `--vcf` flag does, [a documented difference](../divergences/output-format-vcf-fields.md)). One line per input record: the record's own columns, with every allele's consequences comma-joined into the INFO field named by `--vcf_info_field` (`CSQ` by default). An earlier entry under that name is replaced, unless `--keep_csq` keeps it.

## CSQ field structure

The Format list is Ensembl VEP's fixed prefix, `Allele` to `Existing_variation`, followed by the active flag fields and the built-in plugin fields. With the flags `--vcf` implies, the list is the one in the header below. Its 23 fields, numbered by position:

1. `Allele`: the ALT allele of the entry; the one field that is never empty.
2. `Consequence`: the SO consequence terms, joined by `&`.
3. `IMPACT`: `HIGH`, `MODERATE`, `LOW` or `MODIFIER`.
4. `SYMBOL`: the gene symbol (`--symbol`).
5. `Gene`: the Ensembl gene stable ID.
6. `Feature_type`: `Transcript`; `RegulatoryFeature` or `MotifFeature` on the rows `--regulatory` adds ([Regulatory and motif features](../regulatory.md)).
7. `Feature`: the Ensembl transcript stable ID, or the regulatory or motif feature's.
8. `BIOTYPE`: the transcript biotype (`--biotype`), or a regulatory feature's type whatever `--biotype` says.
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

Under `--fields` the Format list is that list alone, in its own order; a name nothing produces is an empty field. A plugin key listed twice, the plugin named twice on the command line or two plugins sharing a key, is one field, where Ensembl VEP 116.3 lists one per instance ([a key listed twice](../divergences/duplicate-output-column.md)).

Positions, peptides and codons follow the rules described for the [default format](vep-and-tab.md#differences-a-reader-will-notice). Under `--shift_3prime 1` the `EXON` and `INTRON` candidates are the features the moved allele overlaps, and the numbering follows the input's position, as in Ensembl VEP; a move across an exon boundary leaves both empty. An absent value is empty, not `-`, except `Allele`. An intergenic entry carries `Allele`, `Consequence` and `IMPACT` alone, every later field empty:

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

Inside a CSQ value: `,` and `|` become `&`, `;` becomes `%3B`, whitespace becomes `_`. Lists join with `&`. Whitespace in a gene symbol is percent-encoded as `%20` first, so it survives the `_` rule, unless `--no_escape`, which leaves it to that rule. `--no_escape` also keeps the `=` of an HGVSp value, written `%3D` otherwise.

## Header

The input's meta lines are kept, minus any earlier definition of the INFO field and any `##VEP` or `##vep-rs` line. A `##fileformat` line is added when the input lacks one. Then come the `##VEP=` line, the `##vep-rs=<version>` line, the field's `##INFO` definition, one `##<KEY>=<description>` line per field of each loaded built-in plugin, the command line and the column line. The `##INFO` definition is one line, however wide it renders here.

```
##VEP="v116.3" API="v116" time="[TIME]" cache="[PATH]/116" ensembl=116.0d85231 ensembl-compara=116.4038708 ensembl-funcgen=116.90049ea ensembl-io=116.6afb5dc ensembl-variation=116.db7394a 1000genomes="phase3" COSMIC="98" ...
##vep-rs=1.0.0
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">
##VEP-command-line='vep --assembly GRCh37 ...'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
```

The `##VEP=` line carries the tokens Ensembl VEP writes on a cache run, in its order (`OutputFactory/VCF.pm` lines 234-263). Its keys `ensembl=` to `ensembl-variation=` state the Ensembl API module versions of the Ensembl VEP 116.3 image (`ensemblorg/ensembl-vep:release_116.3`): they are the reference release's modules, none of which runs in vep-rs, whose own version the `##vep-rs=` line gives. The `time` token is UTC, where Ensembl VEP prints the process's local time (`Utils.pm` lines 467-468). The `##INFO=<ID=<short>_<field>,...>` lines for `--custom` sources follow the plugin lines, one per key of each source ([Custom annotation](../custom-annotation.md#output)).

## Records from a VCF input

A VCF input's line is written back as Ensembl VEP stores it, each column the bytes it was read from, sample columns included. A trailing empty column is dropped, a structural record's QUAL, FILTER and INFO read `.` where absent, empty or `0`, and a QUAL or FILTER absent from any record prints `.`. A record of any other input format gets [a line built for it](#records-from-other-input-formats).

## Records from other input formats

A record of another input format (Ensembl default, region, HGVS, SPDI, ID) has no VCF line of its own, so one is built as Ensembl VEP's `to_VCF_record` builds it (`ensembl-variation`, `VariationFeature.pm` lines 2340-2444, called from `OutputFactory/VCF.pm` line 337 at release 116.3). Its columns are CHROM, POS, the record's name as ID (`chr_start_alleles` for a record its ID column cannot name), REF, ALT, `.` for QUAL and FILTER, and the annotation alone as INFO.

A balanced change keeps its position and alleles (`21 30282169 multi C A,T`). An insertion or deletion takes the base before it as an anchor, from `--fasta` when the run has one and `N` otherwise, and its POS moves one base back. So the deletion `TA/-` at 30348410 is `21 30348409 NTA N` without a FASTA. A record on the reverse strand has its alleles reverse-complemented, in the ALT and REF columns and in the `Allele` field of its CSQ entries, where the default and tab outputs print the allele as written.

A structural class word (`DEL`, `DUP`, `INS`, `TDUP`, `INV` and the rest) is written as `<DEL>`, `<DUP>`, `<INS>` or `<TDUP>`, or as its class term in angle brackets (`<inversion>`), one base before its span with `END=` in INFO. The `GRCh38-input-ensembl`, `GRCh38-input-hgvs` and `GRCh38-input-spdi` corpora compare these lines with Ensembl VEP 116.3's. A COSMIC or HGMD identifier resolved from the cache, whose stored allele is the word `COSMIC_MUTATION` or `HGMD_MUTATION`, takes that word in angle brackets as its ALT ([Identifier input](../divergences/id-input.md)).

The header of such a run has no `##contig` line, as Ensembl VEP 116.3 writes none, so htslib refuses to write the file as BCF (`Contig '21' is not defined in the header`). Adding the contig lines from the reference's FASTA index first lets `bcftools view -Ob` and `bcftools index` read it:

```shell
bcftools reheader --fai reference.fa.fai out.vcf -o with-contigs.vcf
bcftools view -Ob -o out.bcf with-contigs.vcf
```

## Keeping an earlier annotation (`--keep_csq`)

An input that Ensembl VEP or vep-rs annotated before carries a `CSQ=` entry and the header lines that describe it. Without `--keep_csq` the first entry under the field's name is removed from each line and its `##INFO` and `##VEP` header lines are dropped, so the output describes one annotation. With `--keep_csq` every input header line stays and the new entry is appended under the same name. The header then carries two `##INFO=<ID=CSQ,...>` lines (the old Format list and the new one) and two `##VEP=` lines, and a line's INFO column carries two `CSQ=` keys, the old entry first. An empty entry under the name is removed alone, where Ensembl VEP's substitution takes the next field too ([the divergence page](../divergences/empty-csq-entry.md)).

The doubled header and INFO key are Ensembl VEP's own output, reproduced as is; a reader that keeps one value per INFO key sees only one of them. `--vcf_info_field ANN` names the new entry `ANN` and the new header line `##INFO=<ID=ANN,...>`. An earlier `CSQ=` entry is then another key, kept beside it whether or not `--keep_csq` is given. The `keep-csq` and `keep-csq-ann` sets of the `GRCh38-writers` corpus run both on an input Ensembl VEP 116.3 annotated first.

## Where next

- [The default format](vep-and-tab.md#differences-a-reader-will-notice): the rules behind the positions, peptides and codons.
- [Custom annotation](../custom-annotation.md#output): the keys and header lines a `--custom` source adds.
- [Regulatory and motif features](../regulatory.md): the rows `--regulatory` adds.
- [Output formats](../output-formats.md): the other formats and the rules every format follows.
