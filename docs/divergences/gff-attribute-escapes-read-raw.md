# A GFF3 attribute value holding an escape, read as the transcript source or a promoter window

**The defect in one sentence.** Ensembl VEP 116.2's `--gff` and `--custom` GFF promoter sources read every attribute value of the ninth column as written, so a value carrying a GFF3 escape prints the escape in place of the character: `Name=A%2CB` prints the symbol `A%2CB` and `ID=PROM%2c1%3E` the `feature_id` `PROM%2c1%3E`.

**Why this is a defect.** GFF3 reserves `;`, `=`, `,`, `&` and `%` in the ninth column and has a value carry them as `%XX` escapes; Ensembl's own files write `%2C` and `%3B` in gene names. The escaping is the URL convention: any byte may be written as `%XX`, the hexadecimal digits in either case, and a reader decodes every such escape, `%3E` for `>` as much as `%2c` for `,`. The parser Ensembl VEP reads a GFF3 with decodes them in its `get_attribute_by_name` and `get_type` accessors, through its `decode_string`. Its `get_attributes`, the accessor both sources call, stores each value as written. Ensembl's database and cache carry the character, so the `--gff` source prints a symbol no other source prints.

**What vep-rs does instead.** vep-rs decodes every `%XX` escape in every attribute value it reads from a GFF3, in either letter case and in one pass, as the parser's own `decode_string` reads one. So the symbol, identifiers, biotype and support level of a cache built from the file, and the `feature_id` and `associated_gene` of a promoter window a `--custom` GFF source overlays, are the text the file encodes. Ensembl VEP's tests on a value, whether it is set and where its list splits on `,`, run on the written value here too; the decoding follows them. The pair split is the parser's, so `Name=a=b` reads `a` on both engines, and only the decoding differs.

Both engines decode a GTF value the same way, so the GTF form of the shape reaches no class.

**Where it is in Ensembl VEP's source** (`ensembl-io` at release/116, `ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/IO/Parser/GFF3.pm`, lines 190-198, `get_attributes`: `split(';', ...)` then `split('=', $attr)` into two variables (line 194), no decoding. Lines 174-183, `get_attribute_by_name`: `decode_string($value)`; lines 161-164, `get_type`: decoded likewise. `Bio/EnsEMBL/IO/Parser/GXF.pm`, lines 303-307, `decode_string`: `uri_unescape`.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, line 289: `attributes => $parser->get_attributes`, the hash every later read takes a value from.
- `BaseGXF.pm`, lines 632-703, `_add_identifiers`: line 638 the transcript line's `gene_id`, lines 662-666 its `gene_name` or `gene`, line 670 the gene line's `name` or `Name`, line 674 `transcript_support_level` or `tsl`, lines 683-685 the `tag` list. Line 464, `_create_transcript`: `transcript_id`; line 481: `version` or `transcript_version`; line 723, `_add_translation`: the CDS line's `protein_id`.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/GFF.pm`, line 277, `_record_get_id`: `ID`, `Name`, `id` or `name`; line 248, `_record_get_parent_ids`: `Parent` or `parent`; lines 304-335, `_record_get_biotype`: `biotype`, `transcript_type` or `transcript_biotype`.
- `GFF.pm`, lines 198-201, `_create_records`: a `gff_type=gencode_promoter` source's `feature_id` is the hash's `ID` and its `associated_gene` the hash's `gene_id`, the same hash of `BaseGXF.pm`, line 289, so a promoter window's two fields print as written.

**What correcting it would change.** Every value read through `get_attributes` would carry the characters the file encodes, as the values read through `get_attribute_by_name` do. The symbol, identifier, biotype or support level of a record whose value holds an escape, and the `feature_id` and `associated_gene` of a promoter window, would print as vep-rs prints them. The human release files hold no escape in a value the source reads, so their rows would not change.

**The class it produces.** `gff_attribute_escapes_read_raw`: a `SYMBOL`, `Gene`, `Feature`, `BIOTYPE`, `TSL` or `ENSP` value (the JSON `gene_symbol`, `gene_id`, `transcript_id`, `biotype`, `tsl` and `protein_id`) of a transcript row that Ensembl VEP's `--gff` source prints with a `%XX` escape. It also covers a promoter window's `feature_id` or `associated_gene` (the `<short name>_feature_id` and `<short name>_associated_gene` keys) that the `--custom` GFF source prints with one. Its text decoded, and in the VCF format re-spelt as the `CSQ` writer spells `,` and `;`, must equal vep-rs's. The class is documented only: no committed fixture carries an escape in a printed value, and no committed flag set runs Ensembl VEP over one that does.

An identifier whose `gene:` or `transcript:` prefix is written with an escaped colon loses the prefix once decoded; Ensembl VEP strips the prefix off the written text and prints the escaped prefix whole. The class covers that row with the prefix removed from the decoded text, the one place its rule reads the text with a prefix taken off.

**How often the shape occurs.** Tier 3: the shape is a GFF3 line whose read attribute value holds a `%XX` escape. The 964 feature lines of the committed GTF and the 58 records of the committed promoter GFF carry none in a read value. Of the committed flag sets, 14 run a cache built from a GTF or GFF3, one overlays that promoter file, and none runs Ensembl VEP over a file with the shape. The Saccharomyces cerevisiae R64-1-1 GFF3 names three genes with a comma (`ADE5%2C7`, `ARG5%2C6`, `DUR1%2C2`) and the Arabidopsis thaliana TAIR10 GFF3 115 genes with a semicolon (`CYCA2%3B3` among them). The human, mouse and fly files hold escapes in `description` alone, a value neither engine prints.

## 1. A TMPRSS3 gene line whose `Name` holds an escaped comma

**In short.** The GRCh38 release-116 GFF3 lines of the gene TMPRSS3 and its transcript `ENST00001113783`, with the gene's `Name=TMPRSS3` rewritten `Name=A%2CB`, the GFF3 spelling of the symbol `A,B`, make a two-gene file with the gene PWP2 and its transcript `ENST00000456705` as a control. On the record `tmprss3_codon6`, `21 42395401 G>A`, a `missense_variant` on that transcript, vep-rs writes `SYMBOL=A,B`. Ensembl VEP 116.2 reading the file through `--gff` writes `SYMBOL=A%2CB` in every format, and on the yeast release file's ARG5,6 gene the same way, `ARG5%2C6` where its own cache carries `ARG5,6`. A promoter window overlaid through `--custom`, `ID=PROM%2c1%3E;gene_id=ENSG%3B1`, prints `PROM%2c1%3E` and `ENSG%3B1` on Ensembl VEP in every format, `PROM,1>` and `ENSG;1` on vep-rs in the tab and JSON formats.

### 1.1 What a reader needs to know first

A GFF3 line's ninth column is a list of `tag=value` pairs separated by `;`. The format reserves `;`, `=`, `,`, `&` and `%` inside a value and spells each as a percent escape, `%2C` for `,`, so a parser splits the column first and decodes each value after. Any character may be escaped so, `%3E` for `>`, in either letter case. Ensembl VEP's `--gff` source takes a transcript's symbol from the transcript line's `gene_name` or `gene`, else from the gene line's `name` or `Name`, where an Ensembl GFF3 writes it. A `--custom` GFF source of `gff_type=gencode_promoter` overlays promoter windows, each with a `feature_id` from its `ID` and an `associated_gene` from its `gene_id`, read from the same hash.

### 1.2 The record

The gene line, as in the release file but for its `Name` (its `description` and `logic_name` left out here), with its transcript line below it (the three exon and two CDS lines omitted here):

```
21	ensembl_havana	gene	42371840	42396122	.	-	.	ID=gene:ENSG00000160183;Name=A%2CB;biotype=protein_coding;gene_id=ENSG00000160183;version=19
21	havana_tagene	mRNA	42393803	42396122	.	-	.	ID=transcript:ENST00001113783;Parent=gene:ENSG00000160183;Name=TMPRSS3-224;biotype=protein_coding;tag=gencode_basic;transcript_id=ENST00001113783;version=1
```

The promoter window, in a second file, `promoters.gff.gz`, bgzipped and indexed:

```
21	.	promoter_window	42395390	42395410	.	+	.	ID=PROM%2c1%3E;gene_id=ENSG%3B1
```

The input, one record on each gene:

```
21	42395401	tmprss3_codon6	G	A	.	.	.
21	44107424	pwp2_codon3	T	C	.	.	.
```

Both engines run `--symbol --biotype --tsl --protein --protein_version --transcript_version --gene_version --hgvs` with `--custom file=promoters.gff.gz,format=gff,gff_type=gencode_promoter`. Ensembl VEP reads the file, bgzipped and indexed as `tmprss3_escaped_name.gff3.gz`, through `--gff` with a FASTA of chromosome 21; vep-rs reads a cache built from the file and that FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.2:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183.19	ENST00001113783.1	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=A%2CB;BIOTYPE=protein_coding;ENSP=ENSP00000783588.1;SOURCE=tmprss3_escaped_name.gff3.gz;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu;GENCODE_PROMOTER=21:42395390-42395410;GENCODE_PROMOTER_associated_gene=ENSG%3B1;GENCODE_PROMOTER_feature_id=PROM%2c1%3E;GENCODE_PROMOTER_type=promoter_window
pwp2_codon3	21:44107424	C	ENSG00000241945.10	ENST00000456705.1	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;BIOTYPE=protein_coding;TSL=3;ENSP=ENSP00000411749.1;SOURCE=tmprss3_escaped_name.gff3.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

vep-rs:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183.19	ENST00001113783.1	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=A,B;BIOTYPE=protein_coding;ENSP=ENSP00000783588.1;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu;GENCODE_PROMOTER=21:42395390-42395410;GENCODE_PROMOTER_associated_gene=ENSG%3B1;GENCODE_PROMOTER_feature_id=PROM,1>;GENCODE_PROMOTER_type=promoter_window
pwp2_codon3	21:44107424	C	ENSG00000241945.10	ENST00000456705.1	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;BIOTYPE=protein_coding;TSL=3;ENSP=ENSP00000411749.1;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

The tab output carries the same two symbols in its `SYMBOL` column, and the JSON output `"gene_symbol": "A%2CB"` against `"gene_symbol": "A,B"`. In the VCF output Ensembl VEP's `CSQ` entry carries `A%2CB` and vep-rs's `A&B`, the spelling both engines give a `,` inside a `CSQ` value. Beyond the symbol and the window's two fields, the rows differ only in what the class [`file_cache_writes_no_source_field`](file-cache-writes-no-source-field.md) covers, the `SOURCE` item above: the consequence, the positions, the codons and the protein are the same. The rows of the control gene PWP2 differ only in what that class covers, in every format.

The window's two fields print `PROM%2c1%3E` and `ENSG%3B1` on Ensembl VEP in every format, the lower-case escape and the `%3E` as written. vep-rs prints `PROM,1>` in the default, tab and JSON formats and `PROM&1>` in its `CSQ` entry, since a `>` is re-spelt nowhere. It prints `ENSG;1` in the tab and JSON formats. In the default format it prints `ENSG%3B1`, since the `Extra` column re-spells a `;` inside a value that way on both engines, and in its `CSQ` entry the same, the spelling every `CSQ` value gives a `;`. So the `associated_gene` values of the default and VCF formats agree byte for byte while the `feature_id` values differ in all four.

### 1.4 Why Ensembl VEP writes that

The source reads the ninth column through `get_attributes`, which splits on `;` and `=` and stores `A%2CB` as written; `_add_identifiers` finds no `gene_name` on the transcript line and takes the gene line's `Name` from that hash. The decoding step of the same parser lives in `get_attribute_by_name`, an accessor the source never calls, so the escape reaches the output. The promoter source's `_create_records` takes the window's `ID` and `gene_id` from the same hash, so its two fields reach the output as written too, the lower-case `%2c` and the `%3E` with them: no escape of any shape is decoded.

### 1.5 Why vep-rs's row is the correct one

The file spells the symbol `A,B` in the one way GFF3 allows a comma to be written, and the parser's own decoding reads it back as `A,B`. A cache built from the file names the gene as the file names it, which is also how Ensembl's database, its cache and every other Ensembl VEP source name the same gene. On the yeast release file, `ADE5,7`, `ARG5,6` and `DUR1,2` are the names the Ensembl cache carries; the `--gff` source alone prints them with the escape.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gff3 tmprss3_escaped_name.gff3 --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
vep -i record.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 --symbol \
  --custom file=promoters.gff.gz,format=gff,gff_type=gencode_promoter
```

vep-rs writes `SYMBOL=A,B` and `GENCODE_PROMOTER_feature_id=PROM,1>`. Ensembl VEP 116.2 with `--gff tmprss3_escaped_name.gff3.gz --symbol` and the same `--custom` source, both files bgzipped and indexed, writes `SYMBOL=A%2CB` and `GENCODE_PROMOTER_feature_id=PROM%2c1%3E`. Restoring `Name=TMPRSS3`, `ID=PROM1` and `gene_id=ENSG1` makes the engines agree, and so does the same comma in a GTF, `gene_name "A%2CB"`, which both engines read as `A,B`.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/gff_attribute_escapes_read_raw.toml` holds the two files, the record, both engines' symbols and both engines' window fields. The class covers the `SYMBOL`, `Gene`, `Feature`, `BIOTYPE`, `TSL` and `ENSP` values, and their JSON keys, of a transcript whose GFF3 line, or whose gene line, holds a `%XX` escape in the value the source prints. It also covers the `feature_id` and `associated_gene` of a promoter window whose `ID` or `gene_id` holds one. No committed dataset reaches the shape, so the class is documented only.
