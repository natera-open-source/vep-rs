# A GTF attribute key written twice on one line, read as the transcript source

**The defect in one sentence.** Ensembl VEP 116.3's `--gtf` source reads a key written twice on one line as a list of both values and, wherever it expects one value, prints the list joined by commas, prints the list's memory address, or loses the record. A transcript line's `gene_name "TMPRSS3"; gene_name "OTHER"` is the symbol `TMPRSS3,OTHER`; a second `transcript_version` makes the `Feature` `ENST00001113783.ARRAY(0xaaaafe585258)`; a second `gene_id` or `transcript_id` makes the transcript's records `intergenic_variant`, and on an exon line drops that exon. vep-rs refuses to build a cache from a file that repeats a key read as one value, naming the line and the key, and writes no row.

**Why this is a defect.** A GTF gives a gene, a transcript, a biotype, a version or a protein one value per line, and repeats a key only where the value is a list, as the Ensembl files repeat `tag`. The parser keeps every value of a repeated key, as the `tag` list needs, but the source reads the list where it expects text. Ensembl VEP joins a list printed as a field, prints its address when the list is appended to an identifier, and never matches it against the plain identifier on the transcript's other lines, so the record or an exon is lost.

None of the three outputs is a value the file carries. The joined values and the address are written without a warning; the lost record alone draws one, naming the parent that was not found.

**What vep-rs does instead.** `vep-cache-builder` refuses the line: `ERROR: [cache-builder] Failed to parse GTF: <file> line <n>: GTF attribute <key> is repeated; Ensembl VEP reads a repeated key as a list`, exit 1, and no cache is written from the file. The keys refused are the twenty the readers take as one value: `gene_id`, `transcript_id`, `exon_id`, `protein_id`, `version`, `protein_version`, `transcript_version`, `transcript_biotype`, `transcript_type`, `biotype`, `gene_name`, `gene`, `name`, `Name`, `transcript_support_level`, `tsl`, `ccdsid`, `ccds_id`, `dbxref` and `Dbxref`. `tag` is read as the list it is, and a repeated key no reader takes keeps its last value. A key written bare and then with a value, `gene_id; gene_id "G"`, is one value to both engines.

**Where it is in Ensembl VEP's source** (`ensembl-io` at release/116, `ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/IO/Parser/GTF.pm`, lines 116-136, `get_attributes`: line 128 tests whether the key already holds a defined value, lines 129-131 gather the earlier value and the new one into a list and store its reference, line 133 stores a first value. A key whose first piece has no value is undefined, so a second piece replaces it.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/GTF.pm`, lines 156-158, `_record_get_id`: a record's identifier is its `<type>_id` attribute, the list reference when the key is repeated. Lines 133-141, `_record_get_parent_ids`: a transcript's parents are its `gene_id` split on `,`, an exon's or a CDS line's its `transcript_id`, so the split of a reference's text names a parent no record has.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, line 464: `transcript_id` as the stable identifier; line 481: `version` or `transcript_version`; line 638: `gene_id`; lines 662-666: `gene_name` or `gene` as the symbol; line 674: `transcript_support_level` or `tsl`, written `tsl` plus the value; line 723: `protein_id`. `File/GTF.pm`, lines 176-182, `_record_get_biotype`: `transcript_biotype`, `transcript_type` or `biotype`, else the source column.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1474: `Feature` gains `.` and the transcript's version, a reference's text when it is a list; lines 1569-1574: `ENSP` gains the translation's version the same way; lines 1618-1621: `TSL` prints the digits after `tsl`, which `tslARRAY(0x...)` has none of. `Bio/EnsEMBL/VEP/OutputFactory/VEP_output.pm`, lines 104 and 117, through `Bio/EnsEMBL/VEP/Utils.pm` lines 185-193, `convert_arrayref`: a list in a column or an `Extra` field prints its values joined by commas; the JSON writer prints the list itself.

**What correcting it would change.** Were the source to read one value per key it takes as one, as every Ensembl GTF writes one, a line repeating such a key would print one of its values or refuse the file, and no record or exon would be lost. The release files repeat `tag` alone, so their rows would not change.

**The class it produces.** `gtf_repeated_attribute_key`: every row Ensembl VEP's `--gtf` source writes from a file whose gene, transcript, exon, CDS or stop codon line repeats one of the twenty keys above, against vep-rs's refusal to build a cache from the file. Ensembl VEP prints a repeated symbol, biotype or protein identifier as the values joined by commas (a JSON list; `&`-joined in a VCF `CSQ` value). It prints a repeated version, or a repeated protein identifier under `--protein_version`, as the list's address, and no `TSL` for a repeated support level. A transcript whose gene or transcript line repeats `gene_id` or `transcript_id` is written as `intergenic_variant`, and an exon line repeating `transcript_id` loses that exon.

A repeated `exon_id`, `protein_version` or CCDS key, which the source never reads, leaves its rows as they are; vep-rs refuses the line all the same, since the value it reads from the key is not single. The class is documented only: no committed GTF or GFF3 repeats a key the readers take as one value, and no committed flag set runs Ensembl VEP over a file that does.

**How often the shape occurs.** Tier 3: the shape is a GTF line of a type the source reads that writes one of the twenty keys twice. The 964 feature lines of the committed GTF carry none, 14 committed flag sets run a cache built from a GTF or GFF3, and no committed flag set runs Ensembl VEP over a file with the shape. The GRCh38 release-116 and GRCh37 release-87 GTFs, 11,248,794 and 2,613,765 feature lines, repeat `tag` alone.

## 1. A TMPRSS3 transcript line carrying `gene_name` twice

**In short.** The GRCh38 release-116 GTF lines of the gene TMPRSS3 and its transcript `ENST00001113783`, with a second `gene_name "OTHER"` written after the transcript line's `gene_name "TMPRSS3"`, make a two-gene file with the gene PWP2 and its transcript `ENST00000456705` as a control. On the record `tmprss3_codon6`, `21 42395401 G>A`, a `missense_variant` on that transcript, Ensembl VEP 116.3 reading the file through `--gtf` writes `SYMBOL=TMPRSS3,OTHER` in every format. vep-rs refuses to build a cache from the file, naming line 4 and the key, so it writes no row.

### 1.1 What a reader needs to know first

A GTF line's ninth column is a list of `key "value"` pairs separated by `;`. Ensembl VEP's `--gtf` source takes a transcript's symbol from the transcript line's `gene_name`, its biotype from `transcript_biotype`, its version from `transcript_version` and its identifiers from `gene_id` and `transcript_id`, one value each; `tag`, the one key an Ensembl GTF repeats, is read as a list. A GTF read through `--gtf` is bgzip-compressed and tabix-indexed, and vep-rs reads a cache built from the plain file.

### 1.2 The record

The transcript line, as in the release file but for the second `gene_name`, with its gene line above it (the three exon, two CDS, stop codon and UTR lines omitted here):

```
21	ensembl_havana	gene	42371840	42396122	.	-	.	gene_id "ENSG00000160183"; gene_version "19"; gene_name "TMPRSS3"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
21	havana_tagene	transcript	42393803	42396122	.	-	.	gene_id "ENSG00000160183"; gene_version "19"; transcript_id "ENST00001113783"; transcript_version "1"; gene_name "TMPRSS3"; gene_name "OTHER"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TMPRSS3-224"; transcript_source "havana_tagene"; transcript_biotype "protein_coding"; tag "gencode_basic";
```

The input, one record on each gene:

```
21	42395401	tmprss3_codon6	G	A	.	.	.
21	44107424	pwp2_codon3	T	C	.	.	.
```

A third record, `tmprss3_codon60`, `21 42395000 C>T` inside the transcript's third exon, is run beside them for the exon-line shape of the table below. Both engines run `--symbol --biotype --tsl --protein --protein_version --transcript_version --gene_version --hgvs`. Ensembl VEP reads the file, bgzipped and indexed as `tmprss3_repeated_gene_name.gtf.gz`, through `--gtf` with a FASTA of chromosome 21; vep-rs builds a cache from the file and that FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.3:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783.1	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3,OTHER;BIOTYPE=protein_coding;ENSP=ENSP00000783588.1;SOURCE=tmprss3_repeated_gene_name.gtf.gz;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu
pwp2_codon3	21:44107424	C	ENSG00000241945	ENST00000456705.1	Transcript	missense_variant	52	8	3	F/S	tTc/tCc	-	IMPACT=MODERATE;STRAND=1;SYMBOL=PWP2;BIOTYPE=protein_coding;TSL=3;ENSP=ENSP00000411749.1;SOURCE=tmprss3_repeated_gene_name.gtf.gz;HGVSc=ENST00000456705.1:c.8T>C;HGVSp=ENSP00000411749.1:p.Phe3Ser
```

vep-rs, from `vep-cache-builder`, with no cache written and no row:

```
ERROR: [cache-builder] Failed to parse GTF: tmprss3_repeated_gene_name.gtf line 4: GTF attribute gene_name is repeated; Ensembl VEP reads a repeated key as a list
```

The tab output carries `TMPRSS3,OTHER` in its `SYMBOL` column, the JSON output the list `"gene_symbol":["TMPRSS3","OTHER"]`, and the VCF output's `CSQ` entry `TMPRSS3&OTHER`, the spelling a `,` inside a `CSQ` value takes. The rows of the control gene PWP2 are the control file's rows in every format. The same file with another key written twice, each measured on the same records:

| Key written twice | Line | Ensembl VEP 116.3 |
| --- | --- | --- |
| `transcript_biotype` | transcript | `BIOTYPE=protein_coding,lncRNA`, a JSON list, `protein_coding&lncRNA` in the VCF |
| `transcript_version` | transcript | `Feature` and the HGVSc reference name `ENST00001113783.ARRAY(0xaaaafe585258)`, another address in the next run |
| `protein_id` | CDS | `ENSP=ARRAY(0xaaaadaff3988).1` and `HGVSp=ARRAY(0xaaaadaff3988).1:p.Pro6Leu`; without `--protein_version`, `ENSP=ENSP00000783588,OTHER` and an HGVSp reference name that is again an address |
| `transcript_support_level` | transcript | no `TSL` field, where the control gene's single `"3"` prints `TSL=3` |
| `gene_id` or `transcript_id` | gene or transcript | the TMPRSS3 record reads `intergenic_variant` with no transcript row; the warnings file names the parent not found |
| `transcript_id` | exon | the exon is lost: the record inside it, `tmprss3_codon60`, reads `downstream_gene_variant` with `DISTANCE=324` in place of `missense_variant` |
| `exon_id`, `ccds_id`, `protein_version` | exon, transcript, CDS | the control file's rows |

vep-rs refuses each file at its line, naming the key.

### 1.4 Why Ensembl VEP writes that

The source reads the ninth column through `get_attributes`, which finds `gene_name` already defined when the second piece arrives and stores the two values as a list reference. `_add_identifiers` takes the symbol from that hash, and the output writer joins a list with commas, so the symbol is `TMPRSS3,OTHER`. A version is appended to an identifier with `.`, where Perl writes a reference as its address.

A record whose `gene_id` or `transcript_id` is a list is filed under the address, which the plain identifier on the transcript's other lines never names, so the records that depend on it are pruned. The warnings file says `Parent entries with the following IDs were not found or skipped due to invalid types`, naming the address or the identifier, and the output carries no trace.

### 1.5 Why vep-rs's refusal is the correct outcome

The file gives the transcript two symbols, two biotypes or two identifiers where the format and every reader of it expect one, and no reading can tell which the file means. Printing both joined, printing a memory address or dropping the transcript each writes a row no file carries, and the dropped transcript changes the consequence of every record inside it. Refusing the file at the line that repeats the key, with the key named, leaves the correction to the file's author.

### 1.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gtf tmprss3_repeated_gene_name.gtf --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
```

vep-rs exits 1 with the `ERROR: [cache-builder]` line above and writes no cache. Ensembl VEP 116.3 with `--gtf tmprss3_repeated_gene_name.gtf.gz --symbol`, the file bgzipped and indexed, writes `SYMBOL=TMPRSS3,OTHER`. Removing the second `gene_name` makes the engines agree, as does writing the two values on one `tag` key, which both engines read as a list.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/gtf_repeated_attribute_key.toml` holds the file, the records, Ensembl VEP's rows and vep-rs's refusal. The class covers every row of a transcript whose gene, transcript, exon, CDS or stop codon line repeats one of the twenty keys the readers take as one value, in every format, against the refusal. No committed dataset reaches the shape, so the class is documented only.

## 2. A last value padded with blanks at the end of its line

**In short.** The same TMPRSS3 lines with the transcript line's `gene_name "TMPRSS3"` moved to the end of the line, written without its closing `;` and followed by three blanks. Ensembl VEP 116.3 reading the file through `--gtf` keeps the blanks in the symbol and writes `SYMBOL=TMPRSS3   `; vep-rs reads the line without the blanks that pad it and writes `SYMBOL=TMPRSS3`. A GFF3 whose gene line ends `;Name=TMPRSS3   ` differs the same way. A file whose lines end in `\r\n` reads as one whose lines end in `\n` on both engines, in every format, so the line end itself makes no difference. This is a second mechanism, with its own entry, `gxf_last_value_trailing_blanks`, and its own reproduction.

### 2.1 What a reader needs to know first

A GTF value is the quoted text of its pair, and the pair ends at the next `;`; the last pair of an Ensembl line ends with one too. A line may be padded with blanks before its line end, and a file written on another system may end its lines in `\r\n`. Both engines read a tabix-indexed file through the same library, which hands them each line without its `\n` and without a `\r` before it, so neither engine sees the `\r`. vep-rs reads the plain file line by line and takes each line without its line end and without the blanks that pad either end of it.

### 2.2 The record

The transcript line, with `gene_name "TMPRSS3"` moved to the end, unterminated and followed by three blanks (shown here as `   `):

```
21	havana_tagene	transcript	42393803	42396122	.	-	.	gene_id "ENSG00000160183"; gene_version "19"; transcript_id "ENST00001113783"; transcript_version "1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TMPRSS3-224"; transcript_source "havana_tagene"; transcript_biotype "protein_coding"; tag "gencode_basic"; gene_name "TMPRSS3"   
```

The GFF3 form, the gene line of the same gene with its `Name` pair moved to the end and followed by three blanks:

```
21	ensembl_havana	gene	42371840	42396122	.	-	.	ID=gene:ENSG00000160183;biotype=protein_coding;description=transmembrane serine protease 3 [Source:HGNC Symbol%3BAcc:HGNC:11877];gene_id=ENSG00000160183;logic_name=ensembl_havana_gene_homo_sapiens;version=19;Name=TMPRSS3   
```

The input and the flags are those of section 1.2.

### 2.3 What each engine writes

Ensembl VEP 116.3, from the GTF, the three blanks inside the `SYMBOL` value:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783.1	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3   ;BIOTYPE=protein_coding;ENSP=ENSP00000783588.1;SOURCE=tmprss3_last_value_blanks.gtf.gz;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu
```

vep-rs:

```
tmprss3_codon6	21:42395401	A	ENSG00000160183	ENST00001113783.1	Transcript	missense_variant	249	17	6	P/L	cCg/cTg	-	IMPACT=MODERATE;STRAND=-1;SYMBOL=TMPRSS3;BIOTYPE=protein_coding;ENSP=ENSP00000783588.1;HGVSc=ENST00001113783.1:c.17C>T;HGVSp=ENSP00000783588.1:p.Pro6Leu
```

Ensembl VEP's tab output carries `TMPRSS3   ` in its `SYMBOL` column, its JSON output `"gene_symbol":"TMPRSS3   "`, and its VCF output `TMPRSS3%20%20%20`, each blank escaped as the `CSQ` writer escapes a blank in a symbol. The GFF3 form gives the same symbols in the four formats, Ensembl VEP's with the blanks and vep-rs's without. The same files with every line ending in `\r\n` write, on both engines and in all four formats, the rows the files with `\n` line ends write.

### 2.4 Why Ensembl VEP writes that

The library that reads a tabix-indexed line strips its `\n` and a `\r` before it and nothing else, so the line reaches the parser with its padding blanks. `get_attributes` splits the column on `;`, so the last pair runs to the end of the line, and its value is everything after the key with the quotes removed: `TMPRSS3` and the three blanks. The symbol is printed as read, the VCF writer escaping each blank as `%20`.

### 2.5 Why vep-rs's row is the correct one

The file names the gene `TMPRSS3`: the blanks stand outside the quotes, after the value's closing `"`, where the format puts nothing but the `;` and the separator before the next key. A pair written with its `;` reads `TMPRSS3` on both engines, so the blanks print only because the author left the `;` off the last pair, and the symbol `TMPRSS3   ` names no gene. vep-rs reads a padded line as the unpadded line, the reading every pair but the last gets on both engines.

### 2.6 Try it yourself

```
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gtf tmprss3_last_value_blanks.gtf --genome-fasta <corpus>/reference.fa.gz \
  --output-dir cache
vep -i record.vcf -o out.txt --offline --json_cache cache \
  --fasta <corpus>/reference.fa.gz --assembly GRCh38 --symbol
```

vep-rs writes `SYMBOL=TMPRSS3`. Ensembl VEP 116.3 with `--gtf tmprss3_last_value_blanks.gtf.gz --symbol`, the file bgzipped and indexed, writes `SYMBOL=TMPRSS3   `. Closing the pair with `;` makes the engines agree, and so does deleting the blanks. A file whose lines end in `\r\n` is indexed and read by both engines as the same file with `\n` line ends. A line padded with blanks before its first column, or a comment line opening with blanks before its `#`, cannot be indexed by tabix, so Ensembl VEP cannot read such a file. vep-rs reads the first as the unpadded line and skips the second.

### 2.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/gxf_last_value_trailing_blanks.toml` holds the GTF, the record and both engines' `SYMBOL` values as the VCF output spells them, `TMPRSS3%20%20%20` against `TMPRSS3`. The class, `gxf_last_value_trailing_blanks`, covers a `SYMBOL`, `BIOTYPE` or `ENSP` value, and its JSON key, of a transcript whose GTF or GFF3 line ends in a pair that carries a value, closes without `;` and is followed by blanks. Ensembl VEP's value carries the blanks; its text with the trailing blanks removed must equal vep-rs's.

A support level is not covered: Ensembl VEP reads its digits out of the padded value, so both engines print `TSL=3` for a padded `"3"`. No committed dataset reaches the shape, so the class is documented only.
