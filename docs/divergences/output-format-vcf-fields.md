# The CSQ fields of a VCF output selected with `--output_format vcf`

**The defect in one sentence.** Ensembl VEP 116.2 switches a VCF output's `SYMBOL`, `BIOTYPE`, `EXON` and `INTRON` fields on from the `--vcf` flag, not from the output format. A VCF output selected with `--output_format vcf` therefore names the four fields in its `CSQ` layout and leaves each of them empty on every entry. vep-rs switches the fields on from the output format and fills them under `--output_format vcf` as under `--vcf`.

**Why this is a defect.** `--vcf` and `--output_format vcf` select the same writer, and the writer's layout names the four fields whichever of the two selected it. A field the header declares and no entry fills tells a reader that a transcript has no symbol, no biotype and no exon or intron number, which the cache contradicts. Ensembl VEP fills the same fields for the same record under `--vcf`. Nothing but the spelling of the format choice differs between the two runs.

**What vep-rs does instead.** vep-rs switches the fields on from the output format. A VCF output fills `SYMBOL`, `BIOTYPE`, `EXON` and `INTRON`, and carries the `SYMBOL_SOURCE` and `HGNC_ID` columns, under `--output_format vcf` as under `--vcf`.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`):

- `Config.pm`, lines 500-508: the option set of the `vcf` flag sets `output_format` to `vcf` and switches on `symbol`, `biotype` and `numbers`. No option set is keyed on the value of `output_format`.
- `OutputFactory.pm`, lines 229-238: the `output_format` value is lower-cased and names the writer's class, so `--output_format vcf` reaches the VCF writer without that set.
- `OutputFactory/VCF.pm`, lines 100-118: the layout's fixed columns, `SYMBOL`, `BIOTYPE`, `EXON` and `INTRON` among them. Lines 434-447, `fields`, append the columns of the flags that are on.
- `Constants.pm`, line 85: the `symbol` flag's columns are `SYMBOL`, `SYMBOL_SOURCE` and `HGNC_ID`.
- `OutputFactory.pm`, lines 1481, 1520 and 1580: the exon and intron numbers, the symbol and the biotype are written only when `numbers`, `symbol` and `biotype` are on.

**What correcting it would change.** With the three flags switched on for every VCF output, Ensembl VEP's `--output_format vcf` would write the fields and the layout of its `--vcf`, and the two engines would agree.

**The class it produces.** `output_format_vcf_fields_left_empty`: the `SYMBOL`, `BIOTYPE`, `EXON` and `INTRON` fields of every `CSQ` entry, and the `SYMBOL_SOURCE` and `HGNC_ID` columns of the layout, in a VCF output selected with `--output_format vcf` on a run that leaves `--symbol`, `--biotype` or `--numbers` off. Under `--fields` it covers the listed columns of a flag left off. The reference's fields are empty where the candidate's hold the cache's values. The class is documented only: no comparator codes it, since no committed flag set reaches it.

**How often the shape occurs.** Tier 3: the shape is an `--output_format` value that names the VCF writer on a run without one of the three flags. Of the 392 committed flag sets, 1 passes `--output_format`, with the value `ensembl`, and none carries the shape.

## 1. One record's entries with and without their fields

**In short.** The record `21 45573402 C T` lies in the first exon of ENST00000859469, a `protein_coding` transcript of SLC19A1, and in the one intron of ENST00000763268, a `lncRNA`. Under `--output_format vcf`, Ensembl VEP 116.2 writes both entries with empty `SYMBOL`, `BIOTYPE`, `EXON` and `INTRON` fields. vep-rs writes `SLC19A1`, `protein_coding` and `1/7` on the first, and `lncRNA` and `1/1` on the second. Under `--vcf` both engines write vep-rs's entries.

### 1.1 What a reader needs to know first

A VCF output carries a record's annotation in its `CSQ` INFO field, one entry per allele and feature, the entry's fields joined by `|`. The output's `##INFO` header line names the fields in their order, the layout. Three flags fill four of the fields: `--symbol` the gene symbol, `--biotype` the transcript's biotype, and `--numbers` the exon and intron numbers. `--symbol` also adds two columns to the layout, `SYMBOL_SOURCE` and `HGNC_ID`.

The VCF writer is selected in two ways, by `--vcf` or by `--output_format vcf`. Ensembl VEP groups flags into option sets, in which one flag switches others on, and its `--vcf` switches on `--symbol`, `--biotype` and `--numbers` ([Output formats](../output-formats.md)).

### 1.2 The record

The record has the position and alleles of the record `ps_45573402` of the `GRCh38-pick-summary` golden corpus, under an identifier of its own, and is read against that corpus's cache (`tests/golden/116/GRCh38-pick-summary/json_cache`):

```
21	45573402	utr_and_intron	C	T	.	.	.
```

Position `45573402` lies in exon 1 of 7 of ENST00000859469, in that transcript's 5' untranslated region, and in the one intron of ENST00000763268. It also lies upstream of three more transcripts, within 350 bases of each. Both engines take one flag and no format flag:

```
--output_format vcf
```

### 1.3 What each engine writes

The `##INFO` line that states the layout. Ensembl VEP 116.2:

```
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS">
```

vep-rs:

```
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">
```

The `CSQ` entries of ENST00000859469 and ENST00000763268. Ensembl VEP 116.2:

```
T|5_prime_UTR_variant|MODIFIER||ENSG00000173638|Transcript|ENST00000859469||||||163|||||||-1|
T|intron_variant&non_coding_transcript_variant|MODIFIER||ENSG00000299414|Transcript|ENST00000763268|||||||||||||1|
```

vep-rs:

```
T|5_prime_UTR_variant|MODIFIER|SLC19A1|ENSG00000173638|Transcript|ENST00000859469|protein_coding|1/7||||163|||||||-1||HGNC|HGNC:10937
T|intron_variant&non_coding_transcript_variant|MODIFIER||ENSG00000299414|Transcript|ENST00000763268|lncRNA||1/1||||||||||1|||
```

The record's three other entries differ the same way, and every other field agrees. With `--vcf` in place of `--output_format vcf`, both engines write vep-rs's layout and entries, and so they do under `--output_format vcf --symbol --biotype --numbers`. The spelling `--output_format VCF` writes the same output as `vcf` in each engine.

One flag alone fills its own fields in Ensembl VEP. Under `--output_format vcf --symbol` it writes `SLC19A1` and the two added columns and leaves `BIOTYPE`, `EXON` and `INTRON` empty. Under `--output_format vcf --biotype --numbers` it writes those three and leaves `SYMBOL` empty, in the shorter layout.

`--fields` replaces the layout with the list it gives. With `--fields Allele,Consequence` added the engines agree, and with `--fields Allele,Feature,SYMBOL` Ensembl VEP leaves `SYMBOL` empty where vep-rs writes `SLC19A1`.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP reads the command line into a table of flags and then applies its option sets. The set that switches `symbol`, `biotype` and `numbers` on belongs to the `vcf` flag, and it also sets `output_format` to `vcf`. A run that writes `--output_format vcf` itself sets the same format and never fires the set, so the three flags stay off.

The output factory then lower-cases the format and builds the VCF writer. The writer's layout always holds the four fields, since they are among its fixed columns, and it appends `SYMBOL_SOURCE` and `HGNC_ID` only when `symbol` is on. Each entry's symbol, biotype and numbers are computed only under their flags, so each of the four fields is written empty.

Source (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`):

- `Config.pm`, lines 500-508: the set of the `vcf` flag.
- `OutputFactory/VCF.pm`, lines 100-118 and 434-447: the fixed columns and the columns the flags add.
- `OutputFactory.pm`, lines 1481, 1520 and 1580: the three tests that gate the values.

### 1.5 Why vep-rs's row is the correct one

The record's position is inside exon 1 of ENST00000859469 and inside the intron of ENST00000763268, and the cache names the first transcript's gene SLC19A1 and gives both transcripts a biotype. Those are the values vep-rs writes, and Ensembl VEP writes them itself for this record under `--vcf`. The layout promises the four fields in every VCF output. A reader of the `--output_format vcf` output cannot tell an empty `SYMBOL` that means no symbol, as on the lncRNA's entry, from one that means the flag was off. vep-rs writes one VCF output for one record, however the format was chosen.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.vcf --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --output_format vcf
```

vep-rs writes the entries of 1.3 with their fields. Ensembl VEP 116.2 takes the same flag with the cache options of [Reproducing these examples](reproducing.md) and writes them empty; with `--vcf` it writes vep-rs's entries.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/output_format_vcf_fields_left_empty.toml` runs the record of 1.2 on both engines and records the record's whole data line. The class covers the four fields and the two columns in a VCF output selected with `--output_format vcf`. The golden corpora and the format comparator select the VCF output with `--vcf`, under which the engines agree, and no committed flag set passes `--output_format vcf`. So no comparator codes the class.
