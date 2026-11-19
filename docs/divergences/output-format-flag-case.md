# The `--output_format` value read in its given case

**The defect in one sentence.** Ensembl VEP 116.2 lower-cases the `--output_format` value to choose the writer and keeps the spelling as typed for two later readers. So `--output_format VCF` selects the VCF writer and then leaves a gene symbol's whitespace unescaped, and `--output_format JSON` selects the JSON writer and then encodes the `GA4GH_VRS` object into a string. vep-rs reads the value once, lower-cased, so every spelling of a format gives the writer and the rows of the lower-case spelling.

**Why this is a defect.** The value names a writer, and the engine's own test of it says the spelling does not matter: `Vcf`, `VCF` and `vcf` all select the VCF writer. The two readers that then test the spelling describe themselves as tests of the format, one encoding a symbol's spaces for VCF, the other avoiding a double encoding when the output is JSON. Under another spelling of the same word each takes the other format's branch. A VCF output then carries a symbol's whitespace as `_`, the `--no_escape` form, where `vcf` writes `%20`, and a JSON output carries the `GA4GH_VRS` object as a string inside the JSON. Nothing but the case of one word differs.

**What vep-rs does instead.** vep-rs reads the value once, lower-cased, before the option sets and every reader see it, so every spelling of a format gives the writer and the rows of the lower-case spelling.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/Config.pm`, line 81: `output_format` takes any string, and nothing in `Config.pm` tests its value.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 225-238, `new`: a copy of the value is lower-cased (line 229) and looked up in `%FORMAT_MAP` (lines 97-103), which names the writer class; an unknown value is refused.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 150, and `Bio/EnsEMBL/VEP/BaseVEP.pm`, lines 755-764, `add_shortcuts`: the factory copies the parameter as typed into its own `output_format` key.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1527: a symbol's whitespace becomes `%20` only when that key is exactly `vcf` and `--no_escape` is off.
- `Bio/EnsEMBL/VEP/OutputFactory/VCF.pm`, line 403: the VCF writer then turns any whitespace left in a value into `_`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1409: the `GA4GH_VRS` value is encoded into a JSON string unless that key is exactly `json`, a step meant to skip the encoding when the writer is the JSON one.

**What correcting it would change.** The two readers would test the lower-cased value, so every spelling of `vcf` would write `%20` for a symbol's whitespace and every spelling of `json` the `GA4GH_VRS` object, as the lower-case spellings do.

**The class it produces.** `output_format_flag_read_in_its_given_case`: a `SYMBOL` in a VCF output selected with `--output_format` spelled other than `vcf`, whose whitespace Ensembl VEP writes as `_` and vep-rs as `%20`. It also covers a `GA4GH_VRS` in a JSON output selected with `--output_format` spelled other than `json`, which Ensembl VEP writes as a string and vep-rs as an object. The class is documented only: no committed flag set spells the value in another case, so no comparison sets a field aside under it. The one committed `--output_format` run spells `ensembl`. One committed cache, the Arabidopsis TAIR10 dataset's, names a gene with whitespace in its symbol, and Ensembl VEP computes `GA4GH_VRS` in database mode alone.

**How often the shape occurs.** Tier 3: the shape is a spelling of the flag's value that no committed run passes. 1 flag set names the flag, with the value `ensembl`, and 0 corpora carry the shape.

## 1. `--output_format VCF` on a gene whose symbol holds a space

**In short.** The record `4 12007993 chr4_mnv2_12007993_TC_AG TC AG` of the Arabidopsis TAIR10 dataset (`tests/golden/116/TAIR10-arabidopsis/variants.vcf`) overlaps the gene AT4G22890, whose symbol in the cache is `PGR5-LIKE A`. Run with `--symbol --output_format VCF`, Ensembl VEP 116.2 writes `PGR5-LIKE_A` in the `SYMBOL` field of every `CSQ` entry of that gene and vep-rs writes `PGR5-LIKE%20A`. Spelled `vcf`, or given as `--vcf`, both engines write `PGR5-LIKE%20A`; under `--no_escape` both write `PGR5-LIKE_A`. The default, tab and JSON outputs of the record are the same in both engines under either spelling.

### 1.1 What a reader needs to know first

`--output_format` names the output writer; `--vcf` is the shorthand that also switches on `--symbol`, `--biotype` and `--numbers`. In a VCF output every annotation of a record is one `CSQ` entry, its fields joined by `|`, and the writer replaces a value's whitespace by `_` so that the INFO column stays one token. A gene symbol is the one value Ensembl VEP encodes first, its spaces becoming `%20`, unless `--no_escape` is given. vep-rs's `--output_format vcf` fills `BIOTYPE`, `EXON` and `INTRON` as `--vcf` does, where Ensembl VEP's leaves them empty, a documented difference ([The CSQ fields of a VCF output selected with `--output_format vcf`](output-format-vcf-fields.md)) that the rows below carry as well.

### 1.2 The record

```
4	12007993	chr4_mnv2_12007993_TC_AG	TC	AG	.	.	.
```

Both engines run `--symbol --output_format VCF` against the dataset's cache and FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.2, the `CSQ` entry for AT4G22890.1:

```
AG|splice_region_variant&synonymous_variant|LOW|PGR5-LIKE_A|AT4G22890|Transcript|AT4G22890.1||||||484-485|313-314|105|S|TCt/AGt|||1||EntrezGene|
```

vep-rs:

```
AG|splice_region_variant&synonymous_variant|LOW|PGR5-LIKE%20A|AT4G22890|Transcript|AT4G22890.1|protein_coding|3/10||||484-485|313-314|105|S|TCt/AGt|||1||EntrezGene|
```

The record's ten entries differ in the `SYMBOL` field of the five AT4G22890 transcripts the same way, and in the `BIOTYPE` and `EXON` fields vep-rs fills under `--output_format vcf`. Under `--output_format vcf` and under `--vcf` both engines write `PGR5-LIKE%20A`; under `--output_format vcf --no_escape` both write `PGR5-LIKE_A`. The default, tab and JSON outputs are the same in both engines whether the format is spelled `vep`, `VEP`, `tab`, `TAB`, `json` or `Json`.

### 1.4 Why Ensembl VEP writes that

The factory lower-cases a copy of `VCF` to choose the VCF writer and keeps `VCF` in its own `output_format` key. Line 1527 tests that key against `vcf`, a comparison that is case-sensitive, so the `%20` encoding is skipped. The VCF writer then replaces the space with `_` at line 403, as it does for any value under `--no_escape`.

### 1.5 Why vep-rs's row is the correct one

The value selected the VCF writer, so the writer's rule for a symbol applies, and that rule is the one Ensembl VEP applies under `vcf` and under `--vcf`. A flag matched without regard to case for the writer cannot be case-sensitive for the writer's encoding of a field. The `_` form is the output of a different request, `--no_escape`, which neither run made. vep-rs reads the value once and writes the format's encoding whatever the spelling.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.vcf --offline --json_cache <corpus>/json_cache \
  --species arabidopsis_thaliana --assembly TAIR10 \
  --fasta <corpus>/reference.fa.gz --symbol --output_format VCF
```

vep-rs writes `PGR5-LIKE%20A`. Ensembl VEP 116.2 writes `PGR5-LIKE_A`. Spelling the value `vcf` makes the engines agree.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/output_format_flag_read_in_its_given_case.toml` holds this record and run with both engines' `SYMBOL` values. The class covers a `SYMBOL` whose two values differ in `_` against `%20` under `--output_format` spelled other than `vcf`, and a `GA4GH_VRS` written as a string against an object under `--output_format` spelled other than `json`. No committed flag set spells the value in another case, so the class is documented only.

## 2. `--output_format JSON` with `--ga4gh_vrs`

**In short.** Under `--ga4gh_vrs`, every transcript entry of the JSON output carries a `ga4gh_vrs` key. Spelled `--output_format json`, Ensembl VEP 116.2 writes the key's value as a JSON object; spelled `--output_format JSON`, it writes the same object encoded as one JSON string. vep-rs writes the object under either spelling. Ensembl VEP computes the notation in database mode alone and writes no `ga4gh_vrs` key offline, so the shape is reached only there.

### 2.1 What a reader needs to know first

`--ga4gh_vrs` adds the GA4GH Variation Representation of each allele, a nested object with a `location` and a `state`. In the default, tab and VCF outputs the object is written as a JSON string inside the field; in the JSON output it is a value of the record's own JSON. Ensembl VEP's factory encodes the object into a string for every format but the JSON one, which encodes the whole record itself. Ensembl VEP computes the notation through its database, so a `--database` run is needed on its side; vep-rs computes it from the cache and the reference FASTA offline.

### 2.2 The record

```
21	10451779	exon_di_ins_fwd_1	A	ATA	.	.	.
```

The record of the GRCh38 identifiers dataset (`tests/golden/116/GRCh38-identifiers/variants.vcf`). Ensembl VEP runs `--database --db_version 116 --ga4gh_vrs --output_format JSON` and `--output_format json`; vep-rs runs `--offline --ga4gh_vrs` with the dataset's cache and FASTA under both spellings.

### 2.3 What each engine writes

Ensembl VEP 116.2 under `--output_format json`, the first transcript entry's key:

```
"ga4gh_vrs":{"state":{"type":"SequenceState","sequence":"TA"},"type":"Allele","location":{"interval":{"end":10451779,"start":10451779,"type":"SimpleInterval"},"type":"SequenceLocation","sequence_id":"refseq:NC_000021.9"}}
```

and under `--output_format JSON`, the same object as one string:

```
"ga4gh_vrs":"{\"state\":{\"type\":\"SequenceState\",\"sequence\":\"TA\"},\"type\":\"Allele\",\"location\":{\"interval\":{\"end\":10451779,\"start\":\"10451779\",\"type\":\"SimpleInterval\"},\"type\":\"SequenceLocation\",\"sequence_id\":\"refseq:NC_000021.9\"}}"
```

Inside the string the `start` is itself a quoted number: the JSON writer's conversion of number-like text into numbers runs over the record after the string was made and cannot reach inside it. vep-rs writes the object under both spellings, equal to Ensembl VEP's `json` object member for member.

### 2.4 Why Ensembl VEP writes that

The factory's `output_format` key holds `JSON`. Line 1409 encodes the value into a string unless that key is exactly `json`, so the test meant to recognise the JSON writer fails for its own format spelled otherwise. The JSON writer then encodes the string once more as a value of the record.

### 2.5 Why vep-rs's row is the correct one

The step at line 1409 exists to keep the object an object in the JSON output; its own comment says so. The output is the JSON output whatever the spelling that selected it, so the object is the value the format calls for, and a string holding JSON inside JSON is the double encoding the step was written to prevent. vep-rs writes the object under every spelling of the format.

### 2.6 Try it yourself

```
vep -i record.vcf -o out.jsonl --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --fasta <corpus>/reference.fa.gz \
  --ga4gh_vrs --output_format JSON
```

vep-rs writes the `ga4gh_vrs` object. Ensembl VEP 116.2, run with `--database --db_version 116` in place of the cache arguments, writes it as a string. Spelling the value `json` makes the engines agree.

### 2.7 How the comparator treats this class

The reproduction of [section 1.7](#17-how-the-comparator-treats-this-class) reproduces the first shape offline; this one needs Ensembl VEP's database mode and is read from the two JSON outputs. The class covers a `GA4GH_VRS` field written as a string by Ensembl VEP and as an object by vep-rs under `--output_format` spelled other than `json`. No committed flag set spells the value in another case, so the class is documented only.
