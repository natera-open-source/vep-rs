# A regulatory or motif identifier written as a JSON number

**The defect in one sentence.** In the JSON output, Ensembl VEP 116.2 turns every value Perl's number test admits into a number before encoding, exempting six identifier keys, none of them a regulatory or motif key. A `--regulatory_gff` feature or motif whose identifier is all digits is written as `"regulatory_feature_id":12345` or `"motif_feature_id":424242`, a transcription factor named `123` as the number `123`, and an identifier spelt `007` or `1e3` as `7` or `1000`. vep-rs writes each of them as a JSON string holding the text the file or cache carries.

**Why this is a defect.** An identifier is a name, and the JSON output writes names as strings: the same run writes `"transcript_id":"NR_048535.1"`, and the exempt list exists so that an all-digit transcript, gene or variant identifier stays one. A regulatory or motif identifier is the same kind of value under a key the list omits, and a transcription factor name is a name. The conversion loses text: `007` becomes `7` and `1e3` becomes `1000`, so a reader cannot recover the identifier the file wrote, and a reader that expects a string at an identifier key stops on the number. The default, tab and VCF outputs of the same run print the identifier as written.

**What vep-rs does instead.** vep-rs writes every feature identifier (`transcript_id`, `regulatory_feature_id`, `motif_feature_id`) and every transcription factor name as a JSON string holding the text the file or cache carries. It converts the other members as Ensembl VEP converts them, with the rules [the JSON page](../output/json.md#field-details) describes.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/OutputFactory/JSON.pm`, lines 114-121, `%NUMBERIFY_EXEMPT`: `seq_region_name`, `id`, `gene_id`, `gene_symbol`, `transcript_id` and `SYMBOL`. Line 288: `numberify($hash, \%NUMBERIFY_EXEMPT)` runs over the whole record before it is encoded. Lines 349 and 379-383: the `Feature` cell is renamed `<feature type>_id`, so `regulatory_feature_id` or `motif_feature_id`, keys the exempt list does not name.
- `Bio/EnsEMBL/VEP/Utils.pm`, lines 221-247, `numberify`: a hash value under a key the list does not name (line 231), and every element of an array (line 241, with no exemption at all), that `looks_like_number` admits is replaced by itself plus zero. The test admits digit strings, decimals, exponent forms, signed numbers and the infinity and not-a-number words; the addition drops a leading zero and expands an exponent.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/RegFeat.pm`, lines 575-584, `_record_get_id`: a `--regulatory_gff` record's identifier is its `ID` attribute, else `Name`, else `<chr>:<start>-<end>:<type>`. Lines 432-443: the transcription factor names are the `transcription_factor` attribute split on commas. `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1867-1884 and 1900-1943: the identifier becomes the `Feature` cell and the names the `TRANSCRIPTION_FACTORS` list.

**What correcting it would change.** Exempting `regulatory_feature_id`, `motif_feature_id` and the elements of `transcription_factors` would write them as strings in every run, as vep-rs writes them; the default, tab and VCF outputs would not change. The motif's matrix identifier under `motif_name` is converted by both engines and would stay as it is.

**The class it produces.** `regulatory_identifier_numberified`, documented: in the JSON format, a `regulatory_feature_id`, a `motif_feature_id` or an element of `transcription_factors` that Ensembl VEP writes as a number and vep-rs as a string. The shape is a regulatory feature or motif identifier, or a factor name, that Perl's number test admits. Ensembl's own identifiers (`ENSR...`, `ENSM...`) and factor names never are, so no committed file or cache carries the shape, no manifest names the class, and the concordance figures, which score consequence sets, never see the member.

## 1. Three regulatory features and a motif re-identified with digits

**In short.** A `--regulatory_gff` file holds four features over one deletion of the GRCh38 regulatory corpus: an enhancer identified `12345`, a promoter flanking region identified `007`, an open chromatin region identified `1e3` and the corpus's own CTCF binding site `ENSR21_536SRH`. A second file holds a motif identified `424242` whose factors are `HNF4A,123,RXRA`, beside the corpus's motif `ENSM00000574814`. Under `--json`, Ensembl VEP writes the three identifiers as the numbers `12345`, `7` and `1000`, the motif as `424242` and the factor as `123`; vep-rs writes the strings `"12345"`, `"007"`, `"1e3"`, `"424242"` and `"123"`. The two Ensembl-shaped identifiers are strings in both engines, and the default, tab and VCF outputs agree byte for byte.

### 1.1 What a reader needs to know first

Under `--regulatory_gff`, each regulatory feature or motif overlapping a variant gets a row whose `Feature` cell is the record's identifier, taken from its `ID` attribute. In the JSON format that cell becomes the `regulatory_feature_id` or `motif_feature_id` member of an entry under `regulatory_feature_consequences` or `motif_feature_consequences`, and a motif's factor names become the `transcription_factors` list. Both engines then write values that read as numbers as JSON numbers; what differs is which keys are spared.

### 1.2 The record

The regulatory file, the corpus's enhancer and CTCF binding site over the record from `regulatory_features.gff3.gz`, the enhancer re-identified with digits, and two features added over the same bases, identified by a leading-zero and an exponent spelling, compressed with bgzip and indexed with `tabix -p gff`:

```
21	Ensembl	enhancer	5013566	5014331	.	.	.	ID=12345;color=#faca00
21	Ensembl	promoter_flanking_region	5014300	5014400	.	.	.	ID=007;color=#ff0000
21	Ensembl	open_chromatin_region	5014310	5014350	.	.	.	ID=1e3;color=#00ff00
21	Ensembl	CTCF_binding_site	5014322	5014343	.	-	.	ID=ENSR21_536SRH;color=#40e0d0
```

The motif file, the corpus's motif at 5090718 twice, once re-identified with a digit factor name:

```
21	Ensembl	TF_binding_site	5090718	5090731	7.0773697666047	+	.	ID=424242;binding_matrix_id=ENSPFM0290;transcription_factor=HNF4A,123,RXRA
21	Ensembl	TF_binding_site	5090718	5090731	7.0773697666047	+	.	ID=ENSM00000574814;binding_matrix_id=ENSPFM0290;transcription_factor=HNF4A,NR2F1,RXRA
```

The inputs: the corpus record `del_covers_ctcf`, a 23-base deletion at 21:5014322-5014344, run with `--regulatory_gff file=numid.gff3.gz --json`; and a substitution at 21:5090724, `A` to `T`, run with `--regulatory_gff motifs=nummotif.gff3.gz --json`, both with the corpus's FASTA.

### 1.3 What each engine writes

The identifier members, by record:

| Record | Entry | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- | --- |
| `del_covers_ctcf` | enhancer | `"regulatory_feature_id":12345` | `"regulatory_feature_id":"12345"` |
| `del_covers_ctcf` | promoter flanking region | `"regulatory_feature_id":7` | `"regulatory_feature_id":"007"` |
| `del_covers_ctcf` | open chromatin region | `"regulatory_feature_id":1000` | `"regulatory_feature_id":"1e3"` |
| `del_covers_ctcf` | CTCF binding site | `"regulatory_feature_id":"ENSR21_536SRH"` | `"regulatory_feature_id":"ENSR21_536SRH"` |
| 21:5090724 `A/T` | motif `424242` | `"motif_feature_id":424242`, `"transcription_factors":["HNF4A",123,"RXRA"]` | `"motif_feature_id":"424242"`, `"transcription_factors":["HNF4A","123","RXRA"]` |
| 21:5090724 `A/T` | motif `ENSM00000574814` | `"motif_feature_id":"ENSM00000574814"` | `"motif_feature_id":"ENSM00000574814"` |

Ensembl VEP's entry for the motif `424242`, with its members in the order it printed them:

```
{"motif_name":"ENSPFM0290","impact":"MODIFIER","strand":1,"motif_pos":7,"transcription_factors":["HNF4A",123,"RXRA"],"variant_allele":"T","consequence_terms":["TF_binding_site_variant"],"motif_feature_id":424242}
```

vep-rs's entry for the same motif, its keys sorted:

```
{"consequence_terms":["TF_binding_site_variant"],"impact":"MODIFIER","motif_feature_id":"424242","motif_name":"ENSPFM0290","motif_pos":7,"strand":1,"transcription_factors":["HNF4A","123","RXRA"],"variant_allele":"T"}
```

Every other member of every entry agrees. In the default, tab and VCF outputs the two engines' rows are byte for byte the same: the `Feature` cells read `007`, `12345`, `1e3`, `ENSR21_536SRH` and `424242`, and `TRANSCRIPTION_FACTORS=HNF4A,123,RXRA`.

### 1.4 Why Ensembl VEP writes that

Before a record's hash is encoded, `numberify` visits every member. A hash value whose key the exempt list does not name, and every element of a list, is tested with `looks_like_number` and replaced by itself plus zero when the test admits it. The record's identifier sits under `regulatory_feature_id` or `motif_feature_id`, keys the list does not name, so `12345` becomes the number 12345, `007` the number 7 and `1e3` the number 1000; `ENSR21_536SRH` fails the test and stays a string. The factor names sit in a list, where nothing is exempt, so `123` becomes a number beside `HNF4A` and `RXRA`. The encoder then writes each number bare and each string quoted.

### 1.5 Why vep-rs's row is the correct one

The file wrote `007` as the feature's name; `7` is not that name, and no reader of the number can tell `007` from `7` or `1e3` from `1000`. The identifier is text the file carries, and the one JSON value that carries it whole is the string, which is also what both engines' text formats print and what both write for every transcript identifier. vep-rs writes the string for every feature identifier and every factor name; where a member is a measurement, `motif_pos`, `strand` or `motif_score_change`, it writes the number Ensembl VEP writes.

### 1.6 Try it yourself

Write the two files of 1.2 under a `##gff-version 3` line, `bgzip` and `tabix -p gff` each, and run both engines on the records with the regulatory corpus's cache and FASTA; the commands of [Reproducing these examples](reproducing.md#the-commands) give the engines and their caches. vep-rs:

```
vep -i records.vcf -o out.json --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --fasta <corpus>/reference.fa.gz \
  --regulatory_gff file=numid.gff3.gz,motifs=nummotif.gff3.gz --json
```

`grep -o '"regulatory_feature_id":[^,}]*' out.json` prints the four identifiers quoted from vep-rs and three of them bare from Ensembl VEP 116.2. Without `--json` the two engines' outputs are the same.

### 1.7 How the comparator treats this class

No comparator applies the class. The reproduction `scripts/divergences/reproductions/regulatory_identifier_numberified.toml` holds the regulatory file, the record and both engines' entries. No committed regulatory file or cache carries an identifier or a factor name that reads as a number, so no committed flag set meets the shape and no manifest names the class; the concordance figures score consequence sets and never see the member.
