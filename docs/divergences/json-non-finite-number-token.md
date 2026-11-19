# An infinite or not-a-number value written as a bare token in the JSON output

**The defect in one sentence.** In the JSON output format, Ensembl VEP 116.2 turns every value its number test admits into a Perl number before encoding, and the encoder writes an infinite or not-a-number double as the bare token `inf`, `-inf` or `nan`. A plugin or annotation value spelling an infinity (`Inf`, `Infinity`), a not-a-number (`NaN`) or a number past the double's range (`1e400`, or an integer of 310 or more digits) therefore leaves the record's line outside the JSON grammar. vep-rs writes such a value as a JSON string holding the cell's text.

**Why this is a defect.** The JSON grammar has no infinity and no not-a-number: RFC 8259 section 6 permits no numeric value outside its grammar and names those two. The encoder's own documentation says such values "cannot be represented in JSON, and it is an error to pass those in", and advises against numifying a value that merely looks like a number "(including infinities/NaN), without validating". Ensembl VEP numifies on the number test alone.

A reader that follows the grammar refuses the whole line, so every field of the record is lost to it, not the one value. The token also loses the text: `Infinity` and `1e400` both become `inf`, and `-NaN` becomes `nan`. The default, tab and VCF outputs of the same run print the cell as the table wrote it.

**What vep-rs does instead.** vep-rs writes the value as a JSON string holding the cell's text, `"Inf"`, `"NaN"` or `"1e400"`, the text its default, tab and VCF outputs print, and every line it writes is a JSON document. A value whose double is finite is a number in both engines, with the blanks and the forms [the JSON page](../output/json.md#field-details) describes.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2; the `JSON` 4.11 and `JSON::XS` 4.04 modules its image carries):

- `Bio/EnsEMBL/VEP/OutputFactory/JSON.pm`, line 288, `get_all_output_hashes_by_InputBuffer`: `numberify($hash, \%NUMBERIFY_EXEMPT)` runs over every output hash before it is encoded; lines 114-121 exempt the identifier keys alone; lines 195-196 encode the hash with a plain `JSON->new`.
- `Bio/EnsEMBL/VEP/Utils.pm`, lines 221-247, `numberify`: a hash value (line 231) or an array element (line 241) that `looks_like_number` admits is replaced by itself plus zero. Perl's number test admits an infinity and a not-a-number in any letter case with an optional sign, the `1.#INF` and `1.#IND` spellings, and any decimal or exponent form, and the addition makes an infinite double of a number past the double's range.
- `JSON::XS`, the backend `JSON` selects: a floating-point value is written through the C library's `%g` conversion with no check that it is finite. That conversion prints an infinite double `inf` or `-inf` and a not-a-number `nan`, in lower case where Perl's own text for the same values is `Inf`, `-Inf` and `NaN`. The module's documentation of its number mapping (`JSON/XS.pm`, "PERL -> JSON") states the values cannot be represented and are an error to pass in.

**What correcting it would change.** Every line of the JSON output would be a JSON document whatever a plugin's table holds, and the value would need a representation the grammar has: the text as a string, as vep-rs writes it, or `null`.

**The class it produces.** `json_non_finite_number_token`, documented: in the JSON format, a member whose value Ensembl VEP writes as the bare token `inf`, `-inf` or `nan` and vep-rs as a string. The row's plugin or custom-annotation cell spells an infinity, a not-a-number, or a number whose double overflows. No committed table carries such a cell, so no manifest names the class and no comparator applies it; the concordance figures score consequence sets and never see the field. The default, tab and VCF outputs carry no difference.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Eight ClinPred cells on the corpus records, through `--json`

**In short.** A synthetic ClinPred table whose score cells read `Inf`, `-Inf`, `NaN`, `1e400`, `Infinity`, `-1e400`, `0 but true` and ` 0 but true` is run on seven records of the `GRCh38-plugin-clinpred-fixture` corpus, eight record and allele pairs, with `--plugin ClinPred,file=<table> --json`. Ensembl VEP writes `"clinpred":inf` for the `Inf`, `Infinity` and `1e400` cells, `"clinpred":-inf` for `-Inf` and `-1e400`, `"clinpred":nan` for `NaN`, and the five lines carrying them are not JSON. vep-rs writes `"clinpred":"Inf"`, `"clinpred":"-Inf"`, `"clinpred":"NaN"`, `"clinpred":"1e400"`, `"clinpred":"Infinity"` and `"clinpred":"-1e400"`. The two remaining cells agree: `0 but true` is the number `0` on both engines, the one ten-byte spelling Perl's number test admits, and ` 0 but true`, with its leading blank, the string `" 0 but true"` on both.

### 1.1 What a reader needs to know first

In the JSON format both engines write a value that reads as a number as a JSON number, `0.5` for a `0.5` cell and `93` for a `93 ` cell, and anything else as a string. Ensembl VEP decides with Perl's `looks_like_number`, which admits every spelling Perl's own numeric conversion accepts. Besides integers, decimals and exponent forms, those are the words `Inf`, `Infinity` and `NaN` in any letter case with an optional sign, and the ten bytes `0 but true`, Perl's idiom for a true zero. vep-rs reads the same spellings and asks one more question, whether the number is finite.

ClinPred's plugin writes its table's fifth cell as the row's `ClinPred` value without reading it as a number (`ClinPred.pm` lines 127 and 133), so the table decides the spelling the JSON writer meets. Any plugin or custom annotation whose cell reaches the output as text has the same path.

### 1.2 The record

Seven records of the `GRCh38-plugin-clinpred-fixture` corpus (`tests/golden/116/GRCh38-plugin-clinpred-fixture/variants.vcf`), substitutions ClinPred scores, one with two alternate alleles, against a synthetic table written for them whose score column holds the eight spellings:

```
21	30282169	stop_lost_rev	C	A	.	.	.
21	30282172	missense_rev	T	C	.	.	.
21	30282928	missense_rev_two_alts	G	T,A	.	.	.
21	30396087	start_lost_fwd	A	G	.	.	.
21	30396091	missense_fwd	C	T	.	.	.
21	30396101	stop_gained_fwd	C	A	.	.	.
21	30642890	start_lost_fwd_unscored	T	C	.	.	.
```

```
#chr	Start	Ref	Alt	ClinPred_Score
21	30282169	C	A	Inf
21	30282172	T	C	-Inf
21	30282928	G	A	1e400
21	30282928	G	T	NaN
21	30396087	A	G	0 but true
21	30396091	C	T	Infinity
21	30396101	C	A	 0 but true
21	30642890	T	C	-1e400
```

Both engines run `--plugin ClinPred,file=<table> --json` with the corpus's cache; the seventh row's cell opens with one space.

### 1.3 What each engine writes

The `clinpred` member of each transcript entry, by record and allele:

| Record and allele | Cell | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- | --- |
| `stop_lost_rev`, `A` | `Inf` | `"clinpred":inf` | `"clinpred":"Inf"` |
| `missense_rev`, `C` | `-Inf` | `"clinpred":-inf` | `"clinpred":"-Inf"` |
| `missense_rev_two_alts`, `T` | `NaN` | `"clinpred":nan` | `"clinpred":"NaN"` |
| `missense_rev_two_alts`, `A` | `1e400` | `"clinpred":inf` | `"clinpred":"1e400"` |
| `start_lost_fwd`, `G` | `0 but true` | `"clinpred":0` | `"clinpred":0` |
| `missense_fwd`, `T` | `Infinity` | `"clinpred":inf` | `"clinpred":"Infinity"` |
| `stop_gained_fwd`, `A` | ` 0 but true` | `"clinpred":" 0 but true"` | `"clinpred":" 0 but true"` |
| `start_lost_fwd_unscored`, `C` | `-1e400` | `"clinpred":-inf` | `"clinpred":"-1e400"` |

Ensembl VEP's transcript entry for `stop_lost_rev`, with the members in the order it printed them:

```
{"impact":"HIGH","cdna_start":782,"cds_start":764,"amino_acids":"*/L","cds_end":764,"variant_allele":"A","strand":-1,"cdna_end":782,"protein_start":255,"codons":"tGa/tTa","protein_end":255,"transcript_id":"ENST00000340345","consequence_terms":["stop_lost"],"clinpred":inf,"gene_id":"ENSG00000188694"}
```

vep-rs's entry for the same record, its keys sorted:

```
{"amino_acids":"*/L","cdna_end":782,"cdna_start":782,"cds_end":764,"cds_start":764,"clinpred":"Inf","codons":"tGa/tTa","consequence_terms":["stop_lost"],"gene_id":"ENSG00000188694","impact":"HIGH","protein_end":255,"protein_start":255,"strand":-1,"transcript_id":"ENST00000340345","variant_allele":"A"}
```

Every other member of every record agrees between the engines, the version key vep-rs alone writes aside. In the default, tab and VCF outputs the two engines' bodies are byte for byte the same: `ClinPred=Inf`, `ClinPred=-Inf`, `ClinPred=NaN`, `ClinPred=1e400`, `ClinPred=0 but true` and the rest as the table spells them, the VCF `CSQ` entry writing a blank as `_`.

### 1.4 Why Ensembl VEP writes that

Before a record's hash is encoded, `numberify` visits every value and replaces each one `looks_like_number` admits with itself plus zero (`Utils.pm` line 231). Perl's number test admits `Inf`, `-Inf`, `Infinity` and `NaN` in any letter case, so the addition yields the infinite and not-a-number doubles; `1e400` is admitted as an exponent form, and its conversion overflows to the infinite double too. `0 but true` is a special case of the same test, admitted only as those ten bytes, so the leading blank of ` 0 but true` keeps that cell text.

`JSON::XS` then encodes the hash. A string value is quoted; a double is handed to the C library's `%g` conversion, which has no notion of the JSON grammar and prints `inf`, `-inf` or `nan`. The output factory passes a plain `JSON->new` with no option that would refuse or convert such a value, and nothing between the plugin's cell and the printed line reads the text back.

Source (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/OutputFactory/JSON.pm`, lines 114-121: the exempt keys; line 288: the `numberify` call; lines 195-196: the encoder.
- `Bio/EnsEMBL/VEP/Utils.pm`, lines 231 and 241: the `+ 0` on every value `looks_like_number` admits.

### 1.5 Why vep-rs's row is the correct one

The JSON grammar decides. A number in a JSON text is a digit sequence with an optional fraction and exponent, and RFC 8259 section 6 permits no value outside that form, naming infinity and not-a-number as the cases excluded. `inf`, `-inf` and `nan` are not JSON values, and a line holding one is not a JSON text. The cell's text is a fact of the table, so the one value the grammar can carry is that text as a string, which is also what both engines' text formats print. vep-rs writes the string; where the table's text has a finite number, `0 but true` included, it writes the number Ensembl VEP writes.

### 1.6 Try it yourself

Write the table of 1.2 under its header, `bgzip` and `tabix -s 1 -b 2 -e 2` it, and run both engines on the records with the corpus's cache; the commands of [Reproducing these examples](reproducing.md#the-commands) give the engines and their caches. vep-rs:

```
vep -i records.vcf -o out.json --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 \
  --plugin ClinPred,file=clinpred_cells_synthetic.tsv.gz --json
```

Ensembl VEP's lines for five of the seven records are refused by any reader that follows the grammar (`python3 -c 'import json, sys; [json.loads(l) for l in sys.stdin]' < out.json` stops at the first); vep-rs's seven lines all parse. Without `--json` both engines write the cells as the table spells them.

### 1.7 How the comparator treats this class

No comparator applies the class. No committed table carries a cell spelling an infinity, a not-a-number or an overflowing number, so no committed flag set meets the shape and no manifest names the class. The concordance figures score consequence sets and never see a plugin field.
