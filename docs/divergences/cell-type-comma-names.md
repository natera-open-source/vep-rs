# `--cell_type` names with a comma

A difference in the `CELL_TYPE` field under one flag value, and below it a defect of the empty value; nothing is set aside, since no measured dataset runs the flag. Under a `--cell_type` name that holds a comma, Ensembl VEP prints an empty `CELL_TYPE` on every `RegulatoryFeature` row, and vep-rs prints the feature's activity in that epigenome.

`--cell_type` names the epigenomes whose activity each `RegulatoryFeature` row reports, `<epigenome>:<activity>` per name the feature has a record for ([Regulatory and motif features](../regulatory.md)). The names are the cache's: the epigenome short names of the regulatory build with each space written `_`, `A549_(m,_58_y)` or `common_myeloid_progenitor,_CD34-positive_(27-36_y)` (`AnnotationSource/Database/RegFeat.pm` lines 196 to 203). Of the 112 names of the GRCh38 release 116 build, 77 hold a comma, two of them outside the parentheses. The index's count, 385 of 560 epigenome names, sums the five committed caches that list epigenomes: those 112 names in each of five GRCh38 caches.

Ensembl VEP splits the flag's value on every comma (`Config.pm` lines 760 to 762) and the cache's comma-joined `cell_types` line the same way (`CacheDir.pm` line 174), then tests each piece of the value against the pieces of the line (`AnnotationType/RegFeat.pm` lines 77 to 97). `--cell_type A549_(m,_58_y)` is the two pieces `A549_(m` and `_58_y)`, and the cache's line yields the same two pieces, so the check passes. The feature's activity map is keyed by the whole name, so neither piece is found in it (`OutputFactory.pm` lines 1968 to 1978): every `RegulatoryFeature` row prints an empty `CELL_TYPE`, where the feature does carry an activity for `A549_(m,_58_y)`. A name holding a comma can be asked for and never answered.

vep-rs keeps the cache's names as a list and reads the flag's value against it by longest match ([`--cell_type` on the CLI reference](../cli/output-control.md)). At each position the longest name the rest of the value starts with, followed by a comma or the end, is the next name, and that comma is consumed. `A549_(m,_58_y),BJ_(m)` names two epigenomes, and `CELL_TYPE` prints the feature's activity in each, `A549_(m,_58_y):ACTIVE,BJ_(m):INACTIVE`. A value that starts with no name of the list is refused with Ensembl VEP's text, `Cell type "<piece>" unavailable; available cell types are:` and the list on the next line, naming the piece up to the next comma, as Ensembl VEP names it.

For the names without a comma the two engines agree on every row: `--cell_type GM12878_(f),H1_(m),astrocyte` prints the same three activities on both. In the VCF output each engine writes the list joined by `&` and a comma inside a name as `&` too (`OutputFactory/VCF.pm` line 401), so `A549_(m&_58_y):ACTIVE` is the VCF spelling of the first name on both sides.

The record below lies in the enhancer `ENSR21_C6M55` of chromosome 21, active in `A549_(m,_58_y)`:

```
21	46699255	rs868224642	G	C	.	.	.
```

Under `--tab --regulatory --cell_type A549_(m,_58_y)` the two engines write, on its `RegulatoryFeature` row:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `CELL_TYPE` | empty | `A549_(m,_58_y):ACTIVE` |

A golden corpus run under such a value documents each `RegulatoryFeature` row's `CELL_TYPE` in its manifest's `field_divergences` under the class `cell_type_comma_names`. `tools/vep-formats-compare --masks cell_type_comma_names` then checks those cells against the value vep-rs writes in place of Ensembl VEP's empty one, and every other column as usual.

To reproduce: Ensembl VEP 116.3 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--regulatory --cell_type 'A549_(m,_58_y)'`; vep-rs with the same flags over the converted cache.

## An empty `--cell_type` value

The second difference of the flag is a defect of Ensembl VEP 116.3 that vep-rs does not reproduce, a documented class. This section names the rows the class covers; no comparison tool codes the class, since no dataset or corpus of the repository reaches its shape, as [Intended divergences from Ensembl VEP](../intended-divergences.md) lists it. The class carries its tier under that page's rule: tier 1 is a shape the chromosome 21 cuts carry with no flag, tier 2 one that needs a flag or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts are over the cuts' 15,646,481 records, the 152 golden corpora and their 511 flag sets.

**The defect in one sentence.** Under `--regulatory --cell_type ''` or `--cell_type ','`, Ensembl VEP declares no `CELL_TYPE` field yet prints an empty `CELL_TYPE=` on every `RegulatoryFeature` and `MotifFeature` row of the default format and `"cell_type":[]` on each such JSON entry.

**Why this is a defect.** Ensembl VEP splits the value on commas into a list, and an empty list is the flag unset wherever the engine tests its length. Two places do: the option set that switches `--regulatory` on, and the field list that declares the header line and the tab and VCF columns. The output factory copies the list reference itself and tests that, and a reference to an empty list is true. So the field is filled on every regulatory and motif row and printed in the two formats that have no column list, with no header line describing it.

**What vep-rs does instead.** vep-rs reads a value that names no epigenome as the flag absent in every format, so no `CELL_TYPE` key or `cell_type` member is written.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/Config.pm`, lines 760-762: the comma split into a list. Lines 902-906, `_is_flag_active`: `scalar @{$config->{$flag}}` for a list, `0` for the empty one, so the `cell_type` option set at lines 427-432 switches nothing on.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2525-2537, `flag_fields`: the same length test, so `CELL_TYPE` joins no header or column list. Line 217: `cell_type` among the factory's shortcuts, copied by `Bio/EnsEMBL/VEP/BaseVEP.pm`, lines 755-763, `add_shortcuts`, as the array reference. Lines 1879 and 1927: `$hash->{CELL_TYPE} = $self->get_cell_types($rf) if $self->{cell_type}`, true for the empty reference, filling the field with an empty list on every regulatory and motif row.

**What correcting it would change.** The empty value would leave every row without the field, as vep-rs's rows are and as Ensembl VEP's own tab and VCF rows are.

**The class it produces.** `cell_type_empty_value_key`: on a run whose `--cell_type` value splits to no name, the empty `CELL_TYPE` key Ensembl VEP prints in the default format's `Extra` and its `"cell_type": []` on a JSON `regulatory_feature_consequences` or `motif_feature_consequences` entry, which vep-rs leaves out. The tab and VCF formats carry no column on either engine. No row on the committed datasets: every `--cell_type` value a comparison runs names at least one epigenome (`GM12878_(f),H1_(m),astrocyte`, `A549_(m,_58_y)`, `astrocyte`), so the class applies to 0 rows.

**How often the shape occurs.** Tier 3: the shape is a `--cell_type` value that splits to no name, which none of the 511 flag sets passes; 3 of them run the flag, one of those beside `--regulatory_gff`.

### `rs868224642` in the enhancer `ENSR21_C6M55` under `--cell_type ','`

**In short.** The record `21 46699255 rs868224642 G C` over the GRCh38 chromosome 21 cache under `--regulatory --cell_type ','`. Ensembl VEP's default row ends `IMPACT=MODIFIER;BIOTYPE=enhancer;CELL_TYPE=` and each of its four `MotifFeature` rows ends `CELL_TYPE=`, with no `## CELL_TYPE` header line, and its JSON carries `"cell_type":[]` on each of those entries; vep-rs's rows end `BIOTYPE=enhancer` and `TRANSCRIPTION_FACTORS=NRF1`, with no key. `--cell_type ''` gives the same rows. Under `--tab` and `--vcf` the two engines agree, and every other field agrees in every format.
