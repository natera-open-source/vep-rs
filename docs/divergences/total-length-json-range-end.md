# The end key of a position range in JSON under `--total_length`

**The defect in one sentence.** Under `--total_length`, Ensembl VEP 116.3's JSON output writes the first position of a cDNA, coding or protein position range under both the start key and the end key, and writes neither key for a range whose first position is unknown. vep-rs writes the range's last position, with its total, under the end key.

**Why this is a defect.** `cdna_end`, `cds_end` and `protein_end` name the last transcript position a variant covers. For a deletion over cDNA bases 484 to 487, Ensembl VEP's JSON says the variant ends at base 484, and the object holds neither 487 nor the total the flag asks for. For the same record and flag the text formats print `484-487/2674`, and the JSON format without the flag writes `cdna_end` 487. The flag's suffix on the column's text is all that differs between the two JSON runs.

**What vep-rs does instead.** vep-rs splits the column where Ensembl VEP does and keeps the second part. The start key holds the first position as a number, as in Ensembl VEP: that value is correct, and only the end key differs. The end key holds the last position with its total as a string, `"487/2674"`, the form both engines give both keys of a single position under the flag. A range whose first position is unknown (`?-435/2674`) writes the end key alone, as `?-435` does without the flag. A range whose last position is unknown writes the start key alone in both engines.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `OutputFactory.pm`, lines 1737, 1747 and 1751: under `total_length` the cDNA, coding and peptide totals are appended to the text of the position column, after its range, as in `484-487/2674`. Lines 2147, 2153 and 2157 do the same on a structural variant's row.
- `OutputFactory/JSON.pm`, line 335: the column is split at its hyphen into `$s` and `$e`, here `484` and `487/2674`.
- `OutputFactory/JSON.pm`, line 337: the end key takes `$e` only when it is all digits or holds a `?`, else `$s`. An end that carries a total is neither.
- `OutputFactory/JSON.pm`, lines 340-341: a key is deleted unless its value is a number, or a number, a slash and a number. The pattern admits `487/2674`. A `?` that line 337 copied from the start is deleted under both keys.

**What correcting it would change.** With the test of line 337 admitting a number followed by a slash and a total, as the pattern of lines 340-341 does, Ensembl VEP would write `487/2674` under the end key and the two engines would agree.

**The class it produces.** `total_length_json_range_end`: the `cdna_end`, `cds_end` or `protein_end` of a JSON consequence object under `--total_length`, where vep-rs's value is the last position and the total. The class has two forms. In the first, both engines write one number under the start key, Ensembl VEP's end key repeats it and vep-rs's last position is above it. In the second, neither engine writes the start key and Ensembl VEP writes no end key. The golden harness documents each value in the corpus manifest, and the format comparator accepts the two forms as documented, under the class's name.

**How often the shape occurs.** Tier 2: the shape is a record with an allele that spans two or more positions: an insertion, or a change of two or more reference bases. 1,759,673 of the 15,646,481 records of the cuts carry it, 112,464 per million. The transcript decides whether the row prints a range. Of the 511 committed flag sets, 1 passes `--total_length` and compares the JSON format, and 142 of the 152 corpora hold a record of the shape.

## 1. A deletion inside an exon and a deletion across its first base

**In short.** The record `rs11184139` of the `GRCh37-output-flags` corpus deletes four bases of an exon of ENST00000298232, a reverse-strand transcript of TPTE: cDNA bases 484 to 487 of 2,674. Under `--total_length --json`, Ensembl VEP 116.3 writes `cdna_start` 484 and `cdna_end` 484. vep-rs writes `cdna_start` 484 and `cdna_end` `"487/2674"`. The coding and protein keys differ the same way. A second, made-up deletion crosses the exon's first base, so its column reads `?-435/2674`: Ensembl VEP writes no position key for it, and vep-rs writes `cdna_end` `"435/2674"`.

### 1.1 What a reader needs to know first

The text formats print a variant's place on a transcript in three columns, `cDNA_position`, `CDS_position` and `Protein_position`. A column holds one position (`490`) or a range (`484-487`), with `?` for a first or last position that lies outside the transcript's exons. `--total_length` appends a slash and the length of the cDNA, of the coding sequence or of the peptide, as in `484-487/2674` ([Default and tab output](../output/vep-and-tab.md)).

The JSON format splits each column into a start key and an end key: `cdna_start` and `cdna_end`, `cds_start` and `cds_end`, `protein_start` and `protein_end`. A single position fills both keys, and a position written `?` gets no key ([JSON output](../output/json.md)).

### 1.2 The records

Two records, read against the cache of the `GRCh37-output-flags` golden corpus (`tests/golden/116/GRCh37-output-flags/json_cache`). The first is a record of the corpus. The second is made up for this page:

```
21	10970008	rs11184139	CCTTT	C	.	.	.
21	10970060	del_over_exon_start	GGACT	G	.	.	.
```

ENST00000298232 lies on the reverse strand. Its sixth exon is `21:10970009-10970062`, cDNA bases 434 to 487. The first record deletes `21:10970009-10970012`, the exon's last four bases. The second deletes `21:10970061-10970064`: two bases of the intron before the exon, then the exon's first two bases. Both engines take one flag and the JSON format:

```
--total_length --json
```

### 1.3 What each engine writes

The position keys of the transcript's object for the first record. Ensembl VEP 116.3:

```
"cdna_start":484,"cdna_end":484,
"cds_start":116,"cds_end":116,
"protein_start":39,"protein_end":39
```

vep-rs:

```
"cdna_start":484,"cdna_end":"487/2674",
"cds_start":116,"cds_end":"119/1602",
"protein_start":39,"protein_end":"40/533"
```

For the second record Ensembl VEP 116.3 writes none of the six keys. vep-rs writes three:

```
"cdna_end":"435/2674","cds_end":"67/1602","protein_end":"23/533"
```

Every other member of the two objects agrees. In the default format the engines agree on both records. The three columns read `484-487/2674`, `116-119/1602` and `39-40/533` for the first record, and `?-435/2674`, `?-67/1602` and `?-23/533` for the second.

Without `--total_length` both engines write `cdna_end` 487, `cds_end` 119 and `protein_end` 40 for the first record. A single position under the flag is one string under both keys in both engines, the position and the total.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP builds the three columns as text and, under the flag, appends the total after the range. Its JSON writer splits each column at the hyphen. For `484-487/2674` the first part is `484` and the second is `487/2674`.

The writer puts the first part under the start key. It puts the second part under the end key only when that part is all digits or holds a `?`, and the first part otherwise. That fallback gives a single position, which has no second part, one value under both keys. A second part with a total fails the test, so the end key takes `484`.

For `?-435/2674` the first part is `?`, and the fallback copies it under the end key. The writer deletes a key whose value is neither a number nor a number, a slash and a number, so both keys go.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `OutputFactory.pm`, lines 1737, 1747 and 1751: the totals appended to the columns.
- `OutputFactory/JSON.pm`, lines 335-337: the split and the test of the end.
- `OutputFactory/JSON.pm`, lines 340-341: the deletion of a key that holds no position.

### 1.5 Why vep-rs's row is the correct one

A transcript's bases are numbered from its 5' end, which on the reverse strand is each exon's highest coordinate. The five exons before the sixth hold 433 bases, so `21:10970062` is cDNA base 434 and `21:10970009` is base 487. The first record's four bases are cDNA 484 to 487. The coding sequence starts at cDNA base 369, which makes them coding bases 116 to 119, in codons 39 and 40.

So the variant ends at cDNA base 487, coding base 119 and residue 40, the numbers vep-rs writes before each slash. Ensembl VEP prints the same range in its text formats and writes the same three numbers under the end keys without the flag. The totals 2,674, 1,602 and 533 are the ones both engines print in the text columns.

The second record's last base is the exon's second: cDNA base 435, coding base 67, in codon 23. That position is known, and both engines print it after the `?` in the text columns. vep-rs writes it under the end key.

### 1.6 Try it yourself

```
vep -i records.vcf -o out.json --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh37 --total_length --json
```

vep-rs writes the keys of 1.3. Ensembl VEP 116.3 takes the same two flags with the cache options of [Reproducing these examples](reproducing.md) and writes its own.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/total_length_json_range_end.toml` runs the two records of 1.2 on both engines and records their `cdna_end` values. Ensembl VEP's list is `484` alone, since it writes no key for the second record. vep-rs's list is `487/2674` and `435/2674`.

The `total-length` flag set of the `GRCh37-output-flags` corpus names the class in its `masks`. Its manifest documents 104 end keys on 48 transcript objects (48 cDNA, 40 coding and 16 protein keys), each with both engines' values, and `vep-golden-corpus classify` writes those entries from the position columns of the default format. The format comparator accepts an end key under the class only in the two forms the class names. Any other difference on a position key stays a difference.
