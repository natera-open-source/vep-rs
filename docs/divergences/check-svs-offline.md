# `--check_svs` in an offline run

**The defect in one sentence.** Ensembl VEP 116.2 reads the known structural variants that `--check_svs` reports from its variation database alone, although Ensembl publishes the same variants with each release as an indexed VCF. It declares `--check_svs` and `--offline` a pair of flags that cannot be combined, and ends an offline run under the flag before it reads a record, with `Can't use --check_svs and --offline together`. Under `--safe`, which skips the check of such pairs, the offline run declares the `SV` field in its header and fills it on no row. vep-rs writes offline what Ensembl VEP writes with a database connection, wherever Ensembl's file holds the database's mapping.

**Why this is a defect.** Ensembl VEP's documentation states the limit: the flag needs database access. The limit comes from where the data is kept and not from the question asked. The names of the structural variants over a record are fixed by the release, and Ensembl's own file of them holds what the database query returns, apart from three kinds of mapping named below. So a run that may not reach a database cannot have the field, on a record whose answer the release's files hold.

The refusal is no sound guard either. `--safe`, the flag that makes a plugin's failure fatal, switches the check of flag pairs off. An offline run under `--check_svs --safe` then exits 0 with the flag on, the field declared and no row carrying it, and it prints no warning.

**What vep-rs does instead.** vep-rs takes `--check_svs` with a file: Ensembl's `homo_sapiens_structural_variations.vcf.gz` of the matching release, with its index beside it. It reads the file offline and writes the `SV` field ([Known variants](../output/known-variants.md#sv)); [Output control](../cli/output-control.md) describes the flag. The flag without a file is an argument syntax error, exit 2.

The file lacks three kinds of mapping the database holds: a variant on a scaffold or patch, a further location of a variant mapped to several locations, and a variant that starts at the first base of its chromosome. On a record one of them reaches, vep-rs's offline row lacks a name the connected run writes ([Known variants](../output/known-variants.md#sv)).

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/` unless a path says otherwise):

- `Config.pm`, line 125: `check_svs` is declared without a value.
- `vep`, the launcher, line 30: an argument left over after the options ends the run with `Unexpected extra command-line parameter(s)`, so a file after the flag is refused.
- `Config.pm`, line 628: `check_svs` is listed as incompatible with `offline`. Lines 959-964 end the run with `Can't use --check_svs and --offline together` while the configuration is checked, unless the run passes `--safe` (line 961).
- `AnnotationSourceAdaptor.pm`, lines 202-206: the structural-variant source is built only from a variation database adaptor, under the comment `this has no cache equivalent`.
- `AnnotationSource/Database/StructuralVariation.pm`, line 87: the source matches each variant it fetched from the database against the buffer's records.

**What correcting it would change.** A structural-variant source read from the release's VCF would let an offline run of Ensembl VEP write the names of that file, which are the names vep-rs writes. They are the names its database run writes wherever the file holds the database's mapping.

**The class it produces.** `check_svs_refused_offline`: the `SV` field, and the JSON `colocated_structural_variants` list, of every record of an offline run under `--check_svs`. The reference, Ensembl VEP's side of a comparison, ends before its output opens, with exit status 2, or under `--safe` writes its rows without the field, where the candidate, vep-rs's side, carries it. No comparator codes the class. The two committed flag sets that pass the flag are compared with a reference Ensembl VEP wrote with a connection to its database.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a run under `--check_svs`, whatever its records. Of the 392 committed flag sets, 2 pass the flag, both of the `GRCh38-checksvs` golden corpus.

## 1. One substitution under six structural variants

**In short.** The record `21 5047000 T A` lies inside six known structural variants of Ensembl release 116. Ensembl VEP 116.2 declares `--check_svs` and `--offline` a pair that cannot be combined: offline it refuses the flag with `Can't use --check_svs and --offline together` and exit status 2, and writes no output. `--safe` skips that check, and the offline run then exits 0 with the `SV` field declared in its header, no value on the row and nothing on standard error. vep-rs, given the flag with Ensembl's structural-variant file, writes `SV=nsv2770691,nsv3109934,nsv3163927,nsv3167562,nsv533679,nsv6146554`, the six names Ensembl VEP writes for the record with a connection to its database.

### 1.1 What a reader needs to know first

`--check_svs` adds the `SV` field: the names of the known structural variants that overlap the record. Ensembl keeps those variants in its variation database. With each release it also publishes them as one indexed VCF, `homo_sapiens_structural_variations.vcf.gz`.

An Ensembl VEP run reads its annotation from a local cache, from a database, or from both. `--offline` forbids every database connection. The cache carries the known short variants, so `--check_existing` works offline, and it carries no structural variants.

### 1.2 The record

The record has the position and alleles of a record of the `GRCh38-checksvs` golden corpus, under an identifier of its own. It is read against that corpus's cache (`tests/golden/116/GRCh38-checksvs/json_cache`):

```
21	5047000	snv_5047000	T	A	.	.	.
```

Each engine takes the flag in its own form. Ensembl VEP's takes no value. vep-rs's names the corpus's slice of Ensembl's file, the 987 records of the file of release 116 that meet the corpus's records:

```
--check_svs
--check_svs tests/golden/116/GRCh38-checksvs/structural_variants.vcf.gz
```

Every run is offline. Ensembl VEP is run without `--safe` and with it, and vep-rs writes the same row either way.

### 1.3 What each engine writes

Ensembl VEP 116.2, run without `--safe`, writes no output file and this line to standard error, and exits with status 2:

```
ERROR: Can't use --check_svs and --offline together
```

Under `--safe` it skips its check of flag pairs. The run exits 0 with nothing on standard error. Its header declares the field, and the record's one row carries no `SV`, since no structural-variant source is built without a database:

```
## SV : IDs of overlapping structural variants
snv_5047000	21:5047000	A	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

vep-rs writes the record's row with the field, under `--safe` and without it:

```
snv_5047000	21:5047000	A	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;SV=nsv2770691,nsv3109934,nsv3163927,nsv3167562,nsv533679,nsv6146554
```

Run with `--cache` and without `--offline`, so that it may connect to Ensembl's public database, Ensembl VEP writes vep-rs's row. The corpus's `check-svs` flag set records the same `SV` value for the corpus record at that position.

Neither engine takes the other's form of the flag. Given a file after `--check_svs`, Ensembl VEP ends with `ERROR: Unexpected extra command-line parameter(s)` and the file's path, offline or not. Given the flag with no file, vep-rs ends with `a value is required for '--check_svs <FILE>'` and exit status 2.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP checks its flags against a table of pairs that cannot be combined before it opens the input. The table pairs `check_svs` with `offline`, so the run ends there. The check sits inside a test of `--safe` and is skipped under that flag. The pairing states a fact about the source: the one structural-variant source Ensembl VEP has asks the variation database for the variants over each buffer, and the cache has no file it could read in the database's place.

Source (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`):

- `Config.pm`, line 628: the pair. Lines 959-964: the check and the message; line 961: the test of `--safe` around it.
- `AnnotationSourceAdaptor.pm`, lines 202-206: the source is built from a database adaptor alone.

### 1.5 Why vep-rs's row is the correct one

The six variants are records of Ensembl's structural-variant file of release 116, and each spans position 5047000 on chromosome 21. `nsv6146554` runs from 5029000 to 5047000 and `nsv533679` from 5010002 to 5130908. The names are what the flag asks for, and Ensembl VEP writes the same six for the same record whenever it can reach its database. Its offline run under `--safe` declares the field and leaves it off the row, which no reading of the record supports. vep-rs reads them from Ensembl's own file of that database's variants, so its offline row is the row of Ensembl VEP's connected run.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 \
  --check_svs <corpus>/structural_variants.vcf.gz
```

vep-rs writes its row of 1.3. Ensembl VEP 116.2 takes `--check_svs` alone, with the cache options of [Reproducing these examples](reproducing.md). With `--offline` it ends with the message of 1.3, and with `--offline --safe` it writes its row without the field. With `--cache` and a connection to Ensembl's public database it writes vep-rs's row.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/check_svs_refused_offline.toml` runs the record of 1.2 on both engines offline under `--safe` and records the `SV` value of each engine's row: none for Ensembl VEP, the six names for vep-rs. Before that run it runs Ensembl VEP twice more and stops unless each run reads as 1.3 states. Without `--safe` that is exit status 2, the message and no output. With it, that is exit status 0, nothing on standard error and the header line of the field.

The two flag sets of the `GRCh38-checksvs` golden corpus pass the flag, and their reference outputs were written by Ensembl VEP with a connection to its database. Compared with them, vep-rs's rows agree on every record of the corpus. The differences in the field's text need inputs that corpus does not hold: [`check_svs_mate_chromosome`](check-svs-mate-chromosome.md) a breakend, and the three kinds of mapping Ensembl's file lacks ([Known variants](../output/known-variants.md#sv)) a record they reach. No comparator codes this class.
