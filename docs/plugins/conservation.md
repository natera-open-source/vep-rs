# Conservation: GERP scores from a bigWig

Conservation writes the GERP conservation score of a variant's position, read
from a bigWig of Ensembl Compara's per-base scores. It is a record-level
plugin: one `Conservation` value per input record, printed on every row of the
record, intergenic rows included. This page covers its parameters, the data
file and how the value is formed.

## Parameters and data {#conservation}

| Property      | Value                                                                                                                                                                                          |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh38 from Ensembl's file; any build with a bigWig of your own                                                                                                                                 |
| Matches       | position: the record's span, whatever its alleles                                                                                                                                              |
| Output fields | `Conservation`                                                                                                                                                                                 |
| Data source   | Ensembl Compara conservation scores, [GRCh38](https://ftp.ensembl.org/pub/release-116/compara/conservation_scores/92_mammals.gerp_conservation_score/) (`gerp_conservation_scores.homo_sapiens.GRCh38.bw`) |

The plugin takes the bigWig as its first parameter, by position or as `file=`,
and the word `MAX` as an optional last one. The parameters are all positional
or all `key=value`, so `MAX` goes with the positional spelling:

```bash
GERP=/path/to/gerp_conservation_scores.homo_sapiens.GRCh38.bw
--plugin Conservation,$GERP
--plugin Conservation,$GERP,MAX
--plugin Conservation,file=$GERP
```

Without `MAX` the value is the average score over the bases the record
affects; with it, the largest. The header line names the method:
`## Conservation : The average conservation score for this site`, or
`maximum`. `MAX` is read from the last parameter alone, as Ensembl VEP's
plugin reads it, so `MAX` anywhere but last is refused. That plugin
reads a leading `MAX` as the file or drops the requested maximum without a
word ([a documented divergence](../divergences/plugin-unknown-mode-word-refused.md)). Another positional parameter
after the bigWig is not read, as that plugin reads the first alone, and is
noted in the warnings file.

The bigWig is read as it is downloaded; it needs no index. The GRCh38 file is
about 9.6 GB. Ensembl places no restriction on access to or use of its data, so
the file may be copied and redistributed. Ensembl publishes no GRCh37 bigWig,
so on that build the plugin reads a bigWig of per-base scores you provide.

### Two forms that are not provided

Two forms of Ensembl VEP's plugin are not provided. `--plugin
Conservation,database,...` reads the scores from the Ensembl Compara database,
which vep-rs does not connect to, and is refused by name. A first parameter
that is not a local file (`mammals`, or an FTP address) makes Ensembl VEP's
plugin look the file up on Ensembl's FTP site at run time. vep-rs refuses it
and asks for the path of a local bigWig.

### Two invocations that are refused

Two malformed invocations are refused when the plugin loads. Each is a warning,
or the run's error under `--safe`, and the run goes on without the plugin.
`--plugin Conservation` with no parameter is refused with
`Failed to instantiate plugin Conservation: no bigWig file given; use Conservation,<file.bw> or Conservation,<file.bw>,MAX`.
Ensembl VEP's plugin loads with a warning and writes the header line, its unset
method read as `average`, but never a value. Under `--tab` the
`Conservation` column therefore holds `-` on every row (`Conservation.pm`, lines 77-80,
140 and 153; VEP_plugins at release/116).

A first parameter naming a local file that is not a bigWig is refused with
`Failed to instantiate plugin Conservation: <path> is not a readable bigWig: <reason>`,
the reason being the reader's. Ensembl VEP's plugin opens the file on every
row inside `run` (line 195). The library that reads bigWigs ends the
process there with `<path> is not a big wig file`. The plugin's own
warning for an unopened file (lines 197-200) never gets to run. The run ends
with an error and annotates no row.

## The value

A record's chromosome reaches the file's chromosome under its own name or
with a leading `chr` removed, on either side. An input `chr21` reads the file's
`21`, and an input `21` the file's `chr21`. The bases read are the record's
span, or for an insertion the two bases it falls between. Every base the file
covers within the span contributes the value of its interval once; a base the
file does not cover contributes nothing. A record whose span the file covers
nowhere, or whose chromosome the file lacks, gets no `Conservation` field.

Under `--minimal`, each allele of a multi-allelic record the flag splits reads
its own trimmed span, the one its row's alleles describe, as it does when it is
the record's only alternate. Ensembl VEP reads the record's span, the one its
`Location` prints, which its rejoin restores before the plugin runs
([a documented class](../divergences/minimal-plugin-sees-the-rejoined-record.md)).

| Covered bases | Value                                                      |
| ------------- | ---------------------------------------------------------- |
| one           | its score                                                  |
| several       | their average, or under `MAX` the largest of them          |
| none          | no field                                                   |

The value is printed with three decimals (`-1.250`, `3.080`), as Ensembl VEP's
plugin prints it. A structural variant (a symbolic allele such as `<DEL>`) gets
no value: Ensembl VEP hands the plugin sequence variants alone, and vep-rs does
the same.

## Where `MAX` differs from Ensembl VEP

**In short.** Under `MAX`, Ensembl VEP's Conservation plugin sorts a span's
covered scores as text and reports the last, so a span whose scores are all
negative can report a value below its numeric maximum: `-5.2` sorts after
`-2.17`. vep-rs reports the numeric maximum, so wherever the engines differ
Ensembl VEP's value lies below vep-rs's, which lies below zero, the comparator's
rule. On a set of ClinVar GRCh38 chromosome 21 records of
two or more bases against the full GERP bigWig, every record that differs under
`MAX` is a span of all-negative scores, and the average differs on no row. The
class is `conservation_max_text_order`, and
[Intended divergences from Ensembl VEP](../intended-divergences.md) lists it;
the concordance figures never read the field.

### 1.1 What a reader needs to know first

A span of several bases has one score per covered base; the default reports
their average and `MAX` the largest. GERP scores run from about -12, strongly
accelerated, to about 6, strongly constrained; most bases of a genome score
near or below zero. A text sort orders `-5.2`, `-2.17` and `3.1` by their
characters: the minus sign before a digit, then `2` before `5`. Among negative
numbers with the same count of integer digits the one with the larger magnitude
sorts last. But `-12.2` sorts before `-3.5`, whose `3` follows `1`.

### 1.2 The record

The ClinVar GRCh38 record `21 44293093 1005579 CT AC`, a two-base substitution
on chromosome 21, against the Ensembl 116 GRCh38 cache and the Ensembl Compara
GERP bigWig for GRCh38. The bigWig's values over the two bases are `-5.2` and
`-2.17`. Both engines run `--plugin Conservation,<bigWig>,MAX`.

### 1.3 What each engine writes

The record's row on `ENST00000291582`. Ensembl VEP:

```
1005579	21:44293093-44293094	AC	ENSG00000160224	ENST00000291582	Transcript	missense_variant	1327-1328	1196-1197	399	A/D	gCT/gAC	-	IMPACT=MODERATE;STRAND=1;Conservation=-5.200
```

vep-rs:

```
1005579	21:44293093-44293094	AC	ENSG00000160224	ENST00000291582	Transcript	missense_variant	1327-1328	1196-1197	399	A/D	gCT/gAC	-	IMPACT=MODERATE;STRAND=1;Conservation=-2.170
```

Every other column agrees, as it does on the record's other rows.

### 1.4 Why Ensembl VEP writes that

`Conservation.pm`, lines 242-243, sort the covered values with a bare `sort` and
report the last element. Perl's `sort` without a comparison block compares its elements
as strings. The scores are therefore ordered by their characters, and the last
element is the one that sorts last as text, not the largest number. With `-5.2` and
`-2.17` the text order is `-2.17`, `-5.2`, and `-5.2` is reported.

When the covered scores include a non-negative one the text order and the numeric order
agree on the last element for the magnitudes GERP produces. That is why the
difference shows on all-negative spans alone.

### 1.5 Why vep-rs's row is the correct one

`MAX` asks for the largest score, and `-2.17` is larger than `-5.2`. The
plugin's own average over the same two bases is `-3.685`, which lies between
the two values, as the average of two numbers must. Ensembl VEP's `MAX` lies
below its own average, which no maximum can.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin Conservation,gerp_conservation_scores.homo_sapiens.GRCh38.bw,MAX
```

`record.vcf` holds the record of 1.2. `<cache>` is the directory of the
Ensembl 116 GRCh38 cache, the one holding `info.json` and `transcripts/`
([Cache setup](../cache-setup.md)). The bigWig is the GRCh38 file of the Data
source row.

| Mode | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `MAX` | `Conservation=-5.200` | `Conservation=-2.170` |
| the average (no mode) | `Conservation=-3.685` | `Conservation=-3.685` |

The two per-base values can be read from the bigWig with
`bigWigToBedGraph <bigWig> out.bedgraph -chrom=21 -start=44293092 -end=44293094`.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) applies
`conservation_max_text_order` when a comparison case names it. A
`Conservation` value (the JSON `conservation`) that differs is set aside when
both values are numbers and the reference's lies below the candidate's, which
lies below zero, counted once per row. Any other pair of values stays a
difference. The `GRCh38-plugin-conservation` golden corpus keeps records whose
`MAX` agrees between the engines (every multi-base record's largest score at or
above zero). The class therefore applies to none of its rows, and its
comparison case under `MAX` names no class.

## Where next

- [Plugins](../plugins.md): every Ensembl VEP plugin with its support status
  and how one is enabled.
- [The record a plugin reads for a `--minimal` split allele](../divergences/minimal-plugin-sees-the-rejoined-record.md):
  the class behind the span a split allele reads.
- [Intended divergences from Ensembl VEP](../intended-divergences.md): the
  index of every documented difference between the engines.
