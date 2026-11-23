# The record a JSON run's transcript-assembly mismatch warning does not name

**The defect in one sentence.** In a `--json` run on a RefSeq or merged cache, Ensembl VEP writes the transcript-assembly mismatch warning without the record's name, as `WARNING: Transcript-assembly mismatch in `. It writes that line once for the run, or once for each child process under `--fork`, however many records draw it. vep-rs writes `WARNING: Transcript-assembly mismatch in <Uploaded_variation>` once for each such record, in every output format.

**Why this is a defect.** The warning tells a user which record has an alternate allele that an edited model reads as its own bases, so that the record's rows can be read with that in mind. In the default, tab and VCF formats Ensembl VEP names the record. In a JSON run the same record on the same model draws a line that names nothing. A warning whose text a process wrote before is dropped, so the second and every later record that process annotates draws no line at all. The interpreter's `Use of uninitialized value` warning beside the line says the name was read from a place that does not hold it.

**What vep-rs does instead.** vep-rs writes the named line for each record that draws the warning, in a JSON run as in the other three formats. The warnings file of one input is the same whatever format the run writes.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1726-1727, `TranscriptVariationAllele_to_output_hash`: the warning's text is `Transcript-assembly mismatch in ` joined to the `Uploaded_variation` of the row's hash.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 347-349, `get_all_output_hashes_by_VariationFeature`: the default, tab and VCF formats build every row of a record on the record's own hash, whose `Uploaded_variation` lines 894-895 set.
- `Bio/EnsEMBL/VEP/OutputFactory/JSON.pm`, line 320, `add_VariationFeatureOverlapAllele_info`: the JSON format builds every row on an empty hash, so the row holds no `Uploaded_variation` when the warning is drawn.
- `Bio/EnsEMBL/VEP/BaseVEP.pm`, lines 809-813, `warning_msg`: a warning whose text the process wrote before is dropped. Every unnamed line has one text, so a process writes it once: the run itself without `--fork`, each child under it.

**What correcting it would change.** The JSON format would hand each row the record's hash, as the other formats do, or the warning would read the name from the record itself. A JSON run would then write one named line a record, the lines its default run writes, and the interpreter warning would go.

**The class it produces.** `json_mismatch_warning_names_no_record`: in the warnings file of a `--json` run, Ensembl VEP's line `WARNING: Transcript-assembly mismatch in `, written once for the run or for each forked child, where vep-rs writes a named line for each record that draws the warning. The warnings mode of the format comparator sets those lines aside under the class, and section 1.7 says what each comparison reads.

**How often the shape occurs.** Tier 2: the shape is a record that draws the warning, in a run that writes JSON from a RefSeq or merged cache. The cuts hold such records: the two of section 1 are records of the gnomAD v4.1 chromosome 21 cut. 19 committed flag sets pass `--refseq` or `--merged` and compare the JSON format. The JSON reference warnings of 6 of them, in 2 corpora, hold the unnamed line: five of the six flag sets of `GRCh38-refseq` and the one flag set of `GRCh37-refseq`. The other 13 draw no warning in any format: one runs under `--use_given_ref`, and no record of the other twelve draws it.

## 1. Two records of the GRCh38 RefSeq dataset under `--json`

**In short.** The records `rs922296917` and `rs1247639677` of the GRCh38 RefSeq dataset (`tests/golden/116/GRCh38-refseq/variants.vcf`) each draw the warning on an edited model. Under `--refseq --json`, Ensembl VEP 116.3 writes one line to the warnings file, `WARNING: Transcript-assembly mismatch in `, with a blank after `in` and no name. vep-rs writes two lines, one naming each record. In the default format both engines write those two named lines. The two JSON objects are the same in both engines, apart from vep-rs's own `vep_rs_version` key.

### 1.1 What a reader needs to know first

A RefSeq cache carries models whose sequence was edited to match the RefSeq record where the assembly differs. A record draws one `Transcript-assembly mismatch in <record>` warning where its alternate allele, the empty sequence for a deletion, equals the bases of such a model's first edit, as [RefSeq and merged caches](../output/refseq-and-merged.md#the-transcript-assembly-mismatch-warning) describes. Each of the two records below is a deletion with two rows on models whose first edit deletes bases, so both sides of that comparison are empty. The warning goes to the warnings file, the file either engine writes beside the output, and to standard error unless the run is `--quiet`.

Ensembl VEP builds each output row as a hash and draws the warning while it fills the hash, so what the line can name depends on what the hash holds at that point. The name it joins is the row's `Uploaded_variation`, which the default format prints as its first column. A JSON object carries the same name as its `id`.

### 1.2 The records

```
21	45505371	rs922296917	CGGCCCTCCGGGCCCCCCTG	C	.	.	.
21	45505388	rs1247639677	CTGGGCCCCCT	C	.	.	.
```

Both engines run `--refseq` against the dataset's RefSeq cache and FASTA, once with `--json` and once in the default format. Ensembl VEP runs without `--fork` here, so that nothing is relayed; section 1.3 gives the forked form.

### 1.3 What each engine writes

Ensembl VEP 116.3's warnings file of the JSON run is one line, which ends with the blank after `in`:

```
WARNING: Transcript-assembly mismatch in 
```

Its standard error carries, four times, the interpreter's `Use of uninitialized value in concatenation (.) or string` at `OutputFactory.pm` line 1727: once for each of the four rows that draw the warning. Its warnings file of the default run, and vep-rs's warnings file of either run:

```
WARNING: Transcript-assembly mismatch in rs922296917
WARNING: Transcript-assembly mismatch in rs1247639677
```

Under `--fork` the parent relays what a child wrote to standard error. The JSON run's file then opens with the first interpreter warning behind the child's process id, [`fork_warning_pid_prefix`](fork-warning-pid-prefix.md), then holds the unnamed line and the child's other interpreter warnings. The dataset's committed `expected/json.warnings.log`, written under `--fork 1`, holds that block. Each child drops the repeats of its own warnings alone. Under `--fork 1 --buffer_size 1` the two records go to two children, and the file holds two such blocks, each with the unnamed line.

### 1.4 Why Ensembl VEP writes that

In the default, tab and VCF formats a record's rows start as copies of the record's hash, which holds `Uploaded_variation`, so the warning names the record. The JSON format keeps the record's fields in the object itself and builds each row on an empty hash. When the row of an edited model is filled, the warning joins a key the hash lacks: Perl reads the missing value as an empty text and warns about it, and the line ends after `in `.

The second record draws the same text. The routine that writes a warning keeps every text its process wrote and returns before writing one it has seen, so the second line is dropped where one process annotates both records. The same routine lets through the two named lines of a default run, since their texts differ.

### 1.5 Why vep-rs's lines are the correct ones

The condition belongs to the record and the model: each record's alternate allele equals the bases of an edited model's first edit, whatever format the rows are written in and however many processes annotate them. Two records meet it, so two lines are due, and each must name its record to be of use. Ensembl VEP writes exactly those two lines for this input in three of its four formats. vep-rs writes them in all four.

### 1.6 Try it yourself

```
vep -i records.vcf -o out.json --json --offline --cache \
  --dir_cache <corpus>/cache_root --cache_version 116 \
  --assembly GRCh38 --refseq --fasta <corpus>/reference.fa.gz
cat out.json_warnings.txt
```

vep-rs prints the two named lines. Ensembl VEP 116.3 prints the one unnamed line, and the two named lines when `--json` is dropped.

### 1.7 How the comparisons treat this class

The reproduction `scripts/divergences/reproductions/json_mismatch_warning_names_no_record.toml` holds the two records and the run, with the lines of each engine's warnings file that carry the warning.

The golden harness, the test that holds vep-rs to the committed Ensembl VEP outputs of each corpus, compares the warnings file of each run's default-format output with the warnings the corpus's manifest lists, the named lines, and reads the warnings file of no JSON run. The committed `expected/json.warnings.log` files keep Ensembl VEP's text as written. A test of the `vep` binary, `a_json_run_names_each_record_that_warns`, holds vep-rs's JSON run of these two records to the two named lines.

The warnings mode of the format comparator (`tools/vep-formats-compare`) applies the class where `--masks` names it, and to this shape alone: the reference's file holds the unnamed line and no named one, and vep-rs's file holds a named one. It sets aside the reference's unnamed lines and vep-rs's named lines and counts them under `masks_applied`. A file of vep-rs's that names no record leaves the reference's line a difference. The class cannot tell which records vep-rs names, since the reference's line names none: the test above holds that.

Where the reference run was forked, the class also sets aside the interpreter warnings of the statement that wrote the unnamed line. It reads that statement from the line directly before the unnamed line, so an interpreter warning of any other statement stays a difference. The first of those warnings stands behind the child's process id, and the class takes it once [`fork_warning_pid_prefix`](fork-warning-pid-prefix.md), named beside it, has removed that prefix.
