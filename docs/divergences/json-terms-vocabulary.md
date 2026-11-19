# A JSON output under `--terms display` or `--terms NCBI`

**The defect in one sentence.** Under `--terms display` or `--terms NCBI`, Ensembl VEP 116.2's JSON writer ranks a record's terms by looking each printed term up in a table keyed by SO term. A record with two or more terms, one of them a label of the vocabulary, finds no entry for the label, and the run ends there, with nothing of that record's buffer written. vep-rs ranks each consequence by its SO term and writes the record, its terms printed in the vocabulary asked for.

**Why this is a defect.** `display` and `NCBI` are two of the three values `--terms` accepts, and the default, tab and VCF outputs write a run under either. The JSON writer needs a rank for one member, `most_severe_consequence`. Ensembl VEP's own `--most_severe` finds that rank under any vocabulary: it ranks the consequences themselves and prints the first one's label. The JSON writer ranks the printed text instead, which is an SO term under the default vocabulary alone. The run dies on the first such record, and its output holds only the buffers written before that record's.

**What vep-rs does instead.** vep-rs ranks every consequence by its SO term and prints it in the vocabulary asked for. A JSON output under `display` or `NCBI` carries each entry's `consequence_terms` and the record's `most_severe_consequence` as labels.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`; `ensembl-variation` at release/116):

- `Config.pm`, line 598: `terms` takes `SO`, `display` or `NCBI`, and lines 945-955 refuse any other value.
- `OutputFactory.pm`, lines 1340-1341: an entry's terms are printed through the vocabulary's method, `display_term` or `NCBI_term`, and a consequence without a label there prints its SO term.
- `OutputFactory/JSON.pm`, line 374: the JSON writer collects every entry's printed terms.
- `OutputFactory/JSON.pm`, lines 387-388: it copies `%OVERLAP_CONSEQUENCES` and sorts the collected terms by `$all_cons{$a}->rank`. The table is keyed by SO term (`ensembl-variation`, `Bio/EnsEMBL/Variation/Utils/Constants.pm`, line 653), so a label has no entry there.
- `OutputFactory.pm`, lines 431-447, `summary_only`: under `--most_severe` the consequences are sorted by their own rank, and the first is printed through the vocabulary's method.

**What correcting it would change.** Ranking the consequences before they are printed, as `summary_only` does, would let the JSON writer finish. Its objects would carry the labels vep-rs writes.

**The class it produces.** `json_terms_label_without_rank`: a JSON output under `--terms display` or `--terms NCBI`, without `--most_severe`, whose input holds a record with two or more terms, one of them a label of the vocabulary. The reference stops before that record's buffer, with a non-zero exit status, where the candidate has one object per record. The class is documented only: Ensembl VEP writes no whole JSON output to compare under either vocabulary, and no committed flag set asks for one.

**How often the shape occurs.** Tier 3: the shape is a JSON output under one of the two vocabularies, without `--most_severe`. Of the 392 committed flag sets, 2 pass `--terms` with one of them, `terms-display` and `terms-ncbi` of the `GRCh37-output-flags` golden corpus, and both compare the default, tab and VCF formats alone.

## 1. A record with six terms

**In short.** The record `21 45573402 C T` has five transcript entries carrying six terms, the most severe a `5_prime_UTR_variant`. Under `--terms display --json`, Ensembl VEP 116.2 ends with the message `Can't call method "rank" on an undefined value` and a non-zero exit status, and writes no object. vep-rs writes the record's object, with `"most_severe_consequence":"5PRIME_UTR"`. In the tab format under the same vocabulary both engines print the same labels, and under `--most_severe` both print `5PRIME_UTR`.

### 1.1 What a reader needs to know first

`--terms` chooses the vocabulary a run prints its consequence terms in. `SO`, the default, prints Sequence Ontology terms such as `5_prime_UTR_variant`. `display` prints Ensembl's display labels, such as `5PRIME_UTR`, and `NCBI` prints labels such as `untranslated_5`. A consequence without a label in the chosen vocabulary prints its SO term: under `NCBI`, `upstream_gene_variant` and `intergenic_variant` do.

A JSON object carries each entry's terms in its `consequence_terms` array. It also carries one `most_severe_consequence`, the term of the highest severity among all the record's entries ([JSON output](../output/json.md)). Severity is a rank every consequence has, whatever vocabulary prints it.

### 1.2 The record

The record has the position and alleles of the record `ps_45573402` of the `GRCh38-pick-summary` golden corpus, under an identifier of its own, and is read against that corpus's cache (`tests/golden/116/GRCh38-pick-summary/json_cache`):

```
21	45573402	utr_and_intron	C	T	.	.	.
```

It lies in the 5' untranslated region of ENST00000859469, in the one intron of the non-coding ENST00000763268, and upstream of three more transcripts. Both engines take the JSON format and one flag:

```
--terms display
```

### 1.3 What each engine writes

Ensembl VEP 116.2 leaves its output file empty and writes this line to standard error:

```
Can't call method "rank" on an undefined value at /opt/vep/src/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory/JSON.pm line 388, <__ANONIO__> line 3.
```

Run with `--fork 1`, as the reproduction runs it, the line `Died in forked process` and a process number follow, and the exit status is 255. Without `--fork` the line stands alone and the status is another non-zero value.

vep-rs writes the record's object:

```
{"allele_string":"C/T","assembly_name":"GRCh38","end":45573402,"id":"utr_and_intron","input":"21\t45573402\tutr_and_intron\tC\tT\t.\t.\t.","most_severe_consequence":"5PRIME_UTR","seq_region_name":"21","start":45573402,"strand":1,"transcript_consequences":[{"consequence_terms":["UPSTREAM"],"gene_id":"ENSG00000173638","impact":"MODIFIER","strand":-1,"transcript_id":"ENST00000651099","variant_allele":"T"},{"consequence_terms":["UPSTREAM"],"gene_id":"ENSG00000299414","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000763267","variant_allele":"T"},{"consequence_terms":["INTRONIC","WITHIN_NON_CODING_GENE"],"gene_id":"ENSG00000299414","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000763268","variant_allele":"T"},{"consequence_terms":["UPSTREAM"],"gene_id":"ENSG00000299414","impact":"MODIFIER","strand":1,"transcript_id":"ENST00000763269","variant_allele":"T"},{"cdna_end":163,"cdna_start":163,"consequence_terms":["5PRIME_UTR"],"gene_id":"ENSG00000173638","impact":"MODIFIER","strand":-1,"transcript_id":"ENST00000859469","variant_allele":"T"}],"vep_rs_version":"0.6.0"}
```

Under `--terms NCBI` Ensembl VEP exits the same way, and vep-rs writes `"most_severe_consequence":"untranslated_5"`. Without `--terms` both engines write `"most_severe_consequence":"5_prime_UTR_variant"`. In the tab format under `--terms display`, both print `UPSTREAM` on three rows, `INTRONIC,WITHIN_NON_CODING_GENE` on one and `5PRIME_UTR` on one, and with `--most_severe` added both print the one row `5PRIME_UTR`. With `--most_severe` added to the JSON run both engines write `"most_severe_consequence":"5PRIME_UTR"`, since each record then has one term; with `--summary` Ensembl VEP ends as above.

The exit needs neither two different labels nor a transcript, as other records read against the same cache show. `21 21410367 C T`, whose three entries each carry `missense_variant`, and `21 12410040 G A,C,T`, three intergenic entries, end a `display` run the same way. A record with one term is written by both engines: `21 11909772 T C` gives `"most_severe_consequence":"INTERGENIC"`. So is a record whose terms all lack a label in the vocabulary: under `NCBI`, both engines write `21 12410040 G A,C,T` with `intergenic_variant`, and `21 16640433 G A` with the `non_coding_transcript_exon_variant` of its sixteen entries.

A record ahead of such a record is written only when it falls in an earlier buffer. With `21 11909772 T C` placed before the record of 1.2, Ensembl VEP writes nothing at the default buffer size and that one object under `--buffer_size 1`.

### 1.4 Why Ensembl VEP writes that

Each entry's terms reach the JSON writer already printed in the vocabulary asked for. The writer gathers them into one list for the record. To find the most severe it copies the table of all consequences, which is keyed by SO term, and sorts the list by the rank of each term's entry in that table. `UPSTREAM` and `5PRIME_UTR` are not keys of the table, so the lookup of either gives an undefined value, and asking that value for its rank ends the run.

A list of one term is never sorted, since Perl's `sort` compares nothing for a single element, so a record with one term passes. So does a record whose every term printed as its SO term, which the table holds.

Source (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`):

- `OutputFactory.pm`, lines 1340-1341: the terms are printed through the vocabulary's method.
- `OutputFactory/JSON.pm`, line 374: the printed terms are collected.
- `OutputFactory/JSON.pm`, lines 387-388: the table is copied and the terms are sorted by `->rank`.

### 1.5 Why vep-rs's row is the correct one

The record's most severe consequence is the 5' UTR variant on ENST00000859469. Both engines say so in the default vocabulary, `5_prime_UTR_variant`, and both say so under `--terms display --most_severe`, in the tab format and in JSON, `5PRIME_UTR`: the same consequence, printed as its label. vep-rs's JSON object states that label, and each entry's `consequence_terms` are the labels both engines print on the entry's row in the tab format. `--terms display` is a value Ensembl VEP accepts, and it changes how a term is printed, never which consequences a record has. A run that writes the record in three formats has no reason to end on it in the fourth.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.jsonl --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --json --terms display
```

vep-rs writes the object of 1.3. Ensembl VEP 116.2 takes the same two flags with the cache options of [Reproducing these examples](reproducing.md) and ends with the message of 1.3, its output empty.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/json_terms_label_without_rank.toml` runs the record of 1.2 on both engines. It records the `most_severe_consequence` of vep-rs's object, and for Ensembl VEP no value and the exit status 255. The class covers a JSON output under `--terms display` or `--terms NCBI`. Ensembl VEP's output stops before the first record of the shape, so a comparison has no whole reference to read, and the two committed flag sets that pass `--terms` compare the default, tab and VCF formats. No comparator codes the class.
