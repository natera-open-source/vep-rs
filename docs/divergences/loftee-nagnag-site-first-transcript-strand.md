# LoFTEE's `NAGNAG_SITE` flag on an acceptor shared by transcripts on both strands

**The defect in one sentence.** The LoF plugin that Ensembl VEP runs for LoFTEE reads the nine bases around a splice-acceptor variant once per record, on the strand of the first transcript whose acceptor row it tests, and keeps them for the record's other transcripts. A transcript on the other strand that shares the acceptor base is tested in the first transcript's orientation, so its `NAGNAG_SITE` flag follows the order of the transcripts. vep-rs reads the nine bases on each row's own strand.

**Why this is a defect.** `NAGNAG_SITE` says that a second acceptor dinucleotide stands three bases from the transcript's own, where splicing could shift by one codon. That is a statement about one transcript's strand, and the routine orients the bases to the transcript it is called for. The stored copy defeats the orientation for every later transcript of the record. The motif `AG.AG` is not its own reverse complement, so a transcript's flag changes with which acceptor row of the record is tested first, though neither the record nor the transcript changes.

**What vep-rs does instead.** For each `splice_acceptor_variant` row vep-rs reads the nine reference bases around the variant, orients them to the strand of the row's transcript and tests the motif. Nothing is kept from one row to the next, so a row's flag is the same whatever other transcripts the record overlaps.

**Where it is in Ensembl VEP's source** (`LoF.pm` of the LOFTEE repository's `grch38` branch, tag `v1.0.4_GRCh38`; the `master` branch, the plugin for GRCh37, holds the same routine at lines 498-516):

- `LoF.pm`, line 147: the plugin takes the record's variation feature from the row it is given, one object for every transcript row of the record.
- Lines 306-307: on a row whose consequences include `splice_acceptor_variant`, `NAGNAG_SITE` is pushed when `check_nagnag_variant` returns true.
- Lines 489-493, in `check_nagnag_variant`: unless the variation feature already holds a `splice_context_seq_cache`, line 490 reads the nine bases, line 491 reverse-complements them when this row's transcript is on the reverse strand, and line 492 stores them on the variation feature.
- Line 494: every call tests the stored sequence. Lines 497-498: a sequence of nine bases is matched against `AG.AG`, anywhere in it.
- Lines 396-399, in `intron_motif_end`, the neighbouring check: the intron sequence it stores is kept on the transcript, one per intron, so that check reads each transcript's own strand.

**What correcting it would change.** Storing the nine bases as the genome gives them and orienting them at each call, or storing one copy per strand, would test every row on its own strand. Both rows of the record below would then read as vep-rs writes them.

**The class it produces.** `loftee_nagnag_site_first_transcript_strand`: on a `splice_acceptor_variant` row of a record that has a `splice_acceptor_variant` row on the other strand, a `LoF_flags` (JSON `lof_flags`) whose two values differ in the `NAGNAG_SITE` item alone. One engine carries the item and the other does not, with every other item equal and in the same order. The field is set aside and counted once per row. A row whose `LoF`, `LoF_filter` or `LoF_info` also differs, or whose flags differ in another item, stays a difference. The format comparator codes the class, reading the strands from the rows' `STRAND`, and a comparison case names it under `masks`.

**How often the shape occurs.** Tier 3: the shape is a base that is an acceptor base of protein-coding transcripts on both strands. None of the committed caches' 116,165 transcripts has one, and 10 of the 511 committed flag sets run the plugin, each on a corpus without one. No committed corpus row falls under the class. The rows of 1.3 are those of a fabricated model, committed as a test fixture.

## 1. Two acceptors shared by a forward and a reverse transcript

**In short.** In a fabricated model, the forward-strand transcript `FWDT1` and the reverse-strand transcript `REVT1` share two pairs of acceptor bases. A substitution at `toy:400` lies where the motif reads on the forward strand alone, and one at `toy:699` where it reads on the reverse strand alone. Ensembl VEP 116.3 flags both rows at `toy:400` and neither at `toy:699`. vep-rs flags `FWDT1` at `toy:400` and `REVT1` at `toy:699`.

### 1.1 What a reader needs to know first

LoFTEE calls a splice-acceptor variant `HC`, high confidence, unless a filter applies, and adds flags that qualify the call without changing it ([What LoFTEE annotates](../plugins/loftee.md#what-loftee-annotates)). `NAGNAG_SITE` is one: the acceptor `AG` has a second `AG` three bases away, so the spliceosome may use the neighbour and the transcript lose or gain one codon. `NON_CAN_SPLICE` is another: the intron does not begin with `GT` and end with `AG`.

Both are read on the transcript's strand. A transcript on the reverse strand reads the reverse complement of the genome, so its acceptor `AG` stands in the genome as `CT`. Two genes on opposite strands can overlap so that one base ends an intron of each. The two last intron bases of a forward-strand transcript lie before its exon's start, and those of a reverse-strand transcript after its exon's end.

### 1.2 The record

The model is two protein-coding genes on a fabricated sequence, `toy`, of 1,000 bases:

| Transcript | Strand | Exons | Coding sequence |
| --- | --- | --- | --- |
| `FWDT1` | forward | 101-200, 401-500, 701-800 | 150-200, 401-500, 701-750 |
| `REVT1` | reverse | 850-950, 620-698, 320-398 | 850-900, 620-698, 349-398 |

Bases 399-400 end the first intron of `FWDT1` and the second of `REVT1`. Bases 699-700 end the second intron of `FWDT1` and the first of `REVT1`. The two records are one substitution on each pair:

```
toy	400	acceptor_forward_motif	G	C	.	.	.
toy	699	acceptor_reverse_motif	C	G	.	.	.
```

The nine bases `toy:396-404` are `AGCAGGCGG`. They hold the motif at 396-400, and their reverse complement, `CCGCCTGCT`, does not. The nine bases `toy:695-703` are `GGGCCTGCT`. They do not hold the motif, and their reverse complement, `AGCAGGCCC`, holds it.

Both engines take the same plugin flag. `use_gerp_end_trunc:0` selects the rule that reads no GERP file, and `run_splice_predictions:0` switches off the splice-prediction layer, which vep-rs does not run:

```
--plugin LoF,loftee_path:<loftee>,use_gerp_end_trunc:0,run_splice_predictions:0
```

### 1.3 What each engine writes

Each record has two rows. In the default format, with both engines reading the model through `--gff`, Ensembl VEP 116.3 writes:

```
acceptor_forward_motif	toy:400	C	FWDG1	FWDT1	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SOURCE=two_strand_acceptor.gff3.gz;LoF=HC;LoF_flags=NAGNAG_SITE;LoF_info=INTRON_SIZE:200
acceptor_forward_motif	toy:400	C	REVG1	REVT1	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1;SOURCE=two_strand_acceptor.gff3.gz;LoF=HC;LoF_flags=NON_CAN_SPLICE,NAGNAG_SITE;LoF_info=INTRON_SIZE:221
acceptor_reverse_motif	toy:699	G	FWDG1	FWDT1	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SOURCE=two_strand_acceptor.gff3.gz;LoF=HC;LoF_flags=NON_CAN_SPLICE;LoF_info=INTRON_SIZE:200
acceptor_reverse_motif	toy:699	G	REVG1	REVT1	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1;SOURCE=two_strand_acceptor.gff3.gz;LoF=HC;LoF_info=INTRON_SIZE:151
```

and vep-rs:

```
acceptor_forward_motif	toy:400	C	FWDG1	FWDT1	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SOURCE=two_strand_acceptor.gff3.gz;LoF=HC;LoF_flags=NAGNAG_SITE;LoF_info=INTRON_SIZE:200
acceptor_forward_motif	toy:400	C	REVG1	REVT1	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1;SOURCE=two_strand_acceptor.gff3.gz;LoF=HC;LoF_flags=NON_CAN_SPLICE;LoF_info=INTRON_SIZE:221
acceptor_reverse_motif	toy:699	G	FWDG1	FWDT1	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SOURCE=two_strand_acceptor.gff3.gz;LoF=HC;LoF_flags=NON_CAN_SPLICE;LoF_info=INTRON_SIZE:200
acceptor_reverse_motif	toy:699	G	REVG1	REVT1	Transcript	splice_acceptor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1;SOURCE=two_strand_acceptor.gff3.gz;LoF=HC;LoF_flags=NAGNAG_SITE;LoF_info=INTRON_SIZE:151
```

The rows of `FWDT1` are the same from both engines. The rows of `REVT1` differ in `NAGNAG_SITE` alone: Ensembl VEP writes it at `toy:400` and vep-rs at `toy:699`. `NON_CAN_SPLICE` stands on the same two rows in both outputs. Each shared pair of bases is the canonical acceptor of one transcript alone, `FWDT1` at 399-400 and `REVT1` at 699-700. Ensembl VEP's rows are the same under each of eight values of `PERL_HASH_SEED`.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP hands the plugin one row at a time, and every row of a record carries the same variation feature. The plugin's routine for the flag runs on a splice-acceptor row and looks on that object for a stored sequence. On the record's first such row there is none: the routine reads the nine bases around the variant, reverse-complements them if the row's transcript is on the reverse strand, and stores the result on the variation feature. On every later acceptor row of the record it finds the stored sequence and tests that.

The row of `FWDT1` is the first acceptor row of both records. At `toy:400` the stored sequence is the forward reading, which holds the motif, so the row of `REVT1` is flagged too. At `toy:699` the stored sequence is again the forward reading, which does not hold it, so the row of `REVT1` is not flagged.

Source (LOFTEE `grch38` branch): `LoF.pm`, line 147, the variation feature; lines 489-493, the stored sequence and its orientation; line 494, the test of the stored copy; lines 306-307, the flag.

### 1.5 Why vep-rs's row is the correct one

Read each transcript's acceptor on its own strand. `REVT1` reads the reverse complement of the genome. At `toy:400` its nine bases are `CCGCCTGCT`, which hold no `AG` at all. `NAGNAG_SITE` does not describe that transcript there, and vep-rs writes none.

At `toy:699` the nine bases of `REVT1` are `AGCAGGCCC`: the intron's last two bases are `AG`, and a second `AG` stands three bases before them. That is the tandem acceptor the flag names, on the transcript's own strand, and vep-rs writes it.

Ensembl VEP's plugin agrees wherever the transcript's row is the record's first acceptor row. It flags `FWDT1` at `toy:400` and not at `toy:699`, each on that transcript's own reading. Its neighbouring check, `NON_CAN_SPLICE`, keeps its sequence per transcript and reads both transcripts on their own strands, as its rows above show.

### 1.6 Try it yourself

```
vep -i variants.vcf -o out.txt --offline --json_cache <cache> \
  --gff two_strand_acceptor.gff3.gz --fasta two_strand_acceptor.fa \
  --species homo_sapiens --assembly GRCh38 \
  --plugin LoF,use_gerp_end_trunc:0,run_splice_predictions:0
```

The files are those of `crates/vep-cli/tests/fixtures/loftee-acceptor-on-both-strands/`, and `<cache>` is any GRCh38 cache, beside which the model is read: `tests/golden/116/GRCh38-plugin-loftee/json_cache` holds no transcript of `toy`. vep-rs writes its four rows of 1.3. Ensembl VEP 116.3 takes `--gff` and `--fasta` with no cache option, the same `--plugin` flag with `loftee_path:` added, and `--dir_plugins` naming a checkout of the LOFTEE `grch38` branch ([Reproducing the plugin concordance](../plugins/reproducing-concordance.md)). The fixture's `ensembl_vep.txt` holds the rows it wrote.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/loftee_nagnag_site_first_transcript_strand.toml` writes the sequence and the model, builds a cache from them for vep-rs, runs both records on both engines and records the four rows of each. Ensembl VEP reads the model file there and vep-rs the cache, so Ensembl VEP's rows also carry the `SOURCE` field of [`file_cache_writes_no_source_field`](file-cache-writes-no-source-field.md), which the reproduction's compared columns leave out.

The format comparator sets aside the `LoF_flags` field of a row of the class and counts the row once. It reads a record's rows as it reads them for `loftee_first_alt_donor_filter`, which for a default or tab row under an identifier needs the input file, `--input`. On the fixture, a comparison that names the class sets aside the two rows of `REVT1` in each of the default, tab, VCF and JSON formats and reports no difference. Without the class it reports those two rows in each format.
