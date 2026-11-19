# LoFTEE on GRCh37: a GERP walk across a whole exon without the conservation database

**The defect in one sentence.** On GRCh37 the LoF plugin that Ensembl VEP runs for LoFTEE reads the GERP total of each whole exon between a stop-gained or frameshift variant and the transcript's end from a table of its conservation database. Run without `conservation_file`, the plugin's default, it holds the word `false` in place of that database, and the lookup dies. The row carries none of the plugin's fields, under a `Plugin 'LoF' went wrong` warning. vep-rs sums the per-base GERP file across each whole exon and writes `LoF` and `LoF_info` on the row.

**Why this is a defect.** The module is written to run without the database. It defaults `conservation_file` to `false`, it skips its PhyloCSF lookup under that value, the basic usage line of its documentation names no such file, and that documentation says what the plugin does when none is given. The GERP walk alone calls the database without that test. So one run writes a call for a stop gain whose walk meets no whole exon, as in a single-exon transcript or where the stop codon's exon alone follows the variant's, and nothing for a stop gain with a whole exon after it. The plugin's own `grch38` branch sums such an exon from its GERP file.

**What vep-rs does instead.** On GRCh37 vep-rs sums the per-base GERP file across each whole exon, as it does across the variant's own exon and the exon holding the stop codon, and writes `LoF` and `LoF_info` on the row. It reads no conservation database, and it refuses a `conservation_file` that names one ([LoFTEE's parameters](../plugins/loftee.md#loftee)).

**Where it is in Ensembl VEP's source** (`LoF.pm` and `gerp_dist.pl` of the LOFTEE repository's `master` branch, the plugin for GRCh37; `ensembl-vep` at release/116.2):

- `LoF.pm`, lines 67-72: a parameter is stored when it splits on `:` into exactly two fields. Line 119 defaults `conservation_file` to `'false'`, line 120 sets `conservation_database` to the string `'false'`, and lines 121-128 replace it with a database handle only when the key holds another string.
- `LoF.pm`, line 258: a stop-gained or frameshift row with a coding position calls `get_gerp_weighted_dist` with the GERP source and `conservation_database`.
- `gerp_dist.pl`, lines 26-66: the walk over the transcript's exons. Lines 46-57 sum the stop codon's exon and the variant's exon from the GERP source, and lines 58-60 pass every other exon to `get_exon_gerp` with the database.
- `gerp_dist.pl`, line 129, in `get_exon_gerp`: `$gerp_db->prepare(...)`, a method call on the string `false`.
- `LoF.pm`, line 288: the PhyloCSF lookup runs only when `conservation_file` is not `false`, the test the walk lacks.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2225-2245, `run_plugins`: a plugin's `run` is called inside an `eval`, and an exception leaves the row without the plugin's fields and writes the warning. `Bio/EnsEMBL/VEP/BaseVEP.pm`, lines 809-813: a warning already written is not written again.

**What correcting it would change.** Summing a whole exon from the GERP source the run already reads, as the walk does for the exons at its two ends, would write the row. With a per-base file the row would be the one vep-rs writes.

**The class it produces.** `loftee_gerp_walk_without_conservation_database`: on a GRCh37 run of the plugin with no conservation database, the stop-gained and frameshift rows whose walk meets a whole exon. The reference row carries no `LoF`, `LoF_filter`, `LoF_flags` or `LoF_info`. The candidate row carries `LoF` and `LoF_info`, with `LoF_filter` and `LoF_flags` where a filter or a flag applies. The class is documented only: no committed flag set runs LoFTEE on GRCh37, so no comparator codes it.

**How often the shape occurs.** Tier 3: the shape is a GRCh37 run of the LoF plugin whose `conservation_file` is unset or the word `false`. Of the 392 committed flag sets, 4 load the plugin, the four the index counts as reaching the entry. All four run GRCh38, where the plugin's `grch38` branch sums every exon from its GERP file, so none carries the shape. On a run of the shape, the rows it reaches are those of a stop-gained or frameshift variant with a whole exon after its own.

## 1. A stop gain in the first of three exons

**In short.** `21 37442543 C>T` on GRCh37 is a stop gain in the first exon of CBR1, on four transcripts. On the three that have three exons the walk crosses the whole second exon. Ensembl VEP 116.2, running the plugin's `master` branch, writes no LoF field on those three rows, and vep-rs writes `LoF=HC` with a five-item `LoF_info`. On the fourth transcript, of two exons, both engines write the same `LoF=HC` row.

### 1.1 What a reader needs to know first

LoFTEE calls a stop-gained or frameshift variant `HC`, high confidence, unless one of its filters applies, and then `LC` ([What LoFTEE annotates](../plugins/loftee.md#what-loftee-annotates)). The `END_TRUNC` filter is for a variant near the end of the coding sequence. To judge it on GRCh37, the plugin walks the transcript's exons from the variant's exon onward and adds two totals: `BP_DIST`, each exon's span in bases, and `GERP_DIST`, each exon's sum of per-base GERP conservation scores. `LoF_info` prints both, with the variant's `PERCENTILE` in the coding sequence, `DIST_FROM_LAST_EXON` and the `50_BP_RULE` verdict.

The walk has three kinds of exon. The variant's own exon counts from the variant to the exon's end. The exon holding the stop codon counts up to the stop codon. Any other exon the walk visits is a whole exon. The plugin sums the first two kinds from the per-base GERP file that `gerp_file` names. It reads a whole exon's total from the `gerp_exons` table of the database that `conservation_file` names, the database that also holds the PhyloCSF data.

### 1.2 The record

The record has the position and alleles of the record `fwd_first_exon_stop_gt101` of the `GRCh37-plugin-nmd` golden corpus, under an identifier of its own. It is read against the cache of the `GRCh37` golden corpus (`tests/golden/116/GRCh37/json_cache`):

```
21	37442543	first_exon_stop	C	T	.	.	.
```

It turns codon 44 of CBR1 into a stop on four coding transcripts. ENST00000290349 has the exons `21:37442239-37442702`, `21:37443248-37443355` and `21:37444744-37445464`, and its stop codon ends at `37445180`. ENST00000399191 and ENST00000530908 have the same second exon and a third of their own. ENST00000439427 has two exons, `21:37442313-37442702` and `21:37443248-37443881`, the second holding its stop codon.

The GERP table has the layout of the plugin's own GRCh37 file, a chromosome, a position and a score, bgzipped and indexed with `tabix`. It holds one row for each base of `21:37442001-37446000`, and of `21:31586001-31589000` for the two records of 1.3 that lie in another gene. Each row scores `1`, so that a GERP total is the number of bases summed. The stretch under this record opens:

```
21	37442001	1
21	37442002	1
21	37442003	1
```

Both engines take the GRCh37 assembly and one plugin flag, with no `conservation_file`. `human_ancestor_fa:false` switches the ancestral-allele filter off, and `run_splice_predictions:0` the splice-prediction layer, which vep-rs does not run:

```
--plugin LoF,loftee_path:<loftee>,gerp_file:GERP_scores.txt.gz,human_ancestor_fa:false,run_splice_predictions:0
```

### 1.3 What each engine writes

The record has nine rows, four of them `stop_gained`. Those four in the default format, from Ensembl VEP 116.2:

```
first_exon_stop	21:37442543	T	ENSG00000159228	ENST00000290349	Transcript	stop_gained	305	130	44	Q/*	Cag/Tag	-	IMPACT=HIGH;STRAND=1
first_exon_stop	21:37442543	T	ENSG00000159228	ENST00000399191	Transcript	stop_gained	211	130	44	Q/*	Cag/Tag	-	IMPACT=HIGH;STRAND=1
first_exon_stop	21:37442543	T	ENSG00000159228	ENST00000439427	Transcript	stop_gained	231	130	44	Q/*	Cag/Tag	-	IMPACT=HIGH;STRAND=1;LoF=HC;LoF_info=GERP_DIST:408,BP_DIST:406,PERCENTILE:0.242085661080074,DIST_FROM_LAST_EXON:159,50_BP_RULE:PASS
first_exon_stop	21:37442543	T	ENSG00000159228	ENST00000530908	Transcript	stop_gained	305	130	44	Q/*	Cag/Tag	-	IMPACT=HIGH;STRAND=1
```

and from vep-rs:

```
first_exon_stop	21:37442543	T	ENSG00000159228	ENST00000290349	Transcript	stop_gained	305	130	44	Q/*	Cag/Tag	-	IMPACT=HIGH;STRAND=1;LoF=HC;LoF_info=GERP_DIST:705,BP_DIST:702,PERCENTILE:0.155875299760192,DIST_FROM_LAST_EXON:266,50_BP_RULE:PASS
first_exon_stop	21:37442543	T	ENSG00000159228	ENST00000399191	Transcript	stop_gained	211	130	44	Q/*	Cag/Tag	-	IMPACT=HIGH;STRAND=1;LoF=HC;LoF_info=GERP_DIST:540,BP_DIST:537,PERCENTILE:0.194319880418535,DIST_FROM_LAST_EXON:266,50_BP_RULE:PASS
first_exon_stop	21:37442543	T	ENSG00000159228	ENST00000439427	Transcript	stop_gained	231	130	44	Q/*	Cag/Tag	-	IMPACT=HIGH;STRAND=1;LoF=HC;LoF_info=GERP_DIST:408,BP_DIST:406,PERCENTILE:0.242085661080074,DIST_FROM_LAST_EXON:159,50_BP_RULE:PASS
first_exon_stop	21:37442543	T	ENSG00000159228	ENST00000530908	Transcript	stop_gained	305	130	44	Q/*	Cag/Tag	-	IMPACT=HIGH;STRAND=1;LoF=HC;LoF_info=GERP_DIST:393,BP_DIST:390,PERCENTILE:0.24904214559387,DIST_FROM_LAST_EXON:266,50_BP_RULE:PASS
```

The row of the two-exon ENST00000439427 is the same from both engines, and so are the record's five other rows. Ensembl VEP exits with status 0 and writes one line to its warnings file, which opens `WARNING: <process id> :` under `--fork 1`, as the reproduction runs it. The line closes on the number of input lines read, 5 for the reproduction's two header lines and three records:

```
WARNING: Plugin 'LoF' went wrong: Can't locate object method "prepare" via package "false" (perhaps you forgot to load "false"?) at /work/extra/loftee_grch37/gerp_dist.pl line 129, <__ANONIO__> line 5.
```

The reproduction's two other records lie in CLDN8: a frameshift, `21 31588080 AC>A`, and a stop gain, `21 31588082 G>T`. ENST00000286809 holds the gene's whole coding sequence in its first exon and has a second, untranslated exon, which the walk visits as a whole exon. Ensembl VEP writes no LoF field on its two rows. vep-rs writes `LoF=HC` on both, with `GERP_DIST:1498` on the frameshift and `GERP_DIST:1499` on the stop gain. On the single-exon ENST00000399899 both engines write the same two `LoF=HC` rows, flagged `SINGLE_EXON`.

Three more CBR1 records against the same cache and table show where the difference starts and stops. A one-base deletion in the first exon, `21 37442513 GT>G`, is a frameshift on the same four transcripts and is written like the stop gain: no LoF field from Ensembl VEP on the three-exon transcripts, `LoF=HC` from vep-rs, and one shared row on ENST00000439427. A stop gain in the second exon, `21 37443313 C>T`, and one in the last exon of ENST00000290349, `21 37444779 A>T`, walk no whole exon. Both engines write every row of those two records alike, a row filtered `END_TRUNC` on ENST00000530908 among them.

The spelling of the default does not matter. `conservation_file:false`, a bare `conservation_file:` and `conservation_file:a:b`, the last two dropped by the plugin's two-field rule, give the rows above on both engines. So does `conservation_file:False`, by a neighbouring line: the plugin opens an empty database of that name, and the lookup ends at line 130 on the table the file lacks.

Given a database, the plugin completes the walk. An SQLite file with an empty `phylocsf_data` table and a `gerp_exons` table holding `108` for the second exon of each three-exon coding transcript, the total of that exon's 108 bases, was passed as `conservation_file`. Ensembl VEP then writes `LoF=HC` on all four rows, with the five items vep-rs writes on each, followed by `PHYLOCSF_TOO_SHORT`, an item of the PhyloCSF lookup that naming a database switches on.

### 1.4 Why Ensembl VEP writes that

At start-up the plugin stores each parameter given as one `key:value` pair, defaults `conservation_file` to the word `false`, and opens a database only when the key holds another string. Otherwise the variable that should hold the database handle keeps the string `false`.

For each stop-gained or frameshift row the plugin walks the transcript's exons from the variant's exon onward. It sums the variant's exon and the stop codon's exon from the GERP file, through `tabix`. For any other exon it asks the database for the stored total. With a string in the handle's place, Perl looks for a `prepare` method in a package named `false`, finds none and raises an exception.

Ensembl VEP calls each plugin on each row inside an `eval`. The exception discards what the plugin had computed for the row, the row is written without the plugin's fields, and the message goes to the warnings file once, since a warning already written is not repeated. A row whose walk never reaches the lookup is untouched, which is why the two-exon transcript of the same record carries its call.

Source (LOFTEE `master` branch; `ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/`):

- `LoF.pm`, lines 119-128: the default, the string and the database handle; line 258: the walk's call.
- `gerp_dist.pl`, lines 58-60: a whole exon goes to `get_exon_gerp`; line 129: the method call on the string.
- `OutputFactory.pm`, lines 2225-2245: the `eval` and the warning; `BaseVEP.pm`, lines 809-813: a warning is written once.

### 1.5 Why vep-rs's row is the correct one

Follow the walk on ENST00000290349. The variant's exon counts from `37442543` to `37442702`, 160 bases. The second exon, `37443248` to `37443355`, is a whole exon of 108. The last exon counts from `37444744` to the stop codon's end at `37445180`, 437. Every base scores `1`, so the GERP total is 160 + 108 + 437 = 705, the `GERP_DIST` vep-rs prints. `BP_DIST` adds each span's end minus its start: 159 + 107 + 436 = 702.

The call follows from the coordinates. `END_TRUNC` applies on GRCh37 when `DIST_FROM_LAST_EXON` is at most 50 and `GERP_DIST` at most 180. The last exon's coding part spans 436, so `DIST_FROM_LAST_EXON` is 702 - 436 = 266, and the filter cannot apply whatever the GERP total. Measured as the plugin measures it, the stop gain lies 266 exon bases before the last exon, and no other filter applies, so the call is `HC`.

Ensembl VEP's own plugin agrees wherever it can finish. It writes that call and the same five items when a database supplies the second exon's total (1.3). It writes the two-exon transcript's row, whose walk needs no stored total, with no database at all. The row vep-rs writes on the three-exon transcripts is the plugin's rule, evaluated from the per-base scores the run already reads.

### 1.6 Try it yourself

```
{ seq 31586001 31589000; seq 37442001 37446000; } \
  | awk -v OFS='\t' '{ print 21, $1, 1 }' | bgzip > GERP_scores.txt.gz
tabix -s 1 -b 2 -e 2 GERP_scores.txt.gz

vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh37 \
  --plugin LoF,gerp_file:GERP_scores.txt.gz,human_ancestor_fa:false,run_splice_predictions:0
```

`<corpus>` is `tests/golden/116/GRCh37`, and `record.vcf` holds the record of 1.2. vep-rs writes the four rows of 1.3. Ensembl VEP 116.2 takes the same `--plugin` flag with `loftee_path:` added, `--dir_plugins` naming a checkout of the LOFTEE `master` branch and the cache options of [Reproducing these examples](reproducing.md). It leaves the three rows without a LoF field.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/loftee_gerp_walk_without_conservation_database.toml` builds the table of 1.2, runs the record and the two CLDN8 records of 1.3 on both engines and records their eight stop-gained and frameshift rows. The class covers the stop-gained and frameshift rows of a GRCh37 LoFTEE run without a conservation database whose walk meets a whole exon. No committed flag set runs LoFTEE on GRCh37, so no comparator codes the class.
