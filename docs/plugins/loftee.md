# Consequence filter: LoFTEE

LoFTEE's data files, command and parameters for a loss-of-function confidence
call on stop-gained, frameshift and essential-splice rows, and its differences
from Ensembl VEP's LoF plugin.

## Parameters and data {#loftee}

| Property      | Value                                                                                                                                                                  |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh38, verified on real data; the GRCh37 bundle is read, but its GERP table carries no licence text, so GRCh37 verdicts are unverified                                |
| Matches       | consequence filter over the transcript model, `--fasta` splice motifs, GERP                                                                                            |
| Output fields | `LoF`, `LoF_filter`, `LoF_flags`, `LoF_info`                                                                                                                           |
| Data source   | [LoFTEE data](https://personal.broadinstitute.org/konradk/loftee_data/), per assembly: `GRCh38/` (GERP bigWig) or `GRCh37/` (GERP tabix TSV), plus human_ancestor.fa |

**Data files.** Both upstream bundles are published under
<https://personal.broadinstitute.org/konradk/loftee_data/>, one directory per
assembly; the files below are what each assembly needs.

| Assembly | Directory | Files                                                                                                | Parameters                           |
| -------- | --------- | ---------------------------------------------------------------------------------------------------- | ------------------------------------ |
| GRCh38   | `GRCh38/` | `human_ancestor.fa.gz` with its `.fai` and `.gzi`; `gerp_conservation_scores.homo_sapiens.GRCh38.bw` | `human_ancestor_fa=`, `gerp_bigwig=` |
| GRCh37   | `GRCh37/` | `human_ancestor.fa.gz` with its `.fai` and `.gzi`; `GERP_scores.final.sorted.txt.gz` with its `.tbi` | `human_ancestor_fa=`, `gerp_tabix=`, `run_splice_predictions=0` |

Every file is used as shipped. vep-rs reads the bgzipped `human_ancestor.fa.gz`
through the `.fai` and `.gzi` beside it. It reads the GERP bigWig on GRCh38 and
the tabix-indexed GERP TSV with its `.tbi` on GRCh37. The PhyloCSF databases of
the same directories (`loftee.sql.gz` on GRCh38, `phylocsf_gerp.sql.gz` on
GRCh37) feed `conservation_file`, which vep-rs refuses, and are not needed.

A plain `human_ancestor.fa` with its own `.fai` works the same (the bundle's
`.gzi` belongs to the compressed file alone):

```bash
gzip -dc human_ancestor.fa.gz > human_ancestor.fa
samtools faidx human_ancestor.fa
```

### A command for each assembly

`$LOFTEE` below is the unpacked bundle for the run's assembly,
`/path/to/loftee/grch38` or `/path/to/loftee/grch37`:

```bash
# GRCh38: upstream ships GERP as a bigWig
--plugin LoFTEE,human_ancestor_fa=$LOFTEE/human_ancestor.fa,\
gerp_bigwig=$LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw

# GRCh37: upstream ships GERP as a tabix per-base TSV instead, and a run
# names run_splice_predictions=0 or is refused
--plugin LoFTEE,human_ancestor_fa=$LOFTEE/human_ancestor.fa,\
gerp_tabix=$LOFTEE/GERP_scores.final.sorted.txt.gz,\
run_splice_predictions=0
```

`--fasta`, the pipeline reference genome, is required for the splice-motif
checks: the `GC_TO_GT_DONOR` filter and the `NON_CAN_SPLICE` and `NAGNAG_SITE`
flags. Without it those checks are skipped.

LoFTEE targets parity with the `grch38` branch defaults of Ensembl VEP's LoF
plugin, with the splice-prediction layer off. On GRCh37 the matching Ensembl VEP
run is the `master` branch with `run_splice_predictions:0`, which vep-rs
requires there ([LoFTEE on GRCh37](loftee-grch37.md)). LoFTEE answers to both
names: `--plugin LoF,...`, the Ensembl VEP plugin's name, and
`--plugin LoFTEE,...`, so an Ensembl VEP argument string loads it unchanged.

### What each parameter means

| Parameter | Read on | Meaning |
| --- | --- | --- |
| `human_ancestor_fa` | both | Indexed FASTA, plain (`.fa` + `.fai`) or bgzipped (`.fa.gz` + `.fai` + `.gzi`, the form the upstream bundle ships). Enables the `ANC_ALLELE` filter (SNP-only). |
| `gerp_bigwig` | both | GERP conservation bigWig, as shipped in the upstream `grch38` bundle. Enables the default GERP-weighted `END_TRUNC` filter. |
| `gerp_tabix` | both | GERP as a tabix-indexed per-base TSV, as shipped in the upstream `master`/GRCh37 bundle (`GERP_scores.final.sorted.txt.gz` + `.tbi`). Same filter, different on-disk format. Alias `gerp_scores`. Alias `gerp_file`, the `master` branch's own key. |
| `gerp` | both | Assembly-agnostic spelling of the GERP source; classifies the file by its extension. |
| `loftee_path` | both | The LoFTEE directory. When no GERP source is named, the GRCh38 default is `gerp_conservation_scores.homo_sapiens.GRCh38.bw` appended to this value as given (`LoF.pm` lines 83 and 130; `/vep/loftee/` when absent), so the value ends with a slash; on GRCh37 the default is the `master` branch's fixed `/vep/loftee/GERP_scores.final.sorted.txt.gz`. An unreadable default is refused at plugin initialization like any GERP path. |
| `use_gerp_end_trunc` | GRCh38 | Default `true`. With no GERP source, `use_gerp_end_trunc=false` selects the 5%-percentile `END_TRUNC` fallback, which needs no data file. Ensembl VEP's plugin reads the word `false` as on and switches the GERP rule off under `0` alone ([section 3](#3-use_gerp_end_truncfalse-keeps-the-gerp-rule)). On GRCh37 the key is logged and not read: the `master` branch never reads it and has no percentile rule (`LoF.pm` lines 79 and 253), so a GRCh37 run needs its GERP file. |
| `gerp_end_trunc_cutoff` | GRCh38 | Default `-58`, the `grch38` branch's. On GRCh37 the key is logged and not read: the `master` branch compares with a fixed `+180` (`LoF.pm` line 268), a cutoff of the opposite sign. |
| `min_intron_size` | both | Default `15`. Tunable, as in Ensembl VEP's plugin. |
| `filter_position` | GRCh38 | Default `0.05`, the percentile rule's cutoff. Tunable, as in Ensembl VEP's plugin. On GRCh37 the key is logged and not read, the `master` branch having no percentile rule. |
| `check_complete_cds` | both | Default `false`. Under `true`, `INCOMPLETE_CDS` filters the transcripts the cache marks `cds_start_NF` or `cds_end_NF`; Ensembl VEP's plugin filters every multi-exon transcript ([section 2](#2-check_complete_cds-filters-every-multi-exon-transcript)). |
| `skip_lof_info` | GRCh37 | Any value but `0` or an empty one, true by Perl's string truth, drops `LoF_info` from every row and leaves `LoF`, `LoF_filter` and `LoF_flags`, as the `master` branch does (`LoF.pm` lines 243 and 339). On GRCh38, whose branch has no such key, the word is logged as unknown and changes nothing. |
| `tabix_path` | GRCh37 | Accepted and unused: the `master` branch names its `tabix` program with it (`LoF.pm` lines 130 and 141), and vep-rs reads the GERP file itself. On GRCh38, whose branch has no such key, the word is logged as unknown. |
| `run_splice_predictions` | both | `0` switches the splice-prediction layer off and is accepted silently, so an Ensembl VEP `LoF` argument string is reused as given. A value true by Perl's string truth asks for the layer, which vep-rs does not implement, and is refused at plugin initialization. An empty value leaves the branch's default, as an absent key does: the layer off on GRCh38 and on on GRCh37, where a run whose key's last value is not `0` is refused the same way ([LoFTEE on GRCh37](loftee-grch37.md)). |
| `conservation_file` | both | `false` in any letter case, or an empty value, leaves the PhyloCSF lookup off, as the plugin's default does. Any other value names the database Ensembl VEP's plugin reads for the PhyloCSF items ([What LoFTEE annotates](#what-loftee-annotates)), which vep-rs does not compute, and is refused at plugin initialization: the plugin is dropped under `LoFTEE: conservation_file=<value> asks for the PhyloCSF items of LoF_flags and LoF_info, which vep-rs does not compute`, and under `--safe` the run stops with it, as for a true `run_splice_predictions`. On GRCh37 the `master` branch also takes whole-exon GERP totals from this database, which vep-rs sums from the GERP file ([LoFTEE on GRCh37](loftee-grch37.md#what-the-master-branch-does-besides)). |
| `apply_all`, `fast_length_calculation`, `debug` | neither | Accepted with a warning naming the key, and not implemented. Under `apply_all` set to any word but `false`, Ensembl VEP's plugin writes `LoF` with an empty value, and its `LoF_info` items, on every protein-coding row (`LoF.pm` lines 131 and 227); vep-rs writes the `LoF` fields on loss-of-function rows alone. |

### How a parameter word is read

Every parameter is a named `key=value`, read at its first `=`, and optional;
LoFTEE's native `key:value` form is accepted alongside. A word holding a colon
is read as `key:value` whatever `=` it holds, as Ensembl VEP's LoF plugin reads
it. A `loftee_path` whose directory name holds `=` is therefore read whole.

A colon ending the word is dropped first, as Ensembl VEP's LoF plugin drops the
empty field it produces. So `key:value:` reads as `key:value` and `key:` is
dropped. A word still holding a second colon (`key:value:x`, `key::value`, a
URL's `://`) is dropped and logged. Each field keeps its blanks: a padded key is
unknown, a padded `0` true, as Ensembl VEP's LoF plugin reads them (`LoF.pm`
lines 67-72).

A URL (`gerp_bigwig:https://...`) is such a word; vep-rs reads no remote data
file ([What vep-rs does not
do](../getting-started/what-vep-rs-does-not-do.md#plugins)).

### Which GERP file is read, and what a failed read does

Name the GERP source with `gerp_bigwig`, `gerp_tabix` or `gerp_file`, or let
`loftee_path` supply Ensembl VEP's default. Each key keeps its last value, as
Ensembl VEP's LoF plugin keeps one entry per key. A value that is the word
`false` in any letter case, or empty, names no file under that key. The same
word `false` switches `human_ancestor_fa` off.

When several keys name a file, the key Ensembl VEP's LoF plugin reads for the
assembly wins: `gerp_bigwig` on GRCh38 and `gerp_file` on GRCh37. The last-named
of the other keys comes next. Ensembl VEP's LoF plugin never reads those other
keys: given one alone, it reads its default, vep-rs the named file.

An unreadable GERP path, named or derived, is a hard error at plugin
initialization rather than a silent degradation. A chromosome the GERP file
lacks, bigWig or tabix TSV, counts 0 toward `GERP_DIST`, as Ensembl VEP's LoF
plugin counts it. One warning per run names the file and the chromosome.

A read failure the bigWig library or the tabix reader returns fails the row,
which keeps no LoF fields, under LoFTEE's `went wrong` warning. Ensembl VEP's
LoF plugin dies the same way on a bigWig path it cannot open ([its `master`
branch differs](loftee-grch37.md#what-the-master-branch-does-besides)). A
compressed bigWig block that cannot be inflated stops the run, as it stops
Ensembl VEP's.

## What LoFTEE annotates

LoFTEE is not a score lookup. It filters the consequences vep-rs has already
assigned against the transcript model, the `--fasta` splice motifs and GERP, and
emits a loss-of-function confidence call.

LoFTEE only annotates `stop_gained`, `frameshift_variant`,
`splice_acceptor_variant` and `splice_donor_variant` on protein-coding
transcripts; everything else gets no `LoF` field.

The MaxEntScan + logistic-regression + de-novo-donor-SVM splice-prediction layer
and the `OS` confidence class, off by default on the `grch38` branch, are not
implemented.

The PhyloCSF items below are not computed. None is a filter, so none changes the
`HC` or `LC` call.

| Item | Field | Written by Ensembl VEP's plugin (`LoF.pm` lines 274-285) |
| --- | --- | --- |
| `ANN_ORF`, `MAX_ORF` | `LoF_info` | On a stop-gained or frameshift row with an exon number whose exon the database holds: the score of the exon's annotated reading frame and the highest score of any. |
| `PHYLOCSF_TOO_SHORT` | `LoF_info` | On such a row whose exon the database lacks. |
| `PHYLOCSF_WEAK`, `PHYLOCSF_UNLIKELY_ORF` | `LoF_flags` | When the annotated frame's score is negative: `PHYLOCSF_UNLIKELY_ORF` if the highest score is positive, `PHYLOCSF_WEAK` otherwise. |

## Where LoFTEE differs from Ensembl VEP's LoF plugin

Five behaviours of the LoF plugin Ensembl VEP 116.2 runs for LoFTEE, the
`grch38` branch, are not reproduced. Each is a documented class whose section
names the rows a comparison sets aside. The format comparator codes those of
sections 1, 4 and 5 and no tool the other two. No committed row falls under any,
as [Intended divergences from Ensembl VEP](../intended-divergences.md) lists
them. Its [How to read the
counts](../intended-divergences.md#how-to-read-the-counts) defines the format
comparator and the comparison case the sections name.

The concordance figures score consequence sets and never read a plugin field, so
nothing is set aside in them. Every row below is in the default output format,
against the GRCh38 LoFTEE golden corpus's cache, reference and data files
(`tests/golden/116/GRCh38-plugin-loftee/`).

Each class carries its tier under the rule of [Intended divergences from Ensembl
VEP](../intended-divergences.md). Tier 1 is a shape the chromosome 21 population
cuts carry with no flag. Tier 2 is one that needs a flag or plugin on ordinary
data. Tier 3 is one ordinary data does not carry. The counts are over the cuts'
15,646,481 records, the 111 golden corpora and their 392 flag sets.

### 1. The GC-to-GT donor filter judged by a record's first alternate allele

**The defect in one sentence.** Ensembl VEP's LoF plugin decides the
`GC_TO_GT_DONOR` filter from the first alternate allele of the record for every
allele, so on a multi-allelic splice-donor record a later allele is filtered, or
spared, according to its sibling's base.

**Why this is a defect.** The filter marks a donor substitution that turns a
`GC` intron start into the canonical `GT`. That change leaves splicing intact
and so is not a loss of function. Whether an allele does that depends on its own
base: `C` to `T` does, `C` to `A` does not. The plugin reads the base from the
record's allele list, where only the first alternate is reached; on a one-allele
record that is the allele being annotated.

**What vep-rs does instead.** vep-rs tests each allele's own base against the
intron start. An allele is filtered when it makes `GT` and stays `HC` when it
does not, whatever its siblings are.

**Where it is in Ensembl VEP's source** (`LoF.pm`, the `grch38` branch;
`ensembl-variation` at release/116):

- `LoF.pm`, line 300, `run`: `GC_TO_GT_DONOR` is pushed for a
  `splice_donor_variant` row when `intron_motif_start_GC_to_GT` holds for the
  transcript variation. Every allele of the record shares that object.
- `LoF.pm`, line 363, `intron_motif_start_GC_to_GT`: `my ($ref, $alt)` takes the
  first two alleles of the transcript variation's list, the reference and the
  first alternate. Line 367 tests `$ref eq 'C'` and `$alt eq 'T'` against the
  intron's first two bases.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlap.pm`, lines 363-372: the
  list is the reference allele followed by the alternate alleles in record
  order.

**What correcting it would change.** Each allele's `LoF` call would follow its
own base, as Ensembl VEP's does when the record is split into one line per
allele.

#### The donor filter's class and how often its shape occurs

**The class it produces.** `loftee_first_alt_donor_filter`: on a
`splice_donor_variant` row of a multi-allelic record, a `LoF_filter` (JSON
`lof_filter`) whose two values differ by the `GC_TO_GT_DONOR` item alone. It is
set aside with the row's `LoF` (`lof`), which reads `LC` beside the item and
`HC` without it, counted once per row. Either engine may be the side with the
item; any other differing item stays a difference.

No row on the committed corpora: the LoFTEE golden corpus holds no multi-allelic
record (0 of 25). The class therefore applies to 0 rows of its 16 compared
formats, and no comparison case loads the plugin.

**How often the shape occurs.** Tier 2: 3,506 of the population cuts' records
(224 per million) are multi-allelic substitutions whose alleles differ in being
`T`, or `A` on the reverse strand; the transcript decides whether a `GC` donor
is reached. 56 corpora carry the shape and 4 flag sets run the plugin.

#### 1.1 A second alternate allele at a GC donor

**In short.** The record `21 46126243 C T,A`, the LoFTEE golden corpus's `C T`
donor record with a second alternate, on the forward-strand transcripts
`ENST00001028357` and `ENST00001028364`, whose intron opens `GC` at 46126242.
Both engines write `LoF=LC;LoF_filter=GC_TO_GT_DONOR` on the `T` rows, which
make the donor `GT`. On the `A` rows, which make it `GA`, Ensembl VEP writes the
same filter, read from the first alternate `T`. vep-rs writes `LoF=HC` with no
filter, the call Ensembl VEP gives `C A` on a line of its own.

Each of the four donor rows carries
`LoF_flags=NON_CAN_SPLICE;LoF_info=INTRON_SIZE:261` on both engines. Written
`C A,T`, Ensembl VEP spares the `T` rows the filter instead, and vep-rs's rows
do not change.

### 2. `check_complete_cds` filters every multi-exon transcript

**The defect in one sentence.** Under `check_complete_cds:true`, Ensembl VEP's
LoF plugin tests whether a transcript's incomplete-CDS attributes are defined,
not whether any exist, and the lookup always returns a list, so every
stop-gained or frameshift row of a multi-exon transcript is filtered
`INCOMPLETE_CDS`.

**Why this is a defect.** The plugin's own documentation describes the parameter
as a check of the Ensembl "Complete CDS" annotation. The filter is therefore
meant for a transcript whose CDS start or end was not found, the `cds_start_NF`
and `cds_end_NF` attributes. The core API's `get_all_Attributes` returns a
reference to the matching attributes, an empty list when there are none, and a
reference is always defined. The test therefore never distinguishes a complete
CDS from an incomplete one.

**What vep-rs does instead.** Under `check_complete_cds=true` vep-rs filters
`INCOMPLETE_CDS` only on a transcript the cache marks `cds_start_NF` or
`cds_end_NF`, and leaves a complete-CDS transcript's call `HC`.

**Where it is in Ensembl VEP's source** (`LoF.pm`, the `grch38` branch;
`ensembl` at release/116):

- `LoF.pm`, lines 268-270: `push(@filters, 'INCOMPLETE_CDS') if
  (check_incomplete_cds($tv))` when `lc($self->{check_complete_cds}) eq
  'true'`; lines 423-431, `check_incomplete_cds`: `defined($start_annotation)
  || defined($end_annotation)` over the two `get_all_Attributes` results.
- `Bio/EnsEMBL/Transcript.pm`, lines 1232-1252, `get_all_Attributes`: a code
  returns `\@results`, the list of matching attributes, empty or not.

**What correcting it would change.** Only the rows of a transcript carrying
`cds_start_NF` or `cds_end_NF` (`FLAGS=cds_start_NF` or `cds_end_NF` on the
row) would be filtered `INCOMPLETE_CDS`, as vep-rs's rows are.

**The class it produces.** `loftee_incomplete_cds_every_row`: under
`check_complete_cds:true`, on a stop-gained or frameshift row whose `FLAGS`
names neither `cds_start_NF` nor `cds_end_NF`, a `LoF_filter` carrying
`INCOMPLETE_CDS` in the reference alone and the `LoF` that follows it, `LC` in
the reference and `HC` in the candidate when no other filter applies. The JSON
`lof_filter` and `lof` are read alike. No committed row: the LoFTEE golden
corpus's four runs pass no `check_complete_cds`, and no comparison case loads
the plugin.

**How often the shape occurs.** Tier 3: the shape is the parameter
`check_complete_cds:true` itself, which none of the 392 flag sets passes; 4 of
them run the plugin. Every multi-exon stop-gained or frameshift row then carries
the filter, so the defect needs no particular record.

#### 2.1 Four stop-gained rows in UMODL1

**In short.** The GRCh38 record `21 42076175 4454616 G T`
(`tests/golden/116/GRCh38-plugin-loftee/variants.vcf`) is a stop gain on four
UMODL1 transcripts, none marked `cds_start_NF` or `cds_end_NF`. Run with
`check_complete_cds:true`, Ensembl VEP 116.2 writes
`LoF=LC;LoF_filter=INCOMPLETE_CDS` on all four; vep-rs writes `LoF=HC` with the
same `LoF_info`, as both engines do without the parameter. Over the corpus's 25
records Ensembl VEP filters 53 rows and vep-rs 7. The 7 are the frameshift rows
of transcripts flagged `cds_end_NF`, where the two engines agree.

### 3. `use_gerp_end_trunc:false` keeps the GERP rule

**The defect in one sentence.** Ensembl VEP's LoF plugin tests
`use_gerp_end_trunc` for Perl truth, so the word `false`, the plugin's own
spelling of an off switch, keeps the GERP-weighted `END_TRUNC` rule on, and only
the literal `0` selects the percentile rule.

**Why this is a defect.** The plugin stores its switches as the strings `true`
and `false`: `use_gerp_end_trunc` defaults to `'true'`, `check_complete_cds` and
`human_ancestor_fa` to `'false'`. Its documentation says a flag "not provided
(or set to 'false')" is off. Its own `check_complete_cds` test compares the
string with `'true'`; the `use_gerp_end_trunc` test asks whether the string is
true in Perl, which every non-empty string but `0` is. The documented off value
is therefore indistinguishable from the default.

**What vep-rs does instead.** On GRCh38 vep-rs reads `use_gerp_end_trunc=false`
(the word in any letter case, or `0`) as off and applies the percentile rule,
with no GERP data needed. Every other value keeps the GERP-weighted rule. On
GRCh37 neither engine reads the key.

**Where it is in Ensembl VEP's source** (`LoF.pm`, the `grch38` branch):

- `LoF.pm`, lines 67-72: every `key:value` parameter stored as the string
  given. Line 79: `use_gerp_end_trunc` defaults to `'true'`. Lines 77 and 268:
  `check_complete_cds` defaults to `'false'` and is tested with `lc(...) eq
  'true'`.
- `LoF.pm`, lines 240-258: `unless ($self->{use_gerp_end_trunc})` selects the
  percentile rule (line 241), else the GERP rule, whose `GERP_DIST`, `BP_DIST`,
  `DIST_FROM_LAST_EXON` and `50_BP_RULE` items and `END_TRUNC` test follow
  (lines 246-254).

**What correcting it would change.** Under `use_gerp_end_trunc:false` every
stop-gained or frameshift row would carry the percentile rule's `END_TRUNC` and
no GERP items in `LoF_info`, as vep-rs's rows do.

**The class it produces.** `loftee_use_gerp_end_trunc_false`: under
`use_gerp_end_trunc:false`, on a stop-gained or frameshift row, a `LoF_info`
whose reference value carries the GERP items the candidate's lacks, and a
`LoF_filter` `END_TRUNC` and `LoF` that follow each engine's own rule. The items
are `GERP_DIST` and `BP_DIST`, with `DIST_FROM_LAST_EXON` and `50_BP_RULE` on a
row with an exon number and the `LoF_flags` flag `NO_EXON_NUMBER` without one. The JSON `lof_info`, `lof_flags`,
`lof_filter` and `lof` are read alike. No committed row: the LoFTEE golden
corpus's four runs pass no `use_gerp_end_trunc`, and no comparison case loads
the plugin.

**How often the shape occurs.** Tier 3: the shape is the parameter
`use_gerp_end_trunc` spelt with an off word Perl reads as true, which none of
the 392 flag sets passes. 4 of them run the plugin.

#### 3.1 A stop gain in the single exon of CLDN8

**In short.** The GRCh38 record `21 30215764 64605 G T`
(`tests/golden/116/GRCh38-plugin-loftee/variants.vcf`) is a stop gain in the
single exon of CLDN8, 23.9 percent into the coding sequence. Run with
`use_gerp_end_trunc:false`, Ensembl VEP 116.2 writes the GERP rule's row
unchanged,
`LoF=LC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.238938053097345,GERP_DIST:-204.02089612931,BP_DIST:516,DIST_FROM_LAST_EXON:-333,50_BP_RULE:FAIL`;
vep-rs writes
`LoF=HC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.238938053097345`, the
percentile being under 0.95. Of the corpus's 83 stop-gained and frameshift rows,
the 76 that carry `LoF` differ this way, three in the `LoF` call. Under
`use_gerp_end_trunc:0` both engines select the percentile rule.

### 4. The ancestral-allele filter withheld from every allele of a multi-allelic record

**The defect in one sentence.** Under `human_ancestor_fa`, Ensembl VEP's LoF
plugin tests the length of the record's whole allele string before it compares
an allele with the ancestral base. On a multi-allelic record of one-base alleles
no allele is therefore filtered `ANC_ALLELE`, whatever base the ancestral
sequence states.

**Why this is a defect.** The filter marks a loss-of-function allele that is the
ancestral state, a property of one allele and one base of the ancestral
sequence. The comparison itself reads the allele being annotated.

Before it, a guard meant to skip indels tests the record's allele string. It
skips a string holding `-`, which every indel carries, and a string longer than
three characters. A multi-base substitution carries the longer string, and so
does `C/T/A`, the string of a two-allele substitution record. The same `C>T`
allele is therefore filtered when its record carries it alone and spared when a
sibling allele shares the line.

**What vep-rs does instead.** vep-rs judges each allele of the record on its
own. A one-base substitution whose base equals the ancestral base at its
position is filtered `ANC_ALLELE`, whatever its siblings are. That is the call
Ensembl VEP gives the allele on a line of its own.

**Where it is in Ensembl VEP's source** (`LoF.pm`, the `grch38` branch;
`ensembl-vep` at release/116.2):

- `LoF.pm`, lines 312-313, `run`: `ANC_ALLELE` is pushed when
  `check_for_ancestral_allele` holds for the allele, under any
  `human_ancestor_fa` but the word `false` in any letter case. vep-rs reads that
  off word the same way.
- `LoF.pm`, line 531, `check_for_ancestral_allele`: `return 0` when the
  variation feature's allele string matches `/-/` or its length is not 3. Lines
  535-544: the allele's own sequence is compared with the one base
  `samtools faidx` reads from the ancestral FASTA at the record's position.
- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, line 342: the variation feature's allele
  string is the REF and every ALT joined by `/`. Every allele of the record
  shares that one string.

**What correcting it would change.** Each allele's `LoF_filter` would carry `ANC_ALLELE` by its own base, as it does when the record is split into one line per allele.

#### The ancestral-allele filter's class and how often its shape occurs

**The class it produces.** `loftee_anc_allele_multiallelic_record`: on a
stop-gained, frameshift or essential-splice row of a multi-allelic record run
with `human_ancestor_fa`, a `LoF_filter` (JSON `lof_filter`) whose two values
differ by the `ANC_ALLELE` item alone. It is set aside with the row's `LoF`
(`lof`), which reads `LC` beside the item and `HC` without it, counted once per
row. vep-rs is the side with the item, and any other differing item stays a
difference.

No row on the committed corpora: the LoFTEE golden corpus holds no multi-allelic
record (0 of 25) and no `ANC_ALLELE` row. No comparison case loads the plugin.

Against the `ensemblorg/ensembl-vep` image, which ships no `samtools`, the
plugin applies the filter to no record at all. A comparison with that image's
output therefore never reaches the class.

**How often the shape occurs.** Tier 2: 3,892 of the population cuts' records
(249 per million) are multi-allelic records of one-base alleles, whose allele
string is longer than three characters. The transcript and the ancestral
sequence decide whether an allele is a loss of function equal to the ancestral
base. 61 corpora carry the shape and 4 flag sets run the plugin.

#### 4.1 A stop gain whose base is the ancestral one, alone and beside a sibling

**In short.** At `21:32314600` the GRCh38 reference reads `C` and the LoFTEE
golden corpus's ancestral FASTA reads `T`. `C>T` is a stop gain (`R/*`,
`Cga/Tga`) on the forward-strand transcript `ENST00000339944`. The one-allele
record `21 32314600 C T` writes `LoF=LC;LoF_filter=END_TRUNC,ANC_ALLELE` on that
row in both engines. The two-allele records `21 32314600 C T,A` and `C A,T`
write the same `T` row with `LoF_filter=END_TRUNC` alone in Ensembl VEP and
`END_TRUNC,ANC_ALLELE` in vep-rs.

`END_TRUNC` keeps the call `LC` on both sides; at a site where `ANC_ALLELE` is
the only filter the call itself differs, `HC` against `LC`. The `A` rows are
`synonymous_variant` on both engines and carry no `LoF` field.

### 5. GERP totals read from the bigWig's zoom levels

**The defect in one sentence.** Ensembl VEP's LoF plugin totals each exon
interval's GERP score through one summary bin of the bigWig, served from a zoom
level whenever the file carries one coarse enough, so on the upstream GRCh38
bigWig `GERP_DIST` approximates the per-base sum.

**Why this is a defect.** The plugin's rule is a GERP-weighted distance: the
per-base GERP scores from the variant to the stop codon, summed over the exon
intervals between them. The sum is compared with the cutoff `-58` for the
`END_TRUNC` filter.

A zoom level stores one total per fixed span of bases. A bin that only partly
overlaps an interval contributes that whole-span total in proportion to the
overlap, not the scores of the bases inside it. The total therefore depends on
where the span boundaries fall, and two files holding the same per-base values
give two different totals. The plugin's expression, `length * binMean`, equals
the per-base sum only when no zoom level intervenes.

**What vep-rs does instead.** vep-rs reads the per-base values of each interval
and sums them, a base the file leaves unscored counted at the mean of the scored
ones. `GERP_DIST` is therefore the exact total whatever zoom levels the file
carries, and `END_TRUNC` follows that total. On a bigWig without zoom levels the
two engines' rows are identical, so the plugin's expression and vep-rs's sum
agree on what the total is. The zoom level is what moves it.

#### Where Ensembl VEP's source reads a zoom level

**Where it is in Ensembl VEP's source** (`gerp_dist.pl` and `LoF.pm`, the
`grch38` branch; `Bio::DB::BigWig` as the `ensemblorg/ensembl-vep` image ships
it):

- `gerp_dist.pl`, lines 70-79, `get_interval_gerp`: one
  `features(-type => 'summary')` call per exon interval. The total is
  `$c->length * binMean($c->score())` over the features returned, one for the
  interval.
- `Bio/DB/BigWig.pm`, lines 1076-1080, `score`: the first element of
  `statistical_summary(1)`, lines 1049-1059, one `bigWigSummaryArrayExtended`
  bin over the interval; lines 878-882, `binMean`: `sumData / validCount`. The
  library serves that bin from the coarsest zoom level whose span is at most
  half the interval's length. It serves the bin from the full data when no level
  qualifies. On a file whose smallest level spans 10 bases a 19-base interval is
  exact and a 20-base one is not. On the upstream file, whose smallest level
  spans 40, the switch is at 80 bases.
- `LoF.pm`, line 246: `get_gerp_weighted_dist` sums the intervals' totals. Line
  247 writes `GERP_DIST`. Line 254 pushes `END_TRUNC` when the distance from the
  last exon is at most 50 and the total is at most `gerp_end_trunc_cutoff`.

**What correcting it would change.** Every `GERP_DIST` written from a bigWig
with zoom levels would equal the per-base sum, as vep-rs's does. A borderline
`END_TRUNC` call would follow the exact total.

#### The GERP total's class and how often its shape occurs

**The class it produces.** `loftee_gerp_zoom_level_summary`: on a stop-gained or
frameshift row run with a `gerp_bigwig` carrying a zoom level, a `LoF_info`
(JSON `lof_info`) whose two values differ in the `GERP_DIST` item alone, set
aside and counted once per row. When the two totals straddle the cutoff with the
distance from the last exon at most 50, the row's `LoF_filter` (`lof_filter`)
and `LoF` (`lof`) are set aside with it: `END_TRUNC` stands beside the lower
total, `LC` beside the filter and `HC` without it. Any other differing item
stays a difference.

No committed row: the LoFTEE golden corpus's bigWig is cut without zoom levels,
so the 197 LoF rows of its default output are byte-identical. The class applies
to 0 rows of its 16 compared formats, and no comparison case loads the plugin.

**How often the shape occurs.** Tier 3: the shape is a GERP bigWig carrying a
zoom level, which none of the 392 flag sets passes. 4 of them run the plugin,
each with the corpus's zoom-free file.

The upstream GRCh38 bundle's `gerp_conservation_scores.homo_sapiens.GRCh38.bw`,
the file [Parameters and data](#loftee) names, carries ten zoom levels, the
smallest spanning 40 bases. Under it every exon interval of 80 bases or more is
read from a zoom level. With the corpus's values rebuilt into a bigWig carrying
zoom levels, Ensembl VEP 116.2's `GERP_DIST` differs from vep-rs's on every
stop-gained and frameshift row with a `LoF` field, 76 of 83, and `LoF`,
`LoF_filter` and `LoF_flags` on none: no record's two totals straddle the
cutoff.

#### 5.1 A stop gain in the single exon of CLDN8

**In short.** The GRCh38 record `21 30215764 64605 G T`
(`tests/golden/116/GRCh38-plugin-loftee/variants.vcf`) is a stop gain in the
single exon of CLDN8 on the reverse strand. The plugin's one GERP interval is
therefore the 517 bases from the stop codon to the variant,
`21:30215248-30215764`. On the corpus bigWig, which carries no zoom level, both
engines write
`LoF=LC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.238938053097345,GERP_DIST:-204.02089612931,BP_DIST:516,DIST_FROM_LAST_EXON:-333,50_BP_RULE:FAIL`.
With the same per-base values rebuilt into a bigWig carrying zoom levels of 10,
40, 160, 640 and 2,560 bases, Ensembl VEP 116.2 writes
`GERP_DIST:-283.796467024088` and vep-rs `GERP_DIST:-204.02089612931`, the row
otherwise unchanged; the reproduction
`scripts/divergences/reproductions/loftee_gerp_zoom_level_summary.toml` rebuilds
the file and runs both engines.

The library serves the 517-base interval from the 160-base level as four bins.
The first bin spans `21:30215146-30215305`, and 58 of its bases lie in the
interval. Those 58 scores sum to `29.62`, and the bin's whole-span total of
`-221.87` enters at 58/160, `-80.43`. The last bin spans `21:30215626-30215785`,
139 of its bases lie in the interval and sum to `-38.28`, and its whole-span
total of `-9.21` enters at 139/160, `-8.01`. The two middle bins lie inside the
interval and enter whole.

vep-rs's `-204.02089612931` is the sum of the 517 per-base values, the figure
the plugin's own expression yields on the zoom-free file. The upstream bundle's
file reads `-282.329269599915` for the same interval, whose per-base values
there also sum to `-204.02089612931`; its bins fall on other boundaries.

### A sixth behaviour, without a class

**The strand of the `NAGNAG_SITE` context.** `LoF.pm`, lines 489-494, reads the
nine bases around a splice-acceptor variant once per record, oriented to the
strand of the first transcript that asks. It keeps that sequence for the
record's other transcripts. A second protein-coding transcript whose acceptor is
the same base on the other strand is therefore tested in the first transcript's
orientation, and `AG.AG` is not its own reverse complement, so its flag follows
the order of the transcripts. vep-rs orients the nine bases to each transcript's
own strand.

The committed caches hold no acceptor base shared by protein-coding transcripts on both strands, so the behaviour has no class.

## Where next

- [LoFTEE on GRCh37](loftee-grch37.md): the key a GRCh37 run passes.
- [Reproducing the plugin concordance](reproducing-concordance.md): how to score
  the plugin against Ensembl VEP.
- [Intended divergences from Ensembl VEP](../intended-divergences.md): the index
  of every class.
