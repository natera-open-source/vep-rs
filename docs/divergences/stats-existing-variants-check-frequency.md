# The statistics report's count of existing variants under `--check_frequency`

**The defect in one sentence.** Under `--check_frequency`, the `Novel / existing variants` row of Ensembl VEP's statistics report counts every variant the run processes as existing, whether or not the cache holds a known variant at its position. vep-rs counts the variants that have a co-located known variant, the count both engines print for the same records without the flag.

**Why this is a defect.** The row says how many of the processed variants are already known: a variant is existing when a known variant of the cache matches it, and novel otherwise. A frequency check decides which variants a run writes and makes no variant known. Ensembl VEP's frequency routine leaves an empty list of known variants on a variant that had none, and the report counts the list's presence, not its members. So the row reads no novel variant and every variant existing, whatever the input, and disagrees with the `Existing_variation` column of the same run.

**What vep-rs does instead.** vep-rs counts a processed variant as existing when one of its alleles has a co-located known variant, with `--check_frequency` and without it. Its row equals Ensembl VEP's row for the same records under `--check_existing` alone.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm`, lines 316-317, `frequency_check_buffer`: the frequency data of every variant in the buffer is read.
- The same file, line 195, `get_frequency_data`: a `foreach` over the variant's list of known variants. Perl creates the list, empty, on a variant that has none.
- `Bio/EnsEMBL/VEP/Stats.pm`, line 244, `log_VariationFeature`: the count rises when the variant holds the list, and a list is true whether or not it is empty.

**What correcting it would change.** The difference vanishes. `OutputFactory.pm` line 1046 counts the list's members before it writes the known variants. The same count at `Stats.pm` line 244 gives the row of the run without the frequency check.

**The class it produces.** `stats_existing_variants_check_frequency`, the name comparisons count this difference under: the `Novel / existing variants` row of the statistics report of a run under `--check_frequency`, or under `--filter_common`, which switches it on. The concordance figures and the output-format comparisons never read the statistics report, so nothing is set aside in them. Five committed runs of the repository's test corpora, in two corpora, run with the flag.

Set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Two records, one of them known

**In short.** On `21 30107044 . T C`, which has a known variant at its position, and `21 33030010 . A G`, which has none, Ensembl VEP 116.3 prints `0 (0.0) / 2 (100.0)` under `--check_frequency` and `1 (50.0) / 1 (50.0)` under `--check_existing` alone. vep-rs prints `1 (50.0) / 1 (50.0)` in both runs. Every row of the annotation output and every other row of the report agree between the engines.

### 1.1 What a reader needs to know first

`--stats_text` writes a plain-text statistics report, `<output>_summary.txt`, beside the output; the HTML report carries the same table. Its `[General statistics]` section holds the row `Novel / existing variants`, written `novel (percent) / existing (percent)` over the processed variants, or `-` when the run counted no existing variant. `Variants processed` counts the variants the run writes, so a variant a filter drops is in neither figure.

`--check_existing` looks each input variant up among the cache's known variants and prints the matches in the `Existing_variation` column. Ensembl VEP reads the cache from one base before the variant to one base after it. A known variant matches when it shares an allele with the input variant, or when it lists no alleles and has the same start and end (`AnnotationType/Variation.pm`, lines 146-207).

`--check_frequency` switches that lookup on, then drops or keeps alleles by the frequency a known variant records for them; `--filter_common` is a preset of it. A variant with no known variant has no frequency, and an `exclude` filter keeps it.

### 1.2 The records

Both are records of the `GRCh38-freqfilter` test corpus (`tests/golden/116/GRCh38-freqfilter/variants.vcf`), annotated against the corpus's slice of the release-116 GRCh38 cache. The cache's known-variant file holds `rs1254937786` at position 30107044 of chromosome 21, sharing the first record's alternate allele, and the filter below keeps that record.

The file holds no variant from 33030009 to 33030011, the window read for a variant at 33030010 (`VariationTabix.pm`, line 229). The corpus's copy of the file shows it:

```
tabix <corpus>/json_cache/variations/21/all_vars.gz 21:33030009-33030011
```

The command prints no row. Both engines run with `--stats_text` and the corpus's own frequency flags: `--check_frequency --freq_pop AF --freq_freq 0.01 --freq_gt_lt gt --freq_filter exclude`.

### 1.3 What each engine writes

The `[General statistics]` section of each report. Ensembl VEP 116.3:

```
[General statistics]
Lines of input read	2
Variants processed	2
Variants filtered out	0
Novel / existing variants	0 (0.0) / 2 (100.0)
Overlapped genes	3
Overlapped transcripts	20
Overlapped regulatory features	-
```

vep-rs:

```
[General statistics]
Lines of input read	2
Variants processed	2
Variants filtered out	0
Novel / existing variants	1 (50.0) / 1 (50.0)
Overlapped genes	3
Overlapped transcripts	20
Overlapped regulatory features	-
```

Both outputs print `rs1254937786` in the `Existing_variation` column of the first record's rows and `-` in the second's, and the two outputs are equal row for row.

### 1.4 Why Ensembl VEP writes that

Under `--check_frequency` the known-variant source runs `frequency_check_buffer` over the buffer of input variants once the lookup is done, and that routine reads the frequency data of every variant. The read walks the variant's list of known variants (line 195 below). To walk a list that is not there, Perl first creates it, empty. The loop body never runs and no frequency is stored, but the variant leaves the routine holding an empty list.

The statistics logger counts a variant as existing when its list of known variants is a true value (line 244 below), and in Perl a reference to a list is true whatever the list holds. Every variant the routine keeps therefore counts, and the row's second figure equals `Variants processed`.

```
foreach my $ex(@{$vf->{existing}}) {            # BaseCacheVariation.pm:195
$stats->{existing}++ if $vf->{existing};        # Stats.pm:244
```

Without the flag the routine does not run. Ensembl VEP reads the cache through the `Bio::DB::HTS::Tabix` module, and on that route a variant with no match holds no list, so the logger passes over it.

The source keeps a second route, through the `tabix` program, for a run in which the module does not load. That route hands every variant with no match an empty list, with the flag or without it (`VariationTabix.pm`, line 196). Release 116.3 loads the module at start-up and stops without it (`Bio/EnsEMBL/IO/TabixParser.pm`, line 41), so no run takes the second route.

The writers are unaffected. They test the list's length before they read it, so the `Existing_variation` column prints `-` for the second record.

Source (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/Config.pm`, lines 441-453: `--filter_common` sets `check_frequency`, which sets `check_existing`.
- `Bio/EnsEMBL/VEP/AnnotationSource/Cache/VariationTabix.pm`, lines 216-244 and 127, and `Bio/EnsEMBL/VEP/AnnotationType/Variation.pm`, line 96: the lookup through the module gives a variant a list only when a cache row matches it, and the frequency check runs after the lookup, under the flag alone.
- `Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm`, lines 316-317 and 195: every variant's frequency data is read, by a `foreach` over its list of known variants.
- `Bio/EnsEMBL/VEP/Stats.pm`, lines 244 and 726-735: the count, and the row written from it.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1046: the writer's test of the list's length.

### 1.5 Why vep-rs's row is the correct one

The row's two figures split the processed variants by whether a known variant matches them, and the cache answers that for each record without either engine. It holds a matching variant at the first record's position and no variant within a base of the second's. One of the two is existing and one is novel: `1 (50.0) / 1 (50.0)`.

Ensembl VEP prints that row for the same records under `--check_existing` alone. Under the frequency check its own `Existing_variation` column names a known variant for one record of the two, while its report counts two. The filter dropped nothing here, so both runs processed the same variants.

### 1.6 Try it yourself

Run a VCF of the two records through both engines, writing the text report. vep-rs's command is below, `<corpus>` being the corpus's directory; Ensembl VEP takes the same flags with the cache options of [Reproducing these examples](reproducing.md):

```
vep -i two_records.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --stats_text \
  --check_frequency --freq_pop AF --freq_freq 0.01 --freq_gt_lt gt \
  --freq_filter exclude
```

The `Novel / existing variants` row of each run:

| Run | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `--check_existing`, both records | `1 (50.0) / 1 (50.0)` | the same |
| `--check_frequency`, both records | `0 (0.0) / 2 (100.0)` | `1 (50.0) / 1 (50.0)` |
| `--check_frequency`, the second record alone | `0 (0.0) / 1 (100.0)` | `-` |
| `--check_existing`, the corpus's 15 records | `1 (6.7) / 14 (93.3)` | the same |
| `--check_frequency`, the corpus's 15 records | `0 (0.0) / 11 (100.0)` | `1 (9.1) / 10 (90.9)` |

The corpus run under `--check_frequency` drops four records and processes eleven. Ten of the eleven carry an identifier in the `Existing_variation` column of both engines' outputs.

Under `--freq_filter include` the two rows agree. A variant with no known variant has no frequency and is dropped, so every variant such a run writes is existing: at a threshold of 0.0001 the corpus run prints `0 (0.0) / 12 (100.0)` in both engines.

### 1.7 How the comparator treats this class

No committed run of the test corpora keeps the statistics report of a run under `--check_frequency`. The five that run with the flag, in the `GRCh38-freqfilter` and `GRCh38-freqfilter-include` corpora, compare the annotation outputs alone, and `crates/vep-cli/tests/stats_existing_variants.rs` holds vep-rs's rows of this page.

The repository's comparison of two statistics reports, `tools/vep-stats-compare`, holds the row to this class. The rule applies where Ensembl VEP's report names `--check_frequency` or `--filter_common` in its `Command line options` row and its `Novel / existing variants` row counts every one of its `Variants processed` as existing. vep-rs's row then agrees when it is `-`, or a split whose two counts sum to that total, each with its percentage of the total at one decimal. A split counts at least one existing variant, since a run that counts none writes `-`. Any other value is a difference, every other row and section is compared as written, and the comparison's report lists the row under `documented` with the class's name.
