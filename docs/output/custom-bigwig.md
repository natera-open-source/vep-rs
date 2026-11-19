# bigWig sources and conservation shortcuts

A bigWig file annotates each variant with the values of the intervals it matches, as a `--custom` source or through a conservation shortcut like `--phyloP`.

A `--custom` source in `format=bigwig` ([Custom annotation](../output-formats.md#custom-annotation)) is read through the bigWig's own index. This page states:

- how its records match and are named
- how the summary statistics and the span rule run
- what the three conservation shortcuts `--bigwig`, `--phyloP` and `--phastCons` expand to, with the two flags that name their files
- three behaviours of Ensembl VEP's reader that vep-rs does not reproduce

## bigWig sources

A bigWig source (`format=bigwig`) is read as Ensembl VEP reads it. A record is an interval of the file clipped to the variant's query window, with its 0-based start read as 1-based. The window runs from one base before the variant's start to one base past its end.

Under `type=exact` a record matches when its bounds equal the variant's. A one-base interval at a single-base variant so matches, and an interval reaching past the variant's end does not. Under `type=overlap`, `within` and `surrounding`, with `distance`, `overlap_cutoff`, `reciprocal` and `fields=PC`, the record is matched on the same bounds by the rules every format shares.

A record is named by its value, stringified as Ensembl VEP stringifies a number. Under `coords` it is named by the file's chromosome name with the record's 1-based start and end. Ensembl VEP 116.2 names a record that way where its span rule, below, lists it, but elsewhere with the chromosome name less its `chr` and the record's 0-based start (`File/BigWig.pm:420-427`); vep-rs prints the one form on both paths. The difference is documented on [The name of a bigWig record under `coords`](../divergences/bigwig-coords-record-name.md).

The file's chromosomes are reached by their names with a leading `chr` removed, through the run's chromosome synonyms and the `chr` rule. An input `21` so reads a file named `chr21`. `summary_stats` runs over the matched values: `min`, `max`, `mean`, `sum` and `count`, with a value of 0 a minimum.

A variant spanning more than 65,536 bases takes Ensembl VEP's span rule under `type=overlap` with `summary_stats` of `min` and `max` alone and no other matching option. The extrema are exact over the whole span. The record list holds the intervals of the span's first 32,768 bases, and `...` closes it where the span continues. Every variant is annotated once on that path; Ensembl VEP 116.2 lists each record of a shorter variant twice there. The difference is documented on [A bigWig record listed twice under `min` and `max`](../divergences/bigwig-min-max-record-listed-twice.md).

## Conservation shortcuts

`--bigwig <file>` stands for `--custom <file>,,bigwig,exact`. `--phyloP <N>` and `--phastCons <N>`, each repeatable, stand for one `exact` bigWig source each. Each is read from the directory `--ucsc_data_root` names, laid out as UCSC's goldenPath tree with `--ucsc_assembly` as the assembly:

```text
<root><assembly>/phyloP<N>way/<assembly>.phyloP<N>way.bw
<root><assembly>/phastCons<N>way/<assembly>.phastCons<N>way.bw
```

The root is joined to the assembly as given, so it is named with its trailing slash. `--phyloP` or `--phastCons` without `--ucsc_assembly` exits 1 with Ensembl VEP's text, `You must set --ucsc_assembly to use --phyloP`.

The default root is Ensembl VEP's remote `http://hgdownload.cse.ucsc.edu/goldenpath/`. vep-rs reads local files only, so a shortcut under the default root, or under any `--ucsc_data_root` holding a URL scheme, exits 1 with `ERROR: [custom] --ucsc_data_root must name a local directory laid out as <root>/<assembly>/phyloP<N>way/<assembly>.phyloP<N>way.bw; the same source is --custom file=<local>.bw,format=bigwig,type=exact`. A local mirror of the UCSC tree, or the bigWig itself through `--custom`, is the route.

The `--phyloP <N>` source's short name, its output key, is `phyloP<N>way`. Ensembl VEP's option set spells it `phlyoP<N>way` (`Config.pm:555`), a misspelling vep-rs does not reproduce; `--phastCons <N>` is `phastCons<N>way` in both. A run that needs the misspelt key names the source through `--custom ...,short_name=phlyoP<N>way,...`. The difference is documented on [The key of a `--phyloP` source](../divergences/phylop-short-name.md).

## A bigWig record listed twice after a neighbouring variant's capped walk

Three behaviours of Ensembl VEP 116.2's bigWig reader are not reproduced; this section's is a documented class. This section names the rows a comparison would set aside. No comparison tool codes the class, since no committed row falls under it, as [Intended divergences from Ensembl VEP](../intended-divergences.md) lists it. Its [How to read the counts](../intended-divergences.md#how-to-read-the-counts) defines the terms.

The class carries its tier under that page's rule:

- tier 1 is a shape the chromosome 21 population cuts carry with no flag
- tier 2 one that needs a flag or plugin on ordinary data
- tier 3 one ordinary data does not carry

The counts are over the cuts' 15,646,481 records, the 111 golden corpora and their 392 flag sets.

**The defect in one sentence.** For a bigWig source of `type=overlap` or `within` without `summary_stats`, a variant whose window holds more than `num_records` matching intervals stops its walk early and leaves the unread intervals queued in the reader, so the next variant is annotated against those leftovers before its own. A variant at the end of the stopped variant's span lists its record twice.

**Why this is a defect.** The reader's `seek` is documented as fetching the requested region's data and caching it, and the walk reads that cache as the region's records. The cache is a queue that `seek` appends to and only the walk consumes; nothing empties it, so after a walk that stopped at the cap the queue holds another variant's intervals ahead of the requested region's. A variant's record list should come from its own window alone, and the tabix readers of the VCF, BED and GFF formats do exactly that, opening a fresh iterator at every `seek`.

**What vep-rs does instead.** vep-rs reads each variant's own window, or a batch's merged windows cut back to each variant's, and keeps no reader state between variants; the cap lists the window's first `num_records` intervals and `...`.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `ensembl-io` at release/116):

- `Bio/EnsEMBL/IO/Parser/BigWig.pm`, lines 52-73, `seek`: the query's intervals are pushed onto `cache->{features}`, which `Bio/EnsEMBL/IO/BigFileParser.pm` creates once, line 57, and only ever shifts from, lines 350-358, `read_block`.
- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, line 337, `annotate_InputBuffer`: `last if !defined $stats && $record_count > $self->{num_records}` ends the walk with the queue unread. Lines 321-322 are the next variant's `seek` and `next`, and the walk's condition at lines 333-335, `get_start <= $vf_end + 1`, admits the leftovers.

**What correcting it would change.** Each variant would list the records of its own window once, as vep-rs does and as Ensembl VEP does when no walk reaches the cap.

**The class it produces.** `bigwig_overlap_capped_walk_carryover`: on a bigWig source of `type=overlap` or `within` without `summary_stats`, the `<short_name>` list of a variant that follows one whose list closed with `...`, where the reference repeats a value the candidate lists once. The committed `scores.bw` windows hold at most 4 intervals against the cap of 50, and the `--bigwig`, `--phyloP` and `--phastCons` shortcuts are `exact` sources, one record per variant, so it applies to no committed row or comparison case.

**How often the shape occurs.** Tier 2: 36,859 of the population cuts' records (2,356 per million) span more than 50 bases, which a per-base bigWig fills past the record cap; the following record inside the span shows it. 41 corpora carry the shape. No committed flag set names a bigWig source of `type=overlap` or `within` without `summary_stats`: the committed sources are `exact`, or carry `summary_stats`.

### A substitution at the last base of a 100-base deletion's span

**In short.** A bigWig `dense` holds 150 one-base intervals at `21:5028900-5029049` with the values `0.001` to `0.150`; the input is `21 5028900 N <DEL>` with `END=5029000`, then `21 5029000 A G` and `21 5029010 A G`, under `--custom file=dense.bw,short_name=dense,format=bigwig,type=overlap`. The deletion's walk stops after its 51st match, listing 50 values and `...` and leaving the rest of its window queued. The substitution's `seek` adds its own intervals behind them; its walk reads the leftovers, matches the queued interval at `5029000`, then matches its own copy of it. Ensembl VEP 116.2 prints `dense=0.101000003516674,0.101000003516674`, vep-rs `dense=0.101000003516674`; the deletion's row and the row at `5029010`, `dense=0.111000001430511`, agree on both engines.

## Where next

- [Custom annotation](../output-formats.md#custom-annotation): the `--custom` string, the matching rules every format shares and the keys a source adds.
- [The name of a bigWig record under `coords`](../divergences/bigwig-coords-record-name.md): the two forms Ensembl VEP 116.2 gives a record's name.
- [A bigWig record listed twice under `min` and `max`](../divergences/bigwig-min-max-record-listed-twice.md): the records a shorter variant lists twice under the span rule.
- [The key of a `--phyloP` source](../divergences/phylop-short-name.md): the misspelt key and the way to name it.
