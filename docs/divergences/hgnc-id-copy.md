# The `HGNC_ID` of a gene named outside HGNC

A difference in the text of the `HGNC_ID` field on a transcript whose gene takes its symbol from another source (EntrezGene, RFAM) and shares it with an HGNC gene elsewhere on the chromosome; nothing is set aside, since the comparators never see the field.

A gene named by HGNC carries its own HGNC id. A gene named by another source has none, and Ensembl VEP's `merge_features` (`AnnotationType/Transcript.pm`) lends it one: it records the id of every transcript carrying one under the transcript's gene symbol (line 290), then gives every transcript of that symbol the recorded id (line 305).

The routine runs twice. The cache dumper runs it over each 1 Mb region's transcripts as it writes the region's file (`Pipeline/DumpVEP/Dumper/Core.pm` line 109), so a transcript stored in two region files, one spanning the boundary, can carry an id in one file and none in the other. The annotation run then runs it again over the features fetched for every buffer (`AnnotationSource.pm` lines 139-143), so a same-named HGNC gene anywhere in the regions a buffer's records brought in lends its id too.

A transcript's `HGNC_ID` in Ensembl VEP therefore depends on which region file the record's position read it from, and on the other records of the buffer. The release-116 GRCh38 cache carries 134 transcripts, of four gene pairs, whose two region-file copies differ, and four far-apart pairs of one symbol that reach the buffer-time copy.

vep-rs's cache builder lends the id once, as the cache is built, from the transcripts overlapping the region the transcript starts in, and stores each transcript in that region's shard alone. Nothing runs at annotation time, so a transcript prints one `HGNC_ID` whatever the input. The two engines agree on every HGNC-named gene, and on a gene named elsewhere whose same-named HGNC gene lies in every region the transcript spans.

The first example is NMRAL2P on chromosome 3. ENST00000824369, a transcript of the EntrezGene-named gene ENSG00000290862, spans 185959973 to 186012561, across the boundary at 186000000. The HGNC-named NMRAL2P gene ENSG00000171658 lies at 185960039 to 185979837, in the lower region alone. Under `--symbol`, for one record at a time:

| Record | Ensembl VEP | vep-rs, native cache |
| --- | --- | --- |
| `3 186010000 A C` | `SYMBOL_SOURCE=EntrezGene`, no `HGNC_ID` | `HGNC_ID=HGNC:52332` |
| `3 185970000 A C` | `HGNC_ID=HGNC:52332` | `HGNC_ID=HGNC:52332` |

The second is SNORD81 on chromosome 1. ENST00000363064, the RFAM-named gene ENSG00000199934 at 85592280, shares its symbol with the HGNC-named ENSG00000310740 at 173864146, 88 Mb away, which carries `HGNC:10101`:

| Input | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `1 85592300 A C` alone | no `HGNC_ID` | no `HGNC_ID` |
| the same record with `1 173864200 A C` in one buffer | `HGNC_ID=HGNC:10101` | no `HGNC_ID` |

The Ensembl VEP columns are the rows of a run over the release-116 GRCh38 cache, and the transcript counts above come from a scan of that cache's region files; the vep-rs columns follow from the builder's rule above.

A transcript's identifiers are a fact of the transcript. A value that changes with the record's position or the file's other records describes the run rather than the gene, so the buffer-time copy is a defect of the lookup and not a reading vep-rs reproduces.

A cache converted from Ensembl VEP's region files holds both copies of a boundary transcript, one per shard. The loader reads a chromosome's shards in order of their start position and keeps the first copy of each transcript, so a converted cache prints the lowest region file's copy, whatever the input, with the `HGNC_ID` the dumper lent over that region's transcripts. The dumper writes a transcript into every region file its gene overlaps, so that file is the region the transcript starts in unless its gene begins in an earlier one. In the first case the dumper and the builder lend over the same region; in the second the earlier region's lending decides the column.

To reproduce in Ensembl VEP 116.2: the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--symbol`, each input above as its own file; the second example's two records go in one file under the default `--buffer_size`. In vep-rs: the same flags over a native cache built from the release-116 core dumps, or over a converted cache for the second example.
