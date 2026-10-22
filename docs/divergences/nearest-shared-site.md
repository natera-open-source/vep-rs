# The transcripts of a shared start site under `--nearest`

A difference in the text of the `NEAREST` field wherever several protein-coding transcripts start at the winning site; nothing is set aside, since no measured dataset runs the flag.

Under `--nearest`, Ensembl VEP keeps every protein-coding transcript's start site in an interval tree, one point per transcript. The first run on a cache directory writes the table (`AnnotationSource/Cache/Transcript.pm` lines 275-321). Each line is a transcript's start on the forward strand or its end on the reverse (`AnnotationType/Transcript.pm` lines 747-751), inserted as a point (`TranscriptTree.pm` line 165). For each record the tree is asked for one object below the record's end and one above its start (`TranscriptTree.pm` lines 233-236). The nearer of the two is kept, both on a tie (lines 242-253), and the value each object carries is printed once (`AnnotationType/Transcript.pm` lines 707-711).

A site several transcripts share is several points at one coordinate, and each search returns one of them. Nothing in Ensembl VEP's source chooses which. The object depends on the order of the table file's lines, which `readdir` gives (`Cache/Transcript.pm` line 297), and on the direction of the search that found the site. So one record can print one transcript of a site and the next record another, and a record both searches land on can print two of them.

vep-rs lists every transcript starting at a winning site, in stable-id order. Under `gene` and `symbol` each value is printed once, as `get_nearest` prints it, so the two engines agree wherever the transcripts of the site belong to one gene. Where transcripts of several genes share a site, Ensembl VEP prints the gene of the one object and vep-rs each gene once.

The records below sit around the start site of GRIK1 at 29940022, which seven protein-coding transcripts share (the gene is on the reverse strand, so the site is each transcript's end):

```
21	30000000	intergenic_del25	GCCCACTTTTACTTACTCACGTCTTT	G	.	.	.
21	29940020	two_below_site	T	A	.	.	.
```

Under `--nearest transcript` the two engines write, on every row of each record:

| Record | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `intergenic_del25` | `ENST00000399914` | `ENST00000389124,ENST00000389125,ENST00000399907,ENST00000399909,ENST00000399913,ENST00000399914,ENST00001131791` |
| `two_below_site` | `ENST00000389125` | the same seven |

The deletion lies above the site, which the search below its end returns; the substitution lies two bases below it, which the search above its start returns. Each search returns a different transcript of the site. The same records on the same cache, with a table file whose lines fall in another order, print `ENST00001131791` for the deletion. At the start site of POTED at 13609777, which two transcripts share, a record on the site prints `ENST00000620442` and a record two bases below it `ENST00000299443`; vep-rs prints `ENST00000299443,ENST00000620442` for both.

Under `--nearest gene` and `--nearest symbol` every row of these records agrees, `ENSG00000171189` and `GRIK1`. Over the 47 records of `tests/golden/116/GRCh38-identifiers` under `--nearest transcript`, 326 of 808 rows differ, each one id against a list containing it; under `gene` and `symbol` every row agrees. Chromosome 21 of the release-116 GRCh38 cache holds one site shared by three genes, 34073509, where Ensembl VEP prints `ENSG00000243927` and `MRPS6` for a record on the site.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--nearest transcript`, the cache directory writable for the table its first run writes; vep-rs with the same flag over the same cache.
