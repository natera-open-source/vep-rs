# The order of NearestExonJB's exon items

**The defect in one sentence.** Ensembl VEP's NearestExonJB plugin files a row's exons in a hash keyed by exon and builds the value by walking the hash's keys, so two or more items come in Perl's hash order, which changes from one run to the next. vep-rs lists the items in transcript order, so every run of an input writes the same value.

**Why this is a defect.** The record and the transcript fix the items: each names an exon, its distance from the variant, the exon end that distance is measured to, and the exon's length. Their order in Ensembl VEP's output comes from neither. The same input, flags and cache list one exon first on one run and the other first on another, so two runs of one command differ byte for byte. Within a single run, two transcripts that share both exons can list them in opposite orders.

**What vep-rs does instead.** vep-rs lists the items in transcript order, the exon the transcript reads first standing first, so every run of an input writes the same value. The items themselves are the ones Ensembl VEP writes. The plugin's parameters and its field are described under [NearestExonJB](../plugins/transcript-model.md#nearestexonjb).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `NearestExonJB.pm`, lines 219-240: each exon of the row is filed in the hash `%dists` under its stable id, with the shorter of its two distances to the variant, the exon end that distance is measured to, and its length.
- Lines 267-271: `foreach my $exon (keys %dists)` pushes one item for every exon at the smallest distance, and line 274 joins the items with `,`.
- Lines 245-263: the parameter `intronic=1` gives a variant inside an intron, which overlaps no exon, the two exons beside that intron (lines 190-217). A second walk over `keys %dists` pushes an item for both, whatever their distances.

**What correcting it would change.** Walking the exons in the transcript's order, or sorting the hash's keys, would write the items in one order on every run.

**The class it produces.** `nearest_exon_jb_item_order`: a `NearestExonJB` value (JSON `nearestexonjb`) of two or more items. Ensembl VEP's items are vep-rs's, in an order Perl's hash order decided for that run. A value with an item changed, missing or added is not the class and stays a difference.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md). The trigger is a value of two or more items, which the transcript decides and the record alone does not, so the rate is counted from the output rows of `--plugin NearestExonJB` over the cuts. Two of the cuts' 15,646,481 records have such a row, 0.13 per million, on nine rows. Both are the record of section 1, ClinVar variation 2444492, in the ClinVar cuts of the two assemblies, with seven of the nine rows on GRCh38 and two on GRCh37; no record of the three population cuts has one.

Without `intronic=1` a value has two items only where two overlapped exons lie equally near. With `intronic=1`, which the count above does not use, every row of an intronic variant lists the two exons beside its intron. Two committed flag sets run the plugin, both in the `GRCh38-plugin-computed` corpus, one of them with `intronic=1`. That corpus is the one whose committed reference outputs hold such a value: 7 of 961 values under the default run and 362 of 1,060 under `intronic=1`.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion whose ends lie equally near two exons

**In short.** ClinVar variation 2444492 deletes the 1,783 bases `21:46398011-46399793`. On `ENST00000703224` the deletion starts three bases before a 117-base exon and ends three bases before the last base of a 207-base exon, so both exons lie at the smallest distance, 3, and the value has two items. Ensembl VEP 116.2 writes the 207-base exon first when Perl's hash order is fixed, and either order when it is not. vep-rs writes the 117-base exon first, the transcript's order, on every run.

### 1.1 What a reader needs to know first

NearestExonJB reports the exon boundary nearest a variant that overlaps exons of a transcript. For each overlapped exon it measures the variant's start to the exon's start and the variant's end to the exon's end, and keeps the shorter of the two. An item reads `ExonID|distance|start/end|length`, where `start` or `end` names the exon end the distance is measured to. When several exons lie at the smallest distance, the value carries one item for each, joined by `,`.

A Perl hash holds its keys in no order of its own. `keys` returns them in an order set by the hash function and its seed, and Perl draws a new seed for each process unless the environment variable `PERL_HASH_SEED` fixes one. Perl also perturbs the order hash by hash as keys are stored, unless `PERL_PERTURB_KEYS=0` switches that off. So two hashes that hold the same keys can list them differently inside one run, and the order is fixed only with both variables set.

### 1.2 The record

The record is ClinVar variation 2444492 on GRCh38, one of the ClinVar GRCh38 dataset's chromosome 21 records. It is also the 34th record of the `GRCh38-plugin-computed` golden corpus (`tests/golden/116/GRCh38-plugin-computed/variants.vcf`). Its reference allele is 1,784 bases long and its alternate allele is the first of them, so the line is not printed here. The deleted bases are `21:46398011-46399793`.

`ENST00000703224` lies on the forward strand, and three of its exons overlap the deletion:

| Exon | Span | Length | Shorter distance |
| --- | --- | --- | --- |
| `ENSE00003618058` | 46398014-46398130 | 117 | 3, start to start |
| `ENSE00003541937` | 46398235-46398255 | 21 | 224, start to start |
| `ENSE00003551868` | 46399590-46399796 | 207 | 3, end to end |

Both engines run `--plugin NearestExonJB` in the default format.

### 1.3 What each engine writes

The row on `ENST00000703224`. Ensembl VEP 116.2, run with `PERL_HASH_SEED=0` and `PERL_PERTURB_KEYS=0`:

```
2444492	21:46398011-46399793	-	ENSG00000160299	ENST00000703224	Transcript	splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,3_prime_UTR_variant,intron_variant,NMD_transcript_variant	?-5034	-	-	-	-	-	IMPACT=HIGH;STRAND=1;NearestExonJB=ENSE00003551868|3|end|207,ENSE00003618058|3|start|117
```

vep-rs:

```
2444492	21:46398011-46399793	-	ENSG00000160299	ENST00000703224	Transcript	splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,3_prime_UTR_variant,intron_variant,NMD_transcript_variant	?-5034	-	-	-	-	-	IMPACT=HIGH;STRAND=1;NearestExonJB=ENSE00003618058|3|start|117,ENSE00003551868|3|end|207
```

The rows agree in every other column. With the two environment variables unset, ten runs of Ensembl VEP wrote the 207-base exon first on seven and the 117-base exon first on three; vep-rs wrote the 117-base exon first every time.

The record's rows on its six other transcripts carry two items as well, and each of them changed its order between runs. `ENST00000995542` shares both exons with `ENST00000703224`, and in five of the ten runs the two rows listed them in opposite orders.

### 1.4 Why Ensembl VEP writes that

The plugin files each overlapped exon in a hash under the exon's stable id, with the shorter of its two distances. It then walks `keys` of that hash and pushes an item for every exon whose distance is the smallest. The 117-base and the 207-base exon both lie at distance 3, so both are pushed, in the order `keys` returns them; the 21-base exon, at 224, is not.

That order is the hash's own. It depends on where the seeded hash function places each key and on the perturbation Perl applies as keys are stored. The seed is drawn per process and the perturbation per hash, so the order moves from run to run and, inside a run, from one hash to another. Nothing sorts the keys, and the order in which the plugin visited the exons is not kept.

Source (`VEP_plugins` at release/116):

- `NearestExonJB.pm`, lines 219-240: the store into `%dists` under the exon's stable id.
- Lines 267-271: the walk over `keys %dists` that builds the value.

### 1.5 Why vep-rs's row is the correct one

The record's coordinates fix both items. The deletion starts at 46398011, three bases before the 117-base exon's first base at 46398014. It ends at 46399793, three bases before the 207-base exon's last base at 46399796. Both engines print those two items and no other.

Only the order differs, and Ensembl VEP's own runs of this one input disagree about it, so neither of its orders can be read from the record or the transcript. vep-rs's order is the transcript's: on the forward strand the 117-base exon precedes the 207-base one. It is the same on every run, so two runs of one input can be compared byte for byte.

### 1.6 Try it yourself

The command as vep-rs takes it:

```
vep -i <corpus>/variants.vcf -o out.txt --offline \
  --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --plugin NearestExonJB
```

`<corpus>` is `tests/golden/116/GRCh38-plugin-computed`. vep-rs writes `NearestExonJB=ENSE00003618058|3|start|117,ENSE00003551868|3|end|207` on the row of `ENST00000703224` at `21:46398011-46399793`. Ensembl VEP 116.2 takes the same plugin flag with the cache options of [Reproducing these examples](reproducing.md). With `PERL_HASH_SEED=0` and `PERL_PERTURB_KEYS=0` in its environment it writes the two items in the other order; without them, repeated runs write both orders.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies `nearest_exon_jb_item_order` when a comparison case, one run of both engines over one input, names it. A `NearestExonJB` whose two values carry the same items, each as often, in different orders is set aside and counted once per row under the class, whatever the format. A value with an item changed, missing or added, a value one engine alone writes, and every other field stay differences.

A case may name the comparator's built-in class `order_only` instead. That class reads `NearestExonJB`, like the `DOMAINS` and `FREQS` columns, as a multiset of its items and counts the same rows under its own name; the golden harness compares the `GRCh38-plugin-computed` corpus under it. Where a case names both, `order_only` takes the row.
