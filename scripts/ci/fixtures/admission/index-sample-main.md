# Intended divergences from Ensembl VEP

A sample of the index in its two-table shape: defects named by letter with a page each, differences named by their first-cell text with a page or a section of another page. The tables that follow the two entry tables carry plain first cells, which name no entry.

## The defects

| Defect | In one sentence | Classes | Page |
| --- | --- | --- | --- |
| [A](divergences/defect-a.md) | A later differing region erases the splice-region term | 1, 2 | [Defect A](divergences/defect-a.md) |
| [B](divergences/defect-b.md) | `start_lost` and `start_retained_variant` are both written for one allele | 3, 4 | [Defect B](divergences/defect-b.md) |
| [none](divergences/cnv-tr.md) | A representation difference, not a defect: `<CNV:TR>` | 12 | [`<CNV:TR>`](divergences/cnv-tr.md) |
| [M](divergences/defect-m.md) | `*` as the only ALT is annotated from `SVTYPE` | 20 | [Defect M](divergences/defect-m.md) |

## Differences outside the classes

| Difference | In one sentence | Set aside | Page |
| --- | --- | --- | --- |
| [JSON escaping](divergences/json-escaping.md) | Ensembl VEP's JSON percent-encodes the `=` of a synonymous HGVSp notation; vep-rs's JSON is unescaped | nothing: the comparators never see the field | [Unescaped JSON under `--output_format json`](divergences/json-escaping.md) |
| [A bare `M`](divergences/bare-m-chromosome.md) | Under an Ensembl cache's synonyms a record on `M` overlaps nothing in Ensembl VEP; vep-rs annotates it against the mitochondrial transcripts | nothing | [A bare `M` under a cache whose synonyms name `chrM`](divergences/bare-m-chromosome.md) |
| [A chromosome the cache lacks](cache-setup.md#chromosome-names) | Ensembl VEP admits a record whose name is a synonym of a cache chromosome and writes `intergenic_variant`; vep-rs skips the record with the not-found warning | nothing: no measured dataset spells a chromosome that way | [Chromosome names](cache-setup.md#chromosome-names) |
| [SpliceAI on a reverse-strand indel](plugins/spliceai.md#reverse-strand-indels) | Ensembl VEP 116.2's SpliceAI plugin complements a reverse-strand indel's alleles twice; vep-rs complements them once | nothing: no measured run loads the plugin | [Reverse-strand indels](plugins/spliceai.md#reverse-strand-indels) |

## Where each defect is in Ensembl VEP's source

| Defect | What goes wrong in Ensembl VEP | Where | Classes |
| --- | --- | --- | --- |
| A | The splice-region verdict is stored by assignment | `BaseTranscriptVariationAllele.pm` line 215 | 1, 2 |
| B | `start_lost` and `start_retained_variant` are computed by two independent tests | `Utils/VariationEffect.pm` lines 851-963 | 3, 4 |

## The classes

| # | Class | Defect |
| --- | --- | --- |
| 1 | Splice-region term lost to an assignment order | A |
| 3 | `start_lost` and `start_retained_variant` emitted together | B |

The engines, inputs, caches and commands every walkthrough uses are on [Reproducing these examples](divergences/reproducing.md).
