# Intended divergences from Ensembl VEP

A sample of the index in its regrouped shape: one table per part of a run, every entry named by an identifier in its first cell, the defects folded into the tables as entries whose page is a defect page.

## How the entries are ordered {#differences-outside-the-classes}

The entries are grouped by the part of a run they touch, one section each below.

## SNP and indel calling

| Entry | In one sentence | Tier and count | Verified against | Page |
| --- | --- | --- | --- | --- |
| [`start_cooccurrence_swap`](divergences/defect-b.md) | `start_lost` and `start_retained_variant` are both written for one allele | 1 · no flag needed | 116.2 | [Defect B](divergences/defect-b.md), classes 3 and 4 |
| [`splice_family_swap`](divergences/defect-a.md) | A later differing region erases the splice-region term | 1 · no flag needed | 116.2 | [Defect A](divergences/defect-a.md), classes 1 and 2 |

## Structural-variant calling

| Entry | In one sentence | Tier and count | Verified against | Page |
| --- | --- | --- | --- | --- |
| [`reference_skipped_record`](divergences/defect-c.md) | A record marked skipped is dropped before annotation | 1 · no flag needed | 116.2 | [Defect C](divergences/defect-c.md), class 5 |

## Plugins

### SpliceAI

| Entry | In one sentence | Tier and count | Verified against | Page |
| --- | --- | --- | --- | --- |
| [`spliceai_reverse_strand_indel`](plugins/spliceai.md#reverse-strand-indels) | Ensembl VEP 116.2's SpliceAI plugin complements a reverse-strand indel's alleles twice; vep-rs complements them once | 3 · 9 flag sets | 116.2 | [Reverse-strand indels](plugins/spliceai.md#reverse-strand-indels) |

## The classes

| # | Class | Defect |
| --- | --- | --- |
| 3 | `start_lost` and `start_retained_variant` emitted together | B |
