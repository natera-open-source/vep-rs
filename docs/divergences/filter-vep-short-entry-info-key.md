# The field a short consequence entry lacks in `filter_vep`

**The defect in one sentence.** Ensembl VEP's `filter_vep` builds each consequence entry's hash over a copy of the VCF line's INFO keys and sets a header field undefined only when the hash lacks the key. So an entry shorter than the `Format` list reads, for each field past its last `|`, the INFO value of the same name: the input's `AF` is read as the consequence field `AF`. vep-rs's `vep filter` answers such a field from the entry alone, so it is undefined whatever INFO holds.

**Why this is a defect.** The step that completes an entry's fields is written to set its missing fields undefined, and does so for every field whose name is not also an INFO key. The reading is inconsistent on one input: an entry whose `AF` position is present and empty reads undefined, while an entry with no `AF` position reads the INFO value. And the value substituted is a quantity of another definition: the input's `AF` is the allele frequency its own `##INFO` line defines, where the consequence field `AF` is the 1000 Genomes frequency Ensembl VEP writes under `--af`. `AF > 0.1` passes an entry that carries no frequency annotation at all.

**What vep-rs does instead.** `vep filter` reads a consequence field the entry lacks as undefined, whatever the line's INFO column holds. An INFO key is a field of every entry on its line, read when the filter names a key outside the `Format` list; a name in the `Format` list is answered from the entry alone. So `AF exists` fails a short entry beside `AF=0.5`, and `AF nex` passes it.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `filter_vep`, lines 252-257: every INFO item but the literal `CSQ` is keyed into the line's hash by the text before its first `=`.
- `filter_vep`, line 389: `my %data = %{$main_data || {}};`, each entry's hash begins as a copy of the line's hash, INFO keys included.
- `filter_vep`, line 395: `$data{$headers->[$_]} = $split[$_] for 0..$#split;`, the entry's own fields overwrite the copy, as many as the entry has.
- `filter_vep`, lines 413-416: `$data{$_} = undef unless exists($data{$_});` under the comment `# set missing fields to undef`, so a field whose name is an INFO key is never set undefined.
- `filter_vep`, lines 364-366: the `Format` list is the last `##INFO=<ID=CSQ` header line's, so an earlier annotation's shorter entries are read under the later, longer list.
- `OutputFactory/VCF.pm`, lines 215-217 and 266-270: under `--keep_csq` the writer keeps the input's `##INFO=<ID=CSQ` line and appends its own, which is how one file comes to hold entries of two lengths.

**What correcting it would change.** A field past an entry's last `|` would read as undefined, as an empty field between two `|` does, and only a name outside the `Format` list would be answered from the INFO column, as `vep filter` reads them.

**The class it produces.** `filter_vep_short_entry_reads_info_key`: a line `filter_vep` passes or fails on an INFO value read as a consequence field the entry lacks. The class is documented only: no committed tool run carries the shape, and the tool runs are compared byte for byte (tier 3 of [Intended divergences from Ensembl VEP](../intended-divergences.md)).

## 1. A re-annotation under more flags beside an INFO key of a field's name

**In short.** The header declares `AF` and two `CSQ` lines, the first with `Format: Allele|Consequence|IMPACT` and the second with `Format: Allele|Consequence|IMPACT|AF`, as a `--keep_csq` run adding `--af` leaves them. Both engines read every entry under the second list. Only `s1`, a short entry beside an INFO `AF`, differs: `filter_vep` evaluates the input's `AF=0.5` as the consequence field `AF`, where `vep filter` reads that field as undefined. The records are:

```
21  100  s1  A  G  .  .  AF=0.5;CSQ=G|missense_variant|MODERATE
21  200  s2  C  T  .  .  AF=0.5;CSQ=T|missense_variant|MODERATE|0.01
21  300  s3  G  A  .  .  AF=0.5;CSQ=A|missense_variant|MODERATE|
21  400  s4  T  C  .  .  DP=9;CSQ=C|missense_variant|MODERATE
```

`s1` is a short entry beside an INFO `AF`; `s2` writes its own `AF`; `s3` writes an empty `AF`; `s4` is a short entry on a line with no INFO `AF`. The records each filter passes:

```
filter       Ensembl VEP   vep-rs
AF > 0.1     s1            none
AF exists    s1 s2         s2
AF nex       s3 s4         s1 s3 s4
AF < 0.1     s2            s2
AF is 0.5    s1            none
```

Only `s1` differs: `filter_vep` evaluates the input's `AF=0.5` as the consequence field `AF`, so `AF > 0.1`, `AF exists` and `AF is 0.5` pass it and `AF nex` fails it. On `s3`, whose `AF` position is present and empty, both engines read undefined, and so on `s4`, short on a line without the key; `s2` reads its own `0.01` on both, the INFO value ignored. Under `AF exists`, `--count` prints `2` in `filter_vep` and `1` in `vep filter`, and `--soft_filter` marks `s1` `filter_vep_pass` there and `filter_vep_fail` here.

To reproduce: the `ensemblorg/ensembl-vep:release_116.2` image's `filter_vep` and `vep filter` on the records above, under a header declaring `AF` and the two `CSQ` formats quoted, with `--filter` as quoted and `--no_missing_consequence_warning`; neither tool needs a cache.
