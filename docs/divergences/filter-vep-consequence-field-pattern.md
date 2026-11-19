# The consequence entries `filter_vep` reads from a VCF line

**The defect in one sentence.** Ensembl VEP's `filter_vep` finds the consequence entries of a VCF line with one pattern over the whole line: the field's name, `=`, then one or more characters up to the next `;`, blank or line end. So an empty entry followed by anything, another INFO item or the FORMAT column of a line with samples, takes that text as its one entry, and a key that merely ends in the name, `OLD_CSQ`, is read as the field too. vep-rs's `vep filter` reads an empty entry as no entry and takes entries from the named key alone.

**Why this is a defect.** The tool's own INFO loop, just above the pattern, splits the column on `;` and `=` and keys `OLD_CSQ` apart from `CSQ`, and its help names "the CSQ INFO key". The same pattern reads an empty entry at the end of a line as no entry at all. So what an empty entry means depends on where it sits, and another key's value is evaluated as consequence entries. Under `--only_matched` the rewrite built on the same pattern deletes the swallowed item or column, or empties the wrong key.

**What vep-rs does instead.** `vep filter` takes the entries from every INFO item whose key is exactly the field's name and whose value holds an entry, splits each value on `,`, and under `--only_matched` rewrites each such item to its own passing entries. A line with no such item is evaluated on its columns and INFO keys and counted among the lines without the field, wherever an empty entry sits.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `filter_vep`, line 261: `while($line =~ m/($vcf_info_field)\=(.+?)(\;|$|\s)/g)`, unanchored, the value one or more characters. Lines 262-265 split each match on `,` and parse every chunk gathered so far again, so the entries before a second match are evaluated twice.
- `filter_vep`, lines 267-271: a line the pattern matches nowhere is evaluated on its own hash and counted as lacking the field, which is how an empty entry at the end of a line reads.
- `filter_vep`, line 296: under `--only_matched`, `s/$vcf_info_field\=(.+?)(\;|\s|$)/.../` rewrites the first match in the line.
- `filter_vep`, lines 252-258: the INFO loop that keys every item but the literal `CSQ` by the text before its first `=`.

**What correcting it would change.** An empty entry would read as no entry wherever it sits, and only the named key's value would be filtered or rewritten, as `vep filter` reads them.

**The class it produces.** `filter_vep_consequence_field_pattern`: a line `filter_vep` passes or rewrites through an empty entry's neighbour or through another key ending in the field's name. The class is documented only: no committed tool run carries either shape, and the tool runs are compared byte for byte (tier 3 of [Intended divergences from Ensembl VEP](../intended-divergences.md)).

## 1. An empty entry before another INFO item

**In short.** On `21 100 r1 A G . . CSQ=;DP=4`, `filter_vep` reads one entry, `;DP=4`, so its `Allele` field is defined and `--filter "Allele exists"` passes the line. `vep filter` reads no entry and the filter fails it. The sibling `21 200 r2 C T . . DP=4;CSQ=` reads as a line without the field on both engines.

Under `--soft_filter --only_matched --filter "IMPACT is HIGH"` nothing on the first line passes. `filter_vep` has one entry to remove, the swallowed text, and writes `CSQ=`, losing `DP=4`. `vep filter` has no entry to remove and leaves `CSQ=;DP=4` as it is.

```
engine       record  FILTER           INFO written
Ensembl VEP  r1      filter_vep_fail  CSQ=
vep-rs       r1      filter_vep_fail  CSQ=;DP=4
```

The end-of-run warning counts the lines without the field: one in `filter_vep` (r2), two in `vep filter` (r1 and r2).

On a line with sample columns an empty entry at the end of INFO is followed by a tab and the FORMAT column, which the pattern takes as the entry's text. On `21 100 r1 A G . . DP=4;CSQ= GT 0/1`, `--filter "Allele exists"` passes the line in `filter_vep` and fails it in `vep filter`. Under `--soft_filter --only_matched --filter "IMPACT is HIGH"`, `filter_vep` writes the line without its `GT` column, `DP=4;CSQ=` then `0/1`; `vep filter` writes its columns as they are.

## 2. A key ending in the field's name

**In short.** The record is `21 100 r1 A G . . OLD_CSQ=G|stop_gained|HIGH|S1|;CSQ=G|intron_variant|MODIFIER|S1|`, under a header whose `CSQ` format is `Allele|Consequence|IMPACT|SYMBOL|DISTANCE`. `--filter "IMPACT is HIGH"` passes the line in `filter_vep`, which reads `OLD_CSQ`'s value as an entry, and fails it in `vep filter`, whose one entry has `IMPACT` `MODIFIER`. `--filter "OLD_CSQ match stop_gained"` passes the line on both: both key `OLD_CSQ` as an INFO field of its own.

Under `--only_matched --filter "IMPACT is MODIFIER"` the one real entry passes, so `vep filter` writes the line unchanged. `filter_vep` evaluates three hashes over two chunks; the one that passes is the third, whose raw chunk is undefined, so the field is rewritten to nothing. Its rewrite lands on the first `CSQ=` in the line, inside `OLD_CSQ`:

```
engine       INFO written
Ensembl VEP  OLD_CSQ=;CSQ=G|intron_variant|MODIFIER|S1|
vep-rs       OLD_CSQ=G|stop_gained|HIGH|S1|;CSQ=G|intron_variant|MODIFIER|S1|
```

With the keys in the other order, `CSQ=G|intron_variant|MODIFIER|S1|;OLD_CSQ=G|stop_gained|HIGH|S1|`, `--only_matched --filter "IMPACT is HIGH"` makes `filter_vep` pass the line on `OLD_CSQ`'s entry and write the real field emptied, `CSQ=;OLD_CSQ=G|stop_gained|HIGH|S1|`. `vep filter` writes no line. A key the name is a prefix of, `CSQX`, is read as the field by neither engine, since `=` must follow the name.

To reproduce: the `ensemblorg/ensembl-vep:release_116.2` image's `filter_vep` and `vep filter` on the records above, under a header declaring `DP`, `OLD_CSQ` and the `CSQ` format quoted, with the flags quoted; neither tool needs a cache.
