# Abbreviated fields, `is_child` and the regex dialect of `--transcript_filter`

A difference in which transcripts a `--transcript_filter` keeps when its field is an abbreviation, and in which filters run at all; nothing is set aside, since no measured dataset runs the flag.

## An abbreviated field on the first transcript

Ensembl VEP looks a filter's field up on each transcript in turn (`FilterSet.pm` lines 533-613, `get_input`). It tries the exact key, then a synonym recorded earlier in the run, then the field in lower case, in upper case and with a leading `_` (lines 552-565). Last it tries a unique case-insensitive prefix of the transcript's keys (lines 567-572). The exact-key branch reads the value, and the lower-case, upper-case and leading-`_` branches read it as they record the synonym.

The prefix branch records the synonym and reads nothing, so the first transcript a run evaluates against an abbreviated field has no value: a unit with a value fails it (lines 481-483) and `nex` passes it. Every later transcript reads the value through the recorded synonym (lines 546-548).

vep-rs resolves an abbreviation once, when the filter is parsed, over a fixed list of the cache's transcript keys, and reads every transcript's value, the first included. An exact, lower-case, upper-case or `_`-prefixed field name reads the same value on both engines from the first transcript on.

The record below is the `tests/golden/116/GRCh38-transcript-filter` corpus's `tf_17608184`, within reach of 41 transcripts; the corpus's 32 records are used for the whole-corpus counts.

```
21	17608184	tf_17608184	A	G	.	.	.
```

| Filter | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `stable match ENST` | 40 rows | 41 rows, ENST00000626994 `upstream_gene_variant` among them |
| `stable nex` | the ENST00000626994 row alone | one `intergenic_variant` row |
| `stable_id match ENST` | 41 rows | 41 rows |
| `en > 0`, the 32 records | 951 rows | 952 rows, `tf_7130369` on ENST00000623962 among them |

ENST00000626994 is the first transcript Ensembl VEP evaluates for `tf_17608184`, and ENST00000623962 (`downstream_gene_variant`) the first for `tf_7130369`. Under `stable_id match ENST` and `end > 0`, the exact keys, both engines write every row.

## An ambiguous prefix on a transcript without a symbol

Ensembl VEP's prefix search runs over the keys the transcript object holds. The cache dumper deletes `display_xref`, `external_db`, `external_display_name`, `external_name` and `external_status` only when they are defined (`Dumper/Core.pm` lines 155-172). A transcript without a gene symbol therefore keeps all five, undefined. On such a transcript the prefix `e` matches the four `external_` keys and `end` (line 568), and the run dies there with `Multiple data fields found matching query field "e"` and the five keys (lines 574-580). On a transcript with a symbol `e` names `end` alone.

vep-rs's key list leaves the five out, so `e` names `end` on every transcript. Over the corpus's 32 records, `e > 40000000` ends Ensembl VEP's run at ENST00000623962, the first transcript without a symbol, and writes 226 rows in vep-rs, the rows `end > 40000000` writes.

## The regex dialect

`match`, `matches`, `regex`, `re` and `nre` run their pattern through Perl's regex engine in Ensembl VEP (lines 714-715), case-insensitive and unanchored. vep-rs compiles it with the `regex` crate, which lacks back-references, look-around and the `\Z` and `\K` escapes. So `stable_id match ENST\Z` runs in Ensembl VEP, where no stable id ends in `ENST` and every record is intergenic. vep-rs refuses it before the input is opened, `ERROR: [filter] pattern "ENST\Z" is not supported: unrecognized escape sequence`, exit 1. The flag's row on [Filtering, picking and transcript sets](../cli/filtering-picking.md) and the `filter` component on [Compatibility with Ensembl VEP](../cli/compatibility.md#exit-status-and-errors) state the dialect.

## `is_child`

`is_child` tests the pattern as a regex first (line 726) and otherwise asks an ontology adaptor for the parent term's descendants (lines 731-747). The `vep` command builds its filter set without one (`AnnotationSource.pm` lines 87-90); only `filter_vep` sets an adaptor (`filter_vep` lines 127-131). Under `--transcript_filter` the operator therefore passes a transcript through the regex alone. It dies on the first transcript the regex does not match, `No ontology adaptor set` (line 736), exit 255, with the header written and no row.

vep-rs matches the regex first and otherwise walks its built-in Sequence Ontology table of consequence terms, so a transcript the regex does not match is kept when its value is a descendant of the parent term and dropped otherwise. `biotype is_child transcript_variant` ends Ensembl VEP's run at the first transcript; in vep-rs no biotype is a descendant of `transcript_variant`, so every transcript is dropped and each of the corpus's 35 alleles prints an `intergenic_variant` row.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--transcript_filter` with the filter quoted; vep-rs with the same flags over the same cache.
