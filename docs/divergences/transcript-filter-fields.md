# Abbreviated fields, `is_child` and the regex dialect of `--transcript_filter`

A difference in which transcripts a `--transcript_filter` keeps when its field is an abbreviation, or a key a GTF or GFF3 transcript source holds differently from the cache, and in which filters run at all. Ensembl VEP's first transcript has no value under an abbreviated field, and its run dies at the prefix `e` on a transcript without a symbol or where an `is_child` regex does not match. vep-rs reads every transcript's value, reads `e` as `end`, answers `is_child` from its built-in table and alone refuses a pattern the `regex` crate lacks. Nothing is set aside, since no measured dataset runs the flag. The same filter grammar runs in `vep filter`, Ensembl VEP's `filter_vep`, which the closing section covers.

## An abbreviated field on the first transcript

Ensembl VEP looks a filter's field up on each transcript in turn (`FilterSet.pm` lines 533-613, `get_input`). It tries the exact key, then a synonym recorded earlier in the run, then the field in lower case, in upper case and with a leading `_` (lines 552-565). Last it tries a unique case-insensitive prefix of the transcript's keys (lines 567-572). The exact-key branch reads the value, and the lower-case, upper-case and leading-`_` branches read it as they record the synonym.

The prefix branch records the synonym and reads nothing, so the first transcript a run evaluates against an abbreviated field has no value: a unit with a value fails it (lines 481-483) and `nex` passes it. Every later transcript reads the value through the recorded synonym (lines 546-548).

vep-rs resolves an abbreviation once, when the filter is parsed, over a fixed list of the cache's transcript keys, and reads every transcript's value, the first included. An exact, lower-case, upper-case or `_`-prefixed field name reads the same value on both engines from the first transcript on.

The prefix branch reads the written name as a regular expression anchored at its start (`/^$field/i`, line 568), and vep-rs reads it the same way, so a name holding a metacharacter, `bio.ype` or `biotype$`, resolves to `biotype` on both engines.

vep-rs's fixed key list also holds `_source_cache`, `_gene_symbol_id` and `_bam_edit_status` on every cache. Ensembl VEP's transcript carries `_source_cache` on a merged cache (`Database/Transcript.pm` line 278) or from a GTF or GFF3 source (`BaseGXF.pm` line 635). It carries `_gene_symbol_id` on a RefSeq transcript whose gene has one EntrezGene entry (`Database/Transcript.pm` line 563). It carries `_bam_edit_status` on a model edited against the cache's alignment file (`AnnotationType/Transcript.pm` lines 601 and 621). A plain Ensembl cache's transcript holds none of the three.

So `_s` names `_swissprot` alone there on such a cache, and two keys here, where the run is refused as ambiguous. `_gene_symbol_` likewise names `_gene_symbol_source` alone on a transcript without the identifier, and two keys here. A filter naming any of the three keys in full reads the same value on both engines, undefined where the transcript lacks the key.

A `match` or `nre` unit with no operand (`biotype match`) matches nothing here, so `match` keeps no transcript and `nre` every one. Ensembl VEP 116.2 compiles the missing operand under an `uninitialized value` warning (lines 714-715) and keeps every transcript under `match` and none under `nre`.

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

## A GTF or GFF3 as the transcript source

Ensembl VEP hands its `--gtf` and `--gff` sources the filter set its cache takes (`AnnotationSourceAdaptor.pm` line 401). Each transcript the file builds is evaluated against it (`AnnotationType/Transcript.pm` line 116). vep-rs reads such a file through a cache built from it ([A GTF or GFF3 as the transcript model](../cache-from-gtf-or-gff3.md)), and `--transcript_filter` and `--exclude_predicted` run on that cache as on any other.

On a key both engines fill, the two keep the same transcripts. Over a set of ClinVar GRCh38 chromosome 21 records against the release 116 GTF, `biotype is protein_coding` writes the same rows on both engines, on the same transcripts. Over a subset of those records against the RefSeq GRCh38 GFF3, `--exclude_predicted` writes the same rows on both, every `XM_` and `XR_` model dropped.

The transcript Ensembl VEP builds from a GTF or GFF3 (`BaseGXF.pm`, `_create_transcript` and `_add_identifiers`) holds `stable_id`, `version`, `biotype`, `start`, `end`, `strand`, `dbID`, `slice`, `_trans_exon_array`, `_gene_stable_id`, `_gene_version`, `_gene_symbol` and `_source_cache`. It holds `is_canonical` when tagged and `attributes` when its line carries a support level or a tag. When coding it holds `_protein`, `translation`, `cdna_coding_start`, `cdna_coding_end` and `_variation_effect_feature_cache`. It holds no `source`, `_ccds`, `_gene_symbol_source`, `_gene_hgnc_id`, cross-reference or phenotype key, and no `_gene` or `_vep_lazy_loaded`.

The cache vep-rs builds from the file holds the cache's keys, so a filter on a key the two hold differently reads differently. The table shows the three a filter is likely to name, over that subset against the release 116 GTF. Unfiltered, both engines write the same rows there; a record whose every transcript is dropped writes one `intergenic_variant` row.

| Filter | Ensembl VEP `--gtf` | vep-rs, the cache built from the GTF |
| --- | --- | --- |
| `source is ensembl_havana` | one `intergenic_variant` row per record: the key is absent | the rows of the transcripts whose source column is `ensembl_havana`, and one `intergenic_variant` row for each record with none |
| `ccds exists` | one `intergenic_variant` row per record: the key is absent | every row: every transcript holds the cache's `-` placeholder |
| `_source_cache exists` | every row: the key holds the file's name | one `intergenic_variant` row per record: the key is undefined |
| `stable_id match ENST` | every row | every row |

`source` is the file's source column, which Ensembl VEP reads for the biotype of a GTF transcript without one but never stores. The `-` under `ccds` is what the Ensembl cache itself stores for a transcript without a CCDS, so `ccds exists` keeps every transcript of any cache on both engines; only Ensembl VEP's GTF and GFF sources lack the key. The keys no GTF or GFF3 carries (`_gene_symbol_source`, `_gene_symbol_id`, `_gene_hgnc_id`, `_refseq`, `_swissprot`, `_trembl`, `_uniparc`, `_uniprot_isoform`, `_gene_phenotype`, `_bam_edit_status`) are undefined on both engines: `exists` drops every transcript and `nex` keeps every one.

Five bookkeeping keys differ as well. `dbID` is a counter Ensembl VEP assigns to each transcript of the source and is absent from a cache the builder writes, so over the same records `dbID exists` writes every row there and one `intergenic_variant` row per record here. `_gene`, `_vep_lazy_loaded`, `attributes` and `_variation_effect_feature_cache` are defined on every vep-rs cache transcript. The transcript Ensembl VEP builds from the file never holds the first two, holds `attributes` only when its line carries a support level or a tag, and holds `_variation_effect_feature_cache` only when it is coding.

An abbreviated field resolves over that smaller key set in Ensembl VEP and over the cache's list in vep-rs ([the first section](#an-abbreviated-field-on-the-first-transcript)). `_s`, `_t` and `d` each name one key on the file's transcript (`_source_cache`, `_trans_exon_array`, `dbID`) and two in vep-rs's list, where the run is refused with `Multiple data fields found matching query field` before the input is opened. Over the same records `_s exists`, `_t exists` and `d exists` write every row in Ensembl VEP but the ones the first-transcript reads of the first section drop. `_u` names no key on the file's transcript and two in vep-rs's list (`_uniparc`, `_uniprot_isoform`): `_u exists` writes one `intergenic_variant` row per record there and is refused here.

## The regex dialect

`match`, `matches`, `regex`, `re` and `nre` run their pattern through Perl's regex engine in Ensembl VEP (lines 714-715), case-insensitive and unanchored. vep-rs reads the pattern as Perl reads it and compiles it with the `regex` crate: the case fold, the switch to Unicode rules, the brace quantifiers, the octal and control escapes and the flag groups follow Perl.

The crate lacks back-references, look-around, atomic, branch-reset and conditional groups, recursion, code blocks, backtracking control verbs, the `\b{wb}` bound types, the `(?l)` locale flag, the `\X` escape and a `\N{name}` character. vep-rs also refuses a literal run of more than 64 letters `s` under Unicode rules, whose regex grows with the square of the run, and groups nested more than 100 deep. It refuses a `\N{U+}` escape of several characters in a negated class, before a `-` or as the end of a range too, which Perl reads with a warning. A `-` after such an escape opens a range from the letter `N` of the escape in Perl, refused as inverted below `N`.

So `stable_id match ENST\X` runs in Ensembl VEP, where `\X` is one extended grapheme cluster and every stable id with a character after `ENST` passes. vep-rs refuses it before the input is opened, `ERROR: [filter] pattern "ENST\X" is not supported: unrecognized escape sequence`, exit 1. The flag's row on [Filtering, picking and transcript sets](../cli/filtering-picking.md) and the `filter` component on [Compatibility with Ensembl VEP](../cli/compatibility.md#exit-status-and-errors) state the dialect.

Four readings differ beyond the refusals. An alternation of literal words matches a value's `ß` in Perl through a word that ends inside the fold of `ß`: `^(?:ss|s)$` matches `sß` there, where vep-rs compares the whole fold, `sss`, with each word and matches neither. A `u`, `a` or `aa` charset flag rules its own group alone in both engines, but Perl's optimizer merges some adjacent literals across that boundary and folds the merged run whole. So `^(?u:s)s$` matches `ß` there, where vep-rs reads each letter under the charset of its own scope and matches `ss` and `sS` alone.

Perl also keeps a `\N{U+}` escape of several characters as a run of its own, so `^\N{U+73.73}s$` matches `ßS` there and not `Sß`, where vep-rs folds the whole run and matches both. A possessive quantifier, `s++`, runs as the greedy `s+`, so `s++s$` matches `sss` in vep-rs and never in Perl.

## `is_child`

`is_child` tests the pattern as a regex first (line 726) and otherwise asks an ontology adaptor for the parent term's descendants (lines 731-747). The `vep` command builds its filter set without one (`AnnotationSource.pm` lines 87-90); only `filter_vep` sets an adaptor (`filter_vep` lines 127-131). Under `--transcript_filter` the operator therefore passes a transcript through the regex alone. It dies on the first transcript the regex does not match, `No ontology adaptor set` (line 736), exit 255, with the header written and no row.

vep-rs matches the regex first and otherwise walks its built-in Sequence Ontology table of consequence terms, so a transcript the regex does not match is kept when its value is a descendant of the parent term and dropped otherwise. `biotype is_child transcript_variant` ends Ensembl VEP's run at the first transcript; in vep-rs no biotype is a descendant of `transcript_variant`, so every transcript is dropped and each of the corpus's 35 alleles prints an `intergenic_variant` row.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--transcript_filter` with the filter quoted; vep-rs with the same flags over the same cache.

## `vep filter`

`vep filter` runs the same grammar over the lines of an annotated output, as Ensembl VEP's `filter_vep` does ([its page](../cli/filter.md)). Five readings of `filter_vep` are corrected there rather than reproduced; everything else, the cleaning of `-` and empty values, the hash of a line without the consequence field, the `--start`, `--limit` and `--test` arithmetic and the `--soft_filter` column, is as `filter_vep` has it.

**An abbreviated field on the first line.** `filter_vep` resolves a field name against each line's own keys, through the synonym cache described above (`FilterSet.pm` lines 533-586). The prefix branch records the synonym and reads nothing (lines 570-571), so the first entry it evaluates under an abbreviated field has no value, in `filter_vep` as in `--transcript_filter`. `vep filter` reads the value on that first evaluation too.

The `tests/golden/116/GRCh38-filter-tool` corpus shows the difference. On its `annotated.tab.gz`, `Gen match ENSG` (the prefix of `Gene`) writes 160 rows in `vep filter` and 159 in `filter_vep`, which drops the first row. On its `annotated.vcf` under `--only_matched`, the same filter keeps every entry of the first line in `vep filter`, and `filter_vep` keeps all but the first. A field written as a key (`Gene match ENSG`) reads the same on both from the first line on.

**An INFO flag.** `filter_vep` splits each INFO item on `=` (`filter_vep` lines 253-257), so a flag key without `=`, `DB`, is a key with no value: undefined, which `DB` alone and `DB exists` fail and `DB nex` passes on the very line that carries it. `vep filter` reads a flag as defined and empty: `DB` passes the lines that carry the flag and `DB nex` the lines that do not.

The correction is the flag's alone. An empty value, `AF=`, is defined and empty on a line without a consequence entry on both engines, and a word without `=` in a tab line's `Extra` column is undefined on both, like `KEY=` there. A key with a value is undefined on a consequence entry when the value is `-` or empty, on both.

**A value with a second `=`.** The same split keeps the text between the first and the second `=` alone (lines 253 and 402, for an INFO item and an `Extra` pair). So `URL=a=b` reads as `a` in `filter_vep`. An HGVSp of a synonymous change in a tab output written without escaping, `HGVSp=ENSP00000369497.3:p.Leu5=`, reads as `ENSP00000369497.3:p.Leu5`. `vep filter` reads everything after the key's first `=`: `URL is a=b` and `HGVSp match Leu5=` pass in `vep filter` and fail in `filter_vep`.

One more difference is observable through `--list` alone. `filter_vep` classes a header line by the number of `#` characters anywhere in it (lines 358-361): exactly two is a meta line, any other count the column line. `vep filter` classes a line by its prefix, `##` a meta line and a single `#` the column line, so an `##INFO` line whose description holds a `#` is listed as a field by `vep filter` and taken as the column line by `filter_vep`. The data lines are filtered the same, since the real column line comes last on both.

The `is_child` operator in `vep filter` answers from the built-in Sequence Ontology table as `--transcript_filter` does, with or without `--ontology`. A parent the table lacks is refused before the input is read, with `ERROR: [filter] no SO term named <term> in the consequence ontology`, where `filter_vep` dies on the first entry its regex does not match. The database connection flags `filter_vep` takes for `--ontology` are refused by name.

**`--only_matched` on a line with two items under the field's name.** Both engines pass a `--keep_csq` output's line when any entry of either item passes. `vep filter` rewrites each item to its own passing entries. `filter_vep` parses its chunk list again at each match, so its hashes hold the first item's entries twice (lines 263-265). It pairs each passing hash with the chunk at its index and rewrites the first item alone, only when the passing hashes and the chunks differ in number (lines 290, 294 and 296). So an entry of the second item can land in the first, a missing chunk is an empty entry, and a line with as many passing hashes as chunks stands as written.

**A line starting with a blank.** `filter_vep` takes the INFO column as the eighth blank-separated word, and Perl's `split /\s+/` keeps an empty first field for a leading blank (line 252). So on such a line it reads the FILTER column's text as the INFO items, while its pattern still finds the consequence entries anywhere in the line. `vep filter` skips the leading blank and reads the INFO column. A leading blank on the first data line makes `filter_vep` read the file as a tab output (lines 425-429), where `vep filter` reads VCF. A blank inside INFO ends the column on both engines; only `filter_vep`'s pattern finds a consequence item past it.

One Ensembl VEP defect in the reading of a VCF line is not reproduced either. `filter_vep` finds a line's consequence entries with one pattern over the whole line. So an empty entry followed by another INFO item or by the FORMAT column takes that text as its entry, and a key ending in the field's name is read as the field. `vep filter` reads every INFO item whose key is exactly the field's name and whose value holds an entry. Both shapes are on [The consequence entries `filter_vep` reads from a VCF line](filter-vep-consequence-field-pattern.md).

A second Ensembl VEP defect in the same reading is not reproduced. `filter_vep` builds each entry's hash over the line's INFO keys and leaves a field the entry lacks at the INFO value of the same name, so a `--keep_csq` output's shorter entries read the input's `AF` as the consequence field `AF`. `vep filter` reads such a field as undefined, on [The field a short consequence entry lacks in `filter_vep`](filter-vep-short-entry-info-key.md).

A third, in the grammar `--transcript_filter` and both tools share, is the hash-order reading of an `in` list mixing a number and other text, on [An `in` list mixing a number and other text](transcript-filter-in-list-numeric-by-hash-order.md). Write the list's entries in one shape, `HGNC:2559,HGNC:1132` or `2559,1132`, and the engines keep the same rows.

Three more lie in the reading of a unit and of an `in` list file. A unit whose first word is `0`, which Ensembl VEP drops, is on [A filter unit whose first word is `0`](transcript-filter-zero-field.md). A list file whose name carries a mode character is on [An `in` list file whose name carries a mode character](transcript-filter-list-file-two-argument-open.md), and a list that names a directory on [An `in` list that names a directory](transcript-filter-list-directory.md).
