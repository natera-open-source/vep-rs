# Index of the fixture page

Prose before any section, with a `key_in_prose` that names nothing.

## Coordinates {#coordinates}

Opening prose of the first section, with [`linked_here`](x.md) and `bare_here`.

| Entry | In one sentence | Tier and count | Verified against | Page |
| --- | --- | --- | --- | --- |
| [`first_entry`](divergences/first.md) | A sentence with an escaped \| pipe | 1 · 7.3 per million · no flag needed | 116.2 | [Page](divergences/first.md) |
| [`second_entry`](divergences/second.md) | Another | 2 · 156.8 per million · 0 flag sets, 28 corpora | 116.2 | [Page](divergences/second.md) |
   [`third_entry`](divergences/third.md) | No outer pipes, indented | 3 · none on the cuts · 9 flag sets, 2 corpora | 116.1 | [Page](divergences/third.md)

Prose after the table.

#### A deeper heading outside the malformed-input section is body text

| Entry | In one sentence | Tier and count | Verified against | Page |
|:---|---:|:-:|---|---|
| [`fourth_entry`](x.md) | Under the same section | tier 2 | 116.2 | [P](x.md) |

## Plugins

### Several plugins

  | Entry | In one sentence | Tier and count | Verified against | Page |
  | --- | --- | --- | --- | --- |
  | [`plugin_entry`](x.md) | s | 3 · none on the cuts · 1 flag set, 0 corpora | 116.2 | [P](x.md) |

### CADD {#cadd}

| Entry | In one sentence | Tier and count | Verified against | Page |
| --- | --- | --- | --- | --- |
| [`cadd_entry`](x.md) | s | 2 · 98.6 per million by output rows · 15 flag sets, 20 corpora | 116.2 | [P](x.md) |
### A heading ends a table

    | Entry | In one sentence | Tier and count | Verified against | Page |
    | [`not_a_table`](x.md) | indented four blanks | 3 | 116.2 | [P](x.md) |

| Entry | In one sentence | Tier and count | Page |
| --- | --- | --- | --- |
| [`other_header`](x.md) | s | 3 | [P](x.md) |

## Differences only malformed input reaches {#malformed}

The opening names [`first_entry`](divergences/first.md), `second_entry`, [`unknown_entry`](x.md) and `broken_entry`; a `Not_A_Key` and [`open_link`] are no keys.

### Coordinates

A line that names `inside_category` is not the opening.

| Entry | In one sentence | Tier and count | Verified against | Page |
| --- | --- | --- | --- | --- |
| [`broken_entry`](x.md) | s | 3 · none on the cuts · 1 flag set, 0 corpora | 116.2 | [P](x.md) |

### Plugins

#### CADD

| Entry | In one sentence | Tier and count | Verified against | Page |
| --- | --- | --- | --- | --- |
| [`broken_cadd_entry`](x.md) | s | 3 · none on the cuts · 0 flag sets, 0 corpora | 116.2 | [P](x.md) |

## How to read the counts

Prose with `a_key` after the tables.
