# A PubMed or synonyms list the cache cut inside an id

**The defect in one sentence.** The Ensembl caches of release 116 end a known variant's `pubmed` and `var_synonyms` lists at 1,024 characters, in the middle of an id as readily as after a separator, and Ensembl VEP 116.3 prints the cut piece as an id. vep-rs writes the list without the item the cache cut.

**Why this is a defect.** An identifier is whole or it is not that identifier. The first seven digits of an eight-digit PubMed id are a number PubMed may have given to another article, and the first characters of a ClinVar accession name no record. The output marks the piece in no way: it stands last in the list, separated like every id before it. A reader who joins the printed ids to their sources attaches a citation or a record the database never linked to the variant.

**What vep-rs does instead.** vep-rs writes the list without the item the cache cut, in every output format. A stored value of exactly 1,024 characters that ends in no separator is read as the text before its last separator. A list at the cap whose last item cannot be told from a whole one loses that item too. In a synonyms list a source whose only id is the cut one goes with that id. A value ending in a separator, a shorter value and a longer one are written as stored.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, under `Bio/EnsEMBL/VEP/`):

- `Pipeline/DumpVEP/InitDump.pm`, lines 194-201 and 257-269: the dump's two queries build each list with `GROUP_CONCAT` and set no length limit, so the database server's limit applies.
- `Pipeline/DumpVEP/Dumper/Variation.pm`, lines 252-253: the two strings are written into the cache's columns as the queries returned them.
- `OutputFactory.pm`, lines 1126, 1090 and 1141: the stored `pubmed` value is split on `,` into `PUBMED`, and the stored `var_synonyms` strings are joined with `--`.
- `OutputFactory/JSON.pm`, lines 452-457 and 469: the same stored values are split into the `var_synonyms` object and the `pubmed` array.

**What correcting it would change.** A dump whose queries lifted the server's limit would store every list whole, both engines would print the same lists from its caches, and the class would vanish. Against the caches of release 116 no writer can recover a cut item, because the cache does not hold the rest of it.

**The class it produces.** `colocated_list_cut_at_cache_cap`: a `PUBMED` or `VAR_SYNONYMS` value, or a JSON `pubmed` or `var_synonyms` member, read from a stored list of exactly 1,024 characters that ends in no separator. Ensembl VEP 116.3 writes the stored list and vep-rs writes it without its last item.

**How often the shape occurs.** The table counts the stored lists of exactly 1,024 characters over every variation file of the human caches of release 116, 1,487,415,922 known variant rows on GRCh38 and 1,129,087,916 on GRCh37. No stored value of either field is longer than 1,024 characters in either build, and none has 1,023 or 1,025.

| Stored lists of exactly 1,024 characters | GRCh38 | GRCh37 |
| --- | --- | --- |
| `pubmed`, all | 1,122 | 752 |
| `pubmed` ending in `,`, every id whole | 43 | 26 |
| `pubmed` ending inside an item | 1,079 | 726 |
| of those, a last item shorter than every other item | 1,007 | 685 |
| of those, a last item within the length range of the other items | 72 | 41 |
| `var_synonyms`, all | 25 | 6 |
| `var_synonyms` ending in `--`, every id whole | 3 | 0 |
| `var_synonyms` ending inside an id that follows a `,` | 19 | 6 |
| `var_synonyms` ending inside the only id of its source group | 3 | 0 |

## 1. A PubMed list of 113 ids and seven digits of another

**In short.** A known variant's stored `pubmed` value holds 113 ids and the first seven digits of a 114th, 1,024 characters in all. Ensembl VEP 116.3 prints the seven digits as a 114th id in every output format, and vep-rs writes the 113 whole ids. A `var_synonyms` value is read the same way: the cut id goes and the ids before it stay, and a source whose only id is the cut one goes with it. A list of 1,024 characters whose last item is as long as a whole id loses that item as well, because nothing in the stored value tells such an item from the first characters of a longer id.

### 1.1 What a reader needs to know first

`--pubmed` adds the PubMed ids of the known variants a row lists, and `--var_synonyms` adds their synonyms; each switches `--check_existing` on. Both fields are read from the cache's variation files. A `pubmed` cell is the variant's ids joined by `,`. A `var_synonyms` cell is groups joined by `--`, each a source name, `::` and that source's ids joined by `,`. [Co-located variant keys](../output/known-variants.md#co-located-variant-keys) describes how the fields are printed.

The cache dump builds each list in the database, and the database ends a list at 1,024 characters. A stored value of exactly that length either ends in a separator, and then every item in it is whole, or ends inside an item. Every stored value of the two fields is ASCII in both builds, so characters and bytes count the same; vep-rs counts bytes.

An item is the smallest piece the engines print as one identifier: a PubMed id, or one id of a synonyms group.

### 1.2 The record

Every identifier on this page is invented, and the page prints no value of a release cache. The input is one substitution:

```
21	500	.	A	G	.	.	.
```

The cache holds one known variant at the position, `rs9920000001` with the alleles `A/G`. Its `pubmed` cell is the ids `99000001` to `99000113` joined by commas, then a comma and `9900011`:

```
99000001,99000002,...,99000113,9900011
```

The `...` stands for the 110 ids `99000003` to `99000112`. The 113 ids and their commas fill 1,017 characters, and the seven digits bring the value to 1,024. Both engines run with `--check_existing --pubmed --var_synonyms`. The record and its value are those of the test `pubmed_cut_inside_an_id_ends_before_it` in `crates/vep-cli/tests/colocated_list_fields.rs`.

The GRCh38 cache of release 116 stores a list of this shape for the known variant starting at chromosome 16 position 27344882. Its `pubmed` value has 1,024 characters and 115 items: 111 of eight characters, three of seven and a last one of one character.

### 1.3 What each engine writes

Ensembl VEP 116.3 writes the stored value whole. In the default format it is the `PUBMED` pair of the Extra column, and the `--tab` cell holds the same list:

```
PUBMED=99000001,99000002,...,99000113,9900011
```

vep-rs writes the 113 ids:

```
PUBMED=99000001,99000002,...,99000113
```

The VCF format joins the same items with `&` in both engines. In JSON the list is an array of numbers, 114 of them in Ensembl VEP's line, the first below, and 113 in vep-rs's, the second:

```
"pubmed":[99000001,99000002,...,99000113,9900011]
"pubmed":[99000001,99000002,...,99000113]
```

The other stored values of the test file, each of 1,024 characters, are written as follows. A `var_synonyms` value there is one `UniProt` id padded to the length the value needs, then the 70 `ClinVar` ids `RCV999000001` to `RCV999000070`, then the ending the row names.

| Stored value | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `pubmed`: 112 ids, `9900113`, then `99000114` | 114 items, the last `99000114` | 113 items, the last `9900113` |
| `pubmed` of the record, beside a second known variant storing `99000201` | 115 items, `9900011` the 114th | 114 items |
| `var_synonyms` ending `,RCV999` | the stored text; JSON lists `RCV999` as a 71st `ClinVar` id | the text before `,RCV999`; JSON lists 70 ids |
| `var_synonyms` ending `,RCV999000999` | the stored text, 71 `ClinVar` ids | the text before `,RCV999000999`, 70 ids |
| `var_synonyms` ending `--OMIM::6059` | the stored text; JSON lists `6059` under `OMIM` | the text before `--OMIM::6059`; JSON has no `OMIM` key |
| `pubmed` ending in `,` after 114 ids | the 114 ids | the 114 ids |

### 1.4 Why Ensembl VEP writes that

Two steps produce the row, and neither checks the other. The cache dump asks the database for each variant's ids as one concatenated string and sets no limit on its length. The database server's limit for such a string then applies, and the string ends at its 1,024th character. The dump writes the string into the cache column as it was returned.

The annotator reads the stored value as a list of whole items. For the default, tab and VCF formats it splits `pubmed` on `,` and joins the pieces again, and it pushes the `var_synonyms` string whole. For JSON it splits both values on their separators. No step tests the length of the value or of an item, so the piece after the last separator is printed as an id.

Source (`ensembl-vep` at release/116.3, under `Bio/EnsEMBL/VEP/`):

- `Pipeline/DumpVEP/InitDump.pm`, lines 194-201: the PubMed query, `GROUP_CONCAT(p.pmid)` for each variant; lines 257-269: the synonyms query, its groups joined by `--`.
- `Pipeline/DumpVEP/Dumper/Variation.pm`, lines 252-253: the two strings written into the columns.
- `OutputFactory.pm`, line 1126: `split(',', $ex->{pubmed})` pushed onto `PUBMED`; line 1090: the synonyms string pushed whole; line 1141: the join with `--`.
- `OutputFactory/VEP_output.pm`, line 117, `OutputFactory/Tab.pm`, line 128, and `OutputFactory/VCF.pm`, lines 391 and 401: the list joined with `,`, or with `&` in a VCF.
- `OutputFactory/JSON.pm`, lines 452-457: the synonyms split on `--`, `::` and `,`; line 469: `pubmed` split on `,`.

### 1.5 Why vep-rs's row is the correct one

The stored value says where it was cut. It has exactly 1,024 characters, the length no value of the field exceeds, and it ends in no separator. Its last item therefore ends where the limit fell, not where an id ends. Here that item has seven digits where the 113 before it have eight, so `9900011` is the start of an id whose last digit the cache does not hold.

The caches show the same on their own lists. Of the GRCh38 `pubmed` values of 1,024 characters that end inside an item, 1,007 of 1,079 have a last item shorter than every other item of the list; on GRCh37, 685 of 726 do. The lists just under the cap show no such item: of the `pubmed` lists of 974 to 1,023 characters, the fifty lengths below 1,024, none of 72 on GRCh38 and none of 42 on GRCh37 has one. Those 72 lists lie under the cap and are not the 72 of the table, which stand at it. A list of eight-digit ids reaches 1,017 characters after 113 ids, so the cut falls after seven digits of the 114th.

A list at the cap whose last item cannot be told from a whole one loses that item too. Such an item is a whole id or the first digits of a longer one, and the stored value does not say which. On GRCh38 72 `pubmed` lists end this way, and on GRCh37 41. vep-rs prints no id it cannot know to be whole, and on those lists the price is one id that may have been whole.

A cut id that is the only id of its source group takes its source with it, in the row formats and in JSON alike. vep-rs writes no source name without an id and no source with an empty list: Ensembl VEP never prints a source without an id from a whole list, and a reader of the field expects the pair. Three GRCh38 `var_synonyms` values have this shape and no GRCh37 value does. A value left with no whole id writes no field; no cache of release 116 stores one.

The cache holds no more than 1,024 characters of a list, so neither engine can print the rest. vep-rs leaves the cache as it is: a converted cache carries the stored value whole, and the reading is applied as a row is written. A cache `vep-cache-builder` installs stores such a list at its whole items within the length. A reader who needs the stored value reads the `pubmed` or `var_synonyms` column of the variant's row in the cache's variation file. A reader who needs the whole list reads the source database: Ensembl's variation database, or the database a synonym's source name refers to.

### 1.6 Try it yourself

The test file builds every stored value of this page in a cache of its own and holds vep-rs to the rows above:

```
cargo test -p vep-cli --test colocated_list_fields
```

This command prints the record's `pubmed` value, 1,024 characters and a line end:

```
printf '%s,9900011\n' "$(seq -s, 99000001 99000113)"
```

For Ensembl VEP 116.3, copy a GRCh38 cache of release 116 and replace `21/all_vars.gz` with a file of one row. The row is tab-separated in the column order that `variation_cols` gives in `info.txt`. It holds `21` under `chr`, `rs9920000001` under `variation_name`, `500` under `start`, `A/G` under `allele_string`, the value under `pubmed` and `.` in every other cell. Compress the file with `bgzip` and index it with `tabix -C -s 1 -b 5 -e 5`. Then run the record with the cache options of [Reproducing these examples](reproducing.md):

```
vep -i record.vcf -o out.txt --offline --cache --dir_cache <copy> \
  --cache_version 116 --species homo_sapiens --assembly GRCh38 \
  --format vcf --check_existing --pubmed --var_synonyms
```

| Field of the record's row | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `PUBMED` items | 114, the last `9900011` | 113, the last `99000113` |
| `PUBMED` length | 1,024 characters | 1,016 characters |
| JSON `pubmed` elements | 114 | 113 |

The corpus `tests/golden/116/GRCh38-colocated-lists` holds the control. The known variant of its first record stores 113 ids in 1,024 characters ending in a comma: 106 ids of eight digits and seven of nine, 911 digits in all, and 113 commas. Both engines write the 113 ids.

### 1.7 How the comparator treats this class

The format comparator carries the class `colocated_list_cut_at_cache_cap`. It sets a pair aside only when two things hold, which JSON reads from its object in the ways the later paragraphs of this section give. Ensembl VEP's list equals vep-rs's list plus one last item, which in a synonyms list may be a last source group of one id. And Ensembl VEP's printed value shows a stored list of exactly 1,024 characters that ends in no separator.

A list one character shorter or longer, a pair whose earlier item differs, a missing value and a wrong value all stay differences. The one exception is a JSON list whose last id is printed as a number with no decimal part, which the end of this section describes. The class applies only where a comparison case names it among its `masks`, the classes the case sets aside, and no dataset run names it.

In the default, tab and VCF formats a row prints the stored lists of all its known variants joined into one value. The separator that joins them, `,` in `PUBMED` (`&` in a VCF) and `--` in `VAR_SYNONYMS`, is also a separator inside a list, so the value does not show where one list ends. There the comparator takes as one stored list each stretch of exactly 1,024 characters that starts at the value's start or after any such separator and ends at the value's end or before one. A `%3B`, which the default and VCF formats print in place of a `;`, counts as that one character.

vep-rs's value must then equal Ensembl VEP's with that stretch's cut last item taken out and everything before and after it unchanged, a later stretch being read the same way. In JSON it reads one known variant at a time.

In JSON, whose object keys Ensembl VEP writes in hash order, no group reads as the last. So there the one id vep-rs lacks may stand in any group, as that group's last id, and where it is the group's only id the group goes with it. Such a pair is set aside only when the stored value rebuilt from Ensembl VEP's object has exactly 1,024 characters. In the row formats that id must stand in the last group of its stretch.

Ensembl VEP's JSON prints an id that Perl's number test admits as a number (`Utils.pm`, line 241, called from `OutputFactory/JSON.pm`, line 288). A cut piece of digits and a closing point, `999001.`, is therefore printed as `999001`, one character short of the stored piece. One of the three GRCh38 values that end inside the only id of a source group closes on such a point.

Where the id vep-rs lacks is printed as a number with no decimal part, the object does not show whether a point closed the stored piece. The comparator then rebuilds the value twice, as the object gives it and with a point after that number, and sets the pair aside when either text has exactly 1,024 characters. The row formats print the stored text, so on the same record the comparator still requires a stored list of exactly 1,024 characters there.

Two other stored values print the same object as that cut one, though every id in them is whole: a value of 1,023 characters that ends in such a number, and a `var_synonyms` value of 1,024 characters that ends in a comma directly after it. If vep-rs left that last id out of either, JSON would set the pair aside all the same. No stored value of either build has 1,023 characters, and no stored `var_synonyms` value at the cap ends in a comma.

The comparator adds nothing else to a rebuilt value. A cut piece that closes on a zero after its point, `999001.0`, is printed as `999001`, two characters short, and a value that ends in it stays a difference in JSON. An id elsewhere in the list that closes on a zero, `999001.0010`, is printed as `999001.001` and leaves the rebuilt value one character short. Such a pair is set aside only where one of the two texts still has exactly 1,024 characters.
