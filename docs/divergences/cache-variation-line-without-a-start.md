# Known variants behind a cache variation line with no start

**The defect in one sentence.** Under `--check_existing`, a line of a cache's `all_vars.gz` whose start cell holds no number ends Ensembl VEP 116.3's lookup of every record that reads the line first, so those records are written with no known variant. vep-rs passes over the line and names the known variants the file holds.

**Why this is a defect.** The line that cannot be read is not the record's. The known variant at the record's position stands further on in the same file, well-formed, and Ensembl VEP names the known variants of that file wherever a lookup does not start at the line. The run exits 0 and writes no warnings file; two lines of the tabix library on standard error are the only sign, and a run under `--fork` does not print them. A reader of the output takes the empty `Existing_variation` cell for a variant no database knows.

**What vep-rs does instead.** vep-rs reads the file through its own index reader. A line whose start cell no integer leads names no position, so it is no known variant: vep-rs passes over it and reads on. Each record is written with the known variants at its position, and the rows are those of the same file without the line.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, under `Bio/EnsEMBL/VEP/`, and the libraries its installer builds):

- `AnnotationSource/Cache/VariationTabix.pm`, lines 224-229: one tabix query for each record, from one base before its start to one base after its end. Line 232: the loop reads lines until the iterator returns none, and no line tests why it stopped.
- `OutputFactory.pm`, line 1046: every field of a row that comes from a known variant is added only when that loop found one.
- `Bio/DB/HTS.xs` of Bio::DB::HTS 2.11, lines 1453-1456: any negative return of the library's iterator comes back as no line.
- `Runner.pm`, lines 658-662: a forked child closes its standard error, so the library's message is lost under `--fork`.
- `INSTALL.pl`, line 104, and `travisci/get_dependencies.sh`, line 10, which the release's Docker image, `ensemblorg/ensembl-vep:release_116.3`, is built with: both fetch htslib 1.9.
- htslib 1.9, `tbx.c`, lines 95-96: the start cell is read as a number, and a cell with no digit fails the parse. Line 161 prints the message. Lines 174-175: `tbx_readrec` hands on the fields of the line it failed to parse.
- htslib 1.9, `hts.c`, lines 2643-2644: the iterator ends the read at a line whose sequence is not the query's.

**What correcting it would change.** A reader that passed over a line it cannot parse would return the lines after it, and both engines would name the same known variants. From htslib 1.10 the library's index builder passes over such a line; its reader, in every release, does not. A cache dump that wrote no such line would never meet the difference.

**The class it produces.** `cache_variation_line_without_a_start`: the `Existing_variation` cell, and the fields read from a known variant, on the rows of a record whose lookup reads a line with no start before its own lines. A lookup reads the file from the start of the record's index bin, so these are the records of the bin whose stretch of the file begins at the line. Ensembl VEP 116.3 writes no known variant and vep-rs writes the ones the file holds. The example and the cache counted below hold the line as the first of a file; a line further into a file was not run.

**How often the shape occurs.** The table counts the lines whose start cell is not a whole number over the variation files of three caches of release 116. The third is the cache Ensembl Genomes publishes ready-made for the mosquito *Anopheles funestus*, assembly AfunF3, cache version 63.

| Cache | Files | Lines | Lines with no start |
| --- | --- | --- | --- |
| Human, GRCh38 | 1,904 | 1,487,415,922 | 0 |
| Human, GRCh37 | 699 | 1,129,087,916 | 0 |
| *Anopheles funestus*, AfunF3 | 3 | 14,489,964 | 1 |

The one line is the first of the file of chromosome 2, which the cache keeps in the directory named by the accession its `chr_synonyms.txt` pairs with `2`. The index bin that begins at it covers the chromosome's first 16,384 bases and holds 862 known variants, none of which Ensembl VEP 116.3 names. The human counts are of each assembly's cache of the Ensembl transcript set; the RefSeq set's own copies of the files hold the same lines. No other cache was counted.

## 1. A first line whose start cell is `.`

**In short.** The first line of a chromosome's variation file has one cell too many, so its start cell is `.`. Two records lie in the index bin that begins at that line and one lies in the next bin. Ensembl VEP 116.3 writes the first two with no known variant and names the third's. vep-rs names all three.

### 1.1 What a reader needs to know first

`--check_existing` adds the known variants at a record's position to its rows, in the `Existing_variation` column. A cache converted with `convert_cache.pl` holds them in one `all_vars.gz` for each chromosome: a bgzip file of tab-separated lines sorted by start, with a tabix index on the start column. The columns are those `variation_cols` names in the cache's `info.txt`, the chromosome first and the start fifth.

The index files the lines into bins of 16,384 bases and records, for each bin, the stretch of the file that holds its lines. A lookup for a position seeks to the start of its bin's stretch and reads line by line until a line starts past the position. Every line read on the way is parsed for its sequence and start.

The `tabix` of htslib 1.10 or later builds an index past a line it cannot parse: it prints a message and files nothing for the line. A bin's stretch begins where the line filed before it ended, so the stretch of the next line's bin begins at the unparsed line.

### 1.2 The record

Every identifier in this example is invented. The variation file of chromosome 21 holds five lines; the first eight cells of each are shown, and every cell after them is `.`:

```
21	var_fi21	var_first	.	.	84	.	A/G
21	var_first	.	.	84	.	A/G	.
21	var_second	.	.	6052	.	G/C	.
21	var_second	.	.	6052	.	G/C	.
21	var_third	.	.	16419	.	T/C	.
```

The second cell of the first line, `var_fi21`, is the first six characters of a name followed by the chromosome's name. A whole record follows it, one cell to the right of its place, so the line's fifth cell is `.` where the others hold a start. The second known variant stands on two lines, as the cache below writes many of its variants. The input is three substitutions, each at the start of one known variant with its alleles:

```
21	84	at_the_line	A	G	.	.	.
21	6052	past_the_line	G	C	.	.	.
21	16419	past_the_bin	T	C	.	.	.
```

Both engines run with `--check_existing`. The records and the file are those of the tests in `crates/vep-cli/tests/known_variants_past_an_unreadable_line.rs`.

The *Anopheles funestus* cache above has this shape. The first line of its file for chromosome 2 has sixteen cells where every other line has fifteen, and its start cell is `.`. The cache's own index files nothing for the line and begins its first bin's stretch at it. For records at positions 84, 6052, 11759 and 15468 of that chromosome, each the start of a known variant of the file, Ensembl VEP 116.3 names no known variant and vep-rs names the file's variant. For a record at 16419, in the next bin, both engines name it.

### 1.3 What each engine writes

The `Existing_variation` cell of each record's row:

| Record | Location | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `at_the_line` | `21:84` | `-` | `var_first` |
| `past_the_line` | `21:6052` | `-` | `var_second` |
| `past_the_bin` | `21:16419` | `var_third` | `var_third` |

Every other cell of the three rows is the same in both engines: the made-up variants hold a name and alleles alone, so no field read from a known variant differs beside the name. Ensembl VEP exits 0 and writes no warnings file. Run without `--fork`, its standard error holds two lines for each of the first two records. The second quotes the file's first line whole and is cut here at `...`:

```
[E::get_intv] Failed to parse TBX_GENERIC, was wrong -p [type] used?
The offending line was: "21	var_fi21	var_first	.	.	84	...
```

### 1.4 Why Ensembl VEP writes that

Ensembl VEP asks the tabix library for the lines from one base before a record's start to one base after its end, and reads what the library's iterator returns until it returns nothing. It makes one such query for each record.

For the first two records the query lies in the first bin, whose stretch begins at the file's first line. The library parses that line, finds no number in its start cell, prints its message and fails the parse. The reader of htslib 1.9 does not test that failure: it hands the iterator the line's fields, among them a sequence number the failed parse never set. The number is whatever the memory held, so the source alone does not fix what follows. In every run made in the release's Docker image the iterator then ended the read, as it does at a line of another sequence, where a sorted file holds nothing more for the query.

The Perl wrapper returns nothing for any negative result, and Ensembl VEP's loop ends as it does when a position has no known variant. The well-formed lines at the record's position, a few lines further on, are never read. From htslib 1.10 the reader prints the message and returns an error for such a line, which the wrapper returns as nothing in the same way.

The third record's query lies in the second bin, whose stretch begins at the line of `var_third`, so the first line is never read. Source: the lines listed under "Where it is in Ensembl VEP's source".

### 1.5 Why vep-rs's row is the correct one

Ensembl VEP's own rule says which known variants a record has. They are the file's lines from one base before the record's start to one base after its end (`VariationTabix.pm`, line 229). A line counts when it is not marked as failed, the cache's flag for a variant Ensembl's quality checks set aside (`AnnotationType/Variation.pm`, lines 224-227). It must also share an alternate allele with the record at its position (`AnnotationType/Variation.pm`, lines 146-207).

The record at 6052 substitutes `C` for `G`. The file's third and fourth lines start at 6052, are not failed and hold the alleles `G/C`. The rule keeps and matches them, and a name is listed once (`VariationTabix.pm`, line 237). The record has one known variant, `var_second`. The same reading gives `var_first` for the record at 84, from the second line, and `var_third` for the record at 16419.

The first line adds nothing to any record. Its start cell is `.`, which Ensembl VEP's own parser of a variation line reads as no start (`AnnotationSource/Cache/BaseCacheVariation.pm`, line 124), so it names no position a record could share. Passing over it loses nothing. Ending the read at it loses every known variant of the bin.

Each engine also shows the correct row itself. Ensembl VEP names `var_third` from the same file, by the same code, where the lookup does not start at the line. And vep-rs writes the same three rows from the file without its first line.

### 1.6 Try it yourself

The test file builds the cache of this page and holds vep-rs to the rows above, with the line and without it:

```
cargo test -p vep-cli --test known_variants_past_an_unreadable_line
```

For Ensembl VEP 116.3, copy a GRCh38 cache of release 116 and replace `21/all_vars.gz` with a file of the five lines of section 1.2. Lines two to five are tab-separated in the column order that `variation_cols` gives in `info.txt`: `21` under `chr`, the name under `variation_name`, the position under `start`, the alleles under `allele_string` and `.` in every other cell. The first line is the second with `var_fi21` and a tab put before its name.

Compress the file with `bgzip` and index it with `tabix -C -s 1 -b 5 -e 5` of htslib 1.10 or later, which prints the message once and writes the index; `-C` writes a `.csi`, the index the caches carry. The `tabix` of htslib 1.9 files a line it cannot parse under a sequence no line names, and a lookup then never reads it. Then run the three records, under a VCF header, with the cache options of [Reproducing these examples](reproducing.md):

```
vep -i records.vcf -o out.txt --offline --cache --dir_cache <copy> \
  --cache_version 116 --species homo_sapiens --assembly GRCh38 \
  --format vcf --check_existing
```

For vep-rs, make a cache directory of three parts: the same `all_vars.gz` and its `.csi` under `variations/21/`, a file `transcripts/21/1-1000000.json` holding `[]`, and an `info.json` holding `assembly` `GRCh38`, `var_type` `tabix` and the same `variation_cols` as a list. Then run:

```
vep -i records.vcf -o out.txt --offline --json_cache <directory> \
  --assembly GRCh38 --check_existing
```

### 1.7 How the comparators treat this class

No comparator sets these rows aside. The comparators are the project's tools that compare the two engines' outputs, and each sets aside the rows of some documented differences under a class name. A row on which Ensembl VEP names no known variant and vep-rs names one looks, from the two rows alone, like a known variant vep-rs made up, and only the variation file tells the two apart. Such a row therefore stays a difference in every comparison and is traced to the file.

The human caches of release 116 hold no line of this shape, so no comparison over human records meets the class.
