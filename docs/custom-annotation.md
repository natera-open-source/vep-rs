# Custom annotation

`--custom` names a tabix-indexed VCF, BED or GFF file, or a bigWig, whose records annotate the input variants under a short name of the source's choosing. The string is read as Ensembl VEP reads it, in either of its two forms:

```text
--custom file=<file>,short_name=<name>,format=<vcf|bed|bigwig|gff|gtf>,
         type=<overlap|within|surrounding|exact>,coords=<0|1>,fields=<A%B>,
         overlap_cutoff=<pct>,distance=<bp>,same_type=<0|1>,reciprocal=<0|1>,
         num_records=<N|all>,summary_stats=<min%max%mean%count%sum>,
         gff_type=<transcript|gencode_promoter>
--custom <file>,<short_name>,<format>,<type>,<coords>,<field>,<field>,...
```

The keyed form is written on one line; it is broken here to fit the page.

## Defaults

Every option takes Ensembl VEP's default when the string leaves it out.

| Option | Default |
| --- | --- |
| `type` | `overlap` |
| `coords` | off |
| `overlap_cutoff` | off |
| `same_type` | off |
| `reciprocal` | off |
| `distance` | unset |
| `num_records` | 50 |
| `gff_type` | `transcript` |

The short name is the file's last path segment, or `GENCODE_PROMOTER` for a `gencode_promoter` GFF.

## One source per chromosome

A file holding `###CHR###` stands for one source per chromosome whose file exists. The chromosomes are the reference FASTA's sequence names when the run has one, else `1` to `22`, `X`, `Y` and `MT`. A pattern that matches no file exits 1 with `No files with pattern <file> were found`, whatever its position among the sources.

Ensembl VEP 116.3 dies with that text only when the pattern is the run's first `--custom` source. After an earlier source it keeps the sources it has and drops the pattern without a key, a header line or a warning. The test guarding the message counts every source the run has built so far rather than the pattern's own (`AnnotationSourceAdaptor.pm` lines 330 and 451-452). No committed run names such a pattern (tier 3 of [Intended divergences from Ensembl VEP](intended-divergences.md); 19 flag sets run a `--custom` source).

## Strings and files that are refused

A malformed string exits 1 before any output is written, with Ensembl VEP's text under the `custom` component. The string is quoted on a `LINE:` line where Ensembl VEP does. Each of these is refused so:

- an unknown key
- a missing file or format
- an invalid type
- an unknown format
- an unsupported summary statistic

Every source's file is opened before the input is read, under `--show_cache_info` too. A file exits 1 the same way, as Ensembl VEP fails before its output opens, when it:

- is not there
- has no index
- is not bgzipped
- has a `##fileformat` line that does not name major version 4

Whether a file is bgzipped is read from its bytes, never its name. A bgzipped and indexed source named without `.gz` or `.bgz` is read, where Ensembl VEP refuses it by name alone, the documented class [`tabix_source_named_without_gz`](divergences/tabix-source-named-without-gz.md).

## Sources that are not read

A remote file (`http://`, `https://`, `ftp://` or any other URL scheme) is not read and is refused with an error that names the route. Ensembl VEP reads it over the network; vep-rs reads local files only, so copy the file to a local path.

## How a source is read

A VCF, BED or GFF source is read record by record against each input record, as Ensembl VEP reads it. The records a batch can reach are read once. A window one base either side of every input record is merged into as few reads as the file's layout allows. Each input record then takes the records of its own window in file order.

A `--custom` run reads its cache as every run does; Ensembl VEP accepts a `--custom` source as a run's only annotation source. Files are read locally through their tabix index (`.tbi` or `.csi`); a BAM is not read, and a remote file is refused as above.

## Chromosomes

A source lists the chromosomes its index names (a bigWig, its chromosome list). An input chromosome is matched to them as to the cache's, through the synonyms and the `chr` rule of [Chromosome names](cache-setup.md#chromosome-names).

An input record on a chromosome the cache lacks is admitted when a source, or a `--regulatory_gff` file, lists it. It is written as `intergenic_variant` with the sources' values: Ensembl VEP admits a chromosome any annotation source holds (`BaseRunner.pm` lines 444-460). A bare `M` is renamed `MT` when a source or the cache lists `MT` and nothing lists `M` or `chrM` (`Parser.pm` lines 794-796).

## Matching

A BED record spans `chromStart + 1` to `chromEnd`. A structural VCF record (its ALT symbolic without `*`, a bracketed breakend, or not bases alone beside an `SVTYPE`) spans:

- `POS + 1` to `END`
- or to `POS + |SVLEN|` with `SVLEN`
- or `POS + 1` alone with neither

A record whose ALT names `INS` ends at `POS + length(REF)`. Against those spans, and against a structural input record whatever the record, the four types read as Ensembl VEP's.

| Type | Asks for |
| --- | --- |
| `overlap` | a shared base |
| `within` | a record inside the input record |
| `surrounding` | a record holding the input record |
| `exact` | the same start and the same end |

`distance` asks each end of the record to lie within that many bases of the input record's end. `overlap_cutoff` asks the shared bases to be at least that percentage of the record, or of the input record under `reciprocal=1`. `PC` among the `fields` reports that percentage to three decimals.

A sequence VCF record, its ALT bases alone, spans `POS` to `POS + length(REF) - 1`. It matches an `overlap` source on a shared base. Under `exact`, each alternate allele of the record is compared with each of the input record's. Both pairs are trimmed to their minimal form and position, and the record attaches to the alleles whose forms agree. A sequence record never matches a `within` or `surrounding` source, as in Ensembl VEP, whose VCF reader has branches for `overlap` and `exact` alone.

`same_type` admits a record only when its structural class is the input record's, so it matches structural records alone.

## Values

A record is named by its first `ID` (VCF) or its name column (BED), or by `chromosome:start-end` under `coords=1` or when it has no identifier.

The `fields` of a VCF source are its INFO keys, a key without a value reading `1`. A comma-separated value with one entry per alternate allele gives each allele its own entry under `exact`. Under any other type it rejoins the entries in allele order. A value with one entry more drops its first. `--custom_multi_allelic`, and a file whose header says `##source=ClinVar`, keep every value whole on every allele.

A VCF source with `fields` adds the record's `FILTER` entries, split at `;`, unless `--custom_suppress_filter` is given and `fields` does not name `FILTER`. A record's score is its `QUAL` (VCF) or its fifth column (BED).

## Records and statistics per input record

An input record keeps at most `num_records` matching records and, where more matched, one marker record named `...`. The marker belongs to the whole input record and to no allele, which decides [where a row lists it](#the-rows-a-record-prints-on). `num_records=all` keeps every one.

A `num_records` with a fraction or a minus sign (`2.5`, `-1`) is read as its integer part, at least 0. The source writes one warning naming the value and the reading. Ensembl VEP compares the count with the text as a number ([Defect N, class 26](divergences/defect-n.md#26-a---custom-record-limit-that-is-not-a-whole-number)).

`summary_stats` computes `min`, `max`, `mean`, `count` and `sum` over the matching records' scores:

- `count` counts the records
- `min` and `max` keep a score of 0
- `mean` is the sum over the count
- a record without a score counts but adds no value

A `COSV` record can be named like one already kept under the same short name, by this source or an earlier one. Such a record lends the kept one the fields it lacks and is not kept twice. It counts, adds no score and takes no `PC`.

Sources that share a short name list their records under its keys. The name's statistics come from the sources that asked for some. A source that asked for none leaves them as they stand, as in Ensembl VEP: the keys carry the asking source's figures, a `count` of 0 where it asked for `count` and matched nothing. Where two sources ask, vep-rs prints one source's own figures ([the shared statistics](#2-the-summary-statistics-of-two-custom-sources-sharing-a-short-name)).

The records and statistics belong to the whole input record, so every allele of a multi-allelic input record shares them. An `exact` VCF record belongs to the allele it matched, read as the input gives it. A 3' shift of the allele a transcript row prints leaves the record on its own allele, where Ensembl VEP releases before 116.3 attach it by the rotated allele the row prints.

## GFF and GTF transcript models

A GFF or GTF transcript model (`format=gtf`, or `format=gff` in its default `transcript` mode) is a transcript source read beside the cache, as `--gff <file>` and `--gtf <file>` name it. Its file is read whole, so it needs no index and no bgzip, and its transcripts are annotated with the cache's. [A file beside a cache in one run](cache-from-gtf-or-gff3.md#a-file-beside-a-cache-in-one-run) states what the run needs and writes.

## GFF promoter windows

`format=gff,gff_type=gencode_promoter` reads a tabix-indexed GFF of GENCODE promoter windows as an overlay. It is the one GFF mode whose records annotate the input. A record spans its fourth to its fifth column and matches by the rules of [Matching](#matching); `exact` is refused, as Ensembl VEP refuses it for any GFF source, and so is an index `tabix` built past a start that is no number, [`gxf_malformed_coordinate_line_refused`](divergences/gxf-malformed-coordinate-line-refused.md).

Each matching `promoter_window` record becomes one record named `chromosome:start-end`, a leading `chr` dropped, with three fields: `type` (`promoter_window`), `feature_id` (its `ID` attribute) and `associated_gene` (its `gene_id` attribute, else `.`), every `%XX` escape decoded. The short name defaults to `GENCODE_PROMOTER`, the keys `GENCODE_PROMOTER`, `GENCODE_PROMOTER_type`, `GENCODE_PROMOTER_feature_id` and `GENCODE_PROMOTER_associated_gene`. A window without an `ID` leaves `feature_id` undefined: the text formats print the empty entry, the JSON format `null`. `fields` in the string changes nothing, and the records carry no score, so `summary_stats` yields `count` alone.

A matching record of any other type (an `enhancer`, a `gene`) yields no record. It yields one warning per type and file, written as every warning is. The warning reads `Ignoring '<type>' feature_type from <file> GFF/GTF file. This feature_type is not supported in VEP.`. A GFF in a `gff_type` Ensembl VEP does not know opens and annotates nothing. Ensembl VEP reads a GFF source offline only with `--fasta`; vep-rs needs none for a promoter source.

## bigWig sources

A bigWig source's matching rules, statistics and span rule, and the `--bigwig`, `--phyloP` and `--phastCons` shortcuts, with `--ucsc_assembly`, `--ucsc_data_root` and their local mirror: [bigWig sources and conservation shortcuts](output/custom-bigwig.md).

## Output

Every source adds keys to the row, as Ensembl VEP's output factory adds them.

| Key | Value |
| --- | --- |
| `<short_name>` | the names of the records that print on the row |
| `<short_name>_<field>` | the values of a named field, one entry per record carrying it |
| `<short_name>_<stat>` | a statistic as a single value |

A list joins with `,` in the default and tab formats and with `&` inside a VCF `CSQ` entry. A `FILTER` list joins its entries with `;`, which the default format and the VCF escape as `%3B`. A number prints as Ensembl VEP prints one.

The default format sorts the keys among the other `Extra` keys after the flag fields. The tab, VCF and Parquet formats add one column per header key after the plugin columns, source by source in the order given.

`--custom` also adds the `SOURCE` column. It names the cache a merged cache's transcript came from, `Ensembl` or `RefSeq`, and the short name of the GFF or GTF source a transcript was read from beside the cache. It is empty on a transcript of a single-source cache, a cache built from a GTF or GFF included. Ensembl VEP's `--gtf` and `--gff` sources stamp their transcripts with the file's name, which a cache built from that file does not carry.

A source named `SOURCE` lists its record names in that column, under the one key; Ensembl VEP's tab and VCF writers list the column twice, vep-rs once ([the class](divergences/duplicate-output-column.md)).

### Header lines

The headers describe each key.

| Key | Header text |
| --- | --- |
| `<short_name>` | the file, its directory replaced by `[PATH]` |
| `<short_name>_<field>` | `<Description>. <field> field from <file>` |
| `<short_name>_PC` | `Percentage of the annotation record covered by the input variant from <file>` |
| `<short_name>_<stat>` | `<stat> data from <file>` |
| a key that one source declares as a field and as a statistic | both texts: on one line, the statistic's first, joined by `; `, in the tab and VCF headers, and on a line for each declaration in the default format ([the key declared twice](divergences/custom-key-named-as-a-plain-value-field.md#36-the-description-of-a-key-one-source-declares-twice)) |
| a key that two sources of one short name declare, one as a field and one as a statistic | both texts, each with the files of the sources that declare it: on one line, the statistic's first, in the tab and VCF headers, and on separate lines in the default format ([the two declarations](divergences/pc-description-text.md#22-a-key-the-two-sources-declare-apart-as-a-record-field-and-as-a-statistic)) |

A field's description is taken from a VCF source's own `##INFO` line, else it is the field's name with its first letter upper-cased. Under `reciprocal=1` the `PC` text reads `Percentage of the input variant covered by the annotation record`. Ensembl VEP's own two texts name the opposite share ([The description of a custom source's `PC` field](divergences/pc-description-text.md)).

Two sources under one short name that declare a key alike share one line, which lists both files. Ensembl VEP appends the second file to the first custom header line unless the key holds a lowercase `o` ([the repeated short name's header lines](divergences/pc-description-text.md#2-the-header-line-that-receives-a-repeated-custom-short-names-file)). The VCF writes each key as an `##INFO` line after the plugin lines.

A VCF source with `fields` carries `<short_name>_FILTER` on every record. The default format and the JSON print it whether or not `fields` names `FILTER`. The tab, VCF and Parquet columns come from the header, so they carry it only when `fields` names it.

### The rows a record prints on

A record matched to the whole input record prints on every row of that input record, listed before the row's allele records. So does the `...` marker of a `num_records` limit, under an `exact` source too: a row lists the marker before its own allele's records, and a row of an allele with no kept record lists the marker alone. An `exact` VCF record prints on the rows of the allele it matched. It never prints on the reference allele's row of `--process_ref_homs` nor on a `--summary` or `--most_severe` row.

A source that matched no record of an input record prints nothing on its text rows, statistics included, unless another source of its short name lists a record there. The JSON lists the statistics object it asked for beside the records of any other source, whatever its short name.

### The JSON format

In the JSON format `custom_annotations` is an object on the input record, keyed by short name. Each value lists the statistics object first, when statistics were asked for, then one object per record of each source under the short name. A record's object holds its `name`, the `allele` an `exact` match named, its `fields` and its `score`. A value that reads as a number is a number and a `FILTER` list an array. An input record with no matched record has no `custom_annotations` key, and the keys appear in no consequence object.

## Where `--custom` differs from Ensembl VEP

Two defects of Ensembl VEP 116.3's custom-source reader are classes documented below. Six more are on other pages, beside the dropped `###CHR###` pattern and the `PC` description stated above:

- three bigWig behaviours ([bigWig sources and conservation shortcuts](output/custom-bigwig.md))
- the header line of a repeated short name ([the `PC` description page](divergences/pc-description-text.md#2-the-header-line-that-receives-a-repeated-custom-short-names-file))
- a record limit that is not a whole number ([Defect N](divergences/defect-n.md#26-a---custom-record-limit-that-is-not-a-whole-number))
- [a VCF source row whose `END` lies below its position](divergences/tabix-row-with-an-end-below-its-start.md)

Each is a documented class: the section names the rows a comparison sets aside. The class applies to no committed row, as [Intended divergences from Ensembl VEP](intended-divergences.md) lists it. That page's [How to read the counts](intended-divergences.md#how-to-read-the-counts) defines the terms. Every row below is in the default output format, against the `GRCh38-custom` golden corpus's cache and records (`tests/golden/116/GRCh38-custom/`) unless the section names another.

Each class carries its tier under that page's rule:

- tier 1 is a shape the chromosome 21 cuts carry with no flag
- tier 2 one that needs a flag or plugin on ordinary data
- tier 3 one ordinary data does not carry

The counts are over the cuts' 15,646,481 records, the 152 golden corpora and their 511 flag sets.

### 1. An exact custom record repeated once per row beside a whole-variant record of the same name

**The defect in one sentence.** When a short name's records at a variant include a record attached to the whole variant and an `exact` record attached to one allele, Ensembl VEP 116.3 lists the allele's record once for each row of its own allele, on every row of the variant.

**Why this is a defect.** The whole-variant records are stored in the variant's output hash as a list, and each row's hash is a shallow copy of that hash, so every row shares the one list. Each row then pushes its allele's records onto the shared list, and the rows are written only after all of them are built. So a variant of one alternate with nine rows prints the allele's record nine times on each row. With several alternates each allele's record is pushed once for each row of that allele, and every row of the variant prints them all, the rows of the other alleles included.

The record matched the variant once; the number of copies counts the rows of its allele, not its matches. The `FILTER` and field lists a whole-variant record carries repeat with it. Beside a `...` marker, which carries no field, the name list alone repeats and the field lists agree with vep-rs's.

**What vep-rs does instead.** vep-rs keeps the records of each source once and decides per row which of them print, so an allele's record prints once on each row of its allele, beside the whole-variant records.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1004-1014, `VariationFeature_to_output_hash`: the records without an `allele` are pushed onto the variant hash's `<short_name>` list by `_add_custom_annotations_to_hash`, lines 1297-1312.
- The same file, line 391, `get_all_VariationFeatureOverlapAllele_output_hashes`: `my %copy = %$hash;` copies the hash one level deep, so every row's hash holds the same list reference. Lines 1416-1435, `VariationFeatureOverlapAllele_to_output_hash`, push the row's allele records onto it, and lines 266-271 write the rows after every hash is built.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/VCF.pm`, lines 380-387 and 247-263, with `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 571-574: an `exact` source's record matches by span and has no `allele` where its ALT is symbolic without `*`, a breakend, or not bases alone beside an `SVTYPE`, or the input record is structural.
- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 436-440, `annotate_VariationFeature`: the `...` marker of a `num_records` limit is pushed without an `allele` as well, so it stands as a whole-variant record beside the allele records the limit kept. It carries no field, and `OutputFactory.pm` lines 1300-1304 put a field's list on the variant's hash only for a record that has the field. So beside a marker each row builds its own field lists, and the name list alone repeats.

**What correcting it would change.** Each row would list the allele's record once, as vep-rs does and as Ensembl VEP does when the variant's hash holds no list under the short name.

**The class it produces.** `custom_exact_source_allele_record_accumulation`: on a default, tab or VCF row of a sequence variant, a `<short_name>`, `<short_name>_FILTER` or `<short_name>_<field>` list whose reference value holds the whole-variant records and then every allele-matched record of the variant, each repeated once for each row of its own allele. The candidate's row holds the whole-variant records and its own allele's records once, so on a row of another allele the reference lists copies the candidate's row lacks. Beside a `...` marker the `<short_name>` list alone differs. In the shape of 1.2 the flag's list comes first. The JSON `custom_annotations` object belongs to the variant and is unaffected.

The committed `exact` sources hold no structural record, no committed `exact` run lists a `...` marker and no run names a short name twice, so the class applies to no committed row or comparison case.

**How often the shape occurs.** Tier 2: the shape is a record without an allele beside a sequence record of the same `exact` VCF source at the variant. That record is a structural one, matched by span, or the `...` marker of a `num_records` limit the variant's matches exceed. An ordinary record of several alternates reaches the marker route under such a limit, where the source holds a record for each allele.

The index's cell counts that route. Its rate counts the records of two or more sequence alternates: 6,347 of the population cuts' records (405.6 per million), and 80 corpora hold one. Its flag sets are those that set `num_records` on an `exact` VCF source: none of the 511, the one committed limit standing on a `type=overlap` source. A record of one alternate reaches the marker only where the source holds more records for its allele than the limit keeps. None of the 96 records of the committed custom VCF sources is structural, and 6 flag sets run a VCF source.

#### 1.1 A deletion record and a substitution record of one `exact` source at `21:10542415`

**In short.** An `exact` VCF source `mix` holds `21 10542414 svrec N <DEL>` with `SVTYPE=DEL;END=10542415;TAG=sv` and `21 10542415 seqrec A G` with `TAG=seq`, both with `PASS` in their `FILTER` column; the input `21 10542415 A G`, an allele of a `GRCh38-custom` corpus record, has nine transcript rows. The structural record matches by span and belongs to the whole variant; the substitution matches the allele `G`. Ensembl VEP 116.3 prints on all nine rows `mix=svrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec`, with ten `PASS` entries in `mix_FILTER` and `sv,seq,seq,...` in `mix_TAG`; vep-rs prints `mix=svrec,seqrec;mix_FILTER=PASS,PASS;mix_TAG=sv,seq` on each, one entry per record matched.

#### 1.2 An `exact` record under a short name that spells a list field

**In short.** A second shape of the class needs no whole-variant record: the short name spells a list Ensembl VEP keeps in the variant's hash, `ZYG` under `--individual_zyg` (`OutputFactory.pm` line 987), `NEAREST` under `--nearest` (line 1017) or `SV` under `--check_svs` (line 963). Each row pushes its allele's record of an `exact` source onto that list. vep-rs writes the list, then the record once.

| Short name, flag and input | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| Example 1. `ZYG`, `--individual_zyg all`: `21 10542415 A G` with one heterozygous sample `S1`, nine output rows, and a source that holds `seqrec` alone | `ZYG=S1:HET,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec` | `ZYG=S1:HET,seqrec` |
| Example 2. `ZYG`, `--individual_zyg all`: line 111 of `tests/golden/116/GRCh38-individual-zyg/variants.vcf`, one alternate allele, eight samples and two output rows, and a made-up record `made_up_1` on that allele with `AF=0.25` and `PASS`, read with `fields=AF` | the eight zygosities, then `made_up_1,made_up_1` | the eight zygosities, then `made_up_1` |
| Example 3. `NEAREST`, `--nearest symbol`: line 42 of `tests/golden/116/GRCh38-custom/variants.vcf`, one alternate allele and 16 output rows, and the record `multi_7` of the corpus's `custom/sites.vcf.gz` | `NEAREST=SON`, then `multi_7` 16 times | `NEAREST=SON,multi_7` |
| Example 2 with its source read as `type=overlap`, whose record belongs to the whole variant and is pushed once | the eight zygosities, then `made_up_1` | the same |
| Example 2 or 3 with its source under the short name `sites` | `sites=made_up_1`, or `sites=multi_7` | the same |
| Example 2 or 3 with its flag off | `ZYG=made_up_1`, or `NEAREST=multi_7` | the same |
| The `ZYG_AF` and `ZYG_FILTER` lists of example 2, which each row starts for itself | `0.25` and `PASS` | the same |
| `SV`, `--check_svs` | not run: Ensembl VEP refuses the flag without a database, and line 963 sets the list as line 987 does | the names `--check_svs` lists, then the record once, as a unit test holds it |

Each example's cell stands on every row of the variant, the flag's list first; a one-row variant shows no difference. Ensembl VEP's tab and VCF outputs carry the key's column twice ([`duplicate_output_column_collapsed`](divergences/duplicate-output-column.md)). No committed flag set names a source `ZYG`, `NEAREST` or `SV`. The engines order the zygosities differently ([`individual_zyg_list_order`](divergences/individual-zyg-order.md)).

### 2. The summary statistics of two custom sources sharing a short name

**The defect in one sentence.** When two `--custom` sources share a short name and ask for `summary_stats`, Ensembl VEP 116.3 keeps one running statistics hash per short name: `min`, `max` and `sum` pool both sources' records, `count` is the last source's, and `mean` divides the pooled sum by that count.

**Why this is a defect.** Each source walks a variant's records on its own and closes its statistics with its own record count. The hash it closes is keyed by short name alone, so the second source finds the first's `min`, `max` and `sum`, extends them, overwrites `count` with its own and divides the pooled sum by it. For scores 4 and 9, one record each, the keys print `min 4`, `max 9`, `sum 13`, `count 1` and `mean 13`: a mean no record carries, neither source's and not the pooled 6.5. When `sum` is not asked for it is deleted after each source and the second source's `mean` is its own; `min` and `max` pool regardless.

**What vep-rs does instead.** vep-rs closes each source's statistics over its own records and prints, under the shared keys, the figures of the last source that asked for statistics and matched a record at the variant. `<short_name>` lists every source's records, a `COSV` record named like an earlier source's once.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 423-433, `annotate_VariationFeature`: the hash is `$vf->{_custom_annotations_stats}->{short_name}`, read and extended by every source of that name, `min` and `max` at lines 426-431 and `sum` at line 432.
- The same file, lines 344-352, `annotate_InputBuffer`: `count` is set to this source's `$record_count` at line 346, `mean` to the hash's `sum` over that count at lines 347-348, and `sum` is deleted at line 350 unless asked for.

**What correcting it would change.** Each source's keys would describe one population of records, its own as vep-rs prints, or both sources' records pooled with their joint count.

**The class it produces.** `custom_shared_short_name_statistics`: on a run naming one short name on two or more sources with `summary_stats`, the `<short_name>_min`, `_max`, `_mean`, `_count` and `_sum` fields of a row where more than one source matched a record. The reference pools across the sources and the candidate prints the own figures of the last matching source that asked for them; the class is applied from a comparison case's documentation. No committed flag set or comparison case names a short name twice, so it applies to no committed row.

**How often the shape occurs.** Tier 3: the shape is two sources sharing a short name, each asking summary statistics, which none of the 511 flag sets does; 19 of them run a `--custom` source.

#### 2.1 Two BED sources named `regions`, scores 4 and 9, at `21:5380280`

**In short.** Two BED sources `regA.bed.gz` (score 4) and `regB.bed.gz` (score 9), each with one record at `21:5380279-5380280` and each given as `short_name=regions,format=bed,type=overlap,summary_stats=min%max%count%sum%mean`, meet the `GRCh38-custom` corpus record `21 5380280 rs1444031311 A C,G`. Ensembl VEP 116.3 prints `regions=regA,regB;regions_count=1;regions_max=9;regions_mean=13;regions_min=4;regions_sum=13`; vep-rs prints `regions=regA,regB;regions_count=1;regions_max=9;regions_mean=9;regions_min=9;regions_sum=9`, the second source's own figures. Neither 4 nor 13 describes the second source's record, and 13 is the mean of nothing: the two records' sum over one record's count. Without `sum` among the statistics Ensembl VEP prints `regions_mean=9;regions_min=4`.

## Where next

- [A GTF or GFF3 as the transcript model](cache-from-gtf-or-gff3.md): how to build a cache from a GTF or GFF3 file, or read the file beside a cache.
- [Chromosome names](cache-setup.md#chromosome-names): the synonyms and the `chr` rule that match an input chromosome to a source's.
- [Intended divergences from Ensembl VEP](intended-divergences.md): every documented class with its tier, and how to read the counts.
