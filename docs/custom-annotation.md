# Custom annotation

`--custom` names a tabix-indexed VCF, BED or GFF file, or a bigWig, whose records annotate the input variants under a short name of the source's choosing. The string is read as Ensembl VEP reads it, in either of its two forms:

```text
--custom file=<file>,short_name=<name>,format=<vcf|bed|bigwig|gff>,
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

Ensembl VEP 116.2 dies with that text only when the pattern is the run's first `--custom` source. After an earlier source it keeps the sources it has and drops the pattern without a key, a header line or a warning. The test guarding the message counts every source the run has built so far rather than the pattern's own (`AnnotationSourceAdaptor.pm` lines 330 and 451-452). No committed run names such a pattern (tier 3 of [Intended divergences from Ensembl VEP](intended-divergences.md); 14 flag sets run a `--custom` source).

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

Two shapes of source are not read, each refused with an error that names the route:

- A remote file (`http://`, `https://`, `ftp://` or any other URL scheme). Ensembl VEP reads it over the network; vep-rs reads local files only, so copy the file to a local path.
- A GFF or GTF transcript model (`format=gtf`, or `format=gff` in its default `transcript` mode), which Ensembl VEP reads in place of its cache. vep-rs reads transcripts from its cache alone. Build it from the file and the genome with `vep-cache-builder --gtf <file>` (or `--gff3 <file>`) `--genome-fasta <genome.fa>`, then annotate against it. [A GTF or GFF3 as the transcript model](cache-from-gtf-or-gff3.md) describes the steps. A GFF in `gencode_promoter` mode is an overlay source, [GFF promoter windows](#gff-promoter-windows).

## How a source is read

A VCF, BED or GFF source is read record by record against each input record, as Ensembl VEP reads it. The records a batch can reach are read once. A window one base either side of every input record is merged into as few reads as the file's layout allows. Each input record then takes the records of its own window in file order.

A `--custom` run reads its cache as every run does; Ensembl VEP accepts a `--custom` source as a run's only annotation source. Files are read locally through their tabix index (`.tbi` or `.csi`); a BAM is not read, and a remote file is refused as above.

## Chromosomes

A source lists the chromosomes its index names (a bigWig, its chromosome list). An input chromosome is matched to them as to the cache's, through the synonyms and the `chr` rule of [Chromosome names](cache-setup.md#chromosome-names).

An input record on a chromosome the cache lacks is admitted when a source, or a `--regulatory_gff` file, lists it. It is written as `intergenic_variant` with the sources' values: Ensembl VEP admits a chromosome any annotation source holds (`BaseRunner.pm` lines 444-460). A bare `M` is renamed `MT` when a source or the cache lists `MT` and nothing lists `M` or `chrM` (`Parser.pm` lines 794-796).

## Matching

A BED record spans `chromStart + 1` to `chromEnd`. A VCF record naming a structural variant (an `SVTYPE`, a symbolic or a breakend allele) spans:

- `POS + 1` to `END`
- or to `POS + |SVLEN|` when `SVLEN` is given
- or `POS + 1` alone with neither

A record whose ALT names `INS` ends at `POS + length(REF)`. Against those spans, and against a structural input record whatever the record, the four types read as Ensembl VEP's.

| Type | Asks for |
| --- | --- |
| `overlap` | a shared base |
| `within` | a record inside the input record |
| `surrounding` | a record holding the input record |
| `exact` | the same start and the same end |

`distance` asks each end of the record to lie within that many bases of the input record's end. `overlap_cutoff` asks the shared bases to be at least that percentage of the record, or of the input record under `reciprocal=1`. `PC` among the `fields` reports that percentage to three decimals.

A sequence VCF record, its ALT bases alone, spans `POS` to `POS + length(REF) - 1`. It matches an `overlap` source on a shared base. Under `exact`, each alternate allele of the record is compared with each alternate allele of the input record. Both pairs are trimmed to their minimal form and position, and the record attaches to the alleles whose forms agree. A sequence record never matches a `within` or `surrounding` source, as in Ensembl VEP, whose VCF reader has branches for `overlap` and `exact` alone.

`same_type` admits a record only when its structural class is the input record's, so it matches structural records alone.

## Values

A record is named by its first `ID` (VCF) or its name column (BED), or by `chromosome:start-end` under `coords=1` or when it has no identifier.

The `fields` of a VCF source are its INFO keys, a key without a value reading `1`. A comma-separated value with one entry per alternate allele gives each allele its own entry under `exact`. Under any other type it rejoins the entries in allele order. A value with one entry more drops its first. `--custom_multi_allelic`, and a file whose header says `##source=ClinVar`, keep every value whole on every allele.

A VCF source with `fields` adds the record's `FILTER` entries, split at `;`, unless `--custom_suppress_filter` is given and `fields` does not name `FILTER`. A record's score is its `QUAL` (VCF) or its fifth column (BED).

## Records and statistics per input record

An input record keeps at most `num_records` matching records, then one named `...`. `num_records=all` keeps every one. A `num_records` with a fraction or a minus sign (`2.5`, `-1`) is read as its integer part, at least 0. The source writes one warning naming the value and the reading. Ensembl VEP compares the count with the text as a number ([Defect N, class 26](divergences/defect-n.md#26-a---custom-record-limit-that-is-not-a-whole-number)).

`summary_stats` computes `min`, `max`, `mean`, `count` and `sum` over the matching records' scores:

- `count` counts the records
- `min` and `max` keep a score of 0
- `mean` is the sum over the count
- a record without a score counts but adds no value

A `COSV` record can be named like one already kept under the same short name, by this source or an earlier one. Such a record lends the kept one the fields it lacks and is not kept twice. It counts, adds no score and takes no `PC`.

The records and statistics belong to the whole input record, so every allele of a multi-allelic input record shares them. An `exact` VCF record belongs to the allele it matched, read as the input gives it. A 3' shift of the allele a transcript row prints leaves the record on its own allele ([section 3 below](#3-an-exact-custom-record-on-a-row-whose-allele---shift_3prime-1-rotated)).

## GFF promoter windows

`format=gff,gff_type=gencode_promoter` reads a tabix-indexed GFF of GENCODE promoter windows as an overlay. It is the one GFF mode annotating beside the cache rather than standing in for it. A record spans its fourth to its fifth column and matches by the rules of [Matching](#matching); `exact` is refused, as Ensembl VEP refuses it for any GFF source, and so is an index `tabix` built past a start that is no number, [`gxf_malformed_coordinate_line_refused`](divergences/gxf-malformed-coordinate-line-refused.md).

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

`--custom` also adds the `SOURCE` column, which names the cache a merged cache's transcript came from, `Ensembl` or `RefSeq`, and is empty on a transcript of a single-source cache, a cache built from a GTF or GFF included. Ensembl VEP's `--gtf` and `--gff` sources stamp their transcripts with the file's name, which a cache built from that file does not carry.

A source named `SOURCE` lists its record names in that column, under the one key; Ensembl VEP's tab and VCF writers list the column twice, vep-rs once ([the class](divergences/duplicate-output-column.md)).

### Header lines

The headers describe each key.

| Key | Header text |
| --- | --- |
| `<short_name>` | the file, its directory replaced by `[PATH]` |
| `<short_name>_<field>` | `<Description>. <field> field from <file>` |
| `<short_name>_PC` | `Percentage of the annotation record covered by the input variant from <file>` |
| `<short_name>_<stat>` | `<stat> data from <file>` |

A field's description is taken from a VCF source's own `##INFO` line, else it is the field's name with its first letter upper-cased. Under `reciprocal=1` the `PC` text reads `Percentage of the input variant covered by the annotation record`. Ensembl VEP's own two texts name the opposite share ([The description of a custom source's `PC` field](divergences/pc-description-text.md)).

Two sources under one short name share one line, which lists both files. Ensembl VEP appends the second file to the first custom header line unless the key holds a lowercase `o` ([the repeated short name's header lines](divergences/pc-description-text.md#2-the-header-line-that-receives-a-repeated-custom-short-names-file)). The VCF writes each key as an `##INFO` line after the plugin lines.

A VCF source with `fields` carries `<short_name>_FILTER` on every record. The default format and the JSON print it whether or not `fields` names `FILTER`. The tab, VCF and Parquet columns come from the header, so they carry it only when `fields` names it.

### The rows a record prints on

A record matched to the whole input record prints on every row of that input record, listed before the row's allele records. An `exact` VCF record prints on the rows of the allele it matched. It never prints on the reference allele's row of `--process_ref_homs` nor on a `--summary` or `--most_severe` row.

A source that matched no record of an input record prints nothing on its text rows, statistics included. The JSON lists the statistics object it asked for beside another source's records.

### The JSON format

In the JSON format `custom_annotations` is an object on the input record, keyed by short name. Each value lists the statistics object first, when statistics were asked for, then one object per record of each source under the short name. A record's object holds its `name`, the `allele` an `exact` match named, its `fields` and its `score`. A value that reads as a number is a number and a `FILTER` list an array. An input record with no matched record has no `custom_annotations` key, and the keys appear in no consequence object.

## Where `--custom` differs from Ensembl VEP

Three defects of Ensembl VEP 116.2's custom-source reader are classes documented below. Six more are on other pages, beside the dropped `###CHR###` pattern and the `PC` description stated above:

- three bigWig behaviours ([bigWig sources and conservation shortcuts](output/custom-bigwig.md))
- the header line of a repeated short name ([the `PC` description page](divergences/pc-description-text.md#2-the-header-line-that-receives-a-repeated-custom-short-names-file))
- a record limit that is not a whole number ([Defect N](divergences/defect-n.md#26-a---custom-record-limit-that-is-not-a-whole-number))
- [a VCF source row whose `END` lies below its position](divergences/tabix-row-with-an-end-below-its-start.md)

Each is a documented class: the section names the rows a comparison sets aside. The class applies to no committed row, as [Intended divergences from Ensembl VEP](intended-divergences.md) lists it. That page's [How to read the counts](intended-divergences.md#how-to-read-the-counts) defines the terms. Every row below is in the default output format, against the `GRCh38-custom` golden corpus's cache and records (`tests/golden/116/GRCh38-custom/`) unless the section names another.

Each class carries its tier under that page's rule:

- tier 1 is a shape the chromosome 21 population cuts carry with no flag
- tier 2 one that needs a flag or plugin on ordinary data
- tier 3 one ordinary data does not carry

The counts are over the cuts' 15,646,481 records, the 111 golden corpora and their 392 flag sets.

### 1. An exact custom record repeated once per row beside a whole-variant record of the same name

**The defect in one sentence.** When a short name's records at a variant include both a record attached to the whole variant and an `exact` record attached to one allele, Ensembl VEP 116.2 lists the allele's record once per row of the variant, on every row.

**Why this is a defect.** The whole-variant records are stored in the variant's output hash as a list, and each row's hash is a shallow copy of that hash, so every row shares the one list. Each row then pushes its allele's records onto the shared list, and the rows are written only after all of them are built, so a variant with nine rows prints the allele's record nine times on each row. The record matched the variant once; the number of copies states the number of transcripts, not of matches, and the `FILTER` and field lists repeat with it.

**What vep-rs does instead.** vep-rs keeps the records of each source once and decides per row which of them print, so an allele's record prints once on each row of its allele, beside the whole-variant records.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1004-1014, `VariationFeature_to_output_hash`: the records without an `allele` are pushed onto the variant hash's `<short_name>` list by `_add_custom_annotations_to_hash`, lines 1297-1312.
- The same file, line 391, `get_all_VariationFeatureOverlapAllele_output_hashes`: `my %copy = %$hash;` copies the hash one level deep, so every row's hash holds the same list reference. Lines 1416-1426, `VariationFeatureOverlapAllele_to_output_hash`, push the row's allele records onto it, and lines 266-271 write the rows after every hash is built.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/VCF.pm`, lines 380-387 and 247-263, with `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 571-574: a structural record of an `exact` source (a symbolic allele, a breakend or an `SVTYPE`) matches by span and is built without an `allele`. One `exact` source can therefore hold both kinds of record.

**What correcting it would change.** Each row would list the allele's record once, as vep-rs does and as Ensembl VEP does when the short name has no whole-variant record at the variant.

**The class it produces.** `custom_exact_source_allele_record_accumulation`: on a default, tab or VCF row of a sequence variant, a `<short_name>`, `<short_name>_FILTER` or `<short_name>_<field>` list whose reference value is the candidate's with each allele-matched record repeated, the repeat count being the variant's row count. The JSON `custom_annotations` object belongs to the variant and is unaffected. The committed `exact` sources hold no structural record and no run names a short name twice, so the class applies to no committed row or comparison case.

**How often the shape occurs.** Tier 3: the shape is a structural record of an `exact` VCF source, matched by span without an allele, beside a sequence record of the same source at the variant. None of the 96 records of the committed custom VCF sources is structural, and 6 flag sets run a VCF source.

#### 1.1 A deletion record and a substitution record of one `exact` source at `21:10542415`

**In short.** An `exact` VCF source `mix` holds `21 10542414 svrec N <DEL>` with `SVTYPE=DEL;END=10542415;TAG=sv` and `21 10542415 seqrec A G` with `TAG=seq`, both with `PASS` in their `FILTER` column; the input `21 10542415 A G`, an allele of a `GRCh38-custom` corpus record, has nine transcript rows. The structural record matches by span and belongs to the whole variant; the substitution matches the allele `G`. Ensembl VEP 116.2 prints on all nine rows `mix=svrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec,seqrec`, with ten `PASS` entries in `mix_FILTER` and `sv,seq,seq,...` in `mix_TAG`; vep-rs prints `mix=svrec,seqrec;mix_FILTER=PASS,PASS;mix_TAG=sv,seq` on each, one entry per record matched.

### 2. The summary statistics of two custom sources sharing a short name

**The defect in one sentence.** When two `--custom` sources share a short name and ask for `summary_stats`, Ensembl VEP 116.2 keeps one running statistics hash per short name: `min`, `max` and `sum` pool both sources' records, `count` is the last source's, and `mean` divides the pooled sum by that count.

**Why this is a defect.** Each source walks a variant's records on its own and closes its statistics with its own record count. The hash it closes is keyed by short name alone, so the second source finds the first's `min`, `max` and `sum`, extends them, overwrites `count` with its own and divides the pooled sum by it. For scores 4 and 9, one record each, the keys print `min 4`, `max 9`, `sum 13`, `count 1` and `mean 13`: a mean no record carries, neither source's and not the pooled 6.5. When `sum` is not asked for it is deleted after each source and the second source's `mean` is its own; `min` and `max` pool regardless.

**What vep-rs does instead.** vep-rs closes each source's statistics over its own records and prints, under the shared keys, the figures of the last source that matched a record at the variant, while `<short_name>` lists every source's records, a `COSV` record named like an earlier source's once.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 423-433, `annotate_VariationFeature`: the hash is `$vf->{_custom_annotations_stats}->{short_name}`, read and extended by every source of that name, `min` and `max` at lines 426-431 and `sum` at line 432.
- The same file, lines 344-352, `annotate_InputBuffer`: `count` is set to this source's `$record_count` at line 346, `mean` to the hash's `sum` over that count at lines 347-348, and `sum` is deleted at line 350 unless asked for.

**What correcting it would change.** Each source's keys would describe one population of records, its own as vep-rs prints, or both sources' records pooled with their joint count.

**The class it produces.** `custom_shared_short_name_statistics`: on a run naming one short name on two or more sources with `summary_stats`, the `<short_name>_min`, `_max`, `_mean`, `_count` and `_sum` fields of a row where more than one source matched a record. The reference pools across the sources and the candidate prints the last matching source's own figures; the class is applied from a comparison case's documentation. No committed flag set or comparison case names a short name twice, so it applies to no committed row.

**How often the shape occurs.** Tier 3: the shape is two sources sharing a short name, each asking summary statistics, which none of the 392 flag sets does; 14 of them run a `--custom` source.

#### 2.1 Two BED sources named `regions`, scores 4 and 9, at `21:5380280`

**In short.** Two BED sources `regA.bed.gz` (score 4) and `regB.bed.gz` (score 9), each with one record at `21:5380279-5380280` and each given as `short_name=regions,format=bed,type=overlap,summary_stats=min%max%count%sum%mean`, meet the `GRCh38-custom` corpus record `21 5380280 rs1444031311 A C,G`. Ensembl VEP 116.2 prints `regions=regA,regB;regions_count=1;regions_max=9;regions_mean=13;regions_min=4;regions_sum=13`; vep-rs prints `regions=regA,regB;regions_count=1;regions_max=9;regions_mean=9;regions_min=9;regions_sum=9`, the second source's own figures. Neither 4 nor 13 describes the second source's record, and 13 is the mean of nothing: the two records' sum over one record's count. Without `sum` among the statistics Ensembl VEP prints `regions_mean=9;regions_min=4`.

### 3. An exact custom record on a row whose allele `--shift_3prime 1` rotated

**The defect in one sentence.** Under `--shift_3prime 1`, Ensembl VEP 116.2 attaches a row's `exact` custom records by the rotated allele the row prints, so a shifted coding row loses its record, or takes another alternate allele's record whose spelling equals the rotation.

**Why this is a defect.** An `exact` record is matched against the input variant before any shifting and is keyed on the allele as given. The 3' shift rotates an inserted or deleted sequence, `GGGCC` becoming `CGGGC` after a shift of 14, on the transcripts it shifts alone, so one variant carries its record on its non-coding and upstream rows and none on its coding rows. Where one alternate allele's rotation spells another alternate allele, the first allele's rows print the second's record, a `Benign` beside a variant whose own record says `Pathogenic`. Ensembl VEP 116.3 corrects this in `OutputFactory.pm` as "Keep custom annotations on 3' shifted alleles", keying the records on the pre-shift allele the row's shift hash keeps.

**What vep-rs does instead.** vep-rs attaches an `exact` record to the rows of the allele it matched, read as the input gives it, whatever spelling the row's `Allele` column prints after the shift, as Ensembl VEP 116.3 does.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1416-1426, `VariationFeatureOverlapAllele_to_output_hash`: the row's records are those whose `allele` equals `$hash->{Allele}`; line 1347 sets that allele to the row's shifted sequence, and lines 1349-1352 to the genomic shift's under `--shift_genomic 1`.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/VCF.pm`, lines 405-421 and 230-236: an `exact` record is matched against the variant's own alleles and keyed on the matched input allele.
- Ensembl VEP 116.3, the same sub of `OutputFactory.pm`: the records are selected by the shift hash's `alt_orig_allele_string` when the row has one.

**What correcting it would change.** The class vanishes, as it does against Ensembl VEP 116.3: every row of an allele carries the allele's own record.

**The class it produces.** `custom_shifted_allele_records`: under `--shift_3prime 1` with an `exact` VCF source, a `<short_name>`, `<short_name>_<field>` or `<short_name>_FILTER` field the candidate alone fills, or both fill with different records, on a row whose `Allele` is not the alternate allele the input gives it; applied from a comparison case's documentation. No committed flag set or comparison case runs `--shift_3prime 1` together with a `--custom` source, the three `--shift_3prime` flag sets and five cases carrying no source and the twelve cases with a custom source no shift, so it applies to no committed row.

**How often the shape occurs.** Tier 2: 1,203,766 of the population cuts' records (76,935 per million) are indels whose inserted or deleted bases number at least two and are not all one letter, so that a 3' shift can rotate the spelling; the genome decides the shift. 89 corpora carry the shape, and no committed flag set runs `--shift_3prime 1` with an `exact` VCF source.

#### 3.1 A three-base insertion in a coding exon, `21 13620283 G GTCA`

**In short.** The `GRCh38-identifiers` corpus record `cds_tri_ins_fwd_1`, `21 13620283 G GTCA`, meets an `exact` source `ClinVar` whose record `rec_tca`, `G GTCA` with `PASS` in `FILTER`, carries `CLNSIG=Pathogenic`. Under `--shift_3prime 1` the two coding rows, `ENST00000299443` and `ENST00000620442`, print the rotated allele `CAT`: Ensembl VEP 116.2 prints `IMPACT=MODERATE;STRAND=1` and no `ClinVar` key on them, while the unrotated upstream row, allele `TCA`, carries the record. vep-rs prints `ClinVar=rec_tca;ClinVar_CLNSIG=Pathogenic;ClinVar_FILTER=PASS` on all three. With `--shift_genomic 1` added, the input `21 10605702 T TAT,TTA` with records `T TAT` `Pathogenic` and `T TTA` `Benign` shows the second shape. The first allele's eight rotated rows, allele `TA`, print the second allele's `ClinVar=rec_ta;ClinVar_CLNSIG=Benign` in Ensembl VEP 116.2 and their own `rec_at`, `Pathogenic`, in vep-rs.

## Where next

- [A GTF or GFF3 as the transcript model](cache-from-gtf-or-gff3.md): how to build a cache from a GTF or GFF3 file and annotate against it.
- [Chromosome names](cache-setup.md#chromosome-names): the synonyms and the `chr` rule that match an input chromosome to a source's.
- [Intended divergences from Ensembl VEP](intended-divergences.md): every documented class with its tier, and how to read the counts.
