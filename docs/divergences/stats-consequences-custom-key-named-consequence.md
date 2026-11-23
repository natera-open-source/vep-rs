# The statistics report's consequence counts under a `--custom` source keyed `Consequence`

**The defect in one sentence.** Under a `--custom` VCF source of `type=exact` whose key is `Consequence`, the `[Consequences (all)]` section of Ensembl VEP's statistics report counts the name of each matching source record as a consequence type, once for each output row of the allele the record names. The name is the ID of the record's line. vep-rs counts the consequence terms of the run's rows and nothing else.

**Why this is a defect.** The section counts the consequence types of the run's rows. A source record's name identifies a line of the user's annotation file and says nothing of a variant's effect. Ensembl VEP's own report shows the count is not meant. It sorts the section's lines by consequence rank, has no rank for the name and warns of an undefined value. It reads the same counter by consequence type for `[Coding consequences]`. Its `[Consequences (most severe)]` section never counts the name, and the same run written as JSON does not count it at all. Every distinct name among the source's matches adds a line of its own.

**What vep-rs does instead.** vep-rs counts each row's consequence terms, whatever the source's key and whatever the output format. Outside `[VEP run statistics]`, its report equals the report of the same run with the source under any other key.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1341: the row's consequence terms are assigned to the row's `Consequence` list.
- The same file, lines 1426-1435, which call the routine of lines 1297-1312: line 1301 pushes the name of each source record on the row's allele onto the list under the source's key. Under the key `Consequence` that is the list of line 1341.
- `Bio/EnsEMBL/VEP/Stats.pm`, line 282, `log_VariationFeatureOverlapAllele`: every item of the list is counted as a consequence type.

**What correcting it would change.** A correction in Ensembl VEP would remove the difference. Line 1337 of `OutputFactory.pm` holds the row's consequences as objects before any name is pushed. A count of those objects at `Stats.pm` line 282 gives the section vep-rs writes, in every output format.

**The class it produces.** `stats_consequences_custom_key_named_consequence` is the name comparisons count this difference under. Its shape is a line of `[Consequences (all)]` whose name is no consequence type, in the report of a default, tab or VCF run under such a source. The concordance figures and the output-format comparisons never read the statistics report, so nothing is set aside in them. No committed run of the repository's test corpora keys a source `Consequence`.

The `Consequence` cell of the annotation output agrees. In those three formats both engines print the name after the row's terms in the cell of the named allele's rows, and in JSON neither lists it among an entry's `consequence_terms`. Ensembl VEP's tab and VCF outputs also list the key `Consequence` a second time, as one more column and one more `CSQ` field: the documented difference [`duplicate_output_column_collapsed`](duplicate-output-column.md).

Set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. An input record with no transcript in reach

**In short.** On the input record `21 5380280 rs1444031311 A C,G`, whose first alternate the source record `exact_2` names, Ensembl VEP 116.3 writes the line `exact_2 1` above `intergenic_variant 2` in `[Consequences (all)]` of a default, tab or VCF run. In a JSON run it writes the term's line alone. vep-rs writes the term's line alone in all four formats. Every other line of the two reports agrees outside `[VEP run statistics]`, the section of each run's own command line, times and paths.

### 1.1 What a reader needs to know first

`--stats_text` writes a plain-text statistics report, `<output>_summary.txt`, beside the output; the HTML report draws the same sections as charts from the same counts. `[Consequences (most severe)]` counts each input variant once, under its most severe consequence. `[Consequences (all)]` counts every consequence term of every output row, so a row with two terms adds one to each. The report's first section, `[VEP run statistics]`, holds the run's command line in its `Command line options` row.

`--custom` annotates each variant from a file the user names, the source. The source's `short_name` is its key: the names of the matching source records print under it, and each field the source asks for prints under the key, `_` and the field's name. Ensembl VEP's documentation takes any string without commas as a short name. A VCF source record's name is the ID column of its line. A VCF source of `type=exact` matches a source record to the input allele it shares, and the record's name prints on that allele's rows alone.

Ensembl VEP keeps the values of one output row in one table keyed by name, so a source key that spells one of the fixed columns shares that column's cell. Under the key `Consequence` the cell holds the row's terms, then the names of the allele's source records. vep-rs prints the same cell: the column is the place the user's key names. The report's count is another matter. It states how many consequences of each type the run found, and a source record's name is none.

### 1.2 The input record

The input record is line 9 of the `GRCh38-custom` test corpus's input (`tests/golden/116/GRCh38-custom/variants.vcf`), annotated against the release-116 GRCh38 cache:

```
21	5380280	rs1444031311	A	C,G	.	.	.
```

The corpus's `custom/sites.vcf.gz` holds one source record at that position. It is named `exact_2` and its alternate allele is the input's first. Both engines run with `--stats_text` and the corpus's file as a source keyed `Consequence`; section 1.6 gives the flag whole. No committed run of the corpus uses that key: the runs of this page are made for it.

### 1.3 What each engine writes

Three cells of the input record's two rows, the same in both engines' tab outputs:

| Location | Allele | Consequence |
| --- | --- | --- |
| `21:5380280` | `C` | `intergenic_variant,exact_2` |
| `21:5380280` | `G` | `intergenic_variant` |

The two consequence sections of the report, as Ensembl VEP 116.3 writes them in the default, tab and VCF formats:

```
[Consequences (most severe)]
intergenic_variant	1

[Consequences (all)]
exact_2	1
intergenic_variant	2
```

Ensembl VEP 116.3 in the JSON format, and vep-rs in all four formats:

```
[Consequences (most severe)]
intergenic_variant	1

[Consequences (all)]
intergenic_variant	2
```

Each of the three Ensembl VEP runs that count the name also writes one line to standard error, broken here to fit the page; the JSON run writes none:

```
Use of uninitialized value in numeric comparison (<=>) at
/opt/vep/src/ensembl-vep/modules/Bio/EnsEMBL/VEP/Stats.pm line 782,
<__ANONIO__> line 3.
```

### 1.4 Why Ensembl VEP writes that

Ensembl VEP builds one table of values for each output row. It first assigns the row's consequence terms to the entry `Consequence`. It then adds the custom annotations of the row's allele through one routine, which pushes each source record's name onto the list under the source's key. Under the key `Consequence` that list is the one the terms were assigned to, and the writers print it as the column's cell.

The statistics logger receives the same table once the row is built and counts every item of the entry. It makes no test of what an item is.

```
# OutputFactory.pm:1341, the row's terms
$hash->{Consequence} = [map {$_->$term_method || $_->SO_term} @ocs];
# OutputFactory.pm:1301, in the routine called at lines 1426-1435
push @{$hash->{$custom_name}}, $annot->{name};
# Stats.pm:282, the count
$stats->{consequences}->{$_}++ for @{$hash->{Consequence}};
```

The report sorts the section's lines by the severity rank of each consequence type. The name has no rank, so the comparison reads an undefined value and Perl warns. An undefined rank compares as zero, below every rank a consequence has, and the name's line is written first.

The JSON writer takes the custom annotations off the variant before a row's table is built and writes them under the `custom_annotations` member of the variant's JSON object. No name reaches a row's `Consequence` entry, so the report of a JSON run counts the terms alone.

Source (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1341, 1426-1435 and 1297-1312: the terms, then the push of the allele's source records under the source's key.
- The same file, line 401: each row's table is handed to the statistics logger.
- `Bio/EnsEMBL/VEP/Stats.pm`, line 282: the count of every item.
- The same file, lines 511-513, 535-543 and 782: the section, headed `Consequence type`, and the sort of its lines by consequence rank.
- The same file, lines 471-478: `[Coding consequences]` is read from the same counter, one consequence type at a time.
- The same file, line 241: the most severe consequence is read from the variant, never from a row's table.
- `Bio/EnsEMBL/VEP/OutputFactory/JSON.pm`, lines 250-251: the custom annotations are taken off the variant.

### 1.5 Why vep-rs's report is the correct one

The section's counts follow from the input record and the cache, without either engine. The cache's region file for bases 5,000,001 to 6,000,000 of chromosome 21 holds 71 transcripts. None lies within 5,000 bases of position 5380280, the reach of the upstream and downstream terms; the nearest starts 68,307 bases away. The run asks for no regulatory features. Each alternate therefore has one row with no feature, whose one term is `intergenic_variant`: two counts of that term and no other consequence type.

`exact_2` is the third column, the ID, of a line of `sites.vcf.gz` and identifies that line. Ensembl VEP's report agrees wherever it does not read a row's table. The most severe section, read from the variant, counts no name. The JSON run of the same input record counts none either, though the record's consequence terms are the same in every format.

### 1.6 Try it yourself

Write the input record of section 1.2 to `record.vcf` under a VCF header and run both engines with the text report. vep-rs's command is below, `<corpus>` being the corpus's directory; Ensembl VEP takes the same flags with the cache options of [Reproducing these examples](reproducing.md):

```
SOURCE="file=<corpus>/custom/sites.vcf.gz,short_name=Consequence"
SOURCE="$SOURCE,format=vcf,type=exact,fields=AF%CLASS"
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --stats_text --tab \
  --custom "$SOURCE"
```

The lines of `[Consequences (all)]` in each run:

| Run | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `--tab`, `--vcf` or the default format | `exact_2 1`, `intergenic_variant 2` | `intergenic_variant 2` |
| `--json` | `intergenic_variant 2` | the same |
| `--tab`, the source keyed `sites` | `intergenic_variant 2` | the same |

### 1.7 How the comparator treats this class

`crates/vep-cli/tests/statistics_report_under_a_custom_key_named_consequence.rs` holds vep-rs's sections of this page, since no committed run of the test corpora uses the key.

The repository's comparison of two statistics reports, `tools/vep-stats-compare`, holds such a line to this class. The rule applies where the `Command line options` row of Ensembl VEP's report names a `--custom` source keyed `Consequence`. A line of `[Consequences (all)]` that Ensembl VEP's report alone holds, with an integer as its count, is then set aside when its name is no consequence type in any vocabulary `--terms` chooses: `SO`, `display` or `NCBI`. The comparison's report lists the line under `documented` with the class's name.

Ensembl VEP masks the directories in the `Command line options` row. Where a second source is given by position under a directory, the mask takes its `--custom` word and the value before it; the rule then reads no key there, and the line is reported.

Of the source, the rule reads the key alone, neither its format nor its type. Under a source of another kind Ensembl VEP's section holds no name, so there is no line to set aside. Every other difference is reported: a consequence type that one report lacks, a count that differs, a name that vep-rs's report holds and a name in any other section are compared as written.

## 2. An input record whose named allele has nine rows

**In short.** On the input record `21 10542415 multi_eq_ref_second A G,A`, the source record `exact_3` names the first alternate, which overlaps nine transcripts. Ensembl VEP 116.3 writes `exact_3 9` in `[Consequences (all)]` of a tab run: one count for each of the allele's nine rows. vep-rs counts the rows' terms alone, and every other line of the two reports agrees outside `[VEP run statistics]`.

The input record is line 16 of the same corpus input, and the source and the flags are those of section 1. Its second alternate equals the reference, and both engines write no row for it. `[Consequences (all)]` as Ensembl VEP 116.3 writes it, the name's line first:

```
[Consequences (all)]
exact_3	9
missense_variant	7
intron_variant	2
NMD_transcript_variant	2
```

vep-rs writes the three terms' lines with the same counts and no line for the name. In both outputs the nine rows of the first alternate end their `Consequence` cell with `exact_3`: seven with `missense_variant`, two of those with `NMD_transcript_variant` as well, and two with `intron_variant`. The push of section 1.4 is made once for each row, so one input record adds nine to the name's count.

## 3. The edges of the class

The name is counted only where a source record is attached to one allele and its name stands in a row's table when the logger reads it. In the JSON run of section 1 the names are taken off before any table is built. The three cases below are read from Ensembl VEP's source at release/116.3 and were not run; in each, the two engines' sections agree.

- A source of `type=overlap`, and any BED, GFF, GTF or bigWig source. Its records belong to the whole input record. Their names are pushed before the rows are built (`OutputFactory.pm`, lines 1004-1014), and line 1341 then assigns each row's terms over them, so in both engines the cell and the section hold the terms alone. Only an exact match of a VCF source names an allele (`AnnotationSource/File/VCF.pm`, lines 228-246 and 405-421).
- A structural variant as the input record or as the source record. Its match is handed to a routine that names no allele (`AnnotationSource/File/VCF.pm`, lines 380-387, and `AnnotationSource/File.pm`, lines 566-575). The rows of a structural variant take no allele's source records (`OutputFactory.pm`, lines 2032-2061).
- A run under `--summary` or `--most_severe`, which logs no row (`OutputFactory.pm`, lines 387 and 426-459).

One shape inside the class is left as a difference; it too is read from the source and was not run. A source record whose name spells a consequence type, `missense_variant` say, raises that type's count in Ensembl VEP's section. It also raises the type's count in `[Coding consequences]`, the section of the coding types (`Stats.pm`, lines 471-478). The comparison cannot tell such a count from a miscount, so it reports it.
