# The description of a custom source's `PC` field

Two header-line differences of a `--custom` source: the `PC` description below, and the line that receives a repeated short name's file ([section 2](#2-the-header-line-that-receives-a-repeated-custom-short-names-file)). In the first, Ensembl VEP's header describes the field, under either setting of `reciprocal`, with a sentence that does not name the share the value is, and vep-rs's names that share. In the second, Ensembl VEP appends the file to the first custom header line unless the key holds a lowercase `o`, and vep-rs appends it to the key's own line, on the text of the declaration the second source repeats.

**The defect in one sentence.** Under either setting of `reciprocal`, Ensembl VEP's header describes the `PC` field of a `--custom` source with a sentence that does not name the share the value is. Without `reciprocal` it announces the share of the input variant covered while the value is the share of the annotation record covered; with `reciprocal=1` it announces a minimum reciprocal overlap while the value is the share of the input variant covered.

**Why this is a defect.** The header line is the one place an output explains what a number means. The numbers themselves are right in both settings and agree between the engines; the two sentences describing them are each attached to the wrong setting. A reader who trusts the header reads a record 66.667 percent covered by the variant as a variant 66.667 percent covered by the record, and reads a one-sided share as a minimum of two shares that Ensembl VEP never computes.

**What vep-rs does instead.** vep-rs's header names the share the value is: `Percentage of the annotation record covered by the input variant` without `reciprocal`, and `Percentage of the input variant covered by the annotation record` under `reciprocal=1`. The values are Ensembl VEP's.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/AnnotationSourceAdaptor.pm`, lines 395-397: `overlap_def` is set to `Percentage of minimum reciprocal overlap between input variant and reference variant` when `reciprocal` is true and to `Percentage of input variant covered by reference variant` otherwise.
- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 536-564: the share computed. `$length` at line 538 is the input variant's length, and line 554 divides the overlap by it. Lines 556-563 replace that share, when `reciprocal` is false, by the overlap divided by the annotation record's length, under the comment `when reciprocal = 0, check percentage of reference (annotation) variant covered`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2485-2489: the `<key>_PC` header line prints `overlap_def` followed by ` from <file>`.

**What correcting it would change.** The class vanishes: each setting's header would name the share its value is, as vep-rs's does, and no value would move.

**The class it produces.** `pc_description_text`, applied by the golden harness and by the format comparator to the reference's header lines before the comparison. The concordance figures score consequence sets and never read a header line, so nothing is set aside in them.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. The `PC` header names the other share

**In short.** A `--custom` source asked for the `PC` field prints, for each record it matches, a percentage of overlap. Ensembl VEP computes the share of the annotation record the variant covers, or under `reciprocal=1` the share of the variant the record covers, and both engines print the same numbers. Ensembl VEP's header line describes neither setting's share: the default's sentence names the share `reciprocal=1` computes, and the `reciprocal=1` sentence names a minimum reciprocal overlap that no setting computes. vep-rs's header describes the share printed. On a two-base deletion matched by a three-base record and a one-base record, the values `66.667,100.000` and `100.000,50.000` show which share each setting reports.

### 1.1 What a reader needs to know first

A `--custom` source of `type=overlap` matches an annotation record to a variant when their spans overlap. Its `fields=PC` option adds a field, `<short_name>_PC`, holding one percentage per matched record, and its `reciprocal=1` option switches which span the percentage is taken over. Each `--custom` source adds a header line per field, `## <key> : <description>` in the default and tab formats and `##INFO=<ID=<key>,...,Description="<description>">` in the VCF, and the `PC` line's description is the sentence this page is about. The JSON format carries the values and no header.

A percentage of overlap has two readings: the fraction of the annotation record that the variant covers, and the fraction of the variant that the record covers. They differ whenever the two spans differ in length.

### 1.2 The record

The record is `21 5028997 . TCA T` from `variants.vcf` of the `GRCh38-custom` golden corpus (`tests/golden/116/GRCh38-custom/`), a deletion of the two bases 5028998-5028999:

```
21	5028997	.	TCA	T	.	.	.
```

The corpus's `custom/regions.bed.gz` holds two records that overlap it. In BED's half-open coordinates and then in the one-based coordinates both engines compute with:

```
21	5028996	5028999	exact_1	1	+
21	5028998	5028999	overlap_1	2	+
```

`exact_1` spans 5028997-5028999, three bases, and covers both deleted bases; `overlap_1` spans 5028999, one base, and covers one of the two.

### 1.3 What each engine writes

The source without `reciprocal`, `--custom file=<corpus>/custom/regions.bed.gz,short_name=regions,format=bed,type=overlap,fields=PC`. Ensembl VEP's header line and row:

```
## regions_PC : Percentage of input variant covered by reference variant from [PATH]/regions.bed.gz
21_5028997_TCA/T	21:5028998-5028999	-	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;regions=exact_1,overlap_1;regions_PC=66.667,100.000
```

vep-rs:

```
## regions_PC : Percentage of the annotation record covered by the input variant from [PATH]/regions.bed.gz
21_5028997_TCA/T	21:5028998-5028999	-	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;regions=exact_1,overlap_1;regions_PC=66.667,100.000
```

The same source with `reciprocal=1`, the `bed-overlap-stats` flag set of the corpus's manifest, which also asks for the summary statistics. Ensembl VEP:

```
## regions_PC : Percentage of minimum reciprocal overlap between input variant and reference variant from [PATH]/regions.bed.gz
21_5028997_TCA/T	21:5028998-5028999	-	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;regions=exact_1,overlap_1;regions_PC=100.000,50.000;regions_count=2;regions_max=2;regions_mean=1.5;regions_min=1;regions_sum=3
```

vep-rs:

```
## regions_PC : Percentage of the input variant covered by the annotation record from [PATH]/regions.bed.gz
21_5028997_TCA/T	21:5028998-5028999	-	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;regions=exact_1,overlap_1;regions_PC=100.000,50.000;regions_count=2;regions_max=2;regions_mean=1.5;regions_min=1;regions_sum=3
```

In both settings the rows agree in every value and the header lines differ in the sentence before ` from`. The VCF's `##INFO=<ID=regions_PC,...,Description="...">` line carries the same two sentences.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP computes the overlap of the variant and the record as a count of shared bases and divides it by one of the two lengths. Without `reciprocal` it divides by the record's length: `exact_1` shares two of its three bases with the deletion, 66.667 percent, and `overlap_1` shares its one base, 100 percent. With `reciprocal=1` it divides by the variant's length of two: `exact_1` covers both deleted bases, 100 percent, and `overlap_1` one of them, 50 percent. The code comment above the default branch states the share it takes: the percentage of the reference, that is the annotation, record covered.

The header's sentence is chosen elsewhere, where the source's options are parsed, by the same `reciprocal` switch, and the two sentences are attached the other way round. The default sentence announces the share of the input variant covered, which is the share the `reciprocal=1` branch computes. The `reciprocal=1` sentence announces a minimum reciprocal overlap, a smaller-of-two-shares measure that no branch computes: for `exact_1` it would read 66.667 where 100.000 is printed.

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 536-564: `$length` (the variant's) at line 538, the shared bases over it at line 554, and the `unless ($reciprocal)` branch at lines 556-563 that divides by the record's length instead.
- `Bio/EnsEMBL/VEP/AnnotationSourceAdaptor.pm`, lines 395-397: the two sentences, chosen by `reciprocal`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2485-2489: the `<key>_PC` header line.

### 1.5 Why vep-rs's row is the correct one

The value is fixed by the arithmetic and is the same in both engines; the question is which sentence describes it. Without `reciprocal` the printed `66.667` is two bases of `exact_1`'s three, a share of the annotation record, and vep-rs's header says so. Under `reciprocal=1` the printed `100.000` is two bases of the deletion's two, a share of the input variant, and vep-rs's header says so. The option's own name agrees: `reciprocal` turns the measure round, from the record's side to the variant's, and the header should turn with it.

Ensembl VEP's default sentence, `Percentage of input variant covered by reference variant`, would read `100.000` for `exact_1`, not the `66.667` printed, and its `reciprocal=1` sentence names a minimum of the two shares that would read `66.667` where `100.000` is printed. Neither sentence can describe the number under it.

### 1.6 Try it yourself

Run the record of 1.2 through both engines with the corpus's `regions.bed.gz` as the source, once without `reciprocal` and once with `reciprocal=1`; vep-rs's command is below, and Ensembl VEP takes the same `--custom` string with the cache and FASTA options of [Reproducing these examples](reproducing.md):

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --fasta <corpus>/reference.fa.gz --species homo_sapiens --assembly GRCh38 \
  --custom file=<corpus>/custom/regions.bed.gz,short_name=regions,format=bed,type=overlap,fields=PC
```

| Setting | `regions_PC` (both engines) | Ensembl VEP's header | vep-rs's header |
| --- | --- | --- | --- |
| default | `66.667,100.000` | `Percentage of input variant covered by reference variant` | `Percentage of the annotation record covered by the input variant` |
| `reciprocal=1` | `100.000,50.000` | `Percentage of minimum reciprocal overlap between input variant and reference variant` | `Percentage of the input variant covered by the annotation record` |

Checking the default row against the spans of 1.2: two of `exact_1`'s three bases is 66.667 percent, and both deleted bases out of two would be 100 percent, so the default value is the record's share, as vep-rs's header reads.

### 1.7 How the comparator treats this class

The golden harness (`crates/vep-cli/tests/common/mod.rs`, `masked_reference_text`) and the format comparator (`tools/vep-formats-compare`, `reference_text_masked`) apply the class `pc_description_text` when a corpus flag set or a comparison case names it in `masks`. Under the class a reference header line that describes a `<key>_PC` field with either of Ensembl VEP's two sentences is read with vep-rs's sentence for the same setting. The rewrite covers the `## <key>_PC : ` shape and the VCF `##INFO` shape and runs before the header lines are compared. The comparator counts each rewritten line under the class in its report. A row's values compare as any other value, and a header line describing any other key is untouched.

The `bed-overlap-stats` flag set of the `GRCh38-custom` corpus names the class, and so does the comparison case that asks a BED source for `PC`. The class rewrites nothing in the JSON format, which carries no header, and a run that does not name it compares the header line as written.

## 2. The header line that receives a repeated custom short name's file

A second header-line difference of a custom source, a documented class, as [Intended divergences from Ensembl VEP](../intended-divergences.md) lists it. The class is documented only: no committed dataset reaches its shape. The class carries its tier under that page's rule: tier 1 is a shape the chromosome 21 cuts carry with no flag, tier 2 one that needs a flag or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts are over the cuts' 15,646,481 records, the 152 golden corpora and their 511 flag sets.

**The defect in one sentence.** When two `--custom` sources share a short name, Ensembl VEP 116.3 appends the second source's file to the first custom header line of the output instead of the shared key's own line, unless the key happens to contain a lowercase `o`.

**Why this is a defect.** The header is the one place an output names the file a key's values came from. `get_custom_headers` indexes the header keys by position but keeps only the entries that match the regular expression `/o/`, a test on a letter, not on a key. A key without a lowercase `o` has no index, and `|| 0` and `undef / 2` both read as line 0, so the file lands on whichever key was listed first. Every field or statistic key of the repeated name adds a `Use of uninitialized value` warning. The same run with a key holding an `o` lists both files on the key's own line, the shape the code intends.

On the key's own line Ensembl VEP appends the file to the text that line has, whatever the second source declares the key as. The test that finds a listed key reads the key alone. A record field and a statistic of one name share a key, so a statistic that the second source asks for joins the first source's field text as one more file, and the statistic is described in no header ([section 2.2](#22-a-key-the-two-sources-declare-apart-as-a-record-field-and-as-a-statistic)).

**What vep-rs does instead.** vep-rs finds the key's own header line, whatever its letters, and appends `,<file>` there, for the short name and for each of its field and statistic keys; the VCF `##INFO` descriptions follow the same lines. A field's file joins the field's text and a statistic's the statistic's. A declaration the key's texts lack is described by a text of its own.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 2452, `get_custom_headers`: `%pos` is built from the flattened `[key, description]` pairs, keeping only the entries that match `/o/`.
- The same file, line 2458: the short name's line is `($pos{$short_name} || 0) / 2`; lines 2483 and 2500: a field's or statistic's line is `$pos{$sub_id} / 2`, where an undefined position reads as 0.
- The same file, line 2451: the keys and texts listed so far are read once a source, before the source's own are pushed. Lines 2482 and 2499 test a field's and a statistic's key against that reading. Where the key is found, lines 2484 and 2501 append the file to the later of the key's listed texts, whatever it declares, and the text of the second source's declaration (lines 2493 and 2503-2506) is not pushed.
- `Bio/EnsEMBL/VEP/OutputFactory/VCF.pm`, line 276: the same pairs become the `##INFO` lines, so the VCF header carries the misplaced file too.

**What correcting it would change.** The second file would join the shared key's own line in every run, as it does for a key holding an `o`, a declaration of another kind would bring its own text, and no row would change.

**The class it produces.** `custom_repeated_short_name_header_lines`: on a run naming one short name on two or more sources, a `,<file>` that the reference's `## <key> : <text>` and `##INFO` header lines carry on another key's line; the rows agree. Where the sources declare a key apart, the reference lacks the second declaration's text and, for a key holding an `o`, carries the file on the text of the other declaration. The class is documented only: no committed golden corpus and no comparison case names a short name twice.

**How often the shape occurs.** Tier 3: the shape is two sources sharing a short name that holds no lowercase `o`, which none of the 511 flag sets does; 19 of them run a `--custom` source. The shape of section 2.2 needs two sources of one short name as well, and a record field whose name is a statistic's (`count`, `min`, `max`, `mean` or `sum`).

### 2.1 A short name `SITES` on two VCF sources beside a BED source named `cons`

**In short.** A run names a BED source `cons`, then two VCF sources `sites.vcf.gz` and `clinvar_source.vcf.gz` under one short name `SITES` with `type=exact,fields=AF`. Ensembl VEP 116.3 prints `## cons : [PATH]/regions.bed.gz,[PATH]/clinvar_source.vcf.gz,[PATH]/clinvar_source.vcf.gz`, once for the `SITES` key and once for `SITES_AF`, while `## SITES : [PATH]/sites.vcf.gz` and `## SITES_AF` name one file, and warns `Use of uninitialized value $pos{"SITES_AF"} in division (/)`. vep-rs prints `## SITES : [PATH]/sites.vcf.gz,[PATH]/clinvar_source.vcf.gz`, the same two files on `## SITES_AF`, and `cons` with its BED alone. The rows agree in the default, tab and VCF formats, and with the shared name spelt `solo`, which holds an `o`, both engines print the same header lines.

### 2.2 A key the two sources declare apart, as a record field and as a statistic

**In short.** Two VCF sources named `S` each match one record at `21:5380280`. The first names the record field `count`. The second asks for the statistic `count`, the number of records it matched. Both declare the key `S_count`. Ensembl VEP 116.3 describes it with the first source's field text and appends the second file, over a cell that holds the second source's statistic, `1`. In its tab and VCF headers vep-rs describes the statistic and then the field, each with its own file: `count data from [PATH]/b.vcf.gz; The first file's made-up field. count field from [PATH]/a.vcf.gz`. The rows agree.

**The files and the command.** The input record is `21 5380280 rs1444031311 A C,G` of the `GRCh38-custom` corpus. Two made-up site files hold one record each at that position: `exact_2` with `count=7` in `a.vcf.gz`, and `second_1` with `count=9` in `b.vcf.gz`. Each record names the allele `C`, and neither has a `QUAL`. The first file's header describes `count` as `The first file's made-up field`, and the second file's as `The second file's made-up field`.

```
--custom file=a.vcf.gz,short_name=S,format=vcf,type=overlap,fields=count
--custom file=b.vcf.gz,short_name=S,format=vcf,type=overlap,\
summary_stats=count
```

**What each engine writes.** `A` and `B` stand for `[PATH]/a.vcf.gz` and `[PATH]/b.vcf.gz`, and `F` for the field's text, `The first file's made-up field. count field from A`.

| Header | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Tab | `## S_count : F,B` | `## S_count : count data from B; F` |
| VCF | one `##INFO` line of `S_count`, with `F,B` | one line, with `count data from B; F` |
| Default | `## S_count : F,B` | `## S_count : F`, then `## S_count : count data from B` |

vep-rs lists a key once in its tab and VCF headers, so its two texts share a line there. Its default header prints the texts on separate lines, in the order the sources declare them. The entry's reproduction record writes the two files and runs the command in the tab and default formats.

Both engines write `## S : A,B`: the short name's line is the first custom line and `S_count` holds an `o`, so the lowercase-`o` rule of this section moves no file in this run. On the rows of both alleles, `C` and `G`, both write `exact_2,second_1` under `S` and `1` under `S_count`. The JSON record holds both the statistic and the `count` of the first file's record, `7`, in both engines.

**Why Ensembl VEP writes that.** It lists the first source's key with the field's text. Before it reads the second source it takes the keys listed so far, the first source's. It tests the second source's statistic key against them, finds it, and appends the file to the listed text. The branch that pushes the statistic's own text is not taken. A source's own keys are not among those it is tested against, so a source that declares a key twice lists it twice, and a later source's file joins the later listing.

**Why vep-rs's lines are the correct ones.** The second source matched one record, so its count is `1`, and under `type=overlap` the statistic replaces the list of the field's values under the key ([A key that spells a source's statistic](custom-key-named-as-a-plain-value-field.md#3-a-key-that-spells-a-sources-statistic)). The cell therefore holds a count of records. Ensembl VEP's line calls it the values of a record field and names `b.vcf.gz` as a file of that field, which the second source does not name.

A text for each declaration, each with the files of the sources that make it, is true of the cell whichever value it holds. On the one line of the tab and VCF headers the statistic's text stands first, the order of the two values on a row that holds both. Under `type=overlap` a record belongs to the whole input record, and the statistic replaces its value. Under `type=exact` a record belongs to the allele it names, `C`, and its value is listed after the statistic: vep-rs's row of `C` reads `1,7`. Ensembl VEP ends that run before its first row, the documented difference of [that page](custom-key-named-as-a-plain-value-field.md#3-a-key-that-spells-a-sources-statistic).

**The other orders.** The same two files under other declarations; `F2` stands for the second file's field text, `The second file's made-up field. count field from B`, and a `/` parts two lines of the VCF and default headers:

| First source, second source | Ensembl VEP, VCF and default | vep-rs, tab and VCF |
| --- | --- | --- |
| field; field and statistic | `F,B,B` | `count data from B; F,B` |
| statistic; field | `count data from A,B` | `count data from A; F2` |
| statistic; field and statistic | `count data from A,B,B` | `count data from A,B; F2` |
| field and statistic; field | `F` / `count data from A,B` | `count data from A; F,B` |
| field and statistic; both | `F` / `count data from A,B,B` | `count data from A,B; F,B` |
| field and statistic; statistic | `F` / `count data from A,B` | `count data from A,B; F` |

vep-rs's default header prints its texts on separate lines, in the order the sources declare them. In the last row the engines print the same two texts with the same files: the second source repeats the statistic, whose text takes its file. The default headers agree there, and the VCF header holds the two texts on two lines in Ensembl VEP and on one in vep-rs. Where the first source declares the key twice, Ensembl VEP's tab header prints the later text on both of the key's lines ([the description of a key one source declares twice](custom-key-named-as-a-plain-value-field.md#36-the-description-of-a-key-one-source-declares-twice)).

A statistic that needs a score shows why one text cannot stand alone. Give each record the field `min` in place of `count`, and ask for `summary_stats=min`. A VCF record's score is its `QUAL`. With a `QUAL` on the records the cell holds the asking source's statistic. With none it holds the field's values: `7` under the command above, and `7,9` where the second source names the field as well.

vep-rs's line carries both texts in either run. The key `S_min` holds no `o`, so Ensembl VEP appends the second file to the `## S` line a second time, for the statistic the second source asks for: `## S : A,B,B`. Its `## S_min` line keeps the first source's text alone.
