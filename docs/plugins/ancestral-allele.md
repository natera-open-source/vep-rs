# AncestralAllele: the ancestral allele from Ensembl's EPO alignments

AncestralAllele writes the ancestral allele at a variant's position, read from
the FASTA of ancestral sequences Ensembl derives from its primate EPO
alignments. It is a record-level plugin: one `AA` value per input record,
printed on every row of the record, intergenic rows included, or one per allele
of a record `--minimal` splits. This page covers its one parameter, the data
file, how the special values arise and the two places it differs from Ensembl
VEP's plugin.

## Parameters and data {#ancestralallele}

| Property      | Value                                                                                                                                                                                                                                    |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                                                                                                                                           |
| Matches       | position: the record's span, whatever its alleles                                                                                                                                                                                        |
| Output fields | `AA`                                                                                                                                                                                                                                     |
| Data source   | Ensembl ancestral alleles, [GRCh38](https://ftp.ensembl.org/pub/release-116/fasta/ancestral_alleles/) (`homo_sapiens_ancestor_GRCh38.tar.gz`) and [GRCh37](https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/) (`homo_sapiens_ancestor_GRCh37_e71.tar.bz2`) |

The plugin takes the FASTA as its only parameter, by position or as `file=`.
A positional parameter after it is not read, as Ensembl VEP's plugin reads
the first alone, and is noted in the warnings file:

```bash
--plugin AncestralAllele,/path/to/homo_sapiens_ancestor_GRCh38.fa.gz
--plugin AncestralAllele,file=/path/to/homo_sapiens_ancestor_GRCh38.fa.gz
```

The FASTA is indexed: a plain `.fa` with its `.fai`, or the bgzipped `.fa.gz`
with its `.fai` and `.gzi`. The bgzipped form is the one Ensembl VEP's own
plugin reads through htslib. A gzip file that is not BGZF is refused, since its index offsets cannot
be reached. The header line names the file as it was given:
`## AA : Ancestral allele from <file>`.

**Preparing the data.** Ensembl ships one FASTA per chromosome inside a tar
archive. Concatenate them into one bgzipped file and index it; `samtools faidx`
writes the `.fai` and the `.gzi` together:

```bash
tar xfz homo_sapiens_ancestor_GRCh38.tar.gz
cat homo_sapiens_ancestor_GRCh38/*.fa | bgzip -c \
    > homo_sapiens_ancestor_GRCh38.fa.gz
samtools faidx homo_sapiens_ancestor_GRCh38.fa.gz
```

The GRCh38 archive is about 860 MB and the bgzipped FASTA about 830 MB.
Ensembl places no restriction on access to or use of its data, so the file may
be copied and redistributed.

## The value

The records of the FASTA are named as Ensembl writes them,
`ANCESTOR_for_chromosome:GRCh38:21:1:46709983:1`. The chromosome is the third
colon-separated field. A record's chromosome reaches it under its own name,
with a `chr` prefix removed, or with `chr` or `CHR` added. The value is the
FASTA's bases over the record's span, in the case the file gives them, which
carries Ensembl's confidence coding:

| Characters | Meaning                                                                   |
| ---------- | ------------------------------------------------------------------------- |
| `ACGT`     | high-confidence call: the ancestral state is supported by two sequences   |
| `acgt`     | low-confidence call: supported by one sequence                            |
| `N`        | no call: the ancestral state is supported by no other sequence            |
| `-`        | the extant species carries an insertion at this position                  |
| `.`        | no coverage in the alignment                                              |

A deletion or a multi-base substitution reports the bases of its whole span. A
record given on the reverse strand (an Ensembl-format or region input) reports
their reverse complement. Three values stand in where the file cannot answer, as
Ensembl VEP's plugin writes them:

| Value            | When                                                                                                        |
| ---------------- | ----------------------------------------------------------------------------------------------------------- |
| `-`              | an insertion: its span is empty, so no ancestral base exists for it                                         |
| `?`              | the record's chromosome has no record in the FASTA                                                          |
| the empty string | the record's span lies wholly past the end of its chromosome's FASTA record; a span that runs past the end reports the bases up to it |

The empty string prints as `AA=` in the default format, an empty cell in the
tab format and `""` in JSON. Ensembl VEP writes the same from a bgzipped FASTA;
from an uncompressed one it writes the record's last base, the class of
[section 2](#2-a-span-past-the-end-of-an-uncompressed-fastas-record).

A structural variant (a symbolic allele such as `<DEL>`) gets no `AA`: Ensembl
VEP hands the plugin sequence variants alone, and vep-rs does the same.

## Where AncestralAllele differs from Ensembl VEP

Two behaviours of Ensembl VEP 116.3's plugin are not reproduced: one under
`--minimal` alone, one on an uncompressed FASTA whose record ends before a
variant's span. On a bgzipped FASTA without the flag the two engines write the
same `AA` on every row. Each is a documented class: its section names the rows
a comparison would set aside. No comparison tool codes the class, since no
committed row falls under it, as
[Intended divergences from Ensembl VEP](../intended-divergences.md) lists it.
Its [How to read the counts](../intended-divergences.md#how-to-read-the-counts)
defines the terms.

The classes carry their tier under that page's rule:

- tier 1 is a shape the chromosome 21 cuts carry with no flag;
- tier 2 is one that needs a flag or plugin on ordinary data;
- tier 3 is one ordinary data does not carry.

The counts are over the cuts' 15,646,481 records, the 152 golden corpora and
their 511 flag sets.

### 1. The ancestral allele printed for a later allele of a `--minimal` split record

**The defect in one sentence.** Under `--minimal`, Ensembl VEP prints the first
alternate allele's `AA` on every allele's rows of a split multi-allelic record,
so a later allele carries the ancestral allele of a span that is not its own.

**Why this is a defect.** The plugin's own rule values an allele by the span it
reads: the bases over it, or `-` when the span is empty. Ensembl VEP computes
every split allele's value that way on its own trimmed span and caches it on
the allele's feature. It then rejoins the record by pointing each later
allele's overlaps at the first allele's feature, where the first value is
cached, so the later values are never read. The later allele's `AA` matches
neither its own span nor the record's whole span, which both engines print
without the flag.

**What vep-rs does instead.** vep-rs values every split allele over its own
trimmed span. A split record's rows therefore agree with Ensembl VEP's on the
first allele and differ on a later allele whenever the two spans read differently.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3,
`VEP_plugins` at release/116):

- `AncestralAllele.pm`, lines 116-124, 148 and 151, `run`: the value is read
  over the variation feature's own `start` and `end` (`-` when `start` exceeds
  `end`) and cached on that feature; a feature already carrying a value returns
  it.
- `Bio/EnsEMBL/VEP/InputBuffer.pm`, lines 638-709, `split_variants`: a
  multi-allelic record one of whose alleles trims becomes one feature per allele
  on its trimmed span. The first is marked `first`, the others `merge_with`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2298 and 2312,
  `rejoin_variants_in_InputBuffer`: each split feature's rows are built once, so
  each caches its own value. Lines 2300-2302 restore the first feature's
  allele string and span. Line 2323 points every later allele's overlaps at the
  first feature, so the output pass reads its cached value for every allele.

**What correcting it would change.** Every allele's rows would carry the
ancestral allele of their own trimmed span, as vep-rs's rows do. Ensembl
VEP's rows do the same for a bi-allelic record, which `--minimal` trims without
a rejoin.

#### The class for a split record and how often it occurs

**The class it produces.** `minimal_plugin_value_follows_first_allele`: an `AA`
(the VCF `CSQ` field `AA`, the JSON `aa`) that differs on a row of a `MINIMISED`
record whose allele is not the record's first alternate allele. Ensembl VEP's
value there equals its own `AA` on the first allele's rows. A first-allele row stays
a difference. No row is set aside on the published datasets: the plugin's
corpus has two multi-allelic records and no `--minimal` run, and no
`--minimal` run carries the plugin.

AncestralAllele is the one plugin that caches
its value on the variation feature, so the class covers it alone; a plugin that
reads the rejoined record's span or alleles without caching is the sibling class
[`minimal_plugin_sees_the_rejoined_record`](../divergences/minimal-plugin-sees-the-rejoined-record.md).

**How often the shape occurs.** Tier 2: 1,422 of the cuts' records
(91 per million) are multi-allelic records one of whose alleles trims beyond the
anchor, so that `--minimal` splits and rejoins them. The ancestral sequence
decides the values. 22 corpora carry the shape and no committed flag set runs
the plugin under `--minimal`.

#### 1.1 A two-allele indel whose second allele trims to an insertion

**In short.** The record `21 42750000 anc_multi_indel TG T,TGT` of the GRCh38
AncestralAllele corpus (`tests/golden/116/GRCh38-plugin-ancestral/variants.vcf`),
run with the corpus's flags and `--minimal`, is `G/-/GT` at 21:42750001 once
anchored. Its `GT` allele trims to the insertion of `T` after 21:42750001, an
empty span, which the plugin's own rule values `-`. Ensembl VEP 116.3 writes
`AA=G` on all 120 rows of the record; vep-rs writes `AA=G` on the 60 rows of the
`-` allele and `AA=-` on the 60 rows of the `T` allele, every other column
agreeing. Without `--minimal` both engines write `AA=G` on every row, the
record's one span being 21:42750001.

On `21 42600001 TC CC,T`, whose later allele is the deletion of the `C` at
21:42600002, Ensembl VEP 116.3 writes `AA=T` on that allele's rows, the first
allele's base at 21:42600001. The allele's own span reads `C`, which vep-rs
writes, and the record's whole span `TC`, which both engines write on every row
without `--minimal`.

### 2. A span past the end of an uncompressed FASTA's record

**The defect in one sentence.** From an uncompressed ancestral FASTA, Ensembl
VEP 116.3 reads a span that lies wholly past the end of the chromosome's record
as that record's last base, one base whatever the span's length, where the
bgzipped form of the same file gives the empty string.

**Why this is a defect.** The record ends where it ends, and a position past it
has no ancestral base in the file. The plugin's own rule for the bgzipped form,
and for a span straddling the end, is to report the bases the file holds over
the span and nothing more. The last base of the record is the ancestral state of
another position, and a one-base value for a three-base deletion describes no
span at all. The value depends on how the file was compressed, not on what it
says.

**What vep-rs does instead.** vep-rs reads a plain `.fa` and a bgzipped
`.fa.gz` alike. A span wholly past the record's end is the empty string, and a
span that runs past the end reports the bases up to it. These are the values
the bgzipped form gives in both engines.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, BioPerl
1.6.924, the release Ensembl VEP 116.3 installs):

- `AncestralAllele.pm`, line 95: a path ending `.gz` is read through
  `Bio::DB::HTS::Faidx`, any other through `Bio::DB::Fasta`. Line 137 asks
  Faidx for `chromosome:start-end`, which htslib cuts at the record's end, empty
  when wholly past it. Line 140 asks `Bio::DB::Fasta` for
  `seq($name, $start => $end)`.
- `Bio/DB/Fasta.pm`, lines 280-303, `subseq`: the file offsets of the first and
  last residue come from `_calc_offset`. Line 292 reads
  `$filestop - $filestart + 1` bytes between them.
- `Bio/DB/IndexedBase.pm`, lines 772-781, `_calc_offset`: line 779 raises a
  residue index below the first residue to the first. Line 780 lowers one at
  or past the record's length to the last residue. A span wholly past the end
  has both its ends lowered to the last residue, so the read is one byte, the
  record's last base. A span straddling the end is cut there.

**What correcting it would change.** The uncompressed path would report the
bases the record holds over the span, the empty string when it holds none, as
the Faidx path and vep-rs do. The value would be the same whatever the
file's compression.

#### The class for an uncompressed FASTA and how often it occurs

**The class it produces.** `ancestral_allele_uncompressed_fasta_clamp`: an `AA`
(the VCF `CSQ` field `AA`, the JSON `aa`) on a record whose span lies wholly
past the end of its chromosome's record in an uncompressed FASTA. Ensembl VEP's
value is the record's last base and vep-rs's the empty string.

**How often the shape occurs.** Tier 3: the shape is an uncompressed FASTA whose
chromosome record ends before a variant's span, which no committed flag set
carries. Both of the 2 flag sets that run the plugin name a bgzipped FASTA, and
no committed record lies past the end of its record.

#### A worked example on a FASTA cut short

**In short.** The corpus's chromosome 21 record, cut to its first 14,439,400
bases and written uncompressed, ends on the `C` at 21:14439400. Run against it,
a three-base deletion 301 kb past the end and a substitution 4.8 Mb past it get
`AA=C` from Ensembl VEP 116.3 on every row and `AA=` from vep-rs. The corpus's
`anc_intergenic_snv` at the last base gets `AA=C` from both, and the same cut
bgzipped gives vep-rs's values in both engines.

#### 2.1 What a reader needs to know first

The plugin picks its reader by the file's suffix. A `.gz` path is read through
htslib's indexed access. Any other is read through the pure-Perl indexed
reader, which writes its own index beside the file on first use. vep-rs reads
both forms through the `.fai` the file is indexed with.

The three records below use the corpus's cache and flags with the cut record
in place of the committed FASTA. The compressed control is the same cut
compressed with bgzip and indexed.

#### 2.2 The record

The corpus's chromosome 21 record cut to its first 14,439,400 bases, the
committed FASTA's first 144,395 lines written uncompressed as `cut.fa`, with the
`.fai` line:

```
ANCESTOR_for_chromosome:GRCh38:21:1:46709983:1	14439400	48	100	101
```

The record ends on the `C` at 21:14439400. The input is the corpus's
`anc_intergenic_snv` at that last base, a three-base deletion 301 kb past the
end, and the corpus's `rs937007124` 4.8 Mb past it:

```
21	14439400	anc_intergenic_snv	C	T	.	.	.
21	14740399	anc_three_base_deletion	CTAG	C	.	.	.
21	19268555	rs937007124	C	T	.	.	.
```

Both engines run `--plugin AncestralAllele,cut.fa`.

#### 2.3 What each engine writes

Ensembl VEP 116.3 writes `AA=C` on every row of all three records. Its rows for
the deletion (the first of its twenty-seven) and the later substitution:

```
anc_three_base_deletion	21:14740400-14740402	-	ENSG00000232884	ENST00000660333	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;AA=C
rs937007124	21:19268555	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;AA=C
```

vep-rs writes `AA=C` on `anc_intergenic_snv`, whose span is the record's last
base, and the empty string on the other two:

```
anc_three_base_deletion	21:14740400-14740402	-	ENSG00000232884	ENST00000660333	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;AA=
rs937007124	21:19268555	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;AA=
```

On the two records past the end the tab output carries `C` against the empty
cell, the VCF output `C` against an empty `AA` field, and the JSON output
`"aa":"C"` against `"aa":""`. Neither engine warns. On the same cut compressed
with bgzip, both engines write `AA=C` on the first record and `AA=` on the other
two.

#### 2.4 Why Ensembl VEP writes that

For the deletion the plugin asks the pure-Perl reader for residues 14,740,400
to 14,740,402 of a record 14,439,400 residues long. The reader lowers each index
at or past the length to the last residue. The first and last residue of the
span are then both the record's last. The read is one byte, and the value is
`C`, the ancestral base at 21:14439400. The substitution at 21:19268555 reads
the same byte.

From the compressed form the plugin asks htslib for the same regions,
which cuts them at the record's end and returns nothing.

#### 2.5 Why vep-rs's row is the correct one

The file has no base at 21:14740400 or 21:19268555, so no ancestral allele can be
read for the two records. The empty string says so, as the plugin itself
says it from the compressed form. The `C` is the ancestral state of 21:14439400,
a position neither record touches, and a one-base value cannot describe a
three-base span. The record cut short is a fact about the file; the value
reported for it should not change with the file's compression.

#### 2.6 Try it yourself

```
zcat <corpus>/ancestral/homo_sapiens_ancestor_21.slice.fa.gz | head -n 144395 > cut.fa
printf 'ANCESTOR_for_chromosome:GRCh38:21:1:46709983:1\t14439400\t48\t100\t101\n' > cut.fa.fai
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --dir_plugins <corpus> --plugin AncestralAllele,cut.fa
```

`<corpus>` is the corpus directory, `tests/golden/116/GRCh38-plugin-ancestral`,
and `record.vcf` holds the three records of 2.2.

vep-rs writes `AA=` on the two records past the end. Ensembl VEP 116.3 writes
`AA=C`. Compressing `cut.fa` with `bgzip` and indexing it, then naming
`cut.fa.gz`, makes the engines agree.

#### 2.7 How the comparator treats this class

The reproduction
`scripts/divergences/reproductions/ancestral_allele_uncompressed_fasta_clamp.toml`
holds the cut, the three records and both engines' rows. The class is
documented only: no committed run reaches the shape.
Ensembl's per-chromosome dumps run to the chromosome's full length, so an
intact dump has no span past its end; the committed FASTA is bgzipped.

## Where next

- [Plugins](../plugins.md): every Ensembl VEP plugin with its support status
  and how one is enabled.
- [The record a plugin reads for a `--minimal` split allele](../divergences/minimal-plugin-sees-the-rejoined-record.md):
  the sibling class, for a plugin that reads the rejoined record without
  caching.
- [Intended divergences from Ensembl VEP](../intended-divergences.md): the
  index of every documented difference between the engines.
