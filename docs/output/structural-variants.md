# How structural variants are written

This page describes how a structural variant is recognised, spanned, typed and written in every output format. Each rule follows Ensembl VEP 116 unless the text names a difference. The examples are VCF records written as `REF ALT` with their INFO keys where they matter.

## Routing

### What makes a record structural

A record is a structural variant when two conditions hold, as Ensembl VEP routes it. First, REF and its ALTs together are not a run of `ACGT`. Second, either INFO/SVTYPE is set to a value Perl reads as true (not empty, not `0`), or an ALT is written in a structural form: an angle or bracket token without `*`, or a dot-bounded breakend. So `N G` with `SVTYPE=DEL;END=...` is a deletion over END, `<*>` with `SVTYPE=DEL` is a deletion, and `N G` alone is a sequence variant.

Lowercase alleles of a sequence record are upper-cased once the anchor is trimmed (`a g` prints Allele `G`, JSON `allele_string` `A/G`), as Ensembl VEP's `validate_vf` does. The exception is a record the `--chr` list or the chromosome test refuses and `--dont_skip` keeps, which both engines annotate as written.

### The type of a symbolic allele

The record's type is read from its ALTs joined by `/` when that string opens with a VCF 4.4 form in upper case: `<DEL`, `<INS`, `<DUP`, `<INV`, `<CNV`, or `<CN` followed by an optional `=` and a digit (`<CN0`, `<CN=5`). A lowercase `<del>` is outside the forms, so SVTYPE names its type. Without SVTYPE it has no Sequence Ontology term, because the table's keys are upper case, and is an unsupported type.

The type is read through the regexes of Ensembl VEP's `get_SO_term`, so their edges hold. A copy number needs a digit: `<CN=X>` is the unsupported type `CN=X`, and with `SVTYPE=DEL` it is a deletion. A list's repeated parts collapse only when each part loses one bracket a side: `<<DEL>>,<DEL>` is the unsupported `DEL/DEL`, while a lone `<<DEL>>` loses every bracket and is a `deletion`, and `<<CN0>>` a `copy_number_variation`. A mobile element takes one colon before its name: `<INS:ME::ALU>` is a `mobile_element_insertion`, `<INS:ME:ALU>` an `Alu_insertion`. A trailing colon stays in the type: `<DEL:>` is unsupported, `<DEL:X>` a `deletion`.

### gVCF alleles

A gVCF variant site (`A G,<NON_REF>`) is its real alleles: `<NON_REF>` beside another allele names no variant and is dropped from the ALT list. A record whose ALTs are `<NON_REF>` alone, one or several, is a reference block with one intergenic row per ALT. Under a true INFO/SVTYPE every ALT, the block included, is an allele of the structural record, as Ensembl VEP reads it: `G,<NON_REF>` with `SVTYPE=DEL` is two deletion alleles.

### Records that name no allele

A trailing empty ALT field (`G,`, `<DEL>,`) names no allele, as Ensembl VEP's `split` drops it. A record whose ALT column is `.` is dropped with no warning whatever INFO/SVTYPE says, and kept as a non-variant record under `--allow_non_variant`. A record whose every ALT is `*` (`N *`, `*,*`) names no allele whatever INFO/SVTYPE says. It gets no row in the default and tab outputs, its VCF line is carried without CSQ, and its JSON object holds the record's `id`, coordinates, `allele_string` (`N/*`) and `"most_severe_consequence":"?"`, as Ensembl VEP writes the record without `SVTYPE`.

## Spans

A symbolic record starts at POS+1. It ends at POS+|SVLEN| when INFO carries SVLEN, else at INFO/END, else at its start. A tandem repeat or tandem duplication ends at POS plus the largest signed SVLEN, a trailing empty value being no value, so `SVLEN=-300,` ends it before its start and the record is skipped. SVLEN is read as Perl numifies it: `300abc` is 300, and `abc`, `.` or an empty value is 0. A repeated INFO key reads its last value.

`SVLEN=0`, a non-numeric SVLEN or an END at POS leaves the pair `POS+1..POS`, printed `chr:POS-POS+1`, for every class. Such a record gets the positional terms of a sequence insertion, with `feature_elongation` for a duplication or an insertion. An END more than one base before POS (`END=0` included), or one that is not a coordinate, skips the record at warning level, as Ensembl VEP 116 skips it.

Coordinates print as numbers. Ensembl VEP keeps a sign or a fraction it was given (`21:30425651-+30425899`, `21:30425651-30425950.7`); vep-rs prints 30425899 and reads a fractional SVLEN by its integer part. POS is read as the integer its text denotes, in any spelling that denotes one exactly (`3e7`, `30425949.0`, `30425949.`, a leading or trailing blank), on both paths, as Ensembl VEP's structural path numifies it. A POS that is not a number, or a negative one, skips the record at warning level (`Invalid POS`). A structural record whose first base POS + 1 passes 2^62 - 1 fails its line.

A trailing blank at the end of a data line is kept as part of the last INFO value, as Ensembl VEP keeps it: `END2=30425949 ` names the mate `N[21:30425949 [`, and the VCF output carries `END2=30425949 ;CSQ=`.

## Multi-allelic records

A multi-allelic record is one variant, as in VEP without `--minimal`: no per-allele minimisation, and every allele's rows in one record.

Under `--minimal` every bi-allelic sequence record is trimmed from both ends and marked `MINIMISED`, as `minimise_alleles` trims it (`Parser.pm` lines 964-1019). A multi-allelic record one of whose alleles trims is split as `split_variants` splits it (`InputBuffer.pm` lines 638-718). Each allele is annotated on its own trimmed span and numbered by its written position; a repeated ALT writes its own rows, and an allele equal to the reference once trimmed writes none.

The rows of a split record are rejoined as Ensembl VEP rejoins them (`OutputFactory.pm` lines 2272-2354). Location and the JSON `start` and `end` keep the record's span, `UPLOADED_ALLELE` prints the record's chopped string and `REF_ALLELE` the first allele's trimmed reference. The rows sort by allele: an allele whose trimmed string sorts before the letters of a transcript id (`-/A`, `A/G`) comes ahead of the first allele's rows. A record none of whose alleles trims is annotated as it stands. The position columns of a split allele follow the format ([Output formats](../output-formats.md#where-each-rule-lives)).

A record with more than one ALT and a symbolic allele among them has one class for the whole record, read from its ALTs joined by `/` as Ensembl VEP 116 reads them. Identical symbolic ALTs collapse to their type: `<INS>,<INS>` is an `insertion`, `<DEL>,<DEL>` a `deletion`. A DEL beside a DUP in either order, `<DEL>,<DUP:TANDEM>` included, is a `copy_number_variation`, and so is `<CN0>,<CN2>`. A sequence ALT beside a symbolic one (`G,<INV>` with `SVTYPE=DEL`) takes the type INFO/SVTYPE names.

Every ALT then gets its own row for each transcript, labelled with the class term as its Allele. So the rows of one transcript repeat per ALT and differ only in `ALLELE_NUM`, which numbers the parts of Ensembl VEP's allele string from 1. A breakpoint record's reference base is part 1, so a bare `<BND>` prints 2; a record with parsed breakends prints 0 for its own breakend and numbers its mates from 1. A list VEP has no type for (`<INS>,<DEL>`, `<INV>,<DUP>`) keeps one row per ALT under each ALT's own class and the unsupported-type treatment of the [next section](#skipped-records).

A `*` or `.` beside another ALT under a true INFO/SVTYPE is a part like any other: it writes its own row of the record's class, numbered at its written place, as Ensembl VEP writes it. So `C G,*` and `C G,.` with `SVTYPE=DEL` each write two deletion rows per transcript, `ALLELE_NUM` 1 and 2.

## Skipped records

Ensembl VEP 116 drops a record it marks skipped before annotation and writes nothing for it in any format; its `--dont_skip` restores it. Four kinds of record are marked so:

- a record whose end minus start exceeds `--max_sv_size` bases (default 10,000,000; `--max_sv_size -1` lifts the limit in both engines). A sequence record is measured too, on its minimised alleles, so a `<CNV:TR>` expanded with a FASTA is measured as the sequence record it becomes;
- a symbolic ALT its type table lacks (`<CPX>`, `<NON_REF>`);
- a multi-allelic symbolic list it has no type for (`<INS>,<DEL>`, `<INV>,<DUP>`);
- a `<DEL>` with neither `END` nor `SVLEN`, or whose start exceeds its end (`deletion looks incomplete`).

vep-rs annotates the first three kinds, a [documented divergence](../divergences/defect-c.md#5-a-record-ensembl-vep-drops-before-annotation), and drops the fourth as Ensembl VEP drops it:

- the default, tab and Parquet outputs list their rows: an oversize or unsupported record against every transcript it overlaps, a list with no type with one row per ALT under each ALT's own class;
- the VCF output carries the line of an oversize or unsupported record without CSQ;
- the JSON output writes an unsupported record as an `input`-only object and omits an oversize one.

A `<DEL>` with neither `END` nor `SVLEN`, or whose start exceeds its end, is warned `deletion looks incomplete` and written in no format, as Ensembl VEP writes it in none. `--dont_skip` keeps it on both engines, read as a deletion of the base after POS ([Defect C](../divergences/defect-c.md)). A `SVLEN` of `NaN` is the exception: vep-rs reads it as 0 and drops the record, where Ensembl VEP annotates it with the end `NaN` ([Defect K](../divergences/defect-k.md)). The golden corpora record the rows only vep-rs writes for the other three kinds under the manifest class `reference_skipped_record`.

A `<NON_REF>` gVCF reference-confidence block is not a variant allele. It is annotated against no transcript and appears in the default and tab formats as one `intergenic_variant` row spanning the block, with a warning per record.

## Breakends

### Points and spans

A breakend in bracket notation is a point record whose Location is the base after POS. A dot-notation breakend (`N.`, `.N`) is the same point, with or without INFO/SVTYPE. A record with a parsed breakend is one breakpoint whatever its other ALTs: `<BND>,<DEL>` with `SVTYPE=BND` writes every allele as a `chromosome_breakpoint` at the point, its END unread.

A `<BND>` spans POS+1 to POS+|SVLEN| when INFO carries SVLEN and is otherwise the single base at POS+1. INFO/END never ends a breakpoint's span, as in Ensembl VEP 116; a breakpoint whose END is `0` or empty is skipped, as Ensembl VEP skips it. A `<BND>` names a mate only through `INFO/CHR2` written with a value and `END2`, else `END`; POS2 is never read. Its alleles are evaluated over the span.

### Mates

The mate allele of a `<BND>` is written `N[CHR2:END2[` with both values as written (`N[chr21:00030425700[`), whatever the reference base; the record's own breakend takes the reference base (`A.`). An END2 that is not a run of digits (`-7`, `abc`) names no mate, and the record's allele string still carries the text (`N/N[21:-7[`). Whitespace around a bracket, the chromosome or the colon is read through, as Ensembl VEP's breakend regex reads it: `N[ 21 : 30425700 [` names the mate `21:30425700`. An INFO value ends at a second `=`, so `END2=30425700=5` is 30425700.

A dot-bounded ALT beside a bracket ALT (`.N,N[21:30425700[`) lends the own breakend its base and writes no row of its own. A `<BND>` that names no mate is one `chromosome_breakpoint` allele over its span or point, with the region terms and no `feature_truncation`.

### Rows on each side

A paired breakend has two alleles in the output. Its local end is printed as a single breakend (`.A`) and annotated against the transcripts at the record's own position. The bracket allele (`]21:27307377]A`) gets a row for a transcript only when the mate coordinate lies within the up- or downstream distance of that transcript. Ensembl VEP's `_close_to_feature` admits the mate within 5,000 bases whatever `--distance` says, and above that distance writes no row for the record at all ([the divergence page](../divergences/breakend-mate-beyond-5000.md)).

On the mate chromosome the bracket allele's rows describe the mate breakend at the mate coordinate. A mate inside the transcript earns `feature_truncation` with the region term; a mate outside it earns `upstream_gene_variant` or `downstream_gene_variant` with the distance. The `EXON` and `INTRON` fields are those of the mate coordinate. A transcript the local rows already name is not repeated on the mate side.

The record's own breakend is written `N.` or `.N` (the REF base, with the dot on the side of the bracket mate) whether or not `INFO/MATEID` is present. `chromosome_breakpoint` is the Allele of a record without a bracket mate, a single breakend (`.N`) or a bare `<BND>` with neither `CHR2` nor a span, as `class_SO_term` is Ensembl VEP's Allele for an allele without a breakend string.

### Intergenic breakends

A breakend record none of whose breakends lies within the up- or downstream distance of a transcript is intergenic, as Ensembl VEP 116 writes it: one `intergenic_variant` row per allele (`N.` and the bracket mate, or the single `chromosome_breakpoint`), Feature `-`, `IMPACT=MODIFIER`, in every format. When either breakend reaches a transcript the record gets no intergenic row, and an allele of it that reached nothing gets no row. `NMD_transcript_variant` on a breakend row follows the record's own position, so a mate inside an NMD transcript earns `feature_truncation` without it.
