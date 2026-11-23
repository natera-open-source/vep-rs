# How structural variants are written

vep-rs recognises, spans, types and writes structural variants by these rules in every output format, following Ensembl VEP 116 unless a difference is named.

The examples are VCF records written as `REF ALT` with their INFO keys where they matter.

## Routing

### What makes a record structural

A record is a structural variant when two conditions hold, as Ensembl VEP routes it. First, REF and its ALTs together are not a run of `ACGT`. Second, either INFO/SVTYPE is set to a value Perl reads as true (not empty, not `0`), or an ALT is written in a structural form: an angle or bracket token without `*`, or a dot-bounded breakend. So `N G` with `SVTYPE=DEL;END=...` is a deletion over END, `<*>` with `SVTYPE=DEL` is a deletion, and `N G` alone is a sequence variant.

Lowercase alleles of a sequence record are upper-cased once the anchor is trimmed (`a g` prints Allele `G`, JSON `allele_string` `A/G`), as Ensembl VEP's `validate_vf` does. The exception is a record the `--chr` list or the chromosome test refuses and `--dont_skip` keeps, which both engines annotate as written.

### The type of a symbolic allele

The record's type is read from its ALTs joined by `/` when that string opens with a VCF 4.4 form in upper case: `<DEL`, `<INS`, `<DUP`, `<INV`, `<CNV`, or `<CN` followed by an optional `=` and a digit (`<CN0`, `<CN=5`). A lowercase `<del>` is outside the forms, so SVTYPE names its type. Without SVTYPE it has no Sequence Ontology term, because the table's keys are upper case, and is an unsupported type.

The type is read through the regexes of Ensembl VEP's `get_SO_term`, so their edges hold:

- A copy number needs a digit: `<CN=X>` is the unsupported type `CN=X`, and with `SVTYPE=DEL` it is a deletion.
- A list's repeated parts collapse only when each part loses one bracket a side: `<<DEL>>,<DEL>` is the unsupported `DEL/DEL`. A lone `<<DEL>>` loses every bracket and is a `deletion`, and `<<CN0>>` a `copy_number_variation`.
- A mobile element takes one colon before its name: `<INS:ME::ALU>` is a `mobile_element_insertion`, `<INS:ME:ALU>` an `Alu_insertion`.
- A trailing colon stays in the type: `<DEL:>` is unsupported, `<DEL:X>` a `deletion`.

### gVCF alleles

A gVCF variant site (`A G,<NON_REF>`) is its real alleles: `<NON_REF>` beside another allele names no variant and is dropped from the ALT list. A record whose ALTs are `<NON_REF>` alone, one or several, is a reference block with one intergenic row per ALT. Under a true INFO/SVTYPE every ALT, the block included, is an allele of the structural record, as Ensembl VEP reads it. `G,<NON_REF>` with `SVTYPE=DEL` is two deletion alleles.

### Records that name no allele

A trailing empty ALT field (`G,`, `<DEL>,`) names no allele, as Ensembl VEP's `split` drops it. A record whose ALT column is `.` is dropped whatever INFO/SVTYPE says, and kept as a non-variant record under `--allow_non_variant`. A record whose every ALT is `*` (`N *`, `*,*`) names no allele whatever INFO/SVTYPE says. It gets no row in the default and tab outputs, its VCF line is carried without CSQ, and its JSON object holds the record's `id`, coordinates, `allele_string` (`N/*`) and `"most_severe_consequence":"?"`, as Ensembl VEP writes the record without `SVTYPE`.

| Record | vep-rs | Ensembl VEP |
| --- | --- | --- |
| ALT `.` under a true INFO/SVTYPE | no allele: dropped, or a non-variant record under `--allow_non_variant`; either way one warning names the record's line | a structural variant, annotated from INFO alone with or without `--allow_non_variant` ([Defect J](../divergences/defect-j.md)) |

## Spans

A symbolic record starts at POS+1. It ends at:

1. POS+|SVLEN| when INFO carries SVLEN;
2. else INFO/END;
3. else its start.

A tandem repeat or tandem duplication ends at POS plus the largest signed SVLEN, a trailing empty value being no value. So `SVLEN=-300,` ends it before its start and the record is skipped. SVLEN is read as Perl numifies it: `300abc` is 300, and `abc`, `.` or an empty value is 0. A repeated INFO key reads its last value.

`SVLEN=0`, a non-numeric SVLEN or an END at POS leaves the pair `POS+1..POS`, printed `chr:POS-POS+1`, for every class. Such a record gets the positional terms of a sequence insertion, with `feature_elongation` for a duplication or an insertion. An END more than one base before POS (`END=0` included), or one that is not a coordinate, skips the record at warning level, as Ensembl VEP 116 skips it.

Coordinates print as numbers. Ensembl VEP keeps a sign or a fraction it was given (`21:30425651-+30425899`, `21:30425651-30425950.7`); vep-rs prints 30425899 and reads a fractional SVLEN by its integer part. POS is read as the integer its text denotes, in any spelling that denotes one exactly (`3e7`, `30425949.0`, `30425949.`, a leading or trailing blank), on both paths, as Ensembl VEP's structural path numifies it. A POS that is not a number, or a negative one, skips the record at warning level (`Invalid POS`). A structural record whose first base POS + 1 passes 2^62 - 1 fails its line.

A trailing blank at the end of a data line is kept as part of the last INFO value, as Ensembl VEP keeps it. `END2=30425949 ` names the mate `N[21:30425949 [`, and the VCF output carries `END2=30425949 ;CSQ=`.

## Multi-allelic records

A multi-allelic record is one variant, as in Ensembl VEP without `--minimal`: no per-allele minimisation, and every allele's rows in one record. A sequence record whose first ALT is empty (`C ,G`) is the exception on the Ensembl VEP side. It is minimised as a one-allele deletion and its other allele is dropped, where vep-rs writes both alleles' rows ([Defect O, section 25](../divergences/defect-o.md#25-a-sequence-record-whose-first-alt-is-empty)).

### Trimming and splitting under `--minimal`

Under `--minimal` every bi-allelic sequence record is trimmed from both ends and marked `MINIMISED`, as `minimise_alleles` trims it (`Parser.pm` lines 964-1019). A multi-allelic record one of whose alleles trims is split as `split_variants` splits it (`InputBuffer.pm` lines 638-718). Each allele is annotated on its own trimmed span and numbered by its written position. A repeated ALT writes its own rows, and an allele equal to the reference once trimmed writes none.

The rows of a split record are rejoined as Ensembl VEP rejoins them (`OutputFactory.pm` lines 2281-2363). Location and the JSON `start` and `end` keep the record's span, `UPLOADED_ALLELE` prints the record's chopped string and `REF_ALLELE` the first allele's trimmed reference. The rows sort by allele: a later allele's rows are keyed by its trimmed string, `_` and the feature's key ([Row order](vep-and-tab.md#row-order)). An allele whose string sorts before the letters of a stable id (`-/A`, `A/G`) comes ahead of the first allele's transcript and regulatory rows. One opening with `-` comes ahead of the first allele's motif rows too, whose key is the motif's database id.

A record none of whose alleles trims is annotated as it stands.

### A symbolic allele among the ALTs

A record with more than one ALT and a symbolic allele among them has one class for the whole record, read from its ALTs joined by `/` as Ensembl VEP 116 reads them. Identical symbolic ALTs collapse to their type: `<INS>,<INS>` is an `insertion`, `<DEL>,<DEL>` a `deletion`. A DEL beside a DUP in either order, `<DEL>,<DUP:TANDEM>` included, is a `copy_number_variation`, and so is `<CN0>,<CN2>`. A sequence ALT beside a symbolic one (`G,<INV>` with `SVTYPE=DEL`) takes the type INFO/SVTYPE names.

Every ALT gets its own row for each transcript, its Allele the class term. One transcript's rows differ only in `ALLELE_NUM`, which numbers the parts of Ensembl VEP's allele string from 1, and, for a [stated `<CNV:TR>` length](../divergences/cnv-tr-symbolic-loss-named-a-gain.md), in Consequence and IMPACT. A breakpoint record's reference base is part 1, so a bare `<BND>` prints 2. A record with parsed breakends prints 0 for its own breakend and numbers its mates from 1. A list Ensembl VEP has no type for (`<INS>,<DEL>`, `<INV>,<DUP>`) keeps one row per ALT under each ALT's own class and the unsupported-type treatment of the [next section](#skipped-records).

A `*` or `.` beside another ALT under a true INFO/SVTYPE is a part like any other. It writes its own row of the record's class, numbered at its written place, as Ensembl VEP writes it. So `C G,*` and `C G,.` with `SVTYPE=DEL` each write two deletion rows per transcript, `ALLELE_NUM` 1 and 2.

## Skipped records

Ensembl VEP 116 drops a record it marks skipped before annotation and writes nothing for it in any format; its `--dont_skip` restores it. Four kinds of record are marked so:

- a record whose end minus start exceeds `--max_sv_size` bases (default 10,000,000; `--max_sv_size -1` lifts the limit in both engines). A sequence record is measured too, on its minimised alleles, so a `<CNV:TR>` expanded with a FASTA is measured as the sequence record it becomes;
- a symbolic ALT its type table lacks (`<CPX>`, `<NON_REF>`);
- a multi-allelic symbolic list it has no type for (`<INS>,<DEL>`, `<INV>,<DUP>`);
- a `<DEL>` with neither `END` nor `SVLEN`, or whose start exceeds its end (`deletion looks incomplete`).

vep-rs annotates the first three kinds, a [documented divergence](../divergences/defect-c.md#5-a-record-ensembl-vep-drops-before-annotation), and drops the fourth as Ensembl VEP drops it:

- the default, tab and Parquet outputs list their rows: an oversize or unsupported record against every transcript it overlaps (`<NON_REF>` excepted), a list with no type with one row per ALT under each ALT's own class;
- the VCF output carries the line of an oversize or unsupported record without CSQ;
- the JSON output writes an unsupported record as an `input`-only object and omits an oversize one.

A `<DEL>` with neither `END` nor `SVLEN`, or whose start exceeds its end, is warned `deletion looks incomplete` and written in no format, as Ensembl VEP writes it in none. `--dont_skip` keeps it on both engines, read as a deletion of the base after POS ([Defect C](../divergences/defect-c.md)). A `SVLEN` of `NaN` is the exception: vep-rs reads it as 0 and drops the record, where Ensembl VEP annotates it with the end `NaN` ([Defect K](../divergences/defect-k.md)). The golden corpora record the rows only vep-rs writes for the other three kinds under the manifest class `reference_skipped_record`.

A `<NON_REF>` gVCF reference-confidence block, not a variant allele, is annotated against no transcript. The default and tab formats write one `intergenic_variant` row per ALT spanning the block, with a warning per record.

## Breakends

### Points and spans

A breakend in bracket notation is a point record whose Location is the base after POS. A dot-notation breakend (`N.`, `.N`) is the same point, with or without INFO/SVTYPE. A record with a parsed breakend is one breakpoint whatever its other ALTs. `<BND>,<DEL>` with `SVTYPE=BND` writes every allele as a `chromosome_breakpoint` at the point, its END unread.

A `<BND>` spans POS+1 to POS+|SVLEN| when INFO carries SVLEN and is otherwise the single base at POS+1. INFO/END never ends a breakpoint's span, as in Ensembl VEP 116. A breakpoint whose END is `0` or empty is skipped, as Ensembl VEP skips it. A `<BND>` names a mate only through `INFO/CHR2` written with a value and `END2`, else `END`; POS2 is never read. Its alleles are evaluated over the span.

### Mates

The mate allele of a `<BND>` is written `N[CHR2:END2[` with both values as written (`N[chr21:00030425700[`), whatever the reference base. The record's own breakend takes the reference base (`A.`). An END2 that is not a run of digits (`-7`, `abc`) names no mate, and the record's allele string still carries the text (`N/N[21:-7[`).

Whitespace around a bracket, the chromosome or the colon is read through, as Ensembl VEP's breakend regex reads it. `N[ 21 : 30425700 [` names the mate `21:30425700`. An INFO value ends at a second `=`, so `END2=30425700=5` is 30425700.

A dot-bounded ALT beside a bracket ALT (`.N,N[21:30425700[`) lends the own breakend its base and writes no row of its own. A `<BND>` that names no mate is one `chromosome_breakpoint` allele over its span or point, with the region terms and no `feature_truncation`.

### Rows on each side

A paired breakend has two alleles in the output. Its local end is printed as a single breakend (`.A`) and annotated against the transcripts at the record's own position. The bracket allele (`]21:27307377]A`) gets a row for a transcript only when the mate coordinate lies within the up- or downstream distance of that transcript. Ensembl VEP's `_close_to_feature` admits the mate within 5,000 bases whatever `--distance` says, and above that distance writes no row for the record at all ([the divergence page](../divergences/breakend-mate-beyond-5000.md)).

On the mate chromosome the bracket allele's rows describe the mate breakend at the mate coordinate. A mate inside the transcript earns `feature_truncation` with the region term; a mate outside it earns `upstream_gene_variant` or `downstream_gene_variant` with the distance. The `EXON` and `INTRON` fields are those of the mate coordinate. A transcript the local rows already name is not repeated on the mate side.

The record's own breakend is written `N.` or `.N` (the REF base, with the dot on the side of the bracket mate) whether or not `INFO/MATEID` is present. `chromosome_breakpoint` is the Allele of a record without a bracket mate, a single breakend (`.N`) or a bare `<BND>` with neither `CHR2` nor a span, as `class_SO_term` is Ensembl VEP's Allele for an allele without a breakend string.

### Intergenic breakends

A breakend record none of whose breakends lies within the up- or downstream distance of a transcript is intergenic, as Ensembl VEP 116 writes it. Every format has one `intergenic_variant` row per allele (`N.` and the bracket mate, or the single `chromosome_breakpoint`), Feature `-`, `IMPACT=MODIFIER`. When either breakend reaches a transcript the record gets no intergenic row, and an allele of it that reached nothing gets no row. `NMD_transcript_variant` on a breakend row follows the record's own position, so a mate inside an NMD transcript earns `feature_truncation` without it.

## UTRAnnotator on a symbolic allele: the text `<DEL>` scored as a five-base deletion

Of the built-in plugins that annotate a structural variant's rows, UTRAnnotator
is the one whose Ensembl VEP module reads the allele's text as sequence. The
behaviour is a documented class the comparison tools set aside: this section
names its rows,
as [Intended divergences from Ensembl VEP](../intended-divergences.md) lists it.
Its [How to read the counts](../intended-divergences.md#how-to-read-the-counts)
defines the terms.

The class carries its tier under that page's rule: tier 1 is a shape the
chromosome 21 cuts carry with no flag, tier 2 one that needs a flag
or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts are
over the cuts' 15,646,481 records, the 152 golden corpora and their 511 flag
sets. A cut is the chromosome 21 part of a public population dataset, a golden
corpus is a test input kept in the repository with Ensembl VEP's output for it,
and a flag set is one set of options a corpus is run with.

**The defect in one sentence.** Ensembl VEP's UTRAnnotator plugin runs on
structural variants, whose allele string is the raw ALT text, so a `<DEL>`,
`<DUP>` or `<INV>` in a 5' UTR is spliced in as a five-base deletion, the text's
length, whatever span the record's `END` states.

**Why this is a defect.** The plugin's five analyses splice the record's alleles
into the 5' UTR sequence: the reference allele's length is removed and the
alternate allele's bases inserted. For a structural variant the reference allele
is the symbolic text, `<DEL>` or `<DUP>`, and the alternate allele is empty,
because a structural allele has no sequence method. The text's five characters
stand in for five bases: a 40-base deletion is scored as a 5-base one, five being
two modulo three, and a duplication or inversion as a deletion.

The plugin reads the record's `END` for one thing only: the share of a 5' UTR
exon the variant overlaps, which its `max_overlap` parameter tests. The uORF
effects the row then reports describe a change no allele of the record makes.

**What vep-rs does instead.** vep-rs annotates a structural allele from what
its record states:

| Allele | What vep-rs annotates |
| --- | --- |
| a symbolic deletion whose stated span lies inside one exon of the 5' UTR | the sequence deletion of that span |
| the loss alternate (`<DEL>`, `<CN0>`, `<CN=0>`) of a record with several alternates that writes `END`, its span inside one exon of the 5' UTR | the deletion of the record's span |
| an inversion that states one place inside one exon of the 5' UTR | its span replaced by its reverse complement; a span equal to its own reverse complement reports nothing |
| a plain insertion, `<INS>`, that states its length at one point | the frameshift that length decides |
| every other structural allele | `-` in both fields |

The plugin's analyses read the number of removed bases and the reference UTR,
never a letter of the deleted allele, so a deletion's row carries the two fields
the same deletion gets when its record spells the bases out.

A plain insertion, `<INS>`, states a length and no bases. Where the length is
no multiple of three and the point lies inside a uORF or an overlapping ORF, the
plugin's frameshift analysis fires whatever is inserted, and vep-rs writes
`5_prime_UTR_uORF_frameshift_variant`. Each effect carries the six keys that
read no inserted base: `Evidence`, `KozakContext`, `KozakStrength`,
`ref_StartDistanceToCDS`, `ref_type` and `ref_type_length`. It carries
`alt_type=uORF` where a stop codon past the point lands in the start codon's
frame once the insertion moves it, which holds for any inserted bases.

The insertion's row never carries `alt_type_length`, nor `alt_type` where no
such stop codon settles it. It does not say whether the inserted bases bring a
start codon, or split a start or a stop codon: an insertion written with its
bases reports such an effect in the frameshift's place.

The three `Existing_` counts of the reference UTR stand on every row, and both
engines agree on them.

`POS` and `SVLEN` state the deletion's span, or `POS` and `END` where the
record writes no `SVLEN`, and the reference holds every base of it. Ensembl VEP
annotates the same span when the record writes it with explicit alleles.

An insertion's record is read the same way: `SVLEN` bases from the base after
`POS`. A record that writes no `END`, `END` equal to `POS`, or `END` equal to
`POS` plus `SVLEN` states that one point. Any other `END` states an interval
beside the length, and vep-rs leaves `-`.

For the explicit spelling the two engines' rows differ only in ways three other
documented classes of this plugin cover. Ensembl VEP searches a
[window sized by the reference allele](../divergences/utrannotator-gained-start-window-comparison.md)
for a gained start codon. It reports
[another of the analyses](../divergences/utrannotator-consequence-precedence.md)
when several find an effect. It writes the annotation's `key=value` pairs in
[hash order](../divergences/utrannotator-annotation-pair-order.md).

vep-rs leaves `-` in both fields where the record states no change of the
table above. A
deletion or an inversion that ends outside the 5' UTR, or in another of its
exons, changes bases the spliced UTR does not hold. A record marked
`IMPRECISE`, or with `CIPOS` or `CIEND` other than `0,0`, admits other places,
at which another analysis can fire. In a record with several alternates vep-rs
takes the loss alternate's span from `END` alone, never from a list of
lengths, whose span is the first alternate's. Every other alternate keeps `-`.

A mobile element's insertion keeps `-`: its stated length need not count the
target-site duplication the element brings. An insertion under any other text,
such as `<INS:MT>`, keeps `-` too, because vep-rs takes a stated length as the
count of inserted bases for the plain `<INS>` alone. An insertion that states
no length gives nothing to analyse. One whose length is a multiple of three
shifts no frame, and one whose point lies in no reading frame has none to
shift.

vep-rs analyses no duplication. A plain duplication's record does not say where
the copy lies. A tandem duplication's record does, and its row keeps `-` all
the same.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116,
`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `UTRAnnotator.pm`, lines 80-82, `variant_feature_types`: structural variants
  are accepted; lines 151-153, `run`: `$ref` is the first allele of the allele
  string and `$alt` is `''` unless the allele object can `variation_feature_seq`,
  which only `Bio/EnsEMBL/Variation/VariationFeatureOverlapAllele.pm`, line 276,
  defines.
- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, line 572: a structural variant's allele
  string is the ALT text.
- `UTRAnnotator.pm`, lines 1158-1164 and 1297, `get_ref_coding` and
  `get_allele_exon_pos`: the text is the reference allele and its length the
  span; line 1197, `mut_utr_sequence`: that many bases are removed; lines 889 and
  923, `uFrameshift`: the length difference, five, decides the frame.

**What correcting it would change.** A structural deletion inside one exon of
the 5' UTR, a loss alternate among several included, would get the uORF effect
of the span its record states, as vep-rs writes it, and an inversion that of
its reverse complement. A plain insertion
would get at most the frameshift its length decides. Every other structural
allele would get none, and the five-base reading would vanish.

**The class it produces.** `utrannotator_structural_allele_text`: a
`5UTR_consequence` or `5UTR_annotation` (JSON `5utr_consequence`,
`5utr_annotation`) that differs on a row whose `Allele` is a structural class
word, not a sequence, as the rules below state. The class word is any term
a structural record's class is named by but a breakpoint's, `Alu_insertion` or
`copy_number_variation` as much as `deletion`. The ALT text's own length sets
the bases removed, twelve for `<INS:ME:ALU>`, and a record with two symbolic
alternates is read by its first text on both of its rows.

The format comparator's rule (`tools/vep-formats-compare`,
`structural_allele_text_row`) sets both fields aside on such a row where
Ensembl VEP's `5UTR_consequence` is one of the plugin's five terms and vep-rs's
row fills neither field. The plugin class table
(`scripts/concordance/plugin_class_table.tsv`) declares the class, so
`vep-plugin-pairs score` sets the same rows aside. Each row is judged against one
row of vep-rs under the same Location, `Allele` and Feature. Rows that share
all three are paired in file order, and both engines write a record's rows in
the order of its alternates, so each of the two rows of a record with two
alternates meets the row of its own alternate.

Where vep-rs fills the two fields, four more rules read the row's own columns
in the default and tab formats:

| Row's `Allele` | The row's own columns | vep-rs's term, where Ensembl VEP's is one of the five or both its fields are empty | The record's alternates, where the comparison is handed the record |
| --- | --- | --- | --- |
| a deletion's class word (`structural_deletion_span_row`) | the Location spans the distance between the row's two cDNA positions; the Consequence names `5_prime_UTR_variant` beside, at most, `feature_truncation` and `NMD_transcript_variant` | one a deletion can fire: any but the gained stop, and the frameshift only on a span that is no multiple of three | not read |
| `copy_number_variation` | the Location spans that distance; the Consequence names `5_prime_UTR_variant` beside nothing but `NMD_transcript_variant` | one a deletion can fire, as on a deletion's row | two or more, the row's own with the text `<DEL>`, `<CN0>` or `<CN=0>` |
| `inversion` | the Location spans that distance; the Consequence names `5_prime_UTR_variant` beside nothing but `NMD_transcript_variant` | the gained start, the lost start, the gained stop or the lost stop; the frameshift is refused | `<INV>` alone |
| `insertion` (`structural_insertion_length_row`) | the Location's span is no multiple of three; the Consequence names `5_prime_UTR_variant` | the frameshift, each effect holding the six keys read from the reference UTR, with or without `alt_type=uORF`, and no `alt_type_length` | not read |

The comparison is handed a row's record, the input line of `--input`, when it
also names a class that reads records, such as
`utrannotator_gained_start_window`. The row's own alternate is its place among
Ensembl VEP's rows of the record on its transcript, which must be as many as
the alternates. Without a record that states its alternates, the row's own
columns decide.

On such a row the class sets aside each of the two fields that differs from
Ensembl VEP's or that Ensembl VEP leaves `-`. A difference UTRAnnotator's
pair-order, gained-start window or precedence class explains is counted under
that class.

The class takes nothing else. A sequence allele's row, a row whose `Allele` is
a breakend's own string or a word outside those terms and a difference in the
three `Existing_` counts each stay a difference. So does a row vep-rs fills
outside those rules, and any filled row of the VCF or JSON format, which
names no Location.

A `chromosome_breakpoint` row stays a difference too. For a breakpoint record
Ensembl VEP puts the REF column before the allele string it builds from the ALT
text, unless that string opens or ends with a dot. The plugin takes the part
that stands first as its reference allele: the REF column, or the first ALT
text where the string was left as it stood. The word does not tell the two
cases apart, so the class takes neither. The
`GRCh38-plugin-utr-annotator` and `GRCh37-plugin-utr-annotator` corpora hold
substitutions, insertions and deletions alone, so no row of their four
comparison cases is set aside.

**How often the shape occurs.** Tier 2: 756 of the cuts' records (48
per million) carry a symbolic allele; the rows differ where it falls in a uORF
of a 5' UTR. 29 corpora carry the shape and 4 flag sets run the plugin.

### A 40-base `<DEL>` in a 5' UTR

**In short.** The record `21 10524588 . N <DEL> . . END=10524628` deletes 40
bases inside one exon of the 5' UTR of `ENST00000427445`, from three bases past
the `ATG` of a uORF. Ensembl VEP 116.3 removes five bases there, the length of
the text `<DEL>`, finds the uORF's frame shifted by two and writes
`5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;5UTR_annotation=alt_type_length=39:KozakContext=AGAATGT:KozakStrength=Moderate:ref_type_length=51:ref_type=uORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=207`.
vep-rs removes the 40 bases, finds the frame shifted by one and writes the same
term with
`5UTR_annotation=Evidence=False:KozakContext=AGAATGT:KozakStrength=Moderate:alt_type=uORF:alt_type_length=66:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51`.

Written with explicit alleles, the 41 reference bases from 10524588 as REF and
the first of them as ALT, the same deletion gets
`5UTR_annotation=alt_type_length=66:KozakContext=AGAATGT:KozakStrength=Moderate:ref_type_length=51:ref_type=uORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=207`
from Ensembl VEP: vep-rs's pairs in Ensembl VEP's own order. vep-rs writes the
frameshift on all nine rows of the symbolic record, each with the `Existing_*`
counts Ensembl VEP writes
(`Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=6` on that
transcript), and the row's `Allele` is `deletion` on both engines. A real
five-base deletion at the same start, `GGACCG>G`, gives Ensembl VEP the
`5UTR_annotation` it writes for the symbolic record on every transcript.

A constructed record, `21 10524588 . N <INS> . . END=10524588;SVLEN=4`, inserts
four unstated bases at the same place. Ensembl VEP 116.3 writes for it, on all
nine rows, the five plugin fields it writes for the `<DEL>` record. vep-rs
writes the frameshift with
`5UTR_annotation=Evidence=False:KozakContext=AGAATGT:KozakStrength=Moderate:alt_type=uORF:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51`
on `ENST00000427445`. `alt_type=uORF` stands because a stop codon further on
ends the shifted frame whatever is inserted.

### An inverted start codon, and a loss beside a gain

Two constructed records, one for each:

| Constructed record | What it states | Ensembl VEP 116.3, under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0` | vep-rs |
| --- | --- | --- | --- |
| `21 14210904 utr_inv_start_codon N <INV> . . SVTYPE=INV;END=14210907;SVLEN=3`, on `ENST00000400211`, one of six transcripts whose rows differ the same way | The `ATG` 48 bases from the cap of this reverse-strand transcript becomes `CAT`: `CCTCAATGAGTT` becomes `CCTCACATAGTT`. The `TAG` that opens there, 60 bases before the coding sequence, ends in frame the ORF that begins 87 bases before it. `21 14210905 CAT ATG`, the same change with explicit alleles, gets vep-rs's pairs from Ensembl VEP. | `5UTR_annotation=alt_type_length=48:KozakContext=GTGATGA:KozakStrength=Moderate:ref_type_length=1215:ref_type=InFrame_oORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=87;5UTR_consequence=5_prime_UTR_uORF_frameshift_variant` | `5UTR_annotation=Evidence=False:KozakContext=GTGATGA:KozakStrength=Moderate:newSTOPDistanceToCDS=60:ref_StartDistanceToCDS=87:ref_type=InFrame_oORF;5UTR_consequence=5_prime_UTR_uORF_stop_codon_gain_variant` |
| `21 10524588 utr_loss_beside_gain N <CN0>,<CN2> . . SVTYPE=CNV;END=10524590`, on `ENST00000427445`, one of nine such transcripts: the loss alternate's row | `END` states a two-base span, where Ensembl VEP removes five bases, the length of the text `<CN0>`. `21 10524588 GGA G`, the same deletion with explicit alleles, gets vep-rs's pairs from Ensembl VEP. | `5UTR_annotation=alt_type_length=39:KozakContext=AGAATGT:KozakStrength=Moderate:ref_type_length=51:ref_type=uORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=207;5UTR_consequence=5_prime_UTR_uORF_frameshift_variant` | `5UTR_annotation=Evidence=False:KozakContext=AGAATGT:KozakStrength=Moderate:alt_type=uORF:alt_type_length=42:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51;5UTR_consequence=5_prime_UTR_uORF_frameshift_variant` |
| the same record, the other alternate's row | no sequence | the same two fields | `5UTR_annotation=-;5UTR_consequence=-` |

## Where next

- [Output formats](../output-formats.md#where-each-rule-lives): a split allele's position columns follow the format.
