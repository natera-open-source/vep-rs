# `--regulatory`: the warning, the SPDI chromosome name, the GFF index and five defects

Three differences under `--regulatory`, each outside the flags the measured datasets run with, then five defects of Ensembl VEP 116.3's regulatory file reader and breakend lookup that vep-rs does not reproduce; nothing is set aside in the concordance figures. Ensembl VEP writes no message against vep-rs's warning on a cache without regulatory data, the cache's chromosome name against the FASTA's in a regulatory-only record's `--spdi` rows, and its library's message against an `ERROR` for a `--regulatory_gff` file without an index. It misses a promoter past its seek margin, prints GFF3 escapes, drops a row sharing an identifier, and repeats a breakend's row or adds one for another chromosome's mate; vep-rs writes each row once, decoded, and none for that mate.

## A cache without regulatory data

Ensembl VEP builds its regulatory source only when the cache's `info.txt` says `regulatory 1` (`CacheDir.pm` lines 168 to 176). A `--regulatory` run on a cache dumped without regulatory features therefore writes transcript and intergenic rows alone, with no message: the flag is accepted and silently does nothing.

vep-rs writes the same rows and one line to the warnings file, `--regulatory: this cache holds no regulatory data; no RegulatoryFeature or MotifFeature row is written`, so a run that expected regulatory rows and got none can tell the cache from the flag. The line is written once per run, whether `info.json` says `regulatory: false` or the `regulatory/` directory is missing, and `--everything` reaches it the same way, since the set switches `--regulatory` on (`Config.pm` line 362).

## The chromosome name of `--spdi` on a regulatory-only record

When the FASTA names a chromosome unlike the cache, a regulatory-only record's rows under `--spdi` name it as the cache does in Ensembl VEP and as the FASTA does here, `21:30234967:A:G` against `chr21:30234967:A:G`. The regulatory feature hands the variant its slice there (`AnnotationType/RegFeat.pm` line 126), a slice vep-rs never has.

A record that overlaps a transcript takes the transcript's slice first on both engines, so its rows agree; the difference is confined to a record whose only features are regulatory, and to a FASTA whose names differ from the cache's.

## The index of a `--regulatory_gff` file

Both engines read a `--regulatory_gff` file through tabix, so each file is bgzip-compressed and tabix-indexed. The regulatory-feature file Ensembl distributes comes without an index, which `tabix -p gff <file>` writes. Ensembl VEP opens the file through its tabix library (`AnnotationSource/File/RegFeat.pm` lines 204 to 209), which ends the run with the library's own message when the index is missing. vep-rs ends the run before any file is read, with `ERROR: [cli] --regulatory_gff <key>=<file> has no tabix index`, naming the key and the file as the other refusals of the flag do.

## The comparator

A golden corpus whose cache holds no regulatory tree compares the two engines with the mask class `regulatory_feature_rows`, which the test runner and `tools/vep-formats-compare --masks` both apply. The reference's `RegulatoryFeature` and `MotifFeature` rows are set aside, and the `most_severe_consequence` of each JSON record that lost such a row is recomputed from the rows that remain. `--keep-only regulatory_feature_rows` compares those rows alone. Of the committed corpora, the class is carried by `everything` flag sets alone, each of a corpus whose species directory is not `homo_sapiens` and whose converted cache holds no regulatory tree, whatever the cache of its reference run held.

## Where the regulatory file reader and the breakend lookup differ from Ensembl VEP

Five behaviours of Ensembl VEP 116.3 are not reproduced. Each is a documented class, which no comparison tool applies: no committed dataset carries a row of it. [Intended divergences from Ensembl VEP](../intended-divergences.md) lists it; its [How to read the counts](../intended-divergences.md#how-to-read-the-counts) defines the terms. Every row below is in the default output format unless the section names another, against the Ensembl 116 GRCh38 cache of chromosome 21 and, for the `--regulatory_gff` examples, a regulatory GFF3 of one or two lines written for the example.

Each class carries its tier under that page's rule: tier 1 is a shape the chromosome 21 cuts carry with no flag, tier 2 one that needs a flag or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts are over the cuts' 15,646,481 records, the 152 golden corpora and their 511 flag sets. The regulatory and motif rows themselves, their terms and their fields are otherwise the same on both engines ([Regulatory and motif features](../regulatory.md)). The two field values that differ have their own pages: [the motif score of a multi-base allele on a reverse-strand motif](motif-score-reverse-strand.md) and [`--cell_type` names with a comma](cell-type-comma-names.md), which also carries [the empty `--cell_type` value](cell-type-comma-names.md#an-empty---cell_type-value).

### 1. A promoter whose extension reaches past the seek margin of the next 1 Mb region

**The defect in one sentence.** Under `--regulatory_gff --extended_promoters`, Ensembl VEP reads each 1 Mb region of the file with a fixed 10,000-base margin, so a promoter whose core lies further than that beyond the region's edge is missed for a variant in its extension.

**Why this is a defect.** The margin's own comment says it exists "so a promoter whose core lies in a neighbouring region but whose extended bounds reach into this one is still found". A fixed 10,000 bases honours that only for an extension of at most that length, and the `extended_start` and `extended_end` attributes carry no such bound. The rows a variant gets then depend on the other variants in its buffer: one of them inside the core's region makes that region's read find the promoter for every buffered variant.

**What vep-rs does instead.** vep-rs reads each chromosome of the file whole into one index, so a widened promoter is found for every variant its extension covers, whatever the distance from its core and whatever else the buffer holds.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/RegFeat.pm`, lines 110-113: the margin constant and its comment. Lines 301-311, `_get_regfeats_by_coords`: the tabix seek of a 1 Mb region padded by the margin on each side, and the read loop ending at the first record that starts beyond the padded end.
- The same file, lines 360-365, `_record_to_regfeat`: the widening to `extended_start` and `extended_end`, applied only to the records the padded seek returned. `Bio/EnsEMBL/VEP/AnnotationSource.pm`, lines 190-200, and `Bio/EnsEMBL/VEP/AnnotationType/RegFeat.pm`, lines 165-167, `up_down_size`: the regions read are those of the buffer's variants, with no flank for a regulatory source.

**What correcting it would change.** Every variant in a promoter's extension would get the promoter's row, as vep-rs's rows do and as Ensembl VEP's do when another variant of the buffer lies in the core's region.

**The class it produces.** `promoter_extension_beyond_seek_margin`: a `RegulatoryFeature` row of `BIOTYPE` `promoter` (the JSON `biotype`) that vep-rs alone writes under `--regulatory_gff --extended_promoters`, for a promoter whose core lies more than 10,000 bases beyond the edge of the variant's 1 Mb region. The outputs carry no promoter coordinates, so a manifest's candidate-only tuples apply it. No row on the committed datasets: the regulatory file of the GRCh38 regulatory GFF dataset extends each of its 492 promoters by at most 1,000 bases, and none crosses a 1 Mb boundary when widened. The class is documented only.

**How often the shape occurs.** Tier 3: the shape is a promoter whose extension crosses a 1 Mb region boundary by more than the margin. None of the 27,530 records of the committed regulatory and promoter GFF files does, and 1 flag set runs `--extended_promoters`.

#### 1.1 A promoter at `21:21020000-21020500` extended back to `21:20990000`

**In short.** A regulatory file holding `21 Ensembl promoter 21020000 21020500 . . . ID=FAR_BIN;extended_start=20990000;extended_end=21020500` beside the GRCh38 chromosome 21 cache, and the record `21 20995000 vc_far_bin_ext A C`, under `--regulatory_gff file=<that file> --extended_promoters`. The variant lies inside the widened promoter, so vep-rs writes `FAR_BIN RegulatoryFeature regulatory_region_variant ... BIOTYPE=promoter` for it; Ensembl VEP writes its two transcript rows alone. Its region, `21:20000001-21000000`, is read up to `21:21010000`, and the core starts at `21:21020000`. With a second record at `21:21500000` in the input, Ensembl VEP writes the promoter row for the first: that record's region is read, and the promoter found there is tested against every buffered variant.

A core at `21:21000100`, inside the margin, is found by both engines; one at `21:21010001`, one base past it, by vep-rs alone.

### 2. A `--regulatory_gff` attribute value holding a GFF3 escape

**The defect in one sentence.** Ensembl VEP prints a `--regulatory_gff` record's identifier, matrix id and factor names as the file writes them, so a value holding a GFF3 escape (`%3B` for `;`, `%2C` for `,`) prints the escape in place of the character.

**Why this is a defect.** GFF3 reserves `;`, `=`, `,` and `&` in the ninth column and has a value carry them as `%XX` escapes, which the reader decodes. The parser Ensembl VEP reads the file with decodes them in `get_attribute_by_name`, `get_source` and `get_type`, through a `decode_string` whose own description reads "Return the decoded string: %2C will be ,". Its `get_attributes`, the accessor the regulatory source calls, splits the raw column on `;` and `=` and stores each value as written. That split also ends a value at a second `=`, so `ID=a=b` names the feature `a`, in both engines.

**What vep-rs does instead.** vep-rs decodes every `%XX` escape, in either letter case, in every attribute value it prints, so the identifier, `MOTIF_NAME` and `TRANSCRIPTION_FACTORS` carry the characters the file encodes. The `transcription_factor` list is split on the commas written, as Ensembl VEP splits it, each name decoded after. A promoter's `extended_start` or `extended_end` holding an escape is read, under `--extended_promoters`, as the number it encodes. Ensembl VEP compares the written text as the number 0, so an escaped `extended_start` widens the promoter to the chromosome's start and an escaped `extended_end` leaves its end as written.

**Where it is in Ensembl VEP's source** (`ensembl-io` at release/116, `ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/IO/Parser/GFF3.pm`, lines 190-198, `get_attributes`: `split(';', ...)` then `split('=', $attr)`, no decoding. Lines 174-183, `get_attribute_by_name`: `decode_string($value)`; lines 151-164, `get_source` and `get_type`: decoded likewise. `Bio/EnsEMBL/IO/Parser/GXF.pm`, lines 303-307, `decode_string`: `uri_unescape`.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/RegFeat.pm`, line 353: `$parser->get_attributes`. Lines 372, 438 and 441: the `ID` or `Name`, `transcription_factor` and `binding_matrix_id` values read from that hash into the feature.

**What correcting it would change.** A `Feature`, `MOTIF_NAME` or `TRANSCRIPTION_FACTORS` value would print the characters the file encodes, as vep-rs's does; the files Ensembl distributes hold no escape, so their rows would not change.

**The class it produces.** `regulatory_gff_attribute_escapes`: a `Feature`, `MOTIF_NAME` or `TRANSCRIPTION_FACTORS` value (the JSON `regulatory_feature_id`, `motif_feature_id`, `motif_name` and `transcription_factors`) of a `--regulatory_gff` row that Ensembl VEP prints with a `%XX` escape. Its text decoded, and in the VCF format re-spelt as the `CSQ` writer spells `,` and `;`, must equal vep-rs's. It also covers the span of a promoter whose extended bound holds an escape. No row on the committed datasets: the attribute columns of the GRCh38 regulatory GFF dataset's three files hold no `%` (0 of 5,201, 13,577 and 3,493 records). The class is documented only.

**How often the shape occurs.** Tier 3: the shape is an attribute value holding a `%XX` escape, which vep-rs decodes. None of the 27,530 records of the committed regulatory and promoter GFF files has one, and 6 flag sets run `--regulatory_gff`.

#### 2.1 An identifier written `ESC%3BAPE%2C1`

**In short.** A regulatory file holding `21 Ensembl CTCF_binding_site 5014322 5014343 . - . ID=ESC%3BAPE%2C1` beside the GRCh38 chromosome 21 cache, and the record `chr21 5014335 escaped_id G C` under `--regulatory_gff file=<that file>`. Ensembl VEP's tab row names the feature `ESC%3BAPE%2C1`, and its VCF entry the same; vep-rs's tab row names it `ESC;APE,1`, the string the escapes encode, and its VCF entry `ESC%3BAPE&1`, the spelling both engines give a `;` and a `,` inside a `CSQ` value. Every other column of the row agrees. A line with `ID=a=b` names the feature `a` in both engines, the text between the first and the second `=`.

### 3. One identifier in both the regulatory file and the EMAR file of `--regulatory_gff`

**The defect in one sentence.** When a record of the `file=` source and a record of the `emars=` source share an identifier, Ensembl VEP keeps a variant's row for the EMAR record alone and drops the regulatory record's row with no warning.

**Why this is a defect.** Ensembl VEP builds the EMAR source as its own instance because, in its own words, "EMARs are a separate file of RegulatoryFeatures with their own identifier namespace" and "these are genuinely distinct features". Yet every regulatory source adds its overlaps to one hash per variant keyed by the feature's stable id, so the source read last overwrites the other's entry. The identifier check that guards that hash describes the loss ("a collision silently discards one of the annotations - last write wins") but keeps its seen set per source, so a collision across the two files draws no warning.

**What vep-rs does instead.** vep-rs writes one row per record that overlaps the variant, so both features print, each with its own `BIOTYPE`. Two records of one file sharing an identifier are one feature on both engines, the later record's, with Ensembl VEP's warning ([Regulatory features from GFF files](../regulatory.md#regulatory-features-from-gff-files)).

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSourceAdaptor.pm`, lines 272-275 and 284-300, `get_all_regulatory_gff`: the comment quoted above, then the `file` and `emars` instances built in that order.
- `Bio/EnsEMBL/Variation/VariationFeature.pm`, line 918, `add_RegulatoryFeatureVariation`: `$self->{regulatory_feature_variations}->{$rfv->regulatory_feature_stable_id} = $rfv`, one hash for every source; line 710: the rows written from its keys.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/RegFeat.pm`, lines 463-467 and 478-498, `_check_id_collision`: the description quoted above, and a seen set in `$self->{_seen_ids}`, the source's own, so the warning "identifiers must be unique or annotations will be lost" is drawn within one file alone.

**What correcting it would change.** A variant inside two features of one identifier from the two files would get both rows, as vep-rs writes them.

**The class it produces.** `regulatory_duplicate_identifier_rows`: a `RegulatoryFeature` row (a JSON `regulatory_feature_consequences` entry) vep-rs alone writes under `--regulatory_gff` with both `file=` and `emars=`, whose record, allele and `Feature` equal those of an Ensembl VEP row of another `BIOTYPE`. No row on the committed datasets: the identifiers of the GRCh38 regulatory GFF dataset's three files are unique within each file (5,201, 13,577 and 3,493) and disjoint across them (`ENSR21_`, `ENSR0` and `ENSM0` prefixes). The class is documented only.

**How often the shape occurs.** Tier 3: the shape is an EMAR record carrying an identifier of the regulatory file's shape. None of the 3,493 records of the committed EMAR file does, and 6 flag sets run `--regulatory_gff`.

#### 3.1 An enhancer and an EMAR both named `ENSR21_85ZF5`

**In short.** The regulatory file's `21 Ensembl enhancer 5011320 5011550 . . . ID=ENSR21_85ZF5` and an EMAR file's `21 Ensembl EMAR 5011300 5011600 . . . ID=ENSR21_85ZF5`, with the record `21 5011400 v_collision A C` under `--regulatory_gff file=<regulatory>,emars=<EMAR>` over the GRCh38 chromosome 21 cache. Ensembl VEP writes one `RegulatoryFeature` row, `ENSR21_85ZF5 ... BIOTYPE=EMAR`, and no warning; vep-rs writes two, `BIOTYPE=enhancer` then `BIOTYPE=EMAR`. The variant lies inside both records; the EMAR row alone survives because the EMAR source is built after the regulatory one and writes the shared key last. With the EMAR file alone, both engines write the one EMAR row.

### 4. A breakend whose two ends lie in one regulatory feature

**The defect in one sentence.** Under `--regulatory`, Ensembl VEP writes each regulatory or motif row of a paired breakend once per listing of the variant at a position inside the feature, four times when its own position and its mate both lie in the feature.

**Why this is a defect.** The buffer's lookup tree holds a breakend variant once at its own span and once at each mate position, so a query over a feature's span returns the variant once per listing it overlaps. Each returned copy is then pushed once for its own span and once for each breakend position inside the span. A transcript overlap built from those listings is stored under one key per transcript, so the repeats collapse to one row. A regulatory overlap is pushed onto an array, one per listing, so the same row prints as many times as the variant was listed.

**What vep-rs does instead.** vep-rs writes a feature's row once per allele of the record, however many of the record's positions lie inside the feature.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/InputBuffer.pm`, lines 380-390, `interval_tree`: the variant inserted at its own span and at each breakend position. Lines 320-343, `get_overlapping_vfs`: each fetched copy pushed for its own span and again for each breakend position inside the query.
- `Bio/EnsEMBL/VEP/AnnotationType/RegFeat.pm`, lines 125-136, `annotate_InputBuffer`: one `StructuralVariationOverlap` per listing, pushed onto `$vf->{regulation_structural_variations}->{$type}`. `Bio/EnsEMBL/Variation/StructuralVariationFeature.pm`, lines 495-502: that array written as it stands.
- `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, lines 144-151, and `StructuralVariationFeature.pm`, lines 708-728, `add_TranscriptStructuralVariation`: the transcript overlap of the same listings stored under one key per transcript, so it prints once.

**What correcting it would change.** A breakend's regulatory and motif rows would print once per allele and feature, as vep-rs's do and as Ensembl VEP's own transcript rows do.

**The class it produces.** `regulatory_breakend_listing_duplicates`: under `--regulatory`, a `RegulatoryFeature` or `MotifFeature` row of a breakend record (a JSON `regulatory_feature_consequences` or `motif_feature_consequences` entry) that Ensembl VEP writes more than once, byte for byte. The golden harness and the format comparator compare every copy as an entry of its own, so they report the repeats as entries vep-rs lacks. No row on the committed datasets: the two breakends of the GRCh38 regulatory dataset's structural input (`bnd_mate_far`, `bnd_mate_near`) have their mates outside every regulatory feature. The class is documented only.

**How often the shape occurs.** Tier 3: the shape is a bracket breakend whose mate is on its own chromosome; the regulatory build decides whether one feature spans both ends. None of the cuts' 15,646,481 records is one, 10 corpora carry such breakends, and 60 flag sets run `--regulatory`, six of them through `--regulatory_gff`, whose option set switches it on.

#### 4.1 A breakend at `21:5040000` whose mate at `21:5040250` lies in the same feature

**In short.** The record `21 5039999 bnd_both_in_rf N N[21:5040250[ . . SVTYPE=BND` over the GRCh38 chromosome 21 cache under `--regulatory`: both ends lie in `ENSR21_9422C` (`21:5039957-5040269`). Ensembl VEP writes eight `RegulatoryFeature` rows, the `N.` and `N[21:5040250[` rows alternating four times, each `regulatory_region_variant`; vep-rs writes the two rows once. The variant sits in the lookup tree at both positions, and each fetched copy is pushed for its own span and again for its mate, four listings. The two intergenic rows agree. A breakend whose mate alone lies in the feature is listed there once, so both engines write the feature's row once per allele; the bracket allele's term then differs as [Defect E](defect-e.md) describes, Ensembl VEP writing `feature_truncation` and vep-rs `regulatory_region_variant`.

### 5. A breakend whose mate on another chromosome selects a feature of its own chromosome

**The defect in one sentence.** Under `--regulatory`, Ensembl VEP lists a paired breakend for a regulatory or motif feature of the record's own chromosome whenever the number of its mate's position falls in the feature's span, whatever chromosome the mate names, writing the own allele's `intergenic_variant` row for that feature.

**Why this is a defect.** The regulatory annotation asks the buffer for the variants overlapping a feature's span and names no chromosome, and the buffer's chromosome test passes every variant when the query names none. A breakend sits in the buffer's tree at its mate's coordinate as well as its own, so a mate at `1:5040100` is matched by a chromosome 21 feature spanning position 5040100. The transcript annotation received the chromosome-aware lookup in release 116.1; the regulatory one did not. The row carries `intergenic_variant` on a `RegulatoryFeature` row, the fall-through term, and comes and goes with the other records of the buffer, since the feature is loaded only when the buffer's span on that chromosome reaches it.

**What vep-rs does instead.** vep-rs looks a mate's position up against the features of the mate's chromosome alone, so a mate on another chromosome adds no row for the record's own chromosome, and a record's rows do not depend on its neighbours in the buffer.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationType/RegFeat.pm`, line 125, `annotate_InputBuffer`: `get_overlapping_vfs($rf->{start}, $rf->{end})`, two arguments and no chromosome; `AnnotationType/Transcript.pm`, line 108, passes the transcript's chromosome.
- `Bio/EnsEMBL/VEP/InputBuffer.pm`, lines 339-341: a variant is listed when a mate's position lies in the query and the chromosome test passes; line 351: that test returns true when the query names no chromosome. Lines 387-389: the tree holds the variant at each mate's position.
- `Bio/EnsEMBL/Variation/StructuralVariationOverlap.pm`, lines 90-102 and 144-158: the record's own allele is kept for a feature of its own chromosome within 5,000 bases, the other-chromosome mate is not; `BaseVariationFeatureOverlapAllele.pm`, line 285: the default term.
- `Bio/EnsEMBL/VEP/AnnotationSource.pm`, lines 251-253 and 306-309: the mate's region is fetched on the mate's chromosome, and a feature of the record's chromosome is kept only when the buffer's span there reaches it.

**What correcting it would change.** The regulatory lookup would pass the feature's chromosome as the transcript lookup does, so a mate on another chromosome would select no feature of the record's own chromosome, as in vep-rs.

**The class it produces.** `regulatory_other_chromosome_mate_row`: under `--regulatory`, a `RegulatoryFeature` or `MotifFeature` row of a paired breakend's own allele (a JSON `regulatory_feature_consequences` or `motif_feature_consequences` entry) that Ensembl VEP alone writes, with `intergenic_variant`. The feature is on the record's chromosome, and the position number of a mate on another chromosome falls in it. The class is documented only. No committed comparison carries such a row: none of the 6 corpora whose breakends name a mate on another chromosome runs `--regulatory`.

**How often the shape occurs.** Tier 3: the shape is a bracket breakend whose mate is on another chromosome; the regulatory build and the buffer's span decide whether the mate's number lands in a feature within 5,000 bases of the record. None of the cuts' 15,646,481 records is one, 6 corpora carry such breakends, and 60 flag sets run `--regulatory`, six of them through `--regulatory_gff`, whose option set switches it on.

#### 5.1 A breakend at `21:5038000` whose mate is `1:5040100`

**In short.** The records `21 5037999 bnd_x_mate_number N N[1:5040100[ . . SVTYPE=BND` and `21 5040499 anchor_del N <DEL> . . END=5040600;SVTYPE=DEL` in one buffer over the GRCh38 chromosome 21 cache under `--regulatory`. The feature `ENSR21_9422C` spans `21:5039957-5040269`: the mate's number falls in it, and the deletion brings the buffer's chromosome 21 span up to it. Ensembl VEP writes `bnd_x_mate_number 21:5038000 N. - ENSR21_9422C RegulatoryFeature intergenic_variant` before the record's two intergenic rows; vep-rs writes the two intergenic rows alone. The breakend alone, or with its mate number outside every feature (`N[1:5038500[`), writes no feature row in either engine; with the mate on chromosome 21 both write the row, the bracket allele's term differing as [Defect E](defect-e.md) describes.
