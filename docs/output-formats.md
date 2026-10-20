# Output Formats

vep-rs supports five output formats, selected via `--output_format` or the shorthand flags `--vcf`, `--json`, `--tab`:

- `vep` (default): VEP's default tab-delimited format with an `Extra` column; `ensembl` is Ensembl VEP's other name for it
- `tab`: VEP's `--tab` format, every Extra key as its own named column
- `vcf`: the input VCF with a CSQ INFO field appended
- `json`: line-delimited JSON, one object per input record
- `parquet`: a Hive-partitioned Parquet directory produced via a `duckdb` subprocess (see [Parquet Format](#parquet-format) below); vep-rs's own, with no Ensembl VEP counterpart

The `--output_format` value is matched in any case (`VCF`, `Json`), as Ensembl VEP lower-cases it; a value outside the six is an argument error listing them. `--vcf`, `--json` and `--tab` refuse one another with Ensembl VEP's message (`Can't use --json and --vcf together`), and `--vcf` with `--most_severe` or `--summary` is refused the same way ([CLI reference](cli-reference.md#incompatible-flags-configpm615-632)).

### Compressed output

`--compress_output gzip` writes the output as one gzip member and `--compress_output bgzip` as a BGZF stream (gzip members of at most 64 KiB ending in the empty end-of-file member, the form `bgzip`, `tabix` and `samtools` read), to the output file under the name given, or to standard output under `-o STDOUT`. Ensembl VEP pipes its output through the external `gzip -c` or `bgzip -c`; vep-rs encodes in process, so the decompressed text is the same output and the compressed bytes are the encoder's. The flag applies to the default, `tab`, `vcf` and `json` formats and is refused with `parquet`, whose files carry their own compression. The golden corpus `GRCh38-writers` runs both encoders in every text format and compares the decompressed output with Ensembl VEP's.

### Headers

Every text format opens with a header block, and `--no_headers` suppresses it whole, the column line included, in the default, `tab` and `vcf` formats (`json` has none); the Parquet TSV intermediate keeps its column line, which DuckDB reads as the schema, whatever the flag. The default and `tab` headers are Ensembl VEP's `##` block: the version line, the time, the cache path, the API line, one `## <source> version <value>` line per source version of the cache, the column descriptions, the command line and the column line. Two kinds of line Ensembl VEP prints there have no counterpart: the Perl API component versions (`## ensembl version 116.0d85231`, `## ensembl-variation version ...`), since no Perl API runs here, and the `## Connected to <database>` line of a database run, since vep-rs reads a local cache only; the `## <short name> : [PATH]/<file>` lines of `--custom` sources are not written, since custom annotation sources are not supported. The VCF header's `##INFO` definition and plugin lines follow the same rule, below.

## Parity with Ensembl VEP

With the default flag set and `--tab`, `--vcf` or `--json`, the default, `tab`, `vcf` and `json` outputs match what Ensembl VEP 116.2 writes field for field: every column of every row and, for JSON, every key of every object compared as text (VEP writes JSON keys in Perl hash order, so the JSON is not byte-identical and key order is not compared). Header lines that carry a time, a cache path or the command line are compared by shape, and VEP's Perl API component version lines have no counterpart. Every format also names the engine that wrote it, which Ensembl VEP's output does not: the line `## vep-rs <version>` second in the default and tab headers, the meta line `##vep-rs=<version>` after `##VEP=` in the VCF header, and the key `vep_rs_version` in every JSON object (see [Engine identity](#engine-identity)); the comparators set aside exactly one such line per header, and the key on every object, and compare everything else. The contract is enforced by the golden corpora under `tests/golden/<release>/<assembly>/` (`crates/vep-cli/tests/golden.rs` and `format_parity.rs`): every consequence-set combination observed on the benchmark datasets, annotated with the corpus's own pruned cache and compared with VEP's committed output in all four formats, plus the Parquet round trip. Consequence terms themselves are governed by the concordance measurement ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)). The corpus manifest records the (Location, Allele, Feature) keys whose terms differ between the engines (117 of the 17,556 GRCh37 rows, 4 of the 2,818 GRCh37-hgvs rows, 323 of the 20,623 GRCh38 rows and 2 of the 977 GRCh38-release116 rows, none of them `unexplained_residual`) and the rows only one engine emits, each under its divergence class; for a recorded key the parity tests check only that the row's consequence set is one of the recorded sets, compare none of its other columns, and skip the record's `most_severe_consequence` in JSON. Outside the tested flag set there are known differences: `--regulatory` adds its (always empty) fields to the headers and columns and warns that regulatory annotation is not implemented, `--output_format vcf` switches on the same fields as `--vcf` (VEP ties them to the flag alone).

The golden-corpus flag set is `--offline --json_cache <dir> --species homo_sapiens --assembly <assembly> --buffer_size 5000 --no_stats --quiet` (or, for a corpus whose cache sits in Ensembl VEP's `cache_root/homo_sapiens/<release>_<assembly>` layout, `--offline --cache --dir_cache <root> --cache_version <release>` in place of `--json_cache <dir>`) plus one format flag and the corpus manifest's own flags; under it the default, `--tab`, `--vcf` and `--json` data rows match Ensembl VEP 116.2 apart from the consequence-term divergences the manifests document, and header lines that carry a time, a path, the command line or Perl API component versions match in shape. `--vcf` also switches on `--symbol`, `--biotype` and `--numbers`, as in VEP, and `--json` switches on `--no_escape`, as VEP's `rest` option set does; `--output_format json` switches it on here as well, where VEP leaves that output escaped ([intended divergences](intended-divergences.md#unescaped-json-under---output_format-json)). `--everything` switches on exactly Ensembl VEP's set ([CLI reference](cli-reference.md#option-sets-configpm351-593)), so the header keys agree; the regulatory columns it implies through `--regulatory` are not produced, and no RegulatoryFeature or MotifFeature row is written.

The rules the formatters follow are VEP's own, cited to their Perl sources in `crates/vep-io/src/output/fields.rs` (the single field table shared by every format), `crates/vep-effects/src/display.rs` (positions, codons and peptides) and the formatters beside them. The ones a reader of the output most often needs:

- **Row order.** Rows of one input record are transcript-major (sorted by transcript stable id) and allele-minor (ALT order); intergenic rows come last. A record's rows are contiguous in the default and tab formats and become one line (VCF) or one object (JSON).
- **Uploaded_variation.** The first `;`-element of a VCF record's ID (an Ensembl-format ID column is taken whole), unless it is `.`, empty or `0` (the values VEP reads as no ID) or the `chr_start_alleles` name VEP generates for an empty ID column, so `rs2;.` names the record `rs2` while `.;rs2`, `;rs2`, `0` and `0;rs2` do not name it. Such a record is named from its own input line, as Ensembl VEP 116 names it: `CHROM_POS_REF/ALT1/ALT2` with POS and every allele exactly as the VCF line spells them, case included (`21_10795228_C/G/T`; `21_30425650_a/g` for the lowercase `21 30425650 . a g`; `21_47612512_C/CTG` for the insertion `21 47612512 . C CTG`, whose variant starts at 47612513; `21_33867341_C/<CN2>` for the symbolic record `21 33867341 . C <CN2>`, whose variant starts at 33867342), and `chr_start_alleles` for an Ensembl-format line, start the lower of its two coordinates, except that a line whose alleles are written in lower case is named at the start and from the alleles its record carries after Ensembl VEP's own minimisation (`21_30348410_ta/-` for `21 30348409 30348411 cta/c`), the name Ensembl VEP builds before it upper-cases them. An HGVS or SPDI input is named by its notation as written, whitespace removed (`21:g.30282169C>A`, `NC_000021.9:10538884:A:C`), and a region input by its line as written (`21:25585733-25585733:1/G`), in both.
- **Positions.** `cDNA_position`, `CDS_position` and `Protein_position` are `start-end` ranges, a single value when both ends coincide, and `?` for an end that falls in an intron or outside the transcript, and for the CDS and protein columns also in a UTR (`?-335`, `493-?`). An insertion prints the two flanking bases (`1161-1162`). `cDNA_position` appears only for exon-overlapping rows; the CDS and protein columns, `Amino_acids` and `Codons` only for coding rows.
- **Amino_acids and Codons.** `ref/alt` peptides, collapsed to one when unchanged (`P` for a synonymous change); a partial codon appends `X` (`M/X`, `D/DX`); an allele whose codon window is empty prints `-` (`-/S` for an in-frame insertion, `M/-` for a whole-codon deletion). Codons are the reference codon window in lower case with the allele's own bases in upper case, on both sides (`gTGCTGGgc/ggc` for a deletion, `-/TCT` for an insertion). Structural variants carry positions but no codons.
- **DISTANCE** is the shortest gap between the variant and the transcript, the smallest of the four gaps between the record's start or end and the transcript's start or end (so the pair `POS+1..POS` of a zero-length structural variant measures from whichever flank is nearer, and a `<CNV:TR>` contraction read under the cap measures from the run it is located on, not from the units it removes), printed on up- and downstream rows only. **FLAGS** lists every `cds_` transcript attribute (`cds_start_NF,cds_end_NF`). **OverlapBP** and **OverlapPC** appear on every structural variant row that overlaps its feature (not on breakends) in the default Extra column and in JSON (`bp_overlap`, `percentage_overlap`); `--overlaps` adds them as columns to the `tab` and `vcf` outputs, as in Ensembl VEP, and to Parquet. They are reported on transcript rows only: Ensembl VEP also reports them on regulatory-feature and motif-feature rows, and regulatory annotation is not implemented.
- **Structural routing.** A record is a structural variant when REF and its ALTs together are not a run of `ACGT` and either INFO/SVTYPE is set (a value Perl reads as true: not empty, not `0`) or an ALT is written in a structural form (an angle or bracket token without `*`, a dot-bounded breakend), as Ensembl VEP routes it: `N G` with `SVTYPE=DEL;END=...` is a deletion over END, `<*>` with `SVTYPE=DEL` a deletion, `N G` alone a sequence variant. The record's type is read from its ALTs joined by `/` when the string opens with a VCF 4.4 form in upper case: `<DEL`, `<INS`, `<DUP`, `<INV`, `<CNV`, or `<CN` followed by an optional `=` and a digit (`<CN0`, `<CN=5`); a lowercase `<del>` is outside the forms, so SVTYPE names its type, and without SVTYPE it has no Sequence Ontology term (the table's keys are upper case) and is an unsupported type. A gVCF variant site (`A G,<NON_REF>`) is its real alleles: `<NON_REF>` beside another allele names no variant and is dropped from the ALT list; a record whose ALTs are `<NON_REF>` alone, one or several, is a reference block with one intergenic row per ALT; under a true INFO/SVTYPE every ALT, the block included, is an allele of the structural record, as Ensembl VEP reads it (`G,<NON_REF>` with `SVTYPE=DEL` is two deletion alleles). The type is read through `get_SO_term`'s own regexes, so their edges hold: a copy number needs a digit (`<CN=X>` is the unsupported type `CN=X`; with `SVTYPE=DEL` it is a deletion), a list's repeated parts collapse only when each part loses one bracket a side (`<<DEL>>,<DEL>` is the unsupported `DEL/DEL`, while a lone `<<DEL>>` loses every bracket and is a `deletion` and `<<CN0>>` a `copy_number_variation`), a mobile element takes one colon before its name (`<INS:ME::ALU>` is a `mobile_element_insertion`, `<INS:ME:ALU>` an `Alu_insertion`), and a trailing colon stays in the type (`<DEL:>` is unsupported, `<DEL:X>` a `deletion`). A trailing empty ALT field (`G,`, `<DEL>,`) names no allele, as Ensembl VEP's `split` drops it. A record whose ALT column is `.` is skipped at warning level whatever INFO/SVTYPE says. A record whose every ALT is `*` (`N *`, `*,*`) names no allele whatever INFO/SVTYPE says: no row in the default and tab outputs, the VCF line carried without CSQ, and a JSON object of the record's `id`, coordinates, `allele_string` (`N/*`) and `"most_severe_consequence":"?"`, as Ensembl VEP writes the record without `SVTYPE`. Lowercase alleles of a sequence record are upper-cased once the anchor is trimmed (`a g` prints Allele `G`, JSON `allele_string` `A/G`), as Ensembl VEP's `validate_vf` does.
- **Structural spans.** A symbolic record starts at POS+1 and ends at POS+|SVLEN|, else at INFO/END, else at its start; a tandem repeat or tandem duplication reads the largest SVLEN. SVLEN is read as Perl numifies it (`300abc` is 300, `abc`, `.` or an empty value is 0), and a repeated INFO key reads its last value. `SVLEN=0`, a non-numeric SVLEN or an END at POS leaves the pair `POS+1..POS`, printed `chr:POS-POS+1`, for every class: the positional terms of a sequence insertion, with `feature_elongation` for a duplication or an insertion. An END more than one base before POS (`END=0` included) or not a coordinate skips the record at warning level, as Ensembl VEP 116 skips it. Coordinates print as numbers: Ensembl VEP keeps a sign or a fraction it was given (`21:30425651-+30425899`, `21:30425651-30425950.7`), where vep-rs prints 30425899 and reads a fractional SVLEN by its integer part. POS is read as the integer its text denotes, in any spelling that denotes one exactly (`3e7`, `30425949.0`, `30425949.`, a leading or trailing blank), on both paths, as Ensembl VEP's structural path numifies it; a POS that is not a number, or a negative one, skips the record at warning level (`Invalid POS`), and a structural record whose first base POS + 1 passes 2^62 - 1 fails its line. A trailing blank at the end of a data line is kept as part of the last INFO value, as Ensembl VEP keeps it (`END2=30425949 ` names the mate `N[21:30425949 [`, and the VCF output carries `END2=30425949 ;CSQ=`).
- **Multi-allelic records** are one variant, as in VEP without `--minimal`: no per-allele minimisation, every allele's rows in one record.
- **Multi-allelic structural variants.** A record with more than one ALT and a symbolic allele among them has one class for the whole record, read from its ALTs joined by `/` as Ensembl VEP 116 reads them: identical symbolic ALTs collapse to their type (`<INS>,<INS>` is an `insertion`, `<DEL>,<DEL>` a `deletion`), a DEL beside a DUP in either order, `<DEL>,<DUP:TANDEM>` included, is a `copy_number_variation`, `<CN0>,<CN2>` is a `copy_number_variation`, and a sequence ALT beside a symbolic one (`G,<INV>` with `SVTYPE=DEL`) takes the type `INFO/SVTYPE` names. Every ALT then gets its own row for each transcript, labelled with the class term as its Allele, so the rows of one transcript repeat per ALT and differ only in `ALLELE_NUM`, which numbers the parts of Ensembl VEP's allele string from 1: a breakpoint record's reference base is part 1, so a bare `<BND>` prints 2, and a record with parsed breakends prints 0 for its own breakend and numbers its mates from 1. A list VEP has no type for (`<INS>,<DEL>`, `<INV>,<DUP>`) keeps one row per ALT under each ALT's own class and the unsupported-type treatment below.
- **Records Ensembl VEP skips.** Ensembl VEP 116 drops a record it marks skipped before annotation and writes nothing for it in any format (its `--dont_skip` restores it): a structural variant spanning more than `--max_sv_size` bases (default 10,000,000; `--max_sv_size -1` lifts the limit in both engines), a symbolic ALT its type table lacks (`<CPX>`, `<NON_REF>`), a multi-allelic symbolic list it has no type for (`<INS>,<DEL>`, `<INV>,<DUP>`), and a `<DEL>` with neither `END` nor `SVLEN`. vep-rs annotates every one of them, a documented divergence ([intended-divergences.md](intended-divergences.md#5-a-record-ensembl-vep-drops-before-annotation)): the default, tab and Parquet outputs list their rows (an oversize or unsupported record against every transcript it overlaps; a list with no type with one row per ALT under each ALT's own class); the VCF output carries the line of an oversize or unsupported record without CSQ; the JSON output writes an unsupported record as an `input`-only object and omits an oversize one; a `<DEL>` with neither `END` nor `SVLEN` is read as a deletion of the base after POS and written in every format. The golden corpora record the rows only vep-rs writes for such a record under the manifest class `reference_skipped_record`. A `<NON_REF>` gVCF reference-confidence block is not a variant allele: it is annotated against no transcript and appears in the default and tab formats as one `intergenic_variant` row spanning the block, with a warning per record.
- **Breakends.** A breakend in bracket notation is a point record whose Location is the base after POS; a dot-notation breakend (`N.`, `.N`) is the same point with or without INFO/SVTYPE. A record with a parsed breakend is one breakpoint whatever its other ALTs: `<BND>,<DEL>` with `SVTYPE=BND` writes every allele as a `chromosome_breakpoint` at the point, its END unread. A `<BND>` spans POS+1 to POS+|SVLEN| when INFO carries SVLEN and is otherwise the single base at POS+1 (INFO/END never ends a breakpoint's span, as in Ensembl VEP 116; a breakpoint whose END is `0` or empty is skipped, as Ensembl VEP skips it); it names a mate only through `INFO/CHR2` written with a value and `END2`, else `END` (POS2 is never read), and its alleles are evaluated over the span. The mate allele is written `N[CHR2:END2[` with both values as written (`N[chr21:00030425700[`), whatever the reference base, and the record's own breakend takes the reference base (`A.`); an END2 that is not a run of digits (`-7`, `abc`) names no mate, and the record's allele string still carries the text (`N/N[21:-7[`); whitespace around a bracket, the chromosome or the colon is read through, as Ensembl VEP's breakend regex reads it (`N[ 21 : 30425700 [` names the mate `21:30425700`), and an INFO value ends at a second `=` (`END2=30425700=5` is 30425700). A dot-bounded ALT beside a bracket ALT (`.N,N[21:30425700[`) lends the own breakend its base and writes no row of its own. A `<BND>` that names no mate is one `chromosome_breakpoint` allele over its span or point, with the region terms and no `feature_truncation`. A paired breakend has two alleles in the output: its local end, printed as a single breakend (`.A`) and annotated against the transcripts at the record's own position, and the bracket allele (`]21:27307377]A`), which gets a row for a transcript only when the mate coordinate lies within the up- or downstream distance of that transcript, as VEP's `_close_to_feature` admits it. On the mate chromosome the bracket allele's rows describe the mate breakend at the mate coordinate: `feature_truncation` with the region term when the mate lies inside the transcript, `upstream_gene_variant` or `downstream_gene_variant` with the distance when it lies outside, and the `EXON` and `INTRON` of that coordinate. A transcript the local rows already name is not repeated on the mate side. The record's own breakend is written `N.` or `.N` (REF base, dot on the side of the bracket mate) whether or not `INFO/MATEID` is present; `chromosome_breakpoint` is the Allele of a record without a bracket mate, a single breakend (`.N`) or a bare `<BND>` with neither `CHR2` nor a span, as `class_SO_term` is Ensembl VEP's Allele for an allele without a breakend string. A breakend record none of whose breakends lies within the up- or downstream distance of a transcript is intergenic, as Ensembl VEP 116 writes it: one `intergenic_variant` row per allele (`N.` and the bracket mate, or the single `chromosome_breakpoint`), Feature `-`, `IMPACT=MODIFIER`, in every format; when either breakend reaches a transcript the record gets no intergenic row, and an allele of it that reached nothing gets no row. `NMD_transcript_variant` on a breakend row follows the record's own position, so a mate inside an NMD transcript earns `feature_truncation` without it.
- **ENSP.** Under `--protein` every transcript row of a transcript with a translation carries its ENSP, structural-variant rows (`transcript_ablation`, `feature_truncation`, breakends) included.
- **The row filters** `--no_intergenic` and `--coding_only` run before the pick family, where Ensembl VEP's `get_all_VariationFeatureOverlapAlleles` applies them: `--no_intergenic` removes the intergenic row of every allele that overlaps no transcript; `--coding_only` keeps the rows whose allele affects the transcript's coding sequence, as Ensembl VEP's `within_cds` judges it: a CDS position, or a coding consequence term where the position columns print `-` (a deletion whose ends both lie outside the coding sequence keeps its row by the start, stop or coding-sequence term its inner span earns), or a variant covering a translated transcript whole; every other transcript row and the intergenic rows go. A record left with no row is written as an unannotated record: no row in the default and tab formats, its VCF line without a CSQ entry, and in JSON the record's identity, coordinates and allele string with `"most_severe_consequence":"?"` and no consequence list. A pick flag then ranks the rows the filters leave, so `--coding_only --pick` keeps one coding row per record and none for a record without one.
- **The pick family** acts on a whole input record, as Ensembl VEP's `filter_VariationFeatureOverlapAlleles` ranks every overlap allele of a variation feature together: `--pick` keeps one row for the record whatever its allele count, `--per_gene` one row per gene over the record, `--flag_pick` marks one row of the record with `PICK=1`; `--pick_allele`, `--pick_allele_gene` and `--flag_pick_allele` group the rows by the Allele column first, and `--flag_pick_allele_gene` by allele then gene. A structural record collapses `--pick_allele` to `--pick`, `--pick_allele_gene` to `--per_gene`, `--flag_pick_allele` to `--flag_pick` and `--flag_pick_allele_gene` to one flagged row per gene over the whole record. The ranking follows `--pick_order` (default `mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq`): each criterion in turn, a lower value better, the rows tied on it compared on the next; `rank` is the lowest consequence rank of the row, `length` the translateable-sequence length of a translated transcript and the spliced length otherwise, longer better; `ensembl` and `refseq` read the source cache of a merged cache and tie on every other cache; a shorter list falls through to the first row in Ensembl VEP's order, transcript-major by stable id and allele-minor, so a full tie keeps the lexicographically first stable id. A record overlapping no transcript has one intergenic row per allele, which `--pick` keeps for the first allele alone and the flagging variants mark like a transcript row. The `PICK` column exists under `--flag_pick`, `--flag_pick_allele` and `--flag_pick_allele_gene` alone. Ensembl VEP emits the gene buckets of `--per_gene`, the allele buckets of `--pick_allele` and both of `--pick_allele_gene` in Perl hash order, which differs from run to run; vep-rs keeps the row order above, and the golden runs of those three sets compare the rows without their order (`"unordered": true` in the corpus manifest).
- **`--most_severe` and `--summary`** squash the rows left after the pick family into one row per record: `Uploaded_variation`, `Location` and `Consequence` (the one lowest-ranked term, or every term once in rank order), `-` in every other fixed column, `Existing_variation` included even under `--check_existing`, and an Extra or tab cells holding the variant-level keys alone (`VARIANT_CLASS`, `MINIMISED`, `NEAREST`, `SV`) or `-`; no IMPACT, STRAND, ALLELE_NUM or plugin key. `--most_severe --per_gene` summarises the per-gene picks, which can differ from the unpicked set when a gene's picked transcript is not its most severe. The JSON object files the row under `intergenic_consequences` as `{"consequence_terms": [...]}` and keeps `most_severe_consequence`, the variant-level keys and `colocated_variants` at the top level. `--vcf` is refused with either flag, as Ensembl VEP refuses it, and each refuses the per-feature flags of Ensembl VEP's incompatibility list (`biotype`, `symbol`, `canonical`, `pick`, `flag_pick`, ...) with its message.
- **NEAREST.** Under `--nearest transcript|gene|symbol` every row of a record, intergenic rows and the `--summary` row included, carries the identifier(s) of the protein-coding transcription start site(s) nearest the record, found as Ensembl VEP's `TranscriptTree::nearest` finds them: the start site of a protein-coding transcript is its start on the forward strand and its end on the reverse; the two candidates are the highest site at or below the record's end and the lowest site at or above its start plus two (a site exactly one base above the start is the interval containing the start in Ensembl VEP's tree and is never a candidate); each is measured to the nearer end of the record and both are reported when they tie, the site at or below the end first; a site on an end of the record is at distance zero. The identifier is the transcript stable id, the gene stable id, or the gene symbol falling back to the gene stable id; values repeat once. Where several transcripts start at the winning coordinate, Ensembl VEP reports the one its interval tree holds first, which depends on the order its cache directory was listed in, and vep-rs reports every transcript at that coordinate in stable-id order; `gene` and `symbol` agree with Ensembl VEP whenever those transcripts belong to one gene. The site table is built from the run's own transcript set, so `--gencode_basic` or `--gencode_primary` narrows it; Ensembl VEP writes its table (`transcript_gene_tss.txt`) under the filters of the first `--nearest` run on a cache directory and reuses that file afterwards. A chromosome with no protein-coding transcript prints `NEAREST=` with no value in the default format, an empty tab cell and `"nearest": []`, as Ensembl VEP does. The JSON key is the top-level `nearest` array; the VCF CSQ field joins values with `&`.
- **HGVSg.** Under `--hgvsg` (which needs a reference FASTA) every row of a sequence variant, intergenic rows included, carries the HGVS genomic notation of its allele, as Ensembl VEP's `hgvs_genomic` writes it: the chromosome as the input names it, `g.` numbering, or `m.` on the mitochondrial chromosome, and the reference bases read from the FASTA; an insertion or deletion is described at its most 3' position along the chromosome and an insertion repeating the bases beside it as a duplication (`21:g.46925315_46925323dup` for `21 46925297 C CCCCCCTGGG`); a same-length change that is the reverse complement is an inversion, and any other change is clipped to the bases that differ (`21:g.34775851T>C` for `21 34775850 ATG ACG`). Each allele of a multi-allelic record is trimmed against the reference on its own to place and class it; an insertion or deletion is then written from the trimmed bases, any other allele with the bases the record's working allele string holds for it, at the trimmed span, as Ensembl VEP writes it (`21:g.34775850delinsGT` for the `CGT` allele of `21 34775849 CAT C,CGT`, whose working allele string is `AT/-/GT`). An allele outside `ACGT`, an allele equal to the reference bases and a structural variant carry no value. The JSON key is `hgvsg`.
- **SPDI.** Under `--spdi` every row of a sequence variant carries its allele in SPDI notation, as Ensembl VEP's `spdi_genomic` writes it from the record as the parser leaves it: the chromosome, the 0-based position, the deleted bases and the inserted bases joined by `:` (`21:10538884:A:C` for `21 10538885 A C`, `21:10578406::A` for an insertion, `21:20000325:T:` for a one-base deletion, `21:30000000:25:` for a deletion of more than 20 bases, written as their count). A multi-allelic record's alleles are the anchor-trimmed pair at the record's span, so `21 14144642 AT CT,A` prints `21:14144641:AT:CT` and `21:14144641:AT:A`; an allele given on the reverse strand is complemented; nothing is shifted. The chromosome is named as the cache names it on every row overlapping a transcript, and on an intergenic row as the reference FASTA names it (Ensembl VEP needs a FASTA for those rows and stops without one; vep-rs then prints the cache's name). Ensembl VEP keeps the record's notations by allele and looks them up by the row's `Allele`, so under `--shift_3prime 1 --shift_genomic 1` a row whose inserted bases the shift rotated carries no `SPDI`, unless the rotation lands on another allele of the record, in which case the row carries that allele's notation (`21 13620283 G GTCA,GCAT`: the `TCA` allele moved 10 bases prints `CAT` as its Allele and `21:13620283::CAT` as its `SPDI`). A structural variant and an allele outside `ACGTN` carry none. The JSON key is `spdi`.
- **GA4GH_VRS.** Under `--ga4gh_vrs` (refused with `--vcf`) every row that has an SPDI, whether or not `--spdi` prints it, also carries the GA4GH Variation Representation Specification `Allele` Ensembl VEP derives from that SPDI: a `SequenceLocation` on `refseq:<NC_ accession>` with a `SimpleInterval` from the SPDI's position to that position plus the deleted length, and a `SequenceState` holding the inserted bases; no digest. Ensembl VEP builds it only from an SPDI named by an `NC_` accession, which it has with a database, or offline on an intergenic row alone when the FASTA record carries that name; vep-rs takes the accession from the chromosome synonyms table (`--synonyms` or the cache's `chr_synonyms.txt`, the first `NC_` synonym of the cache's chromosome name) and writes the object Ensembl VEP writes with a database, or nothing for a chromosome with no such synonym. The text formats carry the object JSON-encoded on one line, `start` a string and `end` a number as Ensembl VEP's encoder writes them, with the keys in a fixed order where Ensembl VEP's follow its hash order; the JSON format nests the object under `ga4gh_vrs`, both bounds numbers.

## Engine identity

The Ensembl-shaped header lines stay as Ensembl VEP writes them, so a parser written for VEP output reads vep-rs output unchanged; each format adds one item naming the vep-rs version (the release, `0.5.0` in the examples below, read from the crate version `vep --version` prints):

| format | where | form |
| --- | --- | --- |
| default, `--tab` | the second header line, after the `## ENSEMBL VARIANT EFFECT PREDICTOR` banner | `## vep-rs 0.5.0` |
| `--vcf` | the meta line after `##VEP=`, in the `##key=value` form VCF 4.3 section 1.4 defines (a bare `## vep-rs` line is not valid VCF) | `##vep-rs=0.5.0` |
| `--json` | a key of every object, an unannotated record's input-only object included | `"vep_rs_version": "0.5.0"` |
| Parquet | the footer's key-value metadata | `vep_rs_version` |

`--no_headers` drops the default, tab and VCF lines with the rest of the header. Re-annotating a VCF that vep-rs wrote replaces its `##vep-rs=` line along with the earlier `##VEP` lines, as a re-annotation by Ensembl VEP replaces its own. The golden-corpus comparators (`crates/vep-cli/tests/common/mod.rs`) set aside exactly one identity line or key of the version under test and report a missing, repeated or differently versioned one as a header difference.

## VEP Default Format

One line per (allele, transcript) consequence, plus one intergenic line per allele that overlaps no feature.

### Header

```
## ENSEMBL VARIANT EFFECT PREDICTOR v116.2
## vep-rs 0.5.0
## Output produced at [TIME]
## Using cache in [PATH]/116
## Using API version 116, DB version ?
## 1000genomes version phase3
## COSMIC version 98
## ClinVar version 202306
## ...one line per source version recorded in the cache's info.json...
## Column descriptions:
## Uploaded_variation : Identifier of uploaded variant
## Location : Location of variant in standard coordinate format (chr:start or chr:start-end)
## ...one line per column...
## Extra column keys:
## IMPACT : Subjective impact classification of consequence type
## DISTANCE : Shortest distance from variant to transcript
## STRAND : Strand of the feature (1/-1)
## FLAGS : Transcript quality flags
## VEP command-line: vep --assembly GRCh37 --force_overwrite --input_file [PATH]/input.vcf --json_cache [PATH]/116 --offline --output_file [PATH]/out.txt
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
```

The source-version lines come from the cache's `info.json` (`source_versions`, written by `scripts/data/storable_to_json.pl` and `vep-cache-builder`, or by `scripts/golden/prune_json_cache.py --perl-info` from a VEP cache's `info.txt`). The `Using cache in` line prints the cache directory as given; `--mask_header_cache_path` prints it as `[PATH]/<leaf>`, its last component alone, here, in the `--tab` header, in the VCF `cache=` token and in the Parquet footer. The `Extra column keys` block lists every key the active flags can produce, then one `## <KEY> : <description>` line per field of each loaded built-in plugin, with the plugin's own text, as Ensembl VEP lists its plugins' header info. `--no_headers` suppresses the whole block, column line included. `--fields` changes nothing here: the default format keeps its thirteen columns and its `Extra` column whatever the list names, as in Ensembl VEP.

### Columns

| #   | Column             | Description                                                                                                                                                                                              |
| --- | ------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Uploaded_variation | Input ID, else the record's own line as `CHROM_POS_REF/ALT1/ALT2` (POS and alleles as written)                                                                                                           |
| 2   | Location           | `chr:pos` or `chr:start-end` (`chr:end-start` for an insertion)                                                                                                                                          |
| 3   | Allele             | The ALT allele of this row; for a structural variant the SO term of its class (`deletion`, `duplication`, `Alu_insertion`), except that a breakend keeps its literal ALT (`N.`, `A[21:27348902[`, `N[21:00030425700[` for a mate named by `CHR2`; `chromosome_breakpoint` for a record without a parsed mate, a single breakend or a bare `<BND>`), every ALT of a multi-allelic structural variant prints the record's class term, and an unsupported `<CPX>` prints `CPX`. A sequence allele whose reference base is ambiguous (`N G`) prints `G` with `coding_sequence_variant` and no amino acids, since no peptide is read from an ambiguous reference |
| 4   | Gene               | Ensembl gene stable ID                                                                                                                                                                                   |
| 5   | Feature            | Ensembl transcript stable ID                                                                                                                                                                             |
| 6   | Feature_type       | `Transcript`, or `-` on an intergenic row; VEP's `RegulatoryFeature` and `MotifFeature` rows are not produced                                                                                            |
| 7   | Consequence        | SO consequence term(s), comma-separated                                                                                                                                                                  |
| 8   | cDNA_position      | Position or range in the cDNA; `?` for an unmappable end                                                                                                                                                 |
| 9   | CDS_position       | Position or range in the coding sequence                                                                                                                                                                 |
| 10  | Protein_position   | Position or range in the protein                                                                                                                                                                         |
| 11  | Amino_acids        | Reference/variant peptides (`G/A`), one value when unchanged                                                                                                                                             |
| 12  | Codons             | Reference/variant codon windows, allele bases in upper case (`gGa/gCa`)                                                                                                                                  |
| 13  | Existing_variation | Co-located known variant IDs, comma-separated                                                                                                                                                            |
| 14  | Extra              | `KEY=VALUE` pairs joined by `;`, in VEP's key order (see below); `-` when empty                                                                                                                           |

### Extra Field

Keys appear in VEP's order: `IMPACT`, `DISTANCE`, `STRAND`, `FLAGS` (preceded by `ALLELE_NUM` under `--allele_number` and `REF_ALLELE` under `--show_ref_allele`), then the fields of the active flags in VEP's flag order (`VARIANT_CLASS`, `SPDI`, `GA4GH_VRS`, `CLINVAR_SOMATIC_CLASSIFICATION`, `SYMBOL`, `SYMBOL_SOURCE`, `HGNC_ID`, `BIOTYPE`, `CANONICAL`, `MANE`, `TSL`, `APPRIS`, `CCDS`, `ENSP`, `NEAREST`, `SIFT`, `PolyPhen`, `EXON`, `INTRON`, `DOMAINS`, `HGVSc`, `HGVSp`, `HGVS_OFFSET`, allele frequencies, `MAX_AF`, `FREQS`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `PUBMED`, ..., `SHIFT_LENGTH`, `VAR_SYNONYMS`), then any other key alphabetically (the structural-variant `OverlapBP`/`OverlapPC` and plugin keys). Inside a value `;` is escaped as `%3B` and, unless `--no_escape`, the `=` of an HGVSp value as `%3D`; lists join with `,`. `--terms display` or `--terms NCBI` renders every consequence term in that vocabulary, repeats included on a feature row (`SPLICE_SITE,SPLICE_SITE,INTRONIC`) and each label once on a `--summary` row, and leaves the ranks, the impacts and the most severe term on the Sequence Ontology term; `DISTANCE` is printed only when the rendered list names `upstream_gene_variant` or `downstream_gene_variant`, so the display vocabulary, whose labels are `UPSTREAM` and `DOWNSTREAM`, never prints it, as Ensembl VEP does not. `--total_length` appends `/<total>` to `cDNA_position` (the transcript's exon length sum), `CDS_position` (the coding sequence length) and `Protein_position` (the peptide length) wherever the column applies, a `-` placeholder included (`-/1200`). An intergenic row carries `IMPACT=MODIFIER` and no `STRAND`, `DISTANCE` or `FLAGS`; record-level keys such as `VARIANT_CLASS`, `NEAREST`, `CLIN_SIG` and variant-level plugin keys join it under their flags.

`HGVS_OFFSET` is the number of bases the HGVS notations moved an insertion or deletion 3' along the transcript, negative on the reverse strand, present only when something moved; under `--shift_hgvs 0` nothing moves and the key never appears. `SHIFT_LENGTH` (`--shift_length`) is the number of bases `--shift_3prime 1` moved the allele before its terms and positions were computed, written on every transcript row, `0` included, and on no intergenic row; without `--shift_3prime 1` the column is listed and left empty. Under `--shift_3prime 1 --shift_genomic 1` the `Location` column prints the moved position (a transcript row along its transcript's strand, an intergenic row towards higher coordinates) and `Allele` the inserted bases rotated along the move; `Uploaded_variation` keeps the input's name.

The co-located variant keys follow Ensembl VEP's `add_colocated_variant_info`. A row lists the known variants matched on its own allele or on their span alone, in Ensembl VEP's order (germline before somatic, then dbSNP, HGMD, COSMIC, other). `SOMATIC` and `PHENO` carry one `0` or `1` per listed variant and are left out when every value is `0`. `PUBMED` lists every id of every variant, duplicates kept. `VAR_SYNONYMS` (under `--var_synonyms`, which switches `--check_existing` on) joins each variant's `source::id,id` string with `--`. `CLINVAR_SOMATIC_CLASSIFICATION` (under `--clinvar_somatic_classification`) is the cache's classification string as written, `phenotype(classification)` entries joined by `;`, the last listed variant's; Ensembl VEP builds that string from ClinVar's somatic records when it dumps a cache or reads its database, and vep-rs reads a cache only, so a classification the cache lacks is never recovered. `CLIN_SIG` reports, by default (`--clin_sig_allele 1`), the ClinVar entries recorded for the row's own allele when any listed variant carries one, `;`-joined across variants, each variant's terms `,`-joined; where a variant's ClinVar reference allele differs from the input's, its entry alleles are read reverse-complemented first. When no listed variant carries an entry for the row's allele, it lists the terms of the variants read before the first one with allele-specific entries, duplicates kept; under `--clin_sig_allele 0` it lists every listed variant's terms, duplicates kept.

The frequency keys follow Ensembl VEP's `add_colocated_frequency_data`. `--af` adds `AF`, `--af_1kg` the 1000 Genomes columns (`AF`, `AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF`, `SAS_AF`), `--af_gnomade` and its older spelling `--af_gnomad` (one or the other, never both) the gnomAD exome columns (`gnomADe_AF` and nine populations), `--af_gnomadg` the gnomAD genome columns (`gnomADg_AF` and ten populations, `gnomADg_AMI_AF` included), each flag switching `--check_existing` on. A row reports, for each listed variant, the frequency of the record allele matched to the row's own allele, as the cache wrote it (`0.005244`, `6.875e-05`), except `AF`, which prints with four decimals (`0.0068`); the values of every listed variant merge into one `,`-joined list, each value once. A biallelic record whose `AF` stores one allele yields the remainder `1 - total` for the other (`0.9932`); a variant matched on its span alone is read only under `--af` or `--max_af` and contributes only in that case, the `AF` entry of the row's own allele, stored or remainder. `--max_af` adds `MAX_AF`, the largest frequency over every population column of every group whether its flag is set or not, leaving out the combined populations (`AF`, `gnomADe_AF`, `gnomADg_AF`), and `MAX_AF_POPS`, the populations reaching it, in the order the groups sort (`af`, `af_1kg`, `af_gnomad`, `af_gnomade`, `af_gnomadg`) and each group lists its populations; a variant whose populations are all `0` writes `MAX_AF=0` with every zero population listed, and across variants an equal maximum appends its populations as they are. `--af_gnomad` on a cache whose declared variation columns lack `gnomADe` ends the run with `ERROR: [cache] gnomAD data is not available in this cache`, as Ensembl VEP ends it; the check runs before the FASTA and the input open and under `--show_cache_info` too, reads that spelling alone and is skipped under `--everything` (which switches `--af_gnomade` on and is refused beside `--af_gnomad` unless `--safe`), so `--af_gnomade` never checks. Under `--old_maf` the `AFR_AF`, `AMR_AF` and `EUR_AF` columns print empty on a release-116 cache, as they do in Ensembl VEP, whose rewrite of those columns leaves no allele to match.

`--check_frequency` (which `--filter_common` switches on, with `--freq_pop 1KG_ALL --freq_freq 0.01 --freq_gt_lt gt --freq_filter exclude`) filters the input before any transcript is read, as Ensembl VEP's variation source filters its buffer: for each record, the population `--freq_pop` names is read from every listed variant (the column of that name, else the column named by the second of its `_`-separated parts, so `gnomADe_AFR` falls back to the 1000 Genomes `AFR`; a variant on the other strand is read through its reverse complement, and the one allele a variant names no frequency for takes the remainder `1 - total`), each alternate allele takes its own value, else the value of the record allele it was matched to, else `NA`, and passes when the value is at or above `--freq_freq` under `gt` or at or below it under `lt`; a value of `0` is no value. Under `exclude` a record whose every allele passes is dropped, one with some passing alleles keeps the failing ones alone, and one with no data or no passing allele is kept whole; under `include` a record is kept only when every allele passes, so one without data is dropped, as is a non-variant record written through under `--allow_non_variant` unless a known variant lying exactly at its span gives it data: a null-allele variant (`COSMIC_MUTATION`, `HGMD_MUTATION`), whose one allele takes the remainder `1`, except under `--exclude_null_alleles`, which matches none, or any variant with data under `--no_check_alleles`; the record is then kept, as Ensembl VEP keeps it, since it has no allele to fail. A non-variant record with a true INFO/SVTYPE (not empty, not `0`) is dropped under `include` whatever lies at its span, as it is in Ensembl VEP, which reads it as a structural variant and co-locates no known variant with one. A record that loses alleles is renumbered (`ALLELE_NUM`) and its allele string rewritten to the kept alleles (the JSON `allele_string` and `AMBIGUITY` read the rewritten one, and so does `UPLOADED_ALLELE` unless the run is `--minimal` or the alleles differ in length with no `-` among them, when Ensembl VEP prints the uploaded spelling). `FREQS` lists, on every row of a record that had data, `POPULATION:ALLELE:frequency` for every alternate allele as uploaded, dropped alleles included, the population as `--freq_pop` was given in upper case (`AF:C:0.006789`, `GNOMADE_AFR:-:0.0001206`, `AF:A:NA`), in the alleles' order where Ensembl VEP writes them in hash order; the JSON format has no `FREQS`. `--freq_pop` accepts `1KG_ALL`, a declared variation column in any case (`af`, `AFR`, `gnomade_nfe`) and `1KG_` before one, and refuses every other name, the gnomAD genome populations included, with `ERROR: [cache] Invalid population <name> specified` before the input opens; on a cache without variation data the filter is inactive, as Ensembl VEP builds no variation source there. `1KG_ALL` reads the cache's `AF` column where Ensembl VEP 116.2 finds no data and filters nothing ([Intended divergences](intended-divergences.md#--filter_common-and---freq_pop-1kg_all-read-the-1000-genomes-af-column)). The records dropped whole are counted in the run's closing line as `filtered out`.

**SV** (under `--check_svs <file>`) lists the names of the known structural variants overlapping the record, sorted and unique, on every row of the record (the JSON key is `colocated_structural_variants`, a list at the record's top level). Ensembl VEP reads them from its variation database and refuses the flag offline; vep-rs takes the flag with a file, Ensembl's `homo_sapiens_structural_variations.vcf.gz` of the matching release with its `.csi` or `.tbi` beside it, and reads it offline (the flag is accepted under `--offline`). That file is Ensembl's dump of the same database query, and vep-rs reads from it what the query returns: supporting-evidence records (`Parent`) are left out, as the database query leaves them out; records of the `probe` class are left out, as the query leaves the class out (the dump carries none); failed and somatic structural variants are absent from the dump, as they are from the query; and a record's span is read as the dump writes it, `POS + 1` to `END`, which is the database's span on every chromosome 21 record of release 116. The dump lacks three kinds of mapping the database holds: the dump covers the 25 chromosomes (1 to 22, X, Y, MT) and no scaffold or patch, so a variant on an unplaced or unlocalised scaffold gets no name where the database holds a mapping (release 116 holds 27,305 top-level features on 360 other regions, 4,903 of them on 75 `GL` and `KI` scaffolds, 2,772 on the 45 such scaffolds the release-116 cache carries); a structural variant mapped to several locations is in the dump at one location per dump piece, so a further mapping the database holds may be missing (on chromosome 21 of release 116, 713 of the 97,616 database mappings, which lie under the same name's dump mapping except at 590 of the chromosome's 46,709,983 positions, where the database names a structural variant the dump lacks); and a structural variant starting at the first base of its chromosome is absent from the dump (six of release 116 on the chromosomes 4, 9, 12, X and Y, two of them spanning the whole of X and of Y, so every variant on those two chromosomes lacks one name the database reports). A record overlaps the input when their spans meet, an insertion between two bases overlapping the records that span both, on the input's own chromosome; a breakend allele's mate position is read as well, on the mate's chromosome, and a record's alleles share the names (the union over the record), as Ensembl VEP keeps them on the record. Ensembl VEP matches the records of a buffer to its variants without a chromosome check, so a buffer spanning two chromosomes can match across them there.

### Example

Two records of the GRCh37 golden corpus as Ensembl VEP 116.2 wrote them (`tests/golden/116/GRCh37/expected/default.txt.gz`): `1 213061271 806345 G C`, a missense change reported against four transcripts, and the ID-less multi-allelic site `21 10795228 . C G,T`, which overlaps no transcript:

```
806345	1:213061271	C	ENSG00000162769	ENST00000366971	Transcript	missense_variant	1433	1235	412	G/A	gGa/gCa	-	IMPACT=MODERATE;STRAND=1
806345	1:213061271	C	ENSG00000162769	ENST00000419102	Transcript	missense_variant	631	632	211	G/A	gGa/gCa	-	IMPACT=MODERATE;STRAND=1;FLAGS=cds_start_NF
806345	1:213061271	C	ENSG00000162769	ENST00000474693	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2364;STRAND=1
806345	1:213061271	C	ENSG00000162769	ENST00000483790	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
21_10795228_C/G/T	21:10795228	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
21_10795228_C/G/T	21:10795228	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

## Tab Format

Selected with `--tab` or `--output_format tab`. The 13 default columns followed by one column per Extra key the active flags can produce, in the Extra key order above, then one per built-in plugin field; `-` for an absent value, lists joined by `,`. The header repeats the default format's block with a `## <column> : <description>` line for every column (a plugin column takes the plugin's own description) and no `Extra column keys` block, and its column line names every column.

### Choosing the columns (`--fields`)

`--fields` names the columns, comma-separated, and the list stands in for the whole default set, defaults included, in the order given: `--fields Location,Allele,Feature,Consequence,SYMBOL` writes exactly those five columns, each with its description line. Any field a flag or a plugin produces is a valid name (a field whose flag is off prints `-`, so `SYMBOL` needs `--symbol`), and a name nothing produces prints `-` in every row under a `## <name> : ?` description line, as Ensembl VEP writes it. Nothing is trimmed and a trailing comma opens no column. An empty list writes the tab header block and no data line, as Ensembl VEP's output loop does (every row is the empty string, which ends it); the VCF output's entries are then empty (`CSQ=,,`) on both engines. That loop also ends at a row whose whole line is `0`, so `--tab --fields SIFT` makes Ensembl VEP stop at the first SIFT score of `0` and lose every later row; vep-rs writes every row. The list applies to the `tab` and `vcf` outputs; the default and JSON formats ignore it, as Ensembl VEP's do, and `--output_format parquet` refuses it, since the Parquet column set is fixed (`ERROR: [output] --fields applies to --tab and --vcf output; the Parquet column set is fixed`). The `GRCh38-writers` corpus runs a list that reorders the columns and names both a field whose flag is off and a field nothing produces, in every format.

## VCF Format

Selected with `--vcf` or `--output_format vcf`, either of which also switches on `--symbol`, `--biotype` and `--numbers` (in VEP only the `--vcf` flag does). One line per input record: the record's own columns, with every allele's consequences comma-joined into the INFO field named by `--vcf_info_field` (`CSQ` by default); an earlier entry under that name is replaced, unless `--keep_csq` keeps it. A VCF input's line is written back as the bytes it was read from, sample columns included.

A record of another input format (Ensembl default, region, HGVS, SPDI, ID) has no VCF line of its own, so one is built as Ensembl VEP builds it (`to_VCF_record`): CHROM, POS, the record's name as ID (`chr_start_alleles` for a record its ID column cannot name), REF, ALT, `.` for QUAL and FILTER, and the annotation alone as INFO. A balanced change keeps its position and alleles (`21 30282169 multi C A,T`). An insertion or deletion takes the base before it as an anchor, from `--fasta` when the run has one and `N` otherwise, and its POS moves one base back: the deletion `TA/-` at 30348410 is `21 30348409 NTA N` without a FASTA. A record on the reverse strand has its alleles reverse-complemented, in the ALT and REF columns and in the `Allele` subfield of its CSQ entries, where the default and tab outputs print the allele as written. A structural class word (`DEL`, `DUP`, `INS`, `TDUP`, `INV`, ...) is written as `<DEL>`, `<DUP>`, `<INS>` or `<TDUP>`, or as its class term in angle brackets (`<inversion>`), one base before its span with `END=` in INFO. The `GRCh38-input-ensembl`, `GRCh38-input-hgvs` and `GRCh38-input-spdi` corpora compare these lines with Ensembl VEP 116.2's.

### Header

The input's meta lines are kept (minus any earlier definition of the INFO field and any `##VEP` or `##vep-rs` line), a `##fileformat` line is added when the input lacks one, then the `##VEP=` line, the `##vep-rs=<version>` line, the field's `##INFO` definition, one `##<KEY>=<description>` line per field of each loaded built-in plugin, the command line and the column line:

```
##VEP="v116.2" API="v116" time="[TIME]" cache="[PATH]/116" 1000genomes="phase3" COSMIC="98" ...
##vep-rs=0.5.0
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">
##VEP-command-line='vep --assembly GRCh37 ...'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
```

Ensembl VEP's `##VEP=` line also carries its Perl API component versions (`ensembl=116.0d85231 ...`), which have no counterpart here, and its `##INFO=<ID=<short>_<field>,...>` lines for `--custom` sources are not written, since custom annotation sources are not supported.

### Keeping an earlier annotation (`--keep_csq`)

An input that Ensembl VEP or vep-rs annotated before carries a `CSQ=` entry and the header lines that describe it. Without `--keep_csq` the first entry under the field's name is removed from each line and its `##INFO` and `##VEP` header lines are dropped, so the output describes one annotation. With `--keep_csq` every input header line stays and the new entry is appended under the same name, so the header carries two `##INFO=<ID=CSQ,...>` lines (the old Format list and the new one) and two `##VEP=` lines, and a line's INFO column carries two `CSQ=` keys, the old entry first: this is Ensembl VEP's own output, reproduced as is, and a reader that keeps one value per INFO key sees only one of them. `--vcf_info_field ANN` names the new entry `ANN` and the new header line `##INFO=<ID=ANN,...>`; an earlier `CSQ=` entry is then another key, kept beside it whether or not `--keep_csq` is given. The `keep-csq` and `keep-csq-ann` sets of the `GRCh38-writers` corpus run both on an input Ensembl VEP 116.2 annotated first.

### CSQ Field Structure

The Format list is VEP's fixed prefix (`Allele` to `Existing_variation`) followed by the active flag fields and the built-in plugin fields, or the `--fields` list alone in its own order (a name nothing produces is an empty subfield). With the `--vcf` implied flags the default list is the one above.

The CSQ value of record 806345 above, from the corpus's `expected/vcf.vcf.gz`; the file joins its four entries with `,` on one line:

```
CSQ=C|missense_variant|MODERATE|FLVCR1|ENSG00000162769|Transcript|ENST00000366971|protein_coding|6/10||||1433|1235|412|G/A|gGa/gCa|||1||HGNC|24682,
    C|missense_variant|MODERATE|FLVCR1|ENSG00000162769|Transcript|ENST00000419102|protein_coding|5/9||||631|632|211|G/A|gGa/gCa|||1|cds_start_NF|HGNC|24682,
    C|downstream_gene_variant|MODIFIER|FLVCR1|ENSG00000162769|Transcript|ENST00000474693|processed_transcript|||||||||||2364|1||HGNC|24682,
    C|intron_variant&non_coding_transcript_variant|MODIFIER|FLVCR1|ENSG00000162769|Transcript|ENST00000483790|processed_transcript||2/5||||||||||1||HGNC|24682
```

`EXON` and `INTRON` are `n/total` or `a-b/total` over the exons and introns the variant strictly overlaps; under `--shift_3prime 1` the candidates are the features the moved allele overlaps and the numbering follows the input's position, as in Ensembl VEP, so a move across an exon boundary leaves both empty. An absent value is empty (not `-`), except `Allele`. An intergenic entry carries `Allele`, `Consequence` and `IMPACT` alone (`G|intergenic_variant|MODIFIER||||...`).

### Encoding

Inside a CSQ value: `,` and `|` become `&`, `;` becomes `%3B`, whitespace becomes `_`. Lists join with `&`. Whitespace in a gene symbol is percent-encoded as `%20` first, so it survives the `_` rule, unless `--no_escape`, which leaves it to that rule; `--no_escape` also keeps the `=` of an HGVSp value, written `%3D` otherwise.

## JSON Format

Selected with `--json` or `--output_format json`. JSON Lines, one object per input record (a multi-allelic record is one object whose consequence arrays carry every allele). `--json` implies `--no_escape`.

### Schema

The two corpus records above as vep-rs writes them, pretty-printed with the keys in reading order (vep-rs itself writes the keys of an object sorted; VEP's key order is Perl hash order and varies from object to object). Ensembl VEP 116.2's committed objects (`expected/json.jsonl.gz`) carry the same keys and values except `vep_rs_version`, which only vep-rs writes. The `input` value is the record's own VCF line; the first one's INFO column is abridged here.

```json
{
  "id": "806345",
  "seq_region_name": "1",
  "start": 213061271,
  "end": 213061271,
  "strand": 1,
  "allele_string": "G/C",
  "assembly_name": "GRCh37",
  "input": "1\t213061271\t806345\tG\tC\t.\t.\tAF_EXAC=1e-05;ALLELEID=794562;...",
  "transcript_consequences": [
    {
      "variant_allele": "C",
      "gene_id": "ENSG00000162769",
      "transcript_id": "ENST00000366971",
      "consequence_terms": ["missense_variant"],
      "cdna_start": 1433,
      "cdna_end": 1433,
      "cds_start": 1235,
      "cds_end": 1235,
      "protein_start": 412,
      "protein_end": 412,
      "amino_acids": "G/A",
      "codons": "gGa/gCa",
      "impact": "MODERATE",
      "strand": 1
    },
    {
      "variant_allele": "C",
      "gene_id": "ENSG00000162769",
      "transcript_id": "ENST00000419102",
      "consequence_terms": ["missense_variant"],
      "cdna_start": 631,
      "cdna_end": 631,
      "cds_start": 632,
      "cds_end": 632,
      "protein_start": 211,
      "protein_end": 211,
      "amino_acids": "G/A",
      "codons": "gGa/gCa",
      "impact": "MODERATE",
      "strand": 1,
      "flags": ["cds_start_NF"]
    },
    {
      "variant_allele": "C",
      "gene_id": "ENSG00000162769",
      "transcript_id": "ENST00000474693",
      "consequence_terms": ["downstream_gene_variant"],
      "impact": "MODIFIER",
      "distance": 2364,
      "strand": 1
    },
    {
      "variant_allele": "C",
      "gene_id": "ENSG00000162769",
      "transcript_id": "ENST00000483790",
      "consequence_terms": ["intron_variant", "non_coding_transcript_variant"],
      "impact": "MODIFIER",
      "strand": 1
    }
  ],
  "most_severe_consequence": "missense_variant",
  "vep_rs_version": "0.5.0"
}
```

```json
{
  "id": ".",
  "seq_region_name": "21",
  "start": 10795228,
  "end": 10795228,
  "strand": 1,
  "allele_string": "C/G/T",
  "assembly_name": "GRCh37",
  "input": "21\t10795228\t.\tC\tG,T\t100\tPASS\tAC=5,8;AF=0.000998403,0.00159744;AN=5008;NS=2504;DP=64221;EAS_AF=0,0;AMR_AF=0.0029,0;AFR_AF=0.0023,0.0008;EUR_AF=0,0;SAS_AF=0,0.0072;AA=.|||;VT=SNP;MULTI_ALLELIC",
  "intergenic_consequences": [
    { "variant_allele": "G", "consequence_terms": ["intergenic_variant"], "impact": "MODIFIER" },
    { "variant_allele": "T", "consequence_terms": ["intergenic_variant"], "impact": "MODIFIER" }
  ],
  "most_severe_consequence": "intergenic_variant",
  "vep_rs_version": "0.5.0"
}
```

An allele gets an `intergenic_consequences` entry only when it overlaps no transcript, so a record never carries both kinds of entry for the same allele.

### Field Details

- `id` is VEP's `variation_name`: the first `;`-joined ID of a VCF record or the whole ID column of an Ensembl-format line, the literal `.` when that is `.`, and the generated `chr_start_alleles` when it is absent, empty or `0` (an ID Perl reads as false). An HGVS or SPDI input's `id` is its notation and a region input's its line as written, as in VEP.
- `allele_string` is the record's working allele string (anchor-trimmed for a multi-allelic record, `REF/N[chr:pos[` for a paired breakend while a single breakend `N.` stays bare, the symbolic ALT for other structural variants). A `*` keeps its written place (`N/*/G`, `A/G/*/C/GA`) and counts as an allele of one base in the trim (`G/C/*` for `AG AC,*`), a `<NON_REF>` beside a bracket ALT stays in a breakpoint's string (`N/N[21:30425700[/<NON_REF>`), and a sequence record's alleles are upper case (`A/G` for `a g`), as Ensembl VEP writes them.
- Column renames: `Consequence` to `consequence_terms` (array), `Gene` to `gene_id`, `Allele` to `variant_allele`, `Feature` to `transcript_id`, `SYMBOL` to `gene_symbol`, `OverlapBP` to `bp_overlap`, `OverlapPC` to `percentage_overlap`, `ENSP` to `protein_id` (`ENSP00000286808.3` under `--protein_version`), `RefSeq` to `refseq_transcript_ids`; `Feature_type` is dropped; `YES` flags become `1`.
- Positions split into `*_start` and `*_end`, each written only when defined: `390-?` gives `cds_start` 390 alone and `?-335` gives `cdna_end` 335 alone, as Ensembl VEP 116 writes them. A single position fills both keys. Under `--total_length` the split is Ensembl VEP's own: a single position with its total is the string `"428/1200"` at both ends, and a range with its total (`428-430/1200`) prints its start, `428`, at both ends, because the `/1200` tail of the end is not a plain number. SIFT and PolyPhen split into `*_prediction` and `*_score`.
- Numeric-looking values are numbers except identifier keys; a zero fraction prints as an integer (`100`, not `100.0`).
- `gene_symbol`, `gene_symbol_source`, `biotype`, `variant_class` and the other flag-gated keys appear only under their flags. Absent or `-` values are omitted. `--fields` does not apply: every key the flags produce is written, as in Ensembl VEP. `--shift_genomic 1` adds a `location` to every transcript entry (the `Location` column, moved under `--shift_3prime 1`) and to the intergenic entry of an insertion or deletion `--shift_3prime 1` analysed.
- The keys VEP derives from the variant alone sit on the record object, lower-cased, and on no consequence entry: `variant_class` (`--variant_class`), `ambiguity` (`--ambiguity`), `minimised` (`--minimal`) and `check_ref` (`--dont_skip`). `uploaded_allele`, `allele_num` and `ref_allele` sit on every consequence entry of a sequence variant; a structural variant's entries carry `allele_num` alone, as Ensembl VEP writes them.
- `--terms display` or `--terms NCBI` renders `consequence_terms` and `most_severe_consequence` in that vocabulary; Ensembl VEP exits on a record with two or more distinct terms there, so this JSON is one it does not produce.
- Plugin fields, built-in or dylib, sit inside each consequence object; none is written at the record level.
- `vep_rs_version` names the vep-rs release that wrote the object, on every object; Ensembl VEP writes no such key ([Engine identity](#engine-identity)).

### Key set

The object's keys are Ensembl VEP's with one exception: `custom_annotations`, the record-level object Ensembl VEP writes for `--custom` sources, is never written, since custom annotation sources are not supported. A `colocated_variants` entry is written as Ensembl VEP's JSON factory writes it, whatever the flags: every set cache column of the known variant (`seq_region_name`, `id`, `start`, `end`, `strand`, `allele_string`, `clin_sig` and `pubmed` as lists, `somatic` and `phenotype_or_disease` when `1`, `clin_sig_allele`, `clin_sig_ref_allele`, `clinvar_somatic_classification` and `var_synonyms` as `{source: [ids]}`) and `frequencies`, the matched allele's values nested under each alternate allele with lower-cased population keys. The consequence objects carry none of the co-located variant or frequency keys of the text formats.

## Parquet Format

Selected with `--output_format parquet`. Writes a TSV intermediate (one row per allele x transcript consequence, in VEP's row order) at `<output>.tsv.tmp` and post-processes it to Parquet via the `duckdb` CLI; the intermediate is deleted on success unless the hidden `--keep_intermediate` is given. The output is a directory, so `-o STDOUT` is refused (`ERROR: [output] --output_format parquet needs a file path`); `--fields` and `--compress_output` are refused with it too (the column set is fixed, and the files carry their own compression), and `--no_headers` leaves the intermediate's column line in place, since DuckDB reads the schema from it.

**Requires `duckdb` 1.2 or newer on PATH.** Install via `brew install duckdb` (macOS) or the CLI archive from <https://duckdb.org/docs/installation/> (Linux: the `duckdb_cli-linux-*.zip` release asset). vep-rs checks for it after the FASTA, plugins and variation cache are loaded and before the first variant is read, and fails there if `duckdb` is missing or too old, so no annotation is spent on a run it cannot finalize.

### Directory layout

A Hive-partitioned directory keyed by `chrom`, with the partition column also stored inside the files:

```
out.parquet/
├── chrom=1/data_0.parquet
├── chrom=2/data_0.parquet
└── chrom=X/data_0.parquet
```

Read it with `read_parquet('out.parquet/**/*.parquet', hive_partitioning=false)` to keep `chrom` as the stored VARCHAR (Hive type inference would otherwise cast `chrom=21` to a number).

### Schema

Five key columns, then exactly the default `--tab` columns in `--tab` order:

| Column                                       | Type     | Notes                                         |
| -------------------------------------------- | -------- | --------------------------------------------- |
| `chrom`                                      | VARCHAR  | input chromosome name                         |
| `pos`, `end`                                 | BIGINT   | the record's VEP start and end                |
| `ref`, `alt`                                 | VARCHAR  | this row's reference and alternate allele     |
| `Uploaded_variation` ... `Existing_variation` | VARCHAR  | the 13 default columns                        |
| `IMPACT`, `FLAGS`, flag and plugin fields    | VARCHAR  | `-` is stored as NULL                         |
| `DISTANCE`                                   | INTEGER  |                                               |
| `STRAND`                                     | TINYINT  |                                               |

### Row shape (`--parquet_shape`)

- `nested` (default): one row per (record, allele); `Uploaded_variation`, `Location`, `Allele` and `Existing_variation` stay scalar and every per-consequence column is a LIST in the flat row order. Two input records at one site stay two rows when they differ in ID, END or another per-record column; identical ID-less duplicate lines merge into one row. Query with `UNNEST` or `list_contains`:

  ```sql
  SELECT chrom, pos, ref, alt
  FROM read_parquet('out.parquet/**/*.parquet', hive_partitioning=false)
  WHERE list_contains(Feature, 'ENST00000366667')
    AND list_contains(IMPACT, 'HIGH');
  ```

- `flat`: one row per (allele x transcript consequence); every column scalar, so predicates push down directly:

  ```sql
  SELECT chrom, pos, ref, alt, Consequence, IMPACT
  FROM read_parquet('out.parquet/**/*.parquet', hive_partitioning=false)
  WHERE Feature = 'ENST00000366667' AND IMPACT = 'HIGH';
  ```

### Writer options

Chosen for the reader, not the writer: Parquet v2 pages, ZSTD level 9, row groups of 122,880 rows sorted by (chrom, pos, ref, alt) so the `pos` statistics are tight, a dictionary on every column (DuckDB writes a Bloom filter only for dictionary-encoded columns, so the dictionary size limit is set to the row-group size; the deprecated `DICTIONARY_COMPRESSION_RATIO_THRESHOLD` option is also passed and DuckDB 1.2 and later ignore it), Bloom filters at a 0.1% false-positive rate, and the run's identity in the footer's key-value metadata (`vep_rs_version`, `vep_api_version`, `assembly`, `cache`, `command_line`, every `source_<key>` of the cache). In the nested shape a LIST column carries a Bloom filter only when its elements are dictionary-encoded, which needs their distinct count in the row group to fit that limit (DuckDB 1.2 wrote no Bloom filter for nested types at all), so a high-cardinality column such as `Feature` usually has none; the scalar columns, whose distinct count cannot exceed the row count, keep theirs. `parquet_round_trips_to_the_tab_output` in `crates/vep-cli/tests/format_parity.rs` asserts, on every corpus, that both shapes round-trip to the `--tab` rows and, on the flat shape, that every populated column of every row group carries a dictionary, a Bloom filter and statistics and that the footer holds `vep_rs_version`, `assembly` and `command_line`; the page version, compression, row-group size, sort order and false-positive rate are not asserted.

### Round-trip to `--tab`

`scripts/adapters/parquet_to_vep_tab.py <dir|file> [--output tab.txt]` projects the tab columns back out (NULL to `-`, integers to text, nested rows unnested). Rows come out in Parquet order; compare against a `--tab` file as a multiset of lines.

### Debug flags

- `--keep_intermediate` (hidden): preserves the TSV intermediate after DuckDB finalize for inspection.
- `--parquet_row_group_size N` (hidden): rows per row group, rounded up to a multiple of 2,048 (DuckDB writes whole vectors; the dictionary size limit follows the rounded value, otherwise a row group with a distinct value in every row loses its dictionary and its Bloom filter). `parquet_row_groups_prune_on_bloom_filters` writes a corpus at 2,048-row groups and checks with `parquet_bloom_probe` that an absent `Feature` inside the column's value range is excluded by every row group and a single-row `Feature` by every row group but its own.

## Plugin Data in Output

Built-in plugin annotations appear in all output formats:

| Format  | Variant-level plugins (CADD, gnomADc)                        | Transcript-level plugins (REVEL, SpliceAI)   |
| ------- | ------------------------------------------------------------ | -------------------------------------------- |
| VEP     | Extra column key (`CADD_PHRED`)                              | Extra column key (`REVEL_score`)             |
| Tab     | Named column                                                 | Named column                                 |
| VCF     | Additional fields in CSQ                                     | Additional fields in CSQ                     |
| JSON    | Fields repeated in every `transcript_consequences` object    | Fields in the `transcript_consequences` object |
| Parquet | Additional columns                                           | Additional columns                           |

Variant-level plugin data is repeated on every transcript consequence row of that variant, as in Perl VEP; no format writes it at the record level.

A C-ABI dylib plugin (`--plugin <path>`, or `lib<Name>.so` / `lib<Name>.dylib` under `--dir_plugins`) runs per transcript consequence after the built-in plugins and returns a JSON object whose entries become fields of that consequence. Those fields appear in the default format's Extra column (alphabetically after the flag fields) and inside the JSON consequence objects only. The `tab`, `vcf` and Parquet column sets are fixed from the fields the built-in plugins declare before the first record is written, and the FFI contract has no header symbol through which a dylib could declare its own, so a dylib field has no column there and no `##` header line.
