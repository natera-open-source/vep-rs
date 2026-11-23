# Output formats

vep-rs supports five output formats, selected via `--output_format` or the shorthand flags `--vcf`, `--json`, `--tab`:

- `vep` (default): Ensembl VEP's default tab-delimited format with an `Extra` column; `ensembl` is Ensembl VEP's other name for it
- `tab`: Ensembl VEP's `--tab` format, every Extra key as its own named column
- `vcf`: the input VCF with a CSQ INFO field appended
- `json`: line-delimited JSON, one object per input record
- `parquet`: a directory of typed Parquet parts, one row per (allele, feature) pair with the `--tab` columns typed, written natively (see [Parquet Format](output/parquet.md)); vep-rs's own, with no Ensembl VEP counterpart

Four pages describe the formats, the default and tab formats sharing one. A fifth describes the records every format treats specially, structural variants. A sixth describes the known variants a record is matched against: co-located variants, their frequencies and known structural variants.

- [VEP default and tab](output/vep-and-tab.md): the header, the 14 columns, the Extra keys and their order, the `--tab` columns, and the rules behind the values a reader notices first.
- [VCF](output/vcf.md): the header lines, the 23 CSQ fields by position, and the CSQ encoding.
- [JSON](output/json.md): the record object and its consequence objects, key by key, with the renames from the column names.
- [Parquet](output/parquet.md): the directory of parts, the typed schema and its inverse to `--tab`, the ordering precondition, partitioning, overwriting and the reader recipes.
- [How structural variants are written](output/structural-variants.md): routing, spans, multi-allelic records, the records Ensembl VEP skips, and breakends.
- [Known variants](output/known-variants.md): the co-located variant keys, the frequency keys, the frequency filter and the `SV` column.

The `--output_format` value is matched in any case (`VCF`, `Json`) and read once, lower-cased, by every part of the run; Ensembl VEP lower-cases it to choose the writer alone and reads the spelling as typed in two places, a documented intended divergence, [`output_format_flag_read_in_its_given_case`](divergences/output-format-flag-case.md). A value outside the six is an argument error listing them. `--vcf`, `--json` and `--tab` refuse one another with Ensembl VEP's message (`Can't use --json and --vcf together`). `--vcf` with `--most_severe` or `--summary` is refused the same way ([CLI reference](cli-reference.md#incompatible-flags-configpm615-632)).

## Compressed output

`--compress_output gzip` writes the output as one gzip member and `--compress_output bgzip` as a BGZF stream: gzip members of at most 64 KiB ending in the empty end-of-file member, the form `bgzip`, `tabix` and `samtools` read. Either goes to the output file under the name given, or to standard output under `-o STDOUT`. Ensembl VEP pipes its output through the external `gzip -c` or `bgzip -c`. vep-rs encodes in process, so the decompressed text is the same output and the compressed bytes are the encoder's.

The BGZF encoder deflates its blocks on worker threads and writes them in order. The hidden `--compression_threads` sets how many. `0`, the default, is the number of available logical CPUs capped at 16, and `1` deflates on the writing thread, as the one-member gzip encoder does. Compressed bytes depend on the text alone; two runs of a command differ only in the header's time.

The flag applies to the default, `tab`, `vcf` and `json` formats and is refused with `parquet`, whose parts are ZSTD-compressed by the writer itself. The golden corpus `GRCh38-writers` runs both encoders in every text format and compares the decompressed output with Ensembl VEP's.

## Control bytes in values

A value read from an input file can carry a control byte: an identifier or a chromosome name with a tab, a line end or a NUL in it, a plugin data field with one. Written as read, a tab or a line end would split the column or the row it sits in.

The default and tab formats write every byte from `0x00` to `0x1F`, and `0x7F`, as its percent-encoding. That is `%09` for a tab, `%0A` for a line end and `%00` for a NUL. The `%3B` the Extra column writes for a `;` in a value is the same rule. Inside a VCF `CSQ` value the whitespace rule comes first, so a tab is `_`, and every other control byte is percent-encoded the same way. The JSON output escapes them as JSON does.

The columns a VCF line copies from the input, the ID among them, stay the input's bytes, as every format's header lines do. The line built for a record of another input format has each column encoded like any value. Ensembl VEP writes the bytes as read ([the divergence page](divergences/control-bytes-in-output.md)).

## Parity with Ensembl VEP

With the default flag set and `--tab`, `--vcf` or `--json`, vep-rs matches Ensembl VEP exactly except for the documented divergences: the default, `tab`, `vcf` and `json` outputs are compared with what Ensembl VEP 116.3 writes field for field. Every column of every row is compared as text and, for JSON, every key of every object. Ensembl VEP writes JSON keys in Perl hash order, so the JSON is not byte-identical and key order is not compared. Header lines that carry a time, a cache path or the command line are compared by shape. Ensembl VEP's Perl API component version lines have no counterpart.

Every format also names the engine that wrote it, which Ensembl VEP's output does not. The comparators set that one line or key aside and compare everything else ([Engine identity](#engine-identity)).

The contract is enforced by the golden corpora under `tests/golden/<release>/<assembly>/` (`crates/vep-cli/tests/golden.rs` and `format_parity.rs`). Each corpus holds every consequence-set combination observed on the benchmark datasets. It is annotated with its own pruned cache and compared with Ensembl VEP's committed output in all four formats, plus the Parquet round trip. Consequence terms themselves are governed by the concordance measurement ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)).

The corpus manifest records the (Location, Allele, Feature) keys whose terms differ between the engines, and the rows only one engine emits, each under its divergence class. The differing keys of four corpora, none `unexplained_residual`:

| Corpus | Differing keys | Rows |
| --- | ---: | ---: |
| GRCh37 | 447 | 17,556 |
| GRCh37-hgvs | 9 | 2,818 |
| GRCh38 | 1,344 | 20,623 |
| GRCh38-release116 | 3 | 977 |

For a recorded key the parity tests check only that the row's consequence set is one of the recorded sets. They compare none of its other columns and skip the record's `most_severe_consequence` in JSON.

The tested flag set is the one the golden corpora run: `--offline --json_cache <dir> --species <species> --assembly <assembly> --buffer_size 5000 --no_stats --quiet` plus one format flag and the corpus manifest's own flags. The species and assembly are the manifest's (`homo_sapiens` when it names no species). A corpus whose cache sits in Ensembl VEP's `cache_root/<species>/<release>_<assembly>` layout takes `--offline --cache --dir_cache <root> --cache_version <release>` in place of `--json_cache <dir>`.

The `GRCh38-regulatory` corpus records two documented differences. Under `--regulatory` ([Regulatory and motif features](regulatory.md)), a multi-base allele on a reverse-strand motif scores unlike Ensembl VEP ([the divergence page](divergences/motif-score-reverse-strand.md)). Under the same flag, `--cell_type` reads an epigenome name holding a comma whole ([its page](divergences/cell-type-comma-names.md)).

### What the format flags and `--everything` switch on

`--vcf` also switches on `--symbol`, `--biotype` and `--numbers`, as in Ensembl VEP, and `--json` switches on `--no_escape`, as Ensembl VEP's `rest` option set does. `--output_format json` switches it on here as well, where Ensembl VEP leaves that output escaped ([Unescaped JSON under `--output_format json`](divergences/json-escaping.md)). `--output_format vcf` switches on the same fields as `--vcf`, where Ensembl VEP ties them to the flag alone ([The CSQ fields of a VCF output selected with `--output_format vcf`](divergences/output-format-vcf-fields.md)). `--everything` switches on exactly Ensembl VEP's set ([CLI reference](cli-reference.md#option-sets-configpm351-593)), `--regulatory` among them. So its rows include the regulatory and motif rows of a cache that carries regulatory features.

### Where each rule lives

The rules the formatters follow are Ensembl VEP's own, cited to their Perl sources in `crates/vep-io/src/output/fields.rs` (the single field table shared by every format), `crates/vep-effects/src/display.rs` (positions, codons and peptides) and the formatters beside them. The rules a reader of the output most often needs, and the page that states each:

- **Row order**, transcript, regulatory-feature and motif-feature rows, each kind feature-major and allele-minor, intergenic rows last: [VEP default and tab](output/vep-and-tab.md#row-order).
- **Uploaded_variation**, the first `;`-element of the ID, else `CHROM_POS_REF/ALT1/ALT2` as written: [VEP default and tab](output/vep-and-tab.md#uploaded_variation).
- **Positions**, `start-end` ranges with `?` for an unmappable end: [VEP default and tab](output/vep-and-tab.md#positions).
- **Amino_acids and Codons**, `ref/alt` peptides and codon windows: [VEP default and tab](output/vep-and-tab.md#amino_acids-and-codons).
- **DISTANCE**, the shortest gap between variant and transcript: [VEP default and tab](output/vep-and-tab.md#distance).
- **FLAGS, OverlapBP and OverlapPC**: [VEP default and tab](output/vep-and-tab.md#flags-overlapbp-and-overlappc).
- **ENSP** under `--protein`, structural-variant rows included: [VEP default and tab](output/vep-and-tab.md#ensp).
- **The row filters** `--no_intergenic` and `--coding_only`, applied before the pick family: [VEP default and tab](output/vep-and-tab.md#the-row-filters).
- **The pick family**, `--pick`, `--per_gene`, `--flag_pick` and their allele variants, ranked by `--pick_order` over a whole record: [VEP default and tab](output/vep-and-tab.md#the-pick-family).
- **`--most_severe` and `--summary`**, one row per record after the pick family: [VEP default and tab](output/vep-and-tab.md#most_severe-and-summary).
- **NEAREST**, the protein-coding transcription start sites nearest the record under `--nearest`: [VEP default and tab](output/vep-and-tab.md#nearest).
- **HGVSg**, the HGVS genomic notation of a sequence variant's allele under `--hgvsg`: [VEP default and tab](output/vep-and-tab.md#hgvsg).
- **SPDI**, the allele in SPDI notation under `--spdi`: [VEP default and tab](output/vep-and-tab.md#spdi).
- **GA4GH_VRS**, the GA4GH Variation Representation `Allele` derived from the SPDI under `--ga4gh_vrs`: [VEP default and tab](output/vep-and-tab.md#ga4gh_vrs).

#### Rules stated outside the default and tab page

- **Co-located variant keys** under `--check_existing`, `Existing_variation`, `CLIN_SIG`, `SOMATIC`, `PHENO` and `PUBMED`, with `VAR_SYNONYMS` and `CLINVAR_SOMATIC_CLASSIFICATION` under their own flags: [Known variants](output/known-variants.md#co-located-variant-keys).
- **Frequency keys**, `AF`, the 1000 Genomes and gnomAD columns, `MAX_AF` and `MAX_AF_POPS`: [Known variants](output/known-variants.md#frequency-keys).
- **Filtering by frequency**, `--check_frequency` and `--filter_common`, with the `FREQS` key: [Known variants](output/known-variants.md#filtering-by-frequency).
- **SV**, the known structural variants overlapping the record under `--check_svs`: [Known variants](output/known-variants.md#sv).
- **Structural routing**, what makes a record structural and how its type is read: [How structural variants are written](output/structural-variants.md#routing).
- **Structural spans**, where a symbolic record starts and ends: [How structural variants are written](output/structural-variants.md#spans).
- **Multi-allelic records**, sequence and structural: [How structural variants are written](output/structural-variants.md#multi-allelic-records).
- **`--minimal`**, which trims every bi-allelic record from both ends and splits a multi-allelic record whose alleles trim: [How structural variants are written](output/structural-variants.md#multi-allelic-records). A split allele's positions follow the format: [below](#a-split-allele-in-each-format).
- **`--minimal` under `--coding_only`**, the span a split allele's rows are judged over: [below](#a-split-allele-in-each-format).
- **`--minimal` on a structural record**: the record is annotated as written, where Ensembl VEP trims two symbolic ALTs against each other ([a symbolic record with two or more ALTs](divergences/minimal-symbolic-alleles.md)).
- **Records Ensembl VEP skips**, which vep-rs annotates: [How structural variants are written](output/structural-variants.md#skipped-records).
- **Breakends**, points, spans, mates and the rows on each side: [How structural variants are written](output/structural-variants.md#breakends).
- **Per-sample copies and the reference row**, `--individual`, `--individual_zyg` and `--process_ref_homs`: [below](#per-sample-copies-and-the-reference-row).
- **The transcript's reference allele**, `--use_transcript_ref` with `GIVEN_REF` and `USED_REF`: [below](#the-transcripts-reference-allele).
- **RefSeq and merged cache columns**, `REFSEQ_MATCH`, `SOURCE`, `REFSEQ_OFFSET` and `BAM_EDIT`, with the HGVSc of an edited RefSeq model and the term of a variant inside the bases its alignment deleted: [RefSeq and merged cache columns](output/refseq-and-merged.md).

### A split allele in each format

Under `--minimal`, a split allele's `cDNA_position`, `CDS_position` and `Protein_position` follow the format. The default and tab outputs print the record span's values, in the columns the allele's own span admits, as Ensembl VEP recomputes them after its rejoin (`OutputFactory.pm` lines 294-299). The VCF and JSON outputs print the allele's own, as their factories recompute them before it (`OutputFactory/VCF.pm` lines 305-308, `OutputFactory/JSON.pm` lines 217-220). Codons, Amino_acids and the terms are the trimmed allele's in every format. The HGVS notations shift as a bi-allelic record's do ([the HGVSp of a split allele](divergences/minimal-split-hgvsp.md)).

Under `--coding_only`, a split allele's rows are judged over the span the output format rejoins on. That span is the record's in the default and tab outputs and the allele's own in the VCF and JSON outputs. Ensembl VEP's default and tab factories clear the cached CDS coordinates after the rejoin (`OutputFactory.pm` lines 297-323). Its VCF and JSON factories clear them before it (`OutputFactory/VCF.pm` lines 305-308, `OutputFactory/JSON.pm` lines 217-220). Every row of a split record prints the record's own `REF_ALLELE` when the first allele's feature kept no row before the rejoin (`OutputFactory.pm` line 1357 at line 2298).

## Per-sample copies and the reference row

### One copy per individual (`--individual`)

Under `--individual all|<list>`, each record is written once per individual named, as Ensembl VEP copies the record's feature per individual. `all` is every sample column in header order; a list is read in its order.

The default and tab formats repeat `Uploaded_variation` and `Location` for each copy's rows. The VCF output writes one line per copy, the input record with that copy's `CSQ`. The JSON output writes one object per copy. Every row and object carries `IND` and `ZYG` (`ind`, `zyg` in JSON). `ZYG` is `HET` when the genotype holds more than one distinct allele, else `HOM`, with `REF` appended for an individual carrying no alternate allele.

A copy carries the alleles its genotype names:

- `./1` and a haploid `1` annotate the one allele as `HOM`;
- `1/2` annotates both alleles as `HET`;
- `*` is no alternate allele;
- `./.` writes no copy.

An individual carrying no alternate allele writes no default or tab row; its VCF line has no `CSQ` and its JSON object no consequence (`most_severe_consequence` `?`), unless `--process_ref_homs` gives it the reference row.

Under `--most_severe` or `--summary` such a copy writes the record's one row with `-` in every allele and feature column, as Ensembl VEP does. The copy's JSON object has no consequence entry, where Ensembl VEP's JSON writer dies ([a JSON summary row with no consequence](divergences/json-summary-without-term.md)).

A record whose every named individual has a missing genotype is written in no format. A structural variant is written once, without `IND` or `ZYG`, since Ensembl VEP builds it before the per-individual step.

#### `UPLOADED_ALLELE` and `--minimal` on a copy

`UPLOADED_ALLELE` on a copy is the record's allele string as uploaded when its alleles differ in length, and otherwise the copy's own string. That string is `C/C` for a reference-homozygous copy under `--process_ref_homs` and `C/T` for the carrier of one allele of `C/T,G`. Ensembl VEP prints the feature's string where it kept no original. Two per-individual readings differ from Ensembl VEP's: [an individual carrying two alternate alleles](divergences/individual-allele-order.md) and [an individual's copy of a record trimmed past the anchor chop](divergences/individual-trimmable-indel.md), an indel such as `ACAC/AC` with or without `--minimal` and a same-length record such as `GA/GT` under it.

`--individual` runs with `--minimal` under `--safe`, which skips the incompatibility check in both engines. The record is read under `--minimal` before it is copied. So a copy of a multi-allelic record keeps the record's chopped alleles, with no `MINIMISED`, unless it names two alleles that trim. Such a copy is split in both engines, in Ensembl VEP in Perl's hash order. That order decides its `REF_ALLELE`, `VARIANT_CLASS` and `UPLOADED_ALLELE` ([an individual's copy of a `--minimal` record under `--safe`](divergences/safe-individual-minimal.md)).

#### How the sample columns are read

The sample columns are read as ensembl-io reads them. A genotype bit is a Perl array index: `-1` is the last allele and `+1` is 1. A bit past the record's alleles is an empty allele, which makes its copy heterozygous with no alternate row (`3/0`) or vanishes when it trails (`0/3`).

A header repeating a sample name has one sample of that name, reading its last column. A list may name a sample twice and gets two copies. A header whose ninth column is not `FORMAT` names its samples from that column and yields no genotype, so no copy.

A genotype naming the `.` of a non-variant record is annotated as an allele by Ensembl VEP ([Defect J](divergences/defect-j.md)) and as none here. Under `--individual`, `--buffer_size` bounds the copies a batch holds, as Ensembl VEP's buffer counts them.

### One zygosity list per record (`--individual_zyg`)

Under `--individual_zyg all|<list>` each record is written once, as without the flag, and every row and object carries one `ZYG` list. The list holds a `sample:ZYG` entry per named sample whose genotype is called, with the zygosity read as under `--individual`. A sample with a missing genotype has no entry, a record no named sample has a genotype for carries no `ZYG`, and there is no `IND`.

The default and tab formats join the entries with `,`, the VCF's `CSQ` with `&`, and the JSON `zyg` is an array. The entries come in the order the samples are named: the header's column order under `all` and the list's order under a list. A name listed twice has one entry. Ensembl VEP prints them in Perl's hash order ([the ZYG list order](divergences/individual-zyg-order.md)).

A run naming one sample, by a list of one name or by `all` over a single-sample VCF, reads a record that sample carries no alternate allele of (`*` counting as none) as the reference twice. Ensembl VEP then annotates the reference allele alone, so the record is its reference row, as under `--process_ref_homs`, with the one entry.

The reading above of a record trimmed past the anchor chop, an indel or under `--minimal` a same-length record such as `GA/GT`, applies to the list too. Ensembl VEP reads a `0/0` sample of it as `HOM` and keeps the alternate allele's rows, where vep-rs reads `HOMREF` and annotates the reference allele. A `1/1` sample of a deletion whose chopped alternate equals the minimised reference, such as `ATT/AT`, Ensembl VEP reads as `HOMREF`, annotating the reference allele, where vep-rs reads `HOM` and keeps the deletion's rows ([a record trimmed past the chop, per individual](divergences/individual-trimmable-indel.md)).

A structural variant carries no `ZYG`. `--individual_zyg` implies `--allow_non_variant` and cannot be combined with `--individual`. A record is read whole under the flag, its sample columns included, and a batch holds `--buffer_size` records. So the memory a batch takes grows with the sample count; a smaller `--buffer_size` bounds it on a many-sample VCF.

### The reference row (`--process_ref_homs`)

Under `--process_ref_homs`, every feature of every sequence record gains the reference allele's row first, `ALLELE_NUM` 0, its `Allele` the reference. The row carries the consequences the predicates give an allele equal to the reference:

- `synonymous_variant` in a coding sequence;
- `stop_retained_variant` at a stop codon;
- the single reference peptide in `Amino_acids`;
- the reference codon twice in `Codons`.

A reference allele of two or more bases carries no intron or splice-site term, since Ensembl VEP reads those over the bases where the allele differs from the reference. A row left with no term prints `intergenic_variant`. An intergenic record gains an `intergenic_variant` row for the reference too. Under `--coding_only` the reference row stays where its feature has an alternate allele in the coding sequence, since Ensembl VEP's `affects_cds` asks the alternates alone. It goes with a reference-homozygous copy, a non-variant record and a `--minimal` split feature whose allele equals the reference or is `*` or `<DEL:*>`.

Under `--minimal`, a multi-allelic record that splits writes one reference row per alternate allele, before that allele's rows and under its key. The row's allele is the allele's own trimmed reference over its trimmed span, `-` for an insertion's. An allele equal to the reference once trimmed writes its reference row alone, as does a `*` or `<DEL:*>` allele. The record's one intergenic row for the reference comes from its first intergenic allele.

The `REF_ALLELE` of every row of a split record under `--process_ref_homs` is the first allele's trimmed reference, `-` when that allele equals the reference, or the record's own reference where `--coding_only` drops that equal allele's feature. The record's own reference is printed whenever the first allele's feature kept no row before the rejoin, as under `--coding_only` for a first allele outside the coding sequence. One VCF shape differs under `--lookup_ref`: [a repeated alternate allele the looked-up reference equals](divergences/minimal-lookup-ref-repeated-allele.md).

The reference row has no `HGVSc` or `HGVSp`; under `--hgvsg` it prints the alternate allele's `HGVSg` on a bi-allelic record, as Ensembl VEP does. Under `--spdi` it carries no `SPDI`, except where the record's allele string lists the reference again after the first allele. A reference-homozygous individual's copy does in Ensembl VEP, unless `--lookup_ref` rewrites it to the looked-up reference followed by the written one. The row then prints the reference described against itself (`21:99:A:A`; an insertion's empty reference prints as `21:30284295::-`).

## The transcript's reference allele

Under `--use_transcript_ref`, each row is read against the reference allele Ensembl VEP takes from the transcript rather than from the input. Every format carries `GIVEN_REF`, the input's reference allele, and `USED_REF`, the row's. The JSON keys are `given_ref` and `used_ref`, and a `-` under them is kept, not omitted. The row's reference allele is the transcript's spliced cDNA over the variant where the variant lies wholly inside the exons, and `-` for an insertion. Elsewhere, in an intron, across an exon boundary or beyond the transcript, it is the input's allele: Ensembl VEP reads no genomic sequence for it.

The alleles are compared against that reference. An alternate allele equal to the transcript's sequence takes the terms an allele equal to the reference takes. It loses the intron and splice terms a multi-base allele earns only where it differs. Under `--process_ref_homs` the reference row's `Allele` is the transcript's allele, while `REF_ALLELE`, `UPLOADED_ALLELE` and `Uploaded_variation` keep the input's.

The HGVS notations shift a deletion towards the transcript's 3' end along the bases the transcript names. So a deletion whose input misspells the deleted bases shifts as the transcript reads it.

The flag needs `--fasta`, which supplies a non-coding transcript's cDNA; a coding transcript's comes from the cache. `--use_given_ref` given beside it switches the flag off. A BAM-edited cache, a converted RefSeq or merged cache whose metadata names the alignment file its models were corrected against, switches the flag on by itself with an `INFO:` line, as in Ensembl VEP, and adds `BAM_EDIT` ([Cache conversion](cache-conversion.md#refseq-and-merged-caches)). An edited transcript's allele is cut from its edited cDNA, coding or not. `--bam`, which edits the models at run time, is refused.

### Where `USED_REF` follows a deletion's shift

`USED_REF` follows a deletion's 3' shift only where Ensembl VEP rotates the reference allele along with it, through its codon and protein paths. That is a row whose shifted span has both protein coordinates, under `--shift_3prime 1` or `--hgvs` (`TranscriptVariationAllele.pm` lines 674, 813-820 and 1771-1786). There the deleted bases are rotated by the shift towards the transcript's 3' end, so a deletion given as `CA` and shifted one base prints `USED_REF=AC`. A non-coding row, a UTR row and a deletion whose shifted span leaves the coding sequence keep the bases as the input wrote them.

Under `--process_ref_homs` the codon path rotates nothing, because the reference row reads the codon before any alternate allele hands it a shift (`codon` rotates only on the call that memoises it, lines 802-804); `--hgvs` still rotates.

## Engine identity

The Ensembl-shaped header lines stay as Ensembl VEP writes them, so a parser written for Ensembl VEP output reads vep-rs output unchanged. Each format adds one item naming the vep-rs version: the release, `1.0.0` in the examples below, read from the crate version `vep --version` prints.

| Format           | Where                                                                                                                                | Form                         |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------- |
| default, `--tab` | the second header line, after the `## ENSEMBL VARIANT EFFECT PREDICTOR` banner                                                       | `## vep-rs 1.0.0`            |
| `--vcf`          | the meta line after `##VEP=`, in the `##key=value` form VCF 4.3 section 1.4 defines (a bare `## vep-rs` line is not valid VCF)        | `##vep-rs=1.0.0`             |
| `--json`         | a key of every object, an unannotated record's input-only object included                                                            | `"vep_rs_version": "1.0.0"`  |
| Parquet          | the key-value metadata of every part's footer, and the manifest                                                                      | `vep_rs_version`             |

`--no_headers` drops the default, tab and VCF lines with the rest of the header. Re-annotating a VCF that vep-rs wrote replaces its `##vep-rs=` line along with the earlier `##VEP` lines, as a re-annotation by Ensembl VEP replaces its own. The golden-corpus comparators (`crates/vep-cli/tests/common/mod.rs`) set aside exactly one identity line or key of the version under test. They report a missing, repeated or differently versioned one as a header difference.

## Plugin data in output

Built-in plugin annotations appear in all output formats:

| Format  | Variant-level plugins (CADD, gnomADc)                     | Transcript-level plugins (REVEL, SpliceAI)     |
| ------- | --------------------------------------------------------- | ---------------------------------------------- |
| VEP     | Extra column key (`CADD_PHRED`)                           | Extra column key (`REVEL`)                     |
| Tab     | Named column                                              | Named column                                   |
| VCF     | Additional fields in CSQ                                  | Additional fields in CSQ                       |
| JSON    | Fields repeated in every `transcript_consequences` object | Fields in the `transcript_consequences` object |
| Parquet | Named column, typed as a string                           | Named column, typed as a string                |

Variant-level plugin data is repeated on every row of that variant whose kind the plugin declares, the regulatory, motif and intergenic rows of a `--regulatory` run included for a plugin declaring every kind, as in Ensembl VEP. No format writes it at the record level. A transcript-level plugin writes on transcript rows alone ([Rows a plugin writes on](plugins.md#rows-a-plugin-writes-on)).

A C-ABI dylib plugin (`--plugin <path>`, or `lib<Name>.so` / `lib<Name>.dylib` under `--dir_plugins`) runs after the built-in plugins on each row whose kind it lists in `vep_plugin_feature_types`. Without that symbol it runs on transcript rows alone, Ensembl VEP's default ([Writing a plugin](developers/writing-a-plugin.md#vep_plugin_feature_types-optional)). It returns a JSON object whose entries become fields of that row.

Those fields appear in the default format's Extra column (alphabetically after the flag fields) and inside the JSON consequence objects. Every key the plugin declares in `vep_plugin_header_info` is also a column of `--tab`, an entry of the VCF CSQ Format and a Parquet column, each with its header line. An undeclared key reaches the Extra column and the JSON objects alone ([Writing a plugin](developers/writing-a-plugin.md#what-run-returns-and-where-it-lands)). A key named as one of the thirteen default columns (`Codons`) replaces that column's value in every format and never reaches Extra. In JSON, a key named as a column the consequence object leaves out joins the object on the rows the plugin wrote it ([JSON format](output/json.md#field-details)).

## Custom annotation

`--custom` sources, their argument grammar, the way their records match the input, the values they carry and the keys, columns and header lines they add to each format are described on their own page, [Custom annotation](custom-annotation.md).
