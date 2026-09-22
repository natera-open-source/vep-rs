# Output Formats

vep-rs supports five output formats, selected via `--output_format` or the shorthand flags `--vcf`, `--json`, `--tab`:

- `vep` (default): VEP's default tab-delimited format with an `Extra` column
- `tab`: VEP's `--tab` format, every Extra key as its own named column
- `vcf`: the input VCF with a CSQ INFO field appended
- `json`: line-delimited JSON, one object per input record
- `parquet`: a Hive-partitioned Parquet directory produced via a `duckdb` subprocess (see [Parquet Format](#parquet-format) below)

## Parity with Ensembl VEP

With the default flag set and `--tab`, `--vcf` or `--json`, the default, `tab`, `vcf` and `json` outputs match what Ensembl VEP 115.2 writes field for field: every column of every row and, for JSON, every key of every object compared as text (VEP writes JSON keys in Perl hash order, so the JSON is not byte-identical and key order is not compared). Header lines that carry a time, a cache path or the command line are compared by shape, and VEP's Perl API component version lines have no counterpart. The contract is enforced by the golden corpora under `tests/golden/<release>/<assembly>/` (`crates/vep-cli/tests/golden.rs` and `format_parity.rs`): every consequence-set combination observed on the benchmark datasets, annotated with the corpus's own pruned cache and compared with VEP's committed output in all four formats, plus the Parquet round trip. Consequence terms themselves are governed by the concordance measurement ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)). The corpus manifest records the (Location, Allele, Feature) keys whose terms differ between the engines (146 of the 20,111 GRCh37 rows and 406 of the 29,252 GRCh38 rows; 16 and 40 of them are classified `unexplained_residual`) and the rows only one engine emits; for a recorded key the parity tests check only that the row's consequence set is one of the recorded sets, compare none of its other columns, and skip the record's `most_severe_consequence` in JSON. Outside the tested flag set there are known differences: `--regulatory` adds its (always empty) fields to the headers and columns and warns that regulatory annotation is not implemented, `--overlaps` does not exist, `--output_format vcf` switches on the same fields as `--vcf` (VEP ties them to the flag alone), and under `--check_existing` a row none of whose co-located variants is somatic or phenotype-linked prints `SOMATIC=0` and `PHENO=0` where VEP prints neither key.

The rules the formatters follow are VEP's own, cited to their Perl sources in `crates/vep-io/src/output/fields.rs` (the single field table shared by every format), `crates/vep-effects/src/display.rs` (positions, codons and peptides) and the formatters beside them. The ones a reader of the output most often needs:

- **Row order.** Rows of one input record are transcript-major (sorted by transcript stable id) and allele-minor (ALT order); intergenic rows come last. A record's rows are contiguous in the default and tab formats and become one line (VCF) or one object (JSON).
- **Uploaded_variation.** The input ID (the first when several are `;`-joined). An ID-less record is named `chr_start_alleles`, where the alleles are the raw REF and every ALT joined by `/` and `start` is the variant's start after VEP's anchor trimming (`21_10795228_C/G/T`; `21_47612513_C/CTG` for the insertion `21 47612512 . C CTG`); a symbolic structural variant is named by its ALT as written (`21_33867342_<CN2>`).
- **Positions.** `cDNA_position`, `CDS_position` and `Protein_position` are `start-end` ranges, a single value when both ends coincide, and `?` for an end that falls in an intron or outside the transcript, and for the CDS and protein columns also in a UTR (`?-335`, `493-?`). An insertion prints the two flanking bases (`1161-1162`). `cDNA_position` appears only for exon-overlapping rows; the CDS and protein columns, `Amino_acids` and `Codons` only for coding rows.
- **Amino_acids and Codons.** `ref/alt` peptides, collapsed to one when unchanged (`P` for a synonymous change); a partial codon appends `X` (`M/X`, `D/DX`); an allele whose codon window is empty prints `-` (`-/S` for an in-frame insertion, `M/-` for a whole-codon deletion). Codons are the reference codon window in lower case with the allele's own bases in upper case, on both sides (`gTGCTGGgc/ggc` for a deletion, `-/TCT` for an insertion). Structural variants carry positions but no codons.
- **DISTANCE** is the shortest gap between the variant and the transcript, printed on up- and downstream rows only. **FLAGS** lists every `cds_` transcript attribute (`cds_start_NF,cds_end_NF`). **OverlapBP** and **OverlapPC** appear on every structural variant row that overlaps its feature (not on breakends) in the default Extra column and in JSON (`bp_overlap`, `percentage_overlap`), and never in `tab`, `vcf` or Parquet: VEP's `--overlaps` flag, which adds them as columns there, is not implemented.
- **Multi-allelic records** are one variant, as in VEP without `--minimal`: no per-allele minimisation, every allele's rows in one record.
- **Unsupported structural variants.** A symbolic ALT VEP's type table lacks (`<CPX>`, `<NON_REF>`) and a structural variant wider than `--max_sv_size` (default 10,000,000 bases) are annotated in the default and tab formats, but the VCF output carries their line without CSQ and the JSON output writes the unsupported type as an `input`-only object and omits the oversize one, exactly as VEP does.

## VEP Default Format

One line per (allele, transcript) consequence, plus one intergenic line per allele that overlaps no feature.

### Header

```
## ENSEMBL VARIANT EFFECT PREDICTOR v115.2
## Output produced at [TIME]
## Using cache in [PATH]/115
## Using API version 115, DB version ?
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
## VEP command-line: vep --assembly GRCh37 --force_overwrite --input_file [PATH]/input.vcf --json_cache [PATH]/115 --offline --output_file [PATH]/out.txt
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
```

The source-version lines come from the cache's `info.json` (`source_versions`, written by `scripts/data/storable_to_json.pl` and `vep-cache-builder`, or by `scripts/golden/prune_json_cache.py --perl-info` from a VEP cache's `info.txt`). The `Extra column keys` block lists every key the active flags can produce (plugin keys are not listed). `--no_headers` suppresses the whole block, column line included.

### Columns

| #   | Column             | Description                                                                                                                                                                                              |
| --- | ------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Uploaded_variation | Input ID, else `chr_start_alleles` (raw REF and ALTs joined by `/`)                                                                                                                                      |
| 2   | Location           | `chr:pos` or `chr:start-end` (`chr:end-start` for an insertion)                                                                                                                                          |
| 3   | Allele             | The ALT allele of this row; for a structural variant the SO term of its class (`deletion`, `duplication`, `Alu_insertion`), except that a breakend keeps its literal ALT (`N.`, `A[21:27348902[`) and an unsupported `<CPX>` prints `CPX` |
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

Keys appear in VEP's order: `IMPACT`, `DISTANCE`, `STRAND`, `FLAGS` (preceded by `ALLELE_NUM` under `--allele_number` and `REF_ALLELE` under `--show_ref_allele`), then the fields of the active flags in VEP's flag order (`VARIANT_CLASS`, `SYMBOL`, `SYMBOL_SOURCE`, `HGNC_ID`, `BIOTYPE`, `CANONICAL`, `MANE`, `TSL`, `APPRIS`, `CCDS`, `ENSP`, `SIFT`, `PolyPhen`, `EXON`, `INTRON`, `DOMAINS`, `HGVSc`, `HGVSp`, allele frequencies, `MAX_AF`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `PUBMED`, ...), then any other key alphabetically (the structural-variant `OverlapBP`/`OverlapPC` and plugin keys). Inside a value `;` is escaped as `%3B` and, unless `--no_escape`, the `=` of an HGVSp value as `%3D`; lists join with `,`. An intergenic row carries `IMPACT=MODIFIER` and no `STRAND`, `DISTANCE` or `FLAGS`; record-level keys such as `VARIANT_CLASS`, `CLIN_SIG` and variant-level plugin keys join it under their flags.

### Example

Two records of the GRCh37 golden corpus as Ensembl VEP 115.2 wrote them (`tests/golden/115/GRCh37/expected/default.txt.gz`): `1 213061271 806345 G C`, a missense change reported against four transcripts, and the ID-less multi-allelic site `21 10795228 . C G,T`, which overlaps no transcript:

```
806345	1:213061271	C	ENSG00000162769	ENST00000366971	Transcript	missense_variant	1433	1235	412	G/A	gGa/gCa	-	IMPACT=MODERATE;STRAND=1
806345	1:213061271	C	ENSG00000162769	ENST00000419102	Transcript	missense_variant	631	632	211	G/A	gGa/gCa	-	IMPACT=MODERATE;STRAND=1;FLAGS=cds_start_NF
806345	1:213061271	C	ENSG00000162769	ENST00000474693	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2364;STRAND=1
806345	1:213061271	C	ENSG00000162769	ENST00000483790	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
21_10795228_C/G/T	21:10795228	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
21_10795228_C/G/T	21:10795228	T	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

## Tab Format

Selected with `--tab` or `--output_format tab`. The 13 default columns followed by one column per Extra key the active flags can produce, in the Extra key order above, then one per built-in plugin field; `-` for an absent value, lists joined by `,`. The header repeats the default format's block with a `## <column> : <description>` line for every column and no `Extra column keys` block, and its column line names every column.

## VCF Format

Selected with `--vcf` or `--output_format vcf`, either of which also switches on `--symbol`, `--biotype` and `--numbers` (in VEP only the `--vcf` flag does). One line per input record: the record's own columns, with every allele's consequences comma-joined into the CSQ INFO field; an earlier CSQ entry is replaced.

### Header

The input's meta lines are kept (minus any earlier CSQ definition or `##VEP` line), a `##fileformat` line is added when the input lacks one, then:

```
##VEP="v115.2" API="v115" time="[TIME]" cache="[PATH]/115" 1000genomes="phase3" COSMIC="98" ...
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">
##VEP-command-line='vep --assembly GRCh37 ...'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
```

### CSQ Field Structure

The Format list is VEP's fixed prefix (`Allele` to `Existing_variation`) followed by the active flag fields and the built-in plugin fields. With the `--vcf` implied flags the default list is the one above.

The CSQ value of record 806345 above, from the corpus's `expected/vcf.vcf.gz`; the file joins its four entries with `,` on one line:

```
CSQ=C|missense_variant|MODERATE|FLVCR1|ENSG00000162769|Transcript|ENST00000366971|protein_coding|6/10||||1433|1235|412|G/A|gGa/gCa|||1||HGNC|24682,
    C|missense_variant|MODERATE|FLVCR1|ENSG00000162769|Transcript|ENST00000419102|protein_coding|5/9||||631|632|211|G/A|gGa/gCa|||1|cds_start_NF|HGNC|24682,
    C|downstream_gene_variant|MODIFIER|FLVCR1|ENSG00000162769|Transcript|ENST00000474693|processed_transcript|||||||||||2364|1||HGNC|24682,
    C|intron_variant&non_coding_transcript_variant|MODIFIER|FLVCR1|ENSG00000162769|Transcript|ENST00000483790|processed_transcript||2/5||||||||||1||HGNC|24682
```

`EXON` and `INTRON` are `n/total` or `a-b/total` over the exons and introns the variant strictly overlaps. An absent value is empty (not `-`), except `Allele`. An intergenic entry carries `Allele`, `Consequence` and `IMPACT` alone (`G|intergenic_variant|MODIFIER||||...`).

### Encoding

Inside a CSQ value: `,` and `|` become `&`, `;` becomes `%3B`, whitespace becomes `_`. Lists join with `&`.

## JSON Format

Selected with `--json` or `--output_format json`. JSON Lines, one object per input record (a multi-allelic record is one object whose consequence arrays carry every allele). `--json` implies `--no_escape`.

### Schema

The two corpus records above as Ensembl VEP 115.2 wrote them (`expected/json.jsonl.gz`), pretty-printed with the keys in the order vep-rs emits them; VEP's key order is Perl hash order and varies from object to object. The `input` value is the record's own VCF line; the first one's INFO column is abridged here.

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
  "most_severe_consequence": "missense_variant"
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
  "most_severe_consequence": "intergenic_variant"
}
```

An allele gets an `intergenic_consequences` entry only when it overlaps no transcript, so a record never carries both kinds of entry for the same allele.

### Field Details

- `id` is VEP's `variation_name`: the input ID, the literal `.` for a VCF record without one, `chr_start_alleles` for an ensembl-format line without an ID. VEP names a region or HGVS input by the input string itself; vep-rs names those `chr_start_alleles` too.
- `allele_string` is the record's working allele string (anchor-trimmed for a multi-allelic record, `REF/N[chr:pos[` for a paired breakend while a single breakend `N.` stays bare, the symbolic ALT for other structural variants).
- Column renames: `Consequence` to `consequence_terms` (array), `Gene` to `gene_id`, `Allele` to `variant_allele`, `Feature` to `transcript_id`, `SYMBOL` to `gene_symbol`, `OverlapBP` to `bp_overlap`, `OverlapPC` to `percentage_overlap`, `ENSP` to `protein_id`, `RefSeq` to `refseq_transcript_ids`; `Feature_type` is dropped; `YES` flags become `1`.
- Positions split into `*_start` and `*_end`; a `?` end takes the start's value (`390-?` gives `cds_start` 390 and `cds_end` 390) and a `?` start is omitted (`?-335` gives `cdna_end` 335 alone). SIFT and PolyPhen split into `*_prediction` and `*_score`.
- Numeric-looking values are numbers except identifier keys; a zero fraction prints as an integer (`100`, not `100.0`).
- `gene_symbol`, `gene_symbol_source`, `biotype`, `variant_class` and the other flag-gated keys appear only under their flags. Absent or `-` values are omitted.
- Plugin fields, built-in or dylib, sit inside each consequence object; none is written at the record level.

## Parquet Format

Selected with `--output_format parquet`. Writes a TSV intermediate (one row per allele x transcript consequence, in VEP's row order) and post-processes it to Parquet via the `duckdb` CLI. **Requires `duckdb` 1.2 or newer on PATH**; the runner fails at startup otherwise.

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

Five key columns, then exactly the `--tab` columns in `--tab` order:

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

- `nested` (default): one row per (record, allele); `Uploaded_variation`, `Location`, `Allele` and `Existing_variation` stay scalar and every per-consequence column is a LIST in the flat row order. Two input records at one site stay two rows when they differ in ID, END or another per-record column; identical ID-less duplicate lines merge into one row.
- `flat`: one row per (allele x transcript consequence); every column scalar.

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
