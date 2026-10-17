# Output formats

vep-rs supports five output formats, selected via `--output_format` or the shorthand flags `--vcf`, `--json`, `--tab`:

- `vep` (default): VEP's default tab-delimited format with an `Extra` column
- `tab`: VEP's `--tab` format, every Extra key as its own named column
- `vcf`: the input VCF with a CSQ INFO field appended
- `json`: line-delimited JSON, one object per input record
- `parquet`: a Hive-partitioned Parquet directory produced via a `duckdb` subprocess (see [Parquet Format](output/parquet.md))

Four pages describe the formats, the default and tab formats sharing one. A fifth describes the records every format treats specially, structural variants:

| Page                                                                   | What it covers                                                                                     |
| ---------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| [VEP default and tab](output/vep-and-tab.md)                           | The header, the 14 columns, the Extra keys and their order, the `--tab` columns, and the rules behind the values a reader notices first |
| [VCF](output/vcf.md)                                                   | The header lines, the 23 CSQ fields by position, and the CSQ encoding                             |
| [JSON](output/json.md)                                                 | The record object and its consequence objects, key by key, with the renames from the column names  |
| [Parquet](output/parquet.md)                                           | The `duckdb` requirement, the directory layout, the schema, the two row shapes and the writer options |
| [How structural variants are written](output/structural-variants.md)   | Routing, spans, multi-allelic records, the records Ensembl VEP skips, and breakends                 |

## Parity with Ensembl VEP

With the default flag set and `--tab`, `--vcf` or `--json`, the default, `tab`, `vcf` and `json` outputs match what Ensembl VEP 116.2 writes field for field. Every column of every row is compared as text and, for JSON, every key of every object. VEP writes JSON keys in Perl hash order, so the JSON is not byte-identical and key order is not compared. Header lines that carry a time, a cache path or the command line are compared by shape; VEP's Perl API component version lines have no counterpart.

The contract is enforced by the golden corpora under `tests/golden/<release>/<assembly>/` (`crates/vep-cli/tests/golden.rs` and `format_parity.rs`). Each corpus holds every consequence-set combination observed on the benchmark datasets, annotated with the corpus's own pruned cache and compared with VEP's committed output in all four formats, plus the Parquet round trip. Consequence terms themselves are governed by the concordance measurement ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)).

The corpus manifest records the (Location, Allele, Feature) keys whose terms differ between the engines, and the rows only one engine emits, each under its divergence class. The differing keys number 117 of the 17,556 GRCh37 rows, 4 of the 2,818 GRCh37-hgvs rows, 323 of the 20,623 GRCh38 rows and 2 of the 977 GRCh38-release116 rows, none of them `unexplained_residual`. For a recorded key the parity tests check only that the row's consequence set is one of the recorded sets, compare none of its other columns, and skip the record's `most_severe_consequence` in JSON.

Outside the tested flag set there are known differences. `--regulatory` adds its always-empty fields to the headers and columns and warns that regulatory annotation is not implemented. `--overlaps` does not exist. `--output_format vcf` switches on the same fields as `--vcf`, where VEP ties them to the flag alone. Under `--check_existing` a row none of whose co-located variants is somatic or phenotype-linked prints `SOMATIC=0` and `PHENO=0`, where VEP prints neither key.

### Where each rule lives

The rules the formatters follow are VEP's own, cited to their Perl sources in `crates/vep-io/src/output/fields.rs` (the single field table shared by every format), `crates/vep-effects/src/display.rs` (positions, codons and peptides) and the formatters beside them. The rules a reader of the output most often needs, and the page that states each:

- **Row order**, transcript-major and allele-minor, intergenic rows last: [VEP default and tab](output/vep-and-tab.md#row-order).
- **Uploaded_variation**, the first `;`-element of the ID, else `CHROM_POS_REF/ALT1/ALT2` as written: [VEP default and tab](output/vep-and-tab.md#uploaded_variation).
- **Positions**, `start-end` ranges with `?` for an unmappable end: [VEP default and tab](output/vep-and-tab.md#positions).
- **Amino_acids and Codons**, `ref/alt` peptides and codon windows: [VEP default and tab](output/vep-and-tab.md#amino_acids-and-codons).
- **DISTANCE**, the shortest gap between variant and transcript: [VEP default and tab](output/vep-and-tab.md#distance).
- **FLAGS, OverlapBP and OverlapPC**: [VEP default and tab](output/vep-and-tab.md#flags-overlapbp-and-overlappc).
- **ENSP** under `--protein`, structural-variant rows included: [VEP default and tab](output/vep-and-tab.md#ensp).
- **Structural routing**, what makes a record structural and how its type is read: [How structural variants are written](output/structural-variants.md#routing).
- **Structural spans**, where a symbolic record starts and ends: [How structural variants are written](output/structural-variants.md#spans).
- **Multi-allelic records**, sequence and structural: [How structural variants are written](output/structural-variants.md#multi-allelic-records).
- **Records Ensembl VEP skips**, which vep-rs annotates: [How structural variants are written](output/structural-variants.md#skipped-records).
- **Breakends**, points, spans, mates and the rows on each side: [How structural variants are written](output/structural-variants.md#breakends).

## Plugin data in output

Built-in plugin annotations appear in all output formats:

| Format  | Variant-level plugins (CADD, gnomADc)                     | Transcript-level plugins (REVEL, SpliceAI)     |
| ------- | --------------------------------------------------------- | ---------------------------------------------- |
| VEP     | Extra column key (`CADD_PHRED`)                           | Extra column key (`REVEL_score`)               |
| Tab     | Named column                                              | Named column                                   |
| VCF     | Additional fields in CSQ                                  | Additional fields in CSQ                       |
| JSON    | Fields repeated in every `transcript_consequences` object | Fields in the `transcript_consequences` object |
| Parquet | Additional columns                                        | Additional columns                             |

Variant-level plugin data is repeated on every transcript consequence row of that variant, as in Perl VEP; no format writes it at the record level.

A C-ABI dylib plugin (`--plugin <path>`, or `lib<Name>.so` / `lib<Name>.dylib` under `--dir_plugins`) runs per transcript consequence after the built-in plugins and returns a JSON object whose entries become fields of that consequence. Those fields appear in the default format's Extra column (alphabetically after the flag fields) and inside the JSON consequence objects only. The `tab`, `vcf` and Parquet column sets are fixed from the fields the built-in plugins declare before the first record is written. The FFI contract has no header symbol through which a dylib could declare its own, so a dylib field has no column there and no `##` header line.
