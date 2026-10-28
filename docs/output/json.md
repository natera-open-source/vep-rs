# JSON format

This page describes the JSON output, one object per input record, and the keys of each object inside it. Selected with `--json` or `--output_format json`. The output is JSON Lines; a multi-allelic record is one object whose consequence arrays carry every allele. `--json` implies `--no_escape`.

## Schema

The two records of the [default format example](vep-and-tab.md#example) as vep-rs writes them, pretty-printed with the keys in reading order. vep-rs itself writes the keys of an object sorted; VEP's key order is Perl hash order and varies from object to object. Ensembl VEP 116.2's committed objects (`expected/json.jsonl.gz`) carry the same keys and values except `vep_rs_version`, which only vep-rs writes. The `input` value is the record's own VCF line; the first one's INFO column is abridged here.

```json
{
  "id": "806345",
  "seq_region_name": "1",
  "start": 213061271,
  "end": 213061271,
  "strand": 1,
  "allele_string": "G/C",
  "assembly_name": "GRCh37",
  "input": "1\t213061271\t806345\tG\tC\t.\t.\tAF_EXAC=1e-05;...",
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
    {
      "variant_allele": "G",
      "consequence_terms": ["intergenic_variant"],
      "impact": "MODIFIER"
    },
    {
      "variant_allele": "T",
      "consequence_terms": ["intergenic_variant"],
      "impact": "MODIFIER"
    }
  ],
  "most_severe_consequence": "intergenic_variant",
  "vep_rs_version": "0.5.0"
}
```

An allele gets an `intergenic_consequences` entry only when it overlaps no transcript, so a record never carries both kinds of entry for the same allele. Under `--regulatory` ([Regulatory and motif features](../regulatory.md)), an allele inside a regulatory feature or a binding motif adds a `regulatory_feature_consequences` or `motif_feature_consequences` entry, and keeps its intergenic entry when it overlaps no transcript.

### The record object

| Field                     | Type    | Notes                                                                                           |
| ------------------------- | ------- | ----------------------------------------------------------------------------------------------- |
| `id`                      | string  | VEP's `variation_name`; see Field Details                                                        |
| `seq_region_name`         | string  | The chromosome as the input names it                                                             |
| `start`, `end`            | integer | The record's VEP start and end                                                                  |
| `strand`                  | integer | The record's strand                                                                             |
| `allele_string`           | string  | The record's working allele string; see Field Details                                            |
| `assembly_name`           | string  | The assembly, as `--assembly` names it, else as the cache's `info.json` names it; see Field Details |
| `input`                   | string  | The record's own input line as Ensembl VEP stores it: trailing empty columns dropped and, for a structural record, QUAL, FILTER and INFO filled with `.` where absent, empty or `0` |
| `variant_class`           | string  | Under `--variant_class`                                                                          |
| `transcript_consequences` | array   | One object per (allele, transcript); absent when no allele overlaps a transcript                 |
| `intergenic_consequences` | array   | One object per allele that overlaps no transcript                                                |
| `regulatory_feature_consequences` | array | Under `--regulatory`, one object per (allele, regulatory feature); its id key is `regulatory_feature_id`, and `biotype` is the feature's type whatever `--biotype` says |
| `motif_feature_consequences` | array | Under `--regulatory`, one object per (allele, binding motif); its id key is `motif_feature_id`, with `motif_name`, `motif_pos` and `transcription_factors` (an array, empty when the matrix names none) beside `strand` |
| `most_severe_consequence` | string  | The most severe term over every entry; `?` for a record that [names no allele](structural-variants.md#records-that-name-no-allele) |
| `colocated_variants`      | array   | Under `--check_existing`, one object per co-located variant                                      |
| `vep_rs_version`          | string  | The vep-rs release that wrote the object; on every object, see Field Details                    |

### A `transcript_consequences` object

| Field                                                                | Type             | Notes                                                                                    |
| -------------------------------------------------------------------- | ---------------- | ---------------------------------------------------------------------------------------- |
| `variant_allele`                                                     | string           | The `Allele` column                                                                      |
| `gene_id`                                                            | string           | `Gene`                                                                                   |
| `transcript_id`                                                      | string           | `Feature`                                                                                |
| `consequence_terms`                                                  | array of strings | `Consequence`                                                                            |
| `cdna_start`, `cdna_end`                                             | integer          | `cDNA_position` split in two; each key only when its end is defined                      |
| `cds_start`, `cds_end`                                               | integer          | `CDS_position` split in two                                                              |
| `protein_start`, `protein_end`                                       | integer          | `Protein_position` split in two                                                          |
| `amino_acids`                                                        | string           | `Amino_acids`                                                                            |
| `codons`                                                             | string           | `Codons`                                                                                 |
| `impact`                                                             | string           | `IMPACT`                                                                                 |
| `strand`                                                             | integer          | `STRAND`                                                                                 |
| `distance`                                                           | integer          | `DISTANCE`; up- and downstream entries only                                              |
| `flags`                                                              | array of strings | `FLAGS`; absent when the transcript has none                                             |
| `gene_symbol`, `gene_symbol_source`, `biotype`                       | string           | `SYMBOL`, `SYMBOL_SOURCE`, `BIOTYPE`, each under its flag                                |
| `protein_id`                                                         | string           | `ENSP` under `--protein`; `ENSP00000286808.3` under `--protein_version`                  |
| `refseq_transcript_ids`                                              | array            | `RefSeq`                                                                                 |
| `sift_prediction`, `sift_score`, `polyphen_prediction`, `polyphen_score` | string, number | `SIFT` and `PolyPhen`, each split into a prediction and a score                           |
| `bp_overlap`, `percentage_overlap`                                   | integer; number  | `OverlapBP`, `OverlapPC`; structural-variant entries that overlap their feature           |
| Other flag fields                                                    | per field        | The lower-cased key of the field, under its flag; `YES` flags become `1`                 |
| Plugin fields                                                        | per plugin       | Built-in or dylib, inside the object; none is written at the record level                |

### An `intergenic_consequences` object

| Field               | Type             | Notes                                                 |
| ------------------- | ---------------- | ----------------------------------------------------- |
| `variant_allele`    | string           | The allele that overlaps no transcript                |
| `consequence_terms` | array of strings | `["intergenic_variant"]`                              |
| `impact`            | string           | `MODIFIER`                                            |

## Field details

`id` is VEP's `variation_name`: the first `;`-joined ID of a VCF record or the whole ID column of an Ensembl-format line. It is the literal `.` when that is `.`, and the generated `chr_start_alleles` when the ID is absent, empty or `0` (an ID Perl reads as false). An HGVS or SPDI input's `id` is its notation, and a region input's is its line as written, as in VEP. An identifier-list input's (`--format id`) is the identifier as written, whitespace removed, on every record the identifier maps to.

`allele_string` is the record's working allele string: anchor-trimmed for a multi-allelic record, `REF/N[chr:pos[` for a paired breakend while a single breakend `N.` stays bare, and the symbolic ALT for other structural variants. A `*` keeps its written place (`N/*/G`, `A/G/*/C/GA`) and counts as an allele of one base in the trim (`G/C/*` for `AG AC,*`). A `<NON_REF>` beside a bracket ALT stays in a breakpoint's string (`N/N[21:30425700[/<NON_REF>`). A sequence record's alleles are upper case (`A/G` for `a g`), as Ensembl VEP writes them.

Column renames: `Consequence` to `consequence_terms` (an array), `Gene` to `gene_id`, `Allele` to `variant_allele`, `Feature` to `transcript_id`, `SYMBOL` to `gene_symbol`, `OverlapBP` to `bp_overlap`, `OverlapPC` to `percentage_overlap`, `ENSP` to `protein_id` and `RefSeq` to `refseq_transcript_ids`. `Feature_type` is dropped. `YES` flags become `1`.

Positions split into `*_start` and `*_end`, each written only when defined: `390-?` gives `cds_start` 390 alone and `?-335` gives `cdna_end` 335 alone, as Ensembl VEP 116 writes them. A single position fills both keys. Under `--total_length` the split is Ensembl VEP's own: a single position with its total is the string `"428/1200"` at both ends. A range with its total (`428-430/1200`) prints its start, `428`, at both ends, because the `/1200` tail of the end is not a plain number. SIFT and PolyPhen split into `*_prediction` and `*_score`.

Numeric-looking values are numbers, except under the identifier keys; a zero fraction prints as an integer (`100`, not `100.0`). `gene_symbol`, `gene_symbol_source`, `biotype`, `variant_class` and the other flag-gated keys appear only under their flags. Absent or `-` values are omitted, except a `-` under `variant_allele`, `given_ref` and `used_ref`, where it is an empty allele. `--fields` does not apply: every key the flags produce is written, as in Ensembl VEP.

`--shift_genomic 1` adds a `location` to every transcript entry (the `Location` column, moved under `--shift_3prime 1`) and to the intergenic entry of an insertion or deletion `--shift_3prime 1` analysed. Plugin fields, built-in or dylib, sit inside each consequence object; none is written at the record level.

The keys VEP derives from the variant alone sit on the record object, lower-cased, and on no consequence entry: `variant_class` (`--variant_class`), `ambiguity` (`--ambiguity`), `minimised` (`--minimal`) and `check_ref` (`--dont_skip`). `uploaded_allele`, `allele_num` and `ref_allele` sit on every consequence entry of a sequence variant; a structural variant's entries carry `allele_num` alone, as Ensembl VEP writes them.

`assembly_name` is the `--assembly` value, else the assembly the cache's `info.json` names, under `--json_cache` and `--full_cache_dir` alike. Ensembl VEP reads the assembly from the cache directory's name, which `--full_cache_dir` bypasses, and its fallback to a `cache_assembly` parameter is declared and never set (`OutputFactory/JSON.pm` lines 168 and 242), so under `--full_cache_dir` with neither `--assembly` nor `--grch37` it prints `"assembly_name":null`. The same unset parameter leaves its pseudoautosomal redirect of `Y` reads unarmed and its `--fasta_dir` filter passing every file naming the species ([Reference FASTA](../cache-setup.md#reference-fasta)).

`--terms display` or `--terms NCBI` renders `consequence_terms` and `most_severe_consequence` in that vocabulary. Ensembl VEP exits on a record with two or more distinct terms there, so this JSON is one it does not produce.

`vep_rs_version` names the vep-rs release that wrote the object, on every object; Ensembl VEP writes no such key ([Engine identity](../output-formats.md#engine-identity)).

## Key set

The object's keys are Ensembl VEP's plus `vep_rs_version`, with one exception: `custom_annotations`, the record-level object Ensembl VEP writes for `--custom` sources, is never written, since custom annotation sources are not supported.

A `colocated_variants` entry is written as Ensembl VEP's JSON factory writes it, whatever the flags: every set cache column of the known variant, and `frequencies`. The columns are `seq_region_name`, `id`, `start`, `end`, `strand`, `allele_string`, `clin_sig` and `pubmed` as lists, `somatic` and `phenotype_or_disease` when `1`, `clin_sig_allele`, `clin_sig_ref_allele`, `clinvar_somatic_classification` and `var_synonyms` as `{source: [ids]}`. `frequencies` nests the matched allele's values under each alternate allele with lower-cased population keys. The consequence objects carry none of the co-located variant or frequency keys of the text formats.
