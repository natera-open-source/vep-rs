# Filtering, picking and transcript sets

These flags reduce the consequences reported per variant, allele or gene, flag the ones a pick would drop, and shape the allele columns. The transcript-set flags of Ensembl VEP are accepted and select nothing: the transcript set is the `--json_cache` contents.

## Filtering and picking

| Flag | Default | What it does |
| --- | --- | --- |
| [`--pick`](#opt_pick) | off | Picks one consequence per variant |
| [`--pick_allele`](#opt_pick_allele) | off | Picks one consequence per allele |
| [`--per_gene`](#opt_per_gene) | off | Picks one consequence per gene |
| [`--pick_allele_gene`](#opt_pick_allele_gene) | off | Picks one consequence per allele per gene |
| [`--most_severe`](#opt_most_severe) | off | Reports only the single most severe consequence |
| [`--summary`](#opt_summary) | off | Reports a summary of all consequences |
| [`--flag_pick`](#opt_flag_pick) | off | Flags, but keeps, the consequences a pick would drop |
| [`--flag_pick_allele`](#opt_flag_pick_allele) | off | Flags the non-picked consequences per allele |
| [`--flag_pick_allele_gene`](#opt_flag_pick_allele_gene) | off | Flags the non-picked consequences per allele per gene |
| [`--pick_order`](#opt_pick_order) | none | No effect; the pick order is fixed to Ensembl VEP's default |
| [`--allele_number`](#opt_allele_number) | off | Adds the allele number from the input |
| [`--show_ref_allele`](#opt_show_ref_allele) | off | Shows the reference allele in the output |
| [`--minimal`](#opt_minimal) | off | Adds the `MINIMISED` column |
| [`--total_length`](#opt_total_length) | off | No effect |

## Transcript sets

| Flag | Default | What it does |
| --- | --- | --- |
| [`--refseq`](#opt_refseq) | off | Selects no transcript set; suffixes the header's cache path with `_refseq` |
| [`--merged`](#opt_merged) | off | Selects no transcript set; suffixes the header's cache path with `_merged` |
| [`--gencode_basic`](#opt_gencode_basic) | off | No effect |

## Accepted for compatibility

| Flag | What Ensembl VEP does with it |
| --- | --- |
| `--pick_order` | Sets the order of the criteria the pick flags apply |
| `--total_length` | Writes cDNA, CDS and protein positions as position/length |
| `--refseq` | Reads the RefSeq transcript cache |
| `--merged` | Reads the merged Ensembl and RefSeq cache |
| `--gencode_basic` | Limits the analysis to transcripts in the GENCODE basic set |

## Options

### `--pick` {#opt_pick}

**`--pick`**, boolean, default off. Ensembl VEP: same flag and default.

Picks one consequence per variant.

### `--pick_allele` {#opt_pick_allele}

**`--pick_allele`**, boolean, default off. Ensembl VEP: same flag and default.

Picks one consequence per allele.

### `--per_gene` {#opt_per_gene}

**`--per_gene`**, boolean, default off. Ensembl VEP: same flag and default.

Picks one consequence per gene.

### `--pick_allele_gene` {#opt_pick_allele_gene}

**`--pick_allele_gene`**, boolean, default off. Ensembl VEP: same flag and default.

Picks one consequence per allele per gene.

### `--most_severe` {#opt_most_severe}

**`--most_severe`**, boolean, default off. Ensembl VEP: same flag and default.

Reports only the single most severe consequence.

### `--summary` {#opt_summary}

**`--summary`**, boolean, default off. Ensembl VEP: same flag and default.

Reports a summary of all consequences.

### `--flag_pick` {#opt_flag_pick}

**`--flag_pick`**, boolean, default off. Ensembl VEP: same flag and default.

Flags, but keeps, the consequences a pick would drop.

### `--flag_pick_allele` {#opt_flag_pick_allele}

**`--flag_pick_allele`**, boolean, default off. Ensembl VEP: same flag and default.

Flags the non-picked consequences per allele.

### `--flag_pick_allele_gene` {#opt_flag_pick_allele_gene}

**`--flag_pick_allele_gene`**, boolean, default off. Ensembl VEP: same flag and default.

Flags the non-picked consequences per allele per gene.

### `--pick_order` {#opt_pick_order}

**`--pick_order <list>`**, string, no default. Ensembl VEP: accepted, no effect.

Not implemented. The pick order is fixed to Ensembl VEP's default: `mane_select`, `mane_plus_clinical`, `canonical`, `appris`, `tsl`, `biotype`, `ccds`, `rank`. The final `length` tie-break is not applied.

### `--allele_number` {#opt_allele_number}

**`--allele_number`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the allele number from the input.

### `--show_ref_allele` {#opt_show_ref_allele}

**`--show_ref_allele`**, boolean, default off. Ensembl VEP: same flag and default.

Shows the reference allele in the output.

### `--minimal` {#opt_minimal}

**`--minimal`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the `MINIMISED` column. Bi-allelic indels are always minimised, as Ensembl VEP does regardless of the flag; per-allele minimisation of multi-allelic records is not implemented.

### `--total_length` {#opt_total_length}

**`--total_length`**, boolean, default off. Ensembl VEP: accepted, no effect.

Not implemented.

### `--refseq` {#opt_refseq}

**`--refseq`**, boolean, default off. Ensembl VEP: accepted, no effect.

Selects no transcript set: the transcript set is the `--json_cache` contents. It only suffixes the species directory of the header's cache path, `homo_sapiens_refseq`, when `--json_cache` is absent.

### `--merged` {#opt_merged}

**`--merged`**, boolean, default off. Ensembl VEP: accepted, no effect.

Selects no transcript set: the transcript set is the `--json_cache` contents. It only suffixes the species directory of the header's cache path, `homo_sapiens_merged` (`homo_sapiens_refseq_merged` with `--refseq`), when `--json_cache` is absent.

### `--gencode_basic` {#opt_gencode_basic}

**`--gencode_basic`**, boolean, default off. Ensembl VEP: accepted, no effect.

Not implemented: the transcript set is the `--json_cache` contents.
