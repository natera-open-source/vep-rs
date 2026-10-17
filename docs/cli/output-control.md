# Output control

These flags control which annotation fields appear in the output. `--everything` switches on Ensembl VEP's `--everything` set plus `--allele_number` and `--hgvsg`.

## Summary

| Flag | Default | What it does |
| --- | --- | --- |
| [`--everything`, `-e`](#opt_everything) | off | Switches on most of the flags in this table and the allele frequency flags |
| [`--symbol`](#opt_symbol) | off | Gene symbol, for example `BRCA1` |
| [`--biotype`](#opt_biotype) | off | Transcript biotype |
| [`--canonical`](#opt_canonical) | off | Flags canonical transcripts |
| [`--mane`](#opt_mane) | off | MANE Select and MANE Plus Clinical flags |
| [`--mane_select`](#opt_mane_select) | off | MANE Select flag only |
| [`--tsl`](#opt_tsl) | off | Transcript Support Level |
| [`--appris`](#opt_appris) | off | APPRIS isoform annotation |
| [`--ccds`](#opt_ccds) | off | CCDS transcript ID |
| [`--protein`](#opt_protein) | off | Ensembl protein ID |
| [`--protein_version`](#opt_protein_version) | off | Appends the translation version to the Ensembl protein ID |
| [`--uniprot`](#opt_uniprot) | off | UniProt cross-references `SWISSPROT` and `TREMBL` |
| [`--xref_refseq`](#opt_xref_refseq) | off | RefSeq transcript cross-references |
| [`--hgvs`](#opt_hgvs) | off | HGVS coding and protein notations |
| [`--hgvsg`](#opt_hgvsg) | off | Adds the `HGVSg` column, which stays empty |
| [`--sift`](#opt_sift) | none | SIFT predictions: `p`, `s` or `b` |
| [`--polyphen`](#opt_polyphen) | none | PolyPhen predictions: `p`, `s` or `b` |
| [`--numbers`](#opt_numbers) | off | Exon and intron numbers |
| [`--domains`](#opt_domains) | off | Adds the `DOMAINS` column, which is never populated |
| [`--regulatory`](#opt_regulatory) | off | No effect; regulatory-feature annotation is not implemented, and the flag warns |
| [`--variant_class`](#opt_variant_class) | off | Variant class, a Sequence Ontology term |
| [`--check_existing`](#opt_check_existing) | off | Checks for co-located known variants |
| [`--coding_only`](#opt_coding_only) | off | No effect |
| [`--no_intergenic`](#opt_no_intergenic) | off | No effect |

## Accepted for compatibility

`--coding_only` and `--no_intergenic` are accepted and have no effect. Two regulatory flags, `--regulatory_gff` and `--extended_promoters`, are refused with an error and are described under [Refused flags](compatibility.md#refused-flags).

| Flag | What Ensembl VEP does with it |
| --- | --- |
| `--hgvsg` | Adds genomic HGVS notation computed from the input coordinates |
| `--regulatory` | Reports overlaps with regulatory regions and transcription factor binding sites |
| `--coding_only` | Keeps only the consequences in the coding regions of transcripts |
| `--no_intergenic` | Drops intergenic consequences from the output |

## Options

### `--everything`, `-e` {#opt_everything}

**`--everything`**, boolean, default off. Ensembl VEP: same flag and default.

Switches on every flag in the summary table except `--xref_refseq`, `--coding_only` and `--no_intergenic`, plus every [allele frequency flag](allele-frequencies.md), `--mirna`, `--allele_number`, `--sift b` and `--polyphen b`. It implies `--regulatory`, so the not-implemented warning is logged on every `--everything` run.

Ensembl VEP's `--everything` does not switch on `--allele_number` or `--hgvsg`, so the `ALLELE_NUM` and `HGVSg` header keys differ there.

### `--symbol` {#opt_symbol}

**`--symbol`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the gene symbol, for example `BRCA1`.

### `--biotype` {#opt_biotype}

**`--biotype`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the transcript biotype.

### `--canonical` {#opt_canonical}

**`--canonical`**, boolean, default off. Ensembl VEP: same flag and default.

Flags canonical transcripts.

### `--mane` {#opt_mane}

**`--mane`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the MANE Select and MANE Plus Clinical flags.

### `--mane_select` {#opt_mane_select}

**`--mane_select`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the MANE Select flag only.

### `--tsl` {#opt_tsl}

**`--tsl`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the Transcript Support Level.

### `--appris` {#opt_appris}

**`--appris`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the APPRIS isoform annotation.

### `--ccds` {#opt_ccds}

**`--ccds`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the CCDS transcript ID.

### `--protein` {#opt_protein}

**`--protein`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the Ensembl protein ID.

### `--protein_version` {#opt_protein_version}

**`--protein_version`**, boolean, default off. Ensembl VEP: same flag and default.

Appends the translation version to the Ensembl protein ID (`ENSP00000286808.3`), unless it already carries one. Prints nothing without `--protein`. This is Ensembl VEP 116's flag.

### `--uniprot` {#opt_uniprot}

**`--uniprot`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the UniProt cross-references `SWISSPROT` and `TREMBL`. The `UNIPARC` and `UNIPROT_ISOFORM` columns are added but never computed, so they stay empty.

### `--xref_refseq` {#opt_xref_refseq}

**`--xref_refseq`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the RefSeq transcript cross-references.

### `--hgvs` {#opt_hgvs}

**`--hgvs`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the HGVS coding and protein notations, and `HGVS_OFFSET` when an insertion or deletion was shifted 3' along the transcript for them.

As in Ensembl VEP 116, a deletion that runs over either end of the transcript is described over the bases inside it (the first line below). Also as in Ensembl VEP 116, a frameshift whose first changed residue is the stop codon reads as an extension, the count excluding the replaced stop (the second line below, not `p.Ter124IlefsTer15`).

```
ENST00000290239.6:c.-49_77+1248del
p.Ter124IleextTer14
```

### `--hgvsg` {#opt_hgvsg}

**`--hgvsg`**, boolean, default off. Ensembl VEP: accepted, no effect.

HGVSg is not computed: the `HGVSg` column is added but always empty.

### `--sift` {#opt_sift}

**`--sift <p|s|b>`**, string, no default. Ensembl VEP: same flag and default.

SIFT predictions: `p` the prediction, `s` the score, `b` both.

### `--polyphen` {#opt_polyphen}

**`--polyphen <p|s|b>`**, string, no default. Ensembl VEP: same flag and default.

PolyPhen predictions: `p` the prediction, `s` the score, `b` both.

### `--numbers` {#opt_numbers}

**`--numbers`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the exon and intron numbers.

### `--domains` {#opt_domains}

**`--domains`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the `DOMAINS` column, which is never populated; see [Accepted without effect](compatibility.md#accepted-without-effect).

### `--regulatory` {#opt_regulatory}

**`--regulatory`**, boolean, default off. Ensembl VEP: accepted, no effect.

Regulatory-feature annotation is not implemented; the flag is accepted and warns. The six regulatory consequence terms are listed under [Consequence terms never assigned](compatibility.md#consequence-terms-never-assigned).

### `--variant_class` {#opt_variant_class}

**`--variant_class`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the variant class, a Sequence Ontology term.

### `--check_existing` {#opt_check_existing}

**`--check_existing`**, boolean, default off. Ensembl VEP: same flag and default.

Checks for co-located known variants.

### `--coding_only` {#opt_coding_only}

**`--coding_only`**, boolean, default off. Ensembl VEP: accepted, no effect.

Not implemented.

### `--no_intergenic` {#opt_no_intergenic}

**`--no_intergenic`**, boolean, default off. Ensembl VEP: accepted, no effect.

Not implemented.
