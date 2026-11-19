# `vep recode`: the warnings file, a repeated input's warnings and the fields a flag switches on

**In short.** Three readings of Ensembl VEP's `variant_recoder` are corrected in [`vep recode`](../cli/recode.md) rather than reproduced. `--warning_file` writes the warnings to the file it names, where `variant_recoder` accepts the flag and writes no file. An input text repeated in two buffer fills that each raise a warning carries both warnings, where `variant_recoder` keeps the last fill's alone ([below](#the-warnings-of-a-repeated-input)). A field flag switches on what computes its field, so `--ga4gh_vrs` switches `spdi` on, where `variant_recoder` appends the key alone. The four committed recoder cases name `--ga4gh_vrs` under the default field list, where the two agree, and none names `--warning_file` or repeats a line (tier 3 of [Intended divergences from Ensembl VEP](../intended-divergences.md)).

## The warnings file

The warnings raised while an input is read or annotated are carried in the output, under `warnings` on the next object written, and `--warning_file` writes them to a file as well. `variant_recoder` accepts `--warning_file` too, "file to write warnings to" in its own option list (`variant_recoder`, line 78), but writes nothing to it. Its initialisation redirects every warning into memory before one is raised (`Bio/EnsEMBL/VEP/VariantRecoder.pm`, lines 171-181; `Runner.pm`, lines 295-300), and the named file is opened only when no redirection exists (`BaseVEP.pm`, lines 878-900).

On a two-record input whose first record names a chromosome the cache lacks, both engines write one object carrying their own skip warning under `warnings` beside the second record's allele; `vep recode --warning_file w.txt` writes that warning to `w.txt`, where `variant_recoder` leaves no file.

## The warnings of a repeated input

**The defect in one sentence.** `variant_recoder` keys its results by the input text and assigns each buffer fill's warning list to that text's object, so an input given twice, each occurrence in a fill that raises a warning, keeps the second fill's warnings and loses the first's.

**Why this is a defect.** Both warnings were raised and belong to the output: the recoder's own rule puts a fill's warnings on the object it writes. The repeated text shares one object across fills, so the second assignment overwrites the first's list instead of adding to it; a repeated input with a warning in one fill alone, or two distinct inputs each with a warning, keep every warning.

**What vep-rs does instead.** `vep recode` adds each fill's warnings to the object of the input it belongs to, so a repeated input's object lists every warning raised across its fills, in input order.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/VariantRecoder.pm`, line 468: the results hash is keyed by the input text, and `||=` reuses the object an earlier fill created for the same text. Line 485: `$results->{$line_id}->{warnings} = [...]`, a plain assignment of the current fill's list; line 486: the warning store is emptied for the next fill.
- `Bio/EnsEMBL/VEP/Runner.pm`, lines 295-300, `internalise_warnings`: each fill starts with an empty warning string, so `warnings` holds that fill's alone. `Bio/EnsEMBL/VEP/InputBuffer.pm`, lines 211-212, with `VariantRecoder.pm`, lines 96-109: the buffer holds one input by default, so two occurrences are two fills.

**What correcting it would change.** The assignment would append to the list the earlier fill stored, and a repeated input's object would carry every warning, as `vep recode` writes it.

**The class it produces.** `recode_warnings_every_buffer_fill`: the `warnings` list of a recoder output object whose input text occurs in two or more buffer fills that each raise a warning; `vep recode` lists every fill's, `variant_recoder` the last fill's alone, and the other keys of the object agree. The class is documented only: no committed recoder input repeats a line.

### 1. One identifier twice, each after an identifier the cache lacks

**In short.** The identifiers `rs_nope1`, `rs886039732`, `rs_nope2`, `rs886039732`, one per line, under `--fields hgvsg` over the `tests/golden/116/GRCh38-recoder` corpus's cache and FASTA. With the default buffer size of 1, each fill holds the known variant and the `No variant found with ID` warning of the line before it. `vep recode` writes one object: `hgvsg` `NC_000021.9:g.32694272_32694276del`, `input` `rs886039732`, and `warnings` naming `rs_nope1` then `rs_nope2`. `variant_recoder` reads identifiers from the Ensembl database alone. The same shape as VCF lines, `21 32694271 . ATTTGT A` twice with a line naming a chromosome the cache lacks before each, run through its module offline, gives one object whose `warnings` names the second skipped line alone; every other key equals `vep recode`'s.

## The fields a flag switches on

`--ga4gh_vrs` adds a GA4GH VRS Allele object for each `spdi`, which it switches on. In `variant_recoder` the flag appends the key without switching `spdi` on (`VariantRecoder.pm`, lines 120-122 and 148-150; `--mane_select` does switch its three fields on, lines 141-146), so `--fields hgvsg --ga4gh_vrs` writes `hgvsg` alone there and `hgvsg`, `spdi` and `ga4gh_vrs` here. Likewise `vcf_string`, `var_synonyms` or `mane_select` named in `--fields` without its flag is computed here and not there, where those three are attached under their flags alone (lines 131-146); in `variant_recoder` a field name switches on the hgvs, spdi, id and `ga4gh_vrs` computations alone (lines 120-126).

## Try it yourself

The input is two HGVS notations, one per line. `vep recode` reads the `tests/golden/116/GRCh38-recoder` corpus's cache and FASTA with its `chr_synonyms.txt`, which names chromosome 21's `NC_000021.9` accession; `variant_recoder`, which has no cache or FASTA option, reads the Ensembl database (`--db_version 116`):

```
ZZ:g.100A>G
NC_000021.9:g.33553711C>A
```

| Run | `variant_recoder` | `vep recode` |
| --- | --- | --- |
| `--warning_file w.txt` | one object, `A` and `warnings`; no `w.txt` | the same object shape; `w.txt` holds the warning |
| `--fields hgvsg --ga4gh_vrs`, the second line alone | `A` holds `hgvsg` and `input` | `A` holds `hgvsg`, `input`, `spdi` and `ga4gh_vrs` |

Each engine's `warnings` entry is its own text for the unparsable first line: `Unable to parse HGVS notation 'ZZ:g.100A>G'` from `variant_recoder`, and from `vep recode` the same sentence after `line 1 skipped (ZZ:g.100A>G):` and before `The reference FASTA has no sequence named ZZ`. In the second run `vep recode`'s `spdi` is `NC_000021.9:33553710:C:A` and its `ga4gh_vrs` one Allele object on `refseq:NC_000021.9`.
