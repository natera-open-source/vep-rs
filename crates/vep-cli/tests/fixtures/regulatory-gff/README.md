# Regulatory GFF test files

The three files are the records of chromosome 21 whose start lies between
5,000,000 and 11,000,000 of the Ensembl release 116 GRCh38 regulation
annotation, cut unchanged from the files the Ensembl site distributes under
`regulation/homo_sapiens/GRCh38/annotation/`:

| File                          | Source file                                           | Records |
| ----------------------------- | ----------------------------------------------------- | ------- |
| `regulatory_features.gff3.gz` | `Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz` | 516     |
| `motif_features.gff3.gz`      | `Homo_sapiens.GRCh38.motif_features.v116.gff3.gz`      | 2,200   |
| `EMARs.gff.gz`                | `Homo_sapiens.GRCh38.EMARs.v116.gff.gz`                | 279     |

Each is bgzip-compressed with its `tabix -p gff` index beside it, the form
`--regulatory_gff` reads. The regulatory file carries 47 promoters with
`extended_start` and `extended_end`, 193 CTCF binding sites on a strand, 265
enhancers and 11 open chromatin regions; the motif file names 50 binding
matrices by `binding_matrix_id`, each motif with its `transcription_factor`
list and a strand; the EMAR file holds `EMAR` records with an `ID` alone.

Regenerate the window from the source files with

```bash
for f in regulatory_features.v116.gff3 motif_features.v116.gff3 EMARs.v116.gff; do
  out=$(echo "$f" | sed 's/\.v116//')
  zcat "Homo_sapiens.GRCh38.$f.gz" \
    | awk -F'\t' '/^#/ || ($1=="21" && $4>=5000000 && $4<11000000)' > "$out"
  bgzip -f "$out" && tabix -f -p gff "$out.gz"
done
```
