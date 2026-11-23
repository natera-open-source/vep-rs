# No GWAS annotation at position 1 or 2 of a chromosome

**The defect in one sentence.** Ensembl VEP's GWAS plugin reads the catalog rows from two bases before a record's start, so a record at position 1 or 2 asks for a window that starts at 0 or -1 and the read dies. The record carries no `GWAS_` field and a `Plugin 'GWAS' went wrong` warning.

**Why this is a defect.** The window exists so that a row written one or two bases before an indel's start is still read; a record at the first two positions of a chromosome has fewer than two bases before it, and the correct window is the bases that exist. The plugin passes the arithmetic result on: a start of 0 is false in Perl and refused as no start at all, and a start of -1 is an unparsable region for the tabix reader. Either way the exception is caught by the plugin runner and reported, and a catalog row at position 1 or 2 can annotate nothing.

**What vep-rs does instead.** vep-rs clamps the window at position 1 and matches the record against the rows it holds, so a catalog row at position 1 or 2 annotates its record as a row anywhere else does.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-variation` at release/116, `ensembl-vep` at release/116.3):

- `GWAS.pm`, line 176, `run`: `get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})`.
- `Bio/EnsEMBL/Variation/Utils/BaseVepTabixPlugin.pm`, line 354, `get_data`: `die("ERROR: No start specified\n") unless $s`, which a start of 0 fails; lines 596-602, `_get_data_hts`: the region is formatted `chr:start-end`, and `21:-1-1` is refused by the reader with `You must specify a region in the format chr, chr:start or chr:start-end`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2252-2254: the exception from the plugin's `run` becomes the warning `Plugin 'GWAS' went wrong: <reason>` and the row is written without the plugin's fields.

**What correcting it would change.** A record at position 1 or 2 would be matched against the rows from position 1 to its end, as vep-rs matches it, and the warning would go.

**The class it produces.** `gwas_window_below_position_one`: a `GWAS_` field (the JSON `gwas_` keys) that vep-rs alone fills on a record whose start is 1 or 2, and in `--warnings` mode the reference's `Plugin 'GWAS' went wrong` lines on such records. The class is documented only. The full associations download has no row at position 1 or 2 of any chromosome. The committed records whose start is 1 or 2, a mitochondrial deletion at position 1 of the fly corpus and the deletion and insertion written at position 0 of the VCF edge-case corpus, run under no GWAS flag set.

## 1. Three records at the start of chromosome 21

**In short.** With the `GRCh38-plugin-gwas` corpus's catalog slice carrying one row moved to `CHR_POS` 2 (`rs2823265-G`, study `GCST90002367`), the records `21 1 pos1 C G`, `21 2 pos2 C G` and `21 3 pos3 C G` are run. Ensembl VEP writes no `GWAS_` field on any of the three. Its warnings file carries `Plugin 'GWAS' went wrong: You must specify a region in the format chr, chr:start or chr:start-end` for the record at 1 and `Plugin 'GWAS' went wrong: ERROR: No start specified` for the record at 2, whose processed file holds the row. vep-rs writes the row's `GWAS_` fields on the record at 2 and no field on the other two.

vep-rs writes `GWAS_accessions=http://purl.obolibrary.org/obo/OBA_2045260;GWAS_associated_gene=;GWAS_beta_coef=;GWAS_odds_ratio=;GWAS_p_value=1E-37;GWAS_pmid=32888493;GWAS_risk_allele=G;GWAS_study=GCST90002367` on the record at 2 and no field on the other two, whose positions the row does not name. The commands of [Reproducing these examples](reproducing.md#reproducing-these-examples) give the engines and their caches; the catalog slice is `tests/golden/116/GRCh38-plugin-gwas/plugin_data/gwas_catalog_slice.tsv`.
