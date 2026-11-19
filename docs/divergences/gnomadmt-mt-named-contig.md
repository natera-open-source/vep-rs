# A gnomADMt data file whose contig is named `MT`

**The defect in one sentence.** Ensembl VEP's gnomADMt plugin rewrites every mitochondrial record's chromosome to `M` before it queries the data file, and the base class it reads the file through tries that name and `chrM` alone. Against a file whose contig is named `MT` or `chrMT` the plugin finds nothing and writes no field and no warning.

**Why this is a defect.** The plugin's documentation offers the file's fields for a record on the mitochondrial chromosome, and a tabix-indexed sites file states its own contig name in its records and its header. The rewrite exists because the gnomAD file names the contig `chrM`, which the base class reaches from `M` by adding `chr`; it is applied to the record and never to the file's names. A file written with the Ensembl spelling `MT`, the name the cache and the output use, is therefore unreachable, and a silent empty result is indistinguishable from a site the file lacks.

**What vep-rs does instead.** vep-rs queries the file under the mitochondrial chromosome's four spellings, `MT`, `chrMT`, `M` and `chrM`, and writes the record the file holds at the position, whichever name the file uses.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-variation` at release/116):

- `gnomADMt.pm`, lines 256-261, `run`: a `chr` prefix is stripped from the record's chromosome and `MT` becomes `M`; line 267 reads the file with that name.
- `Bio/EnsEMBL/Variation/Utils/BaseVepTabixPlugin.pm`, lines 693-707, `_get_source_chr_name`: a name the file lacks is tried without its `chr` prefix when it has one, else with `chr` added, and nothing else. Lines 596-602 then query the mapped name, so `M` against a file naming `MT` queries `M` and yields no row.

**What correcting it would change.** A record on `MT` would carry the file's values against a file named `MT` or `chrMT` as it does against the gnomAD file, which names the contig `chrM`.

**The class it produces.** `gnomadmt_mt_named_contig`: a `gnomAD_<field>` (JSON `gnomad_<field>`) that vep-rs alone fills on a mitochondrial record whose data file names the contig `MT` or `chrMT`. The class is documented only: the committed slice and every comparison case name the contig `chrM`, as the gnomAD file does.

## 1. The committed slice re-headered to `MT`

**In short.** The `GRCh38-plugin-gnomadmt` corpus's gnomAD slice, with its `##contig` line and the first column of every record rewritten from `chrM` to `MT`, is run on the record `MT 3 . T C` with the default fields. Ensembl VEP writes every row of the record with `IMPACT=MODIFIER` and no `gnomAD_` field, and nothing on standard error. Against the slice as committed the same run writes `gnomAD_hap_AC_het=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0` and the four other haplogroup fields. vep-rs writes those five fields against both files.

vep-rs writes those five fields on every row against both files, byte for byte the rows the corpus's expected output holds for `MT_3_T/C`. The commands of [Reproducing these examples](reproducing.md#reproducing-these-examples) give the engines and their caches; the slice is `tests/golden/116/GRCh38-plugin-gnomadmt/data/gnomad.v3.1.chrM.slice.vcf.bgz`, its header and first column rewritten before `bgzip` and `tabix`.
