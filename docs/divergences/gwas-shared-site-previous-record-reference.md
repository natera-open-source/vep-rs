# A GWAS record's match that depends on the record before it at the same site

**The defect in one sentence.** Ensembl VEP's GWAS plugin lends a catalog row the reference allele of the record that first reads it, and the lent allele stays in the rows it caches for the site. The next record at the same site, with another reference string as given, is matched against the first record's reference and gains or loses its `GWAS_` fields with the order of the input.

**Why this is a defect.** The catalog states a risk allele and no reference base, so the plugin takes the reference from the record being annotated, and only while the row's reference is still empty or `N`. The rows it reads for a site come from a one-entry cache keyed by the site, and the lent allele is written into the cached rows themselves, so the second record at the site finds a reference already filled and keeps it.

Two records at one site whose reference strings as given differ, an Ensembl-format record on the reverse strand beside one on the forward strand, are therefore matched against one reference, the first record's, and the result changes when the two are swapped.

**What vep-rs does instead.** vep-rs lends each record its own reference allele and keeps no state between records, so a record's `GWAS_` fields are the same wherever it stands in the input and equal Ensembl VEP's for the same record alone at its site.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-variation` at release/116):

- `GWAS.pm`, line 176, `run`: the rows are read through `get_data` for the record's chromosome, `start - 2` and end. Line 180 writes the record's reference allele string into each row whose reference is empty or `N`, and lines 182-195 match the record's reference and alleles against the row's and return the first row matched.
- `Bio/EnsEMBL/Variation/Utils/BaseVepTabixPlugin.pm`, line 114: the results cache holds one entry; line 361: its key is the chromosome, start and end queried; lines 366-368: a query with the same key returns the cached rows themselves, the ones lines 558-569 stored.
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 1176-1196, `get_matched_variant_alleles`: the record's reference and alleles are reverse-complemented when its strand differs from the row's, which defaults to forward, while the row's lent reference is compared as given. A reverse-strand record alone at its site therefore matches nothing in either engine.

**What correcting it would change.** Every record would be matched against its own reference, as it is when it stands alone at its site, and the order of the input would change no `GWAS_` field.

**The class it produces.** `gwas_shared_site_previous_record_reference`: a `GWAS_` field (the JSON `gwas_` keys) that differs in either direction on a record sharing its chromosome, start and end with the record before it in the input, the two records' reference alleles as given differing. The shape is a reverse-strand record of the default input format beside a forward-strand one, or a VCF record whose REF contradicts its predecessor's. The class is documented only: no committed comparison carries such a pair, since every committed GWAS input is a VCF whose records carry the genome's reference.

## 1. One site, two records, two orders

**In short.** The Ensembl-format records `21 34300618 34300618 T/C - rev` and `21 34300618 34300618 A/G + fwd` name the same substitution, `A>G` at `21:34300618`, where the `GRCh38-plugin-gwas` corpus's catalog slice has `rs2834439-G` of study `GCST005275`. Given alone, the forward record carries `GWAS_study=GCST005275;GWAS_risk_allele=G` and its six other `GWAS_` fields on each of its four transcript rows in both engines, and the reverse record carries no `GWAS_` field in either. Given together, Ensembl VEP writes the fields on none of the eight rows when the reverse record stands first and on all eight when the forward one does; vep-rs writes them on the forward record's four rows in either order.

Given reverse then forward, Ensembl VEP writes no `GWAS_` field on any of the eight rows: the forward record is matched against the lent `T`. Given forward then reverse, it writes the eight fields on all eight rows: the reverse record is matched against the lent `A`. vep-rs writes the forward record's four rows with the fields and the reverse record's four without them in either order.

**Try it yourself.** Run the two records as one Ensembl-format input in each order against the `GRCh38-plugin-gwas` corpus's cache and catalog slice (`tests/golden/116/GRCh38-plugin-gwas/`), with `--format ensembl --plugin GWAS,file=<slice>`, and compare the `GWAS_study` field of the eight rows. The commands of [Reproducing these examples](reproducing.md#reproducing-these-examples) give the engines and their caches.
