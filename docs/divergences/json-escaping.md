# Unescaped JSON under `--output_format json`

A difference in the text of one field, not in any tuple: the comparators score consequence sets and never see it, and nothing is set aside.

Ensembl VEP percent-encodes the `=` of a synonymous HGVSp notation as `%3D` in every output unless `--no_escape` is set (`OutputFactory.pm`, `TranscriptVariationAllele_to_output_hash`: `$hgvs_p =~ s/\=/\%3D/g if $hgvs_p && !$self->{no_escape}`). It sets `--no_escape` through the `rest` option set that `--json` switches on (`Config.pm` lines 484-498). `--output_format json` is not in that set, so Ensembl VEP's JSON written through it carries `%3D` while its JSON written through `--json` carries `=`. vep-rs reads the resolved output format: `--json` and `--output_format json` both imply `--no_escape`, so its JSON is unescaped either way. A text output stays escaped unless `--no_escape` is given.

The record is `rs148468217`, `21 34655545 C T` on GRCh37, synonymous on ENST00000290200 (`tests/golden/116/GRCh37-output-flags/variants.vcf`, with that corpus's cache and the FASTA its manifest names, the GRCh37 chromosome 21 FASTA of the sibling `GRCh37-hgvs` corpus). Both engines run with `--hgvs --fasta <reference.fa> --output_format json`, and the `hgvsp` of the ENST00000290200 entry is:

```
Ensembl VEP 116.2   "hgvsp":"ENSP00000290200.2:p.Asp215%3D"
vep-rs              "hgvsp":"ENSP00000290200.2:p.Asp215="
```

With `--json` in place of `--output_format json` both engines write `p.Asp215=`. In the default, tab and VCF outputs both write `p.Asp215%3D` unless `--no_escape` is given, under which both write `p.Asp215=`. The `no-escape` flag set of the same corpus compares those. The golden corpora pass `--json` to both engines, so no JSON comparison reads the difference.
