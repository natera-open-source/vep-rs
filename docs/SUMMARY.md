# Summary

[vep-rs](index.md)

# Getting started

- [Installation](getting-started/installation.md)
- [Quick start](getting-started/quick-start.md)
- [Cache setup](cache-setup.md)
- [Coming from Ensembl VEP](getting-started/coming-from-ensembl-vep.md)

# Reference

- [CLI reference](cli-reference.md)
  - [Input and output](cli/input-output.md)
  - [Cache, species and runtime](cli/cache-species-runtime.md)
  - [Output control](cli/output-control.md)
  - [Allele frequencies and co-located variants](cli/allele-frequencies.md)
  - [Filtering, picking and transcript sets](cli/filtering-picking.md)
  - [Plugins](cli/plugins.md)
  - [Statistics, verbosity and miscellaneous](cli/statistics-verbosity-misc.md)
  - [Compatibility with Ensembl VEP](cli/compatibility.md)
- [Output formats](output-formats.md)
  - [VEP default and tab](output/vep-and-tab.md)
  - [VCF](output/vcf.md)
  - [JSON](output/json.md)
  - [Parquet](output/parquet.md)
  - [How structural variants are written](output/structural-variants.md)
- [Plugins](plugins.md)
  - [Annotation-store plugins](plugins/annotation-store.md)
  - [LoFTEE](plugins/loftee.md)
  - [Flat-file plugins](plugins/flat-file.md)
  - [SpliceAI](plugins/spliceai.md)
  - [Preparing data files](plugins/data-files.md)

# Concordance with Ensembl VEP

- [Intended divergences](intended-divergences.md)
  - [Defect A: a later differing region erases the splice-region term](divergences/defect-a.md)
  - [Defect B: `start_lost` and `start_retained_variant` are both written for one allele](divergences/defect-b.md)
  - [Defect C: a record marked skipped is dropped before annotation](divergences/defect-c.md)
  - [Defect P: the tests after the chromosome-aware overlap compare chromosome names](divergences/defect-p.md)
  - [Defect E: a breakend's mate-side row evaluated at the record's own coordinate](divergences/defect-e.md)
  - [A representation difference: `<CNV:TR>`](divergences/cnv-tr.md)
  - [Defect F: a breakpoint record without a mate gives its reference allele a row](divergences/defect-f.md)
  - [Defect G: a coordinate past the 32-bit bound turns the buffer intergenic](divergences/defect-g.md)
  - [Defect H: a gVCF variant site is dropped with its real allele](divergences/defect-h.md)
  - [Defect I: a dot-notation breakend without `SVTYPE` is written as a pair](divergences/defect-i.md)
  - [Defect J: the missing ALT value is annotated as an allele](divergences/defect-j.md)
  - [Defect K: a non-integer coordinate text is annotated](divergences/defect-k.md)
  - [Defect L: an empty or unchanged `<CNV:TR>` expansion is dropped without a row](divergences/defect-l.md)
  - [Defect M: `*` as the only ALT is annotated from `SVTYPE`](divergences/defect-m.md)
  - [Defect N: a POS that is not a number is numified on the structural path](divergences/defect-n.md)
  - [Defect O: a `<CNV:TR>` allele with `RN=0` drops a real allele of the record](divergences/defect-o.md)
  - [Reproducing the examples](divergences/reproducing.md)
- [Divergence classes against 116.2](concordance-provenance/116.2-divergence-classes.md)

# Benchmarks

- [Current release](benchmarks/current-release.md)
- [The figures the paper reports](published-figures.md)

# Developers

- [Writing a plugin](developers/writing-a-plugin.md)
- [Contributing](contributing.md)

# Project

- [Citing vep-rs](project/citing.md)
- [Changelog](changelog.md)
