# ClinVar records whose significance holds several classifications

The inputs from which `vep-cache-builder` builds the variation file of two known variants of
chromosome 21, `rs74315447` (34370639) and `rs1805128` (34449382), on GRCh38. They are cut from
Ensembl's files of release 116:

- `homo_sapiens-chr21.vcf`: the header of the release's variation VCF of chromosome 21 and its two
  lines of the two variants. The builder reads the file as BGZF, so a test or a script compresses
  it first.
- `phenotype_feature.txt`: the 18 ClinVar rows of the two variants in the variation database's
  `phenotype_feature` table, in the table's order.
- `phenotype_feature_attrib.txt`: the `clinvar_clin_sig`, `risk_allele` and `clinvar_ref_allele`
  rows of those 18 records, the three attributes a build reads for a germline record.
- `source.txt`, `seq_region.txt`, `attrib_type.txt`: the rows of those tables that the rows above
  name.
- `phenotype.txt`: empty. A build reads a phenotype's description for a somatic record alone.

Three of the 18 records carry a `clinvar_clin_sig` value of two classifications joined by a
semicolon and a blank: `uncertain significance; risk factor` once for `rs74315447`, and for
`rs1805128` `benign/likely benign; other` once and `likely benign; risk factor` twice.

`tests/clinvar_significance_of_several_classifications_on_a_native_cache.rs` builds the variation
file from these inputs, and the reproduction of the documented difference
`clin_sig_allele_value_cut_at_semicolon` builds a cache from them.
