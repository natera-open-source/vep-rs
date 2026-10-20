# Records `--dont_skip` cannot keep

A difference in which records are written under `--dont_skip`, on input shapes no measured dataset has; nothing is set aside.

`--dont_skip` returns a record whatever `validate_vf` returned for it (`Parser.pm` line 189). Four refused shapes have no record vep-rs can hold, and vep-rs skips them under the flag too, with Ensembl VEP's warning where it writes one:

- a coordinate that is not a number. Ensembl VEP writes the text as the record's position, `21:abc`, and annotates it as `intergenic_variant`, or ends the run with an exception when the record reaches a transcript. A VCF POS or an Ensembl-format start or end may also be a number vep-rs holds no coordinate for: negative (class 21, [Defect N](defect-n.md#21-a-non-numeric-pos-numified-by-the-structural-paths-increment)), fractional or in exponent notation (class 18, [Defect K](defect-k.md#18-a-non-integer-end-svlen-or-pos-text-annotated-and-printed)). Ensembl VEP reads such a number and annotates the record as written; vep-rs skips it with the parser's own refusal, with or without the flag;
- a start past the end plus one. Ensembl VEP's run ends with `start arg must be less than or equal to end arg + 1` on the first transcript the record reaches;
- a line with too few columns to name a variant, or one the parser refuses for its shape. vep-rs skips each such line with the parser's own text, and reads an Ensembl-format, region or HGVS line opening with `#` as a comment, skipped without a word. Ensembl VEP's reading depends on the format:
  - a VCF line of fewer than five columns is carried as a feature with no allele, its missing ALT an empty string, and written in no format but the VCF line and the JSON object;
  - an Ensembl-format line is read as a record of whatever columns are there, or dies;
  - a region line its regex does not match is passed over without a word, and when it opens a buffer fill it ends the run with the records after it unread (`Parser/Region.pm` line 171, `Parser.pm` lines 186-187, `Runner.pm` line 424);
  - an HGVS or SPDI notation its parser cannot read is warned without a skipped-variant record (`Parser/HGVS.pm` line 159, `Parser/SPDI.pm` line 183);
- a region, SPDI or HGVS transcript record whose span exceeds `--max_sv_size`, refused at `validate_vf`'s size test, where Ensembl VEP carries the feature without consequences (`Parser.pm` lines 491-501). The HGVS transcript record is warned about, whatever `--chr` lists, since its line names a transcript or a gene and not the chromosome the list is read against.

A record that fails one of the other checks is kept and annotated as it stands in both engines. Those checks are a chromosome the cache lacks, an invalid allele string, a reference allele `--check_ref` rejects, a reference allele whose length is not the span's, an insertion whose coordinates are not start = end + 1, and a chromosome `--chr` leaves out.
