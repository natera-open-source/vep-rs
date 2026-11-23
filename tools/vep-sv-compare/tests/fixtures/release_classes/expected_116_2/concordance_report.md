# SV Concordance Report

## Overall

- **Reference**: Ensembl VEP 116.2; masks applied: `filter_cnv_tr_expansion_divergences`, `filter_registry_swap_pairs`, `filter_breakend_mate_local_read_pairs`, `filter_reference_skipped_records`, `filter_bnd_own_rows_lost`, `filter_bnd_synonym_mate_rows`, `filter_bnd_bare_reference_allele_rows`, `filter_gvcf_variant_sites`, `filter_dot_breakend_pairs`, `filter_missing_alt_annotated`, `filter_non_integer_end_rows`, `filter_cnv_tr_empty_expansion`, `filter_star_alone_svtype_annotated`, `filter_non_numeric_pos_numified`, `filter_cnv_tr_rn_zero_allele_dropped`
- **116.2 classes**: `reference_skipped_record` 3 vep-rs rows on 2 records set aside; `bnd_own_chromosome_rows_lost` 1 set aside; `bnd_synonym_mate_rows_lost` 3 set aside; `fasta_named_slice_intergenic` 1 records counted, none set aside; `bnd_bare_reference_allele_row` 2 reference-allele rows set aside; `bnd_mate_32bit_buffer_intergenic` 1 records counted, none set aside; `gvcf_variant_site_dropped` 2 vep-rs rows on 1 sites set aside; `dot_breakend_without_svtype_pair` 2 reference and 2 vep-rs rows set aside; `missing_alt_annotated` 1 reference rows set aside; `non_integer_end_annotated` 1 reference and 2 vep-rs rows set aside; `cnv_tr_empty_expansion_dropped` 2 vep-rs rows set aside; `star_alone_svtype_annotated` 2 reference rows set aside; `non_numeric_pos_numified` 2 reference rows set aside; `cnv_tr_rn_zero_allele_dropped` 1 vep-rs rows set aside; `scored_engine_skip_route_residue` 0 scored-engine rows at Locations the reference wrote nothing for counted, none set aside
- **F1 (adjusted)**: 0.7857  *(excluded 22 Rust-extra + 14 Perl-extra transcript divergences)*
- **F1 (raw)**: 0.3438
- **Total input variants**: 25
- **Perl tuples**: 27
- **Rust tuples**: 37
- **Intersection**: 11
- **Only Perl**: 16
- **Only Rust**: 26

SV adjusted F1 additionally excludes the `<CNV:TR>` literal-expansion class (3 of 9 pairs GRCH37, 03_tandem_repeats only). Raw F1 unaffected.

## Per-File Metrics

| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |
|------|-------|-------------|-------------|-----------|--------|--------|
| 01_breakends | 6 | 9 | 14 | 8 | 0.6957 | 0.9412 |
| 02_one_sided | 12 | 11 | 11 | 1 | 0.0909 | 0.4000 |
| 03_tandem_repeats | 6 | 5 | 10 | 1 | 0.1333 | 0.5000 |
| 04_small_variants | 1 | 2 | 2 | 1 | 0.5000 | 1.0000 |

## Per-Variant-Type Metrics

| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |
|--------------|-------|-------------|-------------|-----------|-----------|------|
| BND | 7 | 11 | 16 | 3 | 8 | 0.5926 ** |
| Multi_allelic | 2 | 2 | 5 | 2 | 5 | 0.0000 |
| RefOnly | 1 | 1 | 0 | 1 | 0 | 0.0000 |
| SNV | 1 | 1 | 0 | 1 | 0 | 0.0000 |
| Small_INS | 1 | 2 | 2 | 1 | 1 | 0.5000 ** |
| Spanning | 1 | 2 | 0 | 2 | 0 | 0.0000 |
| Symbolic_CNV | 4 | 2 | 6 | 2 | 6 | 0.0000 |
| Symbolic_DEL | 8 | 6 | 7 | 5 | 6 | 0.1538 ** |

