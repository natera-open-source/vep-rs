# SV Concordance Report

## Overall

- **Reference**: Ensembl VEP 116.2; masks applied: `filter_cnv_tr_expansion_divergences`, `filter_registry_swap_pairs`, `filter_breakend_mate_local_read_pairs`, `filter_reference_skipped_records`, `filter_bnd_own_rows_lost`, `filter_bnd_synonym_mate_rows`, `filter_bnd_bare_reference_allele_rows`, `filter_gvcf_variant_sites`, `filter_dot_breakend_pairs`, `filter_missing_alt_annotated`, `filter_non_integer_end_rows`, `filter_cnv_tr_empty_expansion`, `filter_star_alone_svtype_annotated`, `filter_non_numeric_pos_numified`, `filter_cnv_tr_rn_zero_allele_dropped`
- **116.2 classes**: `reference_skipped_record` 6 vep-rs rows on 4 records set aside; `bnd_own_chromosome_rows_lost` 0 set aside; `bnd_synonym_mate_rows_lost` 0 set aside; `fasta_named_slice_intergenic` 0 records counted, none set aside; `bnd_bare_reference_allele_row` 0 reference-allele rows set aside; `bnd_mate_32bit_buffer_intergenic` 0 records counted, none set aside; `gvcf_variant_site_dropped` 0 vep-rs rows on 0 sites set aside; `dot_breakend_without_svtype_pair` 0 reference and 0 vep-rs rows set aside; `missing_alt_annotated` 0 reference rows set aside; `non_integer_end_annotated` 0 reference and 0 vep-rs rows set aside; `cnv_tr_empty_expansion_dropped` 0 vep-rs rows set aside; `star_alone_svtype_annotated` 0 reference rows set aside; `non_numeric_pos_numified` 0 reference rows set aside; `cnv_tr_rn_zero_allele_dropped` 0 vep-rs rows set aside; `scored_engine_skip_route_residue` 0 scored-engine rows at Locations the reference wrote nothing for counted, none set aside
- **F1 (adjusted)**: 0.7692  *(excluded 6 Rust-extra + 0 Perl-extra transcript divergences)*
- **F1 (raw)**: 0.5263
- **Total input variants**: 8
- **Perl tuples**: 6
- **Rust tuples**: 13
- **Intersection**: 5
- **Only Perl**: 1
- **Only Rust**: 8

## Per-File Metrics

| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |
|------|-------|-------------|-------------|-----------|--------|--------|
| sv | 8 | 6 | 13 | 5 | 0.5263 | 0.7692 |

## Per-Variant-Type Metrics

| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |
|--------------|-------|-------------|-------------|-----------|-----------|------|
| BND | 2 | 4 | 5 | 0 | 1 | 0.8889 ** |
| Complex_SV | 2 | 0 | 4 | 0 | 4 | 0.0000 |
| NON_REF | 1 | 0 | 1 | 0 | 1 | 0.0000 |
| Symbolic_DEL | 3 | 2 | 3 | 1 | 2 | 0.4000 ** |

