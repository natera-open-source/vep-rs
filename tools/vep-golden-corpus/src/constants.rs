// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The names and bounds every subcommand shares: the manifest's fixed values, the mask
//! classes a `masks` list may name, the keys `classify` writes, the columns of the two file
//! formats and the genotype and record classes.

pub const GENERATOR_VERSION: &str = "2";
pub const RESERVOIR: usize = 64;
pub const FAMILY_RESERVOIR: usize = 8;
pub const SEQUENCE_CLASSES: [&str; 4] = ["snv", "deletion", "insertion", "mnv_or_complex"];
pub const OUTPUT_FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

/// How a reference was generated: Ensembl VEP on its cache, connected to its database, or
/// reading its transcripts from the file given to `--gtf` or `--gff`.
pub const PERL_MODES: [&str; 4] = ["offline", "database", "gtf", "gff"];
pub const DEFAULT_INPUT: &str = "variants.vcf";

/// The species of a manifest that names none; another is written in Ensembl's production
/// spelling, the cache directory's name.
pub const DEFAULT_SPECIES: &str = "homo_sapiens";

/// The `masks` class under which the reference's regulatory rows are left out, and the
/// `Feature_type` values of those rows.
pub const REGULATORY_FEATURE_ROWS: &str = "regulatory_feature_rows";
pub const REGULATORY_FEATURE_TYPES: [&str; 2] = ["RegulatoryFeature", "MotifFeature"];

/// The `masks` class under which the set-valued columns agree as a multiset of their items
/// in any order, and those columns: Ensembl VEP prints DOMAINS in the cache's
/// protein_features order (BaseTranscriptVariation.pm:634-643), FREQS in hash order
/// (OutputFactory.pm:1152-1153) and the ZYG list of `--individual_zyg` in hash order
/// (OutputFactory.pm:983); the default format joins the items with `,` (Utils.pm:185-198).
pub const ORDER_ONLY: &str = "order_only";
pub const ORDER_ONLY_COLUMNS: [&str; 3] = ["DOMAINS", "FREQS", "ZYG"];

/// The classes the golden harness alone gives a meaning, in its warnings, header and column
/// comparisons: the reference's warnings about a plugin module vep-rs cannot load, the
/// `phlyoP<N>way` key Ensembl VEP's `--phyloP` set misspells (Config.pm:555), its two
/// descriptions of a custom source's PC field (AnnotationSourceAdaptor.pm:395-397), and the
/// SOURCE column and source header lines of a reference read from a GTF or GFF. The
/// harness's `MASK_CLASSES` names these four after `REGULATORY_FEATURE_ROWS` and
/// `ORDER_ONLY`, and nothing else.
pub const HARNESS_ONLY_MASK_CLASSES: [&str; 4] = [
    "unloaded_module_plugin",
    "conservation_short_name_spelling",
    "pc_description_text",
    "file_cache_writes_no_source_field",
];

/// The classes of the behaviours vep-rs corrects where Ensembl VEP 116.3 is wrong, each a
/// page or section of docs/intended-divergences.md: a corpus documents their rows under
/// `expected_divergence` in its three lists and names the class in `masks`. The golden
/// harness's `CORRECTED_BEHAVIOUR_CLASSES` names every class of this list that a committed
/// corpus names in `masks`. A class whose defect Ensembl VEP 116.3 or an earlier release
/// corrects is no class of this list.
pub const CORRECTED_BEHAVIOUR_CLASSES: [&str; 107] = [
    "start_cooccurrence_swap",
    "splice_family_swap",
    "gxf_reverse_strand_utr_sequence",
    "conservation_max_text_order",
    "gwas_haplotype_cross_product",
    "gwas_interaction_first_snp_allele_whitespace",
    "alphamissense_two_digit_transcript_version",
    "cadd_structural_allele_match",
    "splicevault_sibling_allele_rows",
    "splicevault_json_later_rows",
    "transcript_annotator_sibling_allele_rows",
    "transcript_annotator_insertion_window",
    "utrannotator_annotation_pair_order",
    "nearest_exon_jb_item_order",
    "utrannotator_gained_start_window",
    "utrannotator_consequence_precedence",
    "utrannotator_structural_allele_text",
    "utrannotator_lost_start_after_insertion",
    "utrannotator_gained_start_remade_in_place",
    "utrannotator_frameshift_at_stop_boundary",
    "utrannotator_lost_stop_written_before_stop_codon",
    "utrannotator_without_reference",
    "loeuf_zero_minimum",
    "maxentscan_overlapping_windows",
    "alphamissense_sibling_allele_rows",
    "revel_sibling_allele_rows",
    "revel_other_reference_amino_acid_rows",
    "gnomadmt_info_value_second_equals",
    "loftee_anc_allele_multiallelic_record",
    "loftee_first_alt_donor_filter",
    "loftee_gerp_zoom_level_summary",
    "loftee_gerp_walk_without_conservation_database",
    "loftee_nagnag_site_first_transcript_strand",
    "nearest_exon_jb_stretch_priming",
    "nearest_exon_jb_shifted_span",
    "nearest_exon_jb_id_less_exon_fold",
    "nearest_start_site_one_base_above",
    "downstream_reverse_strand_rotation",
    "spliceai_chopped_multiallelic",
    "spliceai_multi_entry_cell",
    "plugin_contig_mitochondrial_spelling",
    "tabix_allele_case_insensitive_match",
    "dbnsfp_mitochondrial_record_named_m",
    "dosage_sensitivity_regulatory_row_warnings",
    "hgvs_shift_shared_across_alleles",
    "hgvs_shift_sequence_edge_window",
    "refseq_edited_model_hgvsc_offset",
    "refseq_edited_model_shifted_insertion_bases",
    "refseq_edited_model_seam_insertion_hgvsc",
    "refseq_edited_model_unminimised_allele_hgvsc",
    "refseq_edited_model_hgvs_offset",
    "alternate_codon_past_the_coding_end",
    "organelle_peptide_paths_translation_table",
    "minimal_symbolic_alleles",
    "shift_3prime_shared_predicates",
    "shift_start_codon_test",
    "hgvs_shifted_insertion_clamped_end",
    "shift_intron_term_unshifted_span",
    "hgvs_intron_offset_abutting_insertion",
    "hgvs_shifted_indel_written_span_gate",
    "hgvs_shifted_deletion_clamped_span",
    "hgvs_clipped_reference_against_whole_alternate",
    "hgvs_added_bases_left_at_the_span_end",
    "hgvs_insertion_left_where_the_shift_stops",
    "hgvs_added_bases_written_outside_the_transcript",
    "hgvsp_stop_codon_protein_comparison",
    "hgvsp_stop_named_on_incomplete_coding_sequence",
    "hgvsp_frameshift_count_past_a_stop_read_as_a_residue",
    "hgvsp_frameshift_with_no_changed_residue",
    "hgvsp_frameshift_count_at_the_first_residue",
    "hgvsp_shifted_indel_written_span_gate",
    "hgvsp_shifted_indel_written_span_start_test",
    "hgvsp_start_codon_test_past_the_shifted_allele",
    "hgvsp_intron_end_insertion",
    "hgvsp_intron_start_insertion",
    "refseq_edited_model_deleted_span_codons",
    "refseq_edited_model_equal_span_insertion_term",
    "refseq_edited_model_intron_span_coding_reference",
    "refseq_edited_model_inserted_span_codons",
    "use_transcript_ref_cache_cdna",
    "refseq_failed_edit_model_hgvsc_offset",
    "hgnc_id_copy",
    "downstream_utr_frame_past_coding_end",
    "downstream_start_past_own_codon",
    "plugin_hgvs_shifted_cds_coordinates",
    "plugin_value_follows_stop_codon_term_pair",
    "plugin_value_follows_splice_term_past_the_record_span",
    "plugin_value_on_reference_skipped_record",
    "blosum62_hgvs_shift_shared_across_alleles",
    "structural_deletion_stop_lost_without_stop_codon",
    "stop_codon_completed_by_polyadenylation",
    "stop_codon_overlap_without_stop_codon",
    "stop_codon_behind_last_residue_without_stop_codon",
    "start_codon_overlap_without_start_codon",
    "stop_retained_without_stop_codon_test",
    "stop_codon_overlap_genomic_span",
    "stop_retained_three_base_stop_codon_test",
    "stop_codon_removed_without_stop_codon_test",
    "stop_codon_overlap_insertion_range",
    "stop_codon_settled_before_unknown_bases",
    "frameshift_inside_a_held_stop_codon",
    "splice_term_past_the_record_span",
    "start_codon_removed_without_five_prime_utr",
    "start_codon_kept_without_five_prime_utr",
    "cnv_tr_symbolic_loss_named_a_gain",
    "total_length_json_range_end",
    "allele_number_repeated_or_reference_equal_alt",
];

/// Every class a `masks` list, the manifest's or a flag set's, may name, in the order an
/// error lists them.
pub fn mask_classes() -> Vec<&'static str> {
    let mut classes = vec![REGULATORY_FEATURE_ROWS, ORDER_ONLY];
    classes.extend(HARNESS_ONLY_MASK_CLASSES);
    classes.extend(CORRECTED_BEHAVIOUR_CLASSES);
    classes
}

/// The keys `classify` writes into a run, in the order it writes them.
pub const CLASSIFY_KEYS: [&str; 7] = [
    "divergences",
    "vep_rs_only_tuples",
    "divergence_summary",
    "vep_rs_only_summary",
    "fields_compared",
    "field_divergences",
    "field_divergence_summary",
];

pub const DEFAULT_COLUMNS: [&str; 14] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "Extra",
];
pub const FIXED_COLUMNS: [&str; 8] = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
];

/// The genotype classes a kept sample shows, then the ones a record shows across its kept
/// samples; `GENOTYPE_CLASSES` is both, in this order.
pub const SAMPLE_CLASSES: [&str; 8] = [
    "het",
    "hom_alt",
    "hom_ref",
    "missing",
    "partially_missing",
    "haploid",
    "star",
    "multi_allelic",
];
pub const RECORD_CLASSES: [&str; 3] = ["all_missing", "no_carrier", "one_hom_ref_rest_missing"];
pub const GENOTYPE_CLASSES: [&str; 11] = [
    "het",
    "hom_alt",
    "hom_ref",
    "missing",
    "partially_missing",
    "haploid",
    "star",
    "multi_allelic",
    "all_missing",
    "no_carrier",
    "one_hom_ref_rest_missing",
];

/// The record-class key of `--record-classes`: a record's regulatory rows under
/// `--regulatory`. `RECORD_CLASS_KEYS` lists every key.
pub const REGULATORY_CLASS: &str = "regulatory_class";
pub const RECORD_CLASS_KEYS: [&str; 1] = [REGULATORY_CLASS];

/// The class of a row or a field `classify` finds no documented divergence for.
pub const UNEXPLAINED: &str = "unexplained_residual";

/// The key a `divergences` entry may carry, written by hand: vep-rs's consequence sets on
/// the entry's key in a run that has a reference sequence.
pub const WITH_REFERENCE_SETS: &str = "vep_rs_consequence_sets_with_reference_sequence";

/// The parity matrix `covers` and `flag-set` read when `--matrix` names no other.
pub const DEFAULT_MATRIX: &str = "tests/parity/matrix.tsv";

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn mask_classes_are_the_two_run_classes_the_harness_classes_and_the_corrected_ones() {
        let classes = mask_classes();
        assert_eq!(classes.len(), 2 + 4 + CORRECTED_BEHAVIOUR_CLASSES.len());
        assert_eq!(&classes[..2], ["regulatory_feature_rows", "order_only"]);
        assert_eq!(classes[2..6], HARNESS_ONLY_MASK_CLASSES);
        assert_eq!(classes[6..], CORRECTED_BEHAVIOUR_CLASSES);
        let distinct: HashSet<&str> = classes.iter().copied().collect();
        assert_eq!(distinct.len(), classes.len());
    }

    #[test]
    fn genotype_classes_are_the_sample_classes_then_the_record_classes() {
        let mut both = SAMPLE_CLASSES.to_vec();
        both.extend(RECORD_CLASSES);
        assert_eq!(both, GENOTYPE_CLASSES);
    }
}
