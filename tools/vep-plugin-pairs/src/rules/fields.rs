// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The fields each class reads, the classes that are handed more than a row's two values,
//! and the columns Ensembl VEP's own options write.
//!
//! `data/field_classes.tsv` holds one row per class, the class and the field names it
//! reads joined by `,`, each a full name or a prefix ending in `_`, in the order the rules
//! are consulted. Three header rows hold the classes that read a transcript's sequence
//! (`#transcript_reading`), the classes that read the engine's rows of the record
//! (`#sibling_allele`) and the classes that read the reference bases a row's allele stands
//! in for (`#record_reference`).

use std::collections::{HashMap, HashSet};
use std::sync::LazyLock;

const TABLE: &str = include_str!("../../data/field_classes.tsv");

/// The table of `data/field_classes.tsv`.
pub struct FieldClasses {
    /// The classes in the table's order, each with the field names it reads.
    pub classes: Vec<(&'static str, Vec<&'static str>)>,
    by_class: HashMap<&'static str, usize>,
    transcript_reading: Vec<&'static str>,
    sibling_allele: Vec<&'static str>,
    record_reference: Vec<&'static str>,
}

static FIELD_CLASSES: LazyLock<FieldClasses> = LazyLock::new(|| {
    let mut table = FieldClasses {
        classes: Vec::new(),
        by_class: HashMap::new(),
        transcript_reading: Vec::new(),
        sibling_allele: Vec::new(),
        record_reference: Vec::new(),
    };
    for line in TABLE.lines() {
        let (name, listed) = line.split_once('\t').expect("a tab in every row");
        let listed: Vec<&'static str> = listed.split(',').collect();
        match name {
            "#transcript_reading" => table.transcript_reading = listed,
            "#sibling_allele" => table.sibling_allele = listed,
            "#record_reference" => table.record_reference = listed,
            _ if name.starts_with('#') => {}
            _ => {
                table.by_class.insert(name, table.classes.len());
                table.classes.push((name, listed));
            }
        }
    }
    table
});

/// The table, read once.
pub fn field_classes() -> &'static FieldClasses {
    &FIELD_CLASSES
}

/// The field names `klass` reads; `None` for a name the table does not hold.
pub fn class_fields(klass: &str) -> Option<&'static [&'static str]> {
    let table = field_classes();
    table
        .by_class
        .get(klass)
        .map(|&at| table.classes[at].1.as_slice())
}

/// Whether `klass` reads a transcript's sequence, so that an empty value a row writes is
/// handed to it as written.
pub fn transcript_reading(klass: &str) -> bool {
    field_classes().transcript_reading.contains(&klass)
}

/// Whether `klass` reads the engine's rows of the same input record.
pub fn sibling_allele(klass: &str) -> bool {
    field_classes().sibling_allele.contains(&klass)
}

/// Whether `klass` reads the reference bases a row's allele stands in for.
pub fn record_reference(klass: &str) -> bool {
    field_classes().record_reference.contains(&klass)
}

/// The columns Ensembl VEP's own options write (the FLAG_FIELDS list of Constants.pm at
/// release/116.3 with the default output columns): a field of one of these names, in any
/// letter case, is no plugin's.
#[rustfmt::skip]
pub const CORE_COLUMN_NAMES: [&str; 112] = [
    "Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
    "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
    "Existing_variation", "Extra", "AF", "AFR_AF", "ALLELE_NUM", "AMBIGUITY", "AMR_AF",
    "APPRIS", "BAM_EDIT", "BIOTYPE", "CANONICAL", "CCDS", "CELL_TYPE", "CHECK_REF",
    "CLINVAR_SOMATIC_CLASSIFICATION", "CLIN_SIG", "DISTANCE", "DOMAINS", "EAS_AF", "ENSP",
    "EUR_AF", "EXON", "FLAGS", "FREQS", "GA4GH_VRS", "GENCODE_PRIMARY", "GENE_PHENO",
    "GIVEN_REF", "HGNC_ID", "HGVS_OFFSET", "HGVSc", "HGVSg", "HGVSp", "HIGH_INF_POS", "IMPACT",
    "IND", "INTRON", "MANE", "MANE_PLUS_CLINICAL", "MANE_SELECT", "MAX_AF", "MAX_AF_POPS",
    "MINIMISED", "MOTIF_NAME", "MOTIF_POS", "MOTIF_SCORE_CHANGE", "NEAREST", "OverlapBP",
    "OverlapPC", "PHENO", "PICK", "PUBMED", "PolyPhen", "REFSEQ_MATCH", "REFSEQ_OFFSET",
    "REF_ALLELE", "RefSeq", "SAS_AF", "SHIFT_LENGTH", "SIFT", "SOMATIC", "SOURCE", "SPDI",
    "STRAND", "SV", "SWISSPROT", "SYMBOL", "SYMBOL_SOURCE", "TRANSCRIPTION_FACTORS", "TREMBL",
    "TSL", "UNIPARC", "UNIPROT_ISOFORM", "UPLOADED_ALLELE", "USED_REF", "VARIANT_CLASS",
    "VAR_SYNONYMS", "ZYG", "miRNA", "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF", "gnomADe_EAS_AF", "gnomADe_FIN_AF", "gnomADe_MID_AF", "gnomADe_NFE_AF",
    "gnomADe_REMAINING_AF", "gnomADe_SAS_AF", "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF",
    "gnomADg_AMR_AF", "gnomADg_ASJ_AF", "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF",
    "gnomADg_NFE_AF", "gnomADg_REMAINING_AF", "gnomADg_SAS_AF",
];

/// `CORE_COLUMN_NAMES` in lower case.
pub fn core_columns() -> &'static HashSet<String> {
    static LOWERED: LazyLock<HashSet<String>> = LazyLock::new(|| {
        CORE_COLUMN_NAMES
            .iter()
            .map(|name| name.to_lowercase())
            .collect()
    });
    &LOWERED
}

/// Whether `name`, its letter case set aside, is one of `CORE_COLUMN_NAMES`.
pub fn is_core_column(name: &str) -> bool {
    core_columns().contains(&name.to_lowercase())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_holds_every_class_once_with_its_fields() {
        let table = field_classes();
        assert_eq!(table.classes.len(), 85);
        assert_eq!(table.by_class.len(), 85);
        assert_eq!(
            table.classes[0],
            ("gxf_reverse_strand_utr_sequence", vec!["HGVSp", "hgvsp"])
        );
        assert_eq!(
            table.classes[84].0,
            "plugin_value_follows_stop_codon_term_pair"
        );
        assert_eq!(
            table.classes[82],
            (
                "hgvsp_start_codon_test_past_the_shifted_allele",
                vec!["HGVSp", "hgvsp"]
            )
        );
        assert_eq!(table.classes[79].0, "frameshift_inside_a_held_stop_codon");
        assert_eq!(
            class_fields("frameshift_inside_a_held_stop_codon"),
            class_fields("stop_retained_without_stop_codon_test")
        );
        assert_eq!(
            table.classes[81],
            (
                "start_codon_kept_without_five_prime_utr",
                vec!["IMPACT", "impact", "HGVSp", "hgvsp"]
            )
        );
        assert_eq!(
            table.classes[78],
            (
                "start_codon_removed_without_five_prime_utr",
                vec!["IMPACT", "impact", "HGVSp", "hgvsp"]
            )
        );
        assert_eq!(
            class_fields("plugin_value_follows_stop_codon_term_pair"),
            Some(
                &[
                    "DownstreamProtein",
                    "ProteinLengthChange",
                    "NMD",
                    "LoF",
                    "LoF_info",
                    "LoF_filter",
                    "LoF_flags",
                    "REVEL",
                    "am_",
                    "alphamissense",
                    "ClinPred",
                    "dbNSFP_"
                ][..]
            )
        );
        assert_eq!(
            table.classes[80],
            (
                "hgvsp_shifted_indel_written_span_start_test",
                vec!["HGVSp", "hgvsp"]
            )
        );
        assert_eq!(
            table.classes[77].0,
            "stop_codon_settled_before_unknown_bases"
        );
        assert_eq!(
            table.classes[75],
            ("dbnsfp_mitochondrial_record_named_m", vec!["dbNSFP_"])
        );
        assert_eq!(
            class_fields("gwas_haplotype_cross_product"),
            Some(&["gwas_"][..])
        );
        assert_eq!(
            class_fields("loftee_gerp_zoom_level_summary"),
            Some(&["LoF_info", "LoF", "LoF_filter"][..])
        );
        assert_eq!(class_fields("order_only"), None);
        assert!(table.classes.iter().all(|(_, fields)| !fields.is_empty()));
    }

    #[test]
    fn the_membership_rows_name_classes_of_the_table() {
        let table = field_classes();
        assert_eq!(
            table.transcript_reading,
            [
                "downstream_utr_frame_past_coding_end",
                "downstream_start_past_own_codon",
                "plugin_hgvs_shifted_cds_coordinates"
            ]
        );
        assert_eq!(
            table.sibling_allele,
            [
                "splicevault_sibling_allele_rows",
                "alphamissense_sibling_allele_rows",
                "revel_sibling_allele_rows"
            ]
        );
        assert_eq!(
            table.record_reference,
            [
                "splice_term_past_the_record_span",
                "plugin_value_follows_splice_term_past_the_record_span",
                "stop_codon_overlap_insertion_range",
                "plugin_value_follows_stop_codon_term_pair",
                "stop_codon_completed_by_polyadenylation",
                "stop_codon_overlap_without_stop_codon",
                "stop_codon_behind_last_residue_without_stop_codon",
                "start_codon_overlap_without_start_codon"
            ]
        );
        for klass in table
            .transcript_reading
            .iter()
            .chain(&table.sibling_allele)
            .chain(&table.record_reference)
        {
            assert!(class_fields(klass).is_some(), "{klass}");
        }
        assert!(transcript_reading("downstream_start_past_own_codon"));
        assert!(sibling_allele("revel_sibling_allele_rows") && !sibling_allele("REVEL"));
        assert!(record_reference(
            "plugin_value_follows_splice_term_past_the_record_span"
        ));
        assert!(!record_reference("loeuf_zero_minimum"));
    }

    #[test]
    fn a_core_column_is_named_in_any_letter_case() {
        assert_eq!(core_columns().len(), 112);
        assert!(is_core_column("SYMBOL") && is_core_column("symbol") && is_core_column("HGVSc"));
        assert!(is_core_column("gnomadg_sas_af") && is_core_column("Extra"));
        assert!(!is_core_column("REVEL") && !is_core_column("SpliceAI_pred_SYMBOL"));
    }
}
