// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The JSON floors.

use crate::consts::{
    JSON_CONSEQUENCE_LISTS, JSON_RENAME, POSITION_RE, RECORD_LEVEL_COLUMNS, RECORD_LEVEL_COLUMN_RE,
};
use crate::entry::Entry;
use crate::lists::{field_objects, lowers_to};
use crate::pyfmt;
use std::collections::BTreeSet;

/// Where a JSON floor is counted: the fields it is counted under, or the note that says
/// why the format counts nothing for the column.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StrListOrStr {
    StrList(Vec<String>),
    Str(String),
}

pub const FLOOR_CUSTOM_SOURCE_NOT_APPLICABLE: &str = "a custom source's column, which the JSON format carries once per record under custom_annotations and not on its consequence entries, so the floor is checked in the text formats of the same case";

/// The record field holding a JSON record's custom annotations, an object keyed by source.
pub const CUSTOM_ANNOTATIONS_FIELD: &str = "custom_annotations";

/// The columns Ensembl VEP's own options write (the FLAG_FIELDS list of Constants.pm at
/// release/116.3 with the default output columns), which the JSON factory keys at the top of a
/// consequence entry wherever they are populated: a floor on one of them is never read inside a
/// plugin's object, whose member of the same name (SpliceAI's SYMBOL) is the plugin's copy, not
/// the column.
///
/// Membership only.
pub const CORE_COLUMNS: &[&str] = &[
    "af",
    "afr_af",
    "allele",
    "allele_num",
    "ambiguity",
    "amino_acids",
    "amr_af",
    "appris",
    "bam_edit",
    "biotype",
    "canonical",
    "ccds",
    "cdna_position",
    "cds_position",
    "cell_type",
    "check_ref",
    "clin_sig",
    "clinvar_somatic_classification",
    "codons",
    "consequence",
    "distance",
    "domains",
    "eas_af",
    "ensp",
    "eur_af",
    "existing_variation",
    "exon",
    "extra",
    "feature",
    "feature_type",
    "flags",
    "freqs",
    "ga4gh_vrs",
    "gencode_primary",
    "gene",
    "gene_pheno",
    "given_ref",
    "gnomade_af",
    "gnomade_afr_af",
    "gnomade_amr_af",
    "gnomade_asj_af",
    "gnomade_eas_af",
    "gnomade_fin_af",
    "gnomade_mid_af",
    "gnomade_nfe_af",
    "gnomade_remaining_af",
    "gnomade_sas_af",
    "gnomadg_af",
    "gnomadg_afr_af",
    "gnomadg_ami_af",
    "gnomadg_amr_af",
    "gnomadg_asj_af",
    "gnomadg_eas_af",
    "gnomadg_fin_af",
    "gnomadg_mid_af",
    "gnomadg_nfe_af",
    "gnomadg_remaining_af",
    "gnomadg_sas_af",
    "hgnc_id",
    "hgvs_offset",
    "hgvsc",
    "hgvsg",
    "hgvsp",
    "high_inf_pos",
    "impact",
    "ind",
    "intron",
    "location",
    "mane",
    "mane_plus_clinical",
    "mane_select",
    "max_af",
    "max_af_pops",
    "minimised",
    "mirna",
    "motif_name",
    "motif_pos",
    "motif_score_change",
    "nearest",
    "overlapbp",
    "overlappc",
    "pheno",
    "pick",
    "polyphen",
    "protein_position",
    "pubmed",
    "ref_allele",
    "refseq",
    "refseq_match",
    "refseq_offset",
    "sas_af",
    "shift_length",
    "sift",
    "somatic",
    "source",
    "spdi",
    "strand",
    "sv",
    "swissprot",
    "symbol",
    "symbol_source",
    "transcription_factors",
    "trembl",
    "tsl",
    "uniparc",
    "uniprot_isoform",
    "uploaded_allele",
    "uploaded_variation",
    "used_ref",
    "var_synonyms",
    "variant_class",
    "zyg",
];

/// The JSON fields a floor on the default-format `column` is counted under
/// (JSON.pm:311-385): the name lower-cased and renamed as `%RENAME_KEYS` renames it,
/// `Feature` under the per-kind id, SIFT and PolyPhen under either of their two keys, a
/// position under either of its start and end keys (JSON.pm:340-341 delete the start when
/// it is `?` and keep the end, which the default format prints as `?-N`); `None` for
/// a column with no row count in the JSON format: a record-level column, `Feature_type`,
/// which chooses the consequence list an entry joins and is deleted from the entry
/// (JSON.pm:348-350), and `Extra`, the text formats' container of the keyed fields
/// (VEP_output.pm:114-133).
pub fn json_floor_fields(column: &str) -> Option<Vec<String>> {
    let key = pyfmt::lower(column);
    if RECORD_LEVEL_COLUMNS.contains(&column)
        || RECORD_LEVEL_COLUMN_RE.is_match(column)
        || key == "feature_type"
        || key == "extra"
    {
        return None;
    }
    if key == "feature" {
        return Some(
            JSON_CONSEQUENCE_LISTS
                .iter()
                .filter(|(_, kind)| *kind != "intergenic")
                .map(|(_, kind)| format!("{kind}_id"))
                .collect(),
        );
    }
    // `re.match` holds the pattern to the start of the name.
    if let Some(coordinate) = POSITION_RE
        .captures(column)
        .filter(|found| found.get(0).is_some_and(|whole| whole.start() == 0))
        .and_then(|found| found.get(1))
    {
        let coordinate = pyfmt::lower(coordinate.as_str());
        return Some(vec![
            format!("{coordinate}_start"),
            format!("{coordinate}_end"),
        ]);
    }
    if key == "sift" || key == "polyphen" {
        return Some(vec![format!("{key}_prediction"), format!("{key}_score")]);
    }
    let renamed = JSON_RENAME
        .iter()
        .find(|(name, _)| *name == key)
        .map_or(key.as_str(), |(_, renamed)| renamed);
    Some(vec![renamed.to_string()])
}

/// Where a JSON floor on `column` is counted, read against the parsed entries of both
/// sides: the fields of `json_floor_fields` when any entry carries one of them (`None`
/// for a column the format has no row count for); else, for a column outside
/// `CORE_COLUMNS`, a plugin object's member, `<object>.<key>`, when a consequence entry of
/// either side carries an object (or an array of objects) whose member is the column's key;
/// else the custom-source note when a record's `custom_annotations` of either side names a
/// source the column spells (`<source>` or `<source>_<field>`); else the fields of
/// `json_floor_fields`, which no entry carries, so the floor reads as unmet on both sides.
pub fn json_floor_location(
    column: &str,
    expected: &[Entry],
    actual: &[Entry],
) -> Option<StrListOrStr> {
    let fields = json_floor_fields(column)?;
    let entries = || expected.iter().chain(actual);
    if entries().any(|e| fields.iter().any(|field| e.fields.contains_key(field))) {
        return Some(StrListOrStr::StrList(fields));
    }
    let key = pyfmt::lower(column);
    if !CORE_COLUMNS.contains(&key.as_str()) {
        // A field whose `{`-opening text parsed as no object on one entry is not parsed
        // again on the entries after it; an object field is parsed on every entry, since a
        // plugin may omit the member on some.
        let mut not_objects: BTreeSet<&str> = BTreeSet::new();
        for e in entries().filter(|e| !e.key.ends_with("|record")) {
            for (parent, value) in &e.fields {
                if not_objects.contains(parent.as_str()) || !value.starts_with('{') {
                    continue;
                }
                let objects = field_objects(value);
                if objects.is_empty() {
                    not_objects.insert(parent);
                } else if objects
                    .iter()
                    .filter_map(|object| object.as_dict())
                    .any(|members| members.keys().any(|name| lowers_to(name, &key)))
                {
                    return Some(StrListOrStr::StrList(vec![format!("{parent}.{key}")]));
                }
            }
        }
    }
    let mut sources: BTreeSet<String> = BTreeSet::new();
    for e in entries().filter(|e| e.key.ends_with("|record")) {
        let annotations = e
            .fields
            .get(CUSTOM_ANNOTATIONS_FIELD)
            .map_or("", String::as_str);
        for annotations_by_source in field_objects(annotations) {
            sources.extend(
                annotations_by_source
                    .as_dict()
                    .into_iter()
                    .flatten()
                    .map(|(source, _)| source.clone()),
            );
        }
    }
    let custom_source = sources.iter().any(|source| {
        column == source
            || column
                .strip_prefix(source.as_str())
                .is_some_and(|rest| rest.starts_with('_'))
    });
    if custom_source {
        return Some(StrListOrStr::Str(
            FLOOR_CUSTOM_SOURCE_NOT_APPLICABLE.to_string(),
        ));
    }
    Some(StrListOrStr::StrList(fields))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fields(column: &str) -> Option<Vec<String>> {
        json_floor_fields(column)
    }

    fn names(items: &[&str]) -> Option<Vec<String>> {
        Some(items.iter().map(|item| item.to_string()).collect())
    }

    #[test]
    fn a_floor_is_counted_under_the_json_keys_of_its_column() {
        assert_eq!(fields("IMPACT"), names(&["impact"]));
        assert_eq!(fields("SYMBOL"), names(&["gene_symbol"]));
        assert_eq!(fields("Consequence"), names(&["consequence_terms"]));
        assert_eq!(
            fields("Feature"),
            names(&["transcript_id", "regulatory_feature_id", "motif_feature_id"])
        );
        assert_eq!(fields("cDNA_position"), names(&["cdna_start", "cdna_end"]));
        assert_eq!(
            fields("Protein_POSITION"),
            names(&["protein_start", "protein_end"])
        );
        assert_eq!(fields("SIFT"), names(&["sift_prediction", "sift_score"]));
        assert_eq!(
            fields("PolyPhen"),
            names(&["polyphen_prediction", "polyphen_score"])
        );
        assert_eq!(fields("am_class"), names(&["am_class"]));
        for column in [
            "Existing_variation",
            "Location",
            "AF",
            "gnomADe_AF",
            "Feature_type",
            "FEATURE_TYPE",
            "Extra",
        ] {
            assert_eq!(fields(column), None, "{column}");
        }
        assert_eq!(fields("AFR"), names(&["afr"]));
        assert_eq!(fields("x_position_y"), names(&["x_position_y"]));
    }

    #[test]
    fn the_core_columns_are_lower_case_and_in_byte_order() {
        assert!(CORE_COLUMNS.is_sorted());
        assert!(CORE_COLUMNS
            .iter()
            .all(|column| pyfmt::lower(column) == *column));
    }

    fn plugin_entry() -> Entry {
        Entry {
            key: "rec0|G|ENST1|transcript".to_string(),
            record: Some(0),
            allele: "G".to_string(),
            feature: "ENST1".to_string(),
            consequence: "missense_variant".to_string(),
            fields: [
                (
                    "alphamissense",
                    r#"{"am_class":"likely_benign","am_pathogenicity":0.1261}"#,
                ),
                ("domains", "{x},{y}"),
                ("impact", "MODERATE"),
            ]
            .into_iter()
            .map(|(name, value)| (name.to_string(), value.to_string()))
            .collect(),
            ..Entry::default()
        }
    }

    #[test]
    fn a_floor_is_placed_on_a_direct_key_then_inside_a_plugin_object() {
        let entry = [plugin_entry()];
        let list = |items: &[&str]| names(items).map(StrListOrStr::StrList);
        assert_eq!(
            json_floor_location("IMPACT", &entry, &[]),
            list(&["impact"])
        );
        assert_eq!(
            json_floor_location("am_class", &[], &entry),
            list(&["alphamissense.am_class"])
        );
        assert_eq!(
            json_floor_location("am_genome", &entry, &entry),
            list(&["am_genome"])
        );
        assert_eq!(json_floor_location("Existing_variation", &entry, &[]), None);
        // A member spelt like one of Ensembl VEP's own columns is not the column.
        let mut spliceai = plugin_entry();
        spliceai.fields.insert(
            "spliceai".to_string(),
            r#"{"SYMBOL":"G1","DS_AG":"0.10"}"#.to_string(),
        );
        let spliceai = [spliceai];
        assert_eq!(
            json_floor_location("SYMBOL", &spliceai, &spliceai),
            list(&["gene_symbol"])
        );
        assert_eq!(
            json_floor_location("DS_AG", &spliceai, &spliceai),
            list(&["spliceai.ds_ag"])
        );
        let mut record = plugin_entry();
        record.key = "rec0|record".to_string();
        record.fields = [(
            CUSTOM_ANNOTATIONS_FIELD.to_string(),
            r#"{"sites":[{"name":"x"}]}"#.to_string(),
        )]
        .into_iter()
        .collect();
        let record = [record];
        for column in ["sites", "sites_AC"] {
            assert_eq!(
                json_floor_location(column, &record, &[]),
                Some(StrListOrStr::Str(
                    FLOOR_CUSTOM_SOURCE_NOT_APPLICABLE.to_string()
                )),
                "{column}"
            );
        }
        assert_eq!(
            json_floor_location("sitesX", &record, &[]),
            list(&["sitesx"])
        );
    }

    #[test]
    fn a_member_of_a_plugin_object_is_read_whatever_its_case() {
        use crate::lists::{json_object_member, populated_under};
        let entry = plugin_entry();
        let member = |parent: &str, key: &str| json_object_member(&entry, parent, key);
        assert_eq!(
            member("alphamissense", "AM_CLASS").as_deref(),
            Some("likely_benign")
        );
        assert_eq!(
            member("alphamissense", "am_pathogenicity").as_deref(),
            Some("0.1261")
        );
        assert_eq!(member("alphamissense", "am_genome"), None);
        assert_eq!(member("domains", "x"), None);
        assert_eq!(member("impact", "impact"), None);
        assert_eq!(member("absent", "x"), None);
        let entries = [entry.clone()];
        assert_eq!(populated_under(&entries, &["alphamissense.am_class"]), 1);
        assert_eq!(
            populated_under(
                &entries,
                &["am_class", "alphamissense.am_genome", "nothing.x"]
            ),
            0
        );
    }
}
