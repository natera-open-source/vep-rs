// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parse Ensembl variation VCF files and build variation JSON cache.
//!
//! Reads bgzipped VCF files from Ensembl (germline or somatic), extracts
//! variation records in the form Ensembl VEP's own VCF parser gives a record
//! (`Parser/VCF.pm`: the shared anchor base chopped, `-` for an emptied allele,
//! an insertion's `start` one past its `end`), joins the ClinVar allele and
//! somatic columns from the variation database dumps when they are given, and
//! writes the records as JSON files in the directory structure expected by
//! `vep-cli`'s `json_cache.rs::load_all_variations()`:
//!
//! ```text
//! {output_dir}/variations/{chr}/{region_start}-{region_end}.json
//! ```

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;
use tracing::{debug, info};

use crate::core_dumps::read_table;

/// Region size for binning variations (1 Mb), matching the transcript builder.
const REGION_SIZE: u64 = 1_000_000;

/// First 1-based position of the 1 Mb region containing `pos`.
///
/// Extracted so the binning tests exercise this function rather than re-inlining the
/// same arithmetic on literals: a test that recomputes the formula it is checking
/// cannot fail when the formula changes.
fn region_start(pos: u64) -> u64 {
    ((pos - 1) / REGION_SIZE) * REGION_SIZE + 1
}

/// VCF INFO field name to JSON frequency key mapping.
///
/// The Ensembl VCF stores population frequencies as INFO fields (e.g., `AFR_AF`),
/// while the JSON cache stores them as top-level keys (e.g., `AFR`).
const FREQ_FIELD_MAPPING: &[(&str, &str)] = &[
    ("AF", "AF"),
    ("AFR_AF", "AFR"),
    ("AMR_AF", "AMR"),
    ("EAS_AF", "EAS"),
    ("EUR_AF", "EUR"),
    ("SAS_AF", "SAS"),
    ("AA_AF", "AA"),
    ("EA_AF", "EA"),
    ("gnomADe_AF", "gnomADe"),
    ("gnomADe_AFR_AF", "gnomADe_AFR"),
    ("gnomADe_AMR_AF", "gnomADe_AMR"),
    ("gnomADe_ASJ_AF", "gnomADe_ASJ"),
    ("gnomADe_EAS_AF", "gnomADe_EAS"),
    ("gnomADe_FIN_AF", "gnomADe_FIN"),
    ("gnomADe_MID_AF", "gnomADe_MID"),
    ("gnomADe_NFE_AF", "gnomADe_NFE"),
    ("gnomADe_REMAINING_AF", "gnomADe_REMAINING"),
    ("gnomADe_SAS_AF", "gnomADe_SAS"),
    ("gnomADg_AF", "gnomADg"),
    ("gnomADg_AFR_AF", "gnomADg_AFR"),
    ("gnomADg_AMI_AF", "gnomADg_AMI"),
    ("gnomADg_AMR_AF", "gnomADg_AMR"),
    ("gnomADg_ASJ_AF", "gnomADg_ASJ"),
    ("gnomADg_EAS_AF", "gnomADg_EAS"),
    ("gnomADg_FIN_AF", "gnomADg_FIN"),
    ("gnomADg_MID_AF", "gnomADg_MID"),
    ("gnomADg_NFE_AF", "gnomADg_NFE"),
    ("gnomADg_REMAINING_AF", "gnomADg_REMAINING"),
    ("gnomADg_SAS_AF", "gnomADg_SAS"),
];

/// A variation record parsed from an Ensembl VCF, ready for JSON serialization.
///
/// Field names and types match `json_cache.rs::JsonVariation` so that
/// `load_all_variations()` can deserialize these without changes.
///
/// `start`, `end` and `allele_string` are in the form `Parser/VCF.pm` gives a VCF
/// record: `end` is `POS + len(REF) - 1` and, when the record is an indel whose REF
/// and ALTs share their first base, that base is chopped, `start` moves one past
/// `POS` and an emptied allele reads `-`, so an insertion has `start == end + 1`.
///
/// Frequency fields are serialized as top-level keys (e.g., `"AF"`, `"gnomADe"`)
/// matching the Perl cache JSON format that `JsonVariation` expects.
#[derive(Debug, Clone, Serialize)]
#[allow(non_snake_case)] // Field names match Perl JSON cache format (e.g., "gnomADe_AFR")
pub struct VariationRecord {
    pub variation_name: String,
    pub start: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end: Option<u64>,
    pub allele_string: String,
    pub strand: i8,
    pub failed: u8,
    pub somatic: u8,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub minor_allele: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub minor_allele_freq: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub clin_sig: Option<String>,
    pub phenotype_or_disease: u8,
    /// ClinVar significance per risk allele, `allele:term` pairs joined by `;`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub clin_sig_allele: Option<String>,
    /// The ClinVar reference allele behind `clin_sig_allele`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub clin_sig_ref_allele: Option<String>,
    /// ClinVar somatic classification, `phenotype(classification)` entries joined by `;`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub clinical_impact: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pubmed: Option<String>,

    #[serde(rename = "AF", skip_serializing_if = "Option::is_none")]
    pub af: Option<String>,
    #[serde(rename = "AFR", skip_serializing_if = "Option::is_none")]
    pub afr: Option<String>,
    #[serde(rename = "AMR", skip_serializing_if = "Option::is_none")]
    pub amr: Option<String>,
    #[serde(rename = "EAS", skip_serializing_if = "Option::is_none")]
    pub eas: Option<String>,
    #[serde(rename = "EUR", skip_serializing_if = "Option::is_none")]
    pub eur: Option<String>,
    #[serde(rename = "SAS", skip_serializing_if = "Option::is_none")]
    pub sas: Option<String>,
    #[serde(rename = "AA", skip_serializing_if = "Option::is_none")]
    pub aa: Option<String>,
    #[serde(rename = "EA", skip_serializing_if = "Option::is_none")]
    pub ea: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_AFR: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_AMR: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_ASJ: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_EAS: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_FIN: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_MID: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_NFE: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_REMAINING: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_SAS: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_AFR: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_AMI: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_AMR: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_ASJ: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_EAS: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_FIN: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_MID: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_NFE: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_REMAINING: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_SAS: Option<String>,
}

impl Default for VariationRecord {
    fn default() -> Self {
        Self {
            variation_name: String::new(),
            start: 0,
            end: None,
            allele_string: String::new(),
            strand: 1,
            failed: 0,
            somatic: 0,
            minor_allele: None,
            minor_allele_freq: None,
            clin_sig: None,
            phenotype_or_disease: 0,
            clin_sig_allele: None,
            clin_sig_ref_allele: None,
            clinical_impact: None,
            pubmed: None,
            af: None,
            afr: None,
            amr: None,
            eas: None,
            eur: None,
            sas: None,
            aa: None,
            ea: None,
            gnomADe: None,
            gnomADe_AFR: None,
            gnomADe_AMR: None,
            gnomADe_ASJ: None,
            gnomADe_EAS: None,
            gnomADe_FIN: None,
            gnomADe_MID: None,
            gnomADe_NFE: None,
            gnomADe_REMAINING: None,
            gnomADe_SAS: None,
            gnomADg: None,
            gnomADg_AFR: None,
            gnomADg_AMI: None,
            gnomADg_AMR: None,
            gnomADg_ASJ: None,
            gnomADg_EAS: None,
            gnomADg_FIN: None,
            gnomADg_MID: None,
            gnomADg_NFE: None,
            gnomADg_REMAINING: None,
            gnomADg_SAS: None,
        }
    }
}

impl VariationRecord {
    /// Set a frequency field by its JSON key name (e.g., "AF", "gnomADe_AFR").
    fn set_freq(&mut self, key: &str, value: String) {
        match key {
            "AF" => self.af = Some(value),
            "AFR" => self.afr = Some(value),
            "AMR" => self.amr = Some(value),
            "EAS" => self.eas = Some(value),
            "EUR" => self.eur = Some(value),
            "SAS" => self.sas = Some(value),
            "AA" => self.aa = Some(value),
            "EA" => self.ea = Some(value),
            "gnomADe" => self.gnomADe = Some(value),
            "gnomADe_AFR" => self.gnomADe_AFR = Some(value),
            "gnomADe_AMR" => self.gnomADe_AMR = Some(value),
            "gnomADe_ASJ" => self.gnomADe_ASJ = Some(value),
            "gnomADe_EAS" => self.gnomADe_EAS = Some(value),
            "gnomADe_FIN" => self.gnomADe_FIN = Some(value),
            "gnomADe_MID" => self.gnomADe_MID = Some(value),
            "gnomADe_NFE" => self.gnomADe_NFE = Some(value),
            "gnomADe_REMAINING" => self.gnomADe_REMAINING = Some(value),
            "gnomADe_SAS" => self.gnomADe_SAS = Some(value),
            "gnomADg" => self.gnomADg = Some(value),
            "gnomADg_AFR" => self.gnomADg_AFR = Some(value),
            "gnomADg_AMI" => self.gnomADg_AMI = Some(value),
            "gnomADg_AMR" => self.gnomADg_AMR = Some(value),
            "gnomADg_ASJ" => self.gnomADg_ASJ = Some(value),
            "gnomADg_EAS" => self.gnomADg_EAS = Some(value),
            "gnomADg_FIN" => self.gnomADg_FIN = Some(value),
            "gnomADg_MID" => self.gnomADg_MID = Some(value),
            "gnomADg_NFE" => self.gnomADg_NFE = Some(value),
            "gnomADg_REMAINING" => self.gnomADg_REMAINING = Some(value),
            "gnomADg_SAS" => self.gnomADg_SAS = Some(value),
            _ => {}
        }
    }
}

/// Statistics from a variation cache build.
pub struct VariationBuildStats {
    pub variation_count: usize,
    pub region_file_count: usize,
}

/// One germline ClinVar phenotype feature: the attributes
/// `PhenotypeFeatureAdaptor::get_clinsig_alleles_by_location` returns per feature.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct GermlineClinVar {
    variation_name: String,
    risk_allele: Option<String>,
    clinvar_clin_sig: Option<String>,
    clinvar_ref_allele: Option<String>,
}

/// One somatic ClinVar phenotype feature: the attributes
/// `PhenotypeFeatureAdaptor::get_somatic_clin_impact_by_location` returns per feature.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct SomaticClinVar {
    phenotype: String,
    somatic_clin_sig: Option<String>,
    oncogenic_clin_sig: Option<String>,
}

/// The ClinVar phenotype features of a variation database, keyed by the
/// `chromosome:start-end` span `Database/Variation.pm` looks them up by.
///
/// Ensembl VEP's dumper joins `clin_sig_allele`, `clin_sig_ref_allele` and
/// `clinical_impact` onto each variation feature from the `phenotype_feature`
/// and `phenotype_feature_attrib` tables (`Database/Variation.pm:161-162,
/// 196-221`); the Ensembl variation VCFs carry none of the three, so this type
/// reads the same tables from the database's text dumps.
#[derive(Debug, Default)]
pub struct ClinVarAttribs {
    germline: HashMap<String, Vec<GermlineClinVar>>,
    somatic: HashMap<String, Vec<SomaticClinVar>>,
}

impl ClinVarAttribs {
    /// Read the ClinVar phenotype features from the variation database dumps in `dir`:
    /// `source`, `seq_region`, `attrib_type`, `phenotype`, `phenotype_feature` and
    /// `phenotype_feature_attrib`, each as `.txt.gz` or `.txt`, on the sequence regions
    /// named by `chromosomes` (every region when the list is empty).
    ///
    /// Keeps the features of the `ClinVar` source whatever their object type, as the
    /// adaptor's two location queries do: a germline feature (`DNA_type` `Germline` or
    /// NULL) when it carries a `risk_allele` or a `clinvar_ref_allele` attribute
    /// (`PhenotypeFeatureAdaptor.pm:1253-1254`), a somatic feature (`DNA_type`
    /// `Somatic`) with its phenotype description (`:1303, 1320`). Errors when a table is
    /// missing; the rows are read with the core dump reader, MySQL's escapes decoded.
    pub fn load(dir: &Path, chromosomes: &[String]) -> Result<Self> {
        let mut clinvar_source_id = None;
        for row in read_table(dir, "source", 8)? {
            let row = row?;
            if row.field(1) == "ClinVar" {
                clinvar_source_id = Some(row.field(0).to_string());
            }
        }
        let Some(clinvar_source_id) = clinvar_source_id else {
            info!("No ClinVar source in {}", dir.display());
            return Ok(Self::default());
        };

        let mut seq_region_names: HashMap<String, String> = HashMap::new();
        for row in read_table(dir, "seq_region", 3)? {
            let row = row?;
            if chromosomes.is_empty() || chromosomes.iter().any(|c| c == row.field(1)) {
                seq_region_names.insert(row.field(0).to_string(), row.field(1).to_string());
            }
        }

        let mut attrib_codes: HashMap<String, String> = HashMap::new();
        for row in read_table(dir, "attrib_type", 4)? {
            let row = row?;
            attrib_codes.insert(row.field(0).to_string(), row.field(1).to_string());
        }

        let mut phenotypes: HashMap<String, String> = HashMap::new();
        for row in read_table(dir, "phenotype", 5)? {
            let row = row?;
            phenotypes.insert(
                row.field(0).to_string(),
                row.optional(3).unwrap_or("").to_string(),
            );
        }

        struct Feature {
            key: String,
            object_id: String,
            phenotype_id: String,
            somatic: bool,
            risk_allele: Option<String>,
            clinvar_clin_sig: Option<String>,
            clinvar_ref_allele: Option<String>,
            somatic_clin_sig: Option<String>,
            oncogenic_clin_sig: Option<String>,
        }
        let mut features: HashMap<String, Feature> = HashMap::new();
        for row in read_table(dir, "phenotype_feature", 12)? {
            let row = row?;
            if row.field(2) != clinvar_source_id {
                continue;
            }
            let Some(chr) = seq_region_names.get(row.field(7)) else {
                continue;
            };
            features.insert(
                row.field(0).to_string(),
                Feature {
                    key: format!("{chr}:{}-{}", row.field(8), row.field(9)),
                    object_id: row.field(5).to_string(),
                    phenotype_id: row.field(1).to_string(),
                    somatic: row.optional(11) == Some("Somatic"),
                    risk_allele: None,
                    clinvar_clin_sig: None,
                    clinvar_ref_allele: None,
                    somatic_clin_sig: None,
                    oncogenic_clin_sig: None,
                },
            );
        }

        // A code set twice on a feature keeps the last row, as the adaptor's
        // per-feature hash does (PhenotypeFeatureAdaptor.pm:1273-1276 assigns each
        // `code=value` in turn).
        for row in read_table(dir, "phenotype_feature_attrib", 3)? {
            let row = row?;
            let (Some(code), Some(f)) = (
                attrib_codes.get(row.field(1)),
                features.get_mut(row.field(0)),
            ) else {
                continue;
            };
            let value = Some(row.field(2).to_string());
            match code.as_str() {
                "risk_allele" => f.risk_allele = value,
                "clinvar_clin_sig" => f.clinvar_clin_sig = value,
                "clinvar_ref_allele" => f.clinvar_ref_allele = value,
                "somatic_clin_sig" => f.somatic_clin_sig = value,
                "oncogenic_clin_sig" => f.oncogenic_clin_sig = value,
                _ => {}
            }
        }

        let mut ordered: Vec<(u64, Feature)> = features
            .into_iter()
            .map(|(id, f)| (id.parse::<u64>().unwrap_or(u64::MAX), f))
            .collect();
        ordered.sort_by_key(|(id, _)| *id);

        let mut out = Self::default();
        for (_, f) in ordered {
            if f.somatic {
                out.somatic.entry(f.key).or_default().push(SomaticClinVar {
                    phenotype: phenotypes.get(&f.phenotype_id).cloned().unwrap_or_default(),
                    somatic_clin_sig: f.somatic_clin_sig,
                    oncogenic_clin_sig: f.oncogenic_clin_sig,
                });
            } else if f.risk_allele.is_some() || f.clinvar_ref_allele.is_some() {
                out.germline
                    .entry(f.key)
                    .or_default()
                    .push(GermlineClinVar {
                        variation_name: f.object_id,
                        risk_allele: f.risk_allele,
                        clinvar_clin_sig: f.clinvar_clin_sig,
                        clinvar_ref_allele: f.clinvar_ref_allele,
                    });
            }
        }
        info!(
            "Read {} germline and {} somatic ClinVar spans from {}",
            out.germline.len(),
            out.somatic.len(),
            dir.display()
        );
        Ok(out)
    }

    /// Whether no feature was read.
    pub fn is_empty(&self) -> bool {
        self.germline.is_empty() && self.somatic.is_empty()
    }

    /// Fill `clin_sig_allele`, `clin_sig_ref_allele` and `clinical_impact` on a record
    /// at `chr` as `Database/Variation.pm:196-221` does.
    ///
    /// A germline feature counts when its variation name is the record's: each
    /// `risk_allele:clinvar_clin_sig` pair (spaces in the term read `_`) once, joined
    /// by `;`, and its `clinvar_ref_allele` as the reference allele. Every somatic
    /// feature on the span counts, whatever its name, formatted by
    /// `_format_clinical_impact` (`:252-275`). The germline pairs keep their first-seen
    /// order where Perl's hash keys have none.
    fn apply(&self, var: &mut VariationRecord, chr: &str) {
        let key = format!("{chr}:{}-{}", var.start, var.end.unwrap_or(var.start));

        if let Some(features) = self.germline.get(&key) {
            let mut pairs: Vec<String> = Vec::new();
            let mut ref_allele = None;
            for f in features {
                let Some(clin_sig) = &f.clinvar_clin_sig else {
                    continue;
                };
                if f.variation_name != var.variation_name {
                    continue;
                }
                let pair = format!(
                    "{}:{}",
                    f.risk_allele.as_deref().unwrap_or(""),
                    clin_sig.replace(' ', "_")
                );
                if !pairs.contains(&pair) {
                    pairs.push(pair);
                }
                if let Some(r) = f.clinvar_ref_allele.as_deref().filter(|r| !r.is_empty()) {
                    ref_allele = Some(r.to_string());
                }
            }
            if !pairs.is_empty() {
                var.clin_sig_allele = Some(pairs.join(";"));
            }
            var.clin_sig_ref_allele = ref_allele;
        }

        if let Some(features) = self.somatic.get(&key) {
            let impact = format_clinical_impact(features);
            if !impact.is_empty() {
                var.clinical_impact = Some(impact);
            }
        }
    }
}

/// `Database/Variation.pm::_format_clinical_impact` (`:252-275`): for each somatic
/// feature, `phenotype(somatic_clin_sig)` with the classification's commas read `;`,
/// then `phenotype(oncogenicity:oncogenic_clin_sig)`, every entry joined by `;`.
fn format_clinical_impact(features: &[SomaticClinVar]) -> String {
    let mut entries = Vec::new();
    for f in features {
        if let Some(sig) = f.somatic_clin_sig.as_deref().filter(|s| !s.is_empty()) {
            entries.push(format!("{}({})", f.phenotype, sig.replace(',', ";")));
        }
        if let Some(onc) = f.oncogenic_clin_sig.as_deref().filter(|s| !s.is_empty()) {
            entries.push(format!("{}(oncogenicity:{onc})", f.phenotype));
        }
    }
    entries.join(";")
}

/// Parse an Ensembl variation VCF file into variation records.
///
/// Supports both plain-text and bgzipped (.vcf.gz) VCF files. Records where
/// FILTER is not PASS, ".", or empty are skipped. The `somatic` flag is applied
/// to all records from the file (Ensembl distributes germline and somatic VCFs
/// separately). With `clinvar`, each record gains the ClinVar allele and somatic
/// columns of the phenotype features on its span.
pub fn parse_variation_vcf(
    path: &Path,
    somatic: bool,
    clinvar: Option<&ClinVarAttribs>,
) -> Result<Vec<VariationRecord>> {
    use noodles::vcf;

    let file = std::fs::File::open(path)
        .with_context(|| format!("Failed to open VCF: {}", path.display()))?;

    let reader: Box<dyn BufRead + Send> = if path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz")
    {
        Box::new(std::io::BufReader::new(noodles::bgzf::io::Reader::new(
            file,
        )))
    } else {
        Box::new(std::io::BufReader::new(file))
    };

    let mut vcf_reader = vcf::io::Reader::new(reader);
    let header = vcf_reader
        .read_header()
        .with_context(|| format!("Failed to read VCF header: {}", path.display()))?;

    let mut records = Vec::new();
    let mut record = vcf::Record::default();
    let mut skipped = 0u64;

    loop {
        match vcf_reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => match parse_vcf_record(&record, &header, somatic) {
                Ok(Some(mut var)) => {
                    if let Some(clinvar) = clinvar {
                        clinvar.apply(&mut var, record.reference_sequence_name());
                    }
                    records.push(var);
                }
                Ok(None) => skipped += 1,
                Err(e) => {
                    debug!("Skipping malformed VCF record: {}", e);
                    skipped += 1;
                }
            },
            Err(e) => {
                debug!("Error reading VCF record, skipping: {}", e);
                skipped += 1;
            }
        }
    }

    info!(
        "Parsed {} variations from {} ({} skipped)",
        records.len(),
        path.display(),
        skipped
    );

    Ok(records)
}

/// Parse a single VCF record into a `VariationRecord`, or `None` if it should be skipped.
fn parse_vcf_record(
    record: &noodles::vcf::Record,
    header: &noodles::vcf::Header,
    somatic: bool,
) -> Result<Option<VariationRecord>> {
    let filters = record.filters();
    let filter_str: &str = filters.as_ref();
    if !filter_str.is_empty() && filter_str != "." && filter_str != "PASS" {
        return Ok(None);
    }

    let pos: u64 = record
        .variant_start()
        .context("Missing POS")?
        .context("Invalid POS")?
        .get() as u64;

    let ids = record.ids();
    let ids_str: &str = ids.as_ref();
    let variation_name = if ids_str.is_empty() || ids_str == "." {
        String::new()
    } else {
        ids_str.split(';').next().unwrap_or("").to_string()
    };

    let ref_allele = record.reference_bases().to_string();

    use noodles::vcf::variant::record::AlternateBases;
    let alts: Vec<String> = record
        .alternate_bases()
        .iter()
        .filter_map(|r| r.ok())
        .filter(|alt| *alt != ".")
        .map(|s| s.to_string())
        .collect();

    if alts.is_empty() {
        return Ok(None);
    }

    let (start, end, allele_string) = ensembl_alleles(pos, &ref_allele, &alts);

    let info = record.info();
    let info_raw: &str = info.as_ref();

    let minor_allele = get_info_string_value(record, header, "MA", info_raw);
    let minor_allele_freq =
        get_info_string_value(record, header, "MAF", info_raw).and_then(|s| s.parse::<f64>().ok());

    let clin_sig = get_info_string_value(record, header, "CLIN_SIG", info_raw);
    let phenotype_or_disease = u8::from(info_has_flag(
        record,
        header,
        "E_Phenotype_or_Disease",
        info_raw,
    ));

    let pubmed = get_info_string_value(record, header, "PUBMED", info_raw);

    let mut var = VariationRecord {
        variation_name,
        start,
        end: Some(end),
        allele_string,
        strand: 1, // Ensembl VCFs are always forward-strand
        failed: 0,
        somatic: if somatic { 1 } else { 0 },
        minor_allele,
        minor_allele_freq,
        clin_sig,
        phenotype_or_disease,
        pubmed,
        ..Default::default()
    };

    // A frequency entry names the allele as the allele string spells it
    // (`Dumper/Variation.pm:350, 378`), so the chopped ALTs key the entries.
    let cache_alts: Vec<String> = var
        .allele_string
        .split('/')
        .skip(1)
        .map(String::from)
        .collect();
    for &(vcf_field, json_key) in FREQ_FIELD_MAPPING {
        if let Some(raw_value) = get_info_string_value(record, header, vcf_field, info_raw) {
            let freq_str = format_freq_value(&raw_value, &cache_alts);
            if !freq_str.is_empty() {
                var.set_freq(json_key, freq_str);
            }
        }
    }

    Ok(Some(var))
}

/// The `start`, `end` and allele string Ensembl VEP's VCF parser gives a record
/// (`Parser/VCF.pm:254, 297-334`).
///
/// `end` is `POS + len(REF) - 1`. The record is an indel when any ALT's length
/// differs from REF's (`:297`); then, when REF and every ALT share their first base
/// (with one ALT, `:328`; with several, `:306-308`, an ALT carrying `*` left out of the
/// comparison), that base is chopped from each allele but a `*` one, `start` becomes
/// `POS + 1` and an emptied allele reads `-` (`:309-316, 331-334`), so an insertion
/// ends one base before it starts. Same-length alleles are never chopped, nor is a
/// record whose REF opens with a non-ASCII character: the chop is a one-byte slice, a
/// character boundary only in ASCII.
fn ensembl_alleles(pos: u64, ref_allele: &str, alts: &[String]) -> (u64, u64, String) {
    let end = pos + ref_allele.len() as u64 - 1;
    let is_indel = alts.iter().any(|alt| alt.len() != ref_allele.len());
    let first_base = ref_allele.as_bytes().first().filter(|b| b.is_ascii());
    let shares_first_base = |alt: &str| alt.as_bytes().first() == first_base;
    let carries_star = |allele: &str| allele.contains('*');
    let anchored = is_indel
        && first_base.is_some()
        && match alts {
            [single] => shares_first_base(single),
            _ => alts
                .iter()
                .filter(|alt| !carries_star(alt))
                .all(|alt| shares_first_base(alt)),
        };

    let chop = |allele: &str| -> String {
        if !anchored || carries_star(allele) {
            return allele.to_string();
        }
        match &allele[1..] {
            "" => "-".to_string(),
            rest => rest.to_string(),
        }
    };

    let start = if anchored { pos + 1 } else { pos };
    let mut allele_string = chop(ref_allele);
    for alt in alts {
        allele_string.push('/');
        allele_string.push_str(&chop(alt));
    }
    (start, end, allele_string)
}

/// Whether the INFO column carries the flag `key`, through the typed API first and
/// the raw `;`-separated text when the header does not declare it.
fn info_has_flag(
    record: &noodles::vcf::Record,
    header: &noodles::vcf::Header,
    key: &str,
    info_raw: &str,
) -> bool {
    if let Some(Ok(Some(_))) = record.info().get(header, key) {
        return true;
    }
    info_raw
        .split(';')
        .any(|field| field == key || field.split_once('=').is_some_and(|(k, _)| k == key))
}

/// Format a VCF frequency INFO value into the Perl cache "ALT:freq" format.
///
/// For single-allelic sites, this produces "ALT:0.123".
/// For multi-allelic sites with comma-separated frequencies, it produces
/// "ALT1:freq1,ALT2:freq2".
fn format_freq_value(raw_value: &str, alts: &[String]) -> String {
    let freqs: Vec<&str> = raw_value.split(',').collect();

    let mut parts = Vec::new();
    for (i, freq_str) in freqs.iter().enumerate() {
        let freq_str = freq_str.trim();
        if freq_str.is_empty() || freq_str == "." {
            continue;
        }
        let allele = if i < alts.len() {
            &alts[i]
        } else if alts.len() == 1 {
            &alts[0]
        } else {
            continue;
        };
        parts.push(format!("{}:{}", allele, freq_str));
    }

    parts.join(",")
}

/// Try to extract a string INFO field value. Uses the noodles typed API first,
/// then falls back to raw string parsing.
fn get_info_string_value(
    record: &noodles::vcf::Record,
    header: &noodles::vcf::Header,
    key: &str,
    info_raw: &str,
) -> Option<String> {
    use noodles::vcf::variant::record::info::field::Value;

    if let Some(Ok(Some(value))) = record.info().get(header, key) {
        return match value {
            Value::String(s) => {
                let s = s.to_string();
                if s.is_empty() || s == "." {
                    None
                } else {
                    Some(s)
                }
            }
            Value::Integer(n) => Some(n.to_string()),
            Value::Float(f) => Some(f.to_string()),
            Value::Array(arr) => {
                use noodles::vcf::variant::record::info::field::value::Array;
                match arr {
                    Array::Integer(iter) => {
                        let vals: Vec<String> = iter
                            .iter()
                            .filter_map(|r| r.ok())
                            .filter_map(|v| v.map(|n| n.to_string()))
                            .collect();
                        if vals.is_empty() {
                            None
                        } else {
                            Some(vals.join(","))
                        }
                    }
                    Array::Float(iter) => {
                        let vals: Vec<String> = iter
                            .iter()
                            .filter_map(|r| r.ok())
                            .filter_map(|v| v.map(|f| f.to_string()))
                            .collect();
                        if vals.is_empty() {
                            None
                        } else {
                            Some(vals.join(","))
                        }
                    }
                    Array::String(iter) => {
                        let vals: Vec<String> = iter
                            .iter()
                            .filter_map(|r| r.ok())
                            .filter_map(|v| v.map(|s| s.to_string()))
                            .collect();
                        if vals.is_empty() {
                            None
                        } else {
                            Some(vals.join(","))
                        }
                    }
                    Array::Character(iter) => {
                        let vals: Vec<String> = iter
                            .iter()
                            .filter_map(|r| r.ok())
                            .filter_map(|v| v.map(|c| c.to_string()))
                            .collect();
                        if vals.is_empty() {
                            None
                        } else {
                            Some(vals.join(","))
                        }
                    }
                }
            }
            Value::Flag => Some("1".to_string()),
            Value::Character(c) => Some(c.to_string()),
        };
    }

    get_raw_info_value(info_raw, key).map(|s| s.to_string())
}

/// Parse a single INFO field value from the raw semicolon-delimited INFO string.
fn get_raw_info_value<'a>(info_raw: &'a str, key: &str) -> Option<&'a str> {
    if info_raw == "." || info_raw.is_empty() {
        return None;
    }
    for field in info_raw.split(';') {
        if let Some((k, v)) = field.split_once('=') {
            if k == key && !v.is_empty() && v != "." {
                return Some(v);
            }
        }
    }
    None
}

/// Build variation JSON cache files from parsed variation records.
///
/// Groups variations into 1Mb regions and writes one JSON file per region,
/// matching the directory structure expected by `json_cache.rs::load_all_variations()`:
/// `{output_dir}/variations/{chr}/{region_start}-{region_end}.json`
pub fn build_variation_cache(
    records: &[VariationRecord],
    output_dir: &Path,
    chr: &str,
) -> Result<VariationBuildStats> {
    let mut by_region: HashMap<u64, Vec<&VariationRecord>> = HashMap::new();
    for rec in records {
        by_region
            .entry(region_start(rec.start))
            .or_default()
            .push(rec);
    }

    let variations_dir = output_dir.join("variations").join(chr);
    std::fs::create_dir_all(&variations_dir)
        .with_context(|| format!("Failed to create {}", variations_dir.display()))?;

    let mut region_file_count = 0usize;

    for (region_start, region_records) in &by_region {
        let region_end = region_start + REGION_SIZE - 1;
        let filename = format!("{}-{}.json", region_start, region_end);
        let filepath = variations_dir.join(&filename);

        let json =
            serde_json::to_string(region_records).context("Failed to serialize variations")?;

        std::fs::write(&filepath, json)
            .with_context(|| format!("Failed to write {}", filepath.display()))?;

        region_file_count += 1;
    }

    info!(
        "chr {}: {} variations in {} regions",
        chr,
        records.len(),
        region_file_count
    );

    Ok(VariationBuildStats {
        variation_count: records.len(),
        region_file_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal VCF as bytes (plain text, not bgzipped) for testing.
    fn make_vcf_bytes(info_headers: &str, records: &str) -> Vec<u8> {
        format!(
            "##fileformat=VCFv4.1\n\
             ##INFO=<ID=MA,Number=1,Type=String,Description=\"Minor allele\">\n\
             ##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n\
             ##INFO=<ID=CLIN_SIG,Number=.,Type=String,Description=\"Clinical significance\">\n\
             ##INFO=<ID=PUBMED,Number=.,Type=String,Description=\"PubMed IDs\">\n\
             ##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
             ##INFO=<ID=AFR_AF,Number=A,Type=Float,Description=\"AFR allele frequency\">\n\
             ##INFO=<ID=EUR_AF,Number=A,Type=Float,Description=\"EUR allele frequency\">\n\
             ##INFO=<ID=gnomADe_AF,Number=A,Type=Float,Description=\"gnomAD exome AF\">\n\
             {info_headers}\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             {records}"
        )
        .into_bytes()
    }

    /// Write a VCF to a temp file and parse it.
    fn parse_vcf_str(info_headers: &str, records: &str, somatic: bool) -> Vec<VariationRecord> {
        parse_vcf_str_with(info_headers, records, somatic, None)
    }

    fn parse_vcf_str_with(
        info_headers: &str,
        records: &str,
        somatic: bool,
        clinvar: Option<&ClinVarAttribs>,
    ) -> Vec<VariationRecord> {
        let dir = tempfile::tempdir().unwrap();
        let vcf_path = dir.path().join("test.vcf");
        let data = make_vcf_bytes(info_headers, records);
        std::fs::write(&vcf_path, &data).unwrap();
        parse_variation_vcf(&vcf_path, somatic, clinvar).unwrap()
    }

    /// INFO lines of an Ensembl variation VCF header: the `E_*` evidence flags, a
    /// `CLIN_*` flag and `TSA`.
    const ENSEMBL_INFO_HEADERS: &str = "##INFO=<ID=E_Phenotype_or_Disease,Number=0,Type=Flag,Description=\"Phenotype_or_Disease\">\n\
         ##INFO=<ID=E_Freq,Number=0,Type=Flag,Description=\"Frequency\">\n\
         ##INFO=<ID=CLIN_pathogenic,Number=0,Type=Flag,Description=\"pathogenic\">\n\
         ##INFO=<ID=TSA,Number=1,Type=String,Description=\"Type of sequence alteration\">\n";

    /// Ensembl VEP's VCF parser chops the shared anchor base of an indel, moves
    /// `start` one past `POS`, leaves `end` at `POS + len(REF) - 1` and spells an
    /// emptied allele `-` (`Parser/VCF.pm:254, 309-316, 331-334`): a `TA/T` deletion at
    /// 5030218 is `A/-` at 5030219-5030219, the mapping Ensembl itself reports for
    /// rs1416749598; a `C/CGT` insertion at 100 is `-/GT` with `start` 101 and `end`
    /// 100; the shared base of `AAA/A,AA,AAAA` goes too (`:306-308`); a `G/GT,A`
    /// record whose ALTs disagree on the first base keeps every allele (`:306-308`);
    /// an `A/AT,*` record chops around the `*` allele, which stays as written
    /// (`:306, 315`).
    #[test]
    fn indel_allele_string_is_ensembl_style() {
        let records = parse_vcf_str(
            "",
            "21\t5030218\trs1416749598\tTA\tT\t.\t.\t.\n\
             21\t100\tins1\tC\tCGT\t.\t.\t.\n\
             21\t5030913\trs1232130144\tAAA\tA,AA,AAAA\t.\t.\t.\n\
             21\t300\tmixed1\tG\tGT,A\t.\t.\t.\n\
             21\t400\tmnv1\tAC\tGT\t.\t.\t.\n\
             21\t500\tspan1\tA\tAT,*\t.\t.\t.\n",
            false,
        );
        assert_eq!(records.len(), 6);

        let del = &records[0];
        assert_eq!(del.allele_string, "A/-");
        assert_eq!((del.start, del.end), (5030219, Some(5030219)));

        let ins = &records[1];
        assert_eq!(ins.allele_string, "-/GT");
        assert_eq!((ins.start, ins.end), (101, Some(100)));

        let multi = &records[2];
        assert_eq!(multi.allele_string, "AA/-/A/AAA");
        assert_eq!((multi.start, multi.end), (5030914, Some(5030915)));

        let mixed = &records[3];
        assert_eq!(mixed.allele_string, "G/GT/A");
        assert_eq!((mixed.start, mixed.end), (300, Some(300)));

        let mnv = &records[4];
        assert_eq!(mnv.allele_string, "AC/GT");
        assert_eq!((mnv.start, mnv.end), (400, Some(401)));

        let spanning = &records[5];
        assert_eq!(spanning.allele_string, "-/T/*");
        assert_eq!((spanning.start, spanning.end), (501, Some(500)));
    }

    /// An allele whose first character is not ASCII never anchors the chop: the
    /// `ÄT/Ä` record keeps both alleles whole with `start` at `POS` rather than
    /// panicking on a byte slice through the character.
    #[test]
    fn non_ascii_allele_is_left_unchopped() {
        let records = parse_vcf_str("", "21\t100\trs1\tÄT\tÄ\t.\t.\t.\n", false);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].allele_string, "ÄT/Ä");
        assert_eq!(records[0].start, 100);
    }

    /// A frequency entry is keyed by the chopped allele (`Dumper/Variation.pm:350,
    /// 378`): the `C/CGT` insertion's `AF=0.25` reads `GT:0.25`, not `CGT:0.25`.
    #[test]
    fn frequency_entries_name_the_chopped_allele() {
        let records = parse_vcf_str("", "21\t100\tins1\tC\tCGT\t.\t.\tAF=0.25\n", false);
        assert_eq!(records[0].af, Some("GT:0.25".to_string()));
    }

    /// `phenotype_or_disease` is 1 when the variation's evidence set names
    /// Phenotype_or_Disease (`Database/Variation.pm:173`), which the Ensembl VCFs
    /// write as the `E_Phenotype_or_Disease` flag; a record flagged `CLIN_pathogenic`
    /// without that evidence stays 0, and a record with the evidence but no clinical
    /// significance reads 1.
    #[test]
    fn phenotype_or_disease_from_evidence_flag() {
        let records = parse_vcf_str(
            ENSEMBL_INFO_HEADERS,
            "21\t100\trs1\tA\tG\t.\t.\tTSA=SNV;E_Freq;E_Phenotype_or_Disease\n\
             21\t200\trs2\tA\tG\t.\t.\tTSA=SNV;E_Freq;CLIN_pathogenic\n\
             21\t300\trs3\tA\tG\t.\t.\tTSA=SNV;E_Freq\n",
            false,
        );
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].phenotype_or_disease, 1);
        assert_eq!(records[1].phenotype_or_disease, 0);
        assert_eq!(records[2].phenotype_or_disease, 0);
    }

    /// The evidence flag is read from a VCF whose header does not declare it.
    #[test]
    fn phenotype_or_disease_flag_read_without_header_line() {
        let records = parse_vcf_str(
            "",
            "21\t100\trs1\tA\tG\t.\t.\tE_Phenotype_or_Disease\n",
            false,
        );
        assert_eq!(records[0].phenotype_or_disease, 1);
    }

    /// Write each table's tab-separated rows to `<table>.txt.gz` under `dir`,
    /// gzipped as the Ensembl FTP publishes the dumps, with `\N` for NULL.
    fn write_dump_dir(dir: &Path, tables: &[(&str, &str)]) {
        use std::io::Write;
        for (table, rows) in tables {
            let path = dir.join(format!("{table}.txt.gz"));
            let file = std::fs::File::create(&path).unwrap();
            let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
            gz.write_all(rows.as_bytes()).unwrap();
            gz.finish().unwrap();
        }
    }

    /// Variation dump rows for one ClinVar span: 21:100-100 carries two germline
    /// features of rs1 (risk alleles T and C, significances `pathogenic` and
    /// `uncertain significance`, reference allele A), one germline feature of another
    /// variation on the same span, one germline feature whose risk allele is set twice
    /// (A, then C), and four somatic features (a somatic classification with a comma,
    /// an oncogenicity, a classification on a phenotype whose description carries
    /// MySQL-escaped backslashes, and an oncogenicity on a `Gene`-type feature of the
    /// span, which the adaptor's location query joins like any other).
    const DUMP_TABLES: &[(&str, &str)] = &[
        ("source", "1\tdbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation\n32\tClinVar\t202509\t\\N\t\\N\t\\N\tmixed\tphenotype_feature\n"),
        ("seq_region", "131537\t21\t4\n131538\t22\t4\n"),
        ("attrib_type", "480\tclinvar_clin_sig\tClinVar clinical significance\t\\N\n483\trisk_allele\tRisk allele\t\\N\n577\tsomatic_clin_sig\tClinVar somatic classification\t\\N\n578\toncogenic_clin_sig\tClinVar somatic classification of oncogenicity\t\\N\n583\tclinvar_ref_allele\tClinVar reference allele\t\\N\n"),
        ("phenotype", "7\t\\N\t\\N\tBreast cancer\t\\N\n8\t\\N\t\\N\tLung carcinoma\t\\N\n9\t\\N\t\\N\tNot provided\t\\N\n10\t\\N\t\\N\tCarcinoma\\\\of\\\\lung\t\\N\n"),
        ("phenotype_feature", "\
1001\t7\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tGermline\n\
1002\t9\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\t\\N\n\
1003\t7\t32\t\\N\tVariation\trs9\t1\t131537\t100\t100\t1\tGermline\n\
1004\t7\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tGermline\n\
1005\t8\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tSomatic\n\
1006\t7\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tSomatic\n\
1007\t7\t1\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tGermline\n\
1008\t7\t32\t\\N\tGene\tENSG1\t1\t131537\t100\t100\t1\tGermline\n\
1009\t10\t32\t\\N\tVariation\trs1\t1\t131537\t100\t100\t1\tSomatic\n\
1010\t9\t32\t\\N\tGene\tENSG1\t1\t131537\t100\t100\t1\tSomatic\n"),
        ("phenotype_feature_attrib", "\
1001\t480\tpathogenic\n1001\t483\tT\n1001\t583\tA\n\
1002\t480\tuncertain significance\n1002\t483\tC\n\
1003\t480\tbenign\n1003\t483\tG\n\
1004\t480\tnot provided\n\
1005\t577\tTier I - Strong,Tier II - Potential\n\
1006\t578\tOncogenic\n\
1007\t480\tbenign\n1007\t483\tG\n\
1009\t577\tTier II - Potential\n1004\t483\tA\n1004\t483\tC\n\
1010\t578\tLikely oncogenic\n"),
    ];

    /// `clin_sig_allele` joins each `risk_allele:clinvar_clin_sig` pair of the germline
    /// ClinVar features that carry the record's own name on its span, spaces in the
    /// term as `_`, and `clin_sig_ref_allele` is their `clinvar_ref_allele`
    /// (`Database/Variation.pm:196-214`; the pair grammar is what Ensembl VEP reports
    /// for rs113488022, `T:pathogenic;...;C:uncertain_significance`). A feature of
    /// another variation on the same span and a dbSNP-source feature are left out; a
    /// feature whose risk allele is set twice keeps the last row
    /// (`PhenotypeFeatureAdaptor.pm:1273-1276`), so rs1's `not provided` pairs with C.
    /// `clinical_impact` is `_format_clinical_impact` (`:252-275`) over the somatic
    /// features: `phenotype(classification)` with `,` read `;`, then
    /// `phenotype(oncogenicity:term)`, joined by `;`, the dump's `\\` read as one
    /// backslash. A record on another span gets none.
    #[test]
    fn clinical_impact_from_clinvar_vcf() {
        let dir = tempfile::tempdir().unwrap();
        write_dump_dir(dir.path(), DUMP_TABLES);
        let clinvar = ClinVarAttribs::load(dir.path(), &[]).unwrap();
        assert!(!clinvar.is_empty());

        let records = parse_vcf_str_with(
            ENSEMBL_INFO_HEADERS,
            "21\t100\trs1\tA\tC,T\t.\t.\tTSA=SNV;E_Phenotype_or_Disease;CLIN_pathogenic\n\
             21\t200\trs2\tA\tG\t.\t.\tTSA=SNV\n",
            false,
            Some(&clinvar),
        );
        assert_eq!(records.len(), 2);

        let r = &records[0];
        assert_eq!(
            r.clin_sig_allele,
            Some("T:pathogenic;C:uncertain_significance;C:not_provided".to_string())
        );
        assert_eq!(r.clin_sig_ref_allele, Some("A".to_string()));
        assert_eq!(
            r.clinical_impact,
            Some("Lung carcinoma(Tier I - Strong;Tier II - Potential);Breast cancer(oncogenicity:Oncogenic);Carcinoma\\of\\lung(Tier II - Potential);Not provided(oncogenicity:Likely oncogenic)".to_string())
        );

        let other = &records[1];
        assert_eq!(other.clin_sig_allele, None);
        assert_eq!(other.clin_sig_ref_allele, None);
        assert_eq!(other.clinical_impact, None);

        let json: serde_json::Value = serde_json::to_value(r).unwrap();
        assert_eq!(
            json["clin_sig_allele"],
            "T:pathogenic;C:uncertain_significance;C:not_provided"
        );
        assert_eq!(json["clin_sig_ref_allele"], "A");
        assert_eq!(
            json["clinical_impact"],
            "Lung carcinoma(Tier I - Strong;Tier II - Potential);Breast cancer(oncogenicity:Oncogenic);Carcinoma\\of\\lung(Tier II - Potential);Not provided(oncogenicity:Likely oncogenic)"
        );
        let json_other: serde_json::Value = serde_json::to_value(other).unwrap();
        assert!(json_other.get("clin_sig_allele").is_none());
        assert!(json_other.get("clinical_impact").is_none());
    }

    /// The dump reader accepts the uncompressed `.txt` spelling, keeps nothing from a
    /// chromosome outside the list it is given, and reports an empty table set when
    /// the source table names no ClinVar.
    #[test]
    fn clinvar_dump_reader_accepts_plain_text_and_absent_clinvar() {
        let dir = tempfile::tempdir().unwrap();
        for (table, rows) in DUMP_TABLES {
            std::fs::write(dir.path().join(format!("{table}.txt")), rows).unwrap();
        }
        let clinvar = ClinVarAttribs::load(dir.path(), &[]).unwrap();
        assert!(!clinvar.is_empty());
        let other_chromosome = ClinVarAttribs::load(dir.path(), &["22".to_string()]).unwrap();
        assert!(other_chromosome.is_empty());

        let no_clinvar = tempfile::tempdir().unwrap();
        std::fs::write(
            no_clinvar.path().join("source.txt"),
            "1\tdbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation\n",
        )
        .unwrap();
        assert!(ClinVarAttribs::load(no_clinvar.path(), &[])
            .unwrap()
            .is_empty());
    }

    /// A dump directory without the tables is an error naming the missing table.
    #[test]
    fn clinvar_dump_reader_names_the_missing_table() {
        let dir = tempfile::tempdir().unwrap();
        let err = ClinVarAttribs::load(dir.path(), &[]).unwrap_err();
        assert!(err.to_string().contains("source.txt"), "{err}");
    }

    #[test]
    fn test_parse_snv() {
        let records = parse_vcf_str(
            "",
            "21\t25000100\trs699\tA\tG\t.\tPASS\tMA=G;MAF=0.2345\n",
            false,
        );
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.variation_name, "rs699");
        assert_eq!(r.start, 25000100);
        assert_eq!(r.end, Some(25000100));
        assert_eq!(r.allele_string, "A/G");
        assert_eq!(r.strand, 1);
        assert_eq!(r.failed, 0);
        assert_eq!(r.somatic, 0);
        assert_eq!(r.minor_allele, Some("G".to_string()));
        assert_eq!(r.minor_allele_freq, Some(0.2345));
    }

    #[test]
    fn test_parse_deletion() {
        let records = parse_vcf_str("", "21\t100\trs1\tATG\tA\t.\t.\t.\n", false);
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.start, 101);
        assert_eq!(r.end, Some(102));
        assert_eq!(r.allele_string, "TG/-");
    }

    #[test]
    fn test_parse_multi_allelic() {
        let records = parse_vcf_str("", "21\t200\trs2\tA\tG,T\t.\tPASS\tAF=0.3,0.1\n", false);
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.allele_string, "A/G/T");
        assert_eq!(r.af, Some("G:0.3,T:0.1".to_string()));
    }

    #[test]
    fn test_somatic_flag() {
        let records = parse_vcf_str("", "21\t100\tCOSV1\tA\tG\t.\t.\t.\n", true);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].somatic, 1);
    }

    #[test]
    fn test_filter_skips_non_pass() {
        let records = parse_vcf_str(
            "",
            "21\t100\trs1\tA\tG\t.\tPASS\t.\n\
             21\t200\trs2\tA\tG\t.\tLowQual\t.\n\
             21\t300\trs3\tA\tG\t.\t.\t.\n",
            false,
        );
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].variation_name, "rs1");
        assert_eq!(records[1].variation_name, "rs3");
    }

    #[test]
    fn test_clin_sig_and_pubmed() {
        let records = parse_vcf_str(
            "",
            "21\t100\trs1\tA\tG\t.\tPASS\tCLIN_SIG=benign,likely_benign;PUBMED=12345,67890\n",
            false,
        );
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.clin_sig, Some("benign,likely_benign".to_string()));
        assert_eq!(r.pubmed, Some("12345,67890".to_string()));
        assert_eq!(r.phenotype_or_disease, 0);
    }

    #[test]
    fn test_no_clin_sig_phenotype_zero() {
        let records = parse_vcf_str("", "21\t100\trs1\tA\tG\t.\tPASS\t.\n", false);
        assert_eq!(records[0].phenotype_or_disease, 0);
    }

    #[test]
    fn test_frequency_formatting_single_allele() {
        let result = format_freq_value("0.123", &["G".to_string()]);
        assert_eq!(result, "G:0.123");
    }

    #[test]
    fn test_frequency_formatting_multi_allele() {
        let result = format_freq_value("0.3,0.1", &["G".to_string(), "T".to_string()]);
        assert_eq!(result, "G:0.3,T:0.1");
    }

    #[test]
    fn test_frequency_formatting_missing_dot() {
        let result = format_freq_value(".", &["G".to_string()]);
        assert_eq!(result, "");
    }

    #[test]
    fn test_population_frequencies() {
        let records = parse_vcf_str(
            "",
            "21\t100\trs1\tA\tG\t.\tPASS\tAF=0.5;AFR_AF=0.3;EUR_AF=0.6;gnomADe_AF=0.45\n",
            false,
        );
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.af, Some("G:0.5".to_string()));
        assert_eq!(r.afr, Some("G:0.3".to_string()));
        assert_eq!(r.eur, Some("G:0.6".to_string()));
        assert_eq!(r.gnomADe, Some("G:0.45".to_string()));
    }

    #[test]
    fn test_build_variation_cache() {
        let dir = tempfile::tempdir().unwrap();
        let records = vec![
            VariationRecord {
                variation_name: "rs1".to_string(),
                start: 500_000,
                end: Some(500_000),
                allele_string: "A/G".to_string(),
                strand: 1,
                ..Default::default()
            },
            VariationRecord {
                variation_name: "rs2".to_string(),
                start: 1_500_000,
                end: Some(1_500_000),
                allele_string: "C/T".to_string(),
                strand: 1,
                ..Default::default()
            },
            VariationRecord {
                variation_name: "rs3".to_string(),
                start: 999_999,
                end: Some(999_999),
                allele_string: "G/A".to_string(),
                strand: 1,
                ..Default::default()
            },
        ];

        let stats = build_variation_cache(&records, dir.path(), "21").unwrap();
        assert_eq!(stats.variation_count, 3);
        assert_eq!(stats.region_file_count, 2); // 1-1000000 and 1000001-2000000

        let region1 = dir.path().join("variations/21/1-1000000.json");
        let region2 = dir.path().join("variations/21/1000001-2000000.json");
        assert!(region1.exists());
        assert!(region2.exists());

        let data1: Vec<serde_json::Value> =
            serde_json::from_str(&std::fs::read_to_string(&region1).unwrap()).unwrap();
        assert_eq!(data1.len(), 2);

        let data2: Vec<serde_json::Value> =
            serde_json::from_str(&std::fs::read_to_string(&region2).unwrap()).unwrap();
        assert_eq!(data2.len(), 1);
        assert_eq!(data2[0]["variation_name"], "rs2");
    }

    #[test]
    fn test_region_binning_boundary() {
        // Position 1 -> region 1-1000000
        assert_eq!(region_start(1), 1);
        // Position 1000000 -> region 1-1000000
        assert_eq!(region_start(1_000_000), 1);
        // Position 1000001 -> region 1000001-2000000
        assert_eq!(region_start(1_000_001), 1_000_001);
    }

    #[test]
    fn test_json_roundtrip_compatible_with_loader() {
        // The serialized JSON must deserialize through the serde structure that
        // json_cache.rs uses.
        let rec = VariationRecord {
            variation_name: "rs699".to_string(),
            start: 230710048,
            end: Some(230710048),
            allele_string: "A/G".to_string(),
            strand: 1,
            failed: 0,
            somatic: 0,
            minor_allele: Some("G".to_string()),
            minor_allele_freq: Some(0.2345),
            clin_sig: Some("benign".to_string()),
            phenotype_or_disease: 1,
            pubmed: Some("12345".to_string()),
            af: Some("G:0.5".to_string()),
            afr: Some("G:0.3".to_string()),
            gnomADe: Some("G:0.45".to_string()),
            ..Default::default()
        };

        let json = serde_json::to_string(&rec).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();

        assert_eq!(parsed["variation_name"], "rs699");
        assert_eq!(parsed["start"], 230710048);
        assert_eq!(parsed["end"], 230710048);
        assert_eq!(parsed["allele_string"], "A/G");
        assert_eq!(parsed["strand"], 1);
        assert_eq!(parsed["failed"], 0);
        assert_eq!(parsed["somatic"], 0);
        assert_eq!(parsed["minor_allele"], "G");
        assert_eq!(parsed["minor_allele_freq"], 0.2345);
        assert_eq!(parsed["clin_sig"], "benign");
        assert_eq!(parsed["phenotype_or_disease"], 1);
        assert_eq!(parsed["pubmed"], "12345");
        assert_eq!(parsed["AF"], "G:0.5");
        assert_eq!(parsed["AFR"], "G:0.3");
        assert_eq!(parsed["gnomADe"], "G:0.45");

        assert!(parsed.get("gnomADe_AFR").is_none());
        assert!(parsed.get("gnomADg").is_none());
        assert!(parsed.get("var_synonyms").is_none());
    }

    #[test]
    fn test_empty_variation_name() {
        let records = parse_vcf_str("", "21\t100\t.\tA\tG\t.\tPASS\t.\n", false);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].variation_name, "");
    }

    #[test]
    fn test_multiple_ids_takes_first() {
        let records = parse_vcf_str("", "21\t100\trs1;rs2;rs3\tA\tG\t.\tPASS\t.\n", false);
        assert_eq!(records[0].variation_name, "rs1");
    }

    #[test]
    fn test_monomorphic_site_skipped() {
        let records = parse_vcf_str("", "21\t100\trs1\tA\t.\t.\tPASS\t.\n", false);
        assert!(records.is_empty());
    }

    #[test]
    fn test_raw_info_value() {
        assert_eq!(
            get_raw_info_value("MA=G;MAF=0.5;CLIN_SIG=benign", "MA"),
            Some("G")
        );
        assert_eq!(
            get_raw_info_value("MA=G;MAF=0.5;CLIN_SIG=benign", "MAF"),
            Some("0.5")
        );
        assert_eq!(
            get_raw_info_value("MA=G;MAF=0.5;CLIN_SIG=benign", "PUBMED"),
            None
        );
        assert_eq!(get_raw_info_value(".", "MA"), None);
        assert_eq!(get_raw_info_value("", "MA"), None);
    }

    #[test]
    fn test_build_empty_records() {
        let dir = tempfile::tempdir().unwrap();
        let stats = build_variation_cache(&[], dir.path(), "21").unwrap();
        assert_eq!(stats.variation_count, 0);
        assert_eq!(stats.region_file_count, 0);
    }
}
