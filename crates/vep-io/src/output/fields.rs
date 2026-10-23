// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The VEP output field table shared by every text format.
//!
//! VEP renders one hash of named fields per (variant, allele, feature) and then
//! lays that hash out per format: the default format keeps thirteen fixed
//! columns and folds every other key into `Extra` as `KEY=VALUE` pairs, `--tab`
//! prints every key as its own column, VCF joins the CSQ subfields with `|`, and
//! JSON lowercases the keys. Which keys exist is decided by the command-line
//! flags in a fixed order (`@FLAG_FIELDS` in `Bio/EnsEMBL/VEP/Constants.pm:60-138`,
//! ensembl-vep release/116.2), and the `user` group (IMPACT, DISTANCE, STRAND,
//! FLAGS) is always on. This module is that table: the column names, their
//! descriptions, the flag-to-field order, the value of any field for a row, and
//! the order in which one input record's rows come out.

use std::borrow::Cow;
use std::collections::BTreeSet;

use vep_core::consequence::{Consequence, OwnSpanPositions, TranscriptConsequence};
use vep_core::spdi::VrsAllele;
use vep_core::sv_routing::{perl_split, perl_true};
use vep_core::terms::Vocabulary;
use vep_core::variant::{
    colocated_in_output_order, push_decimal, ColocatedVariant, InputVariant, VariantClass,
};
use vep_core::variation::reverse_comp;

use super::frequency::{self, FrequencyGroups};

/// The vep-rs release every output names beside the Ensembl VEP header lines:
/// `## vep-rs <version>` in the default and tab headers, `##vep-rs=<version>`
/// in the VCF header and `vep_rs_version` in every JSON object. Every crate of
/// the workspace carries the release as its version.
pub const VEP_RS_VERSION: &str = env!("CARGO_PKG_VERSION");

/// The thirteen fixed columns of the default and tab formats, in order.
pub const DEFAULT_OUTPUT_COLS: [&str; 13] = [
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
];

/// Which optional field groups are switched on. Field order across groups is
/// fixed by [`flag_fields`]; a group's fields appear only when its flag is set.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FieldOptions {
    /// `--individual`: the `IND` and `ZYG` columns (Constants.pm:69).
    pub individual: bool,
    /// `--individual_zyg`: the `ZYG` column alone (Constants.pm:70), the list of
    /// every named sample's zygosity.
    pub individual_zyg: bool,
    /// `--process_ref_homs`: every feature writes its reference allele's row, so a
    /// `--minimal` split record's `REF_ALLELE` is its first allele's trimmed
    /// reference whatever that allele equals.
    pub process_ref_homs: bool,
    pub allele_number: bool,
    pub show_ref_allele: bool,
    pub uploaded_allele: bool,
    pub include_pick: bool,
    pub variant_class: bool,
    pub minimal: bool,
    pub clinvar_somatic_classification: bool,
    pub symbol: bool,
    pub biotype: bool,
    pub canonical: bool,
    /// `--mane_select`, which `--mane` implies: `MANE` lists `MANE_Select` and
    /// `MANE_SELECT` carries the accession.
    pub mane_select: bool,
    /// `--mane` alone adds the `MANE_PLUS_CLINICAL` column.
    pub mane: bool,
    /// Set by `--mane` and by no flag of its own: `MANE` lists
    /// `MANE_Plus_Clinical` and `MANE_PLUS_CLINICAL` carries the accession.
    pub mane_plus_clinical: bool,
    pub tsl: bool,
    pub gencode_primary: bool,
    pub appris: bool,
    pub ccds: bool,
    pub protein: bool,
    pub uniprot: bool,
    pub xref_refseq: bool,
    pub gene_phenotype: bool,
    pub sift: bool,
    pub polyphen: bool,
    pub numbers: bool,
    pub domains: bool,
    pub mirna: bool,
    pub hgvs: bool,
    pub hgvsg: bool,
    /// `--spdi`: the SPDI column, the allele's notation on the rows whose
    /// Allele is the allele it describes (OutputFactory.pm:1400-1401).
    pub spdi: bool,
    /// `--ga4gh_vrs`: the GA4GH_VRS column, the allele's VRS object on the
    /// same rows, JSON-encoded in the text formats (:1403-1410).
    pub ga4gh_vrs: bool,
    pub af: bool,
    pub af_1kg: bool,
    /// `--af_gnomad`, the older spelling of `--af_gnomade`, its own flag in
    /// Ensembl VEP with the same fields (Config.pm:132, Constants.pm:120).
    pub af_gnomad: bool,
    pub af_gnomade: bool,
    pub af_gnomadg: bool,
    pub max_af: bool,
    pub check_frequency: bool,
    pub check_existing: bool,
    pub pubmed: bool,
    pub check_svs: bool,
    pub dont_skip: bool,
    pub overlaps: bool,
    /// `--nearest`: the NEAREST field lists each record's nearest protein-coding
    /// transcription start sites ([`InputVariant::nearest`]).
    pub nearest: bool,
    /// `--shift_length`: the SHIFT_LENGTH field is listed; its value is written
    /// only under `shift_3prime` (OutputFactory.pm:1695).
    pub shift_length: bool,
    /// `--use_transcript_ref`: the `GIVEN_REF` and `USED_REF` columns
    /// (Constants.pm:99), and the reference row's Allele is the allele the row
    /// was analysed against.
    pub use_transcript_ref: bool,
    /// `--shift_3prime 1`: insertions and deletions were annotated at their most
    /// 3' position on each transcript, so SHIFT_LENGTH has a value and
    /// `shift_genomic` has a shift to print.
    pub shift_3prime: bool,
    /// `--shift_genomic 1`: Location and Allele print the shifted position and
    /// allele of each row under `shift_3prime` (OutputFactory.pm:1349-1352,
    /// 1702-1706, 1992-1998).
    pub shift_genomic: bool,
    pub var_synonyms: bool,
    pub regulatory: bool,
    /// `--transcript_version`: `Feature` carries the transcript's `.<version>`.
    pub transcript_version: bool,
    /// `--gene_version`: `Gene` carries the gene's `.<version>`.
    pub gene_version: bool,
    pub ambiguity: bool,
    /// VEP's `--no_escape` (implied by `--json`): HGVSp keeps its `=` instead of
    /// the `%3D` the text formats print, and a gene symbol's whitespace stays in
    /// VCF output instead of becoming `%20`.
    pub no_escape: bool,
    /// The vocabulary consequence terms print in (`--terms`).
    pub terms: Vocabulary,
    /// Append `/<total>` to the three position columns (`--total_length`).
    pub total_length: bool,
    /// The format clears the cached positions before it rejoins the alleles of a
    /// `--minimal` split (the VCF and JSON factories, `OutputFactory/VCF.pm:305-308`,
    /// `OutputFactory/JSON.pm:217-220`), so an allele of one prints its own
    /// span's positions ([`TranscriptConsequence::own_span_positions`]); the
    /// default and tab formats rejoin first (`OutputFactory.pm:294-299`) and print
    /// the record's span.
    pub positions_of_own_span: bool,
    /// `--clin_sig_allele`: which ClinVar entries CLIN_SIG reports.
    pub clin_sig_allele: ClinSigAllele,
}

/// Which ClinVar significance CLIN_SIG reports for a co-located variant that
/// carries allele-specific entries (`--clin_sig_allele`, Config.pm:151, on by
/// default at :316).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ClinSigAllele {
    /// `--clin_sig_allele 1`: the entries of the row's own allele
    /// (OutputFactory.pm:1097-1122, 1143-1144).
    #[default]
    PerAllele,
    /// `--clin_sig_allele 0`: the variant's terms whatever the allele (:1125).
    WholeVariant,
}

/// The gnomAD exome fields of `--af_gnomad` and `--af_gnomade` (Constants.pm:120-121).
const GNOMADE_FIELDS: &[&str] = &[
    "gnomADe_AF",
    "gnomADe_AFR_AF",
    "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF",
    "gnomADe_EAS_AF",
    "gnomADe_FIN_AF",
    "gnomADe_MID_AF",
    "gnomADe_NFE_AF",
    "gnomADe_REMAINING_AF",
    "gnomADe_SAS_AF",
];

/// One flag-driven field group: the flag, its fields, and whether VEP fills
/// the fields from the variant alone (`VariationFeature_to_output_hash`), so
/// that they survive on a `--most_severe` or `--summary` row, which carries
/// no feature.
struct FieldGroup {
    flag: &'static str,
    fields: &'static [&'static str],
    vf_level: bool,
}

const fn group(flag: &'static str, fields: &'static [&'static str]) -> FieldGroup {
    FieldGroup {
        flag,
        fields,
        vf_level: false,
    }
}

const fn vf_group(flag: &'static str, fields: &'static [&'static str]) -> FieldGroup {
    FieldGroup {
        flag,
        fields,
        vf_level: true,
    }
}

/// The flag-driven field groups in VEP's order (`@FLAG_FIELDS`). The `user`
/// group is unconditional and comes first among the groups that can fire here.
const FLAG_FIELD_GROUPS: &[FieldGroup] = &[
    vf_group("individual", &["IND", "ZYG"]),
    vf_group("individual_zyg", &["ZYG"]),
    group("allele_number", &["ALLELE_NUM"]),
    group("show_ref_allele", &["REF_ALLELE"]),
    group("uploaded_allele", &["UPLOADED_ALLELE"]),
    group("user", &["IMPACT", "DISTANCE", "STRAND", "FLAGS"]),
    group("include_pick", &["PICK"]),
    vf_group("variant_class", &["VARIANT_CLASS"]),
    vf_group("minimal", &["MINIMISED"]),
    group("spdi", &["SPDI"]),
    group("ga4gh_vrs", &["GA4GH_VRS"]),
    group(
        "clinvar_somatic_classification",
        &["CLINVAR_SOMATIC_CLASSIFICATION"],
    ),
    group("symbol", &["SYMBOL", "SYMBOL_SOURCE", "HGNC_ID"]),
    group("biotype", &["BIOTYPE"]),
    group("canonical", &["CANONICAL"]),
    group("mane_select", &["MANE", "MANE_SELECT"]),
    group("mane", &["MANE", "MANE_SELECT", "MANE_PLUS_CLINICAL"]),
    group("tsl", &["TSL"]),
    group("gencode_primary", &["GENCODE_PRIMARY"]),
    group("appris", &["APPRIS"]),
    group("ccds", &["CCDS"]),
    group("protein", &["ENSP"]),
    group(
        "uniprot",
        &["SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM"],
    ),
    group("xref_refseq", &["RefSeq"]),
    group("use_transcript_ref", &["GIVEN_REF", "USED_REF"]),
    group("gene_phenotype", &["GENE_PHENO"]),
    vf_group("nearest", &["NEAREST"]),
    vf_group("ambiguity", &["AMBIGUITY"]),
    group("sift", &["SIFT"]),
    group("polyphen", &["PolyPhen"]),
    group("numbers", &["EXON", "INTRON"]),
    group("domains", &["DOMAINS"]),
    group("mirna", &["miRNA"]),
    group("hgvs", &["HGVSc", "HGVSp", "HGVS_OFFSET"]),
    group("hgvsg", &["HGVSg"]),
    group("af", &["AF"]),
    group(
        "af_1kg",
        &["AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF"],
    ),
    group("af_gnomad", GNOMADE_FIELDS),
    group("af_gnomade", GNOMADE_FIELDS),
    group(
        "af_gnomadg",
        &[
            "gnomADg_AF",
            "gnomADg_AFR_AF",
            "gnomADg_AMI_AF",
            "gnomADg_AMR_AF",
            "gnomADg_ASJ_AF",
            "gnomADg_EAS_AF",
            "gnomADg_FIN_AF",
            "gnomADg_MID_AF",
            "gnomADg_NFE_AF",
            "gnomADg_REMAINING_AF",
            "gnomADg_SAS_AF",
        ],
    ),
    group("max_af", &["MAX_AF", "MAX_AF_POPS"]),
    group("check_frequency", &["FREQS"]),
    group("check_existing", &["CLIN_SIG", "SOMATIC", "PHENO"]),
    group("pubmed", &["PUBMED"]),
    vf_group("check_svs", &["SV"]),
    vf_group("dont_skip", &["CHECK_REF"]),
    group("overlaps", &["OverlapBP", "OverlapPC"]),
    group("shift_length", &["SHIFT_LENGTH"]),
    group("var_synonyms", &["VAR_SYNONYMS"]),
    // BIOTYPE leads the regulatory group (Constants.pm:136), so `--regulatory`
    // without `--biotype` adds the column; the value is written only under
    // `--biotype` (OutputFactory.pm:1580), so the column prints `-`.
    group(
        "regulatory",
        &[
            "BIOTYPE",
            "MOTIF_NAME",
            "MOTIF_POS",
            "HIGH_INF_POS",
            "MOTIF_SCORE_CHANGE",
            "TRANSCRIPTION_FACTORS",
        ],
    ),
];

impl FieldOptions {
    fn group_active(&self, group: &str) -> bool {
        match group {
            "user" => true,
            "individual" => self.individual,
            "individual_zyg" => self.individual_zyg,
            "allele_number" => self.allele_number,
            "show_ref_allele" => self.show_ref_allele,
            "uploaded_allele" => self.uploaded_allele,
            "include_pick" => self.include_pick,
            "variant_class" => self.variant_class,
            "minimal" => self.minimal,
            "clinvar_somatic_classification" => self.clinvar_somatic_classification,
            "symbol" => self.symbol,
            "biotype" => self.biotype,
            "canonical" => self.canonical,
            "mane_select" => self.mane_select,
            "mane" => self.mane,
            "tsl" => self.tsl,
            "gencode_primary" => self.gencode_primary,
            "appris" => self.appris,
            "ccds" => self.ccds,
            "protein" => self.protein,
            "uniprot" => self.uniprot,
            "xref_refseq" => self.xref_refseq,
            "gene_phenotype" => self.gene_phenotype,
            "ambiguity" => self.ambiguity,
            "sift" => self.sift,
            "polyphen" => self.polyphen,
            "numbers" => self.numbers,
            "domains" => self.domains,
            "mirna" => self.mirna,
            "hgvs" => self.hgvs,
            "hgvsg" => self.hgvsg,
            "spdi" => self.spdi,
            "ga4gh_vrs" => self.ga4gh_vrs,
            "af" => self.af,
            "af_1kg" => self.af_1kg,
            "af_gnomad" => self.af_gnomad,
            "af_gnomade" => self.af_gnomade,
            "af_gnomadg" => self.af_gnomadg,
            "max_af" => self.max_af,
            "check_frequency" => self.check_frequency,
            "check_existing" => self.check_existing,
            "pubmed" => self.pubmed,
            "check_svs" => self.check_svs,
            "dont_skip" => self.dont_skip,
            "overlaps" => self.overlaps,
            "nearest" => self.nearest,
            "shift_length" => self.shift_length,
            "use_transcript_ref" => self.use_transcript_ref,
            "var_synonyms" => self.var_synonyms,
            "regulatory" => self.regulatory,
            _ => false,
        }
    }

    /// Whether a row's `DISTANCE` is printed: VEP sets it only when the rendered
    /// consequence list names `upstream_gene_variant` or `downstream_gene_variant`
    /// (`OutputFactory.pm`, `%DISTANCE_CONS` tested against `Consequence`), so the
    /// display vocabulary, which renders them `UPSTREAM` and `DOWNSTREAM`, never
    /// prints it; the SO and NCBI vocabularies render both terms as themselves.
    pub fn distance_shown(&self, tc: &TranscriptConsequence) -> bool {
        tc.consequences.iter().any(|c| {
            matches!(
                self.terms.term(*c),
                "upstream_gene_variant" | "downstream_gene_variant"
            )
        })
    }

    /// Whether the unconditional `user` group prints as it does under the SO
    /// vocabulary: `DISTANCE` on every row that carries one.
    pub fn user_group_is_plain(&self) -> bool {
        self.terms != Vocabulary::Display
    }

    /// True when any group beyond the unconditional `user` group is on, i.e.
    /// when a row can carry a field that is not IMPACT, DISTANCE, STRAND or FLAGS.
    pub fn any_optional_group(&self) -> bool {
        FLAG_FIELD_GROUPS
            .iter()
            .any(|g| g.flag != "user" && self.group_active(g.flag))
    }

    /// The frequency groups the configuration reports.
    pub fn frequency_groups(&self) -> FrequencyGroups {
        FrequencyGroups {
            af: self.af,
            af_1kg: self.af_1kg,
            af_gnomad: self.af_gnomad,
            af_gnomade: self.af_gnomade,
            af_gnomadg: self.af_gnomadg,
            max_af: self.max_af,
        }
    }
}

/// The row-invariant part of a configuration's `Extra` column, computed once
/// per run so the per-row writer neither walks the flag table nor rebuilds the
/// flag-field list.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtraFieldsPlan {
    pub options: FieldOptions,
    /// [`flag_fields`] of `options`.
    pub flag_fields: Vec<&'static str>,
    /// [`FieldOptions::any_optional_group`] of `options`.
    pub any_optional_group: bool,
}

impl ExtraFieldsPlan {
    pub fn new(options: FieldOptions) -> Self {
        let flag_fields = flag_fields(&options);
        let any_optional_group = options.any_optional_group();
        Self {
            options,
            flag_fields,
            any_optional_group,
        }
    }
}

/// Groups whose flag adds a column without guarding the value: VEP writes
/// `OverlapBP` and `OverlapPC` on every structural row that overlaps its
/// feature and `CHECK_REF` on every failed check, flag or no flag
/// (`OutputFactory.pm` `StructuralVariationOverlapAllele_to_output_hash`,
/// `VariationFeature_to_output_hash`).
const UNGUARDED_GROUPS: [&str; 2] = ["overlaps", "dont_skip"];

/// Whether a configuration produces the field: a flag field when one of its
/// groups is on or unguarded (`UNGUARDED_GROUPS`), every other name always.
/// VEP sets a guarded flag field's key in the row hash only under its flag, so
/// a `--fields` list naming an inactive field prints `-` or an empty subfield
/// for it.
pub fn field_is_produced(field: &str, options: &FieldOptions) -> bool {
    let mut listed = false;
    for group in FLAG_FIELD_GROUPS {
        if group.fields.contains(&field) {
            if UNGUARDED_GROUPS.contains(&group.flag) || options.group_active(group.flag) {
                return true;
            }
            listed = true;
        }
    }
    !listed
}

/// The optional fields a configuration emits, in VEP's order, each name once.
pub fn flag_fields(options: &FieldOptions) -> Vec<&'static str> {
    fields_of_groups(options, |_| true)
}

/// The active fields VEP fills from the variant alone, in VEP's order: the
/// keys a row without a feature (a `--most_severe` or `--summary` row) carries.
pub fn vf_level_flag_fields(options: &FieldOptions) -> Vec<&'static str> {
    fields_of_groups(options, |g| g.vf_level)
}

fn fields_of_groups(
    options: &FieldOptions,
    keep: impl Fn(&FieldGroup) -> bool,
) -> Vec<&'static str> {
    let mut seen: BTreeSet<&str> = BTreeSet::new();
    let mut out = Vec::new();
    for group in FLAG_FIELD_GROUPS {
        if !options.group_active(group.flag) || !keep(group) {
            continue;
        }
        for f in group.fields {
            if seen.insert(f) {
                out.push(*f);
            }
        }
    }
    out
}

/// Whether `VariationFeature_to_output_hash` sets `field` from the record alone
/// (OutputFactory.pm:883-1025 at release/116.2), so that a `--most_severe` or
/// `--summary` row prints it: the identity and location columns, and the field
/// of every `vf_group` whether or not its flag is on.
fn is_vf_level_field(field: &str) -> bool {
    matches!(field, "Uploaded_variation" | "Location")
        || FLAG_FIELD_GROUPS
            .iter()
            .any(|group| group.vf_level && group.fields.contains(&field))
}

/// VEP's description of a field, for the `## Column descriptions:` and
/// `## Extra column keys:` header blocks.
pub fn field_description(field: &str) -> Option<&'static str> {
    FIELD_DESCRIPTIONS
        .iter()
        .find(|(k, _)| *k == field)
        .map(|(_, d)| *d)
}

/// `%FIELD_DESCRIPTIONS` (Constants.pm:142-266), every key with its text.
const FIELD_DESCRIPTIONS: &[(&str, &str)] = &[
    ("Uploaded_variation", "Identifier of uploaded variant"),
    ("ID", "Identifier of uploaded variant"),
    (
        "Location",
        "Location of variant in standard coordinate format (chr:start or chr:start-end)",
    ),
    ("Allele", "The variant allele used to calculate the consequence"),
    ("Gene", "Stable ID of affected gene"),
    ("Feature", "Stable ID of feature"),
    (
        "Feature_type",
        "Type of feature - Transcript, RegulatoryFeature or MotifFeature",
    ),
    ("Consequence", "Consequence type"),
    ("cDNA_position", "Relative position of base pair in cDNA sequence"),
    ("CDS_position", "Relative position of base pair in coding sequence"),
    ("Protein_position", "Relative position of amino acid in protein"),
    ("Amino_acids", "Reference and variant amino acids"),
    ("Codons", "Reference and variant codon sequence"),
    ("Existing_variation", "Identifier(s) of co-located known variants"),
    ("IMPACT", "Subjective impact classification of consequence type"),
    ("CANONICAL", "Indicates if transcript is canonical for this gene"),
    (
        "MANE",
        "MANE (Matched Annotation from NCBI and EMBL-EBI) set(s) the transcript belongs to",
    ),
    (
        "MANE_SELECT",
        "MANE Select (Matched Annotation from NCBI and EMBL-EBI) Transcript",
    ),
    (
        "MANE_PLUS_CLINICAL",
        "MANE Plus Clinical (Matched Annotation from NCBI and EMBL-EBI) Transcript",
    ),
    ("TSL", "Transcript support level"),
    (
        "APPRIS",
        "Annotates alternatively spliced transcripts as primary or alternate based on a range of computational methods",
    ),
    ("CCDS", "Indicates if transcript is a CCDS transcript"),
    ("SYMBOL", "Gene symbol (e.g. HGNC)"),
    ("SYMBOL_SOURCE", "Source of gene symbol"),
    ("SOURCE", "Source of transcript"),
    ("HGNC_ID", "Stable identifer of HGNC gene symbol"),
    ("ENSP", "Protein identifer"),
    ("FLAGS", "Transcript quality flags"),
    ("SWISSPROT", "UniProtKB/Swiss-Prot accession"),
    ("TREMBL", "UniProtKB/TrEMBL accession"),
    ("UNIPARC", "UniParc accession"),
    ("UNIPROT_ISOFORM", "Direct mappings to UniProtKB isoforms"),
    ("miRNA", "SO terms of overlapped miRNA secondary structure feature(s)"),
    ("HGVSc", "HGVS coding sequence name"),
    ("HGVSp", "HGVS protein sequence name"),
    ("HGVSg", "HGVS genomic sequence name"),
    ("SPDI", "Genomic SPDI notation"),
    (
        "GA4GH_VRS",
        "GA4GH Variation Representation Specification",
    ),
    ("SIFT", "SIFT prediction and/or score"),
    ("PolyPhen", "PolyPhen prediction and/or score"),
    ("EXON", "Exon number(s) / total"),
    ("INTRON", "Intron number(s) / total"),
    (
        "DOMAINS",
        "The source and identifer of any overlapping protein domains",
    ),
    (
        "MOTIF_NAME",
        "The stable identifier of a transcription factor binding profile (TFBP) aligned at this position",
    ),
    (
        "MOTIF_POS",
        "The relative position of the variation in the aligned TFBP",
    ),
    (
        "HIGH_INF_POS",
        "A flag indicating if the variant falls in a high information position of the TFBP",
    ),
    (
        "MOTIF_SCORE_CHANGE",
        "The difference in motif score of the reference and variant sequences for the TFBP",
    ),
    (
        "TRANSCRIPTION_FACTORS",
        "List of transcription factors which bind to the transcription factor binding profile",
    ),
    (
        "CELL_TYPE",
        "List of cell types and classifications for regulatory feature",
    ),
    (
        "AF",
        "Frequency of existing variant in 1000 Genomes combined population",
    ),
    (
        "AFR_AF",
        "Frequency of existing variant in 1000 Genomes combined African population",
    ),
    (
        "AMR_AF",
        "Frequency of existing variant in 1000 Genomes combined American population",
    ),
    (
        "ASN_AF",
        "Frequency of existing variant in 1000 Genomes combined Asian population",
    ),
    (
        "EAS_AF",
        "Frequency of existing variant in 1000 Genomes combined East Asian population",
    ),
    (
        "EUR_AF",
        "Frequency of existing variant in 1000 Genomes combined European population",
    ),
    (
        "SAS_AF",
        "Frequency of existing variant in 1000 Genomes combined South Asian population",
    ),
    (
        "AA_AF",
        "Frequency of existing variant in NHLBI-ESP African American population",
    ),
    (
        "EA_AF",
        "Frequency of existing variant in NHLBI-ESP European American population",
    ),
    (
        "ExAC_AF",
        "Frequency of existing variant in ExAC combined population",
    ),
    (
        "ExAC_Adj_AF",
        "Adjusted frequency of existing variant in ExAC combined population",
    ),
    (
        "ExAC_AFR_AF",
        "Frequency of existing variant in ExAC African/American population",
    ),
    (
        "ExAC_AMR_AF",
        "Frequency of existing variant in ExAC American population",
    ),
    (
        "ExAC_EAS_AF",
        "Frequency of existing variant in ExAC East Asian population",
    ),
    (
        "ExAC_FIN_AF",
        "Frequency of existing variant in ExAC Finnish population",
    ),
    (
        "ExAC_NFE_AF",
        "Frequency of existing variant in ExAC Non-Finnish European population",
    ),
    (
        "ExAC_OTH_AF",
        "Frequency of existing variant in ExAC other combined populations",
    ),
    (
        "ExAC_SAS_AF",
        "Frequency of existing variant in ExAC South Asian population",
    ),
    (
        "gnomADe_AF",
        "Frequency of existing variant in gnomAD exomes combined population",
    ),
    (
        "gnomADe_AFR_AF",
        "Frequency of existing variant in gnomAD exomes African/American population",
    ),
    (
        "gnomADe_AMR_AF",
        "Frequency of existing variant in gnomAD exomes American population",
    ),
    (
        "gnomADe_ASJ_AF",
        "Frequency of existing variant in gnomAD exomes Ashkenazi Jewish population",
    ),
    (
        "gnomADe_EAS_AF",
        "Frequency of existing variant in gnomAD exomes East Asian population",
    ),
    (
        "gnomADe_FIN_AF",
        "Frequency of existing variant in gnomAD exomes Finnish population",
    ),
    (
        "gnomADe_MID_AF",
        "Frequency of existing variant in gnomAD exomes Mid-eastern population",
    ),
    (
        "gnomADe_NFE_AF",
        "Frequency of existing variant in gnomAD exomes Non-Finnish European population",
    ),
    (
        "gnomADe_REMAINING_AF",
        "Frequency of existing variant in gnomAD exomes remaining combined populations",
    ),
    (
        "gnomADe_SAS_AF",
        "Frequency of existing variant in gnomAD exomes South Asian population",
    ),
    (
        "gnomADg_AF",
        "Frequency of existing variant in gnomAD genomes combined population",
    ),
    (
        "gnomADg_AFR_AF",
        "Frequency of existing variant in gnomAD genomes African/American population",
    ),
    (
        "gnomADg_AMI_AF",
        "Frequency of existing variant in gnomAD genomes Amish population",
    ),
    (
        "gnomADg_AMR_AF",
        "Frequency of existing variant in gnomAD genomes American population",
    ),
    (
        "gnomADg_ASJ_AF",
        "Frequency of existing variant in gnomAD genomes Ashkenazi Jewish population",
    ),
    (
        "gnomADg_EAS_AF",
        "Frequency of existing variant in gnomAD genomes East Asian population",
    ),
    (
        "gnomADg_FIN_AF",
        "Frequency of existing variant in gnomAD genomes Finnish population",
    ),
    (
        "gnomADg_MID_AF",
        "Frequency of existing variant in gnomAD genomes Mid-eastern population",
    ),
    (
        "gnomADg_NFE_AF",
        "Frequency of existing variant in gnomAD genomes Non-Finnish European population",
    ),
    (
        "gnomADg_REMAINING_AF",
        "Frequency of existing variant in gnomAD genomes remaining combined populations",
    ),
    (
        "gnomADg_SAS_AF",
        "Frequency of existing variant in gnomAD genomes South Asian population",
    ),
    (
        "MAX_AF",
        "Maximum observed allele frequency in 1000 Genomes, ESP and ExAC/gnomAD",
    ),
    (
        "MAX_AF_POPS",
        "Populations in which maximum allele frequency was observed",
    ),
    ("SV", "IDs of overlapping structural variants"),
    ("FREQS", "Frequencies of overlapping variants used in filtering"),
    ("DISTANCE", "Shortest distance from variant to transcript"),
    ("CLIN_SIG", "ClinVar clinical significance of the dbSNP variant"),
    (
        "CLINVAR_SOMATIC_CLASSIFICATION",
        "ClinVar somatic classification of the dbSNP variant",
    ),
    ("BIOTYPE", "Biotype of transcript or regulatory feature"),
    ("PUBMED", "Pubmed ID(s) of publications that cite existing variant"),
    (
        "VAR_SYNONYMS",
        "List of known variation synonyms and their sources",
    ),
    (
        "ALLELE_NUM",
        "Allele number from input; 0 is reference, 1 is first alternate etc",
    ),
    ("REF_ALLELE", "Reference allele"),
    ("UPLOADED_ALLELE", "The variant allele uploaded"),
    ("AMBIGUITY", "Allele ambiguity code"),
    ("STRAND", "Strand of the feature (1/-1)"),
    (
        "PICK",
        "Indicates if this consequence has been picked as the most severe",
    ),
    ("SOMATIC", "Somatic status of existing variant"),
    ("REFSEQ_MATCH", "RefSeq transcript match status"),
    (
        "REFSEQ_OFFSET",
        "HGVS adjustment length required due to mismatch between RefSeq transcript and the reference genome",
    ),
    ("BAM_EDIT", "Indicates success or failure of edit using BAM file"),
    ("IND", "Individual name"),
    ("ZYG", "Zygosity of individual genotype at this locus"),
    ("GIVEN_REF", "Reference allele from input"),
    ("USED_REF", "Reference allele as used to get consequences"),
    ("VARIANT_CLASS", "SO variant class"),
    (
        "PHENO",
        "Indicates if existing variant(s) is associated with a phenotype, disease or trait; multiple values correspond to multiple variants",
    ),
    (
        "GENE_PHENO",
        "Indicates if gene is associated with a phenotype, disease or trait",
    ),
    (
        "NEAREST",
        "Identifier(s) of nearest transcription start site",
    ),
    (
        "MINIMISED",
        "Alleles in this variant have been converted to minimal representation before consequence calculation",
    ),
    (
        "HGVS_OFFSET",
        "Indicates by how many bases the HGVS notations for this variant have been shifted",
    ),
    (
        "OverlapBP",
        "Number of base pairs overlapping with the corresponding structural variation feature",
    ),
    (
        "OverlapPC",
        "Percentage of corresponding structural variation feature overlapped by the given input",
    ),
    (
        "CHECK_REF",
        "Reports variants where the input reference does not match the expected reference",
    ),
    ("GENCODE_PRIMARY", "Reports if transcript is GENCODE primary"),
    (
        "SHIFT_LENGTH",
        "Reports the number of bases the insertion or deletion has been shifted relative to the underlying transcript due to right alignment before consequence calculation",
    ),
];

/// The shift `--shift_genomic 1` prints on a row and the strand it moved along:
/// a transcript row's own shift along its transcript (OutputFactory.pm:1702-1706),
/// an intergenic row's shift along the forward strand (:1992-1998). `None` when
/// the flags are off or the row moved nothing, so Location and Allele print as
/// the input placed them.
fn genomic_shift(
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> Option<(u64, i8)> {
    if !(options.shift_3prime && options.shift_genomic) {
        return None;
    }
    let (shift, strand) = match tc {
        Some(tc) => (tc.shift_length, tc.strand),
        None => (variant.intergenic_shift?, 1),
    };
    (shift > 0).then_some((shift, strand))
}

/// Appends the row's `Allele`: the allele as the record gives it; under
/// `--shift_3prime 1 --shift_genomic 1`, rotated along the row's shift
/// (OutputFactory.pm:1349-1352); under `--shift_3prime 1` alone, the sequence
/// the row's codon computation rotated ([`TranscriptConsequence::shifted_allele`],
/// :1347). A reference row's is the reference allele's sequence (:1347):
/// [`TranscriptConsequence::used_ref`] when the row has one.
fn write_row_allele(
    out: &mut Vec<u8>,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) {
    if variant.reference_allele_row {
        if let Some(used_ref) = tc.and_then(|x| x.used_ref.as_deref()) {
            out.extend_from_slice(used_ref.as_bytes());
            return;
        }
    }
    if let Some((shift, strand)) = genomic_shift(variant, tc, options) {
        out.extend_from_slice(variant.display_allele_shifted(shift, strand).as_bytes());
    } else if let Some(allele) = tc.and_then(|x| x.shifted_allele.as_deref()) {
        out.extend_from_slice(allele.as_bytes());
    } else {
        variant.write_display_allele(out);
    }
}

/// [`write_row_allele`] as a `String`.
fn row_allele(
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> String {
    let mut out = Vec::new();
    write_row_allele(&mut out, variant, tc, options);
    written_text(out)
}

/// The text a field writer appended: `str` slices and ASCII, so UTF-8; a byte
/// that is not is carried through as the replacement character.
fn written_text(bytes: Vec<u8>) -> String {
    String::from_utf8(bytes).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned())
}

/// Rewrites the ASCII control bytes (`0x00` to `0x1F`, the tab and line end
/// among them, and `0x7F`) a writer appended to `out` from `start` on as
/// `%XX`, the percent-encoding of the byte in upper-case hex. A value read
/// from an input file, an identifier or a chromosome name among them, can
/// carry a byte the default and tab formats use as a separator, a tab or a
/// line end, or a byte no text format shows; written as is, it would split the
/// row or the column, so every value of those formats passes through here as
/// it is written. A value with no control byte is left untouched.
pub fn percent_encode_control_bytes(out: &mut Vec<u8>, start: usize) {
    if !out[start..].iter().any(u8::is_ascii_control) {
        return;
    }
    let value = out.split_off(start);
    for &b in &value {
        if b.is_ascii_control() {
            out.extend_from_slice(format!("%{b:02X}").as_bytes());
        } else {
            out.push(b);
        }
    }
}

/// The SPDI and GA4GH VRS allele the row prints: VEP keeps the record's
/// notations by allele and reads them by the row's Allele
/// (OutputFactory.pm:1396-1400), so a row prints its own allele's, a row whose
/// Allele a genomic shift rotated onto a sibling allele of the record prints
/// the sibling's, and a row whose rotated Allele is no allele of the record
/// carries neither column.
fn row_spdi<'a>(
    variant: &'a InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> Option<(&'a str, Option<&'a VrsAllele>)> {
    let allele = row_allele(variant, tc, options);
    if allele == variant.display_allele() {
        return variant
            .spdi
            .as_deref()
            .map(|spdi| (spdi, variant.ga4gh_vrs.as_ref()));
    }
    variant
        .sibling_spdi
        .iter()
        .find(|sibling| sibling.allele == allele)
        .map(|sibling| (sibling.spdi.as_str(), sibling.ga4gh_vrs.as_ref()))
}

/// The GA4GH VRS allele the row prints (`row_spdi`), for the JSON format to
/// nest as an object.
pub fn row_ga4gh_vrs<'a>(
    variant: &'a InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> Option<&'a VrsAllele> {
    row_spdi(variant, tc, options).and_then(|(_, vrs)| vrs)
}

/// One output row: an allele of the record and the transcript consequence it is
/// reported against, or `None` for the record's intergenic row or its `--summary`
/// or `--most_severe` row (the allele carrying [`InputVariant::summary_terms`]).
pub type Row<'a> = (&'a InputVariant, Option<&'a TranscriptConsequence>);

/// The rows one input record produces, in VEP's order: transcript-major (by
/// transcript id) and allele-minor (the reference row of `--process_ref_homs`
/// first, then ALT order) across every allele split from the record, then one
/// intergenic row per allele that overlapped nothing, or the one `--summary` or
/// `--most_severe` row the record's first allele carries.
/// The alleles of a `--minimal` split sort by the keys Ensembl VEP's rejoin
/// stores their transcript variations under ([`InputVariant::rejoin_key_prefix`]
/// before the transcript id), so an allele whose key prefix sorts before the
/// letters of a transcript id prints every one of its rows first.
///
/// An allele whose chromosome had no transcripts at all (no `most_severe_consequence`
/// and no summary row) produces no row, as VEP writes nothing for it.
pub fn record_rows<'a>(record: &[&'a InputVariant]) -> Vec<Row<'a>> {
    let mut rows = Vec::new();
    record_rows_into(record, &mut rows);
    rows
}

/// [`record_rows`] into a caller-owned vector, cleared first, so a writer that
/// renders many records reuses one allocation.
pub fn record_rows_into<'a>(record: &[&'a InputVariant], rows: &mut Vec<Row<'a>>) {
    rows.clear();
    for variant in record {
        for tc in &variant.transcript_consequences {
            rows.push((variant, Some(tc)));
        }
    }
    // Stable, so rows with one transcript and allele keep their consequence order.
    // The reference allele's row of `--process_ref_homs` leads each feature's
    // alleles (BaseVariationFeatureOverlap.pm:363-373 lists it first).
    rows.sort_by(|a, b| {
        transcript_order(a, b)
            .then(b.0.reference_allele_row.cmp(&a.0.reference_allele_row))
            .then(a.0.allele_index.cmp(&b.0.allele_index))
    });
    // An allele with no feature row writes its intergenic row, or the record's
    // `--summary` or `--most_severe` row, which stands even with no term
    // (OutputFactory.pm:454-458).
    for variant in record {
        if variant.transcript_consequences.is_empty()
            && (variant.most_severe_consequence.is_some() || variant.summary_terms.is_some())
        {
            rows.push((variant, None));
        }
    }
}

/// The order of two feature rows' sort keys ([`rejoin_key`]): by the loader's
/// ranks when both rows carry one that orders as the key does
/// ([`orderable_rank`]), else by the keys themselves.
fn transcript_order(a: &Row<'_>, b: &Row<'_>) -> std::cmp::Ordering {
    match (orderable_rank(a), orderable_rank(b)) {
        (Some(rank_a), Some(rank_b)) => rank_a.cmp(&rank_b),
        _ => rejoin_key(a).cmp(rejoin_key(b)),
    }
}

/// A feature row's [`TranscriptConsequence::id_rank`] when it orders the row as
/// its [`rejoin_key`] does: the ranks of one chromosome's list order as the ids
/// do and equal ids share a rank. `None` for a rank of 0, a consequence no
/// loader ranked; for an allele with a rejoin prefix, which sorts before the
/// id; and for a breakend, whose consequences include its mate's rows, ranked
/// on the mate's chromosome.
fn orderable_rank(row: &Row<'_>) -> Option<u32> {
    let (variant, tc) = (row.0, row.1.expect("feature row"));
    let ranked_with_its_neighbours = tc.id_rank != 0
        && variant.rejoin_key_prefix.is_none()
        && variant.variant_class != VariantClass::Translocation;
    ranked_with_its_neighbours.then_some(tc.id_rank)
}

/// The key Ensembl VEP sorts a row's transcript variation by (VariationFeature.pm:666):
/// the transcript id, behind the allele's rejoin prefix when it has one.
fn rejoin_key<'a>(row: &Row<'a>) -> impl Iterator<Item = u8> + 'a {
    let (variant, tc) = (row.0, row.1.expect("feature row"));
    variant
        .rejoin_key_prefix
        .as_deref()
        .unwrap_or("")
        .bytes()
        .chain(tc.transcript_id.bytes())
}

/// Consequence terms of a feature row in the configured vocabulary, joined by
/// `,`, a label every time an SO term renders to it.
pub fn consequence_string(tc: &TranscriptConsequence, options: &FieldOptions) -> String {
    rendered_terms(&tc.consequences, options.terms, false).join(",")
}

/// `terms` in the vocabulary `--terms` names, one label per SO term, or with
/// `once` each rendered label the first time it appears: Ensembl VEP's
/// `summary_only` collects a `--summary` row's terms through a seen-set after
/// rendering them (`OutputFactory.pm:426-458`), so two SO terms sharing a display
/// label print it one time there.
pub fn rendered_terms(
    terms: &[Consequence],
    vocabulary: Vocabulary,
    once: bool,
) -> Vec<&'static str> {
    let mut out = Vec::with_capacity(terms.len());
    for &c in terms {
        let term = vocabulary.term(c);
        if once && out.contains(&term) {
            continue;
        }
        out.push(term);
    }
    out
}

/// The version `--transcript_version` and `--gene_version` append to a stable
/// id as `.<version>` (OutputFactory.pm:1465 and :1469 at release/116.2):
/// `Some` only for a non-empty id with a known, non-zero version that does not
/// already end in `.<digits>`.
pub fn version_suffix(id: &str, version: Option<u32>) -> Option<u32> {
    let version = version.filter(|v| *v != 0)?;
    if id.is_empty() {
        return None;
    }
    let already_versioned = id
        .rsplit_once('.')
        .is_some_and(|(_, tail)| !tail.is_empty() && tail.bytes().all(|b| b.is_ascii_digit()));
    (!already_versioned).then_some(version)
}

/// Appends `id`, then `.<version>` when `flag` is on and [`version_suffix`]
/// yields one.
pub fn write_versioned_id(out: &mut Vec<u8>, id: &str, version: Option<u32>, flag: bool) {
    out.extend_from_slice(id.as_bytes());
    if !flag {
        return;
    }
    if let Some(version) = version_suffix(id, version) {
        out.push(b'.');
        push_decimal(out, u64::from(version));
    }
}

/// Appends a cache identifier field as Ensembl VEP reads it: the cache's `-`
/// stands for no value and sets nothing (OutputFactory.pm:1525-1572 at
/// release/116.2), and a comma-separated list of accessions is split as Perl's
/// `split` splits it ([`perl_split`]: a trailing comma adds no element), which
/// the formats join with their own separator.
fn write_identifier_list(out: &mut Vec<u8>, value: Option<&str>, list_sep: &str) {
    match value {
        None | Some("-") => {}
        Some(v) => push_joined(out, perl_split(v, &[',']), list_sep),
    }
}

/// Appends a cache identifier field that is one value: the `-` sentinel
/// prints nothing.
fn write_identifier(out: &mut Vec<u8>, value: Option<&str>) {
    match value {
        None | Some("-") => {}
        Some(v) => out.extend_from_slice(v.as_bytes()),
    }
}

/// Appends `n` in decimal with a leading `-` when negative.
fn push_signed(out: &mut Vec<u8>, n: i64) {
    if n < 0 {
        out.push(b'-');
    }
    push_decimal(out, n.unsigned_abs());
}

/// Appends `value` when the row carries one.
fn push_opt(out: &mut Vec<u8>, value: Option<&str>) {
    if let Some(v) = value {
        out.extend_from_slice(v.as_bytes());
    }
}

/// Appends `values` joined by `sep`.
fn push_joined<'a>(out: &mut Vec<u8>, values: impl IntoIterator<Item = &'a str>, sep: &str) {
    for (i, v) in values.into_iter().enumerate() {
        if i > 0 {
            out.extend_from_slice(sep.as_bytes());
        }
        out.extend_from_slice(v.as_bytes());
    }
}

/// Appends `terms` in the vocabulary `--terms` names, one label per SO term,
/// `,`-joined: a feature row's Consequence column. A summary row lists each
/// rendered label once and so goes through [`write_field_value`] instead.
pub fn write_terms(out: &mut Vec<u8>, terms: &[Consequence], vocabulary: Vocabulary) {
    push_joined(out, terms.iter().map(|&term| vocabulary.term(term)), ",");
}

/// The value of one named field for one row, or an empty string when the field
/// is absent on this row. `list_sep` joins list-valued fields (`,` in the default
/// and tab formats, `&` inside a VCF CSQ entry). Field names not in the table are
/// looked up in the row's plugin data.
///
/// The bytes [`write_field_value`] appends, collected: every format reads the
/// one field table.
pub fn field_value(
    field: &str,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
    list_sep: &str,
) -> String {
    let mut out = Vec::new();
    write_field_value(&mut out, field, variant, tc, options, list_sep);
    written_text(out)
}

/// Appends [`field_value`] to `out`, nothing when the field is absent on this
/// row.
///
/// A row carrying the record's [`InputVariant::summary_terms`] is Ensembl VEP's
/// `summary_only` row (OutputFactory.pm:426-458): the record's hash with the
/// summary's `Consequence`, so it has the variant's identity, location and
/// variant-level fields and nothing else, no Allele, IMPACT, Existing_variation
/// or plugin key included.
pub fn write_field_value(
    out: &mut Vec<u8>,
    field: &str,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
    list_sep: &str,
) {
    if let (None, Some(terms)) = (tc, &variant.summary_terms) {
        match field {
            "Consequence" => {
                let rendered = rendered_terms(terms, options.terms, true);
                if rendered.is_empty() {
                    out.push(b'-');
                } else {
                    push_joined(out, rendered, list_sep);
                }
            }
            // The record hash places the row at the input span; the
            // `--shift_genomic 1` move is written by the feature and intergenic
            // row builders alone (OutputFactory.pm:1702-1706, 1992-1998).
            "Location" => variant.write_location(out),
            _ if is_vf_level_field(field) => {
                write_row_field_value(out, field, variant, None, options, list_sep)
            }
            _ => {}
        }
        return;
    }
    write_row_field_value(out, field, variant, tc, options, list_sep)
}

/// [`write_field_value`] of a row against a feature or an allele's intergenic row.
fn write_row_field_value(
    out: &mut Vec<u8>,
    field: &str,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
    list_sep: &str,
) {
    match field {
        "Uploaded_variation" => variant.write_uploaded_variation(out),
        "Location" => match genomic_shift(variant, tc, options) {
            Some((shift, strand)) => {
                out.extend_from_slice(variant.location_shifted(shift, strand).as_bytes())
            }
            None => variant.write_location(out),
        },
        "Allele" => write_row_allele(out, variant, tc, options),
        "Consequence" => match tc {
            Some(x) => write_terms(out, &x.consequences, options.terms),
            None => {
                let consequence = variant
                    .most_severe_consequence
                    .unwrap_or(Consequence::IntergenicVariant);
                out.extend_from_slice(options.terms.term(consequence).as_bytes());
            }
        },
        "IMPACT" => {
            let impact = match tc {
                Some(x) => x.impact,
                None => variant
                    .most_severe_consequence
                    .unwrap_or(Consequence::IntergenicVariant)
                    .impact(),
            };
            out.extend_from_slice(impact.as_str().as_bytes());
        }
        "SYMBOL" => write_identifier(out, tc.and_then(|x| x.gene_symbol.as_deref())),
        "Gene" => {
            if let Some(x) = tc {
                write_versioned_id(out, &x.gene_id, x.gene_version, options.gene_version);
            }
        }
        "Feature_type" => push_opt(out, tc.map(|x| x.feature_type.as_str())),
        "Feature" => {
            if let Some(x) = tc {
                write_versioned_id(
                    out,
                    &x.transcript_id,
                    x.transcript_version,
                    options.transcript_version,
                );
            }
        }
        "BIOTYPE" if !options.biotype => {}
        "BIOTYPE" => push_opt(out, tc.and_then(|x| x.biotype.as_deref())),
        "EXON" => push_opt(out, tc.and_then(|x| x.exon.as_deref())),
        "INTRON" => push_opt(out, tc.and_then(|x| x.intron.as_deref())),
        "HGVSc" | "HGVSp" if !options.hgvs => {}
        "HGVSc" => push_opt(out, tc.and_then(|x| x.hgvsc.as_deref())),
        "HGVSp" => {
            if let Some(v) = tc.and_then(|x| x.hgvsp.as_deref()) {
                if options.no_escape {
                    out.extend_from_slice(v.as_bytes());
                } else {
                    push_joined(out, v.split('='), "%3D");
                }
            }
        }
        "cDNA_position" => write_position_column(
            out,
            tc.and_then(|x| match own_span_positions(x, options) {
                Some(own) => own.cdna.as_deref(),
                None => x.cdna_position.as_deref(),
            }),
            tc.and_then(|x| x.cdna_total)
                .filter(|_| options.total_length),
        ),
        "CDS_position" => write_position_column(
            out,
            tc.and_then(|x| match own_span_positions(x, options) {
                Some(own) => own.cds.as_deref(),
                None => x.cds_position.as_deref(),
            }),
            tc.and_then(|x| x.cds_total)
                .filter(|_| options.total_length),
        ),
        "Protein_position" => write_position_column(
            out,
            tc.and_then(|x| match own_span_positions(x, options) {
                Some(own) => own.protein.as_deref(),
                None => x.protein_position.as_deref(),
            }),
            tc.and_then(|x| x.protein_total)
                .filter(|_| options.total_length),
        ),
        "Amino_acids" => push_opt(out, tc.and_then(|x| x.amino_acids.as_deref())),
        "Codons" => push_opt(out, tc.and_then(|x| x.codons.as_deref())),
        "Existing_variation" => push_joined(
            out,
            variant.existing_variation.iter().map(String::as_str),
            list_sep,
        ),
        "ALLELE_NUM" => push_decimal(
            out,
            u64::from(
                variant
                    .allele_num
                    .unwrap_or(variant.allele_index as u32 + 1),
            ),
        ),
        "IND" => push_opt(out, variant.individual.as_deref()),
        // The `--individual_zyg` list, every named sample's `sample:ZYG`
        // (OutputFactory.pm:981-987), joined as the format joins a list; else the
        // copy's own zygosity under `--individual` (:974-976).
        "ZYG" => match &variant.sample_zygosities {
            Some(samples) => out.extend_from_slice(samples.joined(list_sep).as_bytes()),
            None => push_opt(out, variant.zygosity().as_deref()),
        },
        // A structural row carries neither `REF_ALLELE` nor `UPLOADED_ALLELE`
        // (OutputFactory.pm `BaseStructuralVariationOverlapAllele_to_output_hash`).
        "REF_ALLELE" if variant.is_structural => {}
        "REF_ALLELE" => out.extend_from_slice(feature_reference(variant, options).as_bytes()),
        "UPLOADED_ALLELE" => push_opt(out, variant.uploaded_allele(options.minimal).as_deref()),
        "DISTANCE" => {
            if let Some(d) = tc
                .filter(|x| options.distance_shown(x))
                .and_then(|x| x.distance)
            {
                push_decimal(out, d);
            }
        }
        "OverlapBP" => {
            if let Some((bp, _)) = sv_overlap(variant, tc) {
                push_decimal(out, bp);
            }
        }
        "OverlapPC" => {
            if let Some((_, pc)) = sv_overlap(variant, tc) {
                out.extend_from_slice(format!("{pc:.2}").as_bytes());
            }
        }
        "STRAND" => {
            if let Some(x) = tc {
                push_signed(out, i64::from(x.strand));
            }
        }
        "FLAGS" => {
            if let Some(x) = tc {
                push_joined(
                    out,
                    x.flags
                        .iter()
                        .filter(|f| f.starts_with("cds_"))
                        .map(String::as_str),
                    list_sep,
                );
            }
        }
        "PICK" => {
            // An intergenic allele is flagged like a transcript allele
            // (`VariationFeatureOverlapAllele_to_output_hash`, OutputFactory.pm:1363).
            let picked = tc.map_or(variant.intergenic_pick, |x| {
                x.flags.iter().any(|f| f == "PICK")
            });
            if options.include_pick && picked {
                out.push(b'1');
            }
        }
        "VARIANT_CLASS" => out.extend_from_slice(variant.class_so_term().as_bytes()),
        "AMBIGUITY" => {
            if let Some(code) = variant.ambiguity_code() {
                let mut encoded = [0u8; 4];
                out.extend_from_slice(char::from(code).encode_utf8(&mut encoded).as_bytes());
            }
        }
        "MINIMISED" => {
            if variant.minimised {
                out.push(b'1');
            }
        }
        "SYMBOL_SOURCE" => write_identifier(out, tc.and_then(|x| x.gene_symbol_source.as_deref())),
        "HGNC_ID" => write_identifier(out, tc.and_then(|x| x.hgnc_id.as_deref())),
        "CANONICAL" => {
            if tc.is_some_and(|x| x.canonical) {
                out.extend_from_slice(b"YES");
            }
        }
        // The MANE sets the transcript belongs to, each listed under its own
        // flag state (OutputFactory.pm:1587-1601): `--mane_select` alone never
        // lists `MANE_Plus_Clinical`.
        "MANE" => {
            let select = options.mane_select && tc.is_some_and(|x| x.mane_select.is_some());
            let plus_clinical =
                options.mane_plus_clinical && tc.is_some_and(|x| x.mane_plus_clinical.is_some());
            let sets = [
                select.then_some("MANE_Select"),
                plus_clinical.then_some("MANE_Plus_Clinical"),
            ];
            push_joined(out, sets.into_iter().flatten(), list_sep);
        }
        "MANE_SELECT" => push_opt(
            out,
            tc.and_then(|x| x.mane_select.as_deref())
                .filter(|v| perl_true(v)),
        ),
        "MANE_PLUS_CLINICAL" => push_opt(
            out,
            tc.and_then(|x| x.mane_plus_clinical.as_deref())
                .filter(|v| perl_true(v)),
        ),
        "TSL" => {
            if let Some(t) = tc.and_then(|x| x.tsl) {
                push_decimal(out, u64::from(t));
            }
        }
        "GENCODE_PRIMARY" => {
            if tc.is_some_and(|x| x.gencode_primary) {
                out.push(b'1');
            }
        }
        // `s/principal/P/; s/alternative/A/` (OutputFactory.pm:1616-1622): one
        // substitution each, the rest of the value kept, so `principalM` is `PM`.
        "APPRIS" => {
            if let Some(v) = tc
                .and_then(|x| x.appris.as_deref())
                .filter(|v| perl_true(v))
            {
                let graded = v
                    .replacen("principal", "P", 1)
                    .replacen("alternative", "A", 1);
                out.extend_from_slice(graded.as_bytes());
            }
        }
        "CCDS" => write_identifier(out, tc.and_then(|x| x.ccds.as_deref())),
        "ENSP" => write_identifier(out, tc.and_then(|x| x.protein_id.as_deref())),
        "SWISSPROT" => {
            write_identifier_list(out, tc.and_then(|x| x.swissprot.as_deref()), list_sep)
        }
        "TREMBL" => write_identifier_list(out, tc.and_then(|x| x.trembl.as_deref()), list_sep),
        "UNIPARC" => write_identifier_list(out, tc.and_then(|x| x.uniparc.as_deref()), list_sep),
        "UNIPROT_ISOFORM" => {
            write_identifier_list(out, tc.and_then(|x| x.uniprot_isoform.as_deref()), list_sep)
        }
        "HGVSg" => push_opt(out, variant.hgvsg.as_deref()),
        "NEAREST" => {
            if let Some(values) = variant.nearest.as_deref() {
                push_joined(out, values.iter().map(String::as_str), list_sep);
            }
        }
        "SPDI" => push_opt(out, row_spdi(variant, tc, options).map(|(spdi, _)| spdi)),
        "GA4GH_VRS" => {
            if let Some(vrs) = row_ga4gh_vrs(variant, tc, options) {
                out.extend_from_slice(vrs.encoded().as_bytes());
            }
        }
        "CHECK_REF" => {
            if variant.check_ref_failed {
                out.extend_from_slice(b"failed");
            }
        }
        "RefSeq" => write_identifier_list(out, tc.and_then(|x| x.refseq.as_deref()), list_sep),
        "GENE_PHENO" => {
            if tc.is_some_and(|x| x.gene_phenotype) {
                out.push(b'1');
            }
        }
        "miRNA" => {
            if let Some(x) = tc {
                push_joined(out, x.mirna.iter().map(String::as_str), list_sep);
            }
        }
        // `join(':', grep {$_} (label, identifier))` (OutputFactory.pm:1500): a
        // part that is empty or `0` is false to Perl and left out, so a feature
        // with no label prints its identifier alone.
        "DOMAINS" => {
            if let Some(x) = tc {
                for (i, (label, id)) in x.domains.iter().enumerate() {
                    if i > 0 {
                        out.extend_from_slice(list_sep.as_bytes());
                    }
                    if perl_true(label) {
                        out.extend_from_slice(label.as_bytes());
                        if perl_true(id) {
                            out.push(b':');
                        }
                    }
                    if perl_true(id) {
                        out.extend_from_slice(id.as_bytes());
                    }
                }
            }
        }
        "HGVS_OFFSET" => {
            if let Some(o) = tc.and_then(|x| x.hgvs_offset) {
                push_signed(out, o);
            }
        }
        // Set on every sequence variant's transcript row, 0 included, when both
        // flags are on (`TranscriptVariationAllele_to_output_hash`,
        // OutputFactory.pm:1690-1696); an intergenic row carries none, and a
        // structural record's rows none (`TranscriptStructuralVariationAllele_to_output_hash`,
        // :2118-2125, sets no such key).
        "SHIFT_LENGTH" if options.shift_3prime && !variant.is_structural => {
            if let Some(x) = tc {
                push_decimal(out, x.shift_length);
            }
        }
        "SHIFT_LENGTH" => {}
        // Both are the reference allele's object's (OutputFactory.pm:1774-1779),
        // so a row without a transcript carries neither.
        "USED_REF" => push_opt(out, tc.and_then(|x| x.used_ref.as_deref())),
        "GIVEN_REF" => {
            if tc.is_some_and(|x| x.used_ref.is_some()) {
                out.extend_from_slice(String::from_utf8_lossy(&variant.ref_allele).as_bytes());
            }
        }
        "CLIN_SIG"
        | "SOMATIC"
        | "PHENO"
        | "PUBMED"
        | "VAR_SYNONYMS"
        | "CLINVAR_SOMATIC_CLASSIFICATION" => {
            let column = colocated_column(field, variant, options, list_sep);
            out.extend_from_slice(column.as_bytes());
        }
        "SV" => push_joined(
            out,
            variant.overlapping_svs.iter().map(String::as_str),
            list_sep,
        ),
        "FREQS" => push_joined(out, variant.freqs.iter().map(String::as_str), list_sep),
        "SIFT" => push_opt(out, tc.and_then(|x| x.sift.as_deref())),
        "PolyPhen" => push_opt(out, tc.and_then(|x| x.polyphen.as_deref())),
        "MOTIF_NAME"
        | "MOTIF_POS"
        | "HIGH_INF_POS"
        | "MOTIF_SCORE_CHANGE"
        | "TRANSCRIPTION_FACTORS" => {}
        "MAX_AF" => {
            if let Some((max, _)) = frequency::max_af(variant, &options.frequency_groups()) {
                out.extend_from_slice(max.as_bytes());
            }
        }
        "MAX_AF_POPS" => {
            if let Some((_, populations)) = frequency::max_af(variant, &options.frequency_groups())
            {
                push_joined(out, populations, list_sep);
            }
        }
        _ => match frequency::frequency_column(field) {
            Some(column) => {
                let values =
                    frequency::frequency_field(column, variant, &options.frequency_groups());
                push_joined(out, values.iter().map(String::as_str), list_sep);
            }
            None => push_opt(
                out,
                tc.and_then(|x| x.plugin_data.get(field))
                    .or_else(|| variant.plugin_data.get(field))
                    .map(String::as_str),
            ),
        },
    }
}

/// The allele's own span's positions, which the three position columns print
/// in place of their own under [`FieldOptions::positions_of_own_span`]. `None`
/// when the option is off or the row carries none, so the columns print their own.
fn own_span_positions<'a>(
    tc: &'a TranscriptConsequence,
    options: &FieldOptions,
) -> Option<&'a OwnSpanPositions> {
    if options.positions_of_own_span {
        tc.own_span_positions.as_deref()
    } else {
        None
    }
}

/// Appends a position column with the total `--total_length` appends:
/// `428/1200`, or `-/1200` where the position is undefined but the total
/// applies, since Ensembl VEP appends to its `-` placeholder (`OutputFactory.pm`
/// `TranscriptVariationAllele_to_output_hash`); without a total the position
/// alone, nothing when undefined.
fn write_position_column(out: &mut Vec<u8>, position: Option<&str>, total: Option<u32>) {
    match total {
        Some(total) => {
            out.extend_from_slice(position.unwrap_or("-").as_bytes());
            out.push(b'/');
            push_decimal(out, u64::from(total));
        }
        None => push_opt(out, position),
    }
}

/// `write_position_column` as a `String`, empty when undefined without a total.
pub fn position_column(position: Option<&str>, total: Option<u32>) -> String {
    let mut out = Vec::new();
    write_position_column(&mut out, position, total);
    written_text(out)
}

/// The fields `add_colocated_variant_info` (OutputFactory.pm:1041-1167 at
/// release/116.2) and `add_colocated_frequency_data` (:1185-1278) set on a row:
/// the co-located variant columns and the fields of the frequency groups
/// (`%FREQUENCY_KEYS`, :107-113, and `--max_af`). The JSON factory stubs both
/// out (OutputFactory/JSON.pm:493-515) and writes the data under
/// `colocated_variants` instead; a plugin's field is never one of them.
pub fn colocated_field(field: &str) -> bool {
    const FREQUENCY_GROUPS: [&str; 6] = [
        "af",
        "af_1kg",
        "af_gnomad",
        "af_gnomade",
        "af_gnomadg",
        "max_af",
    ];
    matches!(
        field,
        "Existing_variation"
            | "CLIN_SIG"
            | "SOMATIC"
            | "PHENO"
            | "PUBMED"
            | "VAR_SYNONYMS"
            | "CLINVAR_SOMATIC_CLASSIFICATION"
            | "FREQS"
    ) || FLAG_FIELD_GROUPS
        .iter()
        .any(|group| FREQUENCY_GROUPS.contains(&group.flag) && group.fields.contains(&field))
}

/// One co-located variant column of a row, as `add_colocated_variant_info`
/// (OutputFactory.pm:1041-1167 at release/116.2) fills it from the variants
/// that print on the row: every matched variant in output order, a variant
/// matched by allele only on the rows of an allele it matched (:1081-1083).
///
/// `SOMATIC` and `PHENO` list one `0` or `1` per variant and are dropped when
/// every value is `0` (:1129-1138); `PUBMED` lists every id of every variant,
/// duplicates kept (:1126); `VAR_SYNONYMS` joins the variants' synonym strings
/// with `--` (:1090, 1141); `CLINVAR_SOMATIC_CLASSIFICATION` is the last
/// variant's `clinical_impact` string (:1093-1095). `CLIN_SIG` under
/// [`ClinSigAllele::PerAllele`] is the row's allele's own entries, `;`-joined
/// across variants and `,`-joined within one, when any variant carries entries
/// for it (:1097-1122, 1143-1144); otherwise, and under
/// [`ClinSigAllele::WholeVariant`], the terms of the variants read before the
/// first one with allele-specific entries, duplicates kept (:1125). Every
/// entry allele of a variant whose `clin_sig_ref_allele` differs from the
/// reference `REF_ALLELE` prints ([`feature_reference`]) is reverse-complemented
/// first (:1049, 1100-1112).
fn colocated_column(
    field: &str,
    variant: &InputVariant,
    options: &FieldOptions,
    list_sep: &str,
) -> String {
    if variant.colocated_variants.is_empty() {
        return String::new();
    }
    let allele = variant.display_allele();
    let ordered = colocated_on_row(variant, &allele);
    match field {
        "CLIN_SIG" => clin_sig_column(variant, &allele, &ordered, options, list_sep),
        "SOMATIC" => zero_one_list(ordered.iter().map(|cv| cv.somatic), list_sep),
        "PHENO" => zero_one_list(ordered.iter().map(|cv| cv.phenotype), list_sep),
        "PUBMED" if options.pubmed => ordered
            .iter()
            .flat_map(|cv| cv.pubmed.iter().map(String::as_str))
            .collect::<Vec<_>>()
            .join(list_sep),
        "VAR_SYNONYMS" if options.var_synonyms => ordered
            .iter()
            .filter_map(|cv| cv.var_synonyms.as_deref().filter(|s| perl_true(s)))
            .collect::<Vec<_>>()
            .join("--"),
        "CLINVAR_SOMATIC_CLASSIFICATION" if options.clinvar_somatic_classification => ordered
            .iter()
            .filter_map(|cv| cv.clinvar.as_deref()?.clinical_impact.as_deref())
            .next_back()
            .unwrap_or_default()
            .to_string(),
        _ => String::new(),
    }
}

/// The co-located variants that print on the row of `allele`, in output order:
/// every matched variant, a variant matched by allele only on the rows of an
/// allele it matched (OutputFactory.pm:1069-1083).
fn colocated_on_row<'v>(variant: &'v InputVariant, allele: &str) -> Vec<&'v ColocatedVariant> {
    let mut ordered = colocated_in_output_order(&variant.colocated_variants);
    ordered.retain(|cv| cv.applies_to_allele(allele));
    ordered
}

/// `0` or `1` per value joined by `list_sep`, or nothing when every value is
/// false (OutputFactory.pm:1129-1138).
fn zero_one_list(values: impl Iterator<Item = bool> + Clone, list_sep: &str) -> String {
    if !values.clone().any(|value| value) {
        return String::new();
    }
    values
        .map(|value| if value { "1" } else { "0" })
        .collect::<Vec<_>>()
        .join(list_sep)
}

/// The reference allele Ensembl VEP reads as the record feature's
/// (`$vf->ref_allele_string`, VariationFeature.pm:1533-1539), which `REF_ALLELE`
/// prints (OutputFactory.pm:1357) and the ClinVar entry alleles are read against
/// (:1049-1051). Every row of a `--minimal` split record reads the reference
/// cached on the record's feature before the rejoin restored the record's
/// string: the first allele's trimmed reference
/// (`vep_core::minimal::rejoined_reference`), or the record's own when the first
/// feature wrote no row (`split_first_feature_rowless`). Every other row reads
/// the variant's own reference.
fn feature_reference<'v>(variant: &'v InputVariant, options: &FieldOptions) -> Cow<'v, str> {
    match variant.record_allele_string_multi.as_deref() {
        Some(record) if variant.split_first_feature_rowless => {
            Cow::Borrowed(record.split('/').next().unwrap_or_default())
        }
        Some(record) if variant.original_start.is_some() => Cow::Owned(
            vep_core::minimal::rejoined_reference(record, options.process_ref_homs),
        ),
        _ => String::from_utf8_lossy(&variant.ref_allele),
    }
}

/// The `CLIN_SIG` of a row over the variants that print on it, in order
/// (`colocated_column`).
fn clin_sig_column(
    variant: &InputVariant,
    allele: &str,
    ordered: &[&ColocatedVariant],
    options: &FieldOptions,
    list_sep: &str,
) -> String {
    let ref_allele = feature_reference(variant, options);
    let per_allele = options.clin_sig_allele == ClinSigAllele::PerAllele;
    let mut terms: Vec<&str> = Vec::new();
    let mut by_allele: Vec<String> = Vec::new();
    let mut allele_entries_seen = false;
    for cv in ordered {
        let allele_specific = cv
            .clinvar
            .as_deref()
            .filter(|clinvar| per_allele && clinvar.clin_sig_allele.is_some());
        if let Some(clinvar) = allele_specific {
            let flip = clinvar
                .clin_sig_ref_allele
                .as_deref()
                .is_some_and(|ref_entry| perl_true(ref_entry) && ref_entry != ref_allele);
            // Perl's `%cs_hash` joins each entry allele's terms with `,` in entry
            // order and is read at the row's allele alone, so only that allele's
            // terms are gathered.
            let mut row_allele_terms = String::new();
            let mut row_allele_entry_found = false;
            for (entry_allele, term) in clinvar.clin_sig_allele_entries() {
                let key = if flip {
                    reverse_comp(entry_allele)
                } else {
                    entry_allele.to_string()
                };
                if key != allele {
                    continue;
                }
                if !row_allele_terms.is_empty() {
                    row_allele_terms.push(',');
                }
                row_allele_terms.push_str(term);
                row_allele_entry_found = true;
            }
            if row_allele_entry_found && !by_allele.contains(&row_allele_terms) {
                by_allele.push(row_allele_terms);
            }
            allele_entries_seen = true;
        }
        if !allele_entries_seen {
            terms.extend(cv.clin_sig.iter().map(String::as_str));
        }
    }
    if by_allele.is_empty() {
        terms.join(list_sep)
    } else {
        by_allele.join(";")
    }
}

/// `OverlapBP` and `OverlapPC` of a structural variant row: the base pairs of
/// the variant span inside the feature and the percentage of the feature they
/// cover. VEP computes them for every structural variant except breakpoints
/// (ensembl-vep `OutputFactory.pm`) and prints them only when the overlap is positive.
fn sv_overlap(variant: &InputVariant, tc: Option<&TranscriptConsequence>) -> Option<(u64, f64)> {
    if !variant.is_structural || variant.variant_class == VariantClass::Translocation {
        return None;
    }
    // VEP reads the allele's own start and end (`OutputFactory.pm`,
    // `StructuralVariationOverlapAllele_to_output_hash`), so an insertion pair
    // (`start == end + 1`) has an overlap length of 0 and carries neither field.
    tc?.feature_overlap(variant.start, variant.end)
}

/// Whether a field whose value is empty is still written, as `NEAREST=` in the
/// default format and an empty tab column: on a chromosome with no
/// protein-coding transcription start site Ensembl VEP sets `NEAREST` to an
/// empty list, which is a value (`OutputFactory.pm:1017`), and under
/// `--clinvar_somatic_classification` it sets `CLINVAR_SOMATIC_CLASSIFICATION`
/// to every defined `clinical_impact` of the row's variants in turn (:1093),
/// the empty cell of a region file included; the default and tab writers print
/// every key set (`VEP_output.pm:117`, `Tab.pm:128`). Every other empty value
/// stands for an absent key. A `--summary` or `--most_severe` row is the
/// record's hash alone (`summary_only`, OutputFactory.pm:426-458) and carries
/// the variant-level fields and nothing else, as [`write_field_value`] writes
/// it, so only a variant-level field's empty value is written on it; `tc` is
/// the row's transcript consequence, as [`write_field_value`] takes it.
pub fn empty_value_is_written(
    field: &str,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> bool {
    let summary_row = tc.is_none() && variant.summary_terms.is_some();
    if summary_row && !is_vf_level_field(field) {
        return false;
    }
    match field {
        "NEAREST" => variant.nearest.as_ref().is_some_and(Vec::is_empty),
        "CLINVAR_SOMATIC_CLASSIFICATION" => {
            options.clinvar_somatic_classification
                && colocated_on_row(variant, &variant.display_allele())
                    .iter()
                    .filter_map(|cv| cv.clinvar.as_deref())
                    .any(|clinvar| clinvar.clinical_impact.is_some())
        }
        _ => false,
    }
}

/// Whether a `-` value is a value of the field rather than an absent key, for
/// the writers that drop a `-`: the JSON writer keeps it on `Allele`,
/// `USED_REF` and `GIVEN_REF`, where it is an empty allele
/// (`OutputFactory/JSON.pm:327`), and the Parquet writer follows the JSON one.
pub fn dash_is_a_value(field: &str) -> bool {
    matches!(field, "Allele" | "USED_REF" | "GIVEN_REF")
}

/// The keys of a row's `Extra` column in VEP's order: the active flag fields
/// first, then every other key the row carries (the structural variant overlap
/// fields and plugin keys), alphabetically.
pub fn extra_keys<'a>(
    variant: &'a InputVariant,
    tc: Option<&'a TranscriptConsequence>,
    options: &FieldOptions,
) -> Vec<&'a str> {
    extra_keys_with(&flag_fields(options), variant, tc, options)
}

/// [`extra_keys`] given the configuration's [`flag_fields`] already computed.
pub fn extra_keys_with<'a>(
    flagged: &[&'static str],
    variant: &'a InputVariant,
    tc: Option<&'a TranscriptConsequence>,
    options: &FieldOptions,
) -> Vec<&'a str> {
    let mut keys: Vec<&str> = flagged.to_vec();
    let mut extra: BTreeSet<&str> = BTreeSet::new();
    if !options.overlaps && sv_overlap(variant, tc).is_some() {
        extra.insert("OverlapBP");
        extra.insert("OverlapPC");
    }
    for k in variant.plugin_data.keys() {
        if !flagged.contains(&k.as_str()) {
            extra.insert(k.as_str());
        }
    }
    if let Some(tc) = tc {
        for k in tc.plugin_data.keys() {
            if !flagged.contains(&k.as_str()) {
                extra.insert(k.as_str());
            }
        }
    }
    keys.extend(extra);
    keys
}

/// A field a loaded plugin adds to the output: its name and the description the
/// headers print for it.
pub type PluginHeader = (String, String);

/// The names of a plugin header list, in its order.
pub fn plugin_field_names(plugin_headers: &[PluginHeader]) -> Vec<String> {
    plugin_headers
        .iter()
        .map(|(name, _)| name.clone())
        .collect()
}

/// The description a header prints for a column: the field table's, else the
/// plugin's own, else `?` for a name nothing produces (`OutputFactory/Tab.pm`
/// `description_headers`).
fn column_description<'a>(column: &str, plugin_headers: &'a [PluginHeader]) -> &'a str {
    field_description(column)
        .or_else(|| {
            plugin_headers
                .iter()
                .find(|(name, _)| name == column)
                .map(|(_, description)| description.as_str())
        })
        .unwrap_or("?")
}

/// The header lines of the default format, in VEP's order, ending with the
/// column line: the flag fields' keys, then one `## <KEY> : <description>` line
/// per plugin field (`OutputFactory/VEP_output.pm` `description_headers`).
/// `--no_headers` callers write none of them.
pub fn default_format_header(
    options: &FieldOptions,
    plugin_headers: &[PluginHeader],
    cache_dir: &str,
    cache_info: Option<&vep_core::cache_info::CacheInfo>,
    command_line: &str,
) -> Vec<String> {
    let mut lines = common_header_prefix(cache_dir, cache_info);
    lines.push("## Column descriptions:".to_string());
    for col in DEFAULT_OUTPUT_COLS {
        lines.push(format!(
            "## {} : {}",
            col,
            field_description(col).unwrap_or("?")
        ));
    }
    lines.push("## Extra column keys:".to_string());
    for key in flag_fields(options) {
        lines.push(format!(
            "## {} : {}",
            key,
            field_description(key).unwrap_or("?")
        ));
    }
    for (name, description) in plugin_headers {
        lines.push(format!("## {name} : {description}"));
    }
    lines.push(format!("## VEP command-line: {}", command_line));
    lines.push(format!("#{}\tExtra", DEFAULT_OUTPUT_COLS.join("\t")));
    lines
}

/// The header lines of the tab format: every column, flag and plugin fields
/// included, gets a description line, and the column line names them all.
/// `selected_fields` is the `--fields` list, which replaces the columns.
pub fn tab_header(
    options: &FieldOptions,
    plugin_headers: &[PluginHeader],
    selected_fields: Option<&[String]>,
    cache_dir: &str,
    cache_info: Option<&vep_core::cache_info::CacheInfo>,
    command_line: &str,
) -> Vec<String> {
    let cols = tab_columns(
        options,
        &plugin_field_names(plugin_headers),
        selected_fields,
    );
    let mut lines = common_header_prefix(cache_dir, cache_info);
    lines.push("## Column descriptions:".to_string());
    for col in &cols {
        lines.push(format!(
            "## {} : {}",
            col,
            column_description(col, plugin_headers)
        ));
    }
    lines.push(format!("## VEP command-line: {}", command_line));
    lines.push(format!("#{}", cols.join("\t")));
    lines
}

/// The tab format's columns: the thirteen defaults, the flag fields, then plugin
/// fields; or, under `--fields`, that list in its own order and nothing else
/// (`OutputFactory/Tab.pm` `fields`, where the flag's value pre-empts the list).
pub fn tab_columns(
    options: &FieldOptions,
    plugin_fields: &[String],
    selected_fields: Option<&[String]>,
) -> Vec<String> {
    if let Some(selected_fields) = selected_fields {
        return selected_fields.to_vec();
    }
    let mut cols: Vec<String> = DEFAULT_OUTPUT_COLS.iter().map(|s| s.to_string()).collect();
    cols.extend(flag_fields(options).into_iter().map(String::from));
    for p in plugin_fields {
        if !cols.iter().any(|c| c == p) {
            cols.push(p.clone());
        }
    }
    cols
}

/// The `##VEP=` meta line of a VCF output: version, API version, time, the
/// cache directory and every cache source version as `key="value"`, sorted.
pub fn vcf_vep_meta_line(
    cache_dir: &str,
    cache_info: Option<&vep_core::cache_info::CacheInfo>,
) -> String {
    let mut line = format!(
        "##VEP=\"v{}.{}\" API=\"v{}\" time=\"{}\"",
        vep_core::VEP_VERSION,
        vep_core::VEP_SUB_VERSION,
        vep_core::VEP_VERSION,
        timestamp()
    );
    if !cache_dir.is_empty() {
        line.push_str(&format!(" cache=\"{cache_dir}\""));
    }
    if let Some(info) = cache_info {
        for (key, value) in &info.source_versions {
            line.push_str(&format!(" {key}=\"{value}\""));
        }
    }
    line
}

/// The header lines the default and tab formats share, in VEP's order, with
/// the vep-rs identity line second so a reader of the first line still sees
/// the Ensembl VEP banner.
fn common_header_prefix(
    cache_dir: &str,
    cache_info: Option<&vep_core::cache_info::CacheInfo>,
) -> Vec<String> {
    let mut lines = vec![
        format!(
            "## ENSEMBL VARIANT EFFECT PREDICTOR v{}.{}",
            vep_core::VEP_VERSION,
            vep_core::VEP_SUB_VERSION
        ),
        format!("## vep-rs {VEP_RS_VERSION}"),
        format!("## Output produced at {}", timestamp()),
        format!("## Using cache in {}", cache_dir),
        format!(
            "## Using API version {}, DB version ?",
            vep_core::VEP_VERSION
        ),
    ];
    if let Some(info) = cache_info {
        for (key, value) in &info.source_versions {
            lines.push(format!("## {} version {}", key, value));
        }
    }
    lines
}

/// The current UTC time in the form VEP's `Output produced at` line uses.
fn timestamp() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Civil date from days since the epoch (proleptic Gregorian), UTC.
    let days = (secs / 86_400) as i64;
    let rem = secs % 86_400;
    let (h, m, s) = (rem / 3600, (rem % 3600) / 60, rem % 60);
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };
    format!("{y:04}-{mo:02}-{d:02} {h:02}:{m:02}:{s:02}")
}

#[cfg(test)]
pub(crate) mod fixtures {
    //! Rows carrying a value in every column the formats can list, for the
    //! byte-identity tests of this module and of the tab writer.

    use super::*;
    use std::sync::Arc;

    use vep_core::consequence::Impact;
    use vep_core::spdi::{SiblingSpdi, VrsAllele};
    use vep_core::variant::{MatchedAllele, SampleZygosity, Zygosity, ZygosityList};
    use vep_core::variation::ClinVarDetail;

    /// Every column name the formats can ask for: the fixed thirteen, every
    /// flag field of every group, the overlap pair, a plugin key held on the
    /// consequence and one held on the variant, and a name nothing produces.
    pub(crate) fn every_column() -> Vec<&'static str> {
        let mut columns: Vec<&'static str> = DEFAULT_OUTPUT_COLS.to_vec();
        for group in FLAG_FIELD_GROUPS {
            for f in group.fields {
                if !columns.contains(f) {
                    columns.push(f);
                }
            }
        }
        for c in ["ZPlugin", "VPlugin", "Unknown"] {
            if !columns.contains(&c) {
                columns.push(c);
            }
        }
        columns
    }

    /// A consequence with a value in every transcript-side field: both id
    /// versions, positions with totals and an own span, an upstream term beside
    /// the coding ones so `DISTANCE` prints, `=` twice in HGVSp, a `PICK` flag
    /// between two `cds_` ones, accession lists with a trailing comma, a domain
    /// without a label, one without an identifier and one with neither, and a
    /// plugin key.
    pub(crate) fn populated_consequence() -> TranscriptConsequence {
        let mut tc = TranscriptConsequence {
            transcript_id: Arc::from("ENST00000300"),
            gene_id: Arc::from("ENSG00000100"),
            gene_symbol: Some(Arc::from("SOD1")),
            gene_symbol_source: Some(Arc::from("HGNC")),
            hgnc_id: Some(Arc::from("HGNC:11179")),
            consequences: [
                Consequence::MissenseVariant,
                Consequence::SpliceRegionVariant,
                Consequence::UpstreamGeneVariant,
            ]
            .into_iter()
            .collect(),
            impact: Impact::MODERATE,
            biotype: Some(Arc::from("protein_coding")),
            canonical: true,
            cdna_position: Some("428".to_string()),
            cds_position: Some("301".to_string()),
            protein_position: Some("101".to_string()),
            cdna_total: Some(1_200),
            cds_total: Some(900),
            protein_total: Some(300),
            own_span_positions: Some(Box::new(OwnSpanPositions {
                cdna: Some("430-431".to_string()),
                cds: None,
                protein: Some("102".to_string()),
            })),
            amino_acids: Some("A/V".to_string()),
            codons: Some("gCc/gTc".to_string()),
            protein_id: Some(Arc::from("ENSP00000300")),
            distance: Some(2_500),
            strand: -1,
            exon: Some("3/5".to_string()),
            intron: Some("2/4".to_string()),
            hgvsc: Some("ENST00000300.1:c.301C>T".to_string()),
            hgvsp: Some("ENSP00000300.1:p.Ala101=;p.=".to_string()),
            hgvs_offset: Some(-3),
            shift_length: 2,
            used_ref: Some(Arc::from("C")),
            sift: Some("deleterious(0.01)".to_string()),
            polyphen: Some("probably_damaging(0.98)".to_string()),
            domains: vec![
                ("Pfam".to_string(), "PF00080".to_string()),
                ("PANTHER".to_string(), "PTHR10003".to_string()),
                (String::new(), "PS50047".to_string()),
                ("Low_complexity_(Seg)".to_string(), "0".to_string()),
                (String::new(), "0".to_string()),
            ],
            flags: Arc::from([
                "cds_start_NF".to_string(),
                "PICK".to_string(),
                "cds_end_NF".to_string(),
            ]),
            tsl: Some(1),
            mane_select: Some("NM_000454.5".to_string()),
            mane_plus_clinical: Some("NM_000454.4".to_string()),
            appris: Some("principal1".to_string()),
            ccds: Some("CCDS33536.1".to_string()),
            gencode_primary: true,
            swissprot: Some("P00441,Q9".to_string()),
            trembl: Some("A0A0,".to_string()),
            uniparc: Some("UPI000012A3C5".to_string()),
            uniprot_isoform: Some("P00441-1".to_string()),
            refseq: Some("NM_000454.5,NM_000454.4".to_string()),
            transcript_version: Some(4),
            gene_version: Some(12),
            gene_phenotype: true,
            mirna: vec!["miRNA_stem".to_string(), "miRNA_loop".to_string()],
            feature_start: 1_000,
            feature_end: 2_000,
            id_rank: 3,
            ..TranscriptConsequence::default()
        };
        tc.plugin_data
            .insert("ZPlugin".to_string(), "a=b;c".to_string());
        tc
    }

    /// A known variant matched on `matched` (input allele, record allele) with
    /// a frequency for every population column of every group.
    fn known(
        id: &str,
        matched: &[(&str, &str)],
        clinvar: Option<ClinVarDetail>,
    ) -> ColocatedVariant {
        let mut frequency_strings = std::collections::HashMap::new();
        for (i, column) in frequency::FREQUENCY_KEYS
            .iter()
            .flat_map(|(_, columns)| columns.iter())
            .enumerate()
        {
            frequency_strings.insert(
                column.to_string(),
                format!("G:{}", 0.0031 * (i as f64 + 1.0)),
            );
        }
        frequency_strings.insert("gnomADe_NFE".to_string(), "G:1.886e-05".to_string());
        ColocatedVariant {
            id: id.to_string(),
            allele_string: Some("A/G/T".to_string()),
            strand: 1,
            clin_sig: vec!["pathogenic".to_string(), "likely_pathogenic".to_string()],
            somatic: id.starts_with("COS"),
            phenotype: true,
            pubmed: vec!["123".to_string(), "456".to_string()],
            matched_alleles: matched
                .iter()
                .enumerate()
                .map(|(i, (input, record))| MatchedAllele {
                    input_allele: input.to_string(),
                    input_index: i,
                    record_allele: record.to_string(),
                    record_index: i,
                })
                .collect(),
            frequency_strings,
            clinvar: clinvar.map(Box::new),
            var_synonyms: Some(format!("ClinVar::RCV{id}")),
            ..Default::default()
        }
    }

    /// A `chr21` `A/G` allele of a three-allele record with every variant-side
    /// value set: ids, two co-located variants (one with allele-specific ClinVar
    /// entries and a somatic classification), every frequency, HGVSg, SPDI with
    /// its VRS allele and a sibling's, nearest sites, an individual's genotype
    /// and the sample list, overlapping structural variants, the frequency
    /// filter's list, a failed reference check, a `--minimal` rejoin, and a
    /// plugin key.
    pub(crate) fn populated_variant() -> InputVariant {
        let mut v = InputVariant::new("21".into(), 10_500, 10_500, b"A".to_vec(), b"G".to_vec());
        v.original_chr = "chr21".to_string();
        v.id = Some("rs1".to_string());
        v.allele_index = 1;
        v.allele_num = Some(2);
        v.minimised = true;
        v.original_allele_string = Some("A/G/T".to_string());
        v.original_start = Some(10_499);
        v.original_end = Some(10_500);
        v.record_allele_string_multi = Some(Arc::from("A/G/T"));
        v.uploaded_allele_string = Some(Arc::from("CA/CG/CT"));
        v.existing_variation = vec!["rs1".to_string(), "COSV2".to_string()];
        v.hgvsg = Some("chr21:g.10500A>G".to_string());
        v.nearest = Some(vec!["ENST00000010".to_string(), "ENST00000020".to_string()]);
        v.spdi = Some("NC_000021.9:10499:A:G".to_string());
        v.ga4gh_vrs = Some(VrsAllele {
            sequence_id: "NC_000021.9".to_string(),
            start: 10_499,
            end: 10_500,
            sequence: "G".to_string(),
        });
        v.sibling_spdi = vec![SiblingSpdi {
            allele: "T".to_string(),
            spdi: "NC_000021.9:10499:A:T".to_string(),
            ga4gh_vrs: None,
        }];
        v.overlapping_svs = vec!["esv1".to_string(), "nsv2".to_string()];
        v.freqs = vec![
            "G:AF:0.0031".to_string(),
            "G:gnomADe_NFE:1.886e-05".to_string(),
        ];
        v.check_ref_failed = true;
        v.intergenic_shift = Some(1);
        v.intergenic_pick = true;
        v.individual = Some(Arc::from("S1"));
        v.genotype = vec!["A".to_string(), "G".to_string()];
        v.sample_zygosities = Some(Arc::new(ZygosityList::new(vec![
            SampleZygosity {
                sample: Arc::from("S1"),
                zygosity: Zygosity::Het,
            },
            SampleZygosity {
                sample: Arc::from("S2"),
                zygosity: Zygosity::HomRef,
            },
        ])));
        v.most_severe_consequence = Some(Consequence::MissenseVariant);
        v.plugin_data
            .insert("VPlugin".to_string(), "x;y".to_string());
        v.colocated_variants = vec![
            known(
                "rs1",
                &[("G", "G")],
                Some(ClinVarDetail {
                    clin_sig_allele: Some("G:pathogenic;T:benign".to_string()),
                    clin_sig_ref_allele: Some("A".to_string()),
                    clinical_impact: Some("Tier_I_-_Strong(diagnostic)".to_string()),
                }),
            ),
            known("COSV2", &[("G", "G"), ("T", "T")], None),
        ];
        v
    }

    /// A structural deletion over the populated consequence's feature, so its
    /// rows carry the overlap pair and the record's own class term.
    pub(crate) fn structural_variant() -> InputVariant {
        let mut v = InputVariant::new("21".into(), 1_500, 1_600, b"N".to_vec(), b"<DEL>".to_vec());
        v.id = Some("sv1".to_string());
        v.is_structural = true;
        v.variant_class = VariantClass::StructuralDeletion;
        v.record_so_term = Some("deletion".to_string());
        v.sv_end = Some(1_600);
        v.most_severe_consequence = Some(Consequence::TranscriptAblation);
        v.existing_variation = vec!["esv1".to_string()];
        v.overlapping_svs = vec!["esv1".to_string()];
        v
    }

    /// Every field group on, every version and total shown, the display
    /// vocabulary off; every field listed so a new option fails here first.
    pub(crate) fn all_options() -> FieldOptions {
        FieldOptions {
            individual: true,
            individual_zyg: true,
            process_ref_homs: true,
            allele_number: true,
            show_ref_allele: true,
            uploaded_allele: true,
            include_pick: true,
            variant_class: true,
            minimal: true,
            clinvar_somatic_classification: true,
            symbol: true,
            biotype: true,
            canonical: true,
            mane_select: true,
            mane: true,
            mane_plus_clinical: true,
            tsl: true,
            gencode_primary: true,
            appris: true,
            ccds: true,
            protein: true,
            uniprot: true,
            xref_refseq: true,
            gene_phenotype: true,
            sift: true,
            polyphen: true,
            numbers: true,
            domains: true,
            mirna: true,
            hgvs: true,
            hgvsg: true,
            spdi: true,
            ga4gh_vrs: true,
            af: true,
            af_1kg: true,
            af_gnomad: false,
            af_gnomade: true,
            af_gnomadg: true,
            max_af: true,
            check_frequency: true,
            check_existing: true,
            pubmed: true,
            check_svs: true,
            dont_skip: true,
            overlaps: true,
            nearest: true,
            shift_length: true,
            use_transcript_ref: true,
            shift_3prime: true,
            shift_genomic: true,
            var_synonyms: true,
            regulatory: true,
            transcript_version: true,
            gene_version: true,
            ambiguity: true,
            no_escape: false,
            terms: Vocabulary::So,
            total_length: true,
            positions_of_own_span: false,
            clin_sig_allele: ClinSigAllele::PerAllele,
        }
    }
}
#[cfg(test)]
mod tests {
    use super::fixtures::*;
    use super::*;
    use std::sync::Arc;
    use vep_core::variant::{MatchedAllele, SampleZygosity, Zygosity, ZygosityList};
    use vep_core::variation::ClinVarDetail;

    // The field table as owned strings, the oracle the byte writer is held to
    // over every column: `string_field_value` and the helpers it calls build
    // each value as a `String` the way the formats read them; the live
    // `write_field_value` appends the same bytes and `field_value` collects them.

    /// The row's `Allele`: the allele as the record gives it; under
    /// `--shift_3prime 1 --shift_genomic 1`, rotated along the row's shift
    /// (OutputFactory.pm:1349-1352); under `--shift_3prime 1` alone, the sequence
    /// the row's codon computation rotated ([`TranscriptConsequence::shifted_allele`],
    /// :1347). A reference row's is the reference allele's sequence (:1347):
    /// [`TranscriptConsequence::used_ref`] when the row has one.
    fn string_row_allele(
        variant: &InputVariant,
        tc: Option<&TranscriptConsequence>,
        options: &FieldOptions,
    ) -> String {
        if variant.reference_allele_row {
            if let Some(used_ref) = tc.and_then(|x| x.used_ref.as_deref()) {
                return used_ref.to_string();
            }
        }
        if let Some((shift, strand)) = genomic_shift(variant, tc, options) {
            variant.display_allele_shifted(shift, strand)
        } else if let Some(allele) = tc.and_then(|x| x.shifted_allele.as_deref()) {
            allele.to_string()
        } else {
            variant.display_allele()
        }
    }

    /// `id`, with `.<version>` appended when `flag` is on and [`version_suffix`]
    /// yields one.
    fn string_versioned_id(id: &str, version: Option<u32>, flag: bool) -> String {
        if !flag {
            return id.to_string();
        }
        match version_suffix(id, version) {
            Some(version) => format!("{id}.{version}"),
            None => id.to_string(),
        }
    }

    /// A cache identifier field as Ensembl VEP reads it: the cache's `-` stands for
    /// no value and sets nothing (OutputFactory.pm:1525-1572 at release/116.2), and
    /// a comma-separated list of accessions is split as Perl's `split` splits it
    /// ([`perl_split`]: a trailing comma adds no element), which the formats join
    /// with their own separator.
    fn string_identifier_list(value: Option<&str>, list_sep: &str) -> String {
        match value {
            None | Some("-") => String::new(),
            Some(v) => perl_split(v, &[',']).collect::<Vec<_>>().join(list_sep),
        }
    }

    /// A cache identifier field that is one value: the `-` sentinel prints nothing.
    fn string_identifier(value: Option<&str>) -> String {
        match value {
            None | Some("-") => String::new(),
            Some(v) => v.to_string(),
        }
    }

    /// The value of one named field for one row, or an empty string when the field
    /// is absent on this row. `list_sep` joins list-valued fields (`,` in the default
    /// and tab formats, `&` inside a VCF CSQ entry). Field names not in the table are
    /// looked up in the row's plugin data.
    ///
    /// A row carrying the record's [`InputVariant::summary_terms`] is Ensembl VEP's
    /// `summary_only` row (OutputFactory.pm:426-458): the record's hash with the
    /// summary's `Consequence`, so it has the variant's identity, location and
    /// variant-level fields and nothing else, no Allele, IMPACT, Existing_variation
    /// or plugin key included.
    pub fn string_field_value(
        field: &str,
        variant: &InputVariant,
        tc: Option<&TranscriptConsequence>,
        options: &FieldOptions,
        list_sep: &str,
    ) -> String {
        if let (None, Some(terms)) = (tc, &variant.summary_terms) {
            return match field {
                "Consequence" => {
                    let rendered = rendered_terms(terms, options.terms, true);
                    if rendered.is_empty() {
                        "-".to_string()
                    } else {
                        rendered.join(list_sep)
                    }
                }
                // The record hash places the row at the input span; the
                // `--shift_genomic 1` move is written by the feature and intergenic
                // row builders alone (OutputFactory.pm:1702-1706, 1992-1998).
                "Location" => variant.location(),
                _ if is_vf_level_field(field) => {
                    string_row_field_value(field, variant, None, options, list_sep)
                }
                _ => String::new(),
            };
        }
        string_row_field_value(field, variant, tc, options, list_sep)
    }

    /// [`string_field_value`] of a row against a feature or an allele's intergenic row.
    fn string_row_field_value(
        field: &str,
        variant: &InputVariant,
        tc: Option<&TranscriptConsequence>,
        options: &FieldOptions,
        list_sep: &str,
    ) -> String {
        match field {
            "Uploaded_variation" => variant.uploaded_variation(),
            "Location" => match genomic_shift(variant, tc, options) {
                Some((shift, strand)) => variant.location_shifted(shift, strand),
                None => variant.location(),
            },
            "Allele" => string_row_allele(variant, tc, options),
            "Consequence" => match tc {
                Some(x) => consequence_string(x, options),
                None => {
                    let consequence = variant
                        .most_severe_consequence
                        .unwrap_or(Consequence::IntergenicVariant);
                    options.terms.term(consequence).to_string()
                }
            },
            "IMPACT" => tc
                .map(|x| x.impact.to_string())
                .or_else(|| {
                    variant
                        .most_severe_consequence
                        .map(|c| c.impact().to_string())
                })
                .unwrap_or_else(|| Consequence::IntergenicVariant.impact().to_string()),
            "SYMBOL" => string_identifier(tc.and_then(|x| x.gene_symbol.as_deref())),
            "Gene" => tc
                .map(|x| string_versioned_id(&x.gene_id, x.gene_version, options.gene_version))
                .unwrap_or_default(),
            "Feature_type" => tc.map(|x| x.feature_type.to_string()).unwrap_or_default(),
            "Feature" => tc
                .map(|x| {
                    string_versioned_id(
                        &x.transcript_id,
                        x.transcript_version,
                        options.transcript_version,
                    )
                })
                .unwrap_or_default(),
            "BIOTYPE" if !options.biotype => String::new(),
            "BIOTYPE" => tc
                .and_then(|x| x.biotype.as_deref())
                .unwrap_or("")
                .to_string(),
            "EXON" => tc.and_then(|x| x.exon.clone()).unwrap_or_default(),
            "INTRON" => tc.and_then(|x| x.intron.clone()).unwrap_or_default(),
            "HGVSc" if !options.hgvs => String::new(),
            "HGVSp" if !options.hgvs => String::new(),
            "HGVSc" => tc.and_then(|x| x.hgvsc.clone()).unwrap_or_default(),
            "HGVSp" => {
                let v = tc.and_then(|x| x.hgvsp.clone()).unwrap_or_default();
                if options.no_escape {
                    v
                } else {
                    v.replace('=', "%3D")
                }
            }
            "cDNA_position" => string_position_column(
                tc.and_then(|x| match own_span_positions(x, options) {
                    Some(own) => own.cdna.as_deref(),
                    None => x.cdna_position.as_deref(),
                }),
                tc.and_then(|x| x.cdna_total)
                    .filter(|_| options.total_length),
            ),
            "CDS_position" => string_position_column(
                tc.and_then(|x| match own_span_positions(x, options) {
                    Some(own) => own.cds.as_deref(),
                    None => x.cds_position.as_deref(),
                }),
                tc.and_then(|x| x.cds_total)
                    .filter(|_| options.total_length),
            ),
            "Protein_position" => string_position_column(
                tc.and_then(|x| match own_span_positions(x, options) {
                    Some(own) => own.protein.as_deref(),
                    None => x.protein_position.as_deref(),
                }),
                tc.and_then(|x| x.protein_total)
                    .filter(|_| options.total_length),
            ),
            "Amino_acids" => tc.and_then(|x| x.amino_acids.clone()).unwrap_or_default(),
            "Codons" => tc.and_then(|x| x.codons.clone()).unwrap_or_default(),
            "Existing_variation" => variant.existing_variation.join(list_sep),
            "ALLELE_NUM" => variant
                .allele_num
                .unwrap_or(variant.allele_index as u32 + 1)
                .to_string(),
            "IND" => variant
                .individual
                .as_deref()
                .unwrap_or_default()
                .to_string(),
            // The `--individual_zyg` list, every named sample's `sample:ZYG`
            // (OutputFactory.pm:981-987), joined as the format joins a list; else the
            // copy's own zygosity under `--individual` (:974-976).
            "ZYG" => match &variant.sample_zygosities {
                Some(samples) => samples.joined(list_sep).into_owned(),
                None => variant.zygosity().unwrap_or_default(),
            },
            // A structural row carries neither `REF_ALLELE` nor `UPLOADED_ALLELE`
            // (OutputFactory.pm `BaseStructuralVariationOverlapAllele_to_output_hash`).
            "REF_ALLELE" if variant.is_structural => String::new(),
            // Every row of a `--minimal` split prints the record feature's reference:
            // the record's own when the first feature wrote no row before the rejoin
            // (`split_first_feature_rowless`), else the first feature's.
            "REF_ALLELE" => match variant.record_allele_string_multi.as_deref() {
                Some(record) if variant.split_first_feature_rowless => {
                    record.split('/').next().unwrap_or_default().to_string()
                }
                Some(record) if variant.original_start.is_some() => {
                    vep_core::minimal::rejoined_reference(record, options.process_ref_homs)
                }
                _ => String::from_utf8_lossy(&variant.ref_allele).to_string(),
            },
            "UPLOADED_ALLELE" => variant.uploaded_allele(options.minimal).unwrap_or_default(),
            "DISTANCE" => tc
                .filter(|x| options.distance_shown(x))
                .and_then(|x| x.distance.map(|d| d.to_string()))
                .unwrap_or_default(),
            "OverlapBP" => sv_overlap(variant, tc)
                .map(|(bp, _)| bp.to_string())
                .unwrap_or_default(),
            "OverlapPC" => sv_overlap(variant, tc)
                .map(|(_, pc)| format!("{pc:.2}"))
                .unwrap_or_default(),
            "STRAND" => tc.map(|x| x.strand.to_string()).unwrap_or_default(),
            "FLAGS" => tc
                .map(|x| {
                    x.flags
                        .iter()
                        .filter(|f| f.starts_with("cds_"))
                        .cloned()
                        .collect::<Vec<_>>()
                        .join(list_sep)
                })
                .unwrap_or_default(),
            "PICK" => {
                // An intergenic allele is flagged like a transcript allele
                // (`VariationFeatureOverlapAllele_to_output_hash`, OutputFactory.pm:1363).
                let picked = tc.map_or(variant.intergenic_pick, |x| {
                    x.flags.iter().any(|f| f == "PICK")
                });
                if options.include_pick && picked {
                    "1".to_string()
                } else {
                    String::new()
                }
            }
            "VARIANT_CLASS" => variant.class_so_term(),
            "AMBIGUITY" => variant
                .ambiguity_code()
                .map(|code| char::from(code).to_string())
                .unwrap_or_default(),
            "MINIMISED" => {
                if variant.minimised {
                    "1".to_string()
                } else {
                    String::new()
                }
            }
            "SYMBOL_SOURCE" => string_identifier(tc.and_then(|x| x.gene_symbol_source.as_deref())),
            "HGNC_ID" => string_identifier(tc.and_then(|x| x.hgnc_id.as_deref())),
            "CANONICAL" => {
                if tc.map(|x| x.canonical).unwrap_or(false) {
                    "YES".to_string()
                } else {
                    String::new()
                }
            }
            // The MANE sets the transcript belongs to, each listed under its own
            // flag state (OutputFactory.pm:1587-1601): `--mane_select` alone never
            // lists `MANE_Plus_Clinical`.
            "MANE" => {
                let mut sets = Vec::new();
                if options.mane_select && tc.is_some_and(|x| x.mane_select.is_some()) {
                    sets.push("MANE_Select");
                }
                if options.mane_plus_clinical && tc.is_some_and(|x| x.mane_plus_clinical.is_some())
                {
                    sets.push("MANE_Plus_Clinical");
                }
                sets.join(list_sep)
            }
            "MANE_SELECT" => tc
                .and_then(|x| x.mane_select.as_deref())
                .filter(|v| perl_true(v))
                .unwrap_or("")
                .to_string(),
            "MANE_PLUS_CLINICAL" => tc
                .and_then(|x| x.mane_plus_clinical.as_deref())
                .filter(|v| perl_true(v))
                .unwrap_or("")
                .to_string(),
            "TSL" => tc
                .and_then(|x| x.tsl.map(|t| t.to_string()))
                .unwrap_or_default(),
            "GENCODE_PRIMARY" => {
                if tc.is_some_and(|x| x.gencode_primary) {
                    "1".to_string()
                } else {
                    String::new()
                }
            }
            // `s/principal/P/; s/alternative/A/` (OutputFactory.pm:1616-1622): one
            // substitution each, the rest of the value kept, so `principalM` is `PM`.
            "APPRIS" => tc
                .and_then(|x| x.appris.as_deref())
                .filter(|v| perl_true(v))
                .map(|v| {
                    v.replacen("principal", "P", 1)
                        .replacen("alternative", "A", 1)
                })
                .unwrap_or_default(),
            "CCDS" => string_identifier(tc.and_then(|x| x.ccds.as_deref())),
            "ENSP" => string_identifier(tc.and_then(|x| x.protein_id.as_deref())),
            "SWISSPROT" => {
                string_identifier_list(tc.and_then(|x| x.swissprot.as_deref()), list_sep)
            }
            "TREMBL" => string_identifier_list(tc.and_then(|x| x.trembl.as_deref()), list_sep),
            "UNIPARC" => string_identifier_list(tc.and_then(|x| x.uniparc.as_deref()), list_sep),
            "UNIPROT_ISOFORM" => {
                string_identifier_list(tc.and_then(|x| x.uniprot_isoform.as_deref()), list_sep)
            }
            "HGVSg" => variant.hgvsg.clone().unwrap_or_default(),
            "NEAREST" => variant
                .nearest
                .as_deref()
                .map(|values| values.join(list_sep))
                .unwrap_or_default(),
            "SPDI" => row_spdi(variant, tc, options)
                .map(|(spdi, _)| spdi.to_string())
                .unwrap_or_default(),
            "GA4GH_VRS" => row_ga4gh_vrs(variant, tc, options)
                .map(VrsAllele::encoded)
                .unwrap_or_default(),
            "CHECK_REF" => {
                if variant.check_ref_failed {
                    "failed".to_string()
                } else {
                    String::new()
                }
            }
            "RefSeq" => string_identifier_list(tc.and_then(|x| x.refseq.as_deref()), list_sep),
            "GENE_PHENO" => {
                if tc.is_some_and(|x| x.gene_phenotype) {
                    "1".to_string()
                } else {
                    String::new()
                }
            }
            "miRNA" => tc.map(|x| x.mirna.join(list_sep)).unwrap_or_default(),
            // `join(':', grep {$_} (label, identifier))` (OutputFactory.pm:1500): a
            // part that is empty or `0` is false to Perl and left out, so a feature
            // with no label prints its identifier alone.
            "DOMAINS" => tc
                .map(|x| {
                    x.domains
                        .iter()
                        .map(|(label, id)| match (perl_true(label), perl_true(id)) {
                            (true, true) => format!("{label}:{id}"),
                            (false, true) => id.clone(),
                            (true, false) => label.clone(),
                            (false, false) => String::new(),
                        })
                        .collect::<Vec<_>>()
                        .join(list_sep)
                })
                .unwrap_or_default(),
            "HGVS_OFFSET" => tc
                .and_then(|x| x.hgvs_offset.map(|o| o.to_string()))
                .unwrap_or_default(),
            // Set on every sequence variant's transcript row, 0 included, when both
            // flags are on (`TranscriptVariationAllele_to_output_hash`,
            // OutputFactory.pm:1690-1696); an intergenic row carries none, and a
            // structural record's rows none (`TranscriptStructuralVariationAllele_to_output_hash`,
            // :2118-2125, sets no such key).
            "SHIFT_LENGTH" if options.shift_3prime && !variant.is_structural => {
                tc.map(|x| x.shift_length.to_string()).unwrap_or_default()
            }
            "SHIFT_LENGTH" => String::new(),
            // Both are the reference allele's object's (OutputFactory.pm:1774-1779),
            // so a row without a transcript carries neither.
            "USED_REF" => tc
                .and_then(|x| x.used_ref.as_deref())
                .unwrap_or_default()
                .to_string(),
            "GIVEN_REF" => {
                if tc.is_some_and(|x| x.used_ref.is_some()) {
                    String::from_utf8_lossy(&variant.ref_allele).to_string()
                } else {
                    String::new()
                }
            }
            "CLIN_SIG"
            | "SOMATIC"
            | "PHENO"
            | "PUBMED"
            | "VAR_SYNONYMS"
            | "CLINVAR_SOMATIC_CLASSIFICATION" => {
                colocated_column(field, variant, options, list_sep)
            }
            "SV" => variant.overlapping_svs.join(list_sep),
            "FREQS" => variant.freqs.join(list_sep),
            "SIFT" => tc.and_then(|x| x.sift.clone()).unwrap_or_default(),
            "PolyPhen" => tc.and_then(|x| x.polyphen.clone()).unwrap_or_default(),
            "MOTIF_NAME"
            | "MOTIF_POS"
            | "HIGH_INF_POS"
            | "MOTIF_SCORE_CHANGE"
            | "TRANSCRIPTION_FACTORS" => String::new(),
            "MAX_AF" => frequency::max_af(variant, &options.frequency_groups())
                .map(|(max, _)| max)
                .unwrap_or_default(),
            "MAX_AF_POPS" => frequency::max_af(variant, &options.frequency_groups())
                .map(|(_, populations)| populations.join(list_sep))
                .unwrap_or_default(),
            _ => match frequency::frequency_column(field) {
                Some(column) => {
                    frequency::frequency_field(column, variant, &options.frequency_groups())
                        .join(list_sep)
                }
                None => tc
                    .and_then(|x| x.plugin_data.get(field))
                    .or_else(|| variant.plugin_data.get(field))
                    .cloned()
                    .unwrap_or_default(),
            },
        }
    }

    /// A position column with the total `--total_length` appends: `428/1200`, or
    /// `-/1200` where the position is undefined but the total applies, since
    /// Ensembl VEP appends to its `-` placeholder (`OutputFactory.pm`
    /// `TranscriptVariationAllele_to_output_hash`); without a total the position
    /// alone, empty when undefined.
    pub fn string_position_column(position: Option<&str>, total: Option<u32>) -> String {
        match total {
            Some(total) => format!("{}/{total}", position.unwrap_or("-")),
            None => position.unwrap_or_default().to_string(),
        }
    }

    fn assert_writer_matches(
        variant: &InputVariant,
        tc: Option<&TranscriptConsequence>,
        options: &FieldOptions,
        list_sep: &str,
    ) -> Vec<&'static str> {
        let mut populated = Vec::new();
        for column in every_column() {
            let expected = string_field_value(column, variant, tc, options, list_sep);
            let mut out = b"prefix".to_vec();
            write_field_value(&mut out, column, variant, tc, options, list_sep);
            assert_eq!(
                String::from_utf8(out).unwrap(),
                format!("prefix{expected}"),
                "{column} with list separator {list_sep:?}"
            );
            assert_eq!(
                field_value(column, variant, tc, options, list_sep),
                expected,
                "{column}"
            );
            if !expected.is_empty() {
                populated.push(column);
            }
        }
        populated
    }

    /// Every byte from `0x00` to `0x1F` and `0x7F` written from `start` on is
    /// rewritten as its upper-case `%XX`; the bytes before `start`, a `%`
    /// already there, a byte outside ASCII and a value with no control byte
    /// stay as written.
    #[test]
    fn percent_encode_control_bytes_rewrites_the_value_alone() {
        let mut out = b"kept\t".to_vec();
        out.extend_from_slice("a\0b\tc\nd\x1Fe\x7F%20\u{E9}".as_bytes());
        percent_encode_control_bytes(&mut out, 5);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "kept\ta%00b%09c%0Ad%1Fe%7F%20\u{E9}"
        );
        let mut plain = b"a\tplain value".to_vec();
        percent_encode_control_bytes(&mut plain, 2);
        assert_eq!(plain, b"a\tplain value");
        for byte in 0u8..=0xFF {
            let mut alone = vec![byte];
            percent_encode_control_bytes(&mut alone, 0);
            let expected = if byte < 0x20 || byte == 0x7F {
                format!("%{byte:02X}").into_bytes()
            } else {
                vec![byte]
            };
            assert_eq!(alone, expected, "{byte:#04x}");
        }
    }

    /// The byte writer agrees with the string table on every column, with both
    /// list separators, with and without HGVS escaping, under the whole-variant
    /// ClinVar rule, with the own-span positions, in the display vocabulary, on
    /// a bare row, an empty consequence and the default options; and between the
    /// populated sequence row and the structural row every column but the
    /// regulatory ones and the name nothing produces prints a value.
    #[test]
    fn write_field_value_matches_the_string_table_on_every_column() {
        let variant = populated_variant();
        let tc = populated_consequence();
        let structural = structural_variant();
        let options = all_options();
        let mut populated: Vec<&str> = Vec::new();
        for sep in [",", "&"] {
            populated.extend(assert_writer_matches(&variant, Some(&tc), &options, sep));
            populated.extend(assert_writer_matches(&variant, None, &options, sep));
            populated.extend(assert_writer_matches(&structural, Some(&tc), &options, sep));
            populated.extend(assert_writer_matches(&structural, None, &options, sep));
        }
        for column in every_column() {
            let never = matches!(
                column,
                "MOTIF_NAME"
                    | "MOTIF_POS"
                    | "HIGH_INF_POS"
                    | "MOTIF_SCORE_CHANGE"
                    | "TRANSCRIPTION_FACTORS"
                    | "Unknown"
            );
            assert_eq!(populated.contains(&column), !never, "{column}");
        }
        for options in [
            FieldOptions {
                no_escape: true,
                ..all_options()
            },
            FieldOptions {
                clin_sig_allele: ClinSigAllele::WholeVariant,
                ..all_options()
            },
            FieldOptions {
                positions_of_own_span: true,
                ..all_options()
            },
            FieldOptions {
                terms: Vocabulary::Display,
                ..all_options()
            },
            FieldOptions {
                shift_3prime: false,
                shift_genomic: false,
                ..all_options()
            },
            FieldOptions::default(),
        ] {
            assert_writer_matches(&variant, Some(&tc), &options, ",");
            assert_writer_matches(&structural, Some(&tc), &options, ",");
        }
        let bare = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let empty = TranscriptConsequence::default();
        for options in [all_options(), FieldOptions::default()] {
            assert_writer_matches(&bare, Some(&empty), &options, ",");
            assert_writer_matches(&bare, None, &options, ",");
            assert_writer_matches(&variant, Some(&empty), &options, ",");
        }
        // An individual's copy without the `--individual_zyg` list prints its own
        // zygosity.
        let mut own_copy = variant.clone();
        own_copy.sample_zygosities = None;
        assert_writer_matches(&own_copy, Some(&tc), &options, ",");
        assert_eq!(
            field_value("ZYG", &own_copy, Some(&tc), &options, ","),
            "HET"
        );
        // The reference row of `--process_ref_homs` prints its transcript's
        // allele under `--use_transcript_ref`, the record's without a `used_ref`.
        let mut reference = variant.clone();
        reference.reference_allele_row = true;
        for options in [all_options(), FieldOptions::default()] {
            assert_writer_matches(&reference, Some(&tc), &options, ",");
            assert_writer_matches(&reference, Some(&empty), &options, ",");
        }
        assert_eq!(
            field_value("Allele", &reference, Some(&tc), &options, ","),
            "C"
        );
        assert_eq!(
            field_value("Allele", &reference, Some(&empty), &options, ","),
            "G"
        );
        assert_eq!(
            field_value("HGVSp", &variant, Some(&tc), &options, ","),
            "ENSP00000300.1:p.Ala101%3D;p.%3D"
        );
        assert_eq!(
            field_value("cDNA_position", &variant, Some(&tc), &options, ","),
            "428/1200"
        );
        assert_eq!(
            field_value("Feature", &variant, Some(&tc), &options, ","),
            "ENST00000300.4"
        );
        assert_eq!(
            field_value("DOMAINS", &variant, Some(&tc), &options, "&"),
            "Pfam:PF00080&PANTHER:PTHR10003&PS50047&Low_complexity_(Seg)&"
        );
        assert_eq!(
            field_value("TREMBL", &variant, Some(&tc), &options, ","),
            "A0A0"
        );
        assert_eq!(
            field_value("ZYG", &variant, Some(&tc), &options, "&"),
            "S1:HET&S2:HOMREF"
        );
        assert_eq!(
            field_value("CLIN_SIG", &variant, Some(&tc), &options, ","),
            "pathogenic"
        );
        assert_eq!(
            field_value("OverlapPC", &structural, Some(&tc), &options, ","),
            "10.09"
        );
        assert_eq!(
            field_value("STRAND", &bare, Some(&empty), &options, ","),
            "1"
        );
        assert_eq!(field_value("ALLELE_NUM", &bare, None, &options, ","), "1");
    }

    /// The byte writer agrees with the string table on a `--summary` row, whose
    /// columns are the record's identity, location, consequence list and
    /// variant-level fields alone, and on a flagged intergenic row.
    #[test]
    fn write_field_value_matches_the_string_table_on_summary_and_intergenic_rows() {
        let mut summary = populated_variant();
        summary.summary_terms = Some(
            [
                Consequence::MissenseVariant,
                Consequence::IntronVariant,
                Consequence::SpliceRegionVariant,
            ]
            .into_iter()
            .collect(),
        );
        for sep in [",", "&"] {
            assert_writer_matches(&summary, None, &all_options(), sep);
        }
        let display = FieldOptions {
            terms: Vocabulary::Display,
            ..all_options()
        };
        assert_writer_matches(&summary, None, &display, ",");
        assert_eq!(
            field_value("Consequence", &summary, None, &display, "&"),
            "NON_SYNONYMOUS_CODING&INTRONIC&SPLICE_SITE"
        );
        assert_eq!(
            field_value("Location", &summary, None, &all_options(), ","),
            "chr21:10499-10500"
        );
        assert_eq!(
            field_value("Allele", &summary, None, &all_options(), ","),
            ""
        );
        assert_eq!(
            field_value("IMPACT", &summary, None, &all_options(), ","),
            ""
        );
        assert_eq!(
            field_value("ZYG", &summary, None, &all_options(), ","),
            "S1:HET,S2:HOMREF"
        );
        assert_eq!(
            field_value("VPlugin", &summary, None, &all_options(), ","),
            ""
        );

        let mut intergenic = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        intergenic.most_severe_consequence = Some(Consequence::IntergenicVariant);
        intergenic.intergenic_pick = true;
        intergenic.intergenic_shift = Some(2);
        assert_writer_matches(&intergenic, None, &all_options(), ",");
        assert_eq!(
            field_value("PICK", &intergenic, None, &all_options(), ","),
            "1"
        );
        assert_eq!(
            field_value("Location", &intergenic, None, &all_options(), ","),
            "21:102"
        );
    }

    fn ranked(id: &str, rank: u32) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript_id: Arc::from(id),
            id_rank: rank,
            ..TranscriptConsequence::default()
        }
    }

    /// The reference order of the feature rows: the byte order of the rejoin
    /// keys, then the reference row, then the allele index, stable.
    fn key_sorted<'a>(record: &[&'a InputVariant]) -> Vec<Row<'a>> {
        let mut rows: Vec<Row<'a>> = Vec::new();
        for variant in record {
            for tc in &variant.transcript_consequences {
                rows.push((variant, Some(tc)));
            }
        }
        rows.sort_by(|a, b| {
            rejoin_key(a)
                .cmp(rejoin_key(b))
                .then(b.0.reference_allele_row.cmp(&a.0.reference_allele_row))
                .then(a.0.allele_index.cmp(&b.0.allele_index))
        });
        rows
    }

    fn feature_rows<'a>(rows: &[Row<'a>]) -> Vec<(*const TranscriptConsequence, usize)> {
        rows.iter()
            .filter_map(|(v, t)| t.map(|t| (t as *const TranscriptConsequence, v.allele_index)))
            .collect()
    }

    fn assert_rank_order_is_key_order(record: &[&InputVariant]) {
        let mut rows = Vec::new();
        record_rows_into(record, &mut rows);
        let feature_count = rows.iter().filter(|r| r.1.is_some()).count();
        assert_eq!(
            feature_rows(&rows[..feature_count]),
            feature_rows(&key_sorted(record))
        );
        for a in &rows[..feature_count] {
            for b in &rows[..feature_count] {
                assert_eq!(
                    transcript_order(a, b),
                    rejoin_key(a).cmp(rejoin_key(b)),
                    "{} against {}",
                    a.1.unwrap().transcript_id,
                    b.1.unwrap().transcript_id
                );
            }
        }
    }

    /// Rows of one chromosome whose ranks follow the loader's byte order of the
    /// ids, with a transcript shared by two alleles and one id carried by two
    /// consequences of one allele, order as the id strings do; the two rows of
    /// one id and allele keep their input order.
    #[test]
    fn ranked_rows_order_as_their_ids_do() {
        let mut a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let mut b = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"T".to_vec());
        b.allele_index = 1;
        let mut first_copy = ranked("ENST00000200", 2);
        first_copy.distance = Some(1);
        let mut second_copy = ranked("ENST00000200", 2);
        second_copy.distance = Some(2);
        a.transcript_consequences = vec![
            ranked("ENST00000300", 3),
            first_copy,
            ranked("ENST0000002", 1),
            second_copy,
            ranked("ENST00000010", 0),
        ];
        b.transcript_consequences = vec![ranked("ENST00000200", 2), ranked("ENST00000010", 0)];
        assert_rank_order_is_key_order(&[&a, &b]);

        let mut rows = Vec::new();
        record_rows_into(&[&a, &b], &mut rows);
        let distances: Vec<Option<u64>> = rows
            .iter()
            .filter(|(v, t)| v.allele_index == 0 && &*t.unwrap().transcript_id == "ENST00000200")
            .map(|(_, t)| t.unwrap().distance)
            .collect();
        assert_eq!(distances, vec![Some(1), Some(2)]);
    }

    /// A consequence built without a loaded transcript carries rank 0 and
    /// compares by key against every neighbour, whatever their ranks.
    #[test]
    fn unranked_rows_fall_back_to_the_key_order() {
        let mut a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        a.transcript_consequences = vec![
            ranked("ENST00000300", 3),
            ranked("ENST00000250", 0),
            ranked("ENST00000200", 2),
            ranked("ENST00000900", 0),
            ranked("ENST00000010", 1),
        ];
        assert_rank_order_is_key_order(&[&a]);
    }

    /// A breakend's rows, whose mate consequences are ranked on the mate's
    /// chromosome so a rank order can contradict the id order, and the rows of
    /// a `--minimal` allele whose rejoin prefix sorts before the ids compare by
    /// key; the reference row of `--process_ref_homs` still leads its feature.
    #[test]
    fn rows_ranked_elsewhere_fall_back_to_the_key_order() {
        let mut breakend =
            InputVariant::new("21".into(), 100, 100, b"N".to_vec(), b"N[1:500[".to_vec());
        breakend.variant_class = VariantClass::Translocation;
        breakend.is_structural = true;
        breakend.mate_chr = Some("1".to_string());
        breakend.transcript_consequences = vec![
            ranked("ENST00000030", 2),
            ranked("ENST00000900", 1),
            ranked("ENST00000020", 1),
        ];
        assert_rank_order_is_key_order(&[&breakend]);
        let mut rows = Vec::new();
        record_rows_into(&[&breakend], &mut rows);
        let ids: Vec<&str> = rows
            .iter()
            .map(|(_, t)| &*t.unwrap().transcript_id)
            .collect();
        assert_eq!(ids, vec!["ENST00000020", "ENST00000030", "ENST00000900"]);

        let mut trimmed = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        trimmed.rejoin_key_prefix = Some(Arc::from("1:100-100:A/G"));
        trimmed.transcript_consequences =
            vec![ranked("ENST00000900", 3), ranked("ENST00000020", 1)];
        let mut plain = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"T".to_vec());
        plain.allele_index = 1;
        plain.transcript_consequences = vec![ranked("ENST00000500", 2), ranked("ENST00000020", 1)];
        let mut reference = plain.reference_allele_row();
        reference.transcript_consequences = vec![ranked("ENST00000500", 2)];
        assert_rank_order_is_key_order(&[&trimmed, &plain, &reference]);
        let mut rows = Vec::new();
        record_rows_into(&[&trimmed, &plain, &reference], &mut rows);
        let order: Vec<(&str, usize, bool)> = rows
            .iter()
            .map(|(v, t)| {
                (
                    &*t.unwrap().transcript_id,
                    v.allele_index,
                    v.reference_allele_row,
                )
            })
            .collect();
        assert_eq!(
            order,
            vec![
                ("ENST00000020", 0, false),
                ("ENST00000900", 0, false),
                ("ENST00000020", 1, false),
                ("ENST00000500", 0, true),
                ("ENST00000500", 1, false),
            ]
        );
    }

    #[test]
    fn user_group_fields_come_first_and_once() {
        let opts = FieldOptions {
            symbol: true,
            af_1kg: true,
            af: true,
            ..Default::default()
        };
        let f = flag_fields(&opts);
        assert_eq!(&f[..4], &["IMPACT", "DISTANCE", "STRAND", "FLAGS"]);
        assert_eq!(&f[4..7], &["SYMBOL", "SYMBOL_SOURCE", "HGNC_ID"]);
        assert_eq!(f.iter().filter(|x| **x == "AF").count(), 1);
    }

    /// `--af_gnomad` and `--af_gnomade` name the same ten gnomAD exome fields
    /// (Constants.pm:120-121); `flag_fields` lists them once when both groups
    /// are on (OutputFactory.pm:2516-2538 uniquifies; the CLI refuses the pair
    /// as Config.pm:619 does), no frequency group is filled from the variant
    /// alone, and the JSON factory stubs every one of them out.
    #[test]
    fn gnomad_exome_fields_print_under_either_flag_and_once_under_both() {
        let gnomad = FieldOptions {
            af_gnomad: true,
            ..Default::default()
        };
        let gnomade = FieldOptions {
            af_gnomade: true,
            ..Default::default()
        };
        let both = FieldOptions {
            af_gnomad: true,
            af_gnomade: true,
            ..Default::default()
        };
        assert_eq!(&flag_fields(&gnomad)[4..], GNOMADE_FIELDS);
        assert_eq!(flag_fields(&gnomade), flag_fields(&gnomad));
        assert_eq!(flag_fields(&both), flag_fields(&gnomad));
        assert!(!flag_fields(&FieldOptions::default())
            .iter()
            .any(|f| f.starts_with("gnomADe_")));
        for field in GNOMADE_FIELDS {
            assert!(field_is_produced(field, &gnomad), "{field}");
            assert!(field_is_produced(field, &gnomade), "{field}");
            assert!(
                !field_is_produced(field, &FieldOptions::default()),
                "{field}"
            );
            assert!(colocated_field(field), "{field}");
        }
        let every_frequency_group = FieldOptions {
            af: true,
            af_1kg: true,
            af_gnomad: true,
            af_gnomade: true,
            af_gnomadg: true,
            max_af: true,
            ..Default::default()
        };
        assert_eq!(
            vf_level_flag_fields(&every_frequency_group),
            Vec::<&str>::new()
        );

        let mut cv = known("rs1", "A/G", &[("G", "G")], &[], None);
        cv.frequency_strings = [
            ("gnomADe".to_string(), "G:0.1234".to_string()),
            ("gnomADe_NFE".to_string(), "G:0.2".to_string()),
        ]
        .into();
        let alleles = record("A", &["G"], vec![cv]);
        for options in [&gnomad, &gnomade, &both] {
            assert_eq!(
                field_value("gnomADe_AF", &alleles[0], None, options, ","),
                "0.1234"
            );
            assert_eq!(
                field_value("gnomADe_NFE_AF", &alleles[0], None, options, ","),
                "0.2"
            );
            assert_eq!(
                field_value("gnomADe_AFR_AF", &alleles[0], None, options, ","),
                ""
            );
        }
        assert_eq!(
            field_value(
                "gnomADe_AF",
                &alleles[0],
                None,
                &FieldOptions::default(),
                ","
            ),
            ""
        );
    }

    /// Under `--shift_3prime 1 --shift_genomic 1` a transcript row's Location moves
    /// by its shift along the transcript strand and an inserted allele prints
    /// rotated along the move (OutputFactory.pm:1349-1352, 1702-1706); an
    /// intergenic row moves towards higher coordinates by the forward-strand
    /// shift (:1992-1998); a row that moved nothing, and every row without both
    /// flags, prints the input position and allele.
    #[test]
    fn shift_genomic_moves_location_and_allele_for_shifted_rows() {
        let both = FieldOptions {
            shift_3prime: true,
            shift_genomic: true,
            ..Default::default()
        };
        let insertion =
            InputVariant::new("21".into(), 1_001, 1_000, b"-".to_vec(), b"CAG".to_vec());
        let forward = TranscriptConsequence {
            strand: 1,
            shift_length: 4,
            ..Default::default()
        };
        assert_eq!(
            field_value("Location", &insertion, Some(&forward), &both, ","),
            "21:1004-1005"
        );
        assert_eq!(
            field_value("Allele", &insertion, Some(&forward), &both, ","),
            "AGC"
        );
        let reverse = TranscriptConsequence {
            strand: -1,
            shift_length: 1,
            ..Default::default()
        };
        assert_eq!(
            field_value("Location", &insertion, Some(&reverse), &both, ","),
            "21:999-1000"
        );
        assert_eq!(
            field_value("Allele", &insertion, Some(&reverse), &both, ","),
            "GCA"
        );
        let unmoved = TranscriptConsequence::default();
        assert_eq!(
            field_value("Location", &insertion, Some(&unmoved), &both, ","),
            "21:1000-1001"
        );
        assert_eq!(
            field_value("Allele", &insertion, Some(&unmoved), &both, ","),
            "CAG"
        );

        let mut deletion =
            InputVariant::new("21".into(), 2_000, 2_001, b"AT".to_vec(), b"-".to_vec());
        deletion.intergenic_shift = Some(3);
        assert_eq!(
            field_value("Location", &deletion, None, &both, ","),
            "21:2003-2004"
        );
        assert_eq!(field_value("Allele", &deletion, None, &both, ","), "-");

        for options in [
            FieldOptions {
                shift_3prime: true,
                ..Default::default()
            },
            FieldOptions {
                shift_genomic: true,
                ..Default::default()
            },
        ] {
            assert_eq!(
                field_value("Location", &insertion, Some(&forward), &options, ","),
                "21:1000-1001"
            );
            assert_eq!(
                field_value("Allele", &insertion, Some(&forward), &options, ","),
                "CAG"
            );
            assert_eq!(
                field_value("Location", &deletion, None, &options, ","),
                "21:2000-2001"
            );
        }
    }

    /// `--spdi` and `--ga4gh_vrs` list their columns after `MINIMISED` and before
    /// `SYMBOL` (Constants.pm:79-85); the values are the allele's notation and
    /// the encoded object (OutputFactory.pm:1400-1410), on a transcript row and
    /// on an intergenic row alike, and an allele with no object prints an empty
    /// `GA4GH_VRS`. Both are set in the overlap-allele hash alone, never by
    /// `VariationFeature_to_output_hash` (:883-932), so neither is a
    /// variant-level field and a `--summary` or `--most_severe` row carries
    /// neither column.
    #[test]
    fn spdi_and_ga4gh_vrs_print_the_allele_notation_and_the_encoded_object() {
        let both = FieldOptions {
            minimal: true,
            spdi: true,
            ga4gh_vrs: true,
            symbol: true,
            ..Default::default()
        };
        assert_eq!(
            &flag_fields(&both)[4..8],
            &["MINIMISED", "SPDI", "GA4GH_VRS", "SYMBOL"]
        );
        assert!(!flag_fields(&FieldOptions::default()).contains(&"SPDI"));
        let mut snv = InputVariant::new(
            "21".into(),
            10538885,
            10538885,
            b"A".to_vec(),
            b"C".to_vec(),
        );
        snv.spdi = Some("21:10538884:A:C".to_string());
        let tc = TranscriptConsequence::default();
        assert_eq!(
            field_value("SPDI", &snv, Some(&tc), &both, ","),
            "21:10538884:A:C"
        );
        assert_eq!(
            field_value("SPDI", &snv, None, &both, ","),
            "21:10538884:A:C"
        );
        assert_eq!(field_value("GA4GH_VRS", &snv, Some(&tc), &both, ","), "");
        snv.ga4gh_vrs = Some(vep_core::spdi::VrsAllele {
            sequence_id: "NC_000021.9".to_string(),
            start: 10538884,
            end: 10538885,
            sequence: "C".to_string(),
        });
        let encoded = concat!(
            r#"{"location":{"type":"SequenceLocation","interval":{"type":"SimpleInterval","#,
            r#""start":"10538884","end":10538885},"sequence_id":"refseq:NC_000021.9"},"#,
            r#""state":{"type":"SequenceState","sequence":"C"},"type":"Allele"}"#
        );
        assert_eq!(field_value("GA4GH_VRS", &snv, None, &both, ","), encoded);
        let bare = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"C".to_vec());
        assert_eq!(field_value("SPDI", &bare, Some(&tc), &both, ","), "");

        assert_eq!(vf_level_flag_fields(&both), ["MINIMISED"]);
        let mut summary = snv.clone();
        summary.most_severe_consequence = Some(Consequence::MissenseVariant);
        summary.summary_terms = Some(smallvec::smallvec![Consequence::MissenseVariant]);
        assert_eq!(field_value("SPDI", &summary, None, &both, ","), "");
        assert_eq!(field_value("GA4GH_VRS", &summary, None, &both, ","), "");
        assert_eq!(
            field_value("SPDI", &summary, Some(&tc), &both, ","),
            "21:10538884:A:C"
        );
        assert_eq!(
            field_value("GA4GH_VRS", &summary, Some(&tc), &both, ","),
            encoded
        );
    }

    /// Ensembl VEP looks the SPDI up by the row's Allele among the record's
    /// alleles as written (OutputFactory.pm:1400): under `--shift_3prime 1
    /// --shift_genomic 1` a row whose inserted bases a shift rotated
    /// (`21 13620283 G GTCA` moved 10 bases prints `CAT`) carries no `SPDI` and no
    /// `GA4GH_VRS`, unless the rotation lands on a sibling allele of the record
    /// (`21 13620283 G GTCA,GCAT`), whose notations it then prints as Ensembl VEP
    /// does (`21:13620283::CAT`); a row whose rotation returns the bases as given
    /// (`21 10451779 A ATA` moved 6) and every unshifted row carry their own.
    #[test]
    fn spdi_is_absent_on_a_row_whose_allele_a_genomic_shift_rotated() {
        let options = FieldOptions {
            spdi: true,
            ga4gh_vrs: true,
            shift_3prime: true,
            shift_genomic: true,
            ..Default::default()
        };
        let mut tri = InputVariant::new(
            "21".into(),
            13620284,
            13620283,
            b"-".to_vec(),
            b"TCA".to_vec(),
        );
        tri.spdi = Some("21:13620283::TCA".to_string());
        tri.ga4gh_vrs = Some(vep_core::spdi::VrsAllele {
            sequence_id: "NC_000021.9".to_string(),
            start: 13620283,
            end: 13620283,
            sequence: "TCA".to_string(),
        });
        let moved_ten = TranscriptConsequence {
            strand: 1,
            shift_length: 10,
            ..Default::default()
        };
        assert_eq!(
            field_value("Allele", &tri, Some(&moved_ten), &options, ","),
            "CAT"
        );
        assert_eq!(
            field_value("SPDI", &tri, Some(&moved_ten), &options, ","),
            ""
        );
        assert_eq!(
            field_value("GA4GH_VRS", &tri, Some(&moved_ten), &options, ","),
            ""
        );
        let unmoved = TranscriptConsequence::default();
        assert_eq!(
            field_value("SPDI", &tri, Some(&unmoved), &options, ","),
            "21:13620283::TCA"
        );
        tri.sibling_spdi = vec![vep_core::spdi::SiblingSpdi {
            allele: "CAT".to_string(),
            spdi: "21:13620283::CAT".to_string(),
            ga4gh_vrs: Some(vep_core::spdi::VrsAllele {
                sequence_id: "NC_000021.9".to_string(),
                start: 13620283,
                end: 13620283,
                sequence: "CAT".to_string(),
            }),
        }];
        assert_eq!(
            field_value("SPDI", &tri, Some(&moved_ten), &options, ","),
            "21:13620283::CAT"
        );
        assert!(
            field_value("GA4GH_VRS", &tri, Some(&moved_ten), &options, ",")
                .contains("\"sequence\":\"CAT\""),
        );
        assert_eq!(
            field_value("SPDI", &tri, Some(&unmoved), &options, ","),
            "21:13620283::TCA"
        );
        let mut di = InputVariant::new(
            "21".into(),
            10451780,
            10451779,
            b"-".to_vec(),
            b"TA".to_vec(),
        );
        di.spdi = Some("21:10451779::TA".to_string());
        let moved_six = TranscriptConsequence {
            strand: 1,
            shift_length: 6,
            ..Default::default()
        };
        assert_eq!(
            field_value("Allele", &di, Some(&moved_six), &options, ","),
            "TA"
        );
        assert_eq!(
            field_value("SPDI", &di, Some(&moved_six), &options, ","),
            "21:10451779::TA"
        );
    }

    /// A `--summary` or `--most_severe` row is the record's hash
    /// (`VariationFeature_to_output_hash`, OutputFactory.pm:883-932), whose
    /// Location is the input span: the `--shift_genomic 1` move belongs to the
    /// feature and intergenic row builders (:1702-1706, :1992-1998), so the row
    /// of a record whose intergenic allele shifted keeps the input position.
    #[test]
    fn summary_row_location_is_the_input_span_under_shift_genomic() {
        let both = FieldOptions {
            shift_3prime: true,
            shift_genomic: true,
            ..Default::default()
        };
        let mut deletion =
            InputVariant::new("21".into(), 2_000, 2_001, b"AT".to_vec(), b"-".to_vec());
        deletion.intergenic_shift = Some(3);
        deletion.most_severe_consequence = Some(Consequence::IntergenicVariant);
        deletion.summary_terms = Some(smallvec::smallvec![Consequence::IntergenicVariant]);
        assert_eq!(
            field_value("Location", &deletion, None, &both, ","),
            "21:2000-2001"
        );
        assert_eq!(
            field_value("Consequence", &deletion, None, &both, ","),
            "intergenic_variant"
        );
        assert_eq!(field_value("Allele", &deletion, None, &both, ","), "");
    }

    /// `--shift_length` lists SHIFT_LENGTH among the flag fields, and under
    /// `--shift_3prime 1` every transcript row of a sequence variant carries its
    /// value, 0 included (OutputFactory.pm:1690-1696); an intergenic row carries
    /// none, and so does a structural record's transcript row, whose output hash
    /// (:2118-2125) never sets the key.
    #[test]
    fn shift_length_prints_zero_on_unshifted_transcript_rows() {
        let options = FieldOptions {
            shift_length: true,
            shift_3prime: true,
            ..Default::default()
        };
        assert!(flag_fields(&options).contains(&"SHIFT_LENGTH"));
        let snv = InputVariant::new("21".into(), 1_000, 1_000, b"A".to_vec(), b"G".to_vec());
        let unmoved = TranscriptConsequence::default();
        assert_eq!(
            field_value("SHIFT_LENGTH", &snv, Some(&unmoved), &options, ","),
            "0"
        );
        let moved = TranscriptConsequence {
            shift_length: 7,
            ..Default::default()
        };
        assert_eq!(
            field_value("SHIFT_LENGTH", &snv, Some(&moved), &options, ","),
            "7"
        );
        assert_eq!(field_value("SHIFT_LENGTH", &snv, None, &options, ","), "");
        let mut deletion =
            InputVariant::new("21".into(), 1_000, 2_000, b"N".to_vec(), b"<DEL>".to_vec());
        deletion.is_structural = true;
        assert_eq!(
            field_value("SHIFT_LENGTH", &deletion, Some(&unmoved), &options, ","),
            ""
        );
    }

    /// Without `--shift_3prime 1` the SHIFT_LENGTH column is listed under
    /// `--shift_length` but no row sets it (the flag field follows the parameter,
    /// OutputFactory.pm:2516-2528; the value needs both, :1695), and without
    /// `--shift_length` the field is not listed at all.
    #[test]
    fn shift_length_is_absent_without_shift_3prime() {
        let listed = FieldOptions {
            shift_length: true,
            ..Default::default()
        };
        assert!(flag_fields(&listed).contains(&"SHIFT_LENGTH"));
        let snv = InputVariant::new("21".into(), 1_000, 1_000, b"A".to_vec(), b"G".to_vec());
        let moved = TranscriptConsequence {
            shift_length: 7,
            ..Default::default()
        };
        assert_eq!(
            field_value("SHIFT_LENGTH", &snv, Some(&moved), &listed, ","),
            ""
        );
        let unlisted = FieldOptions {
            shift_3prime: true,
            ..Default::default()
        };
        assert!(!flag_fields(&unlisted).contains(&"SHIFT_LENGTH"));
    }

    /// VEP computes a structural row's overlap from the allele's own start and end:
    /// a span inside the transcript carries OverlapBP and OverlapPC, an insertion pair
    /// (`start == end + 1`, the record whose END is its position) has length 0 and
    /// carries neither.
    #[test]
    fn insertion_pair_carries_no_overlap_fields() {
        let tc = TranscriptConsequence {
            transcript_id: std::sync::Arc::from("ENST1"),
            feature_start: 1_000,
            feature_end: 2_000,
            ..Default::default()
        };
        let mut span =
            InputVariant::new("21".into(), 1_500, 1_600, b"N".to_vec(), b"<INS>".to_vec());
        span.is_structural = true;
        span.variant_class = VariantClass::StructuralInsertion;
        assert_eq!(
            sv_overlap(&span, Some(&tc)),
            Some((101, 101.0 * 100.0 / 1_001.0))
        );
        assert!(extra_keys(&span, Some(&tc), &FieldOptions::default()).contains(&"OverlapBP"));

        let mut pair =
            InputVariant::new("21".into(), 1_501, 1_500, b"N".to_vec(), b"<INS>".to_vec());
        pair.is_structural = true;
        pair.variant_class = VariantClass::StructuralInsertion;
        assert_eq!(sv_overlap(&pair, Some(&tc)), None);
        assert!(!extra_keys(&pair, Some(&tc), &FieldOptions::default()).contains(&"OverlapBP"));
        assert_eq!(pair.location(), "21:1500-1501");
    }

    /// `NEAREST` is a record-level key: its value joins the list with the
    /// format's separator on a transcript row, an intergenic row and the
    /// `--summary` row alike, and the empty list is a written value
    /// (`OutputFactory.pm:1017` sets the key whenever the list exists), unlike an
    /// absent key, which prints nothing.
    #[test]
    fn nearest_prints_on_every_row_and_its_empty_list_is_a_written_value() {
        let options = FieldOptions {
            nearest: true,
            ..Default::default()
        };
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(field_value("NEAREST", &v, None, &options, ","), "");
        assert!(
            !empty_value_is_written("NEAREST", &v, None, &options),
            "no list: an absent key"
        );
        v.nearest = Some(vec!["GENE1".into(), "GENE2".into()]);
        assert_eq!(
            field_value("NEAREST", &v, None, &options, ","),
            "GENE1,GENE2"
        );
        let tc = TranscriptConsequence::default();
        assert_eq!(
            field_value("NEAREST", &v, Some(&tc), &options, ","),
            "GENE1,GENE2",
            "a transcript row carries the record-level key"
        );
        assert_eq!(
            field_value("NEAREST", &v, None, &options, "&"),
            "GENE1&GENE2"
        );
        assert!(!empty_value_is_written("NEAREST", &v, None, &options));
        v.summary_terms = Some(smallvec::smallvec![Consequence::MissenseVariant]);
        assert_eq!(
            field_value("NEAREST", &v, None, &options, ","),
            "GENE1,GENE2",
            "the summary row carries the record-level key"
        );
        assert_eq!(field_value("IMPACT", &v, None, &options, ","), "");
        v.summary_terms = None;
        v.nearest = Some(Vec::new());
        assert_eq!(field_value("NEAREST", &v, None, &options, ","), "");
        assert!(empty_value_is_written("NEAREST", &v, None, &options));
        assert!(!empty_value_is_written("VARIANT_CLASS", &v, None, &options));
        assert!(flag_fields(&options).contains(&"NEAREST"));
        assert!(!flag_fields(&FieldOptions::default()).contains(&"NEAREST"));
    }

    /// A `--summary` row is the record's hash alone (OutputFactory.pm:426-458
    /// over `VariationFeature_to_output_hash`): the tab format prints its identity,
    /// location and consequence list, `-` in every allele and feature column, no
    /// IMPACT, STRAND, ALLELE_NUM or Existing_variation, and the variant-level
    /// VARIANT_CLASS and MINIMISED beside them.
    #[test]
    fn summary_row_tab_cells_are_dashes() {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.id = Some("rs1".into());
        v.existing_variation = vec!["rs1".into()];
        v.minimised = true;
        v.most_severe_consequence = Some(Consequence::MissenseVariant);
        v.summary_terms = Some(smallvec::smallvec![
            Consequence::MissenseVariant,
            Consequence::IntronVariant
        ]);
        let options = FieldOptions {
            allele_number: true,
            variant_class: true,
            minimal: true,
            check_existing: true,
            symbol: true,
            ..Default::default()
        };
        let rows = record_rows(&[&v]);
        assert_eq!(rows.len(), 1);
        assert!(rows[0].1.is_none());
        let columns = tab_columns(&options, &[], None);
        let cells: Vec<String> = columns
            .iter()
            .map(|col| {
                let value = field_value(col, &v, None, &options, ",");
                if value.is_empty() {
                    "-".to_string()
                } else {
                    value
                }
            })
            .collect();
        assert_eq!(
            cells,
            [
                "rs1",
                "21:100",
                "-",
                "-",
                "-",
                "-",
                "missense_variant,intron_variant",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "SNV",
                "1",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
            ]
        );
        assert_eq!(
            columns[13..],
            [
                "ALLELE_NUM",
                "IMPACT",
                "DISTANCE",
                "STRAND",
                "FLAGS",
                "VARIANT_CLASS",
                "MINIMISED",
                "SYMBOL",
                "SYMBOL_SOURCE",
                "HGNC_ID",
                "CLIN_SIG",
                "SOMATIC",
                "PHENO"
            ]
        );
        assert_eq!(
            field_value("Consequence", &v, None, &options, "&"),
            "missense_variant&intron_variant"
        );
        v.plugin_data
            .insert("CADD_PHRED".to_string(), "9".to_string());
        assert_eq!(field_value("CADD_PHRED", &v, None, &options, ","), "");
    }

    /// An intergenic row flagged by `--flag_pick` prints `PICK=1` as a transcript
    /// row does (OutputFactory.pm:1363 sets PICK on every overlap allele).
    #[test]
    fn intergenic_row_prints_pick_when_flagged() {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let options = FieldOptions {
            include_pick: true,
            ..Default::default()
        };
        assert_eq!(field_value("PICK", &v, None, &options, ","), "");
        v.intergenic_pick = true;
        assert_eq!(field_value("PICK", &v, None, &options, ","), "1");
        assert_eq!(
            field_value("PICK", &v, None, &FieldOptions::default(), ","),
            "",
            "PICK is printed only under a `--flag_pick*` option"
        );
    }

    #[test]
    fn default_options_emit_only_the_user_group() {
        assert_eq!(
            flag_fields(&FieldOptions::default()),
            vec!["IMPACT", "DISTANCE", "STRAND", "FLAGS"]
        );
        assert!(!FieldOptions::default().any_optional_group());
    }

    #[test]
    fn record_rows_are_transcript_major_allele_minor_intergenic_last() {
        let mut a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        a.allele_index = 0;
        let mut b = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"T".to_vec());
        b.allele_index = 1;
        let tc = |id: &str| TranscriptConsequence {
            transcript_id: std::sync::Arc::from(id),
            ..Default::default()
        };
        a.transcript_consequences = vec![tc("ENST2"), tc("ENST1")];
        b.transcript_consequences = vec![tc("ENST1"), tc("ENST2")];
        let mut c = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"C".to_vec());
        c.allele_index = 2;
        c.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let rows = record_rows(&[&a, &b, &c]);
        let order: Vec<(String, usize)> = rows
            .iter()
            .map(|(v, t)| {
                (
                    t.map(|t| t.transcript_id.to_string()).unwrap_or_default(),
                    v.allele_index,
                )
            })
            .collect();
        assert_eq!(
            order,
            vec![
                ("ENST1".to_string(), 0),
                ("ENST1".to_string(), 1),
                ("ENST2".to_string(), 0),
                ("ENST2".to_string(), 1),
                (String::new(), 2),
            ]
        );
    }

    /// The reference row of `--process_ref_homs` leads each transcript's rows
    /// (BaseVariationFeatureOverlap.pm:363-373 lists the reference allele first),
    /// then the alternates in allele order.
    #[test]
    fn reference_row_leads_each_transcripts_rows() {
        let tc = |id: &str| TranscriptConsequence {
            transcript_id: std::sync::Arc::from(id),
            ..Default::default()
        };
        let mut alt1 = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        alt1.transcript_consequences = vec![tc("ENST2"), tc("ENST1")];
        let mut alt2 = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"T".to_vec());
        alt2.allele_index = 1;
        alt2.transcript_consequences = vec![tc("ENST1")];
        let mut reference = alt1.reference_allele_row();
        reference.transcript_consequences = vec![tc("ENST1"), tc("ENST2")];
        let order: Vec<(String, &[u8])> = record_rows(&[&reference, &alt1, &alt2])
            .iter()
            .map(|(v, tc)| (tc.unwrap().transcript_id.to_string(), v.alt_allele()))
            .collect();
        assert_eq!(
            order,
            [
                ("ENST1".to_string(), b"A".as_slice()),
                ("ENST1".to_string(), b"G"),
                ("ENST1".to_string(), b"T"),
                ("ENST2".to_string(), b"A"),
                ("ENST2".to_string(), b"G"),
            ]
        );
    }

    /// The reusable form clears what an earlier record left behind and orders
    /// rows exactly as `record_rows` does.
    #[test]
    fn record_rows_into_reuses_the_vector_and_matches_record_rows() {
        let tc = |id: &str| TranscriptConsequence {
            transcript_id: std::sync::Arc::from(id),
            ..Default::default()
        };
        let mut a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        a.transcript_consequences = vec![tc("ENST9"), tc("ENST1")];
        let mut b = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"T".to_vec());
        b.allele_index = 1;
        b.transcript_consequences = vec![tc("ENST1")];
        let mut c = InputVariant::new("21".into(), 300, 300, b"C".to_vec(), b"T".to_vec());
        c.most_severe_consequence = Some(Consequence::IntergenicVariant);

        let mut rows = Vec::new();
        record_rows_into(&[&c], &mut rows);
        assert_eq!(rows.len(), 1);
        record_rows_into(&[&a, &b], &mut rows);
        let key = |rows: &[Row<'_>]| -> Vec<(String, usize)> {
            rows.iter()
                .map(|(v, t)| {
                    (
                        t.map(|t| t.transcript_id.to_string()).unwrap_or_default(),
                        v.allele_index,
                    )
                })
                .collect()
        };
        assert_eq!(key(&rows), key(&record_rows(&[&a, &b])));
        assert_eq!(
            key(&rows),
            vec![
                ("ENST1".to_string(), 0),
                ("ENST1".to_string(), 1),
                ("ENST9".to_string(), 0),
            ]
        );
    }

    #[test]
    fn extra_fields_plan_carries_the_options_flag_fields() {
        let options = FieldOptions {
            symbol: true,
            hgvs: true,
            ..Default::default()
        };
        let plan = ExtraFieldsPlan::new(options.clone());
        assert_eq!(plan.flag_fields, flag_fields(&options));
        assert!(plan.any_optional_group);
        assert!(!ExtraFieldsPlan::new(FieldOptions::default()).any_optional_group);
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(
            extra_keys_with(&plan.flag_fields, &v, None, &options),
            extra_keys(&v, None, &options)
        );
    }

    /// `GENCODE_PRIMARY` is the literal `1` on a transcript carrying the
    /// `gencode_primary` attribute and empty otherwise (OutputFactory.pm:1604-1606
    /// at release/116.2: `$hash->{GENCODE_PRIMARY} = 1` on the attribute's
    /// presence); the column exists under `--flag_gencode_primary` alone, and its
    /// header text is Ensembl VEP's (Constants.pm:265).
    #[test]
    fn gencode_primary_column_is_1_on_presence() {
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let primary = TranscriptConsequence {
            gencode_primary: true,
            ..Default::default()
        };
        let other = TranscriptConsequence::default();
        let options = FieldOptions {
            gencode_primary: true,
            ..Default::default()
        };
        assert_eq!(
            field_value("GENCODE_PRIMARY", &v, Some(&primary), &options, ","),
            "1"
        );
        assert_eq!(
            field_value("GENCODE_PRIMARY", &v, Some(&other), &options, ","),
            ""
        );
        assert_eq!(field_value("GENCODE_PRIMARY", &v, None, &options, ","), "");
        assert!(extra_keys(&v, Some(&primary), &options).contains(&"GENCODE_PRIMARY"));
        assert!(
            !extra_keys(&v, Some(&primary), &FieldOptions::default()).contains(&"GENCODE_PRIMARY")
        );
        assert_eq!(
            field_description("GENCODE_PRIMARY"),
            Some("Reports if transcript is GENCODE primary")
        );
    }

    /// A transcript row of the attributes corpus: the NDUFV3 transcript
    /// ENST00000354250 as Ensembl VEP 116.2 reports it for `21 42893293 C T`
    /// under the attribute flags (`MANE=MANE_Select;MANE_SELECT=NM_021075.4;
    /// TSL=1;APPRIS=PM;CCDS=CCDS33572.1;...;SWISSPROT=P56181.185;
    /// UNIPARC=UPI000037848B;UNIPROT_ISOFORM=P56181-2;RefSeq=NM_021075.4`),
    /// with the cache's raw attribute values on the consequence.
    fn ndufv3_consequence() -> TranscriptConsequence {
        TranscriptConsequence {
            transcript_id: std::sync::Arc::from("ENST00000354250"),
            gene_id: std::sync::Arc::from("ENSG00000160194"),
            transcript_version: Some(7),
            gene_version: Some(20),
            gene_symbol: Some("NDUFV3".into()),
            biotype: Some("protein_coding".into()),
            canonical: true,
            mane_select: Some("NM_021075.4".to_string()),
            tsl: Some(1),
            appris: Some("principalM".to_string()),
            ccds: Some("CCDS33572.1".to_string()),
            protein_id: Some("ENSP00000346196".into()),
            swissprot: Some("P56181.185".to_string()),
            trembl: Some("-".to_string()),
            uniparc: Some("UPI000037848B".to_string()),
            uniprot_isoform: Some("P56181-2".to_string()),
            refseq: Some("NM_021075.4".to_string()),
            gene_phenotype: false,
            ..Default::default()
        }
    }

    fn attribute_options() -> FieldOptions {
        FieldOptions {
            symbol: true,
            biotype: true,
            canonical: true,
            mane: true,
            mane_select: true,
            mane_plus_clinical: true,
            tsl: true,
            appris: true,
            ccds: true,
            protein: true,
            uniprot: true,
            xref_refseq: true,
            gene_phenotype: true,
            mirna: true,
            transcript_version: true,
            gene_version: true,
            ..Default::default()
        }
    }

    fn value(field: &str, tc: &TranscriptConsequence, options: &FieldOptions) -> String {
        let v = InputVariant::new(
            "21".into(),
            42_893_293,
            42_893_293,
            b"C".to_vec(),
            b"T".to_vec(),
        );
        field_value(field, &v, Some(tc), options, ",")
    }

    /// `s/principal/P/` and `s/alternative/A/` (OutputFactory.pm:1616-1622 at
    /// release/116.2) abbreviate the grade and keep the rest, so `principalM` is
    /// `PM` (ENST00000354250 above) and `alternativeM` is `AM` (its neighbour
    /// ENST00000340344 on the same record: `APPRIS=AM`); `principal1` is `P1`,
    /// `alternative2` is `A2`, and an empty value sets nothing.
    #[test]
    fn appris_prints_p_and_a_prefixes() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("APPRIS", &tc, &options), "PM");
        for (raw, printed) in [
            ("alternativeM", "AM"),
            ("principal1", "P1"),
            ("alternative2", "A2"),
            ("principal4", "P4"),
        ] {
            tc.appris = Some(raw.to_string());
            assert_eq!(value("APPRIS", &tc, &options), printed, "{raw}");
        }
        tc.appris = Some(String::new());
        assert_eq!(value("APPRIS", &tc, &options), "");
        tc.appris = None;
        assert_eq!(value("APPRIS", &tc, &options), "");
    }

    /// `CANONICAL` is the literal `YES` on the canonical transcript and absent
    /// otherwise (OutputFactory.pm:1577; `CANONICAL=YES` on ENST00000354250 and
    /// no key on ENST00000340344).
    #[test]
    fn canonical_prints_yes() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("CANONICAL", &tc, &options), "YES");
        tc.canonical = false;
        assert_eq!(value("CANONICAL", &tc, &options), "");
    }

    /// `BIOTYPE` is the transcript's biotype under `--biotype`
    /// (OutputFactory.pm:1580; `BIOTYPE=protein_coding` above) and is not a
    /// column without the flag.
    #[test]
    fn biotype_column_under_flag() {
        let tc = ndufv3_consequence();
        assert_eq!(
            value("BIOTYPE", &tc, &attribute_options()),
            "protein_coding"
        );
        assert!(flag_fields(&attribute_options()).contains(&"BIOTYPE"));
        assert!(!flag_fields(&FieldOptions::default()).contains(&"BIOTYPE"));
    }

    /// `BIOTYPE` leads the `regulatory` flag group (Constants.pm:136), so
    /// `--tab --regulatory` without `--biotype` prints 23 columns, BIOTYPE the
    /// 18th; its value is set only under `--biotype` (OutputFactory.pm:1580), so
    /// a transcript row prints `-` there and the default format's Extra carries
    /// no `BIOTYPE=`. Ensembl VEP 116.2 prints that header and `-` on every
    /// transcript row of a record at 21:10541165. Under both flags the column
    /// is listed once, at the `biotype` group's place, and filled.
    #[test]
    fn regulatory_group_lists_biotype_without_its_flag() {
        let regulatory = FieldOptions {
            regulatory: true,
            ..Default::default()
        };
        let columns = tab_columns(&regulatory, &[], None);
        assert_eq!(columns.len(), 23);
        assert_eq!(columns[17], "BIOTYPE");
        assert_eq!(
            &columns[17..],
            [
                "BIOTYPE",
                "MOTIF_NAME",
                "MOTIF_POS",
                "HIGH_INF_POS",
                "MOTIF_SCORE_CHANGE",
                "TRANSCRIPTION_FACTORS",
            ]
        );
        let tc = ndufv3_consequence();
        assert_eq!(value("BIOTYPE", &tc, &regulatory), "");
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert!(!empty_value_is_written(
            "BIOTYPE",
            &v,
            Some(&tc),
            &regulatory
        ));
        v.transcript_consequences = vec![tc];
        let tab = crate::output::tab::TabOutputFormatter::new(
            regulatory.clone(),
            Vec::new(),
            None,
            Vec::new(),
            true,
        );
        let row = crate::OutputFormatter::format_record(&tab, &[&v]).unwrap();
        let cells: Vec<&str> = row[0].split('\t').collect();
        assert_eq!(cells.len(), 23);
        assert_eq!(cells[17], "-");

        let both = FieldOptions {
            biotype: true,
            ..regulatory
        };
        let columns = tab_columns(&both, &[], None);
        assert_eq!(columns.len(), 23);
        assert_eq!(columns[17], "BIOTYPE");
        assert_eq!(columns.iter().filter(|c| *c == "BIOTYPE").count(), 1);
        let tab =
            crate::output::tab::TabOutputFormatter::new(both, Vec::new(), None, Vec::new(), true);
        let row = crate::OutputFormatter::format_record(&tab, &[&v]).unwrap();
        assert_eq!(row[0].split('\t').nth(17), Some("protein_coding"));
    }

    /// The cache writes `-` where a transcript has no CCDS identifier; Ensembl
    /// VEP sets `CCDS` only for a value other than `-` (OutputFactory.pm:1534-1537),
    /// so the default format's Extra column carries no `CCDS=-`.
    #[test]
    fn ccds_dash_sentinel_suppressed() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("CCDS", &tc, &options), "CCDS33572.1");
        tc.ccds = Some("-".to_string());
        assert_eq!(value("CCDS", &tc, &options), "");
        tc.ccds = None;
        assert_eq!(value("CCDS", &tc, &options), "");
    }

    /// `SWISSPROT`, `TREMBL`, `UNIPARC` and `UNIPROT_ISOFORM` each print nothing
    /// for the cache's `-` (OutputFactory.pm:1568-1574 reads `-` as undefined):
    /// the row above carries `SWISSPROT`, `UNIPARC` and `UNIPROT_ISOFORM` and no
    /// `TREMBL`, whose cache value is `-`.
    #[test]
    fn uniprot_dash_sentinel_suppressed() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("TREMBL", &tc, &options), "");
        assert_eq!(value("SWISSPROT", &tc, &options), "P56181.185");
        tc.swissprot = Some("-".to_string());
        tc.uniparc = Some("-".to_string());
        tc.uniprot_isoform = Some("-".to_string());
        for field in ["SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM"] {
            assert_eq!(value(field, &tc, &options), "", "{field}");
        }
    }

    /// `UNIPARC` and `UNIPROT_ISOFORM` carry the cache's values
    /// (`UNIPARC=UPI000037848B;UNIPROT_ISOFORM=P56181-2` above); a
    /// comma-separated accession list is a list in every format, joined by the
    /// format's separator (`TREMBL=L7XE61_HUMAN,L7XCZ9_HUMAN,B7Z313_HUMAN` on the
    /// GRCh37 APP transcript ENST00000346798 for `21 27253471 A G`; `&` inside
    /// a VCF entry).
    #[test]
    fn uniparc_and_isoform_columns_populated() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("UNIPARC", &tc, &options), "UPI000037848B");
        assert_eq!(value("UNIPROT_ISOFORM", &tc, &options), "P56181-2");
        tc.trembl = Some("L7XE61_HUMAN,L7XCZ9_HUMAN,B7Z313_HUMAN".to_string());
        assert_eq!(
            value("TREMBL", &tc, &options),
            "L7XE61_HUMAN,L7XCZ9_HUMAN,B7Z313_HUMAN"
        );
        let v = InputVariant::new(
            "21".into(),
            27_253_471,
            27_253_471,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        assert_eq!(
            field_value("TREMBL", &v, Some(&tc), &options, "&"),
            "L7XE61_HUMAN&L7XCZ9_HUMAN&B7Z313_HUMAN"
        );
        tc.trembl = Some("A,,B,".to_string());
        assert_eq!(value("TREMBL", &tc, &options), "A,,B");
        let fields = flag_fields(&options);
        let at = |name: &str| fields.iter().position(|f| *f == name).unwrap();
        assert!(at("SWISSPROT") < at("TREMBL") && at("TREMBL") < at("UNIPARC"));
        assert!(at("UNIPARC") < at("UNIPROT_ISOFORM"));
    }

    /// `RefSeq` is the comma-separated list of RefSeq transcript ids split into
    /// a list (`[split(',', $tr->{_refseq})]`, OutputFactory.pm:1540-1543;
    /// `RefSeq=NM_198996.4,NM_001302998.2` on the LIPI transcript
    /// ENST00000681601 for `21 14108946 C T`) and absent for the cache's `-`.
    #[test]
    fn refseq_dash_sentinel_suppressed() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        tc.refseq = Some("NM_198996.4,NM_001302998.2".to_string());
        assert_eq!(value("RefSeq", &tc, &options), "NM_198996.4,NM_001302998.2");
        tc.refseq = Some("-".to_string());
        assert_eq!(value("RefSeq", &tc, &options), "");
        tc.refseq = None;
        assert_eq!(value("RefSeq", &tc, &options), "");
    }

    /// `GENE_PHENO` is the literal `1` when the transcript's gene has a
    /// phenotype association count above zero and absent otherwise
    /// (OutputFactory.pm:1586 reads `$tr->{_gene_phenotype}` for truth;
    /// `GENE_PHENO=1` on the COL6A2 transcript ENST00000397763 for
    /// `21 46098139 T C`, no key on NDUFV3's ENST00000354250 above); the
    /// co-located variants play no part.
    #[test]
    fn gene_pheno_from_transcript_count() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("GENE_PHENO", &tc, &options), "");
        tc.gene_phenotype = true;
        assert_eq!(value("GENE_PHENO", &tc, &options), "1");
        let mut v = InputVariant::new("21".into(), 1, 1, b"A".to_vec(), b"G".to_vec());
        v.colocated_variants
            .push(vep_core::variant::ColocatedVariant {
                id: "rs1".to_string(),
                phenotype: true,
                ..Default::default()
            });
        tc.gene_phenotype = false;
        assert_eq!(field_value("GENE_PHENO", &v, Some(&tc), &options, ","), "");
    }

    /// `DOMAINS` joins each feature's label and identifier with `:` and, as
    /// `join(':', grep {$_} ...)` does (OutputFactory.pm:1500), prints the
    /// identifier alone for a feature whose cache entry carries no label or
    /// the label `0`, which Perl reads as false.
    #[test]
    fn domains_without_a_label_print_the_identifier_alone() {
        let options = FieldOptions {
            domains: true,
            ..Default::default()
        };
        let mut tc = ndufv3_consequence();
        tc.domains = vec![
            ("Pfam".to_string(), "PF00080".to_string()),
            (String::new(), "CYTOPLASMIC_DOMAIN".to_string()),
            ("Gene3D".to_string(), "3.30.160.20".to_string()),
            ("0".to_string(), "PF00080".to_string()),
        ];
        assert_eq!(
            value("DOMAINS", &tc, &options),
            "Pfam:PF00080,CYTOPLASMIC_DOMAIN,Gene3D:3.30.160.20,PF00080"
        );
        tc.domains.clear();
        assert_eq!(value("DOMAINS", &tc, &options), "");
    }

    /// `--gene_version` appends `.<version>` to the Gene column
    /// (OutputFactory.pm:1469: `ENSG00000160194.20` above); nothing is appended
    /// without the flag, without a version, or to an id already ending in
    /// `.<digits>`.
    #[test]
    fn gene_version_suffix_on_gene_column() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("Gene", &tc, &options), "ENSG00000160194.20");
        let off = FieldOptions {
            gene_version: false,
            ..attribute_options()
        };
        assert_eq!(value("Gene", &tc, &off), "ENSG00000160194");
        tc.gene_version = None;
        assert_eq!(value("Gene", &tc, &options), "ENSG00000160194");
        tc.gene_version = Some(0);
        assert_eq!(value("Gene", &tc, &options), "ENSG00000160194");
        tc.gene_id = std::sync::Arc::from("ENSG00000160194.3");
        tc.gene_version = Some(20);
        assert_eq!(value("Gene", &tc, &options), "ENSG00000160194.3");
        assert_eq!(version_suffix("NM_001.2", Some(5)), None);
        assert_eq!(version_suffix("NM_001.a", Some(5)), Some(5));
        assert_eq!(version_suffix("", Some(5)), None);
    }

    /// `--transcript_version` appends `.<version>` to the Feature column
    /// (OutputFactory.pm:1465: `ENST00000354250.7` above) under the same rule.
    #[test]
    fn transcript_version_suffix_on_feature() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("Feature", &tc, &options), "ENST00000354250.7");
        let off = FieldOptions {
            transcript_version: false,
            ..attribute_options()
        };
        assert_eq!(value("Feature", &tc, &off), "ENST00000354250");
        tc.transcript_id = std::sync::Arc::from("NM_021075.4");
        assert_eq!(value("Feature", &tc, &options), "NM_021075.4");
    }

    /// Under `--mane`, `MANE` lists `MANE_Select` for a MANE Select transcript
    /// and `MANE_Plus_Clinical` for a MANE Plus Clinical one, in that order when
    /// a transcript carries both, and `MANE_SELECT` and `MANE_PLUS_CLINICAL`
    /// carry the attribute values (OutputFactory.pm:1587-1601;
    /// `MANE=MANE_Select;MANE_SELECT=NM_021075.4` above and
    /// `MANE=MANE_Plus_Clinical;MANE_PLUS_CLINICAL=NM_058174.3` on COL6A2's
    /// ENST00000397763 for `21 46098139 T C`). Inside a VCF entry the list
    /// joins with `&`.
    #[test]
    fn mane_lists_select_then_plus_clinical() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("MANE", &tc, &options), "MANE_Select");
        assert_eq!(value("MANE_SELECT", &tc, &options), "NM_021075.4");
        assert_eq!(value("MANE_PLUS_CLINICAL", &tc, &options), "");
        tc.mane_select = None;
        tc.mane_plus_clinical = Some("NM_058174.3".to_string());
        assert_eq!(value("MANE", &tc, &options), "MANE_Plus_Clinical");
        assert_eq!(value("MANE_SELECT", &tc, &options), "");
        assert_eq!(value("MANE_PLUS_CLINICAL", &tc, &options), "NM_058174.3");
        tc.mane_select = Some("NM_021075.4".to_string());
        assert_eq!(
            value("MANE", &tc, &options),
            "MANE_Select,MANE_Plus_Clinical"
        );
        let v = InputVariant::new("21".into(), 1, 1, b"A".to_vec(), b"G".to_vec());
        assert_eq!(
            field_value("MANE", &v, Some(&tc), &options, "&"),
            "MANE_Select&MANE_Plus_Clinical"
        );
        tc.mane_select = None;
        tc.mane_plus_clinical = None;
        assert_eq!(value("MANE", &tc, &options), "");
    }

    /// `MANE_SELECT` is the `MANE_Select` attribute's value as the cache spells
    /// it, a RefSeq accession with its version (`NM_021075.4` above;
    /// `NR_185847.1` is the attribute a non-coding MANE Select transcript of the
    /// GRCh38 cache carries), never the transcript's own name; an attribute
    /// present with an empty or `0` value still lists the set but sets no
    /// accession (OutputFactory.pm:1588-1592).
    #[test]
    fn mane_select_value_is_attribute_verbatim() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("MANE_SELECT", &tc, &options), "NM_021075.4");
        tc.mane_select = Some("NR_185847.1".to_string());
        assert_eq!(value("MANE_SELECT", &tc, &options), "NR_185847.1");
        tc.mane_select = Some(String::new());
        assert_eq!(value("MANE_SELECT", &tc, &options), "");
        assert_eq!(value("MANE", &tc, &options), "MANE_Select");
        tc.mane_select = Some("0".to_string());
        assert_eq!(value("MANE_SELECT", &tc, &options), "");
        assert_eq!(value("MANE", &tc, &options), "MANE_Select");
    }

    /// `--mane_select` alone sets Ensembl VEP's `mane_select` and not
    /// `mane_plus_clinical` (Config.pm:580-586 ties the pair to `--mane`), so a
    /// MANE Plus Clinical transcript gets no `MANE` list entry and the
    /// `MANE_PLUS_CLINICAL` column is not produced (Constants.pm:88).
    #[test]
    fn mane_select_alone_omits_plus_clinical() {
        let options = FieldOptions {
            mane_select: true,
            ..Default::default()
        };
        let mut tc = ndufv3_consequence();
        tc.mane_plus_clinical = Some("NM_058174.3".to_string());
        assert_eq!(value("MANE", &tc, &options), "MANE_Select");
        tc.mane_select = None;
        assert_eq!(value("MANE", &tc, &options), "");
        assert_eq!(
            flag_fields(&options),
            vec![
                "IMPACT",
                "DISTANCE",
                "STRAND",
                "FLAGS",
                "MANE",
                "MANE_SELECT"
            ]
        );
        assert!(!field_is_produced("MANE_PLUS_CLINICAL", &options));
        let both = FieldOptions {
            mane: true,
            mane_select: true,
            mane_plus_clinical: true,
            ..Default::default()
        };
        assert_eq!(value("MANE", &tc, &both), "MANE_Plus_Clinical");
        assert!(field_is_produced("MANE_PLUS_CLINICAL", &both));
    }

    /// `TSL` prints the digits alone (`TSL=1` for `tsl1`, and for
    /// `tsl1 (assigned to previous version 4)` on the TPTE transcript
    /// ENST00000618007; OutputFactory.pm:1609-1613 keeps `$1` of `m/tsl(\d+)/`
    /// when it is not zero) and nothing for `tslNA`, which the loader leaves
    /// `None`.
    #[test]
    fn tsl_prints_digit_only() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("TSL", &tc, &options), "1");
        tc.tsl = Some(5);
        assert_eq!(value("TSL", &tc, &options), "5");
        tc.tsl = None;
        assert_eq!(value("TSL", &tc, &options), "");
    }

    /// The `miRNA` column joins the structure terms the consequence carries with
    /// the format's separator (`miRNA=miRNA_loop,miRNA_stem` on the MIR6724-2
    /// transcript ENST00000616483 for the deletion `21 8249508 TGCG T`; `&`
    /// inside a VCF entry) and is empty when the consequence carries none.
    #[test]
    fn mirna_column_joins_the_terms() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        assert_eq!(value("miRNA", &tc, &options), "");
        tc.mirna = vec!["miRNA_loop".to_string(), "miRNA_stem".to_string()];
        assert_eq!(value("miRNA", &tc, &options), "miRNA_loop,miRNA_stem");
        let v = InputVariant::new("21".into(), 1, 1, b"A".to_vec(), b"G".to_vec());
        assert_eq!(
            field_value("miRNA", &v, Some(&tc), &options, "&"),
            "miRNA_loop&miRNA_stem"
        );
        assert!(flag_fields(&options).contains(&"miRNA"));
        assert!(!field_is_produced("miRNA", &FieldOptions::default()));
    }

    /// `SYMBOL`, `SYMBOL_SOURCE`, `HGNC_ID` and `ENSP` print nothing for the
    /// cache's `-` (OutputFactory.pm:1525-1530, 1556-1559).
    #[test]
    fn identifier_columns_read_dash_as_absent() {
        let options = attribute_options();
        let mut tc = ndufv3_consequence();
        tc.gene_symbol = Some("-".into());
        tc.gene_symbol_source = Some("-".into());
        tc.hgnc_id = Some("-".into());
        tc.protein_id = Some("-".into());
        for field in ["SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "ENSP"] {
            assert_eq!(value(field, &tc, &options), "", "{field}");
        }
    }

    #[test]
    fn flags_keep_cds_codes_only() {
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            flags: std::sync::Arc::from(["cds_start_NF".to_string(), "gencode_basic".to_string()]),
            ..Default::default()
        };
        assert_eq!(
            field_value("FLAGS", &v, Some(&tc), &FieldOptions::default(), ","),
            "cds_start_NF"
        );
    }

    #[test]
    fn header_ends_with_the_column_line() {
        let lines = default_format_header(&FieldOptions::default(), &[], "/c", None, "vep -i x");
        assert!(lines[0].starts_with("## ENSEMBL VARIANT EFFECT PREDICTOR v"));
        assert_eq!(lines[1], format!("## vep-rs {VEP_RS_VERSION}"));
        assert!(lines.contains(&"## Column descriptions:".to_string()));
        assert!(lines.contains(&"## Extra column keys:".to_string()));
        assert!(lines.contains(&"## FLAGS : Transcript quality flags".to_string()));
        assert_eq!(
            lines.last().unwrap(),
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra"
        );
    }

    fn missense_tc() -> TranscriptConsequence {
        TranscriptConsequence {
            transcript_id: std::sync::Arc::from("ENST1"),
            consequences: smallvec::smallvec![Consequence::MissenseVariant],
            hgvsp: Some("ENSP1:p.Arg97=".to_string()),
            ..Default::default()
        }
    }

    /// `OutputFactory.pm` `TranscriptVariationAllele_to_output_hash` turns the `=`
    /// of an HGVSp notation into `%3D` unless `no_escape` is set.
    #[test]
    fn no_escape_leaves_hgvsp_equals_raw() {
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = missense_tc();
        let escaped = FieldOptions {
            hgvs: true,
            ..Default::default()
        };
        assert_eq!(
            field_value("HGVSp", &v, Some(&tc), &escaped, ","),
            "ENSP1:p.Arg97%3D"
        );
        let raw = FieldOptions {
            hgvs: true,
            no_escape: true,
            ..Default::default()
        };
        assert_eq!(
            field_value("HGVSp", &v, Some(&tc), &raw, ","),
            "ENSP1:p.Arg97="
        );
    }

    /// The `AMBIGUITY` group follows `--ambiguity` and sits after `GENE_PHENO`
    /// in `@FLAG_FIELDS`; its value is the record's IUPAC code, so an indel, whose
    /// allele set has none, carries no key (`OutputFactory.pm`
    /// `VariationFeature_to_output_hash`: the key is set only for a defined code).
    #[test]
    fn ambiguity_key_absent_for_indels() {
        let opts = FieldOptions {
            ambiguity: true,
            gene_phenotype: true,
            sift: true,
            ..Default::default()
        };
        let fields = flag_fields(&opts);
        let at = |name: &str| fields.iter().position(|f| *f == name).unwrap();
        assert!(at("GENE_PHENO") < at("AMBIGUITY") && at("AMBIGUITY") < at("SIFT"));
        assert!(!flag_fields(&FieldOptions::default()).contains(&"AMBIGUITY"));

        let snv = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(field_value("AMBIGUITY", &snv, None, &opts, ","), "R");
        let mut multi = InputVariant::new("21".into(), 100, 100, b"C".to_vec(), b"G".to_vec());
        multi.record_allele_string_multi = Some("C/G/T".into());
        assert_eq!(field_value("AMBIGUITY", &multi, None, &opts, ","), "B");
        let deletion = InputVariant::new("21".into(), 101, 102, b"AT".to_vec(), b"-".to_vec());
        assert_eq!(field_value("AMBIGUITY", &deletion, None, &opts, ","), "");
        assert!(!extra_keys(&deletion, None, &opts).is_empty());
        let tc = missense_tc();
        assert_eq!(
            field_value("AMBIGUITY", &deletion, Some(&tc), &opts, ","),
            ""
        );
    }

    /// `UPLOADED_ALLELE` is the record's allele string as uploaded, `REF/ALT1/ALT2`,
    /// the same on every row of the record (`OutputFactory.pm`
    /// `VariationFeatureOverlapAllele_to_output_hash`: `original_allele_string ||
    /// allele_string`), the group follows `--uploaded_allele`, and a structural row
    /// carries neither it nor `REF_ALLELE`
    /// (`BaseStructuralVariationOverlapAllele_to_output_hash`).
    #[test]
    fn uploaded_allele_prints_the_record_allele_string() {
        let opts = FieldOptions {
            uploaded_allele: true,
            ..Default::default()
        };
        assert_eq!(
            flag_fields(&opts),
            ["UPLOADED_ALLELE", "IMPACT", "DISTANCE", "STRAND", "FLAGS"]
        );
        assert!(!flag_fields(&FieldOptions::default()).contains(&"UPLOADED_ALLELE"));
        let mut a = InputVariant::new("21".into(), 101, 102, b"AG".to_vec(), b"A".to_vec());
        a.uploaded_allele_string = Some("TAG/TA/T".into());
        let mut b = InputVariant::new("21".into(), 101, 102, b"AG".to_vec(), b"-".to_vec());
        b.uploaded_allele_string = Some("TAG/TA/T".into());
        for v in [&a, &b] {
            assert_eq!(
                field_value("UPLOADED_ALLELE", v, None, &opts, ","),
                "TAG/TA/T"
            );
            assert_eq!(
                field_value("UPLOADED_ALLELE", v, Some(&missense_tc()), &opts, ","),
                "TAG/TA/T"
            );
        }
        let plain = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(
            field_value("UPLOADED_ALLELE", &plain, None, &opts, ","),
            "A/G"
        );
        let mut sv = InputVariant::new("21".into(), 100, 200, b"N".to_vec(), b"<DEL>".to_vec());
        sv.is_structural = true;
        sv.uploaded_allele_string = Some("<DEL>".into());
        let both = FieldOptions {
            uploaded_allele: true,
            show_ref_allele: true,
            ..Default::default()
        };
        assert_eq!(field_value("UPLOADED_ALLELE", &sv, None, &both, ","), "");
        assert_eq!(field_value("REF_ALLELE", &sv, None, &both, ","), "");
        assert_eq!(field_value("REF_ALLELE", &plain, None, &both, ","), "A");
    }

    /// Every row of a `--minimal` split prints the reference cached on the record's
    /// feature (OutputFactory.pm:1357 at :2289, VariationFeature.pm:1533-1539):
    /// the first allele's trimmed reference once its feature wrote a row, the
    /// record's own when it wrote none, whether the row filters found that
    /// (`split_first_feature_rowless`) or the allele equals the reference once
    /// trimmed, in which case `--process_ref_homs` builds the reference row that
    /// caches `-`. The byte writer and the string table agree on each.
    #[test]
    fn ref_allele_of_a_split_row_is_the_reference_its_first_feature_cached() {
        let options = FieldOptions {
            show_ref_allele: true,
            ..Default::default()
        };
        let reference_rows = FieldOptions {
            process_ref_homs: true,
            ..options.clone()
        };
        let mut deletion = InputVariant::new(
            "21".into(),
            10_542_394,
            10_542_397,
            b"TACG".to_vec(),
            b"-".to_vec(),
        );
        deletion.allele_index = 1;
        deletion.original_start = Some(10_542_393);
        deletion.original_end = Some(10_542_397);
        deletion.record_allele_string_multi = Some(Arc::from("GTACG/CTACG/G"));
        for (rowless, expected) in [(false, "G"), (true, "GTACG")] {
            deletion.split_first_feature_rowless = rowless;
            for options in [&options, &reference_rows] {
                assert_eq!(
                    field_value("REF_ALLELE", &deletion, Some(&missense_tc()), options, ","),
                    expected,
                    "rowless {rowless}"
                );
                assert_writer_matches(&deletion, Some(&missense_tc()), options, ",");
            }
        }
        let mut equal = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        equal.allele_index = 1;
        equal.original_start = Some(100);
        equal.original_end = Some(100);
        equal.record_allele_string_multi = Some(Arc::from("A/A/G"));
        assert_eq!(
            field_value("REF_ALLELE", &equal, Some(&missense_tc()), &options, ","),
            "A"
        );
        assert_eq!(
            field_value(
                "REF_ALLELE",
                &equal,
                Some(&missense_tc()),
                &reference_rows,
                ","
            ),
            "-"
        );
    }

    /// `REF_ALLELE`, like `ALLELE_NUM` and `UPLOADED_ALLELE`, is set per overlap
    /// allele (`OutputFactory.pm` `VariationFeatureOverlapAllele_to_output_hash`),
    /// so a row built from the variant alone, a `--most_severe` or `--summary`
    /// row, never carries it.
    #[test]
    fn ref_allele_absent_on_summary_rows() {
        let opts = FieldOptions {
            allele_number: true,
            show_ref_allele: true,
            uploaded_allele: true,
            symbol: true,
            ..Default::default()
        };
        assert_eq!(vf_level_flag_fields(&opts), Vec::<&str>::new());
        assert!(flag_fields(&opts).contains(&"REF_ALLELE"));
    }

    /// `--individual` adds `IND` and `ZYG` first among the flag fields
    /// (Constants.pm:69, the first `@FLAG_FIELDS` entry) with Constants.pm:193-194's
    /// descriptions; both are set from the variant alone
    /// (`VariationFeature_to_output_hash`, OutputFactory.pm:969-977), so a summary
    /// row keeps them; neither is listed without the flag, and both are empty on an
    /// allele that is no individual's copy. `ZYG` is `HET` for more than one
    /// distinct genotype allele, else `HOM`, with `REF` appended for a
    /// reference-homozygous copy (:976); the reference row of `--process_ref_homs`
    /// numbers 0 under `ALLELE_NUM` and prints the reference as its `Allele`.
    #[test]
    fn ind_and_zyg_columns_follow_the_individual_flags() {
        let opts = FieldOptions {
            individual: true,
            allele_number: true,
            symbol: true,
            ..Default::default()
        };
        assert_eq!(flag_fields(&opts)[..3], ["IND", "ZYG", "ALLELE_NUM"]);
        assert_eq!(vf_level_flag_fields(&opts), ["IND", "ZYG"]);
        assert_eq!(field_description("IND"), Some("Individual name"));
        assert_eq!(
            field_description("ZYG"),
            Some("Zygosity of individual genotype at this locus")
        );
        assert!(!flag_fields(&FieldOptions::default()).contains(&"IND"));
        assert!(!field_is_produced("ZYG", &FieldOptions::default()));

        let mut het = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        het.individual = Some("S1".into());
        het.genotype = vec!["A".into(), "G".into()];
        assert_eq!(field_value("IND", &het, None, &opts, ","), "S1");
        assert_eq!(field_value("ZYG", &het, None, &opts, ","), "HET");
        let mut hom = het.clone();
        hom.genotype = vec!["G".into(), "G".into()];
        assert_eq!(field_value("ZYG", &hom, None, &opts, ","), "HOM");
        hom.genotype = vec!["G".into()];
        assert_eq!(field_value("ZYG", &hom, None, &opts, ","), "HOM");
        let mut hom_ref = het.clone();
        hom_ref.genotype = vec!["A".into(), "A".into()];
        hom_ref.hom_ref = true;
        assert_eq!(field_value("ZYG", &hom_ref, None, &opts, ","), "HOMREF");
        hom_ref.genotype = vec!["A".into(), "*".into()];
        assert_eq!(field_value("ZYG", &hom_ref, None, &opts, ","), "HETREF");
        let plain = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(field_value("IND", &plain, None, &opts, ","), "");
        assert_eq!(field_value("ZYG", &plain, None, &opts, ","), "");
        let mut summary = het.clone();
        summary.summary_terms = Some(smallvec::smallvec![Consequence::MissenseVariant]);
        assert_eq!(field_value("IND", &summary, None, &opts, ","), "S1");
        assert_eq!(field_value("ZYG", &summary, None, &opts, ","), "HET");
        assert_eq!(field_value("SYMBOL", &summary, None, &opts, ","), "");

        let reference = hom_ref.reference_allele_row();
        assert_eq!(field_value("Allele", &reference, None, &opts, ","), "A");
        assert_eq!(field_value("ALLELE_NUM", &reference, None, &opts, ","), "0");
        assert_eq!(field_value("ZYG", &reference, None, &opts, ","), "HETREF");
        let rows = record_rows(&[&reference, &het]);
        assert!(rows.is_empty(), "no feature, no intergenic verdict, no row");
    }

    /// `--use_transcript_ref` adds `GIVEN_REF` and `USED_REF` after `RefSeq`
    /// (Constants.pm:99) with their descriptions (`:252-253`); both are the
    /// reference allele's object's on a transcript row (OutputFactory.pm:1774-1779),
    /// `USED_REF` the allele the row was analysed against and `GIVEN_REF` the
    /// input's, so a row without a transcript carries neither, and a
    /// reference row's Allele is the allele it was analysed against (`:1347`).
    /// The JSON writer keeps a `-` on both (`OutputFactory/JSON.pm:327`).
    #[test]
    fn given_ref_and_used_ref_columns_follow_the_flag() {
        let opts = FieldOptions {
            use_transcript_ref: true,
            xref_refseq: true,
            gene_phenotype: true,
            ..Default::default()
        };
        assert_eq!(
            flag_fields(&opts)[4..],
            ["RefSeq", "GIVEN_REF", "USED_REF", "GENE_PHENO"]
        );
        assert!(vf_level_flag_fields(&opts).is_empty());
        assert_eq!(
            field_description("GIVEN_REF"),
            Some("Reference allele from input")
        );
        assert_eq!(
            field_description("USED_REF"),
            Some("Reference allele as used to get consequences")
        );
        assert!(!flag_fields(&FieldOptions::default()).contains(&"USED_REF"));
        assert!(!field_is_produced("GIVEN_REF", &FieldOptions::default()));

        let given_a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"C".to_vec());
        let read_against_g = TranscriptConsequence {
            used_ref: Some("G".into()),
            ..Default::default()
        };
        let no_flag = TranscriptConsequence::default();
        assert_eq!(
            field_value("GIVEN_REF", &given_a, Some(&read_against_g), &opts, ","),
            "A"
        );
        assert_eq!(
            field_value("USED_REF", &given_a, Some(&read_against_g), &opts, ","),
            "G"
        );
        assert_eq!(
            field_value("Allele", &given_a, Some(&read_against_g), &opts, ","),
            "C"
        );
        assert_eq!(
            field_value("GIVEN_REF", &given_a, Some(&no_flag), &opts, ","),
            ""
        );
        assert_eq!(
            field_value("USED_REF", &given_a, Some(&no_flag), &opts, ","),
            ""
        );
        assert_eq!(field_value("GIVEN_REF", &given_a, None, &opts, ","), "");
        let mut summary = given_a.clone();
        summary.summary_terms = Some(smallvec::smallvec![Consequence::MissenseVariant]);
        assert_eq!(field_value("USED_REF", &summary, None, &opts, ","), "");

        let insertion = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"TT".to_vec());
        let read_against_dash = TranscriptConsequence {
            used_ref: Some("-".into()),
            ..Default::default()
        };
        assert_eq!(
            field_value(
                "GIVEN_REF",
                &insertion,
                Some(&read_against_dash),
                &opts,
                ","
            ),
            "-"
        );
        assert_eq!(
            field_value("USED_REF", &insertion, Some(&read_against_dash), &opts, ","),
            "-"
        );
        assert!(
            dash_is_a_value("USED_REF")
                && dash_is_a_value("GIVEN_REF")
                && dash_is_a_value("Allele")
        );
        assert!(!dash_is_a_value("REF_ALLELE"));

        let reference_row = given_a.reference_allele_row();
        assert_eq!(
            field_value("Allele", &reference_row, Some(&read_against_g), &opts, ","),
            "G"
        );
        assert_eq!(
            field_value("Allele", &reference_row, Some(&no_flag), &opts, ","),
            "A"
        );
        assert_eq!(
            field_value(
                "GIVEN_REF",
                &reference_row,
                Some(&read_against_g),
                &opts,
                ","
            ),
            "A"
        );
    }

    /// `--individual_zyg` adds `ZYG` alone (Constants.pm:70, the second
    /// `@FLAG_FIELDS` entry, so after `IND` and `ZYG` when both flags are on and
    /// the name listed once), set from the variant alone so a summary row keeps
    /// it. Its value is the record's list of `sample:ZYG` entries
    /// (OutputFactory.pm:981-987), joined as the format joins a list, `,` in the
    /// default and tab formats and `&` in the VCF; an empty list prints nothing; a
    /// variant without the list prints the copy's own zygosity or nothing; the
    /// reference row of a record read for one reference-homozygous sample carries
    /// the list.
    #[test]
    fn zyg_lists_every_sample_under_individual_zyg() {
        let opts = FieldOptions {
            individual_zyg: true,
            allele_number: true,
            ..Default::default()
        };
        assert_eq!(flag_fields(&opts)[..2], ["ZYG", "ALLELE_NUM"]);
        assert_eq!(vf_level_flag_fields(&opts), ["ZYG"]);
        assert!(!field_is_produced("IND", &opts));
        let both = FieldOptions {
            individual: true,
            individual_zyg: true,
            ..Default::default()
        };
        assert_eq!(flag_fields(&both)[..2], ["IND", "ZYG"]);
        assert_eq!(
            flag_fields(&both).iter().filter(|f| **f == "ZYG").count(),
            1
        );

        let entry = |sample: &str, zygosity: Zygosity| SampleZygosity {
            sample: sample.into(),
            zygosity,
        };
        let list = |entries: Vec<SampleZygosity>| Some(Arc::new(ZygosityList::new(entries)));
        let mut record = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        record.sample_zygosities = list(vec![
            entry("S1", Zygosity::Het),
            entry("S2", Zygosity::HomRef),
        ]);
        assert_eq!(
            field_value("ZYG", &record, None, &opts, ","),
            "S1:HET,S2:HOMREF"
        );
        assert_eq!(
            field_value("ZYG", &record, None, &opts, "&"),
            "S1:HET&S2:HOMREF"
        );
        assert_eq!(field_value("IND", &record, None, &opts, ","), "");
        let mut summary = record.clone();
        summary.summary_terms = Some(smallvec::smallvec![Consequence::MissenseVariant]);
        assert_eq!(
            field_value("ZYG", &summary, None, &opts, ","),
            "S1:HET,S2:HOMREF"
        );
        let mut none_called = record.clone();
        none_called.sample_zygosities = list(Vec::new());
        assert_eq!(field_value("ZYG", &none_called, None, &opts, ","), "");
        let plain = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(field_value("ZYG", &plain, None, &opts, ","), "");

        let mut hom_ref = record.clone();
        hom_ref.sample_zygosities = list(vec![entry("S2", Zygosity::HomRef)]);
        hom_ref.hom_ref = true;
        let reference = hom_ref.reference_allele_row();
        assert_eq!(field_value("Allele", &reference, None, &opts, ","), "A");
        assert_eq!(field_value("ALLELE_NUM", &reference, None, &opts, ","), "0");
        assert_eq!(
            field_value("ZYG", &reference, None, &opts, ","),
            "S2:HOMREF"
        );
    }

    /// `VARIANT_CLASS`, `MINIMISED`, `AMBIGUITY` and `CHECK_REF` are set from the
    /// variant alone (`OutputFactory.pm` `VariationFeature_to_output_hash`), so a
    /// summary row carries them and nothing else of the field table.
    #[test]
    fn variant_class_survives_on_summary_rows() {
        let opts = FieldOptions {
            variant_class: true,
            minimal: true,
            ambiguity: true,
            dont_skip: true,
            symbol: true,
            hgvs: true,
            ..Default::default()
        };
        assert_eq!(
            vf_level_flag_fields(&opts),
            ["VARIANT_CLASS", "MINIMISED", "AMBIGUITY", "CHECK_REF"]
        );
        assert_eq!(
            vf_level_flag_fields(&FieldOptions {
                variant_class: true,
                ..Default::default()
            }),
            ["VARIANT_CLASS"]
        );
        let mut summary = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        summary.most_severe_consequence = Some(Consequence::MissenseVariant);
        summary.summary_terms = Some(smallvec::smallvec![Consequence::MissenseVariant]);
        assert_eq!(field_value("AMBIGUITY", &summary, None, &opts, ","), "R");
        assert_eq!(
            field_value("VARIANT_CLASS", &summary, None, &opts, ","),
            "SNV"
        );
        assert_eq!(field_value("SYMBOL", &summary, None, &opts, ","), "");
    }

    /// `SV` is set from the variant alone too (`OutputFactory.pm:962-964`), so a
    /// summary row carries it, and only `--check_svs` lists it.
    #[test]
    fn sv_survives_on_summary_rows() {
        let opts = FieldOptions {
            check_svs: true,
            symbol: true,
            ..Default::default()
        };
        assert_eq!(vf_level_flag_fields(&opts), ["SV"]);
        assert_eq!(
            vf_level_flag_fields(&FieldOptions {
                symbol: true,
                ..Default::default()
            }),
            Vec::<&str>::new()
        );
        let mut summary = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        summary.overlapping_svs = vec!["nsv2770691".into(), "nsv533679".into()];
        summary.most_severe_consequence = Some(Consequence::MissenseVariant);
        summary.summary_terms = Some(smallvec::smallvec![Consequence::MissenseVariant]);
        assert_eq!(
            field_value("SV", &summary, None, &opts, ","),
            "nsv2770691,nsv533679"
        );
        assert_eq!(field_value("SYMBOL", &summary, None, &opts, ","), "");
    }

    /// `OutputFactory.pm` sets `DISTANCE` only when the rendered consequence list
    /// names `upstream_gene_variant` or `downstream_gene_variant`
    /// (`%DISTANCE_CONS`): under `--terms display` the row reads `UPSTREAM` and
    /// carries no distance, under `--terms NCBI` the two terms have no label and
    /// render as themselves, so the distance stays.
    #[test]
    fn distance_follows_the_rendered_terms() {
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            consequences: smallvec::smallvec![Consequence::UpstreamGeneVariant],
            distance: Some(141),
            ..Default::default()
        };
        let with = |terms: Vocabulary| FieldOptions {
            terms,
            ..Default::default()
        };
        assert_eq!(
            field_value("DISTANCE", &v, Some(&tc), &with(Vocabulary::So), ","),
            "141"
        );
        assert_eq!(
            field_value("DISTANCE", &v, Some(&tc), &with(Vocabulary::Ncbi), ","),
            "141"
        );
        assert_eq!(
            field_value("DISTANCE", &v, Some(&tc), &with(Vocabulary::Display), ","),
            ""
        );
        assert!(with(Vocabulary::So).user_group_is_plain());
        assert!(!with(Vocabulary::Display).user_group_is_plain());
    }

    /// `--total_length` appends `/<total>` to each position column that carries a
    /// total, the `-` placeholder included (`OutputFactory.pm`
    /// `TranscriptVariationAllele_to_output_hash` appends to `format_coords`'s
    /// result); without the flag the stored positions print as they are, so the
    /// numeric readers of the positions never see a total.
    #[test]
    fn total_length_renders_the_totals_beside_the_positions() {
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            cdna_position: Some("428".into()),
            cds_position: None,
            protein_position: Some("97-98".into()),
            cdna_total: Some(1200),
            cds_total: Some(306),
            protein_total: Some(101),
            ..Default::default()
        };
        let with = FieldOptions {
            total_length: true,
            ..Default::default()
        };
        assert_eq!(
            field_value("cDNA_position", &v, Some(&tc), &with, ","),
            "428/1200"
        );
        assert_eq!(
            field_value("CDS_position", &v, Some(&tc), &with, ","),
            "-/306"
        );
        assert_eq!(
            field_value("Protein_position", &v, Some(&tc), &with, ","),
            "97-98/101"
        );
        let without = FieldOptions::default();
        assert_eq!(
            field_value("cDNA_position", &v, Some(&tc), &without, ","),
            "428"
        );
        assert_eq!(
            field_value("CDS_position", &v, Some(&tc), &without, ","),
            ""
        );
        assert_eq!(position_column(Some("11"), None), "11");
        assert_eq!(position_column(None, None), "");
    }

    /// A `--summary` row lists each rendered label once (`OutputFactory.pm`
    /// `summary_only` pushes a rendered term unless already seen), so splice
    /// terms that share `SPLICE_SITE` collapse to one under `--terms display`;
    /// a feature row keeps every repeat, and the SO vocabulary, whose labels
    /// are distinct, is unchanged either way.
    #[test]
    fn summary_rows_list_each_rendered_label_once() {
        let terms = [
            Consequence::SpliceDonor5thBaseVariant,
            Consequence::SpliceRegionVariant,
            Consequence::IntronVariant,
        ];
        let tc = TranscriptConsequence {
            consequences: smallvec::SmallVec::from_slice(&terms),
            ..Default::default()
        };
        let display = FieldOptions {
            terms: Vocabulary::Display,
            ..Default::default()
        };
        assert_eq!(
            consequence_string(&tc, &display),
            "SPLICE_SITE,SPLICE_SITE,INTRONIC"
        );
        assert_eq!(
            rendered_terms(&terms, Vocabulary::Display, true),
            ["SPLICE_SITE", "INTRONIC"]
        );
        assert_eq!(
            rendered_terms(&terms, Vocabulary::So, true),
            [
                "splice_donor_5th_base_variant",
                "splice_region_variant",
                "intron_variant"
            ]
        );
        let mut summary = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        summary.most_severe_consequence = Some(terms[0]);
        summary.summary_terms = Some(smallvec::SmallVec::from_slice(&terms));
        assert_eq!(
            field_value("Consequence", &summary, None, &display, "&"),
            "SPLICE_SITE&INTRONIC"
        );
    }

    /// The `Consequence` column renders in the `--terms` vocabulary on feature
    /// rows and intergenic rows alike; the default is the SO term.
    #[test]
    fn consequence_column_renders_in_the_terms_vocabulary() {
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            consequences: smallvec::smallvec![
                Consequence::SpliceRegionVariant,
                Consequence::SynonymousVariant
            ],
            ..Default::default()
        };
        let with = |terms: Vocabulary| FieldOptions {
            terms,
            ..Default::default()
        };
        assert_eq!(
            field_value("Consequence", &v, Some(&tc), &FieldOptions::default(), ","),
            "splice_region_variant,synonymous_variant"
        );
        assert_eq!(
            field_value(
                "Consequence",
                &v,
                Some(&tc),
                &with(Vocabulary::Display),
                ","
            ),
            "SPLICE_SITE,SYNONYMOUS_CODING"
        );
        assert_eq!(
            field_value("Consequence", &v, Some(&tc), &with(Vocabulary::Ncbi), ","),
            "splice_region_variant,cds-synon"
        );
        let mut intergenic = v.clone();
        intergenic.most_severe_consequence = Some(Consequence::IntergenicVariant);
        assert_eq!(
            field_value(
                "Consequence",
                &intergenic,
                None,
                &with(Vocabulary::Display),
                ","
            ),
            "INTERGENIC"
        );
        assert_eq!(
            field_value(
                "Consequence",
                &intergenic,
                None,
                &with(Vocabulary::Ncbi),
                ","
            ),
            "intergenic_variant"
        );
    }

    fn plugin_headers() -> Vec<PluginHeader> {
        vec![
            (
                "CADD_PHRED".to_string(),
                "CADD PHRED-scaled pathogenicity score".to_string(),
            ),
            (
                "CADD_RAW".to_string(),
                "CADD raw pathogenicity score".to_string(),
            ),
        ]
    }

    /// `OutputFactory/Tab.pm` `fields` (release 116.2): the `--fields` value stands
    /// in for the whole list, defaults included, in the order given, so
    /// `Feature,Allele,SYMBOL,Location` is the column line and nothing else
    /// follows, while without it the thirteen defaults, the flag fields and the
    /// plugin fields come in that order.
    #[test]
    fn fields_replace_the_tab_columns_in_the_given_order() {
        let options = FieldOptions {
            symbol: true,
            ..Default::default()
        };
        let plugins = plugin_field_names(&plugin_headers());
        let chosen: Vec<String> = ["Feature", "Allele", "SYMBOL", "Location"]
            .map(str::to_string)
            .to_vec();
        assert_eq!(tab_columns(&options, &plugins, Some(&chosen)), chosen);
        let lines = tab_header(
            &options,
            &plugin_headers(),
            Some(&chosen),
            "/c",
            None,
            "vep -i x",
        );
        assert_eq!(lines.last().unwrap(), "#Feature\tAllele\tSYMBOL\tLocation");
        assert_eq!(
            lines
                .iter()
                .filter(|l| l.starts_with("## ") && l.contains(" : "))
                .count(),
            4
        );
        let default = tab_columns(&options, &plugins, None);
        assert_eq!(&default[..13], &DEFAULT_OUTPUT_COLS);
        assert_eq!(
            &default[13..],
            &[
                "IMPACT",
                "DISTANCE",
                "STRAND",
                "FLAGS",
                "SYMBOL",
                "SYMBOL_SOURCE",
                "HGNC_ID",
                "CADD_PHRED",
                "CADD_RAW"
            ]
        );
    }

    /// A plugin header key listed twice, as the registry hands it when one plugin
    /// is named twice on the command line (Runner.pm:925 and :994 load one instance
    /// per `--plugin`) or two plugins share a key, is one tab column here, where
    /// `OutputFactory/Tab.pm:195` appends every key it is handed and Ensembl VEP
    /// 116.2 prints `LoFtool\tLoFtool` after `FLAGS` under `--plugin LoFtool,<file>
    /// --plugin LoFtool,<file>`. The collapse is stated in `docs/output/vep-and-tab.md`.
    #[test]
    fn a_repeated_plugin_key_is_one_tab_column() {
        let options = FieldOptions::default();
        let repeated: Vec<String> = ["LoFtool", "LoFtool", "CADD_PHRED"]
            .map(str::to_string)
            .to_vec();
        let columns = tab_columns(&options, &repeated, None);
        assert_eq!(
            &columns[13..],
            &[
                "IMPACT",
                "DISTANCE",
                "STRAND",
                "FLAGS",
                "LoFtool",
                "CADD_PHRED"
            ]
        );
        assert_eq!(columns.iter().filter(|c| *c == "LoFtool").count(), 1);
    }

    /// `OutputFactory/Tab.pm` `description_headers` prints `## <name> : ?` for a
    /// field neither the field table nor a plugin describes, and
    /// `output_hash_to_line` prints `-` where the row has no such key; Ensembl VEP
    /// 116.2 under `--tab --fields Location,Allele,NOT_A_FIELD` writes the
    /// `## NOT_A_FIELD : ?` line and `-` in every row.
    #[test]
    fn unknown_field_prints_dash_and_question_mark_description() {
        let chosen: Vec<String> = ["Location", "Allele", "NOT_A_FIELD"]
            .map(str::to_string)
            .to_vec();
        let lines = tab_header(
            &FieldOptions::default(),
            &[],
            Some(&chosen),
            "/c",
            None,
            "vep -i x",
        );
        assert!(
            lines.contains(&"## NOT_A_FIELD : ?".to_string()),
            "{lines:?}"
        );
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(
            field_value("NOT_A_FIELD", &v, None, &FieldOptions::default(), ","),
            ""
        );
        let tab = crate::output::tab::TabOutputFormatter::new(
            FieldOptions::default(),
            Vec::new(),
            Some(chosen),
            Vec::new(),
            true,
        );
        let mut v = v;
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        assert_eq!(
            crate::OutputFormatter::format_record(&tab, &[&v]).unwrap(),
            vec!["21:100\tG\t-"]
        );
    }

    /// `OutputFactory.pm` sets a flag field's key only under its flag (`BIOTYPE`
    /// under `--biotype`, `SYMBOL` under `--symbol`, `AF` under `--af` or
    /// `--af_1kg`), so a `--fields` list naming an inactive field prints `-` in
    /// the tab output and an empty subfield in the VCF; the thirteen default
    /// columns, the `user` group, the overlap and check-ref keys VEP writes
    /// without their flag, and a plugin or unknown name are produced whatever
    /// the flags. Ensembl VEP 116.2 under `--tab --fields Location,BIOTYPE`
    /// prints `-` for BIOTYPE on every row of the writers corpus.
    #[test]
    fn inactive_flag_field_named_by_fields_prints_dash() {
        let off = FieldOptions::default();
        let on = FieldOptions {
            biotype: true,
            af_1kg: true,
            ..Default::default()
        };
        assert!(!field_is_produced("BIOTYPE", &off));
        assert!(field_is_produced("BIOTYPE", &on));
        assert!(!field_is_produced("AF", &off));
        assert!(
            field_is_produced("AF", &on),
            "AF is in the af_1kg group too"
        );
        for always in [
            "Location",
            "Allele",
            "IMPACT",
            "STRAND",
            "OverlapBP",
            "CHECK_REF",
            "CADD_PHRED",
            "NOT_A_FIELD",
        ] {
            assert!(field_is_produced(always, &off), "{always}");
        }
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: std::sync::Arc::from("ENST1"),
            biotype: Some("protein_coding".into()),
            ..Default::default()
        }];
        let chosen: Vec<String> = ["Feature", "BIOTYPE"].map(str::to_string).to_vec();
        let tab = crate::output::tab::TabOutputFormatter::new(
            off.clone(),
            Vec::new(),
            Some(chosen.clone()),
            Vec::new(),
            true,
        );
        assert_eq!(
            crate::OutputFormatter::format_record(&tab, &[&v]).unwrap(),
            vec!["ENST1\t-"]
        );
        let tab = crate::output::tab::TabOutputFormatter::new(
            on.clone(),
            Vec::new(),
            Some(chosen),
            Vec::new(),
            true,
        );
        assert_eq!(
            crate::OutputFormatter::format_record(&tab, &[&v]).unwrap(),
            vec!["ENST1\tprotein_coding"]
        );
    }

    /// `CHECK_REF` prints `failed` on every row of a record whose reference allele
    /// failed `--check_ref` and that `--dont_skip` kept (OutputFactory.pm:1020 sets
    /// the key from `check_ref_failed`), in the `Extra` column as `CHECK_REF=failed`
    /// and in the tab column the `dont_skip` group adds (Constants.pm:130); a record
    /// that passed prints nothing.
    #[test]
    fn check_ref_failed_renders_failed() {
        let options = FieldOptions {
            dont_skip: true,
            ..Default::default()
        };
        let mut failed = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        failed.check_ref_failed = true;
        failed.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let passed = InputVariant::new("21".into(), 200, 200, b"C".to_vec(), b"T".to_vec());
        assert_eq!(
            field_value("CHECK_REF", &failed, None, &options, ","),
            "failed"
        );
        assert_eq!(field_value("CHECK_REF", &passed, None, &options, ","), "");
        assert!(
            extra_keys(&failed, None, &options).contains(&"CHECK_REF"),
            "the dont_skip group lists the key for the Extra column"
        );
        let tab = crate::output::tab::TabOutputFormatter::new(
            options,
            Vec::new(),
            Some(["Location", "CHECK_REF"].map(str::to_string).to_vec()),
            Vec::new(),
            true,
        );
        assert_eq!(
            crate::OutputFormatter::format_record(&tab, &[&failed]).unwrap(),
            vec!["21:100\tfailed"]
        );
    }

    /// `OutputFactory/VEP_output.pm` `description_headers` (release 116.2) lists
    /// the plugin keys after the flag fields' `Extra column keys` as
    /// `## <KEY> : <description>`, the plugin's own text, before the command line.
    #[test]
    fn default_header_lists_plugin_key_descriptions() {
        let lines = default_format_header(
            &FieldOptions::default(),
            &plugin_headers(),
            "/c",
            None,
            "vep -i x",
        );
        let position = |text: &str| lines.iter().position(|l| l == text).unwrap();
        let flags = position("## FLAGS : Transcript quality flags");
        let phred = position("## CADD_PHRED : CADD PHRED-scaled pathogenicity score");
        let raw = position("## CADD_RAW : CADD raw pathogenicity score");
        let command = lines
            .iter()
            .position(|l| l.starts_with("## VEP command-line: "))
            .unwrap();
        assert!(flags < phred && phred < raw && raw < command, "{lines:?}");
        // The default format keeps its fixed column line whatever the plugins add.
        assert!(lines
            .last()
            .unwrap()
            .ends_with("\tExisting_variation\tExtra"));
    }

    /// A known variant matched to the input alleles `matched` (input allele,
    /// record allele), with the ClinVar fields `clinvar`.
    fn known(
        id: &str,
        allele_string: &str,
        matched: &[(&str, &str)],
        clin_sig: &[&str],
        clinvar: Option<ClinVarDetail>,
    ) -> ColocatedVariant {
        ColocatedVariant {
            id: id.to_string(),
            allele_string: Some(allele_string.to_string()),
            strand: 1,
            clin_sig: clin_sig.iter().map(|t| t.to_string()).collect(),
            matched_alleles: matched
                .iter()
                .enumerate()
                .map(|(i, (input, record))| MatchedAllele {
                    input_allele: input.to_string(),
                    input_index: i,
                    record_allele: record.to_string(),
                    record_index: i,
                })
                .collect(),
            clinvar: clinvar.map(Box::new),
            ..Default::default()
        }
    }

    /// The alleles of a record, each an `InputVariant` carrying the record's
    /// whole co-located list, as the matcher leaves them.
    fn record(ref_allele: &str, alts: &[&str], known: Vec<ColocatedVariant>) -> Vec<InputVariant> {
        alts.iter()
            .enumerate()
            .map(|(i, alt)| {
                let mut v = InputVariant::new(
                    "21".into(),
                    10_541_121,
                    10_541_121,
                    ref_allele.as_bytes().to_vec(),
                    alt.as_bytes().to_vec(),
                );
                v.allele_index = i;
                v.colocated_variants = known.clone();
                v
            })
            .collect()
    }

    fn with_check_existing(set: impl FnOnce(&mut FieldOptions)) -> FieldOptions {
        let mut options = FieldOptions {
            check_existing: true,
            ..Default::default()
        };
        set(&mut options);
        options
    }

    /// `rs374325784` of the release-116 GRCh38 cache: `G/A/C/T`, `clin_sig`
    /// `likely_benign`, `clin_sig_allele` `T:likely_benign`, `clin_sig_ref_allele`
    /// `G`, phenotype-linked, one ClinVar synonym group.
    fn rs374325784(matched: &[(&str, &str)]) -> ColocatedVariant {
        let mut cv = known(
            "rs374325784",
            "G/A/C/T",
            matched,
            &["likely_benign"],
            Some(ClinVarDetail {
                clin_sig_allele: Some("T:likely_benign".into()),
                clin_sig_ref_allele: Some("G".into()),
                clinical_impact: None,
            }),
        );
        cv.phenotype = true;
        cv.var_synonyms = Some("ClinVar::RCV000886041,VCV000713837".into());
        cv
    }

    /// `SOMATIC` and `PHENO` are one `0` or `1` per co-located variant in output
    /// order (germline first) and absent when every value is `0`
    /// (OutputFactory.pm:1069-1071, 1129-1138). Ensembl VEP 116.2 prints neither
    /// key for `21_5230609_A/C` of the `GRCh38-colocated` corpus (`rs1979841086`
    /// and `rs373402393`, both germline and unlinked) and `SOMATIC=1;PHENO=1` for
    /// `21_38403526_T/TTCC` (`rs2517045861`, somatic and phenotype-linked).
    #[test]
    fn somatic_and_pheno_are_per_variant_lists_dropped_when_all_zero() {
        let options = with_check_existing(|_| {});
        let unlinked = record(
            "A",
            &["C"],
            vec![
                known("rs1979841086", "dbSNP_novariation", &[], &[], None),
                known("rs373402393", "A/C", &[("C", "C")], &[], None),
            ],
        );
        assert_eq!(
            field_value("SOMATIC", &unlinked[0], None, &options, ","),
            ""
        );
        assert_eq!(field_value("PHENO", &unlinked[0], None, &options, ","), "");
        assert_eq!(
            field_value("Existing_variation", &unlinked[0], None, &options, ","),
            "",
            "the matcher fills Existing_variation; the renderer prints it as given"
        );

        let mut somatic = known("rs2517045861", "-/TCC", &[("TCC", "TCC")], &[], None);
        somatic.somatic = true;
        somatic.phenotype = true;
        let linked = record("-", &["TCC"], vec![somatic.clone()]);
        assert_eq!(field_value("SOMATIC", &linked[0], None, &options, ","), "1");
        assert_eq!(field_value("PHENO", &linked[0], None, &options, ","), "1");

        // A germline phenotype-linked variant beside a somatic unlinked one: the
        // germline variant sorts first whatever the cache order, so the lists read
        // `0,1` and `1,0`; `&` joins them inside a VCF CSQ entry.
        somatic.phenotype = false;
        let mut germline = known("rs1", "-/TCC", &[("TCC", "TCC")], &[], None);
        germline.phenotype = true;
        let both = record("-", &["TCC"], vec![somatic, germline]);
        assert_eq!(field_value("SOMATIC", &both[0], None, &options, ","), "0,1");
        assert_eq!(field_value("PHENO", &both[0], None, &options, "&"), "1&0");
        assert!(!FieldOptions::default().check_existing);
        assert!(flag_fields(&options).ends_with(&["CLIN_SIG", "SOMATIC", "PHENO"]));
    }

    /// `CLIN_SIG` without allele-specific entries lists every term of every
    /// co-located variant, duplicates kept (OutputFactory.pm:1125 pushes each
    /// variant's `split(',', clin_sig)` and nothing deduplicates it).
    #[test]
    fn clin_sig_keeps_duplicates() {
        let options = with_check_existing(|_| {});
        let alleles = record(
            "A",
            &["G"],
            vec![
                known(
                    "rs1",
                    "A/G",
                    &[("G", "G")],
                    &["benign", "likely_benign"],
                    None,
                ),
                known("rs2", "A/G", &[("G", "G")], &["benign"], None),
            ],
        );
        assert_eq!(
            field_value("CLIN_SIG", &alleles[0], None, &options, ","),
            "benign,likely_benign,benign"
        );
        assert_eq!(
            field_value("CLIN_SIG", &alleles[0], None, &options, "&"),
            "benign&likely_benign&benign"
        );
    }

    /// Under the default `--clin_sig_allele 1`, a variant with allele-specific
    /// entries contributes the entries of the row's own allele and suppresses its
    /// record-level terms (OutputFactory.pm:1097-1125, 1143-1144). Ensembl VEP
    /// 116.2 on `21 10541121 G A,T` of the `GRCh38-colocated` corpus writes
    /// `CLIN_SIG=likely_benign;PHENO=1;VAR_SYNONYMS=ClinVar::RCV000886041,VCV000713837`
    /// on the `T` rows and `PHENO=1;VAR_SYNONYMS=ClinVar::RCV000886041,VCV000713837`
    /// on the `A` rows (`rs374325784`, entries `T:likely_benign`).
    #[test]
    fn clin_sig_allele_default_reports_this_alleles_entries() {
        let options = with_check_existing(|o| o.var_synonyms = true);
        let alleles = record(
            "G",
            &["A", "T"],
            vec![rs374325784(&[("A", "A"), ("T", "T")])],
        );
        let (a, t) = (&alleles[0], &alleles[1]);
        assert_eq!(field_value("CLIN_SIG", a, None, &options, ","), "");
        assert_eq!(
            field_value("CLIN_SIG", t, None, &options, ","),
            "likely_benign"
        );
        for row in [a, t] {
            assert_eq!(field_value("PHENO", row, None, &options, ","), "1");
            assert_eq!(
                field_value("VAR_SYNONYMS", row, None, &options, ","),
                "ClinVar::RCV000886041,VCV000713837"
            );
        }
        // A record-level term read before the first allele-specific variant
        // survives when no entry names the row's allele (:1125 runs per variant).
        let mut earlier = known("rs1", "G/A/T", &[("A", "A"), ("T", "T")], &["benign"], None);
        earlier.phenotype = false;
        let alleles = record(
            "G",
            &["A", "T"],
            vec![earlier, rs374325784(&[("A", "A"), ("T", "T")])],
        );
        assert_eq!(
            field_value("CLIN_SIG", &alleles[0], None, &options, ","),
            "benign"
        );
        assert_eq!(
            field_value("CLIN_SIG", &alleles[1], None, &options, ","),
            "likely_benign",
            "the allele's own entries replace the record-level terms (:1144)"
        );
        // Two variants with entries for the row's allele: `;`-joined, each
        // variant's `,`-joined terms once (:1114-1116, 1120, 1144).
        let mut other = known("rs2", "G/A/T", &[("A", "A"), ("T", "T")], &[], None);
        other.clinvar = Some(Box::new(ClinVarDetail {
            clin_sig_allele: Some("T:pathogenic;T:risk_factor;A:benign".into()),
            clin_sig_ref_allele: Some("G".into()),
            clinical_impact: None,
        }));
        let alleles = record(
            "G",
            &["A", "T"],
            vec![rs374325784(&[("A", "A"), ("T", "T")]), other],
        );
        assert_eq!(
            field_value("CLIN_SIG", &alleles[1], None, &options, ","),
            "likely_benign;pathogenic,risk_factor"
        );
        assert_eq!(
            field_value("CLIN_SIG", &alleles[1], None, &options, "&"),
            "likely_benign;pathogenic,risk_factor",
            "the VCF writer encodes the separators itself"
        );
        assert_eq!(
            field_value("CLIN_SIG", &alleles[0], None, &options, ","),
            "benign"
        );
    }

    /// The ClinVar entry alleles of a `--minimal` split record's rows are read
    /// against the reference `REF_ALLELE` prints, the one cached on the record's
    /// feature (`$vf->ref_allele_string`, OutputFactory.pm:1049), not against
    /// each allele's own trimmed reference. Ensembl VEP 116.2 on the GRCh38 cache
    /// under `--check_existing --minimal --show_ref_allele` (rs7278737 at
    /// 21:14109044 is `G/A/T` with `clin_sig` `benign` and the entry `T:benign`
    /// against `G`): `21 14109044 G T,GT` prints `REF_ALLELE=G` and
    /// `CLIN_SIG=benign` on the `G/T` rows and on the `-/T` rows, whose Allele is
    /// `T` as well; `21 14109044 G A,GA` prints `REF_ALLELE=G` and no `CLIN_SIG`
    /// on either allele's rows, the `-/A` rows included, which the insertion's own
    /// reference `-` would read as a match once the entry allele is complemented.
    #[test]
    fn clin_sig_of_a_split_row_reads_the_record_features_reference() {
        let options = with_check_existing(|o| o.show_ref_allele = true);
        let rs7278737 = |matched: &str| {
            known(
                "rs7278737",
                "G/A/T",
                &[(matched, matched)],
                &["benign"],
                Some(ClinVarDetail {
                    clin_sig_allele: Some("T:benign".into()),
                    clin_sig_ref_allele: Some("G".into()),
                    clinical_impact: None,
                }),
            )
        };
        let split = |working: &str, alt: &str, own: usize, cv: ColocatedVariant| {
            let mut v = InputVariant::new(
                "21".into(),
                14_109_044,
                14_109_044,
                b"G".to_vec(),
                alt.as_bytes().to_vec(),
            );
            v.allele_index = own;
            v.record_allele_string_multi = Some(Arc::from(working));
            assert_eq!(
                vep_core::minimal::apply_minimal_split(&mut v, own),
                vep_core::minimal::MinimalSplit::Trimmed
            );
            v.colocated_variants = vec![cv];
            v
        };
        let snv = split("G/T/GT", "T", 0, rs7278737("T"));
        let insertion = split("G/T/GT", "GT", 1, rs7278737("T"));
        assert_eq!(
            (
                insertion.start,
                insertion.end,
                insertion.allele_string.as_str(),
                insertion.display_allele()
            ),
            (14_109_045, 14_109_044, "-/T", "T".to_string())
        );
        for row in [&snv, &insertion] {
            assert_eq!(field_value("REF_ALLELE", row, None, &options, ","), "G");
            assert_eq!(field_value("CLIN_SIG", row, None, &options, ","), "benign");
        }
        let snv = split("G/A/GA", "A", 0, rs7278737("A"));
        let insertion = split("G/A/GA", "GA", 1, rs7278737("A"));
        for row in [&snv, &insertion] {
            assert_eq!(field_value("REF_ALLELE", row, None, &options, ","), "G");
            assert_eq!(field_value("CLIN_SIG", row, None, &options, ","), "");
        }
    }

    /// `--clin_sig_allele 0` ignores the allele-specific entries: every
    /// co-located variant's record-level terms print on every row
    /// (OutputFactory.pm:1097 is skipped, :1125 pushes each variant's terms).
    /// Ensembl VEP 116.2 on `21 10541121 G A` of the `GRCh38-colocated-alt`
    /// corpus under `--clin_sig_allele 0` writes `CLIN_SIG=likely_benign` for the
    /// `A` allele, which the default run leaves without CLIN_SIG.
    #[test]
    fn clin_sig_allele_zero_reports_locus_terms() {
        let options = with_check_existing(|o| o.clin_sig_allele = ClinSigAllele::WholeVariant);
        let alleles = record("G", &["A"], vec![rs374325784(&[("A", "A")])]);
        assert_eq!(
            field_value("CLIN_SIG", &alleles[0], None, &options, ","),
            "likely_benign"
        );
        let default = with_check_existing(|_| {});
        assert_eq!(default.clin_sig_allele, ClinSigAllele::PerAllele);
        assert_eq!(
            field_value("CLIN_SIG", &alleles[0], None, &default, ","),
            ""
        );
    }

    /// When a variant's `clin_sig_ref_allele` differs from the input's reference
    /// allele, every entry allele is reverse-complemented before the row's allele
    /// is looked up (OutputFactory.pm:1049, 1100-1112): an entry `A:pathogenic`
    /// recorded against reference `G` names the `T` row of a `C/T` input. With
    /// the reference alleles equal, or none recorded, the entries are read as
    /// written.
    #[test]
    fn clin_sig_entries_reverse_complemented_when_ref_allele_differs() {
        let options = with_check_existing(|_| {});
        let entries = |ref_allele: Option<&str>| {
            let mut cv = known("rs1", "C/T", &[("T", "T")], &["pathogenic"], None);
            cv.clinvar = Some(Box::new(ClinVarDetail {
                clin_sig_allele: Some("A:pathogenic;C:benign".into()),
                clin_sig_ref_allele: ref_allele.map(str::to_string),
                clinical_impact: None,
            }));
            record("C", &["T"], vec![cv])
        };
        assert_eq!(
            field_value("CLIN_SIG", &entries(Some("G"))[0], None, &options, ","),
            "pathogenic",
            "A reverse-complemented is T, the row's allele"
        );
        assert_eq!(
            field_value("CLIN_SIG", &entries(Some("C"))[0], None, &options, ","),
            "",
            "read as written, A names no row and the record-level term is suppressed"
        );
        assert_eq!(
            field_value("CLIN_SIG", &entries(None)[0], None, &options, ","),
            ""
        );
        // The complement of `C:benign` is `G:benign`, which the `T` row never reads.
        let mut flipped = entries(Some("G"));
        flipped[0].colocated_variants[0]
            .clinvar
            .as_mut()
            .unwrap()
            .clin_sig_allele = Some("C:benign".into());
        assert_eq!(
            field_value("CLIN_SIG", &flipped[0], None, &options, ","),
            ""
        );
    }

    /// An entry without an allele (`:pathogenic`) is keyed by the empty string,
    /// which no row's allele equals, and still marks the variant as carrying
    /// allele-specific entries, so its record-level terms are suppressed too
    /// (OutputFactory.pm:1106-1125).
    #[test]
    fn empty_allele_entry_matches_nothing() {
        let options = with_check_existing(|_| {});
        let mut cv = known("rs1", "C/T", &[("T", "T")], &["pathogenic"], None);
        cv.clinvar = Some(Box::new(ClinVarDetail {
            clin_sig_allele: Some(":pathogenic".into()),
            clin_sig_ref_allele: Some("C".into()),
            clinical_impact: None,
        }));
        let alleles = record("C", &["T"], vec![cv]);
        assert_eq!(
            field_value("CLIN_SIG", &alleles[0], None, &options, ","),
            ""
        );
    }

    /// `CLINVAR_SOMATIC_CLASSIFICATION` is the cache's `clinical_impact` string
    /// as written, the last co-located variant's in output order
    /// (OutputFactory.pm:1093-1095), under `--clinvar_somatic_classification`
    /// alone (Constants.pm:82), described as Ensembl VEP describes it (:240).
    #[test]
    fn clinvar_somatic_classification_prints_cache_string() {
        let classification = |id: &str, text: &str| {
            let mut cv = known(id, "C/T", &[("T", "T")], &[], None);
            cv.clinvar = Some(Box::new(ClinVarDetail {
                clin_sig_allele: None,
                clin_sig_ref_allele: None,
                clinical_impact: Some(text.into()),
            }));
            cv
        };
        let alleles = record(
            "C",
            &["T"],
            vec![
                classification("rs1", "Neoplasm(Tier I - Strong)"),
                classification("rs2", "Carcinoma(oncogenicity:Oncogenic)"),
            ],
        );
        let options = with_check_existing(|o| o.clinvar_somatic_classification = true);
        assert_eq!(
            field_value(
                "CLINVAR_SOMATIC_CLASSIFICATION",
                &alleles[0],
                None,
                &options,
                ","
            ),
            "Carcinoma(oncogenicity:Oncogenic)"
        );
        assert!(flag_fields(&options).contains(&"CLINVAR_SOMATIC_CLASSIFICATION"));
        let without = with_check_existing(|_| {});
        assert_eq!(
            field_value(
                "CLINVAR_SOMATIC_CLASSIFICATION",
                &alleles[0],
                None,
                &without,
                ","
            ),
            ""
        );
        assert!(!flag_fields(&without).contains(&"CLINVAR_SOMATIC_CLASSIFICATION"));
        assert_eq!(
            field_description("CLINVAR_SOMATIC_CLASSIFICATION"),
            Some("ClinVar somatic classification of the dbSNP variant")
        );
    }

    /// A region file's empty `clin_sig_allele` and `clinical_impact` cells are the
    /// defined empty string, and Ensembl VEP 116.2 on such a cache prints, for
    /// `21 10602054 ATA A` co-located with rs149363218 (`clin_sig`
    /// `uncertain_significance`, both cells empty): no `CLIN_SIG` under the default
    /// `--clin_sig_allele 1` (the defined cell sets `clin_sig_allele_exists`,
    /// OutputFactory.pm:1097, 1121, and :1125 skips the record's term),
    /// `CLIN_SIG=uncertain_significance` under `--clin_sig_allele 0`, and under
    /// `--clinvar_somatic_classification` the key with an empty value
    /// (`CLINVAR_SOMATIC_CLASSIFICATION=` in the default format, an empty tab cell:
    /// :1093 sets the key to the defined `''`). A record whose cell is unset
    /// prints the term, and the key only when the flag is set; a `--summary` row
    /// never carries the key.
    #[test]
    fn region_file_empty_clinvar_cells_suppress_clin_sig_and_write_an_empty_key() {
        let region = |clinvar: Option<ClinVarDetail>| {
            let cv = known(
                "rs149363218",
                "ATA/A",
                &[("-", "-")],
                &["uncertain_significance"],
                clinvar,
            );
            record("TA", &["-"], vec![cv])
        };
        let defined_empty = region(Some(ClinVarDetail {
            clin_sig_allele: Some(String::new()),
            clin_sig_ref_allele: None,
            clinical_impact: Some(String::new()),
        }));
        let unset = region(None);
        let default = with_check_existing(|_| {});
        let whole = with_check_existing(|o| o.clin_sig_allele = ClinSigAllele::WholeVariant);
        let somatic = with_check_existing(|o| o.clinvar_somatic_classification = true);
        assert_eq!(
            field_value("CLIN_SIG", &defined_empty[0], None, &default, ","),
            ""
        );
        assert_eq!(
            field_value("CLIN_SIG", &defined_empty[0], None, &whole, ","),
            "uncertain_significance"
        );
        assert_eq!(
            field_value("CLIN_SIG", &unset[0], None, &default, ","),
            "uncertain_significance"
        );
        for (alleles, options, written) in [
            (&defined_empty, &somatic, true),
            (&defined_empty, &default, false),
            (&unset, &somatic, false),
        ] {
            assert_eq!(
                field_value(
                    "CLINVAR_SOMATIC_CLASSIFICATION",
                    &alleles[0],
                    None,
                    options,
                    ","
                ),
                ""
            );
            assert_eq!(
                empty_value_is_written(
                    "CLINVAR_SOMATIC_CLASSIFICATION",
                    &alleles[0],
                    None,
                    options
                ),
                written
            );
        }
        // The `--summary` row is the record's hash alone (OutputFactory.pm:426-458)
        // and carries no co-located key, so the empty value is not written there.
        let mut summary = defined_empty[0].clone();
        summary.summary_terms = Some(smallvec::smallvec![Consequence::MissenseVariant]);
        assert_eq!(
            field_value(
                "CLINVAR_SOMATIC_CLASSIFICATION",
                &summary,
                None,
                &somatic,
                ","
            ),
            ""
        );
        assert!(!empty_value_is_written(
            "CLINVAR_SOMATIC_CLASSIFICATION",
            &summary,
            None,
            &somatic
        ));
        let tc = TranscriptConsequence::default();
        assert!(empty_value_is_written(
            "CLINVAR_SOMATIC_CLASSIFICATION",
            &summary,
            Some(&tc),
            &somatic
        ));
    }

    /// `docs/output/known-variants.md` states the subset: the cache's
    /// classification string as written, which Ensembl VEP builds from ClinVar's
    /// somatic records when it dumps the cache or reads a database and vep-rs
    /// never rebuilds.
    #[test]
    fn documented_subset_is_stated() {
        let page = std::fs::read_to_string(format!(
            "{}/../../docs/output/known-variants.md",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let anchor = "`CLINVAR_SOMATIC_CLASSIFICATION` (under `--clinvar_somatic_classification`)";
        let paragraph = page
            .lines()
            .find(|l| l.contains(anchor))
            .expect("docs/output/known-variants.md describes CLINVAR_SOMATIC_CLASSIFICATION");
        for phrase in [
            "the cache's classification string as written",
            "Ensembl VEP builds that string from ClinVar's somatic records",
            "vep-rs reads a cache only",
        ] {
            assert!(
                paragraph.contains(phrase),
                "the CLINVAR_SOMATIC_CLASSIFICATION paragraph lacks {phrase:?}: {paragraph}"
            );
        }
    }

    /// `PUBMED` lists every id of every co-located variant in output order,
    /// duplicates kept (OutputFactory.pm:1126), under `--pubmed` alone. Ensembl
    /// VEP 116.2 writes `PUBMED=26076356` for `21_10602054_ATA/A` of the
    /// `GRCh38-colocated` corpus (`rs149363218`).
    #[test]
    fn pubmed_keeps_duplicates_across_variants() {
        let cited = |id: &str, ids: &[&str]| {
            let mut cv = known(id, "ATA/A", &[("-", "-")], &[], None);
            cv.pubmed = ids.iter().map(|i| i.to_string()).collect();
            cv
        };
        let options = with_check_existing(|o| o.pubmed = true);
        let one = record("TA", &["-"], vec![cited("rs149363218", &["26076356"])]);
        assert_eq!(
            field_value("PUBMED", &one[0], None, &options, ","),
            "26076356"
        );
        let two = record(
            "TA",
            &["-"],
            vec![
                cited("rs1", &["26076356", "11111111"]),
                cited("rs2", &["26076356"]),
            ],
        );
        assert_eq!(
            field_value("PUBMED", &two[0], None, &options, ","),
            "26076356,11111111,26076356"
        );
        assert_eq!(
            field_value("PUBMED", &two[0], None, &options, "&"),
            "26076356&11111111&26076356"
        );
        let without = with_check_existing(|_| {});
        assert_eq!(field_value("PUBMED", &two[0], None, &without, ","), "");
    }

    /// `VAR_SYNONYMS` joins each co-located variant's synonym string with `--`
    /// (OutputFactory.pm:1090, 1141) under `--var_synonyms` alone (Constants.pm:133).
    /// Ensembl VEP 116.2 writes `VAR_SYNONYMS=ClinVar::RCV000886041,VCV000713837`
    /// for `21_10541121_G/A/T` of the `GRCh38-colocated` corpus (`rs374325784`).
    #[test]
    fn var_synonyms_joined_with_double_dash() {
        let options = with_check_existing(|o| o.var_synonyms = true);
        let one = record("G", &["A"], vec![rs374325784(&[("A", "A")])]);
        assert_eq!(
            field_value("VAR_SYNONYMS", &one[0], None, &options, ","),
            "ClinVar::RCV000886041,VCV000713837"
        );
        let mut other = known("rs2", "G/A", &[("A", "A")], &[], None);
        other.var_synonyms = Some("OMIM::605717.0001--UniProt::VAR_089699".into());
        let two = record("G", &["A"], vec![rs374325784(&[("A", "A")]), other]);
        assert_eq!(
            field_value("VAR_SYNONYMS", &two[0], None, &options, "&"),
            "ClinVar::RCV000886041,VCV000713837--OMIM::605717.0001--UniProt::VAR_089699",
            "the cache string is joined as written; the VCF writer encodes its commas"
        );
        assert!(flag_fields(&options).contains(&"VAR_SYNONYMS"));
        let without = with_check_existing(|_| {});
        assert_eq!(
            field_value("VAR_SYNONYMS", &two[0], None, &without, ","),
            ""
        );
        assert!(!flag_fields(&without).contains(&"VAR_SYNONYMS"));
        assert_eq!(
            field_description("VAR_SYNONYMS"),
            Some("List of known variation synonyms and their sources")
        );
    }

    /// A variant matched on another allele of the record prints none of its
    /// columns on this allele's rows (OutputFactory.pm:1081-1083); one matched on
    /// its span alone prints on every row.
    #[test]
    fn colocated_columns_follow_the_row_allele_gate() {
        let options = with_check_existing(|o| o.pubmed = true);
        let mut on_t = known("rs1", "G/T", &[("T", "T")], &["benign"], None);
        on_t.pubmed = vec!["1".into()];
        on_t.phenotype = true;
        let mut span_only = known("rs2", "dbSNP_novariation", &[], &["pathogenic"], None);
        span_only.somatic = true;
        let alleles = record("G", &["A", "T"], vec![on_t, span_only]);
        let (a, t) = (&alleles[0], &alleles[1]);
        assert_eq!(
            field_value("CLIN_SIG", a, None, &options, ","),
            "pathogenic"
        );
        assert_eq!(field_value("PUBMED", a, None, &options, ","), "");
        assert_eq!(field_value("SOMATIC", a, None, &options, ","), "1");
        assert_eq!(field_value("PHENO", a, None, &options, ","), "");
        assert_eq!(
            field_value("CLIN_SIG", t, None, &options, ","),
            "benign,pathogenic"
        );
        assert_eq!(field_value("PUBMED", t, None, &options, ","), "1");
        assert_eq!(field_value("SOMATIC", t, None, &options, ","), "0,1");
        assert_eq!(field_value("PHENO", t, None, &options, ","), "1,0");
    }

    /// The field groups of `@FLAG_FIELDS` (Constants.pm:60-138) that the
    /// co-located variant flags add, in Ensembl VEP's order.
    #[test]
    fn colocated_flag_groups_sit_in_ensembl_veps_order() {
        let options = FieldOptions {
            variant_class: true,
            clinvar_somatic_classification: true,
            symbol: true,
            max_af: true,
            check_frequency: true,
            check_existing: true,
            pubmed: true,
            check_svs: true,
            overlaps: true,
            var_synonyms: true,
            ..Default::default()
        };
        assert_eq!(
            flag_fields(&options)[4..],
            [
                "VARIANT_CLASS",
                "CLINVAR_SOMATIC_CLASSIFICATION",
                "SYMBOL",
                "SYMBOL_SOURCE",
                "HGNC_ID",
                "MAX_AF",
                "MAX_AF_POPS",
                "FREQS",
                "CLIN_SIG",
                "SOMATIC",
                "PHENO",
                "PUBMED",
                "SV",
                "OverlapBP",
                "OverlapPC",
                "VAR_SYNONYMS",
            ]
        );
        for field in ["SV", "FREQS"] {
            assert!(field_description(field).is_some(), "{field}");
        }
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.overlapping_svs = vec!["esv1".into(), "nsv2".into()];
        assert_eq!(field_value("SV", &v, None, &options, ","), "esv1,nsv2");
        assert_eq!(field_value("SV", &v, None, &options, "&"), "esv1&nsv2");
    }

    /// `FREQS` lists the record's `POPULATION:ALLELE:frequency` triplets as
    /// `add_colocated_variant_info` writes them (OutputFactory.pm:1148-1163),
    /// `,`-joined in the default and tab formats and `&`-joined in the VCF
    /// entry, on every row of the record and empty without the filter's data:
    /// Ensembl VEP 116.2 prints `FREQS=AF:A:NA,AF:G:1` on `21 7759128 T A,G`
    /// under `--check_frequency --freq_pop AF` and `AF:A:NA&AF:G:1` in its CSQ.
    #[test]
    fn freqs_column_lists_population_allele_frequency_triplets() {
        let options = FieldOptions {
            check_frequency: true,
            ..Default::default()
        };
        let mut v = InputVariant::new("21".into(), 7759128, 7759128, b"T".to_vec(), b"A".to_vec());
        assert_eq!(field_value("FREQS", &v, None, &options, ","), "");
        v.freqs = vec!["AF:A:NA".into(), "AF:G:1".into()];
        assert_eq!(
            field_value("FREQS", &v, None, &options, ","),
            "AF:A:NA,AF:G:1"
        );
        assert_eq!(
            field_value("FREQS", &v, None, &options, "&"),
            "AF:A:NA&AF:G:1"
        );
        let tc = TranscriptConsequence::default();
        assert_eq!(
            field_value("FREQS", &v, Some(&tc), &options, ","),
            "AF:A:NA,AF:G:1",
            "a transcript row carries the record's list"
        );
        assert!(flag_fields(&options).contains(&"FREQS"));
        assert!(!flag_fields(&FieldOptions::default()).contains(&"FREQS"));
    }

    /// The fields the JSON factory stubs out of the consequence hash
    /// (OutputFactory/JSON.pm:493-515): every co-located variant and frequency key.
    #[test]
    fn colocated_fields_are_the_stubbed_keys() {
        for field in [
            "Existing_variation",
            "CLIN_SIG",
            "SOMATIC",
            "PHENO",
            "PUBMED",
            "VAR_SYNONYMS",
            "CLINVAR_SOMATIC_CLASSIFICATION",
            "FREQS",
            "AF",
            "AFR_AF",
            "gnomADg_AMI_AF",
            "MAX_AF",
            "MAX_AF_POPS",
        ] {
            assert!(colocated_field(field), "{field}");
        }
        for field in [
            "Allele",
            "IMPACT",
            "SV",
            "VARIANT_CLASS",
            "CADD_RAW",
            "HGVSc",
            "gnomAD_AF",
        ] {
            assert!(
                !colocated_field(field),
                "{field} is a plugin's or the row's own"
            );
        }
    }

    /// `OutputFactory/Tab.pm` `description_headers` describes a plugin column with
    /// the plugin's own text (`%other_descs`), not `?`.
    #[test]
    fn tab_header_describes_plugin_columns_with_their_own_text() {
        let lines = tab_header(
            &FieldOptions::default(),
            &plugin_headers(),
            None,
            "/c",
            None,
            "vep -i x",
        );
        assert!(
            lines.contains(&"## CADD_PHRED : CADD PHRED-scaled pathogenicity score".to_string())
        );
        assert!(lines.contains(&"## CADD_RAW : CADD raw pathogenicity score".to_string()));
        assert!(!lines.iter().any(|l| l.ends_with(" : ?")), "{lines:?}");
        assert!(lines
            .last()
            .unwrap()
            .ends_with("\tFLAGS\tCADD_PHRED\tCADD_RAW"));
    }

    /// `description_headers` (Tab.pm:144-156) describes every `--fields` name
    /// from `%FIELD_DESCRIPTIONS` (Constants.pm:142-266), whether or not the run
    /// produces it, and a name outside that table as `?`.
    #[test]
    fn tab_header_describes_every_field_of_the_constants_table() {
        let fields: Vec<String> = [
            "ID",
            "REFSEQ_MATCH",
            "BAM_EDIT",
            "CELL_TYPE",
            "ExAC_AF",
            "NOT_A_FIELD",
        ]
        .map(String::from)
        .to_vec();
        let lines = tab_header(
            &FieldOptions::default(),
            &[],
            Some(&fields),
            "/c",
            None,
            "vep -i x",
        );
        for expected in [
            "## ID : Identifier of uploaded variant",
            "## REFSEQ_MATCH : RefSeq transcript match status",
            "## BAM_EDIT : Indicates success or failure of edit using BAM file",
            "## CELL_TYPE : List of cell types and classifications for regulatory feature",
            "## ExAC_AF : Frequency of existing variant in ExAC combined population",
            "## NOT_A_FIELD : ?",
        ] {
            assert!(
                lines.contains(&expected.to_string()),
                "{expected}\n{lines:?}"
            );
        }
        assert_eq!(
            lines.last().unwrap(),
            "#ID\tREFSEQ_MATCH\tBAM_EDIT\tCELL_TYPE\tExAC_AF\tNOT_A_FIELD"
        );
    }
}
