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

use std::collections::BTreeSet;

use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::sv_routing::perl_true;
use vep_core::variant::{colocated_in_output_order, ColocatedVariant, InputVariant, VariantClass};
use vep_core::variation::reverse_comp;

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
    pub mane_select: bool,
    pub mane: bool,
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
    pub af: bool,
    pub af_1kg: bool,
    pub af_gnomade: bool,
    pub af_gnomadg: bool,
    pub max_af: bool,
    pub check_frequency: bool,
    pub check_existing: bool,
    pub pubmed: bool,
    pub check_svs: bool,
    pub dont_skip: bool,
    pub overlaps: bool,
    /// `--shift_length`: the SHIFT_LENGTH field is listed; its value is written
    /// only under `shift_3prime` (OutputFactory.pm:1695).
    pub shift_length: bool,
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
    /// VEP's `--no_escape` (implied by `--json`): HGVSp keeps its `=` instead of
    /// the `%3D` the text formats print.
    pub no_escape: bool,
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

/// The flag-driven field groups in VEP's order (`@FLAG_FIELDS`). The `user`
/// group is unconditional and comes first among the groups that can fire here.
const FLAG_FIELD_GROUPS: &[(&str, &[&str])] = &[
    ("allele_number", &["ALLELE_NUM"]),
    ("show_ref_allele", &["REF_ALLELE"]),
    ("uploaded_allele", &["UPLOADED_ALLELE"]),
    ("user", &["IMPACT", "DISTANCE", "STRAND", "FLAGS"]),
    ("include_pick", &["PICK"]),
    ("variant_class", &["VARIANT_CLASS"]),
    ("minimal", &["MINIMISED"]),
    (
        "clinvar_somatic_classification",
        &["CLINVAR_SOMATIC_CLASSIFICATION"],
    ),
    ("symbol", &["SYMBOL", "SYMBOL_SOURCE", "HGNC_ID"]),
    ("biotype", &["BIOTYPE"]),
    ("canonical", &["CANONICAL"]),
    ("mane_select", &["MANE", "MANE_SELECT"]),
    ("mane", &["MANE", "MANE_SELECT", "MANE_PLUS_CLINICAL"]),
    ("tsl", &["TSL"]),
    ("gencode_primary", &["GENCODE_PRIMARY"]),
    ("appris", &["APPRIS"]),
    ("ccds", &["CCDS"]),
    ("protein", &["ENSP"]),
    (
        "uniprot",
        &["SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM"],
    ),
    ("xref_refseq", &["RefSeq"]),
    ("gene_phenotype", &["GENE_PHENO"]),
    ("sift", &["SIFT"]),
    ("polyphen", &["PolyPhen"]),
    ("numbers", &["EXON", "INTRON"]),
    ("domains", &["DOMAINS"]),
    ("mirna", &["miRNA"]),
    ("hgvs", &["HGVSc", "HGVSp", "HGVS_OFFSET"]),
    ("hgvsg", &["HGVSg"]),
    ("af", &["AF"]),
    (
        "af_1kg",
        &["AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF"],
    ),
    (
        "af_gnomade",
        &[
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
        ],
    ),
    (
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
    ("max_af", &["MAX_AF", "MAX_AF_POPS"]),
    ("check_frequency", &["FREQS"]),
    ("check_existing", &["CLIN_SIG", "SOMATIC", "PHENO"]),
    ("pubmed", &["PUBMED"]),
    ("check_svs", &["SV"]),
    ("dont_skip", &["CHECK_REF"]),
    ("overlaps", &["OverlapBP", "OverlapPC"]),
    ("shift_length", &["SHIFT_LENGTH"]),
    ("var_synonyms", &["VAR_SYNONYMS"]),
    (
        "regulatory",
        &[
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
            "sift" => self.sift,
            "polyphen" => self.polyphen,
            "numbers" => self.numbers,
            "domains" => self.domains,
            "mirna" => self.mirna,
            "hgvs" => self.hgvs,
            "hgvsg" => self.hgvsg,
            "af" => self.af,
            "af_1kg" => self.af_1kg,
            "af_gnomade" => self.af_gnomade,
            "af_gnomadg" => self.af_gnomadg,
            "max_af" => self.max_af,
            "check_frequency" => self.check_frequency,
            "check_existing" => self.check_existing,
            "pubmed" => self.pubmed,
            "check_svs" => self.check_svs,
            "dont_skip" => self.dont_skip,
            "overlaps" => self.overlaps,
            "shift_length" => self.shift_length,
            "var_synonyms" => self.var_synonyms,
            "regulatory" => self.regulatory,
            _ => false,
        }
    }

    /// True when any group beyond the unconditional `user` group is on, i.e.
    /// when a row can carry a field that is not IMPACT, DISTANCE, STRAND or FLAGS.
    pub fn any_optional_group(&self) -> bool {
        FLAG_FIELD_GROUPS
            .iter()
            .any(|(g, _)| *g != "user" && self.group_active(g))
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
    for (group, fields) in FLAG_FIELD_GROUPS {
        if fields.contains(&field) {
            if UNGUARDED_GROUPS.contains(group) || options.group_active(group) {
                return true;
            }
            listed = true;
        }
    }
    !listed
}

/// The optional fields a configuration emits, in VEP's order, each name once.
pub fn flag_fields(options: &FieldOptions) -> Vec<&'static str> {
    let mut seen: BTreeSet<&str> = BTreeSet::new();
    let mut out = Vec::new();
    for (group, fields) in FLAG_FIELD_GROUPS {
        if !options.group_active(group) {
            continue;
        }
        for f in *fields {
            if seen.insert(f) {
                out.push(*f);
            }
        }
    }
    out
}

/// VEP's description of a field, for the `## Column descriptions:` and
/// `## Extra column keys:` header blocks.
pub fn field_description(field: &str) -> Option<&'static str> {
    FIELD_DESCRIPTIONS
        .iter()
        .find(|(k, _)| *k == field)
        .map(|(_, d)| *d)
}

const FIELD_DESCRIPTIONS: &[(&str, &str)] = &[
    ("Uploaded_variation", "Identifier of uploaded variant"),
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
    ("UPLOADED_ALLELE", "The variant allele as it was uploaded"),
    ("STRAND", "Strand of the feature (1/-1)"),
    (
        "PICK",
        "Indicates if this consequence has been picked as the most severe",
    ),
    ("SOMATIC", "Somatic status of existing variant"),
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

/// One output row: an allele of the record and the transcript consequence it is
/// reported against, or `None` for the record's intergenic row or its `--summary`
/// or `--most_severe` row (the allele carrying [`InputVariant::summary_terms`]).
pub type Row<'a> = (&'a InputVariant, Option<&'a TranscriptConsequence>);

/// The rows one input record produces, in VEP's order: transcript-major (by
/// transcript id) and allele-minor (ALT order) across every allele split from
/// the record, then one intergenic row per allele that overlapped nothing, or
/// the one `--summary` or `--most_severe` row the record's first allele carries.
///
/// An allele whose chromosome had no transcripts at all (no `most_severe_consequence`)
/// produces no row, as VEP writes nothing for it.
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
    rows.sort_by(|a, b| {
        let (ta, tb) = (a.1.expect("feature row"), b.1.expect("feature row"));
        ta.transcript_id
            .as_ref()
            .cmp(tb.transcript_id.as_ref())
            .then(a.0.allele_index.cmp(&b.0.allele_index))
    });
    for variant in record {
        if variant.transcript_consequences.is_empty() && variant.most_severe_consequence.is_some() {
            rows.push((variant, None));
        }
    }
}

/// Consequence terms of a row joined by `,`.
pub fn consequence_string(tc: &TranscriptConsequence) -> String {
    let mut out = String::new();
    for (i, c) in tc.consequences.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(c.so_term());
    }
    out
}

fn first_frequency(variant: &InputVariant, pop: &str) -> Option<f64> {
    for cv in &variant.colocated_variants {
        if let Some(&f) = cv.frequencies.get(pop) {
            return Some(f);
        }
    }
    None
}

/// Frequencies print as VEP prints them (`%.4g`): four significant digits,
/// trailing zeros trimmed, scientific notation outside [1e-4, 1e4).
pub fn format_freq(f: f64) -> String {
    if f == 0.0 || !f.is_finite() {
        return "0".to_string();
    }
    let abs = f.abs();
    let mut out = if !(0.0001..10_000.0).contains(&abs) {
        format!("{f:.3e}")
    } else {
        let digits_before_decimal = abs.log10().floor() as i32 + 1;
        let decimals = (4 - digits_before_decimal).max(0) as usize;
        format!("{f:.decimals$}")
    };
    if let Some((mantissa, exponent)) = out.split_once('e') {
        let trimmed = mantissa.trim_end_matches('0').trim_end_matches('.');
        out = format!("{trimmed}e{exponent}");
    } else if out.contains('.') {
        out = out.trim_end_matches('0').trim_end_matches('.').to_string();
    }
    out
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
pub fn field_value(
    field: &str,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
    list_sep: &str,
) -> String {
    if let (None, Some(terms)) = (tc, &variant.summary_terms) {
        return match field {
            "Consequence" => terms
                .iter()
                .map(|c| c.so_term())
                .collect::<Vec<_>>()
                .join(list_sep),
            // The keys `VariationFeature_to_output_hash` sets (OutputFactory.pm:883-1025
            // at release/116.2) that the field table carries.
            "Uploaded_variation" | "Location" | "VARIANT_CLASS" | "MINIMISED" | "CHECK_REF" => {
                row_field_value(field, variant, None, options, list_sep)
            }
            _ => String::new(),
        };
    }
    row_field_value(field, variant, tc, options, list_sep)
}

/// [`field_value`] of a row against a feature or an allele's intergenic row.
fn row_field_value(
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
        "Allele" => match genomic_shift(variant, tc, options) {
            Some((shift, strand)) => variant.display_allele_shifted(shift, strand),
            None => variant.display_allele(),
        },
        "Consequence" => tc
            .map(consequence_string)
            .or_else(|| {
                variant
                    .most_severe_consequence
                    .map(|c| c.so_term().to_string())
            })
            .unwrap_or_else(|| Consequence::IntergenicVariant.so_term().to_string()),
        "IMPACT" => tc
            .map(|x| x.impact.to_string())
            .or_else(|| {
                variant
                    .most_severe_consequence
                    .map(|c| c.impact().to_string())
            })
            .unwrap_or_else(|| Consequence::IntergenicVariant.impact().to_string()),
        "SYMBOL" => tc
            .and_then(|x| x.gene_symbol.as_deref())
            .unwrap_or("")
            .to_string(),
        "Gene" => tc.map(|x| x.gene_id.as_ref()).unwrap_or("").to_string(),
        "Feature_type" => tc.map(|x| x.feature_type.to_string()).unwrap_or_default(),
        "Feature" => tc
            .map(|x| x.transcript_id.as_ref())
            .unwrap_or("")
            .to_string(),
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
        "cDNA_position" => tc.and_then(|x| x.cdna_position.clone()).unwrap_or_default(),
        "CDS_position" => tc.and_then(|x| x.cds_position.clone()).unwrap_or_default(),
        "Protein_position" => tc
            .and_then(|x| x.protein_position.clone())
            .unwrap_or_default(),
        "Amino_acids" => tc.and_then(|x| x.amino_acids.clone()).unwrap_or_default(),
        "Codons" => tc.and_then(|x| x.codons.clone()).unwrap_or_default(),
        "Existing_variation" => variant.existing_variation.join(list_sep),
        "ALLELE_NUM" => variant
            .allele_num
            .unwrap_or(variant.allele_index as u32 + 1)
            .to_string(),
        "REF_ALLELE" => String::from_utf8_lossy(&variant.ref_allele).to_string(),
        "UPLOADED_ALLELE" => variant
            .uploaded_allele_string
            .as_deref()
            .unwrap_or(&variant.allele_string)
            .to_string(),
        "DISTANCE" => tc
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
        "VARIANT_CLASS" => variant.variant_class.so_term().to_string(),
        "MINIMISED" => {
            if variant.minimised {
                "1".to_string()
            } else {
                String::new()
            }
        }
        "SYMBOL_SOURCE" => tc
            .and_then(|x| x.gene_symbol_source.as_deref())
            .unwrap_or("")
            .to_string(),
        "HGNC_ID" => tc
            .and_then(|x| x.hgnc_id.as_deref())
            .unwrap_or("")
            .to_string(),
        "CANONICAL" => {
            if tc.map(|x| x.canonical).unwrap_or(false) {
                "YES".to_string()
            } else {
                String::new()
            }
        }
        "MANE" => {
            let mut sets = Vec::new();
            if tc.and_then(|x| x.mane_select.as_ref()).is_some() {
                sets.push("MANE_Select");
            }
            if tc.and_then(|x| x.mane_plus_clinical.as_ref()).is_some() {
                sets.push("MANE_Plus_Clinical");
            }
            sets.join(list_sep)
        }
        "MANE_SELECT" => tc.and_then(|x| x.mane_select.clone()).unwrap_or_default(),
        "MANE_PLUS_CLINICAL" => tc
            .and_then(|x| x.mane_plus_clinical.clone())
            .unwrap_or_default(),
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
        "APPRIS" => tc.and_then(|x| x.appris.clone()).unwrap_or_default(),
        "CCDS" => tc.and_then(|x| x.ccds.clone()).unwrap_or_default(),
        "ENSP" => tc
            .and_then(|x| x.protein_id.as_deref())
            .unwrap_or("")
            .to_string(),
        "SWISSPROT" => tc.and_then(|x| x.swissprot.clone()).unwrap_or_default(),
        "TREMBL" => tc.and_then(|x| x.trembl.clone()).unwrap_or_default(),
        "HGVSg" => variant.hgvsg.clone().unwrap_or_default(),
        "UNIPARC" | "UNIPROT_ISOFORM" | "miRNA" | "CHECK_REF" => String::new(),
        "RefSeq" => tc.and_then(|x| x.refseq.clone()).unwrap_or_default(),
        "GENE_PHENO" => {
            if options.gene_phenotype && variant.colocated_variants.iter().any(|cv| cv.phenotype) {
                "1".to_string()
            } else {
                String::new()
            }
        }
        "DOMAINS" => tc
            .map(|x| {
                x.domains
                    .iter()
                    .map(|(src, id)| format!("{src}:{id}"))
                    .collect::<Vec<_>>()
                    .join(list_sep)
            })
            .unwrap_or_default(),
        "HGVS_OFFSET" => tc
            .and_then(|x| x.hgvs_offset.map(|o| o.to_string()))
            .unwrap_or_default(),
        // Set on every transcript row, 0 included, when both flags are on
        // (OutputFactory.pm:1690-1695); an intergenic row carries none.
        "SHIFT_LENGTH" if options.shift_3prime => {
            tc.map(|x| x.shift_length.to_string()).unwrap_or_default()
        }
        "SHIFT_LENGTH" => String::new(),
        "CLIN_SIG"
        | "SOMATIC"
        | "PHENO"
        | "PUBMED"
        | "VAR_SYNONYMS"
        | "CLINVAR_SOMATIC_CLASSIFICATION" => colocated_column(field, variant, options, list_sep),
        "SV" => variant.overlapping_svs.join(list_sep),
        "SIFT" => tc.and_then(|x| x.sift.clone()).unwrap_or_default(),
        "PolyPhen" => tc.and_then(|x| x.polyphen.clone()).unwrap_or_default(),
        "MOTIF_NAME"
        | "MOTIF_POS"
        | "HIGH_INF_POS"
        | "MOTIF_SCORE_CHANGE"
        | "TRANSCRIPTION_FACTORS" => String::new(),
        "AF" => first_frequency(variant, "AF")
            .map(format_freq)
            .unwrap_or_default(),
        f if f.ends_with("_AF") => {
            // `AFR_AF` -> population `AFR`; `gnomADe_NFE_AF` -> `gnomADe_NFE`; `gnomADe_AF` -> `gnomADe`.
            let pop = &f[..f.len() - 3];
            first_frequency(variant, pop)
                .map(format_freq)
                .unwrap_or_default()
        }
        "MAX_AF" => {
            let mut max_af: Option<f64> = None;
            for cv in &variant.colocated_variants {
                for value in cv.frequencies.values() {
                    if max_af.is_none_or(|m| *value > m) {
                        max_af = Some(*value);
                    }
                }
            }
            max_af.map(format_freq).unwrap_or_default()
        }
        "MAX_AF_POPS" => {
            let mut max_af: Option<f64> = None;
            let mut max_pop = String::new();
            for cv in &variant.colocated_variants {
                for (pop, value) in &cv.frequencies {
                    if max_af.is_none_or(|m| *value > m) {
                        max_af = Some(*value);
                        max_pop = pop.clone();
                    }
                }
            }
            max_pop
        }
        _ => tc
            .and_then(|x| x.plugin_data.get(field))
            .or_else(|| variant.plugin_data.get(field))
            .cloned()
            .unwrap_or_default(),
    }
}

/// The fields `add_colocated_variant_info` (OutputFactory.pm:1041-1167 at
/// release/116.2) and `add_colocated_frequency_data` (:1185-1278) set on a row:
/// the co-located variant columns and the fields of the frequency groups
/// (`%FREQUENCY_KEYS`, :107-113, and `--max_af`). The JSON factory stubs both
/// out (OutputFactory/JSON.pm:493-515) and writes the data under
/// `colocated_variants` instead; a plugin's field is never one of them.
pub fn colocated_field(field: &str) -> bool {
    const FREQUENCY_GROUPS: [&str; 5] = ["af", "af_1kg", "af_gnomade", "af_gnomadg", "max_af"];
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
        .any(|(group, fields)| FREQUENCY_GROUPS.contains(group) && fields.contains(&field))
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
/// input's reference allele is reverse-complemented first (:1049, 1100-1112).
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
    let mut ordered = colocated_in_output_order(&variant.colocated_variants);
    ordered.retain(|cv| cv.applies_to_allele(&allele));
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

/// The `CLIN_SIG` of a row over the variants that print on it, in order
/// (`colocated_column`).
fn clin_sig_column(
    variant: &InputVariant,
    allele: &str,
    ordered: &[&ColocatedVariant],
    options: &FieldOptions,
    list_sep: &str,
) -> String {
    let ref_allele = String::from_utf8_lossy(&variant.ref_allele);
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
mod tests {
    use super::*;
    use vep_core::variant::MatchedAllele;
    use vep_core::variation::ClinVarDetail;

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

    /// `--shift_length` lists SHIFT_LENGTH among the flag fields, and under
    /// `--shift_3prime 1` every transcript row carries its value, 0 included
    /// (OutputFactory.pm:1690-1695); an intergenic row carries none.
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

    /// The output-formats page states the subset: the cache's classification
    /// string as written, which Ensembl VEP builds from ClinVar's somatic
    /// records when it dumps the cache or reads a database and vep-rs never
    /// rebuilds.
    #[test]
    fn documented_subset_is_stated() {
        let page = std::fs::read_to_string(format!(
            "{}/../../docs/output-formats.md",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let paragraph = page
            .lines()
            .find(|l| {
                l.contains(
                    "`CLINVAR_SOMATIC_CLASSIFICATION` (under `--clinvar_somatic_classification`)",
                )
            })
            .expect("docs/output-formats.md describes CLINVAR_SOMATIC_CLASSIFICATION");
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
}
